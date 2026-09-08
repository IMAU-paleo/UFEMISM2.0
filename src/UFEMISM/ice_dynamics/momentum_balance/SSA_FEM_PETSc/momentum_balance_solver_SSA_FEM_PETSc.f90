module momentum_balance_solver_SSA_FEM_PETSc

#include <petsc/finclude/petscsys.h>

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA),
  ! discretised and solved entirely with PETSc: DMPlex + PetscFE for the discretisation
  ! (continuous P1 velocity on the mesh vertices) and PetscSNES for the linear solve
  ! inside a Picard viscosity iteration.
  !
  ! This is the finite-element counterpart of the existing finite-difference-based SSA
  ! solver (momentum_balance_solver_SSA). It is selected with
  !   choice_stress_balance_approximation = 'SSA_FEM_PETSc'
  ! and is added alongside the existing solvers without affecting any of them.
  !
  ! Implementation is staged (see SSA_PetscFE_SNES_implementation_plan.md in the
  ! repository root). Current state: Phase 2 - real, spatially varying coefficients.
  !
  ! Weak form (PETSc convention  residual = integral( f0 . phi + f1 : grad(phi) ) = 0 ),
  ! one 2-component vector field (u, v); N = eta*H, beta and the driving stress
  ! (tau_dx, tau_dy) are supplied per vertex through a PetscFE AUXILIARY field with
  ! 4 components a = [N, beta, tau_dx, tau_dy]:
  !
  !   f1[u,x] = 2 N (2 du/dx + dv/dy)      f1[u,y] = N (du/dy + dv/dx)
  !   f1[v,x] = N (du/dy + dv/dx)          f1[v,y] = 2 N (2 dv/dy + du/dx)
  !   f0[u]   = beta u + tau_dx            f0[v]   = beta v + tau_dy
  !
  ! Within one Picard iteration N is frozen (eta is evaluated from the strain rates
  ! of the previous velocity solution), so each SNES solve is linear and converges
  ! in one Newton step; the analytic Jacobian g0 = beta*I, g3 = d f1 / d grad(u) is
  ! exact for that frozen-N problem. The outer Picard loop (relaxation, velocity
  ! limiting, L2 stop criterion) mirrors momentum_balance_solver_SSA and reuses the
  ! same config knobs (visc_it_nit, visc_it_relax, visc_it_norm_dUV_tol, ...).
  !
  ! Not done yet: the basal-drag term makes the operator positive-definite under
  ! natural boundary conditions, so no essential BCs are imposed (Phase 5); the
  ! sub-grid grounded-fraction scaling of beta and the ice-front back-pressure are
  ! also Phase 5; a single-SNES nonlinear residual with eta = eta(grad u) evaluated
  ! pointwise is Phase 3.
  !
  ! The internal unknown is a nodal (vertex, P1) velocity field; the result is
  ! exposed on the triangles as u_vav_b / v_vav_b, exactly like the existing SSA
  ! solver, so that ice-thickness evolution and all other downstream code are
  ! unaffected.

  use precisions, only: dp
  use iso_c_binding, only: c_bool, c_double, c_int, c_intptr_t, c_ptr, c_funptr, c_funloc, &
    c_null_funptr, c_null_ptr, c_f_pointer
  use petsc, only: PetscErrorF, PETSC_COMM_SELF, PETSC_COMM_WORLD, PETSC_TRUE, PETSC_FALSE, &
    PETSC_NULL_DMLABEL, PETSC_NULL_VEC, tDM, tVec, tMat, tSNES, tKSP, tPC, tPetscObject, &
    tPetscFE, tPetscDS, tDMLabel, tPetscSection, &
    PetscFECreateLagrange, PetscFEDestroy, PetscObjectSetName, DMSetField, DMCreateDS, DMGetDS, &
    DMCreateMatrix, DMCreateGlobalVector, DMCreateLocalVector, &
    DMGlobalToLocalBegin, DMGlobalToLocalEnd, DMGetLocalSection, DMGetLabel, DMDestroy, &
    DMPlexGetDepthStratum, DMLabelGetValue, PetscSectionGetOffset, &
    VecSet, VecDestroy, VecGetArrayRead, VecRestoreArrayRead, VecSetValues, VecAssemblyBegin, &
    VecAssemblyEnd, MatDestroy, INSERT_VALUES, &
    SNESCreate, SNESSetDM, SNESSetType, SNESSetTolerances, SNESGetKSP, SNESSolve, SNESDestroy, &
    SNESGetIterationNumber, SNESNEWTONLS, KSPSetType, KSPGetPC, KSPPREONLY, PCSetType, PCLU
  use mpi_f08, only: MPI_ALLTOALL, MPI_ALLTOALLV, MPI_ALLREDUCE, MPI_COMM_WORLD, MPI_IN_PLACE, &
    MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_MAX, MPI_SUM, MPI_LOR, MPI_LOGICAL
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use parameters, only: grav, ice_density
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use bed_roughness_model_types, only: type_bed_roughness_model
  use reallocate_mod, only: reallocate_bounds
  use constitutive_equation, only: calc_ice_rheology_Glen, calc_effective_viscosity_Glen_2D
  use mesh_zeta, only: vertical_average
  use sliding_laws, only: calc_basal_friction_coefficient
  use petsc_dmplex, only: mesh_to_dmplex, dmplex_upsy_vertex_id_label_name
  use mesh_disc_apply_operators, only: map_a_b_2D, ddx_a_a_2D, ddy_a_a_2D

  implicit none

  private

  public :: type_momentum_balance_solver_SSA_FEM_PETSc

  ! Layout of the 4-component PetscFE auxiliary field
  integer, parameter :: i_N     = 1   ! N = eta * H          [Pa yr m]
  integer, parameter :: i_beta  = 2   ! basal friction coeff [Pa yr m^-1]
  integer, parameter :: i_taudx = 3   ! driving stress, x    [Pa]
  integer, parameter :: i_taudy = 4   ! driving stress, y    [Pa]
  integer, parameter :: n_aux_comp = 4

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SSA_FEM_PETSc

    ! Persistent PETSc objects (rebuilt on remap)
    type(tDM)      :: dm            ! primary DM: the P1 velocity field
    type(tPetscFE) :: fe
    type(tSNES)    :: snes
    type(tMat)     :: jac
    type(tDM)      :: dm_aux        ! clone of dm carrying the auxiliary field
    type(tPetscFE) :: fe_aux
    type(tVec)     :: aux_vec       ! local vector of dm_aux: [N, beta, tau_dx, tau_dy] per vertex
    logical        :: petsc_is_built = .false.

    ! Solution
    real(dp), dimension(:), allocatable :: u_vav_a, v_vav_a   ! [m yr^-1] nodal (vertices), vi1:vi2
    real(dp), dimension(:), allocatable :: u_vav_b, v_vav_b   ! [m yr^-1] exposed result (triangles), ti1:ti2

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_SSA_FEM_PETSc_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_SSA_FEM_PETSc_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_SSA_FEM_PETSc_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_SSA_FEM_PETSc_run
      procedure, public :: set_velocities_to_solver_results   => momentum_balance_solver_SSA_FEM_PETSc_set_velocities
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_SSA_FEM_PETSc_remap

      procedure, public :: get_momentum_balance_solver_name
      procedure, public :: create_restart_file_old            => create_restart_file_old_SSA_FEM_PETSc
      procedure, public :: write_to_restart_file_old          => write_to_restart_file_old_SSA_FEM_PETSc

      procedure, private :: build_petsc_objects
      procedure, private :: destroy_petsc_objects
      procedure, private :: solve_linearised_SSA
      procedure, private :: calc_auxiliary_fields

  end type type_momentum_balance_solver_SSA_FEM_PETSc

  ! Direct C bindings for PETSc routines whose Fortran wrappers are absent from the
  ! installed PETSc package (same work-around as ct_PETSc_SNES_Poisson.f90).
  interface

    subroutine DMPlexSetSNESLocalFEM( dm, has_boundary, ctx, ierr)
      import :: c_bool, c_intptr_t, tDM
      type(tDM),            intent(inout) :: dm
      logical(kind=c_bool), intent(in)    :: has_boundary
      integer(c_intptr_t),  intent(in)    :: ctx
      integer,              intent(out)   :: ierr
    end subroutine DMPlexSetSNESLocalFEM

    integer(c_int) function dm_clone( dm, newdm) bind(C, name='DMClone')
      import :: c_int, c_intptr_t
      integer(c_intptr_t), value       :: dm
      integer(c_intptr_t), intent(out) :: newdm
    end function dm_clone

    integer(c_int) function dm_set_auxiliary_vec( dm, label, value, part, aux) bind(C, name='DMSetAuxiliaryVec')
      import :: c_int, c_intptr_t
      integer(c_intptr_t), value :: dm, label, aux
      integer(c_int),      value :: value, part   ! PetscInt (32-bit in this build)
    end function dm_set_auxiliary_vec

    integer(c_int) function petsc_ds_set_residual( ds, field, f0, f1) bind(C, name='PetscDSSetResidual')
      import :: c_funptr, c_int, c_intptr_t
      integer(c_intptr_t), value :: ds
      integer(c_intptr_t), value :: field
      type(c_funptr),      value :: f0, f1
    end function petsc_ds_set_residual

    integer(c_int) function petsc_ds_set_jacobian( ds, field_test, field_trial, g0, g1, g2, g3) &
      bind(C, name='PetscDSSetJacobian')
      import :: c_funptr, c_int, c_intptr_t
      integer(c_intptr_t), value :: ds
      integer(c_intptr_t), value :: field_test, field_trial
      type(c_funptr),      value :: g0, g1, g2, g3
    end function petsc_ds_set_jacobian

    integer(c_int) function snes_set_jacobian( snes, jacobian, preconditioner, function, ctx) &
      bind(C, name='SNESSetJacobian')
      import :: c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: snes, jacobian, preconditioner
      type(c_funptr),      value :: function
      type(c_ptr),         value :: ctx
    end function snes_set_jacobian

    integer(c_int) function snes_get_converged_reason( snes, reason) bind(C, name='SNESGetConvergedReason')
      import :: c_int, c_intptr_t
      integer(c_intptr_t), value :: snes
      integer(c_int),      intent(out) :: reason
    end function snes_get_converged_reason

  end interface

contains

  ! ===== Main routines =====

  subroutine momentum_balance_solver_SSA_FEM_PETSc_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( self%u_vav_a( self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%v_vav_a( self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%u_vav_b( self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%v_vav_b( self%mesh%ti1:self%mesh%ti2), source = 0._dp)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_allocate

  subroutine momentum_balance_solver_SSA_FEM_PETSc_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (self%petsc_is_built) call self%destroy_petsc_objects()

    if (allocated( self%u_vav_a)) deallocate( self%u_vav_a)
    if (allocated( self%v_vav_a)) deallocate( self%v_vav_a)
    if (allocated( self%u_vav_b)) deallocate( self%u_vav_b)
    if (allocated( self%v_vav_b)) deallocate( self%v_vav_b)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_deallocate

  subroutine momentum_balance_solver_SSA_FEM_PETSc_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,'(A)') '    NOTE: the SSA_FEM_PETSc solver is at Phase 2 - real coefficients, ' // &
      'Picard viscosity iteration, natural boundary conditions only.'

    call self%build_petsc_objects()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_initialise

  subroutine momentum_balance_solver_SSA_FEM_PETSc_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the SSA with a Picard viscosity iteration
    !< around a linear PetscFE / PetscSNES solve.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    class(atype_ice_model_data),                       intent(inout) :: ice
    class(atype_ice_geometry_model_data),              intent(in   ) :: geom
    type(type_bed_roughness_model),                    intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,               intent(in   ) :: BC_prescr_mask_b
    real(dp), dimension(:  ), optional,               intent(in   ) :: BC_prescr_u_b
    real(dp), dimension(:  ), optional,               intent(in   ) :: BC_prescr_v_b
    integer,  dimension(:,:), optional,               intent(in   ) :: BC_prescr_mask_bk
    real(dp), dimension(:,:), optional,               intent(in   ) :: BC_prescr_u_bk
    real(dp), dimension(:,:), optional,               intent(in   ) :: BC_prescr_v_bk

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_run'
    integer                                          :: ierr, it, snes_its
    logical                                          :: grounded_ice_exists, has_converged
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: N_a, beta_a, taudx_a, taudy_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: A_flow_vav_a, u_prev, v_prev
    real(dp), dimension(self%mesh%nz)                :: A_prof
    real(dp)                                         :: L2_uv, uabs, umin, umax
    real(dp)                                         :: res1, res2
    integer                                          :: vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! If there is no grounded ice or no sliding, there is nothing to solve
    grounded_ice_exists = any( geom%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists .or. C%choice_sliding_law == 'no_sliding') then
      self%u_vav_a = 0._dp
      self%v_vav_a = 0._dp
      self%u_vav_b = 0._dp
      self%v_vav_b = 0._dp
      self%n_visc_its = 0
      self%n_Axb_its  = 0
      call finalise_routine( routine_name)
      return
    end if

    ! Vertically averaged flow factor A - velocity-independent, computed once
    call calc_ice_rheology_Glen( self%mesh, ice, geom)
    do vi = self%mesh%vi1, self%mesh%vi2
      A_prof = ice%A_flow( vi,:)
      A_flow_vav_a( vi) = vertical_average( self%mesh%zeta, A_prof)
    end do

    ! The Picard viscosity iteration
    self%n_visc_its = 0
    self%n_Axb_its  = 0
    has_converged   = .false.
    do it = 1, C%visc_it_nit

      u_prev = self%u_vav_a
      v_prev = self%v_vav_a

      ! Freeze the coefficients at the current velocity solution
      call self%calc_auxiliary_fields( ice, geom, bed_roughness, A_flow_vav_a, N_a, beta_a, taudx_a, taudy_a)
      call fill_PETSc_aux_from_mesh_vertices( self%dm, self%dm_aux, self%aux_vec, self%mesh, &
        N_a, beta_a, taudx_a, taudy_a)
      ierr = dm_set_auxiliary_vec( self%dm%v, 0_c_intptr_t, 0_c_int, 0_c_int, self%aux_vec%v)
      CHKERRQ( ierr)

      ! One linear SNES solve for the frozen-coefficient SSA
      call self%solve_linearised_SSA( snes_its)
      self%n_visc_its = it
      self%n_Axb_its  = self%n_Axb_its + snes_its

      if (any( isnan( self%u_vav_a)) .or. any( isnan( self%v_vav_a))) &
        call crash('SSA_FEM_PETSc: NaN in the velocity solution')

      ! Relax and limit for stability
      self%u_vav_a = C%visc_it_relax * self%u_vav_a + (1._dp - C%visc_it_relax) * u_prev
      self%v_vav_a = C%visc_it_relax * self%v_vav_a + (1._dp - C%visc_it_relax) * v_prev
      do vi = self%mesh%vi1, self%mesh%vi2
        uabs = sqrt( self%u_vav_a( vi)**2 + self%v_vav_a( vi)**2)
        if (uabs > C%vel_max) then
          self%u_vav_a( vi) = self%u_vav_a( vi) * C%vel_max / uabs
          self%v_vav_a( vi) = self%v_vav_a( vi) * C%vel_max / uabs
        end if
      end do

      ! L2-norm of the change between successive velocity solutions (as in momentum_balance_solver_SSA)
      res1 = sum( (self%u_vav_a - u_prev)**2 + (self%v_vav_a - v_prev)**2)
      res2 = sum( (self%u_vav_a + u_prev)**2 + (self%v_vav_a + v_prev)**2)
      call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      L2_uv = 2._dp * res1 / max( res2, 1e-8_dp)

      umin = minval( sqrt( self%u_vav_a**2 + self%v_vav_a**2))
      umax = maxval( sqrt( self%u_vav_a**2 + self%v_vav_a**2))
      call MPI_ALLREDUCE( MPI_IN_PLACE, umin, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, umax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (par%primary) write(0,'(A,I3,A,ES10.3,A,ES10.3)') &
        '    SSA_FEM_PETSc visc. iter. ', it, ': L2 = ', L2_uv, ', max speed = ', umax

      if (L2_uv < C%visc_it_norm_dUV_tol) then
        has_converged = .true.
        exit
      end if

    end do

    if (.not. has_converged .and. par%primary) &
      call warning('SSA_FEM_PETSc: viscosity iteration did not converge within {int_01} iterations', &
        int_01 = C%visc_it_nit)

    ! Expose the result on the triangles
    call map_a_b_2D( self%mesh, self%u_vav_a, self%u_vav_b)
    call map_a_b_2D( self%mesh, self%v_vav_a, self%v_vav_b)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_run

  subroutine momentum_balance_solver_SSA_FEM_PETSc_set_velocities( self, ice, vel)
    !< Hand the solver result to the ice velocity model.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in   ) :: self
    class(atype_ice_model_data),                       intent(inout) :: ice
    class(atype_ice_velocity_model),                   intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_set_velocities'
    integer                     :: ti, vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Velocities on the triangles
    do ti = self%mesh%ti1, self%mesh%ti2
      vel%u_3D_b( ti,:) = self%u_vav_b( ti)
      vel%v_3D_b( ti,:) = self%v_vav_b( ti)
    end do

    ! Strain rates on the vertices - still zero until Phase 3 computes them from the FE field
    do vi = self%mesh%vi1, self%mesh%vi2
      vel%du_dx_3D( vi,:) = 0._dp
      vel%du_dy_3D( vi,:) = 0._dp
      vel%dv_dx_3D( vi,:) = 0._dp
      vel%dv_dy_3D( vi,:) = 0._dp
    end do

    ! In the SSA, vertical gradients of u,v, and all gradients of w, are neglected
    vel%du_dz_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dv_dz_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dw_dx_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_set_velocities

  subroutine momentum_balance_solver_SSA_FEM_PETSc_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    type(type_mesh),                                   intent(in   ) :: mesh_old
    type(type_mesh), target,                           intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Phase 2: rebuild everything from scratch on the new mesh (self%mesh has already
    ! been repointed to mesh_new by remap_model). Velocities are reset to zero; a
    ! proper remap of u_vav via the a-grid follows in a later phase.
    if (self%petsc_is_built) call self%destroy_petsc_objects()

    call reallocate_bounds( self%u_vav_a, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%v_vav_a, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%u_vav_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%v_vav_b, mesh_new%ti1, mesh_new%ti2)
    self%u_vav_a = 0._dp; self%v_vav_a = 0._dp
    self%u_vav_b = 0._dp; self%v_vav_b = 0._dp

    call self%build_petsc_objects()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'SSA_FEM_PETSc'
  end function get_momentum_balance_solver_name

  ! ===== PETSc object lifecycle =====

  subroutine build_petsc_objects( self)
    !< Build the DMPlex, the P1 vector PetscFE field, the auxiliary-field DM, the
    !< PetscDS weak form and the SNES.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_petsc_objects'
    type(tPetscDS)              :: ds
    type(tPetscObject)          :: fe_object, fe_aux_object
    type(tKSP)                  :: ksp
    type(tPC)                   :: pc
    integer                     :: ierr
    integer(c_intptr_t)         :: no_context

    call init_routine( routine_name)

    no_context = 0_c_intptr_t

    ! DMPlex from the UFEMISM mesh (creates + distributes; preserves UPSY vertex IDs
    ! in the 'upsy_vertex_id' DMLabel)
    call mesh_to_dmplex( self%mesh, self%dm)

    ! Primary field: one 2-component P1 Lagrange velocity
    PetscCall( PetscFECreateLagrange( PETSC_COMM_SELF, 2, 2, PETSC_TRUE, 1, -1, self%fe, ierr))
    PetscCall( PetscObjectSetName( self%fe, 'velocity', ierr))
    PetscObjectSpecificCast( fe_object, self%fe)
    PetscCall( DMSetField( self%dm, 0, PETSC_NULL_DMLABEL, fe_object, ierr))
    PetscCall( DMCreateDS( self%dm, ierr))
    PetscCall( DMGetDS( self%dm, ds, ierr))

    ! Auxiliary field: [N, beta, tau_dx, tau_dy], P1, on a clone of the primary DM
    ierr = dm_clone( self%dm%v, self%dm_aux%v)
    CHKERRQ( ierr)
    PetscCall( PetscFECreateLagrange( PETSC_COMM_SELF, 2, n_aux_comp, PETSC_TRUE, 1, -1, self%fe_aux, ierr))
    PetscCall( PetscObjectSetName( self%fe_aux, 'SSA_coefficients', ierr))
    PetscObjectSpecificCast( fe_aux_object, self%fe_aux)
    PetscCall( DMSetField( self%dm_aux, 0, PETSC_NULL_DMLABEL, fe_aux_object, ierr))
    PetscCall( DMCreateDS( self%dm_aux, ierr))
    PetscCall( DMCreateLocalVector( self%dm_aux, self%aux_vec, ierr))
    PetscCall( VecSet( self%aux_vec, 0._dp, ierr))
    ierr = dm_set_auxiliary_vec( self%dm%v, 0_c_intptr_t, 0_c_int, 0_c_int, self%aux_vec%v)
    CHKERRQ( ierr)

    ! Weak form: residual (f0, f1) and analytic Jacobian (g0, g3), all reading the aux field
    ierr = petsc_ds_set_residual( ds%v, 0_c_intptr_t, c_funloc( SSA_FEM_PETSc_f0), c_funloc( SSA_FEM_PETSc_f1))
    CHKERRQ( ierr)
    ierr = petsc_ds_set_jacobian( ds%v, 0_c_intptr_t, 0_c_intptr_t, &
      c_funloc( SSA_FEM_PETSc_g0), c_null_funptr, c_null_funptr, c_funloc( SSA_FEM_PETSc_g3))
    CHKERRQ( ierr)

    ! SNES; DMPlex assembles the residual and Jacobian from the PetscDS weak form
    PetscCall( SNESCreate( PETSC_COMM_WORLD, self%snes, ierr))
    PetscCall( SNESSetDM( self%snes, self%dm, ierr))
    PetscCall( DMPlexSetSNESLocalFEM( self%dm, PETSC_FALSE, no_context, ierr))
    PetscCall( DMCreateMatrix( self%dm, self%jac, ierr))
    ierr = snes_set_jacobian( self%snes%v, self%jac%v, self%jac%v, c_null_funptr, c_null_ptr)
    CHKERRQ( ierr)

    ! Solver configuration: Newton line search, direct linear solve (LU -> MUMPS on
    ! more than one rank, as in ct_PETSc_SNES_Poisson). Tuning is Phase 6.
    PetscCall( SNESSetType( self%snes, SNESNEWTONLS, ierr))
    call SNESSetTolerances( self%snes, self%PETSc_abstol, self%PETSc_rtol, 1.0e-12_dp, 20, 1000, ierr)
    CHKERRQ( ierr)
    PetscCall( SNESGetKSP( self%snes, ksp, ierr))
    PetscCall( KSPSetType( ksp, KSPPREONLY, ierr))
    PetscCall( KSPGetPC( ksp, pc, ierr))
    PetscCall( PCSetType( pc, PCLU, ierr))

    self%petsc_is_built = .true.

    call finalise_routine( routine_name)

  end subroutine build_petsc_objects

  subroutine destroy_petsc_objects( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'destroy_petsc_objects'
    integer                     :: ierr

    call init_routine( routine_name)

    PetscCall( VecDestroy( self%aux_vec, ierr))
    PetscCall( PetscFEDestroy( self%fe_aux, ierr))
    PetscCall( DMDestroy( self%dm_aux, ierr))
    PetscCall( MatDestroy( self%jac, ierr))
    PetscCall( SNESDestroy( self%snes, ierr))
    PetscCall( PetscFEDestroy( self%fe, ierr))
    PetscCall( DMDestroy( self%dm, ierr))
    self%petsc_is_built = .false.

    call finalise_routine( routine_name)

  end subroutine destroy_petsc_objects

  subroutine solve_linearised_SSA( self, snes_its)
    !< One linear SNES solve for the SSA with the currently attached (frozen) coefficients.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    integer,                                           intent(  out) :: snes_its

    ! Local variables:
    character(len=*), parameter :: routine_name = 'solve_linearised_SSA'
    type(tVec)                  :: solution
    integer                     :: ierr
    integer(c_int)              :: snes_reason

    call init_routine( routine_name)

    PetscCall( DMCreateGlobalVector( self%dm, solution, ierr))
    PetscCall( VecSet( solution, 0._dp, ierr))
    PetscCall( SNESSolve( self%snes, PETSC_NULL_VEC, solution, ierr))
    PetscCall( SNESGetIterationNumber( self%snes, snes_its, ierr))
    ierr = snes_get_converged_reason( self%snes%v, snes_reason)
    CHKERRQ( ierr)
    if (snes_reason < 0) &
      call crash('SSA_FEM_PETSc: SNES diverged (SNESConvergedReason = {int_01})', int_01 = int( snes_reason))

    call copy_PETSc_solution_to_mesh_vertices_vec2( self%dm, solution, self%mesh, self%u_vav_a, self%v_vav_a)
    PetscCall( VecDestroy( solution, ierr))

    call finalise_routine( routine_name)

  end subroutine solve_linearised_SSA

  ! ===== Coefficient (auxiliary-field) calculation =====

  subroutine calc_auxiliary_fields( self, ice, geom, bed_roughness, A_flow_vav_a, N_a, beta_a, taudx_a, taudy_a)
    !< Compute the frozen SSA coefficients on the mesh vertices from the current
    !< velocity solution: N = eta*H, the basal friction coefficient beta, and the
    !< driving stress (tau_dx, tau_dy).

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in   ) :: self
    class(atype_ice_model_data),                       intent(inout) :: ice
    class(atype_ice_geometry_model_data),              intent(in   ) :: geom
    type(type_bed_roughness_model),                    intent(in   ) :: bed_roughness
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2),  intent(in   ) :: A_flow_vav_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2),  intent(  out) :: N_a, beta_a, taudx_a, taudy_a

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_auxiliary_fields'
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: du_dx_a, du_dy_a, dv_dx_a, dv_dy_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: dHs_dx_a, dHs_dy_a
    real(dp)                                         :: n_glen, eps0, A_min, eta_max, eta
    integer                                         :: vi

    call init_routine( routine_name)

    n_glen = C%Glens_flow_law_exponent
    eps0   = C%Glens_flow_law_epsilon_sq_0

    ! Maximum allowed effective viscosity, for stability (as in momentum_balance_solver_SSA)
    A_min   = 1e-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / n_glen) * eps0**((1._dp - n_glen) / (2._dp * n_glen))

    ! Effective viscosity from Glen's flow law and the strain rates of the current solution
    call ddx_a_a_2D( self%mesh, self%u_vav_a, du_dx_a)
    call ddy_a_a_2D( self%mesh, self%u_vav_a, du_dy_a)
    call ddx_a_a_2D( self%mesh, self%v_vav_a, dv_dx_a)
    call ddy_a_a_2D( self%mesh, self%v_vav_a, dv_dy_a)

    do vi = self%mesh%vi1, self%mesh%vi2
      eta = calc_effective_viscosity_Glen_2D( eps0, du_dx_a( vi), du_dy_a( vi), dv_dx_a( vi), dv_dy_a( vi), &
        A_flow_vav_a( vi))
      eta = min( max( eta, C%visc_eff_min), eta_max)
      N_a( vi) = eta * max( 0.1_dp, geom%Hi( vi))
    end do

    ! Driving stress on the vertices: tau_d = -rho g H grad(Hs)
    call ddx_a_a_2D( self%mesh, geom%Hs, dHs_dx_a)
    call ddy_a_a_2D( self%mesh, geom%Hs, dHs_dy_a)
    do vi = self%mesh%vi1, self%mesh%vi2
      taudx_a( vi) = -ice_density * grav * geom%Hi( vi) * dHs_dx_a( vi)
      taudy_a( vi) = -ice_density * grav * geom%Hi( vi) * dHs_dy_a( vi)
    end do

    ! Basal friction coefficient from the sliding law, evaluated at the current velocity,
    ! scaled by the sub-grid grounded fraction so that friction vanishes under floating ice.
    ! momentum_balance_solver_SSA does this on the b-grid (fraction_gr_b) in
    ! calc_applied_basal_friction_coefficient; here everything is on the a-grid, so we use
    ! the vertex grounded fraction geom%fraction_gr with the same exponent.
    call calc_basal_friction_coefficient( self%mesh, geom, bed_roughness, self%u_vav_a, self%v_vav_a, &
      ice%effective_pressure, ice%till_yield_stress, ice%basal_friction_coefficient)
    do vi = self%mesh%vi1, self%mesh%vi2
      beta_a( vi) = ice%basal_friction_coefficient( vi)
      if (C%do_GL_subgrid_friction) then
        beta_a( vi) = beta_a( vi) * geom%fraction_gr( vi)**C%subgrid_friction_exponent_on_B_grid
      end if
    end do

    call finalise_routine( routine_name)

  end subroutine calc_auxiliary_fields

  ! ===== Solution / coefficient transfer between UFEMISM vertex arrays and PETSc =====

  subroutine copy_PETSc_solution_to_mesh_vertices_vec2( dm, solution, mesh, u_a, v_a)
    !< Scatter a 2-component nodal PETSc solution back onto the UFEMISM vertex
    !< distribution (vi1:vi2), using the 'upsy_vertex_id' DMLabel and
    !< mesh%V_owning_process.

    ! In/output variables:
    type(tDM),                              intent(in   ) :: dm
    type(tVec),                             intent(in   ) :: solution
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: u_a, v_a

    ! Local variables:
    type(tVec)                         :: local_solution
    type(tPetscSection)                :: local_section
    type(tDMLabel)                     :: upsy_vertex_id_label
    real(dp), dimension(:), pointer    :: vals
    integer, dimension(:), allocatable :: send_counts, recv_counts, send_displ, recv_displ, send_pos
    integer, dimension(:), allocatable :: send_vi, recv_vi
    real(dp), dimension(:), allocatable :: send_u, send_v, recv_u, recv_v
    integer, dimension(mesh%vi1:mesh%vi2) :: ncopies
    integer                            :: ierr, point, vstart, vend, vi, dest, off, ip, si, ns, nr, k

    allocate( send_counts( 0:par%n-1), recv_counts( 0:par%n-1), source = 0)
    allocate( send_displ ( 0:par%n-1), recv_displ ( 0:par%n-1), send_pos( 0:par%n-1))

    PetscCall( DMCreateLocalVector( dm, local_solution, ierr))
    PetscCall( DMGlobalToLocalBegin( dm, solution, INSERT_VALUES, local_solution, ierr))
    PetscCall( DMGlobalToLocalEnd(   dm, solution, INSERT_VALUES, local_solution, ierr))
    PetscCall( DMGetLocalSection( dm, local_section, ierr))
    PetscCall( DMGetLabel( dm, dmplex_upsy_vertex_id_label_name, upsy_vertex_id_label, ierr))
    PetscCall( DMPlexGetDepthStratum( dm, 0, vstart, vend, ierr))

    do point = vstart, vend - 1
      PetscCall( DMLabelGetValue( upsy_vertex_id_label, point, vi, ierr))
      if (vi < 1 .or. vi > mesh%nV) call crash('DMPlex vertex lacks a valid UPSY vertex ID')
      dest = mesh%V_owning_process( vi)
      if (dest < 0 .or. dest >= par%n) call crash('UPSY vertex has an invalid owning process')
      send_counts( dest) = send_counts( dest) + 1
    end do

    send_displ( 0) = 0
    do ip = 1, par%n-1
      send_displ( ip) = send_displ( ip-1) + send_counts( ip-1)
    end do
    ns = sum( send_counts)
    send_pos = send_displ
    allocate( send_vi( max( 1, ns)), send_u( max( 1, ns)), send_v( max( 1, ns)))

    PetscCall( VecGetArrayRead( local_solution, vals, ierr))
    do point = vstart, vend - 1
      PetscCall( DMLabelGetValue( upsy_vertex_id_label, point, vi, ierr))
      dest = mesh%V_owning_process( vi)
      PetscCall( PetscSectionGetOffset( local_section, point, off, ierr))
      si = send_pos( dest) + 1
      send_vi( si) = vi
      send_u ( si) = vals( off + 1)
      send_v ( si) = vals( off + 2)
      send_pos( dest) = send_pos( dest) + 1
    end do
    PetscCall( VecRestoreArrayRead( local_solution, vals, ierr))
    PetscCall( VecDestroy( local_solution, ierr))

    call MPI_ALLTOALL( send_counts, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    recv_displ( 0) = 0
    do ip = 1, par%n-1
      recv_displ( ip) = recv_displ( ip-1) + recv_counts( ip-1)
    end do
    nr = sum( recv_counts)
    allocate( recv_vi( max( 1, nr)), recv_u( max( 1, nr)), recv_v( max( 1, nr)))

    call MPI_ALLTOALLV( send_vi, send_counts, send_displ, MPI_INTEGER, &
      recv_vi, recv_counts, recv_displ, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLTOALLV( send_u, send_counts, send_displ, MPI_DOUBLE_PRECISION, &
      recv_u, recv_counts, recv_displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    call MPI_ALLTOALLV( send_v, send_counts, send_displ, MPI_DOUBLE_PRECISION, &
      recv_v, recv_counts, recv_displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    u_a = 0._dp
    v_a = 0._dp
    ncopies = 0
    do k = 1, nr
      vi = recv_vi( k)
      if (vi < mesh%vi1 .or. vi > mesh%vi2) call crash('DMPlex solution was sent to the wrong UPSY process')
      u_a( vi) = u_a( vi) + recv_u( k)
      v_a( vi) = v_a( vi) + recv_v( k)
      ncopies( vi) = ncopies( vi) + 1
    end do
    if (any( ncopies == 0)) call crash('DMPlex distribution omitted an UPSY vertex')
    u_a = u_a / real( ncopies, dp)
    v_a = v_a / real( ncopies, dp)

  end subroutine copy_PETSc_solution_to_mesh_vertices_vec2

  subroutine fill_PETSc_aux_from_mesh_vertices( dm_topo, dm_aux, aux_vec, mesh, c1, c2, c3, c4)
    !< Scatter four UFEMISM vertex arrays (vi1:vi2) into the 4-component local
    !< auxiliary vector of dm_aux. Inverse of copy_PETSc_solution_to_mesh_vertices_vec2:
    !< each rank requests, for its local DMPlex vertices, the coefficient values from
    !< the UFEMISM process that owns that vertex.

    ! In/output variables:
    type(tDM),                              intent(in   ) :: dm_topo   ! for topology + 'upsy_vertex_id' label
    type(tDM),                              intent(in   ) :: dm_aux    ! for the 4-component local section
    type(tVec),                             intent(inout) :: aux_vec
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: c1, c2, c3, c4

    ! Local variables:
    type(tPetscSection)                 :: aux_section
    type(tDMLabel)                      :: upsy_vertex_id_label
    integer, dimension(:), allocatable  :: send_counts, recv_counts, send_displ, recv_displ, send_pos
    integer, dimension(:), allocatable  :: send_counts4, recv_counts4, send_displ4, recv_displ4
    integer, dimension(:), allocatable  :: req_vi, recv_req_vi
    integer, dimension(:), allocatable  :: local_pt, local_slot
    real(dp), dimension(:), allocatable :: reply_vals, recv_reply
    integer, dimension(4)               :: idx4
    real(dp), dimension(4)              :: vals4
    integer :: ierr, point, vstart, vend, vi, dest, ip, si, ns, nr, k, nlv, off

    allocate( send_counts( 0:par%n-1), recv_counts( 0:par%n-1), source = 0)
    allocate( send_displ ( 0:par%n-1), recv_displ ( 0:par%n-1), send_pos( 0:par%n-1))

    PetscCall( DMGetLabel( dm_topo, dmplex_upsy_vertex_id_label_name, upsy_vertex_id_label, ierr))
    PetscCall( DMPlexGetDepthStratum( dm_topo, 0, vstart, vend, ierr))
    PetscCall( DMGetLocalSection( dm_aux, aux_section, ierr))
    nlv = vend - vstart
    allocate( local_pt( max( 1, nlv)), local_slot( max( 1, nlv)))

    ! Pass 1: count local DMPlex vertices per owning UFEMISM process
    k = 0
    do point = vstart, vend - 1
      k = k + 1
      PetscCall( DMLabelGetValue( upsy_vertex_id_label, point, vi, ierr))
      if (vi < 1 .or. vi > mesh%nV) call crash('DMPlex vertex lacks a valid UPSY vertex ID')
      local_pt( k) = point
      dest = mesh%V_owning_process( vi)
      if (dest < 0 .or. dest >= par%n) call crash('UPSY vertex has an invalid owning process')
      send_counts( dest) = send_counts( dest) + 1
    end do

    send_displ( 0) = 0
    do ip = 1, par%n-1
      send_displ( ip) = send_displ( ip-1) + send_counts( ip-1)
    end do
    ns = sum( send_counts)
    send_pos = send_displ
    allocate( req_vi( max( 1, ns)))

    ! Pass 2: pack the requested vertex IDs per owner, remember each vertex' slot
    k = 0
    do point = vstart, vend - 1
      k = k + 1
      PetscCall( DMLabelGetValue( upsy_vertex_id_label, point, vi, ierr))
      dest = mesh%V_owning_process( vi)
      si = send_pos( dest)
      req_vi( si + 1) = vi
      local_slot( k) = si
      send_pos( dest) = send_pos( dest) + 1
    end do

    ! Exchange the request lists
    call MPI_ALLTOALL( send_counts, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    recv_displ( 0) = 0
    do ip = 1, par%n-1
      recv_displ( ip) = recv_displ( ip-1) + recv_counts( ip-1)
    end do
    nr = sum( recv_counts)
    allocate( recv_req_vi( max( 1, nr)), reply_vals( max( 1, 4*nr)))

    call MPI_ALLTOALLV( req_vi, send_counts, send_displ, MPI_INTEGER, &
      recv_req_vi, recv_counts, recv_displ, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Fill the replies with the four coefficients for each requested (locally owned) vertex
    do k = 1, nr
      vi = recv_req_vi( k)
      if (vi < mesh%vi1 .or. vi > mesh%vi2) call crash('aux request sent to the wrong UPSY process')
      reply_vals( 4*k-3) = c1( vi)
      reply_vals( 4*k-2) = c2( vi)
      reply_vals( 4*k-1) = c3( vi)
      reply_vals( 4*k  ) = c4( vi)
    end do

    ! Send the replies back (4 doubles per requested vertex)
    allocate( recv_reply( max( 1, 4*ns)))
    allocate( send_counts4( 0:par%n-1), recv_counts4( 0:par%n-1))
    allocate( send_displ4 ( 0:par%n-1), recv_displ4 ( 0:par%n-1))
    send_counts4 = 4 * send_counts
    recv_counts4 = 4 * recv_counts
    send_displ4  = 4 * send_displ
    recv_displ4  = 4 * recv_displ
    call MPI_ALLTOALLV( reply_vals, recv_counts4, recv_displ4, MPI_DOUBLE_PRECISION, &
      recv_reply, send_counts4, send_displ4, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    ! Write into the local auxiliary vector at each vertex' section offset
    do k = 1, nlv
      point = local_pt( k)
      si    = local_slot( k)
      PetscCall( PetscSectionGetOffset( aux_section, point, off, ierr))
      idx4  = [off, off+1, off+2, off+3]
      vals4 = recv_reply( 4*si+1 : 4*si+4)
      PetscCall( VecSetValues( aux_vec, 4, idx4, vals4, INSERT_VALUES, ierr))
    end do
    PetscCall( VecAssemblyBegin( aux_vec, ierr))
    PetscCall( VecAssemblyEnd( aux_vec, ierr))

  end subroutine fill_PETSc_aux_from_mesh_vertices

  ! ===== PetscDS pointwise weak-form functions (bind(C)) =====
  !
  ! PETSc assembles  residual = integral( f0 . phi + f1 : grad(phi) ).
  ! Field 0 is the 2-vector velocity (u, v); dim = 2, Nc = 2.
  ! One auxiliary field with 4 components: a(1..4) = [N, beta, tau_dx, tau_dy].
  ! Gradient layout u_x[c*dim + d]:  u_x(1)=du/dx u_x(2)=du/dy u_x(3)=dv/dx u_x(4)=dv/dy

  subroutine SSA_FEM_PETSc_f0( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f0) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f0
    real(c_double), value :: time
    real(c_double), pointer :: u_values(:), a_values(:), f0_values(:)

    call c_f_pointer( u, u_values, [2])
    call c_f_pointer( a, a_values, [n_aux_comp])
    call c_f_pointer( f0, f0_values, [2])

    ! f0 = beta * u + tau_d   (basal drag minus the RHS driving stress, moved to the LHS)
    f0_values( 1) = a_values( i_beta) * u_values( 1) + a_values( i_taudx)
    f0_values( 2) = a_values( i_beta) * u_values( 2) + a_values( i_taudy)

  end subroutine SSA_FEM_PETSc_f0

  subroutine SSA_FEM_PETSc_f1( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f1) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f1
    real(c_double), value :: time
    real(c_double), pointer :: u_x_values(:), a_values(:), f1_values(:)
    real(c_double)          :: N, du_dx, du_dy, dv_dx, dv_dy

    call c_f_pointer( u_x, u_x_values, [4])
    call c_f_pointer( a, a_values, [n_aux_comp])
    call c_f_pointer( f1, f1_values, [4])

    N     = a_values( i_N)
    du_dx = u_x_values( 1)
    du_dy = u_x_values( 2)
    dv_dx = u_x_values( 3)
    dv_dy = u_x_values( 4)

    f1_values( 1) = 2._c_double * N * (2._c_double * du_dx + dv_dy)   ! [u, d/dx]
    f1_values( 2) =               N * (du_dy + dv_dx)                 ! [u, d/dy]
    f1_values( 3) =               N * (du_dy + dv_dx)                 ! [v, d/dx]
    f1_values( 4) = 2._c_double * N * (2._c_double * dv_dy + du_dx)   ! [v, d/dy]

  end subroutine SSA_FEM_PETSc_f1

  subroutine SSA_FEM_PETSc_g0( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, u_tshift, x, nconstants, constants, g0) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants
    real(c_double), value :: time, u_tshift
    real(c_double), intent(out) :: g0(*)
    real(c_double), pointer :: a_values(:)

    call c_f_pointer( a, a_values, [n_aux_comp])

    ! d f0_c / d u_c'  =  beta * delta_{c c'}   (2x2, row-major)
    g0( 1:4) = 0._c_double
    g0( 1)   = a_values( i_beta)
    g0( 4)   = a_values( i_beta)

  end subroutine SSA_FEM_PETSc_g0

  subroutine SSA_FEM_PETSc_g3( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, u_tshift, x, nconstants, constants, g3) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants
    real(c_double), value :: time, u_tshift
    real(c_double), intent(out) :: g3(*)
    real(c_double), pointer :: a_values(:)
    real(c_double)          :: N

    call c_f_pointer( a, a_values, [n_aux_comp])
    N = a_values( i_N)

    ! g3[c,c',d,d'] = d f1[c,d] / d(du_c'/dx_d'),  stored at
    ! index0 = ((c*Nc + c')*dim + d)*dim + d'   with Nc = dim = 2  (Fortran = index0 + 1).
    ! f1[u,x] = 2N(2 du/dx + dv/dy) ; f1[u,y] = N(du/dy + dv/dx)
    ! f1[v,x] = N(du/dy + dv/dx)    ; f1[v,y] = 2N(2 dv/dy + du/dx)
    g3( 1:16) = 0._c_double
    g3(  1) = 4._c_double * N   ! d f1[u,x] / d(du/dx)   c,c',d,d' = 0,0,0,0
    g3(  6) = 2._c_double * N   ! d f1[u,x] / d(dv/dy)          = 0,1,0,1
    g3(  4) =               N   ! d f1[u,y] / d(du/dy)          = 0,0,1,1
    g3(  7) =               N   ! d f1[u,y] / d(dv/dx)          = 0,1,1,0
    g3( 10) =               N   ! d f1[v,x] / d(du/dy)          = 1,0,0,1
    g3( 13) =               N   ! d f1[v,x] / d(dv/dx)          = 1,1,0,0
    g3( 11) = 2._c_double * N   ! d f1[v,y] / d(du/dx)          = 1,0,1,0
    g3( 16) = 4._c_double * N   ! d f1[v,y] / d(dv/dy)          = 1,1,1,1

  end subroutine SSA_FEM_PETSc_g3

  ! ===== Restart NetCDF files =====

  subroutine create_restart_file_old_SSA_FEM_PETSc( self)
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    ! Phase 8 will reuse the SSA restart-file layout.
  end subroutine create_restart_file_old_SSA_FEM_PETSc

  subroutine write_to_restart_file_old_SSA_FEM_PETSc( self, time)
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in   ) :: self
    real(dp),                                          intent(in   ) :: time
    ! Phase 8 will reuse the SSA restart-file layout.
  end subroutine write_to_restart_file_old_SSA_FEM_PETSc

end module momentum_balance_solver_SSA_FEM_PETSc
