module momentum_balance_solver_SSA_FEM_PETSc

#include <petsc/finclude/petscsys.h>

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA),
  ! discretised and solved entirely with PETSc: DMPlex + PetscFE for the discretisation
  ! (continuous P1 velocity on the mesh vertices) and PetscSNES for the non-linear solve.
  !
  ! This is the finite-element counterpart of the existing finite-difference-based SSA
  ! solver (momentum_balance_solver_SSA). It is selected with
  !   choice_stress_balance_approximation = 'SSA_FEM_PETSc'
  ! and is added alongside the existing solvers without affecting any of them.
  !
  ! The full derivation of the f0/f1/g0/g3 callbacks below - from the two coupled
  ! SSA PDEs, via the weak form, to the pointwise residual/Jacobian expressions and
  ! their flat storage layout - is in SSA_FEM_PETSc_weak_form_derivation.md in the
  ! repository root. Implementation progress is tracked in
  ! SSA_PetscFE_SNES_implementation_plan.md. Current state: Phase 3+ - BOTH
  ! non-linearities of the SSA are
  ! evaluated pointwise inside the residual (the shear-thinning Glen viscosity and
  ! the velocity-dependent basal friction law), so a single Newton solve replaces
  ! the whole viscosity/friction Picard iteration. Analytic Jacobians are provided.
  !
  ! Weak form (PETSc convention  residual = integral( f0 . phi + f1 : grad(phi) ) = 0 ),
  ! one 2-component vector field (u, v). A 5-component P1 auxiliary field carries the
  ! velocity-independent data  a = [Abar, H, tauc_eff, tau_dx, tau_dy]:
  !
  !   eta   = 1/2 Abar^(-1/n) (eps_eff^2 + eps0)^((1-n)/(2n))   (eps_eff^2 from grad u)
  !   N     = eta * H
  !   f1[c,d] = 2 N D[c,d]     with  D = [ 2ux+vy , (uy+vx)/2 ; (uy+vx)/2 , 2vy+ux ]
  !   beta  = tauc_eff |u|^(1/p-1) (|u|+u_t)^(-1/p)             (Zoet-Iverson, from |u|)
  !   f0[c] = beta u_c - tau_d,c
  !
  ! eps0, n and the Zoet-Iverson parameters (p, u_t, delta_v, beta_max) come from
  ! PetscDSSetConstants. tauc_eff is the till yield stress scaled by the sub-grid
  ! grounded fraction, so friction vanishes under floating ice. The pointwise
  ! sliding relation is SSA_FEM_PETSc_sliding_beta below (currently Zoet-Iverson
  ! only; see the TODO there).
  !
  ! Non-dimensionalisation (pulled forward from Phase 6): the PetscFE field, and
  ! therefore self%sol, hold u_hat = u / velocity_scale (dimensionless), not the
  ! physical velocity. f0/f1/g0/g3 convert u_hat (and grad u_hat) to physical
  ! units, run the physics above unchanged, then divide by stress_scale (f0, f1)
  ! or velocity_scale/stress_scale (g0, g3) so the assembled residual/Jacobian is
  ! O(1) instead of spanning ~1e-8 to ~1e10 in raw SI units. Only the copy-back to
  ! u_vav_a/v_vav_a (in solve_SSA_Newton) converts back to physical m/yr.
  !
  ! Solver/preconditioner (Phase 6): C%SSA_FEM_PETSc_pc_type selects 'lu' (direct,
  ! the default), 'gamg' (algebraic multigrid) or 'bjacobi'; a rigid-body
  ! near-null-space is attached unconditionally (needed for 'gamg', harmless
  ! otherwise). C%SSA_FEM_PETSc_snes_{rtol,abstol,maxits} configure the Newton
  ! solve.
  !
  ! Not done yet: essential boundary conditions and the ice-front back-pressure
  ! (Phase 5).
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
    tPetscFE, tPetscDS, tDMLabel, tPetscSection, tMatNullSpace, &
    PetscFECreateLagrange, PetscFEDestroy, PetscObjectSetName, DMSetField, DMCreateDS, DMGetDS, &
    PetscDSSetConstants, DMCreateMatrix, DMCreateGlobalVector, DMCreateLocalVector, &
    DMGlobalToLocalBegin, DMGlobalToLocalEnd, DMLocalToGlobalBegin, DMLocalToGlobalEnd, &
    DMGetLocalSection, DMGetLabel, DMDestroy, &
    DMPlexGetDepthStratum, DMLabelGetValue, PetscSectionGetOffset, DMGetCoordinates, &
    DMPlexMarkBoundaryFaces, DMCreateLabel, &
    VecSet, VecDestroy, VecGetArrayRead, VecRestoreArrayRead, VecSetValues, VecAssemblyBegin, &
    VecAssemblyEnd, MatDestroy, INSERT_VALUES, MatSetBlockSize, MatNullSpaceCreateRigidBody, &
    MatSetNearNullSpace, MatNullSpaceDestroy, &
    SNESCreate, SNESSetDM, SNESSetType, SNESSetTolerances, SNESGetKSP, SNESSolve, SNESDestroy, &
    SNESGetIterationNumber, SNESNEWTONLS, KSPSetType, KSPGetPC, KSPPREONLY, KSPGMRES, PCSetType, PCLU, &
    KSPSetTolerances, SNESGetLinearSolveIterations
  use mpi_f08, only: MPI_ALLTOALL, MPI_ALLTOALLV, MPI_ALLREDUCE, MPI_COMM_WORLD, MPI_IN_PLACE, &
    MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_MAX, MPI_LOR, MPI_LOGICAL
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
  use constitutive_equation, only: calc_ice_rheology_Glen
  use mesh_zeta, only: vertical_average
  use sliding_laws, only: calc_basal_friction_coefficient
  use petsc_dmplex, only: mesh_to_dmplex_masked, dmplex_upsy_vertex_id_label_name
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all
  use mesh_disc_apply_operators, only: map_a_b_2D, ddx_a_a_2D, ddy_a_a_2D

  implicit none

  private

  public :: type_momentum_balance_solver_SSA_FEM_PETSc

  ! Layout of the 5-component PetscFE auxiliary field
  integer, parameter :: i_Abar  = 1   ! vertically averaged flow factor A         [Pa^-n yr^-1]
  integer, parameter :: i_H     = 2   ! ice thickness (>= 0.1 m)                  [m]
  integer, parameter :: i_tauc  = 3   ! till yield stress * fraction_gr**exponent [Pa]
  integer, parameter :: i_taudx = 4   ! driving stress, x                        [Pa]
  integer, parameter :: i_taudy = 5   ! driving stress, y                        [Pa]
  integer, parameter :: n_aux_comp = 5

  ! Index of the Zoet-Iverson parameters within the PetscDS constants array
  integer, parameter :: ic_eps0 = 1, ic_nglen = 2, ic_ZIp = 3, ic_ZIut = 4, ic_dv = 5, ic_betamax = 6
  integer, parameter :: n_ds_constants = 6

  ! Non-dimensionalisation. The PetscFE velocity field holds the DIMENSIONLESS
  ! u_hat = u / velocity_scale; the residual/Jacobian callbacks convert u_hat to
  ! physical u internally, run the physics exactly as before, then divide the
  ! result by stress_scale (f0, f1) or velocity_scale/stress_scale (g0, g3) - see
  ! the "Non-dimensionalisation" section of SSA_FEM_PETSc_weak_form_derivation.md.
  ! self%sol and the SNES tolerances therefore live in these units too; only the
  ! copy back to u_vav_a/v_vav_a (in solve_SSA_Newton) converts to physical m/yr.
  real(dp), parameter :: velocity_scale = 1.0e3_dp   ! [m yr^-1] ~ typical fast-flow speed
  real(dp), parameter :: stress_scale   = 1.0e5_dp   ! [Pa] ~ typical driving stress

  ! d D[m] / d(grad u)[k], with D and grad u in the layout
  ! (1,2,3,4) = (xx, xy, yx, yy) resp. (du/dx, du/dy, dv/dx, dv/dy).
  ! D = [ 2ux+vy , (uy+vx)/2 , (uy+vx)/2 , 2vy+ux ]
  real(dp), parameter, dimension(4,4) :: dD_dgradu = reshape([ &
    2._dp, 0._dp, 0._dp, 1._dp,   &   ! d D[.]/d(du/dx)
    0._dp, 0.5_dp, 0.5_dp, 0._dp, &   ! d D[.]/d(du/dy)
    0._dp, 0.5_dp, 0.5_dp, 0._dp, &   ! d D[.]/d(dv/dx)
    1._dp, 0._dp, 0._dp, 2._dp ], &   ! d D[.]/d(dv/dy)
    [4,4])

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SSA_FEM_PETSc

    ! Persistent PETSc objects (rebuilt on remap)
    type(tDM)      :: dm            ! primary DM: the P1 velocity field
    type(tPetscFE) :: fe
    type(tSNES)    :: snes
    type(tMat)     :: jac
    type(tVec)     :: sol           ! global solution vector, kept as the SNES initial guess
    type(tDM)      :: dm_aux        ! clone of dm carrying the auxiliary field
    type(tPetscFE) :: fe_aux
    type(tVec)     :: aux_vec       ! local vector of dm_aux: [Abar, H, tauc_eff, tau_dx, tau_dy] per vertex
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
      procedure, private :: solve_SSA_Newton
      procedure, private :: calc_auxiliary_fields
      procedure, private :: calc_ice_covered_triangle_mask

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

    if (par%primary) write(0,'(A)') '    NOTE: the SSA_FEM_PETSc solver solves the fully non-linear SSA ' // &
      'with one Newton solve (pointwise Glen viscosity + sliding law); natural boundary conditions only.'

    ! The basal friction non-linearity is evaluated pointwise in the residual, which
    ! currently only covers the Zoet-Iverson sliding law (see SSA_FEM_PETSc_sliding_beta).
    select case (C%choice_sliding_law)
    case ('Zoet-Iverson', 'no_sliding')
      ! supported
    case default
      call crash('SSA_FEM_PETSc evaluates the basal friction law pointwise in the residual and ' // &
        'currently only implements "Zoet-Iverson" (got "' // trim( C%choice_sliding_law) // '"). ' // &
        'TODO: port the other sliding laws - see SSA_FEM_PETSc_sliding_beta in this module.')
    end select

    ! The PETSc DMPlex/FE/DS/SNES machinery is built on the ice-covered sub-mesh
    ! (see build_petsc_objects), which needs geom - not available here. It is
    ! instead (re)built at the top of every run() call, since the ice mask (and
    ! therefore the sub-mesh) can change every timestep.

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_initialise

  subroutine momentum_balance_solver_SSA_FEM_PETSc_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the non-linear SSA with PetscSNES (Newton on
    !< the pointwise Glen viscosity), inside a light Picard loop that refreshes the
    !< velocity-dependent basal friction coefficient.

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
    character(len=*), parameter                        :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_run'
    integer                                            :: ierr, newton_its, vi
    logical                                            :: grounded_ice_exists
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2)   :: A_flow_vav_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2, n_aux_comp) :: coeffs
    real(dp), dimension(self%mesh%nz)                  :: A_prof
    real(dp)                                           :: uabs, umax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! If there is no grounded ice or no sliding, there is nothing to solve
    grounded_ice_exists = any( geom%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists .or. C%choice_sliding_law == 'no_sliding') then
      self%u_vav_a = 0._dp; self%v_vav_a = 0._dp
      self%u_vav_b = 0._dp; self%v_vav_b = 0._dp
      ! self%sol (and the rest of the PETSc machinery) may not exist yet if this
      ! is the very first run() call and there is no ice anywhere yet.
      if (self%petsc_is_built) then
        PetscCall( VecSet( self%sol, 0._dp, ierr))
      end if
      self%n_visc_its = 0
      self%n_Axb_its  = 0
      call finalise_routine( routine_name)
      return
    end if

    ! The DMPlex/FE/DS/SNES machinery lives on the ice-covered sub-mesh only (so its
    ! own topological boundary is the ice margin, for the ice-front back-pressure BC
    ! to come). Since the ice mask can change every timestep, the whole lot is torn
    ! down and rebuilt from scratch here, every call - see the "rebuild every
    ! timestep" checklist in SSA_PetscFE_SNES_implementation_plan.md.
    if (self%petsc_is_built) call self%destroy_petsc_objects()
    call self%build_petsc_objects( geom)

    ! DOF numbering is not stable across independent DMPlex rebuilds, so a fresh
    ! self%sol always starts at zero (set inside build_petsc_objects) even though
    ! self%u_vav_a/v_vav_a still hold the previous physical solution. Re-seed
    ! self%sol from those to keep warm-starting Newton across timesteps (this
    ! only affects convergence speed, not correctness).
    block
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2, 2) :: sol_seed
      type(tVec) :: sol_local
      sol_seed(:,1) = self%u_vav_a / velocity_scale
      sol_seed(:,2) = self%v_vav_a / velocity_scale
      PetscCall( DMCreateLocalVector( self%dm, sol_local, ierr))
      PetscCall( VecSet( sol_local, 0._dp, ierr))
      call fill_PETSc_aux_from_mesh_vertices( self%dm, self%dm, sol_local, self%mesh, sol_seed)
      PetscCall( DMLocalToGlobalBegin( self%dm, sol_local, INSERT_VALUES, self%sol, ierr))
      PetscCall( DMLocalToGlobalEnd(   self%dm, sol_local, INSERT_VALUES, self%sol, ierr))
      PetscCall( VecDestroy( sol_local, ierr))
    end block

    ! Vertically averaged flow factor A - velocity-independent, computed once
    call calc_ice_rheology_Glen( self%mesh, ice, geom)
    do vi = self%mesh%vi1, self%mesh%vi2
      A_prof = ice%A_flow( vi,:)
      A_flow_vav_a( vi) = vertical_average( self%mesh%zeta, A_prof)
    end do

    ! Assemble the velocity-independent auxiliary coefficients and hand them to PETSc.
    ! Both non-linearities of the SSA - the Glen viscosity AND the basal friction law -
    ! are evaluated pointwise inside the residual, so a single Newton solve suffices;
    ! there is no outer Picard iteration. The solve is warm-started from self%sol.
    call self%calc_auxiliary_fields( ice, geom, bed_roughness, A_flow_vav_a, coeffs)
    call fill_PETSc_aux_from_mesh_vertices( self%dm, self%dm_aux, self%aux_vec, self%mesh, coeffs)
    ierr = dm_set_auxiliary_vec( self%dm%v, 0_c_intptr_t, 0_c_int, 0_c_int, self%aux_vec%v)
    CHKERRQ( ierr)

    ! n_visc_its: Newton (nonlinear) iterations; n_Axb_its: cumulative inner
    ! linear-solve (KSP) iterations, set inside solve_SSA_Newton
    call self%solve_SSA_Newton( newton_its)
    self%n_visc_its = newton_its

    if (any( isnan( self%u_vav_a)) .or. any( isnan( self%v_vav_a))) &
      call crash('SSA_FEM_PETSc: NaN in the velocity solution')

    ! Limit velocities for the exposed result
    do vi = self%mesh%vi1, self%mesh%vi2
      uabs = sqrt( self%u_vav_a( vi)**2 + self%v_vav_a( vi)**2)
      if (uabs > C%vel_max) then
        self%u_vav_a( vi) = self%u_vav_a( vi) * C%vel_max / uabs
        self%v_vav_a( vi) = self%v_vav_a( vi) * C%vel_max / uabs
      end if
    end do

    umax = maxval( sqrt( self%u_vav_a**2 + self%v_vav_a**2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, umax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (par%primary) write(0,'(A,I0,A,I0,A,ES10.3)') '    SSA_FEM_PETSc: Newton its = ', newton_its, &
      ', cumulative KSP its = ', self%n_Axb_its, ', max speed = ', umax

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

    ! Velocities on the vertices
    do vi = self%mesh%vi1, self%mesh%vi2
      vel%u_vav( vi) = self%u_vav_a( vi)
      vel%v_vav( vi) = self%v_vav_a( vi)
    end do

    ! Velocities on the triangles
    do ti = self%mesh%ti1, self%mesh%ti2
      vel%u_3D_b( ti,:) = self%u_vav_b( ti)
      vel%v_3D_b( ti,:) = self%v_vav_b( ti)
    end do

    ! Strain rates on the vertices - still zero until they are computed from the FE field
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

    ! Discard the PETSc machinery built on the old mesh (self%mesh has already
    ! been repointed to mesh_new by remap_model). It is not rebuilt here: since
    ! run() rebuilds it unconditionally from geom every call anyway (the ice mask
    ! can change every timestep regardless of remapping), building it here too
    ! would just be wasted work. Velocities are reset to zero; a proper remap of
    ! u_vav via the a-grid follows in a later phase.
    if (self%petsc_is_built) call self%destroy_petsc_objects()

    call reallocate_bounds( self%u_vav_a, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%v_vav_a, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%u_vav_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%v_vav_b, mesh_new%ti1, mesh_new%ti2)
    self%u_vav_a = 0._dp; self%v_vav_a = 0._dp
    self%u_vav_b = 0._dp; self%v_vav_b = 0._dp

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'SSA_FEM_PETSc'
  end function get_momentum_balance_solver_name

  ! ===== PETSc object lifecycle =====

  subroutine build_petsc_objects( self, geom)
    !< Build the DMPlex, the P1 vector PetscFE field, the auxiliary-field DM, the
    !< PetscDS weak form and the SNES, restricted to the ice-covered sub-mesh (so
    !< that its own topological exterior boundary - the ice margin - can later
    !< carry a natural ice-front back-pressure BC). Called anew every run(), since
    !< the ice mask (and therefore the sub-mesh) can change every timestep.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    class(atype_ice_geometry_model_data),               intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'build_petsc_objects'
    logical, dimension(self%mesh%nTri) :: mask_tri
    type(tPetscDS)                 :: ds
    type(tPetscObject)             :: fe_object, fe_aux_object
    type(tKSP)                     :: ksp
    type(tPC)                      :: pc
    integer                        :: ierr
    integer(c_intptr_t)            :: no_context
    real(dp), dimension(n_ds_constants) :: ds_constants

    call init_routine( routine_name)

    no_context = 0_c_intptr_t

    ! DMPlex from the ice-covered subset of the UFEMISM mesh (creates + distributes;
    ! preserves UPSY vertex IDs in the 'upsy_vertex_id' DMLabel)
    call self%calc_ice_covered_triangle_mask( geom, mask_tri)
    call mesh_to_dmplex_masked( self%mesh, mask_tri, self%dm)

    ! Primary field: one 2-component P1 Lagrange velocity
    PetscCall( PetscFECreateLagrange( PETSC_COMM_SELF, 2, 2, PETSC_TRUE, 1, -1, self%fe, ierr))
    PetscCall( PetscObjectSetName( self%fe, 'velocity', ierr))
    PetscObjectSpecificCast( fe_object, self%fe)
    PetscCall( DMSetField( self%dm, 0, PETSC_NULL_DMLABEL, fe_object, ierr))
    PetscCall( DMCreateDS( self%dm, ierr))
    PetscCall( DMGetDS( self%dm, ds, ierr))

    ! Uniform scalars for the pointwise functions:
    ! [eps_sq_0, Glen n, Zoet-Iverson p, Zoet-Iverson u_t, delta_v, beta_max]
    ds_constants( ic_eps0)    = C%Glens_flow_law_epsilon_sq_0
    ds_constants( ic_nglen)   = C%Glens_flow_law_exponent
    ds_constants( ic_ZIp)     = C%slid_ZI_p
    ds_constants( ic_ZIut)    = C%slid_ZI_ut
    ds_constants( ic_dv)      = C%slid_delta_v
    ds_constants( ic_betamax) = C%slid_beta_max
    PetscCall( PetscDSSetConstants( ds, n_ds_constants, ds_constants, ierr))

    ! Auxiliary field: [Abar, H, tauc_eff, tau_dx, tau_dy], P1, on a clone of the primary DM
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
    PetscCall( DMCreateGlobalVector( self%dm, self%sol, ierr))
    PetscCall( VecSet( self%sol, 0._dp, ierr))

    ! Solver configuration: Newton line search; KSP/PC choice from config.
    PetscCall( SNESSetType( self%snes, SNESNEWTONLS, ierr))
    call SNESSetTolerances( self%snes, C%SSA_FEM_PETSc_snes_abstol, C%SSA_FEM_PETSc_snes_rtol, &
      1.0e-12_dp, C%SSA_FEM_PETSc_snes_maxits, 1000, ierr)
    CHKERRQ( ierr)
    PetscCall( SNESGetKSP( self%snes, ksp, ierr))
    select case (C%SSA_FEM_PETSc_pc_type)
    case ('lu')
      ! Direct solve (MUMPS on more than one rank, as in ct_PETSc_SNES_Poisson).
      ! Slower at scale than an iterative solve, but does not need the near-null-space
      ! below and has been the validated default through Phases 1-3.
      PetscCall( KSPSetType( ksp, KSPPREONLY, ierr))
      PetscCall( KSPGetPC( ksp, pc, ierr))
      PetscCall( PCSetType( pc, PCLU, ierr))
    case ('gamg')
      ! Algebraic multigrid; needs the rigid-body near-null-space (set below) to
      ! coarsen a vector-valued (elasticity-like) operator well.
      PetscCall( KSPSetType( ksp, KSPGMRES, ierr))
      call KSPSetTolerances( ksp, 1.0e-8_dp, 1.0e-12_dp, 1.0e5_dp, 10000, ierr)
      CHKERRQ( ierr)
      PetscCall( KSPGetPC( ksp, pc, ierr))
      PetscCall( PCSetType( pc, 'gamg', ierr))
    case ('bjacobi')
      ! Matches the finite-difference SSA/DIVA solver's own default
      ! (solve_matrix_equation_PETSc): 'gmres' + 'bjacobi'.
      PetscCall( KSPSetType( ksp, KSPGMRES, ierr))
      call KSPSetTolerances( ksp, 1.0e-8_dp, 1.0e-12_dp, 1.0e5_dp, 10000, ierr)
      CHKERRQ( ierr)
      PetscCall( KSPGetPC( ksp, pc, ierr))
      PetscCall( PCSetType( pc, 'bjacobi', ierr))
    case default
      call crash('SSA_FEM_PETSc: unknown SSA_FEM_PETSc_pc_type "' // trim( C%SSA_FEM_PETSc_pc_type) // '"')
    end select

    ! Rigid-body near-null-space (2 translations + 1 rotation in 2D), from the DM's
    ! own coordinates. Harmless for 'lu'/'bjacobi'; needed for 'gamg' to perform
    ! acceptably on this vector-valued, elasticity-like operator.
    block
      type(tVec)          :: coords
      type(tMatNullSpace)  :: nullsp
      PetscCall( MatSetBlockSize( self%jac, 2, ierr))
      PetscCall( DMGetCoordinates( self%dm, coords, ierr))
      PetscCall( MatNullSpaceCreateRigidBody( coords, nullsp, ierr))
      PetscCall( MatSetNearNullSpace( self%jac, nullsp, ierr))
      PetscCall( MatNullSpaceDestroy( nullsp, ierr))
    end block

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
    PetscCall( VecDestroy( self%sol, ierr))
    PetscCall( MatDestroy( self%jac, ierr))
    PetscCall( SNESDestroy( self%snes, ierr))
    PetscCall( PetscFEDestroy( self%fe, ierr))
    PetscCall( DMDestroy( self%dm, ierr))
    self%petsc_is_built = .false.

    call finalise_routine( routine_name)

  end subroutine destroy_petsc_objects

  subroutine solve_SSA_Newton( self, snes_its)
    !< One Newton solve of the non-linear SSA with the currently attached (frozen-beta)
    !< coefficients. Warm-started from, and written back to, self%sol.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    integer,                                           intent(  out) :: snes_its

    ! Local variables:
    character(len=*), parameter :: routine_name = 'solve_SSA_Newton'
    integer                     :: ierr
    integer(c_int)              :: snes_reason

    call init_routine( routine_name)

    PetscCall( SNESSolve( self%snes, PETSC_NULL_VEC, self%sol, ierr))
    PetscCall( SNESGetIterationNumber( self%snes, snes_its, ierr))
    ierr = snes_get_converged_reason( self%snes%v, snes_reason)
    CHKERRQ( ierr)
    if (snes_reason < 0) &
      call crash('SSA_FEM_PETSc: SNES diverged (SNESConvergedReason = {int_01})', int_01 = int( snes_reason))

    ! Cumulative inner linear-solve (KSP) iteration count for this Newton solve
    PetscCall( SNESGetLinearSolveIterations( self%snes, self%n_Axb_its, ierr))

    ! self%sol holds the dimensionless u_hat = u/velocity_scale; convert to physical m/yr
    call copy_PETSc_solution_to_mesh_vertices_vec2( self%dm, self%sol, self%mesh, self%u_vav_a, self%v_vav_a)
    self%u_vav_a = velocity_scale * self%u_vav_a
    self%v_vav_a = velocity_scale * self%v_vav_a

    call finalise_routine( routine_name)

  end subroutine solve_SSA_Newton

  ! ===== Coefficient (auxiliary-field) calculation =====

  subroutine calc_auxiliary_fields( self, ice, geom, bed_roughness, A_flow_vav_a, coeffs)
    !< Compute the auxiliary-field coefficients on the mesh vertices:
    !< [Abar, H, beta, tau_dx, tau_dy]. All are velocity-independent except beta,
    !< which is frozen at the current velocity solution (sliding law), scaled by the
    !< sub-grid grounded fraction.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in   ) :: self
    class(atype_ice_model_data),                       intent(inout) :: ice
    class(atype_ice_geometry_model_data),              intent(in   ) :: geom
    type(type_bed_roughness_model),                    intent(in   ) :: bed_roughness
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2),  intent(in   ) :: A_flow_vav_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2, n_aux_comp), intent(  out) :: coeffs

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_auxiliary_fields'
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: dHs_dx_a, dHs_dy_a
    integer                                         :: vi
    real(dp)                                        :: tauc

    call init_routine( routine_name)

    ! Driving stress on the vertices: tau_d = -rho g H grad(Hs)
    call ddx_a_a_2D( self%mesh, geom%Hs, dHs_dx_a)
    call ddy_a_a_2D( self%mesh, geom%Hs, dHs_dy_a)

    ! Till yield stress from the sliding law (velocity-independent). Its basal_friction
    ! _coefficient output is discarded: the friction non-linearity is evaluated pointwise
    ! in the residual. tauc is scaled by the sub-grid grounded fraction so that friction
    ! vanishes under floating ice (a-grid analogue of calc_applied_basal_friction_coefficient).
    call calc_basal_friction_coefficient( self%mesh, geom, bed_roughness, self%u_vav_a, self%v_vav_a, &
      ice%effective_pressure, ice%till_yield_stress, ice%basal_friction_coefficient)

    do vi = self%mesh%vi1, self%mesh%vi2
      tauc = ice%till_yield_stress( vi)
      if (C%do_GL_subgrid_friction) &
        tauc = tauc * geom%fraction_gr( vi)**C%subgrid_friction_exponent_on_B_grid

      coeffs( vi, i_Abar)  = A_flow_vav_a( vi)
      coeffs( vi, i_H)     = max( 0.1_dp, geom%Hi( vi))
      coeffs( vi, i_tauc)  = tauc
      coeffs( vi, i_taudx) = -ice_density * grav * geom%Hi( vi) * dHs_dx_a( vi)
      coeffs( vi, i_taudy) = -ice_density * grav * geom%Hi( vi) * dHs_dy_a( vi)
    end do

    call finalise_routine( routine_name)

  end subroutine calc_auxiliary_fields

  subroutine calc_ice_covered_triangle_mask( self, geom, mask_tri)
    !< Which triangles have all 3 vertices ice-covered (grounded or floating) - the
    !< sub-mesh build_petsc_objects builds self%dm on via mesh_to_dmplex_masked.
    !< Recomputed every run() call since the ice mask changes as the ice sheet
    !< evolves.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in ) :: self
    class(atype_ice_geometry_model_data),               intent(in ) :: geom
    logical, dimension(self%mesh%nTri),                 intent(out) :: mask_tri

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'calc_ice_covered_triangle_mask'
    logical, dimension(:), allocatable :: mask_grounded_tot, mask_floating_tot, mask_ice_a_full
    integer                            :: ti

    call init_routine( routine_name)

    ! Gather the (distributed-shared-memory) "has ice" vertex masks to full,
    ! replicated arrays - mesh connectivity (mesh%Tri) is itself fully
    ! replicated, so the mask needs to be too (see mesh_to_dmplex_masked).
    allocate( mask_grounded_tot( self%mesh%nV))
    allocate( mask_floating_tot( self%mesh%nV))
    call gather_dist_shared_to_all( self%mesh%pai_V, geom%mask_grounded_ice, mask_grounded_tot)
    call gather_dist_shared_to_all( self%mesh%pai_V, geom%mask_floating_ice, mask_floating_tot)
    allocate( mask_ice_a_full( self%mesh%nV))
    mask_ice_a_full = mask_grounded_tot .or. mask_floating_tot

    ! A triangle counts as ice-covered if all 3 of its vertices do (matches the
    ! Hi > 0 rule the existing graph abstraction uses for the same classification)
    do ti = 1, self%mesh%nTri
      mask_tri( ti) = all( mask_ice_a_full( self%mesh%Tri( ti,:)))
    end do

    call finalise_routine( routine_name)

  end subroutine calc_ice_covered_triangle_mask

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
    ! Vertices with ncopies == 0 are outside the ice-covered sub-mesh (self%dm
    ! only spans the ice-covered triangles, see build_petsc_objects) - they have
    ! no SSA solution and are left at their u_a = v_a = 0 default set above.
    where (ncopies > 0)
      u_a = u_a / real( ncopies, dp)
      v_a = v_a / real( ncopies, dp)
    end where

  end subroutine copy_PETSc_solution_to_mesh_vertices_vec2

  subroutine fill_PETSc_aux_from_mesh_vertices( dm_topo, dm_aux, aux_vec, mesh, coeffs)
    !< Scatter the per-vertex coefficient array coeffs(vi1:vi2, 1:nc) into the
    !< nc-component local vector aux_vec of dm_aux (nc taken from coeffs itself, so
    !< this doubles as both the n_aux_comp-component SSA coefficient scatter and the
    !< 2-component self%sol warm-start re-seed scatter). Inverse of
    !< copy_PETSc_solution_to_mesh_vertices_vec2: each rank requests, for its local
    !< DMPlex vertices, the coefficients from the UFEMISM process that owns that vertex.

    ! In/output variables:
    type(tDM),                              intent(in   ) :: dm_topo   ! for topology + 'upsy_vertex_id' label
    type(tDM),                              intent(in   ) :: dm_aux    ! for the local section
    type(tVec),                             intent(inout) :: aux_vec
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:, :), intent(in) :: coeffs

    ! Local variables:
    type(tPetscSection)                 :: aux_section
    type(tDMLabel)                      :: upsy_vertex_id_label
    integer, dimension(:), allocatable  :: send_counts, recv_counts, send_displ, recv_displ, send_pos
    integer, dimension(:), allocatable  :: send_countsc, recv_countsc, send_displc, recv_displc
    integer, dimension(:), allocatable  :: req_vi, recv_req_vi
    integer, dimension(:), allocatable  :: local_pt, local_slot
    real(dp), dimension(:), allocatable :: reply_vals, recv_reply
    integer, dimension(:), allocatable  :: idxc
    real(dp), dimension(:), allocatable :: valsc
    integer :: ierr, point, vstart, vend, vi, dest, ip, si, ns, nr, k, nlv, off, nc, j

    nc = size( coeffs, 2)
    allocate( idxc( nc), valsc( nc))

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
    allocate( recv_req_vi( max( 1, nr)), reply_vals( max( 1, nc*nr)))

    call MPI_ALLTOALLV( req_vi, send_counts, send_displ, MPI_INTEGER, &
      recv_req_vi, recv_counts, recv_displ, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Fill the replies with the n_aux_comp coefficients for each requested (locally owned) vertex
    do k = 1, nr
      vi = recv_req_vi( k)
      if (vi < mesh%vi1 .or. vi > mesh%vi2) call crash('aux request sent to the wrong UPSY process')
      reply_vals( (k-1)*nc+1 : k*nc) = coeffs( vi, 1:nc)
    end do

    ! Send the replies back (n_aux_comp doubles per requested vertex)
    allocate( recv_reply( max( 1, nc*ns)))
    allocate( send_countsc( 0:par%n-1), recv_countsc( 0:par%n-1))
    allocate( send_displc ( 0:par%n-1), recv_displc ( 0:par%n-1))
    send_countsc = nc * send_counts
    recv_countsc = nc * recv_counts
    send_displc  = nc * send_displ
    recv_displc  = nc * recv_displ
    call MPI_ALLTOALLV( reply_vals, recv_countsc, recv_displc, MPI_DOUBLE_PRECISION, &
      recv_reply, send_countsc, send_displc, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    ! Write into the local auxiliary vector at each vertex' section offset
    do k = 1, nlv
      point = local_pt( k)
      si    = local_slot( k)
      PetscCall( PetscSectionGetOffset( aux_section, point, off, ierr))
      do j = 1, nc
        idxc( j)  = off + j - 1
        valsc( j) = recv_reply( si*nc + j)
      end do
      PetscCall( VecSetValues( aux_vec, nc, idxc, valsc, INSERT_VALUES, ierr))
    end do
    PetscCall( VecAssemblyBegin( aux_vec, ierr))
    PetscCall( VecAssemblyEnd( aux_vec, ierr))

  end subroutine fill_PETSc_aux_from_mesh_vertices

  ! ===== PetscDS pointwise weak-form functions (bind(C)) =====
  !
  ! PETSc assembles  residual = integral( f0 . phi + f1 : grad(phi) ).
  ! Field 0 is the 2-vector velocity (u, v); dim = 2, Nc = 2.
  ! Auxiliary field a(1..5) = [Abar, H, beta, tau_dx, tau_dy].
  ! constants(1..2) = [eps_sq_0, Glen exponent n].
  ! grad u layout u_x(1..4) = du/dx, du/dy, dv/dx, dv/dy.
  ! Strain tensor D(1..4) = (xx, xy, yx, yy) = [2ux+vy, (uy+vx)/2, (uy+vx)/2, 2vy+ux].
  ! f1(1..4) = (u,x),(u,y),(v,x),(v,y) = 2 N D(m).

  subroutine SSA_FEM_PETSc_strain( u_x, eps0, n_glen, Abar, D, eps2, eta)
    !< Shared helper: strain tensor D, regularised effective strain rate^2 and viscosity.
    real(c_double), intent(in)  :: u_x(4), eps0, n_glen, Abar
    real(c_double), intent(out) :: D(4), eps2, eta
    real(c_double) :: du_dx, du_dy, dv_dx, dv_dy

    du_dx = u_x(1); du_dy = u_x(2); dv_dx = u_x(3); dv_dy = u_x(4)
    D(1) = 2._c_double*du_dx + dv_dy
    D(2) = 0.5_c_double*(du_dy + dv_dx)
    D(3) = D(2)
    D(4) = 2._c_double*dv_dy + du_dx
    eps2 = du_dx**2 + dv_dy**2 + du_dx*dv_dy + 0.25_c_double*(du_dy + dv_dx)**2 + eps0
    eta  = 0.5_c_double * Abar**(-1._c_double/n_glen) * eps2**((1._c_double - n_glen)/(2._c_double*n_glen))
  end subroutine SSA_FEM_PETSc_strain

  subroutine SSA_FEM_PETSc_sliding_beta( u, v, tauc, c_values, beta, dbeta_duabs)
    !< Pointwise Zoet & Iverson (2020) basal friction relation: returns the friction
    !< coefficient beta (such that tau_b = beta * u) as a function of |u|, and
    !< d beta / d|u| (for the analytic Jacobian). tauc is the (grounded-fraction-
    !< scaled) till yield stress; the ZI parameters p, u_t, delta_v and beta_max come
    !< from the PetscDS constants array.
    !<
    !< TODO: this duplicates the kernel of
    !<   sliding_laws.f90 :: calc_sliding_law_ZoetIverson.
    !< It should eventually become a shared *pointwise* routine living in
    !< sliding_laws, called both by the array-based finite-difference solvers and
    !< from here, with the other sliding laws (Weertman, Budd, Coulomb, Tsai2015,
    !< Schoof2005) given the same treatment. Until then, the guard in
    !< momentum_balance_solver_SSA_FEM_PETSc_initialise restricts this solver to
    !< the Zoet-Iverson law.
    real(c_double), intent(in)  :: u, v, tauc
    real(c_double), intent(in)  :: c_values(:)
    real(c_double), intent(out) :: beta, dbeta_duabs
    real(c_double) :: uabs, aexp, bexp, ZIp, ZIut, dv, betamax

    ZIp     = c_values( ic_ZIp)
    ZIut    = c_values( ic_ZIut)
    dv      = c_values( ic_dv)
    betamax = c_values( ic_betamax)

    uabs = sqrt( dv**2 + u**2 + v**2)
    aexp = 1._c_double / ZIp - 1._c_double        ! exponent on |u|
    bexp = -1._c_double / ZIp                     ! exponent on (|u| + u_t)
    beta = tauc * uabs**aexp * (uabs + ZIut)**bexp
    if (beta >= betamax) then
      beta        = betamax
      dbeta_duabs = 0._c_double
    else
      dbeta_duabs = beta * (aexp / uabs + bexp / (uabs + ZIut))
    end if

  end subroutine SSA_FEM_PETSc_sliding_beta

  subroutine SSA_FEM_PETSc_f0( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f0) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f0
    real(c_double), value :: time
    real(c_double), pointer :: u_hat_values(:), a_values(:), c_values(:), f0_values(:)
    real(c_double)          :: u_phys(2), beta, dbeta_duabs

    call c_f_pointer( u, u_hat_values, [2])
    call c_f_pointer( a, a_values, [n_aux_comp])
    call c_f_pointer( constants, c_values, [int( nconstants)])
    call c_f_pointer( f0, f0_values, [2])

    ! u_hat is dimensionless (self%sol); convert to physical m/yr for the physics
    u_phys = velocity_scale * u_hat_values(1:2)
    call SSA_FEM_PETSc_sliding_beta( u_phys(1), u_phys(2), a_values( i_tauc), c_values, beta, dbeta_duabs)

    ! Weak form  integral( f1:grad(phi) + f0.phi ) = 0  with f1 = membrane stress M.
    ! Integrating div(M) by parts gives f0 = beta*u - tau_d, where tau_d = -rho g H grad(Hs)
    ! is the (downslope) driving stress as defined in UFEMISM. Divide by stress_scale to
    ! keep the dimensionless residual O(1).
    f0_values( 1) = (beta * u_phys( 1) - a_values( i_taudx)) / stress_scale
    f0_values( 2) = (beta * u_phys( 2) - a_values( i_taudy)) / stress_scale

  end subroutine SSA_FEM_PETSc_f0

  subroutine SSA_FEM_PETSc_f1( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f1) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f1
    real(c_double), value :: time
    real(c_double), pointer :: u_x_hat_values(:), a_values(:), c_values(:), f1_values(:)
    real(c_double)          :: u_x_phys(4), D(4), eps2, eta, N
    integer                 :: m

    call c_f_pointer( u_x, u_x_hat_values, [4])
    call c_f_pointer( a, a_values, [n_aux_comp])
    call c_f_pointer( constants, c_values, [int( nconstants)])
    call c_f_pointer( f1, f1_values, [4])

    ! grad(u_hat) is dimensionless; convert to a physical velocity gradient
    u_x_phys = velocity_scale * u_x_hat_values(1:4)
    call SSA_FEM_PETSc_strain( u_x_phys, c_values(1), c_values(2), a_values( i_Abar), D, eps2, eta)
    N = eta * a_values( i_H)

    do m = 1, 4
      f1_values( m) = 2._c_double * N * D( m) / stress_scale
    end do

  end subroutine SSA_FEM_PETSc_f1

  subroutine SSA_FEM_PETSc_g0( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, u_tshift, x, nconstants, constants, g0) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants
    real(c_double), value :: time, u_tshift
    real(c_double), intent(out) :: g0(*)
    real(c_double), pointer :: u_hat_values(:), a_values(:), c_values(:)
    real(c_double)          :: u_phys(2), beta, dbeta_duabs, uabs, rank1, jac_scale

    call c_f_pointer( u, u_hat_values, [2])
    call c_f_pointer( a, a_values, [n_aux_comp])
    call c_f_pointer( constants, c_values, [int( nconstants)])

    u_phys = velocity_scale * u_hat_values(1:2)
    call SSA_FEM_PETSc_sliding_beta( u_phys(1), u_phys(2), a_values( i_tauc), c_values, beta, dbeta_duabs)

    ! d(beta(|u|) u_c) / d u_c'  =  beta delta_{c c'}  +  (dbeta/d|u| / |u|) u_c u_c'
    ! (evaluated in physical units, then chain-ruled through u = velocity_scale*u_hat
    ! and divided by stress_scale, since g0 = d(f0*stress_scale)/d(u_hat) / stress_scale)
    uabs      = sqrt( c_values( ic_dv)**2 + u_phys(1)**2 + u_phys(2)**2)
    rank1     = dbeta_duabs / uabs
    jac_scale = velocity_scale / stress_scale
    g0( 1) = jac_scale * (beta + rank1 * u_phys(1) * u_phys(1))
    g0( 2) = jac_scale * (       rank1 * u_phys(1) * u_phys(2))
    g0( 3) = jac_scale * (       rank1 * u_phys(2) * u_phys(1))
    g0( 4) = jac_scale * (beta + rank1 * u_phys(2) * u_phys(2))

  end subroutine SSA_FEM_PETSc_g0

  subroutine SSA_FEM_PETSc_g3( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, u_tshift, x, nconstants, constants, g3) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants
    real(c_double), value :: time, u_tshift
    real(c_double), intent(out) :: g3(*)
    real(c_double), pointer :: u_x_hat_values(:), a_values(:), c_values(:)
    real(c_double)          :: u_x_phys(4), D(4), eps2, eta, N, n_glen, p, coef, jac_scale
    integer                 :: ci, cj, di, dj, m, k, idx

    call c_f_pointer( u_x, u_x_hat_values, [4])
    call c_f_pointer( a, a_values, [n_aux_comp])
    call c_f_pointer( constants, c_values, [int( nconstants)])

    u_x_phys = velocity_scale * u_x_hat_values(1:4)
    call SSA_FEM_PETSc_strain( u_x_phys, c_values(1), c_values(2), a_values( i_Abar), D, eps2, eta)
    N      = eta * a_values( i_H)
    n_glen = c_values(2)
    p      = (1._c_double - n_glen) / (2._c_double * n_glen)   ! d(ln eta) / d(ln eps2)
    coef   = 2._c_double * a_values( i_H) * eta * p / eps2     ! rank-1 shear-thinning weight
    jac_scale = velocity_scale / stress_scale

    ! g3[c,c',d,d'] = d f1[c,d] / d(du_c'/dx_d'),  index0 = ((c*2 + c')*2 + d)*2 + d'
    !   f1[c,d]     = 2 N D(m),  m = 2c + d + 1
    !   d D(m)/d u_x(k) = dD_dgradu(m,k),  k = 2c' + d' + 1
    !   d N / d u_x(k)  = coef/(2H) * D(k)   ->   d f1[m]/d u_x[k]
    !                   = 2 N dD_dgradu(m,k) + coef * D(m) * D(k)
    ! (all in physical units; jac_scale converts d(f1_phys)/d(grad u_phys) to
    ! d(f1_hat)/d(grad u_hat), since u_x_phys = velocity_scale * u_x_hat and
    ! f1_hat = f1_phys / stress_scale)
    do ci = 0, 1
      do di = 0, 1
        m = 2*ci + di + 1
        do cj = 0, 1
          do dj = 0, 1
            k   = 2*cj + dj + 1
            idx = ((ci*2 + cj)*2 + di)*2 + dj + 1
            g3( idx) = jac_scale * (2._c_double * N * dD_dgradu( m, k) + coef * D( m) * D( k))
          end do
        end do
      end do
    end do

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
