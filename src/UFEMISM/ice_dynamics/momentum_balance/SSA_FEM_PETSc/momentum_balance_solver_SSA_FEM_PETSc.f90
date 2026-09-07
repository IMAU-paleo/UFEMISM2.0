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
  ! Implementation is staged (see SSA_PetscFE_SNES_implementation_plan.md in the
  ! repository root). Current state: Phase 1 - the full
  !   DMPlex -> PetscFE -> PetscDS -> PetscSNES
  ! pipeline is stood up and solves a CONSTANT-COEFFICIENT linear SSA:
  !
  !   residual = integral( f0 . phi + f1 : grad(phi) ) = 0
  !     f1[u,x] = 2 N (2 du/dx + dv/dy)      f1[u,y] = N (du/dy + dv/dx)
  !     f1[v,x] = N (du/dy + dv/dx)          f1[v,y] = 2 N (2 dv/dy + du/dx)
  !     f0[u]   = beta u + tau_dx            f0[v]   = beta v + tau_dy
  !
  ! with N = eta*H, beta and (tau_dx, tau_dy) frozen to uniform constants. The
  ! basal-drag term beta*I makes the operator positive-definite under natural
  ! (do-nothing) boundary conditions, so no essential BCs are imposed yet.
  !
  ! With constant coefficients, constant forcing and natural BCs the exact solution
  ! is the spatially uniform field  u = -tau_dx/beta,  v = -tau_dy/beta,  which lies
  ! in the P1 space; run() checks the computed nodal field against it.
  !
  ! Not done yet: real viscosity/friction/driving-stress fields (Phase 2), the
  ! non-linear Glen residual (Phase 3), the analytic shear-thinning Jacobian
  ! (Phase 4), boundary conditions (Phase 5), scaling and solver tuning (Phase 6).
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
    PetscDSSetConstants, DMCreateMatrix, DMCreateGlobalVector, DMCreateLocalVector, &
    DMGlobalToLocalBegin, DMGlobalToLocalEnd, DMGetLocalSection, DMGetLabel, DMDestroy, &
    DMPlexGetDepthStratum, DMLabelGetValue, PetscSectionGetOffset, &
    VecSet, VecDestroy, VecGetArrayRead, VecRestoreArrayRead, MatDestroy, INSERT_VALUES, &
    SNESCreate, SNESSetDM, SNESSetType, SNESSetTolerances, SNESGetKSP, SNESSolve, SNESDestroy, &
    SNESGetIterationNumber, SNESNEWTONLS, KSPSetType, KSPGetPC, KSPPREONLY, PCSetType, PCLU
  use mpi_f08, only: MPI_ALLTOALL, MPI_ALLTOALLV, MPI_ALLREDUCE, MPI_COMM_WORLD, MPI_IN_PLACE, &
    MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_MAX
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use bed_roughness_model_types, only: type_bed_roughness_model
  use reallocate_mod, only: reallocate_bounds
  use petsc_dmplex, only: mesh_to_dmplex, dmplex_upsy_vertex_id_label_name
  use mesh_disc_apply_operators, only: map_a_b_2D

  implicit none

  private

  public :: type_momentum_balance_solver_SSA_FEM_PETSc

  ! Layout of the PetscDS constants array (see the weak form in the module header)
  integer, parameter :: i_N     = 1   ! N = eta * H          [Pa yr m]
  integer, parameter :: i_beta  = 2   ! basal friction coeff [Pa yr m^-1]
  integer, parameter :: i_taudx = 3   ! driving stress, x    [Pa]
  integer, parameter :: i_taudy = 4   ! driving stress, y    [Pa]
  integer, parameter :: n_petsc_constants = 4

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SSA_FEM_PETSc

    ! Persistent PETSc objects (rebuilt on remap)
    type(tDM)      :: dm
    type(tPetscFE) :: fe
    type(tSNES)    :: snes
    type(tMat)     :: jac
    logical        :: petsc_is_built = .false.

    ! Frozen constant coefficients for the Phase 1 linear SSA
    real(dp), dimension(n_petsc_constants) :: petsc_constants = 0._dp

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

    if (par%primary) write(0,'(A)') '    NOTE: the SSA_FEM_PETSc solver is at Phase 1 - it solves a ' // &
      'constant-coefficient linear SSA, not yet the real momentum balance.'

    ! Phase 1 placeholder coefficients. Chosen so the closed-form check has a
    ! non-trivial answer: u = -tau_dx/beta = -0.2 m/yr, v = -tau_dy/beta = 0.1 m/yr.
    ! Phase 2 replaces these with per-vertex auxiliary fields (H, eta, beta, grad s).
    self%petsc_constants( i_N)     = 1.0e14_dp
    self%petsc_constants( i_beta)  = 1.0e5_dp
    self%petsc_constants( i_taudx) = 2.0e4_dp
    self%petsc_constants( i_taudy) = -1.0e4_dp

    call self%build_petsc_objects()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_initialise

  subroutine momentum_balance_solver_SSA_FEM_PETSc_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the constant-coefficient linear SSA with PetscFE / PetscSNES.

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
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_run'
    type(tVec)                  :: solution
    integer                     :: ierr, snes_its
    integer(c_int)              :: snes_reason
    real(dp)                    :: u_exact, v_exact, u_dev, v_dev, tol

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Solve
    PetscCall( DMCreateGlobalVector( self%dm, solution, ierr))
    PetscCall( VecSet( solution, 0._dp, ierr))
    PetscCall( SNESSolve( self%snes, PETSC_NULL_VEC, solution, ierr))
    PetscCall( SNESGetIterationNumber( self%snes, snes_its, ierr))
    ierr = snes_get_converged_reason( self%snes%v, snes_reason)
    CHKERRQ( ierr)
    if (snes_reason < 0) then
      call crash('SSA_FEM_PETSc: SNES diverged (SNESConvergedReason = {int_01})', int_01 = int( snes_reason))
    end if

    ! Copy the PETSc solution onto the mesh vertices, then map to the triangles
    call copy_PETSc_solution_to_mesh_vertices_vec2( self%dm, solution, self%mesh, self%u_vav_a, self%v_vav_a)
    PetscCall( VecDestroy( solution, ierr))

    call map_a_b_2D( self%mesh, self%u_vav_a, self%u_vav_b)
    call map_a_b_2D( self%mesh, self%v_vav_a, self%v_vav_b)

    ! Phase 1 closed-form check: the exact solution is the uniform field -tau/beta
    u_exact = -self%petsc_constants( i_taudx) / self%petsc_constants( i_beta)
    v_exact = -self%petsc_constants( i_taudy) / self%petsc_constants( i_beta)
    u_dev = maxval( abs( self%u_vav_a - u_exact))
    v_dev = maxval( abs( self%v_vav_a - v_exact))
    call MPI_ALLREDUCE( MPI_IN_PLACE, u_dev, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, v_dev, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (par%primary) write(0,'(A,I0,A,ES10.3,A,ES10.3)') '    SSA_FEM_PETSc Phase-1 solve: SNES its = ', &
      snes_its, ', max|u-u_exact| = ', u_dev, ', max|v-v_exact| = ', v_dev
    tol = max( 1.0e-8_dp, 1.0e-6_dp * max( abs( u_exact), abs( v_exact)))
    if (u_dev > tol .or. v_dev > tol) then
      call warning('SSA_FEM_PETSc: Phase-1 closed-form check failed - the constant-coefficient ' // &
        'SSA solution is not the expected uniform field -tau/beta.')
    end if

    self%n_visc_its = snes_its
    self%n_Axb_its  = 0

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

    ! Phase 1: rebuild everything from scratch on the new mesh (self%mesh has already
    ! been repointed to mesh_new by remap_model). Velocities are reset to zero; a
    ! proper remap of u_vav via the a-grid follows in a later phase.
    if (self%petsc_is_built) call self%destroy_petsc_objects()

    call reallocate_bounds( self%u_vav_a, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%v_vav_a, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%u_vav_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%v_vav_b, mesh_new%ti1, mesh_new%ti2)

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
    !< Build the DMPlex, the P1 vector PetscFE field, the PetscDS weak form, and the SNES.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_petsc_objects'
    type(tPetscDS)              :: ds
    type(tPetscObject)          :: fe_object
    type(tKSP)                  :: ksp
    type(tPC)                   :: pc
    integer                     :: ierr
    integer(c_intptr_t)         :: no_context

    call init_routine( routine_name)

    no_context = 0_c_intptr_t

    ! DMPlex from the UFEMISM mesh (creates + distributes; preserves UPSY vertex IDs
    ! in the 'upsy_vertex_id' DMLabel)
    call mesh_to_dmplex( self%mesh, self%dm)

    ! One 2-component P1 Lagrange field: the vertically averaged horizontal velocity
    PetscCall( PetscFECreateLagrange( PETSC_COMM_SELF, 2, 2, PETSC_TRUE, 1, -1, self%fe, ierr))
    PetscCall( PetscObjectSetName( self%fe, 'velocity', ierr))
    PetscObjectSpecificCast( fe_object, self%fe)
    PetscCall( DMSetField( self%dm, 0, PETSC_NULL_DMLABEL, fe_object, ierr))
    PetscCall( DMCreateDS( self%dm, ierr))
    PetscCall( DMGetDS( self%dm, ds, ierr))

    ! Weak form: constants, residual (f0, f1) and analytic Jacobian (g0, g3)
    PetscCall( PetscDSSetConstants( ds, n_petsc_constants, self%petsc_constants, ierr))
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

    PetscCall( MatDestroy( self%jac, ierr))
    PetscCall( SNESDestroy( self%snes, ierr))
    PetscCall( PetscFEDestroy( self%fe, ierr))
    PetscCall( DMDestroy( self%dm, ierr))
    self%petsc_is_built = .false.

    call finalise_routine( routine_name)

  end subroutine destroy_petsc_objects

  ! ===== Solution transfer: PETSc global Vec -> UFEMISM vertex arrays =====

  subroutine copy_PETSc_solution_to_mesh_vertices_vec2( dm, solution, mesh, u_a, v_a)
    !< Scatter a 2-component nodal PETSc solution back onto the UFEMISM vertex
    !< distribution (vi1:vi2), using the 'upsy_vertex_id' DMLabel and
    !< mesh%V_owning_process. Generalisation of the scalar routine in
    !< ct_PETSc_SNES_Poisson.f90.

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

    ! Count how many local DMPlex vertices belong to each UFEMISM process
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

    ! Pack (vertex ID, u, v) for each local DMPlex vertex
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

    ! Exchange
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

    ! Average over duplicate copies (there should be exactly one per vertex at overlap 0)
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

  ! ===== PetscDS pointwise weak-form functions (bind(C)) =====
  !
  ! PETSc assembles  residual = integral( f0 . phi + f1 : grad(phi) ).
  ! Field 0 is the 2-vector velocity (u, v); dim = 2, Nc = 2.
  ! Gradient layout u_x[c*dim + d]:  u_x(1)=du/dx u_x(2)=du/dy u_x(3)=dv/dx u_x(4)=dv/dy

  subroutine SSA_FEM_PETSc_f0( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f0) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f0
    real(c_double), value :: time
    real(c_double), pointer :: u_values(:), c_values(:), f0_values(:)

    call c_f_pointer( u, u_values, [2])
    call c_f_pointer( constants, c_values, [int( nconstants)])
    call c_f_pointer( f0, f0_values, [2])

    ! f0 = beta * u + tau_d   (basal drag minus RHS driving stress, moved to the LHS)
    f0_values( 1) = c_values( i_beta) * u_values( 1) + c_values( i_taudx)
    f0_values( 2) = c_values( i_beta) * u_values( 2) + c_values( i_taudy)

  end subroutine SSA_FEM_PETSc_f0

  subroutine SSA_FEM_PETSc_f1( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f1) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f1
    real(c_double), value :: time
    real(c_double), pointer :: u_x_values(:), c_values(:), f1_values(:)
    real(c_double)          :: N, du_dx, du_dy, dv_dx, dv_dy

    call c_f_pointer( u_x, u_x_values, [4])
    call c_f_pointer( constants, c_values, [int( nconstants)])
    call c_f_pointer( f1, f1_values, [4])

    N     = c_values( i_N)
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
    real(c_double), pointer :: c_values(:)

    call c_f_pointer( constants, c_values, [int( nconstants)])

    ! d f0_c / d u_c'  =  beta * delta_{c c'}   (2x2, row-major)
    g0( 1:4) = 0._c_double
    g0( 1)   = c_values( i_beta)
    g0( 4)   = c_values( i_beta)

  end subroutine SSA_FEM_PETSc_g0

  subroutine SSA_FEM_PETSc_g3( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, u_tshift, x, nconstants, constants, g3) bind(C)

    integer(c_intptr_t), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants
    real(c_double), value :: time, u_tshift
    real(c_double), intent(out) :: g3(*)
    real(c_double), pointer :: c_values(:)
    real(c_double)          :: N

    call c_f_pointer( constants, c_values, [int( nconstants)])
    N = c_values( i_N)

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
