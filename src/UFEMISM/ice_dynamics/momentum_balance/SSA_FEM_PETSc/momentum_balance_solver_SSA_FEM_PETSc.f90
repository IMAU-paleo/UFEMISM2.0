module momentum_balance_solver_SSA_FEM_PETSc

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
  ! repository root):
  !
  !   Phase 0 (this file, current state)
  !     - the solver exists, is selectable, allocates/initialises/remaps/deallocates
  !       cleanly, and runs to completion, but does NOT solve a momentum balance yet:
  !       it leaves the ice velocity at zero. Its purpose is to verify the plumbing
  !       (dispatch, config, build, integrated-test run) before any PETSc code is added.
  !   Phase 1+  - DMPlex + PetscFE field + SNES skeleton, then the linear constant-
  !               viscosity SSA, then auxiliary fields, the full non-linear residual,
  !               the analytic Jacobian, boundary conditions, scaling and solver tuning.
  !
  ! The internal unknown will be a nodal (vertex, P1) velocity field; the result is
  ! exposed on the triangles as u_vav_b / v_vav_b, exactly like the existing SSA
  ! solver, so that ice-thickness evolution and all other downstream code are unaffected.

  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

  private

  public :: type_momentum_balance_solver_SSA_FEM_PETSc

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SSA_FEM_PETSc

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

  end type type_momentum_balance_solver_SSA_FEM_PETSc

contains

  ! == Main routines

  subroutine momentum_balance_solver_SSA_FEM_PETSc_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Phase 0: nothing to allocate yet.
    ! Phase 1+ : the persistent DMPlex, PetscFE, PetscDS and SNES handles, plus the
    !            b-grid solution fields u_vav_b / v_vav_b, will be allocated here.

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

    ! Phase 0: nothing to deallocate yet.
    ! Phase 1+ : SNESDestroy / PetscFEDestroy / DMDestroy (main + auxiliary) etc.

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

    if (par%primary) write(0,'(A)') '    NOTE: the SSA_FEM_PETSc solver is a Phase-0 placeholder; ' // &
      'it does not solve a momentum balance yet and leaves the ice velocity at zero.'

    ! Phase 0: nothing to initialise yet.
    ! Phase 1+ : build the DMPlex from self%mesh (petsc_dmplex/mesh_to_dmplex),
    !            attach the P1 velocity PetscFE field, create the PetscDS, and set up
    !            the SNES.

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_initialise

  subroutine momentum_balance_solver_SSA_FEM_PETSc_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the SSA with PetscFE / PetscSNES.
    !< Phase 0: no-op placeholder - the velocity is left at zero.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    class(atype_ice_model_data),                       intent(inout) :: ice
    class(atype_ice_geometry_model_data),              intent(in   ) :: geom
    type(type_bed_roughness_model),                    intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,               intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,               intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,               intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,               intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,               intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,               intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FEM_PETSc_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Phase 0: do not solve anything. Just report zero stability counters so the
    !          rest of the model behaves as for a converged solve.
    ! Phase 1+ : assemble the auxiliary fields, call SNESSolve, and copy the nodal
    !            solution back onto the b-grid.

    self%n_visc_its = 0
    self%n_Axb_its  = 0

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_run

  subroutine momentum_balance_solver_SSA_FEM_PETSc_set_velocities( self, ice, vel)
    !< Hand the solver result to the ice velocity model.
    !< Phase 0: write zeros, so downstream code (secondary velocities, ice-thickness
    !< evolution) has well-defined input.

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
      vel%u_3D_b( ti,:) = 0._dp
      vel%v_3D_b( ti,:) = 0._dp
    end do

    ! Strain rates on the vertices
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

    ! Phase 0: nothing to remap yet.
    ! Phase 1+ : destroy and rebuild the DMPlex / PetscFE / PetscDS / SNES for the
    !            new mesh, and remap u_vav_b / v_vav_b via the a-grid.

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FEM_PETSc_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'SSA_FEM_PETSc'
  end function get_momentum_balance_solver_name

  ! == Restart NetCDF files

  subroutine create_restart_file_old_SSA_FEM_PETSc( self)
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(inout) :: self
    ! Phase 0: no restart file. Phase 8 will reuse the SSA restart-file layout.
  end subroutine create_restart_file_old_SSA_FEM_PETSc

  subroutine write_to_restart_file_old_SSA_FEM_PETSc( self, time)
    class(type_momentum_balance_solver_SSA_FEM_PETSc), intent(in   ) :: self
    real(dp),                                          intent(in   ) :: time
    ! Phase 0: no restart file. Phase 8 will reuse the SSA restart-file layout.
  end subroutine write_to_restart_file_old_SSA_FEM_PETSc

end module momentum_balance_solver_SSA_FEM_PETSc
