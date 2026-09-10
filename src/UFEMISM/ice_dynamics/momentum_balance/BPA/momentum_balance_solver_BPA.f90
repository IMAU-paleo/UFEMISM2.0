module momentum_balance_solver_BPA

  ! Routines for calculating ice velocities using the Blatter-Pattyn Approximation (BPA)

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_LOR, MPI_LOGICAL, MPI_MIN, MPI_MAX, MPI_SUM
  use precisions, only: dp
  use UPSY_main, only: UPSY
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: warning, crash, happy, init_routine, finalise_routine
  use model_configuration, only: C
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use parameters
  use mesh_disc_apply_operators, only: map_a_b_2D, map_a_b_3D, ddx_a_b_2D, ddy_a_b_2D, &
    ddx_b_a_3D, ddy_b_a_3D, calc_3D_gradient_bk_ak, calc_3D_gradient_bk_bks, &
    map_ak_bks, map_bks_ak, calc_3D_gradient_ak_bk, calc_3D_gradient_bks_bk, map_b_a_3D, &
    map_b_a_2D
  use mesh_disc_calc_matrix_operators_3D, only: calc_3D_matrix_operators_mesh
  use mesh_zeta, only: vertical_average
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_utilities, only: find_ti_copy_ISMIP_HOM_periodic
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use netcdf_io_main
  use mpi_distributed_memory, only: gather_to_all
  use constitutive_equation, only: calc_effective_viscosity_Glen_3D_uv_only, calc_ice_rheology_Glen
  use zeta_gradients, only: calc_zeta_gradients
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  use bed_roughness_model_types, only: type_bed_roughness_model
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver

  implicit none

  private

  public :: type_momentum_balance_solver_BPA

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_BPA

      ! Solution
      real(dp), dimension(:,:), allocatable :: u_bk                        ! 3-D horizontal ice velocity [m yr^-1]
      real(dp), dimension(:,:), allocatable :: v_bk

      ! Intermediate data fields
      real(dp), dimension(:,:), allocatable :: du_dx_ak                    ! Strain rates [yr^-1]
      real(dp), dimension(:,:), allocatable :: du_dy_ak
      real(dp), dimension(:,:), allocatable :: du_dz_ak
      real(dp), dimension(:,:), allocatable :: dv_dx_ak
      real(dp), dimension(:,:), allocatable :: dv_dy_ak
      real(dp), dimension(:,:), allocatable :: dv_dz_ak
      real(dp), dimension(:,:), allocatable :: du_dx_bks
      real(dp), dimension(:,:), allocatable :: du_dy_bks
      real(dp), dimension(:,:), allocatable :: du_dz_bks
      real(dp), dimension(:,:), allocatable :: dv_dx_bks
      real(dp), dimension(:,:), allocatable :: dv_dy_bks
      real(dp), dimension(:,:), allocatable :: dv_dz_bks
      real(dp), dimension(:,:), allocatable :: eta_ak                      ! Effective viscosity
      real(dp), dimension(:,:), allocatable :: eta_bks
      real(dp), dimension(:,:), allocatable :: eta_bk
      real(dp), dimension(:,:), allocatable :: deta_dx_bk                  ! Gradients of eta
      real(dp), dimension(:,:), allocatable :: deta_dy_bk
      real(dp), dimension(:,:), allocatable :: deta_dz_bk
      real(dp), dimension(:  ), allocatable :: basal_friction_coefficient_b! Friction coefficient (tau_b = u * beta_b)
      real(dp), dimension(:  ), allocatable :: dh_dx_b                     ! Surface slope
      real(dp), dimension(:  ), allocatable :: dh_dy_b
      real(dp), dimension(:  ), allocatable :: db_dx_b                     ! Basal slope
      real(dp), dimension(:  ), allocatable :: db_dy_b
      real(dp), dimension(:  ), allocatable :: tau_dx_b                    ! Driving stress
      real(dp), dimension(:  ), allocatable :: tau_dy_b
      real(dp), dimension(:,:), allocatable :: u_bk_prev                   ! Previous velocity solution
      real(dp), dimension(:,:), allocatable :: v_bk_prev

      ! Restart file
      character(len=256)                    :: restart_filename

    contains

      ! Procedures for model memory management and operation
      procedure, public :: get_momentum_balance_solver_name
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_BPA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_BPA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_BPA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_BPA_run
      procedure, public :: set_velocities_to_solver_results   => momentum_balance_solver_BPA_set_velocities
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_BPA_remap

      procedure, public :: create_restart_file_old            => create_restart_file_old_BPA
      procedure, public :: write_to_restart_file_old          => write_to_restart_file_old_BPA

      procedure, public :: initialise_BPA_velocities_from_file

      procedure, public :: solve_BPA_linearised
      procedure, public :: assemble_BPA_linearised_matrix_eq
      procedure, public :: calc_BPA_stiffness_matrix_row_free
      procedure, public :: calc_BPA_stiffness_matrix_row_BC_surf
      procedure, public :: calc_BPA_stiffness_matrix_row_BC_base
      procedure, public :: calc_BPA_stiffness_matrix_row_BC_west
      procedure, public :: calc_BPA_stiffness_matrix_row_BC_east
      procedure, public :: calc_BPA_stiffness_matrix_row_BC_south
      procedure, public :: calc_BPA_stiffness_matrix_row_BC_north

      procedure, public :: calc_driving_stress
      procedure, public :: calc_strain_rates
      procedure, public :: calc_effective_viscosity
      procedure, public :: calc_applied_basal_friction_coefficient
      procedure, public :: relax_viscosity_iterations
      procedure, public :: calc_visc_iter_UV_resid
      procedure, public :: apply_velocity_limits

  end type type_momentum_balance_solver_BPA

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine solve_BPA_linearised( self, ice, n_Axb_its, &
      BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
      class(type_momentum_balance_solver_BPA),                         intent(inout) :: self
      class(atype_ice_model_data),                                     intent(in   ) :: ice
      integer,                                                         intent(  out) :: n_Axb_its              ! Number of iterations used in the iterative solver
      integer,  dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_mask_bk      ! Mask of triangles where velocity is prescribed
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_u_bk         ! Prescribed velocities in the x-direction
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_v_bk         ! Prescribed velocities in the y-direction
    end subroutine solve_BPA_linearised

    module subroutine assemble_BPA_linearised_matrix_eq( self, ice, &
      BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk, &
      A_CSR, bb, uv_bkuv)
      class(type_momentum_balance_solver_BPA),                         intent(in   ) :: self
      class(atype_ice_model_data),                                     intent(in   ) :: ice
      integer,  dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_mask_bk      ! Mask of triangles where velocity is prescribed
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_u_bk         ! Prescribed velocities in the x-direction
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), intent(in   ) :: BC_prescr_v_bk         ! Prescribed velocities in the y-direction
      type(type_CSR_matrix_dp),                                        intent(inout) :: A_CSR
      real(dp), dimension(:), allocatable,                             intent(inout) :: bb
      real(dp), dimension(:), allocatable,                             intent(inout) :: uv_bkuv
    end subroutine assemble_BPA_linearised_matrix_eq

    module subroutine calc_BPA_stiffness_matrix_row_free( self, A_CSR, bb, row_tikuv)
      class(type_momentum_balance_solver_BPA), intent(in   ) :: self
      type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
      integer,                                 intent(in   ) :: row_tikuv
    end subroutine calc_BPA_stiffness_matrix_row_free

    module subroutine calc_BPA_stiffness_matrix_row_BC_surf( self, ice, A_CSR, bb, row_tikuv)
      class(type_momentum_balance_solver_BPA), intent(in   ) :: self
      class(atype_ice_model_data),             intent(in   ) :: ice
      type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
      integer,                                 intent(in   ) :: row_tikuv
    end subroutine calc_BPA_stiffness_matrix_row_BC_surf

    module subroutine calc_BPA_stiffness_matrix_row_BC_base( self, ice, A_CSR, bb, row_tikuv)
      class(type_momentum_balance_solver_BPA), intent(in   ) :: self
      class(atype_ice_model_data),             intent(in   ) :: ice
      type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
      integer,                                 intent(in   ) :: row_tikuv
    end subroutine calc_BPA_stiffness_matrix_row_BC_base

    module subroutine calc_BPA_stiffness_matrix_row_BC_west( self, A_CSR, bb, row_tikuv)
      class(type_momentum_balance_solver_BPA), intent(in   ) :: self
      type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
      integer,                                 intent(in   ) :: row_tikuv
    end subroutine calc_BPA_stiffness_matrix_row_BC_west

    module subroutine calc_BPA_stiffness_matrix_row_BC_east( self, A_CSR, bb, row_tikuv)
      class(type_momentum_balance_solver_BPA), intent(in   ) :: self
      type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
      integer,                                 intent(in   ) :: row_tikuv
    end subroutine calc_BPA_stiffness_matrix_row_BC_east

    module subroutine calc_BPA_stiffness_matrix_row_BC_south( self, A_CSR, bb, row_tikuv)
      class(type_momentum_balance_solver_BPA), intent(in   ) :: self
      type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
      integer,                                 intent(in   ) :: row_tikuv
    end subroutine calc_BPA_stiffness_matrix_row_BC_south

    module subroutine calc_BPA_stiffness_matrix_row_BC_north( self, A_CSR, bb, row_tikuv)
      class(type_momentum_balance_solver_BPA), intent(in   ) :: self
      type(type_CSR_matrix_dp),                intent(inout) :: A_CSR
      real(dp), dimension(A_CSR%i1:A_CSR%i2),  intent(inout) :: bb
      integer,                                 intent(in   ) :: row_tikuv
    end subroutine calc_BPA_stiffness_matrix_row_BC_north

  end interface

contains

! == Main routines

  subroutine momentum_balance_solver_BPA_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_BPA_allocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( self%u_bk(                         self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%v_bk(                         self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz  ), source = 0._dp)

    ! Intermediate data fields
    allocate( self%du_dx_ak(                     self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%du_dy_ak(                     self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%du_dz_ak(                     self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%dv_dx_ak(                     self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%dv_dy_ak(                     self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%dv_dz_ak(                     self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%du_dx_bks(                    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz-1), source = 0._dp)
    allocate( self%du_dy_bks(                    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz-1), source = 0._dp)
    allocate( self%du_dz_bks(                    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz-1), source = 0._dp)
    allocate( self%dv_dx_bks(                    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz-1), source = 0._dp)
    allocate( self%dv_dy_bks(                    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz-1), source = 0._dp)
    allocate( self%dv_dz_bks(                    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz-1), source = 0._dp)
    allocate( self%eta_ak(                       self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%eta_bks(                      self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz-1), source = 0._dp)
    allocate( self%eta_bk(                       self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%deta_dx_bk(                   self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%deta_dy_bk(                   self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%deta_dz_bk(                   self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz  ), source = 0._dp)
    allocate( self%basal_friction_coefficient_b( self%mesh%ti1:self%mesh%ti2                 ), source = 0._dp)
    allocate( self%dh_dx_b(                      self%mesh%ti1:self%mesh%ti2                 ), source = 0._dp)
    allocate( self%dh_dy_b(                      self%mesh%ti1:self%mesh%ti2                 ), source = 0._dp)
    allocate( self%db_dx_b(                      self%mesh%ti1:self%mesh%ti2                 ), source = 0._dp)
    allocate( self%db_dy_b(                      self%mesh%ti1:self%mesh%ti2                 ), source = 0._dp)
    allocate( self%tau_dx_b(                     self%mesh%ti1:self%mesh%ti2                 ), source = 0._dp)
    allocate( self%tau_dy_b(                     self%mesh%ti1:self%mesh%ti2                 ), source = 0._dp)
    allocate( self%u_bk_prev(                    self%mesh%nTri,1:self%mesh%nz               ), source = 0._dp)
    allocate( self%v_bk_prev(                    self%mesh%nTri,1:self%mesh%nz               ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_BPA_allocate

  subroutine momentum_balance_solver_BPA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_BPA_deallocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    deallocate( self%u_bk)
    deallocate( self%v_bk)

    ! Intermediate data fields
    deallocate( self%du_dx_ak)
    deallocate( self%du_dy_ak)
    deallocate( self%du_dz_ak)
    deallocate( self%dv_dx_ak)
    deallocate( self%dv_dy_ak)
    deallocate( self%dv_dz_ak)
    deallocate( self%du_dx_bks)
    deallocate( self%du_dy_bks)
    deallocate( self%du_dz_bks)
    deallocate( self%dv_dx_bks)
    deallocate( self%dv_dy_bks)
    deallocate( self%dv_dz_bks)
    deallocate( self%eta_ak)
    deallocate( self%eta_bks)
    deallocate( self%eta_bk)
    deallocate( self%deta_dx_bk)
    deallocate( self%deta_dy_bk)
    deallocate( self%deta_dz_bk)
    deallocate( self%basal_friction_coefficient_b)
    deallocate( self%dh_dx_b)
    deallocate( self%dh_dy_b)
    deallocate( self%db_dx_b)
    deallocate( self%db_dy_b)
    deallocate( self%tau_dx_b)
    deallocate( self%tau_dy_b)
    deallocate( self%u_bk_prev)
    deallocate( self%v_bk_prev)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_BPA_deallocate

  subroutine momentum_balance_solver_BPA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'momentum_balance_solver_BPA_initialise'
    character(len=:), allocatable :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the choice of initial velocities for this model region
    select case (self%region_name())
    case default
      call crash('unknown model region "' // self%region_name() // '"!')
    case ('NAM')
      choice_initial_velocity  = C%choice_initial_velocity_NAM
    case ('EAS')
      choice_initial_velocity  = C%choice_initial_velocity_EAS
    case ('GRL')
      choice_initial_velocity  = C%choice_initial_velocity_GRL
    case ('ANT')
      choice_initial_velocity  = C%choice_initial_velocity_ANT
    end select

    ! Initialise velocities according to the specified method
    select case (choice_initial_velocity)
    case ('zero')
      self%u_bk( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
      self%v_bk( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
    case ('read_from_file')
      call self%initialise_BPA_velocities_from_file()
    case default
      call crash('unknown choice_initial_velocity "' // trim( choice_initial_velocity) // '"!')
    end select

    ! Set tolerances for PETSc matrix solver for the linearised BPA
    self%PETSc_rtol   = C%stress_balance_PETSc_rtol
    self%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_BPA_initialise

  subroutine initialise_BPA_velocities_from_file( self)
    !< Initialise the velocities for the BPA solver from an external NetCDF file

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter                                     :: routine_name = 'initialise_BPA_velocities_from_file'
    character(len=:), allocatable                                   :: filename
    real(dp)                                                        :: timeframe
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz) :: d_loc_3D_b

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the filename and timeframe to read for this model region
    select case (self%region_name())
    case default
      call crash('unknown model region "' // self%region_name() // '"!')
    case ('NAM')
      filename  = C%filename_initial_velocity_NAM
      timeframe = C%timeframe_initial_velocity_NAM
    case ('EAS')
      filename  = C%filename_initial_velocity_EAS
      timeframe = C%timeframe_initial_velocity_EAS
    case ('GRL')
      filename  = C%filename_initial_velocity_GRL
      timeframe = C%timeframe_initial_velocity_GRL
    case ('ANT')
      filename  = C%filename_initial_velocity_ANT
      timeframe = C%timeframe_initial_velocity_ANT
    end select

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    if (index( filename,'_LAST.nc') > 1) then
      call find_last_output_file( filename)
      call find_last_timeframe(   filename, timeframe)
    end if

    ! Write to terminal
    if (par%primary) write(0,*) '   Initialising BPA velocities from file "' // &
      UPSY%stru%colour_string( trim( filename),'light blue') // '"...'

    ! Read velocities from the file
    call read_field_from_mesh_file_dp_3D_b( filename, 'u_bk', d_loc_3D_b, time_to_read = timeframe)
    self%u_bk( self%mesh%ti1:self%mesh%ti2,:) = d_loc_3D_b
    call read_field_from_mesh_file_dp_3D_b( filename, 'v_bk', d_loc_3D_b, time_to_read = timeframe)
    self%v_bk( self%mesh%ti1:self%mesh%ti2,:) = d_loc_3D_b

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BPA_velocities_from_file

  subroutine momentum_balance_solver_BPA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the Blatter-Pattyn Approximation

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom
    type(type_bed_roughness_model),          intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,      intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,      intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,      intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,      intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,      intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,      intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'momentum_balance_solver_BPA_run'
    integer                               :: ierr
    logical                               :: grounded_ice_exists
    integer,  dimension(:,:), allocatable :: BC_prescr_mask_bk_applied
    real(dp), dimension(:,:), allocatable :: BC_prescr_u_bk_applied
    real(dp), dimension(:,:), allocatable :: BC_prescr_v_bk_applied
    integer                               :: viscosity_iteration_i
    logical                               :: has_converged
    real(dp)                              :: resid_UV, resid_UV_prev
    real(dp)                              :: uv_min, uv_max
    real(dp)                              :: visc_it_relax_applied
    real(dp)                              :: Glens_flow_law_epsilon_sq_0_applied
    integer                               :: nit_diverg_consec
    integer                               :: n_Axb_its_visc_it

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, or no sliding, no need to solve the BPA
    grounded_ice_exists = any( geom%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists) then
      self%u_bk( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
      self%v_bk( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions
    allocate( BC_prescr_mask_bk_applied( self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz))
    allocate( BC_prescr_u_bk_applied(    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz))
    allocate( BC_prescr_v_bk_applied(    self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz))
    if (present( BC_prescr_mask_bk) .or. present( BC_prescr_u_bk) .or. present( BC_prescr_v_bk)) then
      ! Safety
      if (.not. (present( BC_prescr_mask_bk) .and. present( BC_prescr_u_bk) .and. present( BC_prescr_v_bk))) then
        call crash('need to provide prescribed u,v fields and mask!')
      end if
      BC_prescr_mask_bk_applied = BC_prescr_mask_bk
      BC_prescr_u_bk_applied    = BC_prescr_u_bk
      BC_prescr_v_bk_applied    = BC_prescr_v_bk
    else
      BC_prescr_mask_bk_applied = 0
      BC_prescr_u_bk_applied    = 0._dp
      BC_prescr_v_bk_applied    = 0._dp
    end if

    ! Calculate zeta gradients
    call calc_zeta_gradients( self%mesh, ice, geom)

    ! Calculate 3-D matrix operators for the current ice geometry
    call calc_3D_matrix_operators_mesh( self%mesh, &
      ice%dzeta_dx_ak, ice%dzeta_dy_ak, ice%dzeta_dx_bk, ice%dzeta_dy_bk, &
      ice%dzeta_dz_bk, ice%dzeta_dz_bks, &
      ice%d2zeta_dx2_bk, ice%d2zeta_dxdy_bk, ice%d2zeta_dy2_bk)

    ! Calculate the driving stress
    call self%calc_driving_stress( geom)

    ! Adaptive relaxation parameter for the viscosity iteration
    resid_UV                            = 1E9_dp
    nit_diverg_consec                   = 0
    visc_it_relax_applied               = C%visc_it_relax
    Glens_flow_law_epsilon_sq_0_applied = C%Glens_flow_law_epsilon_sq_0

    ! Initialise stability info
    self%n_visc_its = 0
    self%n_Axb_its  = 0

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.
    viscosity_iteration: do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the strain rates for the current velocity solution
      call self%calc_strain_rates()

      ! Calculate the effective viscosity for the current velocity solution
      call self%calc_effective_viscosity( ice, geom, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      call self%calc_applied_basal_friction_coefficient( ice, geom, bed_roughness)

      ! Solve the linearised BPA to calculate a new velocity solution
      call self%solve_BPA_linearised( ice, n_Axb_its_visc_it, &
        BC_prescr_mask_bk_applied, BC_prescr_u_bk_applied, BC_prescr_v_bk_applied)

      ! Update stability info
      self%n_Axb_its = self%n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call self%apply_velocity_limits()

      ! Reduce the change between velocity solutions
      call self%relax_viscosity_iterations( visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      call self%calc_visc_iter_UV_resid( resid_UV)

      ! if the viscosity iteration diverges, lower the relaxation parameter
      if (resid_UV > resid_UV_prev) then
        nit_diverg_consec = nit_diverg_consec + 1
      else
        nit_diverg_consec = 0
      end if
      if (nit_diverg_consec > 2) then
        nit_diverg_consec = 0
        visc_it_relax_applied               = visc_it_relax_applied               * 0.9_dp
        Glens_flow_law_epsilon_sq_0_applied = Glens_flow_law_epsilon_sq_0_applied * 1.2_dp
      end if
      if (visc_it_relax_applied <= 0.05_dp .or. Glens_flow_law_epsilon_sq_0_applied >= 1E-5_dp) then
        if (visc_it_relax_applied < 0.05_dp) then
          call crash('viscosity iteration still diverges even with very low relaxation factor!')
        elseif (Glens_flow_law_epsilon_sq_0_applied > 1E-5_dp) then
          call crash('viscosity iteration still diverges even with very high effective strain rate regularisation!')
        end if
      end if

      ! DENK DROM
      uv_min = minval( self%u_bk)
      uv_max = maxval( self%u_bk)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      ! if (par%primary) WRITE(0,*) '    BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (resid_UV < C%visc_it_norm_dUV_tol) then
        has_converged = .true.
      end if

      ! if we've reached the maximum allowed number of iterations without converging, throw a warning
      if (viscosity_iteration_i > C%visc_it_nit) then
        if (par%primary) call warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        exit viscosity_iteration
      end if

    end do viscosity_iteration

    ! Stability info
    self%n_visc_its = viscosity_iteration_i

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_BPA_run

  subroutine momentum_balance_solver_BPA_set_velocities( self, ice, vel)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_velocity_model),         intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_BPA_set_velocities'
    integer                     :: ti, vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Velocities
    do ti = self%mesh%ti1, self%mesh%ti2
      vel%u_3D_b( ti,:) = self%u_bk( ti,:)
      vel%v_3D_b( ti,:) = self%v_bk( ti,:)
    end do

    ! Strain rates
    do vi = self%mesh%vi1, self%mesh%vi2
      vel%du_dx_3D( vi,:) = self%du_dx_ak( vi,:)
      vel%du_dy_3D( vi,:) = self%du_dy_ak( vi,:)
      vel%du_dz_3D( vi,:) = self%du_dz_ak( vi,:)
      vel%dv_dx_3D( vi,:) = self%dv_dx_ak( vi,:)
      vel%dv_dy_3D( vi,:) = self%dv_dy_ak( vi,:)
      vel%dv_dz_3D( vi,:) = self%dv_dz_ak( vi,:)
    end do

    ! In the BPA, gradients of w are neglected
    vel%dw_dx_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_BPA_set_velocities

  subroutine momentum_balance_solver_BPA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self
    type(type_mesh),                         intent(in   ) :: mesh_old
    type(type_mesh), target,                 intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'momentum_balance_solver_BPA_remap'
    real(dp), dimension(:,:), allocatable :: u_ak
    real(dp), dimension(:,:), allocatable :: v_ak

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_ak( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( v_ak( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map velocities from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_3D( mesh_old, self%u_bk, u_ak)
    call map_b_a_3D( mesh_old, self%v_bk, v_ak)

    ! Remap velocities from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, u_ak, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, v_ak, '2nd_order_conservative')

    ! reallocate memory for the velocities on the triangles
    call reallocate_bounds( self%u_bk, mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%v_bk, mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map velocities from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_3D( mesh_new, u_ak, self%u_bk)
    call map_a_b_3D( mesh_new, v_ak, self%v_bk)

    ! reallocate everything else
    ! ==========================

    call reallocate_bounds( self%du_dx_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( self%du_dy_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( self%du_dz_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( self%dv_dx_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( self%dv_dy_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( self%dv_dz_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( self%du_dx_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( self%du_dy_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( self%du_dz_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( self%dv_dx_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( self%dv_dy_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( self%dv_dz_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( self%eta_ak                      , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! Effective viscosity
    call reallocate_bounds( self%eta_bks                     , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( self%eta_bk                      , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )
    call reallocate_bounds( self%deta_dx_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    call reallocate_bounds( self%deta_dy_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    call reallocate_bounds( self%deta_dz_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    call reallocate_bounds( self%basal_friction_coefficient_b, mesh_new%ti1 , mesh_new%ti2               )       ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    call reallocate_bounds( self%dh_dx_b                     , mesh_new%ti1 , mesh_new%ti2               )       ! Surface slope
    call reallocate_bounds( self%dh_dy_b                     , mesh_new%ti1 , mesh_new%ti2               )
    call reallocate_bounds( self%db_dx_b                     , mesh_new%ti1 , mesh_new%ti2               )       ! Basal slope
    call reallocate_bounds( self%db_dy_b                     , mesh_new%ti1 , mesh_new%ti2               )
    call reallocate_bounds( self%tau_dx_b                    , mesh_new%ti1 , mesh_new%ti2               )       ! Driving stress
    call reallocate_bounds( self%tau_dy_b                    , mesh_new%ti1 , mesh_new%ti2               )
    call reallocate_clean ( self%u_bk_prev                   , mesh_new%nTri              , mesh_new%nz  )       ! Velocity solution from previous viscosity iteration
    call reallocate_clean ( self%v_bk_prev                   , mesh_new%nTri              , mesh_new%nz  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_BPA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_BPA), intent(in) :: self
    character(len=:), allocatable                  :: model_name
    model_name = 'BPA'
  end function get_momentum_balance_solver_name

! == Calculate several intermediate terms in the BPA

  subroutine calc_driving_stress( self, geom)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_driving_stress'
    integer                     :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate dh/dx, dh/dy, db/dx, db/dy on the b-grid
    call ddx_a_b_2D( self%mesh, geom%Hs , self%dh_dx_b)
    call ddy_a_b_2D( self%mesh, geom%Hs , self%dh_dy_b)
    call ddx_a_b_2D( self%mesh, geom%Hib, self%db_dx_b)
    call ddy_a_b_2D( self%mesh, geom%Hib, self%db_dy_b)

    ! Calculate the driving stress
    do ti = self%mesh%ti1, self%mesh%ti2
      self%tau_dx_b( ti) = -ice_density * grav * self%dh_dx_b( ti)
      self%tau_dy_b( ti) = -ice_density * grav * self%dh_dy_b( ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress

  subroutine calc_strain_rates( self)
    !< Calculate the strain rates

    ! The velocities u and v are defined on the bk-grid (triangles, regular vertical)
    !
    ! The horizontal stretch/shear strain rates du/dx, du/dy, dv/dx, dv/dy are
    ! calculated on the ak-grid (vertices, regular vertical), and are then mapped
    ! to the bks-grid (triangles, staggered vertical)
    !
    ! The vertical shear strain rates du/dz, dv/dz are calculated on the bks-grid
    ! (triangles, staggered vertical), and are then mapped to the ak-grid (vertices,
    ! regular vertical).

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_strain_rates'

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate horizontal stretch/strain rates on the ak-grid
    call calc_3D_gradient_bk_ak(  self%mesh, self%mesh%M_ddx_bk_ak , self%u_bk, self%du_dx_ak )
    call calc_3D_gradient_bk_ak(  self%mesh, self%mesh%M_ddy_bk_ak , self%u_bk, self%du_dy_ak )
    call calc_3D_gradient_bk_ak(  self%mesh, self%mesh%M_ddx_bk_ak , self%v_bk, self%dv_dx_ak )
    call calc_3D_gradient_bk_ak(  self%mesh, self%mesh%M_ddy_bk_ak , self%v_bk, self%dv_dy_ak )

    ! Calculate vertical shear strain rates on the bks-grid
    call calc_3D_gradient_bk_bks( self%mesh, self%mesh%M_ddz_bk_bks, self%u_bk, self%du_dz_bks)
    call calc_3D_gradient_bk_bks( self%mesh, self%mesh%M_ddz_bk_bks, self%v_bk, self%dv_dz_bks)

    ! Map horizontal stretch/shear strain rates from the ak-grid to the bks-grid
    call map_ak_bks( self%mesh, self%mesh%M_map_ak_bks, self%du_dx_ak, self%du_dx_bks)
    call map_ak_bks( self%mesh, self%mesh%M_map_ak_bks, self%du_dy_ak, self%du_dy_bks)
    call map_ak_bks( self%mesh, self%mesh%M_map_ak_bks, self%dv_dx_ak, self%dv_dx_bks)
    call map_ak_bks( self%mesh, self%mesh%M_map_ak_bks, self%dv_dy_ak, self%dv_dy_bks)

    ! Map vertical shear strain rates from the bks-grid to the ak-grid
    call map_bks_ak( self%mesh, self%mesh%M_map_bks_ak, self%du_dz_bks, self%du_dz_ak)
    call map_bks_ak( self%mesh, self%mesh%M_map_bks_ak, self%dv_dz_bks, self%dv_dz_ak)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_strain_rates

  subroutine calc_effective_viscosity( self, ice, geom, Glens_flow_law_epsilon_sq_0_applied)
    !< Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    ! The effective viscosity eta is calculated separately on both the ak-grid (vertices, regular vertical)
    ! and on the bks-grid (triangles, staggered vertical), using the strain rates calculated in calc_strain_rates.
    !
    ! eta_bk, deta_dx_bk, and deta_dy_bk are calculated from eta_ak

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom
    real(dp),                                intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_effective_viscosity'
    real(dp), dimension(:,:), allocatable ::  A_flow_bks
    integer                               :: vi,ti,k,ks
    real(dp)                              :: A_min, eta_max
    real(dp), dimension(:,:), allocatable :: eta_bk_from_ak, eta_bk_from_bks
    real(dp)                              :: uabs_base, uabs_surf, R_shear

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * &
      (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! allocate memory
    allocate( A_flow_bks( self%mesh%ti1:self%mesh%ti2, 1:self%mesh%nz-1))

    ! == Calculate effective viscosity on the ak-grid

    ! Calculate the effective viscosity eta
    select case (C%choice_flow_law)
    case default
      call crash('unknown choice_flow_law "' // trim( C%choice_flow_law) // '"!')
    case ('Glen')
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      call calc_ice_rheology_Glen( self%mesh, ice, geom)

      ! Calculate effective viscosity
      do vi = self%mesh%vi1, self%mesh%vi2
      do k  = 1, self%mesh%nz
        self%eta_ak( vi,k) = calc_effective_viscosity_Glen_3D_uv_only( Glens_flow_law_epsilon_sq_0_applied, &
          self%du_dx_ak( vi,k), self%du_dy_ak( vi,k), self%du_dz_ak( vi,k), &
          self%dv_dx_ak( vi,k), self%dv_dy_ak( vi,k), self%dv_dz_ak( vi,k), ice%A_flow( vi,k))
      end do
      end do

    end select

    ! Safety
    self%eta_ak = min( max( self%eta_ak, C%visc_eff_min), eta_max)

    ! == Calculate effective viscosity on the bks-grid

    ! Calculate the effective viscosity eta
    select case (C%choice_flow_law)
    case default
      call crash('unknown choice_flow_law "' // trim( C%choice_flow_law) // '"!')
    case ('Glen')
      ! Calculate the effective viscosity according to Glen's flow law

      ! Calculate flow factors: map ice flow factor from the ak-grid to the bks-grid
      call map_ak_bks( self%mesh, self%mesh%M_map_ak_bks, ice%A_flow, A_flow_bks)

      ! Calculate effective viscosity
      do ti = self%mesh%ti1, self%mesh%ti2
      do ks  = 1, self%mesh%nz-1
        self%eta_bks( ti,ks) = calc_effective_viscosity_Glen_3D_uv_only( C%Glens_flow_law_epsilon_sq_0, &
          self%du_dx_bks( ti,ks), self%du_dy_bks( ti,ks), self%du_dz_bks( ti,ks), &
          self%dv_dx_bks( ti,ks), self%dv_dy_bks( ti,ks), self%dv_dz_bks( ti,ks), A_flow_bks( ti,ks))
      end do
      end do

    end select

    ! Safety
    self%eta_bks = min( max( self%eta_bks, C%visc_eff_min), eta_max)

    ! Calculate the horizontal gradients of the effective viscosity from its value on the ak-grid
    call calc_3D_gradient_ak_bk(  self%mesh, self%mesh%M_ddx_ak_bk , self%eta_ak , self%deta_dx_bk)
    call calc_3D_gradient_ak_bk(  self%mesh, self%mesh%M_ddy_ak_bk , self%eta_ak , self%deta_dy_bk)

    ! Calculate the vertical gradients of the effective viscosity from its value on the bks-grid
    call calc_3D_gradient_bks_bk( self%mesh, self%mesh%M_ddz_bks_bk, self%eta_bks, self%deta_dz_bk)

    ! Map the effective viscosity from the ak- and bks-grids to the bk-grid

    allocate( eta_bk_from_ak(  self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz))
    allocate( eta_bk_from_bks( self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz))

    call map_a_b_3D( self%mesh, self%eta_ak, eta_bk_from_ak)
    call calc_3D_gradient_bks_bk( self%mesh, self%mesh%M_map_bks_bk, self%eta_bks, eta_bk_from_bks)

    ! Preliminary experiments suggest that in settings where ice flow is dominated by
    ! vertical shear (e.g. the Halfar dome with no sliding), the solver is only stable
    ! when using eta_bk_from_bks. But in settings with a lot of sliding and little
    ! vertical shear (e.g. ISMIP-HOM C), it needs eta_from_ak instead. The "shear factor"
    ! R_shear serves to provide a crude approximation to which flow mode dominates,
    ! which is then use to calculate a weighted average between the two versions of eta.

    do ti = self%mesh%ti1, self%mesh%ti2

      ! Calculate the shear factor R_shear
      uabs_surf = sqrt( 0.1_dp + self%u_bk( ti,1           )**2 + self%v_bk( ti,1           )**2)
      uabs_base = sqrt( 0.1_dp + self%u_bk( ti,self%mesh%nz)**2 + self%v_bk( ti,self%mesh%nz)**2)
      R_shear = uabs_base / uabs_surf

      ! By the nature of ice flow, uabs_base <= uabs_surf, so 0 <= R_shear <= 1,
      ! with 0 indicating no sliding and therefore full vertical shear, and
      ! 1 indicating full sliding.

      ! Weighted average
      self%eta_bk( ti,:) = R_shear * eta_bk_from_ak( ti,:) + (1._dp - R_shear) * eta_bk_from_bks( ti,:)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_applied_basal_friction_coefficient( self, ice, geom, bed_roughness)
    !< Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    !< and scaled with the sub-grid grounded fraction

    ! This is where the sliding law is called!

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom
    type(type_bed_roughness_model),          intent(in   ) :: bed_roughness

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_applied_basal_friction_coefficient'
    integer                                          :: ti
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: u_base_b, v_base_b
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: u_base_a, v_base_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Copy basal velocities from 3-D fields
    u_base_b = self%u_bk( self%mesh%ti1:self%mesh%ti2,self%mesh%nz)
    v_base_b = self%v_bk( self%mesh%ti1:self%mesh%ti2,self%mesh%nz)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!

    ! Map velocities to the a-grid
    call map_b_a_2D( self%mesh, u_base_b, u_base_a)
    call map_b_a_2D( self%mesh, v_base_b, v_base_a)
    call calc_basal_friction_coefficient( self%mesh, geom, bed_roughness, u_base_a, v_base_a, &
      ice%effective_pressure, ice%till_yield_stress, ice%basal_friction_coefficient)

    ! Map basal friction coefficient beta_b to the b-grid
    call map_a_b_2D( self%mesh, ice%basal_friction_coefficient, self%basal_friction_coefficient_b)

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    if (C%do_GL_subgrid_friction) then
      do ti = self%mesh%ti1, self%mesh%ti2
        self%basal_friction_coefficient_b( ti) = self%basal_friction_coefficient_b( ti) * &
          geom%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_applied_basal_friction_coefficient

! == Some useful tools for improving numerical stability of the viscosity iteration

  subroutine relax_viscosity_iterations( self, visc_it_relax)
    !< Reduce the change between velocity solutions

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self
    real(dp),                                intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=*), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                     :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = self%mesh%ti1, self%mesh%ti2
      self%u_bk( ti,:) = (visc_it_relax * self%u_bk( ti,:)) + ((1._dp - visc_it_relax) * self%u_bk_prev( ti,:))
      self%v_bk( ti,:) = (visc_it_relax * self%v_bk( ti,:)) + ((1._dp - visc_it_relax) * self%v_bk_prev( ti,:))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine calc_visc_iter_UV_resid( self, resid_UV)
    !< Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    real(dp),                                intent(  out) :: resid_UV

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                     :: ierr
    integer                     :: ti,k
    real(dp)                    :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = self%mesh%ti1, self%mesh%ti2
    do k  = 1, self%mesh%nz

      res1 = res1 + (self%u_bk( ti,k) - self%u_bk_prev( ti,k))**2
      res1 = res1 + (self%v_bk( ti,k) - self%v_bk_prev( ti,k))**2

      res2 = res2 + (self%u_bk( ti,k) + self%u_bk_prev( ti,k))**2
      res2 = res2 + (self%v_bk( ti,k) + self%v_bk_prev( ti,k))**2

    end do
    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_visc_iter_UV_resid

  subroutine apply_velocity_limits( self)
    !< Limit velocities for improved stability

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'apply_velocity_limits'
    integer                     :: ti,k
    real(dp)                    :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = self%mesh%ti1, self%mesh%ti2
    do k  = 1, self%mesh%nz

      ! Calculate absolute speed
      uabs = SQRT( self%u_bk( ti,k)**2 + self%v_bk( ti,k)**2)

      ! Reduce velocities if neceBPAry
      if (uabs > C%vel_max) then
        self%u_bk( ti,k) = self%u_bk( ti,k) * C%vel_max / uabs
        self%v_bk( ti,k) = self%v_bk( ti,k) * C%vel_max / uabs
      end if

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

! == Restart NetCDF files

  subroutine create_restart_file_old_BPA( self)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'create_restart_file_old_BPA'
    character(len=:), allocatable :: filename_base
    integer                       :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = trim( C%output_dir) // 'restart_ice_velocity_BPA'
    call generate_filename_XXXXXdotnc( filename_base, self%restart_filename)

    ! Print to terminal
    if (par%primary) WRITE(0,'(A)') '   Creating BPA restart file "' // &
      UPSY%stru%colour_string( trim( self%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( self%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( self%restart_filename, ncid, self%mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( self%restart_filename, ncid)

    ! Add a zeta dimension to the file
    call add_zeta_dimension_to_file( self%restart_filename, ncid, self%mesh%zeta)

    ! Add the velocity fields to the file
    call add_field_mesh_dp_3D_b( self%restart_filename, ncid, 'u_bk', long_name = '3-D horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_3D_b( self%restart_filename, ncid, 'v_bk', long_name = '3-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_old_BPA

  subroutine write_to_restart_file_old_BPA( self, time)

    ! In/output variables:
    class(type_momentum_balance_solver_BPA), intent(in   ) :: self
    real(dp),                                intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_to_restart_file_old_BPA'
    integer                     :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) WRITE(0,'(A)') '   Writing to BPA restart file "' // &
      UPSY%stru%colour_string( trim( self%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( self%restart_filename, ncid)

    ! Write the time to the file
    call write_time_to_file( self%restart_filename, ncid, time)

    ! Write the velocity fields to the file
    call write_to_field_multopt_mesh_dp_3D_b( self%mesh, self%restart_filename, ncid, 'u_bk', self%u_bk)
    call write_to_field_multopt_mesh_dp_3D_b( self%mesh, self%restart_filename, ncid, 'v_bk', self%v_bk)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_old_BPA

end module momentum_balance_solver_BPA
