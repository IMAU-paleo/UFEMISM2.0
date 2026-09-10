module momentum_balance_solver_SSA

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA)

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_LOR, MPI_LOGICAL, MPI_MIN, MPI_MAX
  use UPSY_main, only: UPSY
  use mpi_basic, only: par
  use precisions, only: dp
  use parameters, only: grav, ice_density
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use netcdf_io_main
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, ddx_b_a_2D, ddy_b_a_2D, map_b_a_2D
  use constitutive_equation, only: calc_ice_rheology_Glen, calc_effective_viscosity_Glen_2D
  use mesh_zeta, only: vertical_average
  use mpi_distributed_memory, only: gather_to_all
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D
  use bed_roughness_model_types, only: type_bed_roughness_model
  use momentum_balance_solver_SSADIVA, only: atype_momentum_balance_solver_SSADIVA

  implicit none

  private

  public :: type_momentum_balance_solver_SSA

  type, extends(atype_momentum_balance_solver_SSADIVA) :: type_momentum_balance_solver_SSA

      ! Intermediate data fields
      real(dp), dimension(:), allocatable :: A_flow_vav_a                ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter

    contains

      ! Procedures for model memory management and operation
      procedure, public :: get_momentum_balance_solver_name
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_SSA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_SSA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_SSA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_SSA_run
      procedure, public :: set_velocities_to_solver_results   => momentum_balance_solver_SSA_set_velocities
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_SSA_remap

      procedure, public :: create_restart_file_old            => create_restart_file_old_SSA
      procedure, public :: write_to_restart_file_old          => write_to_restart_file_old_SSA

      procedure, public :: initialise_SSA_velocities_from_file
      procedure, public :: calc_vertically_averaged_flow_parameter
      procedure, public :: calc_effective_viscosity
      procedure, public :: calc_applied_basal_friction_coefficient

  end type type_momentum_balance_solver_SSA

contains

  ! == Main routines

  subroutine momentum_balance_solver_SSA_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_allocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate variables that are shared between the SSA and DIVA solvers
    call self%allocate_shared_SSA_DIVA_variables()

    ! Allocate variables that are specific to the SSA solver

    ! Intermediate data fields
    allocate( self%A_flow_vav_a( self%mesh%vi1:self%mesh%vi2), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_allocate

  subroutine momentum_balance_solver_SSA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_deallocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Deallocate variables that are shared between the SSA and DIVA solvers
    call self%deallocate_shared_SSA_DIVA_variables()

    ! Deallocate variables that are specific to the SSA solver

    ! Intermediate data fields
    deallocate( self%A_flow_vav_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_deallocate

  subroutine momentum_balance_solver_SSA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'momentum_balance_solver_SSA_initialise'
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
    case default
      call crash('unknown choice_initial_velocity "' // trim( choice_initial_velocity) // '"!')
    case ('zero')
      self%u_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
      self%v_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
    case ('read_from_file')
      call initialise_SSA_velocities_from_file( self)
    end select

    ! Set tolerances for PETSc matrix solver for the linearised SSA
    self%PETSc_rtol   = C%stress_balance_PETSc_rtol
    self%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_initialise

  subroutine initialise_SSA_velocities_from_file( self)
    ! Initialise the velocities for the SSA solver from an external NetCDF file

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'initialise_SSA_velocities_from_file'
    character(len=:), allocatable                    :: filename
    real(dp)                                         :: timeframe
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: d_loc_b

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

    ! write to terminal
    if (par%primary) write(0,*) '   Initialising SSA velocities from file "' // &
      UPSY%stru%colour_string( trim( filename),'light blue') // '"...'

    ! Read velocities from the file
    call read_field_from_mesh_file_dp_2D_b( filename, 'u_vav_b', d_loc_b, time_to_read = timeframe)
    self%u_vav_b( self%mesh%ti1:self%mesh%ti2) = d_loc_b
    call read_field_from_mesh_file_dp_2D_b( filename, 'v_vav_b', d_loc_b, time_to_read = timeframe)
    self%v_vav_b( self%mesh%ti1:self%mesh%ti2) = d_loc_b

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SSA_velocities_from_file

  subroutine momentum_balance_solver_SSA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the Shallow Shelf Approximation

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
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
    character(len=*), parameter         :: routine_name = 'momentum_balance_solver_SSA_run'
    logical                             :: grounded_ice_exists
    integer                             :: ierr
    integer,  dimension(:), allocatable :: BC_prescr_mask_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_u_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_v_b_applied
    integer                             :: viscosity_iteration_i
    logical                             :: has_converged
    real(dp)                            :: L2_uv, L2_uv_prev
    real(dp)                            :: uv_min, uv_max
    real(dp)                            :: visc_it_relax_applied
    real(dp)                            :: Glens_flow_law_epsilon_sq_0_applied
    integer                             :: nit_diverg_consec
    integer                             :: n_Axb_its_visc_it

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, or no sliding, no need to solve the SSA
    grounded_ice_exists = any( geom%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists .or. C%choice_sliding_law == 'no_sliding') then
      self%u_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
      self%v_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions
    allocate( BC_prescr_mask_b_applied( self%mesh%ti1:self%mesh%ti2))
    allocate( BC_prescr_u_b_applied(    self%mesh%ti1:self%mesh%ti2))
    allocate( BC_prescr_v_b_applied(    self%mesh%ti1:self%mesh%ti2))
    if (present( BC_prescr_mask_b) .or. present( BC_prescr_u_b) .or. present( BC_prescr_v_b)) then
      ! Safety
      if (.not. (present( BC_prescr_mask_b) .and. present( BC_prescr_u_b) .and. present( BC_prescr_v_b))) then
        call crash('need to provide prescribed u,v fields and mask!')
      end if
      BC_prescr_mask_b_applied = BC_prescr_mask_b
      BC_prescr_u_b_applied    = BC_prescr_u_b
      BC_prescr_v_b_applied    = BC_prescr_v_b
    else
      BC_prescr_mask_b_applied = 0
      BC_prescr_u_b_applied    = 0._dp
      BC_prescr_v_b_applied    = 0._dp
    end if

    ! Calculate the driving stress
    call self%calc_driving_stress( geom)

    ! Adaptive relaxation parameter for the viscosity iteration
    L2_uv                               = 1E9_dp
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
      call self%calc_horizontal_strain_rates()

      ! Calculate the effective viscosity for the current velocity solution
      call self%calc_effective_viscosity( ice, geom, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      call self%calc_applied_basal_friction_coefficient( ice, geom, bed_roughness)

      ! Solve the linearised SSA to calculate a new velocity solution
      call self%solve_SSA_DIVA_linearised( self%basal_friction_coefficient_b, n_Axb_its_visc_it, &
        BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Update stability info
      self%n_Axb_its = self%n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call self%apply_velocity_limits()

      ! Reduce the change between velocity solutions
      call self%relax_viscosity_iterations( visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      L2_uv_prev = L2_uv
      call self%calc_L2_norm_uv( L2_uv)

      ! if the viscosity iteration diverges, lower the relaxation parameter
      if (L2_uv > L2_uv_prev) then
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
      uv_min = min( minval( self%u_vav_b( self%mesh%ti1:self%mesh%ti2)), minval( self%v_vav_b( self%mesh%ti1:self%mesh%ti2)))
      uv_max = max( maxval( self%u_vav_b( self%mesh%ti1:self%mesh%ti2)), maxval( self%v_vav_b( self%mesh%ti1:self%mesh%ti2)))
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      ! if (par%primary) write(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], L2_uv = ', L2_uv

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (L2_uv < C%visc_it_norm_dUV_tol) then
        has_converged = .TRUE.
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

  end subroutine momentum_balance_solver_SSA_run

  subroutine momentum_balance_solver_SSA_set_velocities( self, ice, vel)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(in   ) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_velocity_model),         intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_set_velocities'
    integer                     :: ti, vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Velocities
    do ti = self%mesh%ti1, self%mesh%ti2
      vel%u_3D_b( ti,:) = self%u_vav_b( ti)
      vel%v_3D_b( ti,:) = self%v_vav_b( ti)
    end do

    ! Strain rates
    do vi = self%mesh%vi1, self%mesh%vi2
      vel%du_dx_3D( vi,:) = self%du_dx_a( vi)
      vel%du_dy_3D( vi,:) = self%du_dy_a( vi)
      vel%dv_dx_3D( vi,:) = self%dv_dx_a( vi)
      vel%dv_dy_3D( vi,:) = self%dv_dy_a( vi)
    end do

    ! In the SSA, vertical gradients of u,v, and all gradients of w, are neglected
    vel%du_dz_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dv_dz_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dw_dx_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_set_velocities

  subroutine momentum_balance_solver_SSA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
    type(type_mesh),                         intent(in   ) :: mesh_old
    type(type_mesh), target,                 intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_remap'

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap variables that are shared between the SSA and DIVA solvers
    call self%remap_shared_SSA_DIVA_variables( mesh_old, mesh_new)

    ! Remap variables that are specific to the SSA solver

    call reallocate_bounds( self%A_flow_vav_a, mesh_new%vi1, mesh_new%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SSA), intent(in) :: self
    character(len=:), allocatable                  :: model_name
    model_name = 'SSA'
  end function get_momentum_balance_solver_name

  ! == Calculate several intermediate terms in the SSA

  subroutine calc_vertically_averaged_flow_parameter( self, ice)
    !< Calculate the vertical average of Glen's flow parameter A

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
    class(atype_ice_model_data),             intent(in   ) :: ice

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'calc_vertically_averaged_flow_parameter'
    integer                            :: vi
    real(dp), dimension( self%mesh%nz) :: A_prof

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the vertical average of Glen's flow parameter A
    do vi = self%mesh%vi1, self%mesh%vi2
      A_prof = ice%A_flow( vi,:)
      self%A_flow_vav_a( vi) = vertical_average( self%mesh%zeta, A_prof)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertically_averaged_flow_parameter

  subroutine calc_effective_viscosity( self, ice, geom, Glens_flow_law_epsilon_sq_0_applied)
    !< Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom
    real(dp),                                intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_effective_viscosity'
    real(dp)                    :: A_min, eta_max
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * &
      (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! Calculate the effective viscosity eta
    if (C%choice_flow_law == 'Glen') then
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      call calc_ice_rheology_Glen( self%mesh, ice, geom)
      call self%calc_vertically_averaged_flow_parameter( ice)

      ! Calculate effective viscosity
      do vi = self%mesh%vi1, self%mesh%vi2
        self%eta_vav_a( vi) = calc_effective_viscosity_Glen_2D( Glens_flow_law_epsilon_sq_0_applied, &
          self%du_dx_a( vi), self%du_dy_a( vi), self%dv_dx_a( vi), self%dv_dy_a( vi), self%A_flow_vav_a( vi))
      end do

    else
      call crash('unknown choice_flow_law "' // trim( C%choice_flow_law) // '"!')
    end if

    do vi = self%mesh%vi1, self%mesh%vi2
      ! Safety
      self%eta_vav_a( vi) = min( max( self%eta_vav_a( vi), C%visc_eff_min), eta_max)

      ! Calculate the product term N = eta * H on the a-grid
      self%N_a( vi) = self%eta_vav_a( vi) * max( 0.1_dp, geom%Hi( vi))
    end do

    ! Calculate the product term N and its gradients on the b-grid
    call map_a_b_2D( self%mesh, self%N_a, self%N_b    )
    call ddx_a_b_2D( self%mesh, self%N_a, self%dN_dx_b)
    call ddy_a_b_2D( self%mesh, self%N_a, self%dN_dy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_applied_basal_friction_coefficient( self, ice, geom, bed_roughness)
    !< Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    !< and scaled with the sub-grid grounded fraction

    ! NOTE: this is where the sliding law is called!

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom
    type(type_bed_roughness_model),          intent(in   ) :: bed_roughness

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_applied_basal_friction_coefficient'
    integer                                          :: ti
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: u_vav_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: v_vav_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the basal friction coefficient for the current velocity solution
    ! This is where the sliding law is called!
    call map_b_a_2D( self%mesh, self%u_vav_b, u_vav_a)
    call map_b_a_2D( self%mesh, self%v_vav_b, v_vav_a)
    call calc_basal_friction_coefficient( self%mesh, geom, bed_roughness, u_vav_a, v_vav_a, &
      ice%effective_pressure, ice%till_yield_stress, ice%basal_friction_coefficient)

    ! Map the basal friction coefficient to the b-grid
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

  ! == Restart NetCDF files

  subroutine create_restart_file_old_SSA( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'create_restart_file_old_SSA'
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
    filename_base = trim( C%output_dir) // 'restart_ice_velocity_SSA'
    call generate_filename_XXXXXdotnc( filename_base, self%restart_filename)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating SSA restart file "' // &
      UPSY%stru%colour_string( trim( self%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( self%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( self%restart_filename, ncid, self%mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( self%restart_filename, ncid)

    ! Add the velocity fields to the file
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'u_vav_b', long_name = '2-D horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'v_vav_b', long_name = '2-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_old_SSA

  subroutine write_to_restart_file_old_SSA( self, time)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(in   ) :: self
    real(dp),                                intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_to_restart_file_old_SSA'
    integer                     :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to SSA restart file "' // &
      UPSY%stru%colour_string( trim( self%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( self%restart_filename, ncid)

    ! write the time to the file
    call write_time_to_file( self%restart_filename, ncid, time)

    ! write the velocity fields to the file
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'u_vav_b', self%u_vav_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'v_vav_b', self%v_vav_b)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_old_SSA

end module momentum_balance_solver_SSA
