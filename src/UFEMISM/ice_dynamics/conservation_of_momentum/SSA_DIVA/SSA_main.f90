module SSA_main

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA)

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_LOR, MPI_LOGICAL, MPI_MIN, MPI_MAX
  use mpi_basic, only: par
  use precisions, only: dp
  use parameters, only: grav, ice_density
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_SSA
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist
  use netcdf_io_main
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, ddx_b_a_2D, ddy_b_a_2D, map_b_a_2D
  use constitutive_equation, only: calc_ice_rheology_Glen, calc_effective_viscosity_Glen_2D
  use mesh_zeta, only: vertical_average
  use mpi_distributed_memory, only: gather_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use SSA_DIVA_utilities, only: calc_driving_stress, calc_horizontal_strain_rates, relax_viscosity_iterations, &
    apply_velocity_limits, calc_L2_norm_uv
  use solve_linearised_SSA_DIVA, only: solve_SSA_DIVA_linearised
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

contains

  ! == Main routines

  subroutine initialise_SSA_solver( mesh, SSA, region_name)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(  out) :: SSA
    character(len=3),                   intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SSA_solver'
    character(len=256)             :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    call allocate_SSA_solver( mesh, SSA)

    ! Determine the choice of initial velocities for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
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
      SSA%u_b = 0._dp
      SSA%v_b = 0._dp
    case ('read_from_file')
      call initialise_SSA_velocities_from_file( mesh, SSA, region_name)
    end select

    ! Set tolerances for PETSc matrix solver for the linearised SSA
    SSA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    SSA%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SSA_solver

  subroutine solve_SSA( mesh, ice, bed_roughness, SSA, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Calculate ice velocities by solving the Shallow Shelf Approximation

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_bed_roughness_model),     intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA
    integer,                            intent(  out) :: n_visc_its            ! Number of non-linear viscosity iterations
    integer,                            intent(  out) :: n_Axb_its             ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:), optional,   intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,   intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,   intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_SSA'
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
    grounded_ice_exists = any( ice%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists .or. C%choice_sliding_law == 'no_sliding') then
      SSA%u_b = 0._dp
      SSA%v_b = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions
    allocate( BC_prescr_mask_b_applied( mesh%ti1:mesh%ti2))
    allocate( BC_prescr_u_b_applied(    mesh%ti1:mesh%ti2))
    allocate( BC_prescr_v_b_applied(    mesh%ti1:mesh%ti2))
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
    call calc_driving_stress( mesh, ice, SSA%tau_dx_b, SSA%tau_dy_b)

    ! Adaptive relaxation parameter for the viscosity iteration
    L2_uv                               = 1E9_dp
    nit_diverg_consec                   = 0
    visc_it_relax_applied               = C%visc_it_relax
    Glens_flow_law_epsilon_sq_0_applied = C%Glens_flow_law_epsilon_sq_0

    ! Initialise stability info
    n_visc_its = 0
    n_Axb_its  = 0

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.
    viscosity_iteration: do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the strain rates for the current velocity solution
      call calc_horizontal_strain_rates( mesh, SSA%u_b, SSA%v_b, SSA%du_dx_a, SSA%du_dy_a, SSA%dv_dx_a, SSA%dv_dy_a)

      ! Calculate the effective viscosity for the current velocity solution
      call calc_effective_viscosity( mesh, ice, SSA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      call calc_applied_basal_friction_coefficient( mesh, ice, bed_roughness, SSA)

      ! Solve the linearised SSA to calculate a new velocity solution
      call solve_SSA_DIVA_linearised( mesh, SSA%u_b, SSA%v_b, SSA%N_b, SSA%dN_dx_b, SSA%dN_dy_b, &
      SSA%basal_friction_coefficient_b, SSA%tau_dx_b, SSA%tau_dy_b, SSA%u_b_prev, SSA%v_b_prev, &
        SSA%PETSc_rtol, SSA%PETSc_abstol, n_Axb_its_visc_it, &
        BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Update stability info
      n_Axb_its = n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call apply_velocity_limits( mesh, SSA%u_b, SSA%v_b)

      ! Reduce the change between velocity solutions
      call relax_viscosity_iterations( mesh, SSA%u_b, SSA%v_b, SSA%u_b_prev, SSA%v_b_prev, visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      L2_uv_prev = L2_uv
      call calc_L2_norm_uv( mesh, SSA%u_b, SSA%v_b, SSA%u_b_prev, SSA%v_b_prev, L2_uv)

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
      uv_min = minval( SSA%u_b)
      uv_max = maxval( SSA%u_b)
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
    n_visc_its = viscosity_iteration_i

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA

  subroutine remap_SSA_solver( mesh_old, mesh_new, SSA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh_old
    type(type_mesh),                    intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'remap_SSA_solver'
    real(dp), dimension(:    ), allocatable :: u_a
    real(dp), dimension(:    ), allocatable :: v_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_a( mesh_old%vi1: mesh_old%vi2))
    allocate( v_a( mesh_old%vi1: mesh_old%vi2))

    ! Map velocities from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, SSA%u_b, u_a)
    call map_b_a_2D( mesh_old, SSA%v_b, v_a)

    ! Remap velocities from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, u_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, v_a, '2nd_order_conservative')

    ! reallocate memory for the velocities on the triangles
    call reallocate_bounds( SSA%u_b                         , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( SSA%v_b                         , mesh_new%ti1, mesh_new%ti2)

    ! Map velocities from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, u_a, SSA%u_b)
    call map_a_b_2D( mesh_new, v_a, SSA%v_b)

    ! Clean up after yourself
    deallocate( u_a)
    deallocate( v_a)

    ! reallocate everything else
    ! ==========================

    call reallocate_bounds( SSA%A_flow_vav_a                , mesh_new%vi1, mesh_new%vi2)           ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter
    call reallocate_bounds( SSA%du_dx_a                     , mesh_new%vi1, mesh_new%vi2)           ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( SSA%du_dy_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( SSA%dv_dx_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( SSA%dv_dy_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( SSA%eta_a                       , mesh_new%vi1, mesh_new%vi2)           ! Effective viscosity
    call reallocate_bounds( SSA%N_a                         , mesh_new%vi1, mesh_new%vi2)           ! Product term N = eta * H
    call reallocate_bounds( SSA%N_b                         , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( SSA%dN_dx_b                     , mesh_new%ti1, mesh_new%ti2)           ! Gradients of N
    call reallocate_bounds( SSA%dN_dy_b                     , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( SSA%basal_friction_coefficient_b, mesh_new%ti1, mesh_new%ti2)           ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    call reallocate_bounds( SSA%tau_dx_b                    , mesh_new%ti1, mesh_new%ti2)           ! Driving stress
    call reallocate_bounds( SSA%tau_dy_b                    , mesh_new%ti1, mesh_new%ti2)
    call reallocate_clean ( SSA%u_b_prev                    , mesh_new%nTri             )           ! Velocity solution from previous viscosity iteration
    call reallocate_clean ( SSA%v_b_prev                    , mesh_new%nTri             )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SSA_solver

  ! == Calculate several intermediate terms in the SSA

  subroutine calc_vertically_averaged_flow_parameter( mesh, ice, SSA)
    !< Calculate the vertical average of Glen's flow parameter A

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(in   ) :: ice
    type(type_ice_velocity_solver_SSA),  intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_vertically_averaged_flow_parameter'
    integer                        :: vi
    real(dp), dimension( mesh%nz)  :: A_prof

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the vertical average of Glen's flow parameter A
    do vi = mesh%vi1, mesh%vi2
      A_prof = ice%A_flow( vi,:)
      SSA%A_flow_vav_a( vi) = vertical_average( mesh%zeta, A_prof)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertically_averaged_flow_parameter

  subroutine calc_effective_viscosity( mesh, ice, SSA, Glens_flow_law_epsilon_sq_0_applied)
    !< Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA
    real(dp),                           intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_effective_viscosity'
    real(dp)                       :: A_min, eta_max
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! Calculate the effective viscosity eta
    if (C%choice_flow_law == 'Glen') then
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      call calc_ice_rheology_Glen( mesh, ice)
      call calc_vertically_averaged_flow_parameter( mesh, ice, SSA)

      ! Calculate effective viscosity
      do vi = mesh%vi1, mesh%vi2
        SSA%eta_a( vi) = calc_effective_viscosity_Glen_2D( Glens_flow_law_epsilon_sq_0_applied, &
          SSA%du_dx_a( vi), SSA%du_dy_a( vi), SSA%dv_dx_a( vi), SSA%dv_dy_a( vi), SSA%A_flow_vav_a( vi))
      end do

    else
      call crash('unknown choice_flow_law "' // trim( C%choice_flow_law) // '"!')
    end if

    ! Safety
    SSA%eta_a = min( max( SSA%eta_a, C%visc_eff_min), eta_max)

    ! Calculate the product term N = eta * H on the a-grid
    SSA%N_a = SSA%eta_a * max( 0.1_dp, ice%Hi)

    ! Calculate the product term N and its gradients on the b-grid
    call map_a_b_2D( mesh, SSA%N_a, SSA%N_b    )
    call ddx_a_b_2D( mesh, SSA%N_a, SSA%dN_dx_b)
    call ddy_a_b_2D( mesh, SSA%N_a, SSA%dN_dy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_applied_basal_friction_coefficient( mesh, ice, bed_roughness, SSA)
    !< Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    !< and scaled with the sub-grid grounded fraction

    ! NOTE: this is where the sliding law is called!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_bed_roughness_model),     intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_applied_basal_friction_coefficient'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the basal friction coefficient for the current velocity solution
    ! This is where the sliding law is called!
    call calc_basal_friction_coefficient( mesh, ice, bed_roughness, SSA%u_b, SSA%v_b)

    ! Map the basal friction coefficient to the b-grid
    call map_a_b_2D( mesh, ice%basal_friction_coefficient, SSA%basal_friction_coefficient_b)

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    if (C%do_GL_subgrid_friction) then
      do ti = mesh%ti1, mesh%ti2
        SSA%basal_friction_coefficient_b( ti) = SSA%basal_friction_coefficient_b( ti) * ice%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_applied_basal_friction_coefficient

  ! == Initialisation

  subroutine initialise_SSA_velocities_from_file( mesh, SSA, region_name)
    ! Initialise the velocities for the SSA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA
    character(len=3),                   intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SSA_velocities_from_file'
    real(dp)                       :: dummy1
    character(len=256)             :: filename
    real(dp)                       :: timeframe

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings
    dummy1 = mesh%xmin

    ! Determine the filename and timeframe to read for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
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
      colour_string( trim( filename),'light blue') // '"...'

    ! Read velocities from the file
    if (timeframe == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_mesh_file_dp_2D_b( filename, 'u_b', SSA%u_b)
      call read_field_from_mesh_file_dp_2D_b( filename, 'v_b', SSA%v_b)
    else
      ! Read specified timeframe
      call read_field_from_mesh_file_dp_2D_b( filename, 'u_b', SSA%u_b, time_to_read = timeframe)
      call read_field_from_mesh_file_dp_2D_b( filename, 'v_b', SSA%v_b, time_to_read = timeframe)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SSA_velocities_from_file

  subroutine allocate_SSA_solver( mesh, SSA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(  out) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SSA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( SSA%u_b(                          mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( SSA%v_b(                          mesh%ti1:mesh%ti2), source = 0._dp)

    ! Intermediate data fields
    allocate( SSA%A_flow_vav_a(                 mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SSA%du_dx_a(                      mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SSA%du_dy_a(                      mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SSA%dv_dx_a(                      mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SSA%dv_dy_a(                      mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SSA%eta_a(                        mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SSA%N_a(                          mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( SSA%N_b(                          mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( SSA%dN_dx_b(                      mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( SSA%dN_dy_b(                      mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( SSA%basal_friction_coefficient_b( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( SSA%tau_dx_b(                     mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( SSA%tau_dy_b(                     mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( SSA%u_b_prev(                     mesh%nTri        ), source = 0._dp)
    allocate( SSA%v_b_prev(                     mesh%nTri        ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SSA_solver

  ! == Restart NetCDF files

  subroutine write_to_restart_file_SSA( mesh, SSA, time)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(in   ) :: SSA
    real(dp),                           intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_SSA'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to SSA restart file "' // &
      colour_string( trim( SSA%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( SSA%restart_filename, ncid)

    ! write the time to the file
    call write_time_to_file( SSA%restart_filename, ncid, time)

    ! write the velocity fields to the file
    call write_to_field_multopt_mesh_dp_2D_b( mesh, SSA%restart_filename, ncid, 'u_b', SSA%u_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, SSA%restart_filename, ncid, 'v_b', SSA%v_b)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_SSA

  subroutine create_restart_file_SSA( mesh, SSA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_SSA'
    character(len=256)             :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = trim( C%output_dir) // 'restart_ice_velocity_SSA'
    call generate_filename_XXXXXdotnc( filename_base, SSA%restart_filename)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating SSA restart file "' // &
      colour_string( trim( SSA%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( SSA%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( SSA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( SSA%restart_filename, ncid)

    ! Add the velocity fields to the file
    call add_field_mesh_dp_2D_b( SSA%restart_filename, ncid, 'u_b', long_name = '2-D horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( SSA%restart_filename, ncid, 'v_b', long_name = '2-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_SSA

end module SSA_main
