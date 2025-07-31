module predictor_corrector_scheme

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_INTEGER, &
    MPI_MAX, MPI_SUM
  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, crash
  use model_configuration, only: C
  use region_types, only: type_model_region
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_pc
  use reallocate_mod, only: reallocate_bounds
  use netcdf_io_main
  use time_step_criteria, only: calc_critical_timestep_adv
  use conservation_of_mass_main, only: calc_dHi_dt
  use ice_thickness_safeties, only: alter_ice_thickness
  use ice_geometry_basics, only: ice_surface_elevation
  use masks_mod, only: determine_masks
  use subgrid_grounded_fractions_main, only: calc_grounded_fractions
  use conservation_of_momentum_main, only: solve_stress_balance
  use subgrid_ice_margin, only: calc_effective_thickness

  implicit none

  private

  public :: run_ice_dynamics_model_pc, initialise_pc_scheme, remap_pc_scheme, &
    create_restart_file_pc_scheme, write_to_restart_file_pc_scheme

contains

  subroutine run_ice_dynamics_model_pc( region, dt_max)
    !< Calculate a new next modelled ice thickness

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in   ) :: dt_max

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'run_ice_dynamics_model_pc'
    real(dp)                                             :: dt_crit_adv
    integer                                              :: pc_it
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: Hi_dummy
    integer                                              :: vi, n_guilty, n_tot
    integer                                              :: n_visc_its
    integer                                              :: n_Axb_its
    integer                                              :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Store previous ice model state
    region%ice%t_Hi_prev  = region%ice%t_Hi_next
    region%ice%Hi_prev    = region%ice%Hi_next

    ! == Calculate time step ==
    ! =========================

    ! Store previous time step
    region%ice%pc%dt_n = region%ice%pc%dt_np1

    ! Calculate new time step (Robinson et al., 2020, Eq. 33)
    region%ice%pc%dt_np1 = (C%pc_epsilon / region%ice%pc%eta_np1)**(C%pc_k_I + C%pc_k_p) * &
      (C%pc_epsilon / region%ice%pc%eta_n)**(-C%pc_k_p) * region%ice%pc%dt_n

    ! Limit time step to maximum allowed value
    region%ice%pc%dt_np1 = MIN( region%ice%pc%dt_np1, dt_max)

    ! Limit time step to 1.2 times the previous time step
    region%ice%pc%dt_np1 = MIN( region%ice%pc%dt_np1, C%pc_max_time_step_increase * region%ice%pc%dt_n)

    ! Limit time step to minimum allowed value
    region%ice%pc%dt_np1 = MAX( region%ice%pc%dt_np1, C%dt_ice_min)

    ! Limit time step to critical advective time step
    call calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_adv)

    region%ice%pc%dt_np1 = MIN( region%ice%pc%dt_np1, dt_crit_adv)

    ! == Time step iteration: if, at the end of the PC timestep, the truncation error
    !    turns out to be too large, run it again with a smaller dt, until the truncation
    !    decreases to below the specified tolerance
    ! ==================================================================================

    ! Store thinning rates from previous time step
    region%ice%pc%dHi_dt_Hi_nm1_u_nm1 = region%ice%dHi_dt

    ! Store the previous maximum truncation error eta_n
    region%ice%pc%eta_n = region%ice%pc%eta_np1

    pc_it = 0
    iterate_pc_timestep: do while (pc_it < C%pc_nit_max)

      pc_it = pc_it + 1

      ! Calculate time step ratio
      region%ice%pc%zeta_t = region%ice%pc%dt_np1 / region%ice%pc%dt_n

      ! == Predictor step ==
      ! ====================

      ! Calculate thinning rates for current geometry and velocity
      call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, region%AMB%AMB, region%ice%fraction_margin, &
                        region%ice%mask_noice, region%ice%pc%dt_np1, region%ice%pc%dHi_dt_Hi_n_u_n, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

      ! Calculate predicted ice thickness (Robinson et al., 2020, Eq. 30)
      region%ice%pc%Hi_star_np1 = region%ice%Hi_prev + region%ice%pc%dt_np1 * ((1._dp + region%ice%pc%zeta_t / 2._dp) * &
        region%ice%pc%dHi_dt_Hi_n_u_n - (region%ice%pc%zeta_t / 2._dp) * region%ice%pc%dHi_dt_Hi_nm1_u_nm1)

      ! if so desired, modify the predicted ice thickness field based on user-defined settings
      call alter_ice_thickness( region%mesh, region%ice, region%ice%Hi_prev, region%ice%pc%Hi_star_np1, region%refgeo_PD, region%time)

      ! Adjust the predicted dHi_dt to compensate for thickness modifications
      ! This is just Robinson et al., 2020, Eq 30 above rearranged to retrieve
      ! an updated dHi_dt_Hi_n_u_n from the modified Hi_star_np1. if no ice
      ! thickness modifications were applied, then there will be not change.
      region%ice%pc%dHi_dt_Hi_n_u_n = ((region%ice%pc%Hi_star_np1 - region%ice%Hi_prev) / region%ice%pc%dt_np1 + (region%ice%pc%zeta_t / 2._dp) * region%ice%pc%dHi_dt_Hi_nm1_u_nm1) / (1._dp + region%ice%pc%zeta_t / 2._dp)

      ! == Update step ==
      ! =================

      ! Set model geometry to predicted
      region%ice%Hi = region%ice%pc%Hi_star_np1

      ! Set thinning rates to predicted
      region%ice%dHi_dt = (region%ice%Hi - region%ice%Hi_prev) / region%ice%pc%dt_np1

      ! Set model geometry to predicted
      do vi = region%mesh%vi1, region%mesh%vi2
        ! Basic geometry
        region%ice%Hs ( vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
        region%ice%Hib( vi) = region%ice%Hs(  vi) - region%ice%Hi( vi)
      end do

      ! Update masks
      call determine_masks( region%mesh, region%ice)

      ! DENK DROM : assess whether this is important for the velocitiy computation below
      ! ! Calculate zeta gradients
      ! call calc_zeta_gradients( region%mesh, region%ice)

      ! Update sub-grid grounded fractions
      call calc_grounded_fractions( region%mesh, region%ice)

      ! DENK DROM : assess whether this is important for the velocitiy computation below
      ! ! Calculate the basal mass balance
      ! call run_BMB_model( region%mesh, region%ice, region%ocean, region%refgeo_PD, region%SMB, region%BMB, region%name, region%time)

      ! Calculate ice velocities for the predicted geometry
      call solve_stress_balance( region%mesh, region%ice, region%bed_roughness, &
        region%BMB%BMB, region%name, n_visc_its, n_Axb_its)

      ! Update stability info
      region%ice%dt_ice     = region%ice%pc%dt_np1
      region%ice%n_visc_its = n_visc_its
      region%ice%n_Axb_its  = n_Axb_its

      ! == Corrector step ==
      ! ====================

      ! Set model geometry back to original
      do vi = region%mesh%vi1, region%mesh%vi2
        region%ice%Hi(  vi) = region%ice%Hi_prev( vi)
        region%ice%Hs(  vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
        region%ice%Hib( vi) = region%ice%Hs(  vi) - region%ice%Hi( vi)
      end do

      ! Update masks
      call determine_masks( region%mesh, region%ice)

      ! Update sub-grid grounded fractions
      call calc_grounded_fractions( region%mesh, region%ice)

      ! Update effective ice thickness
      call calc_effective_thickness( region%mesh, region%ice, region%ice%Hi, region%ice%Hi_eff, region%ice%fraction_margin)

      ! Calculate thinning rates for the current ice thickness and predicted velocity
      call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, region%AMB%AMB, region%ice%fraction_margin, &
                        region%ice%mask_noice, region%ice%pc%dt_np1, region%ice%pc%dHi_dt_Hi_star_np1_u_np1, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

      ! Calculate corrected ice thickness (Robinson et al. (2020), Eq. 31)
      region%ice%pc%Hi_np1 = region%ice%Hi_prev + (region%ice%pc%dt_np1 / 2._dp) * (region%ice%pc%dHi_dt_Hi_n_u_n + region%ice%pc%dHi_dt_Hi_star_np1_u_np1)

      ! Save "raw" thinning rates, as applied after the corrector step
      region%ice%dHi_dt_raw = (region%ice%pc%Hi_np1 - region%ice%Hi_prev) / region%ice%pc%dt_np1

      ! if so desired, modify the corrected ice thickness field based on user-defined settings
      call alter_ice_thickness( region%mesh, region%ice, region%ice%Hi_prev, region%ice%pc%Hi_np1, region%refgeo_PD, region%time)

      ! Adjust the predicted dHi_dt to compensate for thickness modifications
      ! This is just Robinson et al., 2020, Eq 31 above rearranged to retrieve
      ! an updated dHi_dt_Hi_star_np1_u_np1 from the modified Hi_np1. if no ice
      ! thickness modifications were applied, then there will be not change.
      region%ice%pc%dHi_dt_Hi_star_np1_u_np1 = (region%ice%pc%Hi_np1 - region%ice%Hi_prev) / (region%ice%pc%dt_np1 / 2._dp) - region%ice%pc%dHi_dt_Hi_n_u_n

      ! Add difference between raw and applied dHi_dt to residual tracker
      region%ice%dHi_dt_residual = region%ice%dHi_dt_raw - &  ! Raw change
                                  (region%ice%pc%Hi_np1 - region%ice%Hi_prev) / region%ice%pc%dt_np1  ! Minus applied change

      ! == Truncation error ==
      ! ======================

      ! Estimate truncation error
      call calc_pc_truncation_error( region%mesh, region%ice, region%ice%pc)

      ! == Error assessment ==
      ! ======================

      ! Initialise unstable vertex count
      n_tot = 0
      n_guilty = 0

      ! Determine number of unstable vertices
      do vi = region%mesh%vi1, region%mesh%vi2
        ! Only consider fully grounded vertices
        if (region%ice%fraction_gr( vi) < 1._dp) CYCLE
        ! if so, add to total vertex count
        n_tot = n_tot + 1
        ! if this vertex's error is larger than tolerance
        if (region%ice%pc%tau_np1( vi) > C%pc_epsilon) then
          ! Add to total guilty vertex count
          n_guilty = n_guilty + 1
          ! Add to this vertex's guilty record
          region%ice%pc%tau_n_guilty( vi) = region%ice%pc%tau_n_guilty( vi) + 1
        end if
      end do

      ! Add up findings from each process domain
      call MPI_ALLREDUCE( MPI_IN_PLACE, n_tot,    1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, n_guilty, 1, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Safety
      if (n_tot == 0) n_tot = 1

      ! Check if largest truncation error is small enough; if so, move on
      if (region%ice%pc%eta_np1 < C%pc_epsilon) then
        exit iterate_pc_timestep

      ! if not, check whether that occurs in a significant amount of vertices; if not,
      ! set the truncation error to almost the tolerance (to allow for growth) and move on
      elseif (100._dp * real( n_guilty,dp) / real(n_tot,dp) < C%pc_guilty_max) then
        ! if (par%primary) call warning('{dp_01}% of vertices are changing rapidly, ignoring for now', dp_01 = 100._dp * real( n_guilty,dp) / real(n_tot,dp))
        region%ice%pc%eta_np1 = .95_dp * C%pc_epsilon
        exit iterate_pc_timestep

      ! if not, re-do the PC timestep
      else
        !if (par%primary) call warning('{dp_01}% of vertices ({int_01}) are changing rapidly (eta = {dp_02}), reducing dt and redoing PC timestep', dp_01 = 100._dp * real( n_guilty,dp) / real(n_tot,dp), int_01 = n_guilty, dp_02 = region%ice%pc%eta_np1)
        region%ice%pc%dt_np1 = region%ice%pc%dt_np1 * 0.8_dp
        ! if the timestep has reached the specified lower limit, stop iterating
        if (region%ice%pc%dt_np1 <= C%dt_ice_min) then
          region%ice%pc%dt_np1 = C%dt_ice_min
          exit iterate_pc_timestep
        end if
      end if

    end do iterate_pc_timestep

    ! == Final quantities
    ! ===================

    ! Set next modelled ice thickness
    region%ice%t_Hi_next = region%ice%t_Hi_prev + region%ice%pc%dt_np1
    region%ice%Hi_next   = region%ice%pc%Hi_np1

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_ice_dynamics_model_pc

  subroutine calc_pc_truncation_error( mesh, ice, pc)
    !< Calculate the truncation error tau in the ice thickness
    !< rate of change (Robinson et al., 2020, Eq. 32)

    ! In- and output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice
    type(type_ice_pc),    intent(inout) :: pc

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_pc_truncation_error'
    integer                        :: vi, ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate truncation error tau (Robinson et al., 2020, Eq. 32)
    do vi = mesh%vi1, mesh%vi2
      pc%tau_np1( vi) = pc%zeta_t * ABS( pc%Hi_np1( vi) - pc%Hi_star_np1( vi)) / ((3._dp * pc%zeta_t + 3._dp) * pc%dt_n)
    end do

    ! Calculate the maximum truncation error eta over grounded ice only
    pc%eta_np1 = C%pc_eta_min
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_grounded_ice( vi) .and. .not. ice%mask_gl_gr( vi) .and. ice%fraction_gr( vi) == 1._dp) then
        pc%eta_np1 = MAX( pc%eta_np1, pc%tau_np1( vi))
      end if
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, pc%eta_np1, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pc_truncation_error

  subroutine initialise_pc_scheme( mesh, pc, region_name)
    !< allocate memory and initialise values for the ice thickness predictor/corrector scheme.

    ! In- and output variables
    type(type_mesh),   intent(in   ) :: mesh
    type(type_ice_pc), intent(  out) :: pc
    character(len=3),  intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_pc_scheme'
    character(len=256)             :: pc_choice_initialise
    character(len=256)             :: filename_pc_initialise
    real(dp)                       :: timeframe_pc_initialise

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    ! ===============

    allocate( pc%dHi_dt_Hi_nm1_u_nm1(      mesh%vi1:mesh%vi2))   ! [m/yr] Thinning rates from previous time step
    allocate( pc%dHi_dt_Hi_n_u_n(          mesh%vi1:mesh%vi2))   ! [m/yr] Thinning rates for current time step with old geometry
    allocate( pc%Hi_star_np1(              mesh%vi1:mesh%vi2))   ! [m]    Predicted ice thickness
    allocate( pc%dHi_dt_Hi_star_np1_u_np1( mesh%vi1:mesh%vi2))   ! [m/yr] Thinning rates for predicted ice thickness and updated velocity
    allocate( pc%Hi_np1(                   mesh%vi1:mesh%vi2))   ! [m]    Corrected ice thickness
    allocate( pc%tau_np1(                  mesh%vi1:mesh%vi2))   ! [m]    Truncation error
    allocate( pc%tau_n_guilty(             mesh%vi1:mesh%vi2))   ! [-]    Number of PC iterations where vertex had truncation errors above the tolerance

    ! Initialise
    ! ==========

    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"!')
    case ('NAM')
      pc_choice_initialise    = C%pc_choice_initialise_NAM
      filename_pc_initialise  = C%filename_pc_initialise_NAM
      timeframe_pc_initialise = C%timeframe_pc_initialise_NAM
    case ('EAS')
      pc_choice_initialise    = C%pc_choice_initialise_EAS
      filename_pc_initialise  = C%filename_pc_initialise_EAS
      timeframe_pc_initialise = C%timeframe_pc_initialise_EAS
    case ('GRL')
      pc_choice_initialise    = C%pc_choice_initialise_GRL
      filename_pc_initialise  = C%filename_pc_initialise_GRL
      timeframe_pc_initialise = C%timeframe_pc_initialise_GRL
    case ('ANT')
      pc_choice_initialise    = C%pc_choice_initialise_ANT
      filename_pc_initialise  = C%filename_pc_initialise_ANT
      timeframe_pc_initialise = C%timeframe_pc_initialise_ANT
    end select

    select case (pc_choice_initialise)
    case default
      call crash('unknown pc_choice_initialise "' // trim( pc_choice_initialise) // '"!')
    case ('zero')
      ! Initialise everything from scratch

      pc%dt_n                     = C%dt_ice_min
      pc%dt_np1                   = C%dt_ice_min
      pc%zeta_t                   = 1._dp
      pc%dHi_dt_Hi_nm1_u_nm1      = 0._dp
      pc%dHi_dt_Hi_n_u_n          = 0._dp
      pc%Hi_star_np1              = 0._dp
      pc%dHi_dt_Hi_star_np1_u_np1 = 0._dp
      pc%Hi_np1                   = 0._dp
      pc%tau_np1                  = C%pc_epsilon
      pc%eta_n                    = C%pc_epsilon
      pc%eta_np1                  = C%pc_epsilon

    case ('read_from_file')
      ! Initialise from a (restart) file
      call initialise_pc_scheme_from_file( pc, filename_pc_initialise, timeframe_pc_initialise)
    end select

    ! Initialise the event counter for errors above tolerance
    pc%tau_n_guilty = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_pc_scheme

  subroutine initialise_pc_scheme_from_file( pc, filename, timeframe)
    !< Initialise values for the ice thickness predictor/corrector scheme from a (restart) file.

    ! In- and output variables
    type(type_ice_pc),  intent(inout) :: pc
    character(len=256), intent(in   ) :: filename
    real(dp),           intent(in   ) :: timeframe

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_pc_scheme_from_file'
    character(len=1024)            :: filename_applied
    real(dp)                       :: timeframe_applied

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    filename_applied  = filename
    timeframe_applied = timeframe
    if (index( filename_applied,'_LAST.nc') > 1) then
      call find_last_output_file( filename_applied)
      call find_last_timeframe(   filename_applied, timeframe_applied)
    end if

    ! write to terminal
    if (par%primary) write(0,*) '   Initialising ice thickness predictor/corrector scheme from file "' // colour_string( trim( filename),'light blue') // '"...'

    ! Read values from the file
    if (timeframe_applied == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_file_0D(         filename_applied, 'dt_n'                    , pc%dt_n                    )
      call read_field_from_file_0D(         filename_applied, 'dt_np1'                  , pc%dt_np1                  )
      call read_field_from_file_0D(         filename_applied, 'zeta_t'                  , pc%zeta_t                  )
      call read_field_from_mesh_file_dp_2D( filename_applied, 'dHi_dt_Hi_nm1_u_nm1'     , pc%dHi_dt_Hi_nm1_u_nm1     )
      call read_field_from_mesh_file_dp_2D( filename_applied, 'dHi_dt_Hi_n_u_n'         , pc%dHi_dt_Hi_n_u_n         )
      call read_field_from_mesh_file_dp_2D( filename_applied, 'Hi_star_np1'             , pc%Hi_star_np1             )
      call read_field_from_mesh_file_dp_2D( filename_applied, 'dHi_dt_Hi_star_np1_u_np1', pc%dHi_dt_Hi_star_np1_u_np1)
      call read_field_from_mesh_file_dp_2D( filename_applied, 'Hi_np1'                  , pc%Hi_np1                  )
      call read_field_from_mesh_file_dp_2D( filename_applied, 'tau_np1'                 , pc%tau_np1                 )
      call read_field_from_file_0D(         filename_applied, 'eta_n'                   , pc%eta_n                   )
      call read_field_from_file_0D(         filename_applied, 'eta_np1'                 , pc%eta_np1                 )
    else
      ! Read specified timeframe
      call read_field_from_file_0D(         filename_applied, 'dt_n'                    , pc%dt_n                    , time_to_read = timeframe_applied)
      call read_field_from_file_0D(         filename_applied, 'dt_np1'                  , pc%dt_np1                  , time_to_read = timeframe_applied)
      call read_field_from_file_0D(         filename_applied, 'zeta_t'                  , pc%zeta_t                  , time_to_read = timeframe_applied)
      call read_field_from_mesh_file_dp_2D( filename_applied, 'dHi_dt_Hi_nm1_u_nm1'     , pc%dHi_dt_Hi_nm1_u_nm1     , time_to_read = timeframe_applied)
      call read_field_from_mesh_file_dp_2D( filename_applied, 'dHi_dt_Hi_n_u_n'         , pc%dHi_dt_Hi_n_u_n         , time_to_read = timeframe_applied)
      call read_field_from_mesh_file_dp_2D( filename_applied, 'Hi_star_np1'             , pc%Hi_star_np1             , time_to_read = timeframe_applied)
      call read_field_from_mesh_file_dp_2D( filename_applied, 'dHi_dt_Hi_star_np1_u_np1', pc%dHi_dt_Hi_star_np1_u_np1, time_to_read = timeframe_applied)
      call read_field_from_mesh_file_dp_2D( filename_applied, 'Hi_np1'                  , pc%Hi_np1                  , time_to_read = timeframe_applied)
      call read_field_from_mesh_file_dp_2D( filename_applied, 'tau_np1'                 , pc%tau_np1                 , time_to_read = timeframe_applied)
      call read_field_from_file_0D(         filename_applied, 'eta_n'                   , pc%eta_n                   , time_to_read = timeframe_applied)
      call read_field_from_file_0D(         filename_applied, 'eta_np1'                 , pc%eta_np1                 , time_to_read = timeframe_applied)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_pc_scheme_from_file

  subroutine write_to_restart_file_pc_scheme( mesh, pc, time)
    !< write to the restart NetCDF file for the ice thickness predictor/corrector scheme

    ! In/output variables:
    type(type_mesh),   intent(in   ) :: mesh
    type(type_ice_pc), intent(in   ) :: pc
    real(dp),          intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_pc_scheme'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to ice dynamics restart file "' // &
      colour_string( trim( pc%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( pc%restart_filename, ncid)

    ! write the time to the file
    call write_time_to_file( pc%restart_filename, ncid, time)

    ! write the data fields to the file
    call write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'dt_n'                    , pc%dt_n                    )
    call write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'dt_np1'                  , pc%dt_np1                  )
    call write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'zeta_t'                  , pc%zeta_t                  )
    call write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'dHi_dt_Hi_nm1_u_nm1'     , pc%dHi_dt_Hi_nm1_u_nm1     )
    call write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'dHi_dt_Hi_n_u_n'         , pc%dHi_dt_Hi_n_u_n         )
    call write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'Hi_star_np1'             , pc%Hi_star_np1             )
    call write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'dHi_dt_Hi_star_np1_u_np1', pc%dHi_dt_Hi_star_np1_u_np1)
    call write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'Hi_np1'                  , pc%Hi_np1                  )
    call write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'tau_np1'                 , pc%tau_np1                 )
    call write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'eta_n'                   , pc%eta_n                   )
    call write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'eta_np1'                 , pc%eta_np1                 )

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_pc_scheme

  subroutine create_restart_file_pc_scheme( mesh, pc)
    !< Create a restart NetCDF file for the ice thickness predictor/corrector scheme
    !< Includes generation of the procedural filename (e.g. "restart_pc_00001.nc")

    ! In/output variables:
    type(type_mesh),   intent(in   ) :: mesh
    type(type_ice_pc), intent(inout) :: pc

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_pc_scheme'
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
    filename_base = trim( C%output_dir) // 'restart_pc_scheme'
    call generate_filename_XXXXXdotnc( filename_base, pc%restart_filename)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating ice dynamics restart file "' // &
      colour_string( trim( pc%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( pc%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( pc%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( pc%restart_filename, ncid)

    ! Add the data fields to the file
    call add_field_dp_0D(      pc%restart_filename, ncid, 'dt_n'                    , long_name = 'Previous time step', units = 'yr')
    call add_field_dp_0D(      pc%restart_filename, ncid, 'dt_np1'                  , long_name = 'Current time step' , units = 'yr')
    call add_field_dp_0D(      pc%restart_filename, ncid, 'zeta_t'                  , long_name = 'Ratio between previous and new time step')
    call add_field_mesh_dp_2D( pc%restart_filename, ncid, 'dHi_dt_Hi_nm1_u_nm1'     , long_name = 'Thinning rates from previous time step', units = 'm/yr')
    call add_field_mesh_dp_2D( pc%restart_filename, ncid, 'dHi_dt_Hi_n_u_n'         , long_name = 'Thinning rates for current time step with old geometry', units = 'm/yr')
    call add_field_mesh_dp_2D( pc%restart_filename, ncid, 'Hi_star_np1'             , long_name = 'Predicted ice thickness', units = 'm')
    call add_field_mesh_dp_2D( pc%restart_filename, ncid, 'dHi_dt_Hi_star_np1_u_np1', long_name = 'Thinning rates for predicted ice thickness and updated velocity', units = 'm/yr')
    call add_field_mesh_dp_2D( pc%restart_filename, ncid, 'Hi_np1'                  , long_name = 'Corrected ice thickness', units = 'm')
    call add_field_mesh_dp_2D( pc%restart_filename, ncid, 'tau_np1'                 , long_name = 'Truncation error', units = 'm')
    call add_field_dp_0D(      pc%restart_filename, ncid, 'eta_n'                   , long_name = 'Previous maximum truncation error', units = 'm')
    call add_field_dp_0D(      pc%restart_filename, ncid, 'eta_np1'                 , long_name = 'Current maximum truncation error', units = 'm')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_pc_scheme

  subroutine remap_pc_scheme( mesh_old, mesh_new, pc)
    !< reallocate memory for the ice thickness predictor/corrector scheme.

    ! In- and output variables
    type(type_mesh),   intent(in   ) :: mesh_old
    type(type_mesh),   intent(in   ) :: mesh_new
    type(type_ice_pc), intent(inout) :: pc

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_pc_scheme'

    ! Add routine to path
    call init_routine( routine_name)

    ! reallocate memory
    call reallocate_bounds( pc%dHi_dt_Hi_nm1_u_nm1     , mesh_new%vi1, mesh_new%vi2)   ! [m/yr] Thinning rates from previous time step
    call reallocate_bounds( pc%dHi_dt_Hi_n_u_n         , mesh_new%vi1, mesh_new%vi2)   ! [m/yr] Thinning rates for current time step with old geometry
    call reallocate_bounds( pc%Hi_star_np1             , mesh_new%vi1, mesh_new%vi2)   ! [m]    Predicted ice thickness
    call reallocate_bounds( pc%dHi_dt_Hi_star_np1_u_np1, mesh_new%vi1, mesh_new%vi2)   ! [m/yr] Thinning rates for predicted ice thickness and updated velocity
    call reallocate_bounds( pc%Hi_np1                  , mesh_new%vi1, mesh_new%vi2)   ! [m]    Corrected ice thickness
    call reallocate_bounds( pc%tau_np1                 , mesh_new%vi1, mesh_new%vi2)   ! [m]    Truncation error
    call reallocate_bounds( pc%tau_n_guilty            , mesh_new%vi1, mesh_new%vi2)   ! [-]    Number of events above error tolerance

    ! Reinitialise everything from scratch
    pc%dt_n                     = C%dt_ice_min
    pc%dt_np1                   = C%dt_ice_min
    pc%zeta_t                   = 1._dp
    pc%dHi_dt_Hi_nm1_u_nm1      = 0._dp
    pc%dHi_dt_Hi_n_u_n          = 0._dp
    pc%Hi_star_np1              = 0._dp
    pc%dHi_dt_Hi_star_np1_u_np1 = 0._dp
    pc%Hi_np1                   = 0._dp
    pc%tau_np1                  = C%pc_epsilon
    pc%tau_n_guilty             = 0
    pc%eta_n                    = C%pc_epsilon
    pc%eta_np1                  = C%pc_epsilon

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_pc_scheme

end module predictor_corrector_scheme
