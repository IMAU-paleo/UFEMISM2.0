module inversion_utilities
  !< Some random utilities for inversions

  ! FIXME: group inversion code more logically!

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use model_configuration, only: C
  use region_types, only: type_model_region
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use SMB_model_types, only: type_SMB_model
  use BMB_model_types, only: type_BMB_model
  use LMB_model_types, only: type_LMB_model
  use AMB_model_types, only: type_AMB_model
  use plane_geometry, only: is_in_polygon
  use netcdf_io_main
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  use mesh_utilities, only: extrapolate_Gaussian, interpolate_to_point_dp_2D
  use conservation_of_mass_main, only: calc_dHi_dt
  use mpi_distributed_memory, only: gather_to_all
  use map_velocities_to_c_grid, only: map_velocities_from_b_to_c_2D

  implicit none

  private

  public :: initialise_dHi_dt_target, &
    MISMIPplus_adapt_flow_factor, BMB_inversion, LMB_inversion

contains

  subroutine initialise_dHi_dt_target( mesh, ice, region_name)
    !< Prescribe a target dHi_dt from a file without a time dimension

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice
    character(len=3),     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_dHi_dt_target'
    character(len=256)             :: filename_dHi_dt_target
    real(dp)                       :: timeframe_dHi_dt_target

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"!')
    case ('NAM')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_NAM
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_NAM
    case ('EAS')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_EAS
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_EAS
    case ('GRL')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_GRL
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_GRL
    case ('ANT')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_ANT
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_ANT
    end select

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    if (index( filename_dHi_dt_target,'_LAST.nc') > 1) then
      call find_last_output_file( filename_dHi_dt_target)
      call find_last_timeframe(   filename_dHi_dt_target, timeframe_dHi_dt_target)
    end if

    ! Print to terminal
    if (par%primary)  write(*,"(A)") '     Initialising target ice rates of change from file "' // colour_string( trim( filename_dHi_dt_target),'light blue') // '"...'

    ! Read dHi_dt from file
    if (timeframe_dHi_dt_target == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_file_2D( filename_dHi_dt_target, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target)
    else
      ! Assume the file has a time dimension, and read the specified timeframe
      call read_field_from_file_2D( filename_dHi_dt_target, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target, time_to_read = timeframe_dHi_dt_target)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_dHi_dt_target

  subroutine MISMIPplus_adapt_flow_factor( mesh, ice)
    !< Automatically adapt the uniform flow factor A to achieve
    ! a steady-state mid-stream grounding-line position at x = 450 km in the MISMIP+ experiment

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'MISMIPplus_adapt_flow_factor'
    real(dp), dimension(2)         :: pp,qq
    integer                        :: ti_in
    real(dp)                       :: TAFp,TAFq,lambda_GL, x_GL
    real(dp)                       :: A_flow_old, f, A_flow_new

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (C%choice_ice_rheology_Glen /= 'uniform') then
      call crash('only works in MISMIP+ geometry with a uniform flow factor!')
    end if

    ! Determine mid-channel grounding-line position
    pp = [mesh%xmin, 0._dp]
    qq = pp
    TAFp = 1._dp
    TAFq = 1._dp
    ti_in = 1
    do while (TAFp * TAFq > 0._dp)
      pp   = qq
      TAFp = TAFq
      qq = pp + [C%maximum_resolution_grounding_line, 0._dp]
      call interpolate_to_point_dp_2D( mesh, ice%TAF, qq, ti_in, TAFq)
    end do

    lambda_GL = TAFp / (TAFp - TAFq)
    x_GL = lambda_GL * qq( 1) + (1._dp - lambda_GL) * pp( 1)

    ! Adjust the flow factor
    f = 2._dp ** ((x_GL - 450E3_dp) / 80000._dp)
    C%uniform_Glens_flow_factor = C%uniform_Glens_flow_factor * f

    if (par%primary) write(0,*) '    MISMIPplus_adapt_flow_factor: x_GL = ', x_GL/1E3, ' km; changed flow factor to ', C%uniform_Glens_flow_factor

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine MISMIPplus_adapt_flow_factor

  subroutine BMB_inversion( region, dt)
    !< Invert the basal mass balance that would keep the ice shelves in equilibrium

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in   ) :: dt

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'BMB_inversion'
    integer                                              :: vi
    integer,  dimension(region%mesh%vi1:region%mesh%vi2) :: extrapolation_mask
    real(dp)                                             :: dt_dummy
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, dHi_dt_dummy, Hi_dummy
    real(dp)                                             :: w

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this inversion is desired
    if (C%do_BMB_inversion .and. &
        region%time >= C%BMB_inversion_t_start .and. &
        region%time <= C%BMB_inversion_t_end) then
      ! Go for it
    else
      ! Finalise routine path
      call finalise_routine( routine_name)
      ! And exit subroutine
      return
    end if

    ! Initialise
    extrapolation_mask = 0

    ! == Equilibrium LMB
    ! ==================

    ! Set dummy mass balance terms to 0
    SMB_dummy    = 0._dp
    BMB_dummy    = 0._dp
    LMB_dummy    = 0._dp
    AMB_dummy    = 0._dp

    ! Copy model time step
    dt_dummy = dt

    ! Use no basal or lateral mass balance to invert BMB values for an ice shelf in equilibrium
    call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, BMB_dummy, LMB_dummy, AMB_dummy, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! Initialise
    region%BMB%BMB_inv = 0._dp

    ! Compute equilibrium LMB
    do vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where BMB does not operate
      if (.not. region%ice%mask_gl_gr( vi) .and. &
          .not. region%ice%mask_floating_ice( vi) .and. &
          .not. region%ice%mask_cf_fl( vi)) cycle

      ! Equilibrium BMB field
      region%BMB%BMB_inv( vi) = -dHi_dt_dummy( vi)

      ! Add to extrapolation seeds
      extrapolation_mask( vi) = 2

    end do

    ! == Calving fronts
    ! =================

    do vi = region%mesh%vi1, region%mesh%vi2

      ! Detect shelf fronts where upstream BMB can be extrapolated into
      if (region%ice%mask_cf_fl( vi) .and. .not. region%ice%mask_gl_fl( vi)) then
        extrapolation_mask( vi) = 1
      end if

    end do

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( region%mesh, extrapolation_mask, region%BMB%BMB_inv, 10000._dp)

    ! == Total BMB
    ! ============

    ! Initialise
    region%BMB%BMB        = 0._dp
    region%BMB%BMB_transition_phase = 0._dp

    if (C%do_BMB_transition_phase) then

      ! Safety
      if (C%BMB_transition_phase_t_start < C%BMB_inversion_t_start .or. C%BMB_transition_phase_t_end > C%BMB_inversion_t_end ) then
        ! If the window of smoothing falls outside window of BMB inversion, crash.
        call crash(' The time window for BMB smoothing does not fall within the time window for BMB inversion. Make sure that "BMB_transition_phase_t_start" >= "BMB_inversion_t_start", and "BMB_transition_phase_t_end" <= "BMB_inversion_t_end".')

      elseif (C%BMB_transition_phase_t_start >= C%BMB_transition_phase_t_end) then
        ! If start and end time of smoothing window is equal or start > end, crash.
        call crash(' "BMB_transition_phase_t_start" is equivalent or larger than "BMB_transition_phase_t_end".')

      end if

      ! Compute smoothing weights for BMB inversion smoothing
      if (region%time < C%BMB_transition_phase_t_start) then
        w = 1.0_dp

      elseif (region%time >= C%BMB_transition_phase_t_start .and. &
        region%time <= C%BMB_transition_phase_t_end) then
        w = 1.0_dp - ((region%time - C%BMB_transition_phase_t_start)/(C%BMB_transition_phase_t_end - C%BMB_transition_phase_t_start))

      elseif (region%time > C%BMB_transition_phase_t_end) then
        w = 0.0_dp

      end if

    end if

    ! Compute total BMB
    do vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where BMB does not operate
      if (.not. region%ice%mask_gl_gr( vi) .and. &
          .not. region%ice%mask_floating_ice( vi) .and. &
          .not. region%ice%mask_cf_fl( vi)) cycle

      if (C%do_BMB_transition_phase) then
        ! If BMB_transition_phase is turned ON, use weight 'w' to compute BMB field
        region%BMB%BMB( vi) = w * region%BMB%BMB_inv( vi) + (1.0_dp - w) * region%BMB%BMB_modelled( vi)

        ! Save smoothed BMB field for diagnostic output
        region%BMB%BMB_transition_phase( vi) = region%BMB%BMB( vi)

      else
        ! If BMB_transition_phase is turned OFF, just apply inverted melt rates
        region%BMB%BMB( vi) = region%BMB%BMB_inv( vi)

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine BMB_inversion

end module inversion_utilities
