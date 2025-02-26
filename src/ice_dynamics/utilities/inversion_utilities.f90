module inversion_utilities
  !< Some random utilities for inversions

  ! FIXME: group inversion code more logically!

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string, warning
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

  public :: MB_inversion, initialise_dHi_dt_target, initialise_uabs_surf_target, &
    MISMIPplus_adapt_flow_factor, SMB_inversion, BMB_inversion, LMB_inversion

contains

  subroutine MB_inversion( mesh, ice, refgeo, SMB, BMB, LMB, AMB, dHi_dt_predicted, Hi_predicted, dt, time, region_name)
    !< Calculate the basal mass balance
    !< Use an inversion based on the computed dHi_dt

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    type(type_reference_geometry),          intent(in   ) :: refgeo
    type(type_SMB_model),                   intent(inout) :: SMB
    type(type_BMB_model),                   intent(inout) :: BMB
    type(type_LMB_model),                   intent(inout) :: LMB
    type(type_AMB_model),                   intent(inout) :: AMB
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: dHi_dt_predicted
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_predicted
    real(dp),                               intent(in   ) :: dt
    real(dp),                               intent(in   ) :: time
    character(len=3)                                      :: region_name

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'MB_inversion'
    integer                                :: vi
    logical                                :: do_BMB_inversion, do_LMB_inversion, do_SMB_inversion, do_SMB_absorb
    integer,  dimension(mesh%vi1:mesh%vi2) :: mask
    real(dp), dimension(mesh%vi1:mesh%vi2) :: previous_field
    real(dp)                               :: value_change
    real(dp), dimension(:,:), allocatable  :: poly_ROI
    real(dp), dimension(2)                 :: p

    ! == Initialisation
    ! =================

    if (C%choice_regions_of_interest == 'Patagonia') then
      ! Compute polygon for reconstruction
      call calc_polygon_Patagonia( poly_ROI)
    else
      ! Create a dummy polygon
      allocate( poly_ROI(1,2))
      poly_ROI(1,1) = 0._dp
      poly_ROI(1,2) = 0._dp
    end if

    ! Add routine to path
    call init_routine( routine_name)

    do_BMB_inversion = .false.
    do_LMB_inversion = .false.
    do_SMB_inversion = .false.
    do_SMB_absorb    = .false.

    ! Check whether we want a BMB inversion
    if (C%do_BMB_inversion .and. &
        time >= C%BMB_inversion_t_start .and. &
        time <= C%BMB_inversion_t_end) then
      do_BMB_inversion = .true.
    end if

    ! Check whether we want a LMB inversion
    if (C%do_LMB_inversion .and. &
        time >= C%LMB_inversion_t_start .and. &
        time <= C%LMB_inversion_t_end) then
      do_LMB_inversion = .true.
    end if

    if (C%do_SMB_removal_icefree_land) then
      do_SMB_inversion = .true.
    end if

    ! Check whether we want a SMB adjustment
    if (C%do_SMB_residual_absorb .and. &
        time >= C%SMB_residual_absorb_t_start .and. &
        time <= C%SMB_residual_absorb_t_end) then
      do_SMB_absorb = .true.
    end if

    ! == BMB: first pass
    ! ==================

    ! Store previous values
    previous_field = BMB%BMB_inv

    ! Initialise extrapolation mask
    mask = 0

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      ! Skip if not desired
      if (.not. do_BMB_inversion) cycle

      ! Floating calving fronts
      if (ice%mask_cf_fl( vi)) then

        ! Just mark for eventual extrapolation
        mask( vi) = 1

      elseif (ice%mask_floating_ice( vi)) then

        ! Basal melt will account for all change here
        BMB%BMB_inv( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

        ! Adjust rate of ice thickness change dHi/dt to compensate the change
        dHi_dt_predicted( vi) = 0._dp

        ! Adjust corrected ice thickness to compensate the change
        Hi_predicted( vi) = ice%Hi_prev( vi)

        ! Use this vertex as seed during extrapolation
        mask( vi) = 2

      elseif (ice%mask_gl_gr( vi)) then

        ! For grounding lines, BMB accounts for melt
        if (dHi_dt_predicted( vi) >= 0._dp) then

          ! Basal melt will account for all change here
          BMB%BMB_inv( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)
          ! Adjust rate of ice thickness change dHi/dt to compensate the change
          dHi_dt_predicted( vi) = 0._dp
          ! Adjust corrected ice thickness to compensate the change
          Hi_predicted( vi) = ice%Hi_prev( vi)

        else
          ! Remove basal melt, but do not add refreezing
          BMB%BMB_inv( vi) = 0._dp
        end if

        ! Ignore this vertex during extrapolation
        mask( vi) = 0

      else
        ! Not a place where basal melt operates
        BMB%BMB_inv( vi) = 0._dp
        ! Ignore this vertex during extrapolation
        mask( vi) = 0
      end if

    end do

    ! == Extrapolate into calving fronts
    ! ==================================

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, mask, BMB%BMB_inv, 16000._dp)

    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_cf_fl( vi)) then

        ! Get x and y coordinates of this vertex
        p = mesh%V( vi,:)

        ! Skip vertices within reconstruction polygon
        if (is_in_polygon(poly_ROI, p)) cycle

        ! Skip if not desired
        if (.not. do_BMB_inversion) cycle

        ! Compute change after extrapolation
        value_change = BMB%BMB_inv( vi) - previous_field( vi)

        ! Adjust rate of ice thickness change dHi/dt to compensate the change
        dHi_dt_predicted( vi) = dHi_dt_predicted( vi) + value_change

        ! Adjust new ice thickness to compensate the change
        Hi_predicted( vi) = ice%Hi_prev( vi) + dHi_dt_predicted( vi) * dt

      end if
    end do

    ! == LMB: remaining positive dHi_dt
    ! =================================

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      ! Skip if not desired
      if (.not. do_LMB_inversion) cycle

      ! For these areas, use dHi_dt to get an "inversion" of equilibrium LMB.
      if (ice%mask_cf_fl( vi) .OR. ice%mask_cf_gr( vi) .OR. ice%mask_icefree_ocean( vi)) then

        if (dHi_dt_predicted( vi) >= 0._dp .and. ice%fraction_margin( vi) < 1._dp) then

          ! Assume that calving accounts for all remaining mass loss here (after first BMB pass)
          LMB%LMB_inv( vi) = LMB%LMB( vi) - dHi_dt_predicted( vi)
          ! Adjust rate of ice thickness change dHi/dt to compensate the change
          dHi_dt_predicted( vi) = 0._dp
          ! Adjust corrected ice thickness to compensate the change
          Hi_predicted( vi) = ice%Hi_prev( vi)

        else
          ! Remove lateral melt, but do not add mass
          LMB%LMB_inv( vi) = 0._dp
        end if

      else
        ! Not a place where lateral melt operates
        LMB%LMB_inv( vi) = 0._dp
      end if

    end do ! vi = mesh%vi1, mesh%vi2

    ! ! == BMB: final pass
    ! ! ==================

    ! do vi = mesh%vi1, mesh%vi2
    !   if (ice%mask_cf_fl( vi)) then

    !     ! Get x and y coordinates of this vertex
    !     p = mesh%V( vi,:)

    !     ! Skip vertices within reconstruction polygon
    !     if (is_in_polygon(poly_ROI, p)) cycle

    !     ! Skip if not desired
    !     if (.not. do_BMB_inversion) cycle

    !     ! BMB will absorb all remaining change after calving did its thing
    !     BMB%BMB( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

    !     ! Adjust rate of ice thickness change dHi/dt to compensate the change
    !     dHi_dt_predicted( vi) = 0._dp

    !     ! Adjust corrected ice thickness to compensate the change
    !     Hi_predicted( vi) = ice%Hi_prev( vi)

    !   end if
    ! end do

    ! == SMB: reference ice-free land areas
    ! =====================================

    ! Store pre-adjustment values
    previous_field = SMB%SMB

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      ! Skip if not desired
      if (.not. do_SMB_inversion) cycle

      ! Skip other areas
      if (refgeo%Hb( vi) < refgeo%SL( vi) .OR. refgeo%Hi( vi) > 0._dp) cycle

      ! For reference ice-free land, use dHi_dt to get an "inversion" of equilibrium SMB.
      SMB%SMB( vi) = min( 0._dp, ice%divQ( vi))

      ! Adjust rate of ice thickness change dHi/dt to compensate the change
      dHi_dt_predicted( vi) = 0._dp

      ! Adjust corrected ice thickness to compensate the change
      Hi_predicted( vi) = ice%Hi_prev( vi)

    end do

    ! == SMB: absorb remaining change
    ! ===============================

    ! Store pre-adjustment values
    previous_field = SMB%SMB

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      if (.not. do_SMB_absorb) cycle

      ! For grounded ice, use dHi_dt to get an "inversion" of equilibrium SMB.
      SMB%SMB( vi) = SMB%SMB( vi) - dHi_dt_predicted( vi)

      ! Adjust rate of ice thickness change dHi/dt to compensate the change
      dHi_dt_predicted( vi) = 0._dp

      ! Adjust corrected ice thickness to compensate the change
      Hi_predicted( vi) = ice%Hi_prev( vi)

    end do

    ! == Assign main fields
    ! =====================

    ! Clean up after yourself
    deallocate( poly_ROI)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine MB_inversion

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

    ! Print to terminal
    if (par%master)  write(*,"(A)") '     Initialising target ice rates of change from file "' // colour_string( trim( filename_dHi_dt_target),'light blue') // '"...'

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

  subroutine initialise_uabs_surf_target( mesh, ice, region_name)
    !< Initialise surface ice velocity data from an external NetCDF file

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice
    character(len=3),     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_uabs_surf_target'
    character(len=256)             :: filename_uabs_surf_target
    real(dp)                       :: timeframe_uabs_surf_target

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename and timeframe for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"!')
    case('NAM')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_NAM
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_NAM
    case ('EAS')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_EAS
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_EAS
    case ('GRL')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_GRL
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_GRL
    case ('ANT')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_ANT
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_ANT
    end select

    ! Print to terminal
    if (par%master)  write(*,"(A)") '     Initialising target surface ice speed from file "' // colour_string( trim( filename_uabs_surf_target),'light blue') // '"...'

    if (timeframe_uabs_surf_target == 1E9_dp) then
      call read_field_from_file_2D( filename_uabs_surf_target, 'uabs_surf', mesh, ice%uabs_surf_target)
    else
      call read_field_from_file_2D( filename_uabs_surf_target, 'uabs_surf', mesh, ice%uabs_surf_target, timeframe_uabs_surf_target)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_uabs_surf_target

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

    if (par%master) write(0,*) '    MISMIPplus_adapt_flow_factor: x_GL = ', x_GL/1E3, ' km; changed flow factor to ', C%uniform_Glens_flow_factor

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine MISMIPplus_adapt_flow_factor

  subroutine SMB_inversion( region, dt)
    !< Invert the surface mass balance that would keep the ice sheet in check

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in   ) :: dt

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'SMB_inversion'
    integer                                              :: vi
    integer,  dimension(region%mesh%vi1:region%mesh%vi2) :: extrapolation_mask
    real(dp)                                             :: dt_dummy
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, dHi_dt_dummy, Hi_dummy

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this inversion is desired
    if (.not. C%do_SMB_removal_icefree_land) then
      ! Finalise routine path
      call finalise_routine( routine_name)
      ! And exit subroutine
      return
    end if

    ! Initialise
    extrapolation_mask = 0

    ! == Equilibrium SMB
    ! ==================

    ! Set dummy mass balance terms to 0
    SMB_dummy    = 0._dp
    BMB_dummy    = 0._dp
    LMB_dummy    = 0._dp
    AMB_dummy    = 0._dp

    ! Copy model time step
    dt_dummy = dt

    ! Use full mass balance to invert SMB values
    call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, AMB_dummy, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! Compute equilibrium LMB
    do vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where land should not necessarily be ice-free
      if (.not. region%ice%mask_icefree_land( vi) .OR. .not. region%refgeo_PD%Hi( vi) == 0._dp) CYCLE

      ! Equilibrium SMB field
      region%SMB%SMB( vi) = region%SMB%SMB(vi) - dHi_dt_dummy( vi)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine SMB_inversion

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
    region%BMB%BMB_smooth = 0._dp

    if (C%do_BMB_smooth_inversion) then
      ! Compute smoothing weights for BMB inversion smoothing

      ! Check whether smoothing time window is within inversion time window
      if (C%BMB_smooth_inversion_t_start < C%BMB_inversion_t_start .or. C%BMB_smooth_inversion_t_end > C%BMB_inversion_t_end ) then 
        ! If the window of smoothing falls outside window of BMB inversion, crash. 
        call crash(' The time window for BMB smoothing does not fall within the time window for BMB inversion. Make sure that "BMB_smooth_inversion_t_start" >= "BMB_inversion_t_start", and "BMB_smooth_inversion_t_end" <= "BMB_inversion_t_end".')

      if (C%BMB_smooth_inversion_t_start >= C%BMB_smooth_inversion_t_end) then
        ! If start and end time of smoothing window is equal or start > end, crash.
        call crash(' "BMB_smooth_inversion_t_start" is equivalent or larger than "BMB_smooth_inversion_t_end".')

      elseif (region%time >= C%BMB_smooth_inversion_t_start .and. &
        region%time <= C%BMB_smooth_inversion_t_end) then
        ! If region time is in smoothing window, compute smoothing weights
        w = 1.0_dp - ((region%time - C%BMB_smooth_inversion_t_start)/(C%BMB_smooth_inversion_t_end - C%BMB_smooth_inversion_t_start))

      else 
        ! If region time is outside smoothing window, do nothing

      end if

    end if

    ! Compute total BMB
    do vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where BMB does not operate
      if (.not. region%ice%mask_gl_gr( vi) .and. &
          .not. region%ice%mask_floating_ice( vi) .and. &
          .not. region%ice%mask_cf_fl( vi)) cycle

      if (C%do_BMB_smooth_inversion) then 
        ! If BMB_smooth_inversion is turned ON

        if (region%time >= C%BMB_smooth_inversion_t_start .and. &
            region%time <= C%BMB_smooth_inversion_t_end) then
          ! If time is in smoothing window: apply smoothing from inverted to modelled melt rates
          region%BMB%BMB( vi) = w * region%BMB%BMB_inv( vi) + (1.0_dp - w) * region%BMB%BMB_modelled( vi)

        elseif (region%time < C%BMB_smooth_inversion_t_start ) then 
          ! If time is before start smoothing window: use inverted melt rates
          region%BMB%BMB( vi) = region%BMB%BMB_inv( vi)

        elseif (region%time > C%BMB_smooth_inversion_t_end ) then
          ! If time is beyond end smoothing window: use modelled melt rates
          region%BMB%BMB( vi) = region%BMB%BMB_modelled( vi)

        end if

        ! Save smoothed BMB field for diagnostic output
        region%BMB%BMB_smooth( vi) = region%BMB%BMB( vi) 

      else 
        ! If BMB_smooth_inversion is turned OFF, just apply inverted melt rates
        region%BMB%BMB( vi) = region%BMB%BMB_inv( vi)

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine BMB_inversion

  subroutine LMB_inversion( region, dt)
    !< Invert the lateral mass balance that would keep the calving front in equilibrium

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in   ) :: dt

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'LMB_inversion'
    integer                                              :: vi, vj, ci, ei
    real(dp), dimension(region%mesh%ei1:region%mesh%ei2) :: u_vav_c, v_vav_c
    real(dp), dimension(region%mesh%nE)                  :: u_vav_c_tot, v_vav_c_tot
    real(dp)                                             :: dt_dummy, calving_rate, calving_ratio, calving_perp, L_c, V_calved
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, dHi_dt_dummy, Hi_dummy, divQ_eff, LMB_trans
    real(dp), dimension(region%mesh%nV)                  :: Hi_tot, fraction_margin_tot
    logical                                              :: found_advancing_calving_front, found_calving_front_neighbour
    logical,  dimension(region%mesh%vi1:region%mesh%vi2) :: mask_advancing_calving_front
    logical,  dimension(region%mesh%nV)                  :: mask_advancing_calving_front_tot, mask_cf_fl_tot, mask_icefree_ocean_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this inversion is desired
    if (C%do_LMB_inversion .and. &
        region%time >= C%LMB_inversion_t_start .and. &
        region%time <= C%LMB_inversion_t_end) then
      ! Go for it
    else
      ! Finalise routine path
      call finalise_routine( routine_name)
      ! And exit subroutine
      return
    end if

    ! == Equilibrium LMB
    ! ==================

    ! Set dummy mass balance terms to 0
    SMB_dummy    = 0._dp
    BMB_dummy    = 0._dp
    LMB_dummy    = 0._dp
    AMB_dummy    = 0._dp

    ! Copy model time step
    dt_dummy = dt

    ! Use no lateral mass balance to invert its value for a calving front in equilibrium
    call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, LMB_dummy, AMB_dummy, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! Initialise
    region%LMB%LMB_inv = 0._dp

    ! Compute equilibrium LMB
    do vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where LMB does not operate
      if (.not. region%ice%mask_cf_fl( vi) .and. &
          .not. region%ice%mask_cf_gr( vi) .and. &
          .not. region%ice%mask_icefree_ocean( vi)) CYCLE

      ! Equilibrium LMB field: let positive values remain, since the final goal is
      ! that the *sum* of the equilibrium and transient LMB be equal to the target rate
      region%LMB%LMB_inv( vi) = -dHi_dt_dummy( vi)

    end do

    ! == Transient LMB
    ! ================

    ! Gather data from all processes
    call gather_to_all(      region%ice%Hi_eff, Hi_tot)
    call gather_to_all(      region%ice%fraction_margin, fraction_margin_tot)

    ! Gather masks from all processes
    call gather_to_all( region%ice%mask_cf_fl, mask_cf_fl_tot)
    call gather_to_all( region%ice%mask_icefree_ocean, mask_icefree_ocean_tot)

    ! Calculate vertically averaged ice velocities on the edges
    call map_velocities_from_b_to_c_2D( region%mesh, region%ice%u_vav_b, region%ice%v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all( u_vav_c, u_vav_c_tot)
    call gather_to_all( v_vav_c, v_vav_c_tot)

    ! Initialise
    LMB_trans = 0._dp

    ! ! Translate imposed transient calving rates into LMB
    ! do vi = region%mesh%vi1, region%mesh%vi2

    !   V_calved = 0._dp

    !   ! Skip vertices where LMB does not operate
    !   if (.not. region%ice%mask_cf_fl( vi) .and. .not. region%ice%mask_icefree_ocean( vi)) CYCLE

    !   calving_rate = 0._dp

    !   if (C%choice_refgeo_PD_ANT == 'idealised' .and. region%time > 10000._dp) then
    !     if (C%choice_refgeo_init_idealised == 'calvmip_circular') then
    !       calving_rate = -300._dp * SIN(2._dp * pi * region%time / 1000._dp)
    !     elseif (C%choice_refgeo_init_idealised == 'calvmip_Thule') then
    !       calving_rate = -750._dp * SIN(2._dp * pi * region%time / 1000._dp)
    !     end if
    !   end if

    !   if (region%ice%uabs_vav( vi) > 0._dp) then
    !     calving_ratio = calving_rate / region%ice%uabs_vav( vi)
    !   end if

    !   ! Loop over all connections of vertex vi
    !   do ci = 1, region%mesh%nC( vi)

    !     ! Connection ci from vertex vi leads through edge ei to vertex vj
    !     ei = region%mesh%VE( vi,ci)
    !     vj = region%mesh%C(  vi,ci)

    !     ! The shared Voronoi cell boundary section between the
    !     ! Voronoi cells of vertices vi and vj has length L_c
    !     L_c = region%mesh%Cw( vi,ci)

    !     ! Calculate calving rate component perpendicular to this shared Voronoi cell boundary section
    !     calving_perp = calving_ratio * ABS( u_vav_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + v_vav_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci))

    !     ! Calving front vertices: check if neighbour is open ocean
    !     if (region%ice%mask_cf_fl( vi) .and. mask_icefree_ocean_tot( vj)) then

    !       ! Volume calved laterally: perpendicular calving rate times area of the ice front face [m^3/yr]
    !       V_calved = V_calved + L_c * calving_perp * Hi_tot( vi)

    !     ! Ice-free ocean vertices: check if neighbour is a fully advanced calving front
    !     elseif (region%ice%mask_icefree_ocean( vi) .and. mask_cf_fl_tot( vj) .and. fraction_margin_tot( vj) >= 1._dp) then

    !       ! Volume calved laterally: perpendicular calving rate times area of the ice front face [m^3/yr]
    !       V_calved = V_calved + L_c * calving_perp * Hi_tot( vj)

    !     end if

    !   end do

    !   ! Translate lateral volume loss into vertical thinning rate [m/yr]
    !   LMB_trans( vi) = V_calved / region%mesh%A( vi)

    ! end do

    ! == Total LMB
    ! ============

    ! Initialise
    region%LMB%LMB = 0._dp

    ! Compute total LMB
    do vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where LMB does not operate
      if (.not. region%ice%mask_cf_fl( vi) .and. &
          .not. region%ice%mask_cf_gr( vi) .and. &
          .not. region%ice%mask_icefree_ocean( vi)) cycle

      ! Final LMB field: now _this_ one should never be positive
      region%LMB%LMB( vi) = MIN( 0._dp, region%LMB%LMB_inv( vi) + LMB_trans( vi))

    end do

    ! ! == Effective ice divergence
    ! ! ===========================

    ! ! Set dummy mass balance terms to 0
    ! SMB_dummy    = 0._dp
    ! BMB_dummy    = 0._dp
    ! LMB_dummy    = 0._dp
    ! AMB_dummy    = 0._dp

    ! ! Copy model time step
    ! dt_dummy = dt

    ! ! Use no mass balance to get an estimate of the effective flux divergence
    ! call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, region%ice%fraction_margin, &
    !                   region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! ! Effective flux divergence
    ! divQ_eff = -dHi_dt_dummy

    ! ! Initialise mask of advancing fronts
    ! mask_advancing_calving_front = .false.

    ! ! Identify advancing calving fronts
    ! do vi = region%mesh%vi1, region%mesh%vi2
    !   if (region%ice%mask_icefree_ocean( vi) .and. divQ_eff(vi) < 0._dp) then
    !     mask_advancing_calving_front( vi) = .true.
    !   end if
    ! end do

    ! ! == Valid areas
    ! ! ==============

    ! ! Set dummy mass balance terms to 0
    ! SMB_dummy    = 0._dp
    ! BMB_dummy    = 0._dp
    ! LMB_dummy    = 0._dp
    ! AMB_dummy    = 0._dp

    ! ! Copy model time step
    ! dt_dummy = dt

    ! ! Use total mass balance to check whether the advancing calving front will survive the incoming lateral mass balance
    ! call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, AMB_dummy, region%ice%fraction_margin, &
    !                   region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! ! Check predicted dHi/dt
    ! do vi = region%mesh%vi1, region%mesh%vi2
    !   if (region%ice%mask_icefree_ocean( vi) .and. dHi_dt_dummy( vi) <= 0._dp) then

    !     ! It will not, so do not consider this point an advancing front
    !     mask_advancing_calving_front( vi) = .false.

    !     ! Apply only equilibrium LMB here
    !     region%LMB%LMB( vi) = MIN( 0._dp, region%LMB%LMB_inv( vi))

    !   end if
    ! end do

    ! ! Gather advancing and floating calving front masks from all processes
    ! call gather_to_all( mask_advancing_calving_front, mask_advancing_calving_front_tot)

    ! ! Identify vertices where LMB will operate
    ! do vi = region%mesh%vi1, region%mesh%vi2

    !   ! Valid calving front vertices
    !   if (region%ice%mask_cf_fl( vi)) then
    !     ! Initialise flag
    !     found_advancing_calving_front = .false.

    !     ! Check for advancing front neighbours
    !     do ci = 1, region%mesh%nC( vi)
    !       vj = region%mesh%C( vi,ci)
    !       if (mask_advancing_calving_front_tot( vj)) then
    !         found_advancing_calving_front = .true.
    !         exit
    !       end if
    !     end do

    !     if (found_advancing_calving_front) then
    !       ! do not apply LMB here, since it will be applied on its advancing neighbour
    !       region%LMB%LMB( vi) = 0._dp
    !     end if

    !   end if

    !   ! Valid ocean vertices
    !   if (region%ice%mask_icefree_ocean( vi)) then
    !     ! Initialise flag
    !     found_calving_front_neighbour = .false.

    !     ! Check for calving front neighbours
    !     do ci = 1, region%mesh%nC( vi)
    !       vj = region%mesh%C( vi,ci)
    !       if (mask_cf_fl_tot( vj)) then
    !         found_calving_front_neighbour = .true.
    !         exit
    !       end if
    !     end do

    !     if (.not. found_calving_front_neighbour) then
    !       ! No calving front neighbours: apply only equilibrium LMB here
    !       region%LMB%LMB( vi) = MIN( 0._dp, region%LMB%LMB_inv( vi))
    !     end if

    !   end if

    ! end do

    ! == DENK DROM
    ! ============

    ! if (C%choice_refgeo_PD_ANT == 'idealised' .and. &
    !      (C%choice_refgeo_init_idealised == 'calvmip_circular' .OR. &
    !       C%choice_refgeo_init_idealised == 'calvmip_Thule')) then

    !   if (region%time >= 6000._dp) then
    !     if (C%choice_regions_of_interest == 'CalvMIP_quarter') then
    !       C%ROI_maximum_resolution_grounding_line = 5000._dp
    !       ! C%ROI_maximum_resolution_calving_front  = 8000._dp
    !       ! C%ROI_maximum_resolution_floating_ice   = 10000._dp
    !       ! C%ROI_maximum_resolution_grounded_ice   = 20000._dp
    !     else
    !       C%maximum_resolution_grounding_line = 5000._dp
    !       ! C%maximum_resolution_calving_front  = 8000._dp
    !       ! C%maximum_resolution_floating_ice   = 10000._dp
    !     end if
    !   end if

    !   if (region%time >= 7000._dp) then
    !     if (C%choice_regions_of_interest == 'CalvMIP_quarter') then
    !       C%ROI_maximum_resolution_grounding_line = 3000._dp
    !       ! C%ROI_maximum_resolution_calving_front  = 5000._dp
    !       ! C%ROI_maximum_resolution_floating_ice   = 8000._dp
    !       ! C%ROI_maximum_resolution_grounded_ice   = 16000._dp
    !     else
    !       C%maximum_resolution_grounding_line = 3000._dp
    !       ! C%maximum_resolution_calving_front  = 5000._dp
    !       ! C%maximum_resolution_floating_ice   = 8000._dp
    !     end if
    !   end if

    !   if (region%time >= 9000._dp) then
    !     C%allow_mesh_updates = .false.
    !     if (C%choice_refgeo_init_idealised == 'calvmip_circular') then
    !       C%calving_threshold_thickness_shelf = 10._dp
    !     end if
    !   else
    !     do vi = region%mesh%vi1, region%mesh%vi2
    !       if (SQRT(region%mesh%V( vi,1)**2 + region%mesh%V( vi,2)**2) < 750000._dp) then
    !         region%LMB%LMB( vi) = 0._dp
    !       else
    !         region%LMB%LMB( vi) = -100._dp
    !       end if
    !     end do
    !   end if

    ! end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine LMB_inversion

end module inversion_utilities
