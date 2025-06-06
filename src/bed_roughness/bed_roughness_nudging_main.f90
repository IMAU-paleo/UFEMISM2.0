module bed_roughness_nudging_main

  ! Contains all the routines for nudging the bed roughness

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use bed_roughness_model_types, only: type_basal_inversion
  use region_types, only: type_model_region
  use bed_roughness_nudging_H_dHdt_flowline, only: initialise_basal_inversion_H_dHdt_flowline, run_basal_inversion_H_dHdt_flowline

  implicit none

  private

  public :: initialise_basal_inversion, run_basal_inversion

contains

  subroutine run_basal_inversion( region)
    ! Run the main basal inversion model

    ! Input variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_basal_inversion'
    integer                       :: vi
    real(dp)                      :: wt_prev, wt_next

    ! Add routine to path
    call init_routine( routine_name)

    ! Only do basal inversion within the specified time window
    if (region%time < C%bed_roughness_nudging_t_start) then
      region%BIV%t_next = C%bed_roughness_nudging_t_start
      call finalise_routine( routine_name)
      return
    end if
    if (region%time == C%bed_roughness_nudging_t_end) then
      region%BIV%t_next = C%end_time_of_run
      call finalise_routine( routine_name)
      return
    end if

    ! If the desired time is beyond the time of the next modelled bed roughness,
    ! run the basal inversion model to calculate a new next modelled bed roughness.
    ! =============================================================================

    if (region%time == region%BIV%t_next) then
      ! Need to calculate new predicted bed roughness

      ! Store previous modelled bed roughness
      region%BIV%generic_bed_roughness_1_prev = region%BIV%generic_bed_roughness_1_next
      region%BIV%generic_bed_roughness_2_prev = region%BIV%generic_bed_roughness_2_next
      region%BIV%t_prev = region%BIV%t_next
      region%BIV%t_next = region%BIV%t_prev + C%bed_roughness_nudging_dt

      ! Run the basal inversion model to calculate a new next modelled bed roughness
      select case (C%choice_bed_roughness_nudging_method)
      case default
        call crash('unknown choice_bed_roughness_nudging_method "' // trim( C%choice_bed_roughness_nudging_method) // '"')
      case ('H_dHdt_flowline')
        ! Run with the specified target geometry
        select case (C%choice_inversion_target_geometry)
        case default
          call crash('unknown choice_inversion_target_geometry "' // trim( C%choice_inversion_target_geometry) // '"')
        case ('init')
          call run_basal_inversion_H_dHdt_flowline( region%mesh, region%grid_smooth, region%ice, region%refgeo_init, region%BIV)
        case ('PD')
          call run_basal_inversion_H_dHdt_flowline( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%BIV)
      end select

      end select

    elseif (region%time > region%BIV%t_next) then
      ! This should not be possible
      call crash('overshot the basal inversion time step')
    else
      ! We're within the current BIV prediction window
    end IF

    ! Interpolate between previous and next modelled bed roughness
    ! to find the bed roughness at the desired time
    ! =================================================================

    ! Calculate time interpolation weights
    wt_prev = (region%BIV%t_next - region%time) / (region%BIV%t_next - region%BIV%t_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled bed roughness to desired time
    do vi = region%mesh%vi1, region%mesh%vi2
      region%BIV%generic_bed_roughness_1( vi) = wt_prev * region%BIV%generic_bed_roughness_1_prev( vi) + wt_next * region%BIV%generic_bed_roughness_1_next( vi)
      region%BIV%generic_bed_roughness_2( vi) = wt_prev * region%BIV%generic_bed_roughness_2_prev( vi) + wt_next * region%BIV%generic_bed_roughness_2_next( vi)
    end do

    ! Update sliding law-specific bed roughness
    ! =========================================

    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('no_sliding')
      call crash('cannot run basal inversion for choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('idealised')
      call crash('cannot run basal inversion for choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by slid_beta_sq
      region%ice%slid_beta_sq = region%BIV%generic_bed_roughness_1
    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle
      region%ice%till_friction_angle = region%BIV%generic_bed_roughness_1
    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle
      region%ice%till_friction_angle = region%BIV%generic_bed_roughness_1
    case ('Tsai2015')
      ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
      region%ice%slid_alpha_sq = region%BIV%generic_bed_roughness_1
      region%ice%slid_beta_sq  = region%BIV%generic_bed_roughness_2
    case ('Schoof2005')
      ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
      region%ice%slid_alpha_sq = region%BIV%generic_bed_roughness_1
      region%ice%slid_beta_sq  = region%BIV%generic_bed_roughness_2
    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
      region%ice%till_friction_angle = region%BIV%generic_bed_roughness_1
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_basal_inversion

  subroutine initialise_basal_inversion( mesh, ice, BIV, region_name)
    ! Initialise the main basal inversion model

    ! Input variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(in   ) :: ice
    type(type_basal_inversion),          intent(  out) :: BIV
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_basal_inversion'

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary) write(0,*) ' Initialising basal inversion model "' // &
      colour_string( trim( C%choice_bed_roughness_nudging_method),'light blue') // '"...'

    ! Allocate memory for main variables
    ! ==================================

    allocate( BIV%generic_bed_roughness_1( mesh%vi1:mesh%vi2))
    allocate( BIV%generic_bed_roughness_2( mesh%vi1:mesh%vi2))

    BIV%generic_bed_roughness_1 = 0._dp
    BIV%generic_bed_roughness_2 = 0._dp

    allocate( BIV%generic_bed_roughness_1_prev( mesh%vi1:mesh%vi2))
    allocate( BIV%generic_bed_roughness_2_prev( mesh%vi1:mesh%vi2))
    allocate( BIV%generic_bed_roughness_1_next( mesh%vi1:mesh%vi2))
    allocate( BIV%generic_bed_roughness_2_next( mesh%vi1:mesh%vi2))

    BIV%generic_bed_roughness_1_prev = 0._dp
    BIV%generic_bed_roughness_2_prev = 0._dp
    BIV%generic_bed_roughness_1_next = 0._dp
    BIV%generic_bed_roughness_2_next = 0._dp

    ! Timeframes
    BIV%t_prev   = C%start_time_of_run
    BIV%t_next   = C%start_time_of_run

    ! Get sliding law-specific bed roughness
    ! ======================================

    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('no_sliding')
      call crash('cannot run basal inversion for choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('idealised')
      call crash('cannot run basal inversion for choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by slid_beta_sq
      BIV%generic_bed_roughness_1      = ice%slid_beta_sq
      BIV%generic_bed_roughness_1_prev = ice%slid_beta_sq
      BIV%generic_bed_roughness_1_next = ice%slid_beta_sq
      BIV%generic_bed_roughness_2      = 0._dp
      BIV%generic_bed_roughness_2_prev = 0._dp
      BIV%generic_bed_roughness_2_next = 0._dp
    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle
      BIV%generic_bed_roughness_1      = ice%till_friction_angle
      BIV%generic_bed_roughness_1_prev = ice%till_friction_angle
      BIV%generic_bed_roughness_1_next = ice%till_friction_angle
      BIV%generic_bed_roughness_2      = 0._dp
      BIV%generic_bed_roughness_2_prev = 0._dp
      BIV%generic_bed_roughness_2_next = 0._dp
    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle
      BIV%generic_bed_roughness_1      = ice%till_friction_angle
      BIV%generic_bed_roughness_1_prev = ice%till_friction_angle
      BIV%generic_bed_roughness_1_next = ice%till_friction_angle
      BIV%generic_bed_roughness_2      = 0._dp
      BIV%generic_bed_roughness_2_prev = 0._dp
      BIV%generic_bed_roughness_2_next = 0._dp
    case ('Tsai2015')
      ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
      BIV%generic_bed_roughness_1      = ice%slid_alpha_sq
      BIV%generic_bed_roughness_1_prev = ice%slid_alpha_sq
      BIV%generic_bed_roughness_1_next = ice%slid_alpha_sq
      BIV%generic_bed_roughness_2      = ice%slid_beta_sq
      BIV%generic_bed_roughness_2_prev = ice%slid_beta_sq
      BIV%generic_bed_roughness_2_next = ice%slid_beta_sq
    case ('Schoof2005')
      ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
      BIV%generic_bed_roughness_1      = ice%slid_alpha_sq
      BIV%generic_bed_roughness_1_prev = ice%slid_alpha_sq
      BIV%generic_bed_roughness_1_next = ice%slid_alpha_sq
      BIV%generic_bed_roughness_2      = ice%slid_beta_sq
      BIV%generic_bed_roughness_2_prev = ice%slid_beta_sq
      BIV%generic_bed_roughness_2_next = ice%slid_beta_sq
    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
      BIV%generic_bed_roughness_1      = ice%till_friction_angle
      BIV%generic_bed_roughness_1_prev = ice%till_friction_angle
      BIV%generic_bed_roughness_1_next = ice%till_friction_angle
      BIV%generic_bed_roughness_2      = 0._dp
      BIV%generic_bed_roughness_2_prev = 0._dp
      BIV%generic_bed_roughness_2_next = 0._dp
    end select

    ! Initialise chosen basal inversion model
    ! =======================================

    select case (C%choice_bed_roughness_nudging_method)
    case default
      call crash('unknown choice_bed_roughness_nudging_method "' // trim( C%choice_bed_roughness_nudging_method) // '"')
    case ('H_dHdt_flowline')
      call initialise_basal_inversion_H_dHdt_flowline( mesh, ice, BIV, region_name)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_basal_inversion

end module bed_roughness_nudging_main
