module bed_roughness_nudging_main

  ! Contains all the routines for nudging the bed roughness

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: crash, warning, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use bed_roughness_model_types, only: type_bed_roughness_model
  use region_types, only: type_model_region
  use bed_roughness_nudging_H_dHdt_flowline, only: initialise_bed_roughness_nudging_H_dHdt_flowline, run_bed_roughness_nudging_H_dHdt_flowline
  use bed_roughness_nudging_H_dHdt_local, only: initialise_bed_roughness_nudging_H_dHdt_local, run_bed_roughness_nudging_H_dHdt_local
  use bed_roughness_nudging_H_u_flowline, only: initialise_bed_roughness_nudging_H_u_flowline, run_bed_roughness_nudging_H_u_flowline

  implicit none

  private

  public :: initialise_bed_roughness_nudging_model, run_bed_roughness_nudging_model

contains

  subroutine run_bed_roughness_nudging_model( region)
    ! Run the main bed roughness nudging model

    ! Input variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_bed_roughness_nudging_model'
    integer                       :: vi
    real(dp)                      :: wt_prev, wt_next

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%do_bed_roughness_nudging) then
      call finalise_routine( routine_name)
      return
    end if

    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('Tsai2015','Schoof2005','no_sliding','idealised')
      call crash('Bed roughness nudging not supported for choice_sliding_law "' // &
        trim(C%choice_sliding_law) // '"')
    case ('Weertman','Coulomb','Budd','Zoet-Iverson')
      ! Nudging should work for these sliding laws
    end select

    ! Only do basal inversion within the specified time window
    if (region%time < C%bed_roughness_nudging_t_start) then
      region%bed_roughness%t_next = C%bed_roughness_nudging_t_start
      call finalise_routine( routine_name)
      return
    end if
    if (region%time == C%bed_roughness_nudging_t_end) then
      region%bed_roughness%t_next = C%end_time_of_run
      call finalise_routine( routine_name)
      return
    end if

    ! If the desired time is beyond the time of the next modelled bed roughness,
    ! run the basal inversion model to calculate a new next modelled bed roughness.
    ! =============================================================================

    if (region%time == region%bed_roughness%t_next) then
      ! Need to calculate new predicted bed roughness

      ! Store previous modelled bed roughness
      region%bed_roughness%generic_bed_roughness_prev = region%bed_roughness%generic_bed_roughness_next
      region%bed_roughness%t_prev = region%bed_roughness%t_next
      region%bed_roughness%t_next = region%bed_roughness%t_prev + C%bed_roughness_nudging_dt

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
          call run_bed_roughness_nudging_H_dHdt_flowline( region%mesh, region%grid_smooth, region%ice, region%refgeo_init, region%bed_roughness)
        case ('PD')
          call run_bed_roughness_nudging_H_dHdt_flowline( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%bed_roughness)
        end select
      case ('H_dHdt_local')
        ! Run with the specified target geometry
        select case (C%choice_inversion_target_geometry)
        case default
          call crash('unknown choice_inversion_target_geometry "' // trim( C%choice_inversion_target_geometry) // '"')
        case ('init')
          call run_bed_roughness_nudging_H_dHdt_local( region%mesh, region%grid_smooth, &
            region%ice, region%refgeo_init, region%bed_roughness%generic_bed_roughness_prev, &
            region%bed_roughness%generic_bed_roughness_next, region%bed_roughness%nudging_H_dHdt_local)
        case ('PD')
          call run_bed_roughness_nudging_H_dHdt_local( region%mesh, region%grid_smooth, &
            region%ice, region%refgeo_PD  , region%bed_roughness%generic_bed_roughness_prev, &
            region%bed_roughness%generic_bed_roughness_next, region%bed_roughness%nudging_H_dHdt_local)
        end select
      case ('H_u_flowline')
        ! Run with the specified target geometry
        select case (C%choice_inversion_target_geometry)
        case default
          call crash('unknown choice_inversion_target_geometry "' // trim( C%choice_inversion_target_geometry) // '"')
        case ('init')
          call run_bed_roughness_nudging_H_u_flowline( region%mesh, region%ice, region%refgeo_init, region%bed_roughness)
        case ('PD')
          call run_bed_roughness_nudging_H_u_flowline( region%mesh, region%ice, region%refgeo_PD, region%bed_roughness)
        end select
      end select

    elseif (region%time > region%bed_roughness%t_next) then
      ! This should not be possible
      call crash('overshot the basal inversion time step')
    else
      ! We're within the current BIV prediction window
    end IF

    ! Interpolate between previous and next modelled bed roughness
    ! to find the bed roughness at the desired time
    ! =================================================================

    ! Calculate time interpolation weights
    wt_prev = (region%bed_roughness%t_next - region%time) / (region%bed_roughness%t_next - region%bed_roughness%t_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled bed roughness to desired time
    do vi = region%mesh%vi1, region%mesh%vi2
      region%bed_roughness%generic_bed_roughness( vi) = wt_prev * region%bed_roughness%generic_bed_roughness_prev( vi) + wt_next * region%bed_roughness%generic_bed_roughness_next( vi)
    end do

    ! Update sliding law-specific bed roughness
    ! =========================================

    select case (C%choice_sliding_law)
    case default
      call crash('unknown/unsupported choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by beta_sq
      region%bed_roughness%beta_sq = region%bed_roughness%generic_bed_roughness
    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle
      region%bed_roughness%till_friction_angle = region%bed_roughness%generic_bed_roughness
    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle
      region%bed_roughness%till_friction_angle = region%bed_roughness%generic_bed_roughness
    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
      region%bed_roughness%till_friction_angle = region%bed_roughness%generic_bed_roughness
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_bed_roughness_nudging_model

  subroutine initialise_bed_roughness_nudging_model( mesh, bed_roughness)
    ! Initialise the main bed roughness nudging model

    ! Input variables:
    type(type_mesh),                intent(in   ) :: mesh
    type(type_bed_roughness_model), intent(inout) :: bed_roughness

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_nudging_model'

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%do_bed_roughness_nudging) then
      call finalise_routine( routine_name)
      return
    end if

    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('Tsai2015','Schoof2005','no_sliding','idealised')
      call crash('Bed roughness nudging not supported for choice_sliding_law "' // &
        trim(C%choice_sliding_law) // '"')
    case ('Weertman','Coulomb','Budd','Zoet-Iverson')
      ! Nudging should work for these sliding laws
    end select

    ! Print to terminal
    if (par%primary) write(0,*) ' Initialising bed roughness nudging model "' // &
      colour_string( trim( C%choice_bed_roughness_nudging_method),'light blue') // '"...'

    ! Allocate memory for main variables
    ! ==================================

    allocate( bed_roughness%generic_bed_roughness     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( bed_roughness%generic_bed_roughness_prev( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( bed_roughness%generic_bed_roughness_next( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Timeframes
    bed_roughness%t_prev   = C%start_time_of_run
    bed_roughness%t_next   = C%start_time_of_run

    ! Get sliding law-specific bed roughness
    ! ======================================

    select case (C%choice_sliding_law)
    case default
      call crash('unknown/unsupported choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by beta_sq
      bed_roughness%generic_bed_roughness      = bed_roughness%beta_sq
      bed_roughness%generic_bed_roughness_prev = bed_roughness%beta_sq
      bed_roughness%generic_bed_roughness_next = bed_roughness%beta_sq
    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle
      bed_roughness%generic_bed_roughness      = bed_roughness%till_friction_angle
      bed_roughness%generic_bed_roughness_prev = bed_roughness%till_friction_angle
      bed_roughness%generic_bed_roughness_next = bed_roughness%till_friction_angle
    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle
      bed_roughness%generic_bed_roughness      = bed_roughness%till_friction_angle
      bed_roughness%generic_bed_roughness_prev = bed_roughness%till_friction_angle
      bed_roughness%generic_bed_roughness_next = bed_roughness%till_friction_angle
    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
      bed_roughness%generic_bed_roughness      = bed_roughness%till_friction_angle
      bed_roughness%generic_bed_roughness_prev = bed_roughness%till_friction_angle
      bed_roughness%generic_bed_roughness_next = bed_roughness%till_friction_angle
    end select

    ! Initialise chosen basal inversion model
    ! =======================================

    select case (C%choice_bed_roughness_nudging_method)
    case default
      call crash('unknown choice_bed_roughness_nudging_method "' // trim( C%choice_bed_roughness_nudging_method) // '"')
    case ('H_dHdt_flowline')
      call initialise_bed_roughness_nudging_H_dHdt_flowline( mesh, bed_roughness%nudging_H_dHdt_flowline)
    case ('H_dHdt_local')
      call initialise_bed_roughness_nudging_H_dHdt_local( mesh, bed_roughness%nudging_H_dHdt_local)
    case ('H_u_flowline')
      call initialise_bed_roughness_nudging_H_u_flowline( mesh, bed_roughness%nudging_H_u_flowline)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_nudging_model

end module bed_roughness_nudging_main
