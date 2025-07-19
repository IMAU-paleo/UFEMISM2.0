module ocean_realistic

  ! Realistic ocean models

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use ocean_model_types                                      , only: type_ocean_model
  use netcdf_io_main
  use ocean_extrapolation                                    , only: extrapolate_ocean_forcing
  use series_utilities
  use ocean_deltaT
  use ocean_GlacialIndex

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine run_ocean_model_realistic( mesh, ice, ocean, time)
    ! Calculate the ocean
    !
    ! Use an realistic ocean scheme

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean
    real(dp),                               intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'run_ocean_model_realistic'

    ! Add routine to path
    call init_routine( routine_name)

    ! Run the chosen realistic ocean model
    select case (C%choice_ocean_model_realistic)
      case ('snapshot')
        ! Apply extrapolation method if required
        select case (C%choice_ocean_extrapolation_method)
          case('initialisation')
            ! nothing to do here
          case default
            call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
        end select

      case ('transient')
        call run_ocean_model_transient(mesh, ocean, time)

      case default
        call crash('unknown choice_ocean_model_realistic "' // trim( C%choice_ocean_model_realistic) // '"')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_realistic

  subroutine initialise_ocean_model_realistic( mesh, ice, ocean, region_name, start_time_of_run)
    ! Initialise the ocean model
    !
    ! Use an realistic ocean scheme

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean
    character(len=3),                       intent(in)    :: region_name
    real(dp),                               intent(in)    :: start_time_of_run

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_ocean_model_realistic'
    character(len=256)                                    :: filename_ocean_snapshot

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  write(*,"(A)") '     Initialising realistic ocean model "' // &
      colour_string( trim( C%choice_ocean_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic ocean model
    select case (C%choice_ocean_model_realistic)
      case ('snapshot')
        ! Read single-time data from external file

        select case (region_name)
          case ('NAM')
            filename_ocean_snapshot = C%filename_ocean_snapshot_NAM
          case ('EAS')
            filename_ocean_snapshot = C%filename_ocean_snapshot_EAS
          case ('GRL')
            filename_ocean_snapshot = C%filename_ocean_snapshot_GRL
          case ('ANT')
            filename_ocean_snapshot = C%filename_ocean_snapshot_ANT
          case default
            call crash('unknown region_name "' // region_name // '"')
        end select

        ! Fill in  main variables
        call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_T_ocean, mesh, C%output_dir, C%z_ocean, ocean%T)
        call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_S_ocean, mesh, C%output_dir, C%z_ocean, ocean%S)

        ! Apply extrapolation method if required
        select case (C%choice_ocean_extrapolation_method)
          case('initialisation')
            call extrapolate_ocean_forcing( mesh, ice, ocean%T)
            call extrapolate_ocean_forcing( mesh, ice, ocean%S)
          case default
            call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
        end select

      case ('transient')
          call initialise_ocean_model_transient(mesh, ice, ocean, region_name, start_time_of_run)

      case default
        call crash('unknown choice_ocean_model_realistic "' // trim( C%choice_ocean_model_realistic) // '"')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_realistic

  subroutine initialise_ocean_model_transient(mesh, ice, ocean, region_name, start_time_of_run)
  ! Initialise the ocean model
    !
    ! Use an realistic ocean scheme

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean
    character(len=3),                       intent(in)    :: region_name
    real(dp),                               intent(in)    :: start_time_of_run

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_ocean_model_transient'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_ocean_model_transient)
      case ('deltaT')
        call initialise_ocean_model_deltaT( mesh, ice, ocean, region_name, start_time_of_run)

      case ('GlacialIndex')
        call initialise_ocean_model_GlacialIndex( mesh, ice, ocean, region_name, start_time_of_run)

      case default
        call crash('unknown choice_ocean_model_transient "' // trim( C%choice_ocean_model_transient) // '"')
    end select

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine initialise_ocean_model_transient

  subroutine run_ocean_model_transient(mesh, ocean, time)
  ! Runs a transient ocean model
    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ocean_model),                 intent(inout) :: ocean
    real(dp),                               intent(in)    :: time

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_transient'

    ! Add routine to path
    CALL init_routine( routine_name)

    select case (C%choice_ocean_model_transient)
      case('deltaT')
        call run_ocean_model_deltaT(mesh, ocean, time)

      case('GlacialIndex')
        call run_ocean_model_GlacialIndex(mesh, ocean, time)

      end select
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine run_ocean_model_transient

end module ocean_realistic
