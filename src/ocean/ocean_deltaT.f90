module ocean_deltaT

  ! deltaT ocean model (i.e., T = T0+dT)

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

  implicit none

contains

! ===== Main routines =====
! =========================

subroutine initialise_ocean_model_deltaT( mesh, ice, ocean, region_name, start_time_of_run)
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
    character(len=256), parameter                         :: routine_name = 'initialise_ocean_model_deltaT'
    character(len=256)                                    :: filename_ocean_snapshot, filename_ocean_dT

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary)  write(*,"(A)") '     Initialising transient ocean model "' // &
            colour_string( trim( C%choice_ocean_model_transient),'light blue') // '"...'
            ! We need the snapshot and the dT to apply to it
            select case (region_name)
              case ('NAM')
                filename_ocean_snapshot = C%filename_ocean_snapshot_NAM
                filename_ocean_dT       = C%filename_ocean_dT_NAM
              case ('EAS')
                filename_ocean_snapshot = C%filename_ocean_snapshot_EAS
                filename_ocean_dT       = C%filename_ocean_dT_EAS
              case ('GRL')
                filename_ocean_snapshot = C%filename_ocean_snapshot_GRL
                filename_ocean_dT       = C%filename_ocean_dT_GRL
              case ('ANT')
                filename_ocean_snapshot = C%filename_ocean_snapshot_ANT
                filename_ocean_dT       = C%filename_ocean_dT_ANT
              case default
                call crash('unknown region_name "' // region_name // '"')
            end select

            ! Allocating timeframe variables; the series itself is allocated in the read function below
            allocate(ocean%deltaT%dT_t0)
            allocate(ocean%deltaT%dT_t1)
            allocate(ocean%deltaT%dT_at_t0)
            allocate(ocean%deltaT%dT_at_t1)
            allocate( ocean%deltaT%T0( mesh%vi1:mesh%vi2,C%nz_ocean))
            allocate( ocean%deltaT%S0( mesh%vi1:mesh%vi2,C%nz_ocean))
            ocean%deltaT%T0 = 0._dp
            ocean%deltaT%S0 = 0._dp

            ! Fill in  main variables
            call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_T_ocean,  mesh, C%output_dir, C%z_ocean, ocean%deltaT%T0)
            call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_S_ocean,  mesh, C%output_dir, C%z_ocean, ocean%deltaT%S0)

            call read_field_from_series_file(   filename_ocean_dT,       field_name_options_dT_ocean, ocean%deltaT%dT_series, ocean%deltaT%dT_series_time)
            call update_timeframes_from_record(ocean%deltaT%dT_series_time, ocean%deltaT%dT_series, ocean%deltaT%dT_t0, ocean%deltaT%dT_t1, ocean%deltaT%dT_at_t0, ocean%deltaT%dT_at_t1, start_time_of_run)

            ! Apply extrapolation method if required
            select case (C%choice_ocean_extrapolation_method)
              case('initialisation')
                call extrapolate_ocean_forcing( mesh, ice, ocean%deltaT%T0)
                call extrapolate_ocean_forcing( mesh, ice, ocean%deltaT%S0)
              case default
                call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
            end select

        call finalise_routine(routine_name)

    end subroutine initialise_ocean_model_deltaT

    subroutine run_ocean_model_deltaT(mesh, ocean, time)
    ! Apply the dT anomalies to the snapshot ocean
    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ocean_model),                 intent(inout) :: ocean
    real(dp),                               intent(in)    :: time

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_deltaT'
    REAL(dp)                                           :: dT_at_time
    INTEGER                                            :: vi, z

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (time < ocean%deltaT%dT_t0 .OR. time > ocean%deltaT%dT_t1) THEN
          !IF (par%primary)  WRITE(0,*) '   Model time is out of the current dT timeframes. Updating timeframes...'
          call update_timeframes_from_record(ocean%deltaT%dT_series_time, ocean%deltaT%dT_series, ocean%deltaT%dT_t0, ocean%deltaT%dT_t1, ocean%deltaT%dT_at_t0, ocean%deltaT%dT_at_t1, time)
        END IF

        ! Interpolate the two timeframes - constant dT over the entire region
        call interpolate_value_from_forcing_record(ocean%deltaT%dT_t0, ocean%deltaT%dT_t1, ocean%deltaT%dT_at_t0, ocean%deltaT%dT_at_t1, time, dT_at_time)

        do vi = mesh%vi1, mesh%vi2
          do z = 1, C%nz_ocean
            ocean%T(vi, z) = ocean%deltaT%T0(vi, z) + dT_at_time
          end do
        end do

    call finalise_routine(routine_name)

    end subroutine run_ocean_model_deltaT




end module ocean_deltaT