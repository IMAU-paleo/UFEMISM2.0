
MODULE global_forcings_main

  ! Global forcings used across regions and models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning, insert_val_into_string_int,insert_val_into_string_dp
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE global_forcing_types                                   , ONLY: type_global_forcing
  USE netcdf_io_main
  USE netcdf_basic
  USE series_utilities

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================


  SUBROUTINE initialise_global_forcings( forcing)
    ! initialise the forcing structure to get d18O, CO2, insolation, etc...

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_global_forcing),         INTENT(OUT)   :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'initialise_global_forcings'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! TODO: checks with what exactly we need to load here to know which global forcings need to be read
    ! e.g., insolation, CO2, d18O, etc...
    ! read and load the insolation data only if needed (i.e., we are using IMAU-ITM)
    !IF (C%choice_SMB_parameterised == 'IMAU-ITM') CALL initialise_insolation_forcing( forcing, mesh)

    ! CO2 record - not yet implemented
    
    ! d18O record - not yet implemented

    ! Sea level
    IF (C%choice_sealevel_model == 'prescribed') THEN
      IF (par%primary) WRITE(0,*) ' Initialising sea level record...'
      CALL initialise_sealevel_record(forcing, C%start_time_of_run)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_global_forcings



  SUBROUTINE update_global_forcings(forcing, time)
  ! Update all records such as sea level, CO2, d18O, etc...

  IMPLICIT NONE

  ! Input/output variables
  TYPE(type_global_forcing),         INTENT(INOUT)   :: forcing
  REAL(dp),                             INTENT(IN)   :: time

  ! Local variables:
  CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'update_global_forcings'

  ! Add routine to path
  CALL init_routine( routine_name)

  IF (C%choice_sealevel_model == 'prescribed') THEN
    CALL update_sealevel_at_model_time(forcing, time)
  END IF

  ! Finalise routine path
  CALL finalise_routine( routine_name)

  END SUBROUTINE update_global_forcings

  
  SUBROUTINE initialise_sealevel_record( forcing, time)
    ! Read the NetCDF file containing the prescribed sea-level curve data.

    IMPLICIT NONE

    TYPE(type_global_forcing),         INTENT(INOUT)   :: forcing
    REAL(dp),                             INTENT(IN)   :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_sealevel_record'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocating timeframe variables; the series itself is allocated in the read function below
    allocate(forcing%sl_t0)
    allocate(forcing%sl_t1)
    allocate(forcing%sl_at_t0)
    allocate(forcing%sl_at_t1)

    select case (C%choice_sealevel_model)
      case ('prescribed')
        call read_field_from_series_file( C%filename_prescribed_sealevel, field_name_options_sealevel, forcing%sea_level_record, forcing%sea_level_time)
        call update_timeframes_from_record(forcing%sea_level_time, forcing%sea_level_record, forcing%sl_t0, forcing%sl_t1, forcing%sl_at_t0,forcing%sl_at_t1, time)
      case default
        call crash('Unknown choice of sea level!')
    end select

     ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_sealevel_record


  SUBROUTINE update_sealevel_at_model_time(forcing, time)
  ! Update the current sea level based on the loaded sea level curve

    IMPLICIT NONE

    TYPE(type_global_forcing),         INTENT(INOUT)   :: forcing
    REAL(dp),                          INTENT(IN   )   :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_sealevel_at_model_time'

    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < forcing%sl_t0 .OR. time > forcing%sl_t1) THEN
      !IF (par%primary)  WRITE(0,*) '   Model time is out of the current sea level timeframes. Updating timeframes...'
      call update_timeframes_from_record(forcing%sea_level_time, forcing%sea_level_record, forcing%sl_t0, forcing%sl_t1, forcing%sl_at_t0, forcing%sl_at_t1, time)
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_sealevel_at_model_time

  SUBROUTINE update_sealevel_in_model(forcing, mesh, ice, time)
  ! Update the current sea level based on the loaded sea level curve

    IMPLICIT NONE

    TYPE(type_global_forcing),         INTENT(IN)      :: forcing
    TYPE(type_mesh),                   INTENT(IN   )   :: mesh
    TYPE(type_ice_model),              INTENT(INOUT)   :: ice
    REAL(dp),                          INTENT(IN   )   :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_sealevel_in_model'
    INTEGER                                            :: vi
    REAL(dp)                                           :: computed_sea_level

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate timeframe interpolation weights (plus safety checks for when the extend beyond the record)
    CALL interpolate_value_from_forcing_record(forcing%sl_t0, forcing%sl_t1, forcing%sl_at_t0, forcing%sl_at_t1, time, computed_sea_level) ! We might be calling this twice with no need, but might be worth keeping it like that to facilitate future impoementations

    do vi = mesh%vi1, mesh%vi2
      ice%SL( vi) = computed_sea_level
    end do
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_sealevel_in_model

END MODULE global_forcings_main