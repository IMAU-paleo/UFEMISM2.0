
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
    if (C%choice_matrix_forcing == 'CO2_direct') then
      call initialise_CO2_record( forcing)
    end if
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

  if (C%choice_matrix_forcing == 'CO2_direct') then
      call update_CO2_at_model_time( forcing, time)
  end if

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

  subroutine initialise_CO2_record( forcing)
    ! Read the CO2 record specified in C%filename_CO2_record. Assumes this is a file with time in yr and CO2 in ppmv

    ! NOTE: assumes time is listed in yr BP (so LGM would be -21000)
    implicit none
    
    ! In/output variables:
!    REAL(dp),                            INTENT(IN)    :: time
    type(type_global_forcing),         intent(inout)   :: forcing

    ! Local variables:
    character(LEN=256), parameter                      :: routine_name = 'initialise_CO2_record'
    integer                                            :: i,ios
    integer                                            :: CO2_record_length

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if     (C%choice_matrix_forcing == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    else
      call crash('should only be called when choice_matrix_forcing = "CO2_direct"!')
    end if
    
    if (par%primary)  write(*,"(A)") ' Initialising CO2 record '

    ! Read CO2 record (time and values) from specified text file
    if (par%primary)  write(0,*) ' Reading CO2 record from ', TRIM(C%filename_CO2_record), '...'

    call read_field_from_series_file( C%filename_CO2_record, field_name_options_CO2, forcing%CO2_record, forcing%CO2_time)

    CO2_record_length = size(forcing%CO2_record)
      if (C%start_time_of_run < forcing%CO2_time(1)) THEN
      !print *, "print value of forcing%CO2_time(first), ", forcing%CO2_time(1)
         call warning(' Model time starts before start of CO2 record; constant extrapolation will be used in that case!')
      end if
      if (C%end_time_of_run > forcing%CO2_time(CO2_record_length)) THEN
      !print *, "value of forcing%CO2_time(last), ", forcing%CO2_time(CO2_record_length)
         call warning(' Model time will reach beyond end of CO2 record; constant extrapolation will be used in that case!')
      end if
    
     !Set the value for the current (starting) model time
    call update_CO2_at_model_time( forcing, C%start_time_of_run)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_CO2_record

  subroutine update_CO2_at_model_time( forcing, time)
    ! Interpolate the data in forcing%CO2 to find the value at the queried time.
    ! If time lies outside the range of forcing%CO2_time, return the first/last value
    !
    ! NOTE: assumes time is listed in yr BP, so LGM would be -21000.0, and 0.0 corresponds to January 1st 1900.
    !
    ! NOTE: calculates average value over the preceding 30 years. For paleo this doesn't matter
    !       in the least, but for the historical period this makes everything more smooth.

    ! In/output variables:
    real(dp),                            intent(in)    :: time
    type(type_global_forcing),           intent(inout) :: forcing

    ! Local variables:
    character(LEN=256), parameter                      :: routine_name = 'update_CO2_at_model_time'
    integer                                            :: ti1, ti2, til, tiu
    real(dp)                                           :: a, b, tl, tu, intCO2, dintCO2, CO2_aux
    real(dp), parameter                                :: dt_smooth = 60._dp
    integer                                            :: CO2_record_length

    ! Add routine to path
    call init_routine( routine_name)

    !IF (par%primary)  WRITE(*,"(A)") '   Updating CO2 at model time...'

    ! Safety
    if     (C%choice_matrix_forcing == 'CO2_direct') then
      ! Observed CO2 is needed for these forcing methods.
    else
      call crash('should only be called when choice_matrix_forcing = "CO2_direct"!')
    end if

    CO2_record_length = size(forcing%CO2_record)

    if     (time < MINVAL( forcing%CO2_time)) then
      ! Model time before start of CO2 record; using constant extrapolation
      forcing%CO2_obs = forcing%CO2_record( 1)
    elseif (time > MAXVAL( forcing%CO2_time)) then
      ! Model time beyond end of CO2 record; using constant extrapolation
      forcing%CO2_obs = forcing%CO2_record( CO2_record_length)
    else

      ! Find range of raw time frames enveloping model time
      ti1 = 1
      do while (forcing%CO2_time( ti1) < time - dt_smooth .AND. ti1 < CO2_record_length)
        ti1 = ti1 + 1
      end do
      ti1 = MAX( 1, ti1 - 1)

      ti2 = 2
      do while (forcing%CO2_time( ti2) < time             .AND. ti2 < CO2_record_length)
        ti2 = ti2 + 1
      end do

      ! Calculate conservatively-remapped time-averaged CO2
      intCO2 = 0._dp
      do til = ti1, ti2 - 1
        tiu = til + 1

        ! Linear interpolation between til and tiu: CO2( t) = a + b*t
        b = (forcing%CO2_record( tiu) - forcing%CO2_record( til)) / (forcing%CO2_time( tiu) - forcing%CO2_time( til))
        a = forcing%CO2_record( til) - b*forcing%CO2_time( til)

        ! Window of overlap between [til,tiu] and [t - dt_smooth, t]
        tl = MAX( forcing%CO2_time( til), time - dt_smooth)
        tu = MIN( forcing%CO2_time( tiu), time            )
        dintCO2 = (tu - tl) * (a + b * (tl + tu) / 2._dp)
        intCO2 = intCO2 + dintCO2
      end do
      forcing%CO2_obs = intCO2 / dt_smooth

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_CO2_at_model_time

END MODULE global_forcings_main