MODULE climate_realistic

  ! Realistic climate models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model, type_global_forcing
  use netcdf_io_main

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model_realistic( mesh, ice, climate, time)
    ! Calculate the climate
    !
    ! Use an realistic climate scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_realistic'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen realistic climate model
    IF     (C%choice_climate_model_realistic == 'snapshot') THEN
      ! Do nothing
    ELSEIF (C%choice_climate_model_realistic == 'climate_matrix') THEN
      ! This is probably where we will update insolation, CO2, etc...
      CALL crash('choice_climate_model_realistic climate_matrix not implemented yet!"')
    ELSE
      CALL crash('unknown choice_climate_model_realistic "' // TRIM( C%choice_climate_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_realistic

  SUBROUTINE initialise_climate_model_realistic( mesh, climate, region_name)
    ! Initialise the climate model
    !
    ! Use a realistic climate scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_realistic'
    CHARACTER(LEN=256)                                    :: filename_climate_snapshot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising realistic climate model "' // &
      colour_string( TRIM( C%choice_climate_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic climate model
    IF (C%choice_climate_model_realistic == 'snapshot') THEN
      ! Read single-time data from external file

      ! Determine which climate model to initialise for this region
      IF     (region_name == 'NAM') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_NAM
      ELSEIF (region_name == 'EAS') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_EAS
      ELSEIF (region_name == 'GRL') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_GRL
      ELSEIF (region_name == 'ANT') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_ANT
      ELSE
        CALL crash('unknown region_name "' // region_name // '"')
      END IF

      CALL read_field_from_file_2D( filename_climate_snapshot, 'Hs', mesh, climate%Hs)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m', mesh, climate%T2m)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh, climate%Precip)

      CALL initialise_global_forcings( mesh, forcing, climate) ! TODO: initialise gobal/solar forcing
      ! TODO: what should be the timeframe_init_insolation?
      CALL get_insolation_at_time( mesh, timeframe_init_insolation, forcing, climate%Q_TOA)

    ELSE
      CALL crash('unknown choice_climate_model_realistic "' // TRIM( C%choice_climate_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_realistic

  SUBROUTINE initialise_global_forcings( mesh, forcing, climate)
    ! initialise the forcing structure to get d18O, CO2, insolation, etc...

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                   INTENT(IN)    :: mesh
    TYPE(type_global_forcing),         INTENT(INOUT) :: forcing
    TYPE(type_climate_model),          INTENT(IN)    :: climate

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'initialise_global_forcings'

    ! Add routine to path
    CALL init_routine( routine_name)

    !TODO: do stuff...
    ! read and load the insolation data
    CALL initialise_insolation_forcing( forcing)

    !TODO: read and load whatever record we need/want: d18O, CO2, GI, etc...

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_global_forcings

  ! == Insolation
  SUBROUTINE initialise_insolation_forcing( forcing)
    ! initialise the insolation series in the forcing structure

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_global_forcing),         INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'initialise_insolation_forcing'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_insolation_forcing == 'none') THEN
      ! No insolation included, likely because we're running an idealised-geometry experiment
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN
      ! Initialise insolation

      ! The times at which we have insolation fields from Laskar, between which we'll interpolate
      ! to find the insolation at model time (assuming that t0 <= model_time <= t1)

      ALLOCATE( forcing%ins_t0)
      ALLOCATE( forcing%wins_t0)
      ALLOCATE( forcing%ins_t1)
      ALLOCATE( forcing%wins_t1)


      IF (par%master) THEN
        ! Give impossible values to timeframes, so that the first call to get_insolation_at_time
        ! is guaranteed to first read two new timeframes from the NetCDF file
        forcing%ins_t0 = C%start_time_of_run - 100._dp
        forcing%ins_t1 = C%start_time_of_run - 90._dp
      END IF ! IF (par%master) THEN
      CALL sync

      IF (par%master) WRITE(0,*) ' Initialising insolation data from ', TRIM(C%filename_insolation), '...'

      ! Inquire into the insolation forcing netcdf file
      ALLOCATE( forcing%ins_nyears)
      ALLOCATE( forcing%ins_nlat)

      IF (par%master) CALL inquire_insolation_file( forcing)
      CALL sync

      ! Insolation
      ALLOCATE(forcing%ins_time       ( forcing%ins_nyears))
      ALLOCATE(forcing%ins_lat        (   forcing%ins_nlat))
      ALLOCATE(forcing%Q_TOA0         (forcing%ins_nlat,12))
      ALLOCATE(forcing%Q_TOA1         (forcing%ins_nlat,12))
      
      ! Read time and latitude data
      IF (par%master) THEN

        CALL read_insolation_file_time_lat( forcing) !TODO

        IF (C%start_time_of_run < forcing%ins_time(1)) THEN
          CALL warning(' Model time starts before start of insolation record; the model will crash lol')
        END IF
        IF (C%end_time_of_run > forcing%ins_time(forcing%ins_nyears)) THEN
          CALL warning(' Model time will reach beyond end of insolation record; constant extrapolation will be used in that case!')
        END IF
      END IF

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_insoltion_forcing

  SUBROUTINE get_insolation_at_time( mesh, time, forcing, Q_TOA)
    ! Get monthly insolation at time t on the regional grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                   INTENT(IN)    :: mesh
    TYPE(type_global_forcing),         INTENT(INOUT) :: forcing
    REAL(dp),                          INTENT(IN)    :: time
    REAL(dp), DIMENSION(:,:),          INTENT(OUT)   :: Q_TOA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'get_insolation_at_time'
    REAL(dp)                                         :: time_applied
    INTEGER                                          :: vi,m,ilat_l,ilat_u
    REAL(dp)                                         :: wt0, wt1, wlat_l, wlat_u
    REAL(dp), DIMENSION(:  ), ALLOCATABLE            ::  Q_TOA_int
    INTEGER                                          :: wQ_TOA_int

    ! Add routine to path
    CALL init_routine( routine_name)

    time_applied = 0._dp

    ! Safety
    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static') THEN
      time_applied = C%static_insolation_time
    ELSEIF (C%choice_insolation_forcing == 'realistic') THEN
      time_applied = time
    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time_applied < forcing%ins_t0 .OR. time_applied > forcing%ins_t1) THEN
      CALL update_insolation_timeframes_from_file( time_applied) ! TODO
    END IF

    ! TODO: Allocate shared memory for timeframe-interpolated lat-month-only insolation
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, Q_TOA_int, wQ_TOA_int)
    

    ! Calculate timeframe interpolation weights
    wt0 = (forcing%ins_t1 - time_applied) / (forcing%ins_t1 - forcing%ins_t0)
    wt1 = 1._dp - wt0

    ! Interpolate the two timeframes
    Q_TOA_int = wt0 * forcing%ins_Q_TOA0 + wt1 * forcing%ins_Q_TOA1 ! TODO: does it need to be in par%master?
    

    ! Map the timeframe-interpolated lat-month-only insolation to the model mesh
    DO vi = mesh%vi1, mesh%vi2

      ilat_l = FLOOR(mesh%lat(vi) + 91)
      ilat_u = ilat_l + 1

      wlat_l = forcing%ins_lat(ilat_u) - mesh%lat(j,i)
      wlat_u = 1._dp - wlat_l

      DO m = 1, 12
        Q_TOA( m,vi) = wlat_l * Q_TOA_int( ilat_l,m) + wlat_u * Q_TOA_int( ilat_u,m)
      END DO
    END DO


    ! TODO: do we need the clean up below?
    CALL deallocate_shared( wQ_TOA_int)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_insolation_at_time

END MODULE climate_realistic
