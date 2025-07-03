MODULE climate_realistic

  ! Realistic climate models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning, insert_val_into_string_int,insert_val_into_string_dp
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model, type_climate_model_snapshot
  USE global_forcing_types                                   , ONLY: type_global_forcing
  USE global_forcings_main
  USE netcdf_io_main
  USE netcdf_basic
  use mpi_distributed_memory, only: distribute_from_primary

  IMPLICIT NONE

  private

  public :: run_climate_model_realistic
  public :: initialise_climate_model_realistic
  public :: initialise_global_forcings
  public :: get_insolation_at_time
  public :: update_CO2_at_model_time
  public :: update_sealevel_at_model_time
  public :: initialise_insolation_forcing
  public :: initialise_CO2_record

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model_realistic( mesh, ice, climate, forcing, time)
    ! Calculate the climate
    !
    ! Use an realistic climate scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    TYPE(type_global_forcing),              INTENT(IN)    :: forcing
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_realistic'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen realistic climate model

    ! Update temperature and precipitation fields based on the mismatch between 
    ! the ice sheet surface elevation in the forcing climate and the model's ice sheet surface elevation
    CALL apply_lapse_rate_geometry_corrections( mesh, ice, climate)

    ! if needed for IMAU-ITM or climate matrix, we need to update insolation
    IF (climate%snapshot%has_insolation) THEN
      CALL get_insolation_at_time( mesh, time, climate%snapshot)
    
      IF (C%choice_climate_model_realistic == 'climate_matrix') THEN
        ! This is probably where we will update insolation, CO2, etc...
        CALL crash('choice_climate_model_realistic climate_matrix not implemented yet!"')
        !CALL get_climate_at_time( mesh, time, forcing, climate)
      !ELSE
      !  CALL crash('unknown choice_climate_model_realistic "' // TRIM( C%choice_climate_model_realistic) // '"')
      END IF
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_realistic

  SUBROUTINE initialise_climate_model_realistic( mesh, ice, climate, forcing, region_name)
    ! Initialise the climate model
    !
    ! Use a realistic climate scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    TYPE(type_global_forcing),              INTENT(IN)    :: forcing
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_realistic'
    CHARACTER(LEN=256)                                    :: filename_climate_snapshot
    LOGICAL                                               :: do_lapse_rates
    REAL(dp)                                              :: timeframe_init_insolation

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '     Initialising realistic climate model "' // &
      colour_string( TRIM( C%choice_climate_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic climate model
    climate%snapshot%has_insolation = .FALSE. 
    IF (C%choice_climate_model_realistic == 'snapshot') THEN
      ! Read single-time data from external file

      ! Determine which climate model to initialise for this region
      IF     (region_name == 'NAM') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_NAM
        climate%snapshot%do_lapse_rates    = C%do_lapse_rate_corrections_NAM
        climate%snapshot%lapse_rate_precip = C%lapse_rate_precip_NAM
        climate%snapshot%lapse_rate_temp   = C%lapse_rate_temp_NAM
        IF (C%choice_SMB_model_NAM == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE. 
        END IF
      ELSEIF (region_name == 'EAS') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_EAS
        climate%snapshot%do_lapse_rates    = C%do_lapse_rate_corrections_EAS
        climate%snapshot%lapse_rate_precip = C%lapse_rate_precip_EAS
        climate%snapshot%lapse_rate_temp   = C%lapse_rate_temp_EAS
        IF (C%choice_SMB_model_EAS == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE. 
        END IF
      ELSEIF (region_name == 'GRL') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_GRL
        climate%snapshot%do_lapse_rates    = C%do_lapse_rate_corrections_GRL
        climate%snapshot%lapse_rate_precip = C%lapse_rate_precip_GRL
        climate%snapshot%lapse_rate_temp   = C%lapse_rate_temp_GRL
        IF (C%choice_SMB_model_GRL == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE. 
        END IF
      ELSEIF (region_name == 'ANT') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_ANT
        climate%snapshot%do_lapse_rates    = C%do_lapse_rate_corrections_ANT
        climate%snapshot%lapse_rate_precip = C%lapse_rate_precip_ANT
        climate%snapshot%lapse_rate_temp   = C%lapse_rate_temp_ANT
        IF (C%choice_SMB_model_ANT == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE. 
        END IF
      ELSE
        CALL crash('unknown region_name "' // region_name // '"')
      END IF

      CALL read_field_from_file_2D( filename_climate_snapshot, 'Hs', mesh, climate%snapshot%Hs)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m', mesh, climate%T2m)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh, climate%Precip)
      

      call apply_lapse_rate_geometry_corrections( mesh, ice, climate)

      ! Initialises the insolation (if needed)
      IF (climate%snapshot%has_insolation) THEN  
        IF (C%choice_insolation_forcing == 'none') THEN
          CALL crash('Chosen climate or SMB model cannot be used with choice_insolation_forcing = "none"!')
        ELSE
          CALL initialise_insolation_forcing( climate%snapshot, mesh)
          IF (C%start_time_of_run < 0._dp) THEN
            timeframe_init_insolation = C%start_time_of_run
          ELSE
            timeframe_init_insolation = 0._dp
          END IF
          CALL get_insolation_at_time( mesh, timeframe_init_insolation, climate%snapshot)
        END IF
      END IF

    ELSE
      CALL crash('unknown choice_climate_model_realistic "' // TRIM( C%choice_climate_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_realistic

  SUBROUTINE apply_lapse_rate_geometry_corrections( mesh, ice, climate)
    ! Applies the lapse rate corrections for temperature and precipitation
    ! to correct for the mismatch between T and P at the forcing's ice surface elevation and the model's ice surface elevation

    IMPLICIT NONE

    TYPE(type_mesh),                       INTENT(IN)    :: mesh
    TYPE(type_ice_model),                  INTENT(IN)    :: ice
    TYPE(type_climate_model),              INTENT(INOUT) :: climate

    ! Local Variables
    CHARACTER(LEN=256), PARAMETER                        :: routine_name = 'apply_lapse_rate_geometry_corrections'
    INTEGER                                              :: vi, m
    REAL(dp)                                             :: deltaH, deltaT, deltaP
    REAL(dp), DIMENSION(:,:), ALLOCATABLE                :: T_inv, T_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     ((C%choice_climate_model_realistic == 'snapshot') .AND. (climate%snapshot%do_lapse_rates .eqv. .TRUE.)) THEN

      allocate( T_inv     (mesh%vi1:mesh%vi2, 12))
      allocate( T_inv_ref (mesh%vi1:mesh%vi2, 12))

      
      do vi = mesh%vi1, mesh%vi2

        ! we only apply corrections where it is not open ocean
        if (ice%mask_icefree_ocean( vi) .eqv. .FALSE.) then
          deltaT  = (ice%Hs( vi) - climate%snapshot%Hs( vi)) * (-1._dp * abs(climate%snapshot%lapse_rate_temp))
          do m = 1, 12
            ! Do corrections - based on Eq. 11 of Albrecht et al. (2020; TC) for PISM
            climate%T2m( vi, m)    = climate%T2m( vi, m)    + deltaT
            

            ! Calculate inversion-layer temperatures
            T_inv_ref( vi, m) = 88.9_dp + 0.67_dp *  climate%T2m( vi, m)
            T_inv(     vi, m) = 88.9_dp + 0.67_dp * (climate%T2m( vi, m) - climate%snapshot%lapse_rate_temp * (ice%Hs( vi) - climate%snapshot%Hs( vi)))
            ! Correct precipitation based on a simple Clausius-Clapeyron method (Jouzel & Merlivat, 1984; Huybrechts, 2002)
            ! Same as implemented in IMAU-ICE
            climate%Precip( vi, m) = climate%Precip( vi, m) * (T_inv_ref( vi, m) / T_inv( vi, m))**2 * EXP(22.47_dp * (T0 / T_inv_ref( vi, m) - T0 / T_inv( vi, m)))
            
          end do ! m
        end if
      end do ! vi

      deallocate(T_inv)
      deallocate(T_inv_ref)
      
    ELSEIF (C%choice_climate_model_realistic == 'climate_matrix') THEN
      ! Not yet implemented! Will likely use the lambda field from Berends et al. (2018)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_lapse_rate_geometry_corrections

  ! == Insolation
  SUBROUTINE initialise_insolation_forcing( snapshot, mesh)
    ! initialise the insolation series in the forcing structure

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_climate_model_snapshot),   INTENT(INOUT) :: snapshot
    TYPE(type_mesh),            INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'initialise_insolation_forcing'
    CHARACTER(LEN=256)                               :: str
    INTEGER                                          :: ncid
    REAL(dp)                                         :: closest_t0

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_insolation_forcing == 'none') THEN
      ! No insolation included, likely because we're running an idealised-geometry experiment
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN
      
      ! Initialise insolation
      ! The times at which we have insolation fields from Laskar, between which we'll interpolate
      ! to find the insolation at model time (assuming that t0 <= model_time <= t1)
      IF (par%primary)   WRITE(0,*) ' Initialising insolation data from ', TRIM(C%filename_insolation), '...'

      ! Memory allocation
      ALLOCATE(snapshot%ins_t0)
      ALLOCATE(snapshot%ins_t1)
      ALLOCATE(snapshot%ins_ti0)
      ALLOCATE(snapshot%ins_ti1)
      ALLOCATE(snapshot%ins_nlat)
      ALLOCATE(snapshot%ins_nlon)
      ALLOCATE(snapshot%ins_lat            (   snapshot%ins_nlat))
      ALLOCATE(snapshot%ins_Q_TOA0         (mesh%vi1:mesh%vi2,12))
      ALLOCATE(snapshot%ins_Q_TOA1         (mesh%vi1:mesh%vi2,12))
      !ALLOCATE(snapshot%lambda             (mesh%vi1:mesh%vi2))
      ALLOCATE(snapshot%Q_TOA              (mesh%vi1:mesh%vi2,12))
      ALLOCATE(snapshot%Albedo             (mesh%vi1:mesh%vi2,12))
      ALLOCATE(snapshot%I_abs              (mesh%vi1:mesh%vi2))
      snapshot%ins_t0     = C%start_time_of_run
      snapshot%ins_t1     = C%start_time_of_run
      snapshot%ins_nlat   = 181
      snapshot%ins_nlon   = 360
      snapshot%ins_lat    = 0._dp
      snapshot%ins_Q_TOA0 = 0._dp
      snapshot%ins_Q_TOA1 = 0._dp
      snapshot%Q_TOA      = 0._dp
      
      ! find the closest timeframe to the start of the run
      call read_field_from_file_0D( C%filename_insolation, field_name_options_time, closest_t0, time_to_read = C%start_time_of_run)

      if (C%start_time_of_run >= closest_t0) then
        if (par%primary) WRITE(0,*) '     start time is after ins_t0, reading one step further for ins_t1...'
        snapshot%ins_t0 = closest_t0
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t1, time_to_read = C%start_time_of_run+1000._dp)
      else
        ! otherwise we read one record before for t1
        if (par%primary) WRITE(0,*) '     start time is before closest ins_t0, reading one step earlier for t0, and using that one for t1...'
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t0, time_to_read = C%start_time_of_run-1000._dp)
        snapshot%ins_t1 = closest_t0
      end if
print *, "values of closest_t0 ", closest_t0, "and snapshot%ins_t1 ", snapshot%ins_t1
      ! Read the fields at ins_t0
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, snapshot%ins_Q_TOA0, time_to_read = snapshot%ins_t0)
      ! Read the fields at ins_t1
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, snapshot%ins_Q_TOA1, time_to_read = snapshot%ins_t1)

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_insolation_forcing

  SUBROUTINE get_insolation_at_time( mesh, time, snapshot)
    ! Get monthly insolation at time t on the regional grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model_snapshot),      INTENT(INOUT) :: snapshot
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'get_insolation_at_time'
    REAL(dp)                                         :: time_applied
    INTEGER                                          :: vi,m 
    REAL(dp)                                         :: wt0, wt1

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
      print *, "value of time in get_insolation_at_time", time_applied
    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time_applied < snapshot%ins_t0 .OR. time_applied > snapshot%ins_t1) THEN
      IF (par%primary)  WRITE(0,*) '   Model time is out of the current insolation timeframes. Updating timeframes...'
      CALL update_insolation_timeframes_from_file( snapshot, time_applied, mesh)
    END IF

    ! Calculate timeframe interpolation weights (plus safety checks for when the extend beyond the record)
    if (snapshot%ins_t1 == snapshot%ins_t0) then
      wt0 = 0._dp
      wt1 = 1._dp
    else
      if (time_applied > snapshot%ins_t1) then
        wt0 = 0._dp
      elseif (time_applied < snapshot%ins_t0) then
        wt0 = 1._dp
      else
        wt0 = (snapshot%ins_t1 - time_applied) / (snapshot%ins_t1 - snapshot%ins_t0)
      end if
      wt1 = 1._dp - wt0
    end if
print *, "value of wt0 ", wt0, "  and wt1 ", wt1
print *, "sum of Q_TOA0 ", sum(snapshot%ins_Q_TOA0), "and Q_TOA1 ", sum(snapshot%ins_Q_TOA1)
    ! Interpolate the two timeframes
    do vi = mesh%vi1, mesh%vi2
      do m = 1, 12
      !print *, "value of Q_TOA0 ", snapshot%ins_Q_TOA0(vi, m), "and Q_TOA1 ", snapshot%ins_Q_TOA1( vi, m)
        snapshot%Q_TOA(vi, m) = wt0 * snapshot%ins_Q_TOA0(vi, m) + wt1 * snapshot%ins_Q_TOA1(vi, m)
      end do
    end do

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_insolation_at_time

  SUBROUTINE update_insolation_timeframes_from_file( snapshot, time, mesh)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.

    ! NOTE: assumes time in forcing file is in kyr

    IMPLICIT NONE

    TYPE(type_mesh),                  INTENT(IN)     :: mesh
    TYPE(type_climate_model_snapshot), INTENT(INOUT)  :: snapshot
    REAL(dp),                         INTENT(IN)     :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_insolation_timeframes_from_file'
    INTEGER                                            :: ti0, ti1, ncid
    CHARACTER(LEN=256)                                 :: str

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN

      ! Update insolation
      ! Find time indices to be read
      !IF (par%primary) THEN

        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t0, time_to_read = time)
        
        ! if the desired time is after t0, we read one record after for t1
        if (time >= snapshot%ins_t0) then
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t1, time_to_read = time+1000._dp)
        else
        ! otherwise we read one record before for t0, and that record becomes t1
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t1, time_to_read = time)
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, snapshot%ins_t0, time_to_read = time-1000._dp)
        end if

      !END IF ! IF (par%primary) THEN
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, snapshot%ins_Q_TOA0, time_to_read = snapshot%ins_t0)
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, snapshot%ins_Q_TOA1, time_to_read = snapshot%ins_t1)

      call warning('insolation timeframes at t = {dp_01} are ins_t0={dp_02} and ins_t1={dp_03}', dp_01 =  time, dp_02 = snapshot%ins_t0, dp_03 = snapshot%ins_t1)

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_insolation_timeframes_from_file
  
! == Prescribed CO2 record
  SUBROUTINE update_CO2_at_model_time( time, forcing)
    ! Interpolate the data in forcing%CO2 to find the value at the queried time.
    ! If time lies outside the range of forcing%CO2_time, return the first/last value
    !
    ! NOTE: assumes time is listed in yr BP, so LGM would be -21000.0, and 0.0 corresponds to January 1st 1900.
    !
    ! NOTE: calculates average value over the preceding 30 years. For paleo this doesn't matter
    !       in the least, but for the historical period this makes everything more smooth.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_global_forcing),         INTENT(INOUT)   :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_CO2_at_model_time'
    INTEGER                                            :: ti1, ti2, til, tiu
    REAL(dp)                                           :: a, b, tl, tu, intCO2, dintCO2, CO2_aux
    REAL(dp), PARAMETER                                :: dt_smooth = 60._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_matrix_forcing == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_matrix_forcing = "CO2_direct"!')
    END IF

    !IF (par%primary) THEN

      IF     (time < MINVAL( forcing%CO2_time)) THEN
        ! Model time before start of CO2 record; using constant extrapolation
        forcing%CO2_obs = forcing%CO2_record( 1)
      ELSEIF (time > MAXVAL( forcing%CO2_time)) THEN
        ! Model time beyond end of CO2 record; using constant extrapolation
        forcing%CO2_obs = forcing%CO2_record( C%CO2_record_length)
      ELSE

        ! Find range of raw time frames enveloping model time
        ti1 = 1
        DO WHILE (forcing%CO2_time( ti1) < time - dt_smooth .AND. ti1 < C%CO2_record_length)
          ti1 = ti1 + 1
        END DO
        ti1 = MAX( 1, ti1 - 1)

        ti2 = 2
        DO WHILE (forcing%CO2_time( ti2) < time             .AND. ti2 < C%CO2_record_length)
          ti2 = ti2 + 1
        END DO

        ! Calculate conservatively-remapped time-averaged CO2
        intCO2 = 0._dp
        DO til = ti1, ti2 - 1
          tiu = til + 1

          ! Linear interpolation between til and tiu: CO2( t) = a + b*t
          b = (forcing%CO2_record( tiu) - forcing%CO2_record( til)) / (forcing%CO2_time( tiu) - forcing%CO2_time( til))
          a = forcing%CO2_record( til) - b*forcing%CO2_time( til)

          ! Window of overlap between [til,tiu] and [t - dt_smooth, t]
          tl = MAX( forcing%CO2_time( til), time - dt_smooth)
          tu = MIN( forcing%CO2_time( tiu), time            )
          dintCO2 = (tu - tl) * (a + b * (tl + tu) / 2._dp)
          intCO2 = intCO2 + dintCO2
        END DO
        forcing%CO2_obs = intCO2 / dt_smooth

      END IF
      !CO2_aux = forcing%CO2_obs
      !call MPI_BCAST(forcing%CO2_obs, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !END IF
    !call distribute_from_primary_dp_1D(forcing%CO2_obs)
    ! not working, check in the code of ice_dynamics_main what is distributing with this, a real? or more like a vector...
    ! compare with the things that I asked to chatgpt before, doing this it calculates for every core. w_CO2 = -0.5 now. 
    print *, "print value of forcing%CO2_obs in climate realistic...", forcing%CO2_obs
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_CO2_at_model_time
  SUBROUTINE initialise_CO2_record( forcing)
    ! Read the CO2 record specified in C%filename_CO2_record. Assumes this is an ASCII text file with at least two columns (time in kyr and CO2 in ppmv)
    ! and the number of rows being equal to C%CO2_record_length

    ! NOTE: assumes time is listed in kyr BP (so LGM would be -21.0); converts to yr after reading!

    IMPLICIT NONE
    
    ! In/output variables:
!    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_global_forcing),         INTENT(INOUT)   :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_CO2_record'
    INTEGER                                            :: i,ios

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_matrix_forcing == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_matrix_forcing = "CO2_direct"!')
    END IF
    
      IF (par%primary)  WRITE(*,"(A)") ' Initialising CO2 record '

    ! Allocate shared memory to take the data
    allocate( forcing%CO2_time(   C%CO2_record_length))
    allocate( forcing%CO2_record( C%CO2_record_length))
    !CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_time,   forcing%wCO2_time  )
    !CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_record, forcing%wCO2_record)
    !CALL allocate_shared_dp_0D(                      forcing%CO2_obs,    forcing%wCO2_obs   )

    ! Read CO2 record (time and values) from specified text file
    IF (par%primary)  WRITE(0,*) ' Reading CO2 record from ', TRIM(C%filename_CO2_record), '...'
!    IF (par%primary) THEN

!      WRITE(0,*) ' Reading CO2 record from ', TRIM(C%filename_CO2_record), '...'
!    END IF
! check this?!
! I added field_name_options_CO2 in netcdf field list
! from the funciton 3rd and 4th are outputs

!! HERE IDK IF IS NEEDED TO CALL THEM INSIDE PRIMARY OR NOT.. CHECK
      call read_field_from_series_file( C%filename_CO2_record, field_name_options_CO2, forcing%CO2_record, forcing%CO2_time)

!      OPEN(   UNIT = 1337, FILE=C%filename_CO2_record, ACTION='READ')

!      DO i = 1, C%CO2_record_length
!        READ( UNIT = 1337, FMT=*, IOSTAT=ios) forcing%CO2_time( i), forcing%CO2_record( i)
!        IF (ios /= 0) THEN
!          CALL crash('length of text file "' // TRIM(C%filename_CO2_record) // '" does not match C%CO2_record_length!')
!        END IF
!      END DO

!      CLOSE( UNIT  = 1337)
!    IF (par%primary) THEN

      IF (C%start_time_of_run/1000._dp < forcing%CO2_time(1)) THEN
         CALL warning(' Model time starts before start of CO2 record; constant extrapolation will be used in that case!')
      END IF
      IF (C%end_time_of_run/1000._dp > forcing%CO2_time(C%CO2_record_length)) THEN
         CALL warning(' Model time will reach beyond end of CO2 record; constant extrapolation will be used in that case!')
      END IF

      ! Convert from kyr to yr
      forcing%CO2_time = forcing%CO2_time * 1000._dp

!    END IF ! IF (par%primary)
    
    IF (par%primary)  WRITE(*,"(A)") '   Updating CO2 at model time...'
    ! Set the value for the current (starting) model time
    CALL update_CO2_at_model_time( C%start_time_of_run, forcing)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_CO2_record


END MODULE climate_realistic
