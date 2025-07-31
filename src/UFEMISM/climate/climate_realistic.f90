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
  use climate_matrix_utilities, only: get_insolation_at_time

  IMPLICIT NONE

  private

  public :: run_climate_model_realistic
  public :: initialise_climate_model_realistic
  public :: initialise_insolation_forcing

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

      CALL read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , mesh, C%output_dir, climate%snapshot%Hs)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , mesh, C%output_dir, climate%T2m)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh, C%output_dir, climate%Precip)


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
      !ALLOCATE(snapshot%lambda             (mesh%vi1:mesh%vi2)) ! initialised in climate matrix routine
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

      ! Read the fields at ins_t0
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, C%output_dir, snapshot%ins_Q_TOA0, time_to_read = snapshot%ins_t0)
      ! Read the fields at ins_t1
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, mesh, C%output_dir, snapshot%ins_Q_TOA1, time_to_read = snapshot%ins_t1)

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_insolation_forcing

END MODULE climate_realistic
