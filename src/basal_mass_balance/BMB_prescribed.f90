MODULE BMB_prescribed

  ! Prescribed BMB forcing

  ! Prescribed BMB forcing. Read in a fixed 2D field of BMB values. At each time step, the applied BMB is determined by the chosen
  ! subgrid melt scheme. So in case of FCMP, for example, BMB_gr will be zero throughout the run.
  !
  ! If the BMB model = 'prescribed_fixed', the subgrid melt schem is only applied during initialisation. Hence, the integrated
  ! metric BMB_total will be constant throughout the run, regardless of GL migration or calving. GL advance will lead to nonzero
  ! BMB_gr, which may be useful to maintain the GL position. 

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, warning, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE netcdf_input                                           , ONLY: read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model_prescribed( mesh, ice, BMB, region_name, time)
    ! Calculate the basal mass balance

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_prescribed'
    REAL(dp)                                              :: dummy_dp
    CHARACTER                                             :: dummy_char
    CHARACTER(LEN=256)                                    :: choice_BMB_prescribed

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy_dp   = mesh%xmin
    dummy_dp   = ice%Hi( mesh%vi1)
    dummy_dp   = BMB%BMB_shelf( mesh%vi1)
    dummy_char = region_name( 1:1)
    dummy_dp   = time

    ! Determine the type of prescribed BMB forcing for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_NAM
      CASE ('EAS')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_EAS
      CASE ('GRL')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_GRL
      CASE ('ANT')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Run the chosen type of prescribed BMB forcing
    SELECT CASE (choice_BMB_prescribed)
      CASE ('BMB_no_time')
        ! BMB only, no time
        CALL run_BMB_model_prescribed_notime( mesh, BMB)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_prescribed "' // TRIM( choice_BMB_prescribed) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_prescribed

  SUBROUTINE initialise_BMB_model_prescribed( mesh, BMB, region_name)
    ! Initialise the BMB model
    !
    ! Prescribed BMB forcing

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_prescribed'
    CHARACTER(LEN=256)                                    :: choice_BMB_prescribed

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine the type of prescribed BMB forcing for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_NAM
      CASE ('EAS')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_EAS
      CASE ('GRL')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_GRL
      CASE ('ANT')
        choice_BMB_prescribed  = C%choice_BMB_prescribed_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Initialised the chosen type of prescribed BMB forcing
    SELECT CASE (choice_BMB_prescribed)
      CASE ('BMB_no_time')
        CALL initialise_BMB_model_prescribed_notime( mesh, BMB, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_prescribed "' // TRIM( choice_BMB_prescribed) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_prescribed

  ! == BMB only, no time
  ! ====================

  SUBROUTINE run_BMB_model_prescribed_notime( mesh, BMB)
    ! Calculate the basal mass balance
    !
    ! Prescribed BMB forcing
    !
    ! BMB only, no time

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_prescribed_notime'
    REAL(dp)                                              :: dummy_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy_dp = mesh%xmin
    dummy_dp = BMB%BMB_shelf( mesh%vi1)

    ! No need to do anything, as the BMB was already read during initialisation

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_prescribed_notime

  SUBROUTINE initialise_BMB_model_prescribed_notime( mesh, BMB, region_name)
    ! Initialise the BMB model
    !
    ! Prescribe BMB from a file without a time dimension
    !
    ! BMB only, no time

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_prescribed_notime'
    CHARACTER(LEN=256)                                    :: filename_BMB_prescribed
    REAL(dp)                                              :: timeframe_BMB_prescribed

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE ('NAM')
        filename_BMB_prescribed  = C%filename_BMB_prescribed_NAM
        timeframe_BMB_prescribed = C%timeframe_BMB_prescribed_NAM
      CASE ('EAS')
        filename_BMB_prescribed  = C%filename_BMB_prescribed_EAS
        timeframe_BMB_prescribed = C%timeframe_BMB_prescribed_EAS
      CASE ('GRL')
        filename_BMB_prescribed  = C%filename_BMB_prescribed_GRL
        timeframe_BMB_prescribed = C%timeframe_BMB_prescribed_GRL
      CASE ('ANT')
        filename_BMB_prescribed  = C%filename_BMB_prescribed_ANT
        timeframe_BMB_prescribed = C%timeframe_BMB_prescribed_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising BMB from file "' // colour_string( TRIM( filename_BMB_prescribed),'light blue') // '"...'

    ! Read BMB from file
    IF (timeframe_BMB_prescribed == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_2D( filename_BMB_prescribed, 'BMB||basal_mass_balance||', mesh, BMB%BMB_shelf)
    ELSE
      ! Assume the file has a time dimension, and read the specified timeframe
      CALL read_field_from_file_2D( filename_BMB_prescribed, 'BMB||basal_mass_balance||', mesh, BMB%BMB_shelf, time_to_read = timeframe_BMB_prescribed)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_prescribed_notime

END MODULE BMB_prescribed
