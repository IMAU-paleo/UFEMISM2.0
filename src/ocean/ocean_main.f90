MODULE ocean_main

  ! The main ocean model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ocean_model( mesh, ice, ocean, region_name, time)
    ! Calculate the ocean

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new ocean
    IF (C%do_asynchronous_ocean) THEN
      ! Asynchronous coupling: do not calculate a new ocean in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next ocean time step
      IF (time == ocean%t_next) THEN
        ! Go on to calculate a new ocean
        ocean%t_next = time + C%dt_ocean
      ELSEIF (time > ocean%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the ocean time step')
      ELSE
        ! It is not yet time to calculate a new ocean
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_ocean) THEN
      ! Synchronous coupling: calculate a new ocean in every model loop
    END IF

    ! Determine which ocean model to run for this region
    IF     (region_name == 'NAM') THEN
      choice_ocean_model = C%choice_ocean_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_ocean_model = C%choice_ocean_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_ocean_model = C%choice_ocean_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_ocean_model = C%choice_ocean_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Run the chosen ocean model
    IF (choice_ocean_model == 'none') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model

  SUBROUTINE initialise_ocean_model( mesh, ocean, region_name)
    ! Initialise the ocean model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),               INTENT(OUT)   :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '  Initialising ocean model...'

    ! Determine which ocean model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_ocean_model = C%choice_ocean_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_ocean_model = C%choice_ocean_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_ocean_model = C%choice_ocean_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_ocean_model = C%choice_ocean_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Allocate memory for main variables
    ALLOCATE( ocean%T( mesh%vi1:mesh%vi1,12))
    ALLOCATE( ocean%S( mesh%vi1:mesh%vi1,12))
    ocean%T = 0._dp
    ocean%S = 0._dp

    ! Set time of next calculation to start time
    ocean%t_next = C%start_time_of_run

    ! Determine which ocean model to initialise
    IF (choice_ocean_model == 'none') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model

END MODULE ocean_main
