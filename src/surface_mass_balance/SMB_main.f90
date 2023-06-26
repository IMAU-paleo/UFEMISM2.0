MODULE SMB_main

  ! The main SMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE SMB_idealised                                          , ONLY: initialise_SMB_model_idealised, run_SMB_model_idealised

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_SMB_model( mesh, ice, climate, SMB, region_name, time)
    ! Calculate the surface mass balance

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(IN)    :: climate
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new SMB
    IF (C%do_asynchronous_SMB) THEN
      ! Asynchronous coupling: do not calculate a new SMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next SMB time step
      IF (time == SMB%t_next) THEN
        ! Go on to calculate a new SMB
        SMB%t_next = time + C%dt_SMB
      ELSEIF (time > SMB%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the SMB time step')
      ELSE
        ! It is not yet time to calculate a new SMB
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_SMB) THEN
      ! Synchronous coupling: calculate a new SMB in every model loop
    END IF

    ! Determine which SMB model to run for this region
    IF     (region_name == 'NAM') THEN
      choice_SMB_model = C%choice_SMB_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_SMB_model = C%choice_SMB_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_SMB_model = C%choice_SMB_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_SMB_model = C%choice_SMB_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Run the chosen SMB model
    IF     (choice_SMB_model == 'uniform') THEN
      SMB%SMB = C%uniform_SMB
    ELSEIF (choice_SMB_model == 'idealised') THEN
      CALL run_SMB_model_idealised( mesh, ice, SMB, time)
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model

  SUBROUTINE initialise_SMB_model( mesh, SMB, region_name)
    ! Initialise the SMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(OUT)   :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '  Initialising surface mass balance model...'

    ! Determine which SMB model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_SMB_model = C%choice_SMB_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_SMB_model = C%choice_SMB_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_SMB_model = C%choice_SMB_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_SMB_model = C%choice_SMB_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Allocate memory for main variables
    ALLOCATE( SMB%SMB( mesh%vi1:mesh%vi2))
    SMB%SMB = 0._dp

    ! Set time of next calculation to start time
    SMB%t_next = C%start_time_of_run

    ! Determine which SMB model to initialise
    IF     (choice_SMB_model == 'uniform') THEN
      ! No need to do anything
    ELSEIF (choice_SMB_model == 'idealised') THEN
      CALL initialise_SMB_model_idealised( mesh, SMB)
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model

END MODULE SMB_main
