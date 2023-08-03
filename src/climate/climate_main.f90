MODULE climate_main

  ! The main climate model module.

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
  USE climate_idealised                                      , ONLY: initialise_climate_model_idealised, run_climate_model_idealised
  USE reallocate_mod                                         , ONLY: reallocate_bounds

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model( mesh, ice, climate, region_name, time)
    ! Calculate the climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new climate
    IF (C%do_asynchronous_climate) THEN
      ! Asynchronous coupling: do not calculate a new climate in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next climate time step
      IF (time == climate%t_next) THEN
        ! Go on to calculate a new climate
        climate%t_next = time + C%dt_climate
      ELSEIF (time > climate%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the climate time step')
      ELSE
        ! It is not yet time to calculate a new climate
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_climate) THEN
      ! Synchronous coupling: calculate a new climate in every model loop
      climate%t_next = time + C%dt_climate
    END IF

    ! Determine which climate model to run for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Run the chosen climate model
    IF     (choice_climate_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_climate_model == 'idealised') THEN
      CALL run_climate_model_idealised( mesh, ice, climate, time)
    ELSE
      CALL crash('unknown choice_climate_model "' // TRIM( choice_climate_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model

  SUBROUTINE initialise_climate_model( mesh, climate, region_name)
    ! Initialise the climate model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(OUT)   :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '  Initialising climate model...'

    ! Determine which climate model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Allocate memory for main variables
    ALLOCATE( climate%T2m(    mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%Precip( mesh%vi1:mesh%vi2,12))
    climate%T2m    = 0._dp
    climate%Precip = 0._dp

    ! Set time of next calculation to start time
    climate%t_next = C%start_time_of_run

    ! Determine which climate model to initialise
    IF     (choice_climate_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_climate_model == 'idealised') THEN
      CALL initialise_climate_model_idealised( mesh, climate)
    ELSE
      CALL crash('unknown choice_climate_model "' // TRIM( choice_climate_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model

  SUBROUTINE write_to_restart_file_climate_model( mesh, climate, region_name, time)
    ! Write to the restart file for the climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(IN)    :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which climate model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Write to the restart file of the chosen climate model
    IF     (choice_climate_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_climate_model == 'idealised') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_climate_model "' // TRIM( choice_climate_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_climate_model

  SUBROUTINE create_restart_file_climate_model( mesh, climate, region_name)
    ! Create the restart file for the climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which climate model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Create the restart file of the chosen climate model
    IF     (choice_climate_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_climate_model == 'idealised') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_climate_model "' // TRIM( choice_climate_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_climate_model

  SUBROUTINE remap_climate_model( mesh_old, mesh_new, climate, region_name)
    ! Remap the climate model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_climate_model'
    CHARACTER(LEN=256)                                    :: choice_climate_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '    Remapping climate model data to the new mesh...'

    ! Determine which climate model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_climate_model = C%choice_climate_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_climate_model = C%choice_climate_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_climate_model = C%choice_climate_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_climate_model = C%choice_climate_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Reallocate memory for main variables
    CALL reallocate_bounds( climate%T2m,    mesh_new%vi1, mesh_new%vi2,12)
    CALL reallocate_bounds( climate%Precip, mesh_new%vi1, mesh_new%vi2,12)

    ! Determine which climate model to remap
    IF     (choice_climate_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_climate_model == 'idealised') THEN
      ! No need to remap anything here
    ELSE
      CALL crash('unknown choice_climate_model "' // TRIM( choice_climate_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_climate_model

END MODULE climate_main
