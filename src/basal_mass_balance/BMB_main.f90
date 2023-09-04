MODULE BMB_main

  ! The main BMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE BMB_idealised                                          , ONLY: initialise_BMB_model_idealised, run_BMB_model_idealised
  USE BMB_parameterised                                      , ONLY: initialise_BMB_model_parameterised, run_BMB_model_parameterised
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE math_utilities                                         , ONLY: is_floating

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model( mesh, ice, ocean, refgeo, BMB, region_name, time)
    ! Calculate the basal mass balance

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new BMB
    IF (C%do_asynchronous_BMB) THEN
      ! Asynchronous coupling: do not calculate a new BMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next BMB time step
      IF (time == BMB%t_next) THEN
        ! Go on to calculate a new BMB
        BMB%t_next = time + C%dt_BMB
      ELSEIF (time > BMB%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the BMB time step')
      ELSE
        ! It is not yet time to calculate a new BMB
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_BMB) THEN
      ! Synchronous coupling: calculate a new BMB in every model loop
      BMB%t_next = time + C%dt_BMB
    END IF

    ! Determine which BMB model to run for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Run the chosen BMB model
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        BMB%BMB = C%uniform_BMB
      CASE ('idealised')
        CALL run_BMB_model_idealised( mesh, ice, BMB, time)
      CASE ('parameterised')
        CALL run_BMB_model_parameterised( mesh, ice, ocean, BMB)
      CASE ('inverted')
        CALL run_BMB_model_inverted( mesh, ice, refgeo, BMB, time)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model

  SUBROUTINE initialise_BMB_model( mesh, BMB, region_name)
    ! Initialise the BMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(OUT)   :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising basal mass balance model...'

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Allocate memory for main variables
    ALLOCATE( BMB%BMB( mesh%vi1:mesh%vi2))
    BMB%BMB = 0._dp

    ! Set time of next calculation to start time
    BMB%t_next = C%start_time_of_run

    ! Determine which BMB model to initialise
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        CALL initialise_BMB_model_idealised( mesh, BMB)
      CASE ('parameterised')
        CALL initialise_BMB_model_parameterised( mesh, BMB)
      CASE ('inverted')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model

  SUBROUTINE write_to_restart_file_BMB_model( mesh, BMB, region_name, time)
    ! Write to the restart file for the BMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(IN)    :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Write to the restart file of the chosen BMB model
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        ! No need to do anything
      CASE ('inverted')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_BMB_model

  SUBROUTINE create_restart_file_BMB_model( mesh, BMB, region_name)
    ! Create the restart file for the BMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Create the restart file of the chosen BMB model
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        ! No need to do anything
      CASE ('inverted')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_BMB_model

  SUBROUTINE remap_BMB_model( mesh_old, mesh_new, BMB, region_name)
    ! Remap the BMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_BMB_model),                   INTENT(OUT)   :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '    Remapping basal mass balance model data to the new mesh...'

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Reallocate memory for main variables
    CALL reallocate_bounds( BMB%BMB, mesh_new%vi1, mesh_new%vi2)

    ! Determine which BMB model to initialise
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        CALL crash('Remapping after mesh update not implemented yet for parameterised BMB')
      CASE ('inverted')
        CALL crash('Remapping after mesh update not implemented yet for parameterised BMB')
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BMB_model

  SUBROUTINE run_BMB_model_inverted( mesh, ice, refgeo, BMB, time)
    ! Calculate the basal mass balance
    !
    ! Use a parameterised BMB scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_inverted'
    INTEGER                                            :: vi
    REAL(dp)                                           :: misfit

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      IF (ice%fraction_gr( vi) == 1._dp .AND. .NOT. is_floating( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)) THEN
        CYCLE
      END IF

      ! Ice thickness misfit
      misfit = ice%Hi( vi) - refgeo%Hi( vi)

      ! Is it improving already?
      IF (ice%dHi_dt( vi)*misfit < 0._dp) THEN
        ! Yes, so leave this vertex alone
        CYCLE
      END IF

      BMB%BMB( vi) =  BMB%BMB( vi) - 1.0_dp * TAN( MIN( 1.0_dp, MAX( -1.0_dp, misfit/5E2_dp)))
      BMB%BMB( vi) = MIN( BMB%BMB( vi),   20._dp)
      BMB%BMB( vi) = MAX( BMB%BMB( vi), -200._dp)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_inverted

END MODULE BMB_main
