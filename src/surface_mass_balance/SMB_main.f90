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
  USE SMB_prescribed                                         , ONLY: initialise_SMB_model_prescribed, run_SMB_model_prescribed
  USE reallocate_mod                                         , ONLY: reallocate_bounds

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
      SMB%t_next = time + C%dt_SMB
    END IF

    ! Determine which SMB model to run for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Run the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        SMB%SMB = C%uniform_SMB
      CASE ('idealised')
        CALL run_SMB_model_idealised( mesh, ice, SMB, time)
      CASE ('prescribed')
        CALL run_SMB_model_prescribed( mesh, ice, SMB, region_name, time)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

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
    IF (par%master)  WRITE(*,"(A)") '   Initialising surface mass balance model...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Allocate memory for main variables
    ALLOCATE( SMB%SMB( mesh%vi1:mesh%vi2))
    SMB%SMB = 0._dp

    IF (C%do_corrections_SMB) THEN
      ! Allocate memory for main variables
      ALLOCATE( SMB%SMB_correction( mesh%vi1:mesh%vi2))
      SMB%SMB_correction = 0._dp
    END IF

    ! Set time of next calculation to start time
    SMB%t_next = C%start_time_of_run

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        SMB%SMB = C%uniform_SMB
      CASE ('idealised')
        CALL initialise_SMB_model_idealised( mesh, SMB)
      CASE ('prescribed')
        CALL initialise_SMB_model_prescribed( mesh, SMB, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model

  SUBROUTINE write_to_restart_file_SMB_model( mesh, SMB, region_name, time)
    ! Write to the restart file for the SMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which SMB model to use for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Write to the restart file of the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_SMB_model

  SUBROUTINE create_restart_file_SMB_model( mesh, SMB, region_name)
    ! Create the restart file for the SMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which SMB model to use for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Create the restart file of the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_SMB_model

  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, SMB, region_name)
    ! Remap the SMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_SMB_model),                   INTENT(OUT)   :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '    Remapping surface mass balance model data to the new mesh...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Reallocate memory for main variables
    CALL reallocate_bounds( SMB%SMB, mesh_new%vi1, mesh_new%vi2)

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        CALL crash('Remapping after mesh update not implemented yet for prescribed SMB')
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SMB_model

! ===== Utilities =====
! =====================

  SUBROUTINE SMB_adjustment( mesh, ice, SMB, dHi_dt_predicted, Hi_predicted, time)
    ! Calculate the basal mass balance
    !
    ! Use an inversion based on the computed dHi_dt

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: dHi_dt_predicted
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi_predicted
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'SMB_adjustment'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Correct land SMB based on remaining dHi_dt at the end of an inversion for a better equilibrium.
    IF (.NOT. C%do_corrections_SMB .OR. &
        time < C%SMB_residual_absorb_t_start .OR. &
        time > C%SMB_residual_absorb_t_end) THEN

      ! Finalise routine path
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    DO vi = mesh%vi1, mesh%vi2
      ! If interior grounded point or ice-free land
      IF (.NOT. ice%mask_icefree_ocean( vi) .AND. &
          .NOT. ice%mask_floating_ice( vi) .AND. &
          .NOT. ice%mask_gl_gr( vi) .AND. &
          .NOT. ice%mask_cf_gr( vi)) THEN

        ! For grounded ice, use dHi_dt to get an "inversion" of equilibrium SMB.
        SMB%SMB( vi) = SMB%SMB( vi) - dHi_dt_predicted( vi)

        ! Adjust rate of ice thickness change dHi/dt to compensate the change
        dHi_dt_predicted( vi) = 0._dp

        ! Adjust corrected ice thickness to compensate the change
        Hi_predicted( vi) = ice%Hi_prev( vi)

      END IF
    END DO ! vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE SMB_adjustment

END MODULE SMB_main
