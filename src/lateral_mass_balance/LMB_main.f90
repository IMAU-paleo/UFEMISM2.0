MODULE LMB_main

  ! The main LMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE LMB_model_types                                        , ONLY: type_LMB_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file, open_existing_netcdf_file_for_writing
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, add_time_dimension_to_file, &
                                                                     add_field_mesh_dp_2D, write_to_field_multopt_mesh_dp_2D, write_time_to_file, write_to_field_multopt_mesh_dp_3D_ocean

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_LMB_model( mesh, ice, LMB, region_name, time)
    ! Calculate the lateral mass balance

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_LMB_model),                   INTENT(INOUT) :: LMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_LMB_model'
    CHARACTER(LEN=256)                                    :: choice_LMB_model
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Synchronous coupling: calculate a new LMB in every model loop
    LMB%t_next = time + C%dt_LMB

    ! Determine which LMB model to run for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_LMB_model = C%choice_LMB_model_NAM
      CASE ('EAS')
        choice_LMB_model = C%choice_LMB_model_EAS
      CASE ('GRL')
        choice_LMB_model = C%choice_LMB_model_GRL
      CASE ('ANT')
        choice_LMB_model = C%choice_LMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Run the chosen LMB model
    SELECT CASE (choice_LMB_model)
      CASE ('uniform')
        LMB%LMB = 0._dp
        DO vi = mesh%vi1, mesh%vi2
          IF (ice%mask_cf_fl( vi)) LMB%LMB( vi) = C%uniform_LMB
        END DO
      CASE ('inverted')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_LMB_model "' // TRIM( choice_LMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_LMB_model

  SUBROUTINE initialise_LMB_model( mesh, LMB, region_name)
    ! Initialise the LMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_LMB_model),                   INTENT(OUT)   :: LMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_LMB_model'
    CHARACTER(LEN=256)                                    :: choice_LMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising lateral mass balance model...'

    ! Determine which LMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_LMB_model = C%choice_LMB_model_NAM
      CASE ('EAS')
        choice_LMB_model = C%choice_LMB_model_EAS
      CASE ('GRL')
        choice_LMB_model = C%choice_LMB_model_GRL
      CASE ('ANT')
        choice_LMB_model = C%choice_LMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Allocate memory for main variables
    ALLOCATE( LMB%LMB( mesh%vi1:mesh%vi2))
    LMB%LMB = 0._dp

    ! Set time of next calculation to start time
    LMB%t_next = C%start_time_of_run

    ! Determine which LMB model to initialise
    SELECT CASE (choice_LMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('inverted')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_LMB_model "' // TRIM( choice_LMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_LMB_model

  SUBROUTINE remap_LMB_model( mesh_old, mesh_new, LMB, region_name)
    ! Remap the LMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_LMB_model),                   INTENT(OUT)   :: LMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_LMB_model'
    CHARACTER(LEN=256)                                    :: choice_LMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '    Remapping lateral mass balance model data to the new mesh...'

    ! Determine which LMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_LMB_model = C%choice_LMB_model_NAM
      CASE ('EAS')
        choice_LMB_model = C%choice_LMB_model_EAS
      CASE ('GRL')
        choice_LMB_model = C%choice_LMB_model_GRL
      CASE ('ANT')
        choice_LMB_model = C%choice_LMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Reallocate memory for main variables
    CALL reallocate_bounds( LMB%LMB, mesh_new%vi1, mesh_new%vi2)

    ! Determine which LMB model to initialise
    SELECT CASE (choice_LMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('inverted')
        CALL crash('Remapping after mesh update not implemented yet for parameterised LMB')
      CASE DEFAULT
        CALL crash('unknown choice_LMB_model "' // TRIM( choice_LMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_LMB_model

! ===== Utilities =====
! =====================

  ! SUBROUTINE LMB_inversion( mesh, ice, SMB, LMB, dHi_dt_predicted, Hi_predicted, dt, time, region_name)
  !   ! Calculate the basal mass balance
  !   !
  !   ! Use an inversion based on the computed dHi_dt

  !   IMPLICIT NONE

  !   ! In/output variables:
  !   TYPE(type_mesh),                        INTENT(IN)    :: mesh
  !   TYPE(type_ice_model),                   INTENT(IN)    :: ice
  !   TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
  !   TYPE(type_LMB_model),                   INTENT(INOUT) :: LMB
  !   REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: dHi_dt_predicted
  !   REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi_predicted
  !   REAL(dp),                               INTENT(IN)    :: dt
  !   REAL(dp),                               INTENT(IN)    :: time
  !   CHARACTER(LEN=3)                                      :: region_name

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'LMB_inversion'
  !   INTEGER                                               :: vi
  !   CHARACTER(LEN=256)                                    :: choice_LMB_model
  !   REAL(dp)                                              :: LMB_previous, LMB_change

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Determine filename for this model region
  !   SELECT CASE (region_name)
  !     CASE ('NAM')
  !       choice_LMB_model  = C%choice_LMB_model_NAM
  !     CASE ('EAS')
  !       choice_LMB_model  = C%choice_LMB_model_EAS
  !     CASE ('GRL')
  !       choice_LMB_model  = C%choice_LMB_model_GRL
  !     CASE ('ANT')
  !       choice_LMB_model  = C%choice_LMB_model_ANT
  !     CASE DEFAULT
  !       CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
  !   END SELECT

  !   ! Invert ocean LMB based on the full dHi_dt at each time step
  !   IF (.NOT. choice_LMB_model == 'inverted' .OR. &
  !       time < C%LMB_inversion_t_start .OR. &
  !       time > C%LMB_inversion_t_end) THEN
  !     ! Finalise routine path
  !     CALL finalise_routine( routine_name)
  !     RETURN
  !   END IF

  !   DO vi = mesh%vi1, mesh%vi2

  !     ! For these areas, use dHi_dt to get an "inversion" of equilibrium LMB.
  !     IF (ice%mask_cf_fl( vi) .OR. ice%mask_cf_gr( vi)) THEN

  !       ! Save previous value to compute the actual adjustment later
  !       LMB_previous = LMB%LMB( vi)

  !       ! Assume that calving accounts for half the mass loss here (other half is melting)
  !       LMB%LMB( vi) = MIN( 0._dp, LMB%LMB( vi) - dHi_dt_predicted( vi)/2._dp)

  !       ! Compute actual change in LMB
  !       LMB_change = LMB%LMB( vi) - LMB_previous

  !       ! Adjust rate of ice thickness change dHi/dt to compensate the change
  !       dHi_dt_predicted( vi) = dHi_dt_predicted( vi) + LMB_change

  !       ! Adjust new ice thickness to compensate the change
  !       Hi_predicted( vi) = ice%Hi_prev( vi) + dHi_dt_predicted( vi) * dt

  !     ELSEIF (ice%mask_cf_fl( vi) .OR. ice%mask_cf_gr( vi)) THEN

  !       ! Save previous value to compute the actual adjustment later
  !       LMB_previous = LMB%LMB( vi)

  !       ! Assume that calving accounts for half the mass loss here (other half is melting)
  !       LMB%LMB( vi) = MIN( 0._dp, LMB%LMB( vi) - dHi_dt_predicted( vi)/2._dp)

  !       ! Compute actual change in LMB
  !       LMB_change = LMB%LMB( vi) - LMB_previous

  !       ! Adjust rate of ice thickness change dHi/dt to compensate the change
  !       dHi_dt_predicted( vi) = dHi_dt_predicted( vi) + LMB_change

  !       ! Adjust new ice thickness to compensate the change
  !       Hi_predicted( vi) = ice%Hi_prev( vi) + dHi_dt_predicted( vi) * dt

  !     ELSE
  !       ! Not a place where lateral loss operates
  !       LMB%LMB( vi) = 0._dp
  !     END IF

  !   END DO ! vi = mesh%vi1, mesh%vi2

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE LMB_inversion

END MODULE LMB_main
