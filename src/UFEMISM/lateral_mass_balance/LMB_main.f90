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
  use netcdf_io_main
  use LMB_GlacialIndex

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
          IF (ice%mask_cf_fl( vi) .OR. ice%mask_cf_gr( vi)) THEN
            LMB%LMB( vi) = C%uniform_LMB
          END IF
        END DO
      CASE ('GlacialIndex')
        CALL run_LMB_model_GlacialIndex(mesh, ice, LMB, time)
      CASE ('inverted')
        ! Nothing to do here
      CASE DEFAULT
        CALL crash('unknown choice_LMB_model "' // TRIM( choice_LMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_LMB_model

  SUBROUTINE initialise_LMB_model( mesh, LMB, region_name, start_time_of_run)
    ! Initialise the LMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_LMB_model),                   INTENT(OUT)   :: LMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: start_time_of_run

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_LMB_model'
    CHARACTER(LEN=256)                                    :: choice_LMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising lateral mass balance model...'

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
      CASE ('GlacialIndex')
        if (par%primary)  write(*,"(A)") '     Initialising LMB model "' // &
            colour_string( trim( choice_LMB_model),'light blue') // '"...'
        CALL initialise_LMB_model_GlacialIndex(mesh, LMB, region_name, start_time_of_run)
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
    IF (par%primary)  WRITE(*,"(A)") '    Remapping lateral mass balance model data to the new mesh...'

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
        ! No need to do anything
      CASE ('GlacialIndex')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_LMB_model "' // TRIM( choice_LMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_LMB_model

END MODULE LMB_main
