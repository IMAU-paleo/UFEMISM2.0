MODULE BMB_laddie

  ! LADDIE model

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE netcdf_input                                           , ONLY: read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model_laddie( mesh, ice, BMB, time)
    ! Calculate the basal mass balance
    !
    ! Call the external LADDIE model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_laddie'
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_runtime
    REAL(dp), DIMENSION(:,:), POINTER                     :: BMB_LADDIE

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE ('NAM')
        filename_BMB_laddie_runtime  = C%filename_BMB_laddie_runtime_NAM
      CASE ('EAS')
        filename_BMB_laddie_runtime  = C%filename_BMB_laddie_runtime_EAS
      CASE ('GRL')
        filename_BMB_laddie_runtime  = C%filename_BMB_laddie_runtime_GRL
      CASE ('ANT')
        filename_BMB_laddie_runtime  = C%filename_BMB_laddie_runtime_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Here, run LADDIE

    ! Read BMB from file
    CALL read_field_from_file_2D( filename_BMB_laddie_runtime, 'BMB', mesh, BMB%BMB)

  END SUBROUTINE



  SUBROUTINE initialise_BMB_model_laddie( mesh, ice, BMB, time)
    ! Initialise the BMB model
    !
    ! Call the external LADDIE model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_laddie'
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_initial
    REAL(dp), DIMENSION(:,:), POINTER                     :: BMB_LADDIE

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE ('NAM')
        filename_BMB_laddie_initial  = C%filename_BMB_laddie_initial_NAM
      CASE ('EAS')
        filename_BMB_laddie_initial  = C%filename_BMB_laddie_initial_EAS
      CASE ('GRL')
        filename_BMB_laddie_initial  = C%filename_BMB_laddie_initial_GRL
      CASE ('ANT')
        filename_BMB_laddie_initial  = C%filename_BMB_laddie_initial_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Here, run LADDIE

    ! Read BMB from file
    CALL read_field_from_file_2D( filename_BMB_laddie_initial, 'BMB', mesh, BMB%BMB)

  END SUBROUTINE


END MODULE BMB_idealised