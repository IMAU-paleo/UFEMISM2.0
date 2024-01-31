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

  SUBROUTINE run_BMB_model_laddie( mesh, BMB)
    ! Calculate the basal mass balance
    !
    ! Call the external LADDIE model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_laddie'
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_runtime
    ! CHARACTER(LEN=256)                                    :: filename_BMB_laddie_config
    ! REAL(dp), DIMENSION(:,:), POINTER                     :: BMB_LADDIE

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Here, run LADDIE

    ! Read BMB from file
    CALL read_field_from_file_2D( C%filename_BMB_laddie_runtime, 'BMB', mesh, BMB%BMB)

    ! Convert to m.i.e./yr
    BMB%BMB = 31557600._dp * BMB%BMB / 918._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE



  SUBROUTINE initialise_BMB_model_laddie( mesh, BMB)
    ! Initialise the BMB model
    !
    ! Call the external LADDIE model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_laddie'
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_initial
    ! REAL(dp), DIMENSION(:,:), POINTER                     :: BMB_LADDIE

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Here, run LADDIE

    ! Read BMB from file
    CALL read_field_from_file_2D( C%filename_BMB_laddie_initial, 'BMB', mesh, BMB%BMB)

    ! Convert to m.i.e./yr
    BMB%BMB = 31557600._dp * BMB%BMB / 918._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE


END MODULE BMB_laddie