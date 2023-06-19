MODULE BMB_main

  ! The main BMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE BMB_types                                              , ONLY: type_BMB_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE initialise_BMB_model( mesh, BMB)
    ! Initialise all data fields of the BMB module

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    TYPE(type_BMB_model),          INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_BMB_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      WRITE(*,"(A)") '  Initialising basal mass balance model...'
    END IF
    CALL sync

    ! Allocate memory for main variables
    ALLOCATE( BMB%BMB( mesh%vi1:mesh%vi1))
    BMB%BMB = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model

END MODULE BMB_main
