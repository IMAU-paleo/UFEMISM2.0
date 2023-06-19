MODULE SMB_main

  ! The main SMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE SMB_types                                              , ONLY: type_SMB_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE initialise_SMB_model( mesh, SMB)
    ! Initialise all data fields of the SMB module

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    TYPE(type_SMB_model),          INTENT(INOUT) :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_SMB_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      WRITE(*,"(A)") '  Initialising surface mass balance model...'
    END IF
    CALL sync

    ! Allocate memory for main variables
    ALLOCATE( SMB%SMB( mesh%vi1:mesh%vi1))
    SMB%SMB = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model

END MODULE SMB_main
