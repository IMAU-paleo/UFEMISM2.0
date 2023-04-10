MODULE unit_tests_ice

  ! Unit tests for different ice velocity solvers.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE ice_model_types                                        , ONLY: type_ice_model

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_ice_unit_tests
    ! Run all unit tests for the different ice velocity solvers

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_all_ice_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all unit tests for the different ice velocity solvers

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_ice_unit_tests

END MODULE unit_tests_ice
