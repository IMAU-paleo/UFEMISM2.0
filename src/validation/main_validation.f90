MODULE main_validation

  ! The main validation module

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE unit_tests_mpi                                         , ONLY: run_all_mpi_distributed_memory_unit_tests
  USE unit_tests_petsc                                       , ONLY: run_all_petsc_unit_tests
  USE unit_tests_mesh                                        , ONLY: run_all_mesh_unit_tests
  USE unit_tests_netcdf                                      , ONLY: run_all_netcdf_unit_tests
  USE unit_tests_ice                                         , ONLY: run_all_ice_unit_tests

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_unit_tests
    ! Run all unit tests

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_all_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all unit tests
    CALL run_all_mpi_distributed_memory_unit_tests
    CALL run_all_petsc_unit_tests
    CALL run_all_mesh_unit_tests
    CALL run_all_netcdf_unit_tests

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_unit_tests

END MODULE main_validation
