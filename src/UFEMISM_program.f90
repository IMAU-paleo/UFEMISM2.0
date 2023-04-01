PROGRAM UFEMISM_program

  ! Hieperdepiep!

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync, initialise_parallelisation, finalise_parallelisation
  USE petsc_basic                                            , ONLY: perr
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, initialise_control_and_resource_tracker
  USE UFEMISM_main_model                                     , ONLY: run_model

  USE unit_tests_mpi
  USE unit_tests_petsc
  USE unit_tests_mesh
  USE unit_tests_netcdf

  IMPLICIT NONE

! ===== Main variables =====
! ==========================

! ===== START =====
! =================

  ! Initialise MPI parallelisation
  CALL initialise_parallelisation
  CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Initialise the control and resource tracker
  CALL initialise_control_and_resource_tracker




  ! Test
!  CALL run_all_mpi_distributed_memory_unit_tests
!  CALL run_all_petsc_unit_tests
!  CALL run_all_mesh_unit_tests
  CALL run_all_netcdf_unit_tests




! ===== FINISH =====
! ==================

  ! Finalise MPI parallelisation
  CALL PetscFinalize( perr)
  CALL finalise_parallelisation

END PROGRAM UFEMISM_program