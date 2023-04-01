MODULE ice_model_main

  ! The main ice-dynamical model module.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE ice_model_types                                        , ONLY: type_ice_model

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

END MODULE ice_model_main
