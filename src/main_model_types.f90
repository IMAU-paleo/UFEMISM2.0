MODULE main_model_types

  ! The derived data types for the main UFEMISM model

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model

  IMPLICIT NONE

  ! == Reference geometries
  ! =======================

CONTAINS

! ===== Main routines =====
! =========================

END MODULE main_model_types