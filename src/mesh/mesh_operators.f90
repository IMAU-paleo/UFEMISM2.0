MODULE mesh_operators

  ! Routines for calculating matrix operators on the mesh.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

END MODULE mesh_operators
