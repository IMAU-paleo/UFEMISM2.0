MODULE UFEMISM_main_model

  ! The main regional ice-sheet model

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE parameters
  USE climate_main                                           , ONLY: run_climate_module

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_model( i)

    IMPLICIT NONE

    ! In/output variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER                                            , INTENT(IN   ) :: i

    WRITE(0,*) '   i = ', i

    CALL run_climate_module( i+1)

  END SUBROUTINE run_model


END MODULE UFEMISM_main_model