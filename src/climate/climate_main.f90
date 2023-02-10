MODULE climate_main
  ! The main climate module

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE parameters

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_module( i)

    IMPLICIT NONE

    ! In/output variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER                                            , INTENT(IN   ) :: i

    WRITE(0,*) '   run_climate_module: i = ', i

  END SUBROUTINE run_climate_module


END MODULE climate_main