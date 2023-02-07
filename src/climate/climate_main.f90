MODULE climate_main
  ! The main climate module

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>

  USE precisions                                             , ONLY: dp

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