PROGRAM UFEMISM_program
  ! Hieperdepiep!

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>

  USE precisions                                             , ONLY: dp
  USE UFEMISM_main_model                                     , ONLY: run_model

# if (defined(DO_SELEN))
  USE SELEN_main_module                                      , ONLY: initialise_SELEN, run_SELEN
# endif

! ===== Main variables =====
! ==========================

  IMPLICIT NONE

! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
  INTEGER                                                            :: i

! ===== START =====
! =================

  WRITE(0,*) 'Hello world!'

  i = 1

  CALL run_model( i)

END PROGRAM UFEMISM_program