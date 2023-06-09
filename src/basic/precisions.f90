MODULE precisions

  ! The different precision KINDs used by UFEMISM

! ===== Preamble =====
! ====================

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  INTEGER, PARAMETER  :: dp  = KIND(1.0D0)  ! KIND of double precision numbers. Reals should be declared as: REAL(dp) :: example

END MODULE precisions