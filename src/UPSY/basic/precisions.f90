module precisions

  ! The different precision KINDs used by UFEMISM

  implicit none

  integer, parameter :: dp  = KIND(1.0D0)  ! KIND of double precision numbers. Reals should be declared as: REAL(dp) :: example
  integer, parameter :: int8 = kind(1_8)   ! KIND of 8-byte integer numbers

end module precisions