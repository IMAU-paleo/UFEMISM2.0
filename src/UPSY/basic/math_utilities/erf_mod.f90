module erf_mod

  ! The error function

  use precisions, only: dp
  use parameters, only: pi

  implicit none

contains

  pure function error_function( X) result( ERR)
    ! Tijn Berends, 2024: not sure where we got this code from,
    ! it was already present in ANICE, doesn't look like something
    ! written at IMAU but there was no reference. Ah, well...

    ! Purpose: Compute error function erf(x)
    ! Input:   x   --- Argument of erf(x)
    ! Output:  ERR --- erf(x)

    ! Input variables:
    real(dp), intent(in)  :: X

    ! Output variables:
    real(dp) :: ERR

    ! Local variables:
    real(dp) :: EPS
    real(dp) :: X2
    real(dp) :: ER
    real(dp) :: R
    real(dp) :: C0
    integer  :: k

    EPS = 1.0E-15_dp
    X2  = X * X
    if (abs(X) < 3.5_dp) then
      ER = 1.0_dp
      R  = 1.0_dp
      do k = 1, 50
        R  = R * X2 / (real(k, dp) + 0.5_dp)
        ER = ER+R
        if (abs(R) < abs(ER) * EPS) then
          C0  = 2.0_dp / sqrt(pi) * X * exp(-X2)
          ERR = C0 * ER
          exit
        end if
      end do
    else
      ER = 1.0_dp
      R  = 1.0_dp
      do k = 1, 12
        R  = -R * (real(k, dp) - 0.5_dp) / X2
        ER = ER + R
        C0  = exp(-X2) / (abs(X) * sqrt(pi))
        ERR = 1.0_dp - C0 * ER
        if (X < 0.0_dp) ERR = -ERR
      end do
    end if

  end function error_function

end module erf_mod
