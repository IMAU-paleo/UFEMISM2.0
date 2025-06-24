module Halfar_SIA_solution

  use precisions, only: dp
  use parameters

  implicit none

  private

  public :: Halfar_dome

contains

  subroutine Halfar_dome( A, n, H0, R0, x, y, t, H)
    ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity
    ! function with dome thickness H0 and margin radius R0 at t0.
    !
    ! Halfar, P.: On the Dynamics of Ice Sheets, Journal of Geophysical Research 86,
    !   11065-11072, 1981

    ! In/output variables:
    real(dp), intent(in   ) :: A          ! [Pa^-3 yr^-1] Glen's flow law parameter
    real(dp), intent(in   ) :: n          ! [ ]           Glen's flow law exponent
    real(dp), intent(in   ) :: H0         ! [m]           Thickness at ice divide at t=0
    real(dp), intent(in   ) :: R0         ! [m]           Ice margin radius at t=0
    real(dp), intent(in   ) :: x          ! [m]           x-coordinate
    real(dp), intent(in   ) :: y          ! [m]           y-coordinate
    real(dp), intent(in   ) :: t          ! [yr]          time
    real(dp), intent(  out) :: H          ! [m]           Ice thickness at [x,y,t]

    ! Local variables
    real(dp) :: Gamma, t0, r, f1, f2, f3, tp

    Gamma = (2._dp / 5._dp) * (A / sec_per_year) * (ice_density * grav)**n
    t0 = 1._dp / ((5._dp * n + 3._dp) * Gamma) * ((2._dp * n + 1._dp)/(n + 1._dp))**n * (R0**(n + 1._dp))/(H0**(2._dp * n  + 1))

    tp = (t * sec_per_year) + t0

    r = SQRT(x**2._dp + y**2._dp)

    f1 = (t0/tp)**(2._dp/(5._dp * n + 3._dp))
    f2 = (t0/tp)**(1._dp/(5._dp * n + 3._dp))
    f3 = (r/R0)

    H = H0 * f1 * MAX( 0._dp, (1._dp - (f2*f3)**((n + 1._dp)/n)))**(n/(2._dp * n + 1._dp))

  end subroutine Halfar_dome

end module Halfar_SIA_solution