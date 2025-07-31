module Bueler_SIA_solution

  use precisions, only: dp
  use parameters

  implicit none

  private

  public :: Bueler_dome

contains

  subroutine Bueler_dome( A, n, H0, R0, lambda, x, y, t, H, M)
    ! Extends the Halfar dome solution to include a positive accumulation rate
    !
    ! Bueler, E., Lingle, C. S., Kallen-Brown, J. A., Covey, D. N., and Bowman, L. N.:
    !   Exact solutions and verification of numerical models for isothermal ice sheets,
    !   Journal of Glaciology 51, 291-306, 2005

    ! In/output variables:
    real(dp), intent(in   ) :: A          ! [Pa^-3 yr^-1] Glen's flow law parameter
    real(dp), intent(in   ) :: n          ! [ ]           Glen's flow law exponent
    real(dp), intent(in   ) :: H0         ! [m]           Thickness at ice divide at t=0
    real(dp), intent(in   ) :: R0         ! [m]           Ice margin radius at t=0
    real(dp), intent(in   ) :: lambda     ! [ ]           Mass balance parameter (lambda = 5.0 gives a nice growing ice sheet)
    real(dp), intent(in   ) :: x          ! [m]           x-coordinate
    real(dp), intent(in   ) :: y          ! [m]           y-coordinate
    real(dp), intent(in   ) :: t          ! [yr]          time
    real(dp), intent(  out) :: H          ! [m]           Ice thickness at [x,y,t]
    real(dp), intent(  out) :: M          ! [m/yr]        Accumulation rate at [x,y,t]

    ! Local variables
    real(dp) :: alpha, beta, Gamma, f1, f2, t0, tp, f3, f4

    alpha = (2._dp - (n+1._dp)*lambda) / ((5._dp*n)+3._dp)
    beta  = (1._dp + ((2._dp*n)+1._dp)*lambda) / ((5._dp*n)+3._dp)
    Gamma = 2._dp/5._dp * (A/sec_per_year) * (ice_density * grav)**n

    f1 = ((2._dp*n)+1)/(n+1._dp)
    f2 = (R0**(n+1._dp))/(H0**((2._dp*n)+1._dp))
    t0 = (beta / Gamma) * (f1**n) * f2

    tp = t * sec_per_year

    f1 = (tp / t0)**(-alpha)
    f2 = (tp / t0)**(-beta)
    f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
    f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
    H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))

    M = (lambda / tp) * H * sec_per_year

  end subroutine Bueler_dome

end module Bueler_SIA_solution