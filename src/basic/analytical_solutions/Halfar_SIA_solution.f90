module Halfar_SIA_solution

  use precisions, only: dp
  use parameters, only: sec_per_year, ice_density, grav

  implicit none

  private

  public :: Halfar

  type type_Halfar_solution
    ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity
    ! function with dome thickness H0 and margin radius R0 at t0.
    !
    ! Halfar, P.: On the Dynamics of Ice Sheets, Journal of Geophysical Research 86,
    !   11065-11072, 1981
    !
    ! See also the UFEMISM 1.0 paper (2021) for a more human-readable version

    ! Extended by Tijn to include expressions for the ice velocities and thinning rates

  contains
    procedure :: H => calc_Halfar_ice_thickness

  end type type_Halfar_solution

  type(type_Halfar_solution) :: Halfar

contains

  subroutine calc_Halfar_ice_thickness( self, A, n, H0, R0, x, y, t, H)
    ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity
    ! function with dome thickness H0 and margin radius R0 at t0.
    !
    ! Halfar, P.: On the Dynamics of Ice Sheets, Journal of Geophysical Research 86,
    !   11065-11072, 1981

    ! In/output variables:
    class(type_Halfar_solution), intent(in   ) :: self
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

    r = sqrt(x**2._dp + y**2._dp)

    f1 = (t0/tp)**(2._dp/(5._dp * n + 3._dp))
    f2 = (t0/tp)**(1._dp/(5._dp * n + 3._dp))
    f3 = (r/R0)

    H = H0 * f1 * max( 0._dp, (1._dp - (f2*f3)**((n + 1._dp)/n)))**(n/(2._dp * n + 1._dp))

  end subroutine calc_Halfar_ice_thickness

  real(dp) function r( x,y)
    real(dp), intent(in) :: x,y
    r = sqrt( x**2 + y**2)
  end function r

  real(dp) function p1( n)
    real(dp), intent(in) :: n
    p1 = -2._dp / (5._dp*n + 3._dp)
  end function p1
  real(dp) function p2( n)
    real(dp), intent(in) :: n
    p2 = -1._dp / (5._dp*n + 3._dp)
  end function p2
  real(dp) function p3( n)
    real(dp), intent(in) :: n
    p3 = (n+1._dp) / n
  end function p3
  real(dp) function p4( n)
    real(dp), intent(in) :: n
    p4 = n / (2._dp*n + 1._dp)
  end function p4

  real(dp) function f1( n,t0,t)
    real(dp), intent(in) :: n,t0,t
    f1 = ((t0 + t * sec_per_year) / t0)**p1( n)
  end function f1
  real(dp) function f2( n,t0,t)
    real(dp), intent(in) :: n,t0,t
    f2 = ((t0 + t * sec_per_year) / t0)**p2( n)
  end function f2
  real(dp) function f3( R0,x,y)
    real(dp), intent(in) :: R0,x,y
    f3 = r(x,y) / R0
  end function f3

  real(dp) function G( n,R0,t0,x,y,t)
    real(dp), intent(in) :: n,R0,t0,x,y,t
    G = 1._dp - (f2( n,t0,t) * f3( R0,x,y))**p3( n)
  end function G

  real(dp) function H_new( n,H0,R0,t0,x,y,t)
    real(dp), intent(in) :: n,H0,R0,t0,x,y,t
    H_new = H0 * f1( n,t0,t) * G( n,R0,t0,x,y,t)**p4( n)
  end function H_new

end module Halfar_SIA_solution