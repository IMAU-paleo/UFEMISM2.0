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
    procedure :: H

  end type type_Halfar_solution

  type(type_Halfar_solution) :: Halfar

contains

  real(dp) function Gamma( A,n)
    real(dp), intent(in) :: A,n
    Gamma = (2._dp / 5._dp) * (A / sec_per_year) * (ice_density * grav)**n
  end function Gamma
  real(dp) function t0( A,n,H0,R0)
    real(dp), intent(in) :: A,n,H0,R0
    t0 = 1._dp / ((5._dp * n + 3._dp) * Gamma( A,n)) * ((2._dp * n + 1._dp)/(n + 1._dp))**n * (R0**(n + 1._dp))/(H0**(2._dp * n  + 1))
  end function t0

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

  real(dp) function f1( A,n,H0,R0,t)
    real(dp), intent(in) :: A,n,H0,R0,t
    f1 = ((t0( A,n,H0,R0) + t * sec_per_year) / t0( A,n,H0,R0))**p1( n)
  end function f1
  real(dp) function f2( A,n,H0,R0,t)
    real(dp), intent(in) :: A,n,H0,R0,t
    f2 = ((t0( A,n,H0,R0) + t * sec_per_year) / t0( A,n,H0,R0))**p2( n)
  end function f2
  real(dp) function f3( R0,x,y)
    real(dp), intent(in) :: R0,x,y
    f3 = r(x,y) / R0
  end function f3

  real(dp) function G( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    G = 1._dp - (min( 1._dp, f2( A,n,H0,R0,t) * f3( R0,x,y)))**p3( n)
  end function G

  real(dp) function H( self,A,n,H0,R0,x,y,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,t
    H = H0 * f1( A,n,H0,R0,t) * G( A,n,H0,R0,x,y,t)**p4( n)
  end function H

end module Halfar_SIA_solution