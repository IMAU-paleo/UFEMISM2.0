module Halfar_SIA_solution
  ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity
  ! function with dome thickness H0 and margin radius R0 at t0.

  ! Halfar, P.: On the Dynamics of Ice Sheets, Journal of Geophysical Research 86,
  !   11065-11072, 1981

  ! See also the UFEMISM 1.0 paper (2021) for a more human-readable version

  ! Extended by Tijn to include expressions for the ice velocities and thinning rates

  use precisions, only: dp
  use parameters, only: sec_per_year, ice_density, grav

  implicit none

  private

  public :: Halfar

  type type_Halfar_solution
    private
  contains
    procedure, public :: H, dH_dt, dH_dx, dH_dy
    procedure, public :: u, v, w
    procedure, public :: u_vav, v_vav
  end type type_Halfar_solution

  type(type_Halfar_solution) :: Halfar

contains

  real(dp) function H( self,A,n,H0,R0,x,y,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,t
    H = H0 * f1( A,n,H0,R0,t) * G( A,n,H0,R0,x,y,t)**p4( n)
  end function H
  real(dp) function dH_dt( self,A,n,H0,R0,x,y,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,t
    if (r(x,y) <= R0) then
      dH_dt = H0 * (&
        df1_dt( A,n,H0,R0,t) *          G( A,n,H0,R0,x,y,t)** p4( n) &
       + f1   ( A,n,H0,R0,t) * p4( n) * G( A,n,H0,R0,x,y,t)**(p4( n) - 1._dp) * dG_dt( A,n,H0,R0,x,y,t))
    else
      dH_dt = 0._dp
    end if
  end function dH_dt
  real(dp) function dH_dx( self,A,n,H0,R0,x,y,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,t
    dH_dx = Q( A,n,H0,R0,x,y,t) * x
  end function dH_dx
  real(dp) function dH_dy( self,A,n,H0,R0,x,y,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,t
    dH_dy = Q( A,n,H0,R0,x,y,t) * y
  end function dH_dy

  real(dp) function u( self,A,n,H0,R0,x,y,z,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,z,t
    if (r( x,y) <= R0) then
      u = -c( A,n) * abs( Q( A,n,H0,R0,x,y,t))**n * r( x,y)**(n-1._dp) * D_m( A,n,H0,R0,x,y,z,t,n+1._dp) * x
    else
      u = 0._dp
    end if
  end function u
  real(dp) function v( self,A,n,H0,R0,x,y,z,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,z,t
    if (r( x,y) <= R0) then
      v = -c( A,n) * abs( Q( A,n,H0,R0,x,y,t))**n * r( x,y)**(n-1._dp) * D_m( A,n,H0,R0,x,y,z,t,n+1._dp) * y
    else
      v = 0._dp
    end if
  end function v
  real(dp) function w( self,A,n,H0,R0,x,y,z,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,z,t
    real(dp)                                :: term1, term2
    if (r( x,y) <= R0) then
      term1 = -(U1( A,n,H0,R0,x,y,t) + V1( A,n,H0,R0,x,y,t)) * (z * Halfar%H( A,n,H0,R0,x,y,t)**(n+1._dp) - D_m( A,n,H0,R0,x,y,z,t,n+2._dp)/(n+2._dp))
      term2 = -(U2( A,n,H0,R0,x,y,t) + V2( A,n,H0,R0,x,y,t)) * (z * Halfar%H( A,n,H0,R0,x,y,t)**n         - D_m( A,n,H0,R0,x,y,z,t,n+1._dp)/(n+1._dp))
      w = term1 + term2
    else
      w = 0._dp
    end if
  end function w

  real(dp) function u_vav( self,A,n,H0,R0,x,y,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,t
    if (r( x,y) <= R0) then
      u_vav = -c( A,n) * abs( Q( A,n,H0,R0,x,y,t))**n * r( x,y)**(n-1._dp) * D_m_vav( A,n,H0,R0,x,y,t,n+1._dp) * x
    else
      u_vav = 0._dp
    end if
  end function u_vav
  real(dp) function v_vav( self,A,n,H0,R0,x,y,t)
    class(type_Halfar_solution), intent(in) :: self
    real(dp),                    intent(in) :: A,n,H0,R0,x,y,t
    if (r( x,y) <= R0) then
      v_vav = -c( A,n) * abs( Q( A,n,H0,R0,x,y,t))**n * r( x,y)**(n-1._dp) * D_m_vav( A,n,H0,R0,x,y,t,n+1._dp) * y
    else
      v_vav = 0._dp
    end if
  end function v_vav

  real(dp) function Halfar_Gamma( A,n)
    real(dp), intent(in) :: A,n
    Halfar_Gamma = (2._dp / 5._dp) * (A / sec_per_year) * (ice_density * grav)**n
  end function Halfar_Gamma
  real(dp) function t0( A,n,H0,R0)
    real(dp), intent(in) :: A,n,H0,R0
    t0 = 1._dp / ((5._dp * n + 3._dp) * Halfar_Gamma( A,n)) * ((2._dp * n + 1._dp)/(n + 1._dp))**n * (R0**(n + 1._dp))/(H0**(2._dp * n  + 1))
  end function t0

  real(dp) function r( x,y)
    real(dp), intent(in) :: x,y
    r = sqrt( x**2 + y**2)
  end function r
  real(dp) function dr_dx( x,y)
    real(dp), intent(in) :: x,y
    dr_dx = x / r( x,y)
  end function dr_dx
  real(dp) function dr_dy( x,y)
    real(dp), intent(in) :: x,y
    dr_dy = y / r( x,y)
  end function dr_dy

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

  real(dp) function df1_dt( A,n,H0,R0,t)
    real(dp), intent(in) :: A,n,H0,R0,t
    df1_dt = p1( n) * sec_per_year / t0( A,n,H0,R0) * ((t0( A,n,H0,R0) + t * sec_per_year) / t0( A,n,H0,R0))**(p1( n)-1._dp)
  end function df1_dt
  real(dp) function df2_dt( A,n,H0,R0,t)
    real(dp), intent(in) :: A,n,H0,R0,t
    df2_dt = p2( n) * sec_per_year / t0( A,n,H0,R0) * ((t0( A,n,H0,R0) + t * sec_per_year) / t0( A,n,H0,R0))**(p2( n)-1._dp)
  end function df2_dt
  real(dp) function df3_dx( R0,x,y)
    real(dp), intent(in) :: R0,x,y
    df3_dx = dr_dx( x,y) / R0
  end function df3_dx
  real(dp) function df3_dy( R0,x,y)
    real(dp), intent(in) :: R0,x,y
    df3_dy = dr_dy( x,y) / R0
  end function df3_dy

  real(dp) function G( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    G = 1._dp - (min( 1._dp, f2( A,n,H0,R0,t) * f3( R0,x,y)))**p3( n)
  end function G
  real(dp) function dG_dx( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    dG_dx = -p3( n) * f2( A,n,H0,R0,t)**p3( n) * f3( R0,x,y)**(p3( n)-1._dp) * df3_dx( R0,x,y)
  end function dG_dx
  real(dp) function dG_dy( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    dG_dy = -p3( n) * f2( A,n,H0,R0,t)**p3( n) * f3( R0,x,y)**(p3( n)-1._dp) * df3_dy( R0,x,y)
  end function dG_dy
  real(dp) function dG_dt( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    dG_dt = -p3( n) * f2( A,n,H0,R0,t)**(p3( n) - 1._dp) * df2_dt( A,n,H0,R0,t) * f3( R0,x,y)**p3( n)
  end function dG_dt

  real(dp) function Q0( A,n,H0,R0,t)
    real(dp), intent(in) :: A,n,H0,R0,t
    Q0 = -H0/R0 * p3( n) * p4( n) * f1( A,n,H0,R0,t) * f2( A,n,H0,R0,t)**p3( n)
  end function Q0
  real(dp) function Q( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    Q = Q0( A,n,H0,R0,t) / r( x,y) * f3( R0,x,y)**(p3( n)-1._dp) * G( A,n,H0,R0,x,y,t)**(p4( n)-1)
  end function Q
  real(dp) function dQ_dx( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    real(dp)             :: term1, term2, term3
    term1 = -1 / r( x,y)**2 * dr_dx( x,y) * f3( R0,x,y)**(p3( n)-1._dp) * G( A,n,H0,R0,x,y,t)**(p4( n)-1._dp)
    term2 =  1 / r( x,y) * (p3( n)-1._dp) * f3( R0,x,y)**(p3( n)-2._dp) * df3_dx( R0,x,y) * G( A,n,H0,R0,x,y,t)**(p4( n)-1._dp)
    term3 =  1 / r( x,y) * f3( R0,x,y)**(p3( n)-1._dp) * (p4( n)-1._dp) * G( A,n,H0,R0,x,y,t)**(p4( n)-2._dp) * dG_dx( A,n,H0,R0,x,y,t)
    dQ_dx = Q0( A,n,H0,R0,t) * (term1 + term2 + term3)
  end function dQ_dx
  real(dp) function dQ_dy( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    real(dp)             :: term1, term2, term3
    term1 = -1 / r( x,y)**2 * dr_dy( x,y) * f3( R0,x,y)**(p3( n)-1._dp) * G( A,n,H0,R0,x,y,t)**(p4( n)-1._dp)
    term2 =  1 / r( x,y) * (p3( n)-1._dp) * f3( R0,x,y)**(p3( n)-2._dp) * df3_dy( R0,x,y) * G( A,n,H0,R0,x,y,t)**(p4( n)-1._dp)
    term3 =  1 / r( x,y) * f3( R0,x,y)**(p3( n)-1._dp) * (p4( n)-1._dp) * G( A,n,H0,R0,x,y,t)**(p4( n)-2._dp) * dG_dy( A,n,H0,R0,x,y,t)
    dQ_dy = Q0( A,n,H0,R0,t) * (term1 + term2 + term3)
  end function dQ_dy

  real(dp) function D_m( A,n,H0,R0,x,y,z,t,m)
    real(dp), intent(in) :: A,n,H0,R0,x,y,z,t,m
    D_m = Halfar%H( A,n,H0,R0,x,y,t)**m - (Halfar%H( A,n,H0,R0,x,y,t) - z)**m
  end function D_m
  real(dp) function int_D_m_dz( A,n,H0,R0,x,y,t,m)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t,m
    int_D_m_dz = (1._dp - 1._dp / (m + 1._dp)) * Halfar%H( A,n,H0,R0,x,y,t)**(m+1._dp)
  end function int_D_m_dz
  real(dp) function D_m_vav( A,n,H0,R0,x,y,t,m)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t,m
    D_m_vav = int_D_m_dz( A,n,H0,R0,x,y,t,m) / Halfar%H( A,n,H0,R0,x,y,t)
  end function D_m_vav

  real(dp) function c( A,n)
    real(dp), intent(in) :: A,n
    c = -2._dp * A / (n+1._dp) * (ice_density * grav)**n
  end function c

  real(dp) function U1( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    real(dp)             :: term1, term2, term3
    term1 = n * Q( A,n,H0,R0,x,y,t)**(n-1._dp) * dQ_dx( A,n,H0,R0,x,y,t) * r( x,y)**(n-1._dp) * x
    term2 = Q( A,n,H0,R0,x,y,t)**n * (n-1._dp) * r( x,y)**(n-2._dp) * dr_dx( x,y) * x
    term3 = Q( A,n,H0,R0,x,y,t)**n * r( x,y)**(n-1._dp)
    U1 = c( A,n) * (term1 + term2 + term3)
  end function U1
  real(dp) function V1( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    real(dp)             :: term1, term2, term3
    term1 = n * Q( A,n,H0,R0,x,y,t)**(n-1._dp) * dQ_dy( A,n,H0,R0,x,y,t) * r( x,y)**(n-1._dp) * y
    term2 = Q( A,n,H0,R0,x,y,t)**n * (n-1._dp) * r( x,y)**(n-2._dp) * dr_dy( x,y) * y
    term3 = Q( A,n,H0,R0,x,y,t)**n * r( x,y)**(n-1._dp)
    V1 = c( A,n) * (term1 + term2 + term3)
  end function V1
  real(dp) function U2( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    U2 = c( A,n) * Q( A,n,H0,R0,x,y,t)**n * r( x,y)**(n-1._dp) * x * (n+1._dp) * Halfar%dH_dx( A,n,H0,R0,x,y,t)
  end function U2
  real(dp) function V2( A,n,H0,R0,x,y,t)
    real(dp), intent(in) :: A,n,H0,R0,x,y,t
    V2 = c( A,n) * Q( A,n,H0,R0,x,y,t)**n * r( x,y)**(n-1._dp) * y * (n+1._dp) * Halfar%dH_dy( A,n,H0,R0,x,y,t)
  end function V2

end module Halfar_SIA_solution