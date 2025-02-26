module kelvin_function

  use precisions, only: dp
  use parameters
  
  implicit none
  
  private
  
  public :: kelvin
  
contains

  ! The Kelvin function (just a wrapper for the Zhang & Jin "klvna" subroutine)
  function kelvin(x) RESULT(kei)
    ! Evaluate the Kelvin function in a given point at a distance x from the load. Currently this is implemented
    ! using the klvna routine in mklvna.f (see that routine for details).

    !IMPLICIT NONE

    ! Input variables:
    ! x is the distance between a considered load and the point in which we want to know the deformation, x is in units Lr
    REAL(dp), INTENT(IN) :: x

    ! Result variables:
    REAL(dp)             :: kei

    ! Local variables:
    REAL(dp)     :: xdp,berdp,beidp,gerdp,geidp,derdp,deidp,herdp,heidp

    xdp = x
    CALL klvna( xdp,berdp,beidp,gerdp,geidp,derdp,deidp,herdp,heidp)
    kei = geidp

  END function kelvin
  SUBROUTINE klvna ( x, ber, bei, ger, gei, der, dei, her, hei )

!*****************************************************************************80
!
!! KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real (dp) X, the argument.
!
!    Output, real (dp) BER, BEI, GER, GEI, DER, DEI, HER, HEI,
!    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
!

!  USE parameters_module, ONLY: pi

!  IMPLICIT NONE

  ! In/output variables:
  REAL(dp), INTENT(IN )   :: x
  REAL(dp), INTENT(OUT)   :: ber
  REAL(dp), INTENT(OUT)   :: bei
  REAL(dp), INTENT(OUT)   :: ger
  REAL(dp), INTENT(OUT)   :: gei
  REAL(dp), INTENT(OUT)   :: der
  REAL(dp), INTENT(OUT)   :: dei
  REAL(dp), INTENT(OUT)   :: her
  REAL(dp), INTENT(OUT)   :: hei

  ! Local variables:
  REAL(dp)   :: cn0
  REAL(dp)   :: cp0
  REAL(dp)   :: cs
  REAL(dp)   :: el
  REAL(dp)   :: eps
  REAL(dp)   :: fac
  REAL(dp)   :: gs

  INTEGER    :: k
  INTEGER    :: Km
  INTEGER    :: m

  REAL(dp)   :: pn0
  REAL(dp)   :: pn1
  REAL(dp)   :: pp0
  REAL(dp)   :: pp1
  REAL(dp)   :: qn0
  REAL(dp)   :: qn1
  REAL(dp)   :: qp0
  REAL(dp)   :: qp1
  REAL(dp)   :: r
  REAL(dp)   :: r0
  REAL(dp)   :: r1
  REAL(dp)   :: rc
  REAL(dp)   :: rs
  REAL(dp)   :: sn0
  REAL(dp)   :: sp0
  REAL(dp)   :: ss
  REAL(dp)   :: x2
  REAL(dp)   :: x4
  REAL(dp)   :: xc1
  REAL(dp)   :: xc2
  REAL(dp)   :: xd
  REAL(dp)   :: xe1
  REAL(dp)   :: xe2
  REAL(dp)   :: xt

  el = 0.5772156649015329_dp
  eps = 1.0D-15

  if ( x == 0.0_dp ) then
    ber = 1.0_dp
    bei = 0.0_dp
    ger = 1.0D+300
    gei = -0.25_dp * pi
    der = 0.0_dp
    dei = 0.0_dp
    her = -1.0D+300
    hei = 0.0_dp
    return
  end if

  x2 = 0.25_dp * x * x
  x4 = x2 * x2

  if ( abs ( x ) < 10.0_dp ) then

    ber = 1.0_dp
    r = 1.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) ** 2 * x4
      ber = ber + r
      if ( abs ( r ) < abs ( ber ) * eps ) then
        exit
      end if
    end do

    bei = x2
    r = x2
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      bei = bei + r
      if ( abs ( r ) < abs ( bei ) * eps ) then
        exit
      end if
    end do

    ger = - ( log ( x / 2.0_dp ) + el ) * ber + 0.25_dp * pi * bei
    r = 1.0_dp
    gs = 0.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m - 1.0_dp ) + 1.0_dp / ( 2.0_dp * m )
      ger = ger + r * gs
      if ( abs ( r * gs ) < abs ( ger ) * eps ) then
        exit
      end if
    end do

    gei = x2 - ( log ( x / 2.0_dp ) + el ) * bei - 0.25_dp * pi * ber
    r = x2
    gs = 1.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m ) + 1.0_dp / ( 2.0_dp * m + 1.0_dp )
      gei = gei + r * gs
      if ( abs ( r * gs ) < abs ( gei ) * eps ) then
        exit
      end if
    end do

    der = -0.25_dp * x * x2
    r = der
    do m = 1, 60
      r = -0.25_dp * r / m / ( m + 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      der = der + r
      if ( abs ( r ) < abs ( der ) * eps ) then
        exit
      end if
    end do

    dei = 0.5_dp * x
    r = dei
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) * x4
      dei = dei + r
      if ( abs ( r ) < abs ( dei ) * eps ) then
        exit
      end if
    end do

    r = -0.25_dp * x * x2
    gs = 1.5_dp
    her = 1.5_dp * r - ber / x &
      - ( log ( x / 2.0_dp ) + el ) * der + 0.25_dp * pi * dei
    do m = 1, 60
      r = -0.25_dp * r / m / ( m + 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2 * m + 1.0_dp ) + 1.0_dp &
        / ( 2 * m + 2.0_dp )
      her = her + r * gs
      if ( abs ( r * gs ) < abs ( her ) * eps ) then
        exit
      end if
    end do

    r = 0.5_dp * x
    gs = 1.0_dp
    hei = 0.5_dp * x - bei / x &
      - ( log ( x / 2.0_dp ) + el ) * dei - 0.25_dp * pi * der
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2 * m - 1.0_dp ) &
        / ( 2 * m + 1.0_dp ) * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m ) + 1.0_dp &
        / ( 2 * m + 1.0_dp )
      hei = hei + r * gs
      if ( abs ( r * gs ) < abs ( hei ) * eps ) then
        return
      end if
    end do

  else

    pp0 = 1.0_dp
    pn0 = 1.0_dp
    qp0 = 0.0_dp
    qn0 = 0.0_dp
    r0 = 1.0_dp

    if ( abs ( x ) < 40.0_dp ) then
      km = 18
    else
      km = 10
    end if

    fac = 1.0_dp
    do k = 1, km
      fac = -fac
      xt = 0.25_dp * k * pi - int ( 0.125_dp * k ) * 2.0_dp * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r0 = 0.125_dp * r0 * ( 2.0_dp * k - 1.0_dp ) ** 2 / k / x
      rc = r0 * cs
      rs = r0 * ss
      pp0 = pp0 + rc
      pn0 = pn0 + fac * rc
      qp0 = qp0 + rs
      qn0 = qn0 + fac * rs
    end do

    xd = x / sqrt (2.0_dp )
    xe1 = exp ( xd )
    xe2 = exp ( - xd )
    xc1 = 1.0_dp / sqrt ( 2.0_dp * pi * x )
    xc2 = sqrt ( 0.5_dp * pi / x )
    cp0 = cos ( xd + 0.125_dp * pi )
    cn0 = cos ( xd - 0.125_dp * pi )
    sp0 = sin ( xd + 0.125_dp * pi )
    sn0 = sin ( xd - 0.125_dp * pi )
    ger = xc2 * xe2 * (  pn0 * cp0 - qn0 * sp0 )
    gei = xc2 * xe2 * ( -pn0 * sp0 - qn0 * cp0 )
    ber = xc1 * xe1 * (  pp0 * cn0 + qp0 * sn0 ) - gei / pi
    bei = xc1 * xe1 * (  pp0 * sn0 - qp0 * cn0 ) + ger / pi
    pp1 = 1.0_dp
    pn1 = 1.0_dp
    qp1 = 0.0_dp
    qn1 = 0.0_dp
    r1 = 1.0_dp
    fac = 1.0_dp

    do k = 1, km
      fac = -fac
      xt = 0.25_dp * k * pi - int ( 0.125_dp * k ) * 2.0_dp * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r1 = 0.125_dp * r1 &
        * ( 4.0_dp - ( 2.0_dp * k - 1.0_dp ) ** 2 ) / k / x
      rc = r1 * cs
      rs = r1 * ss
      pp1 = pp1 + fac * rc
      pn1 = pn1 + rc
      qp1 = qp1 + fac * rs
      qn1 = qn1 + rs
    end do

    her = xc2 * xe2 * ( - pn1 * cn0 + qn1 * sn0 )
    hei = xc2 * xe2 * (   pn1 * sn0 + qn1 * cn0 )
    der = xc1 * xe1 * (   pp1 * cp0 + qp1 * sp0 ) - hei / pi
    dei = xc1 * xe1 * (   pp1 * sp0 - qp1 * cp0 ) + her / pi

  end if

  END SUBROUTINE klvna

end module kelvin_function
