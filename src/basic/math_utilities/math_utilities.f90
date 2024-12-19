module math_utilities

  ! Some generally useful tools, and basic mathematical functions

  use tests_main
  use assertions_basic
  use mpi
  use precisions, only: dp
  use mpi_basic, only: par, cerr, ierr
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use parameters
  use reallocate_mod, only: reallocate

  implicit none

contains

! == Floatation criterion, surface elevation, and thickness above floatation

  pure function is_floating( Hi, Hb, SL) result( isso)
    ! The floatation criterion

    ! Input variables:
    real(dp), intent(in) :: Hi     ! [m] Ice thickness
    real(dp), intent(in) :: Hb     ! [m] Bedrock elevation
    real(dp), intent(in) :: SL     ! [m] Water surface elevation

    ! Output variables:
    logical :: isso   ! Whether or not the ice will float

    isso = .false.
    if (Hi < (SL - Hb) * seawater_density/ice_density) isso = .true.

  end function is_floating

  pure function ice_surface_elevation( Hi, Hb, SL) result( Hs)
    ! The ice surface elevation equation

    ! Input variables:
    real(dp), intent(in) :: Hi     ! [m] Ice thickness
    real(dp), intent(in) :: Hb     ! [m] Bedrock elevation
    real(dp), intent(in) :: SL     ! [m] Water surface elevation

    ! Output variables:
    real(dp) :: Hs     ! [m] Ice surface elevation

    Hs = Hi + max( SL - ice_density / seawater_density * Hi, Hb)

  end function ice_surface_elevation

  pure function thickness_above_floatation( Hi, Hb, SL) result( TAF)
    ! The thickness-above-floatation equation

    ! Input variables:
    real(dp), intent(in) :: Hi     ! [m] Ice thickness
    real(dp), intent(in) :: Hb     ! [m] Bedrock elevation
    real(dp), intent(in) :: SL     ! [m] Water surface elevation

    ! Output variables:
    real(dp) :: TAF    ! [m] Ice thickness above floatation

    TAF = Hi - max( 0._dp, (SL - Hb) * (seawater_density / ice_density))

  end function thickness_above_floatation

  function Hi_from_Hb_Hs_and_SL( Hb, Hs, SL) result( Hi)
    ! Calculate Hi from Hb, Hs, and SL

    ! Input variables:
    real(dp) , intent(in) :: Hb     ! [m] Bedrock elevation
    real(dp) , intent(in) :: Hs     ! [m] Surface elevation
    real(dp) , intent(in) :: SL     ! [m] Water surface elevation

    ! Output variables:
    real(dp) :: Hi     ! [m] Ice thickness

    ! Local variables:
    real(dp) :: Hi_float ! [m] Maximum floating ice thickness
    real(dp) :: Hs_float ! [m] Surface elevation of maximum floating ice thickness

    Hi_float = max( 0._dp, (SL - Hb) * (seawater_density / ice_density))
    Hs_float = Hb + Hi_float

    if (Hs > Hs_float) then
      Hi = Hs - Hb
    else
      Hi = min( Hi_float, (Hs - SL) / (1._dp - (ice_density / seawater_density)) )
    end if

  end function Hi_from_Hb_Hs_and_SL

! == The error function

  pure function error_function( X) result( ERR)
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

! == The oblique stereographic projection

  pure subroutine oblique_sg_projection( lambda, phi, lambda_M_deg, phi_M_deg, beta_deg, x, y)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates to a rectangular coordinate system, with coordinates (x,y).
    !
    ! For more information about M, beta_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    ! For North and South Pole: lambda_M_deg = 0._dp, to generate the correct coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Input variables:
    real(dp), intent(in)  :: lambda         ! [degrees east ] Longitude
    real(dp), intent(in)  :: phi            ! [degrees north] Latitude

    ! Polar stereographic projection parameters
    real(dp), intent(in)  :: lambda_M_deg   ! [degrees east]  Longitude of the pole of the oblique stereographic projection
    real(dp), intent(in)  :: phi_M_deg      ! [degrees north] Latitude  of the pole of the oblique stereographic projection
    real(dp), intent(in)  :: beta_deg       ! [degrees]       Standard parallel     of the oblique stereographic projection

    ! Output variables:
    real(dp), intent(out) :: x              ! [m] x-coordinate
    real(dp), intent(out) :: y              ! [m] y-coordinate

    ! Local variables:
    real(dp) :: alpha_deg                   ! [degrees]
    real(dp) :: phi_P                       ! [radians]
    real(dp) :: lambda_P                    ! [radians]
    real(dp) :: t_P_prime
    real(dp) :: lambda_M, phi_M, alpha

    ! Convert beta to alpha
    alpha_deg = 90._dp - beta_deg

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180._dp) * phi
    lambda_P = (pi / 180._dp) * lambda

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1._dp + cos(alpha)) / (1._dp + cos(phi_P) * cos(phi_M) * cos(lambda_P - lambda_M) + sin(phi_P) * sin(phi_M))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x =  earth_radius * (cos(phi_P) * sin(lambda_P - lambda_M)) * t_P_prime
    y =  earth_radius * (sin(phi_P) * cos(phi_M) - (cos(phi_P) * sin(phi_M)) * cos(lambda_P - lambda_M)) * t_P_prime

  end subroutine oblique_sg_projection

  pure subroutine inverse_oblique_sg_projection( x, y, lambda_M_deg, phi_M_deg, beta_deg, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates to a longitude-latitude coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, alpha_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    ! Input variables:
    real(dp), intent(in)  :: x              ! [m] x-coordinate
    real(dp), intent(in)  :: y              ! [m] y-coordinate

    ! Polar stereographic projection parameters
    real(dp), intent(in)  :: lambda_M_deg   ! [degrees east]  Longitude of the pole of the oblique stereographic projection
    real(dp), intent(in)  :: phi_M_deg      ! [degrees north] Latitude  of the pole of the oblique stereographic projection
    real(dp), intent(in)  :: beta_deg       ! [degrees]       Standard parallel     of the oblique stereographic projection

    ! Output variables:
    real(dp), intent(out) :: lambda_P       ! [degrees east ] Longitude
    real(dp), intent(out) :: phi_P          ! [degrees north] Latitude

    ! Local variables:
    real(dp) :: alpha_deg                   ! [degrees]
    real(dp) :: x_3D_P_prime                ! [m]
    real(dp) :: y_3D_P_prime                ! [m]
    real(dp) :: z_3D_P_prime                ! [m]
    real(dp) :: a
    real(dp) :: t_P
    real(dp) :: x_3D_P                      ! [m]
    real(dp) :: y_3D_P                      ! [m]
    real(dp) :: z_3D_P                      ! [m]
    real(dp) :: lambda_M, phi_M, alpha

    ! Convert beta to alpha
    alpha_deg = 90._dp - beta_deg

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = earth_radius * cos(alpha) * cos(lambda_M) * cos(phi_M) - sin(lambda_M) * x - cos(lambda_M) * sin(phi_M) * y
    y_3D_P_prime = earth_radius * cos(alpha) * sin(lambda_M) * cos(phi_M) + cos(lambda_M) * x - sin(lambda_M) * sin(phi_M) * y
    z_3D_P_prime = earth_radius * cos(alpha) *                 sin(phi_M)                     +                 cos(phi_M) * y

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = cos(lambda_M) * cos(phi_M) * x_3D_P_prime  +  sin(lambda_M) * cos(phi_M) * y_3D_P_prime  +  sin(phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * earth_radius**2 + 2._dp * earth_radius * a) / (earth_radius**2 + 2._dp * earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  earth_radius * cos(lambda_M) * cos(phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  earth_radius * sin(lambda_M) * cos(phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  earth_radius *                 sin(phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    if     (x_3D_P <  0._dp                      ) then
      lambda_P = 180._dp + (180._dp / pi) * atan(y_3D_P / x_3D_P)
    elseif (x_3D_P >  0._dp .and. y_3D_P >= 0._dp) then
      lambda_P =           (180._dp / pi) * atan(y_3D_P / x_3D_P)
    elseif (x_3D_P >  0._dp .and. y_3D_P <  0._dp) then
      lambda_P = 360._dp + (180._dp / pi) * atan(y_3D_P / x_3D_P)
    elseif (x_3D_P == 0._dp .and. y_3D_P >  0._dp) then
      lambda_P =  90._dp
    elseif (x_3D_P == 0._dp .and. y_3D_P <  0._dp) then
      lambda_P = 270._dp
    elseif (x_3D_P == 0._dp .and. y_3D_P == 0._dp) then
      lambda_P =   0._dp
    end if

    ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
    if     (x_3D_P /= 0._dp .or. y_3D_P /= 0._dp) then
      phi_P = (180._dp / pi) * atan(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2))
    elseif (z_3D_P >  0._dp) then
      phi_P =   90._dp
    elseif (z_3D_P <  0._dp) then
      phi_P =  -90._dp
    end if

  end subroutine inverse_oblique_sg_projection

! == Line integrals used in conservative remapping

  pure function line_integral_xdy(   p, q, tol_dist) result( I_pq)
    ! Calculate the line integral x dy from p to q

    ! Input variables:
    real(dp), dimension(2)                             , intent(in)    :: p, q
    real(dp)                                           , intent(in)    :: tol_dist

    ! Output variables:
    real(dp)                                                           :: I_pq

    ! Local variables:
    real(dp)                                                           :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    if (abs( yp-yq) < tol_dist) then
      I_pq = 0._dp
      return
    end if

    dx = q( 1) - p( 1)
    dy = q( 2) - p( 2)

    I_pq = xp*dy - yp*dx + (dx / (2._dp*dy)) * (yq**2 - yp**2)

  end function line_integral_xdy

  pure function line_integral_mxydx( p, q, tol_dist) result( I_pq)
    ! Calculate the line integral -xy dx from p to q

    ! Input variables:
    real(dp), dimension(2), intent(in)    :: p, q
    real(dp)              , intent(in)    :: tol_dist

    ! Output variables:
    real(dp) :: I_pq

    ! Local variables:
    real(dp) :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    if (abs( xp-xq) < tol_dist) then
      I_pq = 0._dp
      return
    end if

    dx = q( 1) - p( 1)
    dy = q( 2) - p( 2)

    I_pq = (1._dp/2._dp * (xp*dy/dx - yp) * (xq**2-xp**2)) - (1._dp/3._dp * dy/dx * (xq**3-xp**3))

  end function line_integral_mxydx

  pure function line_integral_xydy(  p, q, tol_dist) result( I_pq)
    ! Calculate the line integral xy dy from p to q

    ! Input variables:
    real(dp), dimension(2), intent(in)    :: p, q
    real(dp)              , intent(in)    :: tol_dist

    ! Output variables:
    real(dp) :: I_pq

    ! Local variables:
    real(dp) :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    if (abs( yp-yq) < tol_dist) then
      I_pq = 0._dp
      return
    end if

    dx = q( 1 ) - p( 1)
    dy = q( 2 ) - p( 2)

    I_pq = (1._dp/2._dp * (xp - yp*dx/dy) * (yq**2-yp**2)) + (1._dp/3._dp * dx/dy * (yq**3-yp**3))

  end function line_integral_xydy

! == Finding inverses of some small matrices

  pure function calc_determinant_2_by_2( A) result( detA)
    ! Determinant of a 2-by-2 matrix

    ! Input variables:
    real(dp), dimension(2,2), intent(in) :: A

    ! Output variables:
    real(dp) :: detA

    ! Calculate the determinant of A
    detA = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

  end function calc_determinant_2_by_2

  pure function calc_determinant_3_by_3( A) result( detA)
    ! Determinant of a 2-by-2 matrix

    ! Input variables:
    real(dp), dimension(3,3), intent(in) :: A

    ! Output variables:
    real(dp) :: detA

    ! Local variables:
    real(dp), dimension(3,3) :: Ainv

    ! Calculate the minors of A
    Ainv( 1,1) = A( 2,2) * A( 3,3) - A( 2,3) * A( 3,2)
    Ainv( 1,2) = A( 2,1) * A( 3,3) - A( 2,3) * A( 3,1)
    Ainv( 1,3) = A( 2,1) * A( 3,2) - A( 2,2) * A( 3,1)
    Ainv( 2,1) = A( 1,2) * A( 3,3) - A( 1,3) * A( 3,2)
    Ainv( 2,2) = A( 1,1) * A( 3,3) - A( 1,3) * A( 3,1)
    Ainv( 2,3) = A( 1,1) * A( 3,2) - A( 1,2) * A( 3,1)
    Ainv( 3,1) = A( 1,2) * A( 2,3) - A( 1,3) * A( 2,2)
    Ainv( 3,2) = A( 1,1) * A( 2,3) - A( 1,3) * A( 2,1)
    Ainv( 3,3) = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Calculate the determinant of A
    detA = A( 1,1) * Ainv( 1,1) - A( 1,2) * Ainv( 1,2) + A( 1,3) * Ainv( 1,3)

  end function calc_determinant_3_by_3

  pure function calc_determinant_5_by_5( A) result( detA)
    ! Determinant of a 5-by-5 matrix
    !
    ! Source: https://caps.gsfc.nasa.gov/simpson/software/m55inv_f90.txt, accessed 2023-02-07

    use iso_fortran_env, only: real64

    ! Input variables:
    real(dp), dimension(5,5), intent(in) :: A

    ! Output variables:
    real(dp) :: detA    ! Determinant of A

    ! Local variables:

    ! Local variables:
    real(real64) :: A11, A12, A13, A14, A15, A21, A22, A23, A24, &
         A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
         A51, A52, A53, A54, A55

    A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
    A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
    A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
    A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
    A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

    detA = real(                                                          &
       A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+       &
       A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-       &
       A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-       &
       A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+       &
       A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+       &
       A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-       &
       A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-       &
       A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-       &
       A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-       &
       A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+       &
       A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+       &
       A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-       &
       A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-       &
       A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+       &
       A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+       &
       A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+       &
       A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+       &
       A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-       &
       A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-       &
       A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+       &
       A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+       &
       A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-       &
       A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-       &
       A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-       &
       A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-       &
       A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+       &
       A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+       &
       A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-       &
       A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-       &
       A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+       &
       A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+       &
       A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+       &
       A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+       &
       A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-       &
       A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-       &
       A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+       &
       A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+       &
       A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-       &
       A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-       &
       A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55, dp)

  end function calc_determinant_5_by_5

  pure function calc_matrix_inverse_2_by_2( A) result( Ainv)
    ! Direct inversion of a 2-by-2 matrix

    ! Input variables:
    real(dp), dimension(2,2), intent(in) :: A       ! The matrix A to be inverted

    ! Output variables:
    real(dp), dimension(2,2) :: Ainv    ! Inverse of A

    ! Local variables:
    real(dp) :: detA

    ! Calculate the determinant of A
    detA = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Safety
    if (abs( detA) < tiny( detA)) then
      error stop 'calc_matrix_inverse_2_by_2 - matrix is singular to working precision!'
    end if

    ! Calculate the inverse of A
    Ainv( 1,1) =  A( 2,2) / detA
    Ainv( 1,2) = -A( 1,2) / detA
    Ainv( 2,1) = -A( 2,1) / detA
    Ainv( 2,2) =  A( 1,1) / detA

  end function calc_matrix_inverse_2_by_2

  pure function calc_matrix_inverse_3_by_3( A) result( Ainv)
    ! Direct inversion of a 3-by-3 matrix
    !
    ! See: https://metric.ma.ic.ac.uk/metric_public/matrices/inverses/inverses2.html

    ! Input variables:
    real(dp), dimension(3,3), intent(in) :: A       ! The matrix A to be inverted

    ! Output variables:
    real(dp), dimension(3,3) :: Ainv    ! Inverse of A

    ! Local variables:
    real(dp) :: detA

    ! Calculate the minors of A
    Ainv( 1,1) = A( 2,2) * A( 3,3) - A( 2,3) * A( 3,2)
    Ainv( 1,2) = A( 2,1) * A( 3,3) - A( 2,3) * A( 3,1)
    Ainv( 1,3) = A( 2,1) * A( 3,2) - A( 2,2) * A( 3,1)
    Ainv( 2,1) = A( 1,2) * A( 3,3) - A( 1,3) * A( 3,2)
    Ainv( 2,2) = A( 1,1) * A( 3,3) - A( 1,3) * A( 3,1)
    Ainv( 2,3) = A( 1,1) * A( 3,2) - A( 1,2) * A( 3,1)
    Ainv( 3,1) = A( 1,2) * A( 2,3) - A( 1,3) * A( 2,2)
    Ainv( 3,2) = A( 1,1) * A( 2,3) - A( 1,3) * A( 2,1)
    Ainv( 3,3) = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Calculate the determinant of A
    detA = A( 1,1) * Ainv( 1,1) - A( 1,2) * Ainv( 1,2) + A( 1,3) * Ainv( 1,3)

    ! Safety
    if (abs( detA) < tiny( detA)) then
      error stop 'calc_matrix_inverse_3_by_3 - matrix is singular to working precision!'
    end if

    ! Change matrix of minors to get the matrix of cofactors
    Ainv( 1,2) = -Ainv( 1,2)
    Ainv( 2,1) = -Ainv( 2,1)
    Ainv( 2,3) = -Ainv( 2,3)
    Ainv( 3,2) = -Ainv( 3,2)

    ! Transpose matrix of cofactors
    Ainv( 1,2) = Ainv( 1,2) + Ainv( 2,1)
    Ainv( 2,1) = Ainv( 1,2) - Ainv( 2,1)
    Ainv( 1,2) = Ainv( 1,2) - Ainv( 2,1)

    Ainv( 1,3) = Ainv( 1,3) + Ainv( 3,1)
    Ainv( 3,1) = Ainv( 1,3) - Ainv( 3,1)
    Ainv( 1,3) = Ainv( 1,3) - Ainv( 3,1)

    Ainv( 2,3) = Ainv( 2,3) + Ainv( 3,2)
    Ainv( 3,2) = Ainv( 2,3) - Ainv( 3,2)
    Ainv( 2,3) = Ainv( 2,3) - Ainv( 3,2)

    ! Divide by det(A)
    Ainv = Ainv / detA

  end function calc_matrix_inverse_3_by_3

  pure function calc_matrix_inverse_5_by_5( A) result( Ainv)
    ! Direct inversion of a 5-by-5 matrix
    !
    ! Source: https://caps.gsfc.nasa.gov/simpson/software/m55inv_f90.txt, accessed 2023-02-07

    use iso_fortran_env, only: real64

    ! Input variables:
    real(dp), dimension(5,5), intent(in)    :: A       ! The matrix A to be inverted

    ! Output variables:
    real(dp), dimension(5,5) :: Ainv    ! Inverse of A

    ! Local variables:
    real(dp) :: detA    ! Determinant of A

    ! Local variables:
    real(real64) :: A11, A12, A13, A14, A15, A21, A22, A23, A24, &
         A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
         A51, A52, A53, A54, A55
    real(real64), dimension(5,5) :: COFACTOR

    A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
    A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
    A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
    A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
    A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

    detA = real(                                                          &
       A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+       &
       A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-       &
       A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-       &
       A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+       &
       A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+       &
       A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-       &
       A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-       &
       A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-       &
       A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-       &
       A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+       &
       A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+       &
       A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-       &
       A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-       &
       A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+       &
       A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+       &
       A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+       &
       A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+       &
       A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-       &
       A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-       &
       A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+       &
       A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+       &
       A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-       &
       A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-       &
       A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-       &
       A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-       &
       A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+       &
       A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+       &
       A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-       &
       A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-       &
       A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+       &
       A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+       &
       A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+       &
       A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+       &
       A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-       &
       A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-       &
       A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+       &
       A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+       &
       A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-       &
       A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-       &
       A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55, dp)

    ! Safety
    if (abs( detA) < tiny( detA)) then
      error stop 'calc_matrix_inverse_3_by_3 - matrix is singular to working precision!'
    end if

    COFACTOR(1,1) = A25*A34*A43*A52-A24*A35*A43*A52-A25*A33*A44*A52+      &
       A23*A35*A44*A52+A24*A33*A45*A52-A23*A34*A45*A52-A25*A34*A42*A53+   &
       A24*A35*A42*A53+A25*A32*A44*A53-A22*A35*A44*A53-A24*A32*A45*A53+   &
       A22*A34*A45*A53+A25*A33*A42*A54-A23*A35*A42*A54-A25*A32*A43*A54+   &
       A22*A35*A43*A54+A23*A32*A45*A54-A22*A33*A45*A54-A24*A33*A42*A55+   &
       A23*A34*A42*A55+A24*A32*A43*A55-A22*A34*A43*A55-A23*A32*A44*A55+   &
       A22*A33*A44*A55

    COFACTOR(2,1) = -A15*A34*A43*A52+A14*A35*A43*A52+A15*A33*A44*A52-     &
       A13*A35*A44*A52-A14*A33*A45*A52+A13*A34*A45*A52+A15*A34*A42*A53-   &
       A14*A35*A42*A53-A15*A32*A44*A53+A12*A35*A44*A53+A14*A32*A45*A53-   &
       A12*A34*A45*A53-A15*A33*A42*A54+A13*A35*A42*A54+A15*A32*A43*A54-   &
       A12*A35*A43*A54-A13*A32*A45*A54+A12*A33*A45*A54+A14*A33*A42*A55-   &
       A13*A34*A42*A55-A14*A32*A43*A55+A12*A34*A43*A55+A13*A32*A44*A55-   &
       A12*A33*A44*A55

    COFACTOR(3,1) = A15*A24*A43*A52-A14*A25*A43*A52-A15*A23*A44*A52+      &
       A13*A25*A44*A52+A14*A23*A45*A52-A13*A24*A45*A52-A15*A24*A42*A53+   &
       A14*A25*A42*A53+A15*A22*A44*A53-A12*A25*A44*A53-A14*A22*A45*A53+   &
       A12*A24*A45*A53+A15*A23*A42*A54-A13*A25*A42*A54-A15*A22*A43*A54+   &
       A12*A25*A43*A54+A13*A22*A45*A54-A12*A23*A45*A54-A14*A23*A42*A55+   &
       A13*A24*A42*A55+A14*A22*A43*A55-A12*A24*A43*A55-A13*A22*A44*A55+   &
       A12*A23*A44*A55

    COFACTOR(4,1) = -A15*A24*A33*A52+A14*A25*A33*A52+A15*A23*A34*A52-     &
       A13*A25*A34*A52-A14*A23*A35*A52+A13*A24*A35*A52+A15*A24*A32*A53-   &
       A14*A25*A32*A53-A15*A22*A34*A53+A12*A25*A34*A53+A14*A22*A35*A53-   &
       A12*A24*A35*A53-A15*A23*A32*A54+A13*A25*A32*A54+A15*A22*A33*A54-   &
       A12*A25*A33*A54-A13*A22*A35*A54+A12*A23*A35*A54+A14*A23*A32*A55-   &
       A13*A24*A32*A55-A14*A22*A33*A55+A12*A24*A33*A55+A13*A22*A34*A55-   &
       A12*A23*A34*A55

    COFACTOR(5,1) = A15*A24*A33*A42-A14*A25*A33*A42-A15*A23*A34*A42+      &
       A13*A25*A34*A42+A14*A23*A35*A42-A13*A24*A35*A42-A15*A24*A32*A43+   &
       A14*A25*A32*A43+A15*A22*A34*A43-A12*A25*A34*A43-A14*A22*A35*A43+   &
       A12*A24*A35*A43+A15*A23*A32*A44-A13*A25*A32*A44-A15*A22*A33*A44+   &
       A12*A25*A33*A44+A13*A22*A35*A44-A12*A23*A35*A44-A14*A23*A32*A45+   &
       A13*A24*A32*A45+A14*A22*A33*A45-A12*A24*A33*A45-A13*A22*A34*A45+   &
       A12*A23*A34*A45

    COFACTOR(1,2) = -A25*A34*A43*A51+A24*A35*A43*A51+A25*A33*A44*A51-     &
       A23*A35*A44*A51-A24*A33*A45*A51+A23*A34*A45*A51+A25*A34*A41*A53-   &
       A24*A35*A41*A53-A25*A31*A44*A53+A21*A35*A44*A53+A24*A31*A45*A53-   &
       A21*A34*A45*A53-A25*A33*A41*A54+A23*A35*A41*A54+A25*A31*A43*A54-   &
       A21*A35*A43*A54-A23*A31*A45*A54+A21*A33*A45*A54+A24*A33*A41*A55-   &
       A23*A34*A41*A55-A24*A31*A43*A55+A21*A34*A43*A55+A23*A31*A44*A55-   &
       A21*A33*A44*A55

    COFACTOR(2,2) = A15*A34*A43*A51-A14*A35*A43*A51-A15*A33*A44*A51+      &
       A13*A35*A44*A51+A14*A33*A45*A51-A13*A34*A45*A51-A15*A34*A41*A53+   &
       A14*A35*A41*A53+A15*A31*A44*A53-A11*A35*A44*A53-A14*A31*A45*A53+   &
       A11*A34*A45*A53+A15*A33*A41*A54-A13*A35*A41*A54-A15*A31*A43*A54+   &
       A11*A35*A43*A54+A13*A31*A45*A54-A11*A33*A45*A54-A14*A33*A41*A55+   &
       A13*A34*A41*A55+A14*A31*A43*A55-A11*A34*A43*A55-A13*A31*A44*A55+   &
       A11*A33*A44*A55

    COFACTOR(3,2) = -A15*A24*A43*A51+A14*A25*A43*A51+A15*A23*A44*A51-     &
       A13*A25*A44*A51-A14*A23*A45*A51+A13*A24*A45*A51+A15*A24*A41*A53-   &
       A14*A25*A41*A53-A15*A21*A44*A53+A11*A25*A44*A53+A14*A21*A45*A53-   &
       A11*A24*A45*A53-A15*A23*A41*A54+A13*A25*A41*A54+A15*A21*A43*A54-   &
       A11*A25*A43*A54-A13*A21*A45*A54+A11*A23*A45*A54+A14*A23*A41*A55-   &
       A13*A24*A41*A55-A14*A21*A43*A55+A11*A24*A43*A55+A13*A21*A44*A55-   &
       A11*A23*A44*A55

    COFACTOR(4,2) = A15*A24*A33*A51-A14*A25*A33*A51-A15*A23*A34*A51+      &
       A13*A25*A34*A51+A14*A23*A35*A51-A13*A24*A35*A51-A15*A24*A31*A53+   &
       A14*A25*A31*A53+A15*A21*A34*A53-A11*A25*A34*A53-A14*A21*A35*A53+   &
       A11*A24*A35*A53+A15*A23*A31*A54-A13*A25*A31*A54-A15*A21*A33*A54+   &
       A11*A25*A33*A54+A13*A21*A35*A54-A11*A23*A35*A54-A14*A23*A31*A55+   &
       A13*A24*A31*A55+A14*A21*A33*A55-A11*A24*A33*A55-A13*A21*A34*A55+   &
       A11*A23*A34*A55

    COFACTOR(5,2) = -A15*A24*A33*A41+A14*A25*A33*A41+A15*A23*A34*A41-     &
       A13*A25*A34*A41-A14*A23*A35*A41+A13*A24*A35*A41+A15*A24*A31*A43-   &
       A14*A25*A31*A43-A15*A21*A34*A43+A11*A25*A34*A43+A14*A21*A35*A43-   &
       A11*A24*A35*A43-A15*A23*A31*A44+A13*A25*A31*A44+A15*A21*A33*A44-   &
       A11*A25*A33*A44-A13*A21*A35*A44+A11*A23*A35*A44+A14*A23*A31*A45-   &
       A13*A24*A31*A45-A14*A21*A33*A45+A11*A24*A33*A45+A13*A21*A34*A45-   &
       A11*A23*A34*A45

    COFACTOR(1,3) = A25*A34*A42*A51-A24*A35*A42*A51-A25*A32*A44*A51+      &
       A22*A35*A44*A51+A24*A32*A45*A51-A22*A34*A45*A51-A25*A34*A41*A52+   &
       A24*A35*A41*A52+A25*A31*A44*A52-A21*A35*A44*A52-A24*A31*A45*A52+   &
       A21*A34*A45*A52+A25*A32*A41*A54-A22*A35*A41*A54-A25*A31*A42*A54+   &
       A21*A35*A42*A54+A22*A31*A45*A54-A21*A32*A45*A54-A24*A32*A41*A55+   &
       A22*A34*A41*A55+A24*A31*A42*A55-A21*A34*A42*A55-A22*A31*A44*A55+   &
       A21*A32*A44*A55

    COFACTOR(2,3) = -A15*A34*A42*A51+A14*A35*A42*A51+A15*A32*A44*A51-     &
       A12*A35*A44*A51-A14*A32*A45*A51+A12*A34*A45*A51+A15*A34*A41*A52-   &
       A14*A35*A41*A52-A15*A31*A44*A52+A11*A35*A44*A52+A14*A31*A45*A52-   &
       A11*A34*A45*A52-A15*A32*A41*A54+A12*A35*A41*A54+A15*A31*A42*A54-   &
       A11*A35*A42*A54-A12*A31*A45*A54+A11*A32*A45*A54+A14*A32*A41*A55-   &
       A12*A34*A41*A55-A14*A31*A42*A55+A11*A34*A42*A55+A12*A31*A44*A55-   &
       A11*A32*A44*A55

    COFACTOR(3,3) = A15*A24*A42*A51-A14*A25*A42*A51-A15*A22*A44*A51+      &
       A12*A25*A44*A51+A14*A22*A45*A51-A12*A24*A45*A51-A15*A24*A41*A52+   &
       A14*A25*A41*A52+A15*A21*A44*A52-A11*A25*A44*A52-A14*A21*A45*A52+   &
       A11*A24*A45*A52+A15*A22*A41*A54-A12*A25*A41*A54-A15*A21*A42*A54+   &
       A11*A25*A42*A54+A12*A21*A45*A54-A11*A22*A45*A54-A14*A22*A41*A55+   &
       A12*A24*A41*A55+A14*A21*A42*A55-A11*A24*A42*A55-A12*A21*A44*A55+   &
       A11*A22*A44*A55

    COFACTOR(4,3) = -A15*A24*A32*A51+A14*A25*A32*A51+A15*A22*A34*A51-     &
       A12*A25*A34*A51-A14*A22*A35*A51+A12*A24*A35*A51+A15*A24*A31*A52-   &
       A14*A25*A31*A52-A15*A21*A34*A52+A11*A25*A34*A52+A14*A21*A35*A52-   &
       A11*A24*A35*A52-A15*A22*A31*A54+A12*A25*A31*A54+A15*A21*A32*A54-   &
       A11*A25*A32*A54-A12*A21*A35*A54+A11*A22*A35*A54+A14*A22*A31*A55-   &
       A12*A24*A31*A55-A14*A21*A32*A55+A11*A24*A32*A55+A12*A21*A34*A55-   &
       A11*A22*A34*A55

    COFACTOR(5,3) = A15*A24*A32*A41-A14*A25*A32*A41-A15*A22*A34*A41+      &
       A12*A25*A34*A41+A14*A22*A35*A41-A12*A24*A35*A41-A15*A24*A31*A42+   &
       A14*A25*A31*A42+A15*A21*A34*A42-A11*A25*A34*A42-A14*A21*A35*A42+   &
       A11*A24*A35*A42+A15*A22*A31*A44-A12*A25*A31*A44-A15*A21*A32*A44+   &
       A11*A25*A32*A44+A12*A21*A35*A44-A11*A22*A35*A44-A14*A22*A31*A45+   &
       A12*A24*A31*A45+A14*A21*A32*A45-A11*A24*A32*A45-A12*A21*A34*A45+   &
       A11*A22*A34*A45

    COFACTOR(1,4) = -A25*A33*A42*A51+A23*A35*A42*A51+A25*A32*A43*A51-     &
       A22*A35*A43*A51-A23*A32*A45*A51+A22*A33*A45*A51+A25*A33*A41*A52-   &
       A23*A35*A41*A52-A25*A31*A43*A52+A21*A35*A43*A52+A23*A31*A45*A52-   &
       A21*A33*A45*A52-A25*A32*A41*A53+A22*A35*A41*A53+A25*A31*A42*A53-   &
       A21*A35*A42*A53-A22*A31*A45*A53+A21*A32*A45*A53+A23*A32*A41*A55-   &
       A22*A33*A41*A55-A23*A31*A42*A55+A21*A33*A42*A55+A22*A31*A43*A55-   &
       A21*A32*A43*A55

    COFACTOR(2,4) = A15*A33*A42*A51-A13*A35*A42*A51-A15*A32*A43*A51+      &
       A12*A35*A43*A51+A13*A32*A45*A51-A12*A33*A45*A51-A15*A33*A41*A52+   &
       A13*A35*A41*A52+A15*A31*A43*A52-A11*A35*A43*A52-A13*A31*A45*A52+   &
       A11*A33*A45*A52+A15*A32*A41*A53-A12*A35*A41*A53-A15*A31*A42*A53+   &
       A11*A35*A42*A53+A12*A31*A45*A53-A11*A32*A45*A53-A13*A32*A41*A55+   &
       A12*A33*A41*A55+A13*A31*A42*A55-A11*A33*A42*A55-A12*A31*A43*A55+   &
       A11*A32*A43*A55

    COFACTOR(3,4) = -A15*A23*A42*A51+A13*A25*A42*A51+A15*A22*A43*A51-     &
       A12*A25*A43*A51-A13*A22*A45*A51+A12*A23*A45*A51+A15*A23*A41*A52-   &
       A13*A25*A41*A52-A15*A21*A43*A52+A11*A25*A43*A52+A13*A21*A45*A52-   &
       A11*A23*A45*A52-A15*A22*A41*A53+A12*A25*A41*A53+A15*A21*A42*A53-   &
       A11*A25*A42*A53-A12*A21*A45*A53+A11*A22*A45*A53+A13*A22*A41*A55-   &
       A12*A23*A41*A55-A13*A21*A42*A55+A11*A23*A42*A55+A12*A21*A43*A55-   &
       A11*A22*A43*A55

    COFACTOR(4,4) = A15*A23*A32*A51-A13*A25*A32*A51-A15*A22*A33*A51+      &
       A12*A25*A33*A51+A13*A22*A35*A51-A12*A23*A35*A51-A15*A23*A31*A52+   &
       A13*A25*A31*A52+A15*A21*A33*A52-A11*A25*A33*A52-A13*A21*A35*A52+   &
       A11*A23*A35*A52+A15*A22*A31*A53-A12*A25*A31*A53-A15*A21*A32*A53+   &
       A11*A25*A32*A53+A12*A21*A35*A53-A11*A22*A35*A53-A13*A22*A31*A55+   &
       A12*A23*A31*A55+A13*A21*A32*A55-A11*A23*A32*A55-A12*A21*A33*A55+   &
       A11*A22*A33*A55

    COFACTOR(5,4) = -A15*A23*A32*A41+A13*A25*A32*A41+A15*A22*A33*A41-     &
       A12*A25*A33*A41-A13*A22*A35*A41+A12*A23*A35*A41+A15*A23*A31*A42-   &
       A13*A25*A31*A42-A15*A21*A33*A42+A11*A25*A33*A42+A13*A21*A35*A42-   &
       A11*A23*A35*A42-A15*A22*A31*A43+A12*A25*A31*A43+A15*A21*A32*A43-   &
       A11*A25*A32*A43-A12*A21*A35*A43+A11*A22*A35*A43+A13*A22*A31*A45-   &
       A12*A23*A31*A45-A13*A21*A32*A45+A11*A23*A32*A45+A12*A21*A33*A45-   &
       A11*A22*A33*A45

    COFACTOR(1,5) = A24*A33*A42*A51-A23*A34*A42*A51-A24*A32*A43*A51+      &
       A22*A34*A43*A51+A23*A32*A44*A51-A22*A33*A44*A51-A24*A33*A41*A52+   &
       A23*A34*A41*A52+A24*A31*A43*A52-A21*A34*A43*A52-A23*A31*A44*A52+   &
       A21*A33*A44*A52+A24*A32*A41*A53-A22*A34*A41*A53-A24*A31*A42*A53+   &
       A21*A34*A42*A53+A22*A31*A44*A53-A21*A32*A44*A53-A23*A32*A41*A54+   &
       A22*A33*A41*A54+A23*A31*A42*A54-A21*A33*A42*A54-A22*A31*A43*A54+   &
       A21*A32*A43*A54

    COFACTOR(2,5) = -A14*A33*A42*A51+A13*A34*A42*A51+A14*A32*A43*A51-     &
       A12*A34*A43*A51-A13*A32*A44*A51+A12*A33*A44*A51+A14*A33*A41*A52-   &
       A13*A34*A41*A52-A14*A31*A43*A52+A11*A34*A43*A52+A13*A31*A44*A52-   &
       A11*A33*A44*A52-A14*A32*A41*A53+A12*A34*A41*A53+A14*A31*A42*A53-   &
       A11*A34*A42*A53-A12*A31*A44*A53+A11*A32*A44*A53+A13*A32*A41*A54-   &
       A12*A33*A41*A54-A13*A31*A42*A54+A11*A33*A42*A54+A12*A31*A43*A54-   &
       A11*A32*A43*A54

    COFACTOR(3,5) = A14*A23*A42*A51-A13*A24*A42*A51-A14*A22*A43*A51+      &
       A12*A24*A43*A51+A13*A22*A44*A51-A12*A23*A44*A51-A14*A23*A41*A52+   &
       A13*A24*A41*A52+A14*A21*A43*A52-A11*A24*A43*A52-A13*A21*A44*A52+   &
       A11*A23*A44*A52+A14*A22*A41*A53-A12*A24*A41*A53-A14*A21*A42*A53+   &
       A11*A24*A42*A53+A12*A21*A44*A53-A11*A22*A44*A53-A13*A22*A41*A54+   &
       A12*A23*A41*A54+A13*A21*A42*A54-A11*A23*A42*A54-A12*A21*A43*A54+   &
       A11*A22*A43*A54

    COFACTOR(4,5) = -A14*A23*A32*A51+A13*A24*A32*A51+A14*A22*A33*A51-     &
       A12*A24*A33*A51-A13*A22*A34*A51+A12*A23*A34*A51+A14*A23*A31*A52-   &
       A13*A24*A31*A52-A14*A21*A33*A52+A11*A24*A33*A52+A13*A21*A34*A52-   &
       A11*A23*A34*A52-A14*A22*A31*A53+A12*A24*A31*A53+A14*A21*A32*A53-   &
       A11*A24*A32*A53-A12*A21*A34*A53+A11*A22*A34*A53+A13*A22*A31*A54-   &
       A12*A23*A31*A54-A13*A21*A32*A54+A11*A23*A32*A54+A12*A21*A33*A54-   &
       A11*A22*A33*A54

    COFACTOR(5,5) = A14*A23*A32*A41-A13*A24*A32*A41-A14*A22*A33*A41+      &
       A12*A24*A33*A41+A13*A22*A34*A41-A12*A23*A34*A41-A14*A23*A31*A42+   &
       A13*A24*A31*A42+A14*A21*A33*A42-A11*A24*A33*A42-A13*A21*A34*A42+   &
       A11*A23*A34*A42+A14*A22*A31*A43-A12*A24*A31*A43-A14*A21*A32*A43+   &
       A11*A24*A32*A43+A12*A21*A34*A43-A11*A22*A34*A43-A13*A22*A31*A44+   &
       A12*A23*A31*A44+A13*A21*A32*A44-A11*A23*A32*A44-A12*A21*A33*A44+   &
       A11*A22*A33*A44

    AINV = real( transpose( COFACTOR),dp) / detA

  end function calc_matrix_inverse_5_by_5

  pure function solve_Axb_2_by_2( A,b) result( x)
    ! Direct solution of the 2-by-2 matrix equation Ax=b

    ! Input variables:
    real(dp), dimension(2,2), intent(in)    :: A       ! The matrix A
    real(dp), dimension(2  ), intent(in)    :: b       ! The right-hand side b

    ! Output variables:
    real(dp), dimension(2  ) :: x       ! The solution x

    ! Local variables:
    real(dp), dimension(2,2) :: Ainv

    ! Calculate the inverse of A
    Ainv = calc_matrix_inverse_2_by_2( A)

    ! Calculate x
    x( 1) = Ainv( 1,1) * b( 1) + Ainv( 1,2) * b( 2)
    x( 2) = Ainv( 2,1) * b( 1) + Ainv( 2,2) * b( 2)

  end function solve_Axb_2_by_2

! == Sorting

  subroutine quick_n_dirty_sort( f, ii)
    ! Inefficient but simple sorting algorithm
    ! Sort f ascending and return list of new indices

    ! In/output variables:
    real(dp), dimension(:), intent(inout) :: f
    integer,  dimension(:), intent(inout) :: ii

    ! Local variables:
    character(len=256), parameter :: routine_name = 'quick_n_dirty_sort'
    integer                       :: n,i,j

    ! Add routine to path
    call init_routine( routine_name)

    ! Number of elements
    n = size( f,1)

    ! Initialise current indices
    do i = 1, n
      ii( i) = i
    end do

    ! Sort
    do i = 1, n-1
      do j = i+1, n
        if (f( i) > f( j)) then

          ! Now: f( i) = a, f( j) = b
          f(  j) = f(  j) + f(  i)
          ii( j) = ii( j) + ii( i)
          ! Now: f( i) = a, f( j) = a+b
          f(  i) = f(  j) - f(  i)
          ii( i) = ii( j) - ii( i)
          ! Now: f( i) = b, f( j) = a+b
          f(  j) = f(  j) - f(  i)
          ii( j) = ii( j) - ii( i)
          ! Now: f( i) = b, f( j) = a

        end if
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine quick_n_dirty_sort

! == Some basic geometry

  pure function cross2( a,b) result( z)
    ! Vector product z between 2-dimensional vectors a and b

    ! Input variables:
    real(dp), dimension(2), intent(in)    :: a, b

    ! Output variables:
    real(dp) :: z

    z = (a( 1) * b( 2)) - (a( 2) * b( 1))

  end function cross2

  pure subroutine segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
    ! Find out if the line segments [pq] and [rs] intersect. If so, return
    ! the coordinates of the point of intersection

    ! Input variables:
    real(dp), dimension(2), intent(in)    :: p, q, r, s
    real(dp)              , intent(in)    :: tol_dist

    ! Output variables
    real(dp), dimension(2), intent(out)   :: llis
    logical               , intent(out)   :: do_cross

    ! Local variables:
    real(dp), dimension(2,2) :: A
    real(dp), dimension(2)   :: x, b

    ! If pq and rs are colinear, define them as not intersecting
    if ( abs( cross2( [q(1)-p(1), q(2)-p(2)], [s(1)-r(1), s(2)-r(2)] )) < tol_dist) then
      llis = [0._dp, 0._dp]
      do_cross = .false.
      return
    end if

    A( 1,:) = [ (p( 1) - q( 1)), (r( 1) - s( 1))]
    A( 2,:) = [ (p( 2) - q( 2)), (r( 2) - s( 2))]
    b = [ (r( 1) - q( 1)), (r( 2) - q( 2))]

    ! Solve for x
    x = solve_Axb_2_by_2( A,b)

    llis = [q( 1) + x( 1) * (p( 1) - q( 1)), q( 2) + x( 1) * (p( 2) - q( 2))]

    if (x( 1) > 0._dp .and. x( 1) < 1._dp .and. x( 2) > 0._dp .and. x( 2) < 1._dp) then
      do_cross = .true.
    else
      do_cross = .false.
    end if

  end subroutine segment_intersection

  pure function is_in_polygons( poly_mult, p) result( isso)
    ! Check if p is inside any of the polygons in poly_mult

    ! Input variables:
    real(dp), dimension(:,:) , intent(in)    :: poly_mult
    real(dp), dimension(2)   , intent(in)    :: p

    ! Output variables:
    logical :: isso

    ! Local variables:
    integer                                 :: n1,n2,nn
    real(dp), dimension(:,:  ), allocatable :: poly
    logical                                 :: isso_single

    isso = .false.

    n1 = 1
    n2 = 0

    do while (n2 < size( poly_mult,1))

      ! Copy a single polygon from poly_mult
      nn = nint( poly_mult( n1,1))
      n2 = n1 + nn
      allocate( poly( nn,2))
      poly = poly_mult( n1+1:n2,:)
      n1 = n2+1

      ! Check if p is inside this single polygon
      isso_single = is_in_polygon( poly, p)

      ! Clean up after yourself
      deallocate( poly)

      if (isso_single) then
        isso = .true.
        exit
      end if

    end do ! do while (n2 < size( poly_mult,1))

  end function is_in_polygons

  pure function is_in_polygon( Vpoly, p) result( isso)
    ! Use the ray-casting algorithm to check if the point p = [px,py]
    ! lies inside the polygon spanned by poly = [x1,y1; x2,y2; ...; xn,yn]

    ! Input variables:
    real(dp), dimension(:,:) , intent(in)    :: Vpoly
    real(dp), dimension(2)   , intent(in)    :: p

    ! Output variables:
    logical :: isso

    ! Local variables:
    real(dp), dimension(2) :: q,r,s,llis
    real(dp)               :: xmin,xmax,ymin,ymax
    integer                :: n_intersects
    integer                :: vi,vj,n_vertices
    logical                :: do_cross
    real(dp), parameter    :: tol_dist = 1E-5_dp

    isso = .false.

    xmin = minval( Vpoly(:,1))
    xmax = maxval( Vpoly(:,1))
    ymin = minval( Vpoly(:,2))
    ymax = maxval( Vpoly(:,2))

    ! Quick test
    if (p(1) < xmin .or. p(1) > xmax .or. &
        p(2) < ymin .or. p(2) > ymax) then
      isso = .false.
      return
    end if

    ! Define the endpoint of the east-pointing ray
    q = [xmax + (xmax - xmin) / 10._dp, p(2)]

    ! Determine how often the ray intersects the polygon

    n_vertices   = size( Vpoly,1)
    n_intersects = 0

    do vi = 1, n_vertices

      ! Find vertices spanning a polygon line section
      if (vi < n_vertices) then
        vj = vi + 1
      else
        vj = 1
      end if

      ! Define the line section
      r = Vpoly( vi,:)
      s = Vpoly( vj,:)

      ! Determine if the ray intersects the line section
      if ((r(2) < p(2) .and. s(2) < p(2)) .or. (r(2) > p(2) .and. s(2) > p(2))) then
        do_cross = .false.
      else
        call segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
      end if

      if (do_cross) n_intersects = n_intersects + 1

    end do ! do vi = 1, n_vertices

    ! If the number of intersections is odd, p lies inside the polygon
    if (MOD( n_intersects,2) == 1) then
      isso = .true.
    else
      isso = .false.
    end if

  end function is_in_polygon

  pure function lies_on_line_segment( pa, pb, pc, tol_dist) result(isso)
    ! Test if the point pc lies within tol_dist of the line pa-pb

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: pa, pb, pc
    real(dp),               intent(in) :: tol_dist
    logical                            :: isso

    ! Local variables:
    real(dp), dimension(2) :: d, e, d_norm, e_par, e_ort

    d = pb - pa
    e = pc - pa

    d_norm = d / norm2( d)

    e_par = (e( 1) * d_norm( 1) + e( 2) * d_norm( 2)) * d_norm
    e_ort = e - e_par

    isso = .true.
    if (norm2( e_ort) > tol_dist) then
      isso = .false.
      return
    end if

    if ((e( 1) * d( 1) + e( 2)* d ( 2)) > 0._dp) then
      if (norm2( e_par) > (norm2( d))) then
        isso = .false.
        return
      end if
    else
      if (norm2( e_par) > 0._dp) then
        isso = .false.
        return
      end if
    end if

  end function lies_on_line_segment

  pure subroutine line_from_points( p, q, la, lb, lc)
    ! Find a,b,c such that the line ax + by = c passes through p and q

    real(dp), dimension(2), intent(in   ) :: p, q
    real(dp),               intent(  out) :: la, lb, lc

    la = q( 2) - p( 2)
    lb = p( 1) - q( 1)
    lc = la * (p( 1))+ lb * (p( 2))

  end subroutine line_from_points

  pure subroutine perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)
    ! Find a,b,c such that the line ax + by = c describes the perpendicular
    ! bisector to the line [pq]

    real(dp), dimension(2), intent(in   ) :: p, q
    real(dp),               intent(in   ) :: la1, lb1
    real(dp),               intent(  out) :: la2, lb2, lc2
    real(dp)                              :: temp
    real(dp), dimension(2)                :: m

    m = (p+q)/2
    lc2 = -lb1*m(1) + la1*m(2)

    temp = la1
    la2 = -lb1
    lb2 = temp

  end subroutine perpendicular_bisector_from_line

  pure subroutine line_line_intersection( la1, lb1, lc1, la2, lb2, lc2, llis)
    ! Find the intersection llis of the lines la1*x+lb1*y=lc1 and la2*x+lb2*y=lc2

    real(dp),               intent(in   ) :: la1, lb1, lc1, la2, lb2, lc2
    real(dp), dimension(2), intent(  out) :: llis
    real(dp)                              :: d

    d = la1*lb2 - la2*lb1
    if (d == 0) then
      ! The lines are parallel.
      llis = [1E30, 1E30]
    else
      llis = [(lb2*lc1 - lb1*lc2), (la1*lc2 - la2*lc1)]/d
    end if

  end subroutine line_line_intersection

  pure function circumcenter( p, q, r) result( cc)
    ! Find the circumcenter cc of the triangle pqr

    ! Some basic vector operations
    ! Find the circumcenter cc of the triangle pqr
    ! If pqr are colinear, returns [1e30,1e30]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp), dimension(2)             :: cc
    real(dp)                           :: la1,lb1,lc1,le1,lf1,lg1
    real(dp)                           :: la2,lb2,lc2,le2,lf2,lg2

    cc = [0._dp, 0._dp]

    ! Line PQ is represented as ax + by = c, Line QR is represented as ex + fy = g
    call line_from_points( p, q, la1, lb1, lc1)
    call line_from_points( q, r, le1, lf1, lg1)

    ! Converting lines PQ and QR to perpendicular
    ! bisectors. After this, L = ax + by = c
    ! M = ex + fy = g
    call perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)
    call perpendicular_bisector_from_line( q, r, le1, lf1, le2, lf2, lg2)

    ! The point of intersection of L and M gives
    ! the circumcenter
    call line_line_intersection( la2, lb2, lc2, le2, lf2, lg2, cc)

  end function circumcenter

  pure function geometric_center( p, q, r) result( gc)
    ! Calculate the geometric centre of triangle [pqr]

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp), dimension(2)             :: gc

    gc = (p + q + r) / 3._dp

  end function geometric_center

  pure function triangle_area( p, q, r) result( A)
    ! Find the area of the triangle [p,q,r]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: A

    A = abs( cross2( [q(1)-p(1), q(2)-p(2)], [r(1)-p(1), r(2)-p(2)] )) / 2._dp

  end function triangle_area

  pure function is_in_triangle( pa, pb, pc, p) result(isso)
    ! Check if the point p lies inside the triangle abc, or within distance tol of its edges

    real(dp), dimension(2), intent(in) :: pa, pb, pc, p
    logical                            :: isso
    real(dp)                           :: as_x, as_y, s1, s2, s3
    real(dp), parameter                :: tol = 1E-9_dp

    as_x = p( 1) - pa( 1)
    as_y = p( 2) - pa( 2)

    s1 = ((pb( 1) - pa( 1)) * as_y - (pb( 2) - pa( 2)) * as_x)
    s2 = ((pc( 1) - pa( 1)) * as_y - (pc( 2) - pa( 2)) * as_x)
    s3 = ((pc( 1) - pb( 1)) * (p( 2) - pb( 2)) - (pc( 2) - pb( 2)) * (p( 1) - pb( 1)))

    isso = .false.

    if (s1 > -tol .and. s2 < tol .and. s3 > -tol) then
      isso = .true.
      return
    end if

  end function is_in_triangle

  pure function longest_triangle_leg( p, q, r) result( d)
    ! Find the longest leg of the triangle [p,q,r]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: d
    real(dp)                           :: d_pq, d_qr, d_rp

    d_pq = norm2( p-q)
    d_qr = norm2( q-r)
    d_rp = norm2( r-p)
    d = max( max( d_pq, d_qr), d_rp)

  end function longest_triangle_leg

  pure function smallest_triangle_angle( p, q, r) result( alpha)
    ! Find the smallest internal angle of the triangle [p,q,r]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: alpha
    real(dp), dimension(2)             :: pq, qr, rp
    real(dp)                           :: ap, aq, ar

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Internal angles
    ap = acos(-(rp( 1) * pq( 1) + rp( 2) * pq( 2)) / (norm2( rp) * norm2( pq)))
    aq = acos(-(pq( 1) * qr( 1) + pq( 2) * qr( 2)) / (norm2( pq) * norm2( qr)))
    ar = acos(-(rp( 1) * qr( 1) + rp( 2) * qr( 2)) / (norm2( rp) * norm2( qr)))

    ! Smallest internal angle
    alpha = min( min( ap, aq), ar)

  end function smallest_triangle_angle

  !> Calculate the largest internal angle of the triangle [p,q,r]
  pure function largest_triangle_angle( p, q, r) result( alpha)

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: alpha

    ! Local variables:
    real(dp), dimension(2) :: pq, qr, rp
    real(dp)               :: ap, aq, ar

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Internal angles
    ap = acos(-(rp( 1) * pq( 1) + rp( 2) * pq( 2)) / (norm2( rp) * norm2( pq)))
    aq = acos(-(pq( 1) * qr( 1) + pq( 2) * qr( 2)) / (norm2( pq) * norm2( qr)))
    ar = acos(-(rp( 1) * qr( 1) + rp( 2) * qr( 2)) / (norm2( rp) * norm2( qr)))

    ! Largest internal angle
    alpha = max( max( ap, aq), ar)

  end function largest_triangle_angle

  !> Calculate the equiangular skewness of a triangle
  pure function equiangular_skewness( p, q, r) result( skewness)
    ! See: https://en.wikipedia.org/wiki/Types_of_mesh#Equiangular_skew

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: skewness

    ! Local variables:
    real(dp) :: theta_e, theta_max, theta_min

    theta_e = 60._dp * pi / 180._dp

    theta_max = largest_triangle_angle ( p, q, r)
    theta_min = smallest_triangle_angle( p, q, r)

    skewness = max( (theta_max - theta_e  ) / (pi / 2._dp - theta_e), &
                    (theta_e   - theta_min) /               theta_e)

    end function equiangular_skewness

  subroutine crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, tol_dist, pp, qq, is_valid_line)
    ! Crop the line [pq] so that it lies within the specified domain;
    ! if [pq] doesn't pass through the domain at all, return is_valid_line = .false.

    ! In/output variables
    real(dp), dimension(2), intent(in   ) :: p, q
    real(dp),               intent(in   ) :: xmin, xmax, ymin, ymax, tol_dist
    real(dp), dimension(2), intent(  out) :: pp, qq
    logical,                intent(  out) :: is_valid_line

    ! Local variables:
    real(dp), dimension(2) :: sw, se, nw, ne
    logical                :: p_inside
    logical                :: p_on_border
    logical                :: p_outside
    logical                :: q_inside
    logical                :: q_on_border
    logical                :: q_outside
    logical                :: do_cross_w
    logical                :: do_cross_e
    logical                :: do_cross_s
    logical                :: do_cross_n
    integer                :: n_cross
    real(dp), dimension(2) :: llis_w, llis_e, llis_s, llis_n, llis1, llis2

    sw = [xmin,ymin]
    se = [xmax,ymin]
    nw = [xmin,ymax]
    ne = [xmax,ymax]

    ! Determine where p and q are relative to the domain

    p_inside    = .false.
    p_on_border = .false.
    p_outside   = .false.

    IF     (p( 1) > xmin .and. p( 1) < xmax .and. p( 2) > ymin .and. p( 2) < ymax) then
      p_inside    = .true.
    elseif (p (1) < xmin .or.  p( 1) > xmax .or.  p( 2) < ymin .or.  p( 2) > ymax) then
      p_outside   = .true.
    else
      p_on_border = .true.
    end if

    q_inside    = .false.
    q_on_border = .false.
    q_outside   = .false.

    IF     (q( 1) > xmin .and. q( 1) < xmax .and. q( 2) > ymin .and. q( 2) < ymax) then
      q_inside    = .true.
    elseif (q (1) < xmin .or.  q( 1) > xmax .or.  q( 2) < ymin .or.  q( 2) > ymax) then
      q_outside   = .true.
    else
      q_on_border = .true.
    end if

    ! If both of them lie inside the domain, the solution is trivial
    if (p_inside .and. q_inside) then
      pp = p
      qq = q
      is_valid_line = .true.
      return
    end if

    ! If both of them lie on the border, the solution is trivial
    if (p_on_border .and. q_on_border) then
      pp = p
      qq = q
      is_valid_line = .true.
      return
    end if

    ! If one of them lies inside the domain and the other on the border, the solution is trivial
    if ((p_inside .and. q_on_border) .or. (p_on_border .and. q_inside)) then
      pp = p
      qq = q
      is_valid_line = .true.
      return
    end if

    ! If one of them lies inside and the other outside, there must be a single border crossing
    if (p_inside .and. q_outside) then
      ! p lies inside the domain, q lies outside

      ! Possible pq passes through a corner of the domain?
      IF     (lies_on_line_segment( p, q, nw, tol_dist)) then
        ! pq passes through the northwest corner
        pp = p
        qq = nw
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, ne, tol_dist)) then
        ! pq passes through the northeast corner
        pp = p
        qq = ne
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, sw, tol_dist)) then
        ! pq passes through the southwest corner
        pp = p
        qq = sw
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, se, tol_dist)) then
        ! pq passes through the southeast corner
        pp = p
        qq = se
        is_valid_line = .true.
        return
      end if

      ! pq must pass through one of the four borders; determine which one
      call segment_intersection( p, q, sw, nw, llis_w, do_cross_w, tol_dist)
      call segment_intersection( p, q, se, ne, llis_e, do_cross_e, tol_dist)
      call segment_intersection( p, q, sw, se, llis_s, do_cross_s, tol_dist)
      call segment_intersection( p, q, nw, ne, llis_n, do_cross_n, tol_dist)

      IF     (do_cross_w) then
        ! pq crosses the western border
        pp = p
        qq = llis_w
        is_valid_line = .true.
        return
      elseif (do_cross_e) then
        ! pq crosses the eastern border
        pp = p
        qq = llis_e
        is_valid_line = .true.
        return
      elseif (do_cross_s) then
        ! pq crosses the southern border
        pp = p
        qq = llis_s
        is_valid_line = .true.
        return
      elseif (do_cross_n) then
        ! pq crosses the northern border
        pp = p
        qq = llis_n
        is_valid_line = .true.
        return
      end if

      ! This point should not be reachable!
      call crash('crop_line_to_domain - p lies inside, q lies outside, couldnt find exit point of pq!')

    elseif (p_outside .and. q_inside) then
      ! p lies outside the domain, q lies inside

      ! Possible pq passes through a corner of the domain?
      IF     (lies_on_line_segment( p, q, nw, tol_dist)) then
        ! pq passes through the northwest corner
        pp = nw
        qq = q
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, ne, tol_dist)) then
        ! pq passes through the northeast corner
        pp = ne
        qq = q
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, sw, tol_dist)) then
        ! pq passes through the southwest corner
        pp = sw
        qq = q
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, se, tol_dist)) then
        ! pq passes through the southeast corner
        pp = se
        qq = q
        is_valid_line = .true.
        return
      end if

      ! pq must pass through one of the four borders; determine which one
      call segment_intersection( p, q, sw, nw, llis_w, do_cross_w, tol_dist)
      call segment_intersection( p, q, se, ne, llis_e, do_cross_e, tol_dist)
      call segment_intersection( p, q, sw, se, llis_s, do_cross_s, tol_dist)
      call segment_intersection( p, q, nw, ne, llis_n, do_cross_n, tol_dist)

      IF     (do_cross_w) then
        ! pq crosses the western border
        pp = llis_w
        qq = q
        is_valid_line = .true.
        return
      elseif (do_cross_e) then
        ! pq crosses the eastern border
        pp = llis_e
        qq = q
        is_valid_line = .true.
        return
      elseif (do_cross_s) then
        ! pq crosses the southern border
        pp = llis_s
        qq = q
        is_valid_line = .true.
        return
      elseif (do_cross_n) then
        ! pq crosses the northern border
        pp = llis_n
        qq = q
        is_valid_line = .true.
        return
      end if

      ! This point should not be reachable!
      call crash('crop_line_to_domain - p lies outside, q lies inside, couldnt find exit point of pq!')

    end if ! if (p_inside .and. q_outside) then

    ! If both of them lie outside the domain, there might still be a section passing through it
    if (p_outside .and. q_outside) then

      ! pq must pass through either none, or two of the four borders; determine which
      call segment_intersection( p, q, sw, nw, llis_w, do_cross_w, tol_dist)
      call segment_intersection( p, q, se, ne, llis_e, do_cross_e, tol_dist)
      call segment_intersection( p, q, sw, se, llis_s, do_cross_s, tol_dist)
      call segment_intersection( p, q, nw, ne, llis_n, do_cross_n, tol_dist)

      n_cross = 0

      if (do_cross_w) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_w
        else
          llis2 = llis_w
        end if
      end if

      if (do_cross_e) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_e
        else
          llis2 = llis_e
        end if
      end if

      if (do_cross_s) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_s
        else
          llis2 = llis_s
        end if
      end if

      if (do_cross_n) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_n
        else
          llis2 = llis_n
        end if
      end if

      IF     (n_cross == 0) then
        ! pq does not pass through the domain at all
        pp = 0._dp
        qq = 0._dp
        is_valid_line = .false.
        return
      elseif (n_cross == 2) then
        ! pq passes through the domain; crop it

        if (norm2( llis1 - p) < norm2( llis2 - p)) then
          ! the cropped line runs from llis1 to llis2
          pp = llis1
          qq = llis2
          is_valid_line = .true.
          return
        else
          ! the cropped lines runs from llis2 to llis1
          pp = llis2
          qq = llis1
          is_valid_line = .true.
          return
        end if

      else
        ! This should not be possible
        call crash('pq crosses the domain border {int_01} times!', int_01 = n_cross)
      end if

    end if ! if (p_outside .and. q_outside) then

    ! If one of them lies on the border and another outside, it is possible
    ! that the line still passes through the domain, but we neglect that possibility for now...
    if ((p_on_border .and. q_outside) .or. (p_outside .and. q_on_border)) then
      pp = 0._dp
      qq = 0._dp
      is_valid_line = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('crop_line_to_domain - reached the unreachable end!')

  end subroutine crop_line_to_domain

  pure function encroaches_upon( pa, pb, pc, tol_dist) result( isso)
    ! Check if c encroaches upon segment ab

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: pa, pb, pc
    real(dp),               intent(in) :: tol_dist
    logical                            :: isso

    isso = norm2( pc - (pa + pb) / 2._dp) < norm2( pa - pb) / 2._dp + tol_dist

  end function encroaches_upon

  pure function linint_points( x1, x2, f1, f2, f0) result( x0)
    ! Given a function f( x) and points x1, x2 such that f( x1) = f1, f( x2) = f2,
    ! interpolate f linearly to find the point x0 such that f( x0) = f0

    real(dp), intent(in) :: x1, x2, f1, f2, f0
    real(dp)             :: x0
    real(dp)             :: lambda

    ! Safety - if f1 == f2, then f = f0 = f1 everywhere
    if (abs( 1._dp - f1/f2) < 1E-9_dp) then
      x0 = (x1 + x2) / 2._dp
      return
    end if

    lambda = (f2 - f1) / (x2 - x1)
    x0 = x1 + (f0 - f1) / lambda

  end function linint_points

! == Basic array operations

  subroutine permute_2D_int( d, map)
    ! Permute a 2-D array

    ! In/output variables:
    integer,  dimension(:,:  ), allocatable, intent(inout) :: d
    integer,  dimension(2),                  intent(in)    :: map

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'permute_2D_int'
    integer                                 :: i,j,n1,n2
    integer,  dimension(:,:  ), allocatable :: d_temp

    ! Add routine to path
    call init_routine( routine_name)

    if (map( 1) == 1 .and. map( 2) == 2) then
      ! Trivial
      return
    elseif (map( 1) == 2 .and. map( 2) == 1) then
      ! 2-D transpose, as expected
    else
      call crash('invalid permutation!')
    end if

    n1 = size( d,1)
    n2 = size( d,2)

    ! allocate temporary memory
    allocate( d_temp( n1, n2))

    ! Copy data to temporary memory
    d_temp = d

    ! deallocate memory
    deallocate( d)

    ! Reallocate transposed memory
    allocate( d( n2, n1))

    ! Copy and transpose data from temporary memory
    do i = 1, n1
    do j = 1, n2
      d( j,i) = d_temp( i,j)
    end do
    end do

    ! deallocate temporary memory
    deallocate( d_temp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine permute_2D_int

  subroutine permute_2D_dp( d, map)
    ! Permute a 2-D array

    ! In/output variables:
    real(dp), dimension(:,:), allocatable, intent(inout) :: d
    integer,  dimension(2),                intent(in)    :: map

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'permute_2D_dp'
    integer                               :: i,j,n1,n2
    real(dp), dimension(:,:), allocatable :: d_temp

    ! Add routine to path
    call init_routine( routine_name)

    if (map( 1) == 1 .and. map( 2) == 2) then
      ! Trivial
      return
    elseif (map( 1) == 2 .and. map( 2) == 1) then
      ! 2-D transpose, as expected
    else
      call crash('invalid permutation!')
    end if

    n1 = size( d,1)
    n2 = size( d,2)

    ! allocate temporary memory
    allocate( d_temp( n1, n2))

    ! Copy data to temporary memory
    d_temp = d

    ! deallocate memory
    deallocate( d)

    ! Reallocate transposed memory
    allocate( d( n2, n1))

    ! Copy and transpose data from temporary memory
    do i = 1, n1
    do j = 1, n2
      d( j,i) = d_temp( i,j)
    end do
    end do

    ! deallocate temporary memory
    deallocate( d_temp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine permute_2D_dp

  subroutine permute_3D_int( d, map)
    ! Permute a 3-D array

    ! In/output variables:
    integer,  dimension(:,:,:), allocatable, intent(inout) :: d
    integer,  dimension(3),                  intent(in)    :: map

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'permute_3D_int'
    integer                                 :: i,j,k,n1,n2,n3
    integer,  dimension(:,:,:), allocatable :: d_temp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! allocate temporary memory
    allocate( d_temp( n1, n2, n3))

    ! Copy data to temporary memory
    d_temp = d

    ! deallocate memory
    deallocate( d)

    ! Different permutation options
    if (map( 1) == 1 .and. map( 2) == 2 .and. map( 3) == 3) then
      ! [i,j,k] -> [i,j,k] (trivial...)

      ! Reallocate permuted memory
      allocate( d( n1, n2, n3))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( i,j,k) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 1 .and. map( 2) == 3 .and. map( 3) == 2) then
      ! [i,j,k] -> [i,k,j]

      ! Reallocate permuted memory
      allocate( d( n1, n3, n2))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( i,k,j) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 2 .and. map( 2) == 1 .and. map( 3) == 3) then
      ! [i,j,k] -> [j,i,k]

      ! Reallocate permuted memory
      allocate( d( n2, n1, n3))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( j,i,k) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 2 .and. map( 2) == 3 .and. map( 3) == 1) then
      ! [i,j,k] -> [j,k,i]

      ! Reallocate permuted memory
      allocate( d( n2, n3, n1))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( j,k,i) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 3 .and. map( 2) == 1 .and. map( 3) == 2) then
      ! [i,j,k] -> [k,i,j]

      ! Reallocate permuted memory
      allocate( d( n3, n1, n2))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( k,i,j) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 3 .and. map( 2) == 2 .and. map( 3) == 1) then
      ! [i,j,k] -> [k,j,i]

      ! Reallocate permuted memory
      allocate( d( n3, n2, n1))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( k,j,i) = d_temp( i,j,k)
      end do
      end do
      end do

    else
      call crash('invalid permutation!')
    end if

    ! deallocate temporary memory
    deallocate( d_temp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine permute_3D_int

  subroutine permute_3D_dp( d, map)
    ! Permute a 3-D array

    ! In/output variables:
    real(dp), dimension(:,:,:), allocatable, intent(inout) :: d
    integer,  dimension(3),                  intent(in)    :: map

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'permute_3D_dp'
    integer                                 :: i,j,k,n1,n2,n3
    real(dp), dimension(:,:,:), allocatable :: d_temp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! allocate temporary memory
    allocate( d_temp( n1, n2, n3))

    ! Copy data to temporary memory
    d_temp = d

    ! deallocate memory
    deallocate( d)

    ! Different permutation options
    if (map( 1) == 1 .and. map( 2) == 2 .and. map( 3) == 3) then
      ! [i,j,k] -> [i,j,k] (trivial...)

      ! Reallocate permuted memory
      allocate( d( n1, n2, n3))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( i,j,k) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 1 .and. map( 2) == 3 .and. map( 3) == 2) then
      ! [i,j,k] -> [i,k,j]

      ! Reallocate permuted memory
      allocate( d( n1, n3, n2))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( i,k,j) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 2 .and. map( 2) == 1 .and. map( 3) == 3) then
      ! [i,j,k] -> [j,i,k]

      ! Reallocate permuted memory
      allocate( d( n2, n1, n3))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( j,i,k) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 2 .and. map( 2) == 3 .and. map( 3) == 1) then
      ! [i,j,k] -> [j,k,i]

      ! Reallocate permuted memory
      allocate( d( n2, n3, n1))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( j,k,i) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 3 .and. map( 2) == 1 .and. map( 3) == 2) then
      ! [i,j,k] -> [k,i,j]

      ! Reallocate permuted memory
      allocate( d( n3, n1, n2))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( k,i,j) = d_temp( i,j,k)
      end do
      end do
      end do

    elseif (map( 1) == 3 .and. map( 2) == 2 .and. map( 3) == 1) then
      ! [i,j,k] -> [k,j,i]

      ! Reallocate permuted memory
      allocate( d( n3, n2, n1))

      ! Copy and permuted data from temporary memory
      do i = 1, n1
      do j = 1, n2
      do k = 1, n3
        d( k,j,i) = d_temp( i,j,k)
      end do
      end do
      end do

    else
      call crash('invalid permutation!')
    end if

    ! deallocate temporary memory
    deallocate( d_temp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine permute_3D_dp

  subroutine flip_1D_int( d)
    ! Flip a 1-D array

    ! In/output variables:
    integer,  dimension(:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_1D_int'
    integer                       :: i,nx,iopp

    ! Add routine to path
    call init_routine( routine_name)

    nx = size( d,1)

    ! Flip the data
    do i = 1, nx
      iopp = nx + 1 - i
      if (iopp <= i) exit         ! [a  ] [b  ]
      d( i   ) = d( i) + d( iopp) ! [a+b] [b  ]
      d( iopp) = d( i) - d( iopp) ! [a+b] [a  ]
      d( i   ) = d( i) - d( iopp) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_1D_int

  subroutine flip_1D_dp( d)
    ! Flip a 1-D array

    ! In/output variables:
    real(dp), dimension(:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_1D_dp'
    integer                       :: i,nx,iopp

    ! Add routine to path
    call init_routine( routine_name)

    nx = size( d,1)

    ! Flip the data
    do i = 1, nx
      iopp = nx + 1 - i
      if (iopp <= i) exit         ! [a  ] [b  ]
      d( i   ) = d( i) + d( iopp) ! [a+b] [b  ]
      d( iopp) = d( i) - d( iopp) ! [a+b] [a  ]
      d( i   ) = d( i) - d( iopp) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_1D_dp

  subroutine flip_2D_x1_int( d)
    ! Flip a 2-D array along the first dimension

    ! In/output variables:
    integer,  dimension(:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_2D_x1_int'
    integer                       :: i,n1,n2,iopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)

    ! Flip the data
    do i = 1, n1
      iopp = n1 + 1 - i
      if (iopp <= i) exit               ! [a  ] [b  ]
      d( i   ,:) = d( i,:) + d( iopp,:) ! [a+b] [b  ]
      d( iopp,:) = d( i,:) - d( iopp,:) ! [a+b] [a  ]
      d( i   ,:) = d( i,:) - d( iopp,:) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_2D_x1_int

  subroutine flip_2D_x1_dp( d)
    ! Flip a 2-D array along the first dimension

    ! In/output variables:
    real(dp), dimension(:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_2D_x1_dp'
    integer                       :: i,n1,n2,iopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)

    ! Flip the data
    do i = 1, n1
      iopp = n1 + 1 - i
      if (iopp <= i) exit               ! [a  ] [b  ]
      d( i   ,:) = d( i,:) + d( iopp,:) ! [a+b] [b  ]
      d( iopp,:) = d( i,:) - d( iopp,:) ! [a+b] [a  ]
      d( i   ,:) = d( i,:) - d( iopp,:) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_2D_x1_dp

  subroutine flip_2D_x2_int( d)
    ! Flip a 2-D array along the second dimension

    ! In/output variables:
    integer,  dimension(:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_2D_x2_int'
    integer                       :: j,n1,n2,jopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)

    ! Flip the data
    do j = 1, n2
      jopp = n2 + 1 - j
      if (jopp <= j) exit               ! [a  ] [b  ]
      d( :,j   ) = d( :,j) + d( :,jopp) ! [a+b] [b  ]
      d( :,jopp) = d( :,j) - d( :,jopp) ! [a+b] [a  ]
      d( :,j   ) = d( :,j) - d( :,jopp) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_2D_x2_int

  subroutine flip_2D_x2_dp( d)
    ! Flip a 2-D array along the second dimension

    ! In/output variables:
    real(dp), dimension(:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_2D_x2_dp'
    integer                       :: j,n1,n2,jopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)

    ! Flip the data
    do j = 1, n2
      jopp = n2 + 1 - j
      if (jopp <= j) exit               ! [a  ] [b  ]
      d( :,j   ) = d( :,j) + d( :,jopp) ! [a+b] [b  ]
      d( :,jopp) = d( :,j) - d( :,jopp) ! [a+b] [a  ]
      d( :,j   ) = d( :,j) - d( :,jopp) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_2D_x2_dp

  subroutine flip_3D_x1_int( d)
    ! Flip a 3-D array along the first dimension

    ! In/output variables:
    integer,  dimension(:,:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_3D_x1_int'
    integer                       :: i,n1,n2,n3,iopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! Flip the data
    do i = 1, n1
      iopp = n1 + 1 - i
      if (iopp <= i) exit                     ! [a  ] [b  ]
      d( i   ,:,:) = d( i,:,:) + d( iopp,:,:) ! [a+b] [b  ]
      d( iopp,:,:) = d( i,:,:) - d( iopp,:,:) ! [a+b] [a  ]
      d( i   ,:,:) = d( i,:,:) - d( iopp,:,:) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_3D_x1_int

  subroutine flip_3D_x1_dp( d)
    ! Flip a 3-D array along the first dimension

    ! In/output variables:
    real(dp), dimension(:,:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_3D_x1_dp'
    integer                       :: i,n1,n2,n3,iopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! Flip the data
    do i = 1, n1
      iopp = n1 + 1 - i
      if (iopp <= i) exit                     ! [a  ] [b  ]
      d( i   ,:,:) = d( i,:,:) + d( iopp,:,:) ! [a+b] [b  ]
      d( iopp,:,:) = d( i,:,:) - d( iopp,:,:) ! [a+b] [a  ]
      d( i   ,:,:) = d( i,:,:) - d( iopp,:,:) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_3D_x1_dp

  subroutine flip_3D_x2_int( d)
    ! Flip a 3-D array along the second dimension

    ! In/output variables:
    integer,  dimension(:,:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_3D_x2_int'
    integer                       :: j,n1,n2,n3,jopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! Flip the data
    do j = 1, n2
      jopp = n2 + 1 - j
      if (jopp <= j) exit                     ! [a  ] [b  ]
      d( :,j   ,:) = d( :,j,:) + d( :,jopp,:) ! [a+b] [b  ]
      d( :,jopp,:) = d( :,j,:) - d( :,jopp,:) ! [a+b] [a  ]
      d( :,j   ,:) = d( :,j,:) - d( :,jopp,:) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_3D_x2_int

  subroutine flip_3D_x2_dp( d)
    ! Flip a 3-D array along the second dimension

    ! In/output variables:
    real(dp), dimension(:,:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_3D_x2_dp'
    integer                       :: j,n1,n2,n3,jopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! Flip the data
    do j = 1, n2
      jopp = n2 + 1 - j
      if (jopp <= j) exit                     ! [a  ] [b  ]
      d( :,j   ,:) = d( :,j,:) + d( :,jopp,:) ! [a+b] [b  ]
      d( :,jopp,:) = d( :,j,:) - d( :,jopp,:) ! [a+b] [a  ]
      d( :,j   ,:) = d( :,j,:) - d( :,jopp,:) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_3D_x2_dp

  subroutine flip_3D_x3_int( d)
    ! Flip a 3-D array along the third dimension

    ! In/output variables:
    integer,  dimension(:,:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_3D_x3_int'
    integer                       :: k,n1,n2,n3,kopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! Flip the data
    do k = 1, n3
      kopp = n3 + 1 - k
      if (kopp <= k) exit                     ! [a  ] [b  ]
      d( :,:,k   ) = d( :,:,k) + d( :,:,kopp) ! [a+b] [b  ]
      d( :,:,kopp) = d( :,:,k) - d( :,:,kopp) ! [a+b] [a  ]
      d( :,:,k   ) = d( :,:,k) - d( :,:,kopp) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_3D_x3_int

  subroutine flip_3D_x3_dp( d)
    ! Flip a 3-D array along the third dimension

    ! In/output variables:
    real(dp), dimension(:,:,:), intent(inout) :: d

    ! Local variables:
    character(len=256), parameter :: routine_name = 'flip_3D_x3_dp'
    integer                       :: k,n1,n2,n3,kopp

    ! Add routine to path
    call init_routine( routine_name)

    n1 = size( d,1)
    n2 = size( d,2)
    n3 = size( d,3)

    ! Flip the data
    do k = 1, n3
      kopp = n3 + 1 - k
      if (kopp <= k) exit                     ! [a  ] [b  ]
      d( :,:,k   ) = d( :,:,k) + d( :,:,kopp) ! [a+b] [b  ]
      d( :,:,kopp) = d( :,:,k) - d( :,:,kopp) ! [a+b] [a  ]
      d( :,:,k   ) = d( :,:,k) - d( :,:,kopp) ! [b  ] [a  ]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_3D_x3_dp

  ! Remapping of a 1-D variable (2nd-order conservative)
  subroutine remap_cons_2nd_order_1D( z_src, mask_src, d_src, z_dst, mask_dst, d_dst)
    ! 2nd-order conservative remapping of a 1-D variable
    !
    ! Used to remap ocean data from the provided vertical grid to the UFEMISM ocean vertical grid
    !
    ! Both z_src and z_dst can be irregular.
    !
    ! Both the src and dst data have a mask, with 0 indicating grid points where no data is defined.
    !
    ! This subroutine is serial, as it will be applied to single grid cells when remapping 3-D data fields,
    !   with the parallelisation being done by distributing the 2-D grid cells over the processes.

    ! In/output variables:
    real(dp), dimension(:), intent(in   ) :: z_src
    integer,  dimension(:), intent(in   ) :: mask_src
    real(dp), dimension(:), intent(in   ) :: d_src
    real(dp), dimension(:), intent(in   ) :: z_dst
    integer,  dimension(:), intent(in   ) :: mask_dst
    real(dp), dimension(:), intent(  out) :: d_dst

    ! Local variables:
    logical                                 :: all_are_masked
    integer                                 :: nz_src, nz_dst
    integer                                 :: k
    real(dp), dimension(:    ), allocatable :: ddz_src
    integer                                 :: k_src, k_dst
    real(dp)                                :: zl_src, zu_src, zl_dst, zu_dst, z_lo, z_hi, z, d
    real(dp)                                :: dz_overlap, dz_overlap_tot, d_int, d_int_tot
    real(dp)                                :: dist_to_dst, dist_to_dst_min, max_dist
    integer                                 :: k_src_nearest_to_dst

    ! Initialise
    d_dst = 0._dp

    ! sizes
    nz_src = size( z_src,1)
    nz_dst = size( z_dst,1)

    ! Maximum distance on combined grids
    max_dist = maxval([ abs( z_src( nz_src) - z_src( 1)), &
                        abs( z_dst( nz_dst) - z_dst( 1)), &
                        abs( z_src( nz_src) - z_dst( 1)), &
                        abs( z_dst( nz_dst) - z_src( 1))])

    ! Exception for when the entire src field is masked
    all_are_masked = .true.
    do k = 1, nz_src
      if (mask_src( k) == 1) all_are_masked = .false.
    end do
    if (all_are_masked) return

    ! Exception for when the entire dst field is masked
    all_are_masked = .true.
    do k = 1, nz_dst
      if (mask_dst( k) == 1) all_are_masked = .false.
    end do
    if (all_are_masked) return

    ! Calculate derivative d_src/dz (one-sided differencing at the boundary, central differencing everywhere else)
    allocate( ddz_src( nz_src))
    do k = 2, nz_src-1
      ddz_src( k    ) = (d_src( k+1   ) - d_src( k-1     )) / (z_src( k+1   ) - z_src( k-1     ))
    end do
    ddz_src(  1     ) = (d_src( 2     ) - d_src( 1       )) / (z_src( 2     ) - z_src( 1       ))
    ddz_src(  nz_src) = (d_src( nz_src) - d_src( nz_src-1)) / (z_src( nz_src) - z_src( nz_src-1))

    ! Perform conservative remapping by finding regions of overlap
    ! between source and destination grid cells

    do k_dst = 1, nz_dst

      ! Skip masked grid cells
      if (mask_dst( k_dst) == 0) then
        d_dst( k_dst) = 0._dp
        cycle
      end if

      ! Find z range covered by this dst grid cell
      if (k_dst > 1) then
        zl_dst = 0.5_dp * (z_dst( k_dst - 1) + z_dst( k_dst))
      else
        zl_dst = z_dst( 1) - 0.5_dp * (z_dst( 2) - z_dst( 1))
      end if
      if (k_dst < nz_dst) then
        zu_dst = 0.5_dp * (z_dst( k_dst + 1) + z_dst( k_dst))
      else
        zu_dst = z_dst( nz_dst) + 0.5_dp * (z_dst( nz_dst) - z_dst( nz_dst-1))
      end if

      ! Find all overlapping src grid cells
      d_int_tot      = 0._dp
      dz_overlap_tot = 0._dp
      do k_src = 1, nz_src

        ! Skip masked grid cells
        if (mask_src( k_src) == 0) cycle

        ! Find z range covered by this src grid cell
        if (k_src > 1) then
          zl_src = 0.5_dp * (z_src( k_src - 1) + z_src( k_src))
        else
          zl_src = z_src( 1) - 0.5_dp * (z_src( 2) - z_src( 1))
        end if
        if (k_src < nz_src) then
          zu_src = 0.5_dp * (z_src( k_src + 1) + z_src( k_src))
        else
          zu_src = z_src( nz_src) + 0.5_dp * (z_src( nz_src) - z_src( nz_src-1))
        end if

        ! Find region of overlap
        z_lo = max( zl_src, zl_dst)
        z_hi = min( zu_src, zu_dst)
        dz_overlap = max( 0._dp, z_hi - z_lo)

        ! Calculate integral over region of overlap and add to sum
        if (dz_overlap > 0._dp) then
          z = 0.5_dp * (z_lo + z_hi)
          d = d_src( k_src) + ddz_src( k_src) * (z - z_src( k_src))
          d_int = d * dz_overlap

          d_int_tot      = d_int_tot      + d_int
          dz_overlap_tot = dz_overlap_tot + dz_overlap
        end if

      end do ! do k_src = 1, nz_src

      if (dz_overlap_tot > 0._dp) then
        ! Calculate dst value
        d_dst( k_dst) = d_int_tot / dz_overlap_tot
      else
        ! Exception for when no overlapping src grid cells were found; use nearest-neighbour extrapolation

        k_src_nearest_to_dst = 0._dp
        dist_to_dst_min      = max_dist
        do k_src = 1, nz_src
          if (mask_src( k_src) == 1) then
            dist_to_dst = abs( z_src( k_src) - z_dst( k_dst))
            if (dist_to_dst < dist_to_dst_min) then
              dist_to_dst_min      = dist_to_dst
              k_src_nearest_to_dst = k_src
            end if
          end if
        end do

        ! Safety
        if (k_src_nearest_to_dst == 0) then
          write(0,*) '  remap_cons_2nd_order_1D - ERROR: couldnt find nearest neighbour on source grid!'
          call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        end if

        d_dst( k_dst) = d_src( k_src_nearest_to_dst)

      end if ! if (dz_overlap_tot > 0._dp) then

    end do ! do k_dst = 1, nz_dst

    ! Clean up after yourself
    deallocate( ddz_src)

  end subroutine remap_cons_2nd_order_1D

  subroutine interpolate_inside_triangle_dp_2D( pa, pb, pc, fa, fb, fc, p, f, tol_dist)

    ! In/output variables
    real(dp), dimension(2), intent(in   ) :: pa, pb, pc
    real(dp),               intent(in   ) :: fa, fb, fc
    real(dp), dimension(2), intent(in   ) :: p
    real(dp),               intent(  out) :: f
    real(dp),               intent(in   ) :: tol_dist

    ! Local variables
    real(dp) :: Atri_abp, Atri_bcp, Atri_cap, Atri_tot, wa, wb, wc

#if (DO_ASSERTIONS)
    call assert( is_in_triangle( pa, pb, pc, p) .or. &
      lies_on_line_segment( pa, pb, p, tol_dist) .or. &
      lies_on_line_segment( pb, pc, p, tol_dist) .or. &
      lies_on_line_segment( pc, pa, p, tol_dist) .or. &
      norm2( pa - p) <= tol_dist .or. &
      norm2( pb - p) <= tol_dist .or. &
      norm2( pc - p) <= tol_dist, 'p does not lie in triangle')
#endif

    ! If p coincides with a, b, or c, copy f from there
    if (norm2( pa - p) <= tol_dist) then
      f = fa
    elseif (norm2( pb - p) <= tol_dist) then
      f = fb
    elseif (norm2( pc - p) <= tol_dist) then
      f = fc

    ! If p lies on one of the three edges, interpolate between its two vertices
    elseif (lies_on_line_segment( pa, pb, p, tol_dist)) then
      wa = norm2( pb - p) / norm2( pb - pa)
      wb = 1._dp - wa
      f = wa * fa + wb * fb
    elseif (lies_on_line_segment( pb, pc, p, tol_dist)) then
      wb = norm2( pc - p) / norm2( pc - pb)
      wc = 1._dp - wb
      f = wb * fb + wc * fc
    elseif (lies_on_line_segment( pc, pa, p, tol_dist)) then
      wc = norm2( pa - p) / norm2( pa - pc)
      wa = 1._dp - wc
      f = wc * fc + wa * fa

    ! Otherwise, p lies inside the triangle; do a trilinear interpolation
    else
      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_tot = Atri_abp + Atri_bcp + Atri_cap
      wa = Atri_bcp / Atri_tot
      wb = Atri_cap / Atri_tot
      wc = Atri_abp / Atri_tot
      f = wa * fa + wb * fb + wc * fc
    end if

  end subroutine interpolate_inside_triangle_dp_2D

  subroutine interpolate_inside_triangle_dp_3D( pa, pb, pc, fa, fb, fc, p, f, tol_dist)
    ! In/output variables
    real(dp), dimension(2), intent(in   ) :: pa, pb, pc
    real(dp), dimension(:), intent(in   ) :: fa, fb, fc
    real(dp), dimension(2), intent(in   ) :: p
    real(dp), dimension(:), intent(  out) :: f
    real(dp),               intent(in   ) :: tol_dist
    ! Local variables
    real(dp) :: Atri_abp, Atri_bcp, Atri_cap, Atri_tot, wa, wb, wc

#if (DO_ASSERTIONS)
    call assert( is_in_triangle( pa, pb, pc, p) .or. &
      lies_on_line_segment( pa, pb, p, tol_dist) .or. &
      lies_on_line_segment( pb, pc, p, tol_dist) .or. &
      lies_on_line_segment( pc, pa, p, tol_dist) .or. &
      norm2( pa - p) <= tol_dist .or. &
      norm2( pb - p) <= tol_dist .or. &
      norm2( pc - p) <= tol_dist, 'p does not lie in triangle')
    call assert( size( fa) == size( fb) .and. size( fb) == size( fc) .and. size( fc) == size( f), &
      'incorrect input dimensions')
#endif

    ! If p coincides with a, b, or c, copy f from there
    if (norm2( pa - p) <= tol_dist) then
      f = fa
    elseif (norm2( pb - p) <= tol_dist) then
      f = fb
    elseif (norm2( pc - p) <= tol_dist) then
      f = fc

    ! If p lies on one of the three edges, interpolate between its two vertices
    elseif (lies_on_line_segment( pa, pb, p, tol_dist)) then
      wa = norm2( pb - p) / norm2( pb - pa)
      wb = 1._dp - wa
      f = wa * fa + wb * fb
    elseif (lies_on_line_segment( pb, pc, p, tol_dist)) then
      wb = norm2( pc - p) / norm2( pc - pb)
      wc = 1._dp - wb
      f = wb * fb + wc * fc
    elseif (lies_on_line_segment( pc, pa, p, tol_dist)) then
      wc = norm2( pa - p) / norm2( pa - pc)
      wa = 1._dp - wc
      f = wc * fc + wa * fa

    ! Otherwise, p lies inside the triangle; do a trilinear interpolation
    else
      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_tot = Atri_abp + Atri_bcp + Atri_cap
      wa = Atri_bcp / Atri_tot
      wb = Atri_cap / Atri_tot
      wc = Atri_abp / Atri_tot
      f = wa * fa + wb * fb + wc * fc
    end if

  end subroutine interpolate_inside_triangle_dp_3D

end module math_utilities
