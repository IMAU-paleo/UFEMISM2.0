module projections

  ! Projections between spherical and Cartesian coordinates
  !
  ! (adapted from the old Oblimap code by Thomas Reerink; see Reerink et al., 2010)

  use precisions, only: dp
  use parameters, only: earth_radius, pi

  implicit none

  private

  public :: oblique_sg_projection, inverse_oblique_sg_projection

contains

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

end module projections
