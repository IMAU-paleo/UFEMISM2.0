module idealised_geometries

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash
  use model_configuration, only: C
  use ice_geometry_basics, only: ice_surface_elevation
  use parameters, only: pi
  use Halfar_SIA_solution, only: Halfar
  use Bueler_SIA_solution, only: Bueler_dome

  implicit none

  private

  public :: calc_idealised_geometry

contains

  subroutine calc_idealised_geometry( x, y, Hi, Hb, Hs, SL, choice_refgeo_idealised)
    !< Calculate an idealised geometry

    ! In/output variables:
    real(dp),            intent(in   ) :: x,y             ! [m] Coordinates
    real(dp),            intent(  out) :: Hi              ! [m] Ice thickness
    real(dp),            intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp),            intent(  out) :: Hs              ! [m] Surface elevation
    real(dp),            intent(  out) :: SL              ! [m] Sea level
    character(len=1024), intent(in   ) :: choice_refgeo_idealised

    ! Calculated the specified idealised geometry
    select case (choice_refgeo_idealised)
    case default
      call crash('unknown choice_refgeo_idealised "' // trim( choice_refgeo_idealised))
    case ('flatearth')
      call calc_idealised_geometry_flatearth( Hi, Hb, Hs, SL)
    case ('slabonaslope')
      call calc_idealised_geometry_slabonaslope( x, Hi, Hb, Hs, SL)
    case ('Halfar')
      call calc_idealised_geometry_Halfar( x, y, Hi, Hb, Hs, SL)
    case ('Bueler')
      call calc_idealised_geometry_Bueler( x, y, Hi, Hb, Hs, SL)
    case ('SSA_icestream')
      call calc_idealised_geometry_SSA_icestream( x, Hi, Hb, Hs, SL)
    case ('MISMIP_mod')
      call calc_idealised_geometry_MISMIP_mod( x, y, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_A')
      call calc_idealised_geometry_ISMIP_HOM_A( x, y, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_B')
      call calc_idealised_geometry_ISMIP_HOM_B( x, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_C', 'ISMIP-HOM_D')
      call calc_idealised_geometry_ISMIP_HOM_CD( x, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_E')
      call crash('ISMIP-HOM E is not implemented in UFEMISM!')
    case ('ISMIP-HOM_F')
      call calc_idealised_geometry_ISMIP_HOM_F( x, y, Hi, Hb, Hs, SL)
    case ('MISMIP+', 'MISMIPplus')
      call calc_idealised_geometry_MISMIPplus( x, y, Hi, Hb, Hs, SL)
    case ('calvmip_circular')
      call calc_idealised_geometry_CalvMIP_circular( x, y, Hi, Hb, Hs, SL)
    case ('calvmip_Thule')
      call calc_idealised_geometry_CalvMIP_Thule( x, y, Hi, Hb, Hs, SL)
    end select

  end subroutine calc_idealised_geometry

  subroutine calc_idealised_geometry_flatearth( Hi, Hb, Hs, SL)
    !< Calculate a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments

    ! In/output variables:
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    Hi = 0._dp
    Hb = 0._dp
    Hs = 0._dp
    SL = -10000._dp

  end subroutine calc_idealised_geometry_flatearth

  subroutine calc_idealised_geometry_slabonaslope( x, Hi, Hb, Hs, SL)
    !< Calculate a 2,000 m thick slab of ice on a flat, inclined plane

    ! In/output variables:
    real(dp), intent(in   ) :: x               ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_slabonaslope_Hi < 100._dp .or. &
        C%refgeo_idealised_slabonaslope_Hi > 10000._dp) then
      call crash('refgeo_idealised_slabonaslope_Hi has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_slabonaslope_Hi)
    end if
    if ( ABS( C%refgeo_idealised_slabonaslope_dhdx) > 0.1_dp) then
      call crash('refgeo_idealised_slabonaslope_dhdx has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_slabonaslope_dhdx)
    end if
#endif

    Hi = C%refgeo_idealised_slabonaslope_Hi
    Hb = C%refgeo_idealised_slabonaslope_dhdx * x
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_slabonaslope

  subroutine calc_idealised_geometry_Halfar( x, y, Hi, Hb, Hs, SL)
    !< Calculate the Halfar dome solution at t = 0

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%uniform_Glens_flow_factor < 1E-18_dp .OR. &
        C%uniform_Glens_flow_factor > 1E-15_dp) then
      call crash('uniform_flow_factor has unrealistic value of {dp_01}!', dp_01 = C%uniform_Glens_flow_factor)
    end if
    if ( C%Glens_flow_law_exponent < 1._dp .OR. C%Glens_flow_law_exponent > 5._dp) then
      call crash('Glens_flow_law_exponent has unrealistic value of {dp_01}!', dp_01 = C%Glens_flow_law_exponent)
    end if
    if ( C%refgeo_idealised_Halfar_H0 < 100._dp .OR. C%refgeo_idealised_Halfar_H0 > 10000._dp) then
      call crash('refgeo_idealised_Halfar_H0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Halfar_H0)
    end if
    if ( C%refgeo_idealised_Halfar_R0 < 100E3_dp .OR. C%refgeo_idealised_Halfar_R0 > 5000E3_dp) then
      call crash('refgeo_idealised_Halfar_R0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Halfar_R0)
    end if
#endif

    Hi = Halfar%H( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_Halfar_H0, C%refgeo_idealised_Halfar_R0, &
      x, y, 0._dp)
    Hb = 0._dp
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_Halfar

  subroutine calc_idealised_geometry_Bueler( x, y, Hi, Hb, Hs, SL)
    !< Calculate the Bueler dome solution at t = 0

    ! In/output variables:
    real(dp),                       intent(in   ) :: x,y             ! [m] Coordinates
    real(dp),                       intent(  out) :: Hi              ! [m] Ice thickness
    real(dp),                       intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp),                       intent(  out) :: Hs              ! [m] Surface elevation
    real(dp),                       intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp) :: M

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%uniform_Glens_flow_factor < 1E-18_dp .OR. &
        C%uniform_Glens_flow_factor > 1E-15_dp) then
      call crash('uniform_flow_factor has unrealistic value of {dp_01}!', dp_01 = C%uniform_Glens_flow_factor)
    end if
    if ( C%Glens_flow_law_exponent < 1._dp .OR. C%Glens_flow_law_exponent > 5._dp) then
      call crash('Glens_flow_law_exponent has unrealistic value of {dp_01}!', dp_01 = C%Glens_flow_law_exponent)
    end if
    if ( C%refgeo_idealised_Bueler_H0 < 100._dp .OR. C%refgeo_idealised_Bueler_H0 > 10000._dp) then
      call crash('refgeo_idealised_Bueler_H0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_H0)
    end if
    if ( C%refgeo_idealised_Bueler_R0 < 100E3_dp .OR. C%refgeo_idealised_Bueler_R0 > 5000E3_dp) then
      call crash('refgeo_idealised_Bueler_R0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_R0)
    end if
    if ( C%refgeo_idealised_Bueler_lambda < -10._dp .OR. C%refgeo_idealised_Bueler_lambda > 10._dp) then
      call crash('refgeo_idealised_Bueler_lambda has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_lambda)
    end if
#endif

    call Bueler_dome( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_Bueler_H0, C%refgeo_idealised_Bueler_R0, &
      C%refgeo_idealised_Bueler_lambda, x, y, 0._dp, Hi, M)
    Hb = 0._dp
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_Bueler

  subroutine calc_idealised_geometry_SSA_icestream( x, Hi, Hb, Hs, SL)
    !< Calculate a thick slab of ice on a flat, inclined plane

    ! In/output variables:
    real(dp), intent(in   ) :: x               ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_SSA_icestream_Hi < 100._dp .OR. &
        C%refgeo_idealised_SSA_icestream_Hi > 10000._dp) then
      call crash('refgeo_idealised_SSA_icestream_Hi has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_SSA_icestream_Hi)
    end if
    if ( abs( C%refgeo_idealised_SSA_icestream_dhdx) > 0.1_dp) then
      call crash('refgeo_idealised_SSA_icestream_dhdx has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_SSA_icestream_dhdx)
    end if
#endif

    Hi = C%refgeo_idealised_SSA_icestream_Hi
    Hb = C%refgeo_idealised_SSA_icestream_dhdx * x
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_SSA_icestream

  subroutine calc_idealised_geometry_MISMIP_mod( x, y, Hi, Hb, Hs, SL)
    !< Calculate the MISMIP_mod cone-shaped island

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_MISMIP_mod_Hi_init < 0._dp .OR. &
        C%refgeo_idealised_MISMIP_mod_Hi_init > 10000._dp) then
      call crash('refgeo_idealised_MISMIP_mod_Hi_init has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_MISMIP_mod_Hi_init)
    end if
#endif

    if (sqrt( x**2 + y**2) > 900E3_dp) then
      Hi = 0._dp
    else
      Hi = C%refgeo_idealised_MISMIP_mod_Hi_init
    end if
    Hb = 150._dp - 400._dp * sqrt( x**2 + y**2)/ 750000._dp
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_MISMIP_mod

  subroutine calc_idealised_geometry_ISMIP_HOM_A( x, y, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment A (slab on a bumpy slope in both directions)

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
        C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 2000._dp - x * tan( 0.5_dp * pi / 180._dp)
    Hb = Hs - 1000._dp + 500._dp * sin( x * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L) &
                                * sin( y * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L)
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_A

  subroutine calc_idealised_geometry_ISMIP_HOM_B( x, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment B (slab on a bumpy slope in only the x-directions)

    ! In/output variables:
    real(dp), intent(in   ) :: x               ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
        C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 2000._dp - x * TAN( 0.5_dp * pi / 180._dp)
    Hb = Hs - 1000._dp + 500._dp * SIN( x * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L)
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_B

  subroutine calc_idealised_geometry_ISMIP_HOM_CD( x, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment C/D (slab on a slope)

    ! In/output variables:
    real(dp),                       intent(in   ) :: x               ! [m] Coordinates
    real(dp),                       intent(  out) :: Hi              ! [m] Ice thickness
    real(dp),                       intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp),                       intent(  out) :: Hs              ! [m] Surface elevation
    real(dp),                       intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
        C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 2000._dp - x * TAN( 0.1_dp * pi / 180._dp)
    Hb = Hs - 1000._dp
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_CD

  subroutine calc_idealised_geometry_ISMIP_HOM_F( x, y, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment F (slab on another bumpy slope)

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp), parameter :: H0    = 1000._dp
    real(dp), parameter :: a0    = 100._dp
    real(dp), parameter :: sigma = 10000._dp

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
        C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 5000._dp - x * tan( 3._dp * pi / 180._dp)
    Hb = Hs - H0 + a0 * exp( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                + a0 * exp( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2)
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_F

  subroutine calc_idealised_geometry_MISMIPplus( x, y, Hi, Hb, Hs, SL)
    !< Calculate the MISMIP+ geometry

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp)            :: xtilde,Bx,By
    real(dp), parameter :: B0     = -150._dp
    real(dp), parameter :: B2     = -728.8_dp
    real(dp), parameter :: B4     = 343.91_dp
    real(dp), parameter :: B6     = -50.57_dp
    real(dp), parameter :: xbar   = 300000._dp
    real(dp), parameter :: fc     = 4000._dp
    real(dp), parameter :: dc     = 500._dp
    real(dp), parameter :: wc     = 24000._dp
    real(dp), parameter :: zbdeep = -720._dp

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_MISMIPplus_Hi_init < 0._dp .OR. &
        C%refgeo_idealised_MISMIPplus_Hi_init > 10000._dp) then
      call crash('refgeo_idealised_MISMIPplus_Hi_init has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_MISMIPplus_Hi_init)
    end if
#endif

    xtilde = x / xbar
    Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
    By = (dc / (1 + exp(-2._dp*(y - wc)/fc))) + &
        (dc / (1 + exp( 2._dp*(y + wc)/fc)))

    if (x > 640E3_dp) then
      Hi = 0._dp
    else
      Hi = C%refgeo_idealised_MISMIPplus_Hi_init
    end if
    Hb = max( Bx + By, zbdeep)
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_MISMIPplus

  subroutine calc_idealised_geometry_CalvMIP_circular( x, y, Hi, Hb, Hs, SL)
    ! Calculate the geometry for the CalvingMIP circular domain

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp), parameter :: R  = 800E3_dp
    real(dp), parameter :: Bc = 900._dp
    real(dp), parameter :: Bl = -2000._dp
    real(dp), parameter :: Ba = 1100._dp
    real(dp), parameter :: rc = 0._dp
    real(dp)            :: radius, theta

    radius = sqrt(x**2 + y**2)
    theta  = atan2(y,x)

    Hi = 0._dp
    Hb = Bc - (Bc-Bl)*(radius-rc)**2 / (R-rc)**2
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_CalvMIP_circular

  subroutine calc_idealised_geometry_CalvMIP_Thule( x, y, Hi, Hb, Hs, SL)
    ! Calculate the geometry for the CalvingMIP Thule domain

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp), parameter :: R  = 800E3_dp
    real(dp), parameter :: Bc = 900._dp
    real(dp), parameter :: Bl = -2000._dp
    real(dp), parameter :: Ba = 1100._dp
    real(dp), parameter :: rc = 0._dp
    real(dp)            :: radius, theta, l, a

    radius = sqrt(x**2 + y**2)
    theta  = atan2(y,x)
    l = R - cos( 2._dp * theta) * R / 2._dp
    a = Bc - (Bc-Bl)*(radius-rc)**2 / (R-rc)**2

    Hi = 0._dp
    Hb = Ba * cos( 3._dp * pi * radius / l) + a
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_CalvMIP_Thule

end module idealised_geometries
