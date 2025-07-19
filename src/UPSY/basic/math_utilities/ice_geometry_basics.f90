module ice_geometry_basics

  ! Floatation criterion, surface elevation, and thickness above floatation

  use precisions, only: dp
  use parameters

  implicit none

contains

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

end module ice_geometry_basics
