module scalar_types

  ! The different data types used for scalar data

! ===== Preamble =====
! ====================

  use precisions, only: dp

  implicit none

! ===== Types =====
! =================

  type type_regional_scalars
    ! Data fields storing the regional scalar data

    ! Geometry
    real(dp)   :: ice_area                 ! [km^2]     Total ice sheet area
    real(dp)   :: ice_area_PD              ! [km^2]     Total present_day ice sheet area
    real(dp)   :: ice_volume               ! [m SLE]    Total ice sheet volume
    real(dp)   :: ice_volume_PD            ! [m SLE]    Total present-day ice sheet volume
    real(dp)   :: ice_volume_af            ! [m SLE]    Total ice sheet volume above floatation
    real(dp)   :: ice_volume_af_PD         ! [m SLE]    Total present-day ice sheet volume above floatation
    real(dp)   :: sea_level_contribution   ! [m SLE]    Total contribution to global mean sea level

    ! Integrated fluxes
    real(dp)   :: SMB_total                ! [Gt yr^-1] Area-integrated surface mass balance
    real(dp)   :: BMB_total                ! [Gt yr^-1] Area-integrated basal mass balance
    real(dp)   :: gl_flux                  ! [Gt yr^-1] Total flux through the grounding line
    real(dp)   :: cf_flux                  ! [Gt yr^-1] Total flux through the calving front
    real(dp)   :: out_flux                 ! [Gt yr^-1] Total flux exiting the ice sheet/shelves

    ! Integrated mass balance
    real(dp)   :: int_T2m                  ! [K]        Mean near-surface air temperatures
    real(dp)   :: int_snowfall             ! [m yr^-1]  Area-integrated snowfall rates
    real(dp)   :: int_rainfall             ! [m yr^-1]  Area-integrated rainfall rates
    real(dp)   :: int_melt                 ! [m yr^-1]  Area-integrated melt rates
    real(dp)   :: int_refreezing           ! [m yr^-1]  Area-integrated refreezing rates
    real(dp)   :: int_runoff               ! [m yr^-1]  Area-integrated runoff rates

  end type type_regional_scalars

contains

end module scalar_types