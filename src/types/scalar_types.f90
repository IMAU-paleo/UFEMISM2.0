module scalar_types

  use precisions, only: dp

  implicit none

  private

  public :: type_regional_scalars, type_scalar_output_buffer

  type type_scalar_output_buffer
    !< Memory for buffering scalar output (from every model time step) between output writing intervals

    integer :: n_mem         !< Number of timeframes for which memory has been allocated
    integer :: n             !< Number of timeframes that are currently buffered

    real(dp), dimension(:), allocatable :: time

    real(dp), dimension(:), allocatable :: ice_area
    real(dp), dimension(:), allocatable :: ice_volume
    real(dp), dimension(:), allocatable :: ice_volume_af
    real(dp), dimension(:), allocatable :: ice_area_PD
    real(dp), dimension(:), allocatable :: ice_volume_PD
    real(dp), dimension(:), allocatable :: ice_volume_af_PD

    real(dp), dimension(:), allocatable :: SMB_total
    real(dp), dimension(:), allocatable :: SMB_gr
    real(dp), dimension(:), allocatable :: SMB_fl
    real(dp), dimension(:), allocatable :: SMB_land
    real(dp), dimension(:), allocatable :: SMB_ocean

    real(dp), dimension(:), allocatable :: BMB_total
    real(dp), dimension(:), allocatable :: BMB_gr
    real(dp), dimension(:), allocatable :: BMB_fl
    real(dp), dimension(:), allocatable :: BMB_land
    real(dp), dimension(:), allocatable :: BMB_ocean

    real(dp), dimension(:), allocatable :: LMB_total
    real(dp), dimension(:), allocatable :: LMB_gr
    real(dp), dimension(:), allocatable :: LMB_fl

    real(dp), dimension(:), allocatable :: AMB_total
    real(dp), dimension(:), allocatable :: AMB_gr
    real(dp), dimension(:), allocatable :: AMB_fl
    real(dp), dimension(:), allocatable :: AMB_land
    real(dp), dimension(:), allocatable :: AMB_ocean

    real(dp), dimension(:), allocatable :: gl_flux
    real(dp), dimension(:), allocatable :: cf_gr_flux
    real(dp), dimension(:), allocatable :: cf_fl_flux
    real(dp), dimension(:), allocatable :: margin_land_flux
    real(dp), dimension(:), allocatable :: margin_ocean_flux

    real(dp), dimension(:), allocatable :: dt_ice
    integer,  dimension(:), allocatable :: n_visc_its
    integer,  dimension(:), allocatable :: n_Axb_its

  end type type_scalar_output_buffer

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
    real(dp)   :: SMB_gr                   ! [Gt yr^-1] Area-integrated surface mass balance over grounded ice
    real(dp)   :: SMB_fl                   ! [Gt yr^-1] Area-integrated surface mass balance over floating ice
    real(dp)   :: SMB_land                 ! [Gt yr^-1] Area-integrated surface mass balance over ice-free land
    real(dp)   :: SMB_ocean                ! [Gt yr^-1] Area-integrated surface mass balance over ice-free ocean
    real(dp)   :: BMB_total                ! [Gt yr^-1] Area-integrated basal mass balance
    real(dp)   :: BMB_gr                   ! [Gt yr^-1] Area-integrated basal mass balance over grounded ice
    real(dp)   :: BMB_fl                   ! [Gt yr^-1] Area-integrated basal mass balance over floating ice
    real(dp)   :: BMB_land                 ! [Gt yr^-1] Area-integrated basal mass balance over ice-free land
    real(dp)   :: BMB_ocean                ! [Gt yr^-1] Area-integrated basal mass balance over ice-free ocean
    real(dp)   :: LMB_total                ! [Gt yr^-1] Area-integrated lateral mass balance
    real(dp)   :: LMB_gr                   ! [Gt yr^-1] Area-integrated lateral mass balance over grounded ice
    real(dp)   :: LMB_fl                   ! [Gt yr^-1] Area-integrated lateral mass balance over floating ice
    real(dp)   :: AMB_total                ! [Gt yr^-1] Area-integrated additional mass balance from other sources
    real(dp)   :: AMB_gr                   ! [Gt yr^-1] Area-integrated additional mass balance from other sources over grounded ice
    real(dp)   :: AMB_fl                   ! [Gt yr^-1] Area-integrated additional mass balance from other sources over floating ice
    real(dp)   :: AMB_land                 ! [Gt yr^-1] Area-integrated additional mass balance from other sources over ice-free land
    real(dp)   :: AMB_ocean                ! [Gt yr^-1] Area-integrated additional mass balance from other sources over ice-free ocean
    real(dp)   :: gl_flux                  ! [Gt yr^-1] Total flux through grounding lines
    real(dp)   :: cf_gr_flux               ! [Gt yr^-1] Total flux through grounded calving fronts
    real(dp)   :: cf_fl_flux               ! [Gt yr^-1] Total flux through floating calving fronts
    real(dp)   :: margin_land_flux         ! [Gt yr^-1] Total flux exiting grounded ice margins
    real(dp)   :: margin_ocean_flux        ! [Gt yr^-1] Total flux exiting marine ice margins

    ! Integrated mass balance
    real(dp)   :: int_T2m                  ! [K]        Mean near-surface air temperatures
    real(dp)   :: int_snowfall             ! [m yr^-1]  Area-integrated snowfall rates
    real(dp)   :: int_rainfall             ! [m yr^-1]  Area-integrated rainfall rates
    real(dp)   :: int_melt                 ! [m yr^-1]  Area-integrated melt rates
    real(dp)   :: int_refreezing           ! [m yr^-1]  Area-integrated refreezing rates
    real(dp)   :: int_runoff               ! [m yr^-1]  Area-integrated runoff rates

    ! Memory to buffer data before writing to NetCDF
    type(type_scalar_output_buffer) :: buffer

  end type type_regional_scalars

contains

end module scalar_types