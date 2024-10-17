MODULE parameters

  ! Some mathematical and physical constants

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  REAL(dp), PARAMETER :: pi                               = 3.141592653589793_dp
  REAL(dp), PARAMETER :: sec_per_year                     = 31556943.36_dp            ! = 365.2424 * 24 * 3600
  REAL(dp), PARAMETER :: sec_per_day                      = 86400.0_dp                ! = 24 * 3600
  REAL(dp), PARAMETER :: T0                               = 273.16_dp                 ! [K]                 Triple point of water
  REAL(dp), PARAMETER :: Clausius_Clapeyron_gradient      = 8.7E-04_dp                ! [K m^-1]            Clausius Clapeyron gradient
  REAL(dp), PARAMETER :: grav                             = 9.81_dp                   ! [m s^-2]            Acceleration of gravity
  REAL(dp), PARAMETER :: earth_radius                     = 6.371221E6_dp             ! [m] Earth           Radius
  REAL(dp), PARAMETER :: L_fusion                         = 3.335E+5_dp               ! [J kg-1]            Latent heat of fusion
  REAL(dp), PARAMETER :: ice_density                      =  910.0_dp                 ! [kg m^-3]           Ice density
  REAL(dp), PARAMETER :: freshwater_density               = 1000.0_dp                 ! [kg m^-3]           Freshwater density
  REAL(dp), PARAMETER :: seawater_density                 = 1028.0_dp                 ! [kg m^-3]           Seawater density
  REAL(dp), PARAMETER :: earth_density                    = 5511.57_dp                ! [kg m^-3]           Total mean Earth density
  REAL(dp), PARAMETER :: R_gas                            = 8.314_dp                  ! [J mol^-1 K^-1]     Gas constant
  REAL(dp), PARAMETER :: cp_ocean                         = 3.974E3_dp                ! [J kg^-1 K^-1]      Specific heat capacity of ocean water
  REAL(dp), PARAMETER :: ocean_area                       = 3.611E14_dp               ! [m^2]               World ocean area

! ===== LADDIE parameters ====
! ============================

  REAL(dp), PARAMETER :: freezing_lambda_1                = -5.73E-2_dp               ! [K PSU^-1]          Freezing point salinity coefficient
  REAL(dp), PARAMETER :: freezing_lambda_2                =  8.32E-2_dp               ! [K]                 Freezing point offset
  REAL(dp), PARAMETER :: freezing_lambda_3                =  7.61E-4_dp               ! [K m^-1]            Freezing point depth coefficient
  REAL(dp), PARAMETER :: cp_ice                           = 2009.0_dp                 ! [J kg^-1 K^-1]      Specific heat capacity of ice
  REAL(dp), PARAMETER :: Stanton_number                   = 5.9E-4_dp                 ! []                  Effective thermal Stanton number
  REAL(dp), PARAMETER :: Prandtl_number                   = 13.8_dp                   ! []                  Prandtl number
  REAL(dp), PARAMETER :: Schmidt_number                   = 2432.0_dp                 ! []                  Schmidt number
  REAL(dp), PARAMETER :: molecular_viscosity              = 1.95E-6_dp                ! [m^2 s^-1]          Molecular viscosity

END MODULE parameters
