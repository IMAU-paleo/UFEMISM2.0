MODULE climate_model_types

  ! The different data types used in the climate modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_climate_model
    ! The climate model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T2m                     ! [K]      Monthly 2-m air temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Precip                  ! [m.w.e.] Monthly precipitation
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                      ! [m] orography
    
    ! lapse rates for GCM snapshots 
    REAL(dp), DIMENSION(:  ), ALLOCATABLE     :: lambda ! Spatially variable (see Berends et al., 2018)
    REAL(dp)                                  :: lapse_rate_precip  ! single-value per region (precipitation)
    REAL(dp)                                  :: lapse_rate_temp    ! single-value per region (precipitation)

    ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Q_TOA                   ! [W/m2] Monthly mean insolation at the top of the atmosphere (taken from the prescribed insolation solution at orbit_time)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Albedo                  ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: I_abs                   ! Total yearly absorbed insolation, used in the climate matrix for interpolation

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename        ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_climate_model

  TYPE type_global_forcing
    ! Data structure containing model forcing data - CO2 record, d18O record, (global) insolation record

      ! CO2 record
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: CO2_time
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: CO2_record
      REAL(dp),                   ALLOCATABLE     :: CO2_obs
      INTEGER :: wCO2_time, wCO2_record, wCO2_obs

      ! d18O record
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: d18O_time
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: d18O_record
      REAL(dp),                   ALLOCATABLE     :: d18O_obs
      REAL(dp),                   ALLOCATABLE     :: d18O_obs_PD
      INTEGER :: wd18O_time, wd18O_record, wd18O_obs, wd18O_obs_PD

      ! Insolation reconstruction
      !TYPE(type_netcdf_insolation)                :: netcdf_ins
      INTEGER,                    ALLOCATABLE     :: ins_nyears
      INTEGER,                    ALLOCATABLE     :: ins_nlat,ins_nlon
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: ins_time
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: ins_lat
      REAL(dp),                   ALLOCATABLE     :: ins_t0, ins_t1
      INTEGER,                    ALLOCATABLE     :: ins_ti0,ins_ti1
      REAL(dp), DIMENSION(:,:  ), ALLOCATABLE     :: ins_Q_TOA0, ins_Q_TOA1
      
      
      INTEGER :: wins_nyears, wins_nlat, wins_time, wins_lat, wins_t0, wins_t1, wins_Q_TOA0, wins_Q_TOA1
      

  END TYPE type_global_forcing

CONTAINS

END MODULE climate_model_types