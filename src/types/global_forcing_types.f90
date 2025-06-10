MODULE global_forcing_types

  ! The different data types used as global forcings for all model regions

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

TYPE type_global_forcing
    ! Data structure containing model forcing data - CO2 record, d18O record, (global) insolation record

      ! CO2 record
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: CO2_time
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: CO2_record
      REAL(dp),                   ALLOCATABLE     :: CO2_obs
      !INTEGER :: wCO2_time, wCO2_record, wCO2_obs

      ! d18O record
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: d18O_time
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: d18O_record
      REAL(dp),                   ALLOCATABLE     :: d18O_obs
      REAL(dp),                   ALLOCATABLE     :: d18O_obs_PD
      !INTEGER :: wd18O_time, wd18O_record, wd18O_obs, wd18O_obs_PD

      ! Global mean sea level
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: sea_level_time
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: sea_level_record 
      REAL(dp),                   ALLOCATABLE     :: sl_t0, sl_t1, sl_at_t0, sl_at_t1
      

  END TYPE type_global_forcing

  CONTAINS

END MODULE global_forcing_types