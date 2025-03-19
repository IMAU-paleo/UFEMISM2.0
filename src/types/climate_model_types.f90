MODULE climate_model_types

  ! The different data types used in the climate modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================
  TYPE type_climate_snapshot
    ! A single climate snapshot

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! Metadata
    REAL(dp),                   ALLOCATABLE     :: CO2
    REAL(dp),                   ALLOCATABLE     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   ALLOCATABLE     :: orbit_ecc                     ! Orbital parameters
    REAL(dp),                   ALLOCATABLE     :: orbit_obl
    REAL(dp),                   ALLOCATABLE     :: orbit_pre
    REAL(dp),                   ALLOCATABLE     :: sealevel
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre, wsealevel

    ! Climate data
    REAL(dp), DIMENSION(:  ), ALLOCATABLE     :: Hs                            ! Orography (m w.r.t. PD sea level)
    INTEGER, DIMENSION(:   ), ALLOCATABLE     :: mask_ice                      ! Climate snapshot ice (1) no_ice (1) mask
    INTEGER, DIMENSION(:   ), ALLOCATABLE     :: mask_ocean                    ! Climate snapshot ocean (1) land (0) mask
    INTEGER, DIMENSION(:   ), ALLOCATABLE     :: mask_shelf                    ! Climate snapshot shelf (1) no shelf (0) mask
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: Wind_WE                       ! Monthly mean west-east   wind speed (m/s)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: Wind_SN                       ! Monthly mean south-north wind speed (m/s)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: Wind_LR                       ! Monthly mean wind speed in the x-direction (m/s)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: Wind_DU                       ! Monthly mean wind speed in the y-direction (m/s)
    !INTEGER :: wHs, wmask_ice, wmask_ocean, wmask_shelf, wT2m, wPrecip, wHs_ref, wWind_WE, wWind_SN, wWind_LR, wWind_DU

    ! Spatially variable lapse rate for GCM snapshots (see Berends et al., 2018)
    REAL(dp), DIMENSION(:  ), ALLOCATABLE     :: lambda
    INTEGER :: wlambda

    ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2) (taken from the prescribed insolation solution at orbit_time)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
    REAL(dp), DIMENSION(:  ), ALLOCATABLE     :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation
    !INTEGER :: wQ_TOA, wAlbedo, wI_abs

  END TYPE type_climate_snapshot
  
   TYPE type_climate_model_matrix
    ! The "matrix" climate model option: three GCM snapshots (warm, cold, and PI), and a PD reanalysis snapshot to use for bias correction

    ! The three GCM snapshots
    TYPE(type_climate_snapshot)             :: GCM_PI
    TYPE(type_climate_snapshot)             :: GCM_warm
    TYPE(type_climate_snapshot)             :: GCM_cold

    ! The present-day climate
    TYPE(type_climate_snapshot)             :: PD_obs

    ! GCM bias
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: GCM_bias_T2m
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: GCM_bias_Precip
    REAL(dp), DIMENSION(:  ), ALLOCATABLE     :: GCM_bias_Hs
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: GCM_bias_Wind_LR
    REAL(dp), DIMENSION(:,:), ALLOCATABLE     :: GCM_bias_Wind_DU
!    INTEGER :: wGCM_bias_T2m, wGCM_bias_Precip, wGCM_bias_Hs, wGCM_bias_Wind_LR, wGCM_bias_Wind_DU

   ! Total yearly absorbed insolation, used in the climate matrix for interpolation
    REAL(dp), DIMENSION(:  ), ALLOCATABLE     :: I_abs

  END TYPE type_climate_model_matrix
  
  TYPE type_climate_model
    ! The climate model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T2m                     ! [K]      Monthly 2-m air temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Precip                  ! [m.w.e.] Monthly precipitation
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                      ! [m] orography
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Wind_LR                 ! [m/s]    Monthly mean wind speed in the x-direction 
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Wind_DU                 ! [m/s]    Monthly mean wind speed in the y-direction 
    
    ! Spatially variable lapse rate for GCM snapshots (see Berends et al., 2018)
    REAL(dp), DIMENSION(:  ), ALLOCATABLE     :: lambda
    INTEGER :: wlambda

    ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Q_TOA                   ! [W/m2] Monthly mean insolation at the top of the atmosphere (taken from the prescribed insolation solution at orbit_time)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Albedo                  ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: I_abs                   ! Total yearly absorbed insolation, used in the climate matrix for interpolation

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename        ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next
    
    ! Add different climate model options
    TYPE(type_climate_model_matrix)         :: matrix             ! The "matrix"          climate model option: three GCM snapshots (warm, cold, and PI), and a PD reanalysis snapshot to use for bias correction

  END TYPE type_climate_model
  
  TYPE type_global_forcing
    ! Data structure containing model forcing data - CO2 record, d18O record, (global) insolation record

      ! CO2 record
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: CO2_time
      REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: CO2_record
      REAL(dp),                   ALLOCATABLE     :: CO2_obs

! not implemented yet, commented for now
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