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

  END TYPE type_climate_model

CONTAINS

END MODULE climate_model_types