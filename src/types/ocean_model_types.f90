MODULE ocean_model_types

  ! The different data types used in the ocean modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================
TYPE type_ocean_model_transient
    ! Main data fields to compute the ocean model T and S transiently
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T0                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S0                           ! [PSU]             Salinity

    ! deltaT record
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: dT_series_time
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: dT_series 
    REAL(dp),                   ALLOCATABLE     :: dT_t0, dT_t1, dT_at_t0, dT_at_t1

  END TYPE type_ocean_model_transient

  TYPE type_ocean_model
    ! The ocean model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S                           ! [PSU]             Salinity

    ! Secondary data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_draft                     ! [degrees Celsius] Temperature at ice base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_freezing_point            ! [degrees Celsius] Pressure freezing point of water

    TYPE(type_ocean_model_transient)        :: transient

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_ocean_model

CONTAINS

END MODULE ocean_model_types