MODULE climate_types

  ! The different data types used in the climate module

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_climate_model
    ! The climate model data structure.

    ! Main fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T2m                         ! [K]      Monthly 2-m air temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Precp                       ! [m.w.e.] Monthly precipitation

  END TYPE type_climate_model

CONTAINS

END MODULE climate_types