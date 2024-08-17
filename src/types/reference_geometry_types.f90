MODULE reference_geometry_types

  ! The reference geometry data type

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE grid_types                                             , ONLY: type_grid
  USE mesh_types                                             , ONLY: type_mesh

  IMPLICIT NONE

! ===== Types =====
! =================

  ! Data structure containing a reference ice-sheet geometry
  TYPE type_reference_geometry

    ! Data on the raw grid/mesh
    TYPE(type_grid)                         :: grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb_grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs_grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL_grid_raw
    TYPE(type_mesh)                         :: mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb_mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs_mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL_mesh_raw

    ! Data on the model mesh
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi                          ! [m] Ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb                          ! [m] Bedrock elevation [w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                          ! [m] Surface elevation [w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL                          ! [m] Sea surface elevation [w.r.t. PD sea level]

  END TYPE type_reference_geometry

CONTAINS

END MODULE reference_geometry_types