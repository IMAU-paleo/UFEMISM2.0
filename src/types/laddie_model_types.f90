MODULE laddie_model_types

  ! The different data types used in the laddie modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

TYPE type_laddie_model
  ! The laddie model structure

  ! Main data fields
  REAL(dp), DIMENSION(:    ), ALLOCATABLE :: H                           ! [m] Layer thickness
  REAL(dp), DIMENSION(:    ), ALLOCATABLE :: U                           ! [m s^-1] 2D velocity
  REAL(dp), DIMENSION(:    ), ALLOCATABLE :: V                           ! [m s^-1]
  REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T                           ! [degrees Celsius] Temperature
  REAL(dp), DIMENSION(:    ), ALLOCATABLE :: S                           ! [PSU]             Salinity  

  

END TYPE type_laddie_model


CONTAINS

END MODULE laddie_model_types