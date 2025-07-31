MODULE LMB_model_types

  ! The different data types used in the LMB modules

! ===== Preamble =====
! ====================

  USE precisions, ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================
  TYPE type_LMB_model_GlacialIndex
    ! Main data fields to compute the ocean model T and S transiently
    REAL(dp),                   ALLOCATABLE     :: LMB_warm                           ! [m.i.e./yr] Lateral mass balance
    REAL(dp),                   ALLOCATABLE     :: LMB_cold                           ! [m.i.e./yr] Lateral mass balance

    ! Glacial Index record
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: GI_series_time
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: GI_series
    REAL(dp),                   ALLOCATABLE     :: GI_t0, GI_t1, GI_at_t0, GI_at_t1


  END TYPE type_LMB_model_GlacialIndex

  TYPE type_LMB_model
    ! The LMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: LMB                         ! [m.i.e./yr] Lateral mass balance

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

    ! Transient GI-dependent LMB
    TYPE(type_LMB_model_GlacialIndex)       :: GI

  END TYPE type_LMB_model

CONTAINS

END MODULE LMB_model_types