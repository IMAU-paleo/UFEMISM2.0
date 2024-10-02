MODULE BMB_model_types

  ! The different data types used in the BMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE laddie_model_types                                     , ONLY: type_laddie_model

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_BMB_model
    ! The BMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB                         ! [m.i.e./yr] Basal mass balance
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_shelf                   ! [m.i.e./yr] Basal mass balance of floating ice, extrapolated below grounded ice
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_inv                     ! [m.i.e./yr] Inverted basal mass balance
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_ref                     ! [m.i.e./yr] Reference basal mass balance (e.g. inverted or prescribed)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_floating_ice           ! T: floating vertex, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_fl                  ! T: floating grounding line vertex, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_gr                  ! T: grounded grounding line vertex, F: otherwise

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

    ! LADDIE
    TYPE(type_laddie_model)                       :: laddie

  END TYPE type_BMB_model

CONTAINS

END MODULE BMB_model_types
