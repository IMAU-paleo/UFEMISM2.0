MODULE BMB_types

  ! The different data types used in the BMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_BMB_model
    ! The BMB model data structure.

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB                         ! [m.i.e./yr] Basal mass balance

  END TYPE type_BMB_model

CONTAINS

END MODULE BMB_types