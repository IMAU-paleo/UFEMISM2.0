MODULE BMB_model_types

  ! The different data types used in the BMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_BMB_model
    ! The BMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB                         ! [m.i.e./yr] Basal mass balance

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_BMB_model

CONTAINS

END MODULE BMB_model_types