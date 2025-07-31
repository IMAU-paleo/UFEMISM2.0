MODULE AMB_model_types

  ! The different data types used in the AMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_AMB_model
    ! The AMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: AMB                         ! [m.i.e./yr] Basal mass balance

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_AMB_model

CONTAINS

END MODULE AMB_model_types