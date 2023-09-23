MODULE LMB_model_types

  ! The different data types used in the LMB modules

! ===== Preamble =====
! ====================

  USE precisions, ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_LMB_model
    ! The LMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: LMB                         ! [m.i.e./yr] Lateral mass balance

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_LMB_model

CONTAINS

END MODULE LMB_model_types