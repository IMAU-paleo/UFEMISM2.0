MODULE SMB_model_types

  ! The different data types used in the SMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_SMB_model
    ! The SMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SMB                         ! [m.i.e./yr] Surface mass balance

    ! Sub-models
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SMB_correction              ! [m.i.e./yr] Surface mass balance

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_SMB_model

CONTAINS

END MODULE SMB_model_types