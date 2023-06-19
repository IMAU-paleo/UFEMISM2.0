MODULE SMB_types

  ! The different data types used in the SMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_SMB_model
    ! The SMB model data structure.

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SMB                         ! [m.i.e./yr] Surface mass balance

  END TYPE type_SMB_model

CONTAINS

END MODULE SMB_types