MODULE basal_inversion_types

  ! The different data types used in the basal inversion modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_basal_inversion_H_dHdt_flowline
    ! The basal inversion model based on flowline-averaged values of H and dH/dt

    ! Main data fields

    ! Sub-models

    ! Timestepping

  END TYPE type_basal_inversion_H_dHdt_flowline

  TYPE type_basal_inversion
    ! The main basal inversion model

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_1
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_2

    ! Sub-models
    TYPE(type_basal_inversion_H_dHdt_flowline) :: flowline_H_dHdt

    ! Timestepping
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_1_prev
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_2_prev
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_1_next
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_2_next
    REAL(dp)                                :: t_prev, t_next

  END TYPE type_basal_inversion

CONTAINS

END MODULE basal_inversion_types