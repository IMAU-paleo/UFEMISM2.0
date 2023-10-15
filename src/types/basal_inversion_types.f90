MODULE basal_inversion_types

  ! The different data types used in the basal inversion modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  ! == Bed roughness
  ! ================

  TYPE type_basal_inversion_H_dHdt_flowline
    ! The basal inversion model based on flowline-averaged values of H and dH/dt

    ! Main data fields

    ! Sub-models

    ! Timestepping

  END TYPE type_basal_inversion_H_dHdt_flowline

  TYPE type_basal_inversion
    ! The main basal roughness inversion model

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_1
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_2

    ! Timestepping
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_1_prev
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_2_prev
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_1_next
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: generic_bed_roughness_2_next
    REAL(dp)                                :: t_prev, t_next

  END TYPE type_basal_inversion

  ! == Pore water pressure
  ! ======================

  TYPE type_hydrology_inversion
    ! The main basal hydrology inversion model

    ! Timestepping
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: pore_water_fraction_prev
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: pore_water_fraction_next
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: pore_water_fraction_app
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_inverted_point
    REAL(dp)                                :: t_prev, t_next

  END TYPE type_hydrology_inversion

CONTAINS

END MODULE basal_inversion_types