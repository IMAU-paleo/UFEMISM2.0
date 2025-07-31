MODULE BMB_model_types

  ! The different data types used in the BMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE laddie_model_types                                     , ONLY: type_laddie_model
  use laddie_forcing_types, only: type_laddie_forcing
  use reference_geometry_types, only: type_reference_geometry

  IMPLICIT NONE

  private

  public :: type_BMB_model, type_BMB_model_inverted

! ===== Types =====
! =================

  type type_BMB_model_inverted

    real(dp), dimension(:), allocatable :: BMB                             !< [m.i.e./yr] Basal mass balance
    type(type_reference_geometry)       :: target_geometry                 !< The geometry that the BMB inversion should aim to reproduce
    logical,  dimension(:), allocatable :: target_mask_shelf               !< Shelf mask of the target geometry

  end type type_BMB_model_inverted

  TYPE type_BMB_model
    ! The BMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB                         ! [m.i.e./yr] Basal mass balance
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_shelf                   ! [m.i.e./yr] Basal mass balance of floating ice, extrapolated below grounded ice
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_inv                     ! [m.i.e./yr] Inverted basal mass balance
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_ref                     ! [m.i.e./yr] Reference basal mass balance (e.g. inverted or prescribed)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_transition_phase        ! [m.i.e./yr] Basal mass balance transition phase, weighted average between inversion and modelled BMB
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: BMB_modelled                ! [m.i.e./yr] Basal mass balance modelled (needs to be saved for computation BMB_transition_phase)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_floating_ice           ! T: floating vertex, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_fl                  ! T: floating grounding line vertex, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_gr                  ! T: grounded grounding line vertex, F: otherwise

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

    ! LADDIE
    type(type_laddie_model)                       :: laddie
    type(type_laddie_forcing)                     :: forcing

    type(type_BMB_model_inverted) :: inv

  END TYPE type_BMB_model

CONTAINS

END MODULE BMB_model_types
