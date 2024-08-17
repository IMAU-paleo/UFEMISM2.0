MODULE GIA_model_types

  ! The different data types used in the GIA modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE grid_types                                             , ONLY: type_grid

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_GIA_model
    ! The GIA model data structure.

    ! The GIA model grid
    TYPE(type_grid)                         :: grid

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: relative_surface_load_mesh  ! [Pa] Surface load relative to the GIA-equilibrium reference geometry, on the UFEMISM mesh
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: relative_surface_load_grid  ! [Pa] Surface load relative to the GIA-equilibrium reference geometry, on the UFEMISM mesh

    ! Sub-models

    ! Timestepping
    REAL(dp)                                :: t_prev                      ! [yr] Time of the previous state
    REAL(dp)                                :: t_next                      ! [yr] Time of the next state
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb_prev                    ! [m]  The previous state (bedrock deflection relative to the GIA-equilibrium reference geometry)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb_next                    ! [m]  The next state     (idem)

  END TYPE type_GIA_model

CONTAINS

END MODULE GIA_model_types