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
  
 ! ==== Type for ELRA model ==== 
 ! =============================

  TYPE type_ELRA_model
    ! The ELRA model data structure.

    ! The ELRA model grid - I think this is not needed, ELRA use the same square grid from GIA
    !TYPE(type_grid)                         :: grid 
    
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: surface_load_mesh                ! [Pa] Total surface load applied at each mesh point
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb_grid_partial                 ! [m]  Partial bedrock displacement on the grid
    REAL(dp), DIMENSION(: , :), ALLOCATABLE :: dHb_grid_tot                     ! [m]  Total bedrock displacement on the grid
    REAL(dp), DIMENSION(: , :), ALLOCATABLE :: dHb_eq_grid                      ! [m]  Bedrock displacement at equilibrium (computed from relative surface load)
    REAL(dp), DIMENSION(: , :), ALLOCATABLE :: dHb_dt_grid                      ! [m/yr] Bedrock deformation rate on the grid
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb_dt_grid_partial              ! [m/yr] Partial bedrock deformation rate on the grid
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb_dt_mesh                      ! [m/yr] Bedrock deformation rate mapped on the mesh
    REAL(dp), DIMENSION(: , :), ALLOCATABLE :: relative_surface_load_grid_tot   ! [Pa] Relative surface load (difference from equilibrium) on the full grid
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: surface_load_GIAeq               ! [Pa] Reference surface load used for GIA equilibrium computation
    integer                                 :: flex_prof_rad                    ! radius of the flexoral profile
    REAL(dp), DIMENSION(: , :), ALLOCATABLE :: flex_prof_grid                   ! radius of the flexoral profile on GIA grid    
    
  END TYPE type_ELRA_model
 
CONTAINS

END MODULE GIA_model_types