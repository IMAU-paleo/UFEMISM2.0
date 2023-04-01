MODULE ice_model_types

  ! The different data types used in the ice modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE ice_configuration                                      , ONLY: type_CFG_ice

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  TYPE type_ice_velocity_solver_SIA
    ! Data fields needed to solve the Shallow Ice Approximation

    ! Solution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D_b                      ! [m yr^-1] 3-D ice velocity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D_a                  ! [yr^-1] Vertical shear strain rates
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D_a

    ! Intermediate data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: D_3D_b                      ! [m yr^-1] Diffusivity

  END TYPE type_ice_velocity_solver_SIA

  TYPE type_ice_model
    ! The ice dynamics model data structure.

    ! Configuration parameters
    TYPE(type_CFG_ice)                      :: C

    ! Basic ice geometry - ice thickness, bedrock elevation, surface elevation, ice base elevation, sea level elevation, thickness above floatation
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi                          ! [m] Ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb                          ! [m] Bedrock elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hib                         ! [m] Ice base elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                          ! [m] Surface elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL                          ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: TAF                         ! [m] Thickness above flotation

    ! Different masks
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_land                   ! T: land above water level (sea and/or lake), F: land below water level (idem)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_ocean                  ! T: land below sea level, F: land above sea level
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_lake                   ! T: land below lake level, F: land above lake level
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_ice                    ! T: Hi > 0, F: Hi = 0
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_sheet                  ! T: grounded ice, F: floating ice or ice-free
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_shelf                  ! T: floating ice, F: grounded ice or ice-free
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_coast                  ! T: land next to sea/lake, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_margin                 ! T: ice next to ice-free, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_gr                  ! T: grounded ice next to floating ice, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_fl                  ! T: floating ice next to grounded ice, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_cf_gr                  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_cf_fl                  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: mask                        ! Combined mask with integers representing different parts, only used for quick visual inspection of output
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: f_grnd                      ! [0-1] Grounded fractions of vertices
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: f_grnd_b                    ! [0-1] Grounded fractions of triangles
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: basin_ID                    ! The drainage basin to which each vertex belongs

    ! Different velocity solvers
    TYPE(type_ice_velocity_solver_SIA)      :: SIA
!    TYPE(type_ice_velocity_solver_SSA)      :: SSA
!    TYPE(type_ice_velocity_solver_DIVA)     :: DIVA
!    TYPE(type_ice_velocity_solver_BPA)      :: BPA

    ! Ice velocities
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D                        ! [m yr^-1] 3-D ice velocity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: w_3D

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav                       ! [m yr^-1] Vertically averaged ice velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_vav
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_vav_b

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_surf                      ! [m yr^-1] Ice velocity at the surface
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_surf_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_surf_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: w_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_surf_b

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_base                      ! [m yr^-1] Ice velocity at the base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_base_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_base_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: w_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_base_b

    ! Strain rates
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dx_3D                    ! [yr^-1]
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dz_3D

    ! Ice thermodynamical/rheological properties
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti                          ! [K]                 Englacial temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: A_flow_3D                   ! [Pa^-3 y^-1]        Glen's flow law parameter
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti_pmp                      ! [K]                 Pressure melting point temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Cpi                         ! [J kg^-1 K^-1]      Specific heat capacity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ki                          ! [J m^-1 K^-1 yr^-1] Thermal conductivity

  END TYPE type_ice_model

CONTAINS

! ===== Subroutines ======
! ========================

END MODULE ice_model_types