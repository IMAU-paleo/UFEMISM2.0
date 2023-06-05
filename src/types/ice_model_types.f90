MODULE ice_model_types

  ! The different data types used in the ice modules

! ===== Preamble =====
! ====================

  USE petscksp                                               , ONLY: tMat
  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_ice_velocity_solver_SIA
    ! Data fields needed to solve the Shallow Ice Approximation

    ! Solution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D_b                      ! [m yr^-1] 3-D ice velocity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D                    ! [yr^-1] Vertical shear strain rates
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D

    ! Intermediate data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: D_3D_b                      ! [m yr^-1] Diffusivity

  END TYPE type_ice_velocity_solver_SIA

  TYPE type_ice_velocity_solver_SSA
    ! Data fields needed to solve the Shallow Shelf Approximation

    ! Solution
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_b                         ! [m yr^-1] 2-D horizontal ice velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_b

    ! Intermediate data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: A_flow_vav_a                ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: du_dx_a                     ! [yr^-1] 2-D horizontal strain rates
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: du_dy_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dv_dx_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dv_dy_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: eta_a                       ! Effective viscosity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: N_a                         ! Product term N = eta * H
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: N_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dN_dx_b                     ! Gradients of N
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dN_dy_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: beta_b_b                    ! Friction coefficient (tau_b = u * beta_b)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dx_b                    ! Driving stress
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dy_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_b_prev                    ! Velocity solution from previous viscosity iteration
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_b_prev

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised SSA
    REAL(dp)                                :: PETSc_rtol
    REAL(dp)                                :: PETSc_abstol

    ! Restart file
    CHARACTER(LEN=256)                      :: restart_filename

  END TYPE type_ice_velocity_solver_SSA

  TYPE type_ice_model
    ! The ice dynamics model data structure.

    ! === Masks ===
    ! =============

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
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_ice_prev               ! T: ice-covered during previous time step, F: otherwise
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: mask                        ! Combined mask with integers representing different parts, only used for quick visual inspection of output

    ! === Basic geometry ===
    ! ======================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi                          ! [m] Ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb                          ! [m] Bedrock elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                          ! [m] Surface elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hib                         ! [m] Ice base elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: TAF                         ! [m] Thickness above flotation

    ! === Geometry changes ===
    ! ========================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi                         ! [m] Ice thickness difference (w.r.t. to reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb                         ! [m] Bedrock elevation difference (w.r.t. to reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHs                         ! [m] Surface elevation difference (w.r.t. to reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHib                        ! [m] Base elevation difference (w.r.t. to reference)

    ! === Geometry rates of change ===
    ! ================================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt                      ! [m yr^-1] Ice thickness rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb_dt                      ! [m yr^-1] Bedrock elevation rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHs_dt                      ! [m yr^-1] Ice surface elevation rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHib_dt                     ! [m yr^-1] Ice base elevation rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_tplusdt                  ! [m] Predicted ice thickness at next time step

    ! === Terrain-following coordinate zeta gradients ===
    ! ===================================================

    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dx_3D                 ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dt_3D

    ! === Thermodynamics ===
    ! ======================

    ! Ice temperatures
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti                          ! [K] Englacial temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti_pmp                      ! [K] Pressure melting point temperature

    ! Physical quantities
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Cpi                         ! [J kg^-1 K^-1] Specific heat capacity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ki                          ! [J m^-1 K^-1 yr^-1] Thermal conductivity

    ! Heating
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: internal_heating            ! [?] Internal heating
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: frictional_heating          ! [?] Frictional heating

    ! === Ice flow ===
    ! ================

    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: A_flow_3D                   ! [Pa^-3 y^-1] Glen's flow law parameter
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: A_flow_vav                  ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter

    ! === Ice velocities ===
    ! ======================

    ! Velocity solvers
    TYPE(type_ice_velocity_solver_SIA)      :: SIA                         ! Shallow Ice Approximation
    TYPE(type_ice_velocity_solver_SSA)      :: SSA                         ! Shallow Shelf Approximation
    ! TYPE(type_ice_velocity_solver_DIVA)     :: DIVA                        ! Depth-Integrated Viscosity Approximation
    ! TYPE(type_ice_velocity_solver_BPA)      :: BPA                         ! Blatter-Pattyn Approximation

    ! 3-D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D                        ! [m yr^-1] 3-D ice velocity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: w_3D

    ! Vertically integrated
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav                       ! [m yr^-1] Vertically averaged ice velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_vav
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_vav_b

    ! Surface
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_surf                      ! [m yr^-1] Ice velocity at the surface
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_surf_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_surf_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: w_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_surf_b

    ! Basal
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_base                      ! [m yr^-1] Ice velocity at the base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_base_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_base_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: w_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_base_b

    ! === Strain rates ===
    ! ====================

    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dx_3D                    ! [yr^-1]
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dz_3D

    ! == Basal conditions ==
    ! ======================

    ! Basal hydrology
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: pore_water_pressure         ! Basal pore water pressure
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: overburden_pressure         ! Basal overburden pressure
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: effective_pressure          ! Basal effective pressure

    ! Basal roughness / friction
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: phi_fric                    ! Till friction angle (degrees)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_c                       ! Till yield stress tauc   (used when choice_sliding_law = "Coloumb" or "Coulomb_regularised")
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: alpha_sq                    ! Coulomb-law friction coefficient [unitless]         (used when choice_sliding_law =             "Tsai2015", or "Schoof2005")
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: beta_sq                     ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3] (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: beta_b                      ! Basal friction (tau_b = u * beta_b)

    ! Basal sliding
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: friction_coef_1             ! Generic basal friction coefficient 1
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: friction_coef_2             ! Generic basal friction coefficient 2
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: basal_shear_stress          ! Basal shear stress

    ! Geothermal heat
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: geothermal_heat_flux        ! Geothermal heat flux

    ! === Sea level ===
    ! =================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL                          ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)

    ! === Area fractions ===
    ! ======================

    ! Grounded
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: bedrock_cdf                 ! Sub-grid bedrock cumulative density function
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: fraction_gr                 ! [0-1] Grounded area fractions of vertices
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: fraction_gr_b               ! [0-1] Grounded area fractions of triangles

    ! Calving front
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: fraction_cf                 ! [0-1] Ice-covered area fractions of calving fronts

    ! === Extras ===
    ! ==============

    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: basin_ID                    ! The drainage basin to which each vertex belongs

  END TYPE type_ice_model

CONTAINS

END MODULE ice_model_types