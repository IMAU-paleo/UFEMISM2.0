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
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: basal_friction_coefficient_b! Basal friction coefficient (tau_b = u * beta_b)
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

  TYPE type_ice_velocity_solver_DIVA
    ! Data fields needed to solve the Depth-Integrated Viscosity Approximation

    ! Solution
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav_b                     ! [m yr^-1] 2-D horizontal ice velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_base_b                    ! [m yr^-1] 2-D horizontal ice velocity at the ice base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_base_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D_b                      ! [m yr^-1] 3-D horizontal ice velocity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D_b

    ! Intermediate data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: du_dx_a                     ! [yr^-1] 2-D horizontal strain rates
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: du_dy_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dv_dx_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dv_dy_a
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D_a                  ! [yr^-1] 3-D vertical shear strain rates
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D_a
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: eta_3D_a                    ! Effective viscosity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: eta_3D_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: eta_vav_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: N_a                         ! Product term N = eta * H
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: N_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dN_dx_b                     ! Gradients of N
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dN_dy_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: F1_3D_a                     ! F-integrals
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: F2_3D_a
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: F1_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: F2_3D_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: basal_friction_coefficient_b! Basal friction coefficient (tau_b = u * beta_b)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: beta_eff_a                  ! "Effective" friction coefficient (turning the SSA into the DIVA)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: beta_eff_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_bx_b                    ! Basal shear stress
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_by_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dx_b                    ! Driving stress
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dy_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_b_prev                    ! Velocity solution from previous viscosity iteration
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_b_prev

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised SSA
    REAL(dp)                                :: PETSc_rtol
    REAL(dp)                                :: PETSc_abstol

    ! Restart file
    CHARACTER(LEN=256)                      :: restart_filename

  END TYPE type_ice_velocity_solver_DIVA

  TYPE type_ice_velocity_solver_BPA
    ! Data fields needed to solve the Blatter-Pattyn Approximation

    ! Solution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_bk                        ! 3-D horizontal ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_bk

    ! Intermediate data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dx_ak                    ! Strain rates [yr^-1]
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dy_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dx_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dy_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dx_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dy_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dx_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dy_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: eta_ak                      ! Effective viscosity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: eta_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: eta_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: deta_dx_bk                  ! Gradients of eta
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: deta_dy_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: deta_dz_bk
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: basal_friction_coefficient_b! Friction coefficient (tau_b = u * beta_b)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dh_dx_b                     ! Surface slope
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dh_dy_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: db_dx_b                     ! Basal slope
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: db_dy_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dx_b                    ! Driving stress
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dy_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_bk_prev                   ! Previous velocity solution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_bk_prev

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised SSA
    REAL(dp)                                :: PETSc_rtol
    REAL(dp)                                :: PETSc_abstol

    ! Restart file
    CHARACTER(LEN=256)                      :: restart_filename

  END TYPE type_ice_velocity_solver_BPA

  TYPE type_ice_velocity_solver_hybrid
    ! Data fields needed to solve the hybrid DIVA/BPA

    ! Solution
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav_b                     ! Vertically averaged horizontal ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_bk                        ! 3-D horizontal ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_bk

    ! DIVA and BPA solvers
    TYPE(type_ice_velocity_solver_DIVA)     :: DIVA                        ! Depth-Integrated Viscosity Approximation
    TYPE(type_ice_velocity_solver_BPA)      :: BPA                         ! Blatter-Pattyn Approximation

    ! Solving masks
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_DIVA_b                 ! T: solve the DIVA here, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_BPA_b                  ! T: solve the BPA  here, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_3D_from_DIVA_b         ! T: calculate 3-D velocities from the vertically averaged DIVA solution here, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_vav_from_BPA_b         ! T: calculate vertically averaged velocities from the 3-D BPA  solution here, F: otherwise

    ! Intermediate data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_bk_prev                   ! Previous velocity solution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_bk_prev

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised SSA
    REAL(dp)                                :: PETSc_rtol
    REAL(dp)                                :: PETSc_abstol

    ! Restart file
    CHARACTER(LEN=256)                      :: restart_filename

  END TYPE type_ice_velocity_solver_hybrid

  TYPE type_ice_pc
    ! Data fields needed for the predictor/corrector time-stepping scheme

    REAL(dp)                                :: dt_n                         ! [yr]   Previous time step
    REAL(dp)                                :: dt_np1                       ! [yr]   Current  time step
    REAL(dp)                                :: zeta_t                       ! [-]    Ratio between previous and new time step
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt_Hi_nm1_u_nm1          ! [m/yr] Thinning rates from previous time step
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt_Hi_n_u_n              ! [m/yr] Thinning rates for current time step with old geometry
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_star_np1                  ! [m]    Predicted ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt_Hi_star_np1_u_np1     ! [m/yr] Thinning rates for predicted ice thickness and updated velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_np1                       ! [m]    Corrected ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_np1                      ! [m]    Truncation error
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: tau_n_guilty                 ! [-]    Number of PC iterations where vertex had truncation errors above the tolerance
    REAL(dp)                                :: eta_n                        ! [m]    Previous maximum truncation error
    REAL(dp)                                :: eta_np1                      ! [m]    Current  maximum truncation error

    ! Restart file
    CHARACTER(LEN=256)                      :: restart_filename

  END TYPE type_ice_pc

  TYPE type_ice_model
    ! The ice dynamics model data structure.

  ! === Ice-sheet geometry ===
  ! ==========================

    ! Basic geometry
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi                          ! [m] Ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb                          ! [m] Bedrock elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                          ! [m] Surface elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL                          ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hib                         ! [m] Ice base elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: TAF                         ! [m] Thickness above flotation
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_eff                      ! [m] Effective thickness of ice margins
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs_slope                    ! [-] Absolute surface gradient

    ! Geometry changes
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi                         ! [m] Ice thickness difference (w.r.t. reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb                         ! [m] Bedrock elevation difference (w.r.t. reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHs                         ! [m] Surface elevation difference (w.r.t. reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHib                        ! [m] Base elevation difference (w.r.t. reference)

    ! Rates of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt                      ! [m yr^-1] Ice thickness rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb_dt                      ! [m yr^-1] Bedrock elevation rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHs_dt                      ! [m yr^-1] Ice surface elevation rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHib_dt                     ! [m yr^-1] Ice base elevation rate of change
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt_raw                  ! [m yr^-1] Ice thickness rate of change before any ice thickness modifications
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt_residual             ! [m yr^-1] Residual ice thickness rate of change for inversions

    ! Horizontal derivatives
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHib_dx_b                   ! [] Horizontal derivative of ice draft on b-grid
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHib_dy_b                   ! [] Horizontal derivative of ice draft on b-grid

    ! Target quantities
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi_dt_target               ! [m yr^-1] Target ice thickness rate of change for inversions

    ! Masks
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_icefree_land           ! T: ice-free land , F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_icefree_ocean          ! T: ice-free ocean, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_grounded_ice           ! T: grounded ice  , F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_floating_ice           ! T: floating ice  , F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_icefree_land_prev      ! T: ice-free land , F: otherwise (during previous time step)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_icefree_ocean_prev     ! T: ice-free ocean, F: otherwise (during previous time step)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_grounded_ice_prev      ! T: grounded ice  , F: otherwise (during previous time step)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_floating_ice_prev      ! T: floating ice  , F: otherwise (during previous time step)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_margin                 ! T: ice next to ice-free, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_gr                  ! T: grounded ice next to floating ice, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_fl                  ! T: floating ice next to grounded ice, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_cf_gr                  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_cf_fl                  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_coastline              ! T: ice-free land next to ice-free ocean, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_noice                  ! T: no ice is allowed here, F: ice is allowed here
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_ROI                    ! T: located in ROI, F: otherwise (all false when no ROI is specified)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_SGD                    ! T: potential subglacial discharge area, F: otherwise
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: mask                        ! Diagnostic, only meant for quick visual inspection in output
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: basin_ID                    ! The drainage basin to which each vertex belongs

    ! Area fractions
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: fraction_gr                 ! [0-1] Grounded area fractions of vertices
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: fraction_gr_b               ! [0-1] Grounded area fractions of triangles
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: fraction_margin             ! [0-1] Ice-covered area fractions of ice margins

    ! Sub-grid bedrock cumulative density functions (CDFs)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: bedrock_cdf                 ! [-] Sub-grid bedrock cumulative density functions on the a-grid (vertices)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: bedrock_cdf_b               ! [-] Sub-grid bedrock cumulative density functions on the b-grid (triangles)

  ! === Terrain-following coordinate zeta gradients ===
  ! ===================================================

    ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta

    ! On the ak-grid (vertices, vertically regular)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dt_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dx_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dy_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dz_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dx2_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dxdy_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dy2_ak

    ! On the bk-grid (triangles, vertically regular)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dx_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dy_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dz_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dx2_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dxdy_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dy2_bk

    ! On the bks-grid (triangles, vertically staggered)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dx_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dy_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dzeta_dz_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dx2_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dxdy_bks
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: d2zeta_dy2_bks

  ! === Thermodynamics and rheology ===
  ! ===================================

    ! Ice temperatures
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti                          ! [K] Englacial temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti_pmp                      ! [K] Pressure melting point temperature
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Ti_hom                      ! [K] Basal temperature w.r.t. pressure melting point

    ! Physical quantities
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Cpi                         ! [J kg^-1 K^-1] Specific heat capacity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ki                          ! [J m^-1 K^-1 yr^-1] Thermal conductivity

    ! Heating
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: internal_heating            ! [?] Internal heating
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: frictional_heating          ! [?] Frictional heating

    ! Glen's flow law factor
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: A_flow                      ! [Pa^-3 y^-1] Glen's flow law factor

  ! === Ice velocities ===
  ! ======================

    ! Velocity solvers
    TYPE(type_ice_velocity_solver_SIA)      :: SIA                         ! Shallow Ice Approximation
    TYPE(type_ice_velocity_solver_SSA)      :: SSA                         ! Shallow Shelf Approximation
    TYPE(type_ice_velocity_solver_DIVA)     :: DIVA                        ! Depth-Integrated Viscosity Approximation
    TYPE(type_ice_velocity_solver_BPA)      :: BPA                         ! Blatter-Pattyn Approximation
    TYPE(type_ice_velocity_solver_hybrid)   :: hybrid                      ! Hybrid DIVA/BPA

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

  ! == Strain rates ==
  ! ==================

    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dx_3D                    ! [yr^-1]
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dz_3D

  ! == Ice flow regime ==
  ! =====================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: divQ                        ! [m yr^-1] Horizontal ice flux divergence
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: R_shear                     ! [0-1]     uabs_base / uabs_surf (0 = pure vertical shear, viscous flow; 1 = pure sliding, plug flow)

  ! == Basal hydrology ==
  ! =====================

    ! Basal hydrology
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: pore_water_pressure         ! [Pa]  Basal pore water pressure
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: overburden_pressure         ! [Pa]  Basal overburden pressure
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: effective_pressure          ! [Pa]  Basal effective pressure
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: pore_water_likelihood       ! [0-1] Basal pore water likelihood
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: pore_water_fraction         ! [0-1] Fraction of overburden pressure reduced by pore water pressure

  ! == Basal sliding ==
  ! ===================

    ! Basal friction and shear stress
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: till_yield_stress           ! [Pa]               Till yield stress (used when choice_sliding_law = "Coloumb", "Budd", or "Zoet-Iverson")
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: basal_friction_coefficient  ! [Pa yr m^-1]       Basal friction coefficient (basal_shear_stress = u_base * basal_friction_coefficient)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: basal_shear_stress          ! [Pa]               Basal shear stress

  ! == Geothermal heat ==
  ! =====================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: geothermal_heat_flux        ! [J m^-2 yr^-1] Geothermal heat flux

  ! === Ice thickness time stepping ===
  ! ===================================

    ! Ice thickness time-stepping solvers
    TYPE(type_ice_pc)                       :: pc

    ! Time frames and ice thicknesses
    REAL(dp)                                :: t_Hi_prev                   ! [yr] Time of the previous state
    REAL(dp)                                :: t_Hi_next                   ! [yr] Time of the next state
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_prev                     ! [m]  The previous state
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_next                     ! [m]  The next state

  ! === Ice temperature time stepping ===
  ! =====================================

    ! Time frames and ice temperatures
    REAL(dp)                                :: t_Ti_prev                   ! [yr] Time of the previous state
    REAL(dp)                                :: t_Ti_next                   ! [yr] Time of the next state
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti_prev                     ! [m]  The previous state
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti_next                     ! [m]  The next state
    CHARACTER(LEN=256)                      :: thermo_restart_filename     !      Filename for thermodynamics restart file

  ! === Numerical stability info ===
  ! ================================

    real(dp)                                :: dt_ice                      ! [yr] Ice-dynamical time step
    integer                                 :: n_visc_its                  !      Number of non-linear viscosity iterations
    integer                                 :: n_Axb_its                   !      Number of iterations in iterative solver for linearised momentum balance

  END TYPE type_ice_model

CONTAINS

END MODULE ice_model_types
