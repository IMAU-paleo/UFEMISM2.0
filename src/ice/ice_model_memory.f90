MODULE ice_model_memory

  ! Routines for administrating the memory for the ice model data.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE control_resources_and_error_messaging                  , ONLY: init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE allocate_ice_model( mesh, ice)
    ! Allocate ice model variables

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),      INTENT(IN)    :: mesh
    TYPE(type_ice_model), INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER       :: routine_name = 'allocate_ice_model'

    ! Add routine to path
    CALL init_routine( routine_name)

  ! === Ice-sheet geometry ===
  ! ==========================

    ! Basic geometry
    ALLOCATE( ice%Hi                          ( mesh%vi1:mesh%vi2        ))  ! [m] Ice thickness
    ALLOCATE( ice%Hb                          ( mesh%vi1:mesh%vi2        ))  ! [m] Bedrock elevation (w.r.t. PD sea level)
    ALLOCATE( ice%Hs                          ( mesh%vi1:mesh%vi2        ))  ! [m] Surface elevation (w.r.t. PD sea level)
    ALLOCATE( ice%SL                          ( mesh%vi1:mesh%vi2        ))  ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    ALLOCATE( ice%Hib                         ( mesh%vi1:mesh%vi2        ))  ! [m] Ice base elevation (w.r.t. PD sea level)
    ALLOCATE( ice%TAF                         ( mesh%vi1:mesh%vi2        ))  ! [m] Thickness above flotation
    ALLOCATE( ice%Hi_eff                      ( mesh%vi1:mesh%vi2        ))  ! [m] Thickness above flotation

    ice%Hi                          = 0._dp
    ice%Hb                          = 0._dp
    ice%Hs                          = 0._dp
    ice%SL                          = 0._dp
    ice%Hib                         = 0._dp
    ice%TAF                         = 0._dp
    ice%Hi_eff                      = 0._dp

    ! Geometry changes
    ALLOCATE( ice%dHi                         ( mesh%vi1:mesh%vi2        ))  ! [m] Ice thickness difference (w.r.t. to reference)
    ALLOCATE( ice%dHb                         ( mesh%vi1:mesh%vi2        ))  ! [m] Bedrock elevation difference (w.r.t. to reference)
    ALLOCATE( ice%dHs                         ( mesh%vi1:mesh%vi2        ))  ! [m] Surface elevation difference (w.r.t. to reference)
    ALLOCATE( ice%dHib                        ( mesh%vi1:mesh%vi2        ))  ! [m] Base elevation difference (w.r.t. to reference)

    ice%dHi                         = 0._dp
    ice%dHb                         = 0._dp
    ice%dHs                         = 0._dp
    ice%dHib                        = 0._dp

    ! Rates of change
    ALLOCATE( ice%dHi_dt                      ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Ice thickness rate of change
    ALLOCATE( ice%dHb_dt                      ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Bedrock elevation rate of change
    ALLOCATE( ice%dHs_dt                      ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Ice surface elevation rate of change
    ALLOCATE( ice%dHib_dt                     ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Ice base elevation rate of change
    ALLOCATE( ice%dHi_dt_raw                  ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Ice thickness rate of change before any modifications
    ALLOCATE( ice%dHi_dt_target               ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Target ice thickness rate of change for inversions
    ALLOCATE( ice%dHi_dt_residual             ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Residual ice thickness rate of change for inversions

    ice%dHi_dt                      = 0._dp
    ice%dHb_dt                      = 0._dp
    ice%dHs_dt                      = 0._dp
    ice%dHib_dt                     = 0._dp
    ice%dHi_dt_raw                  = 0._dp
    ice%dHi_dt_target               = 0._dp
    ice%dHi_dt_residual             = 0._dp

    ! Masks
    ALLOCATE( ice%mask_icefree_land           ( mesh%vi1:mesh%vi2        ))  ! T: ice-free land , F: otherwise
    ALLOCATE( ice%mask_icefree_ocean          ( mesh%vi1:mesh%vi2        ))  ! T: ice-free ocean, F: otherwise
    ALLOCATE( ice%mask_grounded_ice           ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice  , F: otherwise
    ALLOCATE( ice%mask_floating_ice           ( mesh%vi1:mesh%vi2        ))  ! T: floating ice  , F: otherwise
    ALLOCATE( ice%mask_icefree_land_prev      ( mesh%vi1:mesh%vi2        ))  ! T: ice-free land , F: otherwise (during previous time step)
    ALLOCATE( ice%mask_icefree_ocean_prev     ( mesh%vi1:mesh%vi2        ))  ! T: ice-free ocean, F: otherwise (during previous time step)
    ALLOCATE( ice%mask_grounded_ice_prev      ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice  , F: otherwise (during previous time step)
    ALLOCATE( ice%mask_floating_ice_prev      ( mesh%vi1:mesh%vi2        ))  ! T: floating ice  , F: otherwise (during previous time step)
    ALLOCATE( ice%mask_margin                 ( mesh%vi1:mesh%vi2        ))  ! T: ice next to ice-free, F: otherwise
    ALLOCATE( ice%mask_gl_gr                  ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice next to floating ice, F: otherwise
    ALLOCATE( ice%mask_gl_fl                  ( mesh%vi1:mesh%vi2        ))  ! T: floating ice next to grounded ice, F: otherwise
    ALLOCATE( ice%mask_cf_gr                  ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    ALLOCATE( ice%mask_cf_fl                  ( mesh%vi1:mesh%vi2        ))  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    ALLOCATE( ice%mask_coastline              ( mesh%vi1:mesh%vi2        ))  ! T: ice-free land next to ice-free ocean, F: otherwise
    ALLOCATE( ice%mask_noice                  ( mesh%vi1:mesh%vi2        ))  ! T: no ice is allowed here, F: ice is allowed here
    ALLOCATE( ice%mask                        ( mesh%vi1:mesh%vi2        ))  ! Diagnostic, only meant for quick visual inspection in output
    ALLOCATE( ice%basin_ID                    ( mesh%vi1:mesh%vi2        ))  ! The drainage basin to which each vertex belongs

    ice%mask_icefree_land           = .FALSE.
    ice%mask_icefree_ocean          = .FALSE.
    ice%mask_grounded_ice           = .FALSE.
    ice%mask_floating_ice           = .FALSE.
    ice%mask_icefree_land_prev      = .FALSE.
    ice%mask_icefree_ocean_prev     = .FALSE.
    ice%mask_grounded_ice_prev      = .FALSE.
    ice%mask_floating_ice_prev      = .FALSE.
    ice%mask_margin                 = .FALSE.
    ice%mask_gl_gr                  = .FALSE.
    ice%mask_gl_fl                  = .FALSE.
    ice%mask_cf_gr                  = .FALSE.
    ice%mask_cf_fl                  = .FALSE.
    ice%mask_coastline              = .FALSE.
    ice%mask_noice                  = .FALSE.
    ice%mask                        = 0
    ice%basin_ID                    = 0

    ! Area fractions
    ALLOCATE( ice%fraction_gr                 ( mesh%vi1:mesh%vi2        ))  ! [0-1] Grounded area fractions of vertices
    ALLOCATE( ice%fraction_gr_b               ( mesh%ti1:mesh%ti2        ))  ! [0-1] Grounded area fractions of triangles
    ALLOCATE( ice%fraction_margin             ( mesh%vi1:mesh%vi2        ))  ! [0-1] Ice-covered area fractions of calving fronts

    ice%fraction_gr                 = 0._dp
    ice%fraction_gr_b               = 0._dp
    ice%fraction_margin             = 0._dp

    ! Sub-grid bedrock cumulative density functions (CDFs)
    ALLOCATE( ice%bedrock_cdf                 ( mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins))  ! [-] Sub-grid bedrock cumulative density functions on the a-grid (vertices)
    ALLOCATE( ice%bedrock_cdf_b               ( mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins))  ! [-] Sub-grid bedrock cumulative density functions on the b-grid (triangles)

    ice%bedrock_cdf                 = 0._dp
    ice%bedrock_cdf_b               = 0._dp

  ! === Terrain-following coordinate zeta gradients ===
  ! ===================================================

    ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta

    ! On the ak-grid (vertices, vertically regular)
    ALLOCATE( ice%dzeta_dt_ak                 ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dzeta_dx_ak                 ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dzeta_dy_ak                 ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dzeta_dz_ak                 ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%d2zeta_dx2_ak               ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%d2zeta_dxdy_ak              ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%d2zeta_dy2_ak               ( mesh%vi1:mesh%vi2,mesh%nz))

    ice%dzeta_dt_ak                 = 0._dp
    ice%dzeta_dx_ak                 = 0._dp
    ice%dzeta_dy_ak                 = 0._dp
    ice%dzeta_dz_ak                 = 0._dp
    ice%d2zeta_dx2_ak               = 0._dp
    ice%d2zeta_dxdy_ak              = 0._dp
    ice%d2zeta_dy2_ak               = 0._dp

    ! On the bk-grid (triangles, vertically regular)
    ALLOCATE( ice%dzeta_dx_bk                 ( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%dzeta_dy_bk                 ( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%dzeta_dz_bk                 ( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%d2zeta_dx2_bk               ( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%d2zeta_dxdy_bk              ( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%d2zeta_dy2_bk               ( mesh%ti1:mesh%ti2,mesh%nz))

    ice%dzeta_dx_bk                 = 0._dp
    ice%dzeta_dy_bk                 = 0._dp
    ice%dzeta_dz_bk                 = 0._dp
    ice%d2zeta_dx2_bk               = 0._dp
    ice%d2zeta_dxdy_bk              = 0._dp
    ice%d2zeta_dy2_bk               = 0._dp

    ! On the bks-grid (triangles, vertically staggered)
    ALLOCATE( ice%dzeta_dx_bks                ( mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%dzeta_dy_bks                ( mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%dzeta_dz_bks                ( mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%d2zeta_dx2_bks              ( mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%d2zeta_dxdy_bks             ( mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%d2zeta_dy2_bks              ( mesh%ti1:mesh%ti2,mesh%nz-1))

    ice%dzeta_dx_bks                = 0._dp
    ice%dzeta_dy_bks                = 0._dp
    ice%dzeta_dz_bks                = 0._dp
    ice%d2zeta_dx2_bks              = 0._dp
    ice%d2zeta_dxdy_bks             = 0._dp
    ice%d2zeta_dy2_bks              = 0._dp

  ! === Thermodynamics and rheology ===
  ! ===================================

    ! Ice temperatures
    ALLOCATE( ice%Ti                          ( mesh%vi1:mesh%vi2,mesh%nz))  ! [K] Englacial temperature
    ALLOCATE( ice%Ti_pmp                      ( mesh%vi1:mesh%vi2,mesh%nz))  ! [K] Pressure melting point temperature
    ALLOCATE( ice%Ti_hom                      ( mesh%vi1:mesh%vi2        ))  ! [K] Basal temperature w.r.t. pressure melting point

    ice%Ti                          = 0._dp
    ice%Ti_pmp                      = 0._dp
    ice%Ti_hom                      = 0._dp

    ! Physical quantities
    ALLOCATE( ice%Cpi                         ( mesh%vi1:mesh%vi2,mesh%nz))  ! [J kg^-1 K^-1] Specific heat capacity
    ALLOCATE( ice%Ki                          ( mesh%vi1:mesh%vi2,mesh%nz))  ! [J m^-1 K^-1 yr^-1] Thermal conductivity

    ice%Cpi                         = 0._dp
    ice%Ki                          = 0._dp

    ! Heating
    ALLOCATE( ice%internal_heating            ( mesh%vi1:mesh%vi2,mesh%nz))  ! [?] Internal heating
    ALLOCATE( ice%frictional_heating          ( mesh%vi1:mesh%vi2        ))  ! [?] Frictional heating

    ice%internal_heating            = 0._dp
    ice%frictional_heating          = 0._dp

    ! Glen's flow law factor
    ALLOCATE( ice%A_flow                      ( mesh%vi1:mesh%vi2,mesh%nz))  ! [Pa^-3 y^-1] Glen's flow law factor

    ice%A_flow                      = 0._dp

  ! === Ice velocities ===
  ! ======================

    ! 3-D
    ALLOCATE( ice%u_3D                        ( mesh%vi1:mesh%vi2,mesh%nz))  ! [m yr^-1] 3-D ice velocity
    ALLOCATE( ice%v_3D                        ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%u_3D_b                      ( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%v_3D_b                      ( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%w_3D                        ( mesh%vi1:mesh%vi2,mesh%nz))

    ice%u_3D                        = 0._dp
    ice%v_3D                        = 0._dp
    ice%u_3D_b                      = 0._dp
    ice%v_3D_b                      = 0._dp
    ice%w_3D                        = 0._dp

    ! Vertically integrated
    ALLOCATE( ice%u_vav                       ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Vertically averaged ice velocity
    ALLOCATE( ice%v_vav                       ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%u_vav_b                     ( mesh%ti1:mesh%ti2        ))
    ALLOCATE( ice%v_vav_b                     ( mesh%ti1:mesh%ti2        ))
    ALLOCATE( ice%uabs_vav                    ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%uabs_vav_b                  ( mesh%ti1:mesh%ti2        ))

    ice%u_vav                       = 0._dp
    ice%v_vav                       = 0._dp
    ice%u_vav_b                     = 0._dp
    ice%v_vav_b                     = 0._dp
    ice%uabs_vav                    = 0._dp
    ice%uabs_vav_b                  = 0._dp

    ! Surface
    ALLOCATE( ice%u_surf                      ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Ice velocity at the surface
    ALLOCATE( ice%v_surf                      ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%u_surf_b                    ( mesh%ti1:mesh%ti2        ))
    ALLOCATE( ice%v_surf_b                    ( mesh%ti1:mesh%ti2        ))
    ALLOCATE( ice%w_surf                      ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%uabs_surf                   ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%uabs_surf_b                 ( mesh%ti1:mesh%ti2        ))

    ice%u_surf                      = 0._dp
    ice%v_surf                      = 0._dp
    ice%u_surf_b                    = 0._dp
    ice%v_surf_b                    = 0._dp
    ice%w_surf                      = 0._dp
    ice%uabs_surf                   = 0._dp
    ice%uabs_surf_b                 = 0._dp

    ! Basal
    ALLOCATE( ice%u_base                      ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Ice velocity at the base
    ALLOCATE( ice%v_base                      ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%u_base_b                    ( mesh%ti1:mesh%ti2        ))
    ALLOCATE( ice%v_base_b                    ( mesh%ti1:mesh%ti2        ))
    ALLOCATE( ice%w_base                      ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%uabs_base                   ( mesh%vi1:mesh%vi2        ))
    ALLOCATE( ice%uabs_base_b                 ( mesh%ti1:mesh%ti2        ))

    ice%u_base                      = 0._dp
    ice%v_base                      = 0._dp
    ice%u_base_b                    = 0._dp
    ice%v_base_b                    = 0._dp
    ice%w_base                      = 0._dp
    ice%uabs_base                   = 0._dp
    ice%uabs_base_b                 = 0._dp

  ! == Strain rates ==
  ! ==================

    ALLOCATE( ice%du_dx_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))  ! [yr^-1]
    ALLOCATE( ice%du_dy_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%du_dz_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dv_dx_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dv_dy_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dv_dz_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dw_dx_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dw_dy_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dw_dz_3D                    ( mesh%vi1:mesh%vi2,mesh%nz))

    ice%du_dx_3D                    = 0._dp
    ice%du_dy_3D                    = 0._dp
    ice%du_dz_3D                    = 0._dp
    ice%dv_dx_3D                    = 0._dp
    ice%dv_dy_3D                    = 0._dp
    ice%dv_dz_3D                    = 0._dp
    ice%dw_dx_3D                    = 0._dp
    ice%dw_dy_3D                    = 0._dp
    ice%dw_dz_3D                    = 0._dp

  ! == Ice flow regime ==
  ! =====================

    ALLOCATE( ice%divQ                        ( mesh%vi1:mesh%vi2        ))  ! [m yr^-1] Horizontal ice flux divergence
    ALLOCATE( ice%R_shear                     ( mesh%vi1:mesh%vi2        ))  ! [0-1]     uabs_base / uabs_surf (0 = pure vertical shear, viscous flow; 1 = pure sliding, plug flow)

    ice%divQ                        = 0._dp
    ice%R_shear                     = 0._dp

  ! == Basal hydrology ==
  ! =====================

    ! Basal hydrology
    ALLOCATE( ice%pore_water_pressure         ( mesh%vi1:mesh%vi2        ))  ! [Pa]  Basal pore water pressure
    ALLOCATE( ice%overburden_pressure         ( mesh%vi1:mesh%vi2        ))  ! [Pa]  Basal overburden pressure
    ALLOCATE( ice%effective_pressure          ( mesh%vi1:mesh%vi2        ))  ! [Pa]  Basal effective pressure
    ALLOCATE( ice%pore_water_likelihood       ( mesh%vi1:mesh%vi2        ))  ! [0-1] Basal pore water likelihood
    ALLOCATE( ice%pore_water_fraction         ( mesh%vi1:mesh%vi2        ))  ! [0-1] Fraction of overburden pressure reduced by pore water pressure

    ice%pore_water_pressure         = 0._dp
    ice%overburden_pressure         = 0._dp
    ice%effective_pressure          = 0._dp
    ice%pore_water_likelihood       = 0._dp
    ice%pore_water_fraction         = 0._dp

  ! == Basal sliding ==
  ! ===================

    ! Sliding law coefficients
    ALLOCATE( ice%till_friction_angle         ( mesh%vi1:mesh%vi2        ))  ! [degrees]          Till friction angle (degrees)
    ALLOCATE( ice%till_yield_stress           ( mesh%vi1:mesh%vi2        ))  ! [Pa]               Till yield stress (used when choice_sliding_law = "Coloumb", "Budd", or "Zoet-Iverson")
    ALLOCATE( ice%slid_alpha_sq               ( mesh%vi1:mesh%vi2        ))  ! [-]                Coulomb-law friction coefficient (used when choice_sliding_law = "Tsai2015", or "Schoof2005")
    ALLOCATE( ice%slid_beta_sq                ( mesh%vi1:mesh%vi2        ))  ! [Pa m^−1/m yr^1/m] Power-law friction coefficient (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")

    ice%till_friction_angle         = 0._dp
    ice%till_yield_stress           = 0._dp
    ice%slid_alpha_sq               = 0._dp
    ice%slid_beta_sq                = 0._dp

    ! Basal friction and shear stress
    ALLOCATE( ice%basal_friction_coefficient  ( mesh%vi1:mesh%vi2        ))  ! [Pa yr m^-1]       Effective basal friction coefficient (basal_shear_stress = u_base * basal_friction_coefficient)
    ALLOCATE( ice%basal_shear_stress          ( mesh%vi1:mesh%vi2        ))  ! [Pa]               Basal shear stress

    ice%basal_friction_coefficient  = 0._dp
    ice%basal_shear_stress          = 0._dp

  ! == Geothermal heat ==
  ! =====================

    ALLOCATE( ice%geothermal_heat_flux        ( mesh%vi1:mesh%vi2        ))  ! [J m^-2 yr^-1] Geothermal heat flux

    ice%geothermal_heat_flux        = 0._dp

  ! === Ice thickness time stepping ===
  ! ===================================

    ! Predicted model state at next time step
    ALLOCATE( ice%Hi_prev                     ( mesh%vi1:mesh%vi2        ))  ! [m]  The previous state
    ALLOCATE( ice%Hi_next                     ( mesh%vi1:mesh%vi2        ))  ! [m]  The next state

    ice%Hi_prev                     = 0._dp
    ice%Hi_next                     = 0._dp

  ! === Ice temperature time stepping ===
  ! =====================================

    ! Predicted model state at next time step
    ALLOCATE( ice%Ti_prev                     ( mesh%vi1:mesh%vi2,mesh%nz))  ! [m]  The previous state
    ALLOCATE( ice%Ti_next                     ( mesh%vi1:mesh%vi2,mesh%nz))  ! [m]  The next state

    ice%Ti_prev                     = 0._dp
    ice%Ti_next                     = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_ice_model

END MODULE ice_model_memory
