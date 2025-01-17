module ice_model_memory
  !< Routines for administrating the memory for the ice model data.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model

  implicit none

contains

  subroutine allocate_ice_model( mesh, ice)

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(  out) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_ice_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! === Ice-sheet geometry ===
    ! ==========================

    ! Basic geometry
    allocate( ice%Hi      ( mesh%vi1:mesh%vi2))  ! [m] Ice thickness
    allocate( ice%Hb      ( mesh%vi1:mesh%vi2))  ! [m] Bedrock elevation (w.r.t. PD sea level)
    allocate( ice%Hs      ( mesh%vi1:mesh%vi2))  ! [m] Surface elevation (w.r.t. PD sea level)
    allocate( ice%SL      ( mesh%vi1:mesh%vi2))  ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    allocate( ice%Hib     ( mesh%vi1:mesh%vi2))  ! [m] Ice base elevation (w.r.t. PD sea level)
    allocate( ice%TAF     ( mesh%vi1:mesh%vi2))  ! [m] Thickness above flotation
    allocate( ice%Hi_eff  ( mesh%vi1:mesh%vi2))  ! [m] Thickness above flotation
    allocate( ice%Hs_slope( mesh%vi1:mesh%vi2))  ! [-] Absolute surface gradient

    ice%Hi                          = 0._dp
    ice%Hb                          = 0._dp
    ice%Hs                          = 0._dp
    ice%SL                          = 0._dp
    ice%Hib                         = 0._dp
    ice%TAF                         = 0._dp
    ice%Hi_eff                      = 0._dp
    ice%Hs_slope                    = 0._dp

    ! Geometry changes
    allocate( ice%dHi ( mesh%vi1:mesh%vi2))  ! [m] Ice thickness difference (w.r.t. reference)
    allocate( ice%dHb ( mesh%vi1:mesh%vi2))  ! [m] Bedrock elevation difference (w.r.t. reference)
    allocate( ice%dHs ( mesh%vi1:mesh%vi2))  ! [m] Surface elevation difference (w.r.t. reference)
    allocate( ice%dHib( mesh%vi1:mesh%vi2))  ! [m] Base elevation difference (w.r.t. reference)

    ice%dHi                         = 0._dp
    ice%dHb                         = 0._dp
    ice%dHs                         = 0._dp
    ice%dHib                        = 0._dp

    ! Rates of change
    allocate( ice%dHi_dt         ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Ice thickness rate of change
    allocate( ice%dHb_dt         ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Bedrock elevation rate of change
    allocate( ice%dHs_dt         ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Ice surface elevation rate of change
    allocate( ice%dHib_dt        ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Ice base elevation rate of change
    allocate( ice%dHi_dt_raw     ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Ice thickness rate of change before any imposed modifications
    allocate( ice%dHi_dt_residual( mesh%vi1:mesh%vi2))  ! [m yr^-1] Residual ice thickness rate of change after imposed modifications

    ice%dHi_dt                      = 0._dp
    ice%dHb_dt                      = 0._dp
    ice%dHs_dt                      = 0._dp
    ice%dHib_dt                     = 0._dp
    ice%dHi_dt_raw                  = 0._dp
    ice%dHi_dt_residual             = 0._dp

    ! Horizontal derivatives
    allocate( ice%dHib_dx_b( mesh%ti1:mesh%ti2))  ! [] Horizontal derivative of ice draft on b-grid
    allocate( ice%dHib_dy_b( mesh%ti1:mesh%ti2))  ! [] Horizontal derivative of ice draft on b-grid

    ice%dHib_dx_b                   = 0._dp
    ice%dHib_dy_b                   = 0._dp

    ! Target quantities
    allocate( ice%dHi_dt_target   ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Target ice thickness rate of change for inversions
    allocate( ice%uabs_surf_target( mesh%vi1:mesh%vi2))  ! [m yr^-1] Target ice surface speed for inversions

    ice%dHi_dt_target               = 0._dp
    ice%uabs_surf_target            = 0._dp

    ! Masks
    allocate( ice%mask_icefree_land      ( mesh%vi1:mesh%vi2))  ! T: ice-free land , F: otherwise
    allocate( ice%mask_icefree_ocean     ( mesh%vi1:mesh%vi2))  ! T: ice-free ocean, F: otherwise
    allocate( ice%mask_grounded_ice      ( mesh%vi1:mesh%vi2))  ! T: grounded ice  , F: otherwise
    allocate( ice%mask_floating_ice      ( mesh%vi1:mesh%vi2))  ! T: floating ice  , F: otherwise
    allocate( ice%mask_icefree_land_prev ( mesh%vi1:mesh%vi2))  ! T: ice-free land , F: otherwise (during previous time step)
    allocate( ice%mask_icefree_ocean_prev( mesh%vi1:mesh%vi2))  ! T: ice-free ocean, F: otherwise (during previous time step)
    allocate( ice%mask_grounded_ice_prev ( mesh%vi1:mesh%vi2))  ! T: grounded ice  , F: otherwise (during previous time step)
    allocate( ice%mask_floating_ice_prev ( mesh%vi1:mesh%vi2))  ! T: floating ice  , F: otherwise (during previous time step)
    allocate( ice%mask_margin            ( mesh%vi1:mesh%vi2))  ! T: ice next to ice-free, F: otherwise
    allocate( ice%mask_gl_gr             ( mesh%vi1:mesh%vi2))  ! T: grounded ice next to floating ice, F: otherwise
    allocate( ice%mask_gl_fl             ( mesh%vi1:mesh%vi2))  ! T: floating ice next to grounded ice, F: otherwise
    allocate( ice%mask_cf_gr             ( mesh%vi1:mesh%vi2))  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    allocate( ice%mask_cf_fl             ( mesh%vi1:mesh%vi2))  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    allocate( ice%mask_coastline         ( mesh%vi1:mesh%vi2))  ! T: ice-free land next to ice-free ocean, F: otherwise
    allocate( ice%mask_noice             ( mesh%vi1:mesh%vi2))  ! T: no ice is allowed here, F: ice is allowed here
    allocate( ice%mask                   ( mesh%vi1:mesh%vi2))  ! Diagnostic, only meant for quick visual inspection in output
    allocate( ice%basin_ID               ( mesh%vi1:mesh%vi2))  ! The drainage basin to which each vertex belongs

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
    allocate( ice%fraction_gr    ( mesh%vi1:mesh%vi2))  ! [0-1] Grounded area fractions of vertices
    allocate( ice%fraction_gr_b  ( mesh%ti1:mesh%ti2))  ! [0-1] Grounded area fractions of triangles
    allocate( ice%fraction_margin( mesh%vi1:mesh%vi2))  ! [0-1] Ice-covered area fractions of calving fronts

    ice%fraction_gr                 = 0._dp
    ice%fraction_gr_b               = 0._dp
    ice%fraction_margin             = 0._dp

    ! Sub-grid bedrock cumulative density functions (CDFs)
    allocate( ice%bedrock_cdf  ( mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins))  ! [-] Sub-grid bedrock cumulative density functions on the a-grid (vertices)
    allocate( ice%bedrock_cdf_b( mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins))  ! [-] Sub-grid bedrock cumulative density functions on the b-grid (triangles)

    ice%bedrock_cdf                 = 0._dp
    ice%bedrock_cdf_b               = 0._dp

    ! === Terrain-following coordinate zeta gradients ===
    ! ===================================================

    ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta

    ! On the ak-grid (vertices, vertically regular)
    allocate( ice%dzeta_dt_ak   ( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dzeta_dx_ak   ( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dzeta_dy_ak   ( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dzeta_dz_ak   ( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%d2zeta_dx2_ak ( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%d2zeta_dxdy_ak( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%d2zeta_dy2_ak ( mesh%vi1:mesh%vi2,mesh%nz))

    ice%dzeta_dt_ak                 = 0._dp
    ice%dzeta_dx_ak                 = 0._dp
    ice%dzeta_dy_ak                 = 0._dp
    ice%dzeta_dz_ak                 = 0._dp
    ice%d2zeta_dx2_ak               = 0._dp
    ice%d2zeta_dxdy_ak              = 0._dp
    ice%d2zeta_dy2_ak               = 0._dp

    ! On the bk-grid (triangles, vertically regular)
    allocate( ice%dzeta_dx_bk   ( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( ice%dzeta_dy_bk   ( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( ice%dzeta_dz_bk   ( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( ice%d2zeta_dx2_bk ( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( ice%d2zeta_dxdy_bk( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( ice%d2zeta_dy2_bk ( mesh%ti1:mesh%ti2,mesh%nz))

    ice%dzeta_dx_bk                 = 0._dp
    ice%dzeta_dy_bk                 = 0._dp
    ice%dzeta_dz_bk                 = 0._dp
    ice%d2zeta_dx2_bk               = 0._dp
    ice%d2zeta_dxdy_bk              = 0._dp
    ice%d2zeta_dy2_bk               = 0._dp

    ! On the bks-grid (triangles, vertically staggered)
    allocate( ice%dzeta_dx_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1))
    allocate( ice%dzeta_dy_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1))
    allocate( ice%dzeta_dz_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1))
    allocate( ice%d2zeta_dx2_bks ( mesh%ti1:mesh%ti2,mesh%nz-1))
    allocate( ice%d2zeta_dxdy_bks( mesh%ti1:mesh%ti2,mesh%nz-1))
    allocate( ice%d2zeta_dy2_bks ( mesh%ti1:mesh%ti2,mesh%nz-1))

    ice%dzeta_dx_bks                = 0._dp
    ice%dzeta_dy_bks                = 0._dp
    ice%dzeta_dz_bks                = 0._dp
    ice%d2zeta_dx2_bks              = 0._dp
    ice%d2zeta_dxdy_bks             = 0._dp
    ice%d2zeta_dy2_bks              = 0._dp

    ! === Thermodynamics and rheology ===
    ! ===================================

    ! Ice temperatures
    allocate( ice%Ti    ( mesh%vi1:mesh%vi2,mesh%nz))  ! [K] Englacial temperature
    allocate( ice%Ti_pmp( mesh%vi1:mesh%vi2,mesh%nz))  ! [K] Pressure melting point temperature
    allocate( ice%Ti_hom( mesh%vi1:mesh%vi2        ))  ! [K] Basal temperature w.r.t. pressure melting point

    ice%Ti                          = 0._dp
    ice%Ti_pmp                      = 0._dp
    ice%Ti_hom                      = 0._dp

    ! Physical quantities
    allocate( ice%Cpi( mesh%vi1:mesh%vi2,mesh%nz))  ! [J kg^-1 K^-1] Specific heat capacity
    allocate( ice%Ki ( mesh%vi1:mesh%vi2,mesh%nz))  ! [J m^-1 K^-1 yr^-1] Thermal conductivity

    ice%Cpi                         = 0._dp
    ice%Ki                          = 0._dp

    ! Heating
    allocate( ice%internal_heating  ( mesh%vi1:mesh%vi2,mesh%nz))  ! [?] Internal heating
    allocate( ice%frictional_heating( mesh%vi1:mesh%vi2        ))  ! [?] Frictional heating

    ice%internal_heating            = 0._dp
    ice%frictional_heating          = 0._dp

    ! Glen's flow law factor
    allocate( ice%A_flow( mesh%vi1:mesh%vi2,mesh%nz))  ! [Pa^-3 y^-1] Glen's flow law factor

    ice%A_flow                      = 0._dp

    ! === Ice velocities ===
    ! ======================

    ! 3-D
    allocate( ice%u_3D  ( mesh%vi1:mesh%vi2,mesh%nz))  ! [m yr^-1] 3-D ice velocity
    allocate( ice%v_3D  ( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%u_3D_b( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( ice%v_3D_b( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( ice%w_3D  ( mesh%vi1:mesh%vi2,mesh%nz))

    ice%u_3D                        = 0._dp
    ice%v_3D                        = 0._dp
    ice%u_3D_b                      = 0._dp
    ice%v_3D_b                      = 0._dp
    ice%w_3D                        = 0._dp

    ! Vertically integrated
    allocate( ice%u_vav     ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Vertically averaged ice velocity
    allocate( ice%v_vav     ( mesh%vi1:mesh%vi2))
    allocate( ice%u_vav_b   ( mesh%ti1:mesh%ti2))
    allocate( ice%v_vav_b   ( mesh%ti1:mesh%ti2))
    allocate( ice%uabs_vav  ( mesh%vi1:mesh%vi2))
    allocate( ice%uabs_vav_b( mesh%ti1:mesh%ti2))

    ice%u_vav                       = 0._dp
    ice%v_vav                       = 0._dp
    ice%u_vav_b                     = 0._dp
    ice%v_vav_b                     = 0._dp
    ice%uabs_vav                    = 0._dp
    ice%uabs_vav_b                  = 0._dp

    ! Surface
    allocate( ice%u_surf     ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Ice velocity at the surface
    allocate( ice%v_surf     ( mesh%vi1:mesh%vi2))
    allocate( ice%u_surf_b   ( mesh%ti1:mesh%ti2))
    allocate( ice%v_surf_b   ( mesh%ti1:mesh%ti2))
    allocate( ice%w_surf     ( mesh%vi1:mesh%vi2))
    allocate( ice%uabs_surf  ( mesh%vi1:mesh%vi2))
    allocate( ice%uabs_surf_b( mesh%ti1:mesh%ti2))

    ice%u_surf                      = 0._dp
    ice%v_surf                      = 0._dp
    ice%u_surf_b                    = 0._dp
    ice%v_surf_b                    = 0._dp
    ice%w_surf                      = 0._dp
    ice%uabs_surf                   = 0._dp
    ice%uabs_surf_b                 = 0._dp

    ! Basal
    allocate( ice%u_base     ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Ice velocity at the base
    allocate( ice%v_base     ( mesh%vi1:mesh%vi2))
    allocate( ice%u_base_b   ( mesh%ti1:mesh%ti2))
    allocate( ice%v_base_b   ( mesh%ti1:mesh%ti2))
    allocate( ice%w_base     ( mesh%vi1:mesh%vi2))
    allocate( ice%uabs_base  ( mesh%vi1:mesh%vi2))
    allocate( ice%uabs_base_b( mesh%ti1:mesh%ti2))

    ice%u_base                      = 0._dp
    ice%v_base                      = 0._dp
    ice%u_base_b                    = 0._dp
    ice%v_base_b                    = 0._dp
    ice%w_base                      = 0._dp
    ice%uabs_base                   = 0._dp
    ice%uabs_base_b                 = 0._dp

    ! == Strain rates ==
    ! ==================

    allocate( ice%du_dx_3D( mesh%vi1:mesh%vi2,mesh%nz))  ! [yr^-1]
    allocate( ice%du_dy_3D( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%du_dz_3D( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dv_dx_3D( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dv_dy_3D( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dv_dz_3D( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dw_dx_3D( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dw_dy_3D( mesh%vi1:mesh%vi2,mesh%nz))
    allocate( ice%dw_dz_3D( mesh%vi1:mesh%vi2,mesh%nz))

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

    allocate( ice%divQ   ( mesh%vi1:mesh%vi2))  ! [m yr^-1] Horizontal ice flux divergence
    allocate( ice%R_shear( mesh%vi1:mesh%vi2))  ! [0-1]     uabs_base / uabs_surf (0 = pure vertical shear, viscous flow; 1 = pure sliding, plug flow)

    ice%divQ                        = 0._dp
    ice%R_shear                     = 0._dp

  ! == Basal hydrology ==
  ! =====================

    ! Basal hydrology
    allocate( ice%pore_water_pressure  ( mesh%vi1:mesh%vi2))  ! [Pa]  Basal pore water pressure
    allocate( ice%overburden_pressure  ( mesh%vi1:mesh%vi2))  ! [Pa]  Basal overburden pressure
    allocate( ice%effective_pressure   ( mesh%vi1:mesh%vi2))  ! [Pa]  Basal effective pressure
    allocate( ice%pore_water_likelihood( mesh%vi1:mesh%vi2))  ! [0-1] Basal pore water likelihood
    allocate( ice%pore_water_fraction  ( mesh%vi1:mesh%vi2))  ! [0-1] Fraction of overburden pressure reduced by pore water pressure

    ice%pore_water_pressure         = 0._dp
    ice%overburden_pressure         = 0._dp
    ice%effective_pressure          = 0._dp
    ice%pore_water_likelihood       = 0._dp
    ice%pore_water_fraction         = 0._dp

  ! == Basal sliding ==
  ! ===================

    ! Sliding law coefficients
    allocate( ice%till_friction_angle( mesh%vi1:mesh%vi2))  ! [degrees]          Till friction angle (degrees)
    allocate( ice%bed_roughness      ( mesh%vi1:mesh%vi2))  ! [0-1]              Bed roughness fraction
    allocate( ice%till_yield_stress  ( mesh%vi1:mesh%vi2))  ! [Pa]               Till yield stress (used when choice_sliding_law = "Coloumb", "Budd", or "Zoet-Iverson")
    allocate( ice%slid_alpha_sq      ( mesh%vi1:mesh%vi2))  ! [-]                Coulomb-law friction coefficient (used when choice_sliding_law = "Tsai2015", or "Schoof2005")
    allocate( ice%slid_beta_sq       ( mesh%vi1:mesh%vi2))  ! [Pa m^−1/m yr^1/m] Power-law friction coefficient (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")

    ice%till_friction_angle         = 0._dp
    ice%bed_roughness               = 0._dp
    ice%till_yield_stress           = 0._dp
    ice%slid_alpha_sq               = 0._dp
    ice%slid_beta_sq                = 0._dp

    ! Basal friction and shear stress
    allocate( ice%basal_friction_coefficient( mesh%vi1:mesh%vi2))  ! [Pa yr m^-1]       Effective basal friction coefficient (basal_shear_stress = u_base * basal_friction_coefficient)
    allocate( ice%basal_shear_stress        ( mesh%vi1:mesh%vi2))  ! [Pa]               Basal shear stress

    ice%basal_friction_coefficient  = 0._dp
    ice%basal_shear_stress          = 0._dp

    ! == Geothermal heat ==
    ! =====================

    allocate( ice%geothermal_heat_flux( mesh%vi1:mesh%vi2))  ! [J m^-2 yr^-1] Geothermal heat flux

    ice%geothermal_heat_flux        = 0._dp

    ! === Ice thickness time stepping ===
    ! ===================================

    ! Predicted model state at next time step
    allocate( ice%Hi_prev( mesh%vi1:mesh%vi2))  ! [m]  The previous state
    allocate( ice%Hi_next( mesh%vi1:mesh%vi2))  ! [m]  The next state

    ice%Hi_prev                     = 0._dp
    ice%Hi_next                     = 0._dp

  ! === Ice temperature time stepping ===
  ! =====================================

    ! Predicted model state at next time step
    allocate( ice%Ti_prev( mesh%vi1:mesh%vi2,mesh%nz))  ! [m]  The previous state
    allocate( ice%Ti_next( mesh%vi1:mesh%vi2,mesh%nz))  ! [m]  The next state

    ice%Ti_prev                     = 0._dp
    ice%Ti_next                     = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine allocate_ice_model

end module ice_model_memory
