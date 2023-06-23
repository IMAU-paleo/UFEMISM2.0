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

    ice%Hi                          = 0._dp
    ice%Hb                          = 0._dp
    ice%Hs                          = 0._dp
    ice%SL                          = 0._dp
    ice%Hib                         = 0._dp
    ice%TAF                         = 0._dp

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

    ice%dHi_dt                      = 0._dp
    ice%dHb_dt                      = 0._dp
    ice%dHs_dt                      = 0._dp
    ice%dHib_dt                     = 0._dp

    ! Masks
    ALLOCATE( ice%mask_icefree_land           ( mesh%vi1:mesh%vi2        ))  ! T: ice-free land , F: otherwise
    ALLOCATE( ice%mask_icefree_ocean          ( mesh%vi1:mesh%vi2        ))  ! T: ice-free ocean, F: otherwise
    ALLOCATE( ice%mask_grounded_ice           ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice  , F: otherwise
    ALLOCATE( ice%mask_floating_ice           ( mesh%vi1:mesh%vi2        ))  ! T: floating ice  , F: otherwise
    ALLOCATE( ice%mask_icefree_land_prev      ( mesh%vi1:mesh%vi2        ))  ! T: ice-free land , F: otherwise (during previous time step)
    ALLOCATE( ice%mask_icefree_ocean_prev     ( mesh%vi1:mesh%vi2        ))  ! T: ice-free ocean, F: otherwise (during previous time step)
    ALLOCATE( ice%mask_grounded_ice_prev      ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice  , F: otherwise (during previous time step)
    ALLOCATE( ice%mask_floating_ice_prev      ( mesh%vi1:mesh%vi2        ))  ! T: floating ice  , F: otherwise (during previous time step)
    ALLOCATE( ice%mask_gl_gr                  ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice next to floating ice, F: otherwise
    ALLOCATE( ice%mask_gl_fl                  ( mesh%vi1:mesh%vi2        ))  ! T: floating ice next to grounded ice, F: otherwise
    ALLOCATE( ice%mask_cf_gr                  ( mesh%vi1:mesh%vi2        ))  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    ALLOCATE( ice%mask_cf_fl                  ( mesh%vi1:mesh%vi2        ))  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
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
    ice%mask_gl_gr                  = .FALSE.
    ice%mask_gl_fl                  = .FALSE.
    ice%mask_cf_gr                  = .FALSE.
    ice%mask_cf_fl                  = .FALSE.
    ice%mask                        = 0
    ice%basin_ID                    = 0

    ! Area fractions
    ALLOCATE( ice%bedrock_cdf                 ( mesh%vi1:mesh%vi2, 11    ))  ! Sub-grid bedrock cumulative density functions
    ALLOCATE( ice%fraction_gr                 ( mesh%vi1:mesh%vi2        ))  ! [0-1] Grounded area fractions of vertices
    ALLOCATE( ice%fraction_gr_b               ( mesh%ti1:mesh%ti2        ))  ! [0-1] Grounded area fractions of triangles
    ALLOCATE( ice%fraction_cf                 ( mesh%vi1:mesh%vi2        ))  ! [0-1] Ice-covered area fractions of calving fronts

    ice%bedrock_cdf                 = 0._dp
    ice%fraction_gr                 = 0._dp
    ice%fraction_gr_b               = 0._dp
    ice%fraction_cf                 = 0._dp

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

    ice%Ti                          = 0._dp
    ice%Ti_pmp                      = 0._dp

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
    ALLOCATE( ice%A_flow_3D                   ( mesh%vi1:mesh%vi2,mesh%nz))  ! [Pa^-3 y^-1] Glen's flow law factor

    ice%A_flow_3D                   = 0._dp

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

    ! Strain rates
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

  ! == Basal conditions ==
  ! ======================

    ! Basal hydrology
    ALLOCATE( ice%pore_water_pressure         ( mesh%vi1:mesh%vi2        ))  ! Basal pore water pressure
    ALLOCATE( ice%overburden_pressure         ( mesh%vi1:mesh%vi2        ))  ! Basal overburden pressure
    ALLOCATE( ice%effective_pressure          ( mesh%vi1:mesh%vi2        ))  ! Basal effective pressure

    ice%pore_water_pressure         = 0._dp
    ice%overburden_pressure         = 0._dp
    ice%effective_pressure          = 0._dp

    ! Basal roughness / friction
    ALLOCATE( ice%phi_fric                    ( mesh%vi1:mesh%vi2        ))  ! Till friction angle (degrees)
    ALLOCATE( ice%tau_c                       ( mesh%vi1:mesh%vi2        ))  ! Till yield stress tauc   (used when choice_sliding_law = "Coloumb" or "Coulomb_regularised")
    ALLOCATE( ice%alpha_sq                    ( mesh%vi1:mesh%vi2        ))  ! Coulomb-law friction coefficient [unitless]         (used when choice_sliding_law =             "Tsai2015", or "Schoof2005")
    ALLOCATE( ice%beta_sq                     ( mesh%vi1:mesh%vi2        ))  ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3] (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")
    ALLOCATE( ice%beta_b                      ( mesh%vi1:mesh%vi2        ))  ! Basal friction (tau_b = u * beta_b)

    ice%phi_fric                    = 0._dp
    ice%tau_c                       = 0._dp
    ice%alpha_sq                    = 0._dp
    ice%beta_sq                     = 0._dp
    ice%beta_b                      = 0._dp

    ! Basal sliding
    ALLOCATE( ice%friction_coef_1             ( mesh%vi1:mesh%vi2        ))  ! Generic basal friction coefficient 1
    ALLOCATE( ice%friction_coef_2             ( mesh%vi1:mesh%vi2        ))  ! Generic basal friction coefficient 2
    ALLOCATE( ice%basal_shear_stress          ( mesh%vi1:mesh%vi2        ))  ! Basal shear stress

    ice%friction_coef_1             = 0._dp
    ice%friction_coef_2             = 0._dp
    ice%basal_shear_stress          = 0._dp

    ! Geothermal heat
    ALLOCATE( ice%geothermal_heat_flux        ( mesh%vi1:mesh%vi2        ))  ! Geothermal heat flux

    ice%geothermal_heat_flux        = 0._dp

  ! === Time stepping ===
  ! =====================

    ! Predicted model state at next time step
    ALLOCATE( ice%Hi_prev                     ( mesh%vi1:mesh%vi2        ))  ! [m]  The previous state
    ALLOCATE( ice%Hi_next                     ( mesh%vi1:mesh%vi2        ))  ! [m]  The next state

    ice%Hi_prev                     = 0._dp
    ice%Hi_next                     = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_ice_model

END MODULE ice_model_memory
