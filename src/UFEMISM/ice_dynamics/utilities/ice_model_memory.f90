module ice_model_memory
  !< Routines for administrating the memory for the ice model data.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model

  implicit none

  private

  public :: allocate_ice_model

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
    allocate( ice%Hi      ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Hb      ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Hs      ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%SL      ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Hib     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%TAF     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Hi_eff  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Hs_slope( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Geometry changes
    allocate( ice%dHi ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHb ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHs ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHib( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Rates of change
    allocate( ice%dHi_dt         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHb_dt         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHs_dt         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHib_dt        ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHi_dt_raw     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHi_dt_residual( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Horizontal derivatives
    allocate( ice%dHib_dx_b( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%dHib_dy_b( mesh%ti1:mesh%ti2), source = 0._dp)

    ! Target quantities
    allocate( ice%dHi_dt_target   ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Masks
    allocate( ice%mask_icefree_land      ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_icefree_ocean     ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_grounded_ice      ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_floating_ice      ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_icefree_land_prev ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_icefree_ocean_prev( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_grounded_ice_prev ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_floating_ice_prev ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_margin            ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_gl_gr             ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_gl_fl             ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_cf_gr             ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_cf_fl             ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_coastline         ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_ROI               ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_SGD               ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_noice             ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask                   ( mesh%vi1:mesh%vi2), source = 0)
    allocate( ice%basin_ID               ( mesh%vi1:mesh%vi2), source = 0)

    ! Area fractions
    allocate( ice%fraction_gr    ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%fraction_gr_b  ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%fraction_margin( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Sub-grid bedrock cumulative density functions (CDFs)
    allocate( ice%bedrock_cdf  ( mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins), source = 0._dp)
    allocate( ice%bedrock_cdf_b( mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins), source = 0._dp)

    ! === Terrain-following coordinate zeta gradients ===
    ! ===================================================

    ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta

    ! On the ak-grid (vertices, vertically regular)
    allocate( ice%dzeta_dt_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dx_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dy_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dz_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dx2_ak ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dxdy_ak( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dy2_ak ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! On the bk-grid (triangles, vertically regular)
    allocate( ice%dzeta_dx_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dy_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dz_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dx2_bk ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dxdy_bk( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dy2_bk ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)

    ! On the bks-grid (triangles, vertically staggered)
    allocate( ice%dzeta_dx_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%dzeta_dy_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%dzeta_dz_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%d2zeta_dx2_bks ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%d2zeta_dxdy_bks( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%d2zeta_dy2_bks ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)

    ! === Thermodynamics and rheology ===
    ! ===================================

    ! Ice temperatures
    allocate( ice%Ti    ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ti_pmp( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ti_hom( mesh%vi1:mesh%vi2        ), source = 0._dp)

    ! Physical quantities
    allocate( ice%Cpi( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ki ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! Heating
    allocate( ice%internal_heating  ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%frictional_heating( mesh%vi1:mesh%vi2        ), source = 0._dp)

    ! Glen's flow law factor
    allocate( ice%A_flow( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! === Ice velocities ===
    ! ======================

    ! 3-D
    allocate( ice%u_3D  ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%v_3D  ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%u_3D_b( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%v_3D_b( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%w_3D  ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! Vertically integrated
    allocate( ice%u_vav     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%v_vav     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%u_vav_b   ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%v_vav_b   ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%uabs_vav  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%uabs_vav_b( mesh%ti1:mesh%ti2), source = 0._dp)

    ! Surface
    allocate( ice%u_surf     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%v_surf     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%u_surf_b   ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%v_surf_b   ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%w_surf     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%uabs_surf  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%uabs_surf_b( mesh%ti1:mesh%ti2), source = 0._dp)

    ! Basal
    allocate( ice%u_base     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%v_base     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%u_base_b   ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%v_base_b   ( mesh%ti1:mesh%ti2), source = 0._dp)
    allocate( ice%w_base     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%uabs_base  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%uabs_base_b( mesh%ti1:mesh%ti2), source = 0._dp)

    ! == Strain rates ==
    ! ==================

    allocate( ice%du_dx_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%du_dy_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%du_dz_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dv_dx_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dv_dy_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dv_dz_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dw_dx_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dw_dy_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dw_dz_3D( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! == Ice flow regime ==
    ! =====================

    allocate( ice%divQ   ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%R_shear( mesh%vi1:mesh%vi2), source = 0._dp)

    ! == Basal hydrology ==
    ! =====================

    ! Basal hydrology
    allocate( ice%pore_water_pressure  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%overburden_pressure  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%effective_pressure   ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%pore_water_likelihood( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%pore_water_fraction  ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! == Basal sliding ==
    ! ===================

    ! Basal friction and shear stress
    allocate( ice%till_yield_stress         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%basal_friction_coefficient( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%basal_shear_stress        ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! == Geothermal heat ==
    ! =====================

    allocate( ice%geothermal_heat_flux( mesh%vi1:mesh%vi2), source = 0._dp)

    ! === Ice thickness time stepping ===
    ! ===================================

    ! Predicted model state at next time step
    allocate( ice%Hi_prev( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Hi_next( mesh%vi1:mesh%vi2), source = 0._dp)

    ! === Ice temperature time stepping ===
    ! =====================================

    ! Predicted model state at next time step
    allocate( ice%Ti_prev( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ti_next( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_ice_model

end module ice_model_memory
