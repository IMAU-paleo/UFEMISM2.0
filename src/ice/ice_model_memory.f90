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

subroutine allocate_ice_model( mesh, ice)
    ! Allocate ice model variables

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'allocate_ice_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! === Masks ===
    ! =============

    allocate( ice%mask_land     ( mesh%nV ))
    allocate( ice%mask_ocean    ( mesh%nV ))
    allocate( ice%mask_lake     ( mesh%nV ))
    allocate( ice%mask_ice      ( mesh%nV ))
    allocate( ice%mask_sheet    ( mesh%nV ))
    allocate( ice%mask_shelf    ( mesh%nV ))
    allocate( ice%mask_coast    ( mesh%nV ))
    allocate( ice%mask_margin   ( mesh%nV ))
    allocate( ice%mask_gl_gr    ( mesh%nV ))
    allocate( ice%mask_gl_fl    ( mesh%nV ))
    allocate( ice%mask_cf_gr    ( mesh%nV ))
    allocate( ice%mask_cf_fl    ( mesh%nV ))
    allocate( ice%mask_ice_prev ( mesh%nV ))
    ice%mask_ice_prev = .false.
    allocate( ice%mask          ( mesh%nV_loc ))

    ! === Basic geometry ===
    ! ======================

    allocate( ice%Hi  ( mesh%vi1:mesh%vi2 ))   ! Ice thickness
    allocate( ice%Hb  ( mesh%vi1:mesh%vi2 ))   ! Bedrock height
    allocate( ice%Hs  ( mesh%vi1:mesh%vi2 ))   ! Ice surface height
    allocate( ice%Hib ( mesh%vi1:mesh%vi2 ))   ! Ice base height
    allocate( ice%TAF ( mesh%vi1:mesh%vi2 ))   ! Thickness above floatation

    ! === Geometry changes ===
    ! ========================

    allocate( ice%dHi  ( mesh%vi1:mesh%vi2 ))
    ice%dHi = 0d0
    allocate( ice%dHs  ( mesh%vi1:mesh%vi2 ))
    ice%dHs = 0d0
    allocate( ice%dHb  ( mesh%vi1:mesh%vi2 ))
    ice%dHb = 0d0
    allocate( ice%dHib ( mesh%vi1:mesh%vi2 ))
    ice%dHib = 0d0

    ! Predicted ice thickness at next step
    allocate( ice%Hi_tplusdt ( mesh%vi1:mesh%vi2 ))
    ice%Hi_tplusdt = 0d0

    ! === Geometry rates of change ===
    ! ================================

    allocate( ice%dHi_dt  ( mesh%vi1:mesh%vi2 ))
    ice%dHi_dt = 0d0
    allocate( ice%dHb_dt  ( mesh%vi1:mesh%vi2 ))
    ice%dHb_dt = 0d0
    allocate( ice%dHs_dt  ( mesh%vi1:mesh%vi2 ))
    ice%dHs_dt = 0d0
    allocate( ice%dHib_dt ( mesh%vi1:mesh%vi2 ))
    ice%dHib_dt = 0d0

    ! === Terrain-following coordinate zeta gradients ===
    ! ===================================================

    ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta

    ! On the ak-grid (vertices, vertically regular)
    ALLOCATE( ice%dzeta_dt_ak(    mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dzeta_dx_ak(    mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dzeta_dy_ak(    mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%dzeta_dz_ak(    mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%d2zeta_dx2_ak(  mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%d2zeta_dxdy_ak( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( ice%d2zeta_dy2_ak(  mesh%vi1:mesh%vi2,mesh%nz))
    ice%dzeta_dt_ak    = 0._dp
    ice%dzeta_dx_ak    = 0._dp
    ice%dzeta_dy_ak    = 0._dp
    ice%dzeta_dz_ak    = 0._dp
    ice%d2zeta_dx2_ak  = 0._dp
    ice%d2zeta_dxdy_ak = 0._dp
    ice%d2zeta_dy2_ak  = 0._dp

    ! On the bk-grid (triangles, vertically regular)
    ALLOCATE( ice%dzeta_dx_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%dzeta_dy_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%dzeta_dz_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%d2zeta_dx2_bk(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%d2zeta_dxdy_bk( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( ice%d2zeta_dy2_bk(  mesh%ti1:mesh%ti2,mesh%nz))
    ice%dzeta_dx_bk    = 0._dp
    ice%dzeta_dy_bk    = 0._dp
    ice%dzeta_dz_bk    = 0._dp
    ice%d2zeta_dx2_bk  = 0._dp
    ice%d2zeta_dxdy_bk = 0._dp
    ice%d2zeta_dy2_bk  = 0._dp

    ! On the bks-grid (triangles, vertically staggered)
    ALLOCATE( ice%dzeta_dx_bks(    mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%dzeta_dy_bks(    mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%dzeta_dz_bks(    mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%d2zeta_dx2_bks(  mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%d2zeta_dxdy_bks( mesh%ti1:mesh%ti2,mesh%nz-1))
    ALLOCATE( ice%d2zeta_dy2_bks(  mesh%ti1:mesh%ti2,mesh%nz-1))
    ice%dzeta_dx_bks    = 0._dp
    ice%dzeta_dy_bks    = 0._dp
    ice%dzeta_dz_bks    = 0._dp
    ice%d2zeta_dx2_bks  = 0._dp
    ice%d2zeta_dxdy_bks = 0._dp
    ice%d2zeta_dy2_bks  = 0._dp

    ! === Thermodynamics ===
    ! ======================

    ! Ice temperatures
    allocate( ice%Ti     ( mesh%vi1:mesh%vi2, C%nz ))   ! Ice temperatures
    ice%Ti = 0d0
    allocate( ice%Ti_pmp ( mesh%vi1:mesh%vi2, C%nz ))   ! Ice pressure melting point
    ice%Ti_pmp = 0d0

    ! Physical quantities
    allocate( ice%Cpi ( mesh%vi1:mesh%vi2, C%nz ))   ! Ice specific heat capacity
    ice%Cpi = 0d0
    allocate( ice%Ki  ( mesh%vi1:mesh%vi2, C%nz ))   ! Ice conductivity
    ice%Ki = 0d0

    ! Heating
    allocate( ice%internal_heating   ( mesh%vi1:mesh%vi2, C%nz ))
    ice%internal_heating = 0d0
    allocate( ice%frictional_heating ( mesh%vi1:mesh%vi2       ))
    ice%frictional_heating = 0d0

    ! === Ice flow ===
    ! ================

    allocate( ice%A_flow_3D  ( mesh%vi1:mesh%vi2, C%nz ))   ! Ice flow factor
    ice%A_flow_3D = 0d0
    allocate( ice%A_flow_vav ( mesh%vi1:mesh%vi2       ))   ! Vertically integrated
    ice%A_flow_vav = 0d0

    ! === Ice velocities ===
    ! ======================

    ! 3-D
    allocate( ice%u_3D   ( mesh%vi1:mesh%vi2, C%nz   ))
    ice%u_3D = 0d0
    allocate( ice%v_3D   ( mesh%vi1:mesh%vi2, C%nz   ))
    ice%v_3D = 0d0
    allocate( ice%u_3D_b ( mesh%ti1:mesh%ti2, C%nz ))
    ice%u_3D_b = 0d0
    allocate( ice%v_3D_b ( mesh%ti1:mesh%ti2, C%nz ))
    ice%v_3D_b = 0d0
    allocate( ice%w_3D   ( mesh%vi1:mesh%vi2, C%nz   ))
    ice%w_3D = 0d0

    ! Vertically integrated
    allocate( ice%u_vav      ( mesh%vi1:mesh%vi2   ))
    ice%u_vav = 0d0
    allocate( ice%v_vav      ( mesh%vi1:mesh%vi2   ))
    ice%v_vav = 0d0
    allocate( ice%u_vav_b    ( mesh%ti1:mesh%ti2 ))
    ice%u_vav_b = 0d0
    allocate( ice%v_vav_b    ( mesh%ti1:mesh%ti2 ))
    ice%v_vav_b = 0d0
    allocate( ice%uabs_vav   ( mesh%vi1:mesh%vi2   ))
    ice%uabs_vav = 0d0
    allocate( ice%uabs_vav_b ( mesh%ti1:mesh%ti2 ))
    ice%uabs_vav_b = 0d0

    ! Surface
    allocate( ice%u_surf      ( mesh%vi1:mesh%vi2   ))
    ice%u_surf = 0d0
    allocate( ice%v_surf      ( mesh%vi1:mesh%vi2   ))
    ice%v_surf = 0d0
    allocate( ice%u_surf_b    ( mesh%ti1:mesh%ti2 ))
    ice%u_surf_b = 0d0
    allocate( ice%v_surf_b    ( mesh%ti1:mesh%ti2 ))
    ice%v_surf_b = 0d0
    allocate( ice%w_surf      ( mesh%vi1:mesh%vi2   ))
    ice%w_surf = 0d0
    allocate( ice%uabs_surf   ( mesh%vi1:mesh%vi2   ))
    ice%uabs_surf = 0d0
    allocate( ice%uabs_surf_b ( mesh%ti1:mesh%ti2 ))
    ice%uabs_surf_b = 0d0

    ! Basal
    allocate( ice%u_base      ( mesh%vi1:mesh%vi2   ))
    ice%u_base = 0d0
    allocate( ice%v_base      ( mesh%vi1:mesh%vi2   ))
    ice%v_base = 0d0
    allocate( ice%u_base_b    ( mesh%ti1:mesh%ti2 ))
    ice%u_base_b = 0d0
    allocate( ice%v_base_b    ( mesh%ti1:mesh%ti2 ))
    ice%v_base_b = 0d0
    allocate( ice%w_base      ( mesh%vi1:mesh%vi2   ))
    ice%w_base = 0d0
    allocate( ice%uabs_base   ( mesh%vi1:mesh%vi2   ))
    ice%uabs_base = 0d0
    allocate( ice%uabs_base_b ( mesh%ti1:mesh%ti2 ))
    ice%uabs_base_b = 0d0

    ! === Strain rates ===
    ! ====================

    ! u-component
    allocate( ice%du_dx_3D ( mesh%vi1:mesh%vi2, C%nz ))
    allocate( ice%du_dy_3D ( mesh%vi1:mesh%vi2, C%nz ))
    allocate( ice%du_dz_3D ( mesh%vi1:mesh%vi2, C%nz ))
    ice%du_dx_3D = 0d0
    ice%du_dy_3D = 0d0
    ice%du_dz_3D = 0d0

    ! v-component
    allocate( ice%dv_dx_3D ( mesh%vi1:mesh%vi2, C%nz ))
    allocate( ice%dv_dy_3D ( mesh%vi1:mesh%vi2, C%nz ))
    allocate( ice%dv_dz_3D ( mesh%vi1:mesh%vi2, C%nz ))
    ice%dv_dx_3D = 0d0
    ice%dv_dy_3D = 0d0
    ice%dv_dz_3D = 0d0

    ! w-component
    allocate( ice%dw_dx_3D ( mesh%vi1:mesh%vi2, C%nz ))
    allocate( ice%dw_dy_3D ( mesh%vi1:mesh%vi2, C%nz ))
    allocate( ice%dw_dz_3D ( mesh%vi1:mesh%vi2, C%nz ))
    ice%dw_dx_3D = 0d0
    ice%dw_dy_3D = 0d0
    ice%dw_dz_3D = 0d0

    ! === Basal conditions ===
    ! ========================

    ! Basal hydrology
    allocate( ice%pore_water_pressure ( mesh%vi1:mesh%vi2 ))
    ice%pore_water_pressure = 0d0
    allocate( ice%overburden_pressure ( mesh%vi1:mesh%vi2 ))
    ice%overburden_pressure = 0d0
    allocate( ice%effective_pressure  ( mesh%vi1:mesh%vi2 ))
    ice%effective_pressure = 0d0

    ! Basal sliding
    allocate( ice%friction_coef_1    ( mesh%vi1:mesh%vi2 ))
    ice%friction_coef_1 = 0d0
    allocate( ice%friction_coef_2    ( mesh%vi1:mesh%vi2 ))
    ice%friction_coef_2 = 0d0
    allocate( ice%basal_shear_stress ( mesh%vi1:mesh%vi2 ))
    ice%basal_shear_stress = 0d0

    ! Geothermal heat flux
    allocate( ice%geothermal_heat_flux ( mesh%vi1:mesh%vi2 ))
    ice%geothermal_heat_flux = 0d0

    ! === Sea level ===
    ! =================

    ! Regional sea level
    allocate( ice%SL ( mesh%vi1:mesh%vi2 ))
    ice%SL = 0d0

    ! === Area fractions ===
    ! ======================

    ! Grounded
    allocate( ice%bedrock_cdf   ( mesh%vi1:mesh%vi2, 11 ))
    ice%bedrock_cdf = 0d0
    allocate( ice%fraction_gr   ( mesh%vi1:mesh%vi2     ))
    ice%fraction_gr = 0d0
    allocate( ice%fraction_gr_b ( mesh%ti1:mesh%ti2   ))
    ice%fraction_gr_b = 0d0

    ! Calving front
    allocate( ice%fraction_cf   ( mesh%vi1:mesh%vi2   ))
    ice%fraction_cf = 0d0

    ! === Basins ===
    ! ==============

    ! Basin ID's
    allocate( ice%basin_ID ( mesh%vi1:mesh%vi2 ))
    ice%basin_ID = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_ice_model

END MODULE ice_model_memory
