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

    allocate( ice%Hi  ( mesh%nV_loc ))   ! Ice thickness
    allocate( ice%Hb  ( mesh%nV_loc ))   ! Bedrock height
    allocate( ice%Hs  ( mesh%nV_loc ))   ! Ice surface height
    allocate( ice%Hib ( mesh%nV_loc ))   ! Ice base height
    allocate( ice%TAF ( mesh%nV_loc ))   ! Thickness above floatation

    ! === Geometry changes ===
    ! ========================

    allocate( ice%dHi  ( mesh%nV_loc ))
    ice%dHi = 0d0
    allocate( ice%dHs  ( mesh%nV_loc ))
    ice%dHs = 0d0
    allocate( ice%dHb  ( mesh%nV_loc ))
    ice%dHb = 0d0
    allocate( ice%dHib ( mesh%nV_loc ))
    ice%dHib = 0d0

    ! Predicted ice thickness at next step
    allocate( ice%Hi_tplusdt ( mesh%nV_loc ))
    ice%Hi_tplusdt = 0d0

    ! === Geometry rates of change ===
    ! ================================

    allocate( ice%dHi_dt  ( mesh%nV_loc ))
    ice%dHi_dt = 0d0
    allocate( ice%dHb_dt  ( mesh%nV_loc ))
    ice%dHb_dt = 0d0
    allocate( ice%dHs_dt  ( mesh%nV_loc ))
    ice%dHs_dt = 0d0
    allocate( ice%dHib_dt ( mesh%nV_loc ))
    ice%dHib_dt = 0d0

    ! === Thermodynamics ===
    ! ======================

    ! Ice temperatures
    allocate( ice%Ti     ( mesh%nV_loc, C%nz ))   ! Ice temperatures
    ice%Ti = 0d0
    allocate( ice%Ti_pmp ( mesh%nV_loc, C%nz ))   ! Ice pressure melting point
    ice%Ti_pmp = 0d0

    ! Physical quantities
    allocate( ice%Cpi ( mesh%nV_loc, C%nz ))   ! Ice specific heat capacity
    ice%Cpi = 0d0
    allocate( ice%Ki  ( mesh%nV_loc, C%nz ))   ! Ice conductivity
    ice%Ki = 0d0

    ! Heating
    allocate( ice%internal_heating   ( mesh%nV_loc, C%nz ))
    ice%internal_heating = 0d0
    allocate( ice%frictional_heating ( mesh%nV_loc       ))
    ice%frictional_heating = 0d0

    ! === Ice flow ===
    ! ================

    allocate( ice%A_flow_3D  ( mesh%nV_loc, C%nz ))   ! Ice flow factor
    ice%A_flow_3D = 0d0
    allocate( ice%A_flow_vav ( mesh%nV_loc       ))   ! Vertically integrated
    ice%A_flow_vav = 0d0

    ! === Ice velocities ===
    ! ======================

    ! 3-D
    allocate( ice%u_3D   ( mesh%nV_loc, C%nz   ))
    ice%u_3D = 0d0
    allocate( ice%v_3D   ( mesh%nV_loc, C%nz   ))
    ice%v_3D = 0d0
    allocate( ice%u_3D_b ( mesh%nTri_loc, C%nz ))
    ice%u_3D_b = 0d0
    allocate( ice%v_3D_b ( mesh%nTri_loc, C%nz ))
    ice%v_3D_b = 0d0
    allocate( ice%w_3D   ( mesh%nV_loc, C%nz   ))
    ice%w_3D = 0d0

    ! Vertically integrated
    allocate( ice%u_vav      ( mesh%nV_loc   ))
    ice%u_vav = 0d0
    allocate( ice%v_vav      ( mesh%nV_loc   ))
    ice%v_vav = 0d0
    allocate( ice%u_vav_b    ( mesh%nTri_loc ))
    ice%u_vav_b = 0d0
    allocate( ice%v_vav_b    ( mesh%nTri_loc ))
    ice%v_vav_b = 0d0
    allocate( ice%uabs_vav   ( mesh%nV_loc   ))
    ice%uabs_vav = 0d0
    allocate( ice%uabs_vav_b ( mesh%nTri_loc ))
    ice%uabs_vav_b = 0d0

    ! Surface
    allocate( ice%u_surf      ( mesh%nV_loc   ))
    ice%u_surf = 0d0
    allocate( ice%v_surf      ( mesh%nV_loc   ))
    ice%v_surf = 0d0
    allocate( ice%u_surf_b    ( mesh%nTri_loc ))
    ice%u_surf_b = 0d0
    allocate( ice%v_surf_b    ( mesh%nTri_loc ))
    ice%v_surf_b = 0d0
    allocate( ice%w_surf      ( mesh%nV_loc   ))
    ice%w_surf = 0d0
    allocate( ice%uabs_surf   ( mesh%nV_loc   ))
    ice%uabs_surf = 0d0
    allocate( ice%uabs_surf_b ( mesh%nTri_loc ))
    ice%uabs_surf_b = 0d0

    ! Basal
    allocate( ice%u_base      ( mesh%nV_loc   ))
    ice%u_base = 0d0
    allocate( ice%v_base      ( mesh%nV_loc   ))
    ice%v_base = 0d0
    allocate( ice%u_base_b    ( mesh%nTri_loc ))
    ice%u_base_b = 0d0
    allocate( ice%v_base_b    ( mesh%nTri_loc ))
    ice%v_base_b = 0d0
    allocate( ice%w_base      ( mesh%nV_loc   ))
    ice%w_base = 0d0
    allocate( ice%uabs_base   ( mesh%nV_loc   ))
    ice%uabs_base = 0d0
    allocate( ice%uabs_base_b ( mesh%nTri_loc ))
    ice%uabs_base_b = 0d0

    ! === Strain rates ===
    ! ====================

    ! u-component
    allocate( ice%du_dx_3D ( mesh%nV_loc, C%nz ))
    allocate( ice%du_dy_3D ( mesh%nV_loc, C%nz ))
    allocate( ice%du_dz_3D ( mesh%nV_loc, C%nz ))
    ice%du_dx_3D = 0d0
    ice%du_dy_3D = 0d0
    ice%du_dz_3D = 0d0

    ! v-component
    allocate( ice%dv_dx_3D ( mesh%nV_loc, C%nz ))
    allocate( ice%dv_dy_3D ( mesh%nV_loc, C%nz ))
    allocate( ice%dv_dz_3D ( mesh%nV_loc, C%nz ))
    ice%dv_dx_3D = 0d0
    ice%dv_dy_3D = 0d0
    ice%dv_dz_3D = 0d0

    ! w-component
    allocate( ice%dw_dx_3D ( mesh%nV_loc, C%nz ))
    allocate( ice%dw_dy_3D ( mesh%nV_loc, C%nz ))
    allocate( ice%dw_dz_3D ( mesh%nV_loc, C%nz ))
    ice%dw_dx_3D = 0d0
    ice%dw_dy_3D = 0d0
    ice%dw_dz_3D = 0d0

    ! === Basal conditions ===
    ! ========================

    ! Basal hydrology
    allocate( ice%pore_water_pressure ( mesh%nV_loc ))
    ice%pore_water_pressure = 0d0
    allocate( ice%overburden_pressure ( mesh%nV_loc ))
    ice%overburden_pressure = 0d0
    allocate( ice%effective_pressure  ( mesh%nV_loc ))
    ice%effective_pressure = 0d0

    ! Basal sliding
    allocate( ice%friction_coef_1    ( mesh%nV_loc ))
    ice%friction_coef_1 = 0d0
    allocate( ice%friction_coef_2    ( mesh%nV_loc ))
    ice%friction_coef_2 = 0d0
    allocate( ice%basal_shear_stress ( mesh%nV_loc ))
    ice%basal_shear_stress = 0d0

    ! Geothermal heat flux
    allocate( ice%geothermal_heat_flux ( mesh%nV_loc ))
    ice%geothermal_heat_flux = 0d0

    ! === Sea level ===
    ! =================

    ! Regional sea level
    allocate( ice%SL ( mesh%nV_loc ))
    ice%SL = 0d0

    ! === Area fractions ===
    ! ======================

    ! Grounded
    allocate( ice%fraction_gr   ( mesh%nV_loc   ))
    ice%fraction_gr = 0d0
    allocate( ice%fraction_gr_b ( mesh%nTri_loc ))
    ice%fraction_gr_b = 0d0

    ! Calving front
    allocate( ice%fraction_cf   ( mesh%nV_loc   ))
    ice%fraction_cf = 0d0

    ! === Basins ===
    ! ==============

    ! Basin ID's
    allocate( ice%basin_ID ( mesh%nV_loc ))
    ice%basin_ID = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_ice_model

END MODULE ice_model_memory
