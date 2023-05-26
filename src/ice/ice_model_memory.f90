MODULE ice_model_memory

  ! Routines for administrating the memory for the ice model data.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
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

    ! === Sea level ===
    ! =================

    ! Sea level
    allocate( ice%SL ( mesh%nV_loc ))
    ice%SL = 0d0

    ! === Ice change ===
    ! ==================

    ! Rates of change
    allocate( ice%dHi_dt     ( mesh%nV_loc      ))
    ice%dHi_dt = 0d0
    allocate( ice%dHb_dt     ( mesh%nV_loc      ))
    ice%dHb_dt = 0d0
    allocate( ice%dHs_dt     ( mesh%nV_loc      ))
    ice%dHs_dt = 0d0
    allocate( ice%dHib_dt    ( mesh%nV_loc      ))
    ice%dHib_dt = 0d0

    ! Predicted ice thickness at next step
    allocate( ice%Hi_tplusdt ( mesh%nV_loc      ))
    ice%Hi_tplusdt = 0d0

    ! ! Predictor/corrector method
    ! allocate( ice%pc_tau       ( mesh%nV_loc      ))
    ! ice%pc_tau = 0d0
    ! allocate( ice%pc_fcb       ( mesh%nV_loc      ))
    ! ice%pc_fcb = 0d0
    ! allocate( ice%pc_f1        ( mesh%nV_loc      ))
    ! ice%pc_f1 = 0d0
    ! allocate( ice%pc_f2        ( mesh%nV_loc      ))
    ! ice%pc_f2 = 0d0
    ! allocate( ice%pc_f3        ( mesh%nV_loc      ))
    ! ice%pc_f3 = 0d0
    ! allocate( ice%pc_f4        ( mesh%nV_loc      ))
    ! ice%pc_f4 = 0d0
    ! allocate( ice%Hi_old       ( mesh%nV_loc      ))
    ! ice%Hi_old = 0d0
    ! allocate( ice%Hi_pred      ( mesh%nV_loc      ))
    ! ice%Hi_pred = 0d0
    ! allocate( ice%Hi_corr      ( mesh%nV_loc      ))
    ! ice%Hi_corr = 0d0

    ! ! === Thermodynamics ===
    ! ! ======================

    ! ! Ice temperatures
    ! allocate( ice%Ti_a                 ( mesh%nV_loc, C%nz ))   ! Ice temperatures
    ! ice%Ti_pmp_a = 0d0
    ! allocate( ice%Ti_pmp_a             ( mesh%nV_loc, C%nz ))   ! Ice pressure melting point
    ! ice%Ti_pmp_a = 0d0

    ! ! Physical quantities
    ! allocate( ice%Cpi_a                ( mesh%nV_loc, C%nz ))   ! Ice specific heat capacity
    ! ice%Cpi_a = 0d0
    ! allocate( ice%Ki_a                 ( mesh%nV_loc, C%nz ))   ! Ice conductivity
    ! ice%Ki_a = 0d0

    ! ! Heating
    ! allocate( ice%internal_heating_a   ( mesh%nV_loc, C%nz ))
    ! ice%internal_heating_a = 0d0
    ! allocate( ice%frictional_heating_a ( mesh%nV_loc       ))
    ! ice%frictional_heating_a = 0d0

    ! ! Geothermal heat flux
    ! allocate( ice%geothermal_heat_flux                ( mesh%nV_loc       ))   ! Geothermal heat flux

    ! ! === Ice flow ===
    ! ! ================

    ! allocate( ice%A_flow_3D_a  ( mesh%nV_loc, C%nz ))   ! Ice flow factor
    ! ice%A_flow_3D_a = 0d0
    ! allocate( ice%A_flow_vav_a ( mesh%nV_loc       ))   ! Vertically integrated
    ! ice%A_flow_vav_a = 0d0
    ! ice%GHF_a = 0d0

    ! ! === Ice velocities ===
    ! ! ======================

    ! ! 3-D
    ! allocate( ice%u_3D_a ( mesh%nV_loc, C%nz ))
    ! ice%u_3D_a = 0d0
    ! allocate( ice%v_3D_a ( mesh%nV_loc, C%nz ))
    ! ice%v_3D_a = 0d0
    ! allocate( ice%u_3D_b ( nTri_loc, C%nz ))
    ! ice%u_3D_b = 0d0
    ! allocate( ice%v_3D_b ( nTri_loc, C%nz ))
    ! ice%v_3D_b = 0d0
    ! allocate( ice%w_3D_a ( mesh%nV_loc, C%nz ))
    ! ice%w_3D_a = 0d0

    ! ! Vertically integrated
    ! allocate( ice%u_vav_a    ( mesh%nV_loc ))
    ! ice%u_vav_a = 0d0
    ! allocate( ice%v_vav_a    ( mesh%nV_loc ))
    ! ice%v_vav_a = 0d0
    ! allocate( ice%u_vav_b    ( nTri_loc ))
    ! ice%u_vav_b = 0d0
    ! allocate( ice%v_vav_b    ( nTri_loc ))
    ! ice%v_vav_b = 0d0
    ! allocate( ice%uabs_vav_a ( mesh%nV_loc ))
    ! ice%uabs_vav_a = 0d0
    ! allocate( ice%uabs_vav_b ( nTri_loc ))
    ! ice%uabs_vav_b = 0d0

    ! ! Surface
    ! allocate( ice%u_surf_a    ( mesh%nV_loc ))
    ! ice%u_surf_a = 0d0
    ! allocate( ice%v_surf_a    ( mesh%nV_loc ))
    ! ice%v_surf_a = 0d0
    ! allocate( ice%u_surf_b    ( nTri_loc ))
    ! ice%u_surf_b = 0d0
    ! allocate( ice%v_surf_b    ( nTri_loc ))
    ! ice%v_surf_b = 0d0
    ! allocate( ice%w_surf_a    ( mesh%nV_loc ))
    ! ice%w_surf_a = 0d0
    ! allocate( ice%uabs_surf_a ( mesh%nV_loc ))
    ! ice%uabs_surf_a = 0d0
    ! allocate( ice%uabs_surf_b ( nTri_loc ))
    ! ice%uabs_surf_b = 0d0

    ! ! Basal
    ! allocate( ice%u_base_a    ( mesh%nV_loc ))
    ! ice%u_base_a = 0d0
    ! allocate( ice%v_base_a    ( mesh%nV_loc ))
    ! ice%v_base_a = 0d0
    ! allocate( ice%u_base_b    ( nTri_loc ))
    ! ice%u_base_b = 0d0
    ! allocate( ice%v_base_b    ( nTri_loc ))
    ! ice%v_base_b = 0d0
    ! allocate( ice%w_base_a    ( mesh%nV_loc ))
    ! ice%w_base_a = 0d0
    ! allocate( ice%uabs_base_a ( mesh%nV_loc ))
    ! ice%uabs_base_a = 0d0
    ! allocate( ice%uabs_base_b ( nTri_loc ))
    ! ice%uabs_base_b = 0d0

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

    ! ! === Grounded area fractions ===
    ! ! ===============================

    ! allocate( ice%f_grnd_a ( mesh%nV_loc ))
    ! allocate( ice%f_grnd_b ( nTri_loc ))

    ! ! === Zeta derivatives ===
    ! ! ========================

    ! allocate( ice%dzeta_dt_a ( mesh%nV_loc , C%nz ))
    ! ice%dzeta_dt_a = 0d0
    ! allocate( ice%dzeta_dx_a ( mesh%nV_loc , C%nz ))
    ! ice%dzeta_dx_a = 0d0
    ! allocate( ice%dzeta_dy_a ( mesh%nV_loc , C%nz ))
    ! ice%dzeta_dy_a = 0d0
    ! allocate( ice%dzeta_dz_a ( mesh%nV_loc        ))
    ! ice%dzeta_dz_a = 0d0

    ! ! === Calving ===
    ! ! ===============

    ! allocate( ice%float_margin_frac_a ( mesh%nV_loc ))
    ! ice%float_margin_frac_a = 0d0
    ! allocate( ice%Hi_eff_cf_a         ( mesh%nV_loc ))
    ! ice%Hi_eff_cf_a = 0d0

    ! ! === GIA ===
    ! ! ===========

    ! allocate( ice%dHb_dt_a ( mesh%nV_loc ))
    ! ice%dHb_dt_a = 0d0
    ! allocate( ice%dSL_dt_a ( mesh%nV_loc ))
    ! ice%dSL_dt_a = 0d0

    ! ! === Extras ===
    ! ! ==============

    ! allocate( ice%surf_curv    ( mesh%nV ))
    ! ice%surf_curv = 0d0
    ! allocate( ice%log_velocity ( mesh%nV_loc ))
    ! ice%log_velocity = 0d0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_ice_model

END MODULE ice_model_memory
