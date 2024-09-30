MODULE ice_model_main

  ! The main ice-dynamical model module.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF, &
                                                                     save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D, &
                                                                     save_variable_as_netcdf_logical_1D
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_pc
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE region_types                                           , ONLY: type_model_region
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE LMB_model_types                                        , ONLY: type_LMB_model
  USE AMB_model_types                                        , ONLY: type_AMB_model
  USE GIA_model_types                                        , ONLY: type_GIA_model
  USE ice_model_memory                                       , ONLY: allocate_ice_model
  USE ice_model_utilities                                    , ONLY: determine_masks, calc_bedrock_CDFs, calc_grounded_fractions, calc_zeta_gradients, &
                                                                     calc_mask_noice, alter_ice_thickness, initialise_bedrock_CDFs, MB_inversion, calc_effective_thickness, &
                                                                     initialise_dHi_dt_target, initialise_uabs_surf_target
  USE ice_thickness                                          , ONLY: calc_dHi_dt, apply_mask_noice_direct, apply_ice_thickness_BC_explicit
  USE math_utilities                                         , ONLY: ice_surface_elevation, thickness_above_floatation, Hi_from_Hb_Hs_and_SL, is_floating
  USE geothermal_heat_flux                                   , ONLY: initialise_geothermal_heat_flux
  USE basal_hydrology                                        , ONLY: initialise_basal_hydrology_model
  USE bed_roughness                                          , ONLY: initialise_bed_roughness
  USE ice_velocity_main                                      , ONLY: initialise_velocity_solver, solve_stress_balance, remap_velocity_solver, &
                                                                     create_restart_file_ice_velocity, write_to_restart_file_ice_velocity, &
                                                                     map_velocities_from_b_to_c_2D
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D, distribute_from_master_int_1D, &
                                                                     distribute_from_master_int_2D
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file, open_existing_netcdf_file_for_writing
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, add_time_dimension_to_file, &
                                                                     add_field_dp_0D, add_field_mesh_dp_2D, write_time_to_file, write_to_field_multopt_mesh_dp_2D, &
                                                                     write_to_field_multopt_dp_0D
  USE netcdf_input                                           , ONLY: read_field_from_file_0D, read_field_from_mesh_file_2D
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use remapping_main, only: Atlas, map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D, map_from_mesh_to_mesh_2D
  use mesh_data_smoothing, only: smooth_Gaussian_2D
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp
  USE petsc_basic                                            , ONLY: mat_petsc2CSR
  USE reallocate_mod                                         , ONLY: reallocate_bounds_dp_1D
  USE BMB_main                                               , ONLY: run_BMB_model
  USE mesh_operators                                         , ONLY: ddx_a_a_2D, ddy_a_a_2D
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ice_dynamics_model( region)
    ! Calculate ice geometry at the desired time, and update
    ! velocities, thinning rates, and predicted geometry if necessary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ice_dynamics_model'
    REAL(dp)                                              :: wt_prev, wt_next
    INTEGER                                               :: vi
    REAL(dp)                                              :: dt_max
    REAL(dp), DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: dHs_dx, dHs_dy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the current model time is at or beyond the point
    ! when the target dH/dt should be removed from the
    ! continuity equation, set its field to 0
    IF (region%time >= C%target_dHi_dt_t_end) THEN
      region%ice%dHi_dt_target = 0._dp
    END IF

    ! If the desired time is beyond the time of the next modelled ice thickness,
    ! run the ice dynamics model to calculate ice velocities, thinning rates,
    ! and a new next modelled ice thickness.
    ! ======================================

    IF (region%time == region%ice%t_Hi_next) THEN
      ! Need to calculate new ice velocities, thinning rates, and predicted ice thickness

      ! Start with the maximum allowed ice model time step
      dt_max = C%dt_ice_max

      ! Limit time step during the model start-up phase
      ! FIXME

      ! Run the ice dynamics model to calculate ice velocities, thinning rates,
      ! and a new next modelled ice thickness.
      IF     (C%choice_timestepping == 'direct') THEN
        CALL run_ice_dynamics_model_direct( region, dt_max)
      ELSEIF (C%choice_timestepping == 'pc') THEN
        CALL run_ice_dynamics_model_pc( region, dt_max)
      ELSE
        CALL crash('unknown choice_timestepping "' // TRIM( C%choice_timestepping) // '"!')
      END IF

    ELSEIF (region%time > region%ice%t_Hi_next) THEN
      ! This should not be possible
      CALL crash('overshot the ice dynamics time step')
    ELSE
      ! We're within the current ice dynamics prediction window
    END IF ! IF (region%time == region%ice%t_next) THEN

    ! Interpolate between previous and next modelled ice thickness
    ! to find the geometry at the desired time
    ! ========================================

    ! Calculate time interpolation weights
    wt_prev = (region%ice%t_Hi_next - region%time) / (region%ice%t_Hi_next - region%ice%t_Hi_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled ice thickness to desired time
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%Hi( vi) = wt_prev * region%ice%Hi_prev( vi) + wt_next * region%ice%Hi_next( vi)
    END DO

    ! Calculate all other ice geometry quantities
    ! ===========================================

    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Basic geometry
      region%ice%Hs ( vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
      region%ice%Hib( vi) = region%ice%Hs( vi) - region%ice%Hi( vi)
      region%ice%TAF( vi) = thickness_above_floatation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))

      ! Differences w.r.t. present-day
      region%ice%dHi ( vi)  = region%ice%Hi ( vi) - region%refgeo_PD%Hi ( vi)
      region%ice%dHb ( vi)  = region%ice%Hb ( vi) - region%refgeo_PD%Hb ( vi)
      region%ice%dHs ( vi)  = region%ice%Hs ( vi) - region%refgeo_PD%Hs ( vi)
      region%ice%dHib( vi)  = region%ice%Hib( vi) - (region%refgeo_PD%Hs ( vi) - region%refgeo_PD%Hi( vi))

      ! Rates of change
      region%ice%dHi_dt( vi) = (region%ice%Hi_next( vi) - region%ice%Hi_prev( vi)) / (region%ice%t_Hi_next - region%ice%t_Hi_prev)
      IF (region%ice%TAF( vi) > 0._dp) THEN
        ! Grounded ice
        region%ice%dHs_dt ( vi) = region%ice%dHb_dt( vi) + region%ice%dHi_dt( vi)
        region%ice%dHib_dt( vi) = region%ice%dHb_dt( vi)
      ELSE
        ! Floating ice
        region%ice%dHs_dt ( vi) = region%ice%dHi_dt( vi) * (1._dp - ice_density / seawater_density)
        region%ice%dHib_dt( vi) = region%ice%dHi_dt( vi) *          ice_density / seawater_density
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Update masks
    CALL determine_masks( region%mesh, region%ice)

    ! Calculate new effective thickness
    CALL calc_effective_thickness( region%mesh, region%ice, region%ice%Hi, region%ice%Hi_eff, region%ice%fraction_margin)

    ! Calculate absolute surface gradient
    CALL ddx_a_a_2D( region%mesh, region%ice%Hs, dHs_dx)
    CALL ddy_a_a_2D( region%mesh, region%ice%Hs, dHs_dy)
    region%ice%Hs_slope = SQRT( dHs_dx**2 + dHs_dy**2)

    ! NOTE: as calculating the zeta gradients is quite expensive, only do so when necessary,
    !       i.e. when solving the heat equation or the Blatter-Pattyn stress balance
    ! Calculate zeta gradients
    CALL calc_zeta_gradients( region%mesh, region%ice)

    ! Calculate sub-grid grounded-area fractions
    CALL calc_grounded_fractions( region%mesh, region%ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_model

  SUBROUTINE initialise_ice_dynamics_model( mesh, ice, refgeo_init, refgeo_PD, refgeo_GIAeq, GIA, region_name)
    ! Initialise all data fields of the ice module

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo_init
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo_PD
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo_GIAeq
    TYPE(type_GIA_model),                   INTENT(IN)    :: GIA
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ice_dynamics_model'
    INTEGER                                               :: vi
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: dHs_dx, dHs_dy

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      WRITE(*,"(A)") '   Initialising ice dynamics model...'
    END IF
    CALL sync

    ! === Memory allocation ===
    ! =========================

    ! Allocate memory
    CALL allocate_ice_model( mesh, ice)

    ! === Value initialisation ===
    ! ============================

    ! Sea level
    ! =========

    SELECT CASE (C%choice_sealevel_model)

      CASE ('fixed')
        ! Fixed sea level
        ice%SL = C%fixed_sealevel

      CASE ('prescribed')
        ! Sea-level prescribed from external record file
        CALL crash('Sea level initialisation: prescribed method not implement yet!')
        ! ice%SL = forcing%sealevel_obs

      CASE ('eustatic')
        ! Eustatic sea level
        CALL crash('Sea level initialisation: eustatic method not implement yet!')
        ! ice%SL = C%initial_guess_sealevel

      CASE ('SELEN')
        ! Sea level from SELEN
        CALL crash('Sea level initialisation: SELEN method not implement yet!')
        ! ice%SL = C%initial_guess_sealevel

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_sealevel_model "' // &
                    TRIM( C%choice_sealevel_model) // '"!')

    END SELECT

    ! Initialise ice geometry
    ! =======================

    ! ! DENK DROM
    ! IF (par%master) CALL warning('GIA model isnt finished yet - need to include dHb in ice model initialisation')

    ! Basic geometry
    DO vi = mesh%vi1, mesh%vi2
      ice%Hb( vi) = refgeo_GIAeq%Hb( vi)
      ice%Hs( vi) = refgeo_init%Hs ( vi)
      ice%Hi( vi) = Hi_from_Hb_Hs_and_SL( ice%Hb( vi), ice%Hs( vi), ice%SL( vi))
    END DO

    ! Calculate the no-ice mask
    CALL calc_mask_noice( mesh, ice)

    ! Apply no-ice mask
    CALL apply_mask_noice_direct( mesh, ice%mask_noice, ice%Hi)

    ! Apply boundary conditions at the domain border
    CALL apply_ice_thickness_BC_explicit( mesh, ice%mask_noice, ice%Hb, ice%SL, ice%Hi)

    DO vi = mesh%vi1, mesh%vi2

      ! Derived geometry
      ice%Hs ( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))
      ice%Hib( vi) = ice%Hs( vi) - ice%Hi( vi)
      ice%TAF( vi) = thickness_above_floatation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))

      ! Differences w.r.t. present-day
      ice%dHi ( vi)  = ice%Hi ( vi) - refgeo_PD%Hi ( vi)
      ice%dHb ( vi)  = ice%Hb ( vi) - refgeo_PD%Hb ( vi)
      ice%dHs ( vi)  = ice%Hs ( vi) - refgeo_PD%Hs ( vi)
      ice%dHib( vi)  = ice%Hib( vi) - (refgeo_PD%Hs ( vi) - refgeo_PD%Hi( vi))

      ! Rates of change
      ice%dHi_dt ( vi) = 0._dp
      ice%dHb_dt ( vi) = 0._dp
      ice%dHs_dt ( vi) = 0._dp
      ice%dHib_dt( vi) = 0._dp

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Calculate zeta gradients
    CALL calc_zeta_gradients( mesh, ice)

    ! Model states for ice dynamics model
    ice%t_Hi_prev = C%start_time_of_run
    ice%t_Hi_next = C%start_time_of_run
    ice%Hi_prev   = ice%Hi
    ice%Hi_next   = ice%Hi

    ! Initialise masks
    ! ================

    ! Call it twice so also the "prev" versions are set
    CALL determine_masks( mesh, ice)
    CALL determine_masks( mesh, ice)

    ! Effective ice thickness
    ! =======================

    ! Compute effective thickness at calving fronts
    CALL calc_effective_thickness( mesh, ice, ice%Hi, ice%Hi_eff, ice%fraction_margin)

    ! Surface gradients
    ! =================

    ! Calculate absolute surface gradient
    CALL ddx_a_a_2D( mesh, ice%Hs, dHs_dx)
    CALL ddy_a_a_2D( mesh, ice%Hs, dHs_dy)
    ice%Hs_slope = SQRT( dHs_dx**2 + dHs_dy**2)

    ! Target thinning rates
    ! =====================

    ! Load target dHi_dt for inversions
    IF (C%do_target_dHi_dt) THEN
      CALL initialise_dHi_dt_target(mesh, ice, region_name)
    ELSE
      ice%dHi_dt_target = 0._dp
    END IF

    ! Target surface ice speed
    ! ========================

    ! Load target dHi_dt for inversions
    IF (C%do_target_uabs_surf) THEN
      CALL initialise_uabs_surf_target(mesh, ice, region_name)
    ELSE
      ice%uabs_surf_target = 0._dp
    END IF

    ! Sub-grid fractions
    ! ==================

    ! Initialise bedrock cumulative density functions
    CALL initialise_bedrock_CDFs( mesh, refgeo_PD, ice, region_name)
    ! Initialise sub-grid grounded-area fractions
    CALL calc_grounded_fractions( mesh, ice)

    ! Basal conditions
    ! ================

    ! Allocate and initialise basal conditions
    CALL initialise_geothermal_heat_flux(  mesh, ice)
    CALL initialise_basal_hydrology_model( mesh, ice, region_name)
    CALL initialise_bed_roughness(         mesh, ice, region_name)

    ! Velocities
    ! ==========

    ! Initialise data for the chosen velocity solver(s)
    CALL initialise_velocity_solver( mesh, ice, region_name)

    ! Time stepping
    ! =============

    IF     (C%choice_timestepping == 'direct') THEN
      ! No need to initialise anything here
    ELSEIF (C%choice_timestepping == 'pc') THEN
      CALL initialise_pc_scheme( mesh, ice%pc, region_name)
    ELSE
      CALL crash('unknown choice_timestepping "' // TRIM( C%choice_timestepping) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_dynamics_model

  SUBROUTINE write_to_restart_files_ice_model( mesh, ice, time)
    ! Write to all the restart files for the ice dynamics model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_files_ice_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! First for the velocity solver
    CALL write_to_restart_file_ice_velocity( mesh, ice, time)

    ! Then for the time-stepper
    IF     (C%choice_timestepping == 'direct') THEN
      ! Direct time stepping doesn't require a restart file
    ELSEIF (C%choice_timestepping == 'pc') THEN
      CALL write_to_restart_file_pc_scheme( mesh, ice%pc, time)
    ELSE
      CALL crash('unknown choice_timestepping "' // TRIM( C%choice_timestepping) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_files_ice_model

  SUBROUTINE create_restart_files_ice_model( mesh, ice)
    ! Create all the restart files for the ice dynamics model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_files_ice_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! First for the velocity solver
    CALL create_restart_file_ice_velocity( mesh, ice)

    ! Then for the time-stepper
    IF     (C%choice_timestepping == 'direct') THEN
      ! Direct time stepping doesn't require a restart file
    ELSEIF (C%choice_timestepping == 'pc') THEN
      CALL create_restart_file_pc_scheme( mesh, ice%pc)
    ELSE
      CALL crash('unknown choice_timestepping "' // TRIM( C%choice_timestepping) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_files_ice_model

  SUBROUTINE remap_ice_dynamics_model( mesh_old, mesh_new, ice, refgeo_PD, SMB, BMB, LMB, AMB, GIA, time, region_name)
    ! Remap/reallocate all the data of the ice dynamics model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(INOUT) :: mesh_new
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo_PD
    TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                   INTENT(IN)    :: BMB
    TYPE(type_LMB_model),                   INTENT(IN)    :: LMB
    TYPE(type_AMB_model),                   INTENT(IN)    :: AMB
    TYPE(type_GIA_model),                   INTENT(IN)    :: GIA
    REAL(dp),                               INTENT(IN)    :: time
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_ice_dynamics_model'
    INTEGER                                               :: vi,k
    REAL(dp), DIMENSION(mesh_new%vi1:mesh_new%vi2)        :: dHs_dx, dHs_dy
    REAL(dp)                                              :: Ti_min

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '    Remapping ice model data to the new mesh...'

    ! Remap conserved ice model quantities
    ! ====================================

    ! === Ice-sheet geometry ===
    ! ==========================

    ! Remap basic ice geometry Hi,Hb,Hs,SL
    CALL remap_basic_ice_geometry( mesh_old, mesh_new, refgeo_PD, GIA, ice)

    ! Remap dHi/dt to improve stability of the P/C scheme after mesh updates
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, ice%dHi_dt, '2nd_order_conservative')

    ! === Thermodynamics and rheology ===
    ! ===================================

    ! Save minimum temperature from the entire Ti field, so we can
    ! replace any 0s that pop up after the remapping with that.
    ! Otherwise, divisions by 0 might occur during the computation
    ! of the ice flow factor A.
    Ti_min = minval(ice%Ti)

    ! Use 2nd-order conservative remapping for the ice temperature.
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, ice%Ti, '2nd_order_conservative')

    ! Make sure that no values are smaller than the original minimum
    ice%Ti = MAX( ice%Ti, Ti_min)

    ! Predicted model state at next time step
    CALL reallocate_bounds( ice%Ti_prev, mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [K]  The previous state
    CALL reallocate_bounds( ice%Ti_next, mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [K]  The next state

    ! Re-initialise
    ice%Ti_prev = ice%Ti
    ice%Ti_next = ice%Ti

    ! Reallocate memory for all other data fields
    ! (list copied from ice_model_memory/allocate_ice_model)
    ! ======================================================

    ! === Ice-sheet geometry ===
    ! ==========================

    ! Basic geometry
    ! CALL reallocate_bounds( ice%Hi                          , mesh_new%vi1, mesh_new%vi2         )  ! [m] Ice thickness
    ! CALL reallocate_bounds( ice%Hb                          , mesh_new%vi1, mesh_new%vi2         )  ! [m] Bedrock elevation (w.r.t. PD sea level)
    ! CALL reallocate_bounds( ice%Hs                          , mesh_new%vi1, mesh_new%vi2         )  ! [m] Surface elevation (w.r.t. PD sea level)
    ! CALL reallocate_bounds( ice%SL                          , mesh_new%vi1, mesh_new%vi2         )  ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    CALL reallocate_bounds( ice%Hib                         , mesh_new%vi1, mesh_new%vi2         )  ! [m] Ice base elevation (w.r.t. PD sea level)
    CALL reallocate_bounds( ice%TAF                         , mesh_new%vi1, mesh_new%vi2         )  ! [m] Thickness above flotation
    CALL reallocate_bounds( ice%Hi_eff                      , mesh_new%vi1, mesh_new%vi2         )  ! [m] Effective ice thickness
    CALL reallocate_bounds( ice%Hs_slope                    , mesh_new%vi1, mesh_new%vi2         )  ! [-] Absolute surface gradients

    ! Geometry changes
    CALL reallocate_bounds( ice%dHi                         , mesh_new%vi1, mesh_new%vi2         )  ! [m] Ice thickness difference (w.r.t. reference)
    CALL reallocate_bounds( ice%dHb                         , mesh_new%vi1, mesh_new%vi2         )  ! [m] Bedrock elevation difference (w.r.t. reference)
    CALL reallocate_bounds( ice%dHs                         , mesh_new%vi1, mesh_new%vi2         )  ! [m] Surface elevation difference (w.r.t. reference)
    CALL reallocate_bounds( ice%dHib                        , mesh_new%vi1, mesh_new%vi2         )  ! [m] Base elevation difference (w.r.t. reference)

    ! Rates of change
    ! CALL reallocate_bounds( ice%dHi_dt                      , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Ice thickness rate of change
    CALL reallocate_bounds( ice%dHb_dt                      , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Bedrock elevation rate of change
    CALL reallocate_bounds( ice%dHs_dt                      , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Ice surface elevation rate of change
    CALL reallocate_bounds( ice%dHib_dt                     , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Ice base elevation rate of change
    CALL reallocate_bounds( ice%dHi_dt_raw                  , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Ice thickness rate of change before any imposed modifications
    CALL reallocate_bounds( ice%dHi_dt_residual             , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Residual ice thickness rate of change after imposed modifications

    ! Target quantities
    CALL reallocate_bounds( ice%dHi_dt_target               , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Target ice thickness rate of change for inversions
    CALL reallocate_bounds( ice%uabs_surf_target            , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Target ice surface speed for inversions

    ! Masks
    CALL reallocate_bounds( ice%mask_icefree_land           , mesh_new%vi1, mesh_new%vi2         )  ! T: ice-free land , F: otherwise
    CALL reallocate_bounds( ice%mask_icefree_ocean          , mesh_new%vi1, mesh_new%vi2         )  ! T: ice-free ocean, F: otherwise
    CALL reallocate_bounds( ice%mask_grounded_ice           , mesh_new%vi1, mesh_new%vi2         )  ! T: grounded ice  , F: otherwise
    CALL reallocate_bounds( ice%mask_floating_ice           , mesh_new%vi1, mesh_new%vi2         )  ! T: floating ice  , F: otherwise
    CALL reallocate_bounds( ice%mask_icefree_land_prev      , mesh_new%vi1, mesh_new%vi2         )  ! T: ice-free land , F: otherwise (during previous time step)
    CALL reallocate_bounds( ice%mask_icefree_ocean_prev     , mesh_new%vi1, mesh_new%vi2         )  ! T: ice-free ocean, F: otherwise (during previous time step)
    CALL reallocate_bounds( ice%mask_grounded_ice_prev      , mesh_new%vi1, mesh_new%vi2         )  ! T: grounded ice  , F: otherwise (during previous time step)
    CALL reallocate_bounds( ice%mask_floating_ice_prev      , mesh_new%vi1, mesh_new%vi2         )  ! T: floating ice  , F: otherwise (during previous time step)
    CALL reallocate_bounds( ice%mask_margin                 , mesh_new%vi1, mesh_new%vi2         )  ! T: ice next to ice-free, F: otherwise
    CALL reallocate_bounds( ice%mask_gl_gr                  , mesh_new%vi1, mesh_new%vi2         )  ! T: grounded ice next to floating ice, F: otherwise
    CALL reallocate_bounds( ice%mask_gl_fl                  , mesh_new%vi1, mesh_new%vi2         )  ! T: floating ice next to grounded ice, F: otherwise
    CALL reallocate_bounds( ice%mask_cf_gr                  , mesh_new%vi1, mesh_new%vi2         )  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    CALL reallocate_bounds( ice%mask_cf_fl                  , mesh_new%vi1, mesh_new%vi2         )  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    CALL reallocate_bounds( ice%mask_coastline              , mesh_new%vi1, mesh_new%vi2         )  ! T: ice-free land next to ice-free ocean, F: otherwise
    ! CALL reallocate_bounds( ice%mask_noice                  , mesh_new%vi1, mesh_new%vi2         )  ! T: no ice is allowed here, F: ice is allowed here
    CALL reallocate_bounds( ice%mask                        , mesh_new%vi1, mesh_new%vi2         )  ! Diagnostic, only meant for quick visual inspection in output
    CALL reallocate_bounds( ice%basin_ID                    , mesh_new%vi1, mesh_new%vi2         )  ! The drainage basin to which each vertex belongs

    ! Area fractions
    CALL reallocate_bounds( ice%fraction_gr                 , mesh_new%vi1, mesh_new%vi2         )  ! [0-1] Grounded area fractions of vertices
    CALL reallocate_bounds( ice%fraction_gr_b               , mesh_new%ti1, mesh_new%ti2         )  ! [0-1] Grounded area fractions of triangles
    CALL reallocate_bounds( ice%fraction_margin             , mesh_new%vi1, mesh_new%vi2         )  ! [0-1] Ice-covered area fractions of ice margins

    ! Sub-grid bedrock cumulative density functions (CDFs)
    CALL reallocate_bounds( ice%bedrock_cdf                 , mesh_new%vi1, mesh_new%vi2, C%subgrid_bedrock_cdf_nbins)  ! [-] Sub-grid bedrock cumulative density functions on the a-grid (vertices)
    CALL reallocate_bounds( ice%bedrock_cdf_b               , mesh_new%ti1, mesh_new%ti2, C%subgrid_bedrock_cdf_nbins)  ! [-] Sub-grid bedrock cumulative density functions on the b-grid (triangles)

    ! === Terrain-following coordinate zeta gradients ===
    ! ===================================================

    ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta

    ! On the ak-grid (vertices, vertically regular)
    CALL reallocate_bounds( ice%dzeta_dt_ak                 , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dzeta_dx_ak                 , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dzeta_dy_ak                 , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dzeta_dz_ak                 , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%d2zeta_dx2_ak               , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%d2zeta_dxdy_ak              , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%d2zeta_dy2_ak               , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)

    ! On the bk-grid (triangles, vertically regular)
    CALL reallocate_bounds( ice%dzeta_dx_bk                 , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( ice%dzeta_dy_bk                 , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( ice%dzeta_dz_bk                 , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( ice%d2zeta_dx2_bk               , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( ice%d2zeta_dxdy_bk              , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( ice%d2zeta_dy2_bk               , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! On the bks-grid (triangles, vertically staggered)
    CALL reallocate_bounds( ice%dzeta_dx_bks                , mesh_new%ti1, mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( ice%dzeta_dy_bks                , mesh_new%ti1, mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( ice%dzeta_dz_bks                , mesh_new%ti1, mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( ice%d2zeta_dx2_bks              , mesh_new%ti1, mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( ice%d2zeta_dxdy_bks             , mesh_new%ti1, mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( ice%d2zeta_dy2_bks              , mesh_new%ti1, mesh_new%ti2, mesh_new%nz-1)

    ! === Thermodynamics and rheology ===
    ! ===================================

    ! Ice temperatures
    ! CALL reallocate_bounds( ice%Ti                          , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [K] Englacial temperature
    CALL reallocate_bounds( ice%Ti_pmp                      , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [K] Pressure melting point temperature
    CALL reallocate_bounds( ice%Ti_hom                      , mesh_new%vi1, mesh_new%vi2             )  ! [K] Basal temperature w.r.t. pressure melting point

    ! Physical quantities
    CALL reallocate_bounds( ice%Cpi                         , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [J kg^-1 K^-1] Specific heat capacity
    CALL reallocate_bounds( ice%Ki                          , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [J m^-1 K^-1 yr^-1] Thermal conductivity

    ! Heating
    CALL reallocate_bounds( ice%internal_heating            , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [?] Internal heating
    CALL reallocate_bounds( ice%frictional_heating          , mesh_new%vi1, mesh_new%vi2             )  ! [?] Frictional heating

    ! Glen's flow law factor
    CALL reallocate_bounds( ice%A_flow                      , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [Pa^-3 y^-1] Glen's flow law factor

    ! === Ice velocities ===
    ! ======================

    ! 3-D
    CALL reallocate_bounds( ice%u_3D                        , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [m yr^-1] 3-D ice velocity
    CALL reallocate_bounds( ice%v_3D                        , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%u_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( ice%v_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( ice%w_3D                        , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)

    ! Vertically integrated
    CALL reallocate_bounds( ice%u_vav                       , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Vertically averaged ice velocity
    CALL reallocate_bounds( ice%v_vav                       , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%u_vav_b                     , mesh_new%ti1, mesh_new%ti2         )
    CALL reallocate_bounds( ice%v_vav_b                     , mesh_new%ti1, mesh_new%ti2         )
    CALL reallocate_bounds( ice%uabs_vav                    , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%uabs_vav_b                  , mesh_new%ti1, mesh_new%ti2         )

    ! Surface
    CALL reallocate_bounds( ice%u_surf                      , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Ice velocity at the surface
    CALL reallocate_bounds( ice%v_surf                      , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%u_surf_b                    , mesh_new%ti1, mesh_new%ti2         )
    CALL reallocate_bounds( ice%v_surf_b                    , mesh_new%ti1, mesh_new%ti2         )
    CALL reallocate_bounds( ice%w_surf                      , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%uabs_surf                   , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%uabs_surf_b                 , mesh_new%ti1, mesh_new%ti2         )

    ! Basal
    CALL reallocate_bounds( ice%u_base                      , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Ice velocity at the base
    CALL reallocate_bounds( ice%v_base                      , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%u_base_b                    , mesh_new%ti1, mesh_new%ti2         )
    CALL reallocate_bounds( ice%v_base_b                    , mesh_new%ti1, mesh_new%ti2         )
    CALL reallocate_bounds( ice%w_base                      , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%uabs_base                   , mesh_new%vi1, mesh_new%vi2         )
    CALL reallocate_bounds( ice%uabs_base_b                 , mesh_new%ti1, mesh_new%ti2         )

    ! == Strain rates ==
    ! ==================

    CALL reallocate_bounds( ice%du_dx_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)  ! [yr^-1]
    CALL reallocate_bounds( ice%du_dy_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%du_dz_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dv_dx_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dv_dy_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dv_dz_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dw_dx_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dw_dy_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( ice%dw_dz_3D                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)

    ! == Ice flow regime ==
    ! =====================

    CALL reallocate_bounds( ice%divQ                        , mesh_new%vi1, mesh_new%vi2         )  ! [m yr^-1] Horizontal ice flux divergence
    CALL reallocate_bounds( ice%R_shear                     , mesh_new%vi1, mesh_new%vi2         )  ! [0-1]     uabs_base / uabs_surf (0 = pure vertical shear, viscous flow; 1 = pure sliding, plug flow)

    ! == Basal hydrology ==
    ! =====================

    ! Basal hydrology
    CALL reallocate_bounds( ice%pore_water_pressure         , mesh_new%vi1, mesh_new%vi2         )  ! [Pa]  Basal pore water pressure
    CALL reallocate_bounds( ice%overburden_pressure         , mesh_new%vi1, mesh_new%vi2         )  ! [Pa]  Basal overburden pressure
    CALL reallocate_bounds( ice%effective_pressure          , mesh_new%vi1, mesh_new%vi2         )  ! [Pa]  Basal effective pressure
    CALL reallocate_bounds( ice%pore_water_likelihood       , mesh_new%vi1, mesh_new%vi2         )  ! [0-1] Basal pore water likelihood
    CALL reallocate_bounds( ice%pore_water_fraction         , mesh_new%vi1, mesh_new%vi2         )  ! [0-1] Fraction of overburden pressure reduced by pore water pressure

    ! == Basal sliding ==
    ! ===================

    ! Sliding law coefficients
    CALL reallocate_bounds( ice%till_friction_angle         , mesh_new%vi1, mesh_new%vi2         )  ! [degrees]          Till friction angle (degrees)
    CALL reallocate_bounds( ice%bed_roughness               , mesh_new%vi1, mesh_new%vi2         )  ! [0-1]              Bed roughness fraction
    CALL reallocate_bounds( ice%till_yield_stress           , mesh_new%vi1, mesh_new%vi2         )  ! [Pa]               Till yield stress (used when choice_sliding_law = "Coloumb", "Budd", or "Zoet-Iverson")
    CALL reallocate_bounds( ice%slid_alpha_sq               , mesh_new%vi1, mesh_new%vi2         )  ! [-]                Coulomb-law friction coefficient (used when choice_sliding_law = "Tsai2015", or "Schoof2005")
    CALL reallocate_bounds( ice%slid_beta_sq                , mesh_new%vi1, mesh_new%vi2         )  ! [Pa m^−1/m yr^1/m] Power-law friction coefficient (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")

    ! Basal friction and shear stress
    CALL reallocate_bounds( ice%basal_friction_coefficient  , mesh_new%vi1, mesh_new%vi2         )  ! [Pa yr m^-1]       Effective basal friction coefficient (basal_shear_stress = u_base * basal_friction_coefficient)
    CALL reallocate_bounds( ice%basal_shear_stress          , mesh_new%vi1, mesh_new%vi2         )  ! [Pa]               Basal shear stress

    ! == Geothermal heat ==
    ! =====================

    CALL reallocate_bounds( ice%geothermal_heat_flux        , mesh_new%vi1, mesh_new%vi2         )  ! [J m^-2 yr^-1] Geothermal heat flux

    ! === Ice thickness time stepping ===
    ! ===================================

    ! Predicted model state at next time step
    CALL reallocate_bounds( ice%Hi_prev                     , mesh_new%vi1, mesh_new%vi2         )  ! [m]  The previous state
    CALL reallocate_bounds( ice%Hi_next                     , mesh_new%vi1, mesh_new%vi2         )  ! [m]  The next state

    ! Re-initialise the rest of the ice dynamics model
    ! ================================================

    ! Initialise ice geometry
    ! =======================

    DO vi = mesh_new%vi1, mesh_new%vi2

      ! Basic geometry
      ! ice%Hi ( vi) = refgeo_init%Hi( vi)
      ! ice%Hb ( vi) = refgeo_init%Hb( vi)
      ice%Hs ( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))
      ice%Hib( vi) = ice%Hs( vi) - ice%Hi( vi)
      ice%TAF( vi) = thickness_above_floatation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))

      ! Differences w.r.t. present-day
      ice%dHi ( vi)  = ice%Hi ( vi) - refgeo_PD%Hi ( vi)
      ice%dHb ( vi)  = ice%Hb ( vi) - refgeo_PD%Hb ( vi)
      ice%dHs ( vi)  = ice%Hs ( vi) - refgeo_PD%Hs ( vi)
      ice%dHib( vi)  = ice%Hib( vi) - (refgeo_PD%Hs ( vi) - refgeo_PD%Hi( vi))

      ! Rates of change
      ice%dHi_dt ( vi) = 0._dp
      ice%dHb_dt ( vi) = 0._dp
      ice%dHs_dt ( vi) = 0._dp
      ice%dHib_dt( vi) = 0._dp

    END DO ! DO vi = mesh_new%vi1, mesh_new%vi2

    ! Calculate zeta gradients
    CALL calc_zeta_gradients( mesh_new, ice)

    ! Load target dHi_dt for inversions
    IF (C%do_target_dHi_dt) THEN
      CALL initialise_dHi_dt_target(mesh_new, ice, region_name)
    ELSE
      ice%dHi_dt_target = 0._dp
    END IF

    ! Load target surface ice speed for inversions
    IF (C%do_target_uabs_surf) THEN
      CALL initialise_uabs_surf_target(mesh_new, ice, region_name)
    ELSE
      ice%uabs_surf_target = 0._dp
    END IF

    ! Model states for ice dynamics model
    ice%t_Hi_prev = time
    ice%t_Hi_next = time
    ice%Hi_prev   = ice%Hi
    ice%Hi_next   = ice%Hi

    ! Initialise masks
    ! ================

    ! Calculate the no-ice mask
    CALL calc_mask_noice( mesh_new, ice)

    ! Remove ice bleed into forbidden areas
    CALL apply_mask_noice_direct( mesh_new, ice%mask_noice, ice%Hi)
    CALL apply_mask_noice_direct( mesh_new, ice%mask_noice, ice%dHi_dt)

    ! Call it twice so also the "prev" versions are set
    CALL determine_masks( mesh_new, ice)
    CALL determine_masks( mesh_new, ice)

    ! ! Smooth the ice at the calving front to improve model stability
    ! CALL relax_calving_front_after_mesh_update( mesh_new, ice)

    ! Effective ice thickness
    ! =======================

    ! Calculate new effective thickness
    CALL calc_effective_thickness( mesh_new, ice, ice%Hi, ice%Hi_eff, ice%fraction_margin)

    ! Surface gradients
    ! =================

    ! Calculate absolute surface gradient
    CALL ddx_a_a_2D( mesh_new, ice%Hs, dHs_dx)
    CALL ddy_a_a_2D( mesh_new, ice%Hs, dHs_dy)
    ice%Hs_slope = SQRT( dHs_dx**2 + dHs_dy**2)

    ! Sub-grid fractions
    ! ==================

    IF (C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .OR. &
        C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') THEN
      ! Compute bedrock cumulative density function
      CALL calc_bedrock_CDFs( mesh_new, refgeo_PD, ice)
    END IF
    ! Initialise sub-grid grounded-area fractions
    CALL calc_grounded_fractions( mesh_new, ice)

    ! Basal conditions
    ! ================

    ! Allocate and initialise basal conditions
    CALL initialise_geothermal_heat_flux(  mesh_new, ice)
    CALL initialise_basal_hydrology_model( mesh_new, ice, region_name)

    ! FIXME: something should happen here once we start working on remapping of inverted bed roughness!
    CALL initialise_bed_roughness(         mesh_new, ice, region_name)

    ! Velocities
    ! ==========

    ! Remap data for the chosen velocity solver(s)
    CALL remap_velocity_solver( mesh_old, mesh_new, ice)

    ! Time stepping
    ! =============

    IF     (C%choice_timestepping == 'direct') THEN
      ! No need to initialise anything here
    ELSEIF (C%choice_timestepping == 'pc') THEN
      CALL remap_pc_scheme( mesh_old, mesh_new, ice%pc)
    ELSE
      CALL crash('unknown choice_timestepping "' // TRIM( C%choice_timestepping) // '"!')
    END IF

    ! Relax ice geometry around the calving front
    ! ===========================================

    CALL relax_calving_front( mesh_old, mesh_new, ice, SMB, BMB, LMB, AMB, region_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ice_dynamics_model

  SUBROUTINE remap_basic_ice_geometry( mesh_old, mesh_new, refgeo_PD, GIA, ice)
    ! Remap the basic ice geometry Hi,Hb,Hs,SL.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo_PD
    TYPE(type_GIA_model),                   INTENT(IN)    :: GIA
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_basic_ice_geometry'
    REAL(dp), DIMENSION( mesh_old%nV)                     :: Hi_old_tot
    LOGICAL,  DIMENSION( mesh_old%nV)                     :: mask_floating_ice_tot
    LOGICAL,  DIMENSION( mesh_old%nV)                     :: mask_icefree_ocean_tot
    REAL(dp), DIMENSION( mesh_new%vi1:mesh_new%vi2)       :: Hi_new
    REAL(dp), DIMENSION( mesh_new%vi1:mesh_new%vi2)       :: Hs_new
    INTEGER                                               :: vi
    INTEGER                                               :: mi, mi_used
    LOGICAL                                               :: found_map
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_CSR
    INTEGER                                               :: vi_new, k1, k2, k, vi_old
    INTEGER                                               :: n_ice, n_nonice
    INTEGER                                               :: n_shelf, n_open_ocean
    REAL(dp)                                              :: sum_Hi_shelf

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Basic: remap surface elevation Hs from the old mesh, remap bedrock elevation Hb
    !    from its (presumably high-resolution) source grid, define remapped ice thickness
    !    as the difference between the two. As surface elevation is typically much smoother
    !    then ice thickness, remapping works much better.
    ! =====================================================================================

    ! Remap bedrock from the original high-resolution grid, and add the (very smooth) modelled deformation to it
    ! Remapping of Hb in the refgeo structure has already happened, only need to copy the data
    IF (par%master) CALL warning('GIA model isnt finished yet - need to include dHb in mesh update!')
    CALL reallocate_bounds( ice%Hb, mesh_new%vi1, mesh_new%vi2)  ! [m] Bedrock elevation (w.r.t. PD sea level)
    ice%Hb = refgeo_PD%Hb

    ! Remap sea level
    IF (par%master) CALL warning('sea  model isnt finished yet - need to include dSL in mesh update!')
    CALL reallocate_bounds( ice%SL, mesh_new%vi1, mesh_new%vi2)  ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    ice%SL = 0._dp

    ! Gather global ice thickness and masks
    CALL gather_to_all_dp_1D(      ice%Hi                , Hi_old_tot            )
    CALL gather_to_all_logical_1D( ice%mask_floating_ice , mask_floating_ice_tot )
    CALL gather_to_all_logical_1D( ice%mask_icefree_ocean, mask_icefree_ocean_tot)

    ! First, naively remap ice thickness and surface elevation without any restrictions
    CALL map_from_mesh_to_mesh_2D( mesh_old, mesh_new, ice%Hi, Hi_new, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_2D( mesh_old, mesh_new, ice%Hs, Hs_new, '2nd_order_conservative')

    ! Calculate remapped ice thickness as the difference between new bedrock and remapped surface elevation
    DO vi = mesh_new%vi1, mesh_new%vi2
      IF (Hi_new( vi) > 0._dp) THEN
        IF (Hs_new( vi) <= ice%Hb( vi)) THEN
          Hi_new( vi) = 0._dp
        ELSE
          Hi_new( vi) = Hi_from_Hb_Hs_and_SL( ice%Hb( vi), Hs_new( vi), ice%SL( vi))
        END IF
      ELSE
        Hi_new( vi) = 0._dp
      END IF
    END DO ! DO vi = mesh_new%vi1, mesh_new%vi2

    ! Reallocate no-ice mask
    ! T: no ice is allowed here, F: ice is allowed here
    CALL reallocate_bounds( ice%mask_noice, mesh_new%vi1, mesh_new%vi2)

    ! Apply boundary conditions at the domain border
    CALL calc_mask_noice( mesh_new, ice)
    CALL apply_ice_thickness_BC_explicit( mesh_new, ice%mask_noice, ice%Hb, ice%SL, Hi_new)

    ! == Corrections
    ! ==============

    ! Browse the Atlas to find the map between mesh_old and mesh_new
    found_map = .FALSE.
    mi_used   = 0
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == mesh_old%name .AND. Atlas( mi)%name_dst == mesh_new%name &
          .AND. Atlas( mi)%method == '2nd_order_conservative') THEN
        found_map = .TRUE.
        mi_used  = mi
        EXIT
      END IF
    END DO
    ! Safety
    IF (.NOT. found_map) CALL crash('couldnt find which map was used')

    ! Convert the mapping matrix to CSR format
    CALL mat_petsc2CSR( Atlas( mi_used)%M, M_CSR)

    ! == For those vertices of the new mesh that overlap with both old-mesh ice and old-mesh
    !    non-ice, remove very thin remapped ice
    ! ======================================================================================

    DO vi_new = mesh_new%vi1, mesh_new%vi2

      k1 = M_CSR%ptr( vi_new)
      k2 = M_CSR%ptr( vi_new+1) - 1

      n_ice    = 0
      n_nonice = 0

      DO k = k1, k2

        vi_old = M_CSR%ind( k)

        IF     (Hi_old_tot( vi_old) > 1._dp) THEN
          n_ice = n_ice + 1
        ELSEIF (Hi_old_tot( vi_old) < 1._dp) THEN
          n_nonice = n_nonice + 1
        END IF

      END DO ! DO k = k1, k2

      IF (n_ice > 0 .AND. n_nonice > 0) THEN
        ! This new-mesh vertex overlaps with both old-mesh ice vertices,
        ! and old-mesh non-ice vertices
        IF (Hi_new( vi_new) < 1._dp) THEN
          ! Remove very thin remapped ice

          Hi_new( vi_new) = 0._dp
        END IF
      END IF

    END DO ! DO vi_new = mesh_new%vi1, mesh_new%vi2

    ! == For those vertices of the new mesh that overlap with both old-mesh shelf and old-mesh
    !    open ocean, average only over the contributing old-mesh shelf vertices
    ! ======================================================================================

    DO vi_new = mesh_new%vi1, mesh_new%vi2

      k1 = M_CSR%ptr( vi_new)
      k2 = M_CSR%ptr( vi_new+1) - 1

      n_shelf      = 0
      n_open_ocean = 0
      sum_Hi_shelf = 0._dp

      DO k = k1, k2

        vi_old = M_CSR%ind( k)

        IF     (mask_floating_ice_tot(  vi_old)) THEN
          n_shelf      = n_shelf      + 1
          sum_Hi_shelf = sum_Hi_shelf + Hi_old_tot( vi_old)
        ELSEIF (mask_icefree_ocean_tot( vi_old)) THEN
          n_open_ocean = n_open_ocean + 1
        END IF

      END DO ! DO k = k1, k2

      IF (n_shelf > 0 .AND. n_open_ocean > 0) THEN
        ! This new-mesh vertex overlaps with both old-mesh shelf vertices,
        ! and old-mesh open-ocean vertices
        Hi_new( vi_new) = sum_Hi_shelf / REAL( n_shelf,dp)
      END IF

    END DO ! DO vi_new = mesh_new%vi1, mesh_new%vi2

    ! Recalculate Hs
    CALL reallocate_bounds_dp_1D( ice%Hs, mesh_new%vi1, mesh_new%vi2)
    DO vi = mesh_new%vi1, mesh_new%vi2
      ice%Hs( vi) = ice_surface_elevation( Hi_new( vi), ice%Hb( vi), ice%SL( vi))
    END DO

    ! Move Hi_new to ice%Hi
    CALL reallocate_bounds_dp_1D( ice%Hi, mesh_new%vi1, mesh_new%vi2)
    ice%Hi = Hi_new

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_basic_ice_geometry

  SUBROUTINE relax_calving_front( mesh_old, mesh, ice, SMB, BMB, LMB, AMB, region_name)
    ! Relax ice thickness around the calving front
    !
    ! This routine "steps out of time" for a bit (default dt_relax = 2 yr), where it
    ! lets the ice thickness near the calving front relax for a little bit. To achieve
    ! this, it uses the BC_prescr_mask optional arguments of the velocity solver and the
    ! thickness solver to keep velocities and thickness fixed over parts of the domain.
    ! In this case, it keeps them fixed everywhere except over the open ocean, and the
    ! first 5 shelf vertices inward from the calving front. This allows the little bumps
    ! that appear in the ice thickness after remapping, which seriously slow down the
    ! velocity solution, to relax a little.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                   INTENT(IN)    :: BMB
    TYPE(type_LMB_model),                   INTENT(IN)    :: LMB
    TYPE(type_AMB_model),                   INTENT(IN)    :: AMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'relax_calving_front'
    LOGICAL,  DIMENSION(mesh%nV)                          :: mask_icefree_ocean_tot
    LOGICAL,  DIMENSION(mesh%nV)                          :: mask_floating_ice_tot
    LOGICAL,  DIMENSION(mesh%nV)                          :: mask_cf_fl_tot
    INTEGER,  DIMENSION(mesh%nV)                          :: BC_prescr_mask_tot
    INTEGER,  DIMENSION(mesh%nTri)                        :: BC_prescr_mask_b_tot
    INTEGER,  DIMENSION(mesh%nTri,mesh%nz)                :: BC_prescr_mask_bk_tot
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2)                :: BC_prescr_mask
    INTEGER,  DIMENSION(mesh%ti1:mesh%ti2)                :: BC_prescr_mask_b
    INTEGER,  DIMENSION(mesh%ti1:mesh%ti2,mesh%nz)        :: BC_prescr_mask_bk
    INTEGER,  DIMENSION(mesh%nV)                          :: map
    INTEGER,  DIMENSION(mesh%nV)                          :: stack
    INTEGER                                               :: stackN
    INTEGER,  DIMENSION(mesh%nV)                          :: stack2
    INTEGER                                               :: stack2N
    INTEGER                                               :: vi
    INTEGER                                               :: it
    INTEGER,  PARAMETER                                   :: nV_around_calving_front = 5
    INTEGER                                               :: i
    INTEGER                                               :: ci,vj
    INTEGER                                               :: ti,via, vib, vic
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: BC_prescr_Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: BC_prescr_u_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)                :: BC_prescr_v_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,mesh%nz)        :: BC_prescr_u_bk
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,mesh%nz)        :: BC_prescr_v_bk
    REAL(dp), PARAMETER                                   :: dt_relax = 2._dp   ! [yr] Time to relax the ice around the calving front
    REAL(dp)                                              :: t_pseudo
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: SMB_new
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: BMB_new
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: LMB_new
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: AMB_new
    REAL(dp)                                              :: visc_it_norm_dUV_tol_save
    INTEGER                                               :: visc_it_nit_save
    REAL(dp)                                              :: visc_it_relax_save
    REAL(dp)                                              :: visc_eff_min_save
    REAL(dp)                                              :: vel_max_save
    REAL(dp)                                              :: stress_balance_PETSc_rtol_save
    REAL(dp)                                              :: stress_balance_PETSc_abstol_save
    REAL(dp)                                              :: Glens_flow_law_epsilon_sq_0_save
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: Hi_tplusdt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: divQ

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather global masks
    CALL gather_to_all_logical_1D( ice%mask_icefree_ocean, mask_icefree_ocean_tot)
    CALL gather_to_all_logical_1D( ice%mask_floating_ice , mask_floating_ice_tot )
    CALL gather_to_all_logical_1D( ice%mask_cf_fl        , mask_cf_fl_tot        )

    ! == Create the relaxation mask
    ! =============================

    ! Create a mask over the entire ice-free ocean, and 5 vertices into the shelf
    ! NOTE: 0 = velocities and dH/dt are solved
    !       1 = velocities and H are prescribed at remapped values)

    ! Let the master do this (difficult to parallelise)
    IF (par%master) THEN

      ! Initialise mask
      BC_prescr_mask_tot   = 1

      ! First mark the open ocean on the mask
      DO vi = 1, mesh%nV
        IF (mask_icefree_ocean_tot( vi)) THEN
          BC_prescr_mask_tot( vi) = 0
        END IF
      END DO ! DO vi = 1, mesh%nV

      ! Initialise flood-fill
      map    = 0
      stackN = 0

      ! Initialise with the floating calving front
      DO vi = 1, mesh%nV
        IF (mask_cf_fl_tot( vi)) THEN
          map( vi) = 1
          stackN = stackN + 1
          stack( stackN) = vi
        END IF
      END DO ! DO vi = 1, mesh%nV

      ! Expand into shelf
      DO it = 1, nV_around_calving_front

        ! Initialise the second stack
        stack2N = 0

        ! Go over the entire stack
        DO i = 1, stackN

          ! Take a vertex from the stack
          vi = stack( i)

          ! Mark it on the mask
          BC_prescr_mask_tot( vi) = 0

          ! Mark it on the map
          map( vi) = 2

          ! Add all its un-mapped, floating neighbours to the new stack
          DO ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            IF (map( vj) == 0 .AND. mask_floating_ice_tot( vj)) THEN
              map( vj) = 1
              stack2N = stack2N + 1
              stack2( stack2N) = vj
            END IF
          END DO !  DO ci = 1, mesh%nC( vi)

        END DO ! DO i = 1, stackN

        ! Cycle the stacks
        stack( 1:stack2N) = stack2( 1:stack2N)
        stackN = stack2N

      END DO ! DO it = 1, nV_around_calving_front

      ! Fill in the b-grid masks

      BC_prescr_mask_b_tot  = 0
      BC_prescr_mask_bk_tot = 0

      DO ti = 1, mesh%nTri

        ! The three vertices spanning triangle ti
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        ! Only prescribe velocities on triangles where thickness is prescribed on all three corners
        IF (BC_prescr_mask_tot( via) == 1 .AND. &
            BC_prescr_mask_tot( vib) == 1 .AND. &
            BC_prescr_mask_tot( vic) == 1) THEN
          BC_prescr_mask_b_tot(  ti  ) = 1
          BC_prescr_mask_bk_tot( ti,:) = 1
        END IF

      END DO ! DO ti = 1, mesh%nTri

    END IF ! IF (par%master) THEN

    ! Distribute BC masks to all processes
    CALL distribute_from_master_int_1D( BC_prescr_mask_tot   , BC_prescr_mask   )
    CALL distribute_from_master_int_1D( BC_prescr_mask_b_tot , BC_prescr_mask_b )
    CALL distribute_from_master_int_2D( BC_prescr_mask_bk_tot, BC_prescr_mask_bk)

    ! == Fill in prescribed velocities and thicknesses away from the front
    ! ====================================================================

    BC_prescr_Hi   = ice%Hi
    BC_prescr_u_b  = ice%u_vav_b
    BC_prescr_v_b  = ice%v_vav_b
    BC_prescr_u_bk = ice%u_3D_b
    BC_prescr_v_bk = ice%v_3D_b

    ! == Save proper values of config parameters for the velocity solver
    ! ==================================================================

    visc_it_norm_dUV_tol_save                  = C%visc_it_norm_dUV_tol
    visc_it_nit_save                           = C%visc_it_nit
    visc_it_relax_save                         = C%visc_it_relax
    visc_eff_min_save                          = C%visc_eff_min
    vel_max_save                               = C%vel_max
    stress_balance_PETSc_rtol_save             = C%stress_balance_PETSc_rtol
    stress_balance_PETSc_abstol_save           = C%stress_balance_PETSc_abstol
    Glens_flow_law_epsilon_sq_0_save           = C%Glens_flow_law_epsilon_sq_0

    ! == Set temporary, less strict values of config parameters for the velocity solver
    ! =================================================================================

    C%visc_it_norm_dUV_tol                  = 5E-4_dp                          ! Stop criterion for the viscosity iteration: the L2-norm of successive velocity solutions should be smaller than this number
    C%visc_it_nit                           = 20                               ! Maximum number of effective viscosity iterations
    C%visc_it_relax                         = 0.3_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    C%visc_eff_min                          = 1E5_dp                           ! Minimum value for effective viscosity
    C%vel_max                               = 5000._dp                         ! Velocities are limited to this value
    C%stress_balance_PETSc_rtol             = 1E-3_dp                          ! PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    C%stress_balance_PETSc_abstol           = 1E-2_dp                          ! PETSc solver - stop criterion, absolute difference
    C%Glens_flow_law_epsilon_sq_0           = 1E-6_dp                          ! Normalisation term so that zero strain rates produce a high but finite viscosity

    ! == Remap SMB, BMB, LMB, and AMB to get more stable ice thickness
    ! ================================================================

    CALL map_from_mesh_to_mesh_2D( mesh_old, mesh, SMB%SMB, SMB_new, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_2D( mesh_old, mesh, BMB%BMB, BMB_new, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_2D( mesh_old, mesh, LMB%LMB, LMB_new, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_2D( mesh_old, mesh, AMB%AMB, AMB_new, '2nd_order_conservative')

    ! == Relax the ice thickness for a few time steps
    ! ===============================================

    t_pseudo = 0._dp

    pseudo_time: DO WHILE (t_pseudo < dt_relax)

      ! Update velocity solution around the calving front
      CALL solve_stress_balance( mesh, ice, BMB_new, region_name, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

      ! Calculate dH/dt around the calving front
      CALL calc_dHi_dt( mesh, ice%Hi, ice%Hb, ice%SL, ice%u_vav_b, ice%v_vav_b, SMB_new, BMB_new, LMB_new, AMB_new, ice%fraction_margin, ice%mask_noice, C%dt_ice_min, &
        ice%dHi_dt, Hi_tplusdt, divQ, ice%dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)

      ! Update ice thickness and advance pseudo-time
      ice%Hi = Hi_tplusdt
      t_pseudo = t_pseudo + C%dt_ice_min

      ! Update basic geometry
      DO vi = mesh%vi1, mesh%vi2
        ice%Hs ( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))
        ice%Hib( vi) = ice%Hs( vi) - ice%Hi( vi)
        ice%TAF( vi) = thickness_above_floatation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))
      END DO

    END DO pseudo_time ! DO WHILE (t_pseudo < dt_relax)

    ! == Reinstate proper values of config parameters for the velocity solver
    ! =======================================================================

    C%visc_it_norm_dUV_tol                  = visc_it_norm_dUV_tol_save
    C%visc_it_nit                           = visc_it_nit_save
    C%visc_it_relax                         = visc_it_relax_save
    C%visc_eff_min                          = visc_eff_min_save
    C%vel_max                               = vel_max_save
    C%stress_balance_PETSc_rtol             = stress_balance_PETSc_rtol_save
    C%stress_balance_PETSc_abstol           = stress_balance_PETSc_abstol_save
    C%Glens_flow_law_epsilon_sq_0           = Glens_flow_law_epsilon_sq_0_save

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_calving_front

  SUBROUTINE apply_geometry_relaxation( region)
    ! Relax the initial geometry by running the ice dynamics model
    ! for a certain time without any mass balance terms, inversions
    ! or imposed alterations to the evolution of ice thickness

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_carriage_return

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_geometry_relaxation'
    INTEGER                                               :: vi
    REAL(dp)                                              :: t_pseudo, t_step
    REAL(dp)                                              :: visc_it_norm_dUV_tol_config
    REAL(dp), DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: Hi_new, dHi_dt_new
    REAL(dp), DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, dHi_dt_target_dummy
    CHARACTER(LEN=256)                                    :: t_years, r_time, r_step, r_adv, t_format

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master .AND. C%do_time_display .AND. C%geometry_relaxation_t_years > 0._dp) THEN

      IF (C%geometry_relaxation_t_years <= 999._dp .AND. &
          C%geometry_relaxation_t_years >=   1._dp) THEN
        write(*,"(A,F6.2,A)") '   Stepping out of time to relax geometry for ', C%geometry_relaxation_t_years, ' pseudo years...'
      ELSE
        write(*,"(A)") '   Stepping out of time to relax geometry...'
      END IF

      ! Initialise and display pseudo time
      t_pseudo = 0._dp
      write( r_time,"(F7.3)") t_pseudo
      write( r_step,"(F5.3)") C%geometry_relaxation_t_years / 100._dp
      write( *,"(A)", ADVANCE = TRIM( 'no')) c_carriage_return // &
                                             "     t_pseudo = " // TRIM( r_time) // &
                                             " yr - dt = " // TRIM( r_step) // " yr"
    END IF

    DO WHILE (t_pseudo < C%geometry_relaxation_t_years)

      ! Save user-defined viscosity-iteration tolerance to recover it later
      visc_it_norm_dUV_tol_config = C%visc_it_norm_dUV_tol

      ! Set viscosity-iteration to a high enough value, which will allow the
      ! stress balance solver to move forward when given rough initial conditions
      ! NOTE: this is relevant mostly during the first time steps only, since
      ! after a short while the rough surface gradients relax enough to allow
      ! for smaller and smaller residuals during the viscosity iterations
      C%visc_it_norm_dUV_tol = 1._dp

      ! Default time step: set so the relaxation takes 100 iterations
      t_step = C%geometry_relaxation_t_years / 100._dp
      ! But prevent time step larger than maximum allowed
      t_step = MIN( t_step, C%dt_ice_max)

      ! Ignore any mass balance terms
      SMB_dummy = 0._dp
      BMB_dummy = 0._dp
      LMB_dummy = 0._dp
      AMB_dummy = 0._dp

      ! Ignore any target thinning rates
      dHi_dt_target_dummy = 0._dp

      region%ice%effective_pressure = MAX( 0._dp, ice_density * grav * region%ice%Hi_eff) * region%ice%fraction_gr

      ! Calculate ice velocities for the predicted geometry
      CALL solve_stress_balance( region%mesh, region%ice, BMB_dummy, region%name)

      ! Calculate thinning rates for current geometry and velocity
      CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, region%ice%fraction_margin, &
                        region%ice%mask_noice, t_step, dHi_dt_new, Hi_new, region%ice%divQ, dHi_dt_target_dummy)

      ! Set ice model ice thickness to relaxed field
      DO vi = region%mesh%vi1, region%mesh%vi2
        ! Apply relaxation over ice shelves and grounding lines
        IF (region%ice%mask_floating_ice( vi) .OR. region%ice%mask_gl_gr( vi)) THEN
          region%ice%Hi( vi) = Hi_new( vi)
          region%ice%dHi_dt( vi) = dHi_dt_new( vi)
        ! Also over steep-sloped interior ice sheet points
        ELSEIF (region%ice%mask_grounded_ice( vi) .AND. region%ice%Hs_slope( vi) >= 0.03_dp) THEN
          region%ice%Hi( vi) = Hi_new( vi)
          region%ice%dHi_dt( vi) = dHi_dt_new( vi)
        END IF
      END DO

      ! Apply some specific corrections
      DO vi = region%mesh%vi1, region%mesh%vi2
        ! Don't let grounded ice cross the floatation threshold
        ! IF (region%ice%mask_grounded_ice( vi)) THEN
        !   region%ice%Hi( vi) = MAX( region%ice%Hi( vi), (region%ice%SL( vi) - region%ice%Hb( vi)) * seawater_density/ice_density + .1_dp)
        ! END IF
        ! Remove very thin ice
        IF (region%ice%Hi( vi) < C%Hi_min) THEN
          region%ice%Hi( vi) = 0._dp
          region%ice%dHi_dt( vi) = 0._dp
        END IF
        ! Remove ice absent at PD
        IF (region%refgeo_PD%Hi( vi) == 0._dp) THEN
          region%ice%Hi( vi) = 0._dp
          region%ice%dHi_dt( vi) = 0._dp
        END IF
      END DO

      ! Calculate all other ice geometry quantities
      ! ===========================================

      DO vi = region%mesh%vi1, region%mesh%vi2

        ! Basic geometry
        region%ice%Hs ( vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
        region%ice%Hib( vi) = region%ice%Hs( vi) - region%ice%Hi( vi)
        region%ice%TAF( vi) = thickness_above_floatation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))

        IF (region%ice%TAF( vi) > 0._dp) THEN
          ! Grounded ice
          region%ice%dHs_dt ( vi) = region%ice%dHb_dt( vi) + region%ice%dHi_dt( vi)
          region%ice%dHib_dt( vi) = region%ice%dHb_dt( vi)
        ELSE
          ! Floating ice
          region%ice%dHs_dt ( vi) = region%ice%dHi_dt( vi) * (1._dp - ice_density / seawater_density)
          region%ice%dHib_dt( vi) = region%ice%dHi_dt( vi) *          ice_density / seawater_density
        END IF

      END DO

      ! Update masks
      CALL determine_masks( region%mesh, region%ice)

      ! Calculate new effective thickness
      CALL calc_effective_thickness( region%mesh, region%ice, region%ice%Hi, region%ice%Hi_eff, region%ice%fraction_margin)

      ! NOTE: as calculating the zeta gradients is quite expensive, only do so when necessary,
      !       i.e. when solving the heat equation or the Blatter-Pattyn stress balance
      ! Calculate zeta gradients
      CALL calc_zeta_gradients( region%mesh, region%ice)

      ! Calculate sub-grid grounded-area fractions
      CALL calc_grounded_fractions( region%mesh, region%ice)

      ! Reference geometry
      ! ==================

      region%refgeo_PD%Hi  = region%ice%Hi
      region%refgeo_PD%Hs  = region%ice%Hs
      region%refgeo_PD%Hb  = region%ice%Hb

      ! Differences w.r.t. present-day
      region%ice%dHi  = 0._dp
      region%ice%dHb  = 0._dp
      region%ice%dHs  = 0._dp
      region%ice%dHib = 0._dp

      ! Re-initialise previous and next Hi states
      region%ice%Hi_prev = region%ice%Hi
      region%ice%Hi_next = region%ice%Hi

      ! Advance pesudo time
      ! ===================

      t_pseudo = t_pseudo + t_step

      ! Time display
      IF (par%master .AND. C%do_time_display) THEN
        ! Carriage return flag
        r_adv = "no"
        IF (t_pseudo >= C%geometry_relaxation_t_years) r_adv = "yes"
        ! Current pseudo time
        write( r_time,"(F7.3)") MIN( t_pseudo,C%geometry_relaxation_t_years)
        ! Current time step
        write( r_step,"(F5.3)") t_step
        ! Time display message
        write( *,"(A)", ADVANCE = TRIM( r_adv)) c_carriage_return // &
                                                "     t_pseudo = " // TRIM( r_time) // &
                                                " yr - dt = " // TRIM( r_step) // " yr"
      END IF

      ! Retrieve user-defined viscosity-iteration tolerance
      C%visc_it_norm_dUV_tol = visc_it_norm_dUV_tol_config

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_geometry_relaxation

! ===== Predictor-corrector scheme =====
! ======================================

  SUBROUTINE run_ice_dynamics_model_pc( region, dt_max)
    ! Calculate a new next modelled ice thickness

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region
    REAL(dp),                               INTENT(IN)    :: dt_max

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ice_dynamics_model_pc'
    REAL(dp)                                              :: dt_crit_adv
    INTEGER                                               :: pc_it
    REAL(dp), DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: Hi_dummy, dHi_dt_dummy, LMB_dummy, AMB_dummy
    INTEGER                                               :: vi, n_guilty, n_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Store previous ice model state
    region%ice%t_Hi_prev  = region%ice%t_Hi_next
    region%ice%Hi_prev    = region%ice%Hi_next

    ! == Calculate time step ==
    ! =========================

    ! Store previous time step
    region%ice%pc%dt_n = region%ice%pc%dt_np1

    ! Calculate new time step (Robinson et al., 2020, Eq. 33)
    region%ice%pc%dt_np1 = (C%pc_epsilon / region%ice%pc%eta_np1)**(C%pc_k_I + C%pc_k_p) * &
      (C%pc_epsilon / region%ice%pc%eta_n)**(-C%pc_k_p) * region%ice%pc%dt_n

    ! Limit time step to maximum allowed value
    region%ice%pc%dt_np1 = MIN( region%ice%pc%dt_np1, dt_max)

    ! Limit time step to 1.2 times the previous time step
    region%ice%pc%dt_np1 = MIN( region%ice%pc%dt_np1, C%pc_max_time_step_increase * region%ice%pc%dt_n)

    ! Limit time step to minimum allowed value
    region%ice%pc%dt_np1 = MAX( region%ice%pc%dt_np1, C%dt_ice_min)

    ! Limit time step to critical advective time step
    CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_adv)

    region%ice%pc%dt_np1 = MIN( region%ice%pc%dt_np1, dt_crit_adv)

    ! == Time step iteration: if, at the end of the PC timestep, the truncation error
    !    turns out to be too large, run it again with a smaller dt, until the truncation
    !    decreases to below the specified tolerance
    ! ==================================================================================

    ! Store thinning rates from previous time step
    region%ice%pc%dHi_dt_Hi_nm1_u_nm1 = region%ice%dHi_dt

    ! Store the previous maximum truncation error eta_n
    region%ice%pc%eta_n = region%ice%pc%eta_np1

    pc_it = 0
    iterate_pc_timestep: DO WHILE (pc_it < C%pc_nit_max)

      pc_it = pc_it + 1

      ! Calculate time step ratio
      region%ice%pc%zeta_t = region%ice%pc%dt_np1 / region%ice%pc%dt_n

      ! == Predictor step ==
      ! ====================

      ! Invert a basal mass balance field that keeps the ice shelves in equilibrium
      CALL BMB_inversion( region, region%ice%pc%dt_np1)

      ! Invert a lateral mass balance field that keeps the calving fronts in equilibrium
      CALL LMB_inversion( region, region%ice%pc%dt_np1)

      ! Invert a surface mass balance field that keeps the ice sheet in check
      CALL SMB_inversion( region, region%ice%pc%dt_np1)

      ! Calculate thinning rates for current geometry and velocity
      CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, region%AMB%AMB, region%ice%fraction_margin, &
                        region%ice%mask_noice, region%ice%pc%dt_np1, region%ice%pc%dHi_dt_Hi_n_u_n, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

      ! ! If so desired, alter the computed dH/dt by adjusting dummy mass balance fluxes to get an equilibrium state
      ! CALL MB_inversion( region%mesh, region%ice, region%refgeo_PD, region%SMB, region%BMB, region%LMB, region%AMB, region%ice%pc%dHi_dt_Hi_n_u_n, Hi_dummy, region%ice%pc%dt_np1, region%time, region%name)

      ! Calculate predicted ice thickness (Robinson et al., 2020, Eq. 30)
      region%ice%pc%Hi_star_np1 = region%ice%Hi_prev + region%ice%pc%dt_np1 * ((1._dp + region%ice%pc%zeta_t / 2._dp) * &
        region%ice%pc%dHi_dt_Hi_n_u_n - (region%ice%pc%zeta_t / 2._dp) * region%ice%pc%dHi_dt_Hi_nm1_u_nm1)

      ! If so desired, modify the predicted ice thickness field based on user-defined settings
      CALL alter_ice_thickness( region%mesh, region%ice, region%ice%Hi_prev, region%ice%pc%Hi_star_np1, region%refgeo_PD, region%time)

      ! Adjust the predicted dHi_dt to compensate for thickness modifications
      ! This is just Robinson et al., 2020, Eq 30 above rearranged to retrieve
      ! an updated dHi_dt_Hi_n_u_n from the modified Hi_star_np1. If no ice
      ! thickness modifications were applied, then there will be not change.
      region%ice%pc%dHi_dt_Hi_n_u_n = ((region%ice%pc%Hi_star_np1 - region%ice%Hi_prev) / region%ice%pc%dt_np1 + (region%ice%pc%zeta_t / 2._dp) * region%ice%pc%dHi_dt_Hi_nm1_u_nm1) / (1._dp + region%ice%pc%zeta_t / 2._dp)

      ! == Update step ==
      ! =================

      ! Set model geometry to predicted
      region%ice%Hi = region%ice%pc%Hi_star_np1

      ! Set thinning rates to predicted
      region%ice%dHi_dt = (region%ice%Hi - region%ice%Hi_prev) / region%ice%pc%dt_np1

      ! Set model geometry to predicted
      DO vi = region%mesh%vi1, region%mesh%vi2
        ! Basic geometry
        region%ice%Hs ( vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
        region%ice%Hib( vi) = region%ice%Hs(  vi) - region%ice%Hi( vi)
      END DO

      ! Update masks
      CALL determine_masks( region%mesh, region%ice)

      ! DENK DROM : assess whether this is important for the velocitiy computation below
      ! ! Calculate zeta gradients
      ! CALL calc_zeta_gradients( region%mesh, region%ice)

      ! Update sub-grid grounded fractions
      CALL calc_grounded_fractions( region%mesh, region%ice)

      ! DENK DROM : assess whether this is important for the velocitiy computation below
      ! ! Calculate the basal mass balance
      ! CALL run_BMB_model( region%mesh, region%ice, region%ocean, region%refgeo_PD, region%SMB, region%BMB, region%name, region%time)

      ! Calculate ice velocities for the predicted geometry
      CALL solve_stress_balance( region%mesh, region%ice, region%BMB%BMB, region%name)

      ! == Corrector step ==
      ! ====================

      ! Set model geometry back to original
      DO vi = region%mesh%vi1, region%mesh%vi2
        region%ice%Hi(  vi) = region%ice%Hi_prev( vi)
        region%ice%Hs(  vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
        region%ice%Hib( vi) = region%ice%Hs(  vi) - region%ice%Hi( vi)
      END DO

      ! Update masks
      CALL determine_masks( region%mesh, region%ice)

      ! Update sub-grid grounded fractions
      CALL calc_grounded_fractions( region%mesh, region%ice)

      ! Update effective ice thickness
      CALL calc_effective_thickness( region%mesh, region%ice, region%ice%Hi, region%ice%Hi_eff, region%ice%fraction_margin)

      ! Invert a basal mass balance field that keeps the ice shelves in equilibrium
      CALL BMB_inversion( region, region%ice%pc%dt_np1)

      ! Invert a lateral mass balance field that keeps the calving fronts in equilibrium
      CALL LMB_inversion( region, region%ice%pc%dt_np1)

      ! Invert a surface mass balance field that keeps the ice sheet in check
      CALL SMB_inversion( region, region%ice%pc%dt_np1)

      ! Calculate thinning rates for the current ice thickness and predicted velocity
      CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, region%AMB%AMB, region%ice%fraction_margin, &
                        region%ice%mask_noice, region%ice%pc%dt_np1, region%ice%pc%dHi_dt_Hi_star_np1_u_np1, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

      ! ! If so desired, alter the computed dH/dt by adjusting dummy mass balance fluxes to get an equilibrium state
      ! CALL MB_inversion( region%mesh, region%ice, region%refgeo_PD, region%SMB, region%BMB, region%LMB, region%AMB, region%ice%pc%dHi_dt_Hi_star_np1_u_np1, Hi_dummy, region%ice%pc%dt_np1, region%time, region%name)

      ! Calculate corrected ice thickness (Robinson et al. (2020), Eq. 31)
      region%ice%pc%Hi_np1 = region%ice%Hi_prev + (region%ice%pc%dt_np1 / 2._dp) * (region%ice%pc%dHi_dt_Hi_n_u_n + region%ice%pc%dHi_dt_Hi_star_np1_u_np1)

      ! Save "raw" thinning rates, as applied after the corrector step
      region%ice%dHi_dt_raw = (region%ice%pc%Hi_np1 - region%ice%Hi_prev) / region%ice%pc%dt_np1

      ! If so desired, modify the corrected ice thickness field based on user-defined settings
      CALL alter_ice_thickness( region%mesh, region%ice, region%ice%Hi_prev, region%ice%pc%Hi_np1, region%refgeo_PD, region%time)

      ! Adjust the predicted dHi_dt to compensate for thickness modifications
      ! This is just Robinson et al., 2020, Eq 31 above rearranged to retrieve
      ! an updated dHi_dt_Hi_star_np1_u_np1 from the modified Hi_np1. If no ice
      ! thickness modifications were applied, then there will be not change.
      region%ice%pc%dHi_dt_Hi_star_np1_u_np1 = (region%ice%pc%Hi_np1 - region%ice%Hi_prev) / (region%ice%pc%dt_np1 / 2._dp) - region%ice%pc%dHi_dt_Hi_n_u_n

      ! Add difference between raw and applied dHi_dt to residual tracker
      region%ice%dHi_dt_residual = region%ice%dHi_dt_raw - &  ! Raw change
                                   (region%ice%pc%Hi_np1 - region%ice%Hi_prev) / region%ice%pc%dt_np1  ! Minus applied change

      ! == Truncation error ==
      ! ======================

      ! Estimate truncation error
      CALL calc_pc_truncation_error( region%mesh, region%ice, region%ice%pc)

      ! == Error assessment ==
      ! ======================

      ! Initialise unstable vertex count
      n_tot = 0
      n_guilty = 0

      ! Determine number of unstable vertices
      DO vi = region%mesh%vi1, region%mesh%vi2
        ! Only consider fully grounded vertices
        IF (region%ice%fraction_gr( vi) < 1._dp) CYCLE
        ! If so, add to total vertex count
        n_tot = n_tot + 1
        ! If this vertex's error is larger than tolerance
        IF (region%ice%pc%tau_np1( vi) > C%pc_epsilon) THEN
          ! Add to total guilty vertex count
          n_guilty = n_guilty + 1
          ! Add to this vertex's guilty record
          region%ice%pc%tau_n_guilty( vi) = region%ice%pc%tau_n_guilty( vi) + 1
        END IF
      END DO

      ! Add up findings from each process domain
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_tot,    1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_guilty, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Safety
      IF (n_tot == 0) n_tot = 1

      ! Check if largest truncation error is small enough; if so, move on
      IF (region%ice%pc%eta_np1 < C%pc_epsilon) THEN
        EXIT iterate_pc_timestep

      ! If not, check whether that occurs in a significant amount of vertices; if not,
      ! set the truncation error to almost the tolerance (to allow for growth) and move on
      ELSEIF (100._dp * REAL( n_guilty,dp) / REAL(n_tot,dp) < C%pc_guilty_max) THEN
        ! IF (par%master) CALL warning('{dp_01}% of vertices are changing rapidly, ignoring for now', dp_01 = 100._dp * REAL( n_guilty,dp) / REAL(n_tot,dp))
        region%ice%pc%eta_np1 = .95_dp * C%pc_epsilon
        EXIT iterate_pc_timestep

      ! If not, re-do the PC timestep
      ELSE
        IF (par%master) CALL warning('{dp_01}% of vertices ({int_01}) are changing rapidly (eta = {dp_02}), reducing dt and redoing PC timestep', dp_01 = 100._dp * REAL( n_guilty,dp) / REAL(n_tot,dp), int_01 = n_guilty, dp_02 = region%ice%pc%eta_np1)
        region%ice%pc%dt_np1 = region%ice%pc%dt_np1 * 0.8_dp
        ! If the timestep has reached the specified lower limit, stop iterating
        IF (region%ice%pc%dt_np1 <= C%dt_ice_min) THEN
          region%ice%pc%dt_np1 = C%dt_ice_min
          EXIT iterate_pc_timestep
        END IF
      END IF

    END DO iterate_pc_timestep

    ! == Final quantities
    ! ===================

    ! Set next modelled ice thickness
    region%ice%t_Hi_next = region%ice%t_Hi_prev + region%ice%pc%dt_np1
    region%ice%Hi_next   = region%ice%pc%Hi_np1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_model_pc

  SUBROUTINE calc_pc_truncation_error( mesh, ice, pc)
    ! Calculate the truncation error tau in the ice thickness rate of change (Robinson et al., 2020, Eq. 32)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ice_pc),                   INTENT(INOUT) :: pc

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pc_truncation_error'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate truncation error tau (Robinson et al., 2020, Eq. 32)
    DO vi = mesh%vi1, mesh%vi2
      pc%tau_np1( vi) = pc%zeta_t * ABS( pc%Hi_np1( vi) - pc%Hi_star_np1( vi)) / ((3._dp * pc%zeta_t + 3._dp) * pc%dt_n)
    END DO

    ! Calculate the maximum truncation error eta over grounded ice only
    pc%eta_np1 = C%pc_eta_min
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_grounded_ice( vi) .AND. .NOT. ice%mask_gl_gr( vi) .AND. ice%fraction_gr( vi) == 1._dp) THEN
        pc%eta_np1 = MAX( pc%eta_np1, pc%tau_np1( vi))
      END IF
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, pc%eta_np1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pc_truncation_error

  SUBROUTINE initialise_pc_scheme( mesh, pc, region_name)
    ! Allocate memory and initialise values for the ice thickness predictor/corrector scheme.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    TYPE(type_ice_pc),             INTENT(OUT)   :: pc
    CHARACTER(LEN=3),              INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_pc_scheme'
    CHARACTER(LEN=256)                           :: pc_choice_initialise
    CHARACTER(LEN=256)                           :: filename_pc_initialise
    REAL(dp)                                     :: timeframe_pc_initialise

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ! ===============

    ALLOCATE( pc%dHi_dt_Hi_nm1_u_nm1(      mesh%vi1:mesh%vi2))           ! [m/yr] Thinning rates from previous time step
    ALLOCATE( pc%dHi_dt_Hi_n_u_n(          mesh%vi1:mesh%vi2))           ! [m/yr] Thinning rates for current time step with old geometry
    ALLOCATE( pc%Hi_star_np1(              mesh%vi1:mesh%vi2))           ! [m]    Predicted ice thickness
    ALLOCATE( pc%dHi_dt_Hi_star_np1_u_np1( mesh%vi1:mesh%vi2))           ! [m/yr] Thinning rates for predicted ice thickness and updated velocity
    ALLOCATE( pc%Hi_np1(                   mesh%vi1:mesh%vi2))           ! [m]    Corrected ice thickness
    ALLOCATE( pc%tau_np1(                  mesh%vi1:mesh%vi2))           ! [m]    Truncation error
    ALLOCATE( pc%tau_n_guilty(             mesh%vi1:mesh%vi2))           ! [-]    Number of PC iterations where vertex had truncation errors above the tolerance

    ! Initialise
    ! ==========

    IF     (region_name == 'NAM') THEN
      pc_choice_initialise    = C%pc_choice_initialise_NAM
      filename_pc_initialise  = C%filename_pc_initialise_NAM
      timeframe_pc_initialise = C%timeframe_pc_initialise_NAM
    ELSEIF (region_name == 'EAS') THEN
      pc_choice_initialise    = C%pc_choice_initialise_EAS
      filename_pc_initialise  = C%filename_pc_initialise_EAS
      timeframe_pc_initialise = C%timeframe_pc_initialise_EAS
    ELSEIF (region_name == 'GRL') THEN
      pc_choice_initialise    = C%pc_choice_initialise_GRL
      filename_pc_initialise  = C%filename_pc_initialise_GRL
      timeframe_pc_initialise = C%timeframe_pc_initialise_GRL
    ELSEIF (region_name == 'ANT') THEN
      pc_choice_initialise    = C%pc_choice_initialise_ANT
      filename_pc_initialise  = C%filename_pc_initialise_ANT
      timeframe_pc_initialise = C%timeframe_pc_initialise_ANT
    ELSE
      CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END IF

    IF     (pc_choice_initialise == 'zero') THEN
      ! Initialise everything from scratch

      pc%dt_n                     = C%dt_ice_min
      pc%dt_np1                   = C%dt_ice_min
      pc%zeta_t                   = 1._dp
      pc%dHi_dt_Hi_nm1_u_nm1      = 0._dp
      pc%dHi_dt_Hi_n_u_n          = 0._dp
      pc%Hi_star_np1              = 0._dp
      pc%dHi_dt_Hi_star_np1_u_np1 = 0._dp
      pc%Hi_np1                   = 0._dp
      pc%tau_np1                  = C%pc_epsilon
      pc%eta_n                    = C%pc_epsilon
      pc%eta_np1                  = C%pc_epsilon

    ELSEIF (pc_choice_initialise == 'read_from_file') THEN
      ! Initialise from a (restart) file
      CALL initialise_pc_scheme_from_file( pc, filename_pc_initialise, timeframe_pc_initialise)
    ELSE
      CALL crash('unknown pc_choice_initialise "' // TRIM( pc_choice_initialise) // '"!')
    END IF

    ! Initialise the event counter for errors above tolerance
    pc%tau_n_guilty = 0

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_pc_scheme

  SUBROUTINE initialise_pc_scheme_from_file( pc, filename, timeframe)
    ! Initialise values for the ice thickness predictor/corrector scheme from a (restart) file.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_ice_pc),             INTENT(INOUT) :: pc
    CHARACTER(LEN=256),            INTENT(IN)    :: filename
    REAL(dp),                      INTENT(IN)    :: timeframe

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_pc_scheme_from_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write to terminal
    IF (par%master) WRITE(0,*) '   Initialising ice thickness predictor/corrector scheme from file "' // colour_string( TRIM( filename),'light blue') // '"...'

    ! Read values from the file
    IF (timeframe == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_0D(      filename, 'dt_n'                    , pc%dt_n                    )
      CALL read_field_from_file_0D(      filename, 'dt_np1'                  , pc%dt_np1                  )
      CALL read_field_from_file_0D(      filename, 'zeta_t'                  , pc%zeta_t                  )
      CALL read_field_from_mesh_file_2D( filename, 'dHi_dt_Hi_nm1_u_nm1'     , pc%dHi_dt_Hi_nm1_u_nm1     )
      CALL read_field_from_mesh_file_2D( filename, 'dHi_dt_Hi_n_u_n'         , pc%dHi_dt_Hi_n_u_n         )
      CALL read_field_from_mesh_file_2D( filename, 'Hi_star_np1'             , pc%Hi_star_np1             )
      CALL read_field_from_mesh_file_2D( filename, 'dHi_dt_Hi_star_np1_u_np1', pc%dHi_dt_Hi_star_np1_u_np1)
      CALL read_field_from_mesh_file_2D( filename, 'Hi_np1'                  , pc%Hi_np1                  )
      CALL read_field_from_mesh_file_2D( filename, 'tau_np1'                 , pc%tau_np1                 )
      CALL read_field_from_file_0D(      filename, 'eta_n'                   , pc%eta_n                   )
      CALL read_field_from_file_0D(      filename, 'eta_np1'                 , pc%eta_np1                 )
    ELSE
      ! Read specified timeframe
      CALL read_field_from_file_0D(      filename, 'dt_n'                    , pc%dt_n                    , time_to_read = timeframe)
      CALL read_field_from_file_0D(      filename, 'dt_np1'                  , pc%dt_np1                  , time_to_read = timeframe)
      CALL read_field_from_file_0D(      filename, 'zeta_t'                  , pc%zeta_t                  , time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D( filename, 'dHi_dt_Hi_nm1_u_nm1'     , pc%dHi_dt_Hi_nm1_u_nm1     , time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D( filename, 'dHi_dt_Hi_n_u_n'         , pc%dHi_dt_Hi_n_u_n         , time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D( filename, 'Hi_star_np1'             , pc%Hi_star_np1             , time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D( filename, 'dHi_dt_Hi_star_np1_u_np1', pc%dHi_dt_Hi_star_np1_u_np1, time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D( filename, 'Hi_np1'                  , pc%Hi_np1                  , time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D( filename, 'tau_np1'                 , pc%tau_np1                 , time_to_read = timeframe)
      CALL read_field_from_file_0D(      filename, 'eta_n'                   , pc%eta_n                   , time_to_read = timeframe)
      CALL read_field_from_file_0D(      filename, 'eta_np1'                 , pc%eta_np1                 , time_to_read = timeframe)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_pc_scheme_from_file

  SUBROUTINE write_to_restart_file_pc_scheme( mesh, pc, time)
    ! Write to the restart NetCDF file for the ice thickness predictor/corrector scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_pc),                   INTENT(IN)              :: pc
    REAL(dp),                            INTENT(IN)              :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'write_to_restart_file_pc_scheme'
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Writing to ice dynamics restart file "' // &
      colour_string( TRIM( pc%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( pc%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( pc%restart_filename, ncid, time)

    ! Write the data fields to the file
    CALL write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'dt_n'                    , pc%dt_n                    )
    CALL write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'dt_np1'                  , pc%dt_np1                  )
    CALL write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'zeta_t'                  , pc%zeta_t                  )
    CALL write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'dHi_dt_Hi_nm1_u_nm1'     , pc%dHi_dt_Hi_nm1_u_nm1     )
    CALL write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'dHi_dt_Hi_n_u_n'         , pc%dHi_dt_Hi_n_u_n         )
    CALL write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'Hi_star_np1'             , pc%Hi_star_np1             )
    CALL write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'dHi_dt_Hi_star_np1_u_np1', pc%dHi_dt_Hi_star_np1_u_np1)
    CALL write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'Hi_np1'                  , pc%Hi_np1                  )
    CALL write_to_field_multopt_mesh_dp_2D( mesh, pc%restart_filename, ncid, 'tau_np1'                 , pc%tau_np1                 )
    CALL write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'eta_n'                   , pc%eta_n                   )
    CALL write_to_field_multopt_dp_0D(            pc%restart_filename, ncid, 'eta_np1'                 , pc%eta_np1                 )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_pc_scheme

  SUBROUTINE create_restart_file_pc_scheme( mesh, pc)
    ! Create a restart NetCDF file for the ice thickness predictor/corrector scheme
    ! Includes generation of the procedural filename (e.g. "restart_pc_00001.nc")

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_pc),                   INTENT(INOUT)           :: pc

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'create_restart_file_pc_scheme'
    CHARACTER(LEN=256)                                           :: filename_base
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_pc_scheme'
    CALL generate_filename_XXXXXdotnc( filename_base, pc%restart_filename)

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Creating ice dynamics restart file "' // &
      colour_string( TRIM( pc%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( pc%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( pc%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( pc%restart_filename, ncid)

    ! Add the data fields to the file
    CALL add_field_dp_0D(      pc%restart_filename, ncid, 'dt_n'                    , long_name = 'Previous time step', units = 'yr')
    CALL add_field_dp_0D(      pc%restart_filename, ncid, 'dt_np1'                  , long_name = 'Current time step' , units = 'yr')
    CALL add_field_dp_0D(      pc%restart_filename, ncid, 'zeta_t'                  , long_name = 'Ratio between previous and new time step')
    CALL add_field_mesh_dp_2D( pc%restart_filename, ncid, 'dHi_dt_Hi_nm1_u_nm1'     , long_name = 'Thinning rates from previous time step', units = 'm/yr')
    CALL add_field_mesh_dp_2D( pc%restart_filename, ncid, 'dHi_dt_Hi_n_u_n'         , long_name = 'Thinning rates for current time step with old geometry', units = 'm/yr')
    CALL add_field_mesh_dp_2D( pc%restart_filename, ncid, 'Hi_star_np1'             , long_name = 'Predicted ice thickness', units = 'm')
    CALL add_field_mesh_dp_2D( pc%restart_filename, ncid, 'dHi_dt_Hi_star_np1_u_np1', long_name = 'Thinning rates for predicted ice thickness and updated velocity', units = 'm/yr')
    CALL add_field_mesh_dp_2D( pc%restart_filename, ncid, 'Hi_np1'                  , long_name = 'Corrected ice thickness', units = 'm')
    CALL add_field_mesh_dp_2D( pc%restart_filename, ncid, 'tau_np1'                 , long_name = 'Truncation error', units = 'm')
    CALL add_field_dp_0D(      pc%restart_filename, ncid, 'eta_n'                   , long_name = 'Previous maximum truncation error', units = 'm')
    CALL add_field_dp_0D(      pc%restart_filename, ncid, 'eta_np1'                 , long_name = 'Current maximum truncation error', units = 'm')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_pc_scheme

  SUBROUTINE remap_pc_scheme( mesh_old, mesh_new, pc)
    ! Reallocate memory for the ice thickness predictor/corrector scheme.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh_old
    TYPE(type_mesh),               INTENT(IN)    :: mesh_new
    TYPE(type_ice_pc),             INTENT(INOUT) :: pc

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'remap_pc_scheme'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Reallocate memory
    CALL reallocate_bounds( pc%dHi_dt_Hi_nm1_u_nm1     , mesh_new%vi1, mesh_new%vi2)           ! [m/yr] Thinning rates from previous time step
    CALL reallocate_bounds( pc%dHi_dt_Hi_n_u_n         , mesh_new%vi1, mesh_new%vi2)           ! [m/yr] Thinning rates for current time step with old geometry
    CALL reallocate_bounds( pc%Hi_star_np1             , mesh_new%vi1, mesh_new%vi2)           ! [m]    Predicted ice thickness
    CALL reallocate_bounds( pc%dHi_dt_Hi_star_np1_u_np1, mesh_new%vi1, mesh_new%vi2)           ! [m/yr] Thinning rates for predicted ice thickness and updated velocity
    CALL reallocate_bounds( pc%Hi_np1                  , mesh_new%vi1, mesh_new%vi2)           ! [m]    Corrected ice thickness
    CALL reallocate_bounds( pc%tau_np1                 , mesh_new%vi1, mesh_new%vi2)           ! [m]    Truncation error
    CALL reallocate_bounds( pc%tau_n_guilty            , mesh_new%vi1, mesh_new%vi2)           ! [-]    Number of events above error tolerance

    ! Reinitialise everything from scratch
    pc%dt_n                     = C%dt_ice_min
    pc%dt_np1                   = C%dt_ice_min
    pc%zeta_t                   = 1._dp
    pc%dHi_dt_Hi_nm1_u_nm1      = 0._dp
    pc%dHi_dt_Hi_n_u_n          = 0._dp
    pc%Hi_star_np1              = 0._dp
    pc%dHi_dt_Hi_star_np1_u_np1 = 0._dp
    pc%Hi_np1                   = 0._dp
    pc%tau_np1                  = C%pc_epsilon
    pc%tau_n_guilty             = 0
    pc%eta_n                    = C%pc_epsilon
    pc%eta_np1                  = C%pc_epsilon

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_pc_scheme

! ===== Direct scheme =====
! =========================

  SUBROUTINE run_ice_dynamics_model_direct( region, dt_max)
    ! Calculate a new next modelled ice thickness

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region
    REAL(dp),                               INTENT(IN)    :: dt_max

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ice_dynamics_model_direct'
    REAL(dp)                                              :: dt_crit_SIA, dt_crit_adv, dt
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. (C%choice_stress_balance_approximation == 'SIA' .OR. &
               C%choice_stress_balance_approximation == 'SSA' .OR. &
               C%choice_stress_balance_approximation == 'SIA/SSA')) THEN
      CALL crash('direct timestepping only works for SIA, SSA, or SIA/SSA ice dynamics!')
    END IF

    ! Store previous ice model state
    region%ice%t_Hi_prev  = region%ice%t_Hi_next
    region%ice%Hi_prev    = region%ice%Hi_next

    ! Calculate ice velocities
    CALL solve_stress_balance( region%mesh, region%ice, region%BMB%BMB, region%name)

    ! Calculate time step

    ! Start with the maximum allowed time step
    dt = dt_max

    ! Limit to the SIA critical time step
    IF (C%choice_stress_balance_approximation == 'SIA' .OR. &
        C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      CALL calc_critical_timestep_SIA( region%mesh, region%ice, dt_crit_SIA)
      dt = MIN( dt, dt_crit_SIA)
    END IF

    ! Limit to the advective critical time step
    IF (C%choice_stress_balance_approximation == 'SSA' .OR. &
        C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_adv)
      dt = MIN( dt, dt_crit_adv)
    END IF

    ! Limit to the smallest allowed time step
    dt = MAX( C%dt_ice_min, dt)

    ! Calculate thinning rates and predicted geometry
    CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, region%AMB%AMB, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt, region%ice%dHi_dt, region%ice%Hi_next, region%ice%divQ, region%ice%dHi_dt_target)

    ! If so desired, invert/adjust mass balance fluxes to get an equilibrium state
    CALL MB_inversion( region%mesh, region%ice, region%refgeo_PD, region%SMB, region%BMB, region%LMB, region%AMB, region%ice%dHi_dt, region%ice%Hi_next, dt, region%time, region%name)

    ! Save the "raw" dynamical dH/dt before any alterations
    region%ice%dHi_dt_raw = region%ice%dHi_dt

    ! Modify predicted ice thickness if desired
    CALL alter_ice_thickness( region%mesh, region%ice, region%ice%Hi_prev, region%ice%Hi_next, region%refgeo_PD, region%time)

    ! Compute residual between the "raw" and final thinning rates
    region%ice%dHi_dt_residual = region%ice%dHi_dt_raw - (region%ice%Hi_next - region%ice%Hi_prev) / dt

    ! Set next modelled ice thickness timestamp
    region%ice%t_Hi_next = region%ice%t_Hi_prev + dt

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_model_direct

  SUBROUTINE calc_critical_timestep_SIA( mesh, ice, dt_crit_SIA)
    ! Calculate the critical time step for advective ice flow (CFL criterion)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_SIA'
    REAL(dp), DIMENSION(mesh%nV)                       :: Hi_tot
    INTEGER                                            :: ti, via, vib, vic
    REAL(dp)                                           :: d_ab, d_bc, d_ca, d_min, Hi, D_SIA, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather global ice thickness
    CALL gather_to_all_dp_1D( ice%Hi, Hi_tot)

    ! Initialise time step with maximum allowed value
    dt_crit_SIA = C%dt_ice_max

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate shortest triangle side
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      d_ab = NORM2( mesh%V( vib,:) - mesh%V( via,:))
      d_bc = NORM2( mesh%V( vic,:) - mesh%V( vib,:))
      d_ca = NORM2( mesh%V( via,:) - mesh%V( vic,:))

      d_min = MINVAL([ d_ab, d_bc, d_ca])

      ! Find maximum diffusivity in the vertical column
      D_SIA = MAX( 1E2_dp, MAXVAL( ABS( ice%SIA%D_3D_b( ti,:))))

      ! Calculate critical timestep
      Hi = MAXVAL( [0.1_dp, Hi_tot( via), Hi_tot( vib), Hi_tot( vic)])
      dt = d_min**2 / (6._dp * Hi * D_SIA) * dt_correction_factor
      dt_crit_SIA = MIN( dt_crit_SIA, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_SIA = MIN( C%dt_ice_max, dt_crit_SIA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_SIA

  SUBROUTINE calc_critical_timestep_adv( mesh, ice, dt_crit_adv)
    ! Calculate the critical time step for advective ice flow (CFL criterion)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_adv

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_adv'
    REAL(dp), DIMENSION(mesh%nV)                       :: Hi_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_floating_ice_tot
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2)             :: u_vav_c, v_vav_c
    REAL(dp), DIMENSION(mesh%nE)                       :: u_vav_c_tot, v_vav_c_tot
    INTEGER                                            :: ei, vi, vj
    REAL(dp)                                           :: dist, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather global ice thickness
    CALL gather_to_all_dp_1D( ice%Hi, Hi_tot)

    CALL gather_to_all_logical_1D( ice%mask_floating_ice, mask_floating_ice_tot)

    ! Calculate vertically averaged ice velocities on the edges
    CALL map_velocities_from_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_vav_c, v_vav_c)
    CALL gather_to_all_dp_1D( u_vav_c, u_vav_c_tot)
    CALL gather_to_all_dp_1D( v_vav_c, v_vav_c_tot)

    ! Initialise time step with maximum allowed value
    dt_crit_adv = C%dt_ice_max

    DO ei = mesh%ei1, mesh%ei2

      ! Only check at ice-covered vertices
      vi = mesh%EV( ei,1)
      vj = mesh%EV( ei,2)
      IF (Hi_tot( vi) == 0._dp .OR. Hi_tot( vj) == 0._dp) CYCLE

      IF (C%do_grounded_only_adv_dt) THEN
        ! Only check grounded vertices
        IF (mask_floating_ice_tot( vi) .OR. mask_floating_ice_tot( vj)) CYCLE
      END IF

      dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
      dt = dist / MAX( 0.1_dp, ABS( u_vav_c_tot( ei)) + ABS( v_vav_c_tot( ei))) * dt_correction_factor
      dt_crit_adv = MIN( dt_crit_adv, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN( C%dt_ice_max, dt_crit_adv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_adv

! ===== Inversions =====

  SUBROUTINE SMB_inversion( region, dt)
    ! Invert the surface mass balance that would keep the ice sheet in check

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'SMB_inversion'
    INTEGER                                               :: vi
    INTEGER,  DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: extrapolation_mask
    REAL(dp)                                              :: dt_dummy
    REAL(dp), DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, dHi_dt_dummy, Hi_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this inversion is desired
    IF (.NOT. C%do_SMB_removal_icefree_land) THEN
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! And exit subroutine
      RETURN
    END IF

    ! Initialise
    extrapolation_mask = 0

    ! == Equilibrium SMB
    ! ==================

    ! Set dummy mass balance terms to 0
    SMB_dummy    = 0._dp
    BMB_dummy    = 0._dp
    LMB_dummy    = 0._dp
    AMB_dummy    = 0._dp

    ! Copy model time step
    dt_dummy = dt

    ! Use full mass balance to invert SMB values
    CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, AMB_dummy, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! Compute equilibrium LMB
    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where land should not necessarily be ice-free
      IF (.NOT. region%ice%mask_icefree_land( vi) .OR. .NOT. region%refgeo_PD%Hi( vi) == 0._dp) CYCLE

      ! Equilibrium SMB field
      region%SMB%SMB( vi) = region%SMB%SMB(vi) - dHi_dt_dummy( vi)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE SMB_inversion

  SUBROUTINE BMB_inversion( region, dt)
    ! Invert the basal mass balance that would keep the ice shelves in equilibrium

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'BMB_inversion'
    INTEGER                                               :: vi
    INTEGER,  DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: extrapolation_mask
    REAL(dp)                                              :: dt_dummy
    REAL(dp), DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, dHi_dt_dummy, Hi_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this inversion is desired
    IF (C%do_BMB_inversion .AND. &
        region%time >= C%BMB_inversion_t_start .AND. &
        region%time <= C%BMB_inversion_t_end) THEN
      ! Go for it
    ELSE
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! And exit subroutine
      RETURN
    END IF

    ! Initialise
    extrapolation_mask = 0

    ! == Equilibrium LMB
    ! ==================

    ! Set dummy mass balance terms to 0
    SMB_dummy    = 0._dp
    BMB_dummy    = 0._dp
    LMB_dummy    = 0._dp
    AMB_dummy    = 0._dp

    ! Copy model time step
    dt_dummy = dt

    ! Use no basal or lateral mass balance to invert BMB values for an ice shelf in equilibrium
    CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, BMB_dummy, LMB_dummy, AMB_dummy, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! Initialise
    region%BMB%BMB_inv = 0._dp

    ! Compute equilibrium LMB
    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where BMB does not operate
      IF (.NOT. region%ice%mask_gl_gr( vi) .AND. &
          .NOT. region%ice%mask_floating_ice( vi) .AND. &
          .NOT. region%ice%mask_cf_fl( vi)) CYCLE

      ! Equilibrium BMB field
      region%BMB%BMB_inv( vi) = -dHi_dt_dummy( vi)

      ! Add to extrapolation seeds
      extrapolation_mask( vi) = 2

    END DO

    ! == Calving fronts
    ! =================

    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Detect shelf fronts where upstream BMB can be extrapolated into
      IF (region%ice%mask_cf_fl( vi) .AND. .NOT. region%ice%mask_gl_fl( vi)) THEN
        extrapolation_mask( vi) = 1
      END IF

    END DO

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    CALL extrapolate_Gaussian( region%mesh, extrapolation_mask, region%BMB%BMB_inv, 10000._dp)

    ! == Total BMB
    ! ============

    ! Initialise
    region%BMB%BMB = 0._dp

    ! Compute total BMB
    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where BMB does not operate
      IF (.NOT. region%ice%mask_gl_gr( vi) .AND. &
          .NOT. region%ice%mask_floating_ice( vi) .AND. &
          .NOT. region%ice%mask_cf_fl( vi)) CYCLE

      ! Final BMB field
      region%BMB%BMB( vi) = region%BMB%BMB_inv( vi)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE BMB_inversion

  SUBROUTINE LMB_inversion( region, dt)
    ! Invert the lateral mass balance that would keep the calving front in equilibrium

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'LMB_inversion'
    INTEGER                                               :: vi, vj, ci, ei
    REAL(dp), DIMENSION(region%mesh%ei1:region%mesh%ei2)  :: u_vav_c, v_vav_c
    REAL(dp), DIMENSION(region%mesh%nE)                   :: u_vav_c_tot, v_vav_c_tot
    REAL(dp)                                              :: dt_dummy, calving_rate, calving_ratio, D_x, D_y, D, calving_perp, L_c, V_calved
    REAL(dp), DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, dHi_dt_dummy, Hi_dummy, divQ_eff, LMB_trans
    REAL(dp), DIMENSION(region%mesh%nV)                   :: Hi_tot, fraction_margin_tot
    LOGICAL                                               :: found_advancing_calving_front, found_calving_front_neighbour
    LOGICAL,  DIMENSION(region%mesh%vi1:region%mesh%vi2)  :: mask_advancing_calving_front
    LOGICAL,  DIMENSION(region%mesh%nV)                   :: mask_advancing_calving_front_tot, mask_cf_fl_tot, mask_icefree_ocean_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this inversion is desired
    IF (C%do_LMB_inversion .AND. &
        region%time >= C%LMB_inversion_t_start .AND. &
        region%time <= C%LMB_inversion_t_end) THEN
      ! Go for it
    ELSE
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! And exit subroutine
      RETURN
    END IF

    ! == Equilibrium LMB
    ! ==================

    ! Set dummy mass balance terms to 0
    SMB_dummy    = 0._dp
    BMB_dummy    = 0._dp
    LMB_dummy    = 0._dp
    AMB_dummy    = 0._dp

    ! Copy model time step
    dt_dummy = dt

    ! Use no lateral mass balance to invert its value for a calving front in equilibrium
    CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, LMB_dummy, AMB_dummy, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! Initialise
    region%LMB%LMB_inv = 0._dp

    ! Compute equilibrium LMB
    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where LMB does not operate
      IF (.NOT. region%ice%mask_cf_fl( vi) .AND. &
          .NOT. region%ice%mask_cf_gr( vi) .AND. &
          .NOT. region%ice%mask_icefree_ocean( vi)) CYCLE

      ! Equilibrium LMB field: let positive values remain, since the final goal is
      ! that the *sum* of the equilibrium and transient LMB be equal to the target rate
      region%LMB%LMB_inv( vi) = -dHi_dt_dummy( vi)

    END DO

    ! == Transient LMB
    ! ================

    ! Gather data from all processes
    CALL gather_to_all_dp_1D(      region%ice%Hi_eff, Hi_tot)
    CALL gather_to_all_dp_1D(      region%ice%fraction_margin, fraction_margin_tot)

    ! Gather masks from all processes
    CALL gather_to_all_logical_1D( region%ice%mask_cf_fl, mask_cf_fl_tot)
    CALL gather_to_all_logical_1D( region%ice%mask_icefree_ocean, mask_icefree_ocean_tot)

    ! Calculate vertically averaged ice velocities on the edges
    CALL map_velocities_from_b_to_c_2D( region%mesh, region%ice%u_vav_b, region%ice%v_vav_b, u_vav_c, v_vav_c)
    CALL gather_to_all_dp_1D( u_vav_c, u_vav_c_tot)
    CALL gather_to_all_dp_1D( v_vav_c, v_vav_c_tot)

    ! Initialise
    LMB_trans = 0._dp

    ! ! Translate imposed transient calving rates into LMB
    ! DO vi = region%mesh%vi1, region%mesh%vi2

    !   V_calved = 0._dp

    !   ! Skip vertices where LMB does not operate
    !   IF (.NOT. region%ice%mask_cf_fl( vi) .AND. .NOT. region%ice%mask_icefree_ocean( vi)) CYCLE

    !   calving_rate = 0._dp

    !   IF (C%choice_refgeo_PD_ANT == 'idealised' .AND. region%time > 10000._dp) THEN
    !     IF (C%choice_refgeo_init_idealised == 'calvmip_circular') THEN
    !       calving_rate = -300._dp * SIN(2._dp * pi * region%time / 1000._dp)
    !     ELSEIF (C%choice_refgeo_init_idealised == 'calvmip_Thule') THEN
    !       calving_rate = -750._dp * SIN(2._dp * pi * region%time / 1000._dp)
    !     END IF
    !   END IF

    !   IF (region%ice%uabs_vav( vi) > 0._dp) THEN
    !     calving_ratio = calving_rate / region%ice%uabs_vav( vi)
    !   END IF

    !   ! Loop over all connections of vertex vi
    !   DO ci = 1, region%mesh%nC( vi)

    !     ! Connection ci from vertex vi leads through edge ei to vertex vj
    !     ei = region%mesh%VE( vi,ci)
    !     vj = region%mesh%C(  vi,ci)

    !     ! The shared Voronoi cell boundary section between the
    !     ! Voronoi cells of vertices vi and vj has length L_c
    !     L_c = region%mesh%Cw( vi,ci)

    !     ! Calculate calving rate component perpendicular to this shared Voronoi cell boundary section
    !     D_x = region%mesh%V( vj,1) - region%mesh%V( vi,1)
    !     D_y = region%mesh%V( vj,2) - region%mesh%V( vi,2)
    !     D   = SQRT( D_x**2 + D_y**2)
    !     calving_perp = calving_ratio * ABS( u_vav_c_tot( ei) * D_x/D + v_vav_c_tot( ei) * D_y/D)

    !     ! Calving front vertices: check if neighbour is open ocean
    !     IF (region%ice%mask_cf_fl( vi) .AND. mask_icefree_ocean_tot( vj)) THEN

    !       ! Volume calved laterally: perpendicular calving rate times area of the ice front face [m^3/yr]
    !       V_calved = V_calved + L_c * calving_perp * Hi_tot( vi)

    !     ! Ice-free ocean vertices: check if neighbour is a fully advanced calving front
    !     ELSEIF (region%ice%mask_icefree_ocean( vi) .AND. mask_cf_fl_tot( vj) .AND. fraction_margin_tot( vj) >= 1._dp) THEN

    !       ! Volume calved laterally: perpendicular calving rate times area of the ice front face [m^3/yr]
    !       V_calved = V_calved + L_c * calving_perp * Hi_tot( vj)

    !     END IF

    !   END DO

    !   ! Translate lateral volume loss into vertical thinning rate [m/yr]
    !   LMB_trans( vi) = V_calved / region%mesh%A( vi)

    ! END DO

    ! == Total LMB
    ! ============

    ! Initialise
    region%LMB%LMB = 0._dp

    ! Compute total LMB
    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Skip vertices where LMB does not operate
      IF (.NOT. region%ice%mask_cf_fl( vi) .AND. &
          .NOT. region%ice%mask_cf_gr( vi) .AND. &
          .NOT. region%ice%mask_icefree_ocean( vi)) CYCLE

      ! Final LMB field: now _this_ one should never be positive
      region%LMB%LMB( vi) = MIN( 0._dp, region%LMB%LMB_inv( vi) + LMB_trans( vi))

    END DO

    ! ! == Effective ice divergence
    ! ! ===========================

    ! ! Set dummy mass balance terms to 0
    ! SMB_dummy    = 0._dp
    ! BMB_dummy    = 0._dp
    ! LMB_dummy    = 0._dp
    ! AMB_dummy    = 0._dp

    ! ! Copy model time step
    ! dt_dummy = dt

    ! ! Use no mass balance to get an estimate of the effective flux divergence
    ! CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, SMB_dummy, BMB_dummy, LMB_dummy, AMB_dummy, region%ice%fraction_margin, &
    !                   region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! ! Effective flux divergence
    ! divQ_eff = -dHi_dt_dummy

    ! ! Initialise mask of advancing fronts
    ! mask_advancing_calving_front = .FALSE.

    ! ! Identify advancing calving fronts
    ! DO vi = region%mesh%vi1, region%mesh%vi2
    !   IF (region%ice%mask_icefree_ocean( vi) .AND. divQ_eff(vi) < 0._dp) THEN
    !     mask_advancing_calving_front( vi) = .TRUE.
    !   END IF
    ! END DO

    ! ! == Valid areas
    ! ! ==============

    ! ! Set dummy mass balance terms to 0
    ! SMB_dummy    = 0._dp
    ! BMB_dummy    = 0._dp
    ! LMB_dummy    = 0._dp
    ! AMB_dummy    = 0._dp

    ! ! Copy model time step
    ! dt_dummy = dt

    ! ! Use total mass balance to check whether the advancing calving front will survive the incoming lateral mass balance
    ! CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, AMB_dummy, region%ice%fraction_margin, &
    !                   region%ice%mask_noice, dt_dummy, dHi_dt_dummy, Hi_dummy, region%ice%divQ, region%ice%dHi_dt_target)

    ! ! Check predicted dHi/dt
    ! DO vi = region%mesh%vi1, region%mesh%vi2
    !   IF (region%ice%mask_icefree_ocean( vi) .AND. dHi_dt_dummy( vi) <= 0._dp) THEN

    !     ! It will not, so do not consider this point an advancing front
    !     mask_advancing_calving_front( vi) = .FALSE.

    !     ! Apply only equilibrium LMB here
    !     region%LMB%LMB( vi) = MIN( 0._dp, region%LMB%LMB_inv( vi))

    !   END IF
    ! END DO

    ! ! Gather advancing and floating calving front masks from all processes
    ! CALL gather_to_all_logical_1D( mask_advancing_calving_front, mask_advancing_calving_front_tot)

    ! ! Identify vertices where LMB will operate
    ! DO vi = region%mesh%vi1, region%mesh%vi2

    !   ! Valid calving front vertices
    !   IF (region%ice%mask_cf_fl( vi)) THEN
    !     ! Initialise flag
    !     found_advancing_calving_front = .FALSE.

    !     ! Check for advancing front neighbours
    !     DO ci = 1, region%mesh%nC( vi)
    !       vj = region%mesh%C( vi,ci)
    !       IF (mask_advancing_calving_front_tot( vj)) THEN
    !         found_advancing_calving_front = .TRUE.
    !         EXIT
    !       END IF
    !     END DO

    !     IF (found_advancing_calving_front) THEN
    !       ! Do not apply LMB here, since it will be applied on its advancing neighbour
    !       region%LMB%LMB( vi) = 0._dp
    !     END IF

    !   END IF

    !   ! Valid ocean vertices
    !   IF (region%ice%mask_icefree_ocean( vi)) THEN
    !     ! Initialise flag
    !     found_calving_front_neighbour = .FALSE.

    !     ! Check for calving front neighbours
    !     DO ci = 1, region%mesh%nC( vi)
    !       vj = region%mesh%C( vi,ci)
    !       IF (mask_cf_fl_tot( vj)) THEN
    !         found_calving_front_neighbour = .TRUE.
    !         EXIT
    !       END IF
    !     END DO

    !     IF (.NOT. found_calving_front_neighbour) THEN
    !       ! No calving front neighbours: apply only equilibrium LMB here
    !       region%LMB%LMB( vi) = MIN( 0._dp, region%LMB%LMB_inv( vi))
    !     END IF

    !   END IF

    ! END DO

    ! == DENK DROM
    ! ============

    ! IF (C%choice_refgeo_PD_ANT == 'idealised' .AND. &
    !      (C%choice_refgeo_init_idealised == 'calvmip_circular' .OR. &
    !       C%choice_refgeo_init_idealised == 'calvmip_Thule')) THEN

    !   IF (region%time >= 6000._dp) THEN
    !     IF (C%choice_regions_of_interest == 'CalvMIP_quarter') THEN
    !       C%ROI_maximum_resolution_grounding_line = 5000._dp
    !       ! C%ROI_maximum_resolution_calving_front  = 8000._dp
    !       ! C%ROI_maximum_resolution_floating_ice   = 10000._dp
    !       ! C%ROI_maximum_resolution_grounded_ice   = 20000._dp
    !     ELSE
    !       C%maximum_resolution_grounding_line = 5000._dp
    !       ! C%maximum_resolution_calving_front  = 8000._dp
    !       ! C%maximum_resolution_floating_ice   = 10000._dp
    !     END IF
    !   END IF

    !   IF (region%time >= 7000._dp) THEN
    !     IF (C%choice_regions_of_interest == 'CalvMIP_quarter') THEN
    !       C%ROI_maximum_resolution_grounding_line = 3000._dp
    !       ! C%ROI_maximum_resolution_calving_front  = 5000._dp
    !       ! C%ROI_maximum_resolution_floating_ice   = 8000._dp
    !       ! C%ROI_maximum_resolution_grounded_ice   = 16000._dp
    !     ELSE
    !       C%maximum_resolution_grounding_line = 3000._dp
    !       ! C%maximum_resolution_calving_front  = 5000._dp
    !       ! C%maximum_resolution_floating_ice   = 8000._dp
    !     END IF
    !   END IF

    !   IF (region%time >= 9000._dp) THEN
    !     C%allow_mesh_updates = .FALSE.
    !     IF (C%choice_refgeo_init_idealised == 'calvmip_circular') THEN
    !       C%calving_threshold_thickness_shelf = 10._dp
    !     END IF
    !   ELSE
    !     DO vi = region%mesh%vi1, region%mesh%vi2
    !       IF (SQRT(region%mesh%V( vi,1)**2 + region%mesh%V( vi,2)**2) < 750000._dp) THEN
    !         region%LMB%LMB( vi) = 0._dp
    !       ELSE
    !         region%LMB%LMB( vi) = -100._dp
    !       END IF
    !     END DO
    !   END IF

    ! END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE LMB_inversion

END MODULE ice_model_main
