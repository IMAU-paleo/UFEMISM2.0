MODULE ice_model_main

  ! The main ice-dynamical model module.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF, &
                                                                     save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE scalar_types                                           , ONLY: type_regional_scalars
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_pc
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE region_types                                           , ONLY: type_model_region
  USE ice_model_memory                                       , ONLY: allocate_ice_model
  USE ice_model_utilities                                    , ONLY: determine_masks, calc_bedrock_CDFs, calc_grounded_fractions, calc_zeta_gradients, &
                                                                     calc_mask_noice
  USE ice_model_scalars                                      , ONLY: calc_ice_model_scalars
  USE ice_thickness                                          , ONLY: calc_dHi_dt
  USE math_utilities                                         , ONLY: ice_surface_elevation, thickness_above_floatation
  USE geothermal_heat_flux                                   , ONLY: initialise_geothermal_heat_flux
  USE basal_hydrology                                        , ONLY: initialise_basal_hydrology_model
  USE bed_roughness                                          , ONLY: initialise_bed_roughness
  USE ice_velocity_main                                      , ONLY: initialise_velocity_solver, solve_stress_balance, remap_velocity_solver, &
                                                                     create_restart_file_ice_velocity, write_to_restart_file_ice_velocity, &
                                                                     map_velocities_from_b_to_c_2D
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file, open_existing_netcdf_file_for_writing
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, add_time_dimension_to_file, &
                                                                     add_field_dp_0D, add_field_mesh_dp_2D, write_time_to_file, write_to_field_multopt_mesh_dp_2D, &
                                                                     write_to_field_multopt_dp_0D
  USE netcdf_input                                           , ONLY: read_field_from_file_0D, read_field_from_mesh_file_2D

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

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Calculate the no-ice mask
    CALL calc_mask_noice( region%mesh, region%ice)

    ! NOTE: as calculating the zeta gradients is quite expensive, only do so when necessary,
    !       i.e. when solving the heat equation or the Blatter-Pattyn stress balance
!    ! Calculate zeta gradients
!    CALL calc_zeta_gradients( region%mesh, region%ice)

    ! Calculate sub-grid grounded-area fractions
    CALL calc_grounded_fractions( region%mesh, region%ice)

    ! Calculate ice-sheet integrated values (total volume, area, etc.)
    CALL calc_ice_model_scalars( region%mesh, region%ice, region%refgeo_PD, region%scalars)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_model

  SUBROUTINE initialise_ice_dynamics_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ! Initialise all data fields of the ice module

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    TYPE(type_ice_model),          INTENT(INOUT) :: ice
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo_init
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo_PD
    TYPE(type_regional_scalars),   INTENT(OUT)   :: scalars
    CHARACTER(LEN=3),              INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_ice_dynamics_model'
    INTEGER                                      :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      WRITE(*,"(A)") '  Initialising ice dynamics model...'
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

    DO vi = mesh%vi1, mesh%vi2

      ! Basic geometry
      ice%Hi ( vi) = refgeo_init%Hi( vi)
      ice%Hb ( vi) = refgeo_init%Hb( vi)
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

    ! Calculate the no-ice mask
    CALL calc_mask_noice( mesh, ice)

    ! Sub-grid fractions
    ! ==================

    ! Compute bedrock cumulative density function
    CALL calc_bedrock_CDFs( mesh, refgeo_PD, ice)
    ! Initialise sub-grid grounded-area fractions
    CALL calc_grounded_fractions( mesh, ice)

    ! Basal conditions
    ! ================

    ! Allocate and initialise basal conditions
    CALL initialise_geothermal_heat_flux(  mesh, ice)
    CALL initialise_basal_hydrology_model( mesh, ice)
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

    ! Ice-sheet-wide scalars
    ! ======================

    CALL calc_ice_model_scalars( mesh, ice, refgeo_PD, scalars)

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
    REAL(dp), DIMENSION(:    ), ALLOCATABLE               :: Hi_dummy
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( Hi_dummy( region%mesh%vi1:region%mesh%vi2))

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

    ! Calculate time step ratio
    region%ice%pc%zeta_t = region%ice%pc%dt_np1 / region%ice%pc%dt_n

  ! == Predictor step ==
  ! ====================

    ! Store thinning rates from previous time step
    region%ice%pc%dHi_dt_Hi_nm1_u_nm1 = region%ice%dHi_dt

    ! Calculate thinning rates for current geometry and velocity
    CALL calc_dHi_dt( region%mesh, region%ice%Hi_prev, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, &
      region%ice%mask_noice, region%ice%pc%dt_np1, region%ice%pc%dHi_dt_Hi_n_u_n, Hi_dummy)

    ! Calculate predicted ice thickness (Robinson et al., 2020, Eq. 30)
    region%ice%pc%Hi_star_np1 = region%ice%Hi_prev + region%ice%pc%dt_np1 * ((1._dp + region%ice%pc%zeta_t / 2._dp) * &
      region%ice%pc%dHi_dt_Hi_n_u_n - (region%ice%pc%zeta_t / 2._dp) * region%ice%pc%dHi_dt_Hi_nm1_u_nm1)

  ! == Update step ==
  ! =================

    ! Set model geometry to predicted
    region%ice%Hi = region%ice%pc%Hi_star_np1
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%Hs( vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
    END DO

    ! Calculate ice velocities for the predicted geometry
    CALL solve_stress_balance( region%mesh, region%ice, region%BMB%BMB)

  ! == Corrector step ==
  ! ====================

    ! Set model geometry back to original
    region%ice%Hi = region%ice%Hi_prev
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%Hs( vi) = ice_surface_elevation( region%ice%Hi( vi), region%ice%Hb( vi), region%ice%SL( vi))
    END DO
    ! Update masks
    CALL determine_masks( region%mesh, region%ice)

    ! Calculate thinning rates for the predicted ice thickness and updated velocity
    CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, &
      region%ice%mask_noice, region%ice%pc%dt_np1, region%ice%pc%dHi_dt_Hi_star_np1_u_np1, Hi_dummy)

    ! Calculate corrected ice thickness (Robinson et al. (2020), Eq. 31)
    region%ice%pc%Hi_np1 = region%ice%Hi_prev + (region%ice%pc%dt_np1 / 2._dp) * (region%ice%pc%dHi_dt_Hi_n_u_n + region%ice%pc%dHi_dt_Hi_star_np1_u_np1)

    ! Estimate truncation error
    CALL calc_pc_truncation_error( region%mesh, region%ice%pc)

    ! Set next modelled ice thickness
    region%ice%t_Hi_next = region%ice%t_Hi_prev + region%ice%pc%dt_np1
    region%ice%Hi_next   = region%ice%pc%Hi_np1

    ! Clean up after yourself
    DEALLOCATE( Hi_dummy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_model_pc

  SUBROUTINE calc_pc_truncation_error( mesh, pc)
    ! Calculate the truncation error tau in the ice thickness rate of change (Robinson et al., 2020, Eq. 32)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
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

    ! Store the previous maximum truncation error eta_n
    pc%eta_n = pc%eta_np1

    ! Calculate the maximum truncation error eta
    pc%eta_np1 = MAX( C%pc_eta_min, MAXVAL( pc%tau_np1))
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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_pc_scheme

  SUBROUTINE initialise_pc_scheme_from_file( pc, filename, timeframe)
    ! Initialise values for the ice thickness predictor/corrector scheme from a (restart) file.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_ice_pc),             INTENT(OUT)   :: pc
    CHARACTER(LEN=256),            INTENT(IN)    :: filename
    REAL(dp),                      INTENT(IN)    :: timeframe

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_pc_scheme_from_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write to terminal
    IF (par%master) WRITE(0,*) '   Initialising ice thickness predictor/corrector scheme from file "' // colour_string( TRIM( filename),'light blue') // '"...'

    ! Read velocities from the file
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
    IF (par%master) WRITE(0,'(A)') '  Creating ice dynamics restart file "' // &
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
    CALL solve_stress_balance( region%mesh, region%ice, region%BMB%BMB)

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
    CALL calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, &
      region%ice%mask_noice, dt, region%ice%dHi_dt, region%ice%Hi_next)

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
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2)             :: u_vav_c, v_vav_c
    REAL(dp), DIMENSION(mesh%nE)                       :: u_vav_c_tot, v_vav_c_tot
    INTEGER                                            :: ei, vi, vj
    REAL(dp)                                           :: dist, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather global ice thickness
    CALL gather_to_all_dp_1D( ice%Hi, Hi_tot)

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

      dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
      dt = dist / MAX( 0.1_dp, ABS( u_vav_c_tot( ei)) + ABS( v_vav_c_tot( ei))) * dt_correction_factor
      dt_crit_adv = MIN( dt_crit_adv, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN( C%dt_ice_max, dt_crit_adv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_adv

END MODULE ice_model_main
