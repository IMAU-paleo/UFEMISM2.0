MODULE UFEMISM_main_model

  ! The main regional model

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync, ierr
  USE control_resources_and_error_messaging                  , ONLY: happy, warning, crash, init_routine, finalise_routine, colour_string, str2int, int2str, &
                                                                     insert_val_into_string_dp
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE netcdf_debug                                           , ONLY: write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF, &
                                                                     save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D
  USE region_types                                           , ONLY: type_model_region
  USE ice_model_types                                        , ONLY: type_ice_model
  USE mesh_types                                             , ONLY: type_mesh
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE reference_geometries                                   , ONLY: initialise_reference_geometries_raw, initialise_reference_geometries_on_model_mesh
  USE ice_model_main                                         , ONLY: initialise_ice_dynamics_model, run_ice_dynamics_model, remap_ice_dynamics_model, &
                                                                     create_restart_files_ice_model, write_to_restart_files_ice_model, apply_geometry_relaxation
  USE basal_hydrology                                        , ONLY: run_basal_hydrology_model, initialise_pore_water_fraction_inversion, run_pore_water_fraction_inversion
  USE bed_roughness                                          , ONLY: run_bed_roughness_model
  USE thermodynamics_main                                    , ONLY: initialise_thermodynamics_model, run_thermodynamics_model, &
                                                                     create_restart_file_thermo, write_to_restart_file_thermo
  USE climate_main                                           , ONLY: initialise_climate_model, run_climate_model, remap_climate_model, &
                                                                     create_restart_file_climate_model, write_to_restart_file_climate_model
  USE ocean_main                                             , ONLY: initialise_ocean_model, run_ocean_model, remap_ocean_model, &
                                                                     create_restart_file_ocean_model, write_to_restart_file_ocean_model
  USE SMB_main                                               , ONLY: initialise_SMB_model, run_SMB_model, remap_SMB_model, &
                                                                     create_restart_file_SMB_model, write_to_restart_file_SMB_model
  USE BMB_main                                               , ONLY: initialise_BMB_model, run_BMB_model, remap_BMB_model, &
                                                                     create_restart_file_BMB_model, write_to_restart_file_BMB_model
  USE LMB_main                                               , ONLY: initialise_LMB_model, run_LMB_model, remap_LMB_model
  USE AMB_main                                               , ONLY: initialise_AMB_model, remap_AMB_model
  USE GIA_main                                               , ONLY: initialise_GIA_model, run_GIA_model, remap_GIA_model, &
                                                                     create_restart_file_GIA_model, write_to_restart_file_GIA_model
  USE basal_inversion_main                                   , ONLY: initialise_basal_inversion, run_basal_inversion
  USE netcdf_basic                                           , ONLY: open_existing_netcdf_file_for_reading, close_netcdf_file
  USE netcdf_input                                           , ONLY: setup_mesh_from_file
  USE mesh_creation                                          , ONLY: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success
  USE mesh_operators                                         , ONLY: calc_all_matrix_operators_mesh
  USE grid_basic                                             , ONLY: setup_square_grid
  USE main_regional_output                                   , ONLY: create_main_regional_output_file_mesh,   create_main_regional_output_file_grid, &
                                                                     write_to_main_regional_output_file_mesh, write_to_main_regional_output_file_grid, &
                                                                     create_main_regional_output_file_grid_ROI, write_to_main_regional_output_file_grid_ROI, &
                                                                     create_scalar_regional_output_file, write_to_scalar_regional_output_file
  use mesh_ROI_polygons
  USE math_utilities                                         , ONLY: longest_triangle_leg
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_logical_1D
  use apply_maps, only: clear_all_maps_involving_this_mesh
  USE mesh_memory                                            , ONLY: deallocate_mesh
  USE ice_model_scalars                                      , ONLY: calc_ice_model_scalars

  IMPLICIT NONE

CONTAINS

! ===== Main routine =====
! ========================

  SUBROUTINE run_model_region( region, t_end)
    ! Integrate this model region forward in time until t_end

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(INOUT) :: region
    REAL(dp)                                           , INTENT(IN)    :: t_end    ! [yr]

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: routine_name
    INTEGER                                                            :: ndt_av
    REAL(dp)                                                           :: dt_av
    REAL(dp)                                                           :: mesh_fitness_coefficient

    ! Add routine to path
    routine_name = 'run_model('  //  region%name  //  ')'
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE (0,'(A,A,A,A,A,F9.3,A,F9.3,A)') &
      ' Running model region ', colour_string( region%name, 'light blue'), ' (', colour_string( TRIM( region%long_name), 'light blue'), &
      ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'

    ! Initialise average ice-dynamical time step
    ndt_av = 0
    dt_av  = 0._dp

    ! The main UFEMISM time loop
    main_time_loop: DO WHILE (region%time <= t_end)

      ! == Mesh update code
      IF (C%allow_mesh_updates .AND. &
          region%time > region%time_mesh_was_created + C%dt_mesh_update_min .AND. &
          region%time > C%start_time_of_run) THEN

        ! Calculate the mesh fitness coefficient
        CALL calc_mesh_fitness_coefficient( region%mesh, region%ice, mesh_fitness_coefficient)

        ! If necessary and allowed, perform a mesh update
        IF (mesh_fitness_coefficient < C%minimum_mesh_fitness_coefficient) THEN
          CALL update_mesh( region)
        END IF

      END IF ! IF (C%allow_mesh_updates) THEN

      ! Run the subglacial hydrology model
      CALL run_basal_hydrology_model( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%HIV, region%time)

      ! Run the bed roughness model
      CALL run_bed_roughness_model( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%BIV, region%time)

      ! Run the ice dynamics model to calculate ice geometry at the desired time, and update
      ! velocities, thinning rates, and predicted geometry if necessary
      CALL run_ice_dynamics_model( region)

      ! Calculate ice temperature at the desired time, and update
      ! predicted temperature if necessary
      CALL run_thermodynamics_model( region)

      ! Calculate the climate
      CALL run_climate_model( region%mesh, region%ice, region%climate, region%name, region%time)

      ! Calculate the ocean
      CALL run_ocean_model( region%mesh, region%ice, region%ocean, region%name, region%time)

      ! Calculate the surface mass balance
      CALL run_SMB_model( region%mesh, region%grid_smooth, region%ice, region%climate, region%SMB, region%name, region%time)

      ! Calculate the basal mass balance
      CALL run_BMB_model( region%mesh, region%ice, region%ocean, region%refgeo_PD, region%SMB, region%BMB, region%name, region%time)

      ! Calculate the lateral mass balance
      CALL run_LMB_model( region%mesh, region%ice, region%LMB, region%name, region%time)

      ! Calculate bedrock deformation at the desired time, and update
      ! predicted deformation if necessary
      CALL run_GIA_model( region)

      ! Run the basal inversion model
      IF (C%do_bed_roughness_nudging) THEN
        CALL run_basal_inversion( region)
      END IF

      ! Run the pore water fraction inversion model
      IF (C%do_pore_water_nudging) THEN
        CALL run_pore_water_fraction_inversion( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%HIV, region%time)
      END IF

      ! Calculate ice-sheet integrated values (total volume, area, etc.)
      CALL calc_ice_model_scalars( region%mesh, region%ice, region%SMB, region%BMB, region%LMB, region%refgeo_PD, region%scalars)

      ! Write to the main regional output NetCDF file
      CALL write_to_regional_output_files( region)

      ! Advance this region's time to the time of the next "action"
      CALL advance_region_time_to_time_of_next_action( region, t_end)

      ! Keep track of the average ice-dynamical time step and print it to the terminal
      ndt_av = ndt_av + 1
      dt_av  = dt_av  + (region%ice%t_Hi_next - region%ice%t_Hi_prev)
      IF (par%master .AND. C%do_time_display) THEN
        CALL time_display( region, t_end, dt_av, ndt_av)
      END IF

      ! If we've reached the end of this coupling interval, stop
      IF (region%time == t_end) EXIT main_time_loop

    END DO main_time_loop ! DO WHILE (region%time <= t_end)


    IF (region%time == C%end_time_of_run) THEN
      ! Give all processes time to catch up
      CALL sync
      ! Congrats, you've made it
      IF (par%master) WRITE(0,'(A)') ' Finalising regional simulation...'
      ! Write the final model state to output
      CALL write_to_regional_output_files( region)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_region

! ===== Output =====
! ==================

  SUBROUTINE write_to_regional_output_files( region)
    ! Write to all regional output files

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_regional_output_files'
    INTEGER                                                            :: i
    REAL(dp)                                                           :: t_closest
    LOGICAL                                                            :: do_output_main, do_output_restart, do_output_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine time of next output event
    t_closest = MIN( region%output_t_next, region%output_restart_t_next, region%output_grid_t_next)

    ! Determine actions
    IF     (region%time < t_closest) THEN
      ! It is not yet time to write to output
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (region%time > t_closest) THEN
      ! This should not be possible
      CALL crash('overshot the output time step')
    ELSEIF (region%time == t_closest) THEN
      ! It is time to write to output!
    ELSE
      ! region%time is NaN?
      CALL crash('region%time is apparently NaN')
    END IF

    ! Determine type of next output event
    do_output_main    = .FALSE.
    do_output_restart = .FALSE.
    do_output_grid    = .FALSE.

    ! Update time stamps
    IF (region%time == region%output_t_next) THEN
      region%output_t_next = region%output_t_next + C%dt_output
      do_output_main = .TRUE.
    END IF
    IF (region%time == region%output_restart_t_next) THEN
      region%output_restart_t_next = region%output_restart_t_next + C%dt_output_restart
      do_output_restart = .TRUE.
    END IF
    IF (region%time == region%output_grid_t_next) THEN
      region%output_grid_t_next = region%output_grid_t_next + C%dt_output_grid
      do_output_grid = .TRUE.
    END IF

    ! If needed, create a new set of mesh output files for the current model mesh
    IF (.NOT. region%output_files_match_current_mesh) THEN

      ! Create the main regional output files
      CALL create_main_regional_output_file_mesh( region)

      ! Create the restart files for all the model components
      CALL create_restart_files_ice_model   ( region%mesh, region%ice)
      CALL create_restart_file_thermo       ( region%mesh, region%ice)
      CALL create_restart_file_climate_model( region%mesh, region%climate, region%name)
      CALL create_restart_file_ocean_model  ( region%mesh, region%ocean  , region%name)
      CALL create_restart_file_SMB_model    ( region%mesh, region%SMB    , region%name)
      CALL create_restart_file_BMB_model    ( region%mesh, region%BMB    , region%name)
      CALL create_restart_file_GIA_model    ( region%mesh, region%GIA    , region%name)

      ! Confirm that the current set of mesh output files match the current model mesh
      ! (set to false whenever a new mesh is created,
      ! and set to true whenever a new set of output files is created)

      region%output_files_match_current_mesh = .TRUE.

    END IF ! IF (.NOT. region%output_files_match_current_mesh) THEN

    IF (do_output_main) THEN
      ! Write to the main regional output files
      CALL write_to_main_regional_output_file_mesh( region)

      ! Write to the region-of-interest output files
      DO i = 1, region%nROI
        CALL write_to_main_regional_output_file_grid_ROI( region, region%output_grids_ROI( i), region%output_filenames_grid_ROI( i))
      END DO

      ! Write to scalar regional output file
      CALL write_to_scalar_regional_output_file( region)
    END IF

    IF (do_output_restart) THEN
      ! Write to the restart files for all the model components
      CALL write_to_restart_files_ice_model   ( region%mesh, region%ice                 , region%time)
      CALL write_to_restart_file_thermo       ( region%mesh, region%ice                 , region%time)
      CALL write_to_restart_file_climate_model( region%mesh, region%climate, region%name, region%time)
      CALL write_to_restart_file_ocean_model  ( region%mesh, region%ocean  , region%name, region%time)
      CALL write_to_restart_file_SMB_model    ( region%mesh, region%SMB    , region%name, region%time)
      CALL write_to_restart_file_BMB_model    ( region%mesh, region%BMB    , region%name, region%time)
      CALL write_to_restart_file_GIA_model    ( region%mesh, region%GIA    , region%name, region%time)
    END IF

    IF (do_output_grid) THEN
      ! Write to the gridded regional output file
      CALL write_to_main_regional_output_file_grid( region)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_regional_output_files

! ===== Time stepping =====
! =========================

  SUBROUTINE advance_region_time_to_time_of_next_action( region, t_end)
    ! Advance this region's time to the time of the next "action"
    ! (e.g. running the ice dynamics model, writing to output, running
    ! the SMB model, etc.)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(INOUT) :: region
    REAL(dp)                                           , INTENT(IN)    :: t_end    ! [yr]

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'advance_region_time_to_time_of_next_action'
    REAL(dp)                                                           :: time_of_next_action

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ===== Find time of next action =====
    ! ====================================

    ! Start by setting the time of next action to the end of this coupling interval,
    ! then reduce it over all possible actions to find the first one.

    time_of_next_action = t_end

    ! Ice dynamics
    time_of_next_action = MIN( time_of_next_action, region%ice%t_Hi_next)

    ! Thermodynamics
    time_of_next_action = MIN( time_of_next_action, region%ice%t_Ti_next)

    ! Climate
    time_of_next_action = MIN( time_of_next_action, region%climate%t_next)

    ! Ocean
    time_of_next_action = MIN( time_of_next_action, region%ocean%t_next)

    ! SMB
    time_of_next_action = MIN( time_of_next_action, region%SMB%t_next)

    ! BMB
    time_of_next_action = MIN( time_of_next_action, region%BMB%t_next)

    ! LMB
    time_of_next_action = MIN( time_of_next_action, region%LMB%t_next)

    ! GIA
    time_of_next_action = MIN( time_of_next_action, region%GIA%t_next)

    ! Basal inversion
    IF (C%do_bed_roughness_nudging) THEN
      time_of_next_action = MIN( time_of_next_action, region%BIV%t_next)
    END IF

    ! Hydrology inversion
    IF (C%do_pore_water_nudging) THEN
      time_of_next_action = MIN( time_of_next_action, region%HIV%t_next)
    END IF

    ! Target dHi_dt: make sure we don't overshoot its turnoff time
    IF (C%do_target_dHi_dt) THEN
      IF (C%target_dHi_dt_t_end > region%time) THEN
        time_of_next_action = MIN( time_of_next_action, C%target_dHi_dt_t_end)
      END IF
    END IF

    ! Output
    time_of_next_action = MIN( time_of_next_action, region%output_t_next)
    time_of_next_action = MIN( time_of_next_action, region%output_restart_t_next)
    time_of_next_action = MIN( time_of_next_action, region%output_grid_t_next)

    ! ===== Advance region time =====
    ! ===============================

    region%time = time_of_next_action

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE advance_region_time_to_time_of_next_action

! ===== Model initialisation =====
! ================================

  SUBROUTINE initialise_model_region( region, region_name)
    ! Initialise this model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(OUT)   :: region
    CHARACTER(LEN=3),                                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_model_region'
    CHARACTER(LEN=256)                                                 :: grid_name
    REAL(dp)                                                           :: dx_grid_smooth, dx_grid_output

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set region time
    region%time = C%start_time_of_run

    ! ===== Region name =====
    ! =======================

    ! Which region are we initialising?
    region%name = region_name
    IF     (region%name == 'NAM') THEN
      region%long_name = 'North America'
    ELSEIF (region%name == 'EAS') THEN
      region%long_name = 'Eurasia'
    ELSEIF (region%name == 'GRL') THEN
      region%long_name = 'Greenland'
    ELSEIF (region%name == 'ANT') THEN
      region%long_name = 'Antarctica'
    ELSE
      CALL crash('unknown region "' // TRIM( region%name) // '"')
    END IF

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') ''
    IF (par%master) WRITE(0,'(A)') ' Initialising model region ' // colour_string( region%name,'light blue') // ' (' // &
      colour_string( TRIM( region%long_name),'light blue') // ')...'

    ! ===== Reference geometries =====
    ! ================================

    ! Initialise all the reference geometries on their raw input grids
    CALL initialise_reference_geometries_raw( region%name, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq)

    ! ===== Initial mesh =====
    ! ========================

    ! Set up the first model mesh
    CALL setup_first_mesh( region)

    ! Remap reference geometries from their raw input grids to the model mesh
    CALL initialise_reference_geometries_on_model_mesh( region%name, region%mesh, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq)

    ! ===== Square grid used for smoothing =====
    ! ==========================================

    ! Determine resolution for this region's output grid
    IF     (region%name == 'NAM') THEN
      dx_grid_smooth = C%dx_square_grid_smooth_NAM
    ELSEIF (region%name == 'EAS') THEN
      dx_grid_smooth = C%dx_square_grid_smooth_EAS
    ELSEIF (region%name == 'GRL') THEN
      dx_grid_smooth = C%dx_square_grid_smooth_GRL
    ELSEIF (region%name == 'ANT') THEN
      dx_grid_smooth = C%dx_square_grid_smooth_ANT
    ELSE
      CALL crash('unknown region name "' // region%name // '"!')
    END IF

    ! Create the square smooth grid
    grid_name = 'square_grid_smooth_' // region%name
    CALL setup_square_grid( grid_name, region%mesh%xmin, region%mesh%xmax, region%mesh%ymin, region%mesh%ymax, &
       dx_grid_smooth, region%grid_smooth, &
       lambda_M = region%mesh%lambda_M, phi_M = region%mesh%phi_M, beta_stereo = region%mesh%beta_stereo)

    ! ===== Ice dynamics =====
    ! ========================

    CALL initialise_ice_dynamics_model( region%mesh, region%ice, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%GIA, region%name)

    ! ===== Climate =====
    ! ===================

    CALL initialise_climate_model( region%mesh, region%climate, region%name)

    ! ===== Ocean =====
    ! =================

    CALL initialise_ocean_model( region%mesh, region%ocean, region%name)

    ! ===== Surface mass balance =====
    ! ================================

    CALL initialise_SMB_model( region%mesh, region%SMB, region%name)

    ! ===== Basal mass balance =====
    ! ==============================

    CALL initialise_BMB_model( region%mesh, region%ice, region%BMB, region%name)

    ! ===== Lateral mass balance =====
    ! ================================

    CALL initialise_LMB_model( region%mesh, region%LMB, region%name)

    ! ===== Artificial mass balance =====
    ! ===================================

    CALL initialise_AMB_model( region%mesh, region%AMB)

    ! ===== Run the climate, ocean, SMB, BMB, and LMB models =====
    ! ============================================================

    ! Run the models
    CALL run_climate_model( region%mesh, region%ice, region%climate, region%name, C%start_time_of_run)
    CALL run_ocean_model( region%mesh, region%ice, region%ocean, region%name, C%start_time_of_run)
    CALL run_SMB_model( region%mesh, region%grid_smooth, region%ice, region%climate, region%SMB, region%name, C%start_time_of_run)
    CALL run_BMB_model( region%mesh, region%ice, region%ocean, region%refgeo_PD, region%SMB, region%BMB, region%name, C%start_time_of_run)
    CALL run_LMB_model( region%mesh, region%ice, region%LMB, region%name, region%time)

    ! Reset the timers
    region%climate%t_next = C%start_time_of_run
    region%ocean%t_next   = C%start_time_of_run
    region%SMB%t_next     = C%start_time_of_run
    region%BMB%t_next     = C%start_time_of_run
    region%LMB%t_next     = C%start_time_of_run

    ! ===== Thermodynamics =====
    ! ==========================

    ! NOTE: must be initialised after the climate, ocean, SMB, and BMB have been initialised,
    !       (some of) which are needed to set up the first englacial ice temperature field.

    CALL initialise_thermodynamics_model( region)

    ! ===== Glacial isostatic adjustment =====
    ! ========================================

    CALL initialise_GIA_model( region%mesh, region%GIA, region%name)

    ! ===== Basal inversion =====
    ! ===========================

    IF (C%do_bed_roughness_nudging) THEN
      CALL initialise_basal_inversion( region%mesh, region%ice, region%BIV, region%name)
    END IF

    IF (C%do_pore_water_nudging) THEN
      CALL initialise_pore_water_fraction_inversion( region%mesh, region%ice, region%HIV, region%name)
    END IF

    ! ===== Corrections =====
    ! =======================

    CALL apply_regional_corrections( region)
    CALL apply_geometry_relaxation(  region)

    ! ===== Integrated scalars =====
    ! ==============================

    ! Calculate ice-sheet integrated values (total volume, area, etc.)
    CALL calc_ice_model_scalars( region%mesh, region%ice, region%SMB, region%BMB, region%LMB, region%refgeo_PD, region%scalars)

    ! ===== Regional output =====
    ! ===========================

    ! Determine resolution for this region's output grid
    IF     (region%name == 'NAM') THEN
      dx_grid_output = C%dx_output_grid_NAM
    ELSEIF (region%name == 'EAS') THEN
      dx_grid_output = C%dx_output_grid_EAS
    ELSEIF (region%name == 'GRL') THEN
      dx_grid_output = C%dx_output_grid_GRL
    ELSEIF (region%name == 'ANT') THEN
      dx_grid_output = C%dx_output_grid_ANT
    ELSE
      CALL crash('unknown region name "' // region%name // '"!')
    END IF

    ! Create the square output grid
    grid_name = 'square_grid_output_' // region%name
    CALL setup_square_grid( grid_name, region%mesh%xmin, region%mesh%xmax, region%mesh%ymin, region%mesh%ymax, &
       dx_grid_output, region%output_grid, &
       lambda_M = region%mesh%lambda_M, phi_M = region%mesh%phi_M, beta_stereo = region%mesh%beta_stereo)

    ! Create the main regional output files
    CALL create_main_regional_output_file_mesh( region)
    CALL create_main_regional_output_file_grid( region)

    ! Create the main regional output files for the regions of interest
    CALL setup_ROI_grids_and_output_files( region)

    ! Create the restart files for all the model components
    CALL create_restart_files_ice_model   ( region%mesh, region%ice)
    CALL create_restart_file_thermo       ( region%mesh, region%ice)
    CALL create_restart_file_climate_model( region%mesh, region%climate, region%name)
    CALL create_restart_file_ocean_model  ( region%mesh, region%ocean  , region%name)
    CALL create_restart_file_SMB_model    ( region%mesh, region%SMB    , region%name)
    CALL create_restart_file_BMB_model    ( region%mesh, region%BMB    , region%name)
    CALL create_restart_file_GIA_model    ( region%mesh, region%GIA    , region%name)

    ! Create the scalar regional output file
    CALL create_scalar_regional_output_file( region)

    ! Set output writing time to start of run, so the initial state will be written to output
    IF (C%do_create_netcdf_output) THEN
      region%output_t_next = C%start_time_of_run
      region%output_restart_t_next = C%start_time_of_run
      region%output_grid_t_next = C%start_time_of_run
    ELSE
      region%output_t_next = C%end_time_of_run
      region%output_restart_t_next = C%end_time_of_run
      region%output_grid_t_next = C%end_time_of_run
    END IF

    ! Confirm that the current set of mesh output files match the current model mesh
    ! (set to false whenever a new mesh is created,
    ! and set to true whenever a new set of output files is created)

    region%output_files_match_current_mesh = .TRUE.

    ! ===== Finalisation =====
    ! ========================

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') ' Finished initialising model region ' // colour_string( TRIM( region%long_name),'light blue')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

    ! CALL write_to_main_regional_output_file_mesh( region)
    ! stop ':)'

  END SUBROUTINE initialise_model_region

  SUBROUTINE setup_first_mesh( region)
    ! Set up the first model mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'setup_first_mesh'
    CHARACTER(LEN=256)                                                 :: choice_initial_mesh
    CHARACTER(LEN=256)                                                 :: mesh_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') '   Setting up the first mesh for model region ' // colour_string( region%name,'light blue') // '...'

    ! Get settings from config
    IF     (region%name == 'NAM') THEN
      choice_initial_mesh = C%choice_initial_mesh_NAM
    ELSEIF (region%name == 'EAS') THEN
      choice_initial_mesh = C%choice_initial_mesh_EAS
    ELSEIF (region%name == 'GRL') THEN
      choice_initial_mesh = C%choice_initial_mesh_GRL
    ELSEIF (region%name == 'ANT') THEN
      choice_initial_mesh = C%choice_initial_mesh_ANT
    ELSE
      CALL crash('unknown region name "' // TRIM( region%name) // '"!')
    END IF

    ! Mesh name
    mesh_name = 'model_mesh_' // TRIM( region%name) // '_00001'

    ! Calculate a new mesh based on the initial ice-sheet geometry, or read an existing mesh from a file
    IF     (choice_initial_mesh == 'calc_from_initial_geometry') THEN
      region%time_mesh_was_created = C%start_time_of_run
      CALL setup_first_mesh_from_initial_geometry( mesh_name, region)
    ELSEIF (choice_initial_mesh == 'read_from_file') THEN
      region%time_mesh_was_created = C%start_time_of_run - 2._dp * C%dt_mesh_update_min ! So that we allow a mesh update right at the start
      CALL setup_first_mesh_from_file( mesh_name, region)
    ELSE
      CALL crash('unknown choice_initial_mesh "' // TRIM( choice_initial_mesh) // '"!')
    END IF

    ! Write the mesh creation success message to the terminal
    CALL write_mesh_success( region%mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_first_mesh

  SUBROUTINE setup_first_mesh_from_initial_geometry( mesh_name, region)
    ! Set up the first model mesh based on the initial ice-sheet geometry

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                                  INTENT(IN)    :: mesh_name
    TYPE(type_model_region),                             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'setup_first_mesh_from_initial_geometry'
    REAL(dp)                                                           :: xmin, xmax, ymin, ymax
    REAL(dp)                                                           :: lambda_M, phi_M, beta_stereo

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') '     Creating mesh from initial geometry...'

    ! Determine model domain
    IF     (region%name == 'NAM') THEN
      xmin        = C%xmin_NAM
      xmax        = C%xmax_NAM
      ymin        = C%ymin_NAM
      ymax        = C%ymax_NAM
      lambda_M    = C%lambda_M_NAM
      phi_M       = C%phi_M_NAM
      beta_stereo = C%beta_stereo_NAM
    ELSEIF (region%name == 'EAS') THEN
      xmin        = C%xmin_EAS
      xmax        = C%xmax_EAS
      ymin        = C%ymin_EAS
      ymax        = C%ymax_EAS
      lambda_M    = C%lambda_M_EAS
      phi_M       = C%phi_M_EAS
      beta_stereo = C%beta_stereo_EAS
    ELSEIF (region%name == 'GRL') THEN
      xmin        = C%xmin_GRL
      xmax        = C%xmax_GRL
      ymin        = C%ymin_GRL
      ymax        = C%ymax_GRL
      lambda_M    = C%lambda_M_GRL
      phi_M       = C%phi_M_GRL
      beta_stereo = C%beta_stereo_GRL
    ELSEIF (region%name == 'ANT') THEN
      xmin        = C%xmin_ANT
      xmax        = C%xmax_ANT
      ymin        = C%ymin_ANT
      ymax        = C%ymax_ANT
      lambda_M    = C%lambda_M_ANT
      phi_M       = C%phi_M_ANT
      beta_stereo = C%beta_stereo_ANT
    ELSE
      CALL crash('unknown region name "' // TRIM( region%name) // '"!')
    END IF

    ! Determine if the initial geometry is provided gridded or meshed
    IF (ALLOCATED( region%refgeo_init%grid_raw%x)) THEN
      ! Gridded

      ! Safety
      IF (ALLOCATED( region%refgeo_init%mesh_raw%V)) CALL crash('found boht grid and mesh in region%refgeo_init!')

      ! Create mesh from gridded initial geometry data
      CALL create_mesh_from_gridded_geometry( region%name, mesh_name, &
        region%refgeo_init%grid_raw, &
        region%refgeo_init%Hi_grid_raw, &
        region%refgeo_init%Hb_grid_raw, &
        region%refgeo_init%Hs_grid_raw, &
        region%refgeo_init%SL_grid_raw, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        region%mesh)

    ELSEIF (ALLOCATED( region%refgeo_init%mesh_raw%V)) THEN
      ! Meshed

      ! Safety
      IF (ALLOCATED( region%refgeo_init%grid_raw%x)) CALL crash('found boht grid and mesh in region%refgeo_init!')

      ! Create mesh from meshed initial geometry data
      CALL create_mesh_from_meshed_geometry( region%name, mesh_name, &
        region%refgeo_init%mesh_raw, &
        region%refgeo_init%Hi_mesh_raw, &
        region%refgeo_init%Hb_mesh_raw, &
        region%refgeo_init%Hs_mesh_raw, &
        region%refgeo_init%SL_mesh_raw, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        region%mesh)

    ELSE
      CALL crash('no grid or mesh is found in region%refgeo_init!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_first_mesh_from_initial_geometry

  SUBROUTINE setup_first_mesh_from_file( mesh_name, region)
    ! Set up the first model mesh from an external NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                                  INTENT(IN)    :: mesh_name
    TYPE(type_model_region),                             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'setup_first_mesh'
    REAL(dp)                                                           :: xmin, xmax, ymin, ymax
    REAL(dp)                                                           :: lambda_M, phi_M, beta_stereo
    CHARACTER(LEN=256)                                                 :: filename_initial_mesh
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine model domain
    IF     (region%name == 'NAM') THEN
      xmin        = C%xmin_NAM
      xmax        = C%xmax_NAM
      ymin        = C%ymin_NAM
      ymax        = C%ymax_NAM
      lambda_M    = C%lambda_M_NAM
      phi_M       = C%phi_M_NAM
      beta_stereo = C%beta_stereo_NAM
    ELSEIF (region%name == 'EAS') THEN
      xmin        = C%xmin_EAS
      xmax        = C%xmax_EAS
      ymin        = C%ymin_EAS
      ymax        = C%ymax_EAS
      lambda_M    = C%lambda_M_EAS
      phi_M       = C%phi_M_EAS
      beta_stereo = C%beta_stereo_EAS
    ELSEIF (region%name == 'GRL') THEN
      xmin        = C%xmin_GRL
      xmax        = C%xmax_GRL
      ymin        = C%ymin_GRL
      ymax        = C%ymax_GRL
      lambda_M    = C%lambda_M_GRL
      phi_M       = C%phi_M_GRL
      beta_stereo = C%beta_stereo_GRL
    ELSEIF (region%name == 'ANT') THEN
      xmin        = C%xmin_ANT
      xmax        = C%xmax_ANT
      ymin        = C%ymin_ANT
      ymax        = C%ymax_ANT
      lambda_M    = C%lambda_M_ANT
      phi_M       = C%phi_M_ANT
      beta_stereo = C%beta_stereo_ANT
    ELSE
      CALL crash('unknown region name "' // TRIM( region%name) // '"!')
    END IF

    ! Get filename from config
    IF     (region%name == 'NAM') THEN
      filename_initial_mesh = C%filename_initial_mesh_NAM
    ELSEIF (region%name == 'EAS') THEN
      filename_initial_mesh = C%filename_initial_mesh_EAS
    ELSEIF (region%name == 'GRL') THEN
      filename_initial_mesh = C%filename_initial_mesh_GRL
    ELSEIF (region%name == 'ANT') THEN
      filename_initial_mesh = C%filename_initial_mesh_ANT
    ELSE
      CALL crash('unknown region name "' // TRIM( region%name) // '"!')
    END IF

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') '   Reading mesh from file "' // colour_string( TRIM( filename_initial_mesh),'light blue') // '"...'

    ! Read the mesh from the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename_initial_mesh, ncid)
    CALL setup_mesh_from_file(                  filename_initial_mesh, ncid, region%mesh)
    CALL close_netcdf_file(                                            ncid)

    ! Give the mesh a nice name
    region%mesh%name = mesh_name

    ! Safety - check if the mesh we read from the file matches this region's domain and projection
    IF (region%mesh%xmin        /= xmin       ) CALL crash('expected xmin        = {dp_01}, found {dp_02}', dp_01 = xmin       , dp_02 = region%mesh%xmin       )
    IF (region%mesh%xmax        /= xmax       ) CALL crash('expected xmax        = {dp_01}, found {dp_02}', dp_01 = xmax       , dp_02 = region%mesh%xmax       )
    IF (region%mesh%ymin        /= ymin       ) CALL crash('expected ymin        = {dp_01}, found {dp_02}', dp_01 = ymin       , dp_02 = region%mesh%ymin       )
    IF (region%mesh%ymax        /= ymax       ) CALL crash('expected ymax        = {dp_01}, found {dp_02}', dp_01 = ymax       , dp_02 = region%mesh%ymax       )
    IF (region%mesh%lambda_M    /= lambda_M   ) CALL crash('expected lambda_M    = {dp_01}, found {dp_02}', dp_01 = lambda_M   , dp_02 = region%mesh%lambda_M   )
    IF (region%mesh%phi_M       /= phi_M      ) CALL crash('expected phi_M       = {dp_01}, found {dp_02}', dp_01 = phi_M      , dp_02 = region%mesh%phi_M      )
    IF (region%mesh%beta_stereo /= beta_stereo) CALL crash('expected beta_stereo = {dp_01}, found {dp_02}', dp_01 = beta_stereo, dp_02 = region%mesh%beta_stereo)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_first_mesh_from_file

  SUBROUTINE setup_ROI_grids_and_output_files( region)
    ! Set up the (high-resolution) square grids and output files
    ! for this model region's Regions of Interest

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'setup_ROI_grids_and_output_files'
    REAL(dp)                                      :: dx_ROI
    CHARACTER(LEN=256)                            :: all_names_ROI, name_ROI
    INTEGER                                       :: i
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: poly_ROI
    REAL(dp)                                      :: xmin, xmax, ymin, ymax, xmid, ymid, dx, dy
    CHARACTER(LEN=256)                            :: grid_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no regions of interest are specified, do nothing
    IF (C%choice_regions_of_interest == '') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Determine the resolution for this region's region-of-interest output files
    SELECT CASE (region%name)
      CASE ('NAM')
        dx_ROI = C%dx_output_grid_ROI_NAM
      CASE ('EAS')
        dx_ROI = C%dx_output_grid_ROI_EAS
      CASE ('GRL')
        dx_ROI = C%dx_output_grid_ROI_GRL
      CASE ('ANT')
        dx_ROI = C%dx_output_grid_ROI_ANT
      CASE DEFAULT
        CALL crash('unknown region name "' // TRIM( region%name) // '"!')
    END SELECT

    ! Initialise
    region%nROI = 0

    ! Go over all listed regions of interest
    all_names_ROI = C%choice_regions_of_interest

    DO WHILE (.TRUE.)

      ! == Parse list of input ROIs
      ! ===========================

      ! Get the first region of interest from the list
      i = INDEX( all_names_ROI, '||')
      IF (i == 0) THEN
        ! There is only one left in the list
        name_ROI = TRIM( all_names_ROI)
        all_names_ROI = ''
      ELSE
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_TRIM( all_names_ROI))
      END IF

      ! == Check validity of requested ROIs
      ! ===================================

      ! Check if current region is indeed defined in the model
      SELECT CASE (name_ROI)
        CASE ('')
         ! No region requested: don't need to do anything
         EXIT
        CASE ('PineIsland','Thwaites','Amery','RiiserLarsen','SipleCoast','LarsenC','TransMounts','DotsonCrosson', & ! Antarctica
              'Narsarsuaq','Nuuk','Jakobshavn','NGIS','Qaanaaq', &                                                   ! Greenland
              'Patagonia', &                                                                                         ! Patagonia
              'Tijn_test_ISMIP_HOM_A','CalvMIP_quarter')                                                             ! Idealised
          ! List of known regions of interest: these pass the test
        CASE DEFAULT
          ! Region not found
          CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
      END SELECT

      ! == Calculate ROIs
      ! =================

      ! Calculate the polygon describing the specified region of interest
      SELECT CASE (region%name)
        CASE ('NAM')
          ! North america

          SELECT CASE (name_ROI)
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE ('EAS')
          ! Eurasia

          SELECT CASE (name_ROI)
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE ('GRL')
          ! Greenland

          SELECT CASE (name_ROI)
            CASE ('Narsarsuaq')
              CALL calc_polygon_Narsarsuaq( poly_ROI)
            CASE ('Nuuk')
              CALL calc_polygon_Nuuk( poly_ROI)
            CASE ('Jakobshavn')
              CALL calc_polygon_Jakobshavn( poly_ROI)
            CASE ('NGIS')
              CALL calc_polygon_NGIS( poly_ROI)
            CASE ('Qaanaaq')
              CALL calc_polygon_Qaanaaq( poly_ROI)
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE ('ANT')

          SELECT CASE (name_ROI)
            CASE ('PineIsland')
              CALL calc_polygon_Pine_Island_Glacier( poly_ROI)
            CASE ('Thwaites')
              CALL calc_polygon_Thwaites_Glacier( poly_ROI)
            CASE ('Amery')
              CALL calc_polygon_Amery_ice_shelf( poly_ROI)
            CASE ('RiiserLarsen')
              CALL calc_polygon_Riiser_Larsen_ice_shelf( poly_ROI)
            CASE ('SipleCoast')
              CALL calc_polygon_Siple_Coast( poly_ROI)
            CASE ('LarsenC')
              CALL calc_polygon_Larsen_ice_shelf( poly_ROI)
            CASE ('TransMounts')
              CALL calc_polygon_Transantarctic_Mountains( poly_ROI)
            CASE ('DotsonCrosson')
              CALL calc_polygon_DotsonCrosson_ice_shelf( poly_ROI)
            CASE ('Patagonia')
              CALL calc_polygon_Patagonia( poly_ROI)
            CASE ('Tijn_test_ISMIP_HOM_A')
              CALL calc_polygon_Tijn_test_ISMIP_HOM_A( poly_ROI)
            CASE ('CalvMIP_quarter')
              ! Not interesting; skip
              CYCLE
            CASE DEFAULT
              ! Requested area not in this model domain; skip
              CYCLE
          END SELECT

        CASE DEFAULT
          CALL crash('unknown region name "' // region%name // '"!')
      END SELECT

      ! Add to this region's regions of interest
      region%nROI = region%nROI + 1

      ! Determine the bounds of this region of interest's square domain
      xmin = MINVAL( poly_ROI( :,1))
      xmax = MAXVAL( poly_ROI( :,1))
      ymin = MINVAL( poly_ROI( :,2))
      ymax = MAXVAL( poly_ROI( :,2))

      ! Make it 20% in all directions to fully cover the hi-res area
      xmid = (xmin + xmax) / 2._dp
      ymid = (ymin + ymax) / 2._dp
      dx   = xmax - xmid;
      dy   = ymax - ymid
      xmin = xmid - dx * 1.2_dp
      xmax = xmid + dx * 1.2_dp
      ymin = ymid - dy * 1.2_dp
      ymax = ymid + dy * 1.2_dp

      ! Set up a square grid for this region of interest
      grid_name = 'square_grid_output_' // region%name // '_ROI_' // TRIM( name_ROI)
      CALL setup_square_grid( grid_name, xmin, xmax, ymin, ymax, dx_ROI, region%output_grids_ROI( region%nROI), &
        region%mesh%lambda_M, region%mesh%phi_M, region%mesh%beta_stereo)

      ! Create an output file for this region of interest
      region%output_filenames_grid_ROI( region%nROI) = TRIM( C%output_dir) // 'main_output_' // region%name // '_grid_ROI_' // TRIM( name_ROI) // '.nc'
      CALL create_main_regional_output_file_grid_ROI( region, region%output_grids_ROI( region%nROI), region%output_filenames_grid_ROI( region%nROI))

      ! Clean up after yourself
      DEALLOCATE( poly_ROI)

      ! If no names are left, we are finished
      IF (all_names_ROI == '') EXIT

    END DO ! DO WHILE (.TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_ROI_grids_and_output_files

! ===== Mesh update =====
! =======================

  SUBROUTINE update_mesh( region)
    ! Update the model mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'update_mesh'
    REAL(dp)                                                           :: xmin, xmax, ymin, ymax
    REAL(dp)                                                           :: lambda_M, phi_M, beta_stereo
    INTEGER                                                            :: n_old, stat, n_new
    CHARACTER(LEN=5)                                                   :: n_new_str
    CHARACTER(LEN=256)                                                 :: new_mesh_name
    TYPE(type_mesh)                                                    :: mesh_new
    REAL(dp)                                                           :: tstart, tstop
    CHARACTER(LEN=256)                                                 :: str

    ! Add routine to path
    CALL init_routine( routine_name)

    ! We want to know how long this takes, just to show off...
    tstart = MPI_WTIME()

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') ''
    IF (par%master) WRITE(0,'(A)') '   Creating a new mesh for model region ' // colour_string( region%name,'light blue') // '...'

    ! Determine model domain
    IF     (region%name == 'NAM') THEN
      xmin        = C%xmin_NAM
      xmax        = C%xmax_NAM
      ymin        = C%ymin_NAM
      ymax        = C%ymax_NAM
      lambda_M    = C%lambda_M_NAM
      phi_M       = C%phi_M_NAM
      beta_stereo = C%beta_stereo_NAM
    ELSEIF (region%name == 'EAS') THEN
      xmin        = C%xmin_EAS
      xmax        = C%xmax_EAS
      ymin        = C%ymin_EAS
      ymax        = C%ymax_EAS
      lambda_M    = C%lambda_M_EAS
      phi_M       = C%phi_M_EAS
      beta_stereo = C%beta_stereo_EAS
    ELSEIF (region%name == 'GRL') THEN
      xmin        = C%xmin_GRL
      xmax        = C%xmax_GRL
      ymin        = C%ymin_GRL
      ymax        = C%ymax_GRL
      lambda_M    = C%lambda_M_GRL
      phi_M       = C%phi_M_GRL
      beta_stereo = C%beta_stereo_GRL
    ELSEIF (region%name == 'ANT') THEN
      xmin        = C%xmin_ANT
      xmax        = C%xmax_ANT
      ymin        = C%ymin_ANT
      ymax        = C%ymax_ANT
      lambda_M    = C%lambda_M_ANT
      phi_M       = C%phi_M_ANT
      beta_stereo = C%beta_stereo_ANT
    ELSE
      CALL crash('unknown region name "' // TRIM( region%name) // '"!')
    END IF

    ! Create numbered name for the new mesh
    CALL str2int( region%mesh%name( 16:20), n_old, stat)
    IF (stat /= 0) CALL crash('couldnt read number of old mesh!')
    n_new = n_old + 1
    CALL int2str( n_new, n_new_str)
    new_mesh_name = 'model_mesh_' // TRIM( region%name) // '_' // n_new_str

    ! Create a mesh from the modelled ice geometry
    CALL create_mesh_from_meshed_geometry( region%name, new_mesh_name, &
      region%mesh, &
      region%ice%Hi, &
      region%ice%Hb, &
      region%ice%Hs, &
      region%ice%SL, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
      mesh_new)

    ! Write the mesh creation success message to the terminal
    CALL write_mesh_success( mesh_new)

    ! Remap the reference geometries to the new mesh (deallocation of old data happens in there)
    CALL initialise_reference_geometries_on_model_mesh( region%name, mesh_new, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq)

    ! Remap all the model data from the old mesh to the new mesh
    CALL remap_ice_dynamics_model( region%mesh, mesh_new, region%ice, region%refgeo_PD, region%SMB, region%BMB, region%LMB, region%AMB, region%GIA, region%time, region%name)
    CALL remap_climate_model(      region%mesh, mesh_new,             region%climate, region%name)
    CALL remap_ocean_model(        region%mesh, mesh_new,             region%ocean  , region%name)
    CALL remap_SMB_model(          region%mesh, mesh_new,             region%SMB    , region%name)
    CALL remap_BMB_model(          region%mesh, mesh_new, region%ice, region%BMB    , region%name)
    CALL remap_LMB_model(          region%mesh, mesh_new,             region%LMB    , region%name)
    CALL remap_AMB_model(          region%mesh, mesh_new,             region%AMB                 )
    CALL remap_GIA_model(          region%mesh, mesh_new,             region%GIA                 )

    ! Set all model component timers so that they will all be run right after the mesh update
    region%ice%t_Hi_next  = region%time
    region%ice%t_Ti_next  = region%time
    region%climate%t_next = region%time
    region%ocean%t_next   = region%time
    region%SMB%t_next     = region%time
    region%BMB%t_next     = region%time
    region%LMB%t_next     = region%time
    region%GIA%t_next     = region%time

    ! Throw away the mapping operators involving the old mesh
    CALL clear_all_maps_involving_this_mesh( region%mesh)

    ! Deallocate the old mesh, bind the region%mesh pointers to the new mesh.
    CALL deallocate_mesh( region%mesh)
    region%mesh = mesh_new
    region%time_mesh_was_created = region%time

    ! Confirm that the current set of mesh output files do not match the current model mesh
    ! (set to false whenever a new mesh is created,
    ! and set to true whenever a new set of output files is created)

    region%output_files_match_current_mesh = .FALSE.

    ! We want to know how long this takes, just to show off...
    tstop = MPI_WTIME()

    ! Print to screen
    str = '   Finished the mesh update in {dp_01} seconds'
    CALL insert_val_into_string_dp( str, '{dp_01}', tstop-tstart)
    IF (par%master) WRITE(0,'(A)') str

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_mesh

  SUBROUTINE calc_mesh_fitness_coefficient( mesh, ice, f)
    ! Calculate the fitness coefficient of the current mesh, which is defined as the
    ! fraction of the grounding line that is still covered by a good-resolution mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                                INTENT(IN)    :: ice
    REAL(dp),                                            INTENT(OUT)   :: f

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'calc_mesh_fitness_coefficient'
    INTEGER                                                            :: vi
    INTEGER                                                            :: nV_grounding_line_tot, nV_calving_front_tot, nV_ice_front_tot, nV_coastline_tot
    INTEGER                                                            :: nV_grounding_line_bad, nV_calving_front_bad, nV_ice_front_bad, nV_coastline_bad
    REAL(dp)                                                           :: f_grounding_line, f_calving_front, f_ice_front, f_coastline

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    nV_grounding_line_tot     = 0
    nV_calving_front_tot      = 0
    nV_ice_front_tot          = 0
    nV_coastline_tot          = 0

    nV_grounding_line_bad     = 0
    nV_calving_front_bad      = 0
    nV_ice_front_bad          = 0
    nV_coastline_bad          = 0

    ! Check resolution criteria

    DO vi = mesh%vi1, mesh%vi2

      ! Check resolution criteria for the different lines

      ! Grounding line
      IF (ice%mask_gl_gr( vi) .OR. ice%mask_gl_fl( vi)) THEN
        nV_grounding_line_tot = nV_grounding_line_tot + 1
        IF (mesh%R( vi) > C%maximum_resolution_grounding_line * C%mesh_resolution_tolerance) THEN
          nV_grounding_line_bad = nV_grounding_line_bad + 1
        END IF
      END IF ! IF (ice%mask_gl_gr( vi) .OR. ice%mask_gl_fl( vi)) THEN

      ! Calving front
      IF (ice%mask_cf_gr( vi) .OR. ice%mask_cf_fl( vi)) THEN
        nV_calving_front_tot = nV_calving_front_tot + 1
        IF (mesh%R( vi) > C%maximum_resolution_calving_front * C%mesh_resolution_tolerance) THEN
          nV_calving_front_bad = nV_calving_front_bad + 1
        END IF
      END IF ! IF (ice%mask_gl_gr( vi) .OR. ice%mask_gl_fl( vi)) THEN

      ! Ice front
      IF (ice%mask_margin( vi)) THEN
        nV_ice_front_tot = nV_ice_front_tot + 1
        IF (mesh%R( vi) > C%maximum_resolution_ice_front * C%mesh_resolution_tolerance) THEN
          nV_ice_front_bad = nV_ice_front_bad + 1
        END IF
      END IF ! IF (ice%mask_gl_gr( vi) .OR. ice%mask_gl_fl( vi)) THEN

      ! Coastline
      IF (ice%mask_coastline( vi)) THEN
        nV_coastline_tot = nV_coastline_tot + 1
        IF (mesh%R( vi) > C%maximum_resolution_coastline * C%mesh_resolution_tolerance) THEN
          nV_coastline_bad = nV_coastline_bad + 1
        END IF
      END IF ! IF (ice%mask_gl_gr( vi) .OR. ice%mask_gl_fl( vi)) THEN

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Gather data from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_grounding_line_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_calving_front_tot , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_ice_front_tot     , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_coastline_tot     , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_grounding_line_bad, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_calving_front_bad , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_ice_front_bad     , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nV_coastline_bad     , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate mesh fitness

    IF (nV_grounding_line_tot > 0) THEN
      f_grounding_line = 1._dp - REAL( nV_grounding_line_bad,dp) / REAL( nV_grounding_line_tot,dp)
    ELSE
      f_grounding_line = 1._dp
    END IF

    IF (nV_calving_front_tot > 0) THEN
      f_calving_front = 1._dp - REAL( nV_calving_front_bad,dp) / REAL( nV_calving_front_tot,dp)
    ELSE
      f_calving_front = 1._dp
    END IF

    IF (nV_ice_front_tot > 0) THEN
      f_ice_front = 1._dp - REAL( nV_ice_front_bad,dp) / REAL( nV_ice_front_tot,dp)
    ELSE
      f_ice_front = 1._dp
    END IF

    IF (nV_coastline_tot > 0) THEN
      f_coastline = 1._dp - REAL( nV_coastline_bad,dp) / REAL( nV_coastline_tot,dp)
    ELSE
      f_coastline = 1._dp
    END IF

    f = MINVAL( [f_grounding_Line, f_calving_front, f_ice_front, f_coastline])

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_mesh_fitness_coefficient

! ===== Extras =====
! ==================

  SUBROUTINE time_display( region, t_end, dt_av, ndt_av)
    ! Little time display for the screen

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_carriage_return

    IMPLICIT NONE

    ! Input/Ouput variables
    TYPE(type_model_region),                             INTENT(IN)    :: region
    REAL(dp),                                            INTENT(IN)    :: t_end
    REAL(dp),                                            INTENT(IN)    :: dt_av
    INTEGER,                                             INTENT(IN)    :: ndt_av

    ! Local variables
    REAL(dp)                                                           :: dt_ice
    CHARACTER(LEN=9)                                                   :: r_time, r_step, r_adv

    dt_ice = region%ice%t_Hi_next - region%ice%t_Hi_prev

    IF (region%time < t_end) THEN
      r_adv = "no"
      WRITE( r_time,"(F8.3)") MIN( region%time,t_end) / 1000._dp
      WRITE( r_step,"(F6.3)") MAX( dt_ice, C%dt_ice_min)
      WRITE( *     ,"(A)", ADVANCE = TRIM( r_adv)) c_carriage_return // &
              "   t = " // TRIM( r_time) // " kyr - dt = " // TRIM( r_step) // " yr"
    ELSE
      r_adv = "yes"
      WRITE( r_time,"(F8.3)") MIN( region%time,t_end) / 1000._dp
      WRITE( r_step,"(F6.3)") dt_av / REAL( ndt_av,dp)
      WRITE( *     ,"(A)", ADVANCE = TRIM( r_adv)) c_carriage_return // &
            "   t = " // TRIM( r_time) // " kyr - dt_av = " // TRIM( r_step) // " yr"
    END IF
    IF (region%time == region%output_t_next .OR. &
        region%time == region%output_restart_t_next .OR. &
        region%time == region%output_grid_t_next) THEN
      r_adv = "yes"
      WRITE( *,"(A)", ADVANCE = TRIM( r_adv)) c_carriage_return
    END IF

  END SUBROUTINE time_display

  SUBROUTINE apply_regional_corrections( region)
    ! Apply some regional modifications to a few fields which
    ! require other fields that are usually initialised later
    ! than the fields of interest.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region), INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER          :: routine_name = 'apply_regional_corrections'
    INTEGER                                :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,'(A)') '   Implementing regional corrections...'

    ! === SMB over ice-free land ===
    ! ==============================

    ! If so desired, remove positive SMB over ice-free land
    IF (C%do_SMB_removal_icefree_land) THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        IF (region%ice%mask_icefree_land( vi)) region%SMB%SMB( vi) = MIN( region%SMB%SMB( vi), 0._dp)
      END DO
    END IF

    ! == Target dHi_dt
    ! ================

    ! If so desired, limit target dH/dt to available SMB. If the target dH/dt is positive,
    ! it will remove ice during the spinup; we don't want it removing more ice than the SMB
    ! is providing, to avoid friction inversions going crazy to keep the ice in place.
    IF (C%do_target_dHi_dt .AND. C%do_limit_target_dHi_dt_to_SMB) THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        IF (region%ice%dHi_dt_target( vi) > 0._dp) THEN
          region%ice%dHi_dt_target( vi) = MAX( 0._dp, MIN( region%ice%dHi_dt_target( vi), region%SMB%SMB( vi)))
        END IF
      END DO
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_regional_corrections

END MODULE UFEMISM_main_model
