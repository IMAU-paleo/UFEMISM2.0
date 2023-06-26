MODULE UFEMISM_main_model

  ! The main regional model

! ===== Preamble =====
! ====================

  ! USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: happy, warning, crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE region_types                                           , ONLY: type_model_region
  USE reference_geometries                                   , ONLY: initialise_reference_geometries_raw, initialise_reference_geometries_on_model_mesh
  USE ice_model_main                                         , ONLY: initialise_ice_dynamics_model  , run_ice_dynamics_model
  USE thermodynamics_main                                    , ONLY: initialise_thermodynamics_model, run_thermodynamics_model
  USE climate_main                                           , ONLY: initialise_climate_model       , run_climate_model
  USE ocean_main                                             , ONLY: initialise_ocean_model         , run_ocean_model
  USE SMB_main                                               , ONLY: initialise_SMB_model           , run_SMB_model
  USE BMB_main                                               , ONLY: initialise_BMB_model           , run_BMB_model
  USE GIA_main                                               , ONLY: initialise_GIA_model           , run_GIA_model
  USE netcdf_basic                                           , ONLY: open_existing_netcdf_file_for_reading, close_netcdf_file
  USE netcdf_input                                           , ONLY: setup_mesh_from_file
  USE mesh_creation                                          , ONLY: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success
  USE mesh_operators                                         , ONLY: calc_all_matrix_operators_mesh
  USE grid_basic                                             , ONLY: setup_square_grid
  USE main_regional_output                                   , ONLY: create_main_regional_output_files, write_to_main_regional_output_files

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

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

    ! Add routine to path
    routine_name = 'run_model('  //  region%name  //  ')'
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE (0,'(A,A,A,A,A,F9.3,A,F9.3,A)') &
      '  Running model region ', colour_string( region%name, 'light blue'), ' (', colour_string( TRIM( region%long_name), 'light blue'), &
      ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'

    ! Initialise average ice-dynamical time step
    ndt_av = 0
    dt_av  = 0._dp

    ! The main UFEMISM time loop
    main_time_loop: DO WHILE (region%time <= t_end)

      ! Run the ice dynamics model to calculate ice geometry at the desired time, and update
      ! velocities, thinning rates, and predicted geometry if necessary
      CALL run_ice_dynamics_model( region)

      ! Keep track of the average ice-dynamical time step and print it to the terminal
      ndt_av = ndt_av + 1
      dt_av  = dt_av  + (region%ice%t_Hi_next - region%ice%t_Hi_prev)
      IF (par%master .AND. C%do_time_display) CALL time_display( region, t_end, dt_av, ndt_av)

      ! Calculate ice temperature at the desired time, and update
      ! predicted temperature if necessary
      CALL run_thermodynamics_model( region)

      ! Calculate the climate
      CALL run_climate_model( region%mesh, region%ice, region%climate, region%name, region%time)

      ! Calculate the ocean
      CALL run_ocean_model( region%mesh, region%ice, region%ocean, region%name, region%time)

      ! Calculate the surface mass balance
      CALL run_SMB_model( region%mesh, region%ice, region%climate, region%SMB, region%name, region%time)

      ! Calculate the basal mass balance
      CALL run_BMB_model( region%mesh, region%ice, region%ocean, region%BMB, region%name, region%time)

      ! Calculate bedrock deformation at the desired time, and update
      ! predicted deformation if necessary
      CALL run_GIA_model( region)

      ! Write to the main regional output NetCDF file
      CALL write_to_main_regional_output_files( region)

      ! Advance this region's time to the time of the next "action"
      CALL advance_region_time_to_time_of_next_action( region, t_end)

      ! If we've reached the end of this coupling interval, stop
      IF (region%time == t_end) EXIT main_time_loop

    END DO main_time_loop ! DO WHILE (region%time <= t_end)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_region

  ! == Time stepping
  ! ================

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

  ! == Find time of next action
  ! ===========================

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

    ! GIA
    time_of_next_action = MIN( time_of_next_action, region%GIA%t_next)

    ! Output
    time_of_next_action = MIN( time_of_next_action, region%output_t_next)

  ! == Advance region time
  ! ======================

    region%time = time_of_next_action

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE advance_region_time_to_time_of_next_action

  ! == Model initialisation
  ! =======================

  SUBROUTINE initialise_model_region( region, region_name)
    ! Initialise this model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(OUT)   :: region
    CHARACTER(LEN=3),                                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_model_region'
    CHARACTER(LEN=256)                                                 :: grid_name
    REAL(dp)                                                           :: dx_grid_output

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! ===== Ice dynamics =====
    ! ========================

    CALL initialise_ice_dynamics_model( region%mesh, region%ice, region%refgeo_init, region%refgeo_PD, region%scalars, region%name)

    ! ===== Thermodynamics =====
    ! ==========================

    CALL initialise_thermodynamics_model( region)

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

    CALL initialise_BMB_model( region%mesh, region%BMB, region%name)

    ! ===== Glacial isostatic adjustment =====
    ! ========================================

    CALL initialise_GIA_model( region%mesh, region%GIA, region%name)

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

    ! Create the mesh and grid output files
    CALL create_main_regional_output_files( region)

    ! Set output writing time to stat of run, so the initial state will be written to output
    IF (C%do_create_netcdf_output) THEN
      region%output_t_next = C%start_time_of_run
    ELSE
      region%output_t_next = C%end_time_of_run
    END IF

    ! ===== Finalisation =====
    ! ========================

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') ' Finished initialising model region ' // colour_string( TRIM( region%long_name),'light blue')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

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
    IF (par%master) WRITE(0,'(A)') '  Setting up the first mesh for model region ' // colour_string( region%name,'light blue') // '...'

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
      CALL setup_first_mesh_from_initial_geometry( mesh_name, region)
    ELSEIF (choice_initial_mesh == 'read_from_file') THEN
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
    IF (par%master) WRITE(0,'(A)') '   Creating mesh from initial geometry...'

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

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( region%mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_first_mesh_from_file

  ! == Extras
  ! =========

  SUBROUTINE time_display( region, t_end, dt_av, ndt_av)
    ! Little time display for the screen

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_backspace

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
      WRITE( *     ,"(A)", ADVANCE = TRIM( r_adv)) REPEAT( c_backspace,999) // &
              "   t = " // TRIM( r_time) // " kyr - dt = " // TRIM( r_step) // " yr"
    ELSE
      r_adv = "yes"
      WRITE( r_time,"(F8.3)") MIN( region%time,t_end) / 1000._dp
      WRITE( r_step,"(F6.3)") dt_av / REAL( ndt_av,dp)
      WRITE( *     ,"(A)", ADVANCE = TRIM( r_adv)) REPEAT( c_backspace,999) // &
            "   t = " // TRIM( r_time) // " kyr - dt_av = " // TRIM( r_step) // " yr"
    END IF
    IF (region%time == region%output_t_next) THEN
      r_adv = "no"
      WRITE( *,"(A)", ADVANCE = TRIM( r_adv)) REPEAT( c_backspace,999)
    END IF

  END SUBROUTINE time_display

END MODULE UFEMISM_main_model