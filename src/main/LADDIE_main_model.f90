module LADDIE_main_model

  ! The main regional model

! ===== Preamble =====
! ====================

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_WTIME
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: happy, warning, crash, init_routine, finalise_routine, colour_string, str2int, int2str, &
                                                                     insert_val_into_string_dp
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE region_types                                           , ONLY: type_model_region
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE mesh_types                                             , ONLY: type_mesh
  use grid_types, only: type_grid
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE global_forcing_types                                   , ONLY: type_global_forcing
  use reference_geometries_main, only: initialise_reference_geometry_raw_from_file, &
    remap_reference_geometry_to_mesh
  use ice_dynamics_main, only: initialise_ice_dynamics_model, run_ice_dynamics_model, remap_ice_dynamics_model, &
    create_restart_files_ice_model, write_to_restart_files_ice_model, apply_geometry_relaxation
  use basal_hydrology_main, only: run_basal_hydrology_model
  use bed_roughness_main, only: initialise_bed_roughness_model
  USE thermodynamics_main                                    , ONLY: initialise_thermodynamics_model, run_thermodynamics_model, &
                                                                     create_restart_file_thermo, write_to_restart_file_thermo
  USE global_forcings_main                                   , ONLY: initialise_global_forcings, update_sealevel_in_model, update_sealevel_at_model_time
  USE ocean_main                                             , ONLY: initialise_ocean_model, run_ocean_model, remap_ocean_model, &
                                                                     create_restart_file_ocean_model, write_to_restart_file_ocean_model
  USE BMB_main                                               , ONLY: initialise_BMB_model, run_BMB_model, remap_BMB_model, &
                                                                     create_restart_file_BMB_model, write_to_restart_file_BMB_model
  use netcdf_io_main
  USE mesh_creation_main                                     , ONLY: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success
  USE grid_basic                                             , ONLY: setup_square_grid
  USE mesh_output_files, only: create_main_regional_output_file_mesh, write_to_main_regional_output_file_mesh
  use grid_output_files, only: create_main_regional_output_file_grid, write_to_main_regional_output_file_grid, &
    create_main_regional_output_file_grid_ROI, write_to_main_regional_output_file_grid_ROI
  use scalar_output_files, only: create_scalar_regional_output_file, buffer_scalar_output, write_to_scalar_regional_output_file
  use mesh_ROI_polygons
  use plane_geometry, only: longest_triangle_leg
  use apply_maps, only: clear_all_maps_involving_this_mesh
  USE mesh_memory                                            , ONLY: deallocate_mesh
  use ice_mass_and_fluxes, only: calc_ice_mass_and_fluxes
  use tracer_tracking_model_main, only: initialise_tracer_tracking_model, run_tracer_tracking_model, &
    remap_tracer_tracking_model
  use transects_main, only: initialise_transects, write_to_transect_netcdf_output_files
  use UFEMISM_main_model, only: write_to_regional_output_files, setup_ROI_grids_and_output_files
  use laddie_output, only: create_laddie_output_fields_file 
  use laddie_main, only: initialise_laddie_model

  IMPLICIT NONE

CONTAINS

! ===== Main routine =====
! ========================

  SUBROUTINE run_model_region( region, t_end, forcing)
    ! Integrate this model region forward in time until t_end

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(INOUT) :: region
    REAL(dp)                                           , INTENT(IN)    :: t_end    ! [yr]
    TYPE(type_global_forcing)                          , INTENT(IN)    :: forcing

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: routine_name
    INTEGER                                                            :: ndt_av
    REAL(dp)                                                           :: dt_av
    REAL(dp)                                                           :: mesh_fitness_coefficient
    TYPE(type_global_forcing)                                          :: regional_forcing

    ! Add routine to path
    routine_name = 'run_model('  //  region%name  //  ')'
    CALL init_routine( routine_name)

    IF (par%primary) WRITE(0,*) ''
    IF (par%primary) WRITE (0,'(A,A,A,A,A,F9.3,A,F9.3,A)') &
      ' Running model region ', colour_string( region%name, 'light blue'), ' (', colour_string( TRIM( region%long_name), 'light blue'), &
      ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'

    ! Initialise average ice-dynamical time step
    ndt_av = 0
    dt_av  = 0._dp

    ! we create a copy of the forcing type so the different regional models can run asynchronously
    regional_forcing = forcing

    ! The main UFEMISM time loop
    main_time_loop: DO WHILE (region%time <= t_end)

      ! Run the ice dynamics model to calculate ice geometry at the desired time, and update
      ! velocities, thinning rates, and predicted geometry if necessary
      CALL run_ice_dynamics_model( region)

      ! Calculate the ocean
      CALL run_ocean_model( region%mesh, region%ice, region%ocean, region%name, region%time)

      ! Calculate the basal mass balance
      CALL run_BMB_model( region%mesh, region%ice, region%ocean, region%refgeo_PD, region%SMB, region%BMB, region%name, region%time, is_initial=.FALSE.)

      ! Calculate ice-sheet integrated values (total volume, area, etc.)
      CALL calc_ice_mass_and_fluxes( region%mesh, region%ice, region%SMB, region%BMB, region%LMB, region%refgeo_PD, region%scalars)

      ! Write to the main regional output NetCDF file
      CALL write_to_regional_output_files( region)

      ! If we've reached the end of this coupling interval, stop
      IF (region%time == t_end) EXIT main_time_loop

    END DO main_time_loop ! DO WHILE (region%time <= t_end)

    IF (region%time == C%end_time_of_run) THEN
      ! Give all processes time to catch up
      CALL sync
      ! Congrats, you've made it
      IF (par%primary) WRITE(0,'(A)') ' Finalising regional simulation...'
      ! Write the final model state to output
      CALL write_to_regional_output_files( region)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_region

! ===== Model initialisation =====
! ================================

  SUBROUTINE initialise_model_region( region, region_name, laddie, refgeo, mesh, forcing, start_time_of_run)
    ! Initialise this model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(OUT)   :: region
    CHARACTER(LEN=3),                                    INTENT(IN)    :: region_name
    type(type_laddie_model),                             intent(out)   :: laddie
    type(type_reference_geometry),                       intent(out)   :: refgeo
    type(type_mesh),                                     intent(out)   :: mesh
    TYPE(type_global_forcing)                          , INTENT(IN)    :: forcing
    REAL(dp)                                           , INTENT(IN)    :: start_time_of_run

    ! Local variables:
    character(len=256), parameter                                      :: routine_name = 'initialise_model_region'
    character(len=256)                                                 :: grid_name
    character(len=256)                                                 :: filename_refgeo
    real(dp)                                                           :: timeframe_refgeo
    real(dp)                                                           :: dx_grid_output
    real(dp)                                                           :: time
    type(type_ocean_model)                                             :: ocean
    type(type_ice_model)                                               :: ice
    type(type_grid)                                                    :: output_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Set time
    time = C%start_time_of_run

    ! Print to screen
    if (par%primary) write(0,'(A)') ''
    if (par%primary) write(0,'(A)') ' Initialising model ...'

    ! ===== Raw geometry =====
    ! ========================

    filename_refgeo = C%filename_refgeo_init_ANT
    timeframe_refgeo = C%timeframe_refgeo_init_ANT

    ! Clean up memory if necessary
    if (allocated( refgeo%Hi_grid_raw)) deallocate( refgeo%Hi_grid_raw)
    if (allocated( refgeo%Hb_grid_raw)) deallocate( refgeo%Hb_grid_raw)
    if (allocated( refgeo%Hs_grid_raw)) deallocate( refgeo%Hs_grid_raw)
    if (allocated( refgeo%SL_grid_raw)) deallocate( refgeo%SL_grid_raw)

    ! Initialise the reference geometry on their raw input grids
    call initialise_reference_geometry_raw_from_file( region_name, 'refgeo', refgeo, filename_refgeo, timeframe_refgeo)

    ! ===== Geometry on mesh =====
    ! ============================

    ! Set up the model mesh
    call setup_mesh( region_name, mesh, refgeo)

    ! Remap reference geometries from their raw input grids to the model mesh

    ! deallocate existing memory if needed
    if (allocated( refgeo%Hi)) deallocate( refgeo%Hi)
    if (allocated( refgeo%Hb)) deallocate( refgeo%Hb)
    if (allocated( refgeo%Hs)) deallocate( refgeo%Hs)
    if (allocated( refgeo%SL)) deallocate( refgeo%SL)

    ! allocate memory for reference ice geometry on the model mesh
    allocate( refgeo%Hi( mesh%vi1:mesh%vi2))
    allocate( refgeo%Hb( mesh%vi1:mesh%vi2))
    allocate( refgeo%Hs( mesh%vi1:mesh%vi2))
    allocate( refgeo%SL( mesh%vi1:mesh%vi2))

    call remap_reference_geometry_to_mesh( mesh, refgeo)

    if (par%primary) write(0,'(A)') ' Got geometry ...'

    ! ===== Ice dynamics =====
    ! ========================

    ! CALL initialise_ice_dynamics_model( region%mesh, region%ice, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%GIA, region%name, regional_forcing, time)

    ! ===== Ocean =====
    ! =================

    call initialise_ocean_model( mesh, ice, ocean, region_name, time)

    ! ===== LADDIE =====
    ! ==================

    call initialise_laddie_model( mesh, laddie, ocean, ice, region_name)

    ! ===== Integrated scalars =====
    ! ==============================

    ! Calculate ice-sheet integrated values (total volume, area, etc.)
    ! CALL calc_ice_mass_and_fluxes( mesh, ice, region%SMB, region%BMB, region%LMB, region%refgeo_PD, region%scalars)

    ! ===== Regional output =====
    ! ===========================

    dx_grid_output = C%dx_output_grid_ANT

    ! Create the square output grid
    grid_name = 'square_grid_output_' // region_name
    call setup_square_grid( grid_name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, &
       dx_grid_output, output_grid, &
       lambda_M = region%mesh%lambda_M, phi_M = region%mesh%phi_M, beta_stereo = region%mesh%beta_stereo)

    ! Create the main regional output files
    call create_laddie_output_fields_file( mesh, laddie, region_name)
    !CALL create_main_regional_output_file_grid( region)

    ! Create the main regional output files for the regions of interest
    CALL setup_ROI_grids_and_output_files( region)

    ! Create the restart files for all the model components
    !CALL create_restart_file_BMB_model    ( region%mesh, region%BMB    , region%name)

    ! Initialise the output transects
    call initialise_transects( region)

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

    region%output_files_match_current_mesh             = .TRUE.
    region%BMB%laddie%output_fields_file_matches_current_mesh = .TRUE.

    ! ===== Finalisation =====
    ! ========================

    ! Print to screen
    IF (par%primary) WRITE(0,'(A)') ' Finished initialising model region ' // colour_string( TRIM( region%long_name),'light blue')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_model_region

  subroutine setup_mesh( region_name, mesh, refgeo)
    ! Set up the model mesh

    ! In/output variables:
    character(len=3),                                    intent(in   ) :: region_name
    type(type_mesh),                                     intent(inout) :: mesh
    type(type_reference_geometry),                       intent(in)    :: refgeo

    ! Local variables:
    character(len=256), parameter                                      :: routine_name = 'setup_mesh'
    character(len=256)                                                 :: mesh_name

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to screen
    if (par%primary) write(0,'(A)') '   Setting up the mesh ...'

    ! Mesh name
    mesh_name = 'model_mesh_' // trim( region_name) // '_00001'

    ! Calculate a new mesh based on the initial ice-sheet geometry, or read an existing mesh from a file
    call setup_mesh_from_geometry( region_name, mesh_name, mesh, refgeo)

    ! Write the mesh creation success message to the terminal
    call write_mesh_success( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_mesh

  subroutine setup_mesh_from_geometry( region_name, mesh_name, mesh, refgeo)
    ! Set up the model mesh based on the ice-sheet geometry

    ! In/output variables:
    character(len=3),                                    intent(in   ) :: region_name
    character(len=256),                                  intent(in)    :: mesh_name
    type(type_mesh),                                     intent(inout) :: mesh
    type(type_reference_geometry),                       intent(in)    :: refgeo

    ! Local variables:
    character(len=256), parameter                                      :: routine_name = 'setup_mesh_from_geometry'
    real(dp)                                                           :: xmin, xmax, ymin, ymax
    real(dp)                                                           :: lambda_M, phi_M, beta_stereo

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to screen
    if (par%primary) write(0,'(A)') '     Creating mesh from geometry...'

    ! Determine model domain
    xmin = C%xmin_ANT
    xmax = C%xmax_ANT
    ymin = C%ymin_ANT
    ymax = C%ymax_ANT

    ! Determine if the geometry is provided gridded or meshed
    if (allocated( refgeo%grid_raw%x)) then
      ! Gridded

      ! Safety
      if (allocated( refgeo%mesh_raw%V)) call crash('found both grid and mesh in refgeo!')

      ! Create mesh from gridded geometry
      call create_mesh_from_gridded_geometry( region_name, mesh_name, &
        refgeo%grid_raw, &
        refgeo%Hi_grid_raw, &
        refgeo%Hb_grid_raw, &
        refgeo%Hs_grid_raw, &
        refgeo%SL_grid_raw, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        mesh)

    elseif (allocated( refgeo%mesh_raw%V)) then
      ! Meshed

      ! Safety
      if (allocated( refgeo%grid_raw%x)) call crash('found both grid and mesh in refgeo!')

      ! Create mesh from meshed initial geometry data
      call create_mesh_from_meshed_geometry( region_name, mesh_name, &
        refgeo%mesh_raw, &
        refgeo%Hi_mesh_raw, &
        refgeo%Hb_mesh_raw, &
        refgeo%Hs_mesh_raw, &
        refgeo%SL_mesh_raw, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        mesh)

    else
      call crash('no grid or mesh is found in refgeo!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_mesh_from_geometry

end module LADDIE_main_model
