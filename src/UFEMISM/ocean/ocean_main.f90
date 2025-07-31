MODULE ocean_main

  ! The main ocean model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE ocean_utilities                                        , ONLY: initialise_ocean_vertical_grid, calc_ocean_temperature_at_shelf_base, calc_ocean_freezing_point_at_shelf_base
  USE ocean_realistic                                        , ONLY: initialise_ocean_model_realistic, run_ocean_model_realistic
  USE ocean_idealised                                        , ONLY: initialise_ocean_model_idealised, run_ocean_model_idealised
  use netcdf_io_main

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ocean_model( mesh, ice, ocean, region_name, time)
    ! Calculate the ocean

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new ocean
    IF (C%do_asynchronous_ocean) THEN
      ! Asynchronous coupling: do not calculate a new ocean in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next ocean time step
      IF (time == ocean%t_next) THEN
        ! Go on to calculate a new ocean
        ocean%t_next = time + C%dt_ocean
      ELSEIF (time > ocean%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the ocean time step')
      ELSE
        ! It is not yet time to calculate a new ocean
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_ocean) THEN
      ! Synchronous coupling: calculate a new ocean in every model loop
      ocean%t_next = time + C%dt_ocean
    END IF

    ! Determine which ocean model to run for this region
    IF     (region_name == 'NAM') THEN
      choice_ocean_model = C%choice_ocean_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_ocean_model = C%choice_ocean_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_ocean_model = C%choice_ocean_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_ocean_model = C%choice_ocean_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Run the chosen ocean model
    IF (choice_ocean_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_ocean_model == 'idealised') THEN
      CALL run_ocean_model_idealised( mesh, ice, ocean)
    ELSEIF (choice_ocean_model == 'realistic') THEN
      CALL run_ocean_model_realistic( mesh, ice, ocean, time)
    ELSE
      CALL crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    END IF

    ! Compute secondary variables
    CALL calc_ocean_temperature_at_shelf_base(    mesh, ice, ocean)
    CALL calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model

  SUBROUTINE initialise_ocean_model( mesh, ice, ocean, region_name, start_time_of_run)
    ! Initialise the ocean model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(OUT)   :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: start_time_of_run

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising ocean model...'

    ! Determine which ocean model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_ocean_model = C%choice_ocean_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_ocean_model = C%choice_ocean_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_ocean_model = C%choice_ocean_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_ocean_model = C%choice_ocean_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Initialise vertical grid: C%z_ocean and C%nz_ocean
    CALL initialise_ocean_vertical_grid

    ! Allocate memory for main variables
    ALLOCATE( ocean%T( mesh%vi1:mesh%vi2,C%nz_ocean))
    ALLOCATE( ocean%S( mesh%vi1:mesh%vi2,C%nz_ocean))
    ocean%T = 0._dp
    ocean%S = 0._dp

    ! Allocate memory for secondary variables
    ALLOCATE( ocean%T_draft(          mesh%vi1:mesh%vi2))
    ALLOCATE( ocean%T_freezing_point( mesh%vi1:mesh%vi2))
    ocean%T_draft          = 0._dp
    ocean%T_freezing_point = 0._dp

    ! Set time of next calculation to start time
    ocean%t_next = C%start_time_of_run

    ! Determine which ocean model to initialise
    IF     (choice_ocean_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_ocean_model == 'idealised') THEN
      CALL initialise_ocean_model_idealised( mesh, ocean)
    ELSEIF (choice_ocean_model == 'realistic') THEN
      CALL initialise_ocean_model_realistic( mesh, ice, ocean, region_name, start_time_of_run)
    ELSE
      CALL crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model

  SUBROUTINE write_to_restart_file_ocean_model( mesh, ocean, region_name, time)
    ! Write to the restart file for the ocean model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which ocean model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_ocean_model = C%choice_ocean_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_ocean_model = C%choice_ocean_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_ocean_model = C%choice_ocean_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_ocean_model = C%choice_ocean_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Write to the restart file of the chosen ocean model
    IF     (choice_ocean_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_ocean_model == 'idealised') THEN
      ! No need to do anything
    ELSEIF (choice_ocean_model == 'realistic') THEN
      CALL write_to_restart_file_ocean_model_region( mesh, ocean, region_name, time)
    ELSE
      CALL crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_ocean_model

  SUBROUTINE write_to_restart_file_ocean_model_region( mesh, ocean, region_name, time)
    ! Write to the restart NetCDF file for the ocean model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),        INTENT(IN)    :: mesh
    TYPE(type_ocean_model), INTENT(IN)    :: ocean
    CHARACTER(LEN=3),       INTENT(IN)    :: region_name
    REAL(dp),               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER         :: routine_name = 'write_to_restart_file_ocean_model_region'
    INTEGER                               :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Writing to ocean restart file "' // &
      colour_string( TRIM( ocean%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( ocean%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( ocean%restart_filename, ncid, time)

    ! ! Write the velocity fields to the file
    CALL write_to_field_multopt_mesh_dp_3D_ocean( mesh, ocean%restart_filename, ncid, 'T', ocean%T)
    CALL write_to_field_multopt_mesh_dp_3D_ocean( mesh, ocean%restart_filename, ncid, 'S', ocean%S)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_ocean_model_region

  SUBROUTINE create_restart_file_ocean_model( mesh, ocean, region_name)
    ! Create the restart file for the ocean model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which ocean model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_ocean_model = C%choice_ocean_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_ocean_model = C%choice_ocean_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_ocean_model = C%choice_ocean_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_ocean_model = C%choice_ocean_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Create the restart file of the chosen ocean model
    IF     (choice_ocean_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_ocean_model == 'idealised') THEN
      ! No need to do anything
    ELSEIF (choice_ocean_model == 'realistic') THEN
      CALL create_restart_file_ocean_model_region( mesh, ocean, region_name)
    ELSE
      CALL crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_ocean_model

  SUBROUTINE create_restart_file_ocean_model_region( mesh, ocean, region_name)
    ! Create a restart NetCDF file for the ocean submodel
    ! Includes generation of the procedural filename (e.g. "restart_ocean_00001.nc")

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),        INTENT(IN)    :: mesh
    TYPE(type_ocean_model), INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER         :: routine_name = 'create_restart_file_ocean_model_region'
    CHARACTER(LEN=256)                    :: filename_base
    INTEGER                               :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_ocean_' // region_name
    CALL generate_filename_XXXXXdotnc( filename_base, ocean%restart_filename)

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Creating ocean model restart file "' // &
      colour_string( TRIM( ocean%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( ocean%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( ocean%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( ocean%restart_filename, ncid)

    ! Add a depth dimension to the file
    CALL add_depth_dimension_to_file( ocean%restart_filename, ncid, C%z_ocean)

    ! Add the data fields to the file
    CALL add_field_mesh_dp_3D_ocean( ocean%restart_filename, ncid, 'T' , long_name = 'Ocean temperatures', units = 'degrees C')
    CALL add_field_mesh_dp_3D_ocean( ocean%restart_filename, ncid, 'S' , long_name = 'Ocean salinity', units = 'PSU')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_ocean_model_region

  SUBROUTINE remap_ocean_model( mesh_old, mesh_new, ice, ocean, region_name, time)
    ! Remap the ocean model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping ocean model data to the new mesh...'

    ! Determine which ocean model to initialise for this region
    IF     (region_name == 'NAM') THEN
      choice_ocean_model = C%choice_ocean_model_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_ocean_model = C%choice_ocean_model_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_ocean_model = C%choice_ocean_model_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_ocean_model = C%choice_ocean_model_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Reallocate memory for main variables
    CALL reallocate_bounds( ocean%T, mesh_new%vi1, mesh_new%vi2, C%nz_ocean)
    CALL reallocate_bounds( ocean%S, mesh_new%vi1, mesh_new%vi2, C%nz_ocean)

    ! Reallocate memory for secondary variables
    CALL reallocate_bounds( ocean%T_draft,          mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( ocean%T_freezing_point, mesh_new%vi1, mesh_new%vi2)

    ! Determine which ocean model to remap
    IF     (choice_ocean_model == 'none') THEN
      ! No need to do anything
    ELSEIF (choice_ocean_model == 'idealised') THEN
      CALL initialise_ocean_model_idealised( mesh_new, ocean)
    ELSEIF (choice_ocean_model == 'realistic') THEN
        IF     (C%choice_ocean_model_realistic == 'snapshot') THEN
          CALL initialise_ocean_model_realistic( mesh_new, ice, ocean, region_name, time)
        ELSE
          CALL crash('Remapping after mesh update for realistic ocean is only implemented for a snapshot!')
        END IF
    ELSE
      CALL crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ocean_model

END MODULE ocean_main
