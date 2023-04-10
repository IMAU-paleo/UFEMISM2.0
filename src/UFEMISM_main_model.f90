MODULE UFEMISM_main_model

  ! The main regional model

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE reference_geometries                                   , ONLY: type_reference_geometry, initialise_reference_geometries_raw
  USE ice_model_types                                        , ONLY: type_ice_model
  USE netcdf_basic                                           , ONLY: open_existing_netcdf_file_for_reading, close_netcdf_file
  USE netcdf_input                                           , ONLY: setup_mesh_from_file
  USE mesh_creation                                          , ONLY: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry

  ! DENK DROM
  USE netcdf_basic , ONLY: create_new_netcdf_file_for_writing
  USE netcdf_output, ONLY: setup_mesh_in_netcdf_file

  IMPLICIT NONE

  ! Data structure containing a single UFEMISM model region
  TYPE type_model_region

    ! Metadata
    CHARACTER(LEN=3)                        :: name                        ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                      :: long_name                   ! North America, Eurasia, Greenland, Antarctica

    ! The current time of this particular region.
    REAL(dp)                                :: time

    ! The mesh that all model components define their data on
    TYPE(type_mesh)                         :: mesh

    ! Reference geometries
    TYPE(type_reference_geometry)           :: refgeo_init
    TYPE(type_reference_geometry)           :: refgeo_PD
    TYPE(type_reference_geometry)           :: refgeo_GIAeq

    ! The ice dynamics model
    TYPE(type_ice_model)                    :: ice

  END TYPE type_model_region

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
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_region

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

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Initialise all the reference geometries
    CALL initialise_reference_geometries_raw( region%name, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq)

    ! Set up the first model mesh
    CALL setup_first_mesh( region)

    ! DENK DROM
    CALL crash('whoopsiedaisy!')

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') ''
    IF (par%master) WRITE(0,'(A)') ' Finished initialising model region ' // colour_string( region%long_name,'light blue')

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

    CHARACTER(LEN=256) :: filename
    INTEGER :: ncid

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

    IF     (choice_initial_mesh == 'calc_from_initial_geometry') THEN
      CALL setup_first_mesh_from_initial_geometry( region)
    ELSEIF (choice_initial_mesh == 'read_from_file') THEN
      CALL setup_first_mesh_from_file( region)
    ELSE
      CALL crash('unknown choice_initial_mesh "' // TRIM( choice_initial_mesh) // '"!')
    END IF

    ! DENK DROM
    filename = TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, region%mesh)
    CALL close_netcdf_file( ncid)

    ! DENK DROM
    CALL crash('whoopsiedaisy!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_first_mesh

  SUBROUTINE setup_first_mesh_from_initial_geometry( region)
    ! Set up the first model mesh based on the initial ice-sheet geometry

    IMPLICIT NONE

    ! In/output variables:
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
      CALL create_mesh_from_gridded_geometry( &
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

      ! Create mesh from gridded initial geometry data
      CALL create_mesh_from_meshed_geometry( &
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

  SUBROUTINE setup_first_mesh_from_file( region)
    ! Set up the first model mesh from an external NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'setup_first_mesh'
    CHARACTER(LEN=256)                                                 :: filename_initial_mesh
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

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
    IF (par%master) WRITE(0,'(A)') '   Reading mesh from file "' // colour_string( region%name,'light blue') // '"...'

    ! Read the mest from the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename_initial_mesh, ncid)
    CALL setup_mesh_from_file(                  filename_initial_mesh, ncid, region%mesh)
    CALL close_netcdf_file(                                            ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_first_mesh_from_file

END MODULE UFEMISM_main_model