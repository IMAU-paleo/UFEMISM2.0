MODULE reference_geometries

  ! Contains the routines for setting up the three "reference geometries":
  ! - refgeo_init:   initial, used to initialise the simulation
  ! - refgeo_PD:     present-day, used to calculate sea-level contribution, isotope change, and more
  ! - refgeo_GIA_eq: GIA equilibrium, used for the GIA model

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid, setup_square_grid, distribute_gridded_data_from_master_dp_2D
  USE math_utilities                                         , ONLY: ice_surface_elevation
  USE analytical_solutions                                   , ONLY: Halfar_dome, Bueler_dome
  USE netcdf_basic                                           , ONLY: inquire_xy_grid, inquire_lonlat_grid, inquire_mesh, open_existing_netcdf_file_for_reading, &
                                                                     inquire_var_multopt, close_netcdf_file
  USE netcdf_input                                           , ONLY: read_field_from_xy_file_2D, read_field_from_mesh_file_2D
  USE mesh_remapping                                         , ONLY: map_from_xy_grid_to_mesh_2D, map_from_mesh_to_mesh_2D

  IMPLICIT NONE

  ! Data structure containing a reference ice-sheet geometry
  TYPE type_reference_geometry

    ! Data on the raw grid/mesh
    TYPE(type_grid)                         :: grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb_grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs_grid_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL_grid_raw
    TYPE(type_mesh)                         :: mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi_mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb_mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs_mesh_raw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL_mesh_raw

    ! Data on the model mesh
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi                          ! [m] Ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb                          ! [m] Bedrock elevation [w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                          ! [m] Surface elevation [w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL                          ! [m] Sea surface elevation [w.r.t. PD sea level]

  END TYPE type_reference_geometry

CONTAINS

! ===== Subroutines =====
! =======================

  ! Initialise reference geometries on their raw input grid/mesh
  ! ============================================================

  SUBROUTINE initialise_reference_geometries_raw( region_name, refgeo_init, refgeo_PD, refgeo_GIAeq)
    ! Initialise a reference geometry on the raw grid/mesh for the given set of config choices

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3)                                   , INTENT(IN)    :: region_name
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo_init
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo_PD
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo_GIAeq

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_reference_geometries_raw'
    CHARACTER(LEN=256)                                                 :: choice_refgeo
    CHARACTER(LEN=256)                                                 :: choice_refgeo_idealised
    REAL(dp)                                                           :: dx_refgeo_idealised
    CHARACTER(LEN=256)                                                 :: filename_refgeo
    REAL(dp)                                                           :: timeframe_refgeo

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Initial geometry
    ! ===================

    ! Get config choices for this model region
    IF     (region_name == 'NAM') THEN
      choice_refgeo           = C%choice_refgeo_init_NAM
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_NAM
      timeframe_refgeo        = C%timeframe_refgeo_init_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_refgeo           = C%choice_refgeo_init_EAS
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_EAS
      timeframe_refgeo        = C%timeframe_refgeo_init_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_refgeo           = C%choice_refgeo_init_GRL
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_GRL
      timeframe_refgeo        = C%timeframe_refgeo_init_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_refgeo           = C%choice_refgeo_init_ANT
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_ANT
      timeframe_refgeo        = C%timeframe_refgeo_init_ANT
    ELSE
      CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END IF

    ! Initialise reference geometry
    CALL initialise_reference_geometry_raw( region_name, 'initial', refgeo_init, choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)

    ! == Present-day geometry
    ! =======================

    ! Get config choices for this model region
    IF     (region_name == 'NAM') THEN
      choice_refgeo           = C%choice_refgeo_PD_NAM
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_NAM
      timeframe_refgeo        = C%timeframe_refgeo_PD_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_refgeo           = C%choice_refgeo_PD_EAS
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_EAS
      timeframe_refgeo        = C%timeframe_refgeo_PD_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_refgeo           = C%choice_refgeo_PD_GRL
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_GRL
      timeframe_refgeo        = C%timeframe_refgeo_PD_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_refgeo           = C%choice_refgeo_PD_ANT
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_ANT
      timeframe_refgeo        = C%timeframe_refgeo_PD_ANT
    ELSE
      CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END IF

    ! Initialise reference geometry
    CALL initialise_reference_geometry_raw( region_name, 'present-day', refgeo_PD, choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)

    ! == GIA equilibrium geometry
    ! ===========================

    ! Get config choices for this model region
    IF     (region_name == 'NAM') THEN
      choice_refgeo           = C%choice_refgeo_GIAeq_NAM
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_NAM
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_refgeo           = C%choice_refgeo_GIAeq_EAS
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_EAS
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_refgeo           = C%choice_refgeo_GIAeq_GRL
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_GRL
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_refgeo           = C%choice_refgeo_GIAeq_ANT
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_ANT
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_ANT
    ELSE
      CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END IF

    ! Initialise reference geometry
    CALL initialise_reference_geometry_raw( region_name, 'GIA equilibrium', refgeo_GIAeq, choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometries_raw

  SUBROUTINE initialise_reference_geometry_raw( region_name, refgeo_name, refgeo, choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)
    ! Initialise a reference geometry on the raw grid/mesh for the given set of config choices

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3)                                   , INTENT(IN)    :: region_name
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: refgeo_name
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_refgeo_idealised
    REAL(dp)                                           , INTENT(IN)    :: dx_refgeo_idealised
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename_refgeo
    REAL(dp)                                           , INTENT(IN)    :: timeframe_refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_reference_geometry_raw'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (choice_refgeo == 'idealised') THEN
      CALL initialise_reference_geometry_raw_idealised( region_name, refgeo_name, refgeo, choice_refgeo_idealised, dx_refgeo_idealised)
    ELSEIF (choice_refgeo == 'read_from_file') THEN
      CALL initialise_reference_geometry_raw_from_file( region_name, refgeo_name, refgeo, choice_refgeo, filename_refgeo, timeframe_refgeo)
    ELSE
      CALL crash('unknown choice_refgeo "' // TRIM( choice_refgeo) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_raw

  SUBROUTINE initialise_reference_geometry_raw_idealised( region_name, refgeo_name, refgeo, choice_refgeo_idealised, dx_refgeo_idealised)
    ! Initialise a reference geometry on the raw grid/mesh for the given set of config choices
    !
    ! For the case of an internally-generated idealised geometry

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3)                                   , INTENT(IN)    :: region_name
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: refgeo_name
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_refgeo_idealised
    REAL(dp)                                           , INTENT(IN)    :: dx_refgeo_idealised

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_reference_geometry_raw_idealised'
    REAL(dp)                                                           :: xmin, xmax, ymin, ymax
    CHARACTER(LEN=256)                                                 :: name
    TYPE(type_grid)                                                    :: grid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: Hi, Hb, Hs, SL
    INTEGER                                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') '  Initialising ' // TRIM( refgeo_name) // ' geometry for model region ' // &
      colour_string( region_name,'light blue') // ' from idealised case "' // colour_string( TRIM( choice_refgeo_idealised),'light blue') // '"...'

    ! Get domain size for this model region
    IF     (region_name == 'NAM') THEN
      xmin = C%xmin_NAM
      xmax = C%xmax_NAM
      ymin = C%ymin_NAM
      ymax = C%ymax_NAM
    ELSEIF (region_name == 'EAS') THEN
      xmin = C%xmin_EAS
      xmax = C%xmax_EAS
      ymin = C%ymin_EAS
      ymax = C%ymax_EAS
    ELSEIF (region_name == 'GRL') THEN
      xmin = C%xmin_GRL
      xmax = C%xmax_GRL
      ymin = C%ymin_GRL
      ymax = C%ymax_GRL
    ELSEIF (region_name == 'ANT') THEN
      xmin = C%xmin_ANT
      xmax = C%xmax_ANT
      ymin = C%ymin_ANT
      ymax = C%ymax_ANT
    ELSE
      CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END IF

    ! Set up a square grid to generate the idealised geometry on
    name = 'temp_grid_for_idealised_geometry_lines'
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx_refgeo_idealised, refgeo%grid_raw)

    ! Allocate memory for partial grid data
    ALLOCATE( refgeo%Hi_grid_raw( refgeo%grid_raw%n_loc), source = 0._dp)
    ALLOCATE( refgeo%Hb_grid_raw( refgeo%grid_raw%n_loc), source = 0._dp)
    ALLOCATE( refgeo%Hs_grid_raw( refgeo%grid_raw%n_loc), source = 0._dp)
    ALLOCATE( refgeo%SL_grid_raw( refgeo%grid_raw%n_loc), source = 0._dp)

    ! Allocate memory for the full grid data on the master
    IF (par%master) THEN
      ALLOCATE( Hi( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      ALLOCATE( Hb( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      ALLOCATE( Hs( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      ALLOCATE( SL( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
    END IF

    ! Calculate the idealised geometry on the grid (master only)
    IF (par%master) THEN
      DO i = 1, refgeo%grid_raw%nx
      DO j = 1, refgeo%grid_raw%ny
        CALL calc_idealised_geometry( refgeo%grid_raw%x( i), refgeo%grid_raw%y( j), Hi( i,j), Hb( i,j), Hs( i,j), SL( i,j), choice_refgeo_idealised)
      END DO
      END DO
    END IF

    ! Distribute the data over the processes in vector form
    CALL distribute_gridded_data_from_master_dp_2D( refgeo%grid_raw, Hi, refgeo%Hi_grid_raw)
    CALL distribute_gridded_data_from_master_dp_2D( refgeo%grid_raw, Hb, refgeo%Hb_grid_raw)
    CALL distribute_gridded_data_from_master_dp_2D( refgeo%grid_raw, Hs, refgeo%Hs_grid_raw)
    CALL distribute_gridded_data_from_master_dp_2D( refgeo%grid_raw, SL, refgeo%SL_grid_raw)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( Hi)
      DEALLOCATE( Hb)
      DEALLOCATE( Hs)
      DEALLOCATE( SL)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_raw_idealised

  SUBROUTINE initialise_reference_geometry_raw_from_file( region_name, refgeo_name, refgeo, choice_refgeo, filename_refgeo, timeframe_refgeo)
    ! Initialise a reference geometry on the raw grid/mesh for the given set of config choices!
    !
    ! For the case of a (probably) realistic geometry provided through a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3)                                   , INTENT(IN)    :: region_name
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: refgeo_name
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename_refgeo
    REAL(dp)                                           , INTENT(IN)    :: timeframe_refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_reference_geometry_raw_from_file'
    LOGICAL                                                            :: has_xy_grid, has_lonlat_grid, has_mesh

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to screen
    IF (par%master) WRITE(0,'(A)') '  Initialising ' // TRIM( refgeo_name) // ' geometry for model region ' // &
      colour_string( region_name,'light blue') // ' from file "' // colour_string( TRIM( filename_refgeo),'light blue') // '"...'


    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename_refgeo, has_xy_grid    )
    CALL inquire_lonlat_grid( filename_refgeo, has_lonlat_grid)
    CALL inquire_mesh(        filename_refgeo, has_mesh       )

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename_refgeo) // '" contains both an x/y-grid and a lon/lat-grid!')
    IF (has_xy_grid     .AND. has_mesh       ) CALL crash('file "' // TRIM( filename_refgeo) // '" contains both an x/y-grid and a mesh!')
    IF (has_lonlat_grid .AND. has_mesh       ) CALL crash('file "' // TRIM( filename_refgeo) // '" contains both a lon/lat-grid and a mesh!')

    ! Read the grid'mesh and data from the file
    IF (has_xy_grid) THEN
      ! Read reference ice sheet geometry data from an xy-gridded file
      CALL initialise_reference_geometry_raw_from_file_grid( region_name, refgeo_name, refgeo, choice_refgeo, filename_refgeo, timeframe_refgeo)
    ELSEIF (has_mesh) THEN
      ! Read reference ice sheet geometry data from a meshed file
      CALL initialise_reference_geometry_raw_from_file_mesh( region_name, refgeo_name, refgeo, choice_refgeo, filename_refgeo, timeframe_refgeo)
    ELSE
      CALL crash('can only read reference geometry from gridded or meshed data files!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_raw_from_file

  SUBROUTINE initialise_reference_geometry_raw_from_file_grid( region_name, refgeo_name, refgeo, choice_refgeo, filename_refgeo, timeframe_refgeo)
    ! Initialise a reference geometry on the raw grid/mesh for the given set of config choices!
    !
    ! For the case of a (probably) realistic geometry provided through a gridded NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3)                                   , INTENT(IN)    :: region_name
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: refgeo_name
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename_refgeo
    REAL(dp)                                           , INTENT(IN)    :: timeframe_refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_reference_geometry_raw_from_file_grid'
    INTEGER                                                            :: ncid, id_var
    LOGICAL                                                            :: has_SL

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if a sea level variable exists in the file
    CALL open_existing_netcdf_file_for_reading( filename_refgeo, ncid)
    CALL inquire_var_multopt( filename_refgeo, ncid, 'default_options_SL', id_var)
    has_SL = id_var /= -1
    CALL close_netcdf_file( ncid)

    IF (timeframe_refgeo /= 1E9_dp) THEN
      ! We need to read a specific time frame

      CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_grid_raw, time_to_read = timeframe_refgeo, grid = refgeo%grid_raw)
      CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_grid_raw, time_to_read = timeframe_refgeo)
      CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_grid_raw, time_to_read = timeframe_refgeo)

      ! If the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      IF (has_SL) THEN
        CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_SL', refgeo%SL_grid_raw, time_to_read = timeframe_refgeo)
      ELSE
        ALLOCATE( refgeo%SL_grid_raw( refgeo%grid_raw%n_loc), source = 0._dp)
      END IF

    ELSE !  IF (timeframe_refgeo /= 1E9_dp) THEN
      ! We need to read data from a time-less NetCDF file

      CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_grid_raw, grid = refgeo%grid_raw)
      CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_grid_raw)
      CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_grid_raw)

      ! If the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      IF (has_SL) THEN
        CALL read_field_from_xy_file_2D( filename_refgeo, 'default_options_SL', refgeo%SL_grid_raw)
      ELSE
        ALLOCATE( refgeo%SL_grid_raw( refgeo%grid_raw%n_loc), source = 0._dp)
      END IF

    END IF !  IF (timeframe_refgeo /= 1E9_dp) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_raw_from_file_grid

  SUBROUTINE initialise_reference_geometry_raw_from_file_mesh( region_name, refgeo_name, refgeo, choice_refgeo, filename_refgeo, timeframe_refgeo)
    ! Initialise a reference geometry on the raw grid/mesh for the given set of config choices!
    !
    ! For the case of a (probably) realistic geometry provided through a meshed NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3)                                   , INTENT(IN)    :: region_name
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: refgeo_name
    TYPE(type_reference_geometry)                      , INTENT(OUT)   :: refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_refgeo
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename_refgeo
    REAL(dp)                                           , INTENT(IN)    :: timeframe_refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_reference_geometry_raw_from_file_mesh'
    INTEGER                                                            :: ncid, id_var
    LOGICAL                                                            :: has_SL

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if a sea level variable exists in the file
    CALL open_existing_netcdf_file_for_reading( filename_refgeo, ncid)
    CALL inquire_var_multopt( filename_refgeo, ncid, 'default_options_SL', id_var)
    has_SL = id_var /= -1
    CALL close_netcdf_file( ncid)

    IF (timeframe_refgeo /= 1E9_dp) THEN
      ! We need to read a specific time frame

      CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_mesh_raw, time_to_read = timeframe_refgeo, mesh = refgeo%mesh_raw)
      CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_mesh_raw, time_to_read = timeframe_refgeo)
      CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_mesh_raw, time_to_read = timeframe_refgeo)

      ! If the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      IF (has_SL) THEN
        CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_SL', refgeo%SL_mesh_raw, time_to_read = timeframe_refgeo)
      ELSE
        ALLOCATE( refgeo%SL_mesh_raw( refgeo%mesh_raw%nV_loc), source = 0._dp)
      END IF

    ELSE !  IF (timeframe_refgeo /= 1E9_dp) THEN
      ! We need to read data from a time-less NetCDF file

      CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_mesh_raw, mesh = refgeo%mesh_raw)
      CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_mesh_raw)
      CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_mesh_raw)

      ! If the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      IF (has_SL) THEN
        CALL read_field_from_mesh_file_2D( filename_refgeo, 'default_options_SL', refgeo%SL_mesh_raw)
      ELSE
        ALLOCATE( refgeo%SL_mesh_raw( refgeo%mesh_raw%nV_loc), source = 0._dp)
      END IF

    END IF !  IF (timeframe_refgeo /= 1E9_dp) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_raw_from_file_mesh

  ! Remap reference geometry to the model mesh
  ! ==========================================

  SUBROUTINE remap_reference_geometry_to_mesh( mesh, refgeo)
    ! Remap reference geometry to the model mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh)                                    , INTENT(IN)    :: mesh
    TYPE(type_reference_geometry)                      , INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'remap_reference_geometry_to_mesh'
    CHARACTER(LEN=256)                                                 :: method
    INTEGER                                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Deallocate existing memory if needed
    IF (ALLOCATED( refgeo%Hi)) DEALLOCATE( refgeo%Hi)
    IF (ALLOCATED( refgeo%Hb)) DEALLOCATE( refgeo%Hb)
    IF (ALLOCATED( refgeo%Hs)) DEALLOCATE( refgeo%Hs)
    IF (ALLOCATED( refgeo%SL)) DEALLOCATE( refgeo%SL)

    ! Allocate memory for reference ice geometry on the model mesh
    ALLOCATE( refgeo%Hi( mesh%vi1:mesh%vi2))
    ALLOCATE( refgeo%Hb( mesh%vi1:mesh%vi2))
    ALLOCATE( refgeo%Hs( mesh%vi1:mesh%vi2))
    ALLOCATE( refgeo%SL( mesh%vi1:mesh%vi2))

    ! Determine if the initial geometry is provided gridded or meshed
    IF (ALLOCATED( refgeo%grid_raw%x)) THEN
      ! Gridded

      ! Safety
      IF (ALLOCATED( refgeo%mesh_raw%V)) CALL crash('found boht grid and mesh in refgeo!')

      ! Remap data to the model mesh
      CALL map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, refgeo%Hi_grid_raw, refgeo%Hi)
      CALL map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, refgeo%Hb_grid_raw, refgeo%Hb)
      CALL map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, refgeo%SL_grid_raw, refgeo%SL)

    ELSEIF (ALLOCATED( refgeo%mesh_raw%V)) THEN
      ! Meshed

      ! Safety
      IF (ALLOCATED( refgeo%grid_raw%x)) CALL crash('found boht grid and mesh in refgeo!')

      ! Remap data to the model mesh
      method = '2nd_order_conservative'
      CALL map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, refgeo%Hi_mesh_raw, refgeo%Hi, method)
      CALL map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, refgeo%Hb_mesh_raw, refgeo%Hb, method)
      CALL map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, refgeo%SL_mesh_raw, refgeo%SL, method)

    ELSE
      CALL crash('no grid or mesh is found in refgeo!')
    END IF

    ! Don't remap Hs, but recalculate it after remapping Hi,Hb,SL
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = ice_surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), refgeo%SL( vi))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_reference_geometry_to_mesh

  ! Initialise idealised geometries on the model mesh
  ! =================================================

  SUBROUTINE initialise_reference_geometry_idealised( mesh, choice_refgeo_idealised, refgeo)
    ! Initialise an idealised reference geometry on the model mesh

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh)                                    , INTENT(IN)    :: mesh
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_refgeo_idealised
    TYPE(type_reference_geometry)                      , INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'initialise_reference_geometry_idealised'
    INTEGER                                                            :: vii,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vii = 1, mesh%nV_loc
      vi = mesh%vi1 + vii - 1
      CALL calc_idealised_geometry( mesh%V( vi,1), mesh%V( vi,2), refgeo%Hi( vii), refgeo%Hb( vii), refgeo%Hs( vii), refgeo%SL( vii), choice_refgeo_idealised)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised

  ! Calculate various idealsed geometries
  ! =====================================

  SUBROUTINE calc_idealised_geometry( x, y, Hi, Hb, Hs, SL, choice_refgeo_idealised)
    ! Calculate an idealised geometry

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level
    CHARACTER(LEN=256),             INTENT(IN)    :: choice_refgeo_idealised

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Calculated the specified idealised geometry
    IF     (choice_refgeo_idealised == 'flatearth') THEN
      CALL calc_idealised_geometry_flatearth( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'slabonaslope') THEN
      CALL calc_idealised_geometry_slabonaslope( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'Halfar') THEN
      CALL calc_idealised_geometry_Halfar( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'Bueler') THEN
      CALL calc_idealised_geometry_Bueler( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'SSA_icestream') THEN
      CALL calc_idealised_geometry_SSA_icestream( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'MISMIP_mod') THEN
      CALL calc_idealised_geometry_MISMIP_mod( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'ISMIP-HOM_A') THEN
      CALL calc_idealised_geometry_ISMIP_HOM_A( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'ISMIP-HOM_B') THEN
      CALL calc_idealised_geometry_ISMIP_HOM_B( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'ISMIP-HOM_C' .OR. &
            choice_refgeo_idealised == 'ISMIP-HOM_D') THEN
      CALL calc_idealised_geometry_ISMIP_HOM_CD( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'ISMIP-HOM_E') THEN
      CALL calc_idealised_geometry_ISMIP_HOM_E( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'ISMIP-HOM_F') THEN
      CALL calc_idealised_geometry_ISMIP_HOM_F( x, y, Hi, Hb, Hs, SL)
    ELSEIF (choice_refgeo_idealised == 'MISMIP+' .OR. &
            choice_refgeo_idealised == 'MISMIPplus') THEN
      CALL calc_idealised_geometry_MISMIPplus( x, y, Hi, Hb, Hs, SL)
    ELSE
      CALL crash('unknown choice_refgeo_idealised "' // TRIM( choice_refgeo_idealised) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry

  SUBROUTINE calc_idealised_geometry_flatearth( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_flatearth'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    Hi = 0._dp
    Hb = 0._dp
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_flatearth

  SUBROUTINE calc_idealised_geometry_slabonaslope( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! A 2,000 m thick slab of ice on a flat, inclined plane (10 m per km slope)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_slabonaslope'
    REAL(dp), PARAMETER                           :: dhdx = -0.01_dp

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_slabonaslope_Hi < 100._dp .OR. &
         C%refgeo_idealised_slabonaslope_Hi > 10000._dp) THEN
      CALL crash('refgeo_idealised_slabonaslope_Hi has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_slabonaslope_Hi)
    END IF
    IF ( ABS( C%refgeo_idealised_slabonaslope_dhdx) > 0.1_dp) THEN
      CALL crash('refgeo_idealised_slabonaslope_dhdx has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_slabonaslope_dhdx)
    END IF

    Hi = C%refgeo_idealised_slabonaslope_Hi
    Hb = C%refgeo_idealised_slabonaslope_dhdx * x
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_slabonaslope

  SUBROUTINE calc_idealised_geometry_Halfar( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! The Halfar dome solution at t = 0

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_Halfar'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%uniform_flow_factor < 1E-18_dp .OR. &
         C%uniform_flow_factor > 1E-15_dp) THEN
      CALL crash('uniform_flow_factor has unrealistic value of {dp_01}!', dp_01 = C%uniform_flow_factor)
    END IF
    IF ( C%n_flow < 1._dp .OR. C%n_flow > 5._dp) THEN
      CALL crash('n_flow has unrealistic value of {dp_01}!', dp_01 = C%n_flow)
    END IF
    IF ( C%refgeo_idealised_Halfar_H0 < 100._dp .OR. C%refgeo_idealised_Halfar_H0 > 10000._dp) THEN
      CALL crash('refgeo_idealised_Halfar_H0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Halfar_H0)
    END IF
    IF ( C%refgeo_idealised_Halfar_R0 < 100E3_dp .OR. C%refgeo_idealised_Halfar_R0 > 5000E3_dp) THEN
      CALL crash('refgeo_idealised_Halfar_R0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Halfar_R0)
    END IF

    CALL Halfar_dome( C%uniform_flow_factor, C%n_flow, C%refgeo_idealised_Halfar_H0, C%refgeo_idealised_Halfar_R0, &
      x, y, 0._dp, Hi)
    Hb = 0._dp
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_Halfar

  SUBROUTINE calc_idealised_geometry_Bueler( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! The Bueler dome solution at t = 0

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_Bueler'
    REAL(dp)                                      :: M

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%uniform_flow_factor < 1E-18_dp .OR. &
         C%uniform_flow_factor > 1E-15_dp) THEN
      CALL crash('uniform_flow_factor has unrealistic value of {dp_01}!', dp_01 = C%uniform_flow_factor)
    END IF
    IF ( C%n_flow < 1._dp .OR. C%n_flow > 5._dp) THEN
      CALL crash('n_flow has unrealistic value of {dp_01}!', dp_01 = C%n_flow)
    END IF
    IF ( C%refgeo_idealised_Bueler_H0 < 100._dp .OR. C%refgeo_idealised_Bueler_H0 > 10000._dp) THEN
      CALL crash('refgeo_idealised_Bueler_H0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_H0)
    END IF
    IF ( C%refgeo_idealised_Bueler_R0 < 100E3_dp .OR. C%refgeo_idealised_Bueler_R0 > 5000E3_dp) THEN
      CALL crash('refgeo_idealised_Bueler_R0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_R0)
    END IF
    IF ( C%refgeo_idealised_Bueler_lambda < -10._dp .OR. C%refgeo_idealised_Bueler_lambda > 10._dp) THEN
      CALL crash('refgeo_idealised_Bueler_lambda has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_lambda)
    END IF

    CALL Bueler_dome( C%uniform_flow_factor, C%n_flow, C%refgeo_idealised_Bueler_H0, C%refgeo_idealised_Bueler_R0, &
      C%refgeo_idealised_Bueler_lambda, x, y, 0._dp, Hi, M)
    Hb = 0._dp
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_Bueler

  SUBROUTINE calc_idealised_geometry_SSA_icestream( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! A 2,000 m thick slab of ice on a flat, inclined plane (10 m per km slope)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_SSA_icestream'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_SSA_icestream_Hi < 100._dp .OR. &
         C%refgeo_idealised_SSA_icestream_Hi > 10000._dp) THEN
      CALL crash('refgeo_idealised_SSA_icestream_Hi has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_SSA_icestream_Hi)
    END IF
    IF ( ABS( C%refgeo_idealised_SSA_icestream_dhdx) > 0.1_dp) THEN
      CALL crash('refgeo_idealised_SSA_icestream_dhdx has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_SSA_icestream_dhdx)
    END IF

    Hi = C%refgeo_idealised_SSA_icestream_Hi
    Hb = C%refgeo_idealised_SSA_icestream_dhdx * x
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_SSA_icestream

  SUBROUTINE calc_idealised_geometry_MISMIP_mod( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! The MISMIP_mod cone-shaped island

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_MISMIP_mod'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_MISMIP_mod_Hi_init < 0._dp .OR. &
         C%refgeo_idealised_MISMIP_mod_Hi_init > 10000._dp) THEN
      CALL crash('refgeo_idealised_MISMIP_mod_Hi_init has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_MISMIP_mod_Hi_init)
    END IF

    Hi = C%refgeo_idealised_MISMIP_mod_Hi_init
    Hb = 720._dp - 778.5_dp * SQRT( x**2 + y**2)/ 750000._dp
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_MISMIP_mod

  SUBROUTINE calc_idealised_geometry_ISMIP_HOM_A( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! ISMIP-HOM Experiment A (slab on a bumpy slope in both directions)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_ISMIP_HOM_A'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) THEN
      CALL crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    END IF

    Hs = 2000._dp - x * TAN( 0.5_dp * pi / 180._dp)
    Hb = Hs - 1000._dp + 500._dp * SIN( x * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L) * SIN( y * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L)
    Hi = Hs - Hb
    SL = -10000._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_ISMIP_HOM_A

  SUBROUTINE calc_idealised_geometry_ISMIP_HOM_B( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! ISMIP-HOM Experiment B (slab on a bumpy slope in only the x-directions)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_ISMIP_HOM_B'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) THEN
      CALL crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    END IF

    Hs = 2000._dp - x * TAN( 0.5_dp * pi / 180._dp)
    Hb = Hs - 1000._dp + 500._dp * SIN( x * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L)
    Hi = Hs - Hb
    SL = -10000._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_ISMIP_HOM_B

  SUBROUTINE calc_idealised_geometry_ISMIP_HOM_CD( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! ISMIP-HOM Experiment C/D (slab on a slope)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_ISMIP_HOM_CD'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) THEN
      CALL crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    END IF

    Hs = 2000._dp - x * TAN( 0.5_dp * pi / 180._dp)
    Hb = Hs - 1000._dp
    Hi = Hs - Hb
    SL = -10000._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_ISMIP_HOM_CD

  SUBROUTINE calc_idealised_geometry_ISMIP_HOM_E( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! ISMIP-HOM Experiment E (Haut Glacier d'Arolla)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_ISMIP_HOM_E'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    CALL crash('ISMIP-HOM E is not implemented in UFEMISM!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_ISMIP_HOM_E

  SUBROUTINE calc_idealised_geometry_ISMIP_HOM_F( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! ISMIP-HOM Experiment F (slab on another bumpy slope)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_ISMIP_HOM_F'
    REAL(dp), PARAMETER                           :: H0    = 1000._dp
    REAL(dp), PARAMETER                           :: a0    = 100._dp
    REAL(dp), PARAMETER                           :: sigma = 10000._dp

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) THEN
      CALL crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    END IF

    Hs = 5000._dp - x * TAN( 3._dp * pi / 180._dp)
    Hb = Hs - H0 + a0 * EXP( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * EXP( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2)
    Hi = Hs - Hb
    SL = -10000._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_ISMIP_HOM_F

  SUBROUTINE calc_idealised_geometry_MISMIPplus( x, y, Hi, Hb, Hs, SL)
    ! Calculate an idealised geometry
    !
    ! The MISMIpplus fjord geometry

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                       INTENT(IN)    :: x,y             ! [m] Coordinates
    REAL(dp),                       INTENT(OUT)   :: Hi              ! [m] Ice thickness
    REAL(dp),                       INTENT(OUT)   :: Hb              ! [m] Bedrock elevation
    REAL(dp),                       INTENT(OUT)   :: Hs              ! [m] Surface elevation
    REAL(dp),                       INTENT(OUT)   :: SL              ! [m] Sea level

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_idealised_geometry_MISMIPplus'
    REAL(dp)                                      :: xtilde,Bx,By
    REAL(dp), PARAMETER                           :: B0     = -150._dp
    REAL(dp), PARAMETER                           :: B2     = -728.8_dp
    REAL(dp), PARAMETER                           :: B4     = 343.91_dp
    REAL(dp), PARAMETER                           :: B6     = -50.57_dp
    REAL(dp), PARAMETER                           :: xbar   = 300000._dp
    REAL(dp), PARAMETER                           :: fc     = 4000._dp
    REAL(dp), PARAMETER                           :: dc     = 500._dp
    REAL(dp), PARAMETER                           :: wc     = 24000._dp
    REAL(dp), PARAMETER                           :: zbdeep = -720._dp

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety
    IF ( C%refgeo_idealised_MISMIPplus_Hi_init < 0._dp .OR. &
         C%refgeo_idealised_MISMIPplus_Hi_init > 10000._dp) THEN
      CALL crash('refgeo_idealised_MISMIPplus_Hi_init has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_MISMIPplus_Hi_init)
    END IF

    xtilde = x / xbar
    Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
    By = (dc / (1 + EXP(-2._dp*(y - wc)/fc))) + &
         (dc / (1 + EXP( 2._dp*(y + wc)/fc)))

    Hi = C%refgeo_idealised_MISMIPplus_Hi_init
    Hb = MAX( Bx + By, zbdeep)
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_idealised_geometry_MISMIPplus

END MODULE reference_geometries
