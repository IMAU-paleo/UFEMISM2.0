MODULE netcdf_input_module

! ===== Flexible reading of input files =====
! ===========================================
!
! These routines allow for flexible reading of input files. The three top-level functions
! allow you to read 2-D, 2-D monthly, and 3-D data fields from NetCDF files. The files can
! contain the data on a regular x/y-grid, a regular lon/lat-grid, or a mesh. The routines
! will automatically detect which one of these it is, and select the appropriate
! subroutine. The data will also be automatically mapped to the provided model mesh.

! ===== Preamble =====
! ====================

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
                                             allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                             allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                             allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                             allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                             allocate_shared_int_4D,   allocate_shared_dp_4D, &
                                             allocate_shared_bool_0D,  allocate_shared_bool_1D, &
                                             reallocate_shared_int_0D, reallocate_shared_dp_0D, &
                                             reallocate_shared_int_1D, reallocate_shared_dp_1D, &
                                             reallocate_shared_int_2D, reallocate_shared_dp_2D, &
                                             reallocate_shared_int_3D, reallocate_shared_dp_3D, &
                                             deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_grid, type_grid_lonlat, type_model_region
  USE netcdf,                          ONLY: NF90_MAX_VAR_DIMS
  USE mesh_memory_module,              ONLY: allocate_mesh_primary, allocate_mesh_secondary, deallocate_mesh_all
  USE mesh_help_functions_module,      ONLY: calc_triangle_geometric_centres, find_Voronoi_cell_areas, calc_lat_lon_coordinates, &
                                             find_triangle_areas, find_connection_widths, determine_mesh_resolution, check_mesh, &
                                             find_Voronoi_cell_geometric_centres
  USE mesh_ArakawaC_module,            ONLY: make_Ac_mesh
  USE mesh_operators_module,           ONLY: calc_matrix_operators_mesh_basic
  USE mesh_mapping_module,             ONLY: map_from_xy_grid_to_mesh_2D        , map_from_lonlat_grid_to_mesh_2D        , map_from_mesh_to_mesh_2D, &
                                             map_from_xy_grid_to_mesh_2D_monthly, map_from_lonlat_grid_to_mesh_2D_monthly, map_from_mesh_to_mesh_2D_monthly, &
                                             map_from_xy_grid_to_mesh_3D        , map_from_lonlat_grid_to_mesh_3D        , map_from_mesh_to_mesh_3D
  USE utilities_module,                ONLY: flip_1D_dp, flip_2D_x1_dp, flip_2D_x2_dp, flip_3D_x1_dp, flip_3D_x2_dp, flip_3D_x3_dp, &
                                             permute_2D_dp, permute_3D_dp, permute_2D_int, permute_3D_int, inverse_oblique_sg_projection, &
                                             deallocate_grid, deallocate_grid_lonlat, remap_zeta_grid_dp, remap_zeta_mesh_dp
  USE netcdf_basic_module,             ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, field_name_options_z_ocean, &
                                             field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                                             field_name_options_dim_nV, field_name_options_dim_nTri, field_name_options_dim_nC_mem, &
                                             field_name_options_dim_nAc, field_name_options_dim_two, field_name_options_dim_three, &
                                             field_name_options_dim_six, field_name_options_V, field_name_options_Tri, field_name_options_nC, &
                                             field_name_options_C, field_name_options_niTri, field_name_options_iTri, &
                                             field_name_options_edge_index, field_name_options_Tricc, field_name_options_TriC, &
                                             field_name_options_Tri_edge_index, field_name_options_VAc, field_name_options_Aci, &
                                             field_name_options_iAci, field_name_options_A, field_name_options_R, &
                                             field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                                             field_name_options_SL, field_name_options_Ti, &
                                             open_existing_netcdf_file_for_reading, close_netcdf_file, &
                                             inquire_dim_multiple_options, inquire_var_multiple_options, &
                                             read_var_int_0D, read_var_int_1D, read_var_int_2D, read_var_int_3D, read_var_int_4D, &
                                             read_var_dp_0D , read_var_dp_1D , read_var_dp_2D , read_var_dp_3D , read_var_dp_4D, &
                                             check_x, check_y, check_lon, check_lat, check_mesh_dimensions, check_zeta, check_z_ocean, find_timeframe, &
                                             check_xy_grid_field_int_2D, check_xy_grid_field_dp_2D, check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, &
                                             check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, check_lonlat_grid_field_dp_2D_monthly, check_lonlat_grid_field_dp_3D, &
                                             check_mesh_field_int_2D, check_mesh_field_dp_2D, check_mesh_field_dp_2D_monthly, check_mesh_field_dp_3D, &
                                             check_xy_grid_field_dp_3D_ocean, check_lonlat_grid_field_dp_3D_ocean, check_mesh_field_dp_3D_ocean, &
                                             inquire_xy_grid, inquire_lonlat_grid, inquire_mesh

  IMPLICIT NONE

CONTAINS

  ! ===== Top-level functions =====
  ! ===============================

  ! Read and map to mesh
  SUBROUTINE read_field_from_file_2D(         filename, field_name_options, mesh, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_2D'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid, has_mesh, has_time
    TYPE(type_grid)                                    :: grid_from_file
    TYPE(type_grid_lonlat)                             :: grid_lonlat_from_file
    TYPE(type_mesh)                                    :: mesh_from_file
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid_from_file
    INTEGER                                            :: wd_grid_from_file
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid_lonlat_from_file
    INTEGER                                            :: wd_grid_lonlat_from_file
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_mesh_from_file
    INTEGER                                            :: wd_mesh_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)
    CALL inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    IF (has_xy_grid     .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a mesh!')
    IF (has_lonlat_grid .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      ! Data is provided on an x/y-grid

      ! Read grid and gridded data
      CALL read_field_from_xy_file_2D( filename, field_name_options, region_name, grid_from_file, d_grid_from_file, wd_grid_from_file, time_to_read)

      ! Remap data
      CALL map_from_xy_grid_to_mesh_2D( grid_from_file, mesh, d_grid_from_file, d)

      ! Clean up after yourself
      CALL deallocate_grid(      grid_from_file)
      CALL deallocate_shared( wd_grid_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_2D( filename, field_name_options, region_name, grid_lonlat_from_file, d_grid_lonlat_from_file, wd_grid_lonlat_from_file, time_to_read)

      ! Remap data
      CALL map_from_lonlat_grid_to_mesh_2D( grid_lonlat_from_file, mesh, d_grid_lonlat_from_file, d)

      ! Clean up after yourself
      CALL deallocate_grid_lonlat(    grid_lonlat_from_file)
      CALL deallocate_shared(      wd_grid_lonlat_from_file)

    ELSEIF (has_mesh) THEN
      ! Data is provided on a mesh

      ! Read grid and gridded data
      CALL read_field_from_mesh_file_2D( filename, field_name_options, region_name, mesh_from_file, d_mesh_from_file, wd_mesh_from_file, time_to_read)

      ! Remap data
      CALL map_from_mesh_to_mesh_2D( mesh_from_file, mesh, d_mesh_from_file, d)

      ! Clean up after yourself
      CALL deallocate_mesh_all(    mesh_from_file)
      CALL deallocate_shared(   wd_mesh_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_2D

  SUBROUTINE read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_2D_monthly'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid, has_mesh, has_time
    TYPE(type_grid)                                    :: grid_from_file
    TYPE(type_grid_lonlat)                             :: grid_lonlat_from_file
    TYPE(type_mesh)                                    :: mesh_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_from_file
    INTEGER                                            :: wd_grid_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_lonlat_from_file
    INTEGER                                            :: wd_grid_lonlat_from_file
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_mesh_from_file
    INTEGER                                            :: wd_mesh_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)
    CALL inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    IF (has_xy_grid     .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a mesh!')
    IF (has_lonlat_grid .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      ! Data is provided on an x/y-grid

      ! Read grid and gridded data
      CALL read_field_from_xy_file_2D_monthly( filename, field_name_options, region_name, grid_from_file, d_grid_from_file, wd_grid_from_file, time_to_read)

      ! Remap data
      CALL map_from_xy_grid_to_mesh_2D_monthly( grid_from_file, mesh, d_grid_from_file, d)

      ! Clean up after yourself
      CALL deallocate_grid(      grid_from_file)
      CALL deallocate_shared( wd_grid_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_2D_monthly( filename, field_name_options, region_name, grid_lonlat_from_file, d_grid_lonlat_from_file, wd_grid_lonlat_from_file, time_to_read)

      ! Remap data
      CALL map_from_lonlat_grid_to_mesh_2D_monthly( grid_lonlat_from_file, mesh, d_grid_lonlat_from_file, d)

      ! Clean up after yourself
      CALL deallocate_grid_lonlat(    grid_lonlat_from_file)
      CALL deallocate_shared(      wd_grid_lonlat_from_file)

    ELSEIF (has_mesh) THEN
      ! Data is provided on a mesh

      ! Read grid and gridded data
      CALL read_field_from_mesh_file_2D_monthly( filename, field_name_options, region_name, mesh_from_file, d_mesh_from_file, wd_mesh_from_file, time_to_read)

      ! Remap data
      CALL map_from_mesh_to_mesh_2D_monthly( mesh_from_file, mesh, d_mesh_from_file, d)

      ! Clean up after yourself
      CALL deallocate_mesh_all(    mesh_from_file)
      CALL deallocate_shared(   wd_mesh_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_2D_monthly

  SUBROUTINE read_field_from_file_3D(         filename, field_name_options, mesh, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.
    !
    ! NOTE: only meant to be used for 3-D englacial data, not for 3-D ocean data!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_3D'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid, has_mesh, has_time
    TYPE(type_grid)                                    :: grid_from_file
    TYPE(type_grid_lonlat)                             :: grid_lonlat_from_file
    TYPE(type_mesh)                                    :: mesh_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_from_file
    INTEGER                                            :: wd_grid_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_lonlat_from_file
    INTEGER                                            :: wd_grid_lonlat_from_file
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_mesh_from_file
    INTEGER                                            :: wd_mesh_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)
    CALL inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    IF (has_xy_grid     .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a mesh!')
    IF (has_lonlat_grid .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      ! Data is provided on an x/y-grid

      ! Read grid and gridded data
      CALL read_field_from_xy_file_3D( filename, field_name_options, region_name, grid_from_file, d_grid_from_file, wd_grid_from_file, time_to_read)

      ! Remap data
      CALL map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, d_grid_from_file, d)

      ! Clean up after yourself
      CALL deallocate_grid(      grid_from_file)
      CALL deallocate_shared( wd_grid_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_3D( filename, field_name_options, region_name, grid_lonlat_from_file, d_grid_lonlat_from_file, wd_grid_lonlat_from_file, time_to_read)

      ! Remap data
      CALL map_from_lonlat_grid_to_mesh_3D( grid_lonlat_from_file, mesh, d_grid_lonlat_from_file, d)

      ! Clean up after yourself
      CALL deallocate_grid_lonlat(    grid_lonlat_from_file)
      CALL deallocate_shared(      wd_grid_lonlat_from_file)

    ELSEIF (has_mesh) THEN
      ! Data is provided on a mesh

      ! Read file mesh and data
      CALL read_field_from_mesh_file_3D( filename, field_name_options, region_name, mesh_from_file, d_mesh_from_file, wd_mesh_from_file, time_to_read)

      ! Remap data
      CALL map_from_mesh_to_mesh_3D( mesh_from_file, mesh, d_mesh_from_file, d)

      ! Clean up after yourself
      CALL deallocate_mesh_all(    mesh_from_file)
      CALL deallocate_shared(   wd_mesh_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_3D

  ! ===== Medium-level functions =====
  ! ==================================

  ! Read a field from an x/y-grid file
  SUBROUTINE read_field_from_xy_file_2D(             filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_2D'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, ncid, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%nx, grid%ny, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_2D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%nx, grid%ny, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:) = d_with_time( grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_2D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%ny, grid%nx, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%ny, grid%nx, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2) = d_with_time( :,grid%i1:grid%i2,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_2D_dp( d, wd, map = [2,1])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_2D_x1_dp( d)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_2D_x2_dp( d)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 19)

  END SUBROUTINE read_field_from_xy_file_2D

  SUBROUTINE read_field_from_xy_file_2D_monthly(     filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D monthly data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_2D_monthly'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, ncid, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nx, grid%ny, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nx, grid%ny, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, 12, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%ny, grid%nx, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%ny, grid%nx, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%ny, grid%nx, 12, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_3D_dp( d, wd, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_3D_x1_dp( d)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_3D_x2_dp( d)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 19)

  END SUBROUTINE read_field_from_xy_file_2D_monthly

  SUBROUTINE read_field_from_xy_file_3D(             filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 3-D data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_3D'
    INTEGER                                            :: ncid
    INTEGER                                            :: nzeta_from_file
    REAL(dp), DIMENSION(:    ), POINTER                ::  zeta_from_file
    INTEGER                                            :: wzeta_from_file
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_zeta_from_file
    INTEGER                                            :: wd_zeta_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, ncid, grid, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_zeta_from_file( filename, ncid, nzeta_from_file, zeta_from_file, wzeta_from_file)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nx, grid%ny, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nx, grid%ny, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%ny, grid%nx, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%ny, grid%nx, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%ny, grid%nx, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_3D_dp( d_zeta_from_file, wd_zeta_from_file, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_3D_x1_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_3D_x2_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, C%nz, d, wd)

    ! Remap to the model vertical grid
    CALL remap_zeta_grid_dp( zeta_from_file, d_zeta_from_file, C%zeta, d)

    ! Clean up after yourself
    CALL deallocate_shared(   wzeta_from_file)
    CALL deallocate_shared( wd_zeta_from_file)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 19)

  END SUBROUTINE read_field_from_xy_file_3D

  SUBROUTINE read_field_from_xy_file_3D_ocean(       filename, field_name_options, region_name, grid, d, wd, nz_ocean, z_ocean, wz_ocean, time_to_read)
    ! Read a 3-D ocean data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read
    INTEGER,                             INTENT(OUT)   :: nz_ocean
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: z_ocean
    INTEGER,                             INTENT(OUT)   :: wz_ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_3D_ocean'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, ncid, grid, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_z_ocean_from_file( filename, ncid, nz_ocean, z_ocean, wz_ocean)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz_ocean, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nx, grid%ny, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nz_ocean, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%ny, grid%nx, nz_ocean, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%ny, grid%nx, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%ny, grid%nx, nz_ocean, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_3D_dp( d, wd, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_3D_x1_dp( d)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_3D_x2_dp( d)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 19)

  END SUBROUTINE read_field_from_xy_file_3D_ocean

  ! Read a field from a lon/lat-grid file
  SUBROUTINE read_field_from_lonlat_file_2D(         filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_2D'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, ncid, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%nlon, grid%nlat, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_2D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:) = d_with_time( grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%nlat, grid%nlon, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_2D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%nlat, grid%nlon, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2) = d_with_time( :,grid%i1:grid%i2,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_2D_dp( d, wd, map = [2,1])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_2D_x1_dp( d)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_2D_x2_dp( d)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_2D( filename, grid, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_lonlat_file_2D

  SUBROUTINE read_field_from_lonlat_file_2D_monthly( filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D monthly data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_2D'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, ncid, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlon, grid%nlat, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 12, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlat, grid%nlon, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlat, grid%nlon, 12, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_3D_dp( d, wd, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_3D_x1_dp( d)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_3D_x2_dp( d)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_3D( filename, grid, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_lonlat_file_2D_monthly

  SUBROUTINE read_field_from_lonlat_file_3D(         filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 3-D data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_3D'
    INTEGER                                            :: ncid
    INTEGER                                            :: nzeta_from_file
    REAL(dp), DIMENSION(:    ), POINTER                ::  zeta_from_file
    INTEGER                                            :: wzeta_from_file
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_zeta_from_file
    INTEGER                                            :: wd_zeta_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, ncid, grid, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_zeta_from_file( filename, ncid, nzeta_from_file, zeta_from_file, wzeta_from_file)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlon, grid%nlat, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlat, grid%nlon, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlat, grid%nlon, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_3D_dp( d_zeta_from_file, wd_zeta_from_file, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_3D_x1_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_3D_x2_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_3D( filename, grid, d_zeta_from_file)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, C%nz, d, wd)

    ! Remap to the model vertical grid
    CALL remap_zeta_grid_dp( zeta_from_file, d_zeta_from_file, C%zeta, d)

    ! Clean up after yourself
    CALL deallocate_shared(   wzeta_from_file)
    CALL deallocate_shared( wd_zeta_from_file)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_lonlat_file_3D

  SUBROUTINE read_field_from_lonlat_file_3D_ocean(   filename, field_name_options, region_name, grid, d, wd, nz_ocean, z_ocean, wz_ocean, time_to_read)
    ! Read a 3-D ocean data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read
    INTEGER,                             INTENT(OUT)   :: nz_ocean
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: z_ocean
    INTEGER,                             INTENT(OUT)   :: wz_ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_3D_ocean'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, ncid, grid, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_z_ocean_from_file( filename, ncid, nz_ocean, z_ocean, wz_ocean)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, nz_ocean, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlon, grid%nlat, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, nz_ocean, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, nz_ocean, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, ncid, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlat, grid%nlon, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlat, grid%nlon, nz_ocean, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_3D_dp( d, wd, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_3D_x1_dp( d)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_3D_x2_dp( d)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_3D( filename, grid, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_lonlat_file_3D_ocean

  ! Read a field from a mesh file
  SUBROUTINE read_field_from_mesh_file_2D(           filename, field_name_options, region_name, mesh, d, wd, time_to_read)
    ! Read a 2-D data field from a NetCDF file on a mesh,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_mesh_file_2D'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_mesh_from_file( filename, ncid, mesh, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, d, wd)

    ! Read data from file
    IF (.NOT. PRESENT( time_to_read)) THEN
      CALL read_var_dp_1D( filename, ncid, id_var, d)
    ELSE
      ! Allocate shared memory
      CALL allocate_shared_dp_2D( mesh%nV, 1, d_with_time, wd_with_time)
      ! Find out which timeframe to read
      CALL find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      CALL read_var_dp_2D( filename, ncid, id_var, d_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )
      ! Copy to output memory
      d( mesh%vi1:mesh%vi2) = d_with_time( mesh%vi1:mesh%vi2,1)
      ! Clean up after yourself
      CALL deallocate_shared( wd_with_time)
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_mesh_file_2D

  SUBROUTINE read_field_from_mesh_file_2D_monthly(   filename, field_name_options, region_name, mesh, d, wd, time_to_read)
    ! Read a 2-D monthly data field from a NetCDF file on a mesh,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_mesh_file_2D_monthly'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_mesh_from_file( filename, ncid, mesh, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, 12, d, wd)

    ! Read data from file
    IF (.NOT. PRESENT( time_to_read)) THEN
      CALL read_var_dp_2D( filename, ncid, id_var, d)
    ELSE
      ! Allocate shared memory
      CALL allocate_shared_dp_3D( mesh%nV, 12, 1, d_with_time, wd_with_time)
      ! Find out which timeframe to read
      CALL find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      CALL read_var_dp_3D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, 12, 1 /) )
      ! Copy to output memory
      d( mesh%vi1:mesh%vi2,:) = d_with_time( mesh%vi1:mesh%vi2,:,1)
      ! Clean up after yourself
      CALL deallocate_shared( wd_with_time)
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_mesh_file_2D_monthly

  SUBROUTINE read_field_from_mesh_file_3D(           filename, field_name_options, region_name, mesh, d, wd, time_to_read)
    ! Read a 3-D data field from a NetCDF file on a mesh,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_mesh_file_3D'
    INTEGER                                            :: ncid
    INTEGER                                            :: nzeta_from_file
    REAL(dp), DIMENSION(:    ), POINTER                ::  zeta_from_file
    INTEGER                                            :: wzeta_from_file
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_zeta_from_file
    INTEGER                                            :: wd_zeta_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_mesh_from_file( filename, ncid, mesh, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_zeta_from_file( filename, ncid, nzeta_from_file, zeta_from_file, wzeta_from_file)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

    ! Read data from file
    IF (.NOT. PRESENT( time_to_read)) THEN
      CALL read_var_dp_2D( filename, ncid, id_var, d_zeta_from_file)
    ELSE
      ! Allocate shared memory
      CALL allocate_shared_dp_3D( mesh%nV, nzeta_from_file, 1, d_with_time, wd_with_time)
      ! Find out which timeframe to read
      CALL find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      CALL read_var_dp_3D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, nzeta_from_file, 1 /) )
      ! Copy to output memory
      d_zeta_from_file( mesh%vi1:mesh%vi2,:) = d_with_time( mesh%vi1:mesh%vi2,:,1)
      ! Clean up after yourself
      CALL deallocate_shared( wd_with_time)
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, d, wd)

    ! Remap to the model vertical grid
    CALL remap_zeta_mesh_dp( zeta_from_file, d_zeta_from_file, C%zeta, d)

    ! Clean up after yourself
    CALL deallocate_shared(   wzeta_from_file)
    CALL deallocate_shared( wd_zeta_from_file)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_mesh_file_3D

  SUBROUTINE read_field_from_mesh_file_3D_ocean(     filename, field_name_options, region_name, mesh, d, wd, nz_ocean, z_ocean, wz_ocean, time_to_read)
    ! Read a 3-D ocean data field from a NetCDF file on a mesh,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read
    INTEGER,                             INTENT(OUT)   :: nz_ocean
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: z_ocean
    INTEGER,                             INTENT(OUT)   :: wz_ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_mesh_file_3D_ocean'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_mesh_from_file( filename, ncid, mesh, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_z_ocean_from_file( filename, ncid, nz_ocean, z_ocean, wz_ocean)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, nz_ocean, d, wd)

    ! Read data from file
    IF (.NOT. PRESENT( time_to_read)) THEN
      CALL read_var_dp_2D( filename, ncid, id_var, d)
    ELSE
      ! Allocate shared memory
      CALL allocate_shared_dp_3D( mesh%nV, nz_ocean, 1, d_with_time, wd_with_time)
      ! Find out which timeframe to read
      CALL find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      CALL read_var_dp_3D( filename, ncid, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, nz_ocean, 1 /) )
      ! Copy to output memory
      d( mesh%vi1:mesh%vi2,:) = d_with_time( mesh%vi1:mesh%vi2,:,1)
      ! Clean up after yourself
      CALL deallocate_shared( wd_with_time)
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE read_field_from_mesh_file_3D_ocean

  ! ===== Set up grids/mesh from a NetCDF file =====
  ! ================================================

  SUBROUTINE setup_xy_grid_from_file(     filename, ncid, grid, region_name)
    ! Set up an x/y-grid from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the grid at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_xy_grid_from_file'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var_x, id_var_y
    INTEGER                                            :: i,j,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check grid dimensions and variables for validity
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Allocate memory for the grid size
    CALL allocate_shared_int_0D( grid%nx, grid%wnx)
    CALL allocate_shared_int_0D( grid%ny, grid%wny)
    CALL allocate_shared_int_0D( grid%n , grid%wn )

    ! Inquire x and y dimensions
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_x, id_dim_x, dim_length = grid%nx)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_y, id_dim_y, dim_length = grid%ny)
    grid%n = grid%nx * grid%ny

    ! Allocate memory for x and y
    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)

    ! Inquire x and y variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_x, id_var_x)
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_y, id_var_y)

    ! Read x and y
    CALL read_var_dp_1D(  filename, ncid, id_var_x, grid%x)
    CALL read_var_dp_1D(  filename, ncid, id_var_y, grid%y)

    ! Allocate memory for, and calculate, some secondary grid properties
    CALL allocate_shared_dp_0D(                    grid%xmin       , grid%wxmin       )
    CALL allocate_shared_dp_0D(                    grid%xmax       , grid%wxmax       )
    CALL allocate_shared_dp_0D(                    grid%ymin       , grid%wymin       )
    CALL allocate_shared_dp_0D(                    grid%ymax       , grid%wymax       )
    CALL allocate_shared_dp_0D(                    grid%dx         , grid%wdx         )
    CALL allocate_shared_dp_0D(                    grid%tol_dist   , grid%wtol_dist   )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n       , grid%wij2n       )
    CALL allocate_shared_int_2D( grid%n , 2,       grid%n2ij       , grid%wn2ij       )
    CALL allocate_shared_dp_0D(                    grid%lambda_m   , grid%wlambda_m   )
    CALL allocate_shared_dp_0D(                    grid%phi_m      , grid%wphi_m      )
    CALL allocate_shared_dp_0D(                    grid%beta_stereo, grid%wbeta_stereo)
    CALL allocate_shared_dp_2D(  grid%nx, grid%ny, grid%lon        , grid%wlon        )
    CALL allocate_shared_dp_2D(  grid%nx, grid%ny, grid%lat        , grid%wlat        )

    ! Calculate secondary grid data
    IF (par%master) THEN

      ! Resolution
      grid%dx   = ABS( grid%x( 2) - grid%x( 1))

      ! Safety
      DO i = 1, grid%nx-1
        IF (1._dp - ABS(grid%x( i+1) - grid%x( i)) / grid%dx > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular x-dimension!')
      END DO
      DO j = 1, grid%ny-1
        IF (1._dp - ABS(grid%y( j+1) - grid%y( j)) / grid%dx > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular y-dimension!')
      END DO

      ! Domain size
      grid%xmin = MINVAL( grid%x)
      grid%xmax = MAXVAL( grid%x)
      grid%ymin = MINVAL( grid%y)
      grid%ymax = MAXVAL( grid%y)

      ! Tolerance; points lying within this distance of each other are treated as identical
      grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

      ! Conversion tables for grid-form vs. vector-form data
      n = 0
      DO i = 1, grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, grid%ny
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = grid%ny, 1, -1
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO

    END IF ! IF (par%master) THEN
    CALL sync

    ! Set up parallelisation domains
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! Projection parameters for this region
    IF     (region_name == 'NAM') THEN
      grid%lambda_M     = C%lambda_M_NAM
      grid%phi_M        = C%phi_M_NAM
      grid%beta_stereo  = C%beta_stereo_NAM
    ELSEIF (region_name == 'EAS') THEN
      grid%lambda_M     = C%lambda_M_EAS
      grid%phi_M        = C%phi_M_EAS
      grid%beta_stereo  = C%beta_stereo_EAS
    ELSEIF (region_name == 'GRL') THEN
      grid%lambda_M     = C%lambda_M_GRL
      grid%phi_M        = C%phi_M_GRL
      grid%beta_stereo  = C%beta_stereo_GRL
    ELSEIF (region_name == 'ANT') THEN
      grid%lambda_M     = C%lambda_M_ANT
      grid%phi_M        = C%phi_M_ANT
      grid%beta_stereo  = C%beta_stereo_ANT
    END IF

    ! Lon/lat coordinates
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), grid%lambda_M, grid%phi_M, grid%beta_stereo, grid%lon( i,j), grid%lat( i,j))
    END DO
    END DO
    CALL sync

    ! Give the grid a nice name
    grid%name = 'xy_grid_from_file_"' // TRIM( filename) // '"'

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 18)

  END SUBROUTINE setup_xy_grid_from_file

  SUBROUTINE setup_lonlat_grid_from_file( filename, ncid, grid, region_name)
    ! Set up a lon/lat-grid from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the grid at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_lonlat_grid_from_file'
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: id_var_lon, id_var_lat
    INTEGER                                            :: i,j,n
    REAL(dp)                                           :: dlon

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check grid dimensions and variables for validity
    CALL check_lon( filename, ncid)
    CALL check_lat( filename, ncid)

    ! Allocate memory for the grid size
    CALL allocate_shared_int_0D( grid%nlon, grid%wnlon)
    CALL allocate_shared_int_0D( grid%nlat, grid%wnlat)
    CALL allocate_shared_int_0D( grid%n   , grid%wn   )

    ! Inquire lon and lat dimensions
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = grid%nlon)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = grid%nlat)

    IF (par%master) THEN
      grid%n = grid%nlon * grid%nlat
    END IF

    ! Allocate memory for lon and lat
    CALL allocate_shared_dp_1D( grid%nlon, grid%lon, grid%wlon)
    CALL allocate_shared_dp_1D( grid%nlat, grid%lat, grid%wlat)

    ! Inquire lon and lat variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lon, id_var_lon)
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read lon and lat
    CALL read_var_dp_1D(  filename, ncid, id_var_lon, grid%lon)
    CALL read_var_dp_1D(  filename, ncid, id_var_lat, grid%lat)

    ! Allocate memory for, and calculate, some secondary grid properties
    CALL allocate_shared_dp_0D( grid%lonmin     , grid%wlonmin     )
    CALL allocate_shared_dp_0D( grid%lonmax     , grid%wlonmax     )
    CALL allocate_shared_dp_0D( grid%latmin     , grid%wlatmin     )
    CALL allocate_shared_dp_0D( grid%latmax     , grid%wlatmax     )
    CALL allocate_shared_dp_0D( grid%dlon       , grid%wdlon       )
    CALL allocate_shared_dp_0D( grid%dlat       , grid%wdlat       )
    CALL allocate_shared_dp_0D( grid%lambda_m   , grid%wlambda_m   )
    CALL allocate_shared_dp_0D( grid%phi_m      , grid%wphi_m      )
    CALL allocate_shared_dp_0D( grid%beta_stereo, grid%wbeta_stereo)

    CALL allocate_shared_int_2D( grid%nlon, grid%nlat, grid%ij2n, grid%wij2n)
    CALL allocate_shared_int_2D( grid%n   , 2,         grid%n2ij, grid%wn2ij)

    ! Calculate secondary grid data
    IF (par%master) THEN

      ! Resolution
      grid%dlon = ABS( grid%lon( 2) - grid%lon( 1))
      grid%dlat = ABS( grid%lat( 2) - grid%lat( 1))

      ! Safety
      DO i = 1, grid%nlon - 1
        ! Check for regularity in longitude, but allow for 360-degree jumps
        dlon = MIN( ABS( grid%lon( i+1) - grid%lon( i)), ABS( grid%lon( i+1) + 360._dp - grid%lon( i)))
        IF (ABS( 1._dp - dlon / grid%dlon) > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular longitude dimension!')
      END DO
      DO j = 1, grid%nlat - 1
        IF (ABS( 1._dp - ABS( grid%lat( j+1) - grid%lat( j)) / grid%dlat) > 1E-6_dp) &
          CALL crash('file "' // TRIM( filename) // '" has an irregular latitude dimension!')
      END DO

      ! Domain size
      grid%lonmin = MINVAL( grid%lon)
      grid%lonmax = MAXVAL( grid%lon)
      grid%latmin = MINVAL( grid%lat)
      grid%latmax = MAXVAL( grid%lat)

      ! Conversion tables for grid-form vs. vector-form data
      n = 0
      DO i = 1, grid%nlon
        IF (MOD(i,2) == 1) THEN
          DO j = 1, grid%nlat
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = grid%nlat, 1, -1
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO

    END IF ! IF (par%master) THEN
    CALL sync

    ! Projection parameters for this region
    IF     (region_name == 'NAM') THEN
      grid%lambda_M     = C%lambda_M_NAM
      grid%phi_M        = C%phi_M_NAM
      grid%beta_stereo  = C%beta_stereo_NAM
    ELSEIF (region_name == 'EAS') THEN
      grid%lambda_M     = C%lambda_M_EAS
      grid%phi_M        = C%phi_M_EAS
      grid%beta_stereo  = C%beta_stereo_EAS
    ELSEIF (region_name == 'GRL') THEN
      grid%lambda_M     = C%lambda_M_GRL
      grid%phi_M        = C%phi_M_GRL
      grid%beta_stereo  = C%beta_stereo_GRL
    ELSEIF (region_name == 'ANT') THEN
      grid%lambda_M     = C%lambda_M_ANT
      grid%phi_M        = C%phi_M_ANT
      grid%beta_stereo  = C%beta_stereo_ANT
    END IF

    ! Set up parallelisation domains
    CALL partition_list( grid%nlon, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%nlat, par%i, par%n, grid%j1, grid%j2)

    ! Give the grid a nice name
    grid%name = 'lonlat_grid_from_file_"' // TRIM( filename) // '"'

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 10)

  END SUBROUTINE setup_lonlat_grid_from_file

  SUBROUTINE setup_mesh_from_file(        filename, ncid, mesh, region_name)
    ! Set up a mesh from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the mesh at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_mesh_from_file'
    INTEGER                                            :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    INTEGER                                            :: nV_mem, nTri_mem, nC_mem, n_two, n_three
    INTEGER                                            :: id_var_V, id_var_nC, id_var_C, id_var_niTri, id_var_iTri, id_var_edge_index
    INTEGER                                            :: id_var_Tri, id_var_Tricc, id_var_TriC, id_var_Tri_edge_index
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check mesh dimensions and variables for validity
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_nV    , id_dim_vi   , dim_length = nV_mem  )
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   , dim_length = nTri_mem)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   , dim_length = nC_mem  )
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_two   , id_dim_two  , dim_length = n_two   )
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_three , id_dim_three, dim_length = n_three )

    ! Allocate memory for the mesh
    CALL allocate_mesh_primary( mesh, region_name, nV_mem, nTri_mem, nC_mem)

    ! Inquire mesh variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_V             , id_var_V             )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_nC            , id_var_nC            )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_C             , id_var_C             )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_niTri         , id_var_niTri         )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_iTri          , id_var_iTri          )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_edge_index    , id_var_edge_index    )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_Tri           , id_var_Tri           )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_Tricc         , id_var_Tricc         )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_TriC          , id_var_TriC          )
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_Tri_edge_index, id_var_Tri_edge_index)

    ! Read mesh data
    CALL read_var_dp_2D(  filename, ncid, id_var_V             , mesh%V             )
    CALL read_var_int_1D( filename, ncid, id_var_nC            , mesh%nC            )
    CALL read_var_int_2D( filename, ncid, id_var_C             , mesh%C             )
    CALL read_var_int_1D( filename, ncid, id_var_niTri         , mesh%niTri         )
    CALL read_var_int_2D( filename, ncid, id_var_iTri          , mesh%iTri          )
    CALL read_var_int_1D( filename, ncid, id_var_edge_index    , mesh%edge_index    )
    CALL read_var_int_2D( filename, ncid, id_var_Tri           , mesh%Tri           )
    CALL read_var_dp_2D(  filename, ncid, id_var_Tricc         , mesh%Tricc         )
    CALL read_var_int_2D( filename, ncid, id_var_TriC          , mesh%TriC          )
    CALL read_var_int_1D( filename, ncid, id_var_Tri_edge_index, mesh%Tri_edge_index)

    ! Calculate secondary mesh data

    IF (par%master) THEN
      mesh%nV       = mesh%nV_mem
      mesh%nTri     = mesh%nTri_mem
      mesh%xmin     = MINVAL( mesh%V( :,1))
      mesh%xmax     = MAXVAL( mesh%V( :,1))
      mesh%ymin     = MINVAL( mesh%V( :,2))
      mesh%ymax     = MAXVAL( mesh%V( :,2))
      mesh%tol_dist = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) * tol / 2._dp
    END IF
    CALL sync

    ! Determine vertex and triangle domains
    CALL partition_list( mesh%nV,   par%i, par%n, mesh%vi1, mesh%vi2)
    CALL partition_list( mesh%nTri, par%i, par%n, mesh%ti1, mesh%ti2)

    ! Calculate extra mesh data
    CALL allocate_mesh_secondary(             mesh)    ! Adds  9 MPI windows
    CALL calc_triangle_geometric_centres(     mesh)
    CALL find_Voronoi_cell_areas(             mesh)
    CALL calc_lat_lon_coordinates(            mesh)
    CALL find_triangle_areas(                 mesh)
    CALL find_connection_widths(              mesh)
    CALL make_Ac_mesh(                        mesh)    ! Adds  5 MPI windows
    CALL calc_matrix_operators_mesh_basic(    mesh)    ! Adds 42 MPI windows (6 CSR matrices, 7 windows each)
    CALL determine_mesh_resolution(           mesh)
    CALL find_Voronoi_cell_geometric_centres( mesh)

    CALL check_mesh( mesh)

    ! Give the mesh a nice name
    mesh%name = 'mesh_from_file_"' // TRIM( filename) // '"'

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 107)

  END SUBROUTINE setup_mesh_from_file

  SUBROUTINE setup_zeta_from_file(        filename, ncid, nzeta, zeta, wzeta)
    ! Set up a zeta coordinate from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(OUT)   :: nzeta
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   ::  zeta
    INTEGER,                             INTENT(OUT)   :: wzeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_mesh_from_file'
    INTEGER                                            :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check zeta dimension and variable for validity
    CALL check_zeta( filename, ncid)

    ! Inquire zeta dimension
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_zeta, id_dim_zeta, dim_length = nzeta)

    ! Inquire zeta variable
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_zeta, id_var_zeta)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( nzeta, zeta, wzeta)

    ! Read zeta from file
    CALL read_var_dp_1D( filename, ncid, id_var_zeta, zeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE setup_zeta_from_file

  SUBROUTINE setup_z_ocean_from_file(     filename, ncid, nz_ocean, z_ocean, wz_ocean)
    ! Set up a z_ocean coordinate from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(OUT)   :: nz_ocean
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   ::  z_ocean
    INTEGER,                             INTENT(OUT)   :: wz_ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_z_ocean_from_file'
    INTEGER                                            :: id_dim_z_ocean, id_var_z_ocean

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check z_ocean dimension and variable for validity
    CALL check_z_ocean( filename, ncid)

    ! Inquire z_ocean dimension
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_z_ocean, id_dim_z_ocean, dim_length = nz_ocean)

    ! Inquire z_ocean variable
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_z_ocean, id_var_z_ocean)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( nz_ocean, z_ocean, wz_ocean)

    ! Read z_ocean from file
    CALL read_var_dp_1D( filename, ncid, id_var_z_ocean, z_ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE setup_z_ocean_from_file

  ! ===== Determine indexing and dimension directions =====
  ! =======================================================

  SUBROUTINE determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)
    ! Determine the indexing and dimension directions of a variable in an x/y-grid file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=256),                  INTENT(OUT)   :: indexing
    CHARACTER(LEN=256),                  INTENT(OUT)   :: xdir
    CHARACTER(LEN=256),                  INTENT(OUT)   :: ydir

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_xy_indexing'
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: nx, ny
    REAL(dp), DIMENSION(:    ), POINTER                :: x, y
    INTEGER                                            :: wx, wy
    INTEGER                                            :: id_var_x, id_var_y, id_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the x and y dimensions and variables of this file are valid
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Inquire x and y dimensions
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_x, id_dim_x, dim_length = nx)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_y, id_dim_y, dim_length = ny)

    ! Allocate memory for x and y
    CALL allocate_shared_dp_1D( nx, x, wx)
    CALL allocate_shared_dp_1D( ny, y, wy)

    ! Inquire x and y variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_x, id_var_x)
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_y, id_var_y)

    ! Read x and y
    CALL read_var_dp_1D(  filename, ncid, id_var_x, x)
    CALL read_var_dp_1D(  filename, ncid, id_var_y, y)

    ! Determine directions of x and y
    IF (x( 2) > x( 1)) THEN
      xdir = 'normal'
    ELSE
      xdir = 'reverse'
    END IF
    IF (y( 2) > y( 1)) THEN
      ydir = 'normal'
    ELSE
      ydir = 'reverse'
    END IF

    ! Inquire dimensions of the specified field variable
    CALL inquire_var_multiple_options( filename, ncid, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    IF     (dims_of_var( 1) == id_dim_x .AND. dims_of_var( 2) == id_dim_y) THEN
      indexing = 'xy'
    ELSEIF (dims_of_var( 1) == id_dim_y .AND. dims_of_var( 2) == id_dim_x) THEN
      indexing = 'yx'
    ELSE
      CALL crash('x and y are not the first two dimensions of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wx)
    CALL deallocate_shared( wy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_xy_indexing

  SUBROUTINE determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)
    ! Determine the indexing and dimension directions of a variable in a lon/lat-grid file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=256),                  INTENT(OUT)   :: indexing
    CHARACTER(LEN=256),                  INTENT(OUT)   :: londir
    CHARACTER(LEN=256),                  INTENT(OUT)   :: latdir

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_lonlat_indexing'
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: nlon, nlat
    REAL(dp), DIMENSION(:    ), POINTER                :: lon, lat
    INTEGER                                            :: wlon, wlat
    INTEGER                                            :: id_var_lon, id_var_lat, id_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the lon and lat dimensions and variables of this file are valid
    CALL check_lon( filename, ncid)
    CALL check_lat( filename, ncid)

    ! Inquire lon and lat dimensions
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = nlon)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = nlat)

    ! Allocate memory for lon and lat
    CALL allocate_shared_dp_1D( nlon, lon, wlon)
    CALL allocate_shared_dp_1D( nlat, lat, wlat)

    ! Inquire lon and lon variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lon, id_var_lon)
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read lon and lat
    CALL read_var_dp_1D(  filename, ncid, id_var_lon, lon)
    CALL read_var_dp_1D(  filename, ncid, id_var_lat, lat)

    ! Determine directions of x and y
    IF (lon( 2) > lon( 1)) THEN
      londir = 'normal'
    ELSE
      londir = 'reverse'
    END IF
    IF (lat( 2) > lat( 1)) THEN
      latdir = 'normal'
    ELSE
      latdir = 'reverse'
    END IF

    ! Inquire dimensions of the specified field variable
    CALL inquire_var_multiple_options( filename, ncid, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    IF     (dims_of_var( 1) == id_dim_lon .AND. dims_of_var( 2) == id_dim_lat) THEN
      indexing = 'lonlat'
    ELSEIF (dims_of_var( 1) == id_dim_lat .AND. dims_of_var( 2) == id_dim_lon) THEN
      indexing = 'latlon'
    ELSE
      CALL crash('longitude and latitude are not the first two dimensions of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wlon)
    CALL deallocate_shared( wlat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_lonlat_indexing

  ! ===== Longitude corrections to a lon/lat grid + data field =====
  ! ================================================================

  SUBROUTINE correct_longitude_shifts_and_range_2D( filename, grid_lonlat, d_grid)
    ! Make sure longitude is bounded between 0 and 360, and increases monotonically

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid_lonlat
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'correct_longitude_shifts_and_range_2D'
    INTEGER                                            :: i,j
    LOGICAL                                            :: is_correct
    INTEGER                                            :: n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the grid is already correct, do nothing
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    CALL sync
    IF (is_correct) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Limit values to [0,360]
    IF (par%master) THEN
      DO i = 1, grid_lonlat%nlon
        IF     (grid_lonlat%lon( i) <   0._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) + 360._dp
        ELSEIF (grid_lonlat%lon( i) > 360._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) - 360._dp
        END IF
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Fix shifts
    n = 0
    DO i = 1, grid_lonlat%nlon-1
      IF (grid_lonlat%lon( i) > grid_lonlat%lon( i+1)) THEN
        n = i
        EXIT
      END IF
    END DO

    IF (n > 0) THEN

      ! Fix lon
      IF (par%master) THEN
        grid_lonlat%lon = [grid_lonlat%lon( n+1:grid_lonlat%nlon), grid_lonlat%lon( 1:n)]
      END IF
      CALL sync

      ! Fix data field
      DO j = grid_lonlat%j1, grid_lonlat%j2
        d_grid( :,j) = [d_grid( n+1:grid_lonlat%nlon,j), d_grid( 1:n,j)]
      END DO
      CALL sync

    END IF ! IF (n > 0) THEN

    ! The grid should now be correct
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    IF (.NOT. is_correct) CALL crash('something is seriously wrong with the longitude of file "' // TRIM( filename) // '"!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_longitude_shifts_and_range_2D

  SUBROUTINE correct_longitude_shifts_and_range_3D( filename, grid_lonlat, d_grid)
    ! Make sure longitude is bounded between 0 and 360, and increases monotonically

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid_lonlat
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'correct_longitude_shifts_and_range_3D'
    INTEGER                                            :: i,j,k
    LOGICAL                                            :: is_correct
    INTEGER                                            :: n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the grid is already correct, do nothing
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    CALL sync
    IF (is_correct) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Limit values to [0,360]
    IF (par%master) THEN
      DO i = 1, grid_lonlat%nlon
        IF     (grid_lonlat%lon( i) <   0._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) + 360._dp
        ELSEIF (grid_lonlat%lon( i) > 360._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) - 360._dp
        END IF
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Fix shifts
    n = 0
    DO i = 1, grid_lonlat%nlon-1
      IF (grid_lonlat%lon( i) > grid_lonlat%lon( i+1)) THEN
        n = i
        EXIT
      END IF
    END DO

    IF (n > 0) THEN

      ! Fix lon
      IF (par%master) THEN
        grid_lonlat%lon = [grid_lonlat%lon( n+1:grid_lonlat%nlon), grid_lonlat%lon( 1:n)]
      END IF
      CALL sync

      ! Fix data field
      DO j = grid_lonlat%j1, grid_lonlat%j2
      DO k = 1, SIZE( d_grid,3)
        d_grid( :,j,k) = [d_grid( n+1:grid_lonlat%nlon,j,k), d_grid( 1:n,j,k)]
      END DO
      END DO
      CALL sync

    END IF ! IF (n > 0) THEN

    ! The grid should now be correct
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    IF (.NOT. is_correct) CALL crash('something is seriously wrong with the longitude of file "' // TRIM( filename) // '"!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_longitude_shifts_and_range_3D

END MODULE netcdf_input_module