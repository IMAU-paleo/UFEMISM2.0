MODULE netcdf_input

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

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mpi_distributed_memory                                 , ONLY: distribute_from_master_dp_1D, distribute_from_master_dp_2D
  USE grid_basic                                             , ONLY: type_grid, calc_secondary_grid_data, deallocate_grid, &
                                                                     distribute_gridded_data_from_master_dp_2D, distribute_gridded_data_from_master_dp_3D
  USE grid_lonlat_basic                                      , ONLY: type_grid_lonlat, calc_lonlat_field_to_vector_form_translation_tables, &
                                                                     distribute_lonlat_gridded_data_from_master_dp_2D, deallocate_lonlat_grid, &
                                                                     distribute_lonlat_gridded_data_from_master_dp_3D
  USE math_utilities                                         , ONLY: permute_2D_int, permute_2D_dp, permute_3D_int, permute_3D_dp, &
                                                                     flip_1D_dp, flip_2D_x1_dp, flip_2D_x2_dp, flip_3D_x1_dp, flip_3D_x2_dp, flip_3D_x3_dp
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: allocate_mesh_primary, deallocate_mesh
  USE mesh_utilities                                         , ONLY: check_mesh
  USE mesh_secondary                                         , ONLY: calc_all_secondary_mesh_data
  USE mesh_parallel_creation                                 , ONLY: broadcast_merged_mesh
  USE mesh_remapping                                         , ONLY: map_from_xy_grid_to_mesh_2D, map_from_lonlat_grid_to_mesh_2D, map_from_mesh_to_mesh_2D, &
                                                                     map_from_xy_grid_to_mesh_3D, map_from_lonlat_grid_to_mesh_3D, map_from_mesh_to_mesh_3D

  USE netcdf      , ONLY: NF90_MAX_VAR_DIMS
  USE netcdf_basic, ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, &
                          field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                          field_name_options_dim_nV, field_name_options_dim_nTri, field_name_options_dim_nC_mem, &
                          field_name_options_dim_nE, field_name_options_dim_two, field_name_options_dim_three, &
                          field_name_options_dim_four, field_name_options_V, field_name_options_Tri, field_name_options_nC, &
                          field_name_options_C, field_name_options_niTri, field_name_options_iTri, &
                          field_name_options_VBI, field_name_options_Tricc, field_name_options_TriC, &
                          field_name_options_TriBI, field_name_options_E, field_name_options_VE, field_name_options_EV, &
                          field_name_options_ETri, field_name_options_EBI, field_name_options_A, field_name_options_R, &
                          field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                          field_name_options_SL, field_name_options_Ti, get_first_option_from_list, &
                          open_existing_netcdf_file_for_reading, close_netcdf_file, &
                          inquire_dim_multopt, inquire_var_multopt, &
                          read_var_master_int_0D, read_var_master_int_1D, read_var_master_int_2D, read_var_master_int_3D, read_var_master_int_4D, &
                          read_var_master_dp_0D , read_var_master_dp_1D , read_var_master_dp_2D , read_var_master_dp_3D , read_var_master_dp_4D, &
                          check_x, check_y, check_lon, check_lat, check_mesh_dimensions, check_zeta, check_month, find_timeframe, &
                          check_xy_grid_field_int_2D, check_xy_grid_field_dp_2D, check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, &
                          check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, check_lonlat_grid_field_dp_2D_monthly, check_lonlat_grid_field_dp_3D, &
                          check_mesh_field_int_2D, check_mesh_field_dp_2D, check_mesh_field_dp_2D_monthly, check_mesh_field_dp_3D, &
                          inquire_xy_grid, inquire_lonlat_grid, inquire_mesh

  IMPLICIT NONE

CONTAINS

  ! ===== Top-level functions =====
  ! ===============================

  ! Read and map to mesh
  SUBROUTINE read_field_from_file_2D(         filename, field_name_options, mesh, d_partial, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.
    !
    ! NOTE: memory for d_partial is allocated here!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT)   :: d_partial
    REAL(dp), OPTIONAL,                      INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_file_2D'
    LOGICAL                                                :: file_exists
    LOGICAL                                                :: has_xy_grid, has_lonlat_grid, has_mesh
    TYPE(type_grid)                                        :: grid_from_file
    TYPE(type_grid_lonlat)                                 :: grid_lonlat_from_file
    TYPE(type_mesh)                                        :: mesh_from_file
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: d_grid_vec_partial_from_file
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: d_grid_lonlat_vec_partial_from_file
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: d_mesh_partial_from_file
    CHARACTER(LEN=256), PARAMETER                          :: method_mesh2mesh = '2nd_order_conservative'

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
      CALL read_field_from_xy_file_2D( filename, field_name_options, d_grid_vec_partial_from_file, grid = grid_from_file, time_to_read = time_to_read)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc), source = 0._dp)

      ! Remap data
      CALL map_from_xy_grid_to_mesh_2D( grid_from_file, mesh, d_grid_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      CALL deallocate_grid( grid_from_file)
      DEALLOCATE( d_grid_vec_partial_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_2D( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, grid = grid_lonlat_from_file, time_to_read = time_to_read)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc), source = 0._dp)

      ! Remap data
      CALL map_from_lonlat_grid_to_mesh_2D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      CALL deallocate_lonlat_grid( grid_lonlat_from_file)
      DEALLOCATE( d_grid_lonlat_vec_partial_from_file)

    ELSEIF (has_mesh) THEN
      ! Data is provided on a mesh

      ! Read grid and gridded data
      CALL read_field_from_mesh_file_2D( filename, field_name_options, d_mesh_partial_from_file, mesh = mesh_from_file, time_to_read = time_to_read)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc), source = 0._dp)

      ! Remap data
      CALL map_from_mesh_to_mesh_2D( mesh_from_file, mesh, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      CALL deallocate_mesh( mesh_from_file)
      DEALLOCATE( d_mesh_partial_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_2D

  SUBROUTINE read_field_from_file_2D_monthly( filename, field_name_options, mesh, d_partial, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.
    !
    ! NOTE: memory for d_partial is allocated here!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_partial
    REAL(dp), OPTIONAL,                      INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_file_2D_monthly'
    LOGICAL                                                :: file_exists
    LOGICAL                                                :: has_xy_grid, has_lonlat_grid, has_mesh
    TYPE(type_grid)                                        :: grid_from_file
    TYPE(type_grid_lonlat)                                 :: grid_lonlat_from_file
    TYPE(type_mesh)                                        :: mesh_from_file
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_grid_vec_partial_from_file
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_grid_lonlat_vec_partial_from_file
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_mesh_partial_from_file
    CHARACTER(LEN=256), PARAMETER                          :: method_mesh2mesh = '2nd_order_conservative'

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
      CALL read_field_from_xy_file_2D_monthly( filename, field_name_options, d_grid_vec_partial_from_file, grid = grid_from_file, time_to_read = time_to_read)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc, 12), source = 0._dp)

      ! Remap data
      CALL map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, d_grid_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      CALL deallocate_grid( grid_from_file)
      DEALLOCATE( d_grid_vec_partial_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_2D_monthly( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, grid = grid_lonlat_from_file, time_to_read = time_to_read)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc, 12), source = 0._dp)

      ! Remap data
      CALL map_from_lonlat_grid_to_mesh_3D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      CALL deallocate_lonlat_grid( grid_lonlat_from_file)
      DEALLOCATE( d_grid_lonlat_vec_partial_from_file)

    ELSEIF (has_mesh) THEN
      ! Data is provided on a mesh

      ! Read grid and gridded data
      CALL read_field_from_mesh_file_2D_monthly( filename, field_name_options, d_mesh_partial_from_file, mesh = mesh_from_file, time_to_read = time_to_read)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc, 12), source = 0._dp)

      ! Remap data
      CALL map_from_mesh_to_mesh_3D( mesh_from_file, mesh, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      CALL deallocate_mesh( mesh_from_file)
      DEALLOCATE( d_mesh_partial_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_2D_monthly

  SUBROUTINE read_field_from_file_3D( filename, field_name_options, mesh, d_partial, time_to_read, nzeta, zeta)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.
    !
    ! NOTE: memory for d_partial is allocated here!
    !
    ! NOTE: data is returned on the horizontal model mesh, but on the vertical grid
    !       of the input file!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_partial
    REAL(dp), OPTIONAL,                      INTENT(IN)    :: time_to_read
    INTEGER ,                                INTENT(OUT), OPTIONAL :: nzeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT), OPTIONAL :: zeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_file_3D'
    LOGICAL                                                :: file_exists
    LOGICAL                                                :: has_xy_grid, has_lonlat_grid, has_mesh
    TYPE(type_grid)                                        :: grid_from_file
    TYPE(type_grid_lonlat)                                 :: grid_lonlat_from_file
    TYPE(type_mesh)                                        :: mesh_from_file
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_grid_vec_partial_from_file
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_grid_lonlat_vec_partial_from_file
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_mesh_partial_from_file
    CHARACTER(LEN=256), PARAMETER                          :: method_mesh2mesh = '2nd_order_conservative'
    INTEGER                                                :: nzeta_loc
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: zeta_loc

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
      CALL read_field_from_xy_file_3D( filename, field_name_options, d_grid_vec_partial_from_file, grid = grid_from_file, &
        time_to_read = time_to_read, nzeta = nzeta_loc, zeta = zeta_loc)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc, nzeta_loc), source = 0._dp)

      ! Remap data
      CALL map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, d_grid_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      CALL deallocate_grid( grid_from_file)
      DEALLOCATE( d_grid_vec_partial_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_3D( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, grid = grid_lonlat_from_file, &
        time_to_read = time_to_read, nzeta = nzeta_loc, zeta = zeta_loc)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc, nzeta_loc), source = 0._dp)

      ! Remap data
      CALL map_from_lonlat_grid_to_mesh_3D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      CALL deallocate_lonlat_grid( grid_lonlat_from_file)
      DEALLOCATE( d_grid_lonlat_vec_partial_from_file)

    ELSEIF (has_mesh) THEN
      ! Data is provided on a mesh

      ! Read grid and gridded data
      CALL read_field_from_mesh_file_3D( filename, field_name_options, d_mesh_partial_from_file, mesh = mesh_from_file, &
        time_to_read = time_to_read, nzeta = nzeta_loc, zeta = zeta_loc)

      ! Allocate memory for data on the model mesh
      ALLOCATE( d_partial( mesh%nV_loc, nzeta_loc), source = 0._dp)

      ! Remap data
      CALL map_from_mesh_to_mesh_3D( mesh_from_file, mesh, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      CALL deallocate_mesh( mesh_from_file)
      DEALLOCATE( d_mesh_partial_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    END IF

    ! If so specified, return the zeta read from file as output
    IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN
      ! Safety
      IF (.NOT. PRESENT( nzeta) .OR. .NOT. PRESENT( zeta)) CALL crash('should ask for both nzeta and zeta!')
      nzeta = nzeta_loc
      ALLOCATE( zeta( nzeta))
      zeta = zeta_loc
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_3D

  ! ===== Medium-level functions =====
  ! ==================================

  ! Read data fields from an x/y-grid file
  SUBROUTINE read_field_from_xy_file_2D(             filename, field_name_options, d_grid_vec_partial, grid, time_to_read)
    ! Read a 2-D data field from a NetCDF file on an x/y-grid, and optionally return the grid as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT)   :: d_grid_vec_partial
    TYPE(type_grid),            OPTIONAL,    INTENT(OUT)   :: grid
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_xy_file_2D'
    INTEGER                                                :: ncid
    TYPE(type_grid)                                        :: grid_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    CHARACTER(LEN=256)                                     :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_grid
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE                :: d_grid_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nx, grid_loc%ny))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_2D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nx, grid_loc%ny, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%ny, grid_loc%nx))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_2D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%ny, grid_loc%nx, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Perform necessary corrections to the gridded data
  ! ====================================================

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      IF (par%master) CALL permute_2D_dp( d_grid, map = [2,1])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid_loc%x)
      IF (par%master) CALL flip_2D_x1_dp( d_grid)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid_loc%y)
      IF (par%master) CALL flip_2D_x2_dp( d_grid)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid_loc%n_loc))

    ! Distribute data
    CALL distribute_gridded_data_from_master_dp_2D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_grid)

  ! == If so specified, return the read grid as output
  ! ==================================================

    IF (PRESENT( grid)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_xy_grid_from_file( filename, ncid, grid)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    ! Clean up after yourself
    CALL deallocate_grid( grid_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_xy_file_2D

  SUBROUTINE read_field_from_xy_file_2D_monthly(     filename, field_name_options, d_grid_vec_partial, grid, time_to_read)
    ! Read a 2-D monthly data field from a NetCDF file on an x/y-grid, and optionally return the grid as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_grid_vec_partial
    TYPE(type_grid),            OPTIONAL,    INTENT(OUT)   :: grid
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_xy_file_2D_monthly'
    INTEGER                                                :: ncid
    TYPE(type_grid)                                        :: grid_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    CHARACTER(LEN=256)                                     :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE              :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE              :: d_grid_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the file has a valid month dimension
    CALL check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nx, grid_loc%ny, 12))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nx, grid_loc%ny, 12, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 12, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%ny, grid_loc%nx, 12))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%ny, grid_loc%nx, 12, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 12, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Perform necessary corrections to the gridded data
  ! ====================================================

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      IF (par%master) CALL permute_3D_dp( d_grid, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid_loc%x)
      IF (par%master) CALL flip_3D_x1_dp( d_grid)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid_loc%y)
      IF (par%master) CALL flip_3D_x2_dp( d_grid)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid_loc%n_loc, 12))

    ! Distribute data
    CALL distribute_gridded_data_from_master_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_grid)

  ! == If so specified, return the read grid as output
  ! ==================================================

    IF (PRESENT( grid)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_xy_grid_from_file( filename, ncid, grid)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    ! Clean up after yourself
    CALL deallocate_grid( grid_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_xy_file_2D_monthly

  SUBROUTINE read_field_from_xy_file_3D(             filename, field_name_options, d_grid_vec_partial, grid, time_to_read, nzeta, zeta)
    ! Read a 3-D data field from a NetCDF file on an x/y-grid, and optionally return the grid as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_grid_vec_partial
    TYPE(type_grid),            OPTIONAL,    INTENT(OUT)   :: grid
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read
    INTEGER ,                   OPTIONAL,    INTENT(OUT)   :: nzeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT), OPTIONAL :: zeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_xy_file_3D'
    INTEGER                                                :: ncid
    TYPE(type_grid)                                        :: grid_loc
    INTEGER                                                :: nzeta_loc
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: zeta_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    CHARACTER(LEN=256)                                     :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE              :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE              :: d_grid_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, ncid, grid_loc)

    ! Set up the vertical coordinate zeta from the file
    CALL setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nx, grid_loc%ny, nzeta_loc))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nx, grid_loc%ny, nzeta_loc, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, nzeta_loc, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%ny, grid_loc%nx, nzeta_loc))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%ny, grid_loc%nx, nzeta_loc, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, nzeta_loc, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Perform necessary corrections to the gridded data
  ! ====================================================

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      IF (par%master) CALL permute_3D_dp( d_grid, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid_loc%x)
      IF (par%master) CALL flip_3D_x1_dp( d_grid)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid_loc%y)
      IF (par%master) CALL flip_3D_x2_dp( d_grid)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid_loc%n_loc, nzeta_loc))

    ! Distribute data
    CALL distribute_gridded_data_from_master_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_grid)

  ! == If so specified, return the read grid as output
  ! ==================================================

    IF (PRESENT( grid)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_xy_grid_from_file( filename, ncid, grid)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN
      ! Safety
      IF (.NOT. PRESENT( nzeta) .OR. .NOT. PRESENT( zeta)) CALL crash('should ask for both nzeta and zeta!')
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_zeta_from_file( filename, ncid, nzeta, zeta)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN

    ! Clean up after yourself
    CALL deallocate_grid( grid_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_xy_file_3D

  ! Read data fields from an x/y-grid file
  SUBROUTINE read_field_from_lonlat_file_2D(         filename, field_name_options, d_grid_vec_partial, grid, time_to_read)
    ! Read a 2-D data field from a NetCDF file on a lon/lat-grid, and optionally return the grid as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT)   :: d_grid_vec_partial
    TYPE(type_grid_lonlat),     OPTIONAL,    INTENT(OUT)   :: grid
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_lonlat_file_2D'
    INTEGER                                                :: ncid
    TYPE(type_grid_lonlat)                                 :: grid_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    CHARACTER(LEN=256)                                     :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_grid
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE                :: d_grid_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nlon, grid_loc%nlat))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_2D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nlat, grid_loc%nlon))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_2D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Perform necessary corrections to the gridded data
  ! ====================================================

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      IF (par%master) CALL permute_2D_dp( d_grid, map = [2,1])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      IF (par%master) CALL flip_2D_x1_dp( d_grid)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      IF (par%master) CALL flip_2D_x2_dp( d_grid)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid_loc%n_loc))

    ! Distribute data
    CALL distribute_lonlat_gridded_data_from_master_dp_2D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_grid)

  ! == If so specified, return the read grid as output
  ! ==================================================

    IF (PRESENT( grid)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_lonlat_grid_from_file( filename, ncid, grid)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    ! Clean up after yourself
    CALL deallocate_lonlat_grid( grid_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_2D

  SUBROUTINE read_field_from_lonlat_file_2D_monthly( filename, field_name_options, d_grid_vec_partial, grid, time_to_read)
    ! Read a 2-D monthly data field from a NetCDF file on a lon/lat-grid, and optionally return the grid as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_grid_vec_partial
    TYPE(type_grid_lonlat),     OPTIONAL,    INTENT(OUT)   :: grid
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_lonlat_file_2D_monthly'
    INTEGER                                                :: ncid
    TYPE(type_grid_lonlat)                                 :: grid_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    CHARACTER(LEN=256)                                     :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE              :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE              :: d_grid_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the file has a valid month dimension
    CALL check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nlon, grid_loc%nlat, 12))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, 12, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, 12, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nlat, grid_loc%nlon, 12))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, 12, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, 12, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Perform necessary corrections to the gridded data
  ! ====================================================

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      IF (par%master) CALL permute_3D_dp( d_grid, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      IF (par%master) CALL flip_3D_x1_dp( d_grid)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      IF (par%master) CALL flip_3D_x2_dp( d_grid)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid_loc%n_loc, 12))

    ! Distribute data
    CALL distribute_lonlat_gridded_data_from_master_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_grid)

  ! == If so specified, return the read grid as output
  ! ==================================================

    IF (PRESENT( grid)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_lonlat_grid_from_file( filename, ncid, grid)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    ! Clean up after yourself
    CALL deallocate_lonlat_grid( grid_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_2D_monthly

  SUBROUTINE read_field_from_lonlat_file_3D(         filename, field_name_options, d_grid_vec_partial, grid, time_to_read, nzeta, zeta)
    ! Read a 2-D monthly data field from a NetCDF file on a lon/lat-grid, and optionally return the grid as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_grid_vec_partial
    TYPE(type_grid_lonlat),     OPTIONAL,    INTENT(OUT)   :: grid
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read
    INTEGER ,                   OPTIONAL,    INTENT(OUT)   :: nzeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT), OPTIONAL :: zeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_lonlat_file_3D'
    INTEGER                                                :: ncid
    TYPE(type_grid_lonlat)                                 :: grid_loc
    INTEGER                                                :: nzeta_loc
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: zeta_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    CHARACTER(LEN=256)                                     :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE              :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE              :: d_grid_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, ncid, grid_loc)

    ! Set up the vertical coordinate zeta from the file
    CALL setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nlon, grid_loc%nlat, nzeta_loc))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, nzeta_loc, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, nzeta_loc, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate memory
      IF (par%master) ALLOCATE( d_grid( grid_loc%nlat, grid_loc%nlon, nzeta_loc))

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_master_dp_3D( filename, ncid, id_var, d_grid)
      ELSE
        ! Allocate memory
        IF (par%master) ALLOCATE( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, nzeta_loc, 1))
        ! Find out which timeframe to read
        CALL find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        CALL read_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, nzeta_loc, 1 /) )
        ! Copy to output memory
        IF (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        IF (par%master) DEALLOCATE( d_grid_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Perform necessary corrections to the gridded data
  ! ====================================================

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      IF (par%master) CALL permute_3D_dp( d_grid, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      IF (par%master) CALL flip_3D_x1_dp( d_grid)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      IF (par%master) CALL flip_3D_x2_dp( d_grid)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid_loc%n_loc, nzeta_loc))

    ! Distribute data
    CALL distribute_lonlat_gridded_data_from_master_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_grid)

  ! == If so specified, return the read grid as output
  ! ==================================================

    IF (PRESENT( grid)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_lonlat_grid_from_file( filename, ncid, grid)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN
      ! Safety
      IF (.NOT. PRESENT( nzeta) .OR. .NOT. PRESENT( zeta)) CALL crash('should ask for both nzeta and zeta!')
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_zeta_from_file( filename, ncid, nzeta, zeta)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN

    ! Clean up after yourself
    CALL deallocate_lonlat_grid( grid_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_3D

  ! Read data fields from a mesh file
  SUBROUTINE read_field_from_mesh_file_2D(           filename, field_name_options, d_mesh_partial, mesh, time_to_read)
    ! Read a 2-D data field from a NetCDF file on a mesh, and optionally return the mesh as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT)   :: d_mesh_partial
    TYPE(type_mesh),            OPTIONAL,    INTENT(OUT)   :: mesh
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_mesh_file_2D'
    INTEGER                                                :: ncid
    TYPE(type_mesh)                                        :: mesh_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: d_mesh
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_mesh_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Allocate memory
    IF (par%master) ALLOCATE( d_mesh( mesh_loc%nV))

    ! Read data from file
    IF (.NOT. PRESENT( time_to_read)) THEN
      CALL read_var_master_dp_1D( filename, ncid, id_var, d_mesh)
    ELSE
      ! Allocate memory
      IF (par%master) ALLOCATE( d_mesh_with_time( mesh_loc%nV, 1))
      ! Find out which timeframe to read
      CALL find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      CALL read_var_master_dp_2D( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, ti /), count = (/ mesh_loc%nV, 1 /) )
      ! Copy to output memory
      IF (par%master) d_mesh = d_mesh_with_time( :,1)
      ! Clean up after yourself
      IF (par%master) DEALLOCATE( d_mesh_with_time)
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_mesh_partial( mesh_loc%nV_loc))

    ! Distribute data
    CALL distribute_from_master_dp_1D( d_mesh, d_mesh_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_mesh)

  ! == If so specified, return the read mesh as output
  ! ==================================================

    IF (PRESENT( mesh)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_mesh_from_file( filename, ncid, mesh)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    ! Clean up after yourself
    CALL deallocate_mesh( mesh_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_mesh_file_2D

  SUBROUTINE read_field_from_mesh_file_2D_monthly(   filename, field_name_options, d_mesh_partial, mesh, time_to_read)
    ! Read a 2-D data monthly field from a NetCDF file on a mesh, and optionally return the mesh as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_mesh_partial
    TYPE(type_mesh),            OPTIONAL,    INTENT(OUT)   :: mesh
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_mesh_file_2D_monthly'
    INTEGER                                                :: ncid
    TYPE(type_mesh)                                        :: mesh_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_mesh
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE                :: d_mesh_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the file has a valid month dimension
    CALL check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    CALL check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Allocate memory
    IF (par%master) ALLOCATE( d_mesh( mesh_loc%nV, 12))

    ! Read data from file
    IF (.NOT. PRESENT( time_to_read)) THEN
      CALL read_var_master_dp_2D( filename, ncid, id_var, d_mesh)
    ELSE
      ! Allocate memory
      IF (par%master) ALLOCATE( d_mesh_with_time( mesh_loc%nV, 12, 1))
      ! Find out which timeframe to read
      CALL find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      CALL read_var_master_dp_3D( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, 12, 1 /) )
      ! Copy to output memory
      IF (par%master) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      IF (par%master) DEALLOCATE( d_mesh_with_time)
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_mesh_partial( mesh_loc%nV_loc, 12))

    ! Distribute data
    CALL distribute_from_master_dp_2D( d_mesh, d_mesh_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_mesh)

  ! == If so specified, return the read mesh as output
  ! ==================================================

    IF (PRESENT( mesh)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_mesh_from_file( filename, ncid, mesh)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    ! Clean up after yourself
    CALL deallocate_mesh( mesh_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_mesh_file_2D_monthly

  SUBROUTINE read_field_from_mesh_file_3D(           filename, field_name_options, d_mesh_partial, mesh, time_to_read, nzeta, zeta)
    ! Read a 2-D data monthly field from a NetCDF file on a mesh, and optionally return the mesh as well.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                        INTENT(IN)    :: filename
    CHARACTER(LEN=*),                        INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: d_mesh_partial
    TYPE(type_mesh),            OPTIONAL,    INTENT(OUT)   :: mesh
    REAL(dp),                   OPTIONAL,    INTENT(IN)    :: time_to_read
    INTEGER ,                   OPTIONAL,    INTENT(OUT)   :: nzeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(OUT), OPTIONAL :: zeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'read_field_from_mesh_file_3D'
    INTEGER                                                :: ncid
    TYPE(type_mesh)                                        :: mesh_loc
    INTEGER                                                :: nzeta_loc
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: zeta_loc
    INTEGER                                                :: id_var
    CHARACTER(LEN=256)                                     :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_mesh
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE                :: d_mesh_with_time
    INTEGER                                                :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    CALL setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Set up the vertical coordinate zeta from the file
    CALL setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = PRESENT( time_to_read))

    ! Allocate memory
    IF (par%master) ALLOCATE( d_mesh( mesh_loc%nV, nzeta_loc))

    ! Read data from file
    IF (.NOT. PRESENT( time_to_read)) THEN
      CALL read_var_master_dp_2D( filename, ncid, id_var, d_mesh)
    ELSE
      ! Allocate memory
      IF (par%master) ALLOCATE( d_mesh_with_time( mesh_loc%nV, nzeta_loc, 1))
      ! Find out which timeframe to read
      CALL find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      CALL read_var_master_dp_3D( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, nzeta_loc, 1 /) )
      ! Copy to output memory
      IF (par%master) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      IF (par%master) DEALLOCATE( d_mesh_with_time)
    END IF

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

  ! == Distribute gridded data from the master to all processes in partial vector form
  ! ==================================================================================

    ! Allocate memory
    ALLOCATE( d_mesh_partial( mesh_loc%nV_loc, nzeta_loc))

    ! Distribute data
    CALL distribute_from_master_dp_2D( d_mesh, d_mesh_partial)

    ! Clean up gridded data on the master
    IF (par%master) DEALLOCATE( d_mesh)

  ! == If so specified, return the read mesh as output
  ! ==================================================

    IF (PRESENT( mesh)) THEN
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_mesh_from_file( filename, ncid, mesh)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( grid)) THEN

    IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN
      ! Safety
      IF (.NOT. PRESENT( nzeta) .OR. .NOT. PRESENT( zeta)) CALL crash('should ask for both nzeta and zeta!')
      CALL open_existing_netcdf_file_for_reading( filename, ncid)
      CALL setup_zeta_from_file( filename, ncid, nzeta, zeta)
      CALL close_netcdf_file( ncid)
    END IF ! IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN

    ! Clean up after yourself
    CALL deallocate_mesh( mesh_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_mesh_file_3D

  ! ===== Set up grids/mesh from a NetCDF file =====
  ! ================================================

  SUBROUTINE setup_xy_grid_from_file(     filename, ncid, grid)
    ! Set up an x/y-grid from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_grid),                     INTENT(OUT)   :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_xy_grid_from_file'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var_x, id_var_y

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Give the grid a nice name
    grid%name = 'xy_grid_from_file_"' // TRIM( filename) // '"'

    ! Check grid dimensions and variables for validity
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Inquire x and y dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x, dim_length = grid%nx)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y, dim_length = grid%ny)

    ! Allocate memory for x and y
    ALLOCATE( grid%x( grid%nx))
    ALLOCATE( grid%y( grid%ny))

    ! Inquire x and y variables
    CALL inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    CALL inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Read x and y
    CALL read_var_master_dp_1D(  filename, ncid, id_var_x, grid%x)
    CALL read_var_master_dp_1D(  filename, ncid, id_var_y, grid%y)

    ! Broadcast x and y from the master to the other processes
    CALL MPI_BCAST( grid%x, grid%nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( grid%y, grid%ny, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Calculate secondary grid geometry data
    CALL calc_secondary_grid_data( grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_xy_grid_from_file

  SUBROUTINE setup_lonlat_grid_from_file( filename, ncid, grid)
    ! Set up a lon/lat-grid from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_grid_lonlat),              INTENT(OUT)   :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_lonlat_grid_from_file'
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: id_var_lon, id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Give the grid a nice name
    grid%name = 'lonlat_grid_from_file_"' // TRIM( filename) // '"'

    ! Check grid dimensions and variables for validity
    CALL check_lon( filename, ncid)
    CALL check_lat( filename, ncid)

    ! Inquire lon and lat dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = grid%nlon)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = grid%nlat)

    ! Allocate memory for lon and lat
    ALLOCATE( grid%lon( grid%nlon))
    ALLOCATE( grid%lat( grid%nlat))

    ! Inquire lon and lat variables
    CALL inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    CALL inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read x and y
    CALL read_var_master_dp_1D(  filename, ncid, id_var_lon, grid%lon)
    CALL read_var_master_dp_1D(  filename, ncid, id_var_lat, grid%lat)

    ! Broadcast x and y from the master to the other processes
    CALL MPI_BCAST( grid%lon, grid%nlon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( grid%lat, grid%nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Secondary data
    CALL calc_lonlat_field_to_vector_form_translation_tables( grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_lonlat_grid_from_file

  SUBROUTINE setup_mesh_from_file(        filename, ncid, mesh)
    ! Set up a mesh from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_mesh),                     INTENT(OUT)   :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_mesh_from_file'
    CHARACTER(LEN=256)                                 :: name
    INTEGER                                            :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    INTEGER                                            :: nV_mem, nTri_mem, nC_mem, n_two, n_three
    INTEGER                                            :: id_var_xmin, id_var_xmax, id_var_ymin, id_var_ymax, id_var_tol_dist, id_var_lambda_M, id_var_phi_M, id_var_beta_stereo
    INTEGER                                            :: id_var_V, id_var_nC, id_var_C, id_var_niTri, id_var_iTri, id_var_VBI
    INTEGER                                            :: id_var_Tri, id_var_Tricc, id_var_TriC
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Give the mesh a nice name
    name = 'mesh_from_file_"' // TRIM( filename) // '"'

    ! Check mesh dimensions and variables for validity
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   , dim_length = nV_mem  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   , dim_length = nTri_mem)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   , dim_length = nC_mem  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  , dim_length = n_two   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three, dim_length = n_three )

    ! Allocate memory for the mesh
    IF (par%master) THEN
      CALL allocate_mesh_primary( mesh, name, nV_mem, nTri_mem, nC_mem)
      mesh%nV   = mesh%nV_mem
      mesh%nTri = mesh%nTri_mem
    END IF

  ! == Inquire mesh variables
  ! =========================

    ! Metadata
    CALL inquire_var_multopt( filename, ncid, 'xmin'                           , id_var_xmin          )
    CALL inquire_var_multopt( filename, ncid, 'xmax'                           , id_var_xmax          )
    CALL inquire_var_multopt( filename, ncid, 'ymin'                           , id_var_ymin          )
    CALL inquire_var_multopt( filename, ncid, 'ymax'                           , id_var_ymax          )
    CALL inquire_var_multopt( filename, ncid, 'tol_dist'                       , id_var_tol_dist      )
    CALL inquire_var_multopt( filename, ncid, 'lambda_M'                       , id_var_lambda_M      )
    CALL inquire_var_multopt( filename, ncid, 'phi_M'                          , id_var_phi_M         )
    CALL inquire_var_multopt( filename, ncid, 'beta_stereo'                    , id_var_beta_stereo   )

    ! Vertex data
    CALL inquire_var_multopt( filename, ncid, field_name_options_V             , id_var_V             )
    CALL inquire_var_multopt( filename, ncid, field_name_options_nC            , id_var_nC            )
    CALL inquire_var_multopt( filename, ncid, field_name_options_C             , id_var_C             )
    CALL inquire_var_multopt( filename, ncid, field_name_options_niTri         , id_var_niTri         )
    CALL inquire_var_multopt( filename, ncid, field_name_options_iTri          , id_var_iTri          )
    CALL inquire_var_multopt( filename, ncid, field_name_options_VBI           , id_var_VBI           )

    ! Triangle data
    CALL inquire_var_multopt( filename, ncid, field_name_options_Tri           , id_var_Tri           )
    CALL inquire_var_multopt( filename, ncid, field_name_options_Tricc         , id_var_Tricc         )
    CALL inquire_var_multopt( filename, ncid, field_name_options_TriC          , id_var_TriC          )

  ! == Read mesh data
  ! =================

    ! Metadata
    CALL read_var_master_dp_0D(  filename, ncid, id_var_xmin          , mesh%xmin          )
    CALL read_var_master_dp_0D(  filename, ncid, id_var_xmax          , mesh%xmax          )
    CALL read_var_master_dp_0D(  filename, ncid, id_var_ymin          , mesh%ymin          )
    CALL read_var_master_dp_0D(  filename, ncid, id_var_ymax          , mesh%ymax          )
    CALL read_var_master_dp_0D(  filename, ncid, id_var_tol_dist      , mesh%tol_dist      )
    CALL read_var_master_dp_0D(  filename, ncid, id_var_lambda_M      , mesh%lambda_M      )
    CALL read_var_master_dp_0D(  filename, ncid, id_var_phi_M         , mesh%phi_M         )
    CALL read_var_master_dp_0D(  filename, ncid, id_var_beta_stereo   , mesh%beta_stereo   )

    ! Vertex data
    CALL read_var_master_dp_2D(  filename, ncid, id_var_V             , mesh%V             )
    CALL read_var_master_int_1D( filename, ncid, id_var_nC            , mesh%nC            )
    CALL read_var_master_int_2D( filename, ncid, id_var_C             , mesh%C             )
    CALL read_var_master_int_1D( filename, ncid, id_var_niTri         , mesh%niTri         )
    CALL read_var_master_int_2D( filename, ncid, id_var_iTri          , mesh%iTri          )
    CALL read_var_master_int_1D( filename, ncid, id_var_VBI           , mesh%VBI           )

    ! Triangle data
    CALL read_var_master_int_2D( filename, ncid, id_var_Tri           , mesh%Tri           )
    CALL read_var_master_dp_2D(  filename, ncid, id_var_Tricc         , mesh%Tricc         )
    CALL read_var_master_int_2D( filename, ncid, id_var_TriC          , mesh%TriC          )

    ! Safety - check if the mesh data read from NetCDF makes sense
    IF (par%master) CALL check_mesh( mesh)

    ! Broadcast read mesh from the master to the other processes
    CALL broadcast_merged_mesh( mesh)

    ! Calculate secondary mesh data
    CALL calc_all_secondary_mesh_data( mesh, mesh%lambda_M, mesh%phi_M, mesh%beta_stereo)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_mesh_from_file

  SUBROUTINE setup_zeta_from_file(        filename, ncid, nzeta, zeta)
    ! Set up a zeta coordinate from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(OUT)   :: nzeta
    REAL(dp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)   ::  zeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_zeta_from_file'
    INTEGER                                            :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check zeta dimension and variable for validity
    CALL check_zeta( filename, ncid)

    ! Inquire zeta dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta, dim_length = nzeta)

    ! Inquire zeta variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var_zeta)

    ! Allocate memory
    ALLOCATE( zeta( nzeta))

    ! Read zeta from file
    CALL read_var_master_dp_1D( filename, ncid, id_var_zeta, zeta)

    ! Broadcast zeta from master to all other processes
    CALL MPI_BCAST( zeta, nzeta, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_zeta_from_file

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
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x, y
    INTEGER                                            :: id_var_x, id_var_y, id_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the x and y dimensions and variables of this file are valid
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Inquire x and y dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x, dim_length = nx)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y, dim_length = ny)

    ! Allocate memory for x and y
    ALLOCATE( x( nx))
    ALLOCATE( y( ny))

    ! Inquire x and y variables
    CALL inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    CALL inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Read x and y
    CALL read_var_master_dp_1D(  filename, ncid, id_var_x, x)
    CALL read_var_master_dp_1D(  filename, ncid, id_var_y, y)
    CALL MPI_BCAST( x, nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( y, ny, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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
    CALL inquire_var_multopt( filename, ncid, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    IF     (dims_of_var( 1) == id_dim_x .AND. dims_of_var( 2) == id_dim_y) THEN
      indexing = 'xy'
    ELSEIF (dims_of_var( 1) == id_dim_y .AND. dims_of_var( 2) == id_dim_x) THEN
      indexing = 'yx'
    ELSE
      CALL crash('x and y are not the first two dimensions of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF

    ! Clean up after yourself
    DEALLOCATE( x)
    DEALLOCATE( y)

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
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: lon, lat
    INTEGER                                            :: id_var_lon, id_var_lat, id_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the lon and lat dimensions and variables of this file are valid
    CALL check_lon( filename, ncid)
    CALL check_lat( filename, ncid)

    ! Inquire lon and lat dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = nlon)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = nlat)

    ! Allocate memory for lon and lat
    ALLOCATE( lon( nlon))
    ALLOCATE( lat( nlat))

    ! Inquire lon and lon variables
    CALL inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    CALL inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read lon and lat
    CALL read_var_master_dp_1D(  filename, ncid, id_var_lon, lon)
    CALL read_var_master_dp_1D(  filename, ncid, id_var_lat, lat)
    CALL MPI_BCAST( lon, nlon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( lat, nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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
    CALL inquire_var_multopt( filename, ncid, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    IF     (dims_of_var( 1) == id_dim_lon .AND. dims_of_var( 2) == id_dim_lat) THEN
      indexing = 'lonlat'
    ELSEIF (dims_of_var( 1) == id_dim_lat .AND. dims_of_var( 2) == id_dim_lon) THEN
      indexing = 'latlon'
    ELSE
      CALL crash('longitude and latitude are not the first two dimensions of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF

    ! Clean up after yourself
    DEALLOCATE( lon)
    DEALLOCATE( lat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_lonlat_indexing

END MODULE netcdf_input