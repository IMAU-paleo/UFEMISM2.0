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
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE grid_basic                                             , ONLY: type_grid, type_grid_lonlat, calc_secondary_grid_data, deallocate_grid, &
                                                                     distribute_gridded_data_from_master_dp_2D, distribute_gridded_data_from_master_dp_3D
  USE math_utilities                                         , ONLY: permute_2D_int, permute_2D_dp, permute_3D_int, permute_3D_dp, &
                                                                     flip_1D_dp, flip_2D_x1_dp, flip_2D_x2_dp, flip_3D_x1_dp, flip_3D_x2_dp, flip_3D_x3_dp
  USE mesh_types                                             , ONLY: type_mesh

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
                          inquire_dim_multiple_options, inquire_var_multiple_options, &
                          read_var_master_int_0D, read_var_master_int_1D, read_var_master_int_2D, read_var_master_int_3D, read_var_master_int_4D, &
                          read_var_master_dp_0D , read_var_master_dp_1D , read_var_master_dp_2D , read_var_master_dp_3D , read_var_master_dp_4D, &
                          check_x, check_y, check_lon, check_lat, check_mesh_dimensions, check_zeta, find_timeframe, &
                          check_xy_grid_field_int_2D, check_xy_grid_field_dp_2D, check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, &
                          check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, check_lonlat_grid_field_dp_2D_monthly, check_lonlat_grid_field_dp_3D, &
                          check_mesh_field_int_2D, check_mesh_field_dp_2D, check_mesh_field_dp_2D_monthly, check_mesh_field_dp_3D, &
                          inquire_xy_grid, inquire_lonlat_grid, inquire_mesh

  IMPLICIT NONE

CONTAINS

  ! ===== Top-level functions =====
  ! ===============================

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
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
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

      ! Allocate memory
      ALLOCATE( grid%x( grid_loc%nx))
      ALLOCATE( grid%y( grid_loc%ny))

      ! Copy data
      grid%name     = grid_loc%name
      grid%nx       = grid_loc%nx
      grid%ny       = grid_loc%ny
      grid%dx       = grid_loc%dx
      grid%x        = grid_loc%x
      grid%y        = grid_loc%y

      ! Calculate secondary geometry data
      CALL calc_secondary_grid_data( grid)

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
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

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

      ! Allocate memory
      ALLOCATE( grid%x( grid_loc%nx))
      ALLOCATE( grid%y( grid_loc%ny))

      ! Copy data
      grid%name     = grid_loc%name
      grid%nx       = grid_loc%nx
      grid%ny       = grid_loc%ny
      grid%dx       = grid_loc%dx
      grid%x        = grid_loc%x
      grid%y        = grid_loc%y

      ! Calculate secondary geometry data
      CALL calc_secondary_grid_data( grid)

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
    CALL inquire_var_multiple_options( filename, ncid, field_name_options, id_var, var_name = var_name)
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

      ! Allocate memory
      ALLOCATE( grid%x( grid_loc%nx))
      ALLOCATE( grid%y( grid_loc%ny))

      ! Copy data
      grid%name     = grid_loc%name
      grid%nx       = grid_loc%nx
      grid%ny       = grid_loc%ny
      grid%dx       = grid_loc%dx
      grid%x        = grid_loc%x
      grid%y        = grid_loc%y

      ! Calculate secondary geometry data
      CALL calc_secondary_grid_data( grid)

    END IF ! IF (PRESENT( grid)) THEN

    IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN

      ! Safety
      IF (.NOT. PRESENT( nzeta) .OR. .NOT. PRESENT( zeta)) CALL crash('should ask for both nzeta and zeta!')

      nzeta = nzeta_loc
      ALLOCATE( zeta( nzeta))
      zeta = zeta_loc

    END IF ! IF (PRESENT( nzeta) .OR. PRESENT( zeta)) THEN

    ! Clean up after yourself
    CALL deallocate_grid( grid_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_xy_file_3D

  ! ===== Set up grids/mesh from a NetCDF file =====
  ! ================================================

  SUBROUTINE setup_xy_grid_from_file(     filename, ncid, grid)
    ! Set up an x/y-grid from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the grid at all

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
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_x, id_dim_x, dim_length = grid%nx)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_y, id_dim_y, dim_length = grid%ny)

    ! Allocate memory for x and y
    ALLOCATE( grid%x( grid%nx))
    ALLOCATE( grid%y( grid%ny))

    ! Inquire x and y variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_x, id_var_x)
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_y, id_var_y)

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

  SUBROUTINE setup_lonlat_grid_from_file( filename, ncid, grid, region_name)
    ! Set up a lon/lat-grid from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_grid_lonlat),              INTENT(OUT)   :: grid
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_lonlat_grid_from_file'
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: id_var_lon, id_var_lat
    INTEGER                                            :: i,j,n
    REAL(dp)                                           :: dlon

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

!    ! Check grid dimensions and variables for validity
!    CALL check_lon( filename, ncid)
!    CALL check_lat( filename, ncid)
!
!    ! Inquire lon and lat dimensions
!    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = grid%nlon)
!    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = grid%nlat)
!
!    grid%n = grid%nlon * grid%nlat
!
!    ! Allocate memory for lon and lat
!    ALLOCATE( grid%lon( grid%nlon))
!    ALLOCATE( grid%lat( grid%nlat))
!
!    ! Inquire lon and lat variables
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lon, id_var_lon)
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lat, id_var_lat)
!
!    ! Read lon and lat
!    CALL read_var_dp_1D(  filename, ncid, id_var_lon, grid%lon)
!    CALL read_var_dp_1D(  filename, ncid, id_var_lat, grid%lat)
!
!    ! Calculate secondary grid data
!
!    ! Resolution
!    grid%dlon = ABS( grid%lon( 2) - grid%lon( 1))
!    grid%dlat = ABS( grid%lat( 2) - grid%lat( 1))
!
!    ! Safety
!    DO i = 1, grid%nlon - 1
!      ! Check for regularity in longitude, but allow for 360-degree jumps
!      dlon = MIN( ABS( grid%lon( i+1) - grid%lon( i)), ABS( grid%lon( i+1) + 360._dp - grid%lon( i)))
!      IF (ABS( 1._dp - dlon / grid%dlon) > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular longitude dimension!')
!    END DO
!    DO j = 1, grid%nlat - 1
!      IF (ABS( 1._dp - ABS( grid%lat( j+1) - grid%lat( j)) / grid%dlat) > 1E-6_dp) &
!        CALL crash('file "' // TRIM( filename) // '" has an irregular latitude dimension!')
!    END DO
!
!    ! Domain size
!    grid%lonmin = MINVAL( grid%lon)
!    grid%lonmax = MAXVAL( grid%lon)
!    grid%latmin = MINVAL( grid%lat)
!    grid%latmax = MAXVAL( grid%lat)
!
!    ! Conversion tables for grid-form vs. vector-form data
!    n = 0
!    DO i = 1, grid%nlon
!      IF (MOD(i,2) == 1) THEN
!        DO j = 1, grid%nlat
!          n = n+1
!          grid%ij2n( i,j) = n
!          grid%n2ij( n,:) = [i,j]
!        END DO
!      ELSE
!        DO j = grid%nlat, 1, -1
!          n = n+1
!          grid%ij2n( i,j) = n
!          grid%n2ij( n,:) = [i,j]
!        END DO
!      END IF
!    END DO
!
!    ! Projection parameters for this region
!    IF     (region_name == 'NAM') THEN
!      grid%lambda_M     = C%lambda_M_NAM
!      grid%phi_M        = C%phi_M_NAM
!      grid%beta_stereo  = C%beta_stereo_NAM
!    ELSEIF (region_name == 'EAS') THEN
!      grid%lambda_M     = C%lambda_M_EAS
!      grid%phi_M        = C%phi_M_EAS
!      grid%beta_stereo  = C%beta_stereo_EAS
!    ELSEIF (region_name == 'GRL') THEN
!      grid%lambda_M     = C%lambda_M_GRL
!      grid%phi_M        = C%phi_M_GRL
!      grid%beta_stereo  = C%beta_stereo_GRL
!    ELSEIF (region_name == 'ANT') THEN
!      grid%lambda_M     = C%lambda_M_ANT
!      grid%phi_M        = C%phi_M_ANT
!      grid%beta_stereo  = C%beta_stereo_ANT
!    END IF
!
!    ! Set up parallelisation domains
!    CALL partition_list( grid%nlon, par%i, par%n, grid%i1, grid%i2)
!    CALL partition_list( grid%nlat, par%i, par%n, grid%j1, grid%j2)
!
!    ! Give the grid a nice name
!    grid%name = 'lonlat_grid_from_file_"' // TRIM( filename) // '"'

    ! Finalise routine path
    CALL finalise_routine( routine_name)

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

    ! DENK DROM
    CALL crash('fixme!')

!    ! Check mesh dimensions and variables for validity
!    CALL check_mesh_dimensions( filename, ncid)
!
!    ! Inquire mesh dimensions
!    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_nV    , id_dim_vi   , dim_length = nV_mem  )
!    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   , dim_length = nTri_mem)
!    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   , dim_length = nC_mem  )
!    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_two   , id_dim_two  , dim_length = n_two   )
!    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_dim_three , id_dim_three, dim_length = n_three )
!
!    ! Allocate memory for the mesh
!    CALL allocate_mesh_primary( mesh, region_name, nV_mem, nTri_mem, nC_mem)
!
!    ! Inquire mesh variables
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_V             , id_var_V             )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_nC            , id_var_nC            )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_C             , id_var_C             )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_niTri         , id_var_niTri         )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_iTri          , id_var_iTri          )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_edge_index    , id_var_edge_index    )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_Tri           , id_var_Tri           )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_Tricc         , id_var_Tricc         )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_TriC          , id_var_TriC          )
!    CALL inquire_var_multiple_options( filename, ncid, field_name_options_Tri_edge_index, id_var_Tri_edge_index)
!
!    ! Read mesh data
!    CALL read_var_dp_2D(  filename, ncid, id_var_V             , mesh%V             )
!    CALL read_var_int_1D( filename, ncid, id_var_nC            , mesh%nC            )
!    CALL read_var_int_2D( filename, ncid, id_var_C             , mesh%C             )
!    CALL read_var_int_1D( filename, ncid, id_var_niTri         , mesh%niTri         )
!    CALL read_var_int_2D( filename, ncid, id_var_iTri          , mesh%iTri          )
!    CALL read_var_int_1D( filename, ncid, id_var_edge_index    , mesh%edge_index    )
!    CALL read_var_int_2D( filename, ncid, id_var_Tri           , mesh%Tri           )
!    CALL read_var_dp_2D(  filename, ncid, id_var_Tricc         , mesh%Tricc         )
!    CALL read_var_int_2D( filename, ncid, id_var_TriC          , mesh%TriC          )
!    CALL read_var_int_1D( filename, ncid, id_var_Tri_edge_index, mesh%Tri_edge_index)
!
!    ! Calculate secondary mesh data
!
!    IF (par%master) THEN
!      mesh%nV       = mesh%nV_mem
!      mesh%nTri     = mesh%nTri_mem
!      mesh%xmin     = MINVAL( mesh%V( :,1))
!      mesh%xmax     = MAXVAL( mesh%V( :,1))
!      mesh%ymin     = MINVAL( mesh%V( :,2))
!      mesh%ymax     = MAXVAL( mesh%V( :,2))
!      mesh%tol_dist = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) * tol / 2._dp
!    END IF
!    CALL sync
!
!    ! Determine vertex and triangle domains
!    CALL partition_list( mesh%nV,   par%i, par%n, mesh%vi1, mesh%vi2)
!    CALL partition_list( mesh%nTri, par%i, par%n, mesh%ti1, mesh%ti2)
!
!    ! Calculate extra mesh data
!    CALL allocate_mesh_secondary(             mesh)    ! Adds  9 MPI windows
!    CALL calc_triangle_geometric_centres(     mesh)
!    CALL find_Voronoi_cell_areas(             mesh)
!    CALL calc_lat_lon_coordinates(            mesh)
!    CALL find_triangle_areas(                 mesh)
!    CALL find_connection_widths(              mesh)
!    CALL make_Ac_mesh(                        mesh)    ! Adds  5 MPI windows
!    CALL calc_matrix_operators_mesh_basic(    mesh)    ! Adds 42 MPI windows (6 CSR matrices, 7 windows each)
!    CALL determine_mesh_resolution(           mesh)
!    CALL find_Voronoi_cell_geometric_centres( mesh)
!
!    CALL check_mesh( mesh)
!
!    ! Give the mesh a nice name
!    mesh%name = 'mesh_from_file_"' // TRIM( filename) // '"'

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
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_zeta, id_dim_zeta, dim_length = nzeta)

    ! Inquire zeta variable
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_zeta, id_var_zeta)

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
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_x, id_dim_x, dim_length = nx)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_y, id_dim_y, dim_length = ny)

    ! Allocate memory for x and y
    ALLOCATE( x( nx))
    ALLOCATE( y( ny))

    ! Inquire x and y variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_x, id_var_x)
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_y, id_var_y)

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
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = nlon)
    CALL inquire_dim_multiple_options( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = nlat)

    ! Allocate memory for lon and lat
    ALLOCATE( lon( nlon))
    ALLOCATE( lat( nlat))

    ! Inquire lon and lon variables
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lon, id_var_lon)
    CALL inquire_var_multiple_options( filename, ncid, field_name_options_lat, id_var_lat)

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
    DEALLOCATE( lon)
    DEALLOCATE( lat)

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
    IF (is_correct) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Limit values to [0,360]
    DO i = 1, grid_lonlat%nlon
      IF     (grid_lonlat%lon( i) <   0._dp) THEN
        grid_lonlat%lon( i) = grid_lonlat%lon( i) + 360._dp
      ELSEIF (grid_lonlat%lon( i) > 360._dp) THEN
        grid_lonlat%lon( i) = grid_lonlat%lon( i) - 360._dp
      END IF
    END DO

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
      grid_lonlat%lon = [grid_lonlat%lon( n+1:grid_lonlat%nlon), grid_lonlat%lon( 1:n)]

      ! Fix data field
      DO j = 1, grid_lonlat%nlon
        d_grid( :,j) = [d_grid( n+1:grid_lonlat%nlon,j), d_grid( 1:n,j)]
      END DO

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
    IF (is_correct) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Limit values to [0,360]
    DO i = 1, grid_lonlat%nlon
      IF     (grid_lonlat%lon( i) <   0._dp) THEN
        grid_lonlat%lon( i) = grid_lonlat%lon( i) + 360._dp
      ELSEIF (grid_lonlat%lon( i) > 360._dp) THEN
        grid_lonlat%lon( i) = grid_lonlat%lon( i) - 360._dp
      END IF
    END DO

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
      grid_lonlat%lon = [grid_lonlat%lon( n+1:grid_lonlat%nlon), grid_lonlat%lon( 1:n)]

      ! Fix data field
      DO j = 1, grid_lonlat%nlat
      DO k = 1, SIZE( d_grid,3)
        d_grid( :,j,k) = [d_grid( n+1:grid_lonlat%nlon,j,k), d_grid( 1:n,j,k)]
      END DO
      END DO

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

END MODULE netcdf_input