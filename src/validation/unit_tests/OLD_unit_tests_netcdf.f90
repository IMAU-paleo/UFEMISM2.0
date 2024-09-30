MODULE unit_tests_netcdf

  ! Unit tests for different NetCDF in/output routines.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE grid_basic                                             , ONLY: type_grid, setup_square_grid, check_if_grids_are_identical
  USE grid_lonlat_basic                                      , ONLY: type_grid_lonlat, setup_simple_lonlat_grid, check_if_lonlat_grids_are_identical
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: allocate_mesh_primary, deallocate_mesh
  USE mesh_utilities                                         , ONLY: check_mesh, check_if_meshes_are_identical
  USE mesh_dummy_meshes                                      , ONLY: initialise_dummy_mesh_5
  USE mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_refinement_fun, only: mesh_add_smileyface, mesh_add_UFEMISM_letters
  use mesh_Lloyds_algorithm, only: Lloyds_algorithm_single_iteration
  USE mesh_secondary                                         , ONLY: calc_all_secondary_mesh_data
  USE mesh_operators                                         , ONLY: calc_all_matrix_operators_mesh
  use apply_maps, only: clear_all_maps_involving_this_mesh
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, open_existing_netcdf_file_for_reading, close_netcdf_file
  USE netcdf_output                                          , ONLY: setup_xy_grid_in_netcdf_file, add_field_grid_dp_2D_notime, add_month_dimension_to_file, &
                                                                     write_to_field_multopt_grid_dp_2D_notime, add_field_grid_dp_2D_monthly_notime, &
                                                                     write_to_field_multopt_grid_dp_2D_monthly_notime, add_zeta_dimension_to_file, &
                                                                     add_field_grid_dp_3D_notime, write_to_field_multopt_grid_dp_3D_notime, &
                                                                     setup_lonlat_grid_in_netcdf_file, add_field_lonlat_grid_dp_2D_notime, &
                                                                     write_to_field_multopt_lonlat_grid_dp_2D_notime, add_field_lonlat_grid_dp_2D_monthly_notime, &
                                                                     write_to_field_multopt_lonlat_grid_dp_2D_monthly_notime, add_field_lonlat_grid_dp_3D_notime, &
                                                                     write_to_field_multopt_lonlat_grid_dp_3D_notime, &
                                                                     setup_mesh_in_netcdf_file, add_field_mesh_dp_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_notime, add_field_mesh_dp_2D_monthly_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_monthly_notime, add_field_mesh_dp_3D_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_notime
  USE netcdf_input                                           , ONLY: setup_xy_grid_from_file, setup_lonlat_grid_from_file, setup_mesh_from_file, &
                                                                     read_field_from_xy_file_2D, read_field_from_xy_file_2D_monthly, read_field_from_xy_file_3D, &
                                                                     read_field_from_lonlat_file_2D, read_field_from_lonlat_file_2D_monthly, read_field_from_lonlat_file_3D, &
                                                                     read_field_from_mesh_file_2D, read_field_from_mesh_file_2D_monthly, read_field_from_mesh_file_3D, &
                                                                     read_field_from_file_2D, read_field_from_file_2D_monthly, read_field_from_file_3D, &
                                                                     setup_zeta_from_file

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_netcdf_unit_tests
    ! Run all NetCDF unit tests

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_all_netcdf_unit_tests'
    CHARACTER(LEN=256)                                 :: filename_grid_2D, filename_grid_2D_monthly, filename_grid_3D
    CHARACTER(LEN=256)                                 :: filename_lonlat_grid_2D, filename_lonlat_grid_2D_monthly, filename_lonlat_grid_3D
    TYPE(type_mesh)                                    :: mesh, mesh2
    CHARACTER(LEN=256)                                 :: filename_mesh_2D, filename_mesh_2D_monthly, filename_mesh_3D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Test reading/writing gridded data
    CALL test_netcdf_grid_in_and_output_2D(                filename_grid_2D               )
    CALL test_netcdf_grid_in_and_output_2D_monthly(        filename_grid_2D_monthly       )
    CALL test_netcdf_grid_in_and_output_3D(                filename_grid_3D               )

    ! Test reading/writing lon/lat-gridded data
    CALL test_netcdf_lonlat_grid_in_and_output_2D(         filename_lonlat_grid_2D        )
    CALL test_netcdf_lonlat_grid_in_and_output_2D_monthly( filename_lonlat_grid_2D_monthly)
    CALL test_netcdf_lonlat_grid_in_and_output_3D(         filename_lonlat_grid_3D        )

    ! Test reading/writing meshed data
    CALL create_test_mesh( mesh)
    CALL test_netcdf_mesh_in_and_output_2D(          mesh, filename_mesh_2D               )
    CALL test_netcdf_mesh_in_and_output_2D_monthly(  mesh, filename_mesh_2D_monthly       )
    CALL test_netcdf_mesh_in_and_output_3D(          mesh, filename_mesh_3D               )

    ! Test flexible input reading
    CALL create_test_mesh2( mesh2)
    CALL test_flexible_input_2D(         mesh2, filename_grid_2D        , filename_lonlat_grid_2D        , filename_mesh_2D        )
    CALL test_flexible_input_2D_monthly( mesh2, filename_grid_2D_monthly, filename_lonlat_grid_2D_monthly, filename_mesh_2D_monthly)
    CALL test_flexible_input_3D(         mesh2, filename_grid_3D        , filename_lonlat_grid_3D        , filename_mesh_3D        )
    CALL clear_all_maps_involving_this_mesh( mesh2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_netcdf_unit_tests

! == Test reading/writing gridded data

  SUBROUTINE test_netcdf_grid_in_and_output_2D( filename)
    ! Test the NetCDF routines handling gridded in/output data

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_grid_in_and_output_2D'
    TYPE(type_grid)                                    :: grid
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256)                                 :: name
    REAL(dp)                                           :: dx
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial
    INTEGER                                            :: n,i,j
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    INTEGER                                            :: ncid
    TYPE(type_grid)                                    :: grid_read
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing grid and gridded data to NetCDF
  ! ===============================================

    ! Set up a square grid
    name = 'test_grid'
    dx   = 32E3_dp
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid, lambda_M, phi_M, beta_stereo)

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n1:grid%n2))
    DO n = grid%n1, grid%n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      x = grid%x( i)
      y = grid%y( j)
      CALL test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_grid_vec_partial( n) = d
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    ! Add all the variables
    CALL add_field_grid_dp_2D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_xy_grid_from_file( filename, ncid, grid_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_grid_vec_partial_read( grid_read%n1:grid_read%n2))

    ! Read data from file
    CALL read_field_from_xy_file_2D( filename, 'd', d_grid_vec_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_grids_are_identical( grid, grid_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('grids are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the data read from the file is identical to the original
    DO n = grid%n1, grid%n2
      IF (ABS( d_grid_vec_partial_read( n) - d_grid_vec_partial( n)) > &
          ABS(1E-9_dp * MAX(d_grid_vec_partial_read(n), d_grid_vec_partial( n) ))) found_errors = .TRUE.
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF gridded 2-D in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF gridded 2-D in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_grid_in_and_output_2D

  SUBROUTINE test_netcdf_grid_in_and_output_2D_monthly( filename)
    ! Test the NetCDF routines handling gridded in/output data

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_grid_in_and_output_2D_monthly'
    TYPE(type_grid)                                    :: grid
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256)                                 :: name
    REAL(dp)                                           :: dx
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial
    INTEGER                                            :: n,i,j,m
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    INTEGER                                            :: ncid
    TYPE(type_grid)                                    :: grid_read
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing grid and gridded data to NetCDF
  ! ===============================================

    ! Set up a square grid
    name = 'test_grid'
    dx   = 32E3_dp
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid, lambda_M, phi_M, beta_stereo)

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n1:grid%n2, 12))
    DO n = grid%n1, grid%N2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      x = grid%x( i)
      y = grid%y( j)
      CALL test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      DO m = 1, 12
        d_grid_vec_partial( n,m) = d + REAL( m,dp)
      END DO
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    CALL add_month_dimension_to_file( filename, ncid)

    ! Add all the variables
    CALL add_field_grid_dp_2D_monthly_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_grid_dp_2D_monthly_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_xy_grid_from_file( filename, ncid, grid_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_grid_vec_partial_read( grid_read%n1:grid_read%n2,12))

    ! Read data from file
    CALL read_field_from_xy_file_2D_monthly( filename, 'd', d_grid_vec_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_grids_are_identical( grid, grid_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('grids are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the data read from the file is identical to the original
    DO n = grid%n1, grid%n2
    DO m = 1, 12
      IF (ABS( 1._dp - d_grid_vec_partial_read( n,m) / d_grid_vec_partial( n,m)) > 1E-9_dp) found_errors = .TRUE.
    END DO
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF gridded 2-D monthly in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF gridded 2-D monthly in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_grid_in_and_output_2D_monthly

  SUBROUTINE test_netcdf_grid_in_and_output_3D( filename)
    ! Test the NetCDF routines handling gridded in/output data

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_grid_in_and_output_3D'
    TYPE(type_grid)                                    :: grid
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256)                                 :: name
    REAL(dp)                                           :: dx
    INTEGER                                            :: nzeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: zeta
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial
    INTEGER                                            :: n,i,j,k
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    INTEGER                                            :: ncid
    TYPE(type_grid)                                    :: grid_read
    INTEGER                                            :: nzeta_read
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: zeta_read
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing grid and gridded data to NetCDF
  ! ===============================================

    ! Set up a square grid
    name = 'test_grid'
    dx   = 32E3_dp
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid, lambda_M, phi_M, beta_stereo)

    ! Set up a zeta coordinate
    nzeta = C%nz
    ALLOCATE( zeta( nzeta))
    DO k = 1, nzeta
      zeta( k) = REAL( k-1,dp) / REAL( nzeta-1,dp)
    END DO

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n1:grid%n2, nzeta))
    DO n = grid%n1, grid%n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      x = grid%x( i)
      y = grid%y( j)
      CALL test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      DO k = 1, nzeta
        d_grid_vec_partial( n,k) = d + zeta( k)
      END DO
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    CALL add_zeta_dimension_to_file( filename, ncid, zeta)

    ! Add all the variables
    CALL add_field_grid_dp_3D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_grid_dp_3D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_xy_grid_from_file( filename, ncid, grid_read)
    CALL setup_zeta_from_file( filename, ncid, nzeta_read, zeta_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_grid_vec_partial_read( grid_read%n1:grid_read%n2, nzeta_read))

    ! Read data from file
    CALL read_field_from_xy_file_3D( filename, 'd', d_grid_vec_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_grids_are_identical( grid, grid_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('grids are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the zeta read from the file is identical to the original
    IF (nzeta /= nzeta_read) found_errors = .TRUE.
    DO k = 1, MIN( nzeta, nzeta_read)
      IF (ABS( 1._dp - MAX( 0.001_dp, zeta( k)) / MAX( 0.001_dp, zeta_read( k))) > 1E-9_dp) THEN
        CALL warning('zeta grids are not identical')
        found_errors = .TRUE.
      END IF
    END DO

    ! Check if the data read from the file is identical to the original
    DO n = grid%n1, grid%n2
    DO k = 1, MIN( nzeta, nzeta_read)
      IF ( 1._dp - MAX( 0.001_dp, ABS( d_grid_vec_partial_read( n,k))) / MAX( 0.001_dp, ABS( d_grid_vec_partial( n,k))) > 1E-9_dp) THEN
        found_errors = .TRUE.
      END IF
    END DO
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF gridded 3-D in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF gridded 3-D in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_grid_in_and_output_3D

! == Test reading/writing gridded data

  SUBROUTINE test_netcdf_lonlat_grid_in_and_output_2D( filename)
    ! Test the NetCDF routines handling lon/lat-gridded in/output data

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_lonlat_grid_in_and_output_2D'
    CHARACTER(LEN=256)                                 :: name
    INTEGER                                            :: nlon, nlat
    TYPE(type_grid_lonlat)                             :: grid
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial
    INTEGER                                            :: n,i,j
    REAL(dp)                                           :: lon,lat,d
    INTEGER                                            :: ncid
    TYPE(type_grid_lonlat)                             :: grid_read
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing grid and gridded data to NetCDF
  ! ===============================================

    ! Set up a lon/lat-grid
    name = 'test_lonlat_grid'
    nlon = 360
    nlat = 180
    CALL setup_simple_lonlat_grid( name, nlon, nlat, grid)

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n1:grid%n2))
    DO n = grid%n1, grid%n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      lon = grid%lon( i)
      lat = grid%lat( j)
      CALL test_function_lonlat( lon, lat, d)
      d_grid_vec_partial( n) = d
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_lonlat_grid_in_netcdf_file( filename, ncid, grid)

    ! Add all the variables
    CALL add_field_lonlat_grid_dp_2D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_lonlat_grid_dp_2D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_lonlat_grid_from_file( filename, ncid, grid_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_grid_vec_partial_read( grid_read%n1:grid_read%n2))

    ! Read data from file
    CALL read_field_from_lonlat_file_2D( filename, 'd', d_grid_vec_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_lonlat_grids_are_identical( grid, grid_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('grids are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the data read from the file is identical to the original
    DO n = grid%n1, grid%n2
      IF ( 1._dp - MAX( 0.001_dp, ABS( d_grid_vec_partial_read( n))) / MAX( 0.001_dp, ABS( d_grid_vec_partial( n))) > 1E-7_dp) found_errors = .TRUE.
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF lon/lat-gridded 2-D in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF lon/lat-gridded 2-D in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_lonlat_grid_in_and_output_2D

  SUBROUTINE test_netcdf_lonlat_grid_in_and_output_2D_monthly( filename)
    ! Test the NetCDF routines handling lon/lat-gridded in/output data

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_lonlat_grid_in_and_output_2D_monthly'
    CHARACTER(LEN=256)                                 :: name
    INTEGER                                            :: nlon, nlat
    TYPE(type_grid_lonlat)                             :: grid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial
    INTEGER                                            :: n,i,j,m
    REAL(dp)                                           :: lon,lat,d
    INTEGER                                            :: ncid
    TYPE(type_grid_lonlat)                             :: grid_read
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing grid and gridded data to NetCDF
  ! ===============================================

    ! Set up a lon/lat-grid
    name = 'test_lonlat_grid'
    nlon = 360
    nlat = 180
    CALL setup_simple_lonlat_grid( name, nlon, nlat, grid)

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n1:grid%n2, 12))
    DO n = grid%n1, grid%n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      lon = grid%lon( i)
      lat = grid%lat( j)
      CALL test_function_lonlat( lon, lat, d)
      DO m = 1, 12
        d_grid_vec_partial( n,m) = d + REAL( m,dp)
      END DO
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_lonlat_grid_in_netcdf_file( filename, ncid, grid)
    CALL add_month_dimension_to_file( filename, ncid)

    ! Add all the variables
    CALL add_field_lonlat_grid_dp_2D_monthly_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_lonlat_grid_dp_2D_monthly_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_lonlat_grid_from_file( filename, ncid, grid_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_grid_vec_partial_read( grid_read%n1:grid_read%n2,12))

    ! Read data from file
    CALL read_field_from_lonlat_file_2D_monthly( filename, 'd', d_grid_vec_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_lonlat_grids_are_identical( grid, grid_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('grids are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the data read from the file is identical to the original
    DO n = grid%n1, grid%n2
    DO m = 1, 12
      IF ( 1._dp - MAX( 0.001_dp, ABS( d_grid_vec_partial_read( n,m))) / MAX( 0.001_dp, ABS( d_grid_vec_partial( n,m))) > 1E-7_dp) found_errors = .TRUE.
    END DO
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF lon/lat-gridded 2-D monthly in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF lon/lat-gridded 2-D monthly in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_lonlat_grid_in_and_output_2D_monthly

  SUBROUTINE test_netcdf_lonlat_grid_in_and_output_3D( filename)
    ! Test the NetCDF routines handling lon/lat-gridded in/output data

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_lonlat_grid_in_and_output_3D'
    CHARACTER(LEN=256)                                 :: name
    INTEGER                                            :: nlon, nlat
    TYPE(type_grid_lonlat)                             :: grid
    INTEGER                                            :: nzeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: zeta
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial
    INTEGER                                            :: n,i,j,k
    REAL(dp)                                           :: lon,lat,d
    INTEGER                                            :: ncid
    TYPE(type_grid_lonlat)                             :: grid_read
    INTEGER                                            :: nzeta_read
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: zeta_read
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing grid and gridded data to NetCDF
  ! ===============================================

    ! Set up a lon/lat-grid
    name = 'test_lonlat_grid'
    nlon = 360
    nlat = 180
    CALL setup_simple_lonlat_grid( name, nlon, nlat, grid)

    ! Set up a zeta coordinate
    nzeta = C%nz
    ALLOCATE( zeta( nzeta))
    DO k = 1, nzeta
      zeta( k) = REAL( k-1,dp) / REAL( nzeta-1,dp)
    END DO

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n1:grid%n2, nzeta))
    DO n = grid%n1, grid%n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      lon = grid%lon( i)
      lat = grid%lat( j)
      CALL test_function_lonlat( lon, lat, d)
      DO k = 1, nzeta
        d_grid_vec_partial( n,k) = d + zeta( k)
      END DO
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_lonlat_grid_in_netcdf_file( filename, ncid, grid)
    CALL add_zeta_dimension_to_file( filename, ncid, zeta)

    ! Add all the variables
    CALL add_field_lonlat_grid_dp_3D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_lonlat_grid_dp_3D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_lonlat_grid_from_file( filename, ncid, grid_read)
    CALL setup_zeta_from_file( filename, ncid, nzeta_read, zeta_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_grid_vec_partial_read( grid_read%n1:grid_read%n2,nzeta_read))

    ! Read grid and data from file
    CALL read_field_from_lonlat_file_3D( filename, 'd', d_grid_vec_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_lonlat_grids_are_identical( grid, grid_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('grids are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the zeta read from the file is identical to the original
    IF (nzeta /= nzeta_read) found_errors = .TRUE.
    DO k = 1, MIN( nzeta, nzeta_read)
      IF (1._dp - MAX( 0.001_dp, ABS( zeta( k))) / MAX( 0.001_dp, ABS( zeta_read( k))) > 1E-9_dp) found_errors = .TRUE.
    END DO

    ! Check if the data read from the file is identical to the original
    DO n = grid%n1, grid%n2
    DO k = 1, nzeta
      IF ( 1._dp - MAX( 0.001_dp, ABS( d_grid_vec_partial_read( n,k))) / MAX( 0.001_dp, ABS( d_grid_vec_partial( n,k))) > 1E-7_dp) found_errors = .TRUE.
    END DO
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF lon/lat-gridded 3-D in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF lon/lat-gridded 3-D in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_lonlat_grid_in_and_output_3D

! == Test reading/writing meshed data

  SUBROUTINE create_test_mesh( mesh)
    ! Create a simple test mesh to save time

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(OUT)   :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_test_mesh'
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256)                                 :: name
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Create a nice mesh, with a smileyface and the UFEMISM letters
  ! ================================================================

    ! Allocate memory
    name = 'test_mesh'
    CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Initialise the dummy mesh
    CALL initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh with a uniform 400 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 400E3_dp
    CALL refine_mesh_uniform( mesh, res_max, alpha_min)

    ! Smooth the mesh by applying a single iteration of Lloyd's algorithm
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)

    ! Add a smileyface
    res   = 80E3_dp
    width = 100E3_dp
    CALL mesh_add_smileyface( mesh, res, width)

    ! Add the UFEMISM letters
    res   = 50E3_dp
    width = 40E3_dp
    CALL mesh_add_UFEMISM_letters( mesh, res, width)

    ! Smooth the mesh again
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_test_mesh

  SUBROUTINE test_netcdf_mesh_in_and_output_2D( mesh, filename)
    ! Test the NetCDF routines handling mesh in/output data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_mesh_in_and_output_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial
    INTEGER                                            :: vi
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    INTEGER                                            :: ncid
    TYPE(type_mesh)                                    :: mesh_read
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing mesh and meshed data to NetCDF
  ! ===============================================

    ! Generate a simple test data field
    ALLOCATE( d_mesh_partial( mesh%vi1:mesh%vi2))
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_mesh_partial( vi) = d
    END DO

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd', d_mesh_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading mesh and meshed data from NetCDF
  ! =================================================

    ! Set up the mesh from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_mesh_from_file( filename, ncid, mesh_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_mesh_partial_read( mesh_read%vi1:mesh_read%vi2))

    ! Read grid and data from file
    CALL read_field_from_mesh_file_2D( filename, 'd', d_mesh_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_meshes_are_identical( mesh, mesh_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('meshes are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the data read from the file is identical to the original
    DO vi = mesh%vi1, mesh%vi2
      IF (1._dp - MAX( 0.001_dp, ABS( d_mesh_partial_read( vi))) / MAX( 0.001_dp, ABS( d_mesh_partial( vi))) > 1E-9_dp) found_errors = .TRUE.
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF meshed 2-D in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF meshed 2-D in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_mesh_in_and_output_2D

  SUBROUTINE test_netcdf_mesh_in_and_output_2D_monthly( mesh, filename)
    ! Test the NetCDF routines handling mesh in/output data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_mesh_in_and_output_2D_monthly'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_mesh_partial
    INTEGER                                            :: vi,m
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    INTEGER                                            :: ncid
    TYPE(type_mesh)                                    :: mesh_read
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_mesh_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing mesh and meshed data to NetCDF
  ! ===============================================

    ! Generate a simple test data field
    ALLOCATE( d_mesh_partial( mesh%vi1:mesh%vi2, 12))
    DO vi = mesh%vi1,mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      DO m = 1, 12
        d_mesh_partial( vi,m) = d + REAL( m,dp)
      END DO
    END DO

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL add_month_dimension_to_file( filename, ncid)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_monthly_notime( mesh, filename, ncid, 'd', d_mesh_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_mesh_from_file( filename, ncid, mesh_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_mesh_partial_read( mesh_read%vi1:mesh_read%vi2,12))

    ! Read grid and data from file
    CALL read_field_from_mesh_file_2D_monthly( filename, 'd', d_mesh_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_meshes_are_identical( mesh, mesh_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('meshes are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the data read from the file is identical to the original
    DO vi = mesh%vi1, mesh%vi2
    DO m = 1, 12
      IF ( 1._dp - MAX( 0.001_dp, ABS( d_mesh_partial_read( vi,m))) / MAX( 0.001_dp, ABS( d_mesh_partial( vi,m))) > 1E-9_dp) found_errors = .TRUE.
    END DO
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF meshed 2-D monthly in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF meshed 2-D monthly in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_mesh_in_and_output_2D_monthly

  SUBROUTINE test_netcdf_mesh_in_and_output_3D( mesh, filename)
    ! Test the NetCDF routines handling mesh in/output data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=256)                 , INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_netcdf_mesh_in_and_output_3D'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_mesh_partial
    INTEGER                                            :: vi,k
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    INTEGER                                            :: ncid
    TYPE(type_mesh)                                    :: mesh_read
    INTEGER                                            :: nzeta_read
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: zeta_read
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_mesh_partial_read
    LOGICAL                                            :: found_errors
    LOGICAL                                            :: are_identical

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Test writing mesh and meshed data to NetCDF
  ! ===============================================

    ! Generate a simple test data field
    ALLOCATE( d_mesh_partial( mesh%vi1:mesh%vi2, mesh%nz))
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      DO k = 1, mesh%nz
        d_mesh_partial( vi,k) = d + mesh%zeta( k)
      END DO
    END DO

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL add_zeta_dimension_to_file( filename, ncid, mesh%zeta)

    ! Add all the variables
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_3D_notime( mesh, filename, ncid, 'd', d_mesh_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Set up the grid from the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)
    CALL setup_mesh_from_file( filename, ncid, mesh_read)
    CALL setup_zeta_from_file( filename, ncid, nzeta_read, zeta_read)
    CALL close_netcdf_file( ncid)

    ! Allocate memory for read data
    ALLOCATE( d_mesh_partial_read( mesh_read%vi1:mesh_read%vi2,mesh%nz))

    ! Read grid and data from file
    CALL read_field_from_mesh_file_3D( filename, 'd', d_mesh_partial_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_meshes_are_identical( mesh, mesh_read, are_identical)
    IF (.NOT. are_identical) THEN
      CALL warning('meshes are not identical')
      found_errors = .TRUE.
    END IF

    ! Check if the zeta read from the file is identical to the original
    IF (mesh%nz /= nzeta_read) found_errors = .TRUE.
    DO k = 1, MIN( mesh%nz, nzeta_read)
      IF ( 1._dp - MAX( 0.001_dp, ABS( mesh%zeta( k))) / MAX( 0.001_dp, ABS( zeta_read( k))) > 1E-9_dp) found_errors = .TRUE.
    END DO

    ! Check if the data read from the file is identical to the original
    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, MIN( mesh%nz, nzeta_read)
      IF ( 1._dp - MAX( 0.001_dp, ABS( d_mesh_partial_read( vi,k))) / MAX( 0.001_dp, ABS( d_mesh_partial( vi,k))) > 1E-9_dp) found_errors = .TRUE.
    END DO
    END DO

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all NetCDF meshed 3-D in/output routines')
    ELSE
      IF (par%master) CALL warning('found errors in NetCDF meshed 3-D in/output routines')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_netcdf_mesh_in_and_output_3D

! == Test flexible input reading

  SUBROUTINE create_test_mesh2( mesh)
    ! Create another simple test mesh to save time

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(OUT)   :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_test_mesh2'
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256)                                 :: name
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Create a nice mesh, with a smileyface and the UFEMISM letters
  ! ================================================================

    ! Allocate memory
    name = 'test_mesh_new'
    CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Initialise the dummy mesh
    CALL initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh with a uniform 400 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 300E3_dp
    CALL refine_mesh_uniform( mesh, res_max, alpha_min)

    ! Smooth the mesh by applying a single iteration of Lloyd's algorithm
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)

    ! Add a smileyface
    res   = 80E3_dp
    width = 100E3_dp
    CALL mesh_add_smileyface( mesh, res, width)

    ! Smooth the mesh again
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_test_mesh2

  SUBROUTINE test_flexible_input_2D( mesh, filename_grid, filename_lonlat_grid, filename_mesh)
    ! Test the flexible input reading routines

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=256)                 , INTENT(IN)    :: filename_grid, filename_lonlat_grid, filename_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_flexible_input_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_from_grid, d_from_lonlat_grid, d_from_mesh
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    LOGICAL                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( d_from_grid(        mesh%vi1:mesh%vi2))
    ALLOCATE( d_from_lonlat_grid( mesh%vi1:mesh%vi2))
    ALLOCATE( d_from_mesh(        mesh%vi1:mesh%vi2))

    ! Read data from files
    CALL read_field_from_file_2D( filename_grid       , 'd', mesh, d_from_grid       )
    CALL read_field_from_file_2D( filename_lonlat_grid, 'd', mesh, d_from_lonlat_grid)
    CALL read_field_from_file_2D( filename_mesh       , 'd', mesh, d_from_mesh       )

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_from_grid'       )
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_from_lonlat_grid')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_from_mesh'       )

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_from_grid'       , d_from_grid       )
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_from_lonlat_grid', d_from_lonlat_grid)
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_from_mesh'       , d_from_mesh       )

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Validation
  ! =============

    ! No validation so far - if it doesn't crash then it's good enough

    found_errors = .FALSE.

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all flexible 2-D input routines')
    ELSE
      IF (par%master) CALL warning('found errors in flexible 2-D input routines')
    END IF

    ! Clean up after yourself
    DEALLOCATE( d_from_grid)
    DEALLOCATE( d_from_lonlat_grid)
    DEALLOCATE( d_from_mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_flexible_input_2D

  SUBROUTINE test_flexible_input_2D_monthly( mesh, filename_grid, filename_lonlat_grid, filename_mesh)
    ! Test the flexible input reading routines

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=256)                 , INTENT(IN)    :: filename_grid, filename_lonlat_grid, filename_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_flexible_input_2D_monthly'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_from_grid, d_from_lonlat_grid, d_from_mesh
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    LOGICAL                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( d_from_grid(        mesh%vi1:mesh%vi2,12))
    ALLOCATE( d_from_lonlat_grid( mesh%vi1:mesh%vi2,12))
    ALLOCATE( d_from_mesh(        mesh%vi1:mesh%vi2,12))

    ! Read data from files
    CALL read_field_from_file_2D_monthly( filename_grid       , 'd', mesh, d_from_grid       )
    CALL read_field_from_file_2D_monthly( filename_lonlat_grid, 'd', mesh, d_from_lonlat_grid)
    CALL read_field_from_file_2D_monthly( filename_mesh       , 'd', mesh, d_from_mesh       )

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL add_month_dimension_to_file( filename, ncid)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'd_from_grid'       )
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'd_from_lonlat_grid')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'd_from_mesh'       )

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_monthly_notime( mesh, filename, ncid, 'd_from_grid'       , d_from_grid       )
    CALL write_to_field_multopt_mesh_dp_2D_monthly_notime( mesh, filename, ncid, 'd_from_lonlat_grid', d_from_lonlat_grid)
    CALL write_to_field_multopt_mesh_dp_2D_monthly_notime( mesh, filename, ncid, 'd_from_mesh'       , d_from_mesh       )

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Validation
  ! =============

    ! No validation so far - if it doesn't crash then it's good enough

    found_errors = .FALSE.

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all flexible 2-D monthly input routines')
    ELSE
      IF (par%master) CALL warning('found errors in flexible 2-D monthly input routines')
    END IF

    ! Clean up after yourself
    DEALLOCATE( d_from_grid)
    DEALLOCATE( d_from_lonlat_grid)
    DEALLOCATE( d_from_mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_flexible_input_2D_monthly

  SUBROUTINE test_flexible_input_3D( mesh, filename_grid, filename_lonlat_grid, filename_mesh)
    ! Test the flexible input reading routines

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=256)                 , INTENT(IN)    :: filename_grid, filename_lonlat_grid, filename_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_flexible_input_3D'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_from_grid, d_from_lonlat_grid, d_from_mesh
    INTEGER                                            :: nzeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: zeta
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    LOGICAL                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( d_from_grid(        mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( d_from_lonlat_grid( mesh%vi1:mesh%vi2,mesh%nz))
    ALLOCATE( d_from_mesh(        mesh%vi1:mesh%vi2,mesh%nz))

    ! Read data from files
    CALL read_field_from_file_3D( filename_grid       , 'd', mesh, d_from_grid       , nzeta = nzeta, zeta = zeta)
    CALL read_field_from_file_3D( filename_lonlat_grid, 'd', mesh, d_from_lonlat_grid)
    CALL read_field_from_file_3D( filename_mesh       , 'd', mesh, d_from_mesh       )

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL add_zeta_dimension_to_file( filename, ncid, zeta)

    ! Add all the variables
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'd_from_grid'       )
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'd_from_lonlat_grid')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'd_from_mesh'       )

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_3D_notime( mesh, filename, ncid, 'd_from_grid'       , d_from_grid       )
    CALL write_to_field_multopt_mesh_dp_3D_notime( mesh, filename, ncid, 'd_from_lonlat_grid', d_from_lonlat_grid)
    CALL write_to_field_multopt_mesh_dp_3D_notime( mesh, filename, ncid, 'd_from_mesh'       , d_from_mesh       )

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Validation
  ! =============

    ! No validation so far - if it doesn't crash then it's good enough

    found_errors = .FALSE.

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all flexible 3-D input routines')
    ELSE
      IF (par%master) CALL warning('found errors in flexible 3-D input routines')
    END IF

    ! Clean up after yourself
    DEALLOCATE( d_from_grid)
    DEALLOCATE( d_from_lonlat_grid)
    DEALLOCATE( d_from_mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_flexible_input_3D

! == Test functions for generating dummy data fields

  SUBROUTINE test_function( x, y, xmin, xmax, ymin, ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
    ! A simple test function to validate matrix operators

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                   INTENT(IN)        :: x,y
    REAL(dp),                   INTENT(IN)        :: xmin,xmax,ymin,ymax
    REAL(dp),                   INTENT(OUT)       :: d,ddx,ddy
    REAL(dp),                   INTENT(OUT)       :: d2dx2,d2dxdy,d2dy2

    ! Local variables:
    REAL(dp)                                      :: c1,c2

    c1 = 2._dp * pi / (xmax - xmin)
    c2 = 3._dp * pi / (ymax - ymin)

    d      =            SIN( c1 * (x - xmin)) *            SIN( c2 * (y - ymin))
    ddx    =   c1     * COS( c1 * (x - xmin)) *            SIN( c2 * (y - ymin))
    ddy    =            SIN( c1 * (x - xmin)) *   c2     * COS( c2 * (y - ymin))
    d2dx2  = (-c1)**2 * SIN( c1 * (x - xmin)) *            SIN( c2 * (y - ymin))
    d2dxdy =   c1     * COS( c1 * (x - xmin)) *   c2     * COS( c2 * (y - ymin))
    d2dy2  =            SIN( c1 * (x - xmin)) * (-c2)**2 * SIN( c2 * (y - ymin))

  END SUBROUTINE test_function

  SUBROUTINE test_function_lonlat( lon, lat, d)
    ! A simple test function to validate lon/lat gridding stuff

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                   INTENT(IN)        :: lon,lat
    REAL(dp),                   INTENT(OUT)       :: d

    d = SIN( lon * pi / 180._dp) * COS( lat * pi / 180)

  END SUBROUTINE test_function_lonlat

END MODULE unit_tests_netcdf
