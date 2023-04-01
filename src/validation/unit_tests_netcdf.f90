MODULE unit_tests_netcdf

  ! Unit tests for different NetCDF in/output routines.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE parameters
  USE grid_basic                                             , ONLY: type_grid, setup_square_grid, check_if_grids_are_identical
  USE mesh_types                                             , ONLY: type_mesh
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file
  USE netcdf_output                                          , ONLY: setup_xy_grid_in_netcdf_file, add_field_grid_dp_2D_notime, add_month_dimension_to_file, &
                                                                     write_to_field_multiple_options_grid_dp_2D_notime, add_field_grid_dp_2D_monthly_notime, &
                                                                     write_to_field_multiple_options_grid_dp_2D_monthly_notime, add_zeta_dimension_to_file, &
                                                                     add_field_grid_dp_3D_notime, write_to_field_multiple_options_grid_dp_3D_notime
  USE netcdf_input                                           , ONLY: read_field_from_xy_file_2D, read_field_from_xy_file_2D_monthly, read_field_from_xy_file_3D

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

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all NetCDF unit tests
    CALL test_netcdf_grid_in_and_output_2D(         filename_grid_2D        )
    CALL test_netcdf_grid_in_and_output_2D_monthly( filename_grid_2D_monthly)
    CALL test_netcdf_grid_in_and_output_3D(         filename_grid_3D        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_netcdf_unit_tests

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
    INTEGER                                            :: n,n_glob,i,j
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
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, lambda_M, phi_M, beta_stereo, grid)

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n_loc))
    DO n = 1, grid%n_loc
      n_glob = grid%n1 + n - 1
      i = grid%n2ij( n_glob,1)
      j = grid%n2ij( n_glob,2)
      x = grid%x( i)
      y = grid%y( j)
      CALL test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_grid_vec_partial( n) = d
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    ! Add all the variables
    CALL add_field_grid_dp_2D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multiple_options_grid_dp_2D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Read grid and data from file
    CALL read_field_from_xy_file_2D( filename, 'd', d_grid_vec_partial_read, grid = grid_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_grids_are_identical( grid, grid_read, are_identical)
    found_errors = found_errors .OR. (.NOT. are_identical)

    ! Check if the data read from the file is identical to the original
    DO n = 1, grid%n_loc
      IF (ABS( 1._dp - d_grid_vec_partial_read( n) / d_grid_vec_partial( n)) > 1E-9_dp) found_errors = .TRUE.
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
    INTEGER                                            :: n,n_glob,i,j,m
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
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, lambda_M, phi_M, beta_stereo, grid)

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n_loc, 12))
    DO n = 1, grid%n_loc
      n_glob = grid%n1 + n - 1
      i = grid%n2ij( n_glob,1)
      j = grid%n2ij( n_glob,2)
      x = grid%x( i)
      y = grid%y( j)
      CALL test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      DO m = 1, 12
        d_grid_vec_partial( n,m) = d + REAL( m,dp)
      END DO
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    CALL add_month_dimension_to_file( filename, ncid)

    ! Add all the variables
    CALL add_field_grid_dp_2D_monthly_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multiple_options_grid_dp_2D_monthly_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Read grid and data from file
    CALL read_field_from_xy_file_2D_monthly( filename, 'd', d_grid_vec_partial_read, grid = grid_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_grids_are_identical( grid, grid_read, are_identical)
    found_errors = found_errors .OR. (.NOT. are_identical)

    ! Check if the data read from the file is identical to the original
    DO n = 1, grid%n_loc
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
    INTEGER                                            :: n,n_glob,i,j,k
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
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, lambda_M, phi_M, beta_stereo, grid)

    ! Set up a zeta coordinate
    nzeta = 11
    ALLOCATE( zeta( 11))
    zeta = [0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp, 0.7_dp, 0.8_dp, 0.9_dp, 1.0_dp]

    ! Generate a simple test data field
    ALLOCATE( d_grid_vec_partial( grid%n_loc, nzeta))
    DO n = 1, grid%n_loc
      n_glob = grid%n1 + n - 1
      i = grid%n2ij( n_glob,1)
      j = grid%n2ij( n_glob,2)
      x = grid%x( i)
      y = grid%y( j)
      CALL test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      DO k = 1, nzeta
        d_grid_vec_partial( n,k) = d + zeta( k)
      END DO
    END DO

    ! Create a file and write the grid to it
    filename = TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    CALL add_zeta_dimension_to_file( filename, ncid, zeta)

    ! Add all the variables
    CALL add_field_grid_dp_3D_notime( filename, ncid, 'd')

    ! Write all the variables
    CALL write_to_field_multiple_options_grid_dp_3D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Test reading grid and gridded data from NetCDF
  ! =================================================

    ! Read grid and data from file
    CALL read_field_from_xy_file_3D( filename, 'd', d_grid_vec_partial_read, grid = grid_read, nzeta = nzeta_read, zeta = zeta_read)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Check if the grid read from the file is identical to the original
    CALL check_if_grids_are_identical( grid, grid_read, are_identical)
    found_errors = found_errors .OR. (.NOT. are_identical)

    ! Check if the zeta read from the file is identical to the original
    IF (nzeta /= nzeta_read) found_errors = .TRUE.
    DO k = 1, MIN( nzeta, nzeta_read)
      IF (ABS( 1._dp - zeta( k) / zeta_read( k)) > 1E-9_dp) found_errors = .TRUE.
    END DO

    ! Check if the data read from the file is identical to the original
    DO n = 1, grid%n_loc
    DO k = 1, MIN( nzeta, nzeta_read)
      IF (ABS( 1._dp - d_grid_vec_partial_read( n,k) / d_grid_vec_partial( n,k)) > 1E-9_dp) found_errors = .TRUE.
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

END MODULE unit_tests_netcdf
