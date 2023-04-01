MODULE unit_tests_mesh

  ! Unit tests for different mesh creation / discretisation stuff.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: allocate_mesh_primary, deallocate_mesh
  USE mesh_utilities                                         , ONLY: check_mesh
  USE mesh_creation                                          , ONLY: initialise_dummy_mesh
  USE mesh_refinement                                        , ONLY: refine_mesh_uniform, Lloyds_algorithm_single_iteration, mesh_add_smileyface, &
                                                                     mesh_add_UFEMISM_letters
  USE mesh_parallel_creation                                 , ONLY: merge_submeshes, broadcast_merged_mesh
  USE mesh_secondary                                         , ONLY: calc_all_secondary_mesh_data
  USE mesh_operators                                         , ONLY: calc_all_matrix_operators_mesh
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file
  USE netcdf_output                                          , ONLY: setup_mesh_in_netcdf_file, add_field_mesh_dp_2D_notime, add_field_mesh_dp_2D_b_notime, &
                                                                     add_field_mesh_dp_2D_c_notime, write_to_field_multiple_options_mesh_dp_2D_notime, &
                                                                     write_to_field_multiple_options_mesh_dp_2D_b_notime, write_to_field_multiple_options_mesh_dp_2D_c_notime, &
                                                                     setup_xy_grid_in_netcdf_file, add_field_grid_dp_2D_notime, write_to_field_multiple_options_grid_dp_2D_notime
  USE petsc_basic                                            , ONLY: multiply_CSR_matrix_with_vector_1D
  USE grid_basic                                             , ONLY: type_grid, setup_square_grid, distribute_gridded_data_from_master_dp_2D
  USE mesh_remapping                                         , ONLY: map_from_xy_grid_to_mesh_2D, map_from_mesh_to_xy_grid_2D

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  LOGICAL :: do_write_results_to_netcdf = .TRUE.

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_mesh_unit_tests
    ! Run all mesh unit tests

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_all_mesh_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all mesh unit tests
    CALL test_mesh_creation_basic_single_core
    CALL test_mesh_creation_basic_two_cores
    CALL test_mesh_operators_basic
    CALL test_remapping_grid2mesh
    CALL test_remapping_mesh2grid

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_mesh_unit_tests

  SUBROUTINE test_mesh_creation_basic_single_core
    ! Test creation of a very simple mesh without parallelised mesh generation.

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_mesh_creation_basic_single_core'
    REAL(dp)                                           :: tstart, tcomp
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256), PARAMETER                      :: name = 'test_mesh'
    TYPE(type_mesh)                                    :: mesh
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width
    INTEGER                                            :: vi, ci, vj
    REAL(dp)                                           :: RA, RA_max
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Create a nice mesh, with a smileyface and the UFEMISM letters
  ! ================================================================

    tcomp = 0._dp

    ! Allocate memory
    CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Initialise the dummy mesh
    CALL initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)
    CALL check_mesh( mesh)

    ! Refine the mesh with a uniform 400 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 400E3_dp
    tstart = MPI_WTIME()
    CALL refine_mesh_uniform( mesh, res_max, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Smooth the mesh by applying a single iteration of Lloyd's algorithm
    tstart = MPI_WTIME()
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Add a smileyface
    res   = 80E3_dp
    width = 100E3_dp
    tstart = MPI_WTIME()
    CALL mesh_add_smileyface( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Add the UFEMISM letters
    res   = 50E3_dp
    width = 40E3_dp
    tstart = MPI_WTIME()
    CALL mesh_add_UFEMISM_letters( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Smooth the mesh again
    tstart = MPI_WTIME()
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    tstart = MPI_WTIME()
    CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    tcomp = tcomp + MPI_WTIME() - tstart

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! 2023-03-23: mesh has 6394 vertices
    IF (mesh%nV < 5000 .OR. mesh%nV > 7500) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexepcted amount of vertices! Expected 6394, found {int_01}', int_01 = mesh%nV)
    END IF

    ! 2023-03-23: mesh has 12682 triangles
    IF (mesh%nTri < 10000 .OR. mesh%nTri > 15000) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexepcted amount of triangles! Expected 12682, found {int_01}', int_01 = mesh%nTri)
    END IF

    ! 2023-03-23: mesh has 19075 edges
    IF (mesh%nE < 15000 .OR. mesh%nE > 25000) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexepcted amount of edges! Expected 19075, found {int_01}', int_01 = mesh%nE)
    END IF

    ! 2023-03-23: smallest vertex has a Voronoi cell area of 0.64840E+09
    IF (MINVAL( mesh%A) < 0.3E09) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = 0.64840E+09, found {dp_01}', dp_01 = MINVAL( mesh%A))
    END IF

    ! 2023-03-23: largest vertex has a Voronoi cell area of 0.83562E+11
    IF (MAXVAL( mesh%A) > 1.6E12) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexpectedly large vertices! Expected MAXVAL( mesh%A) = 0.83562E+11, found {dp_01}', dp_01 = MAXVAL( mesh%A))
    END IF

    ! 2023-03-23: mesh has a maximum ratio of Voronoi cell area A of adjacent vertices of 5.719
    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      END DO
    END DO
    IF (RA_max > 7._dp)  THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexpectedly high resolution gradient! Expected RA_max = 0.57190E+01, found {dp_01}', dp_01 = RA_max)
    END IF

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('created a test mesh with {int_01} vertices, {int_02} triangles, and {int_03} edges in {dp_01} seconds', &
        int_01 = mesh%nV, int_02 = mesh%nTri, int_03 = mesh%nE, dp_01 = tcomp)
    ELSE
      IF (par%master) CALL warning('found errors in basic single-core mesh creation')
    END IF

    ! Write the resulting mesh to a NetCDF file
    IF (do_write_results_to_netcdf) THEN
      filename = TRIM( routine_name) // '_output.nc'
      CALL create_new_netcdf_file_for_writing( filename, ncid)
      CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
      CALL close_netcdf_file( ncid)
    END IF

    ! Clean up after yourself
    CALL deallocate_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_creation_basic_single_core

  SUBROUTINE test_mesh_creation_basic_two_cores
    ! Test creation of a very simple mesh without parallelised mesh generation.

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_mesh_creation_basic_two_cores'
    REAL(dp)                                           :: tstart, tcomp
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256), PARAMETER                      :: name = 'test_mesh'
    TYPE(type_mesh)                                    :: mesh
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width
    INTEGER                                            :: vi, ci, vj
    REAL(dp)                                           :: RA, RA_max
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Create a nice mesh on each process, with a smileyface and the UFEMISM letters
  ! ================================================================================

    tcomp = 0._dp

    ! Allocate memory
    CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Safety - must be run on at least two cores!
    IF (par%n < 2) CALL crash('this validation check needs to be run on at least two cores!')

    ! Initialise the dummy mesh
    IF (par%i == 0) THEN
      CALL initialise_dummy_mesh( mesh, xmin, (xmin + xmax) / 2._dp, ymin, ymax)
    ELSE
      CALL initialise_dummy_mesh( mesh, (xmin + xmax) / 2._dp, xmax, ymin, ymax)
    END IF
    CALL check_mesh( mesh)

    ! Refine the mesh with a uniform 400 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 400E3_dp
    tstart = MPI_WTIME()
    CALL refine_mesh_uniform( mesh, res_max, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Smooth the mesh by applying a single iteration of Lloyd's algorithm
    tstart = MPI_WTIME()
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Add a smileyface
    res   = 80E3_dp
    width = 100E3_dp
    tstart = MPI_WTIME()
    CALL mesh_add_smileyface( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Add the UFEMISM letters
    res   = 50E3_dp
    width = 40E3_dp
    tstart = MPI_WTIME()
    CALL mesh_add_UFEMISM_letters( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Smooth the mesh again
    tstart = MPI_WTIME()
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

  ! == Merge the two meshes
  ! =======================

    ! Merge submeshes
    tstart = MPI_WTIME()
    CALL merge_submeshes( mesh, 0, 1, 'east-west')
    tcomp = tcomp + MPI_WTIME() - tstart
    IF (par%master) CALL check_mesh( mesh)

    ! Smooth again
    tstart = MPI_WTIME()
    IF (par%master) CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    IF (par%master) CALL check_mesh( mesh)

    ! Broadcast from Master
    tstart = MPI_WTIME()
    CALL broadcast_merged_mesh( mesh)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    tstart = MPI_WTIME()
    CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    tcomp = tcomp + MPI_WTIME() - tstart

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! 2023-03-23: mesh has 8934 vertices
    IF (mesh%nV < 8000 .OR. mesh%nV > 12000) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexepcted amount of vertices! Expected 8934, found {int_01}', int_01 = mesh%nV)
    END IF

    ! 2023-03-23: mesh has 17737 triangles
    IF (mesh%nTri < 15000 .OR. mesh%nTri > 20000) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexepcted amount of triangles! Expected 17737, found {int_01}', int_01 = mesh%nTri)
    END IF

    ! 2023-03-23: mesh has 26670 edges
    IF (mesh%nE < 22000 .OR. mesh%nE > 30000) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexepcted amount of edges! Expected 19075, found {int_01}', int_01 = mesh%nE)
    END IF

    ! 2023-03-23: smallest vertex has a Voronoi cell area of 0.62425E+09
    IF (MINVAL( mesh%A) < 0.3E09) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = 0.62425E+09, found {dp_01}', dp_01 = MINVAL( mesh%A))
    END IF

    ! 2023-03-23: largest vertex has a Voronoi cell area of 0.74116E+11
    IF (MAXVAL( mesh%A) > 1.6E12) THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexpectedly large vertices! Expected MAXVAL( mesh%A) = 0.83562E+11, found {dp_01}', dp_01 = MAXVAL( mesh%A))
    END IF

    ! 2023-03-23: mesh has a maximum ratio of Voronoi cell area A of adjacent vertices of 0.41585E+01
    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      END DO
    END DO
    IF (RA_max > 7._dp)  THEN
      found_errors = .TRUE.
      CALL warning('mesh has unexpectedly high resolution gradient! Expected RA_max = 0.41585E+01, found {dp_01}', dp_01 = RA_max)
    END IF

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('created a test mesh with {int_01} vertices, {int_02} triangles, and {int_03} edges in {dp_01} seconds', &
        int_01 = mesh%nV, int_02 = mesh%nTri, int_03 = mesh%nE, dp_01 = tcomp)
    ELSE
      IF (par%master) CALL warning('found errors in basic two-core mesh creation')
    END IF

    ! Write the resulting mesh to a NetCDF file
    IF (do_write_results_to_netcdf) THEN
      filename = TRIM( routine_name) // '_output.nc'
      CALL create_new_netcdf_file_for_writing( filename, ncid)
      CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
      CALL close_netcdf_file( ncid)
    END IF

    ! Clean up after yourself
    CALL deallocate_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_creation_basic_two_cores

  SUBROUTINE test_mesh_operators_basic
    ! Test the basic mapping/gradient matrix operators

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_mesh_operators_basic'
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256), PARAMETER                      :: name = 'test_mesh'
    TYPE(type_mesh)                                    :: mesh
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width
    INTEGER                                            :: vi,ti,ei,row
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_a_ex, ddx_a_ex, ddy_a_ex
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_c_ex, ddx_c_ex, ddy_c_ex
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            ::        d_a_b, d_a_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_b_a,        d_b_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_c_a, d_c_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddx_a_a, ddx_a_b, ddx_a_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddx_b_a, ddx_b_b, ddx_b_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddx_c_a, ddx_c_b, ddx_c_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddy_a_a, ddy_a_b, ddy_a_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddy_b_a, ddy_b_b, ddy_b_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddy_c_a, ddy_c_b, ddy_c_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d2dx2_b_b, d2dxdy_b_b, d2dy2_b_b
    REAL(dp)                                           ::               maxerr_d_a_b, maxerr_d_a_c
    REAL(dp)                                           :: maxerr_d_b_a,               maxerr_d_b_c
    REAL(dp)                                           :: maxerr_d_c_a, maxerr_d_c_b
    REAL(dp)                                           :: maxerr_ddx_a_a, maxerr_ddx_a_b, maxerr_ddx_a_c
    REAL(dp)                                           :: maxerr_ddx_b_a, maxerr_ddx_b_b, maxerr_ddx_b_c
    REAL(dp)                                           :: maxerr_ddx_c_a, maxerr_ddx_c_b, maxerr_ddx_c_c
    REAL(dp)                                           :: maxerr_ddy_a_a, maxerr_ddy_a_b, maxerr_ddy_a_c
    REAL(dp)                                           :: maxerr_ddy_b_a, maxerr_ddy_b_b, maxerr_ddy_b_c
    REAL(dp)                                           :: maxerr_ddy_c_a, maxerr_ddy_c_b, maxerr_ddy_c_c
    REAL(dp)                                           :: maxerr_d2dx2_b_b, maxerr_d2dxdy_b_b, maxerr_d2dy2_b_b
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Create a nice mesh, with a smileyface and the UFEMISM letters
  ! ================================================================

    ! Allocate memory
    CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Initialise the dummy mesh
    CALL initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)

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

  ! == Calculate and validate all matrix operators
  ! ==============================================

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( mesh)

    ! Allocate distributed memory
    ALLOCATE( d_a_ex(      mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddx_a_ex(    mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddy_a_ex(    mesh%vi1:mesh%vi2), source = 0._dp)

    ALLOCATE( d_b_ex(      mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddx_b_ex(    mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddy_b_ex(    mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( d2dx2_b_ex(  mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( d2dxdy_b_ex( mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( d2dy2_b_ex(  mesh%ti1:mesh%ti2), source = 0._dp)

    ALLOCATE( d_c_ex(      mesh%ei1:mesh%ei2), source = 0._dp)
    ALLOCATE( ddx_c_ex(    mesh%ei1:mesh%ei2), source = 0._dp)
    ALLOCATE( ddy_c_ex(    mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( d_a_b(       mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( d_a_c(       mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( d_b_a(       mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( d_b_c(       mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( d_c_a(       mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( d_c_b(       mesh%ti1:mesh%ti2), source = 0._dp)

    ALLOCATE( ddx_a_a(     mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddx_a_b(     mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddx_a_c(     mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( ddx_b_a(     mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddx_b_b(     mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddx_b_c(     mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( ddx_c_a(     mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddx_c_b(     mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddx_c_c(     mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( ddy_a_a(     mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddy_a_b(     mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddy_a_c(     mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( ddy_b_a(     mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddy_b_b(     mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddy_b_c(     mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( ddy_c_a(     mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( ddy_c_b(     mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( ddy_c_c(     mesh%ei1:mesh%ei2), source = 0._dp)

    ALLOCATE( d2dx2_b_b(   mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( d2dxdy_b_b(  mesh%ti1:mesh%ti2), source = 0._dp)
    ALLOCATE( d2dy2_b_b(   mesh%ti1:mesh%ti2), source = 0._dp)

    ! Calculate exact solutions

    ! a-grid (vertices)
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      row = mesh%vi2n( vi)
      d_a_ex(   row) = d
      ddx_a_ex( row) = ddx
      ddy_a_ex( row) = ddy
    END DO

    ! b-grid (triangles)
    DO ti = mesh%ti1, mesh%ti2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      row = mesh%ti2n( ti)
      d_b_ex(      row) = d
      ddx_b_ex(    row) = ddx
      ddy_b_ex(    row) = ddy
      d2dx2_b_ex(  row) = d2dx2
      d2dxdy_b_ex( row) = d2dxdy
      d2dy2_b_ex(  row) = d2dy2
    END DO

    ! c-grid (edges)
    DO ei = mesh%ei1, mesh%ei2
      x = mesh%E( ei,1)
      y = mesh%E( ei,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      row = mesh%ei2n( ei)
      d_c_ex(   row) = d
      ddx_c_ex( row) = ddx
      ddy_c_ex( row) = ddy
    END DO

    ! Calculate discretised solutions
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_a_b    , d_a_ex, d_a_b     )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_a_c    , d_a_ex, d_a_c     )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_b_a    , d_b_ex, d_b_a     )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_b_c    , d_b_ex, d_b_c     )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_c_a    , d_c_ex, d_c_a     )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_c_b    , d_c_ex, d_c_b     )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_a    , d_a_ex, ddx_a_a   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_b    , d_a_ex, ddx_a_b   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_c    , d_a_ex, ddx_a_c   )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_a    , d_b_ex, ddx_b_a   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_b    , d_b_ex, ddx_b_b   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_c    , d_b_ex, ddx_b_c   )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_c_a    , d_c_ex, ddx_c_a   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_c_b    , d_c_ex, ddx_c_b   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_c_c    , d_c_ex, ddx_c_c   )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_a    , d_a_ex, ddy_a_a   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_b    , d_a_ex, ddy_a_b   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_c    , d_a_ex, ddy_a_c   )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_a    , d_b_ex, ddy_b_a   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_b    , d_b_ex, ddy_b_b   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_c    , d_b_ex, ddy_b_c   )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_c_a    , d_c_ex, ddy_c_a   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_c_b    , d_c_ex, ddy_c_b   )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_c_c    , d_c_ex, ddy_c_c   )

    CALL multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dx2_b_b , d_b_ex, d2dx2_b_b )
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dxdy_b_b, d_b_ex, d2dxdy_b_b)
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dy2_b_b , d_b_ex, d2dy2_b_b )

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Calculate maximum errors

    ! a-grid (vertices)

    maxerr_d_b_a   = 0._dp
    maxerr_d_c_a   = 0._dp
    maxerr_ddx_a_a = 0._dp
    maxerr_ddx_b_a = 0._dp
    maxerr_ddx_c_a = 0._dp
    maxerr_ddy_a_a = 0._dp
    maxerr_ddy_b_a = 0._dp
    maxerr_ddy_c_a = 0._dp

    DO vi = mesh%vi1, mesh%vi2

      ! Skip border vertices, gradients are not needed there anyway
      IF (mesh%VBI( vi) > 0) CYCLE

      ! Calculate errors
      maxerr_d_b_a   = MAX( maxerr_d_b_a  , ABS( d_b_a(   vi) - d_a_ex(   vi)))
      maxerr_d_c_a   = MAX( maxerr_d_c_a  , ABS( d_c_a(   vi) - d_a_ex(   vi)))

      maxerr_ddx_a_a = MAX( maxerr_ddx_a_a, ABS( ddx_a_a( vi) - ddx_a_ex( vi)))
      maxerr_ddx_b_a = MAX( maxerr_ddx_b_a, ABS( ddx_b_a( vi) - ddx_a_ex( vi)))
      maxerr_ddx_c_a = MAX( maxerr_ddx_c_a, ABS( ddx_c_a( vi) - ddx_a_ex( vi)))

      maxerr_ddy_a_a = MAX( maxerr_ddy_a_a, ABS( ddy_a_a( vi) - ddy_a_ex( vi)))
      maxerr_ddy_b_a = MAX( maxerr_ddy_b_a, ABS( ddy_b_a( vi) - ddy_a_ex( vi)))
      maxerr_ddy_c_a = MAX( maxerr_ddy_c_a, ABS( ddy_c_a( vi) - ddy_a_ex( vi)))

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Find maximum errors across all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_b_a  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_c_a  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_a_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_b_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_c_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_a_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_b_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_c_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr_d_b_a > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_map_b_a')
    END IF
    IF (maxerr_d_c_a > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_map_c_a')
    END IF

    IF (maxerr_ddx_a_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_a_a')
    END IF
    IF (maxerr_ddx_b_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_b_a')
    END IF
    IF (maxerr_ddx_c_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_c_a')
    END IF

    IF (maxerr_ddy_a_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_a_a')
    END IF
    IF (maxerr_ddy_b_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_b_a')
    END IF
    IF (maxerr_ddy_c_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_c_a')
    END IF

    ! b-grid (triangles)

    maxerr_d_a_b   = 0._dp
    maxerr_d_c_b   = 0._dp
    maxerr_ddx_a_b = 0._dp
    maxerr_ddx_b_b = 0._dp
    maxerr_ddx_c_b = 0._dp
    maxerr_ddy_a_b = 0._dp
    maxerr_ddy_b_b = 0._dp
    maxerr_ddy_c_b = 0._dp

    DO ti = mesh%ti1, mesh%ti2

      ! Skip border vertices, gradients are not needed there anyway
      IF (mesh%TriBI( ti) > 0) CYCLE

      ! Calculate errors
      maxerr_d_a_b   = MAX( maxerr_d_a_b  , ABS( d_a_b(   ti) - d_b_ex(   ti)))
      maxerr_d_c_b   = MAX( maxerr_d_c_b  , ABS( d_c_b(   ti) - d_b_ex(   ti)))

      maxerr_ddx_a_b = MAX( maxerr_ddx_a_b, ABS( ddx_a_b( ti) - ddx_b_ex( ti)))
      maxerr_ddx_b_b = MAX( maxerr_ddx_b_b, ABS( ddx_b_b( ti) - ddx_b_ex( ti)))
      maxerr_ddx_c_b = MAX( maxerr_ddx_c_b, ABS( ddx_c_b( ti) - ddx_b_ex( ti)))

      maxerr_ddy_a_b = MAX( maxerr_ddy_a_b, ABS( ddy_a_b( ti) - ddy_b_ex( ti)))
      maxerr_ddy_b_b = MAX( maxerr_ddy_b_b, ABS( ddy_b_b( ti) - ddy_b_ex( ti)))
      maxerr_ddy_c_b = MAX( maxerr_ddy_c_b, ABS( ddy_c_b( ti) - ddy_b_ex( ti)))

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Find maximum errors across all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_a_b  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_c_b  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_a_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_b_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_c_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_a_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_b_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_c_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr_d_a_b > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_map_a_b')
    END IF
    IF (maxerr_d_c_b > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_map_c_b')
    END IF

    IF (maxerr_ddx_a_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_a_b')
    END IF
    IF (maxerr_ddx_b_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_b_b')
    END IF
    IF (maxerr_ddx_c_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_c_b')
    END IF

    IF (maxerr_ddy_a_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_a_b')
    END IF
    IF (maxerr_ddy_b_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_b_b')
    END IF
    IF (maxerr_ddy_c_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_c_b')
    END IF

    ! c-grid (edges)

    maxerr_d_a_c   = 0._dp
    maxerr_d_b_c   = 0._dp
    maxerr_ddx_a_c = 0._dp
    maxerr_ddx_b_c = 0._dp
    maxerr_ddx_c_c = 0._dp
    maxerr_ddy_a_c = 0._dp
    maxerr_ddy_b_c = 0._dp
    maxerr_ddy_c_c = 0._dp

    DO ei = mesh%ei1, mesh%ei2

      ! Skip border vertices, gradients are not needed there anyway
      IF (mesh%EBI( ei) > 0) CYCLE

      ! Calculate errors
      maxerr_d_a_c   = MAX( maxerr_d_a_c  , ABS( d_a_c(   ei) - d_c_ex(   ei)))
      maxerr_d_b_c   = MAX( maxerr_d_b_c  , ABS( d_b_c(   ei) - d_c_ex(   ei)))

      maxerr_ddx_a_c = MAX( maxerr_ddx_a_c, ABS( ddx_a_c( ei) - ddx_c_ex( ei)))
      maxerr_ddx_b_c = MAX( maxerr_ddx_b_c, ABS( ddx_b_c( ei) - ddx_c_ex( ei)))
      maxerr_ddx_c_c = MAX( maxerr_ddx_c_c, ABS( ddx_c_c( ei) - ddx_c_ex( ei)))

      maxerr_ddy_a_c = MAX( maxerr_ddy_a_c, ABS( ddy_a_c( ei) - ddy_c_ex( ei)))
      maxerr_ddy_b_c = MAX( maxerr_ddy_b_c, ABS( ddy_b_c( ei) - ddy_c_ex( ei)))
      maxerr_ddy_c_c = MAX( maxerr_ddy_c_c, ABS( ddy_c_c( ei) - ddy_c_ex( ei)))

    END DO ! DO ei = mesh%ei1, mesh%ei2

    ! Find maximum errors across all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_a_c  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_b_c  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_a_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_b_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_c_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_a_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_b_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_c_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr_d_a_c > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_map_a_c')
    END IF
    IF (maxerr_d_b_c > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_map_b_c')
    END IF

    IF (maxerr_ddx_a_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_a_c')
    END IF
    IF (maxerr_ddx_b_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_b_c')
    END IF
    IF (maxerr_ddx_c_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddx_c_c')
    END IF

    IF (maxerr_ddy_a_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_a_c')
    END IF
    IF (maxerr_ddy_b_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_b_c')
    END IF
    IF (maxerr_ddy_c_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      CALL warning('inaccuracies found in M_ddy_c_c')
    END IF

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all basic mesh matrix operators')
    ELSE
      IF (par%master) CALL warning('found errors in basic mesh matrix operators')
    END IF

    ! Write results to a NetCDF file
    IF (do_write_results_to_netcdf) THEN

      ! Create a file and write the mesh to it
      filename = TRIM( routine_name) // '_output.nc'
      CALL create_new_netcdf_file_for_writing( filename, ncid)
      CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

      ! Add all the variables
      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'd_a_ex'     )
      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_a_ex'   )
      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_a_ex'   )

      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_b_ex'     )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_b_ex'   )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_b_ex'   )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dx2_b_ex' )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dxdy_b_ex')
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dy2_b_ex' )

      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'd_c_ex'     )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_c_ex'   )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_c_ex'   )

      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_a_b'      )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'd_a_c'      )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'd_b_a'      )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'd_b_c'      )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'd_c_a'      )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_c_b'      )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_a_a'    )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_a_b'    )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_a_c'    )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_b_a'    )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_b_b'    )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_b_c'    )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_c_a'    )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_c_b'    )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_c_c'    )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_a_a'    )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_a_b'    )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_a_c'    )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_b_a'    )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_b_b'    )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_b_c'    )

      CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_c_a'    )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_c_b'    )
      CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_c_c'    )

      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dx2_b_b'  )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dxdy_b_b' )
      CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dy2_b_b'  )

      ! Write all the variables
      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_a_ex'     , d_a_ex     )
      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_ex'   , ddx_a_ex   )
      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_ex'   , ddy_a_ex   )

      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_b_ex'     , d_b_ex     )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_ex'   , ddx_b_ex   )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_ex'   , ddy_b_ex   )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_ex' , d2dx2_b_ex )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_ex', d2dxdy_b_ex)
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_ex' , d2dy2_b_ex )

      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_c_ex'     , d_c_ex     )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_c_ex'   , ddx_c_ex   )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_c_ex'   , ddy_c_ex   )

      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_a_b'      , d_a_b      )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_a_c'      , d_a_c      )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_b_a'      , d_b_a      )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_b_c'      , d_b_c      )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_c_a'      , d_c_a      )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_c_b'      , d_c_b      )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_a'    , ddx_a_a    )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_a_b'    , ddx_a_b    )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_a_c'    , ddx_a_c    )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_b_a'    , ddx_b_a    )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_b'    , ddx_b_b    )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_b_c'    , ddx_b_c    )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_c_a'    , ddx_c_a    )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_c_b'    , ddx_c_b    )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_c_c'    , ddx_c_c    )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_a'    , ddy_a_a    )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_a_b'    , ddy_a_b    )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_a_c'    , ddy_a_c    )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_b_a'    , ddy_b_a    )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_b'    , ddy_b_b    )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_b_c'    , ddy_b_c    )

      CALL write_to_field_multiple_options_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_c_a'    , ddy_c_a    )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_c_b'    , ddy_c_b    )
      CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_c_c'    , ddy_c_c    )

      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_b'  , d2dx2_b_b  )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_b' , d2dxdy_b_b )
      CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_b'  , d2dy2_b_b  )

      ! Close the file
      CALL close_netcdf_file( ncid)

    END IF ! IF (do_write_results_to_netcdf) THEN

    ! Clean up after yourself
    DEALLOCATE( d_a_ex     )
    DEALLOCATE( ddx_a_ex   )
    DEALLOCATE( ddy_a_ex   )

    DEALLOCATE( d_b_ex     )
    DEALLOCATE( ddx_b_ex   )
    DEALLOCATE( ddy_b_ex   )
    DEALLOCATE( d2dx2_b_ex )
    DEALLOCATE( d2dxdy_b_ex)
    DEALLOCATE( d2dy2_b_ex )

    DEALLOCATE( d_c_ex     )
    DEALLOCATE( ddx_c_ex   )
    DEALLOCATE( ddy_c_ex   )

    DEALLOCATE( d_a_b      )
    DEALLOCATE( d_a_c      )

    DEALLOCATE( d_b_a      )
    DEALLOCATE( d_b_c      )

    DEALLOCATE( d_c_a      )
    DEALLOCATE( d_c_b      )

    DEALLOCATE( ddx_a_a    )
    DEALLOCATE( ddx_a_b    )
    DEALLOCATE( ddx_a_c    )

    DEALLOCATE( ddx_b_a    )
    DEALLOCATE( ddx_b_b    )
    DEALLOCATE( ddx_b_c    )

    DEALLOCATE( ddx_c_a    )
    DEALLOCATE( ddx_c_b    )
    DEALLOCATE( ddx_c_c    )

    DEALLOCATE( ddy_a_a    )
    DEALLOCATE( ddy_a_b    )
    DEALLOCATE( ddy_a_c    )

    DEALLOCATE( ddy_b_a    )
    DEALLOCATE( ddy_b_b    )
    DEALLOCATE( ddy_b_c    )

    DEALLOCATE( ddy_c_a    )
    DEALLOCATE( ddy_c_b    )
    DEALLOCATE( ddy_c_c    )

    DEALLOCATE( d2dx2_b_b  )
    DEALLOCATE( d2dxdy_b_b )
    DEALLOCATE( d2dy2_b_b  )

    CALL deallocate_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_operators_basic

  SUBROUTINE test_remapping_grid2mesh
    ! Test remapping from a square grid to a mesh

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_remapping_grid2mesh'
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256)                                 :: name
    TYPE(type_mesh)                                    :: mesh
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width
    REAL(dp)                                           :: dx
    TYPE(type_grid)                                    :: grid
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_ex_partial
    INTEGER                                            :: vi, vi_glob
    REAL(dp)                                           :: maxerr
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Create a nice mesh, with a smileyface and the UFEMISM letters
  ! ================================================================

    ! Allocate memory
    name = 'test_mesh'
    CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Initialise the dummy mesh
    CALL initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)

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

  ! == Set up a square grid
  ! =======================

    name = 'test_grid'
    dx   = 32E3_dp
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, lambda_M, phi_M, beta_stereo, grid)

  ! == Calculate, apply, and validate grid-to-mesh remapping operator
  ! =================================================================

    found_errors = .FALSE.

    ! Create a nice gridded data field on the master
    IF (par%master) THEN
      ALLOCATE( d_grid( grid%nx, grid%ny), source = 0._dp)
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        x = grid%x( i)
        y = grid%x( j)
        CALL test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_grid( i,j) = d
      END DO
      END DO
    END IF
    CALL sync

    ! Distribute gridded data
    ALLOCATE( d_grid_vec_partial( grid%n_loc))
    CALL distribute_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
    IF (par%master) DEALLOCATE( d_grid)

    ! Map gridded data to the mesh
    ALLOCATE( d_mesh_partial( mesh%nV_loc))
    CALL map_from_xy_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Calculate exact solution
    ALLOCATE( d_mesh_ex_partial( mesh%nV_loc))
    DO vi = 1, mesh%nV_loc
      vi_glob = mesh%vi1 + vi - 1
      x = mesh%V( vi_glob,1)
      y = mesh%V( vi_glob,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_mesh_ex_partial( vi) = d
    END DO

    ! Calculate error
    maxerr = 0._dp
    DO vi = 1, mesh%nV_loc
      ! Skip border vertices
      vi_glob = mesh%vi1 + vi - 1
      IF (mesh%VBI( vi_glob) > 0) CYCLE
      maxerr = MAX( maxerr, ABS( d_mesh_partial( vi) - d_mesh_ex_partial( vi)))
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr > 0.6E-1_dp) found_errors = .TRUE.

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated grid to mesh remapping')
    ELSE
      IF (par%master) CALL warning('found errors in grid to mesh remapping')
    END IF

    ! Write results to a NetCDF file
    IF (do_write_results_to_netcdf) THEN

      ! Create a file and write the mesh to it
      filename = TRIM( routine_name) // '_output.nc'
      CALL create_new_netcdf_file_for_writing( filename, ncid)
      CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

      ! Add all the variables
      CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh')
      CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh_ex')

      ! Write all the variables
      CALL write_to_field_multiple_options_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh'   , d_mesh_partial   )
      CALL write_to_field_multiple_options_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex_partial)

      ! Close the file
      CALL close_netcdf_file( ncid)

    END IF ! IF (do_write_results_to_netcdf) THEN

    ! Clean up after yourself
    CALL deallocate_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_remapping_grid2mesh

  SUBROUTINE test_remapping_mesh2grid
    ! Test remapping from a mesh to a square grid

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_remapping_mesh2grid'
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256)                                 :: name
    TYPE(type_mesh)                                    :: mesh
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width
    REAL(dp)                                           :: dx
    TYPE(type_grid)                                    :: grid
    INTEGER                                            :: i,j,n,n_glob
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial
    INTEGER                                            :: vi, vi_glob
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial, d_grid_ex_vec_partial
    REAL(dp)                                           :: maxerr
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Create a nice mesh, with a smileyface and the UFEMISM letters
  ! ================================================================

    ! Allocate memory
    name = 'test_mesh'
    CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Initialise the dummy mesh
    CALL initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)

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

  ! == Set up a square grid
  ! =======================

    name = 'test_grid'
    dx   = 32E3_dp
    CALL setup_square_grid( name, xmin, xmax, ymin, ymax, dx, lambda_M, phi_M, beta_stereo, grid)

  ! == Calculate, apply, and validate grid-to-mesh remapping operator
  ! =================================================================

    found_errors = .FALSE.

    ! Create a nice data field on the mesh
    ALLOCATE( d_mesh_partial( mesh%nV_loc), source = 0._dp)
    DO vi = 1, mesh%nV_loc
      vi_glob = mesh%vi1 + vi - 1
      x = mesh%V( vi_glob,1)
      y = mesh%V( vi_glob,2)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_mesh_partial( vi) = d
    END DO

    ! Map meshed data to the grid
    ALLOCATE( d_grid_vec_partial( grid%n_loc))
    CALL map_from_mesh_to_xy_grid_2D( mesh, grid, d_mesh_partial, d_grid_vec_partial)

    ! Calculate exact solution on the grid
    ALLOCATE( d_grid_ex_vec_partial( grid%n_loc))
    DO n = 1, grid%n_loc
      n_glob = grid%n1 + n - 1
      i = grid%n2ij( n_glob,1)
      j = grid%n2ij( n_glob,2)
      x = grid%x( i)
      y = grid%y( j)
      CALL test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_grid_ex_vec_partial( n) = d
    END DO

    ! Calculate error
    maxerr = 0._dp
    DO n = 1, grid%n_loc
      n_glob = grid%n1 + n - 1
      i = grid%n2ij( n_glob,1)
      j = grid%n2ij( n_glob,2)
      ! Skip the border
      IF (i == 1 .OR. i == grid%nx .OR. j == 1 .OR. j == grid%ny) CYCLE
      maxerr = MAX( maxerr, ABS( d_grid_vec_partial( n) - d_grid_ex_vec_partial( n)))
    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr > 0.6E-1_dp) found_errors = .TRUE.

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated mesh to grid remapping')
    ELSE
      IF (par%master) CALL warning('found errors in mesh to grid remapping')
    END IF

    ! Write results to a NetCDF file
    IF (do_write_results_to_netcdf) THEN

      ! Create a file and write the mesh to it
      filename = TRIM( routine_name) // '_output.nc'
      CALL create_new_netcdf_file_for_writing( filename, ncid)
      CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)

      ! Add all the variables
      CALL add_field_grid_dp_2D_notime( filename, ncid, 'd_grid')
      CALL add_field_grid_dp_2D_notime( filename, ncid, 'd_grid_ex')

      ! Write all the variables
      CALL write_to_field_multiple_options_grid_dp_2D_notime( grid, filename, ncid, 'd_grid'   , d_grid_vec_partial   )
      CALL write_to_field_multiple_options_grid_dp_2D_notime( grid, filename, ncid, 'd_grid_ex', d_grid_ex_vec_partial)

      ! Close the file
      CALL close_netcdf_file( ncid)

    END IF ! IF (do_write_results_to_netcdf) THEN

    ! Clean up after yourself
    CALL deallocate_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_remapping_mesh2grid

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

END MODULE unit_tests_mesh
