MODULE unit_tests_mesh

  ! Unit tests for different mesh creation / discretisation stuff.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: allocate_mesh_primary, deallocate_mesh
  USE mesh_utilities                                         , ONLY: check_mesh
  USE mesh_creation                                          , ONLY: initialise_dummy_mesh
  USE mesh_refinement                                        , ONLY: refine_mesh_uniform, Lloyds_algorithm_single_iteration, mesh_add_smileyface, &
                                                                     mesh_add_UFEMISM_letters
  USE mesh_parallel_creation                                 , ONLY: merge_submeshes, broadcast_mesh
  USE mesh_secondary                                         , ONLY: calc_all_secondary_mesh_data
  USE mesh_operators                                         , ONLY: calc_all_matrix_operators_mesh, calc_3D_matrix_operators_mesh, calc_3D_gradient_bk_bk
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file
  USE netcdf_output                                          , ONLY: setup_mesh_in_netcdf_file, add_field_mesh_dp_2D_notime, add_field_mesh_dp_2D_b_notime, &
                                                                     add_field_mesh_dp_2D_c_notime, write_to_field_multopt_mesh_dp_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_b_notime, write_to_field_multopt_mesh_dp_2D_c_notime, &
                                                                     setup_xy_grid_in_netcdf_file, add_field_grid_dp_2D_notime, &
                                                                     write_to_field_multopt_grid_dp_2D_notime, add_zeta_dimension_to_file, &
                                                                     add_field_mesh_dp_3D_b_notime, write_to_field_multopt_mesh_dp_3D_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_b_notime
  USE petsc_basic                                            , ONLY: multiply_CSR_matrix_with_vector_1D
  USE grid_basic                                             , ONLY: type_grid, setup_square_grid, distribute_gridded_data_from_master_dp_2D, &
                                                                     calc_grid_mask_as_polygons
  USE grid_lonlat_basic                                      , ONLY: type_grid_lonlat, setup_simple_lonlat_grid, distribute_lonlat_gridded_data_from_master_dp_2D
  USE mesh_remapping                                         , ONLY: map_from_xy_grid_to_mesh_2D, map_from_mesh_to_xy_grid_2D, map_from_lonlat_grid_to_mesh_2D, &
                                                                     map_from_mesh_to_mesh_2D, clear_all_maps_involving_this_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ice_model_memory                                       , ONLY: allocate_ice_model
  USE ice_model_utilities                                    , ONLY: calc_zeta_gradients

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_mesh_unit_tests
    ! Run all mesh unit tests

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_all_mesh_unit_tests'
    TYPE(type_mesh)                                    :: mesh, mesh2

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all mesh unit tests
    CALL test_mesh_creation_basic_single_core
    CALL test_mesh_creation_basic_two_cores(       mesh)
    CALL test_mesh_creation_basic_two_cores_prime( mesh2)
    CALL test_mesh_operators_basic(                mesh)
    CALL test_mesh_operators_3D(                   mesh)
    CALL test_remapping_grid2mesh(                 mesh)
    CALL test_remapping_mesh2grid(                 mesh)
    CALL test_remapping_lonlat2mesh(               mesh)
    CALL test_remapping_mesh2mesh(                 mesh, mesh2)

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
    ! Expected values for the mesh (last updated: 2023-07-05)
    INTEGER,  PARAMETER                                :: nV_expected     = 4435
    INTEGER,  PARAMETER                                :: nTri_expected   = 8765
    INTEGER,  PARAMETER                                :: nE_expected     = 13199
    REAL(dp), PARAMETER                                :: Amin_expected   = 0.10398E+10_dp
    REAL(dp), PARAMETER                                :: Amax_expected   = 0.84217E+11_dp
    REAL(dp), PARAMETER                                :: RA_max_expected = 0.48151E+01_dp

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

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( mesh)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    IF (mesh%nV < NINT( REAL( nV_expected,dp) * 0.9_dp) .OR. mesh%nV > NINT( REAL( nV_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of vertices! Expected {int_01}, found {int_02}', int_01 = nV_expected, int_02 = mesh%nV)
    END IF

    IF (mesh%nTri < NINT( REAL( nTri_expected,dp) * 0.9_dp) .OR. mesh%nTri > NINT( REAL( nTri_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of triangles! Expected {int_01}, found {int_02}', int_01 = nTri_expected, int_02 = mesh%nTri)
    END IF

    IF (mesh%nE < NINT( REAL( nE_expected,dp) * 0.9_dp) .OR. mesh%nE > NINT( REAL( nE_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of edges! Expected {int_01}, found {int_02}', int_01 = nE_expected, int_02 = mesh%nE)
    END IF

    IF (MINVAL( mesh%A) < Amin_expected * 0.5_dp) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amin_expected, dp_02 = MINVAL( mesh%A))
    END IF

    IF (MAXVAL( mesh%A) > Amax_expected * 2.0_dp) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly small vertices! Expected MAXVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amax_expected, dp_02 = MAXVAL( mesh%A))
    END IF

    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      END DO
    END DO
    IF (RA_max > RA_max_expected * 2.0_dp)  THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly high resolution gradient! Expected RA_max = {dp_01}, found {dp_02}', dp_01 = RA_max_expected, dp_02 = RA_max)
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
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    CALL deallocate_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_creation_basic_single_core

  SUBROUTINE test_mesh_creation_basic_two_cores( mesh)
    ! Test creation of a very simple mesh with parallelised mesh generation.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(OUT)   :: mesh

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
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width
    INTEGER                                            :: vi, ci, vj
    REAL(dp)                                           :: RA, RA_max
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    ! Expected values for the mesh (last updated: 2023-07-05)
    INTEGER,  PARAMETER                                :: nV_expected     = 5873
    INTEGER,  PARAMETER                                :: nTri_expected   = 11619
    INTEGER,  PARAMETER                                :: nE_expected     = 17491
    REAL(dp), PARAMETER                                :: Amin_expected   = 0.11096E+10_dp
    REAL(dp), PARAMETER                                :: Amax_expected   = 0.74507E+11_dp
    REAL(dp), PARAMETER                                :: RA_max_expected = 0.34108E+01_dp

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
    CALL broadcast_mesh( mesh)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    tstart = MPI_WTIME()
    CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    tcomp = tcomp + MPI_WTIME() - tstart

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( mesh)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    IF (mesh%nV < NINT( REAL( nV_expected,dp) * 0.9_dp) .OR. mesh%nV > NINT( REAL( nV_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of vertices! Expected {int_01}, found {int_02}', int_01 = nV_expected, int_02 = mesh%nV)
    END IF

    IF (mesh%nTri < NINT( REAL( nTri_expected,dp) * 0.9_dp) .OR. mesh%nTri > NINT( REAL( nTri_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of triangles! Expected {int_01}, found {int_02}', int_01 = nTri_expected, int_02 = mesh%nTri)
    END IF

    IF (mesh%nE < NINT( REAL( nE_expected,dp) * 0.9_dp) .OR. mesh%nE > NINT( REAL( nE_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of edges! Expected {int_01}, found {int_02}', int_01 = nE_expected, int_02 = mesh%nE)
    END IF

    IF (MINVAL( mesh%A) < Amin_expected * 0.5_dp) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amin_expected, dp_02 = MINVAL( mesh%A))
    END IF

    IF (MAXVAL( mesh%A) > Amax_expected * 2.0_dp) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly small vertices! Expected MAXVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amax_expected, dp_02 = MAXVAL( mesh%A))
    END IF

    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      END DO
    END DO
    IF (RA_max > RA_max_expected * 2.0_dp)  THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly high resolution gradient! Expected RA_max = {dp_01}, found {dp_02}', dp_01 = RA_max_expected, dp_02 = RA_max)
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
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_creation_basic_two_cores

  SUBROUTINE test_mesh_creation_basic_two_cores_prime( mesh)
    ! Test creation of another very simple mesh with parallelised mesh generation.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(OUT)   :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_mesh_creation_basic_two_cores_prime'
    REAL(dp)                                           :: tstart, tcomp
    REAL(dp), PARAMETER                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    REAL(dp), PARAMETER                                :: xmax =  3040E3_dp
    REAL(dp), PARAMETER                                :: ymin = -3040E3_dp
    REAL(dp), PARAMETER                                :: ymax =  3040E3_dp
    REAL(dp), PARAMETER                                :: lambda_M    = 0._dp
    REAL(dp), PARAMETER                                :: phi_M       = -90._dp
    REAL(dp), PARAMETER                                :: beta_stereo = 71._dp
    CHARACTER(LEN=256), PARAMETER                      :: name = 'test_mesh2'
    REAL(dp)                                           :: alpha_min, res_max
    REAL(dp)                                           :: res, width
    INTEGER                                            :: vi, ci, vj
    REAL(dp)                                           :: RA, RA_max
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    ! Expected values for the mesh (last updated: 2023-07-05)
    INTEGER,  PARAMETER                                :: nV_expected     = 3299
    INTEGER,  PARAMETER                                :: nTri_expected   = 6460
    INTEGER,  PARAMETER                                :: nE_expected     = 9758
    REAL(dp), PARAMETER                                :: Amin_expected   = 0.30736E+10_dp
    REAL(dp), PARAMETER                                :: Amax_expected   = 0.53848E+11_dp
    REAL(dp), PARAMETER                                :: RA_max_expected = 0.29651E+01_dp

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

    ! Refine the mesh with a uniform 350 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 350E3_dp
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

    ! Smooth the mesh again
    tstart = MPI_WTIME()
    CALL Lloyds_algorithm_single_iteration( mesh, alpha_min)
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
    CALL broadcast_mesh( mesh)
    tcomp = tcomp + MPI_WTIME() - tstart
    CALL check_mesh( mesh)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    tstart = MPI_WTIME()
    CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    tcomp = tcomp + MPI_WTIME() - tstart

    ! Calculate all matrix operators
    CALL calc_all_matrix_operators_mesh( mesh)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    IF (mesh%nV < NINT( REAL( nV_expected,dp) * 0.9_dp) .OR. mesh%nV > NINT( REAL( nV_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of vertices! Expected {int_01}, found {int_02}', int_01 = nV_expected, int_02 = mesh%nV)
    END IF

    IF (mesh%nTri < NINT( REAL( nTri_expected,dp) * 0.9_dp) .OR. mesh%nTri > NINT( REAL( nTri_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of triangles! Expected {int_01}, found {int_02}', int_01 = nTri_expected, int_02 = mesh%nTri)
    END IF

    IF (mesh%nE < NINT( REAL( nE_expected,dp) * 0.9_dp) .OR. mesh%nE > NINT( REAL( nE_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexepcted amount of edges! Expected {int_01}, found {int_02}', int_01 = nE_expected, int_02 = mesh%nE)
    END IF

    IF (MINVAL( mesh%A) < Amin_expected * 0.5_dp) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amin_expected, dp_02 = MINVAL( mesh%A))
    END IF

    IF (MAXVAL( mesh%A) > Amax_expected * 2.0_dp) THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly small vertices! Expected MAXVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amax_expected, dp_02 = MAXVAL( mesh%A))
    END IF

    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      END DO
    END DO
    IF (RA_max > RA_max_expected * 2.0_dp)  THEN
      found_errors = .TRUE.
      IF (par%master) CALL warning('mesh has unexpectedly high resolution gradient! Expected RA_max = {dp_01}, found {dp_02}', dp_01 = RA_max_expected, dp_02 = RA_max)
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
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_creation_basic_two_cores_prime

  SUBROUTINE test_mesh_operators_basic( mesh)
    ! Test the basic mapping/gradient matrix operators

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_mesh_operators_basic'
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

  ! == Validate all matrix operators
  ! ================================

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

    maxerr_d_a_b      = 0._dp
    maxerr_d_c_b      = 0._dp
    maxerr_ddx_a_b    = 0._dp
    maxerr_ddx_b_b    = 0._dp
    maxerr_ddx_c_b    = 0._dp
    maxerr_ddy_a_b    = 0._dp
    maxerr_ddy_b_b    = 0._dp
    maxerr_ddy_c_b    = 0._dp
    maxerr_d2dx2_b_b  = 0._dp
    maxerr_d2dxdy_b_b = 0._dp
    maxerr_d2dy2_b_b  = 0._dp

    DO ti = mesh%ti1, mesh%ti2

      ! Skip border vertices, gradients are not needed there anyway
      IF (mesh%TriBI( ti) > 0) CYCLE

      ! Calculate errors
      maxerr_d_a_b      = MAX( maxerr_d_a_b     , ABS( d_a_b(      ti) - d_b_ex(      ti)))
      maxerr_d_c_b      = MAX( maxerr_d_c_b     , ABS( d_c_b(      ti) - d_b_ex(      ti)))
      maxerr_ddx_a_b    = MAX( maxerr_ddx_a_b   , ABS( ddx_a_b(    ti) - ddx_b_ex(    ti)))
      maxerr_ddx_b_b    = MAX( maxerr_ddx_b_b   , ABS( ddx_b_b(    ti) - ddx_b_ex(    ti)))
      maxerr_ddx_c_b    = MAX( maxerr_ddx_c_b   , ABS( ddx_c_b(    ti) - ddx_b_ex(    ti)))
      maxerr_ddy_a_b    = MAX( maxerr_ddy_a_b   , ABS( ddy_a_b(    ti) - ddy_b_ex(    ti)))
      maxerr_ddy_b_b    = MAX( maxerr_ddy_b_b   , ABS( ddy_b_b(    ti) - ddy_b_ex(    ti)))
      maxerr_ddy_c_b    = MAX( maxerr_ddy_c_b   , ABS( ddy_c_b(    ti) - ddy_b_ex(    ti)))
      maxerr_d2dx2_b_b  = MAX( maxerr_d2dx2_b_b , ABS( d2dx2_b_b(  ti) - d2dx2_b_ex(  ti)))
      maxerr_d2dxdy_b_b = MAX( maxerr_d2dxdy_b_b, ABS( d2dxdy_b_b( ti) - d2dxdy_b_ex( ti)))
      maxerr_d2dy2_b_b  = MAX( maxerr_d2dy2_b_b , ABS( d2dy2_b_b(  ti) - d2dy2_b_ex(  ti)))

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Find maximum errors across all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_a_b     , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_c_b     , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_a_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_b_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_c_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_a_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_b_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_c_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d2dx2_b_b , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d2dxdy_b_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d2dy2_b_b , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

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

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
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
    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_a_ex'     , d_a_ex     )
    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_ex'   , ddx_a_ex   )
    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_ex'   , ddy_a_ex   )

    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_b_ex'     , d_b_ex     )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_ex'   , ddx_b_ex   )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_ex'   , ddy_b_ex   )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_ex' , d2dx2_b_ex )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_ex', d2dxdy_b_ex)
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_ex' , d2dy2_b_ex )

    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_c_ex'     , d_c_ex     )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_c_ex'   , ddx_c_ex   )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_c_ex'   , ddy_c_ex   )

    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_a_b'      , d_a_b      )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_a_c'      , d_a_c      )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_b_a'      , d_b_a      )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_b_c'      , d_b_c      )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_c_a'      , d_c_a      )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_c_b'      , d_c_b      )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_a'    , ddx_a_a    )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_a_b'    , ddx_a_b    )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_a_c'    , ddx_a_c    )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_b_a'    , ddx_b_a    )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_b'    , ddx_b_b    )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_b_c'    , ddx_b_c    )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_c_a'    , ddx_c_a    )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_c_b'    , ddx_c_b    )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_c_c'    , ddx_c_c    )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_a'    , ddy_a_a    )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_a_b'    , ddy_a_b    )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_a_c'    , ddy_a_c    )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_b_a'    , ddy_b_a    )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_b'    , ddy_b_b    )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_b_c'    , ddy_b_c    )

    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_c_a'    , ddy_c_a    )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_c_b'    , ddy_c_b    )
    CALL write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_c_c'    , ddy_c_c    )

    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_b'  , d2dx2_b_b  )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_b' , d2dxdy_b_b )
    CALL write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_b'  , d2dy2_b_b  )

    ! Close the file
    CALL close_netcdf_file( ncid)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_operators_basic

  SUBROUTINE test_mesh_operators_3D( mesh)
    ! Test the 3-D gradient matrix operators

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_mesh_operators_3D'
    TYPE(type_ice_model)                               :: ice
    INTEGER                                            :: vi,ti,k
    REAL(dp)                                           :: x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: f_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dx_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dy_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dz_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dx2_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dxdy_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dy2_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dz2_bk_ex
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dx_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dy_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dz_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dx2_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dxdy_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dy2_bk
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dz2_bk
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate ice model memory
    CALL allocate_ice_model( mesh, ice)

    ! Set up a nice test geometry

    xmin = mesh%xmin
    xmax = mesh%xmax
    ymin = mesh%ymin
    ymax = mesh%ymax

    DO vi = mesh%vi1, mesh%vi2
    DO k  = 1, mesh%nz

      x    = mesh%V( vi,1)
      y    = mesh%V( vi,2)
      zeta = mesh%zeta( k)

      CALL test_function_3D( x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2)

      ice%Hi( vi) = Hi
      ice%Hb( vi) = Hb
      ice%Hs( vi) = Hs
      ice%SL( vi) = SL

    END DO
    END DO

    ! Calculate the exact solution

    ALLOCATE( f_bk_ex(        mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( df_dx_bk_ex(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( df_dy_bk_ex(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( df_dz_bk_ex(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dx2_bk_ex(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dxdy_bk_ex( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dy2_bk_ex(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dz2_bk_ex(  mesh%ti1:mesh%ti2,mesh%nz))

    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, mesh%nz

      x    = mesh%TriGC( ti,1)
      y    = mesh%TriGC( ti,2)
      zeta = mesh%zeta( k)

      CALL test_function_3D( x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2)

      f_bk_ex(        ti,k) = f
      df_dx_bk_ex(    ti,k) = dfdx
      df_dy_bk_ex(    ti,k) = dfdy
      df_dz_bk_ex(    ti,k) = dfdz
      d2f_dx2_bk_ex(  ti,k) = d2fdx2
      d2f_dxdy_bk_ex( ti,k) = d2fdxdy
      d2f_dy2_bk_ex(  ti,k) = d2fdy2
      d2f_dz2_bk_ex(  ti,k) = d2fdz2

    END DO
    END DO

    ! Calculate gradients of zeta
    CALL calc_zeta_gradients( mesh, ice)

    ! Calculate 3-D gradient operators
    CALL calc_3D_matrix_operators_mesh( mesh, ice)

    ! Calculate gradients of f

    ALLOCATE( df_dx_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( df_dy_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( df_dz_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dx2_bk(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dxdy_bk( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dy2_bk(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dz2_bk(  mesh%ti1:mesh%ti2,mesh%nz))

    CALL calc_3D_gradient_bk_bk( mesh, mesh%M2_ddx_bk_bk   , f_bk_ex, df_dx_bk   )
    CALL calc_3D_gradient_bk_bk( mesh, mesh%M2_ddy_bk_bk   , f_bk_ex, df_dy_bk   )
    CALL calc_3D_gradient_bk_bk( mesh, mesh%M2_ddz_bk_bk   , f_bk_ex, df_dz_bk   )
    CALL calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dx2_bk_bk , f_bk_ex, d2f_dx2_bk )
    CALL calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dxdy_bk_bk, f_bk_ex, d2f_dxdy_bk)
    CALL calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dy2_bk_bk , f_bk_ex, d2f_dy2_bk )
    CALL calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dz2_bk_bk , f_bk_ex, d2f_dz2_bk )

    ! If no errors occurred, we are happy
    found_errors = .FALSE.
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all 3-D mesh matrix operators')
    ELSE
      IF (par%master) CALL warning('found errors in 3-D mesh matrix operators')
    END IF

    ! Write results to a NetCDF file

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
    CALL add_zeta_dimension_to_file( filename, ncid, mesh%zeta)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'Hi')
    CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'Hb')
    CALL add_field_mesh_dp_2D_notime(   filename, ncid, 'Hs')

    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'f_bk_ex'       )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dx_bk_ex'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dy_bk_ex'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dz_bk_ex'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dx2_bk_ex' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dxdy_bk_ex')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dy2_bk_ex' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dz2_bk_ex' )

    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dx_bk'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dy_bk'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dz_bk'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dx2_bk' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dxdy_bk')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dy2_bk' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dz2_bk' )

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'Hi'            , ice%Hi        )
    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'Hb'            , ice%Hb        )
    CALL write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'Hs'            , ice%Hs        )

    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'f_bk_ex'       , f_bk_ex       )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dx_bk_ex'   , df_dx_bk_ex   )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dy_bk_ex'   , df_dy_bk_ex   )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dz_bk_ex'   , df_dz_bk_ex   )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dx2_bk_ex' , d2f_dx2_bk_ex )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dxdy_bk_ex', d2f_dxdy_bk_ex)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dy2_bk_ex' , d2f_dy2_bk_ex )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dz2_bk_ex' , d2f_dz2_bk_ex )

    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dx_bk'      , df_dx_bk      )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dy_bk'      , df_dy_bk      )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dz_bk'      , df_dz_bk      )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dx2_bk'    , d2f_dx2_bk    )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dxdy_bk'   , d2f_dxdy_bk   )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dy2_bk'    , d2f_dy2_bk    )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dz2_bk'    , d2f_dz2_bk    )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    DEALLOCATE( f_bk_ex)
    DEALLOCATE( df_dx_bk_ex)
    DEALLOCATE( df_dy_bk_ex)
    DEALLOCATE( df_dz_bk_ex)
    DEALLOCATE( d2f_dx2_bk_ex)
    DEALLOCATE( d2f_dxdy_bk_ex)
    DEALLOCATE( d2f_dy2_bk_ex)
    DEALLOCATE( d2f_dz2_bk_ex)
    DEALLOCATE( df_dx_bk)
    DEALLOCATE( df_dy_bk)
    DEALLOCATE( df_dz_bk)
    DEALLOCATE( d2f_dx2_bk)
    DEALLOCATE( d2f_dxdy_bk)
    DEALLOCATE( d2f_dy2_bk)
    DEALLOCATE( d2f_dz2_bk)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_mesh_operators_3D

  SUBROUTINE test_remapping_grid2mesh( mesh)
    ! Test remapping from a square grid to a mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_remapping_grid2mesh'
    CHARACTER(LEN=256)                                 :: name
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

  ! == Set up a square grid
  ! =======================

    name = 'test_grid'
    dx   = 32E3_dp
    CALL setup_square_grid( name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, dx, grid, lambda_M = mesh%lambda_M, phi_M = mesh%phi_M, beta_stereo = mesh%beta_stereo)

  ! == Calculate, apply, and validate grid-to-mesh remapping operator
  ! =================================================================

    found_errors = .FALSE.

    ! Create a nice gridded data field on the master
    IF (par%master) THEN
      ALLOCATE( d_grid( grid%nx, grid%ny), source = 0._dp)
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        x = grid%x( i)
        y = grid%y( j)
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

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh_ex')

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh'   , d_mesh_partial   )
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    CALL clear_all_maps_involving_this_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_remapping_grid2mesh

  SUBROUTINE test_remapping_mesh2grid( mesh)
    ! Test remapping from a mesh to a square grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_remapping_mesh2grid'
    CHARACTER(LEN=256)                                 :: name
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

  ! == Set up a square grid
  ! =======================

    name = 'test_grid'
    dx   = 32E3_dp
    CALL setup_square_grid( name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, dx, grid, lambda_M = mesh%lambda_M, phi_M = mesh%phi_M, beta_stereo = mesh%beta_stereo)

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

    ! Create a file and write the mesh to it
    filename = TRIM(C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    ! Add all the variables
    CALL add_field_grid_dp_2D_notime( filename, ncid, 'd_grid')
    CALL add_field_grid_dp_2D_notime( filename, ncid, 'd_grid_ex')

    ! Write all the variables
    CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid'   , d_grid_vec_partial   )
    CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid_ex', d_grid_ex_vec_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    CALL clear_all_maps_involving_this_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_remapping_mesh2grid

  SUBROUTINE test_remapping_lonlat2mesh( mesh)
    ! Test remapping from a lon/lat-grid to a mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_remapping_lonlat2mesh'
    CHARACTER(LEN=256)                                 :: name
    INTEGER                                            :: nlon, nlat
    TYPE(type_grid_lonlat)                             :: grid
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial
    REAL(dp)                                           :: lon,lat,d
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_ex_partial
    INTEGER                                            :: vi, vi_glob
    REAL(dp)                                           :: maxerr
    LOGICAL                                            :: found_errors
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Set up a square grid
  ! =======================

    ! Set up a lon/lat-grid
    name = 'test_lonlat_grid'
    nlon = 360
    nlat = 180
    CALL setup_simple_lonlat_grid( name, nlon, nlat, grid)

  ! == Calculate, apply, and validate grid-to-mesh remapping operator
  ! =================================================================

    found_errors = .FALSE.

    ! Create a nice gridded data field on the master
    IF (par%master) THEN
      ALLOCATE( d_grid( grid%nlon, grid%nlat))
      DO i = 1, grid%nlon
      DO j = 1, grid%nlat
        lon = grid%lon( i)
        lat = grid%lat( j)
        CALL test_function_lonlat( lon, lat, d)
        d_grid( i,j) = d
      END DO
      END DO
    END IF
    CALL sync

    ! Distribute gridded data
    ALLOCATE( d_grid_vec_partial( grid%n_loc))
    CALL distribute_lonlat_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
    IF (par%master) DEALLOCATE( d_grid)

    ! Map gridded data to the mesh
    ALLOCATE( d_mesh_partial( mesh%nV_loc))
    CALL map_from_lonlat_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Calculate exact solution
    ALLOCATE( d_mesh_ex_partial( mesh%nV_loc))
    DO vi = 1, mesh%nV_loc
      vi_glob = mesh%vi1 + vi - 1
      lon = mesh%lon( vi_glob)
      lat = mesh%lat( vi_glob)
      CALL test_function_lonlat( lon, lat, d)
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

    IF (maxerr > 0.5E-4_dp) found_errors = .TRUE.

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated lon/lat to mesh remapping')
    ELSE
      IF (par%master) CALL warning('found errors in lon/lat to mesh remapping')
    END IF

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh_ex')

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh'   , d_mesh_partial   )
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex_partial)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    CALL clear_all_maps_involving_this_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_remapping_lonlat2mesh

  SUBROUTINE test_remapping_mesh2mesh( mesh1, mesh2)
    ! Test remapping from a mesh to another mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh1, mesh2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_remapping_mesh2mesh'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d1_ex, d2_ex
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d12_nn, d12_trilin, d12_cons
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d21_nn, d21_trilin, d21_cons
    INTEGER                                            :: vi
    REAL(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    CHARACTER(LEN=256)                                 :: method
    LOGICAL                                            :: found_errors
    REAL(dp)                                           :: maxerr
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( d1_ex(      mesh1%vi1:mesh1%vi2), source = 0._dp)
    ALLOCATE( d2_ex(      mesh2%vi1:mesh2%vi2), source = 0._dp)

    ALLOCATE( d12_nn(     mesh2%vi1:mesh2%vi2), source = 0._dp)
    ALLOCATE( d12_trilin( mesh2%vi1:mesh2%vi2), source = 0._dp)
    ALLOCATE( d12_cons(   mesh2%vi1:mesh2%vi2), source = 0._dp)

    ALLOCATE( d21_nn(     mesh1%vi1:mesh1%vi2), source = 0._dp)
    ALLOCATE( d21_trilin( mesh1%vi1:mesh1%vi2), source = 0._dp)
    ALLOCATE( d21_cons(   mesh1%vi1:mesh1%vi2), source = 0._dp)

    ! Initialise exact solutions
    DO vi = mesh1%vi1, mesh1%vi2
      x = mesh1%V( vi,1)
      y = mesh1%V( vi,2)
      CALL test_function( x, y, mesh1%xmin, mesh1%xmax, mesh1%ymin, mesh1%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d1_ex( vi) = d
    END DO
    DO vi = mesh2%vi1, mesh2%vi2
      x = mesh2%V( vi,1)
      y = mesh2%V( vi,2)
      CALL test_function( x, y, mesh2%xmin, mesh2%xmax, mesh2%ymin, mesh2%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d2_ex( vi) = d
    END DO

  ! == Remap data
  ! =============

    method = 'nearest_neighbour'
    CALL map_from_mesh_to_mesh_2D( mesh1, mesh2, d1_ex, d12_nn, method)
    CALL map_from_mesh_to_mesh_2D( mesh2, mesh1, d2_ex, d21_nn, method)

    method = 'trilin'
    CALL map_from_mesh_to_mesh_2D( mesh1, mesh2, d1_ex, d12_trilin, method)
    CALL map_from_mesh_to_mesh_2D( mesh2, mesh1, d2_ex, d21_trilin, method)

    method = '2nd_order_conservative'
    CALL map_from_mesh_to_mesh_2D( mesh1, mesh2, d1_ex, d12_cons, method)
    CALL map_from_mesh_to_mesh_2D( mesh2, mesh1, d2_ex, d21_cons, method)

  ! == Validation
  ! =============

    found_errors = .FALSE.

    ! Nearest-neighbour
    ! 2023-06-05: maxerr = 0.11995E+00
    maxerr = 0._dp
    DO vi = mesh2%vi1, mesh2%vi2
      IF (mesh2%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d12_nn( vi) - d2_ex( vi)))
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.11995E+00 * 2.0) THEN
      CALL warning('unexpectedly high local errors detected in mesh-to-mesh nearest-neighbour remapping; expected 0.11995E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    END IF

    ! 2023-06-05: maxerr = 0.16505E+00
    maxerr = 0._dp
    DO vi = mesh1%vi1, mesh1%vi2
      IF (mesh1%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d21_nn( vi) - d1_ex( vi)))
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.16505E+00 * 2.0) THEN
      CALL warning('unexpectedly high local errors detected in mesh-to-mesh nearest-neighbour remapping: expected 0.16505E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    END IF

    ! Trilinear
    ! 2023-06-05: maxerr = 0.16466E-01
    maxerr = 0._dp
    DO vi = mesh2%vi1, mesh2%vi2
      IF (mesh2%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d12_trilin( vi) - d2_ex( vi)))
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.16466E-01 * 2.0) THEN
      CALL warning('unexpectedly high local errors detected in mesh-to-mesh trilinear remapping: expected 0.16466E-01, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    END IF

    ! 2023-06-05: maxerr = 0.18900E-01
    maxerr = 0._dp
    DO vi = mesh1%vi1, mesh1%vi2
      IF (mesh1%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d21_trilin( vi) - d1_ex( vi)))
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.18900E-01 * 2.0) THEN
      CALL warning('unexpectedly high local errors detected in mesh-to-mesh trilinear remapping: expected 0.18900E-01, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    END IF

    ! 2nd-order conservative
    ! 2023-06-05: maxerr = 0.27873E+00
    maxerr = 0._dp
    DO vi = mesh2%vi1, mesh2%vi2
      IF (mesh2%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d12_cons( vi) - d2_ex( vi)))
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.27873E+00 * 2.0) THEN
      CALL warning('unexpectedly high local errors detected in mesh-to-mesh 2nd-order conservative remapping: expected 0.27873E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    END IF

    ! 2023-06-05: maxerr = 0.56043E+00
    maxerr = 0._dp
    DO vi = mesh1%vi1, mesh1%vi2
      IF (mesh1%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d21_cons( vi) - d1_ex( vi)))
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.56043E+00 * 2.0) THEN
      CALL warning('unexpectedly high local errors detected in mesh-to-mesh 2nd-order conservative remapping: expected 0.56043E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    END IF

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated mesh to mesh remapping')
    ELSE
      IF (par%master) CALL warning('found errors in mesh to mesh remapping')
    END IF

    ! Write results to a NetCDF file

    ! Mesh 1
    ! ======

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output1.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh1)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd1_ex')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd21_nn')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd21_trilin')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd21_cons')

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd1_ex'     , d1_ex     )
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd21_nn'    , d21_nn    )
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd21_trilin', d21_trilin)
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd21_cons'  , d21_cons  )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Mesh 2
    ! ======

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output2.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh2)

    ! Add all the variables
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd2_ex')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd12_nn')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd12_trilin')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'd12_cons')

    ! Write all the variables
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd2_ex'     , d2_ex     )
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd12_nn'    , d12_nn    )
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd12_trilin', d12_trilin)
    CALL write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd12_cons'  , d12_cons  )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    DEALLOCATE( d1_ex)
    DEALLOCATE( d2_ex)
    DEALLOCATE( d12_nn)
    DEALLOCATE( d12_trilin)
    DEALLOCATE( d12_cons)
    DEALLOCATE( d21_nn)
    DEALLOCATE( d21_trilin)
    DEALLOCATE( d21_cons)

    ! Clean up after yourself
    CALL clear_all_maps_involving_this_mesh( mesh1)
    CALL clear_all_maps_involving_this_mesh( mesh2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_remapping_mesh2mesh

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

  SUBROUTINE test_function_3D( x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2)
    ! A simple test function to validate the 3-D matrix operators

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                   INTENT(IN)        :: x,y,xmin,xmax,ymin,ymax,zeta
    REAL(dp),                   INTENT(OUT)       :: Hi, Hb, Hs, SL
    REAL(dp),                   INTENT(OUT)       :: f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2

    ! Local variables:
    REAL(dp), PARAMETER                           :: a  = 500._dp  ! Amplitude of bedrock undulations
    REAL(dp), PARAMETER                           :: h0 = 2000._dp ! Uniform surface elevation
    REAL(dp)                                      :: cx,cy,z,cz

    ! Ice-sheet geometry: sort of like ISMIP-HOM but without
    ! the uniform bed slope. So a layer of ice with a flat surface at h = h0,
    ! lying on a bed with undulations of amplitude a.

    cx = 2._dp * pi / (xmax - xmin)
    cy = 3._dp * pi / (ymax - ymin)

    Hb = a * SIN( cx * x) * SIN( cy * y)
    Hs = h0
    Hi = Hs - Hb
    SL = -10000._dp

    ! Actual vertical coordinate
    z  = Hs - zeta * Hi
    cz = 4._dp * pi / (h0 + a)

    ! The function f and its gradients
    f       = (cx * x)**3 + (cy * y)**3 + (cz * z)**3 + (cx * x * cy * y)**2
    dfdx    = 3._dp * cx**3 * x**2 + 2._dp * cx**2 * cy**2 * x * y**2
    dfdy    = 3._dp * cy**3 * y**2 + 2._dp * cx**2 * cy**2 * x**2 * y
    dfdz    = 3._dp * cz**3 * z**2
    d2fdx2  = 6._dp * cx**3 * x + 2._dp * cx**2 * cy**2 * y**2
    d2fdxdy = 4._dp * cx**2 * cy**2 * x * y
    d2fdy2  = 6._dp * cy**3 * y + 2._dp * cx**2 * cy**2 * x**2
    d2fdz2  = 6._dp * cz**3 * z

  END SUBROUTINE test_function_3D

END MODULE unit_tests_mesh
