module unit_tests_mesh_old

  ! Unit tests for different mesh creation / discretisation stuff.

! ===== Preamble =====
! ====================

  use mpi
  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging                  , only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mesh_types                                             , only: type_mesh
  use mesh_memory                                            , only: allocate_mesh_primary, deallocate_mesh
  use mesh_utilities                                         , only: check_mesh
  use mesh_creation                                          , only: initialise_dummy_mesh
  use mesh_refinement                                        , only: refine_mesh_uniform, Lloyds_algorithm_single_iteration, mesh_add_smileyface, &
                                                                     mesh_add_UFEMISM_letters
  use mesh_parallel_creation                                 , only: merge_submeshes, broadcast_mesh
  use mesh_secondary                                         , only: calc_all_secondary_mesh_data
  use mesh_operators                                         , only: calc_all_matrix_operators_mesh, calc_3D_matrix_operators_mesh, calc_3D_gradient_bk_bk
  use netcdf_basic                                           , only: create_new_netcdf_file_for_writing, close_netcdf_file
  use netcdf_output                                          , only: setup_mesh_in_netcdf_file, add_field_mesh_dp_2D_notime, add_field_mesh_dp_2D_b_notime, &
                                                                     add_field_mesh_dp_2D_c_notime, write_to_field_multopt_mesh_dp_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_b_notime, write_to_field_multopt_mesh_dp_2D_c_notime, &
                                                                     setup_xy_grid_in_netcdf_file, add_field_grid_dp_2D_notime, &
                                                                     write_to_field_multopt_grid_dp_2D_notime, add_zeta_dimension_to_file, &
                                                                     add_field_mesh_dp_3D_b_notime, write_to_field_multopt_mesh_dp_3D_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_b_notime
  use petsc_basic                                            , only: multiply_CSR_matrix_with_vector_1D
  use grid_basic                                             , only: type_grid, setup_square_grid, distribute_gridded_data_from_master_dp_2D, &
                                                                     calc_grid_mask_as_polygons
  use grid_lonlat_basic                                      , only: type_grid_lonlat, setup_simple_lonlat_grid, distribute_lonlat_gridded_data_from_master_dp_2D
  use mesh_remapping                                         , only: map_from_xy_grid_to_mesh_2D, map_from_mesh_to_xy_grid_2D, map_from_lonlat_grid_to_mesh_2D, &
                                                                     map_from_mesh_to_mesh_2D, clear_all_maps_involving_this_mesh
  use ice_model_types                                        , only: type_ice_model
  use ice_model_memory                                       , only: allocate_ice_model
  use ice_model_utilities                                    , only: calc_zeta_gradients

  implicit none

contains

! ===== Subroutines =====
! =======================

  subroutine unit_tests_mesh_creation_main( test_name_parent)
    ! Run all mesh unit tests

    implicit none

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_mesh_creation_main'
    character(len=1024), parameter :: test_name_local = 'mesh_creation'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh, mesh2

    ! Add routine to path
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Run all mesh creation unit tests
    call test_mesh_creation_basic_single_core
    call test_mesh_creation_basic_two_cores ( mesh)
    call test_mesh_creation_basic_two_cores_prime( mesh2)
    !    call test_mesh_operators_basic(                mesh)
    !    call test_mesh_operators_3D(                   mesh)
    !    call test_remapping_grid2mesh(                 mesh)
    !    call test_remapping_mesh2grid(                 mesh)
    !    call test_remapping_lonlat2mesh(               mesh)
    !    call test_remapping_mesh2mesh(                 mesh, mesh2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine unit_tests_mesh_creation_main

  subroutine test_mesh_creation_basic_single_core
    ! Test creation of a very simple mesh without parallelised mesh generation.

    implicit none

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_mesh_creation_basic_single_core'
    real(dp)                                           :: tstart, tcomp
    real(dp), parameter                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    real(dp), parameter                                :: xmax =  3040E3_dp
    real(dp), parameter                                :: ymin = -3040E3_dp
    real(dp), parameter                                :: ymax =  3040E3_dp
    real(dp), parameter                                :: lambda_M    = 0._dp
    real(dp), parameter                                :: phi_M       = -90._dp
    real(dp), parameter                                :: beta_stereo = 71._dp
    character(len=1024), parameter                     :: name = 'test_mesh'
    type(type_mesh)                                    :: mesh
    real(dp)                                           :: alpha_min, res_max
    real(dp)                                           :: res, width
    integer                                            :: vi, ci, vj
    real(dp)                                           :: RA, RA_max
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid
    ! Expected values for the mesh (last updated: 2023-07-05)
    integer,  parameter                                :: nV_expected     = 4435
    integer,  parameter                                :: nTri_expected   = 8765
    integer,  parameter                                :: nE_expected     = 13199
    real(dp), parameter                                :: Amin_expected   = 0.10398E+10_dp
    real(dp), parameter                                :: Amax_expected   = 0.84217E+11_dp
    real(dp), parameter                                :: RA_max_expected = 0.48151E+01_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! == Create a nice mesh, with a smileyface and the UFEMISM letters
    ! ================================================================

    tcomp = 0._dp

    ! Allocate memory
    call allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Initialise the dummy mesh
    call initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)
    call check_mesh( mesh)

    ! Refine the mesh with a uniform 400 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 400E3_dp
    tstart = MPI_WTIME()
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Smooth the mesh by applying a single iteration of Lloyd's algorithm
    tstart = MPI_WTIME()
    call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Add a smileyface
    res   = 80E3_dp
    width = 100E3_dp
    tstart = MPI_WTIME()
    call mesh_add_smileyface( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Add the UFEMISM letters
    res   = 50E3_dp
    width = 40E3_dp
    tstart = MPI_WTIME()
    call mesh_add_UFEMISM_letters( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Smooth the mesh again
    tstart = MPI_WTIME()
    call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    tstart = MPI_WTIME()
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    tcomp = tcomp + MPI_WTIME() - tstart

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! == Validation
    ! =============

    found_errors = .FALSE.

    IF (mesh%nV < NINT( real( nV_expected,dp) * 0.9_dp) .OR. mesh%nV > NINT( real( nV_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of vertices! Expected {int_01}, found {int_02}', int_01 = nV_expected, int_02 = mesh%nV)
    end IF

    IF (mesh%nTri < NINT( real( nTri_expected,dp) * 0.9_dp) .OR. mesh%nTri > NINT( real( nTri_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of triangles! Expected {int_01}, found {int_02}', int_01 = nTri_expected, int_02 = mesh%nTri)
    end IF

    IF (mesh%nE < NINT( real( nE_expected,dp) * 0.9_dp) .OR. mesh%nE > NINT( real( nE_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of edges! Expected {int_01}, found {int_02}', int_01 = nE_expected, int_02 = mesh%nE)
    end IF

    IF (MINVAL( mesh%A) < Amin_expected * 0.5_dp) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amin_expected, dp_02 = MINVAL( mesh%A))
    end IF

    IF (MAXVAL( mesh%A) > Amax_expected * 2.0_dp) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly small vertices! Expected MAXVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amax_expected, dp_02 = MAXVAL( mesh%A))
    end IF

    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      end DO
    end DO
    IF (RA_max > RA_max_expected * 2.0_dp)  THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly high resolution gradient! Expected RA_max = {dp_01}, found {dp_02}', dp_01 = RA_max_expected, dp_02 = RA_max)
    end IF

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('created a test mesh with {int_01} vertices, {int_02} triangles, and {int_03} edges in {dp_01} seconds', &
        int_01 = mesh%nV, int_02 = mesh%nTri, int_03 = mesh%nE, dp_01 = tcomp)
    ELSE
      IF (par%master) call warning('found errors in basic single-core mesh creation')
    end IF

    ! Write the resulting mesh to a NetCDF file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call deallocate_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine test_mesh_creation_basic_single_core

  subroutine test_mesh_creation_basic_two_cores( mesh)
    ! Test creation of a very simple mesh with parallelised mesh generation.

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(OUT)   :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_mesh_creation_basic_two_cores'
    real(dp)                                           :: tstart, tcomp
    real(dp), parameter                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    real(dp), parameter                                :: xmax =  3040E3_dp
    real(dp), parameter                                :: ymin = -3040E3_dp
    real(dp), parameter                                :: ymax =  3040E3_dp
    real(dp), parameter                                :: lambda_M    = 0._dp
    real(dp), parameter                                :: phi_M       = -90._dp
    real(dp), parameter                                :: beta_stereo = 71._dp
    character(len=1024), parameter                     :: name = 'test_mesh'
    real(dp)                                           :: alpha_min, res_max
    real(dp)                                           :: res, width
    integer                                            :: vi, ci, vj
    real(dp)                                           :: RA, RA_max
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid
    ! Expected values for the mesh (last updated: 2023-07-05)
    integer,  parameter                                :: nV_expected     = 5873
    integer,  parameter                                :: nTri_expected   = 11619
    integer,  parameter                                :: nE_expected     = 17491
    real(dp), parameter                                :: Amin_expected   = 0.11096E+10_dp
    real(dp), parameter                                :: Amax_expected   = 0.74507E+11_dp
    real(dp), parameter                                :: RA_max_expected = 0.34108E+01_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! == Create a nice mesh on each process, with a smileyface and the UFEMISM letters
    ! ================================================================================

    tcomp = 0._dp

    ! Allocate memory
    call allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Safety - must be run on at least two cores!
    IF (par%n < 2) call crash('this validation check needs to be run on at least two cores!')

    ! Initialise the dummy mesh
    IF (par%i == 0) THEN
      call initialise_dummy_mesh( mesh, xmin, (xmin + xmax) / 2._dp, ymin, ymax)
    ELSE
      call initialise_dummy_mesh( mesh, (xmin + xmax) / 2._dp, xmax, ymin, ymax)
    end IF
    call check_mesh( mesh)

    ! Refine the mesh with a uniform 400 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 400E3_dp
    tstart = MPI_WTIME()
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Smooth the mesh by applying a single iteration of Lloyd's algorithm
    tstart = MPI_WTIME()
    call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Add a smileyface
    res   = 80E3_dp
    width = 100E3_dp
    tstart = MPI_WTIME()
    call mesh_add_smileyface( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Add the UFEMISM letters
    res   = 50E3_dp
    width = 40E3_dp
    tstart = MPI_WTIME()
    call mesh_add_UFEMISM_letters( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Smooth the mesh again
    tstart = MPI_WTIME()
    call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! == Merge the two meshes
    ! =======================

    ! Merge submeshes
    tstart = MPI_WTIME()
    call merge_submeshes( mesh, 0, 1, 'east-west')
    tcomp = tcomp + MPI_WTIME() - tstart
    IF (par%master) call check_mesh( mesh)

    ! Smooth again
    tstart = MPI_WTIME()
    IF (par%master) call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    IF (par%master) call check_mesh( mesh)

    ! Broadcast from Master
    tstart = MPI_WTIME()
    call broadcast_mesh( mesh)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    tstart = MPI_WTIME()
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    tcomp = tcomp + MPI_WTIME() - tstart

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! == Validation
    ! =============

    found_errors = .FALSE.

    IF (mesh%nV < NINT( real( nV_expected,dp) * 0.9_dp) .OR. mesh%nV > NINT( real( nV_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of vertices! Expected {int_01}, found {int_02}', int_01 = nV_expected, int_02 = mesh%nV)
    end IF

    IF (mesh%nTri < NINT( real( nTri_expected,dp) * 0.9_dp) .OR. mesh%nTri > NINT( real( nTri_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of triangles! Expected {int_01}, found {int_02}', int_01 = nTri_expected, int_02 = mesh%nTri)
    end IF

    IF (mesh%nE < NINT( real( nE_expected,dp) * 0.9_dp) .OR. mesh%nE > NINT( real( nE_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of edges! Expected {int_01}, found {int_02}', int_01 = nE_expected, int_02 = mesh%nE)
    end IF

    IF (MINVAL( mesh%A) < Amin_expected * 0.5_dp) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amin_expected, dp_02 = MINVAL( mesh%A))
    end IF

    IF (MAXVAL( mesh%A) > Amax_expected * 2.0_dp) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly small vertices! Expected MAXVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amax_expected, dp_02 = MAXVAL( mesh%A))
    end IF

    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      end DO
    end DO
    IF (RA_max > RA_max_expected * 2.0_dp)  THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly high resolution gradient! Expected RA_max = {dp_01}, found {dp_02}', dp_01 = RA_max_expected, dp_02 = RA_max)
    end IF

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('created a test mesh with {int_01} vertices, {int_02} triangles, and {int_03} edges in {dp_01} seconds', &
        int_01 = mesh%nV, int_02 = mesh%nTri, int_03 = mesh%nE, dp_01 = tcomp)
    ELSE
      IF (par%master) call warning('found errors in basic two-core mesh creation')
    end IF

    ! Write the resulting mesh to a NetCDF file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine test_mesh_creation_basic_two_cores

  subroutine test_mesh_creation_basic_two_cores_prime( mesh)
    ! Test creation of another very simple mesh with parallelised mesh generation.

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(OUT)   :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_mesh_creation_basic_two_cores_prime'
    real(dp)                                           :: tstart, tcomp
    real(dp), parameter                                :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    real(dp), parameter                                :: xmax =  3040E3_dp
    real(dp), parameter                                :: ymin = -3040E3_dp
    real(dp), parameter                                :: ymax =  3040E3_dp
    real(dp), parameter                                :: lambda_M    = 0._dp
    real(dp), parameter                                :: phi_M       = -90._dp
    real(dp), parameter                                :: beta_stereo = 71._dp
    character(len=1024), parameter                     :: name = 'test_mesh2'
    real(dp)                                           :: alpha_min, res_max
    real(dp)                                           :: res, width
    integer                                            :: vi, ci, vj
    real(dp)                                           :: RA, RA_max
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid
    ! Expected values for the mesh (last updated: 2023-07-05)
    integer,  parameter                                :: nV_expected     = 3299
    integer,  parameter                                :: nTri_expected   = 6460
    integer,  parameter                                :: nE_expected     = 9758
    real(dp), parameter                                :: Amin_expected   = 0.30736E+10_dp
    real(dp), parameter                                :: Amax_expected   = 0.53848E+11_dp
    real(dp), parameter                                :: RA_max_expected = 0.29651E+01_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! == Create a nice mesh on each process, with a smileyface and the UFEMISM letters
    ! ================================================================================

    tcomp = 0._dp

    ! Allocate memory
    call allocate_mesh_primary( mesh, name, 1000, 2000, 32)

    ! Safety - must be run on at least two cores!
    IF (par%n < 2) call crash('this validation check needs to be run on at least two cores!')

    ! Initialise the dummy mesh
    IF (par%i == 0) THEN
      call initialise_dummy_mesh( mesh, xmin, (xmin + xmax) / 2._dp, ymin, ymax)
    ELSE
      call initialise_dummy_mesh( mesh, (xmin + xmax) / 2._dp, xmax, ymin, ymax)
    end IF
    call check_mesh( mesh)

    ! Refine the mesh with a uniform 350 km resolution
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 350E3_dp
    tstart = MPI_WTIME()
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Smooth the mesh by applying a single iteration of Lloyd's algorithm
    tstart = MPI_WTIME()
    call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Add a smileyface
    res   = 80E3_dp
    width = 100E3_dp
    tstart = MPI_WTIME()
    call mesh_add_smileyface( mesh, res, width)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Smooth the mesh again
    tstart = MPI_WTIME()
    call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! == Merge the two meshes
    ! =======================

    ! Merge submeshes
    tstart = MPI_WTIME()
    call merge_submeshes( mesh, 0, 1, 'east-west')
    tcomp = tcomp + MPI_WTIME() - tstart
    IF (par%master) call check_mesh( mesh)

    ! Smooth again
    tstart = MPI_WTIME()
    IF (par%master) call Lloyds_algorithm_single_iteration( mesh, alpha_min)
    tcomp = tcomp + MPI_WTIME() - tstart
    IF (par%master) call check_mesh( mesh)

    ! Broadcast from Master
    tstart = MPI_WTIME()
    call broadcast_mesh( mesh)
    tcomp = tcomp + MPI_WTIME() - tstart
    call check_mesh( mesh)

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    tstart = MPI_WTIME()
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    tcomp = tcomp + MPI_WTIME() - tstart

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! == Validation
    ! =============

    found_errors = .FALSE.

    IF (mesh%nV < NINT( real( nV_expected,dp) * 0.9_dp) .OR. mesh%nV > NINT( real( nV_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of vertices! Expected {int_01}, found {int_02}', int_01 = nV_expected, int_02 = mesh%nV)
    end IF

    IF (mesh%nTri < NINT( real( nTri_expected,dp) * 0.9_dp) .OR. mesh%nTri > NINT( real( nTri_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of triangles! Expected {int_01}, found {int_02}', int_01 = nTri_expected, int_02 = mesh%nTri)
    end IF

    IF (mesh%nE < NINT( real( nE_expected,dp) * 0.9_dp) .OR. mesh%nE > NINT( real( nE_expected,dp) * 1.1_dp)) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexepcted amount of edges! Expected {int_01}, found {int_02}', int_01 = nE_expected, int_02 = mesh%nE)
    end IF

    IF (MINVAL( mesh%A) < Amin_expected * 0.5_dp) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly small vertices! Expected MINVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amin_expected, dp_02 = MINVAL( mesh%A))
    end IF

    IF (MAXVAL( mesh%A) > Amax_expected * 2.0_dp) THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly small vertices! Expected MAXVAL( mesh%A) = {dp_01}, found {dp_02}', dp_01 = Amax_expected, dp_02 = MAXVAL( mesh%A))
    end IF

    RA_max = 0._dp
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        RA = MAX( mesh%A( vi) / mesh%A( vj), mesh%A( vj) / mesh%A( vi))
        RA_max = MAX( RA_max, RA)
      end DO
    end DO
    IF (RA_max > RA_max_expected * 2.0_dp)  THEN
      found_errors = .TRUE.
      IF (par%master) call warning('mesh has unexpectedly high resolution gradient! Expected RA_max = {dp_01}, found {dp_02}', dp_01 = RA_max_expected, dp_02 = RA_max)
    end IF

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('created a test mesh with {int_01} vertices, {int_02} triangles, and {int_03} edges in {dp_01} seconds', &
        int_01 = mesh%nV, int_02 = mesh%nTri, int_03 = mesh%nE, dp_01 = tcomp)
    ELSE
      IF (par%master) call warning('found errors in basic two-core mesh creation')
    end IF

    ! Write the resulting mesh to a NetCDF file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine test_mesh_creation_basic_two_cores_prime

  subroutine test_mesh_operators_basic( mesh)
    ! Test the basic mapping/gradient matrix operators

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_mesh_operators_basic'
    integer                                            :: vi,ti,ei,row
    real(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_a_ex, ddx_a_ex, ddy_a_ex
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_b_ex, ddx_b_ex, ddy_b_ex, d2dx2_b_ex, d2dxdy_b_ex, d2dy2_b_ex
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_c_ex, ddx_c_ex, ddy_c_ex
    real(dp), DIMENSION(:    ), ALLOCATABLE            ::        d_a_b, d_a_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_b_a,        d_b_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_c_a, d_c_b
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: ddx_a_a, ddx_a_b, ddx_a_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: ddx_b_a, ddx_b_b, ddx_b_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: ddx_c_a, ddx_c_b, ddx_c_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: ddy_a_a, ddy_a_b, ddy_a_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: ddy_b_a, ddy_b_b, ddy_b_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: ddy_c_a, ddy_c_b, ddy_c_c
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d2dx2_b_b, d2dxdy_b_b, d2dy2_b_b
    real(dp)                                           ::               maxerr_d_a_b, maxerr_d_a_c
    real(dp)                                           :: maxerr_d_b_a,               maxerr_d_b_c
    real(dp)                                           :: maxerr_d_c_a, maxerr_d_c_b
    real(dp)                                           :: maxerr_ddx_a_a, maxerr_ddx_a_b, maxerr_ddx_a_c
    real(dp)                                           :: maxerr_ddx_b_a, maxerr_ddx_b_b, maxerr_ddx_b_c
    real(dp)                                           :: maxerr_ddx_c_a, maxerr_ddx_c_b, maxerr_ddx_c_c
    real(dp)                                           :: maxerr_ddy_a_a, maxerr_ddy_a_b, maxerr_ddy_a_c
    real(dp)                                           :: maxerr_ddy_b_a, maxerr_ddy_b_b, maxerr_ddy_b_c
    real(dp)                                           :: maxerr_ddy_c_a, maxerr_ddy_c_b, maxerr_ddy_c_c
    real(dp)                                           :: maxerr_d2dx2_b_b, maxerr_d2dxdy_b_b, maxerr_d2dy2_b_b
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid

    ! Add routine to path
    call init_routine( routine_name)

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
      call test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      row = mesh%vi2n( vi)
      d_a_ex(   row) = d
      ddx_a_ex( row) = ddx
      ddy_a_ex( row) = ddy
    end DO

    ! b-grid (triangles)
    DO ti = mesh%ti1, mesh%ti2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      call test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      row = mesh%ti2n( ti)
      d_b_ex(      row) = d
      ddx_b_ex(    row) = ddx
      ddy_b_ex(    row) = ddy
      d2dx2_b_ex(  row) = d2dx2
      d2dxdy_b_ex( row) = d2dxdy
      d2dy2_b_ex(  row) = d2dy2
    end DO

    ! c-grid (edges)
    DO ei = mesh%ei1, mesh%ei2
      x = mesh%E( ei,1)
      y = mesh%E( ei,2)
      call test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      row = mesh%ei2n( ei)
      d_c_ex(   row) = d
      ddx_c_ex( row) = ddx
      ddy_c_ex( row) = ddy
    end DO

    ! Calculate discretised solutions
    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_a_b    , d_a_ex, d_a_b     )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_a_c    , d_a_ex, d_a_c     )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_b_a    , d_b_ex, d_b_a     )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_b_c    , d_b_ex, d_b_c     )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_c_a    , d_c_ex, d_c_a     )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_c_b    , d_c_ex, d_c_b     )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_a    , d_a_ex, ddx_a_a   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_b    , d_a_ex, ddx_a_b   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_c    , d_a_ex, ddx_a_c   )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_a    , d_b_ex, ddx_b_a   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_b    , d_b_ex, ddx_b_b   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_c    , d_b_ex, ddx_b_c   )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_c_a    , d_c_ex, ddx_c_a   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_c_b    , d_c_ex, ddx_c_b   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_c_c    , d_c_ex, ddx_c_c   )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_a    , d_a_ex, ddy_a_a   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_b    , d_a_ex, ddy_a_b   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_c    , d_a_ex, ddy_a_c   )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_a    , d_b_ex, ddy_b_a   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_b    , d_b_ex, ddy_b_b   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_c    , d_b_ex, ddy_b_c   )

    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_c_a    , d_c_ex, ddy_c_a   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_c_b    , d_c_ex, ddy_c_b   )
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_c_c    , d_c_ex, ddy_c_c   )

    call multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dx2_b_b , d_b_ex, d2dx2_b_b )
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dxdy_b_b, d_b_ex, d2dxdy_b_b)
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dy2_b_b , d_b_ex, d2dy2_b_b )

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

    end DO ! DO vi = mesh%vi1, mesh%vi2

    ! Find maximum errors across all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_b_a  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_c_a  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_a_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_b_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_c_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_a_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_b_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_c_a, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr_d_b_a > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_map_b_a')
    end IF
    IF (maxerr_d_c_a > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_map_c_a')
    end IF

    IF (maxerr_ddx_a_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_a_a')
    end IF
    IF (maxerr_ddx_b_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_b_a')
    end IF
    IF (maxerr_ddx_c_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_c_a')
    end IF

    IF (maxerr_ddy_a_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_a_a')
    end IF
    IF (maxerr_ddy_b_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_b_a')
    end IF
    IF (maxerr_ddy_c_a > 1.0E-7_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_c_a')
    end IF

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

    end DO ! DO ti = mesh%ti1, mesh%ti2

    ! Find maximum errors across all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_a_b     , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_c_b     , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_a_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_b_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_c_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_a_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_b_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_c_b   , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d2dx2_b_b , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d2dxdy_b_b, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d2dy2_b_b , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr_d_a_b > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_map_a_b')
    end IF
    IF (maxerr_d_c_b > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_map_c_b')
    end IF

    IF (maxerr_ddx_a_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_a_b')
    end IF
    IF (maxerr_ddx_b_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_b_b')
    end IF
    IF (maxerr_ddx_c_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_c_b')
    end IF

    IF (maxerr_ddy_a_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_a_b')
    end IF
    IF (maxerr_ddy_b_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_b_b')
    end IF
    IF (maxerr_ddy_c_b > 1.0E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_c_b')
    end IF

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

    end DO ! DO ei = mesh%ei1, mesh%ei2

    ! Find maximum errors across all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_a_c  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_d_b_c  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_a_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_b_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddx_c_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_a_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_b_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr_ddy_c_c, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr_d_a_c > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_map_a_c')
    end IF
    IF (maxerr_d_b_c > 0.25E-1_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_map_b_c')
    end IF

    IF (maxerr_ddx_a_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_a_c')
    end IF
    IF (maxerr_ddx_b_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_b_c')
    end IF
    IF (maxerr_ddx_c_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddx_c_c')
    end IF

    IF (maxerr_ddy_a_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_a_c')
    end IF
    IF (maxerr_ddy_b_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_b_c')
    end IF
    IF (maxerr_ddy_c_c > 0.25E-6_dp) THEN
      found_errors = .TRUE.
      call warning('inaccuracies found in M_ddy_c_c')
    end IF

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('validated all basic mesh matrix operators')
    ELSE
      IF (par%master) call warning('found errors in basic mesh matrix operators')
    end IF

    ! Write results to a NetCDF file

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'd_a_ex'     )
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_a_ex'   )
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_a_ex'   )

    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_b_ex'     )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_b_ex'   )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_b_ex'   )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dx2_b_ex' )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dxdy_b_ex')
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dy2_b_ex' )

    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'd_c_ex'     )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_c_ex'   )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_c_ex'   )

    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_a_b'      )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'd_a_c'      )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'd_b_a'      )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'd_b_c'      )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'd_c_a'      )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_c_b'      )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_a_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_a_b'    )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_a_c'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_b_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_b_b'    )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_b_c'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddx_c_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddx_c_b'    )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddx_c_c'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_a_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_a_b'    )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_a_c'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_b_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_b_b'    )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_b_c'    )

    call add_field_mesh_dp_2D_notime(   filename, ncid, 'ddy_c_a'    )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'ddy_c_b'    )
    call add_field_mesh_dp_2D_c_notime( filename, ncid, 'ddy_c_c'    )

    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dx2_b_b'  )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dxdy_b_b' )
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd2dy2_b_b'  )

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_a_ex'     , d_a_ex     )
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_ex'   , ddx_a_ex   )
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_ex'   , ddy_a_ex   )

    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_b_ex'     , d_b_ex     )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_ex'   , ddx_b_ex   )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_ex'   , ddy_b_ex   )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_ex' , d2dx2_b_ex )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_ex', d2dxdy_b_ex)
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_ex' , d2dy2_b_ex )

    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_c_ex'     , d_c_ex     )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_c_ex'   , ddx_c_ex   )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_c_ex'   , ddy_c_ex   )

    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_a_b'      , d_a_b      )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_a_c'      , d_a_c      )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_b_a'      , d_b_a      )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'd_b_c'      , d_b_c      )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'd_c_a'      , d_c_a      )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_c_b'      , d_c_b      )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_a_a'    , ddx_a_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_a_b'    , ddx_a_b    )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_a_c'    , ddx_a_c    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_b_a'    , ddx_b_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_b_b'    , ddx_b_b    )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_b_c'    , ddx_b_c    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddx_c_a'    , ddx_c_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddx_c_b'    , ddx_c_b    )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddx_c_c'    , ddx_c_c    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_a_a'    , ddy_a_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_a_b'    , ddy_a_b    )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_a_c'    , ddy_a_c    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_b_a'    , ddy_b_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_b_b'    , ddy_b_b    )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_b_c'    , ddy_b_c    )

    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'ddy_c_a'    , ddy_c_a    )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'ddy_c_b'    , ddy_c_b    )
    call write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, 'ddy_c_c'    , ddy_c_c    )

    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dx2_b_b'  , d2dx2_b_b  )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dxdy_b_b' , d2dxdy_b_b )
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd2dy2_b_b'  , d2dy2_b_b  )

    ! Close the file
    call close_netcdf_file( ncid)

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
    call finalise_routine( routine_name)

  end subroutine test_mesh_operators_basic

  subroutine test_mesh_operators_3D( mesh)
    ! Test the 3-D gradient matrix operators

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_mesh_operators_3D'
    type(type_ice_model)                               :: ice
    integer                                            :: vi,ti,k
    real(dp)                                           :: x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: f_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dx_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dy_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dz_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dx2_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dxdy_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dy2_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dz2_bk_ex
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dx_bk
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dy_bk
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: df_dz_bk
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dx2_bk
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dxdy_bk
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dy2_bk
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d2f_dz2_bk
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate ice model memory
    call allocate_ice_model( mesh, ice)

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

      call test_function_3D( x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2)

      ice%Hi( vi) = Hi
      ice%Hb( vi) = Hb
      ice%Hs( vi) = Hs
      ice%SL( vi) = SL

    end DO
    end DO

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

      call test_function_3D( x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2)

      f_bk_ex(        ti,k) = f
      df_dx_bk_ex(    ti,k) = dfdx
      df_dy_bk_ex(    ti,k) = dfdy
      df_dz_bk_ex(    ti,k) = dfdz
      d2f_dx2_bk_ex(  ti,k) = d2fdx2
      d2f_dxdy_bk_ex( ti,k) = d2fdxdy
      d2f_dy2_bk_ex(  ti,k) = d2fdy2
      d2f_dz2_bk_ex(  ti,k) = d2fdz2

    end DO
    end DO

    ! Calculate gradients of zeta
    call calc_zeta_gradients( mesh, ice)

    ! Calculate 3-D gradient operators
    call calc_3D_matrix_operators_mesh( mesh, ice)

    ! Calculate gradients of f

    ALLOCATE( df_dx_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( df_dy_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( df_dz_bk(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dx2_bk(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dxdy_bk( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dy2_bk(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( d2f_dz2_bk(  mesh%ti1:mesh%ti2,mesh%nz))

    call calc_3D_gradient_bk_bk( mesh, mesh%M2_ddx_bk_bk   , f_bk_ex, df_dx_bk   )
    call calc_3D_gradient_bk_bk( mesh, mesh%M2_ddy_bk_bk   , f_bk_ex, df_dy_bk   )
    call calc_3D_gradient_bk_bk( mesh, mesh%M2_ddz_bk_bk   , f_bk_ex, df_dz_bk   )
    call calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dx2_bk_bk , f_bk_ex, d2f_dx2_bk )
    call calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dxdy_bk_bk, f_bk_ex, d2f_dxdy_bk)
    call calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dy2_bk_bk , f_bk_ex, d2f_dy2_bk )
    call calc_3D_gradient_bk_bk( mesh, mesh%M2_d2dz2_bk_bk , f_bk_ex, d2f_dz2_bk )

    ! If no errors occurred, we are happy
    found_errors = .FALSE.
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('validated all 3-D mesh matrix operators')
    ELSE
      IF (par%master) call warning('found errors in 3-D mesh matrix operators')
    end IF

    ! Write results to a NetCDF file

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_zeta_dimension_to_file( filename, ncid, mesh%zeta)

    ! Add all the variables
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'Hi')
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'Hb')
    call add_field_mesh_dp_2D_notime(   filename, ncid, 'Hs')

    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'f_bk_ex'       )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dx_bk_ex'   )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dy_bk_ex'   )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dz_bk_ex'   )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dx2_bk_ex' )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dxdy_bk_ex')
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dy2_bk_ex' )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dz2_bk_ex' )

    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dx_bk'   )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dy_bk'   )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'df_dz_bk'   )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dx2_bk' )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dxdy_bk')
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dy2_bk' )
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd2f_dz2_bk' )

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'Hi'            , ice%Hi        )
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'Hb'            , ice%Hb        )
    call write_to_field_multopt_mesh_dp_2D_notime(   mesh, filename, ncid, 'Hs'            , ice%Hs        )

    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'f_bk_ex'       , f_bk_ex       )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dx_bk_ex'   , df_dx_bk_ex   )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dy_bk_ex'   , df_dy_bk_ex   )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dz_bk_ex'   , df_dz_bk_ex   )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dx2_bk_ex' , d2f_dx2_bk_ex )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dxdy_bk_ex', d2f_dxdy_bk_ex)
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dy2_bk_ex' , d2f_dy2_bk_ex )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dz2_bk_ex' , d2f_dz2_bk_ex )

    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dx_bk'      , df_dx_bk      )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dy_bk'      , df_dy_bk      )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'df_dz_bk'      , df_dz_bk      )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dx2_bk'    , d2f_dx2_bk    )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dxdy_bk'   , d2f_dxdy_bk   )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dy2_bk'    , d2f_dy2_bk    )
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd2f_dz2_bk'    , d2f_dz2_bk    )

    ! Close the file
    call close_netcdf_file( ncid)

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
    call finalise_routine( routine_name)

  end subroutine test_mesh_operators_3D

  subroutine test_remapping_grid2mesh( mesh)
    ! Test remapping from a square grid to a mesh

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_remapping_grid2mesh'
    character(len=1024)                                :: name
    real(dp)                                           :: dx
    type(type_grid)                                    :: grid
    integer                                            :: i,j
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial
    real(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_ex_partial
    integer                                            :: vi, vi_glob
    real(dp)                                           :: maxerr
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! == Set up a square grid
    ! =======================

    name = 'test_grid'
    dx   = 32E3_dp
    call setup_square_grid( name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, dx, grid, lambda_M = mesh%lambda_M, phi_M = mesh%phi_M, beta_stereo = mesh%beta_stereo)

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
        call test_function( x, y, grid%xmin, grid%xmax, grid%ymin, grid%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_grid( i,j) = d
      end DO
      end DO
    end IF
    call sync

    ! Distribute gridded data
    ALLOCATE( d_grid_vec_partial( grid%n_loc))
    call distribute_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
    IF (par%master) DEALLOCATE( d_grid)

    ! Map gridded data to the mesh
    ALLOCATE( d_mesh_partial( mesh%nV_loc))
    call map_from_xy_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Calculate exact solution
    ALLOCATE( d_mesh_ex_partial( mesh%nV_loc))
    DO vi = 1, mesh%nV_loc
      vi_glob = mesh%vi1 + vi - 1
      x = mesh%V( vi_glob,1)
      y = mesh%V( vi_glob,2)
      call test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_mesh_ex_partial( vi) = d
    end DO

    ! Calculate error
    maxerr = 0._dp
    DO vi = 1, mesh%nV_loc
      ! Skip border vertices
      vi_glob = mesh%vi1 + vi - 1
      IF (mesh%VBI( vi_glob) > 0) CYCLE
      maxerr = MAX( maxerr, ABS( d_mesh_partial( vi) - d_mesh_ex_partial( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr > 0.6E-1_dp) found_errors = .TRUE.

    ! == Validation
    ! =============

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('validated grid to mesh remapping')
    ELSE
      IF (par%master) call warning('found errors in grid to mesh remapping')
    end IF

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh_ex')

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh'   , d_mesh_partial   )
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex_partial)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call clear_all_maps_involving_this_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine test_remapping_grid2mesh

  subroutine test_remapping_mesh2grid( mesh)
    ! Test remapping from a mesh to a square grid

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_remapping_mesh2grid'
    character(len=1024)                                :: name
    real(dp)                                           :: dx
    type(type_grid)                                    :: grid
    integer                                            :: i,j,n,n_glob
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial
    integer                                            :: vi, vi_glob
    real(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial, d_grid_ex_vec_partial
    real(dp)                                           :: maxerr
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! == Set up a square grid
    ! =======================

    name = 'test_grid'
    dx   = 32E3_dp
    call setup_square_grid( name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, dx, grid, lambda_M = mesh%lambda_M, phi_M = mesh%phi_M, beta_stereo = mesh%beta_stereo)

    ! == Calculate, apply, and validate grid-to-mesh remapping operator
    ! =================================================================

    found_errors = .FALSE.

    ! Create a nice data field on the mesh
    ALLOCATE( d_mesh_partial( mesh%nV_loc), source = 0._dp)
    DO vi = 1, mesh%nV_loc
      vi_glob = mesh%vi1 + vi - 1
      x = mesh%V( vi_glob,1)
      y = mesh%V( vi_glob,2)
      call test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_mesh_partial( vi) = d
    end DO

    ! Map meshed data to the grid
    ALLOCATE( d_grid_vec_partial( grid%n_loc))
    call map_from_mesh_to_xy_grid_2D( mesh, grid, d_mesh_partial, d_grid_vec_partial)

    ! Calculate exact solution on the grid
    ALLOCATE( d_grid_ex_vec_partial( grid%n_loc))
    DO n = 1, grid%n_loc
      n_glob = grid%n1 + n - 1
      i = grid%n2ij( n_glob,1)
      j = grid%n2ij( n_glob,2)
      x = grid%x( i)
      y = grid%y( j)
      call test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_grid_ex_vec_partial( n) = d
    end DO

    ! Calculate error
    maxerr = 0._dp
    DO n = 1, grid%n_loc
      n_glob = grid%n1 + n - 1
      i = grid%n2ij( n_glob,1)
      j = grid%n2ij( n_glob,2)
      ! Skip the border
      IF (i == 1 .OR. i == grid%nx .OR. j == 1 .OR. j == grid%ny) CYCLE
      maxerr = MAX( maxerr, ABS( d_grid_vec_partial( n) - d_grid_ex_vec_partial( n)))
    end DO

    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr > 0.6E-1_dp) found_errors = .TRUE.

    ! == Validation
    ! =============

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('validated mesh to grid remapping')
    ELSE
      IF (par%master) call warning('found errors in mesh to grid remapping')
    end IF

    ! Create a file and write the mesh to it
    filename = TRIM(C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    ! Add all the variables
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid')
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid_ex')

    ! Write all the variables
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid'   , d_grid_vec_partial   )
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid_ex', d_grid_ex_vec_partial)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call clear_all_maps_involving_this_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine test_remapping_mesh2grid

  subroutine test_remapping_lonlat2mesh( mesh)
    ! Test remapping from a lon/lat-grid to a mesh

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_remapping_lonlat2mesh'
    character(len=1024)                                :: name
    integer                                            :: nlon, nlat
    type(type_grid_lonlat)                             :: grid
    integer                                            :: i,j
    real(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial
    real(dp)                                           :: lon,lat,d
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_partial
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d_mesh_ex_partial
    integer                                            :: vi, vi_glob
    real(dp)                                           :: maxerr
    logical                                            :: found_errors
    character(len=1024)                                :: filename
    integer                                            :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! == Set up a square grid
    ! =======================

    ! Set up a lon/lat-grid
    name = 'test_lonlat_grid'
    nlon = 360
    nlat = 180
    call setup_simple_lonlat_grid( name, nlon, nlat, grid)

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
        call test_function_lonlat( lon, lat, d)
        d_grid( i,j) = d
      end DO
      end DO
    end IF
    call sync

    ! Distribute gridded data
    ALLOCATE( d_grid_vec_partial( grid%n_loc))
    call distribute_lonlat_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
    IF (par%master) DEALLOCATE( d_grid)

    ! Map gridded data to the mesh
    ALLOCATE( d_mesh_partial( mesh%nV_loc))
    call map_from_lonlat_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Calculate exact solution
    ALLOCATE( d_mesh_ex_partial( mesh%nV_loc))
    DO vi = 1, mesh%nV_loc
      vi_glob = mesh%vi1 + vi - 1
      lon = mesh%lon( vi_glob)
      lat = mesh%lat( vi_glob)
      call test_function_lonlat( lon, lat, d)
      d_mesh_ex_partial( vi) = d
    end DO

    ! Calculate error
    maxerr = 0._dp
    DO vi = 1, mesh%nV_loc
      ! Skip border vertices
      vi_glob = mesh%vi1 + vi - 1
      IF (mesh%VBI( vi_glob) > 0) CYCLE
      maxerr = MAX( maxerr, ABS( d_mesh_partial( vi) - d_mesh_ex_partial( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    IF (maxerr > 0.5E-4_dp) found_errors = .TRUE.

    ! == Validation
    ! =============

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('validated lon/lat to mesh remapping')
    ELSE
      IF (par%master) call warning('found errors in lon/lat to mesh remapping')
    end IF

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add all the variables
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh_ex')

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh'   , d_mesh_partial   )
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex_partial)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call clear_all_maps_involving_this_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine test_remapping_lonlat2mesh

  subroutine test_remapping_mesh2mesh( mesh1, mesh2)
    ! Test remapping from a mesh to another mesh

    implicit none

    ! In/output variables:
    type(type_mesh),                     INTENT(INOUT) :: mesh1, mesh2

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'test_remapping_mesh2mesh'
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d1_ex, d2_ex
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d12_nn, d12_trilin, d12_cons
    real(dp), DIMENSION(:    ), ALLOCATABLE            :: d21_nn, d21_trilin, d21_cons
    integer                                            :: vi
    real(dp)                                           :: x,y,d,ddx,ddy,d2dx2,d2dxdy,d2dy2
    character(len=1024)                                :: method
    logical                                            :: found_errors
    real(dp)                                           :: maxerr
    character(len=1024)                                :: filename
    integer                                            :: ncid

    ! Add routine to path
    call init_routine( routine_name)

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
      call test_function( x, y, mesh1%xmin, mesh1%xmax, mesh1%ymin, mesh1%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d1_ex( vi) = d
    end DO
    DO vi = mesh2%vi1, mesh2%vi2
      x = mesh2%V( vi,1)
      y = mesh2%V( vi,2)
      call test_function( x, y, mesh2%xmin, mesh2%xmax, mesh2%ymin, mesh2%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d2_ex( vi) = d
    end DO

    ! == Remap data
    ! =============

    method = 'nearest_neighbour'
    call map_from_mesh_to_mesh_2D( mesh1, mesh2, d1_ex, d12_nn, method)
    call map_from_mesh_to_mesh_2D( mesh2, mesh1, d2_ex, d21_nn, method)

    method = 'trilin'
    call map_from_mesh_to_mesh_2D( mesh1, mesh2, d1_ex, d12_trilin, method)
    call map_from_mesh_to_mesh_2D( mesh2, mesh1, d2_ex, d21_trilin, method)

    method = '2nd_order_conservative'
    call map_from_mesh_to_mesh_2D( mesh1, mesh2, d1_ex, d12_cons, method)
    call map_from_mesh_to_mesh_2D( mesh2, mesh1, d2_ex, d21_cons, method)

    ! == Validation
    ! =============

    found_errors = .FALSE.

    ! Nearest-neighbour
    ! 2023-06-05: maxerr = 0.11995E+00
    maxerr = 0._dp
    DO vi = mesh2%vi1, mesh2%vi2
      IF (mesh2%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d12_nn( vi) - d2_ex( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.11995E+00 * 2.0) THEN
      call warning('unexpectedly high local errors detected in mesh-to-mesh nearest-neighbour remapping; expected 0.11995E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    end IF

    ! 2023-06-05: maxerr = 0.16505E+00
    maxerr = 0._dp
    DO vi = mesh1%vi1, mesh1%vi2
      IF (mesh1%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d21_nn( vi) - d1_ex( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.16505E+00 * 2.0) THEN
      call warning('unexpectedly high local errors detected in mesh-to-mesh nearest-neighbour remapping: expected 0.16505E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    end IF

    ! Trilinear
    ! 2023-06-05: maxerr = 0.16466E-01
    maxerr = 0._dp
    DO vi = mesh2%vi1, mesh2%vi2
      IF (mesh2%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d12_trilin( vi) - d2_ex( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.16466E-01 * 2.0) THEN
      call warning('unexpectedly high local errors detected in mesh-to-mesh trilinear remapping: expected 0.16466E-01, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    end IF

    ! 2023-06-05: maxerr = 0.18900E-01
    maxerr = 0._dp
    DO vi = mesh1%vi1, mesh1%vi2
      IF (mesh1%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d21_trilin( vi) - d1_ex( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.18900E-01 * 2.0) THEN
      call warning('unexpectedly high local errors detected in mesh-to-mesh trilinear remapping: expected 0.18900E-01, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    end IF

    ! 2nd-order conservative
    ! 2023-06-05: maxerr = 0.27873E+00
    maxerr = 0._dp
    DO vi = mesh2%vi1, mesh2%vi2
      IF (mesh2%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d12_cons( vi) - d2_ex( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.27873E+00 * 2.0) THEN
      call warning('unexpectedly high local errors detected in mesh-to-mesh 2nd-order conservative remapping: expected 0.27873E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    end IF

    ! 2023-06-05: maxerr = 0.56043E+00
    maxerr = 0._dp
    DO vi = mesh1%vi1, mesh1%vi2
      IF (mesh1%VBI( vi) > 0) CYCLE ! Skip border vertices, we know remapping is not always accurate there
      maxerr = MAX( maxerr, ABS( d21_cons( vi) - d1_ex( vi)))
    end DO
    call MPI_ALLREDUCE( MPI_IN_PLACE, maxerr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (maxerr > 0.56043E+00 * 2.0) THEN
      call warning('unexpectedly high local errors detected in mesh-to-mesh 2nd-order conservative remapping: expected 0.56043E+00, found {dp_01}', dp_01 = maxerr)
      found_errors = .TRUE.
    end IF

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) call happy('validated mesh to mesh remapping')
    ELSE
      IF (par%master) call warning('found errors in mesh to mesh remapping')
    end IF

    ! Write results to a NetCDF file

    ! Mesh 1
    ! ======

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output1.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh1)

    ! Add all the variables
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd1_ex')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd21_nn')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd21_trilin')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd21_cons')

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd1_ex'     , d1_ex     )
    call write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd21_nn'    , d21_nn    )
    call write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd21_trilin', d21_trilin)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh1, filename, ncid, 'd21_cons'  , d21_cons  )

    ! Close the file
    call close_netcdf_file( ncid)

    ! Mesh 2
    ! ======

    ! Create a file and write the mesh to it
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output2.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh2)

    ! Add all the variables
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd2_ex')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd12_nn')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd12_trilin')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd12_cons')

    ! Write all the variables
    call write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd2_ex'     , d2_ex     )
    call write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd12_nn'    , d12_nn    )
    call write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd12_trilin', d12_trilin)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd12_cons'  , d12_cons  )

    ! Close the file
    call close_netcdf_file( ncid)

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
    call clear_all_maps_involving_this_mesh( mesh1)
    call clear_all_maps_involving_this_mesh( mesh2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine test_remapping_mesh2mesh

  subroutine test_function( x, y, xmin, xmax, ymin, ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
    ! A simple test function to validate matrix operators

    implicit none

    ! In/output variables:
    real(dp),                   INTENT(IN)        :: x,y
    real(dp),                   INTENT(IN)        :: xmin,xmax,ymin,ymax
    real(dp),                   INTENT(OUT)       :: d,ddx,ddy
    real(dp),                   INTENT(OUT)       :: d2dx2,d2dxdy,d2dy2

    ! Local variables:
    real(dp)                                      :: c1,c2

    c1 = 2._dp * pi / (xmax - xmin)
    c2 = 3._dp * pi / (ymax - ymin)

    d      =            SIN( c1 * (x - xmin)) *            SIN( c2 * (y - ymin))
    ddx    =   c1     * COS( c1 * (x - xmin)) *            SIN( c2 * (y - ymin))
    ddy    =            SIN( c1 * (x - xmin)) *   c2     * COS( c2 * (y - ymin))
    d2dx2  = (-c1)**2 * SIN( c1 * (x - xmin)) *            SIN( c2 * (y - ymin))
    d2dxdy =   c1     * COS( c1 * (x - xmin)) *   c2     * COS( c2 * (y - ymin))
    d2dy2  =            SIN( c1 * (x - xmin)) * (-c2)**2 * SIN( c2 * (y - ymin))

  end subroutine test_function

  subroutine test_function_lonlat( lon, lat, d)
    ! A simple test function to validate lon/lat gridding stuff

    implicit none

    ! In/output variables:
    real(dp),                   INTENT(IN)        :: lon,lat
    real(dp),                   INTENT(OUT)       :: d

    d = SIN( lon * pi / 180._dp) * COS( lat * pi / 180)

  end subroutine test_function_lonlat

  subroutine test_function_3D( x, y, xmin, xmax, ymin, ymax, zeta, Hi, Hb, Hs, SL, f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2)
    ! A simple test function to validate the 3-D matrix operators

    implicit none

    ! In/output variables:
    real(dp),                   INTENT(IN)        :: x,y,xmin,xmax,ymin,ymax,zeta
    real(dp),                   INTENT(OUT)       :: Hi, Hb, Hs, SL
    real(dp),                   INTENT(OUT)       :: f, dfdx, dfdy, dfdz, d2fdx2, d2fdxdy, d2fdy2, d2fdz2

    ! Local variables:
    real(dp), parameter                           :: a  = 500._dp  ! Amplitude of bedrock undulations
    real(dp), parameter                           :: h0 = 2000._dp ! Uniform surface elevation
    real(dp)                                      :: cx,cy,z,cz

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

  end subroutine test_function_3D

end module unit_tests_mesh_old
