PROGRAM UFEMISM_program

  ! Hieperdepiep!

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync, initialise_parallelisation, finalise_parallelisation
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine, initialise_control_and_resource_tracker
  USE UFEMISM_main_model                                     , ONLY: run_model

  USE parameters
  USE mesh_types
  USE mesh_memory
  USE mesh_utilities
  USE mesh_Delaunay
  USE mesh_creation
  USE basic_data_types
  USE math_utilities
  USE mesh_edges
  USE mesh_secondary
  USE netcdf_basic
  USE netcdf_input
  USE netcdf_output

  IMPLICIT NONE

! ===== Main variables =====
! ==========================

  TYPE(type_mesh) :: mesh
  CHARACTER(LEN=256) :: name, filename
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: poly, line_GL
  REAL(dp) :: alpha_min
  REAL(dp) :: tstart, tstop
  TYPE(type_grid) :: grid
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Hi, Hb, Hs, TAF
  INTEGER :: i,j
  INTEGER :: ncid
  REAL(dp) :: lambda_M, phi_M, beta_stereo

! ===== START =====
! =================

  ! Initialise MPI parallelisation
  CALL initialise_parallelisation

  ! Initialise the control and resource tracker
  CALL initialise_control_and_resource_tracker




  ! == READ ICE GEOMETRY ==
  ! =======================

  filename = '/Users/berends/Documents/Datasets/BedMachine_Antarctica/Bedmachine_v1_Antarctica_40km.nc'
  CALL read_field_from_xy_file_2D( filename, 'default_options_Hi', Hi, grid = grid)
  CALL read_field_from_xy_file_2D( filename, 'default_options_Hb', Hb)
  CALL read_field_from_xy_file_2D( filename, 'default_options_Hs', Hs)

  ! Calculate grounding line
  ALLOCATE( TAF( grid%nx, grid%ny))
  DO i = 1, grid%nx
  DO j = 1, grid%ny
    TAF( i,j) = thickness_above_floatation( Hi( i,j), Hb( i,j), 0._dp)
  END DO
  END DO
  CALL calc_grid_contour_as_line( grid, TAF, 0._dp, line_GL)



  ! == GENERATE MESH ==
  ! ===================


  tstart = MPI_WTIME()

  ! DENK DROM
  name = 'test_mesh'
  CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)
  CALL initialise_dummy_mesh( mesh, grid%xmin, grid%xmax, grid%ymin, grid%ymax)
  alpha_min = 25._dp * pi / 180._dp

  ! Move central vertex
  CALL move_vertex( mesh, 5, [1E4_dp, 3.7E4_dp])

  ! Uniform
  CALL refine_mesh_uniform( mesh, 4E5_dp, alpha_min)

  ! Smooth
  CALL Lloyds_algorithm_single_iteration( mesh)

  ! Grounding line
  CALL refine_mesh_line( mesh, line_GL, 5E4_dp, alpha_min)

  ! Smileyface
  CALL mesh_add_smileyface( mesh, 5E4_dp)

  ! UFEMISM letters
  CALL mesh_add_UFEMISM_letters( mesh, 5E4_dp)

  ! Smooth again
  CALL Lloyds_algorithm_single_iteration( mesh)

  ! Calculate all secondary mesh data
  lambda_M    = 0._dp
  phi_M       = -90._dp
  beta_stereo = 71._dp
  CALL calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

  tstop = MPI_WTIME()
  CALL warning('Finished generation mesh of {int_01} vertices in {dp_01} seconds.', int_01 = mesh%nV, dp_01 = tstop - tstart)

  CALL write_mesh_to_text_file( mesh, 'test_mesh.txt')



  ! == WRITE TO NETCDF ==
  ! =====================

  filename = 'test_mesh_netcdf.nc'
  CALL create_new_netcdf_file_for_writing( filename, ncid)
  CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)
  CALL close_netcdf_file( ncid)



! ===== FINISH =====
! ==================

  ! Finalise MPI parallelisation
  CALL finalise_parallelisation

END PROGRAM UFEMISM_program