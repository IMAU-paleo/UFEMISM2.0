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

  IMPLICIT NONE

! ===== Main variables =====
! ==========================

  TYPE(type_mesh) :: mesh
  CHARACTER(LEN=256) :: name
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: poly
  REAL(dp) :: alpha_min
  REAL(dp) :: tstart, tstop

! ===== START =====
! =================

  ! Initialise MPI parallelisation
  CALL initialise_parallelisation

  ! Initialise the control and resource tracker
  CALL initialise_control_and_resource_tracker



  tstart = MPI_WTIME()

  ! DENK DROM
  name = 'test_mesh'
  CALL allocate_mesh_primary( mesh, name, 1000, 2000, 32)
  CALL initialise_dummy_mesh( mesh, -1._dp, 1._dp, -1._dp, 1._dp)
  alpha_min = 25._dp * pi / 180._dp

  ! Move central vertex
  IF (par%master) CALL move_vertex( mesh, 5, [0.01_dp, 0.037_dp])

  ! Uniform
  IF (par%master) CALL refine_mesh_uniform( mesh, 0.1_dp, alpha_min)

  ! Smooth
  IF (par%master) CALL Lloyds_algorithm_single_iteration( mesh)

  ! Smileyface
  CALL mesh_add_smileyface( mesh, 0.005_dp)

  ! UFEMISM letters
  CALL mesh_add_UFEMISM_letters( mesh, 0.005_dp)

  ! Polygon
  ALLOCATE( poly( 4,2))
  poly( 1,:) = [0.80_dp, -0.80_dp]
  poly( 2,:) = [0.80_dp, -0.55_dp]
  poly( 3,:) = [0.55_dp, -0.60_dp]
  poly( 4,:) = [0.60_dp, -0.80_dp]
  IF (par%master) CALL refine_mesh_polygon( mesh, poly, 0.01_dp, alpha_min)
  DEALLOCATE( poly)

  ! Smooth again
  IF (par%master) CALL Lloyds_algorithm_single_iteration( mesh)

  tstop = MPI_WTIME()
  IF (par%master) CALL warning('Finished generation mesh of {int_01} vertices in {dp_01} seconds.', int_01 = mesh%nV, dp_01 = tstop - tstart)

  IF (par%master) CALL write_mesh_to_text_file( mesh, 'test_mesh.txt')




! ===== FINISH =====
! ==================

  ! Finalise MPI parallelisation
  CALL finalise_parallelisation

END PROGRAM UFEMISM_program