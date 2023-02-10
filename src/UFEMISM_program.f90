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
  USE mesh_creation

  IMPLICIT NONE

! ===== Main variables =====
! ==========================

  TYPE(type_mesh) :: mesh
  CHARACTER(LEN=256) :: name
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: line, poly
  INTEGER :: i,n
  REAL(dp) :: r, theta1, theta2, alpha_min
  REAL(dp), DIMENSION(2) :: p
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

  n = 50
  r = 0.75_dp
  ALLOCATE( line( n,4))
  DO i = 1, n
    theta1 = 2._dp * pi * REAL( i-1,dp) / REAL( n-1,dp)
    theta2 = 2._dp * pi * REAL( i  ,dp) / REAL( n-1,dp)
    line( i,:) = [r * COS( theta1), r * SIN( theta1), r * COS( theta2), r * SIN( theta2)]
  END DO
  CALL sync
  IF (par%master) CALL refine_mesh_line( mesh, line, 0.01_dp, alpha_min)
  DEALLOCATE( line)

  n = 30
  r = 0.4_dp
  ALLOCATE( line( n,4))
  DO i = 1, n
    theta1 = -2._dp * pi * REAL( i-1,dp) / REAL( 2*n-1,dp)
    theta2 = -2._dp * pi * REAL( i  ,dp) / REAL( 2*n-1,dp)
    line( i,:) = [r * COS( theta1), r * SIN( theta1), r * COS( theta2), r * SIN( theta2)]
  END DO
  CALL sync
  IF (par%master) CALL refine_mesh_line( mesh, line, 0.01_dp, alpha_min)
  DEALLOCATE( line)

  p = [-0.3_dp, 0.4_dp]
  IF (par%master) CALL refine_mesh_point( mesh, p, 0.01_dp, alpha_min)
  p = [ 0.3_dp, 0.4_dp]
  IF (par%master) CALL refine_mesh_point( mesh, p, 0.01_dp, alpha_min)

  n = 100
  ALLOCATE( line( n,4), source = 0._dp)
  line(  1,:) = [-0.80_dp,  0.20_dp, -0.80_dp, -0.20_dp]
  line(  2,:) = [-0.80_dp, -0.20_dp, -0.65_dp, -0.20_dp]
  line(  3,:) = [-0.65_dp, -0.20_dp, -0.65_dp,  0.20_dp]
  line(  4,:) = [-0.55_dp,  0.20_dp, -0.55_dp, -0.20_dp]
  line(  5,:) = [-0.55_dp,  0.20_dp, -0.40_dp,  0.20_dp]
  line(  6,:) = [-0.55_dp,  0.00_dp, -0.40_dp,  0.00_dp]
  IF (par%master) CALL refine_mesh_line( mesh, line, 0.01_dp, alpha_min)
  line(  7,:) = [-0.30_dp,  0.20_dp, -0.30_dp, -0.20_dp]
  IF (par%master) CALL refine_mesh_line( mesh, line, 0.01_dp, alpha_min)
  line(  8,:) = [-0.30_dp,  0.20_dp, -0.15_dp,  0.20_dp]
  line(  9,:) = [-0.30_dp,  0.00_dp, -0.15_dp,  0.00_dp]
  line( 10,:) = [-0.30_dp, -0.20_dp, -0.15_dp, -0.20_dp]
  line( 11,:) = [-0.05_dp,  0.20_dp, -0.05_dp, -0.20_dp]
  line( 12,:) = [-0.05_dp,  0.20_dp,  0.05_dp,  0.00_dp]
  line( 13,:) = [ 0.05_dp,  0.00_dp,  0.15_dp,  0.20_dp]
  line( 14,:) = [ 0.15_dp,  0.20_dp,  0.15_dp, -0.20_dp]
  line( 15,:) = [ 0.25_dp,  0.20_dp,  0.25_dp, -0.20_dp]
  line( 16,:) = [ 0.35_dp,  0.20_dp,  0.50_dp,  0.20_dp]
  line( 17,:) = [ 0.35_dp,  0.20_dp,  0.35_dp,  0.00_dp]
  IF (par%master) CALL refine_mesh_line( mesh, line, 0.01_dp, alpha_min)
  line( 18,:) = [ 0.35_dp,  0.00_dp,  0.50_dp,  0.00_dp]
  line( 19,:) = [ 0.50_dp,  0.00_dp,  0.50_dp, -0.20_dp]
  line( 20,:) = [ 0.35_dp, -0.20_dp,  0.50_dp, -0.20_dp]
  line( 21,:) = [ 0.60_dp,  0.20_dp,  0.60_dp, -0.20_dp]
  line( 22,:) = [ 0.60_dp,  0.20_dp,  0.70_dp,  0.00_dp]
  line( 23,:) = [ 0.70_dp,  0.00_dp,  0.80_dp,  0.20_dp]
  line( 24,:) = [ 0.80_dp,  0.20_dp,  0.80_dp, -0.20_dp]
  IF (par%master) CALL refine_mesh_line( mesh, line, 0.01_dp, alpha_min)
  DEALLOCATE( line)


  ! Polygon
  ALLOCATE( poly( 4,2))
  poly( 1,:) = [0.80_dp, -0.80_dp]
  poly( 2,:) = [0.80_dp, -0.55_dp]
  poly( 3,:) = [0.55_dp, -0.60_dp]
  poly( 4,:) = [0.60_dp, -0.80_dp]
  IF (par%master) CALL refine_mesh_polygon( mesh, poly, 0.01_dp, alpha_min)
  DEALLOCATE( poly)


  tstop = MPI_WTIME()

  IF (par%master) CALL warning('Finished generation mesh of {int_01} vertices in {dp_01} seconds.', int_01 = mesh%nV, dp_01 = tstop - tstart)


  IF (par%master) CALL write_mesh_to_text_file( mesh, 'test_mesh.txt')




! ===== FINISH =====
! ==================

  ! Finalise MPI parallelisation
  CALL finalise_parallelisation

END PROGRAM UFEMISM_program