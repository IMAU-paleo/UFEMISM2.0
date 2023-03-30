MODULE mpi_distributed_memory

  ! Some routine to work with distributed memory in the MPI parallelised architecture.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine

  IMPLICIT NONE

  INTERFACE allgather_array
    PROCEDURE :: allgather_array_dp
    PROCEDURE :: allgather_array_dp_2d
    PROCEDURE :: allgather_array_int
  END INTERFACE

CONTAINS

! ===== Subroutines =====
! =======================

  ! Partition a list of ntot elements over the n processes
  SUBROUTINE partition_list( ntot, i, n, i1, i2)
    ! Partition a list into parallel ranges (e.g. vertex domains)

    ! In/output variables:
    INTEGER,                    INTENT(IN)        :: ntot, i, n
    INTEGER,                    INTENT(OUT)       :: i1, i2
    integer                                       :: remainder, slice_size

    !it must be calculated exactly as petsc does
    !TODO it would be better to get this ranges from PetsC then ...
    IF (ntot > n*2) THEN
      remainder  = MOD( ntot,n)
      slice_size = ntot/n

      i1 = slice_size*i     + MIN( i  , remainder) + 1
      i2 = slice_size*(i+1) + MIN( i+1, remainder)
    ELSE
      IF (i==0) THEN
        i1 = 1
        i2 = ntot
      ELSE
        i1 = 1
        i2 = 0
      END IF
    END IF

  END SUBROUTINE partition_list

! ===== Gather distributed variables to the Master =====
! ======================================================

  SUBROUTINE gather_to_master_int_1D( d_partial, d_tot)
    ! Gather a distributed 1-D integer variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:      ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:      ), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_int_1D'
    INTEGER                                                            :: n1,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot))
    ELSE
      ALLOCATE( d_tot( 0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    CALL MPI_GATHERV( d_partial, n1, MPI_INTEGER, d_tot, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_int_1D

  SUBROUTINE gather_to_master_int_2D( d_partial, d_tot)
    ! Gather a distributed 2-D integer variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:    ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:,:    ), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_int_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2))
    ELSE
      ALLOCATE( d_tot( 0,0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
      CALL MPI_GATHERV( d_partial( :,j), n1, MPI_INTEGER, d_tot( :,j), counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_int_2D

  SUBROUTINE gather_to_master_int_3D( d_partial, d_tot)
    ! Gather a distributed 3-D integer variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:,:,:  ), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_int_3D'
    INTEGER                                                            :: n1,n2,n3,i,n2_proc,n3_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2, n3))
    ELSE
      ALLOCATE( d_tot( 0,0,0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
      CALL MPI_GATHERV( d_partial( :,j,k), n1, MPI_INTEGER, d_tot( :,j,k), counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_int_3D

  SUBROUTINE gather_to_master_int_4D( d_partial, d_tot)
    ! Gather a distributed 4-D integer variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:,:),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:,:,:,:), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_int_4D'
    INTEGER                                                            :: n1,n2,n3,n4,i,n2_proc,n3_proc,n4_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k,r
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)
    n4 = SIZE( d_partial,4)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n4     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n4_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
        IF (n4_proc /= n4) CALL crash('On master: n4 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n4, int_02 = i, int_03 = n4_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2, n3, n4))
    ELSE
      ALLOCATE( d_tot( 0,0,0,0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
    DO r = 1, n4
      CALL MPI_GATHERV( d_partial( :,j,k,r), n1, MPI_INTEGER, d_tot( :,j,k,r), counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_int_4D

  SUBROUTINE gather_to_master_dp_1D( d_partial, d_tot)
    ! Gather a distributed 1-D dp variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:      ), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_dp_1D'
    INTEGER                                                            :: n1,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot))
    ELSE
      ALLOCATE( d_tot( 0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    CALL MPI_GATHERV( d_partial, n1, MPI_DOUBLE_PRECISION, d_tot, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_dp_1D

  SUBROUTINE gather_to_master_dp_2D( d_partial, d_tot)
    ! Gather a distributed 2-D dp variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_dp_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2))
    ELSE
      ALLOCATE( d_tot( 0,0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
      CALL MPI_GATHERV( d_partial( :,j), n1, MPI_DOUBLE_PRECISION, d_tot( :,j), counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_dp_2D

  SUBROUTINE gather_to_master_dp_3D( d_partial, d_tot)
    ! Gather a distributed 3-D dp variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_dp_3D'
    INTEGER                                                            :: n1,n2,n3,i,n2_proc,n3_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2, n3))
    ELSE
      ALLOCATE( d_tot( 0,0,0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
      CALL MPI_GATHERV( d_partial( :,j,k), n1, MPI_DOUBLE_PRECISION, d_tot( :,j,k), counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_dp_3D

  SUBROUTINE gather_to_master_dp_4D( d_partial, d_tot)
    ! Gather a distributed 4-D dp variable to the master
    ! (e.g. for writing to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:,:),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE,           INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_dp_4D'
    INTEGER                                                            :: n1,n2,n3,n4,i,n2_proc,n3_proc,n4_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k,r
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)
    n4 = SIZE( d_partial,4)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n4     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n4_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
        IF (n4_proc /= n4) CALL crash('On master: n4 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n4, int_02 = i, int_03 = n4_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_GATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety - check if the distribution matches expectations
    CALL partition_list( n_tot, par%i, par%n, i1, i2)
    IF (i2 + 1 - i1 /= n1) CALL crash('process {int_01} should have {int_02} values but has {int_03}!', &
      int_01 = par%i, int_02 = i2 + 1 - i1, int_03 = n1)

    ! Allocate memory for the gathered data only on the master
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2, n3, n4))
    ELSE
      ALLOCATE( d_tot( 0,0,0,0))
    END IF

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
    DO r = 1, n4
      CALL MPI_GATHERV( d_partial( :,j,k,r), n1, MPI_DOUBLE_PRECISION, d_tot( :,j,k,r), counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_dp_4D

! ===== Distribute variables from the Master =====
! ================================================

  SUBROUTINE distribute_from_master_int_1D( d_tot, d_partial)
    ! Distribute a 1-D integer variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:      ),                        INTENT(IN)    :: d_tot

    ! Output variables:
    INTEGER , DIMENSION(:      ), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_int_1D'
    INTEGER                                                            :: n1,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1))

    ! Distribute data from the master
    CALL MPI_SCATTERV( d_tot, counts, displs, MPI_INTEGER, d_partial, n1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_int_1D

  SUBROUTINE distribute_from_master_int_2D( d_tot, d_partial)
    ! Distribute a 2-D integer variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:    ),                        INTENT(IN)    :: d_tot

    ! Output variables:
    INTEGER , DIMENSION(:,:    ), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_int_2D'
    INTEGER                                                            :: n1,n2,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    n2    = SIZE( d_tot,2)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n2   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1,n2))

    ! Distribute data from the master
    DO j = 1, n2
      CALL MPI_SCATTERV( d_tot( :,j), counts, displs, MPI_INTEGER, d_partial( :,j), n1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_int_2D

  SUBROUTINE distribute_from_master_int_3D( d_tot, d_partial)
    ! Distribute a 3-D integer variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_tot

    ! Output variables:
    INTEGER , DIMENSION(:,:,:  ), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_int_3D'
    INTEGER                                                            :: n1,n2,n3,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    n2    = SIZE( d_tot,2)
    n3    = SIZE( d_tot,3)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n2   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n3   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1,n2,n3))

    ! Distribute data from the master
    DO j = 1, n2
    DO k = 1, n3
      CALL MPI_SCATTERV( d_tot( :,j,k), counts, displs, MPI_INTEGER, d_partial( :,j,k), n1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_int_3D

  SUBROUTINE distribute_from_master_int_4D( d_tot, d_partial)
    ! Distribute a 4-D integer variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:,:),                        INTENT(IN)    :: d_tot

    ! Output variables:
    INTEGER , DIMENSION(:,:,:,:), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_int_4D'
    INTEGER                                                            :: n1,n2,n3,n4,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k,r
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    n2    = SIZE( d_tot,2)
    n3    = SIZE( d_tot,3)
    n4    = SIZE( d_tot,3)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n2   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n3   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n4   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1,n2,n3,n4))

    ! Distribute data from the master
    DO j = 1, n2
    DO k = 1, n3
    DO r = 1, n4
      CALL MPI_SCATTERV( d_tot( :,j,k,r), counts, displs, MPI_INTEGER, d_partial( :,j,k,r), n1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_int_4D

  SUBROUTINE distribute_from_master_dp_1D( d_tot, d_partial)
    ! Distribute a 1-D dp variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      ),                        INTENT(IN)    :: d_tot

    ! Output variables:
    REAL(dp), DIMENSION(:      ), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_dp_1D'
    INTEGER                                                            :: n1,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1))

    ! Distribute data from the master
    CALL MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_PRECISION, d_partial, n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_dp_1D

  SUBROUTINE distribute_from_master_dp_2D( d_tot, d_partial)
    ! Distribute a 2-D dp variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    ),                        INTENT(IN)    :: d_tot

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_dp_2D'
    INTEGER                                                            :: n1,n2,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    n2    = SIZE( d_tot,2)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n2   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1,n2))

    ! Distribute data from the master
    DO j = 1, n2
      CALL MPI_SCATTERV( d_tot( :,j), counts, displs, MPI_DOUBLE_PRECISION, d_partial( :,j), n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_dp_2D

  SUBROUTINE distribute_from_master_dp_3D( d_tot, d_partial)
    ! Distribute a 3-D dp variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_tot

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_dp_3D'
    INTEGER                                                            :: n1,n2,n3,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    n2    = SIZE( d_tot,2)
    n3    = SIZE( d_tot,3)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n2   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n3   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1,n2,n3))

    ! Distribute data from the master
    DO j = 1, n2
    DO k = 1, n3
      CALL MPI_SCATTERV( d_tot( :,j,k), counts, displs, MPI_DOUBLE_PRECISION, d_partial( :,j,k), n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_dp_3D

  SUBROUTINE distribute_from_master_dp_4D( d_tot, d_partial)
    ! Distribute a 4-D dp variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:,:),                        INTENT(IN)    :: d_tot

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE,           INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_dp_4D'
    INTEGER                                                            :: n1,n2,n3,n4,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k,r
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total size of array
    n_tot = SIZE( d_tot,1)
    n2    = SIZE( d_tot,2)
    n3    = SIZE( d_tot,3)
    n4    = SIZE( d_tot,3)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n2   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n3   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( n4   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Determine array sizes to be received by each process
    DO i = 0, par%n-1
      CALL partition_list( n_tot, i, par%n, i1, i2)
      counts( i+1) = i2 + 1 - i1
      IF (i == par%i) n1 = i2 + 1 - i1
    END DO

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Allocate memory for the distributed data
    ALLOCATE( d_partial( n1,n2,n3,n4))

    ! Distribute data from the master
    DO j = 1, n2
    DO k = 1, n3
    DO r = 1, n4
      CALL MPI_SCATTERV( d_tot( :,j,k,r), counts, displs, MPI_DOUBLE_PRECISION, d_partial( :,j,k,r), n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    END DO
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_dp_4D

! ===== Gather distributed variables to all processes =====
! =========================================================

  ! Versions where each process has a total array partially filled
  ! (i.e. recv_buf is the same as send_buf)

  SUBROUTINE allgather_array_dp( array, i1_, i2_)
    ! Gather a distributed 1-D dp variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      )                       , INTENT(INOUT) :: array
    INTEGER                                  , OPTIONAL, INTENT(IN)    :: i1_,i2_

    ! Local variables:
    INTEGER                                                            :: i1,i2, n
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Gather sizes that will be sent
    IF (PRESENT(i1_) .AND. PRESENT(i2_)) THEN
      i1 = i1_
      i2 = i2_
    ELSE
      CALL partition_list( SIZE( array), par%i, par%n, i1, i2)
    END IF
    CALL MPI_ALLGATHER( i2-i1+1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Calculate offsets through the sizes
    displs( 1) = 0
    DO n = 2, SIZE( displs)
      displs( n) = displs( n-1) + counts( n-1)
    END DO

    ! Safety
    IF (SUM( counts) /= SIZE( array)) THEN
      CALL crash("sizes dont match")
    END IF

    ! Gather data to all processes
    CALL MPI_ALLGATHERV( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, array, counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr)

  END SUBROUTINE allgather_array_dp

  SUBROUTINE allgather_array_dp_2d(array, i1_, i2_)
    ! Gather a distributed 1-D dp variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    )                       , INTENT(INOUT) :: array
    INTEGER                                  , OPTIONAL, INTENT(IN)    :: i1_,i2_

    ! Local variables:
    INTEGER                                                            :: i1,i2, n
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    IF (PRESENT(i1_) .AND. PRESENT(i2_)) THEN
      DO n = 1, size(array,2)
        CALL allgather_array_dp ( array(:,n), i1_, i2_ )
      END DO
    ELSE
      DO n = 1, size(array,2)
        CALL allgather_array_dp ( array(:,n) )
      END DO
    END IF

  END SUBROUTINE allgather_array_dp_2d

  SUBROUTINE allgather_array_int( array,i1_,i2_)
    ! Gather a distributed 1-D dp variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:      )                       , INTENT(INOUT) :: array
    INTEGER                                  , OPTIONAL, INTENT(IN)    :: i1_,i2_

    ! Local variables:
    INTEGER                                                            :: i1,i2, n
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Gather sizes that will be sent
    IF (PRESENT(i1_) .AND. PRESENT(i2_)) THEN
      i1 = i1_
      i2 = i2_
    ELSE
      CALL partition_list( SIZE( array), par%i, par%n, i1, i2)
    END IF
    CALL MPI_ALLGATHER( i2-i1+1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Calculate offsets through the sizes
    displs( 1) = 0
    DO n = 2, SIZE( displs)
      displs( n) = displs( n-1) + counts( n-1)
    END DO

    ! Safety
    IF (SUM( counts) /= SIZE( array)) THEN
      CALL crash("sizes dont match")
    END IF

    ! Gather data to all processes
    CALL MPI_ALLGATHERV( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, array, counts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  END SUBROUTINE allgather_array_int

  ! Versions where each process has a partial array totally filled
  ! (i.e. recv_buf is different from send_buf)

  SUBROUTINE gather_to_all_int_1D( d_partial, d_tot)
    ! Gather a distributed 1-D integer variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:      ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:      ),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_int_1D'
    INTEGER                                                            :: n1,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to all processes
    CALL MPI_ALLGATHERV( d_partial, n1, MPI_INTEGER, d_tot, counts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_int_1D

  SUBROUTINE gather_to_all_int_2D( d_partial, d_tot)
    ! Gather a distributed 2-D integer variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:    ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:,:    ),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_int_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
      CALL MPI_ALLGATHERV( d_partial( :,j), n1, MPI_INTEGER, d_tot( :,j), counts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_int_2D

  SUBROUTINE gather_to_all_int_3D( d_partial, d_tot)
    ! Gather a distributed 3-D integer variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:,:,:  ),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_int_3D'
    INTEGER                                                            :: n1,n2,n3,i,n2_proc,n3_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
      CALL MPI_ALLGATHERV( d_partial( :,j,k), n1, MPI_INTEGER, d_tot( :,j,k), counts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_int_3D

  SUBROUTINE gather_to_all_int_4D( d_partial, d_tot)
    ! Gather a distributed 4-D integer variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:,:),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:,:,:,:),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_int_4D'
    INTEGER                                                            :: n1,n2,n3,n4,i,n2_proc,n3_proc,n4_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k,q
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)
    n4 = SIZE( d_partial,4)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n4     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n4_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
        IF (n4_proc /= n4) CALL crash('On master: n4 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n4, int_02 = i, int_03 = n4_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
    DO q = 1, n4
      CALL MPI_ALLGATHERV( d_partial( :,j,k,q), n1, MPI_INTEGER, d_tot( :,j,k,q), counts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    END DO
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_int_4D

  SUBROUTINE gather_to_all_dp_1D( d_partial, d_tot)
    ! Gather a distributed 1-D dp variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:      ),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_dp_1D'
    INTEGER                                                            :: n1,i
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to all processes
    CALL MPI_ALLGATHERV( d_partial, n1, MPI_DOUBLE_PRECISION, d_tot, counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_dp_1D

  SUBROUTINE gather_to_all_dp_2D( d_partial, d_tot)
    ! Gather a distributed 2-D dp variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_dp_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
      CALL MPI_ALLGATHERV( d_partial( :,j), n1, MPI_DOUBLE_PRECISION, d_tot( :,j), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_dp_2D

  SUBROUTINE gather_to_all_dp_3D( d_partial, d_tot)
    ! Gather a distributed 3-D dp variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_dp_3D'
    INTEGER                                                            :: n1,n2,n3,i,n2_proc,n3_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
      CALL MPI_ALLGATHERV( d_partial( :,j,k), n1, MPI_DOUBLE_PRECISION, d_tot( :,j,k), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_dp_3D

  SUBROUTINE gather_to_all_dp_4D( d_partial, d_tot)
    ! Gather a distributed 4-D dp variable to all processes

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:,:),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:,:),                        INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_all_dp_4D'
    INTEGER                                                            :: n1,n2,n3,n4,i,n2_proc,n3_proc,n4_proc
    INTEGER                                                            :: i1,i2
    INTEGER                                                            :: n_tot, j,k,q
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)
    n3 = SIZE( d_partial,3)
    n4 = SIZE( d_partial,4)

    ! Safety - check if the secondary dimensions are the same everywhere
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n3     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
        CALL MPI_SEND( n4     , 1, MPI_INTEGER, 0, 0          , MPI_COMM_WORLD            , ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( n4_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('On master: n2 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n2, int_02 = i, int_03 = n2_proc)
        IF (n3_proc /= n3) CALL crash('On master: n3 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n3, int_02 = i, int_03 = n3_proc)
        IF (n4_proc /= n4) CALL crash('On master: n4 == {int_01}, but on process {int_02} its {int_03}!',&
          int_01 =  n4, int_02 = i, int_03 = n4_proc)
      END IF
    END DO

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    DO j = 1, n2
    DO k = 1, n3
    DO q = 1, n4
      CALL MPI_ALLGATHERV( d_partial( :,j,k,q), n1, MPI_DOUBLE_PRECISION, d_tot( :,j,k,q), counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    END DO
    END DO
    END DO

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_dp_4D

! ===== Unit tests =====
! ======================

  SUBROUTINE run_all_mpi_distributed_memory_unit_tests
    ! Run all unit tests for the MPI distributed memory subroutines

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_all_mpi_distributed_memory_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all unit tests for the MPI distributed memory subroutines
    CALL test_gather_to_all

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_mpi_distributed_memory_unit_tests

  SUBROUTINE test_gather_to_all
    ! Test the gather_to_all_TYPE_DIM subroutines

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_gather_to_all_int_1D'
    INTEGER,  DIMENSION(:      ), ALLOCATABLE                          :: int_1D_01, int_1D_02, int_1D_03
    INTEGER,  DIMENSION(:,:    ), ALLOCATABLE                          :: int_2D_01, int_2D_02, int_2D_03
    INTEGER,  DIMENSION(:,:,:  ), ALLOCATABLE                          :: int_3D_01, int_3D_02, int_3D_03
    INTEGER,  DIMENSION(:,:,:,:), ALLOCATABLE                          :: int_4D_01, int_4D_02, int_4D_03
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: dp_1D_01, dp_1D_02, dp_1D_03
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: dp_2D_01, dp_2D_02, dp_2D_03
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE                          :: dp_3D_01, dp_3D_02, dp_3D_03
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE                          :: dp_4D_01, dp_4D_02, dp_4D_03
    LOGICAL                                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - should be run on at least two cores
    IF (par%n < 2) CALL crash('should be run on at least two cores')
    IF (par%i > 1) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    found_errors = .FALSE.

  ! == int_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_1D_01( 2), source = 0)
      ALLOCATE( int_1D_02( 7), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_1D_01( 5), source = 0)
      ALLOCATE( int_1D_02( 7), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    END IF

    ! Fill test data
    int_1D_03 = [1, 2, 3, 4, 5, 6, 7]

    IF     (par%i == 0) THEN
      int_1D_01 = int_1D_03( 1:2)
    ELSEIF (par%i == 1) THEN
      int_1D_01 = int_1D_03( 3:7)
    END IF

    ! Gather data
    CALL gather_to_all_int_1D( int_1D_01, int_1D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( int_1D_02 /= int_1D_03)

    ! Clean up after yourself
    DEALLOCATE( int_1D_01)
    DEALLOCATE( int_1D_02)
    DEALLOCATE( int_1D_03)

  ! == int_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_2D_01( 2,2), source = 0)
      ALLOCATE( int_2D_02( 7,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_2D_01( 5,2), source = 0)
      ALLOCATE( int_2D_02( 7,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    END IF

    ! Fill test data
    int_2D_03( :,1) = [ 1,  2,  3,  4,  5,  6,  7]
    int_2D_03( :,2) = [ 8,  9, 10, 11, 12, 13, 14]

    IF     (par%i == 0) THEN
      int_2D_01 = int_2D_03( 1:2,:)
    ELSEIF (par%i == 1) THEN
      int_2D_01 = int_2D_03( 3:7,:)
    END IF

    ! Gather data
    CALL gather_to_all_int_2D( int_2D_01, int_2D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( int_2D_02 /= int_2D_03)

    ! Clean up after yourself
    DEALLOCATE( int_2D_01)
    DEALLOCATE( int_2D_02)
    DEALLOCATE( int_2D_03)

  ! == int_3D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_3D_01( 2,2,2), source = 0)
      ALLOCATE( int_3D_02( 7,2,2), source = 0)
      ALLOCATE( int_3D_03( 7,2,2), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_3D_01( 5,2,2), source = 0)
      ALLOCATE( int_3D_02( 7,2,2), source = 0)
      ALLOCATE( int_3D_03( 7,2,2), source = 0)
    END IF

    ! Fill test data
    int_3D_03( :,1,1) = [ 1,  2,  3,  4,  5,  6,  7]
    int_3D_03( :,2,1) = [ 8,  9, 10, 11, 12, 13, 14]
    int_3D_03( :,1,2) = [15, 16, 17, 18, 19, 20, 21]
    int_3D_03( :,2,2) = [22, 23, 24, 25, 26, 27, 28]

    IF     (par%i == 0) THEN
      int_3D_01 = int_3D_03( 1:2,:,:)
    ELSEIF (par%i == 1) THEN
      int_3D_01 = int_3D_03( 3:7,:,:)
    END IF

    ! Gather data
    CALL gather_to_all_int_3D( int_3D_01, int_3D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( int_3D_02 /= int_3D_03)

    ! Clean up after yourself
    DEALLOCATE( int_3D_01)
    DEALLOCATE( int_3D_02)
    DEALLOCATE( int_3D_03)

  ! == int_4D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_4D_01( 2,2,2,2), source = 0)
      ALLOCATE( int_4D_02( 7,2,2,2), source = 0)
      ALLOCATE( int_4D_03( 7,2,2,2), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_4D_01( 5,2,2,2), source = 0)
      ALLOCATE( int_4D_02( 7,2,2,2), source = 0)
      ALLOCATE( int_4D_03( 7,2,2,2), source = 0)
    END IF

    ! Fill test data
    int_4D_03( :,1,1,1) = [ 1,  2,  3,  4,  5,  6,  7]
    int_4D_03( :,2,1,1) = [ 8,  9, 10, 11, 12, 13, 14]
    int_4D_03( :,1,2,1) = [15, 16, 17, 18, 19, 20, 21]
    int_4D_03( :,2,2,1) = [22, 23, 24, 25, 26, 27, 28]
    int_4D_03( :,1,1,2) = [29, 30, 31, 32, 33, 34, 35]
    int_4D_03( :,2,1,2) = [36, 37, 38, 39, 40, 41, 42]
    int_4D_03( :,1,2,2) = [43, 44, 45, 46, 47, 48, 49]
    int_4D_03( :,2,2,2) = [50, 51, 52, 53, 54, 55, 56]

    IF     (par%i == 0) THEN
      int_4D_01 = int_4D_03( 1:2,:,:,:)
    ELSEIF (par%i == 1) THEN
      int_4D_01 = int_4D_03( 3:7,:,:,:)
    END IF

    ! Gather data
    CALL gather_to_all_int_4D( int_4D_01, int_4D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( int_4D_02 /= int_4D_03)

    ! Clean up after yourself
    DEALLOCATE( int_4D_01)
    DEALLOCATE( int_4D_02)
    DEALLOCATE( int_4D_03)

  ! == dp_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_1D_01( 2), source = 0._dp)
      ALLOCATE( dp_1D_02( 7), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_1D_01( 5), source = 0._dp)
      ALLOCATE( dp_1D_02( 7), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    END IF

    ! Fill test data
    dp_1D_03 = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp]

    IF     (par%i == 0) THEN
      dp_1D_01 = dp_1D_03( 1:2)
    ELSEIF (par%i == 1) THEN
      dp_1D_01 = dp_1D_03( 3:7)
    END IF

    ! Gather data
    CALL gather_to_all_dp_1D( dp_1D_01, dp_1D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( dp_1D_02 /= dp_1D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_1D_01)
    DEALLOCATE( dp_1D_02)
    DEALLOCATE( dp_1D_03)

  ! == dp_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_2D_01( 2,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_2D_01( 5,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    END IF

    ! Fill test data
    dp_2D_03( :,1) = [ 1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    dp_2D_03( :,2) = [ 8._dp,  9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]

    IF     (par%i == 0) THEN
      dp_2D_01 = dp_2D_03( 1:2,:)
    ELSEIF (par%i == 1) THEN
      dp_2D_01 = dp_2D_03( 3:7,:)
    END IF

    ! Gather data
    CALL gather_to_all_dp_2D( dp_2D_01, dp_2D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( dp_2D_02 /= dp_2D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_2D_01)
    DEALLOCATE( dp_2D_02)
    DEALLOCATE( dp_2D_03)

  ! == dp_3D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_3D_01( 2,2,2), source = 0._dp)
      ALLOCATE( dp_3D_02( 7,2,2), source = 0._dp)
      ALLOCATE( dp_3D_03( 7,2,2), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_3D_01( 5,2,2), source = 0._dp)
      ALLOCATE( dp_3D_02( 7,2,2), source = 0._dp)
      ALLOCATE( dp_3D_03( 7,2,2), source = 0._dp)
    END IF

    ! Fill test data
    dp_3D_03( :,1,1) = [ 1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    dp_3D_03( :,2,1) = [ 8._dp,  9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]
    dp_3D_03( :,1,2) = [15._dp, 16._dp, 17._dp, 18._dp, 19._dp, 20._dp, 21._dp]
    dp_3D_03( :,2,2) = [22._dp, 23._dp, 24._dp, 25._dp, 26._dp, 27._dp, 28._dp]

    IF     (par%i == 0) THEN
      dp_3D_01 = dp_3D_03( 1:2,:,:)
    ELSEIF (par%i == 1) THEN
      dp_3D_01 = dp_3D_03( 3:7,:,:)
    END IF

    ! Gather data
    CALL gather_to_all_dp_3D( dp_3D_01, dp_3D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( dp_3D_02 /= dp_3D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_3D_01)
    DEALLOCATE( dp_3D_02)
    DEALLOCATE( dp_3D_03)

  ! == dp_4D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_4D_01( 2,2,2,2), source = 0._dp)
      ALLOCATE( dp_4D_02( 7,2,2,2), source = 0._dp)
      ALLOCATE( dp_4D_03( 7,2,2,2), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_4D_01( 5,2,2,2), source = 0._dp)
      ALLOCATE( dp_4D_02( 7,2,2,2), source = 0._dp)
      ALLOCATE( dp_4D_03( 7,2,2,2), source = 0._dp)
    END IF

    ! Fill test data
    dp_4D_03( :,1,1,1) = [ 1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    dp_4D_03( :,2,1,1) = [ 8._dp,  9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]
    dp_4D_03( :,1,2,1) = [15._dp, 16._dp, 17._dp, 18._dp, 19._dp, 20._dp, 21._dp]
    dp_4D_03( :,2,2,1) = [22._dp, 23._dp, 24._dp, 25._dp, 26._dp, 27._dp, 28._dp]
    dp_4D_03( :,1,1,2) = [29._dp, 30._dp, 31._dp, 32._dp, 33._dp, 34._dp, 35._dp]
    dp_4D_03( :,2,1,2) = [36._dp, 37._dp, 38._dp, 39._dp, 40._dp, 41._dp, 42._dp]
    dp_4D_03( :,1,2,2) = [43._dp, 44._dp, 45._dp, 46._dp, 47._dp, 48._dp, 49._dp]
    dp_4D_03( :,2,2,2) = [50._dp, 51._dp, 52._dp, 53._dp, 54._dp, 55._dp, 56._dp]

    IF     (par%i == 0) THEN
      dp_4D_01 = dp_4D_03( 1:2,:,:,:)
    ELSEIF (par%i == 1) THEN
      dp_4D_01 = dp_4D_03( 3:7,:,:,:)
    END IF

    ! Gather data
    CALL gather_to_all_dp_4D( dp_4D_01, dp_4D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( dp_4D_02 /= dp_4D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_4D_01)
    DEALLOCATE( dp_4D_02)
    DEALLOCATE( dp_4D_03)

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all gather_to_all routines')
    ELSE
      IF (par%master) CALL warning('found errors in gather_to_all routines')
    END IF

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_gather_to_all

END MODULE mpi_distributed_memory
