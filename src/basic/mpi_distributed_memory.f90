MODULE mpi_distributed_memory

  ! Some routine to work with distributed memory in the MPI parallelised architecture.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  ! Partition a list of ntot elements over the n processes
  SUBROUTINE partition_list( ntot, i, n, i1, i2)
    ! Partition a list of ntot elements over the n processes

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                    INTENT(IN)        :: ntot, i, n
    INTEGER,                    INTENT(OUT)       :: i1, i2

    IF (ntot > n*2) THEN
      i1 = MAX(1,    FLOOR(REAL(ntot *  i      / n)) + 1)
      i2 = MIN(ntot, FLOOR(REAL(ntot * (i + 1) / n)))
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

END MODULE mpi_distributed_memory
