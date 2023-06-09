MODULE mpi_distributed_memory

  ! Some routine to work with distributed memory in the MPI parallelised architecture.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string

  IMPLICIT NONE

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
    ! Gather a distributed 1-D integer variable to the Master

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:      ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:      ), optional,              INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_int_1D'
    INTEGER                                                            :: n1,i
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
    IF (par%master) then
      if (present(d_tot)) then
        if( n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')
      else
        CALL crash('d_tot must be present on master process')
      endif
    endif

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    CALL MPI_GATHERV( d_partial, n1, MPI_INTEGER, d_tot, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_int_1D

  SUBROUTINE gather_to_master_int_2D( d_partial, d_tot)
    ! Gather a distributed 2-D integer variable to the Master

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:    ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    INTEGER , DIMENSION(:,:    ), optional,              INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_int_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: j
    integer                                                            :: dummy(1)

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Check sizes
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('n2 = {int_01} on master, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      END IF
    END DO

    IF (par%master) then
      if (.not. present(d_tot)) CALL crash('d_tot must be present on master process')
    endif


    DO j = 1, n2
      if (par%master) then
        CALL gather_to_master_int_1D( d_partial(:,j), d_tot( :, j))
      else
        CALL gather_to_master_int_1D( d_partial(:,j), dummy)
      end if
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_int_2D

  SUBROUTINE gather_to_master_dp_1D( d_partial, d_tot)
    ! Gather a distributed 1-D dp variable to the Master

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:      ), optional,              INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_dp_1D'
    INTEGER                                                            :: n1,i
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
    IF (par%master) then
      if (present(d_tot)) then
        if( n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')
      else
        CALL crash('d_tot must be present on master process')
      endif
    endif

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Gather data to the master
    CALL MPI_GATHERV( d_partial, n1, MPI_DOUBLE_PRECISION, d_tot, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_dp_1D

  SUBROUTINE gather_to_master_dp_2D( d_partial, d_tot)
    ! Gather a distributed 2-D dp variable to the Master

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    ),                        INTENT(IN)    :: d_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ), optional,              INTENT(OUT)   :: d_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_to_master_dp_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: j
    real(dp)                                                           :: dummy(1)

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Check sizes
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('n2 = {int_01} on master, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      END IF
    END DO

    IF (par%master) then
      if (.not. present(d_tot)) CALL crash('d_tot must be present on master process')
    endif

    DO j = 1, n2
      if (par%master) then
        CALL gather_to_master_dp_1D( d_partial(:,j), d_tot( :, j))
      else
        CALL gather_to_master_dp_1D( d_partial(:,j), dummy)
      end if
    END DO


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_master_dp_2D

! ===== Gather distributed variables to all processes =====
! =========================================================

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

    ! Finalise routine path
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
    INTEGER , DIMENSION(:      ), ALLOCATABLE                          :: d_partial_1D
    INTEGER , DIMENSION(:      ), ALLOCATABLE                          :: d_tot_1D
    INTEGER                                                            :: n1_tot,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Check sizes
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('n2 = {int_01} on master, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      END IF
    END DO

    ! Gather 1 column at a time
    ALLOCATE( d_partial_1D( n1))
    n1_tot = SIZE( d_tot,1)
    ALLOCATE( d_tot_1D( n1_tot))

    DO j = 1, n2
      d_partial_1D = d_partial( :,j)
      CALL gather_to_all_int_1D( d_partial_1D, d_tot_1D)
      d_tot( :,j) = d_tot_1D
    END DO

    ! Clean up after yourself
    DEALLOCATE( d_partial_1D)
    DEALLOCATE( d_tot_1D)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_int_2D

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

    ! Finalise routine path
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
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_partial_1D
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_tot_1D
    INTEGER                                                            :: n1_tot,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Check sizes
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('n2 = {int_01} on master, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      END IF
    END DO

    ! Gather 1 column at a time
    ALLOCATE( d_partial_1D( n1))
    n1_tot = SIZE( d_tot,1)
    ALLOCATE( d_tot_1D( n1_tot))

    DO j = 1, n2
      d_partial_1D = d_partial( :,j)
      CALL gather_to_all_dp_1D( d_partial_1D, d_tot_1D)
      d_tot( :,j) = d_tot_1D
    END DO

    ! Clean up after yourself
    DEALLOCATE( d_partial_1D)
    DEALLOCATE( d_tot_1D)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_to_all_dp_2D

! ===== Distribute variables from the Master =====
! ================================================

  SUBROUTINE distribute_from_master_int_1D( d_tot, d_partial)
    ! Distribute a 1-D integer variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:      ), optional,              INTENT(IN)    :: d_tot

    ! Output variables:
    INTEGER , DIMENSION(:      ),                        INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_int_1D'
    INTEGER                                                            :: n1,n_tot,i
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the partial array owned by this process
    n1 = SIZE( d_partial,1)

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (par%master) then
      if( .not. present( d_tot)) CALL crash('d_tot must be present on master')
      if( n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')
    end if

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Scatter data to all the processes
    CALL MPI_SCATTERV( d_tot, counts, displs, MPI_INTEGER, d_partial, n1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_int_1D

  SUBROUTINE distribute_from_master_int_2D( d_tot, d_partial)
    ! Distribute a 2-D integer variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:    ), optional,              INTENT(IN)    :: d_tot

    ! Output variables:
    INTEGER , DIMENSION(:,:    ),                        INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_int_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Check sizes
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('n2 = {int_01} on master, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      END IF
    END DO

    DO j = 1, n2
      if (par%master) then
        CALL distribute_from_master_int_1D( d_tot( :, j), d_partial( : ,j))
      else
        CALL distribute_from_master_int_1D( d_partial=d_partial( : ,j))
      endif
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_int_2D

  SUBROUTINE distribute_from_master_dp_1D( d_tot, d_partial)
    ! Distribute a 1-D dp variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      ), optional,              INTENT(IN)    :: d_tot

    ! Output variables:
    REAL(dp), DIMENSION(:      ),                        INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_dp_1D'
    INTEGER                                                            :: n1,n_tot,i
    INTEGER,  DIMENSION(1:par%n)                                       :: counts, displs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the partial array owned by this process
    n1 = SIZE( d_partial,1)

    ! Determine total size of distributed array
    CALL MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    n_tot = SUM( counts)
    CALL MPI_BCAST( n_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (par%master) then
      if( .not. present( d_tot)) CALL crash('d_tot must be present on master')
      if( n_tot /= SIZE( d_tot,1)) CALL crash('combined sizes of d_partial dont match size of d_tot')
    end if

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    DO i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    END DO

    ! Scatter data to all the processes
    CALL MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_PRECISION, d_partial, n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_dp_1D

  SUBROUTINE distribute_from_master_dp_2D( d_tot, d_partial)
    ! Distribute a 2-D dp variable from the master
    ! (e.g. after reading from to NetCDF)

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    ), optional,              INTENT(IN)    :: d_tot

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ),                        INTENT(OUT)   :: d_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_from_master_dp_2D'
    INTEGER                                                            :: n1,n2,i,n2_proc
    INTEGER                                                            :: j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = SIZE( d_partial,1)
    n2 = SIZE( d_partial,2)

    ! Check sizes
    DO i = 1, par%n-1
      IF (par%i == i) THEN
        CALL MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
      ELSEIF (par%master) THEN
        CALL MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        IF (n2_proc /= n2) CALL crash('n2 = {int_01} on master, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      END IF
    END DO

    ! Distribute 1 column at a time
    DO j = 1, n2
      if (par%master) then
        CALL distribute_from_master_dp_1D( d_tot( :, j), d_partial( : ,j))
      else
        CALL distribute_from_master_dp_1D( d_partial=d_partial( : ,j))
      endif
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_from_master_dp_2D

END MODULE mpi_distributed_memory
