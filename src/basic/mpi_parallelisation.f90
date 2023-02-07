MODULE mpi_parallelisation

  ! A collection of different routines that make parallel programming in UFEMISM a lot easier.

! ===== Preamble =====
! ====================

  USE mpi

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  INTEGER :: cerr, ierr    ! Error flags for MPI routines

  TYPE parallel_info

    INTEGER                       :: i        ! ID of this process (0 = master, >0 = slave)
    INTEGER                       :: n        ! Total number of processes (0 = single-core, >0 = master+slaves)
    LOGICAL                       :: master   ! Whether or not the current process is the master process

  END TYPE parallel_info

  TYPE(parallel_info), SAVE :: par

CONTAINS

! ===== Subroutinea ======
! ========================

  ! Initialise the MPI parallelisation
  SUBROUTINE initialise_parallelisation

    IMPLICIT NONE

    ! MPI Initialisation
    ! ==================

    ! Use MPI to create copies of the program on all the processors, so the model can run in parallel.
    CALL MPI_INIT(ierr)

    ! Get rank of current process and total number of processes
    CALL MPI_COMM_RANK(       MPI_COMM_WORLD, par%i, ierr)
    CALL MPI_COMM_SIZE(       MPI_COMM_WORLD, par%n, ierr)
    par%master = (par%i == 0)

  END SUBROUTINE initialise_parallelisation

  ! Synchronise the different processes
  SUBROUTINE sync
    ! Use MPI_BARRIER to synchronise all the processes

    IMPLICIT NONE

    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

  END SUBROUTINE sync

  ! Partition a list of ntot elements over the n processes
  SUBROUTINE partition_list( ntot, i, n, i1, i2)
    ! Partition a list into parallel ranges (e.g. vertex domains)

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


END MODULE mpi_parallelisation
