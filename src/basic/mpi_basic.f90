MODULE mpi_basic

  ! Some very basic stuff to support the MPI parallelised architecture.

! ===== Preamble =====
! ====================

  USE mpi

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  INTEGER :: cerr, ierr                    ! Error flags for MPI routines
  INTEGER :: MPI_status( MPI_STATUS_SIZE)  ! Status flag for MPI_RECV

  TYPE parallel_info

    INTEGER                       :: i        ! ID of this process (0 = master, >0 = slave)
    INTEGER                       :: n        ! Total number of processes (0 = single-core, >0 = master+slaves)
    LOGICAL                       :: master   ! Whether or not the current process is the master process

  END TYPE parallel_info

  TYPE(parallel_info), SAVE :: par

CONTAINS

! ===== Very basic stuff ======
! =============================

  ! Initialise the MPI parallelisation
  SUBROUTINE initialise_parallelisation

    IMPLICIT NONE

    ! Use MPI to create copies of the program on all the processors, so the model can run in parallel.
    CALL MPI_INIT( ierr)

    ! Get rank of current process and total number of processes
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, par%i, ierr)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, par%n, ierr)
    par%master = (par%i == 0)

  END SUBROUTINE initialise_parallelisation

  ! Finalise the MPI parallelisation
  SUBROUTINE finalise_parallelisation

    IMPLICIT NONE

    CALL MPI_FINALIZE( ierr)

  END SUBROUTINE finalise_parallelisation

  ! Synchronise the different processes
  SUBROUTINE sync
    ! Use MPI_BARRIER to synchronise all the processes

    IMPLICIT NONE

    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

  END SUBROUTINE sync

END MODULE mpi_basic
