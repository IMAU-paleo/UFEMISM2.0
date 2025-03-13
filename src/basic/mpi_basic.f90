module mpi_basic

  ! Some very basic stuff to support the MPI parallelised architecture.

  use mpi_f08, only: MPI_COMM, MPI_INIT, MPI_COMM_SIZE, MPI_COMM_RANK, MPI_COMM_WORLD, &
    MPI_COMM_SPLIT_TYPE, MPI_COMM_TYPE_SHARED, MPI_BARRIER, MPI_INFO_NULL

  implicit none

  private

  public :: par, sync, initialise_parallelisation

  type parallel_info

    integer :: n        ! Global number of processes
    integer :: i        ! Global ID of this process
    logical :: primary  ! Whether or not the current process is the primary process (shorthand for par%i == 0)

  end type parallel_info

  type(parallel_info), save :: par

contains

  subroutine initialise_parallelisation

    ! Local variables:
    integer :: ierr

    ! Use MPI to create copies of the program on all the processors, so the model can run in parallel.
    call MPI_INIT( ierr)

    ! Get rank of current process and total number of processes
    call MPI_COMM_RANK( MPI_COMM_WORLD, par%i, ierr)
    call MPI_COMM_SIZE( MPI_COMM_WORLD, par%n, ierr)
    par%primary = (par%i == 0)

  end subroutine initialise_parallelisation

  subroutine sync

    ! Local variables:
    integer :: ierr

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

  end subroutine sync

end module mpi_basic
