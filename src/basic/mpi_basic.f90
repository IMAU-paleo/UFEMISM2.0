module mpi_basic

  ! Some very basic stuff to support the MPI parallelised architecture.

  use mpi

  implicit none

  private

  public :: par, sync, initialise_parallelisation

  type parallel_info

    integer :: n        ! Total number of processes (0 = single-core, >0 = master+slaves)
    integer :: i        ! Global ID of this process (0 = master, >0 = slave)
    logical :: master   ! Whether or not the current process is the master process

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
    par%master = (par%i == 0)

  end subroutine initialise_parallelisation

  subroutine sync

    ! Local variables:
    integer :: ierr

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

  end subroutine sync

end module mpi_basic
