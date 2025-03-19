module mpi_basic

  ! Some very basic stuff to support the MPI parallelised architecture.

  use mpi_f08, only: MPI_COMM, MPI_INIT, MPI_COMM_SIZE, MPI_COMM_RANK, MPI_COMM_WORLD, &
    MPI_COMM_SPLIT_TYPE, MPI_COMM_TYPE_SHARED, MPI_BARRIER, MPI_INFO_NULL, MPI_COMM_SPLIT, &
    MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_SEND, MPI_RECV, MPI_STATUS, &
    MPI_ANY_TAG

  implicit none

  private

  public :: par, sync, initialise_parallelisation, sync_node

  type parallel_info

    integer        :: n                       ! Global number of processes
    integer        :: i                       ! Global ID of this process
    logical        :: primary                 ! Whether or not the current process is the primary process (shorthand for par%i == 0)

    integer        :: n_nodes                 ! Total number of shared-memory nodes
    integer        :: node_ID                 ! ID of the node this process is part of
    type(MPI_COMM) :: mpi_comm_node           ! MPI communicator for all processes within this node
    integer        :: n_node                  ! Number of processes in the node that this process is part of
    integer        :: i_node                  ! ID of this process within mpi_comm_node
    logical        :: node_primary            ! Whether or not the current process is the master of this node (shorthand for par%i_node == 0)
    type(MPI_COMM) :: mpi_comm_node_primaries ! MPI communicator for all node master processes
    type(MPI_COMM) :: mpi_comm_secondaries    ! MPI communicator for all non-node-master processes (unused, but needed to complete the MPI_SPLIT_COMM call)

  end type parallel_info

  type(parallel_info), save :: par

contains

  subroutine initialise_parallelisation( UFEMISM_program_input_argument)

    ! In/output variables:
    character(len=*), intent(in) :: UFEMISM_program_input_argument

    ! Local variables:
    integer :: ierr, i, n

    ! Use MPI to create copies of the program on all the processors, so the model can run in parallel.
    call MPI_INIT( ierr)

    ! Get global number and rank of processes
    call MPI_COMM_SIZE( MPI_COMM_WORLD, par%n, ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, par%i, ierr)
    par%primary = (par%i == 0)

    if (UFEMISM_program_input_argument == 'unit_tests_multinode') then
      ! In this case, we're actually running on 7 processes on the same
      ! machine; regardless, we'll "pretend" they are on 3 separate nodes,
      ! and split the global communicator accordingly, so we can test
      ! the hybrid distributed/shared memory code.

      ! Safety
      if (.not. par%n == 7) stop 'multi-node unit tests should be run on 4 processes!'

      if (par%i == 0 .or. par%i == 1) then
        call MPI_COMM_SPLIT( MPI_COMM_WORLD, 0, par%i, par%mpi_comm_node, ierr)
      elseif (par%i == 2 .or. par%i == 3 .or. par%i == 4) then
        call MPI_COMM_SPLIT( MPI_COMM_WORLD, 1, par%i, par%mpi_comm_node, ierr)
      else
        call MPI_COMM_SPLIT( MPI_COMM_WORLD, 2, par%i, par%mpi_comm_node, ierr)
      end if

    else
      ! Split global communicator into communicators per shared-memory node
      call MPI_COMM_SPLIT_TYPE( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, par%i, MPI_INFO_NULL, par%mpi_comm_node, ierr)
    end if

    ! Get node number and rank of processes
    call MPI_COMM_SIZE( par%mpi_comm_node, par%n_node, ierr)
    call MPI_COMM_RANK( par%mpi_comm_node, par%i_node, ierr)
    par%node_primary = (par%i_node == 0)

    ! Determine number of nodes
    call determine_number_of_nodes_and_node_IDs

    ! Create communicator for primaries only
    if (par%node_primary) then
      call MPI_COMM_SPLIT( MPI_COMM_WORLD, 0, par%i, par%mpi_comm_node_primaries, ierr)

      ! Safety
      call MPI_COMM_SIZE( par%mpi_comm_node_primaries, n, ierr)
      if (n /= par%n_nodes) stop 'number of node primaries should be equal to n_nodes!'
      call MPI_COMM_RANK( par%mpi_comm_node_primaries, i, ierr)
      if (i /= par%node_ID) stop 'rank in mpi_comm_node_primaries should be equal to node_ID!'
    else
      call MPI_COMM_SPLIT( MPI_COMM_WORLD, 1, par%i, par%mpi_comm_secondaries, ierr)
    end if

    if (UFEMISM_program_input_argument == 'unit_tests_multinode') call print_parallelisation_info

  end subroutine initialise_parallelisation

  subroutine determine_number_of_nodes_and_node_IDs

    ! Let each process, in order, tell the primary what their rank within
    ! their node communicator is. Since they were split using their global rank as the ley,
    ! they are in the same order, so the primary will receive numbers like 0-1-2-3-0-1-2-3-0-1-2-3.
    ! So, the number of nodes is equal to the number of 0's (or more simply, the number of node primaries),
    ! and the node ID's can be determined by increasing the count every time a process reports
    ! a rank of 0 within their node communicator

    ! Local variables:
    integer                       :: ierr, i
    integer, dimension(0:par%n-1) :: i_node, node_ID
    type(MPI_STATUS)              :: recv_status

    ! Number of nodes = number of processes that are a node primary
    if (par%node_primary) then
      par%n_nodes = 1
    else
      par%n_nodes = 0
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, par%n_nodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)

    if (par%primary) then

      ! Receive rank within node from all global processes
      i_node( 0) = par%i_node
      do i = 1, par%n-1
        call MPI_RECV( i_node( i), 1, MPI_INTEGER, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
      end do

      ! Determine ID of each process' node
      node_ID(0) = 0
      do i = 1, par%n-1
        if (i_node( i) == 0) then
          node_ID( i) = node_ID( i-1) + 1
        else
          node_ID( i) = node_ID( i-1)
        end if
      end do
      par%node_ID = node_ID(0)

      ! Send IDs of each process' nodes to them
      do i = 1, par%n-1
        call MPI_SEND( node_ID( i), 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
      end do

    else

      ! Send this process' rank within node to the global primary
      call MPI_SEND( par%i_node, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)

      ! Receive this process' node ID from the global primary
      call MPI_RECV( par%node_ID, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)

    end if

  end subroutine determine_number_of_nodes_and_node_IDs

  subroutine print_parallelisation_info

    ! Local variables:
    integer :: i

    if (par%primary) write(0,*) ''
    call sync

    do i = 0, par%n-1
      if (i == par%i) then
        write(0,'(&
          A,I1,A,I1,&
          A,L1,&
          A,I1,A,I1,&
          A,L1,&
          A,I1,A,I1)') &
          ' Process ', par%i, '/', par%n, &
          ': primary = ', par%primary, &
          ', process ', par%i_node, '/', par%n_node, &
          ' (node primary = ', par%node_primary, &
          ') on node ', par%node_ID, '/', par%n_nodes
      end if
      call sync
    end do

  end subroutine print_parallelisation_info

  subroutine sync

    ! Local variables:
    integer :: ierr

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

  end subroutine sync

  subroutine sync_node

    ! Local variables:
    integer :: ierr

    call MPI_BARRIER( par%mpi_comm_node, ierr)

  end subroutine sync_node

end module mpi_basic
