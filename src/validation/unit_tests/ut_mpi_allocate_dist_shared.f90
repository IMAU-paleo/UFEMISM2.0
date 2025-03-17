module ut_mpi_allocate_dist_shared

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD

  implicit none

  private

  public :: test_allocate_dist_shared

contains

  subroutine test_allocate_dist_shared( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_allocate_dist_shared'
    character(len=1024), parameter :: test_name_local = 'allocate_dist_shared'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_allocate_dist_shared_logical_1D( test_name)
    call test_allocate_dist_shared_logical_2D( test_name)
    call test_allocate_dist_shared_logical_3D( test_name)

    call test_allocate_dist_shared_int_1D( test_name)
    call test_allocate_dist_shared_int_2D( test_name)
    call test_allocate_dist_shared_int_3D( test_name)

    call test_allocate_dist_shared_dp_1D( test_name)
    call test_allocate_dist_shared_dp_2D( test_name)
    call test_allocate_dist_shared_dp_3D( test_name)

    call test_allocate_dist_shared_complex_1D( test_name)
    call test_allocate_dist_shared_complex_2D( test_name)
    call test_allocate_dist_shared_complex_3D( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared

  subroutine test_allocate_dist_shared_logical_1D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_allocate_dist_shared_logical_1D'
    character(len=1024), parameter :: test_name_local = 'logical_1D'
    character(len=1024)            :: test_name
    logical, dimension(:), pointer :: d => null()
    type(MPI_WIN)                  :: w
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15) = .true.
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13) .eqv. .true.
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14) .eqv. .true.
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15) .eqv. .true.
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_logical_1D

  subroutine test_allocate_dist_shared_logical_2D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'test_allocate_dist_shared_logical_2D'
    character(len=1024), parameter   :: test_name_local = 'logical_2D'
    character(len=1024)              :: test_name
    logical, dimension(:,:), pointer :: d => null()
    type(MPI_WIN)                    :: w
    logical                          :: test_result
    integer                          :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = .true.
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3) .eqv. .true.
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4) .eqv. .true.
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5) .eqv. .true.
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_logical_2D

  subroutine test_allocate_dist_shared_logical_3D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_allocate_dist_shared_logical_3D'
    character(len=1024), parameter     :: test_name_local = 'logical_3D'
    character(len=1024)                :: test_name
    logical, dimension(:,:,:), pointer :: d => null()
    type(MPI_WIN)                      :: w
    logical                            :: test_result
    integer                            :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10, 5)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20, 6)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30, 7)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = .true.
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3,2) .eqv. .true.
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4,3) .eqv. .true.
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5,4) .eqv. .true.
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_logical_3D

  subroutine test_allocate_dist_shared_int_1D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_allocate_dist_shared_int_1D'
    character(len=1024), parameter :: test_name_local = 'int_1D'
    character(len=1024)            :: test_name
    integer, dimension(:), pointer :: d => null()
    type(MPI_WIN)                  :: w
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15) = 39
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13) == 37
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14) == 38
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15) == 39
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_int_1D

  subroutine test_allocate_dist_shared_int_2D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'test_allocate_dist_shared_int_2D'
    character(len=1024), parameter   :: test_name_local = 'int_2D'
    character(len=1024)              :: test_name
    integer, dimension(:,:), pointer :: d => null()
    type(MPI_WIN)                    :: w
    logical                          :: test_result
    integer                          :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = 39
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3) == 37
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4) == 38
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5) == 39
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_int_2D

  subroutine test_allocate_dist_shared_int_3D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_allocate_dist_shared_int_3D'
    character(len=1024), parameter     :: test_name_local = 'int_3D'
    character(len=1024)                :: test_name
    integer, dimension(:,:,:), pointer :: d => null()
    type(MPI_WIN)                      :: w
    logical                            :: test_result
    integer                            :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10, 5)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20, 6)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30, 7)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = 39
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3,2) == 37
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4,3) == 38
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5,4) == 39
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_int_3D

  subroutine test_allocate_dist_shared_dp_1D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_allocate_dist_shared_dp_1D'
    character(len=1024), parameter  :: test_name_local = 'dp_1D'
    character(len=1024)             :: test_name
    real(dp), dimension(:), pointer :: d => null()
    type(MPI_WIN)                   :: w
    logical                         :: test_result
    integer                         :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15) = 39._dp
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13) == 37._dp
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14) == 38._dp
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15) == 39._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_dp_1D

  subroutine test_allocate_dist_shared_dp_2D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'test_allocate_dist_shared_dp_2D'
    character(len=1024), parameter    :: test_name_local = 'dp_2D'
    character(len=1024)               :: test_name
    real(dp), dimension(:,:), pointer :: d => null()
    type(MPI_WIN)                     :: w
    logical                           :: test_result
    integer                           :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = 39._dp
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3) == 37._dp
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4) == 38._dp
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5) == 39._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_dp_2D

  subroutine test_allocate_dist_shared_dp_3D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_allocate_dist_shared_dp_3D'
    character(len=1024), parameter      :: test_name_local = 'dp_3D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:,:,:), pointer :: d => null()
    type(MPI_WIN)                       :: w
    logical                             :: test_result
    integer                             :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10, 5)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20, 6)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30, 7)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = 39._dp
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3,2) == 37._dp
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4,3) == 38._dp
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5,4) == 39._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_dp_3D

  subroutine test_allocate_dist_shared_complex_1D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'test_allocate_dist_shared_complex_1D'
    character(len=1024), parameter    :: test_name_local = 'complex_1D'
    character(len=1024)               :: test_name
    complex*16, dimension(:), pointer :: d => null()
    type(MPI_WIN)                     :: w
    logical                           :: test_result
    integer                           :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14) = complex( 14._dp, 38_dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15) = complex( 15._dp, 39_dp)
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13) == complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14) == complex( 14._dp, 38_dp)
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15) == complex( 15._dp, 39_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_complex_1D

  subroutine test_allocate_dist_shared_complex_2D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_allocate_dist_shared_complex_2D'
    character(len=1024), parameter      :: test_name_local = 'complex_2D'
    character(len=1024)                 :: test_name
    complex*16, dimension(:,:), pointer :: d => null()
    type(MPI_WIN)                       :: w
    logical                             :: test_result
    integer                             :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = complex( 14._dp, 38_dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = complex( 15._dp, 39_dp)
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3) == complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4) == complex( 14._dp, 38_dp)
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5) == complex( 15._dp, 39_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_complex_2D

  subroutine test_allocate_dist_shared_complex_3D( test_name_parent)
    ! Test the allocate_dist_shared subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_allocate_dist_shared_complex_3D'
    character(len=1024), parameter        :: test_name_local = 'complex_3D'
    character(len=1024)                   :: test_name
    complex*16, dimension(:,:,:), pointer :: d => null()
    type(MPI_WIN)                         :: w
    logical                               :: test_result
    integer                               :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20, 10, 5)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30, 20, 6)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 30, 7)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = complex( 14._dp, 38_dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = complex( 15._dp, 39_dp)
    call sync

    ! ! Test if other processes in the nodes can see this
    test_result = .true.
    if (par%node_ID == 0 .and. par%i_node == 1) test_result = d(13,3,2) == complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%i_node == 1) test_result = d(14,4,3) == complex( 14._dp, 38_dp)
    if (par%node_ID == 2 .and. par%i_node == 1) test_result = d(15,5,4) == complex( 15._dp, 39_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_complex_3D

end module ut_mpi_allocate_dist_shared
