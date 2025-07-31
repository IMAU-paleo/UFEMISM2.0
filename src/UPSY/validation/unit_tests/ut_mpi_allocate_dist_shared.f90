module ut_mpi_allocate_dist_shared

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_distributed_shared_memory
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD

  implicit none

  private

  public :: test_allocate_dist_shared, setup_simple_parallel_array_info

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
    character(len=1024), parameter             :: routine_name = 'test_allocate_dist_shared_logical_1D'
    character(len=1024), parameter             :: test_name_local = 'logical_1D'
    character(len=1024)                        :: test_name
    type(type_par_arr_info)                    :: pai
    logical, dimension(:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                              :: wd_nih
    logical                                    :: test_result
    integer                                    :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih)
    d_nih( pai%i1_nih:pai%i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = .true.
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13) .eqv. .true.
    if (par%node_ID == 1) test_result = d_nih( 37) .eqv. .true.
    if (par%node_ID == 2) test_result = d_nih( 72) .eqv. .true.
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter               :: routine_name = 'test_allocate_dist_shared_logical_2D'
    character(len=1024), parameter               :: test_name_local = 'logical_2D'
    character(len=1024)                          :: test_name
    type(type_par_arr_info)                      :: pai
    integer                                      :: nz
    logical, dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                :: wd_nih
    logical                                      :: test_result
    integer                                      :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = .true.
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1) .eqv. .true.
    if (par%node_ID == 1) test_result = d_nih( 37,2) .eqv. .true.
    if (par%node_ID == 2) test_result = d_nih( 72,3) .eqv. .true.
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter                 :: routine_name = 'test_allocate_dist_shared_logical_3D'
    character(len=1024), parameter                 :: test_name_local = 'logical_3D'
    character(len=1024)                            :: test_name
    type(type_par_arr_info)                        :: pai
    integer                                        :: nz, nl
    logical, dimension(:,:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                  :: wd_nih
    logical                                        :: test_result
    integer                                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3
    nl = 5

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz, nl)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz, 1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = .true.
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1,2) .eqv. .true.
    if (par%node_ID == 1) test_result = d_nih( 37,2,3) .eqv. .true.
    if (par%node_ID == 2) test_result = d_nih( 72,3,5) .eqv. .true.
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter             :: routine_name = 'test_allocate_dist_shared_int_1D'
    character(len=1024), parameter             :: test_name_local = 'int_1D'
    character(len=1024)                        :: test_name
    type(type_par_arr_info)                    :: pai
    integer, dimension(:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                              :: wd_nih
    logical                                    :: test_result
    integer                                    :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih)
    d_nih( pai%i1_nih:pai%i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = 42
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = 42
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = 42
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13) == 42
    if (par%node_ID == 1) test_result = d_nih( 37) == 42
    if (par%node_ID == 2) test_result = d_nih( 72) == 42
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter               :: routine_name = 'test_allocate_dist_shared_int_2D'
    character(len=1024), parameter               :: test_name_local = 'int_2D'
    character(len=1024)                          :: test_name
    type(type_par_arr_info)                      :: pai
    integer                                      :: nz
    integer, dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                :: wd_nih
    logical                                      :: test_result
    integer                                      :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = 42
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = 42
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = 42
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1) == 42
    if (par%node_ID == 1) test_result = d_nih( 37,2) == 42
    if (par%node_ID == 2) test_result = d_nih( 72,3) == 42
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter                 :: routine_name = 'test_allocate_dist_shared_int_3D'
    character(len=1024), parameter                 :: test_name_local = 'int_3D'
    character(len=1024)                            :: test_name
    type(type_par_arr_info)                        :: pai
    integer                                        :: nz, nl
    integer, dimension(:,:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                  :: wd_nih
    logical                                        :: test_result
    integer                                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3
    nl = 5

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz, nl)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz, 1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = 42
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = 42
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = 42
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1,2) == 42
    if (par%node_ID == 1) test_result = d_nih( 37,2,3) == 42
    if (par%node_ID == 2) test_result = d_nih( 72,3,5) == 42
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter              :: routine_name = 'test_allocate_dist_shared_dp_1D'
    character(len=1024), parameter              :: test_name_local = 'dp_1D'
    character(len=1024)                         :: test_name
    type(type_par_arr_info)                     :: pai
    real(dp), dimension(:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                               :: wd_nih
    logical                                     :: test_result
    integer                                     :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih)
    d_nih( pai%i1_nih:pai%i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = 42._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = 42._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = 42._dp
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13) == 42._dp
    if (par%node_ID == 1) test_result = d_nih( 37) == 42._dp
    if (par%node_ID == 2) test_result = d_nih( 72) == 42._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter                :: routine_name = 'test_allocate_dist_shared_dp_2D'
    character(len=1024), parameter                :: test_name_local = 'dp_2D'
    character(len=1024)                           :: test_name
    type(type_par_arr_info)                       :: pai
    integer                                       :: nz
    real(dp), dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                 :: wd_nih
    logical                                       :: test_result
    integer                                       :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = 42._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = 42._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = 42._dp
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1) == 42._dp
    if (par%node_ID == 1) test_result = d_nih( 37,2) == 42._dp
    if (par%node_ID == 2) test_result = d_nih( 72,3) == 42._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter                  :: routine_name = 'test_allocate_dist_shared_dp_3D'
    character(len=1024), parameter                  :: test_name_local = 'dp_3D'
    character(len=1024)                             :: test_name
    type(type_par_arr_info)                         :: pai
    integer                                         :: nz, nl
    real(dp), dimension(:,:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                   :: wd_nih
    logical                                         :: test_result
    integer                                         :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3
    nl = 5

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz, nl)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz, 1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = 42._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = 42._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = 42._dp
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1,2) == 42._dp
    if (par%node_ID == 1) test_result = d_nih( 37,2,3) == 42._dp
    if (par%node_ID == 2) test_result = d_nih( 72,3,5) == 42._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter                :: routine_name = 'test_allocate_dist_shared_complex_1D'
    character(len=1024), parameter                :: test_name_local = 'complex_1D'
    character(len=1024)                           :: test_name
    type(type_par_arr_info)                       :: pai
    complex*16, dimension(:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                 :: wd_nih
    logical                                       :: test_result
    integer                                       :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih)
    d_nih( pai%i1_nih:pai%i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = complex( 13._dp, 37._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = complex( 13._dp, 37._dp)
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13) == complex( 13._dp, 37._dp)
    if (par%node_ID == 1) test_result = d_nih( 37) == complex( 13._dp, 37._dp)
    if (par%node_ID == 2) test_result = d_nih( 72) == complex( 13._dp, 37._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter                  :: routine_name = 'test_allocate_dist_shared_complex_2D'
    character(len=1024), parameter                  :: test_name_local = 'complex_2D'
    character(len=1024)                             :: test_name
    type(type_par_arr_info)                         :: pai
    integer                                         :: nz
    complex*16, dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                   :: wd_nih
    logical                                         :: test_result
    integer                                         :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = complex( 13._dp, 37._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = complex( 13._dp, 37._dp)
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1) == complex( 13._dp, 37._dp)
    if (par%node_ID == 1) test_result = d_nih( 37,2) == complex( 13._dp, 37._dp)
    if (par%node_ID == 2) test_result = d_nih( 72,3) == complex( 13._dp, 37._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

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
    character(len=1024), parameter                    :: routine_name = 'test_allocate_dist_shared_complex_3D'
    character(len=1024), parameter                    :: test_name_local = 'complex_3D'
    character(len=1024)                               :: test_name
    type(type_par_arr_info)                           :: pai
    integer                                           :: nz, nl
    complex*16, dimension(:,:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                     :: wd_nih
    logical                                           :: test_result
    integer                                           :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define sizes for a basic hybrid distributed/shared array including halos
    call setup_simple_parallel_array_info( pai)
    nz = 3
    nl = 5

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, pai%n_nih, nz, nl)
    d_nih( pai%i1_nih:pai%i2_nih, 1:nz, 1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = complex( 13._dp, 37._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = complex( 13._dp, 37._dp)
    call sync

    ! ! Test if other processes in the nodes can see this
    if (par%node_ID == 0) test_result = d_nih( 13,1,2) == complex( 13._dp, 37._dp)
    if (par%node_ID == 1) test_result = d_nih( 37,2,3) == complex( 13._dp, 37._dp)
    if (par%node_ID == 2) test_result = d_nih( 72,3,5) == complex( 13._dp, 37._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_allocate_dist_shared_complex_3D

  subroutine setup_simple_parallel_array_info( pai)

    ! In/output variables:
    type(type_par_arr_info), intent(  out) :: pai

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_simple_parallel_array_info'

    ! Add routine to call stack
    call init_routine( routine_name)

    pai%n = 90

    if (par%node_ID == 0) then
      ! Node 0: 20 elements

      pai%i1_node = 1
      pai%i2_node = 20

      if (par%i == 0) then
        ! Process 0: 12 elements
        pai%i1 = 1
        pai%i2 = 12
      elseif (par%i == 1) then
        ! Process 1: 8 elements
        pai%i1 = 13
        pai%i2 = 20
      end if

    elseif (par%node_ID == 1) then
      ! Node 1: 30 elements

      pai%i1_node = 21
      pai%i2_node = 50

      if (par%i == 2) then
        ! Process 2: 18 elements
        pai%i1 = 21
        pai%i2 = 38
      elseif (par%i == 3) then
        ! Process 3: 9 elements
        pai%i1 = 39
        pai%i2 = 47
      elseif (par%i == 4) then
        ! Process 4: 3 elements
        pai%i1 = 48
        pai%i2 = 50
      end if

    elseif (par%node_ID == 2) then
      ! Node 2: 40 elements

      pai%i1_node = 51
      pai%i2_node = 90

      if (par%i == 5) then
        ! Process 5: 25 elements
        pai%i1 = 51
        pai%i2 = 75
      elseif (par%i == 6) then
        ! Process 6: 15 elements
        pai%i1 = 76
        pai%i2 = 90
      end if

    end if

    ! Node 0: no left border
    if (par%node_ID == 0) then
      pai%i1_hle =  0
      pai%i2_hle = -1
      pai%i1_hli =  0
      pai%i2_hli = -1
    end if

    ! Border between nodes 0 and 1: 4 halos on the left, 5 on the right
    if (par%node_ID == 0) then
      pai%i1_hri = pai%i2_node-3
      pai%i2_hri = pai%i2_node
      pai%i1_hre = pai%i2_node+1
      pai%i2_hre = pai%i2_node+5
    end if
    if (par%node_ID == 1) then
      pai%i1_hle = pai%i1_node-4
      pai%i2_hle = pai%i1_node-1
      pai%i1_hli = pai%i1_node
      pai%i2_hli = pai%i1_node+4
    end if

    ! Border between nodes 1 and 2: 6 halos on the left, 7 on the right
    if (par%node_ID == 1) then
      pai%i1_hri = pai%i2_node-5
      pai%i2_hri = pai%i2_node
      pai%i1_hre = pai%i2_node+1
      pai%i2_hre = pai%i2_node+7
    end if
    if (par%node_ID == 2) then
      pai%i1_hle = pai%i1_node-6
      pai%i2_hle = pai%i1_node-1
      pai%i1_hli = pai%i1_node
      pai%i2_hli = pai%i1_node+6
    end if

    ! Node 2: no right border
    if (par%node_ID == 2) then
      pai%i1_hre =  0
      pai%i2_hre = -1
      pai%i1_hri =  0
      pai%i2_hri = -1
    end if

    if (par%node_ID == 0) then
      pai%i1_nih = 1
      pai%i2_nih = pai%i2_hre
    elseif (par%node_ID == 1) then
      pai%i1_nih = pai%i1_hle
      pai%i2_nih = pai%i2_hre
    elseif (par%node_ID == 2) then
      pai%i1_nih = pai%i1_hle
      pai%i2_nih = pai%i2_node
    end if

    pai%n_loc  = pai%i2      + 1 - pai%i1
    pai%n_node = pai%i2_node + 1 - pai%i1_node
    pai%n_nih  = pai%i2_nih  + 1 - pai%i1_nih
    pai%n_hle  = pai%i2_hle  + 1 - pai%i1_hle
    pai%n_hli  = pai%i2_hli  + 1 - pai%i1_hli
    pai%n_hre  = pai%i2_hre  + 1 - pai%i1_hre
    pai%n_hri  = pai%i2_hri  + 1 - pai%i1_hri

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine setup_simple_parallel_array_info

end module ut_mpi_allocate_dist_shared
