module ut_mpi_gather_dist_shared_to_all

  ! Unit tests for MPI hybrid distributed/shared memory code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_distributed_shared_memory
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use ut_mpi_allocate_dist_shared, only: setup_simple_parallel_array_info

  implicit none

  private

  public :: test_gather_dist_shared_to_all

contains

  subroutine test_gather_dist_shared_to_all( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_gather_dist_shared_to_all'
    character(len=1024), parameter :: test_name_local = 'gather_dist_shared_to_all'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_gather_dist_shared_to_all_logical_1D( test_name)
    call test_gather_dist_shared_to_all_logical_2D( test_name)
    call test_gather_dist_shared_to_all_logical_3D( test_name)

    call test_gather_dist_shared_to_all_int_1D( test_name)
    call test_gather_dist_shared_to_all_int_2D( test_name)
    call test_gather_dist_shared_to_all_int_3D( test_name)

    call test_gather_dist_shared_to_all_dp_1D( test_name)
    call test_gather_dist_shared_to_all_dp_2D( test_name)
    call test_gather_dist_shared_to_all_dp_3D( test_name)

    call test_gather_dist_shared_to_all_complex_1D( test_name)
    call test_gather_dist_shared_to_all_complex_2D( test_name)
    call test_gather_dist_shared_to_all_complex_3D( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all

  subroutine test_gather_dist_shared_to_all_logical_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_gather_dist_shared_to_all_logical_1D'
    character(len=1024), parameter             :: test_name_local = 'logical_1D'
    character(len=1024)                        :: test_name
    type(type_par_arr_info)                    :: pai
    logical, dimension(:), contiguous, pointer :: d_nih => null()
    logical, dimension(:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                              :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = .true.
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13) .eqv. .true.) .and. &
      (d_tot( 37) .eqv. .true.) .and. &
      (d_tot( 72) .eqv. .true.)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_1D

  subroutine test_gather_dist_shared_to_all_logical_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_gather_dist_shared_to_all_logical_2D'
    character(len=1024), parameter               :: test_name_local = 'logical_2D'
    character(len=1024)                          :: test_name
    type(type_par_arr_info)                      :: pai
    integer                                      :: nz
    logical, dimension(:,:), contiguous, pointer :: d_nih => null()
    logical, dimension(:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = .true.
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1) .eqv. .true.) .and. &
      (d_tot( 37,2) .eqv. .true.) .and. &
      (d_tot( 72,3) .eqv. .true.)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_2D

  subroutine test_gather_dist_shared_to_all_logical_3D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_gather_dist_shared_to_all_logical_3D'
    character(len=1024), parameter                 :: test_name_local = 'logical_3D'
    character(len=1024)                            :: test_name
    type(type_par_arr_info)                        :: pai
    integer                                        :: nz, nl
    logical, dimension(:,:,:), contiguous, pointer :: d_nih => null()
    logical, dimension(:,:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                  :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz, nl)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = .true.
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, nl, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1,2) .eqv. .true.) .and. &
      (d_tot( 37,2,3) .eqv. .true.) .and. &
      (d_tot( 72,3,5) .eqv. .true.)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_3D

  subroutine test_gather_dist_shared_to_all_int_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_gather_dist_shared_to_all_int_1D'
    character(len=1024), parameter             :: test_name_local = 'int_1D'
    character(len=1024)                        :: test_name
    type(type_par_arr_info)                    :: pai
    integer, dimension(:), contiguous, pointer :: d_nih => null()
    integer, dimension(:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                              :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = 42
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = 42
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = 42
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13) == 42) .and. &
      (d_tot( 37) == 42) .and. &
      (d_tot( 72) == 42)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_1D

  subroutine test_gather_dist_shared_to_all_int_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_gather_dist_shared_to_all_int_2D'
    character(len=1024), parameter               :: test_name_local = 'int_2D'
    character(len=1024)                          :: test_name
    type(type_par_arr_info)                      :: pai
    integer                                      :: nz
    integer, dimension(:,:), contiguous, pointer :: d_nih => null()
    integer, dimension(:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = 42
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = 42
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = 42
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1) == 42) .and. &
      (d_tot( 37,2) == 42) .and. &
      (d_tot( 72,3) == 42)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_2D

  subroutine test_gather_dist_shared_to_all_int_3D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_gather_dist_shared_to_all_int_3D'
    character(len=1024), parameter                 :: test_name_local = 'int_3D'
    character(len=1024)                            :: test_name
    type(type_par_arr_info)                        :: pai
    integer                                        :: nz, nl
    integer, dimension(:,:,:), contiguous, pointer :: d_nih => null()
    integer, dimension(:,:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                  :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz, nl)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = 42
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = 42
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = 42
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, nl, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1,2) == 42) .and. &
      (d_tot( 37,2,3) == 42) .and. &
      (d_tot( 72,3,5) == 42)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_3D

  subroutine test_gather_dist_shared_to_all_dp_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_gather_dist_shared_to_all_dp_1D'
    character(len=1024), parameter              :: test_name_local = 'dp_1D'
    character(len=1024)                         :: test_name
    type(type_par_arr_info)                     :: pai
    real(dp), dimension(:), contiguous, pointer :: d_nih => null()
    real(dp), dimension(:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                               :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = 42._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = 42._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = 42._dp
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13) == 42._dp) .and. &
      (d_tot( 37) == 42._dp) .and. &
      (d_tot( 72) == 42._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_1D

  subroutine test_gather_dist_shared_to_all_dp_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_gather_dist_shared_to_all_dp_2D'
    character(len=1024), parameter                :: test_name_local = 'dp_2D'
    character(len=1024)                           :: test_name
    type(type_par_arr_info)                       :: pai
    integer                                       :: nz
    real(dp), dimension(:,:), contiguous, pointer :: d_nih => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                 :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = 42._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = 42._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = 42._dp
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1) == 42._dp) .and. &
      (d_tot( 37,2) == 42._dp) .and. &
      (d_tot( 72,3) == 42._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_2D

  subroutine test_gather_dist_shared_to_all_dp_3D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_gather_dist_shared_to_all_dp_3D'
    character(len=1024), parameter                  :: test_name_local = 'dp_3D'
    character(len=1024)                             :: test_name
    type(type_par_arr_info)                         :: pai
    integer                                         :: nz, nl
    real(dp), dimension(:,:,:), contiguous, pointer :: d_nih => null()
    real(dp), dimension(:,:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                   :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz, nl)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = 42._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = 42._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = 42._dp
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, nl, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1,2) == 42._dp) .and. &
      (d_tot( 37,2,3) == 42._dp) .and. &
      (d_tot( 72,3,5) == 42._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_3D

  subroutine test_gather_dist_shared_to_all_complex_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_gather_dist_shared_to_all_complex_1D'
    character(len=1024), parameter                :: test_name_local = 'complex_1D'
    character(len=1024)                           :: test_name
    type(type_par_arr_info)                       :: pai
    complex*16, dimension(:), contiguous, pointer :: d_nih => null()
    complex*16, dimension(:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                 :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37) = complex( 13._dp, 37._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72) = complex( 13._dp, 37._dp)
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13) == complex( 13._dp, 37._dp)) .and. &
      (d_tot( 37) == complex( 13._dp, 37._dp)) .and. &
      (d_tot( 72) == complex( 13._dp, 37._dp))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_1D

  subroutine test_gather_dist_shared_to_all_complex_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_gather_dist_shared_to_all_complex_2D'
    character(len=1024), parameter                  :: test_name_local = 'complex_2D'
    character(len=1024)                             :: test_name
    type(type_par_arr_info)                         :: pai
    integer                                         :: nz
    complex*16, dimension(:,:), contiguous, pointer :: d_nih => null()
    complex*16, dimension(:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                   :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = complex( 13._dp, 37._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = complex( 13._dp, 37._dp)
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1) == complex( 13._dp, 37._dp)) .and. &
      (d_tot( 37,2) == complex( 13._dp, 37._dp)) .and. &
      (d_tot( 72,3) == complex( 13._dp, 37._dp))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_2D

  subroutine test_gather_dist_shared_to_all_complex_3D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                    :: routine_name = 'test_gather_dist_shared_to_all_complex_3D'
    character(len=1024), parameter                    :: test_name_local = 'complex_3D'
    character(len=1024)                               :: test_name
    type(type_par_arr_info)                           :: pai
    integer                                           :: nz, nl
    complex*16, dimension(:,:,:), contiguous, pointer :: d_nih => null()
    complex*16, dimension(:,:,:), contiguous, pointer :: d_tot => null()
    type(MPI_WIN)                                     :: wd_nih, wd_tot
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
    call allocate_dist_shared( d_tot, wd_tot, pai%n, nz, nl)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = complex( 13._dp, 37._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = complex( 13._dp, 37._dp)
    call sync

    ! Gather data
    call gather_dist_shared_to_all( pai, nz, nl, d_nih, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1,2) == complex( 13._dp, 37._dp)) .and. &
      (d_tot( 37,2,3) == complex( 13._dp, 37._dp)) .and. &
      (d_tot( 72,3,5) == complex( 13._dp, 37._dp))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_3D

end module ut_mpi_gather_dist_shared_to_all
