module ut_mpi_gather_dist_shared_to_all

  ! Unit tests for MPI hybrid distributed/shared memory code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, &
    gather_dist_shared_to_all
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use ut_mpi_gather_dist_shared_to_primary, only: simple_nih_sizes

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
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_gather_dist_shared_to_all_logical_1D'
    character(len=1024), parameter             :: test_name_local = 'logical_1D'
    character(len=1024)                        :: test_name
    integer                                    :: n_tot, i1, i2, n
    integer                                    :: i1_node, i2_node, n_node
    integer                                    :: i1_nih, i2_nih, n_nih
    integer                                    :: i1_hle, i2_hle, n_hle
    integer                                    :: i1_hli, i2_hli, n_hli
    integer                                    :: i1_hre, i2_hre, n_hre
    integer                                    :: i1_hri, i2_hri, n_hri
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
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih)
    d_nih( i1_nih:i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72) = .true.
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot)

    ! Evaluate test result and write to output file
    test_result = d_tot( 13) .and. d_tot( 37) .and. d_tot( 72)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_1D

  subroutine test_gather_dist_shared_to_all_logical_2D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_gather_dist_shared_to_all_logical_2D'
    character(len=1024), parameter               :: test_name_local = 'logical_2D'
    character(len=1024)                          :: test_name
    integer                                      :: n_tot, i1, i2, n, nz
    integer                                      :: i1_node, i2_node, n_node
    integer                                      :: i1_nih, i2_nih, n_nih
    integer                                      :: i1_hle, i2_hle, n_hle
    integer                                      :: i1_hli, i2_hli, n_hli
    integer                                      :: i1_hre, i2_hre, n_hre
    integer                                      :: i1_hri, i2_hri, n_hri
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
    nz = 3
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
    d_nih( i1_nih:i2_nih,1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3) = .true.
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, d_tot)

    ! Evaluate test result and write to output file
    test_result = d_tot( 13,1) .and. d_tot( 37,2) .and. d_tot( 72,3)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_2D

  subroutine test_gather_dist_shared_to_all_logical_3D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_gather_dist_shared_to_all_logical_3D'
    character(len=1024), parameter                 :: test_name_local = 'logical_3D'
    character(len=1024)                            :: test_name
    integer                                        :: n_tot, i1, i2, n, nz, nl
    integer                                        :: i1_node, i2_node, n_node
    integer                                        :: i1_nih, i2_nih, n_nih
    integer                                        :: i1_hle, i2_hle, n_hle
    integer                                        :: i1_hli, i2_hli, n_hli
    integer                                        :: i1_hre, i2_hre, n_hre
    integer                                        :: i1_hri, i2_hri, n_hri
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
    nz = 3
    nl = 5
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
    d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1,2) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2,3) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3,5) = .true.
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz, nl)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, nl, d_tot)

    ! Evaluate test result and write to output file
    test_result = d_tot( 13,1,2) .and. d_tot( 37,2,3) .and. d_tot( 72,3,5)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_3D

  subroutine test_gather_dist_shared_to_all_int_1D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_gather_dist_shared_to_all_int_1D'
    character(len=1024), parameter             :: test_name_local = 'int_1D'
    character(len=1024)                        :: test_name
    integer                                    :: n_tot, i1, i2, n
    integer                                    :: i1_node, i2_node, n_node
    integer                                    :: i1_nih, i2_nih, n_nih
    integer                                    :: i1_hle, i2_hle, n_hle
    integer                                    :: i1_hli, i2_hli, n_hli
    integer                                    :: i1_hre, i2_hre, n_hre
    integer                                    :: i1_hri, i2_hri, n_hri
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
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih)
    d_nih( i1_nih:i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13) = 1
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37) = 2
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72) = 3
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13) == 1) .and. &
      (d_tot( 37) == 2) .and. &
      (d_tot( 72) == 3)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_1D

  subroutine test_gather_dist_shared_to_all_int_2D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_gather_dist_shared_to_all_int_2D'
    character(len=1024), parameter               :: test_name_local = 'int_2D'
    character(len=1024)                          :: test_name
    integer                                      :: n_tot, i1, i2, n, nz
    integer                                      :: i1_node, i2_node, n_node
    integer                                      :: i1_nih, i2_nih, n_nih
    integer                                      :: i1_hle, i2_hle, n_hle
    integer                                      :: i1_hli, i2_hli, n_hli
    integer                                      :: i1_hre, i2_hre, n_hre
    integer                                      :: i1_hri, i2_hri, n_hri
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
    nz = 3
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
    d_nih( i1_nih:i2_nih,1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1) = 1
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2) = 2
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3) = 3
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1) == 1) .and. &
      (d_tot( 37,2) == 2) .and. &
      (d_tot( 72,3) == 3)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_2D

  subroutine test_gather_dist_shared_to_all_int_3D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_gather_dist_shared_to_all_int_3D'
    character(len=1024), parameter                 :: test_name_local = 'int_3D'
    character(len=1024)                            :: test_name
    integer                                        :: n_tot, i1, i2, n, nz, nl
    integer                                        :: i1_node, i2_node, n_node
    integer                                        :: i1_nih, i2_nih, n_nih
    integer                                        :: i1_hle, i2_hle, n_hle
    integer                                        :: i1_hli, i2_hli, n_hli
    integer                                        :: i1_hre, i2_hre, n_hre
    integer                                        :: i1_hri, i2_hri, n_hri
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
    nz = 3
    nl = 5
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
    d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1,2) = 1
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2,3) = 2
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3,5) = 3
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz, nl)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, nl, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1,2) == 1) .and. &
      (d_tot( 37,2,3) == 2) .and. &
      (d_tot( 72,3,5) == 3)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_3D

  subroutine test_gather_dist_shared_to_all_dp_1D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_gather_dist_shared_to_all_dp_1D'
    character(len=1024), parameter              :: test_name_local = 'dp_1D'
    character(len=1024)                         :: test_name
    integer                                     :: n_tot, i1, i2, n
    integer                                     :: i1_node, i2_node, n_node
    integer                                     :: i1_nih, i2_nih, n_nih
    integer                                     :: i1_hle, i2_hle, n_hle
    integer                                     :: i1_hli, i2_hli, n_hli
    integer                                     :: i1_hre, i2_hre, n_hre
    integer                                     :: i1_hri, i2_hri, n_hri
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
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih)
    d_nih( i1_nih:i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13) = 1._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37) = 2._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72) = 3._dp
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13) == 1._dp) .and. &
      (d_tot( 37) == 2._dp) .and. &
      (d_tot( 72) == 3._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_1D

  subroutine test_gather_dist_shared_to_all_dp_2D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_gather_dist_shared_to_all_dp_2D'
    character(len=1024), parameter                :: test_name_local = 'dp_2D'
    character(len=1024)                           :: test_name
    integer                                       :: n_tot, i1, i2, n, nz
    integer                                       :: i1_node, i2_node, n_node
    integer                                       :: i1_nih, i2_nih, n_nih
    integer                                       :: i1_hle, i2_hle, n_hle
    integer                                       :: i1_hli, i2_hli, n_hli
    integer                                       :: i1_hre, i2_hre, n_hre
    integer                                       :: i1_hri, i2_hri, n_hri
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
    nz = 3
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
    d_nih( i1_nih:i2_nih,1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1) = 1._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2) = 2._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3) = 3._dp
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1) == 1._dp) .and. &
      (d_tot( 37,2) == 2._dp) .and. &
      (d_tot( 72,3) == 3._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_2D

  subroutine test_gather_dist_shared_to_all_dp_3D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_gather_dist_shared_to_all_dp_3D'
    character(len=1024), parameter                  :: test_name_local = 'dp_3D'
    character(len=1024)                             :: test_name
    integer                                         :: n_tot, i1, i2, n, nz, nl
    integer                                         :: i1_node, i2_node, n_node
    integer                                         :: i1_nih, i2_nih, n_nih
    integer                                         :: i1_hle, i2_hle, n_hle
    integer                                         :: i1_hli, i2_hli, n_hli
    integer                                         :: i1_hre, i2_hre, n_hre
    integer                                         :: i1_hri, i2_hri, n_hri
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
    nz = 3
    nl = 5
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
    d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1,2) = 1._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2,3) = 2._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3,5) = 3._dp
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz, nl)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, nl, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1,2) == 1._dp) .and. &
      (d_tot( 37,2,3) == 2._dp) .and. &
      (d_tot( 72,3,5) == 3._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_3D

  subroutine test_gather_dist_shared_to_all_complex_1D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_gather_dist_shared_to_all_complex_1D'
    character(len=1024), parameter                :: test_name_local = 'complex_1D'
    character(len=1024)                           :: test_name
    integer                                       :: n_tot, i1, i2, n
    integer                                       :: i1_node, i2_node, n_node
    integer                                       :: i1_nih, i2_nih, n_nih
    integer                                       :: i1_hle, i2_hle, n_hle
    integer                                       :: i1_hli, i2_hli, n_hli
    integer                                       :: i1_hre, i2_hre, n_hre
    integer                                       :: i1_hri, i2_hri, n_hri
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
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih)
    d_nih( i1_nih:i2_nih) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13) = complex( 1._dp, 17._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37) = complex( 2._dp, 17._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72) = complex( 3._dp, 17._dp)
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13) == complex( 1._dp, 17._dp)) .and. &
      (d_tot( 37) == complex( 2._dp, 17._dp)) .and. &
      (d_tot( 72) == complex( 3._dp, 17._dp))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_1D

  subroutine test_gather_dist_shared_to_all_complex_2D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_gather_dist_shared_to_all_complex_2D'
    character(len=1024), parameter                  :: test_name_local = 'complex_2D'
    character(len=1024)                             :: test_name
    integer                                         :: n_tot, i1, i2, n, nz
    integer                                         :: i1_node, i2_node, n_node
    integer                                         :: i1_nih, i2_nih, n_nih
    integer                                         :: i1_hle, i2_hle, n_hle
    integer                                         :: i1_hli, i2_hli, n_hli
    integer                                         :: i1_hre, i2_hre, n_hre
    integer                                         :: i1_hri, i2_hri, n_hri
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
    nz = 3
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
    d_nih( i1_nih:i2_nih,1:nz) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1) = complex( 1._dp, 17._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2) = complex( 2._dp, 17._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3) = complex( 3._dp, 17._dp)
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1) == complex( 1._dp, 17._dp)) .and. &
      (d_tot( 37,2) == complex( 2._dp, 17._dp)) .and. &
      (d_tot( 72,3) == complex( 3._dp, 17._dp))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_2D

  subroutine test_gather_dist_shared_to_all_complex_3D( test_name_parent)
    !< Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                    :: routine_name = 'test_gather_dist_shared_to_all_complex_3D'
    character(len=1024), parameter                    :: test_name_local = 'complex_3D'
    character(len=1024)                               :: test_name
    integer                                           :: n_tot, i1, i2, n, nz, nl
    integer                                           :: i1_node, i2_node, n_node
    integer                                           :: i1_nih, i2_nih, n_nih
    integer                                           :: i1_hle, i2_hle, n_hle
    integer                                           :: i1_hli, i2_hli, n_hli
    integer                                           :: i1_hre, i2_hre, n_hre
    integer                                           :: i1_hri, i2_hri, n_hri
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
    nz = 3
    nl = 5
    call simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
      i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
      i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! Allocate node-shared memory including halos
    call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
    d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d_nih(13,1,2) = complex( 1._dp, 17._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih(37,2,3) = complex( 2._dp, 17._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih(72,3,5) = complex( 3._dp, 17._dp)
    call sync

    ! Gather data to the primary
    call allocate_dist_shared( d_tot, wd_tot, n_tot, nz, nl)
    call gather_dist_shared_to_all( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, nz, nl, d_tot)

    ! Evaluate test result and write to output file
    test_result = &
      (d_tot( 13,1,2) == complex( 1._dp, 17._dp)) .and. &
      (d_tot( 37,2,3) == complex( 2._dp, 17._dp)) .and. &
      (d_tot( 72,3,5) == complex( 3._dp, 17._dp))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    call deallocate_dist_shared( d_tot, wd_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_3D

end module ut_mpi_gather_dist_shared_to_all
