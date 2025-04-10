module ut_mpi_gather_dist_shared_to_primary

  ! Unit tests for MPI hybrid distributed/shared memory code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, &
    gather_dist_shared_to_primary
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, &
    MPI_INTEGER, MPI_ALLGATHER

  implicit none

  private

  public :: test_gather_dist_shared_to_primary

contains

  subroutine test_gather_dist_shared_to_primary( test_name_parent)
    ! Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_gather_dist_shared_to_primary'
    character(len=1024), parameter :: test_name_local = 'gather_dist_shared_to_primary'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_gather_dist_shared_to_primary_logical_1D( test_name)
    call test_gather_dist_shared_to_primary_logical_2D( test_name)
    call test_gather_dist_shared_to_primary_logical_3D( test_name)

    call test_gather_dist_shared_to_primary_int_1D( test_name)
    call test_gather_dist_shared_to_primary_int_2D( test_name)
    call test_gather_dist_shared_to_primary_int_3D( test_name)

    call test_gather_dist_shared_to_primary_dp_1D( test_name)
    call test_gather_dist_shared_to_primary_dp_2D( test_name)
    call test_gather_dist_shared_to_primary_dp_3D( test_name)

    call test_gather_dist_shared_to_primary_complex_1D( test_name)
    call test_gather_dist_shared_to_primary_complex_2D( test_name)
    call test_gather_dist_shared_to_primary_complex_3D( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary

  subroutine test_gather_dist_shared_to_primary_logical_1D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_gather_dist_shared_to_primary_logical_1D'
    character(len=1024), parameter     :: test_name_local = 'logical_1D'
    character(len=1024)                :: test_name
    integer                            :: n_tot, i1, i2, n
    integer                            :: i1_node, i2_node, n_node
    integer                            :: i1_nih, i2_nih, n_nih
    integer                            :: i1_hle, i2_hle, n_hle
    integer                            :: i1_hli, i2_hli, n_hli
    integer                            :: i1_hre, i2_hre, n_hre
    integer                            :: i1_hri, i2_hri, n_hri
    logical, dimension(:), pointer     :: d_nih => null()
    type(MPI_WIN)                      :: wd_nih
    logical, dimension(:), allocatable :: d_tot
    logical                            :: test_result

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
    if (par%primary) then
      allocate( d_tot( n_tot), source = .false.)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13) .and. d_tot( 37) .and. d_tot( 72)
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_logical_1D

  subroutine test_gather_dist_shared_to_primary_logical_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_gather_dist_shared_to_primary_logical_2D'
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
    type(MPI_WIN)                                :: wd_nih
    logical, dimension(:,:), allocatable         :: d_tot
    logical                                      :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = .true.
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = .false.)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13,1) .and. d_tot( 37,2) .and. d_tot( 72,3)
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_logical_2D

  subroutine test_gather_dist_shared_to_primary_logical_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_gather_dist_shared_to_primary_logical_3D'
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
    type(MPI_WIN)                                  :: wd_nih
    logical, dimension(:,:,:), allocatable         :: d_tot
    logical                                        :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = .true.
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = .false.)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13,1,2) .and. d_tot( 37,2,3) .and. d_tot( 72,3,5)
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_logical_3D

  subroutine test_gather_dist_shared_to_primary_int_1D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_gather_dist_shared_to_primary_int_1D'
    character(len=1024), parameter     :: test_name_local = 'int_1D'
    character(len=1024)                :: test_name
    integer                            :: n_tot, i1, i2, n
    integer                            :: i1_node, i2_node, n_node
    integer                            :: i1_nih, i2_nih, n_nih
    integer                            :: i1_hle, i2_hle, n_hle
    integer                            :: i1_hli, i2_hli, n_hli
    integer                            :: i1_hre, i2_hre, n_hre
    integer                            :: i1_hri, i2_hri, n_hri
    integer, dimension(:), pointer     :: d_nih => null()
    type(MPI_WIN)                      :: wd_nih
    integer, dimension(:), allocatable :: d_tot
    logical                            :: test_result

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
    if (par%primary) then
      allocate( d_tot( n_tot), source = 0)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13)==1 .and. d_tot( 37)==2 .and. d_tot( 72)==3
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_int_1D

  subroutine test_gather_dist_shared_to_primary_int_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_gather_dist_shared_to_primary_int_2D'
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
    type(MPI_WIN)                                :: wd_nih
    integer, dimension(:,:), allocatable         :: d_tot
    logical                                      :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = 4
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = 5
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = 6
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = 0)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13,1)==4 .and. d_tot( 37,2)==5 .and. d_tot( 72,3)==6
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_int_2D

  subroutine test_gather_dist_shared_to_primary_int_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_gather_dist_shared_to_primary_int_3D'
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
    type(MPI_WIN)                                  :: wd_nih
    integer, dimension(:,:,:), allocatable         :: d_tot
    logical                                        :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = 4
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = 5
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = 6
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = 0)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13,1,2)==4 .and. d_tot( 37,2,3)==5 .and. d_tot( 72,3,5)==6
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_int_3D

  subroutine test_gather_dist_shared_to_primary_dp_1D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_gather_dist_shared_to_primary_dp_1D'
    character(len=1024), parameter      :: test_name_local = 'dp_1D'
    character(len=1024)                 :: test_name
    integer                             :: n_tot, i1, i2, n
    integer                             :: i1_node, i2_node, n_node
    integer                             :: i1_nih, i2_nih, n_nih
    integer                             :: i1_hle, i2_hle, n_hle
    integer                             :: i1_hli, i2_hli, n_hli
    integer                             :: i1_hre, i2_hre, n_hre
    integer                             :: i1_hri, i2_hri, n_hri
    real(dp), dimension(:), pointer     :: d_nih => null()
    type(MPI_WIN)                       :: wd_nih
    real(dp), dimension(:), allocatable :: d_tot
    logical                             :: test_result

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
    if (par%primary) then
      allocate( d_tot( n_tot), source = 0._dp)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13)==1._dp .and. d_tot( 37)==2._dp .and. d_tot( 72)==3._dp
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_dp_1D

  subroutine test_gather_dist_shared_to_primary_dp_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_gather_dist_shared_to_primary_dp_2D'
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
    type(MPI_WIN)                                 :: wd_nih
    real(dp), dimension(:,:), allocatable         :: d_tot
    logical                                       :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = 4._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = 5._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = 6._dp
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = 0._dp)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13,1)==4._dp .and. d_tot( 37,2)==5._dp .and. d_tot( 72,3)==6._dp
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_dp_2D

  subroutine test_gather_dist_shared_to_primary_dp_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_gather_dist_shared_to_primary_dp_3D'
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
    type(MPI_WIN)                                   :: wd_nih
    real(dp), dimension(:,:,:), allocatable         :: d_tot
    logical                                         :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = 4._dp
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = 5._dp
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = 6._dp
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = 0._dp)
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = d_tot( 13,1,2)==4._dp .and. d_tot( 37,2,3)==5._dp .and. d_tot( 72,3,5)==6._dp
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_dp_3D

  subroutine test_gather_dist_shared_to_primary_complex_1D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_gather_dist_shared_to_primary_complex_1D'
    character(len=1024), parameter        :: test_name_local = 'complex_1D'
    character(len=1024)                   :: test_name
    integer                               :: n_tot, i1, i2, n
    integer                               :: i1_node, i2_node, n_node
    integer                               :: i1_nih, i2_nih, n_nih
    integer                               :: i1_hle, i2_hle, n_hle
    integer                               :: i1_hli, i2_hli, n_hli
    integer                               :: i1_hre, i2_hre, n_hre
    integer                               :: i1_hri, i2_hri, n_hri
    complex*16, dimension(:), pointer     :: d_nih => null()
    type(MPI_WIN)                         :: wd_nih
    complex*16, dimension(:), allocatable :: d_tot
    logical                               :: test_result

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
    if (par%primary) then
      allocate( d_tot( n_tot), source = complex( 0._dp, 0._dp))
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = &
        (d_tot( 13) == complex( 1._dp, 17._dp)) .and. &
        (d_tot( 37) == complex( 2._dp, 17._dp)) .and. &
        (d_tot( 72) == complex( 3._dp, 17._dp))
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_complex_1D

  subroutine test_gather_dist_shared_to_primary_complex_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_gather_dist_shared_to_primary_complex_2D'
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
    type(MPI_WIN)                                   :: wd_nih
    complex*16, dimension(:,:), allocatable         :: d_tot
    logical                                         :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1) = complex( 4._dp, 18._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2) = complex( 5._dp, 18._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3) = complex( 6._dp, 18._dp)
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = complex( 0._dp, 0._dp))
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = &
        (d_tot( 13,1) == complex( 4._dp, 18._dp)) .and. &
        (d_tot( 37,2) == complex( 5._dp, 18._dp)) .and. &
        (d_tot( 72,3) == complex( 6._dp, 18._dp))
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_complex_2D

  subroutine test_gather_dist_shared_to_primary_complex_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                    :: routine_name = 'test_gather_dist_shared_to_primary_complex_3D'
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
    type(MPI_WIN)                                     :: wd_nih
    complex*16, dimension(:,:,:), allocatable         :: d_tot
    logical                                           :: test_result

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
    if (par%node_ID == 0 .and. par%node_primary) d_nih( 13,1,2) = complex( 4._dp, 18._dp)
    if (par%node_ID == 1 .and. par%node_primary) d_nih( 37,2,3) = complex( 5._dp, 18._dp)
    if (par%node_ID == 2 .and. par%node_primary) d_nih( 72,3,5) = complex( 6._dp, 18._dp)
    call sync

    ! Gather data to the primary
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = complex( 0._dp, 0._dp))
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call gather_dist_shared_to_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%primary) then
      test_result = &
        (d_tot( 13,1,2) == complex( 4._dp, 18._dp)) .and. &
        (d_tot( 37,2,3) == complex( 5._dp, 18._dp)) .and. &
        (d_tot( 72,3,5) == complex( 6._dp, 18._dp))
      call unit_test( test_result, test_name)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_complex_3D

  subroutine simple_nih_sizes( n_tot, i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
    i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
    i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! In/output variables:
    integer, intent(  out) :: n_tot
    integer, intent(  out) :: i1, i2, n
    integer, intent(  out) :: i1_node, i2_node, n_node
    integer, intent(  out) :: i1_nih, i2_nih, n_nih
    integer, intent(  out) :: i1_hle, i2_hle, n_hle
    integer, intent(  out) :: i1_hli, i2_hli, n_hli
    integer, intent(  out) :: i1_hre, i2_hre, n_hre
    integer, intent(  out) :: i1_hri, i2_hri, n_hri

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'simple_nih_sizes'

    ! Add routine to call stack
    call init_routine( routine_name)

    n_tot = 90

    if (par%node_ID == 0) then
      ! Node 0: 20 elements

      i1_node = 1
      i2_node = 20

      if (par%i == 0) then
        ! Process 0: 12 elements
        i1 = 1
        i2 = 12
        n  = 12
      elseif (par%i == 1) then
        ! Process 1: 8 elements
        i1 = 13
        i2 = 20
        n  = 8
      end if

    elseif (par%node_ID == 1) then
      ! Node 1: 30 elements

      i1_node = 21
      i2_node = 50

      if (par%i == 2) then
        ! Process 2: 18 elements
        i1 = 21
        i2 = 38
        n  = 18
      elseif (par%i == 3) then
        ! Process 3: 9 elements
        i1 = 39
        i2 = 47
        n  = 9
      elseif (par%i == 4) then
        ! Process 4: 3 elements
        i1 = 48
        i2 = 50
        n  = 3
      end if

    elseif (par%node_ID == 2) then
      ! Node 2: 40 elements

      i1_node = 51
      i2_node = 90

      if (par%i == 5) then
        ! Process 5: 25 elements
        i1 = 51
        i2 = 75
        n  = 25
      elseif (par%i == 6) then
        ! Process 6: 15 elements
        i1 = 76
        i2 = 90
        n  = 15
      end if

    end if

    ! Node 0: no left border
    if (par%node_ID == 0) then
      i1_hle =  0
      i2_hle = -1
      i1_hli =  0
      i2_hli = -1
    end if

    ! Border between nodes 0 and 1: 4 halos on the left, 5 on the right
    if (par%node_ID == 0) then
      i1_hri = i2_node-3
      i2_hri = i2_node
      i1_hre = i2_node+1
      i2_hre = i2_node+5
    end if
    if (par%node_ID == 1) then
      i1_hle = i1_node-4
      i2_hle = i1_node-1
      i1_hli = i1_node
      i2_hli = i1_node+4
    end if

    ! Border between nodes 1 and 2: 6 halos on the left, 7 on the right
    if (par%node_ID == 1) then
      i1_hri = i2_node-5
      i2_hri = i2_node
      i1_hre = i2_node+1
      i2_hre = i2_node+7
    end if
    if (par%node_ID == 2) then
      i1_hle = i1_node-6
      i2_hle = i1_node-1
      i1_hli = i1_node
      i2_hli = i1_node+6
    end if

    ! Node 2: no right border
    if (par%node_ID == 2) then
      i1_hre =  0
      i2_hre = -1
      i1_hri =  0
      i2_hri = -1
    end if

    if (par%node_ID == 0) then
      i1_nih = 1
      i2_nih = i2_hre
    elseif (par%node_ID == 1) then
      i1_nih = i1_hle
      i2_nih = i2_hre
    elseif (par%node_ID == 2) then
      i1_nih = i1_hle
      i2_nih = i2_node
    end if

    n_node = i2_node + 1 - i1_node
    n_nih  = i2_nih  + 1 - i1_nih
    n_hle  = i2_hle  + 1 - i1_hle
    n_hli  = i2_hli  + 1 - i1_hli
    n_hre  = i2_hre  + 1 - i1_hre
    n_hri  = i2_hri  + 1 - i1_hri

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine simple_nih_sizes

  subroutine print_simple_nih_sizes( i1, i2, n, i1_node, i2_node, n_node, i1_nih, i2_nih, n_nih, &
    i1_hle, i2_hle, n_hle, i1_hli, i2_hli, n_hli, &
    i1_hre, i2_hre, n_hre, i1_hri, i2_hri, n_hri)

    ! In/output variables:
    integer, intent(in   ) :: i1, i2, n
    integer, intent(in   ) :: i1_node, i2_node, n_node
    integer, intent(in   ) :: i1_nih, i2_nih, n_nih
    integer, intent(in   ) :: i1_hle, i2_hle, n_hle
    integer, intent(in   ) :: i1_hli, i2_hli, n_hli
    integer, intent(in   ) :: i1_hre, i2_hre, n_hre
    integer, intent(in   ) :: i1_hri, i2_hri, n_hri

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_simple_nih_sizes'
    integer, dimension(0:par%n-1)  :: node_all, process_all
    integer, dimension(0:par%n-1)  :: i1_all, i2_all, n_all
    integer, dimension(0:par%n-1)  :: i1_node_all, i2_node_all, n_node_all
    integer, dimension(0:par%n-1)  :: i1_nih_all, i2_nih_all, n_nih_all
    integer, dimension(0:par%n-1)  :: i1_hle_all, i2_hle_all, n_hle_all
    integer, dimension(0:par%n-1)  :: i1_hli_all, i2_hli_all, n_hli_all
    integer, dimension(0:par%n-1)  :: i1_hre_all, i2_hre_all, n_hre_all
    integer, dimension(0:par%n-1)  :: i1_hri_all, i2_hri_all, n_hri_all
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call MPI_ALLGATHER( par%node_ID  , 1, MPI_INTEGER, node_all    , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( par%i        , 1, MPI_INTEGER, process_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( i1     , 1, MPI_INTEGER, i1_all     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( i2     , 1, MPI_INTEGER, i2_all     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n      , 1, MPI_INTEGER, n_all      , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( i1_node, 1, MPI_INTEGER, i1_node_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( i2_node, 1, MPI_INTEGER, i2_node_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_node , 1, MPI_INTEGER, n_node_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( i1_nih , 1, MPI_INTEGER, i1_nih_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( i2_nih , 1, MPI_INTEGER, i2_nih_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_nih  , 1, MPI_INTEGER, n_nih_all  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( i1_hle , 1, MPI_INTEGER, i1_hle_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( i2_hle , 1, MPI_INTEGER, i2_hle_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_hle  , 1, MPI_INTEGER, n_hle_all  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( i1_hli , 1, MPI_INTEGER, i1_hli_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( i2_hli , 1, MPI_INTEGER, i2_hli_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_hli  , 1, MPI_INTEGER, n_hli_all  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( i1_hre , 1, MPI_INTEGER, i1_hre_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( i2_hre , 1, MPI_INTEGER, i2_hre_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_hre  , 1, MPI_INTEGER, n_hre_all  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( i1_hri , 1, MPI_INTEGER, i1_hri_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( i2_hri , 1, MPI_INTEGER, i2_hri_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_hri  , 1, MPI_INTEGER, n_hri_all  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if (par%primary) then
      write(0,'(A,7I8)') 'Node     :', node_all
      write(0,'(A,7I8)') 'Process  :', process_all
      write(0,*) ''
      write(0,'(A,7I8)') 'n       :', n_all
      write(0,'(A,7I8)') 'i1      :', i1_all
      write(0,'(A,7I8)') 'i2      :', i2_all
      write(0,*) ''
      write(0,'(A,7I8)') 'n_node  :', n_node_all
      write(0,'(A,7I8)') 'i1_node :', i1_node_all
      write(0,'(A,7I8)') 'i2_node :', i2_node_all
      write(0,*) ''
      write(0,'(A,7I8)') 'n_nih   :', n_nih_all
      write(0,'(A,7I8)') 'i1_nih  :', i1_nih_all
      write(0,'(A,7I8)') 'i2_nih  :', i2_nih_all
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hle   :', n_hle_all
      write(0,'(A,7I8)') 'i1_hle  :', i1_hle_all
      write(0,'(A,7I8)') 'i2_hle  :', i2_hle_all
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hli   :', n_hli_all
      write(0,'(A,7I8)') 'i1_hli  :', i1_hli_all
      write(0,'(A,7I8)') 'i2_hli  :', i2_hli_all
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hre   :', n_hre_all
      write(0,'(A,7I8)') 'i1_hre  :', i1_hre_all
      write(0,'(A,7I8)') 'i2_hre  :', i2_hre_all
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hri   :', n_hri_all
      write(0,'(A,7I8)') 'i1_hri  :', i1_hri_all
      write(0,'(A,7I8)') 'i2_hri  :', i2_hri_all
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine print_simple_nih_sizes

end module ut_mpi_gather_dist_shared_to_primary
