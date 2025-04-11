module ut_mpi_distribute_dist_shared_from_primary

  ! Unit tests for MPI hybrid distributed/shared memory code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, &
    distribute_dist_shared_from_primary
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use ut_mpi_gather_dist_shared_to_primary, only: simple_nih_sizes

  implicit none

  private

  public :: test_distribute_dist_shared_from_primary

contains

  subroutine test_distribute_dist_shared_from_primary( test_name_parent)
    ! Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_distribute_dist_shared_from_primary'
    character(len=1024), parameter :: test_name_local = 'distribute_dist_shared_from_primary'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_distribute_dist_shared_from_primary_logical_1D( test_name)
    call test_distribute_dist_shared_from_primary_logical_2D( test_name)
    call test_distribute_dist_shared_from_primary_logical_3D( test_name)

    call test_distribute_dist_shared_from_primary_int_1D( test_name)
    call test_distribute_dist_shared_from_primary_int_2D( test_name)
    call test_distribute_dist_shared_from_primary_int_3D( test_name)

    call test_distribute_dist_shared_from_primary_dp_1D( test_name)
    call test_distribute_dist_shared_from_primary_dp_2D( test_name)
    call test_distribute_dist_shared_from_primary_dp_3D( test_name)

    call test_distribute_dist_shared_from_primary_complex_1D( test_name)
    call test_distribute_dist_shared_from_primary_complex_2D( test_name)
    call test_distribute_dist_shared_from_primary_complex_3D( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary

  subroutine test_distribute_dist_shared_from_primary_logical_1D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_distribute_dist_shared_from_primary_logical_1D'
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
    integer                            :: ierr

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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot), source = .false.)
      d_tot(13) = .true.
      d_tot(37) = .true.
      d_tot(72) = .true.
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13)
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37)
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72)
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_logical_1D

  subroutine test_distribute_dist_shared_from_primary_logical_2D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_distribute_dist_shared_from_primary_logical_2D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = .false.)
      d_tot(13,1) = .true.
      d_tot(37,2) = .true.
      d_tot(72,3) = .true.
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1)
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2)
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3)
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_logical_2D

  subroutine test_distribute_dist_shared_from_primary_logical_3D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_distribute_dist_shared_from_primary_logical_3D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = .false.)
      d_tot(13,1,2) = .true.
      d_tot(37,2,3) = .true.
      d_tot(72,3,5) = .true.
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1,2)
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2,3)
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3,5)
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_logical_3D

  subroutine test_distribute_dist_shared_from_primary_int_1D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_distribute_dist_shared_from_primary_int_1D'
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
    integer                            :: ierr

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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot), source = 0)
      d_tot(13) = 1
      d_tot(37) = 2
      d_tot(72) = 3
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13) == 1
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37) == 2
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72) == 3
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_int_1D

  subroutine test_distribute_dist_shared_from_primary_int_2D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_distribute_dist_shared_from_primary_int_2D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = 0)
      d_tot(13,1) = 1
      d_tot(37,2) = 2
      d_tot(72,3) = 3
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1) == 1
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2) == 2
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3) == 3
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_int_2D

  subroutine test_distribute_dist_shared_from_primary_int_3D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'test_distribute_dist_shared_from_primary_int_3D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = 0)
      d_tot(13,1,2) = 1
      d_tot(37,2,3) = 2
      d_tot(72,3,5) = 3
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1,2) == 1
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2,3) == 2
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3,5) == 3
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_int_3D

  subroutine test_distribute_dist_shared_from_primary_dp_1D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_distribute_dist_shared_from_primary_dp_1D'
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
    integer                             :: ierr

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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot), source = 0._dp)
      d_tot(13) = 1._dp
      d_tot(37) = 2._dp
      d_tot(72) = 3._dp
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13) == 1._dp
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37) == 2._dp
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72) == 3._dp
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_dp_1D

  subroutine test_distribute_dist_shared_from_primary_dp_2D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_distribute_dist_shared_from_primary_dp_2D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = 0._dp)
      d_tot(13,1) = 1._dp
      d_tot(37,2) = 2._dp
      d_tot(72,3) = 3._dp
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1) == 1._dp
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2) == 2._dp
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3) == 3._dp
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_dp_2D

  subroutine test_distribute_dist_shared_from_primary_dp_3D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_distribute_dist_shared_from_primary_dp_3D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = 0._dp)
      d_tot(13,1,2) = 1._dp
      d_tot(37,2,3) = 2._dp
      d_tot(72,3,5) = 3._dp
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1,2) == 1._dp
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2,3) == 2._dp
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3,5) == 3._dp
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_dp_3D

  subroutine test_distribute_dist_shared_from_primary_complex_1D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_distribute_dist_shared_from_primary_complex_1D'
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
    integer                               :: ierr

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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot), source = complex( 0._dp, 0._dp))
      d_tot(13) = complex( 1._dp, 17._dp)
      d_tot(37) = complex( 2._dp, 17._dp)
      d_tot(72) = complex( 3._dp, 17._dp)
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, n_tot)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13) == complex( 1._dp, 17._dp)
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37) == complex( 2._dp, 17._dp)
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72) == complex( 3._dp, 17._dp)
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_complex_1D

  subroutine test_distribute_dist_shared_from_primary_complex_2D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'test_distribute_dist_shared_from_primary_complex_2D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz), source = complex( 0._dp, 0._dp))
      d_tot(13,1) = complex( 1._dp, 17._dp)
      d_tot(37,2) = complex( 2._dp, 17._dp)
      d_tot(72,3) = complex( 3._dp, 17._dp)
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1) == complex( 1._dp, 17._dp)
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2) == complex( 2._dp, 17._dp)
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3) == complex( 3._dp, 17._dp)
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_complex_2D

  subroutine test_distribute_dist_shared_from_primary_complex_3D( test_name_parent)
    !< Test the distribute_dist_shared_from_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                    :: routine_name = 'test_distribute_dist_shared_from_primary_complex_3D'
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

    ! Let the primary write some data to the memory
    if (par%primary) then
      allocate( d_tot( n_tot,nz,nl), source = complex( 0._dp, 0._dp))
      d_tot(13,1,2) = complex( 1._dp, 17._dp)
      d_tot(37,2,3) = complex( 2._dp, 17._dp)
      d_tot(72,3,5) = complex( 3._dp, 17._dp)
    end if

    ! Gather data to the primary
    if (par%primary) then
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl, d_tot = d_tot)
    else
      call distribute_dist_shared_from_primary( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
        n_tot, nz, nl)
    end if

    ! Evaluate test result and write to output file
    if (par%node_ID == 0) then
      test_result = d_nih( 13,1,2) == complex( 1._dp, 17._dp)
    elseif (par%node_ID == 1) then
      test_result = d_nih( 37,2,3) == complex( 2._dp, 17._dp)
    elseif (par%node_ID == 2) then
      test_result = d_nih( 72,3,5) == complex( 3._dp, 17._dp)
    end if
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)
    if (par%primary) deallocate( d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_dist_shared_from_primary_complex_3D

end module ut_mpi_distribute_dist_shared_from_primary
