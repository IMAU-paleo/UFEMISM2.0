module ut_halo_exchange

  ! Unit tests for halo exchange code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, basic_halo_exchange

  implicit none

  private

  public :: test_halo_exchange_main

contains

subroutine test_halo_exchange_main( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'test_halo_exchange_main'
  character(len=1024), parameter :: test_name_local = 'halo_exchange'
  character(len=1024)            :: test_name

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  call test_halo_exchange_logical_1D( test_name)
  call test_halo_exchange_logical_2D( test_name)
  call test_halo_exchange_logical_3D( test_name)

  call test_halo_exchange_int_1D( test_name)
  call test_halo_exchange_int_2D( test_name)
  call test_halo_exchange_int_3D( test_name)

  call test_halo_exchange_dp_1D( test_name)
  call test_halo_exchange_dp_2D( test_name)
  call test_halo_exchange_dp_3D( test_name)

  call test_halo_exchange_complex_1D( test_name)
  call test_halo_exchange_complex_2D( test_name)
  call test_halo_exchange_complex_3D( test_name)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_main

subroutine test_halo_exchange_logical_1D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'test_halo_exchange_logical_1D'
  character(len=1024), parameter :: test_name_local = 'logical_1D'
  character(len=1024)            :: test_name
  logical, dimension(:), pointer :: d_nih => null()
  type(MPI_WIN)                  :: wd_nih
  integer                        :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                        :: i1_hle, i2_hle
  integer                        :: i1_hli, i2_hli
  integer                        :: i1_hre, i2_hre
  integer                        :: i1_hri, i2_hri
  integer                        :: i
  logical                        :: test_result
  integer                        :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  call allocate_dist_shared( d_nih, wd_nih, n_nih)
  d_nih( i1_nih:i2_nih) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = .false.
    do i = i1_node, i2_node
      d_nih( i) = .true.
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    test_result = test_result .and. d_nih( i)
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_logical_1D

subroutine test_halo_exchange_logical_2D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter               :: routine_name = 'test_halo_exchange_logical_2D'
  character(len=1024), parameter               :: test_name_local = 'logical_2D'
  character(len=1024)                          :: test_name
  logical, dimension(:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                :: wd_nih
  integer                                      :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                      :: i1_hle, i2_hle
  integer                                      :: i1_hli, i2_hli
  integer                                      :: i1_hre, i2_hre
  integer                                      :: i1_hri, i2_hri
  integer                                      :: nz, i, k
  logical                                      :: test_result
  integer                                      :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
  d_nih( i1_nih:i2_nih,1:nz) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = .false.
    do i = i1_node, i2_node
      do k = 1, nz
        d_nih( i,k) = .true.
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      test_result = test_result .and. d_nih( i,k)
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_logical_2D

subroutine test_halo_exchange_logical_3D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter                 :: routine_name = 'test_halo_exchange_logical_3D'
  character(len=1024), parameter                 :: test_name_local = 'logical_3D'
  character(len=1024)                            :: test_name
  logical, dimension(:,:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                  :: wd_nih
  integer                                        :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                        :: i1_hle, i2_hle
  integer                                        :: i1_hli, i2_hli
  integer                                        :: i1_hre, i2_hre
  integer                                        :: i1_hri, i2_hri
  integer                                        :: nz, nl, i, k, l
  logical                                        :: test_result
  integer                                        :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  nl = 3
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
  d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = .false.
    do i = i1_node, i2_node
      do k = 1, nz
        do l = 1, nl
          d_nih( i,k,l) = .true.
        end do
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz, nl)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      do l = 1, nl
        test_result = test_result .and. d_nih( i,k,l)
      end do
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_logical_3D

subroutine test_halo_exchange_int_1D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'test_halo_exchange_int_1D'
  character(len=1024), parameter :: test_name_local = 'int_1D'
  character(len=1024)            :: test_name
  integer, dimension(:), pointer :: d_nih => null()
  type(MPI_WIN)                  :: wd_nih
  integer                        :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                        :: i1_hle, i2_hle
  integer                        :: i1_hli, i2_hli
  integer                        :: i1_hre, i2_hre
  integer                        :: i1_hri, i2_hri
  integer                        :: i
  logical                        :: test_result
  integer                        :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  call allocate_dist_shared( d_nih, wd_nih, n_nih)
  d_nih( i1_nih:i2_nih) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      d_nih( i) = 2*i + 1
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    test_result = test_result .and. (d_nih( i) == (2*i + 1))
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_int_1D

subroutine test_halo_exchange_int_2D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter               :: routine_name = 'test_halo_exchange_int_2D'
  character(len=1024), parameter               :: test_name_local = 'int_2D'
  character(len=1024)                          :: test_name
  integer, dimension(:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                :: wd_nih
  integer                                      :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                      :: i1_hle, i2_hle
  integer                                      :: i1_hli, i2_hli
  integer                                      :: i1_hre, i2_hre
  integer                                      :: i1_hri, i2_hri
  integer                                      :: nz, i, k
  logical                                      :: test_result
  integer                                      :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
  d_nih( i1_nih:i2_nih,1:nz) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      do k = 1, nz
        d_nih( i,k) = (nz-1) * (2*i + 1) + k
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      test_result = test_result .and. (d_nih( i,k) == (nz-1) * (2*i + 1) + k)
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_int_2D

subroutine test_halo_exchange_int_3D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter                 :: routine_name = 'test_halo_exchange_int_3D'
  character(len=1024), parameter                 :: test_name_local = 'int_3D'
  character(len=1024)                            :: test_name
  integer, dimension(:,:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                  :: wd_nih
  integer                                        :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                        :: i1_hle, i2_hle
  integer                                        :: i1_hli, i2_hli
  integer                                        :: i1_hre, i2_hre
  integer                                        :: i1_hri, i2_hri
  integer                                        :: nz, nl, i, k, l
  logical                                        :: test_result
  integer                                        :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  nl = 3
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
  d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      do k = 1, nz
        do l = 1, nl
          d_nih( i,k,l) = (nl-1) * ((nz-1) * (2*i + 1) + k) + l
        end do
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz, nl)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      do l = 1, nl
        test_result = test_result .and. (d_nih( i,k,l) == (nl-1) * ((nz-1) * (2*i + 1) + k) + l)
      end do
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_int_3D

subroutine test_halo_exchange_dp_1D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter  :: routine_name = 'test_halo_exchange_dp_1D'
  character(len=1024), parameter  :: test_name_local = 'dp_1D'
  character(len=1024)             :: test_name
  real(dp), dimension(:), pointer :: d_nih => null()
  type(MPI_WIN)                   :: wd_nih
  integer                         :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                         :: i1_hle, i2_hle
  integer                         :: i1_hli, i2_hli
  integer                         :: i1_hre, i2_hre
  integer                         :: i1_hri, i2_hri
  integer                         :: i
  logical                         :: test_result
  integer                         :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  call allocate_dist_shared( d_nih, wd_nih, n_nih)
  d_nih( i1_nih:i2_nih) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      d_nih( i) = real( 2*i + 1, dp)
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    test_result = test_result .and. (d_nih( i) == real( 2*i + 1, dp))
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_dp_1D

subroutine test_halo_exchange_dp_2D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter                :: routine_name = 'test_halo_exchange_dp_2D'
  character(len=1024), parameter                :: test_name_local = 'dp_2D'
  character(len=1024)                           :: test_name
  real(dp), dimension(:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                 :: wd_nih
  integer                                       :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                       :: i1_hle, i2_hle
  integer                                       :: i1_hli, i2_hli
  integer                                       :: i1_hre, i2_hre
  integer                                       :: i1_hri, i2_hri
  integer                                       :: nz, i, k
  logical                                       :: test_result
  integer                                       :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
  d_nih( i1_nih:i2_nih,1:nz) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      do k = 1, nz
        d_nih( i,k) = real( (nz-1) * (2*i + 1) + k, dp)
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      test_result = test_result .and. (d_nih( i,k) == real( (nz-1) * (2*i + 1) + k, dp))
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_dp_2D

subroutine test_halo_exchange_dp_3D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter                  :: routine_name = 'test_halo_exchange_dp_3D'
  character(len=1024), parameter                  :: test_name_local = 'dp_3D'
  character(len=1024)                             :: test_name
  real(dp), dimension(:,:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                   :: wd_nih
  integer                                         :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                         :: i1_hle, i2_hle
  integer                                         :: i1_hli, i2_hli
  integer                                         :: i1_hre, i2_hre
  integer                                         :: i1_hri, i2_hri
  integer                                         :: nz, nl, i, k, l
  logical                                         :: test_result
  integer                                         :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  nl = 3
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
  d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      do k = 1, nz
        do l = 1, nl
          d_nih( i,k,l) = real( (nl-1) * ((nz-1) * (2*i + 1) + k) + l, dp)
        end do
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz, nl)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      do l = 1, nl
        test_result = test_result .and. (d_nih( i,k,l) == real( (nl-1) * ((nz-1) * (2*i + 1) + k) + l, dp))
      end do
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_dp_3D

subroutine test_halo_exchange_complex_1D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter    :: routine_name = 'test_halo_exchange_complex_1D'
  character(len=1024), parameter    :: test_name_local = 'complex_1D'
  character(len=1024)               :: test_name
  complex*16, dimension(:), pointer :: d_nih => null()
  type(MPI_WIN)                     :: wd_nih
  integer                           :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                           :: i1_hle, i2_hle
  integer                           :: i1_hli, i2_hli
  integer                           :: i1_hre, i2_hre
  integer                           :: i1_hri, i2_hri
  integer                           :: i
  logical                           :: test_result
  integer                           :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  call allocate_dist_shared( d_nih, wd_nih, n_nih)
  d_nih( i1_nih:i2_nih) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      d_nih( i) = complex( real( 2*i + 1, dp), 1._dp)
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    test_result = test_result .and. (d_nih( i) == complex( real( 2*i + 1, dp), 1._dp))
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_complex_1D

subroutine test_halo_exchange_complex_2D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter                  :: routine_name = 'test_halo_exchange_complex_2D'
  character(len=1024), parameter                  :: test_name_local = 'complex_2D'
  character(len=1024)                             :: test_name
  complex*16, dimension(:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                   :: wd_nih
  integer                                         :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                         :: i1_hle, i2_hle
  integer                                         :: i1_hli, i2_hli
  integer                                         :: i1_hre, i2_hre
  integer                                         :: i1_hri, i2_hri
  integer                                         :: nz, i, k
  logical                                         :: test_result
  integer                                         :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz)
  d_nih( i1_nih:i2_nih,1:nz) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      do k = 1, nz
        d_nih( i,k) = complex( real( (nz-1) * (2*i + 1) + k, dp), 1._dp)
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      test_result = test_result .and. (d_nih( i,k) == &
        complex( real( (nz-1) * (2*i + 1) + k, dp), 1._dp))
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_complex_2D

subroutine test_halo_exchange_complex_3D( test_name_parent)
  ! Test the halo exchange subroutines

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter                    :: routine_name = 'test_halo_exchange_complex_3D'
  character(len=1024), parameter                    :: test_name_local = 'complex_3D'
  character(len=1024)                               :: test_name
  complex*16, dimension(:,:,:), pointer, contiguous :: d_nih => null()
  type(MPI_WIN)                                     :: wd_nih
  integer                                           :: i1_nih, i2_nih, i1_node, i2_node, n_nih
  integer                                           :: i1_hle, i2_hle
  integer                                           :: i1_hli, i2_hli
  integer                                           :: i1_hre, i2_hre
  integer                                           :: i1_hri, i2_hri
  integer                                           :: nz, nl, i, k, l
  logical                                           :: test_result
  integer                                           :: ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Safety - should be run on 7 cores
  call assert( test_eq( par%n, 7), 'should be run on 7 cores')

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  if (par%node_ID == 0) then
    ! Node 0: 20 elements
    i1_node = 1
    i2_node = 20
  elseif (par%node_ID == 1) then
    ! Node 1: 30 elements
    i1_node = 21
    i2_node = 50
  elseif (par%node_ID == 2) then
    ! Node 2: 40 elements
    i1_node = 51
    i2_node = 90
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

  n_nih = i2_nih + 1 - i1_nih

  ! Allocate some node-shared memory
  nz = 2
  nl = 3
  call allocate_dist_shared( d_nih, wd_nih, n_nih, nz, nl)
  d_nih( i1_nih:i2_nih,1:nz,1:nl) => d_nih

  ! Fill in data
  if (par%node_primary) then
    d_nih = -1._dp
    do i = i1_node, i2_node
      do k = 1, nz
        do l = 1, nl
          d_nih( i,k,l) = complex( real( (nl-1) * ((nz-1) * (2*i + 1) + k) + l, dp), 1._dp)
        end do
      end do
    end do
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, nz, nl)

  ! Verify results
  test_result = .true.
  do i = i1_nih, i2_nih
    do k = 1, nz
      do l = 1, nl
        test_result = test_result .and. (d_nih( i,k,l) == &
          complex( real( (nl-1) * ((nz-1) * (2*i + 1) + k) + l, dp), 1._dp))
      end do
    end do
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  call unit_test( test_result, test_name)

  ! Clean up after yourself
  call deallocate_dist_shared( d_nih, wd_nih)

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine test_halo_exchange_complex_3D

end module ut_halo_exchange
