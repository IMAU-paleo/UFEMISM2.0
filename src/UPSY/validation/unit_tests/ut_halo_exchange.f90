module ut_halo_exchange

  ! Unit tests for halo exchange code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use mpi_distributed_shared_memory
  use ut_mpi_allocate_dist_shared, only: setup_simple_parallel_array_info

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
  character(len=1024), parameter             :: routine_name = 'test_halo_exchange_logical_1D'
  character(len=1024), parameter             :: test_name_local = 'logical_1D'
  character(len=1024)                        :: test_name
  type(type_par_arr_info)                    :: pai
  logical, dimension(:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                              :: wd_nih
  integer                                    :: i
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = .false.
    d_nih( pai%i1_node:pai%i2_node) = .true.
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
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
  type(type_par_arr_info)                      :: pai
  integer                                      :: nz
  logical, dimension(:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                :: wd_nih
  integer                                      :: i,k
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = .false.
    d_nih( pai%i1_node:pai%i2_node,:) = .true.
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
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
  type(type_par_arr_info)                        :: pai
  integer                                        :: nz, nl
  logical, dimension(:,:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                  :: wd_nih
  integer                                        :: i,k,l
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = .false.
    d_nih( pai%i1_node:pai%i2_node,:,:) = .true.
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, nl, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    do k = 1, nz
      do l = 1, nz
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
  character(len=1024), parameter             :: routine_name = 'test_halo_exchange_int_1D'
  character(len=1024), parameter             :: test_name_local = 'int_1D'
  character(len=1024)                        :: test_name
  type(type_par_arr_info)                    :: pai
  integer, dimension(:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                              :: wd_nih
  integer                                    :: i
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = 0
    d_nih( pai%i1_node:pai%i2_node) = 1
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    test_result = test_result .and. d_nih( i) == 1
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
  type(type_par_arr_info)                      :: pai
  integer                                      :: nz
  integer, dimension(:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                :: wd_nih
  integer                                      :: i,k
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = 0
    d_nih( pai%i1_node:pai%i2_node,:) = 1
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    do k = 1, nz
      test_result = test_result .and. d_nih( i,k) == 1
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
  type(type_par_arr_info)                        :: pai
  integer                                        :: nz, nl
  integer, dimension(:,:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                  :: wd_nih
  integer                                        :: i,k,l
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = 0
    d_nih( pai%i1_node:pai%i2_node,:,:) = 1
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, nl, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    do k = 1, nz
      do l = 1, nz
        test_result = test_result .and. d_nih( i,k,l) == 1
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
  character(len=1024), parameter             :: routine_name = 'test_halo_exchange_dp_1D'
  character(len=1024), parameter             :: test_name_local = 'dp_1D'
  character(len=1024)                        :: test_name
  type(type_par_arr_info)                    :: pai
  real(dp), dimension(:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                              :: wd_nih
  integer                                    :: i
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = 0._dp
    d_nih( pai%i1_node:pai%i2_node) = 1._dp
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    test_result = test_result .and. d_nih( i) == 1._dp
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
  character(len=1024), parameter               :: routine_name = 'test_halo_exchange_dp_2D'
  character(len=1024), parameter               :: test_name_local = 'dp_2D'
  character(len=1024)                          :: test_name
  type(type_par_arr_info)                      :: pai
  integer                                      :: nz
  real(dp), dimension(:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                :: wd_nih
  integer                                      :: i,k
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = 0._dp
    d_nih( pai%i1_node:pai%i2_node,:) = 1._dp
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    do k = 1, nz
      test_result = test_result .and. d_nih( i,k) == 1._dp
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
  character(len=1024), parameter                 :: routine_name = 'test_halo_exchange_dp_3D'
  character(len=1024), parameter                 :: test_name_local = 'dp_3D'
  character(len=1024)                            :: test_name
  type(type_par_arr_info)                        :: pai
  integer                                        :: nz, nl
  real(dp), dimension(:,:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                  :: wd_nih
  integer                                        :: i,k,l
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = 0._dp
    d_nih( pai%i1_node:pai%i2_node,:,:) = 1._dp
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, nl, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    do k = 1, nz
      do l = 1, nz
        test_result = test_result .and. d_nih( i,k,l) == 1._dp
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
  character(len=1024), parameter             :: routine_name = 'test_halo_exchange_complex_1D'
  character(len=1024), parameter             :: test_name_local = 'complex_1D'
  character(len=1024)                        :: test_name
  type(type_par_arr_info)                    :: pai
  complex*16, dimension(:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                              :: wd_nih
  integer                                    :: i
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = complex( 0._dp, 0._dp)
    d_nih( pai%i1_node:pai%i2_node) = complex( 13._dp, 37._dp)
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    test_result = test_result .and. d_nih( i) == complex( 13._dp, 37._dp)
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
  character(len=1024), parameter               :: routine_name = 'test_halo_exchange_complex_2D'
  character(len=1024), parameter               :: test_name_local = 'complex_2D'
  character(len=1024)                          :: test_name
  type(type_par_arr_info)                      :: pai
  integer                                      :: nz
  complex*16, dimension(:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                :: wd_nih
  integer                                      :: i,k
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = complex( 0._dp, 0._dp)
    d_nih( pai%i1_node:pai%i2_node,:) = complex( 13._dp, 37._dp)
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    do k = 1, nz
      test_result = test_result .and. d_nih( i,k) == complex( 13._dp, 37._dp)
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
  character(len=1024), parameter                 :: routine_name = 'test_halo_exchange_complex_3D'
  character(len=1024), parameter                 :: test_name_local = 'complex_3D'
  character(len=1024)                            :: test_name
  type(type_par_arr_info)                        :: pai
  integer                                        :: nz, nl
  complex*16, dimension(:,:,:), contiguous, pointer :: d_nih => null()
  type(MPI_WIN)                                  :: wd_nih
  integer                                        :: i,k,l
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

  ! Fill in data in the interior only
  if (par%node_primary) then
    d_nih = complex( 0._dp, 0._dp)
    d_nih( pai%i1_node:pai%i2_node,:,:) = complex( 13._dp, 37._dp)
  end if
  call sync_node

  ! Exchange halos
  call basic_halo_exchange( pai, nz, nl, d_nih)

  ! Verify results
  test_result = .true.
  do i = pai%i1_nih, pai%i2_nih
    do k = 1, nz
      do l = 1, nz
        test_result = test_result .and. d_nih( i,k,l) == complex( 13._dp, 37._dp)
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
