module ut_mpi

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use mpi_distributed_memory, only: gather_to_primary, gather_to_all, distribute_from_primary
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, &
    gather_dist_shared_to_primary, gather_dist_shared_to_all
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD

  implicit none

  private

  public :: unit_tests_mpi_distributed_memory_main, unit_tests_mpi_hybrid_distributed_shared_memory_main

contains

  subroutine unit_tests_mpi_distributed_memory_main( test_name_parent)
    ! Run all unit tests for the MPI distributed memory subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_mpi_distributed_memory_main'
    character(len=1024), parameter :: test_name_local = 'mpi'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Run all unit tests for the MPI distributed memory subroutines
    call test_gather_to_primary( test_name)
    call test_gather_to_all( test_name)
    call test_distribute_from_primary( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_mpi_distributed_memory_main

  subroutine unit_tests_mpi_hybrid_distributed_shared_memory_main( test_name_parent)
    ! Run all unit tests for the MPI hybrid distributed/shared memory subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_mpi_hybrid_distributed_shared_memory_main'
    character(len=1024), parameter :: test_name_local = 'mpi'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! ! Run all unit tests for the MPI distributed memory subroutines
    call test_allocate_dist_shared( test_name)
    call test_gather_dist_shared_to_primary( test_name)
    call test_gather_dist_shared_to_all( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_mpi_hybrid_distributed_shared_memory_main

  ! ===== Gather to primary =====
  ! =============================

  subroutine test_gather_to_primary( test_name_parent)
    ! Test the gather_to_primary_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_gather_to_primary'
    character(len=1024), parameter :: test_name_local = 'gather_to_primary'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_gather_to_primary_int_1D( test_name)
    call test_gather_to_primary_int_2D( test_name)
    call test_gather_to_primary_dp_1D(  test_name)
    call test_gather_to_primary_dp_2D(  test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_primary

  subroutine test_gather_to_primary_int_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter       :: routine_name = 'test_gather_to_primary_int_1D'
    character(len=1024), parameter      :: test_name_local = 'int_1D'
    character(len=1024)                 :: test_name
    integer,  dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2))
      allocate( bb( 7))
      allocate( cc( 7))
    elseif (par%i == 1) then
      allocate( aa( 5))
      allocate( bb( 7))
      allocate( cc( 7))
    end if

    ! Fill test data
    cc = [1, 2, 3, 4, 5, 6, 7]

    if (par%i == 0) then
      aa = cc( 1:2)
    elseif (par%i == 1) then
      aa = cc( 3:7)
    end if

    ! Gather data
    call gather_to_primary( aa, bb)

    ! Check results
    if (par%primary) call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_primary_int_1D

  subroutine test_gather_to_primary_int_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter        :: routine_name = 'test_gather_to_primary_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    elseif (par%i == 1) then
      allocate( aa( 5,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    end if

    ! Fill test data
    cc(:,1) = [ 1,  2,  3,  4,  5,  6,  7]
    cc(:,2) = [ 8,  9, 10, 11, 12, 13, 14]

    if (par%i == 0) then
      aa = cc( 1:2,:)
    elseif (par%i == 1) then
      aa = cc( 3:7,:)
    end if

    ! Gather data
    call gather_to_primary( aa, bb)

    ! Check results
    if (par%primary) call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_primary_int_2D

  subroutine test_gather_to_primary_dp_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter       :: routine_name = 'test_gather_to_primary_dp_1D'
    character(len=1024), parameter      :: test_name_local = 'dp_1D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2))
      allocate( bb( 7))
      allocate( cc( 7))
    elseif (par%i == 1) then
      allocate( aa( 5))
      allocate( bb( 7))
      allocate( cc( 7))
    end if

    ! Fill test data
    cc = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp]

    if (par%i == 0) then
      aa = cc( 1:2)
    elseif (par%i == 1) then
      aa = cc( 3:7)
    end if

    ! Gather data
    call gather_to_primary( aa, bb)

    ! Check results
    if (par%primary) call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_primary_dp_1D

  subroutine test_gather_to_primary_dp_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter          :: routine_name = 'test_gather_to_primary_dp_2D'
    character(len=1024), parameter         :: test_name_local = 'dp_2D'
    character(len=1024)                    :: test_name
    real(dp),  dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    elseif (par%i == 1) then
      allocate( aa( 5,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    end if

    ! Fill test data
    cc(:,1) = [ 1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    cc(:,2) = [ 8._dp,  9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]

    if (par%i == 0) then
      aa = cc( 1:2,:)
    elseif (par%i == 1) then
      aa = cc( 3:7,:)
    end if

    ! Gather data
    call gather_to_primary( aa, bb)

    ! Check results
    if (par%primary) call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_primary_dp_2D

  ! ===== Gather to all =====
  ! =========================

  subroutine test_gather_to_all( test_name_parent)
    ! Test the gather_to_all_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_gather_to_all'
    character(len=1024), parameter :: test_name_local = 'gather_to_all'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_gather_to_all_int_1D( test_name)
    call test_gather_to_all_int_2D( test_name)
    call test_gather_to_all_dp_1D(  test_name)
    call test_gather_to_all_dp_2D(  test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_all

  subroutine test_gather_to_all_int_1D( test_name_parent)
    ! Test the gather_to_all_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_gather_to_all_int_1D'
    character(len=1024), parameter     :: test_name_local = 'int_1D'
    character(len=1024)                :: test_name
    integer, dimension(:), allocatable :: aa, bb, cc

    ! Add routine to path
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2))
      allocate( bb( 7))
      allocate( cc( 7))
    elseif (par%i == 1) then
      allocate( aa( 5))
      allocate( bb( 7))
      allocate( cc( 7))
    end if

    ! Fill test data
    cc = [1, 2, 3, 4, 5, 6, 7]

    if (par%i == 0) then
      aa = cc( 1:2)
    elseif (par%i == 1) then
      aa = cc( 3:7)
    end if

    ! Gather data
    call gather_to_all( aa, bb)

    ! Check results
    call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_all_int_1D

  subroutine test_gather_to_all_int_2D( test_name_parent)
    ! Test the gather_to_all_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'test_gather_to_all_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to path
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    elseif (par%i == 1) then
      allocate( aa( 5,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    end if

    ! Fill test data
    cc(:,1) = [1, 2,  3,  4,  5,  6,  7]
    cc(:,1) = [8, 9, 10, 11, 12, 13, 14]

    if (par%i == 0) then
      aa = cc( 1:2,:)
    elseif (par%i == 1) then
      aa = cc( 3:7,:)
    end if

    ! Gather data
    call gather_to_all( aa, bb)

    ! Check results
    call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_all_int_2D

  subroutine test_gather_to_all_dp_1D( test_name_parent)
    ! Test the gather_to_all_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_gather_to_all_dp_1D'
    character(len=1024), parameter      :: test_name_local = 'dp_1D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:), allocatable :: aa, bb, cc

    ! Add routine to path
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2))
      allocate( bb( 7))
      allocate( cc( 7))
    elseif (par%i == 1) then
      allocate( aa( 5))
      allocate( bb( 7))
      allocate( cc( 7))
    end if

    ! Fill test data
    cc = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp]

    if (par%i == 0) then
      aa = cc( 1:2)
    elseif (par%i == 1) then
      aa = cc( 3:7)
    end if

    ! Gather data
    call gather_to_all( aa, bb)

    ! Check results
    call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_all_dp_1D

  subroutine test_gather_to_all_dp_2D( test_name_parent)
    ! Test the gather_to_all_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_gather_to_all_dp_2D'
    character(len=1024), parameter        :: test_name_local = 'dp_2D'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to path
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 2,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    elseif (par%i == 1) then
      allocate( aa( 5,2))
      allocate( bb( 7,2))
      allocate( cc( 7,2))
    end if

    ! Fill test data
    cc(:,1) = [1._dp, 2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    cc(:,1) = [8._dp, 9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]

    if (par%i == 0) then
      aa = cc( 1:2,:)
    elseif (par%i == 1) then
      aa = cc( 3:7,:)
    end if

    ! Gather data
    call gather_to_all( aa, bb)

    ! Check results
    call unit_test( test_eq( bb, cc), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_all_dp_2D

  ! ===== Distribute from primary =====
  ! ===================================

  subroutine test_distribute_from_primary( test_name_parent)
    ! Test the distribute_from_primary_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_distribute_from_primary'
    character(len=1024), parameter :: test_name_local = 'distribute_from_primary'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_distribute_from_primary_int_1D( test_name)
    call test_distribute_from_primary_int_2D( test_name)
    call test_distribute_from_primary_dp_1D( test_name)
    call test_distribute_from_primary_dp_2D( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_primary

  subroutine test_distribute_from_primary_int_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_distribute_from_primary_int_1D'
    character(len=1024), parameter     :: test_name_local = 'int_1D'
    character(len=1024)                :: test_name
    integer, dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 7))
      allocate( bb( 2))
      allocate( cc( 7))
    elseif (par%i == 1) then
      allocate( aa( 7))
      allocate( bb( 5))
      allocate( cc( 7))
    end if

    ! Fill test data
    cc = [1, 2, 3, 4, 5, 6, 7]

    if (par%primary) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_primary( aa, bb)

    ! Check results
    if (par%primary) then
      call unit_test( test_eq(bb( 1:2), cc( 1:2)), test_name)
    elseif (par%i == 1) then
      call unit_test( test_eq(bb( 1:5), cc( 3:7)), test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_primary_int_1D

  subroutine test_distribute_from_primary_int_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'test_distribute_from_primary_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 7,2))
      allocate( bb( 2,2))
      allocate( cc( 7,2))
    elseif (par%i == 1) then
      allocate( aa( 7,2))
      allocate( bb( 5,2))
      allocate( cc( 7,2))
    end if

    ! Fill test data
    cc(:,1) = [1, 2, 3 , 4 , 5 , 6 , 7 ]
    cc(:,2) = [8, 9, 10, 11, 12, 13, 14]

    if (par%primary) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_primary( aa, bb)

    ! Check results
    if (par%primary) then
      call unit_test( test_eq(bb( 1:2,:), cc( 1:2,:)), test_name)
    elseif (par%i == 1) then
      call unit_test( test_eq(bb( 1:5,:), cc( 3:7,:)), test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_primary_int_2D

  subroutine test_distribute_from_primary_dp_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_distribute_from_primary_dp_1D'
    character(len=1024), parameter      :: test_name_local = 'dp_1D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 7))
      allocate( bb( 2))
      allocate( cc( 7))
    elseif (par%i == 1) then
      allocate( aa( 7))
      allocate( bb( 5))
      allocate( cc( 7))
    end if

    ! Fill test data
    cc = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp]

    if (par%primary) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_primary( aa, bb)

    ! Check results
    if (par%primary) then
      call unit_test( test_eq(bb( 1:2), cc( 1:2)), test_name)
    elseif (par%i == 1) then
      call unit_test( test_eq(bb( 1:5), cc( 3:7)), test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_primary_dp_1D

  subroutine test_distribute_from_primary_dp_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_distribute_from_primary_dp_2D'
    character(len=1024), parameter        :: test_name_local = 'dp_2D'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%i == 0) then
      allocate( aa( 7,2))
      allocate( bb( 2,2))
      allocate( cc( 7,2))
    elseif (par%i == 1) then
      allocate( aa( 7,2))
      allocate( bb( 5,2))
      allocate( cc( 7,2))
    end if

    ! Fill test data
    cc(:,1) = [1._dp, 2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    cc(:,2) = [8._dp, 9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]

    if (par%primary) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_primary( aa, bb)

    ! Check results
    if (par%primary) then
      call unit_test( test_eq(bb( 1:2,:), cc( 1:2,:)), test_name)
    elseif (par%i == 1) then
      call unit_test( test_eq(bb( 1:5,:), cc( 3:7,:)), test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_primary_dp_2D

  ! ===== (de)allocate_dist_shared =====
  ! ====================================

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

! ===== gather_dist_shared_to_primary =====
! =========================================

  subroutine test_gather_dist_shared_to_primary( test_name_parent)
    ! Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_gather_dist_shared_to_primary'
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
    logical, dimension(:), pointer     :: d => null()
    logical, dimension(:), allocatable :: d_tot
    type(MPI_WIN)                      :: w
    logical                            :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
      allocate( d_tot( 90), source = .false.)
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

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13) .and. d_tot( 34) .and. d_tot( 65)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_logical_1D

  subroutine test_gather_dist_shared_to_primary_logical_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'test_gather_dist_shared_to_primary_logical_2D'
    character(len=1024), parameter       :: test_name_local = 'logical_2D'
    character(len=1024)                  :: test_name
    logical, dimension(:,:), pointer     :: d => null()
    logical, dimension(:,:), allocatable :: d_tot
    type(MPI_WIN)                        :: w
    logical                              :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10)
      allocate( d_tot( 90,10), source = .false.)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = .true.
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3) .and. d_tot( 34,4) .and. d_tot( 65,5)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_logical_2D

  subroutine test_gather_dist_shared_to_primary_logical_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_gather_dist_shared_to_primary_logical_3D'
    character(len=1024), parameter         :: test_name_local = 'logical_3D'
    character(len=1024)                    :: test_name
    logical, dimension(:,:,:), pointer     :: d => null()
    logical, dimension(:,:,:), allocatable :: d_tot
    type(MPI_WIN)                          :: w
    logical                                :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10, 5)
      allocate( d_tot( 90,10,5), source = .false.)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10, 5)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = .true.
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3,2) .and. d_tot( 34,4,3) .and. d_tot( 65,5,4)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

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
    integer, dimension(:), pointer     :: d => null()
    integer, dimension(:), allocatable :: d_tot
    type(MPI_WIN)                      :: w
    logical                            :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
      allocate( d_tot( 90), source = 0)
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

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13) == 37 .and. d_tot( 34) == 38 .and. d_tot( 65) == 39
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_int_1D

  subroutine test_gather_dist_shared_to_primary_int_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'test_gather_dist_shared_to_primary_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:,:), pointer     :: d => null()
    integer, dimension(:,:), allocatable :: d_tot
    type(MPI_WIN)                        :: w
    logical                              :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10)
      allocate( d_tot( 90,10), source = 0)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = 39
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3) == 37 .and. d_tot( 34,4) == 38 .and. d_tot( 65,5) == 39
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_int_2D

  subroutine test_gather_dist_shared_to_primary_int_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_gather_dist_shared_to_primary_int_3D'
    character(len=1024), parameter         :: test_name_local = 'int_3D'
    character(len=1024)                    :: test_name
    integer, dimension(:,:,:), pointer     :: d => null()
    integer, dimension(:,:,:), allocatable :: d_tot
    type(MPI_WIN)                          :: w
    logical                                :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10, 5)
      allocate( d_tot( 90,10,5), source = 0)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10, 5)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = 39
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3,2) == 37 .and. d_tot( 34,4,3) == 38 .and. d_tot( 65,5,4) == 39
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

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
    real(dp), dimension(:), pointer     :: d => null()
    real(dp), dimension(:), allocatable :: d_tot
    type(MPI_WIN)                       :: w
    logical                             :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
      allocate( d_tot( 90), source = 0._dp)
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

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13) == 37._dp .and. d_tot( 34) == 38._dp .and. d_tot( 65) == 39._dp
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_dp_1D

  subroutine test_gather_dist_shared_to_primary_dp_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_gather_dist_shared_to_primary_dp_2D'
    character(len=1024), parameter        :: test_name_local = 'dp_2D'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), pointer     :: d => null()
    real(dp), dimension(:,:), allocatable :: d_tot
    type(MPI_WIN)                         :: w
    logical                               :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10)
      allocate( d_tot( 90,10), source = 0._dp)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = 39._dp
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3) == 37._dp .and. d_tot( 34,4) == 38._dp .and. d_tot( 65,5) == 39._dp
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_dp_2D

  subroutine test_gather_dist_shared_to_primary_dp_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'test_gather_dist_shared_to_primary_dp_3D'
    character(len=1024), parameter          :: test_name_local = 'dp_3D'
    character(len=1024)                     :: test_name
    real(dp), dimension(:,:,:), pointer     :: d => null()
    real(dp), dimension(:,:,:), allocatable :: d_tot
    type(MPI_WIN)                           :: w
    logical                                 :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10, 5)
      allocate( d_tot( 90,10,5), source = 0._dp)
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10, 5)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = 39._dp
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3,2) == 37._dp .and. d_tot( 34,4,3) == 38._dp .and. d_tot( 65,5,4) == 39._dp
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

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
    complex*16, dimension(:), pointer     :: d => null()
    complex*16, dimension(:), allocatable :: d_tot
    type(MPI_WIN)                         :: w
    logical                               :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20)
      allocate( d_tot( 90), source = complex( 0._dp, 0._dp))
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14) = complex( 14._dp, 38._dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15) = complex( 15._dp, 39._dp)
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13) == complex( 13._dp, 37._dp) .and. d_tot( 34) == complex( 14._dp, 38._dp) .and. d_tot( 65) == complex( 15._dp, 39._dp)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_complex_1D

  subroutine test_gather_dist_shared_to_primary_complex_2D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'test_gather_dist_shared_to_primary_complex_2D'
    character(len=1024), parameter          :: test_name_local = 'complex_2D'
    character(len=1024)                     :: test_name
    complex*16, dimension(:,:), pointer     :: d => null()
    complex*16, dimension(:,:), allocatable :: d_tot
    type(MPI_WIN)                           :: w
    logical                                 :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10)
      allocate( d_tot( 90,10), source = complex( 0._dp, 0._dp))
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = complex( 14._dp, 38._dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = complex( 15._dp, 39._dp)
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3) == complex( 13._dp, 37._dp) .and. d_tot( 34,4) == complex( 14._dp, 38._dp) .and. d_tot( 65,5) == complex( 15._dp, 39._dp)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_complex_2D

  subroutine test_gather_dist_shared_to_primary_complex_3D( test_name_parent)
    !< Test the gather_dist_shared_to_primary subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'test_gather_dist_shared_to_primary_complex_3D'
    character(len=1024), parameter            :: test_name_local = 'complex_3D'
    character(len=1024)                       :: test_name
    complex*16, dimension(:,:,:), pointer     :: d => null()
    complex*16, dimension(:,:,:), allocatable :: d_tot
    type(MPI_WIN)                             :: w
    logical                                   :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate different sizes on different nodes
    if (par%node_ID == 0) then
      call allocate_dist_shared( d, w, 20 ,10, 5)
      allocate( d_tot( 90,10,5), source = complex( 0._dp, 0._dp))
    elseif (par%node_ID == 1) then
      call allocate_dist_shared( d, w, 30 ,10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40 ,10, 5)
    end if

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = complex( 14._dp, 38._dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = complex( 15._dp, 39._dp)
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_primary( d, d_tot)
    if (par%node_ID == 0) then
      test_result = d_tot( 13,3,2) == complex( 13._dp, 37._dp) .and. d_tot( 34,4,3) == complex( 14._dp, 38._dp) .and. d_tot( 65,5,4) == complex( 15._dp, 39._dp)
    end if

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    if (par%node_ID == 0) then
      deallocate( d_tot)
    end if

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_primary_complex_3D

! ===== gather_dist_shared_to_primary =====
! =========================================

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
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_gather_dist_shared_to_all_logical_1D'
    character(len=1024), parameter :: test_name_local = 'logical_1D'
    character(len=1024)            :: test_name
    logical, dimension(:), pointer :: d => null()
    logical, dimension(:), pointer :: d_tot => null()
    type(MPI_WIN)                  :: w, w_tot
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
    call allocate_dist_shared( d_tot, w_tot, 90)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15) = .true.
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13) .and. d_tot( 34) .and. d_tot( 65)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_1D

  subroutine test_gather_dist_shared_to_all_logical_2D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'test_gather_dist_shared_to_all_logical_2D'
    character(len=1024), parameter   :: test_name_local = 'logical_2D'
    character(len=1024)              :: test_name
    logical, dimension(:,:), pointer :: d => null()
    logical, dimension(:,:), pointer :: d_tot => null()
    type(MPI_WIN)                    :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = .true.
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3) .and. d_tot( 34,4) .and. d_tot( 65,5)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_2D

  subroutine test_gather_dist_shared_to_all_logical_3D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_gather_dist_shared_to_all_logical_3D'
    character(len=1024), parameter     :: test_name_local = 'logical_3D'
    character(len=1024)                :: test_name
    logical, dimension(:,:,:), pointer :: d => null()
    logical, dimension(:,:,:), pointer :: d_tot => null()
    type(MPI_WIN)                      :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10, 5)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10, 5)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = .true.
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = .true.
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = .true.
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3,2) .and. d_tot( 34,4,3) .and. d_tot( 65,5,4)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_logical_3D

  subroutine test_gather_dist_shared_to_all_int_1D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_gather_dist_shared_to_all_int_1D'
    character(len=1024), parameter :: test_name_local = 'int_1D'
    character(len=1024)            :: test_name
    integer, dimension(:), pointer :: d => null()
    integer, dimension(:), pointer :: d_tot => null()
    type(MPI_WIN)                  :: w, w_tot
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
    call allocate_dist_shared( d_tot, w_tot, 90)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15) = 39
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13) == 37 .and. d_tot( 34) == 38 .and. d_tot( 65) == 39
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_1D

  subroutine test_gather_dist_shared_to_all_int_2D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'test_gather_dist_shared_to_all_int_2D'
    character(len=1024), parameter   :: test_name_local = 'int_2D'
    character(len=1024)              :: test_name
    integer, dimension(:,:), pointer :: d => null()
    integer, dimension(:,:), pointer :: d_tot => null()
    type(MPI_WIN)                    :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = 39
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3) == 37 .and. d_tot( 34,4) == 38 .and. d_tot( 65,5) == 39
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_2D

  subroutine test_gather_dist_shared_to_all_int_3D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_gather_dist_shared_to_all_int_3D'
    character(len=1024), parameter     :: test_name_local = 'int_3D'
    character(len=1024)                :: test_name
    integer, dimension(:,:,:), pointer :: d => null()
    integer, dimension(:,:,:), pointer :: d_tot => null()
    type(MPI_WIN)                      :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10, 5)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10, 5)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = 37
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = 38
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = 39
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3,2) == 37 .and. d_tot( 34,4,3) == 38 .and. d_tot( 65,5,4) == 39
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_int_3D

  subroutine test_gather_dist_shared_to_all_dp_1D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_gather_dist_shared_to_all_dp_1D'
    character(len=1024), parameter  :: test_name_local = 'dp_1D'
    character(len=1024)             :: test_name
    real(dp), dimension(:), pointer :: d => null()
    real(dp), dimension(:), pointer :: d_tot => null()
    type(MPI_WIN)                   :: w, w_tot
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
    call allocate_dist_shared( d_tot, w_tot, 90)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15) = 39._dp
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13) == 37._dp .and. d_tot( 34) == 38._dp .and. d_tot( 65) == 39._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_1D

  subroutine test_gather_dist_shared_to_all_dp_2D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'test_gather_dist_shared_to_all_dp_2D'
    character(len=1024), parameter    :: test_name_local = 'dp_2D'
    character(len=1024)               :: test_name
    real(dp), dimension(:,:), pointer :: d => null()
    real(dp), dimension(:,:), pointer :: d_tot => null()
    type(MPI_WIN)                     :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = 39._dp
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3) == 37._dp .and. d_tot( 34,4) == 38._dp .and. d_tot( 65,5) == 39._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_2D

  subroutine test_gather_dist_shared_to_all_dp_3D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_gather_dist_shared_to_all_dp_3D'
    character(len=1024), parameter      :: test_name_local = 'dp_3D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:,:,:), pointer :: d => null()
    real(dp), dimension(:,:,:), pointer :: d_tot => null()
    type(MPI_WIN)                       :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10, 5)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10, 5)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = 37._dp
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = 38._dp
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = 39._dp
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3,2) == 37._dp .and. d_tot( 34,4,3) == 38._dp .and. d_tot( 65,5,4) == 39._dp
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_dp_3D

  subroutine test_gather_dist_shared_to_all_complex_1D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'test_gather_dist_shared_to_all_complex_1D'
    character(len=1024), parameter    :: test_name_local = 'complex_1D'
    character(len=1024)               :: test_name
    complex*16, dimension(:), pointer :: d => null()
    complex*16, dimension(:), pointer :: d_tot => null()
    type(MPI_WIN)                     :: w, w_tot
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
    call allocate_dist_shared( d_tot, w_tot, 90)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14) = complex( 14._dp, 38._dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15) = complex( 15._dp, 39._dp)
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13) == complex( 13._dp, 37._dp) .and. d_tot( 34) == complex( 14._dp, 38._dp) .and. d_tot( 65) == complex( 15._dp, 39._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_1D

  subroutine test_gather_dist_shared_to_all_complex_2D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_gather_dist_shared_to_all_complex_2D'
    character(len=1024), parameter      :: test_name_local = 'complex_2D'
    character(len=1024)                 :: test_name
    complex*16, dimension(:,:), pointer :: d => null()
    complex*16, dimension(:,:), pointer :: d_tot => null()
    type(MPI_WIN)                       :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14,4) = complex( 14._dp, 38._dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15,5) = complex( 15._dp, 39._dp)
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3) == complex( 13._dp, 37._dp) .and. d_tot( 34,4) == complex( 14._dp, 38._dp) .and. d_tot( 65,5) == complex( 15._dp, 39._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_2D

  subroutine test_gather_dist_shared_to_all_complex_3D( test_name_parent)
    ! Test the gather_dist_shared_to_all subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_gather_dist_shared_to_all_complex_3D'
    character(len=1024), parameter        :: test_name_local = 'complex_3D'
    character(len=1024)                   :: test_name
    complex*16, dimension(:,:,:), pointer :: d => null()
    complex*16, dimension(:,:,:), pointer :: d_tot => null()
    type(MPI_WIN)                         :: w, w_tot
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
      call allocate_dist_shared( d, w, 30, 10, 5)
    elseif (par%node_ID == 2) then
      call allocate_dist_shared( d, w, 40, 10, 5)
    end if
    call allocate_dist_shared( d_tot, w_tot, 90, 10, 5)

    ! Let the node primaries write some data to the memory
    if (par%node_ID == 0 .and. par%node_primary) d(13,3,2) = complex( 13._dp, 37._dp)
    if (par%node_ID == 1 .and. par%node_primary) d(14,4,3) = complex( 14._dp, 38._dp)
    if (par%node_ID == 2 .and. par%node_primary) d(15,5,4) = complex( 15._dp, 39._dp)
    call sync

    ! Gather data to the primary
    call gather_dist_shared_to_all( d, d_tot)
    test_result = d_tot( 13,3,2) == complex( 13._dp, 37._dp) .and. d_tot( 34,4,3) == complex( 14._dp, 38._dp) .and. d_tot( 65,5,4) == complex( 15._dp, 39._dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    call deallocate_dist_shared( d, w)
    call deallocate_dist_shared( d_tot, w_tot)

    ! Evaluate test result and write to output file
    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_dist_shared_to_all_complex_3D

end module ut_mpi
