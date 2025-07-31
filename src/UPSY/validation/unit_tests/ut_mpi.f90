module ut_mpi

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use mpi_distributed_memory, only: gather_to_primary, gather_to_all, distribute_from_primary
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, &
    gather_dist_shared_to_primary, gather_dist_shared_to_all, distribute_dist_shared_from_primary
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD

  implicit none

  private

  public :: unit_tests_mpi_distributed_memory_main

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

end module ut_mpi
