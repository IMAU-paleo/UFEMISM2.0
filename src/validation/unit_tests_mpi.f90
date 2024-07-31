module unit_tests_mpi

  ! Unit tests for different MPI routines

  use mpi
  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging                  , only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use mpi_distributed_memory                                 , only: gather_to_all_int_1D, gather_to_all_int_2D, gather_to_all_dp_1D, gather_to_all_dp_2D, &
                                                                     gather_to_master_int_1D, gather_to_master_int_2D, gather_to_master_dp_1D, &
                                                                     gather_to_master_dp_2D, distribute_from_master_int_1D, distribute_from_master_int_2D, &
                                                                     distribute_from_master_dp_1D, distribute_from_master_dp_2D
  use assertions_unit_tests, only: ASSERTION, UNIT_TEST, test_eqv, test_neqv, test_eq, test_neq, test_gt, test_lt, test_ge, test_le, test_ge_le, test_tol

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
    character(len=1024), parameter :: test_name_local = 'mpi_distributed_memory'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Run all unit tests for the MPI distributed memory subroutines
    call test_gather_to_master( test_name)
    call test_gather_to_all( test_name)
    call test_distribute_from_master( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_mpi_distributed_memory_main

  ! ===== Gather to master =====
  ! ============================

  subroutine test_gather_to_master( test_name_parent)
    ! Test the gather_to_master_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_gather_to_master'
    character(len=1024), parameter :: test_name_local = 'gather_to_master'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_gather_to_master_int_1D( test_name)
    call test_gather_to_master_int_2D( test_name)
    call test_gather_to_master_dp_1D(  test_name)
    call test_gather_to_master_dp_2D(  test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_master

  subroutine test_gather_to_master_int_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter       :: routine_name = 'test_gather_to_master_int_1D'
    character(len=1024), parameter      :: test_name_local = 'int_1D'
    character(len=1024)                 :: test_name
    integer,  dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    call gather_to_master_int_1D( aa, bb)

    ! Check results
    if (par%master) call test_eq( bb, cc, UNIT_TEST, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_master_int_1D

  subroutine test_gather_to_master_int_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter        :: routine_name = 'test_gather_to_master_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    CALL gather_to_master_int_2D( aa, bb)

    bb(1,1) = 0

    ! Check results
    if (par%master) call test_eq( bb, cc, UNIT_TEST, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_master_int_2D

  subroutine test_gather_to_master_dp_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter       :: routine_name = 'test_gather_to_master_dp_1D'
    character(len=1024), parameter      :: test_name_local = 'dp_1D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    call gather_to_master_dp_1D( aa, bb)

    ! Check results
    if (par%master) call test_eq( bb, cc, UNIT_TEST, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_master_dp_1D

  subroutine test_gather_to_master_dp_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024),parameter          :: routine_name = 'test_gather_to_master_dp_2D'
    character(len=1024), parameter         :: test_name_local = 'dp_2D'
    character(len=1024)                    :: test_name
    real(dp),  dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    CALL gather_to_master_dp_2D( aa, bb)

    ! Check results
    if (par%master) call test_eq( bb, cc, UNIT_TEST, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_master_dp_2D

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
    CALL init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    CALL gather_to_all_int_1D( aa, bb)

    ! Check results
    call test_eq( bb, cc, UNIT_TEST, test_name)

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
    CALL init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    CALL gather_to_all_int_2D( aa, bb)

    ! Check results
    call test_eq( bb, cc, UNIT_TEST, test_name)

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
    CALL init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    CALL gather_to_all_dp_1D( aa, bb)

    ! Check results
    call test_eq( bb, cc, UNIT_TEST, test_name)

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
    CALL init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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
    CALL gather_to_all_dp_2D( aa, bb)

    ! Check results
    call test_eq( bb, cc, UNIT_TEST, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_gather_to_all_dp_2D

  ! ===== Distribute from master =====
  ! ==================================

  subroutine test_distribute_from_master( test_name_parent)
    ! Test the distribute_from_master_TYPE_DIM subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_distribute_from_master'
    character(len=1024), parameter :: test_name_local = 'distribute_from_master'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_distribute_from_master_int_1D( test_name)
    call test_distribute_from_master_int_2D( test_name)
    call test_distribute_from_master_dp_1D( test_name)
    call test_distribute_from_master_dp_2D( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_master

  subroutine test_distribute_from_master_int_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_distribute_from_master_int_1D'
    character(len=1024), parameter     :: test_name_local = 'int_1D'
    character(len=1024)                :: test_name
    integer, dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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

    if (par%master) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_master_int_1D( aa, bb)

    ! Check results
    if (par%master) then
      call test_eq(bb( 1:2), cc( 1:2), UNIT_TEST, test_name)
    elseif (par%i == 1) then
      call test_eq(bb( 1:5), cc( 3:7), UNIT_TEST, test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_master_int_1D

  subroutine test_distribute_from_master_int_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'test_distribute_from_master_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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

    if (par%master) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_master_int_2D( aa, bb)

    ! Check results
    if (par%master) then
      call test_eq(bb( 1:2,:), cc( 1:2,:), UNIT_TEST, test_name)
    elseif (par%i == 1) then
      call test_eq(bb( 1:5,:), cc( 3:7,:), UNIT_TEST, test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_master_int_2D

  subroutine test_distribute_from_master_dp_1D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_distribute_from_master_dp_1D'
    character(len=1024), parameter      :: test_name_local = 'dp_1D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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

    if (par%master) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_master_dp_1D( aa, bb)

    ! Check results
    if (par%master) then
      call test_eq(bb( 1:2), cc( 1:2), UNIT_TEST, test_name)
    elseif (par%i == 1) then
      call test_eq(bb( 1:5), cc( 3:7), UNIT_TEST, test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_master_dp_1D

  subroutine test_distribute_from_master_dp_2D( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_distribute_from_master_dp_2D'
    character(len=1024), parameter        :: test_name_local = 'dp_2D'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), allocatable :: aa, bb, cc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Safety - should be run on two cores
    call test_eq( par%n, 2, ASSERTION, 'should be run on two cores')

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

    if (par%master) then
      aa = cc
    end if

    ! Distribute data
    call distribute_from_master_dp_2D( aa, bb)

    ! Check results
    if (par%master) then
      call test_eq(bb( 1:2,:), cc( 1:2,:), UNIT_TEST, test_name)
    elseif (par%i == 1) then
      call test_eq(bb( 1:5,:), cc( 3:7,:), UNIT_TEST, test_name)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_distribute_from_master_dp_2D

end module unit_tests_mpi
