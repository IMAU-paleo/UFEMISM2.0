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
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD

  implicit none

  private

  public :: test_gather_dist_shared_to_primary

contains

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

end module ut_mpi_gather_dist_shared_to_primary
