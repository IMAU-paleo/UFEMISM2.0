module ut_mpi_multiply_CSR_matrix_vector

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, MPI_WIN
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared

  implicit none

  private

  public :: test_multiply_CSR_matrix_with_vector

contains

  subroutine test_multiply_CSR_matrix_with_vector( test_name_parent)
    ! Test the multiply_CSR_matrix_with_vector subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_multiply_CSR_matrix_with_vector'
    character(len=1024), parameter :: test_name_local = 'multiply_CSR_matrix_with_vector'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_multiply_CSR_matrix_with_vector_1D_dist_dist    ( test_name)
    call test_multiply_CSR_matrix_with_vector_1D_dist_hybrid  ( test_name)
    call test_multiply_CSR_matrix_with_vector_1D_hybrid_dist  ( test_name)
    call test_multiply_CSR_matrix_with_vector_1D_hybrid_hybrid( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_CSR_matrix_with_vector

  subroutine test_multiply_CSR_matrix_with_vector_1D_dist_dist( test_name_parent)
    ! Test the multiply_CSR_matrix_with_vector subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_multiply_CSR_matrix_with_vector_1D_dist_dist'
    character(len=1024), parameter      :: test_name_local = '1D/dist_dist'
    character(len=1024)                 :: test_name
    type(type_sparse_matrix_CSR_dp)     :: AA
    real(dp), dimension(:), allocatable :: xx, yy, yy_correct
    logical                             :: test_result
    integer                             :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Verify that we can solve the following matrix equation:
    !
    !              A             x   =    y
    !
    !    [1,  ,  ,  ,  ,  ,  ]  [1]     [ 1]
    !    [ , 2,  ,  ,  ,  ,  ]  [2]     [ 4]
    !    [ , 3,  ,  , 4,  ,  ]  [3]     [26]
    !    [ ,  , 2,  ,  , 1,  ]  [4]  =  [12]
    !    [ ,  ,  , 5,  ,  , 2]  [5]     [34]
    !    [ , 1,  ,  , 3,  ,  ]  [6]     [17]
    !    [ ,  , 7,  ,  ,  , 1]  [7]     [28]

    call allocate_matrix_CSR_dist( AA, 7, 7, 1, 1, 2)
    allocate( xx        ( AA%j1:AA%j2))
    allocate( yy        ( AA%i1:AA%i2))
    allocate( yy_correct( AA%i1:AA%i2))

    ! Define A, x, and y_correct

    if     (par%i == 0) then
      call add_entry_CSR_dist( AA, 1, 1, 1._dp)
      xx        (1) = 1
      yy_correct(1) = 1
    elseif (par%i == 1) then
      call add_entry_CSR_dist( AA, 2, 2, 2._dp)
      xx        (2) = 2
      yy_correct(2) = 4
    elseif (par%i == 2) then
      call add_entry_CSR_dist( AA, 3, 2, 3._dp)
      call add_entry_CSR_dist( AA, 3, 5, 4._dp)
      xx        (3) = 3
      yy_correct(3) = 26
    elseif (par%i == 3) then
      call add_entry_CSR_dist( AA, 4, 3, 2._dp)
      call add_entry_CSR_dist( AA, 4, 6, 1._dp)
      xx        (4) = 4
      yy_correct(4) = 12
    elseif (par%i == 4) then
      call add_entry_CSR_dist( AA, 5, 4, 5._dp)
      call add_entry_CSR_dist( AA, 5, 7, 2._dp)
      xx        (5) = 5
      yy_correct(5) = 34
    elseif (par%i == 5) then
      call add_entry_CSR_dist( AA, 6, 2, 1._dp)
      call add_entry_CSR_dist( AA, 6, 5, 3._dp)
      xx        (6) = 6
      yy_correct(6) = 17
    elseif (par%i == 6) then
      call add_entry_CSR_dist( AA, 7, 3, 7._dp)
      call add_entry_CSR_dist( AA, 7, 7, 1._dp)
      xx        (7) = 7
      yy_correct(7) = 28
    end if

    call multiply_CSR_matrix_with_vector_1D( AA, xx, yy, xx_is_hybrid = .false., yy_is_hybrid = .false.)

    test_result = test_eq( yy( AA%i1:AA%i2), yy_correct( AA%i1:AA%i2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    call unit_test( test_result, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_CSR_matrix_with_vector_1D_dist_dist

  subroutine test_multiply_CSR_matrix_with_vector_1D_dist_hybrid( test_name_parent)
    ! Test the multiply_CSR_matrix_with_vector subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_multiply_CSR_matrix_with_vector_1D_dist_hybrid'
    character(len=1024), parameter      :: test_name_local = '1D/dist_hybrid'
    character(len=1024)                 :: test_name
    type(type_sparse_matrix_CSR_dp)     :: AA
    real(dp), dimension(:), allocatable :: xx
    real(dp), dimension(:), pointer     :: yy => null()
    real(dp), dimension(:), pointer     :: yy_correct => null()
    type(MPI_WIN)                       :: wyy, wyy_correct
    logical                             :: test_result
    integer                             :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Verify that we can solve the following matrix equation:
    !
    !              A             x   =    y
    !
    !    [1,  ,  ,  ,  ,  ,  ]  [1]     [ 1]
    !    [ , 2,  ,  ,  ,  ,  ]  [2]     [ 4]
    !    [ , 3,  ,  , 4,  ,  ]  [3]     [26]
    !    [ ,  , 2,  ,  , 1,  ]  [4]  =  [12]
    !    [ ,  ,  , 5,  ,  , 2]  [5]     [34]
    !    [ , 1,  ,  , 3,  ,  ]  [6]     [17]
    !    [ ,  , 7,  ,  ,  , 1]  [7]     [28]

    call allocate_matrix_CSR_dist( AA, 7, 7, 1, 1, 2)
    allocate( xx        ( AA%j1:AA%j2))

    call allocate_dist_shared( yy, wyy, AA%m_node)
    yy( AA%i1_node:AA%i2_node) => yy
    call allocate_dist_shared( yy_correct, wyy_correct, AA%m_node)
    yy_correct( AA%i1_node:AA%i2_node) => yy_correct

    ! Define A, x, and y_correct

    if     (par%i == 0) then
      call add_entry_CSR_dist( AA, 1, 1, 1._dp)
      xx        (1) = 1
      yy_correct(1) = 1
    elseif (par%i == 1) then
      call add_entry_CSR_dist( AA, 2, 2, 2._dp)
      xx        (2) = 2
      yy_correct(2) = 4
    elseif (par%i == 2) then
      call add_entry_CSR_dist( AA, 3, 2, 3._dp)
      call add_entry_CSR_dist( AA, 3, 5, 4._dp)
      xx        (3) = 3
      yy_correct(3) = 26
    elseif (par%i == 3) then
      call add_entry_CSR_dist( AA, 4, 3, 2._dp)
      call add_entry_CSR_dist( AA, 4, 6, 1._dp)
      xx        (4) = 4
      yy_correct(4) = 12
    elseif (par%i == 4) then
      call add_entry_CSR_dist( AA, 5, 4, 5._dp)
      call add_entry_CSR_dist( AA, 5, 7, 2._dp)
      xx        (5) = 5
      yy_correct(5) = 34
    elseif (par%i == 5) then
      call add_entry_CSR_dist( AA, 6, 2, 1._dp)
      call add_entry_CSR_dist( AA, 6, 5, 3._dp)
      xx        (6) = 6
      yy_correct(6) = 17
    elseif (par%i == 6) then
      call add_entry_CSR_dist( AA, 7, 3, 7._dp)
      call add_entry_CSR_dist( AA, 7, 7, 1._dp)
      xx        (7) = 7
      yy_correct(7) = 28
    end if

    call multiply_CSR_matrix_with_vector_1D( AA, xx, yy, xx_is_hybrid = .false., yy_is_hybrid = .true.)

    test_result = test_eq( yy( AA%i1:AA%i2), yy_correct( AA%i1:AA%i2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( yy, wyy)
    call deallocate_dist_shared( yy_correct, wyy_correct)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_CSR_matrix_with_vector_1D_dist_hybrid

  subroutine test_multiply_CSR_matrix_with_vector_1D_hybrid_dist( test_name_parent)
    ! Test the multiply_CSR_matrix_with_vector subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_multiply_CSR_matrix_with_vector_1D_hybrid_dist'
    character(len=1024), parameter      :: test_name_local = '1D/hybrid_dist'
    character(len=1024)                 :: test_name
    type(type_sparse_matrix_CSR_dp)     :: AA
    real(dp), dimension(:), pointer     :: xx => null()
    type(MPI_WIN)                       :: wxx
    real(dp), dimension(:), allocatable :: yy, yy_correct
    logical                             :: test_result
    integer                             :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Verify that we can solve the following matrix equation:
    !
    !              A             x   =    y
    !
    !    [1,  ,  ,  ,  ,  ,  ]  [1]     [ 1]
    !    [ , 2,  ,  ,  ,  ,  ]  [2]     [ 4]
    !    [ , 3,  ,  , 4,  ,  ]  [3]     [26]
    !    [ ,  , 2,  ,  , 1,  ]  [4]  =  [12]
    !    [ ,  ,  , 5,  ,  , 2]  [5]     [34]
    !    [ , 1,  ,  , 3,  ,  ]  [6]     [17]
    !    [ ,  , 7,  ,  ,  , 1]  [7]     [28]

    call allocate_matrix_CSR_dist( AA, 7, 7, 1, 1, 2)
    call allocate_dist_shared( xx, wxx, AA%n_node)
    xx( AA%j1_node:AA%j2_node) => xx
    allocate( yy        ( AA%i1:AA%i2))
    allocate( yy_correct( AA%i1:AA%i2))

    ! Define A, x, and y_correct

    if     (par%i == 0) then
      call add_entry_CSR_dist( AA, 1, 1, 1._dp)
      xx        (1) = 1
      yy_correct(1) = 1
    elseif (par%i == 1) then
      call add_entry_CSR_dist( AA, 2, 2, 2._dp)
      xx        (2) = 2
      yy_correct(2) = 4
    elseif (par%i == 2) then
      call add_entry_CSR_dist( AA, 3, 2, 3._dp)
      call add_entry_CSR_dist( AA, 3, 5, 4._dp)
      xx        (3) = 3
      yy_correct(3) = 26
    elseif (par%i == 3) then
      call add_entry_CSR_dist( AA, 4, 3, 2._dp)
      call add_entry_CSR_dist( AA, 4, 6, 1._dp)
      xx        (4) = 4
      yy_correct(4) = 12
    elseif (par%i == 4) then
      call add_entry_CSR_dist( AA, 5, 4, 5._dp)
      call add_entry_CSR_dist( AA, 5, 7, 2._dp)
      xx        (5) = 5
      yy_correct(5) = 34
    elseif (par%i == 5) then
      call add_entry_CSR_dist( AA, 6, 2, 1._dp)
      call add_entry_CSR_dist( AA, 6, 5, 3._dp)
      xx        (6) = 6
      yy_correct(6) = 17
    elseif (par%i == 6) then
      call add_entry_CSR_dist( AA, 7, 3, 7._dp)
      call add_entry_CSR_dist( AA, 7, 7, 1._dp)
      xx        (7) = 7
      yy_correct(7) = 28
    end if

    call multiply_CSR_matrix_with_vector_1D( AA, xx, yy, xx_is_hybrid = .true., yy_is_hybrid = .false.)

    test_result = test_eq( yy( AA%i1:AA%i2), yy_correct( AA%i1:AA%i2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( xx, wxx)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_CSR_matrix_with_vector_1D_hybrid_dist

  subroutine test_multiply_CSR_matrix_with_vector_1D_hybrid_hybrid( test_name_parent)
    ! Test the multiply_CSR_matrix_with_vector subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_multiply_CSR_matrix_with_vector_1D_hybrid_hybrid'
    character(len=1024), parameter      :: test_name_local = '1D/hybrid_hybrid'
    character(len=1024)                 :: test_name
    type(type_sparse_matrix_CSR_dp)     :: AA
    real(dp), dimension(:), pointer     :: xx => null()
    real(dp), dimension(:), pointer     :: yy => null()
    real(dp), dimension(:), pointer     :: yy_correct => null()
    type(MPI_WIN)                       :: wxx, wyy, wyy_correct
    logical                             :: test_result
    integer                             :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Verify that we can solve the following matrix equation:
    !
    !              A             x   =    y
    !
    !    [1,  ,  ,  ,  ,  ,  ]  [1]     [ 1]
    !    [ , 2,  ,  ,  ,  ,  ]  [2]     [ 4]
    !    [ , 3,  ,  , 4,  ,  ]  [3]     [26]
    !    [ ,  , 2,  ,  , 1,  ]  [4]  =  [12]
    !    [ ,  ,  , 5,  ,  , 2]  [5]     [34]
    !    [ , 1,  ,  , 3,  ,  ]  [6]     [17]
    !    [ ,  , 7,  ,  ,  , 1]  [7]     [28]

    call allocate_matrix_CSR_dist( AA, 7, 7, 1, 1, 2)

    call allocate_dist_shared( xx, wxx, AA%n_node)
    xx( AA%j1_node:AA%j2_node) => xx
    call allocate_dist_shared( yy, wyy, AA%m_node)
    yy( AA%i1_node:AA%i2_node) => yy
    call allocate_dist_shared( yy_correct, wyy_correct, AA%m_node)
    yy_correct( AA%i1_node:AA%i2_node) => yy_correct

    ! Define A, x, and y_correct

    if     (par%i == 0) then
      call add_entry_CSR_dist( AA, 1, 1, 1._dp)
      xx        (1) = 1
      yy_correct(1) = 1
    elseif (par%i == 1) then
      call add_entry_CSR_dist( AA, 2, 2, 2._dp)
      xx        (2) = 2
      yy_correct(2) = 4
    elseif (par%i == 2) then
      call add_entry_CSR_dist( AA, 3, 2, 3._dp)
      call add_entry_CSR_dist( AA, 3, 5, 4._dp)
      xx        (3) = 3
      yy_correct(3) = 26
    elseif (par%i == 3) then
      call add_entry_CSR_dist( AA, 4, 3, 2._dp)
      call add_entry_CSR_dist( AA, 4, 6, 1._dp)
      xx        (4) = 4
      yy_correct(4) = 12
    elseif (par%i == 4) then
      call add_entry_CSR_dist( AA, 5, 4, 5._dp)
      call add_entry_CSR_dist( AA, 5, 7, 2._dp)
      xx        (5) = 5
      yy_correct(5) = 34
    elseif (par%i == 5) then
      call add_entry_CSR_dist( AA, 6, 2, 1._dp)
      call add_entry_CSR_dist( AA, 6, 5, 3._dp)
      xx        (6) = 6
      yy_correct(6) = 17
    elseif (par%i == 6) then
      call add_entry_CSR_dist( AA, 7, 3, 7._dp)
      call add_entry_CSR_dist( AA, 7, 7, 1._dp)
      xx        (7) = 7
      yy_correct(7) = 28
    end if

    call multiply_CSR_matrix_with_vector_1D( AA, xx, yy, xx_is_hybrid = .true., yy_is_hybrid = .true.)

    test_result = test_eq( yy( AA%i1:AA%i2), yy_correct( AA%i1:AA%i2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( xx, wxx)
    call deallocate_dist_shared( yy, wyy)
    call deallocate_dist_shared( yy_correct, wyy_correct)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_CSR_matrix_with_vector_1D_hybrid_hybrid

end module ut_mpi_multiply_CSR_matrix_vector
