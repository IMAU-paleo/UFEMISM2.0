module ut_mpi_CSR_matrix_vector_multiplication

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, finalise_matrix_CSR_dist
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D, &
    multiply_CSR_matrix_with_vector_1D_wrapper
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, MPI_WIN
  use mpi_distributed_shared_memory
  use mpi_distributed_memory, only: distribute_from_primary

  implicit none

  private

  public :: test_CSR_matrix_vector_multiplication_main

contains

  subroutine test_CSR_matrix_vector_multiplication_main( test_name_parent)
    ! Test the parallelised CSR matrix algebra subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_CSR_matrix_vector_multiplication_main'
    character(len=1024), parameter :: test_name_local = 'CSR_matrix_vector_multiplication'
    character(len=1024)            :: test_name
    type(type_sparse_matrix_CSR_dp) :: A1, A2
    real(dp), dimension(1)          :: x1, x2, y_correct1, y_correct2

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call initialise_simple_matrix_equation_1( A1, x1, y_correct1)
    call initialise_simple_matrix_equation_2( A2, x2, y_correct2)

    call test_multiply_CSR_matrix_with_vector_1D( test_name, A1, x1, y_correct1, 1)
    call test_multiply_CSR_matrix_with_vector_1D( test_name, A2, x2, y_correct2, 2)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_CSR_matrix_vector_multiplication_main

  subroutine test_multiply_CSR_matrix_with_vector_1D( test_name_parent, &
    A, x, y_correct, test_number)
    ! Test the multiply_CSR_matrix_with_vector subroutines

    ! In/output variables:
    character(len=*),                intent(in) :: test_name_parent
    type(type_sparse_matrix_CSR_dp), intent(in) :: A
    real(dp), dimension(1),          intent(in) :: x, y_correct
    integer,                         intent(in) :: test_number

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_multiply_CSR_matrix_with_vector_1D'
    character(len=1024), parameter              :: test_name_local = 'dp_1D'
    character(len=1024)                         :: test_name
    type(type_par_arr_info)                     :: pai
    real(dp), dimension(1)                      :: y
    real(dp), dimension(:), contiguous, pointer :: x_nih => null()
    real(dp), dimension(:), contiguous, pointer :: y_nih => null()
    real(dp), dimension(:), contiguous, pointer :: y_correct_nih => null()
    type(MPI_WIN)                               :: wx_nih, wy_nih, wy_correct_nih
    logical                                     :: test_result
    integer                                     :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    write( test_name,'(i1)') test_number
    test_name = trim( test_name_parent) // '/' // trim( test_name_local) // '_' // trim( test_name)

    ! Convert x,y to nih arrays
    call initialise_parallel_array_info( pai, &
      x_nih, wx_nih, y_nih, wy_nih, y_correct_nih, wy_correct_nih)

    call dist_to_hybrid( pai, x, x_nih)
    call dist_to_hybrid( pai, y, y_nih)

    ! Exchange halos
    call basic_halo_exchange( pai, x_nih)

    ! == Dist-dist ==
    ! ===============

    ! Perform matrix-vector multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( A, pai, x, pai, y, &
      xx_is_hybrid = .false., yy_is_hybrid = .false.)

    ! Verify result
    test_result = test_tol( y, y_correct, 1e-10_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_dist_dist')

    ! == Dist-hybrid ==
    ! =================

    ! Perform matrix-vector multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( A, pai, x, pai, y_nih, &
      xx_is_hybrid = .false., yy_is_hybrid = .true.)

    ! Convert y back to fully distributed array
    call hybrid_to_dist( pai, y_nih, y)

    ! Verify result
    test_result = test_tol( y, y_correct, 1e-10_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_dist_hybrid')

    ! == Hybrid-dist ==
    ! =================

    ! Perform matrix-vector multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( A, pai, x_nih, pai, y, &
      xx_is_hybrid = .true., yy_is_hybrid = .false.)

    ! Verify result
    test_result = test_tol( y, y_correct, 1e-10_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_hybrid_dist')

    ! == Hybrid-hybrid ==
    ! ===================

    ! Perform matrix-vector multiplication
    call multiply_CSR_matrix_with_vector_1D_wrapper( A, pai, x_nih, pai, y_nih, &
      xx_is_hybrid = .true., yy_is_hybrid = .true.)

    ! Convert y back to fully distributed array
    call hybrid_to_dist( pai, y_nih, y)

    ! Verify result
    test_result = test_tol( y, y_correct, 1e-10_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_hybrid_hybrid')

    ! == Hybrid hybrid, no wrapper ==
    ! ===============================

    ! Perform matrix-vector multiplication
    call multiply_CSR_matrix_with_vector_1D( A, pai, x_nih, pai, y_nih)

    ! Convert y back to fully distributed array
    call hybrid_to_dist( pai, y_nih, y)

    ! Verify result
    test_result = test_tol( y, y_correct, 1e-10_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_hybrid_hybrid_nowrapper')

    ! Clean up after yourself
    call deallocate_dist_shared( x_nih, wx_nih)
    call deallocate_dist_shared( y_nih, wy_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_CSR_matrix_with_vector_1D

  subroutine initialise_simple_matrix_equation_1( A, x, y)

    ! In/output variables:
    type(type_sparse_matrix_CSR_dp), intent(out) :: A
    real(dp), dimension(1),          intent(out) :: x,y

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_simple_matrix_equation_1'

    ! Add routine to call stack
    call init_routine( routine_name)

    !              A            x  =   y
    !     1  2  3  4  5  6  7
    ! 1  [1,  ,  ,  ,  ,  ,  ] [1]   [ 1]
    ! 2  [ , 2, 3,  ,  ,  ,  ] [2]   [13]
    ! 3  [ , 4,  ,  , 1,  ,  ] [3]   [ 9]
    ! 4  [ ,  , 2, 3,  , 4,  ] [4] = [26]
    ! 5  [ ,  ,  , 1, 2,  ,  ] [1]   [ 6]
    ! 6  [ ,  ,  ,  , 3,  , 4] [2]   [15]
    ! 7  [ ,  ,  ,  ,  , 1, 2] [3]   [ 8]

    call allocate_matrix_CSR_dist( A, 7, 7, 1, 1, 3)

    if     (par%i == 0) then
      call add_entry_CSR_dist( A, 1, 1, 1._dp)
      x = 1._dp
      y = 1._dp
    elseif (par%i == 1) then
      call add_entry_CSR_dist( A, 2, 2, 2._dp)
      call add_entry_CSR_dist( A, 2, 3, 3._dp)
      x = 2._dp
      y = 13._dp
    elseif (par%i == 2) then
      call add_entry_CSR_dist( A, 3, 2, 4._dp)
      call add_entry_CSR_dist( A, 3, 5, 1._dp)
      x = 3._dp
      y = 9._dp
    elseif (par%i == 3) then
      call add_entry_CSR_dist( A, 4, 3, 2._dp)
      call add_entry_CSR_dist( A, 4, 4, 3._dp)
      call add_entry_CSR_dist( A, 4, 6, 4._dp)
      x = 4._dp
      y = 26._dp
    elseif (par%i == 4) then
      call add_entry_CSR_dist( A, 5, 4, 1._dp)
      call add_entry_CSR_dist( A, 5, 5, 2._dp)
      x = 1._dp
      y = 6._dp
    elseif (par%i == 5) then
      call add_entry_CSR_dist( A, 6, 5, 3._dp)
      call add_entry_CSR_dist( A, 6, 7, 4._dp)
      x = 2._dp
      y = 15._dp
    elseif (par%i == 6) then
      call add_entry_CSR_dist( A, 7, 6, 1._dp)
      call add_entry_CSR_dist( A, 7, 7, 2._dp)
      x = 3._dp
      y = 8._dp
    end if

    call finalise_matrix_CSR_dist( A)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_simple_matrix_equation_1

  subroutine initialise_simple_matrix_equation_2( A, x, y)

    ! In/output variables:
    type(type_sparse_matrix_CSR_dp), intent(out) :: A
    real(dp), dimension(1),          intent(out) :: x,y

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_simple_matrix_equation_2'

    ! Add routine to call stack
    call init_routine( routine_name)

    !              A            x    =  y
    !     1  2  3  4  5  6  7
    ! 1  [1,  ,  ,  , 5,  ,  ] [1]   [ 6]
    ! 2  [ , 2, 3,  ,  , 6,  ] [2]   [25]
    ! 3  [ , 4,  ,  , 1,  ,  ] [3]   [ 9]
    ! 4  [5,  , 2, 3,  , 4, 5] [4] = [46]
    ! 5  [ , 6,  , 1, 2,  ,  ] [1]   [18]
    ! 6  [5,  ,  ,  , 3,  , 4] [2]   [20]
    ! 7  [ , 6,  ,  ,  , 1, 2] [3]   [20]

    call allocate_matrix_CSR_dist( A, 7, 7, 1, 1, 7)

    if     (par%i == 0) then
      call add_entry_CSR_dist( A, 1, 1, 1._dp)
      call add_entry_CSR_dist( A, 1, 5, 5._dp)
      x = 1._dp
      y = 6._dp
    elseif (par%i == 1) then
      call add_entry_CSR_dist( A, 2, 2, 2._dp)
      call add_entry_CSR_dist( A, 2, 3, 3._dp)
      call add_entry_CSR_dist( A, 2, 6, 6._dp)
      x = 2._dp
      y = 25._dp
    elseif (par%i == 2) then
      call add_entry_CSR_dist( A, 3, 2, 4._dp)
      call add_entry_CSR_dist( A, 3, 5, 1._dp)
      x = 3._dp
      y = 9._dp
    elseif (par%i == 3) then
      call add_entry_CSR_dist( A, 4, 1, 5._dp)
      call add_entry_CSR_dist( A, 4, 3, 2._dp)
      call add_entry_CSR_dist( A, 4, 4, 3._dp)
      call add_entry_CSR_dist( A, 4, 6, 4._dp)
      call add_entry_CSR_dist( A, 4, 7, 5._dp)
      x = 4._dp
      y = 46._dp
    elseif (par%i == 4) then
      call add_entry_CSR_dist( A, 5, 2, 6._dp)
      call add_entry_CSR_dist( A, 5, 4, 1._dp)
      call add_entry_CSR_dist( A, 5, 5, 2._dp)
      x = 1._dp
      y = 18._dp
    elseif (par%i == 5) then
      call add_entry_CSR_dist( A, 6, 1, 5._dp)
      call add_entry_CSR_dist( A, 6, 5, 3._dp)
      call add_entry_CSR_dist( A, 6, 7, 4._dp)
      x = 2._dp
      y = 20._dp
    elseif (par%i == 6) then
      call add_entry_CSR_dist( A, 7, 2, 6._dp)
      call add_entry_CSR_dist( A, 7, 6, 1._dp)
      call add_entry_CSR_dist( A, 7, 7, 2._dp)
      x = 3._dp
      y = 20._dp
    end if

    call finalise_matrix_CSR_dist( A)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_simple_matrix_equation_2

  subroutine initialise_parallel_array_info( pai, &
    x_nih, wx_nih, y_nih, wy_nih, y_correct_nih, wy_correct_nih)

    ! In/output variables:
    type(type_par_arr_info),         intent(  out) :: pai
    real(dp), dimension(:), pointer, intent(  out) :: x_nih
    type(MPI_WIN),                   intent(  out) :: wx_nih
    real(dp), dimension(:), pointer, intent(  out) :: y_nih
    type(MPI_WIN),                   intent(  out) :: wy_nih
    real(dp), dimension(:), pointer, intent(  out) :: y_correct_nih
    type(MPI_WIN),                   intent(  out) :: wy_correct_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_parallel_array_info'

    ! Add routine to call stack
    call init_routine( routine_name)

    pai%n = 7
    pai%i1 = par%i+1
    pai%i2 = par%i+1

    if (par%node_ID == 0) then

      pai%i1_node = 1
      pai%i2_node = 2

      pai%i1_nih  = 1
      pai%i2_nih  = 3

      pai%i1_hle = 0
      pai%i2_hle = -1
      pai%i1_hli = 0
      pai%i2_hli = -1
      pai%i1_hre = 3
      pai%i2_hre = 3
      pai%i1_hri = 2
      pai%i2_hri = 2

    elseif (par%node_ID == 1) then

      pai%i1_node = 3
      pai%i2_node = 5

      pai%i1_nih  = 2
      pai%i2_nih  = 6

      pai%i1_hle = 2
      pai%i2_hle = 2
      pai%i1_hli = 3
      pai%i2_hli = 3
      pai%i1_hre = 6
      pai%i2_hre = 6
      pai%i1_hri = 5
      pai%i2_hri = 5

    elseif (par%node_ID == 2) then

      pai%i1_node = 6
      pai%i2_node = 7

      pai%i1_nih  = 5
      pai%i2_nih  = 7

      pai%i1_hle = 5
      pai%i2_hle = 5
      pai%i1_hli = 6
      pai%i2_hli = 6
      pai%i1_hre = 0
      pai%i2_hre = -1
      pai%i1_hri = 0
      pai%i2_hri = -1

    end if

    pai%n_node = pai%i2_node + 1 - pai%i1_node
    pai%n_nih  = pai%i2_nih  + 1 - pai%i1_nih
    pai%n_hle  = pai%i2_hle  + 1 - pai%i1_hle
    pai%n_hli  = pai%i2_hli  + 1 - pai%i1_hli
    pai%n_hre  = pai%i2_hre  + 1 - pai%i1_hre
    pai%n_hri  = pai%i2_hri  + 1 - pai%i1_hri

    x_nih => null()
    y_nih => null()
    y_correct_nih => null()

    call allocate_dist_shared( x_nih        , wx_nih        , pai%n_nih)
    call allocate_dist_shared( y_nih        , wy_nih        , pai%n_nih)
    call allocate_dist_shared( y_correct_nih, wy_correct_nih, pai%n_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_parallel_array_info

end module ut_mpi_CSR_matrix_vector_multiplication
