module ut_mpi_CSR_matrix_solving

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, finalise_matrix_CSR_dist
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, MPI_WIN
  use mpi_distributed_shared_memory
  use mpi_distributed_memory, only: distribute_from_primary
  use CSR_matrix_solving, only: solve_matrix_equation_CSR_Jacobi_wrapper, solve_matrix_equation_CSR_Jacobi

  implicit none

  private

  public :: test_CSR_matrix_solving_main

contains

  subroutine test_CSR_matrix_solving_main( test_name_parent)
    ! Test the parallelised CSR matrix algebra subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_CSR_matrix_solving_main'
    character(len=1024), parameter  :: test_name_local = 'CSR_matrix_solving'
    character(len=1024)             :: test_name
    type(type_sparse_matrix_CSR_dp) :: A1, A2
    real(dp), dimension(1)          :: b1, b2, x_ex1, x_ex2

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call initialise_simple_matrix_equation( A1, b1, x_ex1, .false.)
    call initialise_simple_matrix_equation( A2, b2, x_ex2, .true.)

    call test_solve_CSR_Jacobi( test_name, A1, b1, x_ex1, 1)
    call test_solve_CSR_Jacobi( test_name, A2, b2, x_ex2, 2)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_CSR_matrix_solving_main

  subroutine test_solve_CSR_Jacobi( test_name_parent, A, b, x_ex, test_number)
    ! Test the solve_CSR_Jacobi code

    ! In/output variables:
    character(len=*),                intent(in) :: test_name_parent
    type(type_sparse_matrix_CSR_dp), intent(in) :: A
    real(dp), dimension(1),          intent(in) :: b, x_ex
    integer,                         intent(in) :: test_number

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_solve_CSR_SOR'
    character(len=1024), parameter              :: test_name_local = 'Jacobi'
    character(len=1024)                         :: test_name
    type(type_par_arr_info)                     :: pai
    real(dp), dimension(1)                      :: x = 0._dp
    real(dp), dimension(:), contiguous, pointer :: b_nih => null()
    real(dp), dimension(:), contiguous, pointer :: x_nih => null()
    type(MPI_WIN)                               :: wb_nih, wx_nih
    integer                                     :: nit
    real(dp)                                    :: tol
    logical                                     :: test_result
    integer                                     :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    write( test_name,'(i1)') test_number
    test_name = trim( test_name_parent) // '/' // trim( test_name_local) // '_' // trim( test_name)

    ! Convert b,x,x_ex to nih arrays
    call initialise_parallel_array_info( pai, b_nih, wb_nih, x_nih, wx_nih)

    call dist_to_hybrid( pai, b, b_nih)
    call dist_to_hybrid( pai, x, x_nih)

    ! Exchange halos
    call basic_halo_exchange( pai, x_nih)

    ! == Dist-dist ==
    ! ===============

    ! Solve matrix equation
    nit = 500
    tol = 1e-7_dp
    x   = 0._dp
    call solve_matrix_equation_CSR_Jacobi_wrapper( A, pai, x, pai, b, nit, tol, &
      xx_is_hybrid = .false., bb_is_hybrid = .false.)

    ! Verify result
    test_result = test_tol( x, x_ex, 1e-5_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_dist_dist')

    ! == Dist-hybrid ==
    ! =================

    ! Solve matrix equation
    nit = 500
    tol = 1e-7_dp
    x   = 0._dp
    call solve_matrix_equation_CSR_Jacobi_wrapper( A, pai, x, pai, b_nih, nit, tol, &
      xx_is_hybrid = .false., bb_is_hybrid = .true.)

    ! Verify result
    test_result = test_tol( x, x_ex, 1e-5_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_dist_hybrid')

    ! == Hybrid-dist ==
    ! =================

    ! Solve matrix equation
    nit = 500
    tol = 1e-7_dp
    x   = 0._dp
    call solve_matrix_equation_CSR_Jacobi_wrapper( A, pai, x_nih, pai, b, nit, tol, &
      xx_is_hybrid = .true., bb_is_hybrid = .false.)

    ! Convert x back to fully distributed array
    call hybrid_to_dist( pai, x_nih, x)

    ! Verify result
    test_result = test_tol( x, x_ex, 1e-5_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_hybrid_dist')

    ! == Hybrid-hybrid ==
    ! ===================

    ! Solve matrix equation
    nit = 500
    tol = 1e-7_dp
    x   = 0._dp
    call solve_matrix_equation_CSR_Jacobi_wrapper( A, pai, x_nih, pai, b_nih, nit, tol, &
      xx_is_hybrid = .true., bb_is_hybrid = .true.)

    ! Convert x back to fully distributed array
    call hybrid_to_dist( pai, x_nih, x)

    ! Verify result
    test_result = test_tol( x, x_ex, 1e-5_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_hybrid_hybrid')

    ! == Hybrid hybrid, no wrapper ==
    ! ===============================

    ! Solve matrix equation
    nit = 500
    tol = 1e-7_dp
    x   = 0._dp
    call solve_matrix_equation_CSR_Jacobi( A, pai, x_nih, pai, b_nih, nit, tol)

    ! Convert x back to fully distributed array
    call hybrid_to_dist( pai, x_nih, x)

    ! Verify result
    test_result = test_tol( x, x_ex, 1e-5_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '_hybrid_hybrid_nowrapper')

    ! Clean up after yourself
    call deallocate_dist_shared( x_nih, wx_nih)
    call deallocate_dist_shared( b_nih, wb_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_solve_CSR_Jacobi

  subroutine initialise_simple_matrix_equation( A, b, x, with_halos)

    ! In/output variables:
    type(type_sparse_matrix_CSR_dp), intent(  out) :: A
    real(dp), dimension(1),          intent(  out) :: b, x
    logical,                         intent(in   ) :: with_halos

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'initialise_simple_matrix_equation'
    type(type_par_arr_info)                     :: pai
    real(dp), dimension(:), contiguous, pointer :: b_nih => null()
    real(dp), dimension(:), contiguous, pointer :: x_nih => null()
    type(MPI_WIN)                               :: wb_nih, wx_nih

    ! Add routine to call stack
    call init_routine( routine_name)

    if (with_halos) then
      call initialise_parallel_array_info( pai, b_nih, wb_nih, x_nih, wx_nih)
      call allocate_matrix_CSR_dist( A, 7, 7, 1, 1, 3, pai_x = pai, pai_y = pai)
    else
      call allocate_matrix_CSR_dist( A, 7, 7, 1, 1, 3)
    end if

    !              A                    x   =   b
    !      1   2   3   4   5   6   7
    ! 1  [ 1,   ,   ,   ,   ,   ,   ] [1.0]   [1.0]
    ! 2  [-1,  2, -1,   ,   ,   ,   ] [3.5]   [1.0]
    ! 3  [  , -1,  2, -1,   ,   ,   ] [5.0]   [1.0]
    ! 4  [  ,   , -1,  2, -1,   ,   ] [5.5] = [1.0]
    ! 5  [  ,   ,   , -1,  2, -1,   ] [5.0]   [1.0]
    ! 6  [  ,   ,   ,   , -1,  2, -1] [3.5]   [1.0]
    ! 7  [  ,   ,   ,   ,   ,   ,  1] [1.0]   [1.0]


    if     (par%i == 0) then
      call add_entry_CSR_dist( A, 1, 1, 1._dp)
      x = 1._dp
      b = 1._dp
    elseif (par%i == 1) then
      call add_entry_CSR_dist( A, 2, 1, -1._dp)
      call add_entry_CSR_dist( A, 2, 2,  2._dp)
      call add_entry_CSR_dist( A, 2, 3, -1._dp)
      x = 3.5_dp
      b = 1._dp
    elseif (par%i == 2) then
      call add_entry_CSR_dist( A, 3, 2, -1._dp)
      call add_entry_CSR_dist( A, 3, 3,  2._dp)
      call add_entry_CSR_dist( A, 3, 4, -1._dp)
      x = 5._dp
      b = 1._dp
    elseif (par%i == 3) then
      call add_entry_CSR_dist( A, 4, 3, -1._dp)
      call add_entry_CSR_dist( A, 4, 4,  2._dp)
      call add_entry_CSR_dist( A, 4, 5, -1._dp)
      x = 5.5_dp
      b = 1._dp
    elseif (par%i == 4) then
      call add_entry_CSR_dist( A, 5, 4, -1._dp)
      call add_entry_CSR_dist( A, 5, 5,  2._dp)
      call add_entry_CSR_dist( A, 5, 6, -1._dp)
      x = 5._dp
      b = 1._dp
    elseif (par%i == 5) then
      call add_entry_CSR_dist( A, 6, 5, -1._dp)
      call add_entry_CSR_dist( A, 6, 6,  2._dp)
      call add_entry_CSR_dist( A, 6, 7, -1._dp)
      x = 3.5_dp
      b = 1._dp
    elseif (par%i == 6) then
      call add_entry_CSR_dist( A, 7, 7, 1._dp)
      x = 1._dp
      b = 1._dp
    end if

    call finalise_matrix_CSR_dist( A)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_simple_matrix_equation

  subroutine initialise_parallel_array_info( pai, b_nih, wb_nih, x_nih, wx_nih)

    ! In/output variables:
    type(type_par_arr_info),         intent(  out) :: pai
    real(dp), dimension(:), pointer, intent(  out) :: b_nih
    type(MPI_WIN),                   intent(  out) :: wb_nih
    real(dp), dimension(:), pointer, intent(  out) :: x_nih
    type(MPI_WIN),                   intent(  out) :: wx_nih

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

    b_nih => null()
    x_nih => null()

    call allocate_dist_shared( b_nih, wb_nih, pai%n_nih)
    call allocate_dist_shared( x_nih, wx_nih, pai%n_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_parallel_array_info

end module ut_mpi_CSR_matrix_solving
