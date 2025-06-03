module ut_mpi_CSR_matrix_algebra

  ! Unit tests for different MPI routines

  use tests_main
  use assertions_basic
  use ut_basic
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_basic, only: par
  use ut_mpi_CSR_matrix_vector_multiplication, only: test_CSR_matrix_vector_multiplication_main
  use ut_mpi_CSR_matrix_solving, only: test_CSR_matrix_solving_main

  implicit none

  private

  public :: test_CSR_matrix_algebra_main

contains

  subroutine test_CSR_matrix_algebra_main( test_name_parent)
    ! Test the parallelised CSR matrix algebra subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_CSR_matrix_algebra_main'
    character(len=1024), parameter :: test_name_local = 'CSR_matrix_algebra_main'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_CSR_matrix_vector_multiplication_main( test_name)
    call test_CSR_matrix_solving_main( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_CSR_matrix_algebra_main

end module ut_mpi_CSR_matrix_algebra
