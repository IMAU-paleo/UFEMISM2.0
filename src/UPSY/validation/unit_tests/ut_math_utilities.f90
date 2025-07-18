module ut_math_utilities

  ! Unit tests for the math utilities.

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use tridiagonal_solver, only: solve_tridiagonal_matrix_equation

  implicit none

  private

  public :: unit_tests_math_utilities_main

contains

  subroutine unit_tests_math_utilities_main( test_name_parent)
    ! Test the math_utilities subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_math_utilities_main'
    character(len=1024), parameter :: test_name_local = 'math_utilities'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_tridiagonal_solver( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_math_utilities_main

  subroutine test_tridiagonal_solver( test_name_parent)
    ! Test the math_utilities subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_tridiagonal_solver'
    character(len=1024), parameter      :: test_name_local = 'tridiagonal_solver'
    character(len=1024)                 :: test_name
    integer                             :: ni,n,i_mid
    real(dp), dimension(:), allocatable :: a,b,c,d,x
    real(dp)                            :: dx
    logical                             :: verified

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Solve the tridiagonal matrix equation representing d2f/dx2 = 2, f(-1) = 1, f(1) = 1
    ! (which has the analytical solution f = x^2) for different resolutions.
    ! Check that f(x=0) = 0

    verified = .true.

    do ni = 3, 12

      n     = 1 + 2**ni
      i_mid = 1 + 2**(ni-1)

      allocate( a( 2:n  ), source = 0._dp)
      allocate( b( 1:n  ), source = 0._dp)
      allocate( c( 1:n-1), source = 0._dp)
      allocate( d( 1:n  ), source = 0._dp)
      allocate( x( 1:n  ), source = 0._dp)

      dx = 2._dp / real( n-1,dp)

      a =  1._dp / dx**2
      b = -2._dp / dx**2
      c =  1._dp / dx**2
      d =  2._dp

      ! Boundary conditions
      b(1) = 1._dp
      c(1) = 0._dp
      d(1) = 1._dp

      a(n) = 0._dp
      b(n) = 1._dp
      d(n) = 1._dp

      ! Solve the equation
      call solve_tridiagonal_matrix_equation( n, a, b, c, d, x)

      ! Check the solution
      verified = verified .and. &
        test_tol( x( 1    ), 1._dp, 1e-12_dp) .and. &
        test_tol( x( i_mid), 0._dp, 1e-12_dp) .and. &
        test_tol( x( n    ), 1._dp, 1e-12_dp)

      ! Deallocate memory to make room for the next round
      deallocate( a,b,c,d,x)

    end do

    call unit_test( verified, trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_tridiagonal_solver

end module ut_math_utilities