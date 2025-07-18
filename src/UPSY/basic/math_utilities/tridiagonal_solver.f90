module tridiagonal_solver

  use precisions, only: dp

  implicit none

  private

  public :: solve_tridiagonal_matrix_equation

contains

  subroutine solve_tridiagonal_matrix_equation( n, a, b, c, &
    d, x)
    ! Solve a tridiagonal matrix equation of the form:
    !
    ! [ b(1) c(1)                          ] [x(1)  ]   [d(1)  ]
    ! [ a(2) b(2) c(2)                     ] [x(2)  ]   [d(2)  ]
    ! [      a(3) b(3) c(3)                ] [x(3)  ]   [d(3)  ]
    ! [           ...  ...  ...            ] [...   ] = [...   ]
    ! [                a(n-1 b(n-1) c(n-1) ] [x(n-1)]   [d(n-1)]
    ! [                      a(n)   b(n)   ] [x(n)  ]   [d(n)  ]
    !
    ! (see: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)

    ! In/output variables:
    integer,                    intent(in   ) :: n
    real(dp), dimension(2:n  ), intent(inout) :: a
    real(dp), dimension(1:n  ), intent(inout) :: b
    real(dp), dimension(1:n-1), intent(inout) :: c
    real(dp), dimension(1:n  ), intent(inout) :: d
    real(dp), dimension(1:n  ), intent(  out) :: x

    ! Local variables:
    integer  :: i
    real(dp) :: w

    do i = 2, n
      w = a( i) / b( i-1)
      b( i) = b( i) - w * c( i-1)
      d( i) = d( i) - w * d( i-1)
    end do

    x = 0._dp

    x( n) = d( n) / b( n)
    do i = n-1, 1, -1
      x( i) = (d( i) - c( i) * x( i+1)) / b( i)
    end do

  end subroutine solve_tridiagonal_matrix_equation

end module tridiagonal_solver