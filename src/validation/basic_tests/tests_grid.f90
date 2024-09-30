module tests_grid

  ! The assertions/unit tests for simple square x/y-grids.

  use assertions_unit_tests_basic, only: ASSERTION, UNIT_TEST, process_test_result
  use precisions, only: dp
  use grid_types, only: type_grid

  implicit none

  private

  public :: test_tol_grid

contains

  subroutine test_tol_grid( grid1, grid2, tol_dist, test_mode, message)
    ! In/output variables:
    type(type_grid),  intent(in   ) :: grid1, grid2
    real(dp),         intent(in   ) :: tol_dist
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = .true.

    test_result = test_result .and. grid1%nx          == grid2%nx
    test_result = test_result .and. grid1%ny          == grid2%ny
    test_result = test_result .and. grid1%n           == grid2%n
    test_result = test_result .and. (grid1%dx          - grid2%dx         ) <= tol_dist
    test_result = test_result .and. (grid1%xmin        - grid2%xmin       ) <= tol_dist
    test_result = test_result .and. (grid1%xmax        - grid2%xmax       ) <= tol_dist
    test_result = test_result .and. (grid1%ymin        - grid2%ymin       ) <= tol_dist
    test_result = test_result .and. (grid1%ymax        - grid2%ymax       ) <= tol_dist
    test_result = test_result .and. (grid1%lambda_M    - grid2%lambda_M   ) <= tol_dist
    test_result = test_result .and. (grid1%phi_M       - grid2%phi_M      ) <= tol_dist
    test_result = test_result .and. (grid1%beta_stereo - grid2%beta_stereo) <= tol_dist

    if (size(grid1%x) == size(grid2%x)) then
      test_result = test_result .and. all(abs(grid1%x - grid2%x) <= tol_dist)
    else
      test_result = .false.
    end if

    if (size(grid1%y) == size(grid2%y)) then
      test_result = test_result .and. all(abs(grid1%y - grid2%y) <= tol_dist)
    else
      test_result = .false.
    end if

    call process_test_result( test_mode, test_result, message)

  end subroutine test_tol_grid

end module tests_grid
