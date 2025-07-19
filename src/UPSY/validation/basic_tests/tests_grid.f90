module tests_grid

  ! Basic tests for simple square x/y-grids.

  use precisions, only: dp
  use grid_types, only: type_grid

  implicit none

  private

  public :: test_tol_grid

contains

  pure function test_tol_grid( grid1, grid2, tol_dist) result( res)
    type(type_grid),  intent(in) :: grid1, grid2
    real(dp),         intent(in) :: tol_dist
    logical :: res

    res = .true.

    res = res .and. grid1%nx          == grid2%nx
    res = res .and. grid1%ny          == grid2%ny
    res = res .and. grid1%n           == grid2%n
    res = res .and. (grid1%dx          - grid2%dx         ) <= tol_dist
    res = res .and. (grid1%xmin        - grid2%xmin       ) <= tol_dist
    res = res .and. (grid1%xmax        - grid2%xmax       ) <= tol_dist
    res = res .and. (grid1%ymin        - grid2%ymin       ) <= tol_dist
    res = res .and. (grid1%ymax        - grid2%ymax       ) <= tol_dist
    res = res .and. (grid1%lambda_M    - grid2%lambda_M   ) <= tol_dist
    res = res .and. (grid1%phi_M       - grid2%phi_M      ) <= tol_dist
    res = res .and. (grid1%beta_stereo - grid2%beta_stereo) <= tol_dist

    if (size(grid1%x) == size(grid2%x)) then
      res = res .and. all(abs(grid1%x - grid2%x) <= tol_dist)
    else
      res = .false.
    end if

    if (size(grid1%y) == size(grid2%y)) then
      res = res .and. all(abs(grid1%y - grid2%y) <= tol_dist)
    else
      res = .false.
    end if

  end function test_tol_grid

end module tests_grid
