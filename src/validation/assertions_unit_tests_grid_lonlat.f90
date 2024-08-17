module assertions_unit_tests_grid_lonlat

  ! The assertions/unit tests for simple lon/lat-grids.

  use assertions_unit_tests_basic, only: ASSERTION, UNIT_TEST, process_test_result
  use precisions, only: dp
  use grid_types, only: type_grid_lonlat

  implicit none

  private

  public :: test_tol_grid_lonlat

contains

subroutine test_tol_grid_lonlat( grid1, grid2, tol_angle_deg, test_mode, message)
  ! In/output variables:
  type(type_grid_lonlat),  intent(in   ) :: grid1, grid2
  real(dp),                intent(in   ) :: tol_angle_deg
  integer,                 intent(in   ) :: test_mode
  character(len=*),        intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = .true.

  test_result = test_result .and. grid1%nlon        == grid2%nlon
  test_result = test_result .and. grid1%nlat        == grid2%nlat
  test_result = test_result .and. grid1%n           == grid2%n
  test_result = test_result .and. (grid1%dlon        - grid2%dlon  ) <= tol_angle_deg
  test_result = test_result .and. (grid1%dlat        - grid2%dlat  ) <= tol_angle_deg
  test_result = test_result .and. (grid1%lonmin      - grid2%lonmin) <= tol_angle_deg
  test_result = test_result .and. (grid1%lonmax      - grid2%lonmax) <= tol_angle_deg
  test_result = test_result .and. (grid1%latmin      - grid2%latmin) <= tol_angle_deg
  test_result = test_result .and. (grid1%latmax      - grid2%latmax) <= tol_angle_deg

  if (size(grid1%lon) == size(grid2%lon)) then
    test_result = test_result .and. all(abs(grid1%lon - grid2%lon) <= tol_angle_deg)
  else
    test_result = .false.
  end if

  if (size(grid1%lat) == size(grid2%lat)) then
    test_result = test_result .and. all(abs(grid1%lat - grid2%lat) <= tol_angle_deg)
  else
    test_result = .false.
  end if

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_grid_lonlat

end module assertions_unit_tests_grid_lonlat
