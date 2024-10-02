module tests_grid_lonlat

  ! Basic tests for simple lon/lat-grids.

  use precisions, only: dp
  use grid_types, only: type_grid_lonlat

  implicit none

  private

  public :: test_tol_grid_lonlat

contains

  pure function test_tol_grid_lonlat( grid1, grid2, tol_angle_deg) result( res)
    type(type_grid_lonlat), intent(in) :: grid1, grid2
    real(dp),               intent(in) :: tol_angle_deg
    logical :: res

    res = .true.

    res = res .and. grid1%nlon        == grid2%nlon
    res = res .and. grid1%nlat        == grid2%nlat
    res = res .and. grid1%n           == grid2%n
    res = res .and. (grid1%dlon        - grid2%dlon  ) <= tol_angle_deg
    res = res .and. (grid1%dlat        - grid2%dlat  ) <= tol_angle_deg
    res = res .and. (grid1%lonmin      - grid2%lonmin) <= tol_angle_deg
    res = res .and. (grid1%lonmax      - grid2%lonmax) <= tol_angle_deg
    res = res .and. (grid1%latmin      - grid2%latmin) <= tol_angle_deg
    res = res .and. (grid1%latmax      - grid2%latmax) <= tol_angle_deg

    if (size(grid1%lon) == size(grid2%lon)) then
      res = res .and. all(abs(grid1%lon - grid2%lon) <= tol_angle_deg)
    else
      res = .false.
    end if

    if (size(grid1%lat) == size(grid2%lat)) then
      res = res .and. all(abs(grid1%lat - grid2%lat) <= tol_angle_deg)
    else
      res = .false.
    end if

  end function test_tol_grid_lonlat

end module tests_grid_lonlat
