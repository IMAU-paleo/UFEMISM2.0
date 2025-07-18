module ut_netcdf

  ! Unit tests for the netcdf i/o.

  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use ut_netcdf_xy_grid, only: unit_tests_netcdf_xy_grid
  use ut_netcdf_mesh, only: unit_tests_netcdf_mesh

  implicit none

  private

  public :: unit_tests_netcdf_main

contains

  subroutine unit_tests_netcdf_main( test_name_parent)
    ! Test the netcdf i/o subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_netcdf_main'
    character(len=1024), parameter :: test_name_local = 'netcdf'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call unit_tests_netcdf_xy_grid( test_name)
    call unit_tests_netcdf_mesh   ( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_main

end module ut_netcdf