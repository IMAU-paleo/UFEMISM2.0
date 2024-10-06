module ut_mesh_remapping

  ! Unit tests for mesh functions - remapping.

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use grid_types, only: type_grid
  use grid_basic, only: setup_square_grid
  use mpi_basic, only: par
  use ut_mesh_remapping_trace_line_grid_start, only: test_trace_line_grid_start
  use ut_mesh_remapping_trace_line_grid_a, only: test_trace_line_grid_a
  use ut_mesh_remapping_trace_line_grid_b, only: test_trace_line_grid_b
  use ut_mesh_remapping_trace_line_grid_cx, only: test_trace_line_grid_cx
  use ut_mesh_remapping_trace_line_grid_cy, only: test_trace_line_grid_cy

  implicit none

  private

  public :: test_remapping

contains

  subroutine test_remapping( test_name_parent)
    ! Test the remapping subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_remapping'
    character(len=1024), parameter :: test_name_local = 'remapping'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remapping

  subroutine test_trace_line_grid( test_name_parent)
    ! Test the trace_line_grid subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid'
    character(len=1024), parameter :: test_name_local = 'trace_line_grid'
    character(len=1024)            :: test_name
    real(dp)                       :: xmin, xmax, ymin, ymax, dx
    character(len=1024)            :: name
    type(type_grid)                :: grid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Initialise square grid
    name = 'test_grid'
    xmin = 0._dp
    xmax = 10._dp
    ymin = 0._dp
    ymax = 10._dp
    dx   = 1._dp
    call setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid)

    call test_trace_line_grid_start( test_name, grid)
    call test_trace_line_grid_a    ( test_name, grid)
    call test_trace_line_grid_b    ( test_name, grid)
    call test_trace_line_grid_cx   ( test_name, grid)
    call test_trace_line_grid_cy   ( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid

end module ut_mesh_remapping