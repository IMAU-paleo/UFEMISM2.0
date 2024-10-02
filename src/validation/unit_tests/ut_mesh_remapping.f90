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
  use ut_mesh_remapping_trace_line_grid, only: test_trace_line_grid

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

end module ut_mesh_remapping