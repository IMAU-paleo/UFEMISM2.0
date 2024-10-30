module ut_mesh

  ! Unit tests for mesh functions.

  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use ut_mesh_Delaunay, only: test_Delaunay
  use ut_mesh_refinement, only: test_refinement
  use ut_mesh_remapping, only: test_remapping

  implicit none

  private

  public :: unit_tests_mesh_main

contains

  subroutine unit_tests_mesh_main( test_name_parent)
    ! Run all unit tests for the MPI distributed memory subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_mesh_main'
    character(len=1024), parameter :: test_name_local = 'mesh'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Run all unit tests for the mesh creation subroutines
    call test_Delaunay  ( test_name)
    call test_refinement( test_name)
    call test_remapping ( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_mesh_main

end module ut_mesh