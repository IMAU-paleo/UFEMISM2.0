module ut_mesh_remapping

  ! Unit tests for mesh functions - remapping.

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use parameters, only: pi
  use grid_types, only: type_grid
  use grid_basic, only: setup_square_grid
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use ut_mesh_remapping_trace_line_grid_start, only: test_trace_line_grid_start
  use ut_mesh_remapping_trace_line_grid_a, only: test_trace_line_grid_a
  use ut_mesh_remapping_trace_line_grid_b, only: test_trace_line_grid_b
  use ut_mesh_remapping_trace_line_grid_cx, only: test_trace_line_grid_cx
  use ut_mesh_remapping_trace_line_grid_cy, only: test_trace_line_grid_cy
  use ut_mesh_remapping_trace_line_grid, only: test_trace_line_grid
  use ut_mesh_remapping_trace_line_tri_start, only: test_trace_line_tri_start
  use ut_mesh_remapping_trace_line_tri_ti, only: test_trace_line_tri_ti
  use ut_mesh_remapping_trace_line_tri_vi, only: test_trace_line_tri_vi
  use ut_mesh_remapping_trace_line_tri_ei, only: test_trace_line_tri_ei
  use ut_mesh_remapping_trace_line_tri, only: test_trace_line_tri
  use ut_mesh_remapping_trace_line_Vor_start, only: test_trace_line_Vor_start
  use ut_mesh_remapping_trace_line_Vor_vi, only: test_trace_line_Vor_vi
  use ut_mesh_remapping_trace_line_Vor_ti, only: test_trace_line_Vor_ti
  use ut_mesh_remapping_trace_line_Vor_ei, only: test_trace_line_Vor_ei
  use ut_mesh_remapping_trace_line_Vor, only: test_trace_line_Vor
  use ut_mesh_remapping_mesh_to_mesh, only: test_remapping_mesh_to_mesh

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

    call test_trace_line_grid_main  ( test_name)
    call test_trace_line_tri_main   ( test_name)
    call test_trace_line_Vor_main   ( test_name)
    call test_remapping_mesh_to_mesh( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remapping

  subroutine test_trace_line_grid_main( test_name_parent)
    ! Test the trace_line_grid subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_main'
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
    xmax = 10._dp * pi
    ymin = 0._dp
    ymax = 10._dp * pi
    dx   = pi
    call setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid)

    call test_trace_line_grid_start( test_name, grid)
    call test_trace_line_grid_a    ( test_name, grid)
    call test_trace_line_grid_b    ( test_name, grid)
    call test_trace_line_grid_cx   ( test_name, grid)
    call test_trace_line_grid_cy   ( test_name, grid)
    call test_trace_line_grid      ( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_main

  subroutine test_trace_line_tri_main( test_name_parent)
    ! Test the trace_line_tri subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_main'
    character(len=1024), parameter :: test_name_local = 'trace_line_tri'
    character(len=1024)            :: test_name
    real(dp)                       :: xmin, xmax, ymin, ymax, alpha_min, res_max
    character(len=1024)            :: name
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    name = 'test_mesh'
    xmin = 0._dp
    xmax = pi
    ymin = 0._dp
    ymax = pi
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 20._dp

    call allocate_mesh_primary( mesh, name, 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    ! Run unit tests on this mesh
    call test_trace_line_tri_start( test_name, mesh)
    call test_trace_line_tri_ti   ( test_name, mesh)
    call test_trace_line_tri_vi   ( test_name, mesh)
    call test_trace_line_tri_ei   ( test_name, mesh)
    call test_trace_line_tri      ( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_main

  subroutine test_trace_line_Vor_main( test_name_parent)
    ! Test the trace_line_Vor subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_main'
    character(len=1024), parameter :: test_name_local = 'trace_line_Vor'
    character(len=1024)            :: test_name
    real(dp)                       :: xmin, xmax, ymin, ymax, alpha_min, res_max
    character(len=1024)            :: name
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    name = 'test_mesh'
    xmin = 0._dp
    xmax = pi
    ymin = 0._dp
    ymax = pi
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 20._dp

    call allocate_mesh_primary( mesh, name, 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    ! Run unit tests on this mesh
    call test_trace_line_Vor_start( test_name, mesh)
    call test_trace_line_Vor_vi   ( test_name, mesh)
    call test_trace_line_Vor_ti   ( test_name, mesh)
    call test_trace_line_Vor_ei   ( test_name, mesh)
    call test_trace_line_Vor      ( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_main

end module ut_mesh_remapping