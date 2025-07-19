module ut_mesh_Delaunay

  ! Unit tests for mesh functions - Delaunay triangulation.

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use parameters
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use split_border_edges, only: split_border_edge
  use split_edges, only: split_edge
  use split_triangles, only: split_triangle
  use move_vertices, only: move_vertex

  implicit none

  private

  public :: test_Delaunay

contains

  subroutine test_Delaunay( test_name_parent)
    ! Test the Delaunay triangulation subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_Delaunay'
    character(len=1024), parameter :: test_name_local = 'Delaunay'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_split_triangle   ( test_name)
    call test_split_edge       ( test_name)
    call test_split_border_edge( test_name)
    call test_move_vertex      ( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_Delaunay

  subroutine test_split_triangle( test_name_parent)
    ! Test the split_triangle subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_split_triangle'
    character(len=1024), parameter :: test_name_local = 'split_triangle'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), dimension(2)         :: p_new

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southern triangle
    p_new = [0._dp, -1._dp]
    call split_triangle( mesh, 1, p_new)

    ! Check if there is indeed a new vertex in the location we wanted
    call unit_test( test_eq ( mesh%nV, 6), trim(test_name)//'/new_vertex')
    call unit_test( test_tol( mesh%V(6,:), p_new, 1E-10_dp), trim(test_name)//'/new_vertex_location')

    ! Check if there is indeed one new triangle
    call unit_test( test_eq ( mesh%nTri, 5), trim(test_name)//'/new_triangle')

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_split_triangle

  subroutine test_split_edge( test_name_parent)
    ! Test the split_edge subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_split_edge'
    character(len=1024), parameter :: test_name_local = 'split_edge'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), dimension(2)         :: p_new

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southwestern edge
    p_new = (mesh%V(1,:) + mesh%V(5,:)) / 2._dp
    call split_edge( mesh, 1, 5, p_new)

    ! Check if there is indeed a new vertex in the location we wanted
    call unit_test( test_eq ( mesh%nV, 6), trim(test_name)//'/new_vertex')
    call unit_test( test_tol( mesh%V(6,:), p_new, 1E-10_dp), trim(test_name)//'/new_vertex_location')

    ! Check if there are indeed two new triangles
    call unit_test( test_eq ( mesh%nTri, 6), trim(test_name)//'/new_triangles')

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_split_edge

  subroutine test_split_border_edge( test_name_parent)
    ! Test the split_border_edge subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'split_border_edge'
    character(len=1024), parameter :: test_name_local = 'split_border_edge'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), dimension(2)         :: p_new

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southern border edge
    p_new = (mesh%V(1,:) + mesh%V(2,:)) / 2._dp
    call split_border_edge( mesh, 1, 2, p_new)

    ! Check if there is indeed a new vertex in the location we wanted
    call unit_test( test_eq ( mesh%nV, 6), trim(test_name)//'/new_vertex')
    call unit_test( test_tol( mesh%V(6,:), p_new, 1E-10_dp), trim(test_name)//'/new_vertex_location')

    ! Check if there is indeed one new triangle
    call unit_test( test_eq ( mesh%nTri, 5), trim(test_name)//'/new_triangle')

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_split_border_edge

  subroutine test_move_vertex( test_name_parent)
    ! Test the move_vertex subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_move_vertex'
    character(len=1024), parameter :: test_name_local = 'move_vertex'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), dimension(2)         :: p_new

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Move the central vertex
    p_new = [(mesh%xmax - mesh%xmin) / 53._dp, (mesh%ymax - mesh%ymin) / 43.2_dp]
    call move_vertex( mesh, 5, p_new)

    ! Check if the vertex is in the correct location
    call unit_test( test_tol( mesh%V(5,:), p_new, 1E-10_dp), trim(test_name)//'/new_vertex_location')

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_move_vertex

end module ut_mesh_Delaunay