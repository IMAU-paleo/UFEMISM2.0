module unit_tests_mesh

  ! Unit tests for mesh functions.

  use precisions, only: dp
  use parameters
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use assertions_unit_tests
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_creation, only: initialise_dummy_mesh
  use mesh_Delaunay, only: split_triangle, split_edge, split_border_edge, move_vertex
  use mesh_refinement, only: refine_mesh_uniform

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

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_mesh_main

  ! ===== Delaunay triangulation =====
  ! ==================================

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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southern triangle
    p_new = [0._dp, -1._dp]
    call split_triangle( mesh, 1, p_new)

    ! Check if there is indeed a new vertex in the location we wanted
    call test_eq ( mesh%nV    , 6    , UNIT_TEST, trim(test_name)//'/new_vertex')
    call test_tol( mesh%V(6,:), p_new, 1E-10_dp, UNIT_TEST, trim(test_name)//'/new_vertex_location')

    ! Check if there is indeed one new triangle
    call test_eq ( mesh%nTri  , 5    , UNIT_TEST, trim(test_name)//'/new_triangle')

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southwestern edge
    p_new = (mesh%V(1,:) + mesh%V(5,:)) / 2._dp
    call split_edge( mesh, 1, 5, p_new)

    ! Check if there is indeed a new vertex in the location we wanted
    call test_eq ( mesh%nV    , 6    , UNIT_TEST, trim(test_name)//'/new_vertex')
    call test_tol( mesh%V(6,:), p_new, 1E-10_dp, UNIT_TEST, trim(test_name)//'/new_vertex_location')

    ! Check if there are indeed two new triangles
    call test_eq ( mesh%nTri  , 6    , UNIT_TEST, trim(test_name)//'/new_triangles')

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southern border edge
    p_new = (mesh%V(1,:) + mesh%V(2,:)) / 2._dp
    call split_border_edge( mesh, 1, 2, p_new)

    ! Check if there is indeed a new vertex in the location we wanted
    call test_eq ( mesh%nV    , 6    , UNIT_TEST, trim(test_name)//'/new_vertex')
    call test_tol( mesh%V(6,:), p_new, 1E-10_dp, UNIT_TEST, trim(test_name)//'/new_vertex_location')

    ! Check if there is indeed one new triangle
    call test_eq ( mesh%nTri  , 5    , UNIT_TEST, trim(test_name)//'/new_triangle')

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Move the central vertex
    p_new = [(mesh%xmax - mesh%xmin) / 53._dp, (mesh%ymax - mesh%ymin) / 43.2_dp]
    call move_vertex( mesh, 5, p_new)

    ! Check if the vertex is in the correct location
    call test_tol( mesh%V(5,:), p_new, 1E-10_dp, UNIT_TEST, trim(test_name)//'/new_vertex_location')

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_move_vertex

  ! ===== Refinement =====
  ! ======================

  subroutine test_refinement( test_name_parent)
    ! Test the mesh refinement subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_refinement'
    character(len=1024), parameter :: test_name_local = 'refinement'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_refine_mesh_uniform( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refinement

  subroutine test_refine_mesh_uniform( test_name_parent)
    ! Test the refine_mesh_uniform subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_refine_mesh_uniform'
    character(len=1024), parameter :: test_name_local = 'refine_mesh_uniform'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), parameter            :: xmin = -3040E3_dp  ! Just use the standard Antarctica domain; doesn't really matter here...
    real(dp), parameter            :: xmax =  3040E3_dp
    real(dp), parameter            :: ymin = -3040E3_dp
    real(dp), parameter            :: ymax =  3040E3_dp
    real(dp)                       :: res_max, alpha_min

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max   = 400E3_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_uniform

end module unit_tests_mesh