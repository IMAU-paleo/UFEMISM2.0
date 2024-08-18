module unit_tests_mesh

  ! Unit tests for mesh functions.

  use mpi
  use precisions, only: dp
  use mpi_basic, only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use assertions_unit_tests, only: ASSERTION, UNIT_TEST, test_eqv, test_neqv, test_eq, test_neq, test_gt, test_lt, test_ge, test_le, test_ge_le, test_tol, test_eq_permute
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_creation, only: initialise_dummy_mesh
  use mesh_Delaunay, only: split_triangle, split_edge, split_border_edge

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
    call test_Delaunay( test_name)

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

    call test_split_triangle( test_name)
    call test_split_edge    ( test_name)

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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200, 300)

    ! Initialise dummy mesh
    call initialise_dummy_mesh( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southern triangle
    p_new = [0._dp, -1._dp]
    call split_triangle( mesh, 1, p_new)

    ! Check results

    ! Vertex data
    call test_eq        ( mesh%nV         , 6                             , UNIT_TEST, trim(test_name)//'/nV')
    call test_tol       ( mesh%V(6,:)     , [0._dp, -1._dp]     , 1E-10_dp, UNIT_TEST, trim(test_name)//'/V')
    call test_eq        ( mesh%nC(6)      , 3                             , UNIT_TEST, trim(test_name)//'/nC')
    call test_eq        ( mesh%C(6,1:3)   , [2,5,1]                       , UNIT_TEST, trim(test_name)//'/C')
    call test_eq        ( mesh%niTri(6)   , 2                             , UNIT_TEST, trim(test_name)//'/niTri')
    call test_eq        ( mesh%iTri(6,1:2), [5,1]                         , UNIT_TEST, trim(test_name)//'/iTri')

    ! Triangle data
    call test_eq        ( mesh%nTri       , 5                             , UNIT_TEST, trim(test_name)//'/nTri')

    call test_eq_permute( mesh%Tri(1,:)   , [1,6,5]                       , UNIT_TEST, trim(test_name)//'/Tri1')
    call test_tol       ( mesh%Tricc(1,:) , [-0.5_dp, -0.5_dp]  , 1E-2_dp , UNIT_TEST, trim(test_name)//'/Tricc1')
    call test_eq_permute( mesh%TriC(1,:)  , [5,4,0]                       , UNIT_TEST, trim(test_name)//'/TriC1')

    call test_eq_permute( mesh%Tri(5,:)   , [6,2,5]                       , UNIT_TEST, trim(test_name)//'/Tri5')
    call test_tol       ( mesh%Tricc(5,:) , [0.5_dp, -0.5_dp]   , 1E-2_dp , UNIT_TEST, trim(test_name)//'/Tricc5')
    call test_eq_permute( mesh%TriC(5,:)  , [2,1,0]                       , UNIT_TEST, trim(test_name)//'/TriC5')

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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 100, 200, 300)

    ! Initialise dummy mesh
    call initialise_dummy_mesh( mesh, -1._dp, 1._dp, -1._dp, 1._dp)

    ! Split the southwestern edge
    p_new = (mesh%V(1,:) + mesh%V(5,:)) / 2._dp
    call split_edge( mesh, 1, 5, p_new)

    ! Check results

    ! Vertex data
    call test_eq        ( mesh%nV         , 6                             , UNIT_TEST, trim(test_name)//'/nV')
    call test_tol       ( mesh%V(6,:)     , p_new               , 1E-10_dp, UNIT_TEST, trim(test_name)//'/V')
    call test_eq        ( mesh%nC(6)      , 4                             , UNIT_TEST, trim(test_name)//'/nC')
    call test_eq_permute( mesh%C(6,1:4)   , [1,2,5,4]                     , UNIT_TEST, trim(test_name)//'/C')
    call test_eq        ( mesh%niTri(6)   , 4                             , UNIT_TEST, trim(test_name)//'/niTri')
    ! call test_eq_permute( mesh%iTri(6,1:4), [1,5,4,6]                     , UNIT_TEST, trim(test_name)//'/iTri')

    ! ! Triangle data
    ! call test_eq        ( mesh%nTri       , 5                             , UNIT_TEST, trim(test_name)//'/nTri')

    ! call test_eq_permute( mesh%Tri(1,:)   , [1,6,5]                       , UNIT_TEST, trim(test_name)//'/Tri1')
    ! call test_tol       ( mesh%Tricc(1,:) , [-0.5_dp, -0.5_dp]  , 1E-2_dp , UNIT_TEST, trim(test_name)//'/Tricc1')
    ! call test_eq_permute( mesh%TriC(1,:)  , [5,4,0]                       , UNIT_TEST, trim(test_name)//'/TriC1')

    ! call test_eq_permute( mesh%Tri(5,:)   , [6,2,5]                       , UNIT_TEST, trim(test_name)//'/Tri5')
    ! call test_tol       ( mesh%Tricc(5,:) , [0.5_dp, -0.5_dp]   , 1E-2_dp , UNIT_TEST, trim(test_name)//'/Tricc5')
    ! call test_eq_permute( mesh%TriC(5,:)  , [2,1,0]                       , UNIT_TEST, trim(test_name)//'/TriC5')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_split_edge

end module unit_tests_mesh