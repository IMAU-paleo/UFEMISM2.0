module unit_tests_mesh_refinement

  ! Unit tests for mesh functions.

  use precisions, only: dp
  use parameters
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use assertions_unit_tests
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform, refine_mesh_point, refine_mesh_line, refine_mesh_polygon
  use mesh_utilities, only: find_containing_triangle, write_mesh_to_text_file
  use math_utilities, only: longest_triangle_leg, is_in_polygon

  implicit none

  private

  public :: test_refinement

contains

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
    call test_refine_mesh_point  ( test_name)
    call test_refine_mesh_line   ( test_name)
    call test_refine_mesh_polygon( test_name)

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
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    real(dp), parameter            :: res_max = 0.05_dp
    real(dp), parameter            :: alpha_min = 0.4363_dp
    integer                        :: vi,ci,vj
    real(dp)                       :: dist, dist_min, dist_max

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    C%mesh_resolution_tolerance = 1._dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

    ! Check if the resolution* is indeed uniformly below, but not too far below, the specified tolerance
    ! *defined here as the nearest-neighbour distance

    dist_min = hypot( mesh%xmax - mesh%xmin, mesh%ymax - mesh%ymin)
    dist_max = 0._dp
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        dist = hypot( mesh%V( vj,1) - mesh%V( vi,1), mesh%V( vj,2) - mesh%V( vi,2))
        dist_min = min( dist_min, dist)
        dist_max = max( dist_max, dist)
      end do
    end do

    call test_gt( dist_min, res_max / 5._dp, UNIT_TEST, trim(test_name)//'/minimum_resolution')
    call test_lt( dist_max, res_max, UNIT_TEST, trim(test_name)//'/maximum_resolution')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_uniform

  subroutine test_refine_mesh_point( test_name_parent)
    ! Test the refine_mesh_point subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'test_refine_mesh_point'
    character(len=1024), parameter    :: test_name_local = 'refine_mesh_point'
    character(len=1024)               :: test_name
    type(type_mesh)                   :: mesh
    real(dp), parameter               :: xmin = -1._dp
    real(dp), parameter               :: xmax =  1._dp
    real(dp), parameter               :: ymin = -1._dp
    real(dp), parameter               :: ymax =  1._dp
    real(dp), dimension(2), parameter :: POI =  [0.23_dp, -0.129_dp]
    real(dp), parameter               :: res_max = 0.05_dp
    real(dp), parameter               :: alpha_min = 0.4363_dp
    integer                           :: ti
    real(dp), dimension(2)            :: p, q, r
    real(dp)                          :: longest_leg

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    C%mesh_resolution_tolerance = 1._dp
    call refine_mesh_point( mesh, POI, res_Max, alpha_min)

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

    ! Check if the resolution* at the POI is indeed below, but not too far below, the specified tolerance
    ! *defined here as the longest leg of the triangle containing the POI

    ti = 1
    call find_containing_triangle( mesh, POI, ti)
    p = mesh%V( mesh%Tri( ti,1),:)
    q = mesh%V( mesh%Tri( ti,2),:)
    r = mesh%V( mesh%Tri( ti,3),:)
    longest_leg = longest_triangle_leg( p, q, r)

    call test_lt( longest_leg, res_max, UNIT_TEST, trim(test_name)//'/maximum_resolution')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_point

  subroutine test_refine_mesh_line( test_name_parent)
    ! Test the refine_mesh_line subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_refine_mesh_line'
    character(len=1024), parameter :: test_name_local = 'refine_mesh_line'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    real(dp), dimension(1,4)       :: p_line
    real(dp), parameter            :: res_max = 0.05_dp
    real(dp)                       :: width = 0.03_dp
    real(dp), parameter            :: alpha_min = 0.4363_dp
    integer                        :: n,i
    real(dp), dimension(2)         :: pp
    integer                        :: ti
    real(dp), dimension(2)         :: p, q, r
    real(dp)                       :: longest_leg

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    C%mesh_resolution_tolerance = 1._dp
    p_line( 1,:) = [-0.43_dp, -0.129_dp, 0.67_dp, 0.85_dp]
    width = res_max / 2._dp
    call refine_mesh_line( mesh, p_line, res_max, width, alpha_min)

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

    ! Check if the resolution* along the line is indeed below, but not too far below, the specified tolerance
    ! *defined here as the longest leg of the triangle containing the POI

    n = 100
    ti = 1
    longest_leg = 0._dp
    do i = 1, n
      ! Walk along the line in small steps
      pp = [p_line(1,1) + (p_line(1,3) - p_line(1,1)) * (real(i-1,dp) / real(n-1,dp)), &
            p_line(1,2) + (p_line(1,4) - p_line(1,2)) * (real(i-1,dp) / real(n-1,dp))]
      ! Find the longest triangle leg along the line
      call find_containing_triangle( mesh, pp, ti)
      p = mesh%V( mesh%Tri( ti,1),:)
      q = mesh%V( mesh%Tri( ti,2),:)
      r = mesh%V( mesh%Tri( ti,3),:)
      longest_leg = max( longest_leg, longest_triangle_leg( p, q, r))
    end do

    call test_lt( longest_leg, res_max, UNIT_TEST, trim(test_name)//'/maximum_resolution')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_line

  subroutine test_refine_mesh_polygon( test_name_parent)
    ! Test the refine_mesh_polygon subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_refine_mesh_polygon'
    character(len=1024), parameter :: test_name_local = 'refine_mesh_polygon'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    real(dp), dimension(5,2)       :: poly
    real(dp), parameter            :: res_max = 0.05_dp
    real(dp), parameter            :: alpha_min = 0.4363_dp
    integer                        :: n,i,j
    real(dp), dimension(2)         :: pp
    integer                        :: ti
    real(dp), dimension(2)         :: p, q, r
    real(dp)                       :: longest_leg

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000, 32)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    C%mesh_resolution_tolerance = 1._dp
    poly( 1,:) = [-0.43_dp, -0.129_dp]
    poly( 2,:) = [0.27_dp, -0.45_dp]
    poly( 3,:) = [0.67_dp, 0.16_dp]
    poly( 4,:) = [0.54_dp, 0.27_dp]
    poly( 5,:) = [-0.02_dp, 0.1_dp]
    call refine_mesh_polygon( mesh, poly, res_max, alpha_min)

    ! Check if the mesh is still self-consistent
    call test_mesh_is_self_consistent( mesh, UNIT_TEST, trim(test_name)//'/mesh_self_consistency')

    ! Check if the resolution* inside the polygon is indeed below, but not too far below, the specified tolerance
    ! *defined here as the longest leg of the triangles

    n = 100
    ti = 1
    longest_leg = 0._dp
    do i = 1, n
    do j = 1, n
      pp = [mesh%xmin + (mesh%xmax - mesh%xmin) * (real(i-1,dp) / real(n-1,dp)), &
            mesh%ymin + (mesh%ymax - mesh%ymin) * (real(j-1,dp) / real(n-1,dp))]
      if (is_in_polygon( poly, pp)) then
        call find_containing_triangle( mesh, pp, ti)
        p = mesh%V( mesh%Tri( ti,1),:)
        q = mesh%V( mesh%Tri( ti,2),:)
        r = mesh%V( mesh%Tri( ti,3),:)
        longest_leg = max( longest_leg, longest_triangle_leg( p, q, r))
      end if
    end do
    end do

    call test_lt( longest_leg, res_max, UNIT_TEST, trim(test_name)//'/maximum_resolution')

    call write_mesh_to_text_file( mesh, trim(C%output_dir)//'/test_mesh.txt')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_polygon

end module unit_tests_mesh_refinement