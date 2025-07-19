module ut_mesh_refinement

  ! Unit tests for mesh functions.

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use parameters
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform, refine_mesh_point, refine_mesh_line, refine_mesh_polygon
  use mesh_refinement_basic_ROI, only: refine_mesh_line_ROI, refine_mesh_polygon_ROI
  use mesh_utilities, only: find_containing_triangle, calc_smallest_internal_angle_mesh, calc_mean_skewness
  use plane_geometry, only: longest_triangle_leg, is_in_polygon
  use mesh_Lloyds_algorithm, only: Lloyds_algorithm_single_iteration
  use mesh_contiguous_domains, only: enforce_contiguous_process_domains
  use mesh_secondary, only: calc_triangle_geometric_centres
  use mpi_basic, only: par

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

    call test_refine_mesh_uniform    ( test_name)
    call test_refine_mesh_point      ( test_name)
    call test_refine_mesh_line       ( test_name)
    call test_refine_mesh_polygon    ( test_name)
    call test_refine_mesh_line_ROI   ( test_name)
    call test_refine_mesh_polygon_ROI( test_name)
    call test_Lloyds_algorithm       ( test_name)
    call test_contiguous_domains     ( test_name)

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
    real(dp)                       :: smallest_internal_angle
    integer                        :: vi,ci,vj
    real(dp)                       :: dist, dist_min, dist_max

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    call refine_mesh_uniform( mesh, res_max, alpha_min)

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Check if the smallest internal angle is indeed larger than alpha_min
    call calc_smallest_internal_angle_mesh( mesh, smallest_internal_angle)
    call unit_test( test_ge( smallest_internal_angle, alpha_min), trim(test_name)//'/smallest_internal_angle')

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

    call unit_test( test_ge( dist_min, res_max / 5._dp), trim(test_name)//'/minimum_resolution')
    call unit_test( test_le( dist_max, res_max), trim(test_name)//'/maximum_resolution')

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
    real(dp)                          :: smallest_internal_angle
    integer                           :: ti
    real(dp), dimension(2)            :: p, q, r
    real(dp)                          :: longest_leg

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    call refine_mesh_point( mesh, POI, res_Max, alpha_min)

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Check if the smallest internal angle is indeed larger than alpha_min
    call calc_smallest_internal_angle_mesh( mesh, smallest_internal_angle)
    call unit_test( test_ge( smallest_internal_angle, alpha_min), trim(test_name)//'/smallest_internal_angle')

    ! Check if the resolution* at the POI is indeed below, but not too far below, the specified tolerance
    ! *defined here as the longest leg of the triangle containing the POI

    ti = 1
    call find_containing_triangle( mesh, POI, ti)
    p = mesh%V( mesh%Tri( ti,1),:)
    q = mesh%V( mesh%Tri( ti,2),:)
    r = mesh%V( mesh%Tri( ti,3),:)
    longest_leg = longest_triangle_leg( p, q, r)

    call unit_test( test_le( longest_leg, res_max), trim(test_name)//'/maximum_resolution')

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
    real(dp)                       :: smallest_internal_angle
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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    p_line( 1,:) = [-0.43_dp, -0.129_dp, 0.67_dp, 0.85_dp]
    width = res_max / 2._dp
    call refine_mesh_line( mesh, p_line, res_max, width, alpha_min)

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Check if the smallest internal angle is indeed larger than alpha_min
    call calc_smallest_internal_angle_mesh( mesh, smallest_internal_angle)
    call unit_test( test_ge( smallest_internal_angle, alpha_min), trim(test_name)//'/smallest_internal_angle')

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

    call unit_test( test_le( longest_leg, res_max), trim(test_name)//'/maximum_resolution')

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
    real(dp)                       :: smallest_internal_angle
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
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    poly( 1,:) = [-0.43_dp, -0.129_dp]
    poly( 2,:) = [0.27_dp, -0.45_dp]
    poly( 3,:) = [0.67_dp, 0.16_dp]
    poly( 4,:) = [0.54_dp, 0.27_dp]
    poly( 5,:) = [-0.02_dp, 0.1_dp]
    call refine_mesh_polygon( mesh, poly, res_max, alpha_min)

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Check if the smallest internal angle is indeed larger than alpha_min
    call calc_smallest_internal_angle_mesh( mesh, smallest_internal_angle)
    call unit_test( test_ge( smallest_internal_angle, alpha_min), trim(test_name)//'/smallest_internal_angle')

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

    call unit_test( test_le( longest_leg, res_max), trim(test_name)//'/maximum_resolution')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_polygon

  subroutine test_refine_mesh_line_ROI( test_name_parent)
    ! Test the refine_mesh_line_ROI subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'test_refine_mesh_line_ROI'
    character(len=1024), parameter        :: test_name_local = 'refine_mesh_line_ROI'
    character(len=1024)                   :: test_name
    type(type_mesh)                       :: mesh
    real(dp), parameter                   :: xmin = -1._dp
    real(dp), parameter                   :: xmax =  1._dp
    real(dp), parameter                   :: ymin = -1._dp
    real(dp), parameter                   :: ymax =  1._dp
    real(dp), dimension(:,:), allocatable :: p_line
    integer                               :: n_line
    real(dp), dimension(2)                :: p_start, p_end, pp, qq
    real(dp)                              :: w1, w2
    real(dp), dimension(4,2)              :: poly_ROI
    real(dp), parameter                   :: res_max = 0.05_dp
    real(dp)                              :: width
    real(dp), parameter                   :: alpha_min = 0.4363_dp
    real(dp)                              :: smallest_internal_angle
    integer                               :: n,i
    integer                               :: ti
    real(dp), dimension(2)                :: p, q, r
    real(dp)                              :: longest_leg_in, longest_leg_out

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    n_line = 100
    allocate( p_line(n_line,4))
    p_start = [-0.5_dp, 0._dp]
    p_end   = [ 0.5_dp, 0._dp]
    do i = 1, n_line
      w1 = real( i-1,dp) / real( n_line,dp)
      w2 = real( i  ,dp) / real( n_line,dp)
      pp = (1._dp - w1) * p_start + w1 * p_end
      qq = (1._dp - w2) * p_start + w2 * p_end
      p_line( i,:) = [pp(1), pp(2), qq(1), qq(2)]
    end do

    poly_ROI( 1,:) = [0._dp, -1._dp]
    poly_ROI( 2,:) = [1._dp, -1._dp]
    poly_ROI( 3,:) = [1._dp,  1._dp]
    poly_ROI( 4,:) = [0._dp,  1._dp]

    width = res_max / 2._dp

    call refine_mesh_line_ROI( mesh, p_line, res_max, width, alpha_min, poly_ROI)

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Check if the smallest internal angle is indeed larger than alpha_min
    call calc_smallest_internal_angle_mesh( mesh, smallest_internal_angle)
    call unit_test( test_ge( smallest_internal_angle, alpha_min), trim(test_name)//'/smallest_internal_angle')

    ! Check if the resolution* along the part of the line
    ! outside of the ROI is indeed below the specified tolerance,
    ! and the part of the line outside of the ROI is indeed
    ! above it.

    n = 100
    ti = 1
    longest_leg_in  = 0._dp
    longest_leg_out = 0._dp
    do i = 1, n
      ! Walk along the line in small steps
      pp = [p_line(1,1) + (p_line(1,3) - p_line(1,1)) * (real(i-1,dp) / real(n-1,dp)), &
            p_line(1,2) + (p_line(1,4) - p_line(1,2)) * (real(i-1,dp) / real(n-1,dp))]
      ! Find the longest triangle leg along the line
      call find_containing_triangle( mesh, pp, ti)
      p = mesh%V( mesh%Tri( ti,1),:)
      q = mesh%V( mesh%Tri( ti,2),:)
      r = mesh%V( mesh%Tri( ti,3),:)
      if (is_in_polygon( poly_ROI, pp)) then
        longest_leg_in = max( longest_leg_in, longest_triangle_leg( p, q, r))
      else
        longest_leg_out = max( longest_leg_out, longest_triangle_leg( p, q, r))
      end if
    end do

    call unit_test( test_le( longest_leg_in , res_max), trim(test_name)//'/maximum_resolution_in')
    call unit_test( test_ge( longest_leg_out, res_max), trim(test_name)//'/maximum_resolution_out')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_line_ROI

  subroutine test_refine_mesh_polygon_ROI( test_name_parent)
    ! Test the refine_mesh_polygon_ROI subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_refine_mesh_polygon_ROI'
    character(len=1024), parameter :: test_name_local = 'refine_mesh_polygon_ROI'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    real(dp), dimension(5,2)       :: poly
    real(dp), dimension(4,2)       :: poly_ROI
    real(dp), parameter            :: res_max = 0.05_dp
    real(dp), parameter            :: alpha_min = 0.4363_dp
    real(dp)                       :: smallest_internal_angle
    integer                        :: n,i,j
    real(dp), dimension(2)         :: pp
    integer                        :: ti
    real(dp), dimension(2)         :: p, q, r
    real(dp)                       :: longest_leg_in, longest_leg_out

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    poly( 1,:) = [-0.43_dp, -0.129_dp]
    poly( 2,:) = [0.27_dp, -0.45_dp]
    poly( 3,:) = [0.67_dp, 0.16_dp]
    poly( 4,:) = [0.54_dp, 0.27_dp]
    poly( 5,:) = [-0.02_dp, 0.1_dp]

    poly_ROI( 1,:) = [0._dp, -1._dp]
    poly_ROI( 2,:) = [1._dp, -1._dp]
    poly_ROI( 3,:) = [1._dp,  1._dp]
    poly_ROI( 4,:) = [0._dp,  1._dp]

    call refine_mesh_polygon_ROI( mesh, poly, res_max, alpha_min, poly_ROI)

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Check if the smallest internal angle is indeed larger than alpha_min
    call calc_smallest_internal_angle_mesh( mesh, smallest_internal_angle)
    call unit_test( test_ge( smallest_internal_angle, alpha_min), trim(test_name)//'/smallest_internal_angle')

    ! Check if the resolution* inside the polygon is indeed below, but not too far below, the specified tolerance
    ! *defined here as the longest leg of the triangles

    n = 100
    ti = 1
    longest_leg_in  = 0._dp
    longest_leg_out = 0._dp
    do i = 1, n
    do j = 1, n
      pp = [mesh%xmin + (mesh%xmax - mesh%xmin) * (real(i-1,dp) / real(n-1,dp)), &
            mesh%ymin + (mesh%ymax - mesh%ymin) * (real(j-1,dp) / real(n-1,dp))]
      if (is_in_polygon( poly, pp)) then
        call find_containing_triangle( mesh, pp, ti)
        p = mesh%V( mesh%Tri( ti,1),:)
        q = mesh%V( mesh%Tri( ti,2),:)
        r = mesh%V( mesh%Tri( ti,3),:)
        if (is_in_polygon( poly_ROI, pp)) then
          longest_leg_in = max( longest_leg_in, longest_triangle_leg( p, q, r))
        else
          longest_leg_out = max( longest_leg_out, longest_triangle_leg( p, q, r))
        end if
      end if
    end do
    end do

    call unit_test( test_le( longest_leg_in , res_max), trim(test_name)//'/maximum_resolution_in')
    call unit_test( test_ge( longest_leg_out, res_max), trim(test_name)//'/maximum_resolution_out')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_refine_mesh_polygon_ROI

  subroutine test_Lloyds_algorithm( test_name_parent)
    ! Test the Lloyds_algorithm_single_iteration subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'test_Lloyds_algorithm'
    character(len=1024), parameter    :: test_name_local = 'test_Lloyds_algorithm'
    character(len=1024)               :: test_name
    type(type_mesh)                   :: mesh
    real(dp), parameter               :: xmin = -1._dp
    real(dp), parameter               :: xmax =  1._dp
    real(dp), parameter               :: ymin = -1._dp
    real(dp), parameter               :: ymax =  1._dp
    real(dp), dimension(2), parameter :: POI =  [0._dp, 0._dp]
    real(dp), parameter               :: res_max = 0.02_dp
    real(dp), parameter               :: alpha_min = 0.4363_dp
    real(dp)                          :: smallest_internal_angle
    integer,  parameter               :: nit_Lloyds_algorithm = 5
    integer                           :: it_Lloyds_algorithm
    character(len=3)                  :: it_Lloyds_algorithm_str
    real(dp)                          :: mean_skewness, mean_skewness_prev

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh around a point, to generate strong resolution gradients
    call refine_mesh_point( mesh, POI, res_max, alpha_min)

    ! Check if the refined mesh is self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency_before')

    ! Check if the smallest internal angle is indeed larger than alpha_min
    call calc_smallest_internal_angle_mesh( mesh, smallest_internal_angle)
    call unit_test( test_ge( smallest_internal_angle, alpha_min), trim(test_name)//'/smallest_internal_angle_before')

    ! Apply a few iterations of Lloyds algorithm, and check if the mean skewness does indeed decrease

    call calc_mean_skewness( mesh, mean_skewness)

    do it_Lloyds_algorithm = 1, nit_Lloyds_algorithm

      write( it_Lloyds_algorithm_str,'(i3)') it_Lloyds_algorithm

      call Lloyds_algorithm_single_iteration( mesh, alpha_min)

      call unit_test( test_mesh_is_self_consistent( mesh), &
        trim(test_name)//'/mesh_self_consistency_it_'//trim(adjustl(it_Lloyds_algorithm_str)))

      mean_skewness_prev = mean_skewness
      call calc_mean_skewness( mesh, mean_skewness)
      call unit_test( test_le( mean_skewness, mean_skewness_prev), &
        trim(test_name)//'/mean_skewness_it_'//trim(adjustl(it_Lloyds_algorithm_str)))

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_Lloyds_algorithm

  subroutine test_contiguous_domains( test_name_parent)
    ! Test the enforce_contiguous_process_domains subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_contiguous_domains'
    character(len=1024), parameter :: test_name_local = 'contiguous_domains'
    character(len=1024)            :: test_name
    type(type_mesh)                :: mesh
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    real(dp), parameter            :: res_max = 0.05_dp
    real(dp), parameter            :: alpha_min = 0.4363_dp
    logical                        :: vertices_are_sorted, triangles_are_sorted
    integer                        :: vi, ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(test_name)//'_mesh', 1000, 2000)

    ! Initialise dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    call refine_mesh_uniform( mesh, res_max, alpha_min)

    ! Reorder the vertices and triangles to enforce contiguous process domains
    call enforce_contiguous_process_domains( mesh)

    ! Calculate triangle geometric centres (needed to see if they are sorted)
    call calc_triangle_geometric_centres( mesh)

    ! Check if the mesh is still self-consistent
    call unit_test( test_mesh_is_self_consistent( mesh), trim(test_name)//'/mesh_self_consistency')

    ! Check if the vertices and triangles really have been reordered properly
    vertices_are_sorted = .true.
    do vi = 2, mesh%nV
      vertices_are_sorted = vertices_are_sorted .and. mesh%V( vi,1) >= mesh%V( vi-1,1) - mesh%tol_dist
    end do
    call unit_test( vertices_are_sorted, trim(test_name)//'/vertices_are_sorted')

    triangles_are_sorted = .true.
    do ti = 2, mesh%nTri
      triangles_are_sorted = triangles_are_sorted .and. mesh%Trigc( ti,1) >= mesh%Trigc( ti-1,1) - mesh%tol_dist
    end do
    call unit_test( triangles_are_sorted, trim(test_name)//'/triangles_are_sorted')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_contiguous_domains

end module ut_mesh_refinement