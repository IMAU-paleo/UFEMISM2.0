module mesh_refinement_basic_ROI

  ! The basic mesh refinement routines: uniform (dimensionless), point (0D), line (1D), and poly( 2D)
  ! ...but only within a specific Region Of Interest.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use grid_basic, only: poly2line
  use plane_geometry, only: is_in_polygon, crop_line_to_domain, smallest_triangle_angle, &
    is_in_triangle, segment_intersection, lies_on_line_segment, longest_triangle_leg, circumcenter
  use mesh_utilities, only: find_containing_triangle, add_triangle_to_refinement_stack_last, &
    remove_triangle_from_refinement_stack
  use mesh_memory, only: extend_mesh_primary, crop_mesh_primary
  use split_triangles, only: split_triangle
  use mesh_refinement_basic, only: refine_mesh_split_encroaching_triangles

  implicit none

  private

  public :: refine_mesh_line_ROI, refine_mesh_polygon_ROI

contains

  subroutine refine_mesh_line_ROI( mesh, p_line_tot, res_max, width, alpha_min, poly_ROI)
    ! Refine a mesh based on a 1-D line criterion
    !
    ! Only consider line segments lying inside the region of interest
    ! described by the provided polygon

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh          ! The mesh that should be refined
    real(dp), dimension(:,:  ), intent(in)        :: p_line_tot    ! Collection of line segments
    real(dp),                   intent(in)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    real(dp),                   intent(in)        :: width         ! Maximum allowed resolution for triangles crossed by any of these line segments
    real(dp),                   intent(in)        :: alpha_min     ! minimum allowed internal triangle angle
    real(dp), dimension(:,:  ), intent(in)        :: poly_ROI      ! Polygon describing the region of interest

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'refine_mesh_line_ROI'
    integer                                       :: nl_tot
    real(dp), dimension(:,:  ), allocatable       :: p_line
    integer                                       :: nl
    integer                                       :: ti,li
    real(dp), dimension(2)                        :: pp,qq,pp2,qq2,pp_cropped,qq_cropped,dd
    logical                                       :: is_valid_line
    integer                                       :: tip, tiq, via, vib, vic, tip2, tiq2
    real(dp), dimension(2)                        :: va,vb,vc,llis
    logical                                       :: do_cross, crosses
    integer                                       :: li_min, li_max
    real(dp)                                      :: longest_leg, smallest_angle
    logical                                       :: meets_resolution_criterion
    logical                                       :: meets_geometry_criterion
    real(dp), dimension(2)                        :: p_new

    ! Add routine to path
    call init_routine( routine_name)

    nl_tot = size( p_line_tot,1)

    ! If no line is provided, do nothing
    if (nl_tot == 0) then
      call finalise_routine( routine_name)
      return
    end if

    ! Safety
    if (size( p_line_tot,2) /= 4) call crash('line must be an n-by-4 array!')

    ! == Only consider line segments lying inside the region of interest
    !    described by the provided polygon
    ! ====================================

    allocate( p_line( nl_tot,4), source = 0._dp)
    nl = 0

    do li = 1, nl_tot

      ! Line endpoints
      pp = [p_line_tot( li,1), p_line_tot( li,2)]
      qq = [p_line_tot( li,3), p_line_tot( li,4)]

      ! If this line lies inside the region of interest, consider it
      if (is_in_polygon( poly_ROI, pp) .or. is_in_polygon( poly_ROI, qq)) then
        nl = nl + 1
        p_line( nl,:) = p_line_tot( li,:)
      end if

    end do ! do li = 1, nl_tot

    ! == Initialise triangle-line overlap ranges ==
    ! =============================================

    ! The array Tri_li has nTri-by-2 elements, and describes the range of line segments
    ! crossing through each triangle. E.g., if Tri_li( 5,:) = [175,179], that means that
    ! triangle 5 is crosses by line segments 175 to 179.
    ! Since the lines (grounding line, calving front, etc.) are constructed in such a way
    ! that the vast majority of line segments are contiguous, generally a triangle will
    ! only overlap with a small range of them.

    mesh%Tri_li = 0
    mesh%Tri_li( 1:mesh%nTri,1) = nl + 1
    mesh%Tri_li( 1:mesh%nTri,2) = 0

    tip = 1
    tiq = 1

    do li = 1, nl

      ! Line endpoints
      pp = [p_line( li,1), p_line( li,2)]
      qq = [p_line( li,3), p_line( li,4)]

      ! Crop line to mesh domain
      call crop_line_to_domain( pp, qq, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
      pp = pp_cropped
      qq = qq_cropped

      ! If the line segment lies outside of the mesh domain altogether, skip it
      if (.not. is_valid_line) cycle

      ! If they coincide, this line is invalid - skip
      if (norm2( pp - qq) < mesh%tol_dist) cycle

      ! Find the triangles containing p and q
      tip = tiq
      call find_containing_triangle( mesh, pp, tip)
      tiq = tip
      call find_containing_triangle( mesh, qq, tiq)

      ! Add li to the overlap lists of tip and tiq
      mesh%Tri_li( tip,1) = min( mesh%Tri_li( tip,1),li)
      mesh%Tri_li( tip,2) = max( mesh%Tri_li( tip,2),li)
      mesh%Tri_li( tiq,1) = min( mesh%Tri_li( tiq,1),li)
      mesh%Tri_li( tiq,2) = max( mesh%Tri_li( tiq,2),li)

      ! If they both lie inside the same triangle, no need to trace
      if (tip == tiq) cycle

      ! p and q lie in different triangles - perform a pseudo-trace
      do while (norm2( pp - qq) > res_max)

        ! Move sideways to include line width
        pp2 = pp - [(pp( 2) - qq( 2)), qq( 1) - pp( 1)]
        qq2 = pp + [(pp( 2) - qq( 2)), qq( 1) - pp( 1)]

        call crop_line_to_domain( pp2, qq2, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
        pp2 = pp_cropped
        qq2 = qq_cropped

        if (is_valid_line) then

          tip2 = tiq
          call find_containing_triangle( mesh, pp, tip2)
          tiq2 = tip2
          call find_containing_triangle( mesh, qq, tiq2)

          do while (norm2 ( pp2 - qq2) > res_max)

            ! Move pp in the direction of qq
            dd = qq2 - pp2
            dd = dd / norm2( dd)
            dd = dd * res_max
            pp2 = pp2 + dd

            ! Find the triangle that now contains pp
            call find_containing_triangle( mesh, pp2, tip2)

            ! Add li to the overlap lists of tip
            mesh%Tri_li( tip2,1) = min( mesh%Tri_li( tip2,1),li)
            mesh%Tri_li( tip2,2) = max( mesh%Tri_li( tip2,2),li)

            ! If we've reached tiq, the trace is done
            if (tip2 == tiq2) exit

          end do ! do while (norm2 ( pp2 - qq2) > res_max)

        end if ! if (is_valid_line) then

        ! Move along pq

        ! Move pp in the direction of qq
        dd = qq - pp
        dd = dd / norm2( dd)
        dd = dd * res_max
        pp = pp + dd

        ! Find the triangle that now contains pp
        call find_containing_triangle( mesh, pp, tip)

        ! Add li to the overlap lists of tip
        mesh%Tri_li( tip,1) = min( mesh%Tri_li( tip,1),li)
        mesh%Tri_li( tip,2) = max( mesh%Tri_li( tip,2),li)

        ! If we've reached tiq, the trace is done
        if (tip == tiq) exit

      end do ! do while (tip /= tiq)

    end do ! do li = 1, nl

    ! == Iteratively refine the mesh ==
    ! =================================

    ! If a triangle overlaps with any of the line segments, check if it is small enough.
    ! If not, split it. The new triangles will inherit the old one's line segment range.
    ! Update the (now reduced) overlap range for the new triangles.

    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0

    ! Mark which triangles need to be refined right now
    do ti = 1, mesh%nTri
      if (mesh%Tri_li( ti,2) == 0) then
        ! This triangle does not overlap with any line segments, so it doesn't need refining
        cycle
      else
        ! Mark this triangle for refinement
        call add_triangle_to_refinement_stack_last( mesh, ti)
      end if
    end do

    ! Keep refining until all triangles match the criterion
    do while (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      if (mesh%nV > mesh%nV_mem - 10 .or. mesh%nTri > mesh%nTri_mem - 10) then
        call extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      end if

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      call remove_triangle_from_refinement_stack( mesh, ti)

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      va  = mesh%V( via,:)
      vb  = mesh%V( vib,:)
      vc  = mesh%V( vic,:)

      ! Check if it meets the geometry criterion

      smallest_angle = smallest_triangle_angle( va, vb, vc)
      meets_geometry_criterion = smallest_angle >= alpha_min

      if (.not. meets_geometry_criterion) then
        ! This triangle should be split anyway; no need to check
        ! the resolution criterion, so we can save some time)

        meets_resolution_criterion = .true.

      else
        ! This triangle meets the geometry criterion; check if
        ! it also meets the resolution criterion

        ! Initial guess for this triangle's line overlap range (taken from their parent)
        li_min = mesh%Tri_li( ti,1)
        li_max = mesh%Tri_li( ti,2)

        ! If that's already zero, skip
        if (li_min == nl+1 .or. li_max == 0) then
          ! No line overlap anyway

          meets_resolution_criterion = .true.

        else  ! if (li_min == nl+1 .or. li_max == 0) then
          ! Recalculate triangle overlap range

          mesh%Tri_li( ti,:) = [nl+1,0]

          do li = li_min, li_max

            ! Line endpoints
            pp = [p_line( li,1), p_line( li,2)]
            qq = [p_line( li,3), p_line( li,4)]

            ! Crop line to mesh domain
            call crop_line_to_domain( pp, qq, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
            pp = pp_cropped
            qq = qq_cropped

            ! If the line segment lies outside of the mesh domain altogether, skip it
            if (.not. is_valid_line) cycle

            ! If they coincide, this line is invalid - skip
            if (norm2( pp - qq) < mesh%tol_dist) cycle

            ! Check if this line segment crosses triangle ti
            crosses = .false.
            crosses = crosses .or. is_in_triangle( va, vb, vc, pp)
            crosses = crosses .or. is_in_triangle( va, vb, vc, qq)
            call segment_intersection( pp, qq, va, vb, llis, do_cross, mesh%tol_dist)
            crosses = crosses .or. do_cross
            call segment_intersection( pp, qq, vb, vc, llis, do_cross, mesh%tol_dist)
            crosses = crosses .or. do_cross
            call segment_intersection( pp, qq, vc, va, llis, do_cross, mesh%tol_dist)
            crosses = crosses .or. do_cross
            crosses = crosses .or. lies_on_line_segment( pp, qq, va, width)
            crosses = crosses .or. lies_on_line_segment( pp, qq, vb, width)
            crosses = crosses .or. lies_on_line_segment( pp, qq, vc, width)

            ! If so, add it to the overlap range
            if (crosses) then
              mesh%Tri_li( ti,1) = min( mesh%Tri_li( ti,1),li)
              mesh%Tri_li( ti,2) = max( mesh%Tri_li( ti,2),li)
            end if

          end do ! do li = li_min, li_max

          ! If this triangle overlaps with any line segments,
          ! check if it meets the resolution criterion

          meets_resolution_criterion = .true.
          if (mesh%Tri_li( ti,1) <= mesh%Tri_li( ti,2)) then
            longest_leg = longest_triangle_leg( va, vb, vc)
            if (longest_leg > res_max * mesh%resolution_tolerance) meets_resolution_criterion = .false.
          end if

        end if ! if (li_min == nl+1 .or. li_max == 0) then

      end if ! if (.not. meets_geometry_criterion) then

      ! If either of the two criteria is not met, split the triangle
      if (.not. meets_geometry_criterion .or. .not. meets_resolution_criterion) then
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        call split_triangle( mesh, ti, p_new)
      end if

    end do ! do while (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    call refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    call crop_mesh_primary( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine refine_mesh_line_ROI

  subroutine refine_mesh_polygon_ROI( mesh, poly, res_max, alpha_min, poly_ROI, poly_mult_not)
    ! Refine a mesh based on a 2-D polygon criterion
    !
    ! Only consider triangles lying inside the region of interest
    ! described by the provided polygon

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh          ! The mesh that should be refined
    real(dp), dimension(:,:),           intent(in)    :: poly          ! Polygon of interest
    real(dp),                           intent(in)    :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    real(dp),                           intent(in)    :: alpha_min     ! minimum allowed internal triangle angle
    real(dp), dimension(:,:  ),         intent(in)    :: poly_ROI      ! Polygon describing the region of interest
    real(dp), dimension(:,:), optional, intent(in)    :: poly_mult_not  ! Exclude triangles overlapping with these set of polygons

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'refine_mesh_polygon_ROI'
    real(dp), dimension(:,:  ), allocatable           :: p_line
    integer                                           :: ti, via, vib, vic
    real(dp), dimension(2)                            :: va, vb, vc
    logical                                           :: has_any_overlap, is_in_poly_not
    real(dp)                                          :: longest_leg, smallest_angle
    logical                                           :: meets_resolution_criterion
    logical                                           :: meets_geometry_criterion
    real(dp), dimension(2)                            :: p_new
    real(dp), dimension(:,:  ), allocatable           :: poly_not
    integer                                           :: n1,nn,n2

    ! Add routine to path
    call init_routine( routine_name)

    ! First refine the mesh along the polygon perimeter by treating it as a 1-D line
    call poly2line( poly, p_line)
    call refine_mesh_line_ROI( mesh, p_line, res_max, res_max, alpha_min, poly_ROI)

    ! Initialise the refinement stack with all triangles lying (partly) inside the polygon
    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0
    has_any_overlap = .false.
    do ti = 1, mesh%nTri

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      va  = mesh%V( via,:)
      vb  = mesh%V( vib,:)
      vc  = mesh%V( vic,:)

      ! Check that none of the vertices lie within the
      ! secondary do-not-refine-here polygons
      n1 = 1
      n2 = 0
      is_in_poly_not = .false.

      do while (present(poly_mult_not) .and. n2 < size( poly_mult_not,1))

        ! Copy a single polygon from poly_mult
        nn = NinT( poly_mult_not( n1,1))
        n2 = n1 + nn
        allocate( poly_not( nn,2))
        poly_not = poly_mult_not( n1+1:n2,:)
        n1 = n2+1

        if (is_in_polygon( poly_not, va) .and. &
            is_in_polygon( poly_not, vb) .and. &
            is_in_polygon( poly_not, vc)) then
          ! Triangle lies within the no-remesh zone
          is_in_poly_not = .true.
        end if

        ! Clean up after yourself
        deallocate( poly_not)

        ! If we already know, stop checking
        if (is_in_poly_not) exit

      end do ! do while (n2 < size( poly_mult_not,1))

      ! Triangle overlaps with the no-remesh zone, so skip it
      if (is_in_poly_not) cycle

      ! If this triangle lies (partly) inside the polygon, mark it for refinement
      if ((is_in_polygon( poly    , va) .or. is_in_polygon( poly    , vb) .or. is_in_polygon( poly    , vc)) .and. &
          (is_in_polygon( poly_ROI, va) .or. is_in_polygon( poly_ROI, vb) .or. is_in_polygon( poly_ROI, vc))) then
        has_any_overlap = .true.
        call add_triangle_to_refinement_stack_last( mesh, ti)
      end if

    end do ! do ti = 1, mesh%nTri

    ! Safety: if the polygon doesn't overlap with the mesh
    !         at all, no refinement is needed (or indeed possible).
    if (.not. has_any_overlap) then
      call finalise_routine( routine_name)
      return
    end if

    ! Keep refining until all triangles match the criterion
    do while (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      if (mesh%nV > mesh%nV_mem - 10 .or. mesh%nTri > mesh%nTri_mem - 10) then
        call extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      end if

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      call remove_triangle_from_refinement_stack( mesh, ti)

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      va  = mesh%V( via,:)
      vb  = mesh%V( vib,:)
      vc  = mesh%V( vic,:)

      ! Check if it meets the geometry criterion
      smallest_angle = smallest_triangle_angle( va, vb, vc)
      meets_geometry_criterion = smallest_angle >= alpha_min

      ! If it does not meet the geometry criterion, it will be split anyway -
      ! no need to check the resolution criterion
      if (.not. meets_geometry_criterion) then
        meets_resolution_criterion = .true.
      else
        ! Check if the triangle meets the resolution criterion

        ! Check if the triangle lies inside the polygon
        if ((is_in_polygon( poly    , va) .or. is_in_polygon( poly    , vb) .or. is_in_polygon( poly    , vc)) .and. &
            (is_in_polygon( poly_ROI, va) .or. is_in_polygon( poly_ROI, vb) .or. is_in_polygon( poly_ROI, vc))) then
          ! The triangle lies inside the polygon

          longest_leg = longest_triangle_leg( va, vb, vc)
          meets_resolution_criterion = longest_leg <= res_max * mesh%resolution_tolerance

        else
          ! The triangle does not lie inside the polygon
          meets_resolution_criterion = .true.
        end if

      end if ! if (.not. meets_geometry_criterion) then

      ! If either of the two criteria is not met, split the triangle
      if (.not. meets_geometry_criterion .or. .not. meets_resolution_criterion) then
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        call split_triangle( mesh, ti, p_new)
      end if

    end do ! do while (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    call refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    call crop_mesh_primary( mesh)

    ! Clean up after yourself
    deallocate( p_line)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine refine_mesh_polygon_ROI

end module mesh_refinement_basic_ROI
