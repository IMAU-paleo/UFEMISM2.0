MODULE mesh_refinement

  ! Routines used to refine a mesh.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate
  USE math_utilities                                         , ONLY: segment_intersection, is_in_triangle, longest_triangle_leg, smallest_triangle_angle, &
                                                                     circumcenter, lies_on_line_segment, crop_line_to_domain, geometric_center, is_in_polygon, &
                                                                     cross2, quick_n_dirty_sort
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: extend_mesh_primary, crop_mesh_primary
  USE mesh_utilities                                         , ONLY: update_triangle_circumcenter, find_containing_triangle, add_triangle_to_refinement_stack_last, &
                                                                     remove_triangle_from_refinement_stack
  USE mesh_Delaunay                                          , ONLY: split_triangle, move_vertex
  USE grid_basic                                             , ONLY: poly2line

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Mesh refinement based on different criteria

  SUBROUTINE refine_mesh_uniform( mesh, res_max, alpha_min)
    ! Refine a mesh to a uniform maximum resolution

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_uniform'
    INTEGER                                       :: ti, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the refinement stack with all triangles
    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0
    DO ti = 1, mesh%nTri
      CALL add_triangle_to_refinement_stack_last( mesh, ti)
    END DO

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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

      ! Check if it meets the resolution criterion
      longest_leg = longest_triangle_leg( va, vb, vc)
      meets_resolution_criterion = longest_leg <= res_max * C%mesh_resolution_tolerance

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_resolution_criterion) THEN
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    CALL refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_uniform

  SUBROUTINE refine_mesh_point( mesh, POI, res_max, alpha_min)
    ! Refine a mesh based on a 0-D point criterion

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: POI           ! [m] x,y-coordinates of point of interest where the mesh should be refined
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_point'
    INTEGER                                       :: ti_in
    INTEGER                                       :: ti, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety: if the POI lies outside the mesh domain,
    !         no refinement is needed (or indeed possible).
    IF (POI( 1) < mesh%xmin .OR. POI( 1) > mesh%xmax .OR. POI( 2) < mesh%ymin .OR. POI( 2) > mesh%ymax) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Initialise the refinement stack with the triangle containing the point
    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0
    ti_in = 1
    CALL find_containing_triangle( mesh, POI, ti_in)
    CALL add_triangle_to_refinement_stack_last( mesh, ti_in)

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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

      ! Check if it meets the resolution criterion
      IF (ti_in == ti) THEN
        ! This triangle contains the point
        longest_leg = longest_triangle_leg( va, vb, vc)
        meets_resolution_criterion = longest_leg <= res_max * C%mesh_resolution_tolerance
      ELSE
        ! This triangle does not contain the point
        meets_resolution_criterion = .TRUE.
      END IF

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_resolution_criterion) THEN
        ! Split triangle ti at its circumcenter

        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)

        ! Find out again which triangle contains the point, and add it to the stack
        CALL find_containing_triangle( mesh, POI, ti_in)
        CALL add_triangle_to_refinement_stack_last( mesh, ti_in)

      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    CALL refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_point

  SUBROUTINE refine_mesh_line( mesh, p_line, res_max, width, alpha_min)
    ! Refine a mesh based on a 1-D line criterion

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line        ! Collection of line segments
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: width         ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_line'
    INTEGER                                       :: nl
    INTEGER                                       :: ti,li
    REAL(dp), DIMENSION(2)                        :: pp,qq,pp2,qq2,pp_cropped,qq_cropped,dd
    LOGICAL                                       :: is_valid_line
    INTEGER                                       :: tip, tiq, via, vib, vic, tip2, tiq2
    REAL(dp), DIMENSION(2)                        :: va,vb,vc,llis
    LOGICAL                                       :: do_cross, crosses
    INTEGER                                       :: li_min, li_max
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the line is empty, the solution is trivial
    IF (SIZE( p_line,1) == 0 .AND. SIZE( p_line,2) == 0) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    nl = SIZE( p_line,1)
    IF (SIZE( p_line,2) /= 4) CALL crash('line must be an n-by-4 array!')

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

    DO li = 1, nl

      ! Line endpoints
      pp = [p_line( li,1), p_line( li,2)]
      qq = [p_line( li,3), p_line( li,4)]

      ! Crop line to mesh domain
      CALL crop_line_to_domain( pp, qq, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
      pp = pp_cropped
      qq = qq_cropped

      ! If the line segment lies outside of the mesh domain altogether, skip it
      IF (.NOT. is_valid_line) CYCLE

      ! If they coincide, this line is invalid - skip
      IF (NORM2( pp - qq) < mesh%tol_dist) CYCLE

      ! Find the triangles containing p and q
      tip = tiq
      CALL find_containing_triangle( mesh, pp, tip)
      tiq = tip
      CALL find_containing_triangle( mesh, qq, tiq)

      ! Add li to the overlap lists of tip and tiq
      mesh%Tri_li( tip,1) = MIN( mesh%Tri_li( tip,1),li)
      mesh%Tri_li( tip,2) = MAX( mesh%Tri_li( tip,2),li)
      mesh%Tri_li( tiq,1) = MIN( mesh%Tri_li( tiq,1),li)
      mesh%Tri_li( tiq,2) = MAX( mesh%Tri_li( tiq,2),li)

      ! If they both lie inside the same triangle, no need to trace
      IF (tip == tiq) CYCLE

      ! p and q lie in different triangles - perform a pseudo-trace
      DO WHILE (NORM2( pp - qq) > res_max)

        ! Move sideways to include line width
        pp2 = pp - [(pp( 2) - qq( 2)), qq( 1) - pp( 1)]
        qq2 = pp + [(pp( 2) - qq( 2)), qq( 1) - pp( 1)]

        CALL crop_line_to_domain( pp2, qq2, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
        pp2 = pp_cropped
        qq2 = qq_cropped

        IF (is_valid_line) THEN

          tip2 = tiq
          CALL find_containing_triangle( mesh, pp, tip2)
          tiq2 = tip2
          CALL find_containing_triangle( mesh, qq, tiq2)

          DO WHILE (NORM2 ( pp2 - qq2) > res_max)

            ! Move pp in the direction of qq
            dd = qq2 - pp2
            dd = dd / NORM2( dd)
            dd = dd * res_max
            pp2 = pp2 + dd

            ! Find the triangle that now contains pp
            CALL find_containing_triangle( mesh, pp2, tip2)

            ! Add li to the overlap lists of tip
            mesh%Tri_li( tip2,1) = MIN( mesh%Tri_li( tip2,1),li)
            mesh%Tri_li( tip2,2) = MAX( mesh%Tri_li( tip2,2),li)

            ! If we've reached tiq, the trace is done
            IF (tip2 == tiq2) EXIT

          END DO ! DO WHILE (NORM2 ( pp2 - qq2) > res_max)

        END IF ! IF (is_valid_line) THEN

        ! Move along pq

        ! Move pp in the direction of qq
        dd = qq - pp
        dd = dd / NORM2( dd)
        dd = dd * res_max
        pp = pp + dd

        ! Find the triangle that now contains pp
        CALL find_containing_triangle( mesh, pp, tip)

        ! Add li to the overlap lists of tip
        mesh%Tri_li( tip,1) = MIN( mesh%Tri_li( tip,1),li)
        mesh%Tri_li( tip,2) = MAX( mesh%Tri_li( tip,2),li)

        ! If we've reached tiq, the trace is done
        IF (tip == tiq) EXIT

      END DO ! DO WHILE (tip /= tiq)

    END DO ! DO li = 1, nl

  ! == Iteratively refine the mesh ==
  ! =================================

    ! If a triangle overlaps with any of the line segments, check if it is small enough.
    ! If not, split it. The new triangles will inherit the old one's line segment range.
    ! Update the (now reduced) overlap range for the new triangles.

    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0

    ! Mark which triangles need to be refined right now
    DO ti = 1, mesh%nTri
      IF (mesh%Tri_li( ti,2) == 0) THEN
        ! This triangle does not overlap with any line segments, so it doesn't need refining
        CYCLE
      ELSE
        ! Mark this triangle for refinement
        CALL add_triangle_to_refinement_stack_last( mesh, ti)
      END IF
    END DO

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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

      IF (.NOT. meets_geometry_criterion) THEN
        ! This triangle should be split anyway; no need to check
        ! the resolution criterion, so we can save some time)

        meets_resolution_criterion = .TRUE.

      ELSE
        ! This triangle meets the geometry criterion; check if
        ! it also meets the resolution criterion

        ! Initial guess for this triangle's line overlap range (taken from their parent)
        li_min = mesh%Tri_li( ti,1)
        li_max = mesh%Tri_li( ti,2)

        ! If that's already zero, skip
        IF (li_min == nl+1 .OR. li_max == 0) THEN
          ! No line overlap anyway

          meets_resolution_criterion = .TRUE.

        ELSE  ! IF (li_min == nl+1 .OR. li_max == 0) THEN
          ! Recalculate triangle overlap range

          mesh%Tri_li( ti,:) = [nl+1,0]

          DO li = li_min, li_max

            ! Line endpoints
            pp = [p_line( li,1), p_line( li,2)]
            qq = [p_line( li,3), p_line( li,4)]

            ! Crop line to mesh domain
            CALL crop_line_to_domain( pp, qq, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
            pp = pp_cropped
            qq = qq_cropped

            ! If the line segment lies outside of the mesh domain altogether, skip it
            IF (.NOT. is_valid_line) CYCLE

            ! If they coincide, this line is invalid - skip
            IF (NORM2( pp - qq) < mesh%tol_dist) CYCLE

            ! Check if this line segment crosses triangle ti
            crosses = .FALSE.
            crosses = crosses .OR. is_in_triangle( va, vb, vc, pp)
            crosses = crosses .OR. is_in_triangle( va, vb, vc, qq)
            CALL segment_intersection( pp, qq, va, vb, llis, do_cross, mesh%tol_dist)
            crosses = crosses .OR. do_cross
            CALL segment_intersection( pp, qq, vb, vc, llis, do_cross, mesh%tol_dist)
            crosses = crosses .OR. do_cross
            CALL segment_intersection( pp, qq, vc, va, llis, do_cross, mesh%tol_dist)
            crosses = crosses .OR. do_cross
            crosses = crosses .OR. lies_on_line_segment( pp, qq, va, width)
            crosses = crosses .OR. lies_on_line_segment( pp, qq, vb, width)
            crosses = crosses .OR. lies_on_line_segment( pp, qq, vc, width)

            ! If so, add it to the overlap range
            IF (crosses) THEN
              mesh%Tri_li( ti,1) = MIN( mesh%Tri_li( ti,1),li)
              mesh%Tri_li( ti,2) = MAX( mesh%Tri_li( ti,2),li)
            END IF

          END DO ! DO li = li_min, li_max

          ! If this triangle overlaps with any line segments,
          ! check if it meets the resolution criterion

          meets_resolution_criterion = .TRUE.
          IF (mesh%Tri_li( ti,1) <= mesh%Tri_li( ti,2)) THEN
            longest_leg = longest_triangle_leg( va, vb, vc)
            IF (longest_leg > res_max * C%mesh_resolution_tolerance) meets_resolution_criterion = .FALSE.
          END IF

        END IF ! IF (li_min == nl+1 .OR. li_max == 0) THEN

      END IF ! IF (.NOT. meets_geometry_criterion) THEN

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_resolution_criterion) THEN
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    CALL refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_line

  SUBROUTINE refine_mesh_polygon( mesh, poly, res_max, alpha_min)
    ! Refine a mesh based on a 2-D polygon criterion

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: poly          ! Polygon of interest
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_polygon'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line
    INTEGER                                       :: ti, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    LOGICAL                                       :: has_any_overlap
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the polygon is empty, the solution is trivial
    IF (SIZE( poly,1) == 0 .AND. SIZE( poly,2) == 0) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! First refine the mesh along the polygon perimeter by treating it as a 1-D line
    CALL poly2line( poly, p_line)
    CALL refine_mesh_line( mesh, p_line, res_max, res_max, alpha_min)

    ! Initialise the refinement stack with all triangles lying (partly) inside the polygon
    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0
    has_any_overlap = .FALSE.
    DO ti = 1, mesh%nTri

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      va  = mesh%V( via,:)
      vb  = mesh%V( vib,:)
      vc  = mesh%V( vic,:)

      ! If this triangle lies (partly) inside the polygon, mark it for refinement
      IF (is_in_polygon( poly, va) .OR. &
          is_in_polygon( poly, vb) .OR. &
          is_in_polygon( poly, vc)) THEN
        has_any_overlap = .TRUE.
        CALL add_triangle_to_refinement_stack_last( mesh, ti)
      END IF

    END DO ! DO ti = 1, mesh%nTri

    ! Safety: if the polygon doesn't overlap with the mesh
    !         at all, no refinement is needed (or indeed possible).
    IF (.NOT. has_any_overlap) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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
      IF (.NOT. meets_geometry_criterion) THEN
        meets_resolution_criterion = .TRUE.
      ELSE
        ! Check if the triangle meets the resolution criterion

        ! Check if the triangle lies inside the polygon
        IF (is_in_polygon( poly, va) .OR. &
            is_in_polygon( poly, vb) .OR. &
            is_in_polygon( poly, vc)) THEN
          ! The triangle lies inside the polygon

          longest_leg = longest_triangle_leg( va, vb, vc)
          meets_resolution_criterion = longest_leg <= res_max * C%mesh_resolution_tolerance

        ELSE
          ! The triangle does not lie inside the polygon
          meets_resolution_criterion = .TRUE.
        END IF

      END IF ! IF (.NOT. meets_geometry_criterion) THEN

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_resolution_criterion) THEN
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    CALL refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Clean up after yourself
    DEALLOCATE( p_line)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_polygon

  SUBROUTINE refine_mesh_split_encroaching_triangles( mesh, alpha_min)
    ! As a last step in any mesh refinement, make sure no triangles exist anymore
    ! that encroach on the domain border (i.e. have their circumcenter lying
    ! outside the domain)
    !
    ! Only check triangles on the domain border

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_split_encroaching_triangles'
    INTEGER                                       :: ti, n, vi, n_vertices_on_border
    INTEGER                                       :: via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    REAL(dp)                                      :: smallest_angle
    LOGICAL                                       :: meets_geometry_criterion
    LOGICAL                                       :: meets_encroachment_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the refinement stack with all border triangles
    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0
    DO ti = 1, mesh%nTri

      ! Count the number of vertices in this triangle on the domain border
      n_vertices_on_border = 0
      DO n = 1, 3
        vi = mesh%Tri( ti,n)
        IF (mesh%VBI( vi) > 0) THEN
          n_vertices_on_border = n_vertices_on_border + 1
        END IF
      END DO

      IF (n_vertices_on_border >= 2) THEN
        ! This triangle lies on the domain border
        CALL add_triangle_to_refinement_stack_last( mesh, ti)
      END IF

    END DO ! DO ti = 1, mesh%nTri

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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

      ! Check if it meets the encroachment criterion
      meets_encroachment_criterion = .TRUE.
      IF (mesh%Tricc( ti,1) < mesh%xmin .OR. mesh%Tricc( ti,1) > mesh%xmax .OR. &
          mesh%Tricc( ti,2) < mesh%ymin .OR. mesh%Tricc( ti,2) > mesh%ymax) THEN
        meets_encroachment_criterion = .FALSE.
      END IF

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_encroachment_criterion) THEN
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Safety - check if there are no more encroaching triangles left
      meets_encroachment_criterion = .TRUE.
    DO ti = 1, mesh%nTri
      IF (mesh%Tricc( ti,1) < mesh%xmin .OR. mesh%Tricc( ti,1) > mesh%xmax .OR. &
          mesh%Tricc( ti,2) < mesh%ymin .OR. mesh%Tricc( ti,2) > mesh%ymax) THEN
        meets_encroachment_criterion = .FALSE.
      END IF
    END DO
    IF (.NOT. meets_encroachment_criterion) CALL crash('did not remove all encroaching triangles!')

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_split_encroaching_triangles

  SUBROUTINE refine_mesh_split_encroaching_triangles_all( mesh, alpha_min)
    ! As a last step in any mesh refinement, make sure no triangles exist anymore
    ! that encroach on the domain border (i.e. have their circumcenter lying
    ! outside the domain)
    !
    ! Check all triangles

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_split_encroaching_triangles'
    INTEGER                                       :: ti
    INTEGER                                       :: via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    REAL(dp)                                      :: smallest_angle
    LOGICAL                                       :: meets_geometry_criterion
    LOGICAL                                       :: meets_encroachment_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the refinement stack with all border triangles
    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0
    DO ti = 1, mesh%nTri
      CALL add_triangle_to_refinement_stack_last( mesh, ti)
    END DO ! DO ti = 1, mesh%nTri

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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

      ! Check if it meets the encroachment criterion
      meets_encroachment_criterion = .TRUE.
      IF (mesh%Tricc( ti,1) < mesh%xmin .OR. mesh%Tricc( ti,1) > mesh%xmax .OR. &
          mesh%Tricc( ti,2) < mesh%ymin .OR. mesh%Tricc( ti,2) > mesh%ymax) THEN
        meets_encroachment_criterion = .FALSE.
      END IF

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_encroachment_criterion) THEN
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Safety - check if there are no more encroaching triangles left
      meets_encroachment_criterion = .TRUE.
    DO ti = 1, mesh%nTri
      IF (mesh%Tricc( ti,1) < mesh%xmin .OR. mesh%Tricc( ti,1) > mesh%xmax .OR. &
          mesh%Tricc( ti,2) < mesh%ymin .OR. mesh%Tricc( ti,2) > mesh%ymax) THEN
        meets_encroachment_criterion = .FALSE.
      END IF
    END DO
    IF (.NOT. meets_encroachment_criterion) CALL crash('did not remove all encroaching triangles!')

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_split_encroaching_triangles_all

! == Regions of interest

  SUBROUTINE refine_mesh_line_ROI( mesh, p_line_tot, res_max, width, alpha_min, poly_ROI)
    ! Refine a mesh based on a 1-D line criterion
    !
    ! Only consider line segments lying inside the region of interest
    ! described by the provided polygon

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: p_line_tot    ! Collection of line segments
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: width         ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_ROI      ! Polygon describing the region of interest

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_line_ROI'
    INTEGER                                       :: nl_tot
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line
    INTEGER                                       :: nl
    INTEGER                                       :: ti,li
    REAL(dp), DIMENSION(2)                        :: pp,qq,pp2,qq2,pp_cropped,qq_cropped,dd
    LOGICAL                                       :: is_valid_line
    INTEGER                                       :: tip, tiq, via, vib, vic, tip2, tiq2
    REAL(dp), DIMENSION(2)                        :: va,vb,vc,llis
    LOGICAL                                       :: do_cross, crosses
    INTEGER                                       :: li_min, li_max
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    nl_tot = SIZE( p_line_tot,1)

    ! If no line is provided, do nothing
    IF (nl_tot == 0) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (SIZE( p_line_tot,2) /= 4) CALL crash('line must be an n-by-4 array!')

  ! == Only consider line segments lying inside the region of interest
  !    described by the provided polygon
  ! ====================================

    ALLOCATE( p_line( nl_tot,4), source = 0._dp)
    nl = 0

    DO li = 1, nl_tot

      ! Line endpoints
      pp = [p_line_tot( li,1), p_line_tot( li,2)]
      qq = [p_line_tot( li,3), p_line_tot( li,4)]

      ! If this line lies inside the region of interest, consider it
      IF (is_in_polygon( poly_ROI, pp) .OR. is_in_polygon( poly_ROI, qq)) THEN
        nl = nl + 1
        p_line( nl,:) = p_line_tot( li,:)
      END IF

    END DO ! DO li = 1, nl_tot

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

    DO li = 1, nl

      ! Line endpoints
      pp = [p_line( li,1), p_line( li,2)]
      qq = [p_line( li,3), p_line( li,4)]

      ! Crop line to mesh domain
      CALL crop_line_to_domain( pp, qq, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
      pp = pp_cropped
      qq = qq_cropped

      ! If the line segment lies outside of the mesh domain altogether, skip it
      IF (.NOT. is_valid_line) CYCLE

      ! If they coincide, this line is invalid - skip
      IF (NORM2( pp - qq) < mesh%tol_dist) CYCLE

      ! Find the triangles containing p and q
      tip = tiq
      CALL find_containing_triangle( mesh, pp, tip)
      tiq = tip
      CALL find_containing_triangle( mesh, qq, tiq)

      ! Add li to the overlap lists of tip and tiq
      mesh%Tri_li( tip,1) = MIN( mesh%Tri_li( tip,1),li)
      mesh%Tri_li( tip,2) = MAX( mesh%Tri_li( tip,2),li)
      mesh%Tri_li( tiq,1) = MIN( mesh%Tri_li( tiq,1),li)
      mesh%Tri_li( tiq,2) = MAX( mesh%Tri_li( tiq,2),li)

      ! If they both lie inside the same triangle, no need to trace
      IF (tip == tiq) CYCLE

      ! p and q lie in different triangles - perform a pseudo-trace
      DO WHILE (NORM2( pp - qq) > res_max)

        ! Move sideways to include line width
        pp2 = pp - [(pp( 2) - qq( 2)), qq( 1) - pp( 1)]
        qq2 = pp + [(pp( 2) - qq( 2)), qq( 1) - pp( 1)]

        CALL crop_line_to_domain( pp2, qq2, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
        pp2 = pp_cropped
        qq2 = qq_cropped

        IF (is_valid_line) THEN

          tip2 = tiq
          CALL find_containing_triangle( mesh, pp, tip2)
          tiq2 = tip2
          CALL find_containing_triangle( mesh, qq, tiq2)

          DO WHILE (NORM2 ( pp2 - qq2) > res_max)

            ! Move pp in the direction of qq
            dd = qq2 - pp2
            dd = dd / NORM2( dd)
            dd = dd * res_max
            pp2 = pp2 + dd

            ! Find the triangle that now contains pp
            CALL find_containing_triangle( mesh, pp2, tip2)

            ! Add li to the overlap lists of tip
            mesh%Tri_li( tip2,1) = MIN( mesh%Tri_li( tip2,1),li)
            mesh%Tri_li( tip2,2) = MAX( mesh%Tri_li( tip2,2),li)

            ! If we've reached tiq, the trace is done
            IF (tip2 == tiq2) EXIT

          END DO ! DO WHILE (NORM2 ( pp2 - qq2) > res_max)

        END IF ! IF (is_valid_line) THEN

        ! Move along pq

        ! Move pp in the direction of qq
        dd = qq - pp
        dd = dd / NORM2( dd)
        dd = dd * res_max
        pp = pp + dd

        ! Find the triangle that now contains pp
        CALL find_containing_triangle( mesh, pp, tip)

        ! Add li to the overlap lists of tip
        mesh%Tri_li( tip,1) = MIN( mesh%Tri_li( tip,1),li)
        mesh%Tri_li( tip,2) = MAX( mesh%Tri_li( tip,2),li)

        ! If we've reached tiq, the trace is done
        IF (tip == tiq) EXIT

      END DO ! DO WHILE (tip /= tiq)

    END DO ! DO li = 1, nl

  ! == Iteratively refine the mesh ==
  ! =================================

    ! If a triangle overlaps with any of the line segments, check if it is small enough.
    ! If not, split it. The new triangles will inherit the old one's line segment range.
    ! Update the (now reduced) overlap range for the new triangles.

    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0

    ! Mark which triangles need to be refined right now
    DO ti = 1, mesh%nTri
      IF (mesh%Tri_li( ti,2) == 0) THEN
        ! This triangle does not overlap with any line segments, so it doesn't need refining
        CYCLE
      ELSE
        ! Mark this triangle for refinement
        CALL add_triangle_to_refinement_stack_last( mesh, ti)
      END IF
    END DO

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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

      IF (.NOT. meets_geometry_criterion) THEN
        ! This triangle should be split anyway; no need to check
        ! the resolution criterion, so we can save some time)

        meets_resolution_criterion = .TRUE.

      ELSE
        ! This triangle meets the geometry criterion; check if
        ! it also meets the resolution criterion

        ! Initial guess for this triangle's line overlap range (taken from their parent)
        li_min = mesh%Tri_li( ti,1)
        li_max = mesh%Tri_li( ti,2)

        ! If that's already zero, skip
        IF (li_min == nl+1 .OR. li_max == 0) THEN
          ! No line overlap anyway

          meets_resolution_criterion = .TRUE.

        ELSE  ! IF (li_min == nl+1 .OR. li_max == 0) THEN
          ! Recalculate triangle overlap range

          mesh%Tri_li( ti,:) = [nl+1,0]

          DO li = li_min, li_max

            ! Line endpoints
            pp = [p_line( li,1), p_line( li,2)]
            qq = [p_line( li,3), p_line( li,4)]

            ! Crop line to mesh domain
            CALL crop_line_to_domain( pp, qq, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp_cropped, qq_cropped, is_valid_line)
            pp = pp_cropped
            qq = qq_cropped

            ! If the line segment lies outside of the mesh domain altogether, skip it
            IF (.NOT. is_valid_line) CYCLE

            ! If they coincide, this line is invalid - skip
            IF (NORM2( pp - qq) < mesh%tol_dist) CYCLE

            ! Check if this line segment crosses triangle ti
            crosses = .FALSE.
            crosses = crosses .OR. is_in_triangle( va, vb, vc, pp)
            crosses = crosses .OR. is_in_triangle( va, vb, vc, qq)
            CALL segment_intersection( pp, qq, va, vb, llis, do_cross, mesh%tol_dist)
            crosses = crosses .OR. do_cross
            CALL segment_intersection( pp, qq, vb, vc, llis, do_cross, mesh%tol_dist)
            crosses = crosses .OR. do_cross
            CALL segment_intersection( pp, qq, vc, va, llis, do_cross, mesh%tol_dist)
            crosses = crosses .OR. do_cross
            crosses = crosses .OR. lies_on_line_segment( pp, qq, va, width)
            crosses = crosses .OR. lies_on_line_segment( pp, qq, vb, width)
            crosses = crosses .OR. lies_on_line_segment( pp, qq, vc, width)

            ! If so, add it to the overlap range
            IF (crosses) THEN
              mesh%Tri_li( ti,1) = MIN( mesh%Tri_li( ti,1),li)
              mesh%Tri_li( ti,2) = MAX( mesh%Tri_li( ti,2),li)
            END IF

          END DO ! DO li = li_min, li_max

          ! If this triangle overlaps with any line segments,
          ! check if it meets the resolution criterion

          meets_resolution_criterion = .TRUE.
          IF (mesh%Tri_li( ti,1) <= mesh%Tri_li( ti,2)) THEN
            longest_leg = longest_triangle_leg( va, vb, vc)
            IF (longest_leg > res_max * C%mesh_resolution_tolerance) meets_resolution_criterion = .FALSE.
          END IF

        END IF ! IF (li_min == nl+1 .OR. li_max == 0) THEN

      END IF ! IF (.NOT. meets_geometry_criterion) THEN

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_resolution_criterion) THEN
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    CALL refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_line_ROI

  SUBROUTINE refine_mesh_polygon_ROI( mesh, poly, res_max, alpha_min, poly_ROI)
    ! Refine a mesh based on a 2-D polygon criterion
    !
    ! Only consider triangles lying inside the region of interest
    ! described by the provided polygon

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: poly          ! Polygon of interest
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: poly_ROI      ! Polygon describing the region of interest

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_polygon_ROI'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: p_line
    INTEGER                                       :: ti, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    LOGICAL                                       :: has_any_overlap
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! First refine the mesh along the polygon perimeter by treating it as a 1-D line
    CALL poly2line( poly, p_line)
    CALL refine_mesh_line_ROI( mesh, p_line, res_max, res_max, alpha_min, poly_ROI)

    ! Initialise the refinement stack with all triangles lying (partly) inside the polygon
    mesh%refinement_stackN = 0
    mesh%refinement_map    = 0
    has_any_overlap = .FALSE.
    DO ti = 1, mesh%nTri

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      va  = mesh%V( via,:)
      vb  = mesh%V( vib,:)
      vc  = mesh%V( vic,:)

      ! If this triangle lies (partly) inside the polygon, mark it for refinement
      IF ((is_in_polygon( poly    , va) .OR. is_in_polygon( poly    , vb) .OR. is_in_polygon( poly    , vc)) .AND. &
          (is_in_polygon( poly_ROI, va) .OR. is_in_polygon( poly_ROI, vb) .OR. is_in_polygon( poly_ROI, vc))) THEN
        has_any_overlap = .TRUE.
        CALL add_triangle_to_refinement_stack_last( mesh, ti)
      END IF

    END DO ! DO ti = 1, mesh%nTri

    ! Safety: if the polygon doesn't overlap with the mesh
    !         at all, no refinement is needed (or indeed possible).
    IF (.NOT. has_any_overlap) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Keep refining until all triangles match the criterion
    DO WHILE (mesh%refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF

      ! Take the first triangle in the stack
      ti = mesh%refinement_stack( 1)
      CALL remove_triangle_from_refinement_stack( mesh, ti)

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
      IF (.NOT. meets_geometry_criterion) THEN
        meets_resolution_criterion = .TRUE.
      ELSE
        ! Check if the triangle meets the resolution criterion

        ! Check if the triangle lies inside the polygon
        IF ((is_in_polygon( poly    , va) .OR. is_in_polygon( poly    , vb) .OR. is_in_polygon( poly    , vc)) .AND. &
            (is_in_polygon( poly_ROI, va) .OR. is_in_polygon( poly_ROI, vb) .OR. is_in_polygon( poly_ROI, vc))) THEN
          ! The triangle lies inside the polygon

          longest_leg = longest_triangle_leg( va, vb, vc)
          meets_resolution_criterion = longest_leg <= res_max * C%mesh_resolution_tolerance

        ELSE
          ! The triangle does not lie inside the polygon
          meets_resolution_criterion = .TRUE.
        END IF

      END IF ! IF (.NOT. meets_geometry_criterion) THEN

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_resolution_criterion) THEN
        ! Split triangle ti at its circumcenter
        p_new = circumcenter( va, vb, vc)
        CALL split_triangle( mesh, ti, p_new)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Final step to ensure a nice clean mesh
    CALL refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Clean up after yourself
    DEALLOCATE( p_line)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_polygon_ROI

! == Lloyd's algorithm for "smoothing" a mesh

  SUBROUTINE Lloyds_algorithm_single_iteration( mesh, alpha_min)
    ! Lloyd's algorithm: move all vertices to the geometric centers of their Voronoi cells, and update the triangulation.
    ! This "smooths" the mesh, reducing resolution gradients and widening internal angles, thus making it more
    ! suitable for numerical methods (particularly the SSA).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'Lloyds_algorithm_single_iteration'
    INTEGER                                       :: vi, ci, cip1
    REAL(dp)                                      :: VorTriA, sumVorTriA
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc, VorGC

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Move all non-border vertices to their Voronoi cell geometric centre
    DO vi = 1, mesh%nV

      ! Leave border vertices where they are
      IF (mesh%VBI( vi) > 0) CYCLE

      ! Find the geometric centre of this vertex' Voronoi cell
      VorGC      = 0._dp
      sumVorTriA = 0._dp

      DO ci = 1, mesh%nC( vi)

        cip1 = ci + 1
        IF (cip1 > mesh%nC( vi)) cip1 = 1

        pa = mesh%V( vi,:)
        pb = mesh%V( mesh%C( vi,ci  ),:)
        pc = mesh%V( mesh%C( vi,cip1),:)

        VorTriA = cross2( pb - pa, pc - pa)

        VorGC = VorGC + VorTriA * (pa + pb + pc) / 3._dp
        sumVorTriA   = sumVorTriA   + VorTriA

      END DO ! DO ci = 1, mesh%nC( vi)

      VorGC = VorGC / sumVorTriA

      ! Move the vertex
      CALL move_vertex( mesh, vi, VorGC)

    END DO

    ! Final step to ensure a nice clean mesh
    CALL refine_mesh_split_encroaching_triangles_all( mesh, alpha_min)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE Lloyds_algorithm_single_iteration

! == Enforce contiguous process domains

  SUBROUTINE enforce_contiguous_process_domains( mesh)
    ! Shuffle vertices, triangles, and edges so that the vertices owned by each process form
    ! a contiguous domain with short borders, to minimise the data volumes for halo exchanges.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'enforce_contiguous_process_domains'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Shuffle vertices
    CALL enforce_contiguous_process_domains_vertices( mesh)

    ! Shuffle triangles
    CALL enforce_contiguous_process_domains_triangles( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE enforce_contiguous_process_domains

  SUBROUTINE enforce_contiguous_process_domains_vertices( mesh)
    ! Shuffle vertices so that the vertices owned by each process form
    ! a contiguous domain with short borders, to minimise the data volumes for halo exchanges.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'enforce_contiguous_process_domains_vertices'
    REAL(dp), DIMENSION(mesh%nV)                  :: xx
    INTEGER,  DIMENSION(mesh%nV)                  :: vi_new2vi_old, vi_old2vi_new
    INTEGER                                       :: vi_old, vi_new
    REAL(dp), DIMENSION(mesh%nV,2)                :: V_old
    INTEGER,  DIMENSION(mesh%nV)                  :: nC_old
    INTEGER,  DIMENSION(mesh%nV,mesh%nC_mem)      :: C_old
    INTEGER,  DIMENSION(mesh%nV)                  :: niTri_old
    INTEGER,  DIMENSION(mesh%nV,mesh%nC_mem)      :: iTri_old
    INTEGER,  DIMENSION(mesh%nV)                  :: VBI_old
    INTEGER                                       :: ci,ti,n,vj_old,vj_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Sort vertices by x-coordinate
    xx = mesh%V( :,1)
    CALL quick_n_dirty_sort( xx, vi_new2vi_old)

    ! Calculate translation table in the opposite direction
    DO vi_new = 1, mesh%nV
      vi_old = vi_new2vi_old( vi_new)
      vi_old2vi_new( vi_old) = vi_new
    END DO

    ! Shuffle vertex data: V, nC, C, niTri, iTri, VBI
    ! ===============================================

    V_old     = mesh%V
    nC_old    = mesh%nC
    C_old     = mesh%C
    niTri_old = mesh%niTri
    iTri_old  = mesh%iTri
    VBI_old   = mesh%VBI

    mesh%V     = 0._dp
    mesh%nC    = 0
    mesh%C     = 0
    mesh%niTri = 0
    mesh%iTri  = 0
    mesh%VBI   = 0

    DO vi_new = 1, mesh%nV

      ! This new vertex corresponds to this old vertex
      vi_old = vi_new2vi_old( vi_new)

      ! V
      mesh%V( vi_new,:) = V_old( vi_old,:)

      ! nC
      mesh%nC( vi_new) = nC_old( vi_old)

      ! C
      DO ci = 1, mesh%nC( vi_new)
        vj_old = C_old( vi_old,ci)
        vj_new = vi_old2vi_new( vj_old)
        mesh%C( vi_new,ci) = vj_new
      END DO

      ! niTri
      mesh%niTri( vi_new) = niTri_old( vi_old)

      ! iTri
      mesh%iTri( vi_new,:) = iTri_old( vi_old,:)

      ! VBI
      mesh%VBI( vi_new) = VBI_old( vi_old)

    END DO

    ! Shuffle triangle data: Tri
    ! ==========================

    DO ti = 1, mesh%nTri
      DO n = 1, 3
        vi_old = mesh%Tri( ti,n)
        vi_new = vi_old2vi_new( vi_old)
        mesh%Tri( ti,n) = vi_new
      END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE enforce_contiguous_process_domains_vertices

  SUBROUTINE enforce_contiguous_process_domains_triangles( mesh)
    ! Shuffle triangles so that the triangles owned by each process form
    ! a contiguous domain with short borders, to minimise the data volumes for halo exchanges.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'enforce_contiguous_process_domains_triangles'
    REAL(dp), DIMENSION(mesh%nTri)                :: xx
    INTEGER                                       :: ti,via,vib,vic
    REAL(dp), DIMENSION(2)                        :: va,vb,vc,gc
    INTEGER,  DIMENSION(mesh%nTri)                :: ti_new2ti_old, ti_old2ti_new
    INTEGER                                       :: ti_old, ti_new
    INTEGER,  DIMENSION(mesh%nV,mesh%nC_mem)      :: iTri_old
    INTEGER,  DIMENSION(mesh%nTri,3)              :: Tri_old
    INTEGER,  DIMENSION(mesh%nTri,3)              :: TriC_old
    REAL(dp), DIMENSION(mesh%nTri,2)              :: Tricc_old
    INTEGER                                       :: vi,iti,n,tj_old,tj_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate triangle geometric centres
    DO ti = 1, mesh%nTri

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      va = mesh%V( via,:)
      vb = mesh%V( vib,:)
      vc = mesh%V( vic,:)

      gc = geometric_center( va, vb, vc)

      xx( ti) = gc( 1)

    END DO

    ! Sort triangles by x-coordinate
    CALL quick_n_dirty_sort( xx, ti_new2ti_old)

    ! Calculate translation table in the opposite direction
    DO ti_new = 1, mesh%nTri
      ti_old = ti_new2ti_old( ti_new)
      ti_old2ti_new( ti_old) = ti_new
    END DO

    ! Shuffle vertex data: iTri
    ! =========================

    iTri_old  = mesh%iTri

    mesh%iTri  = 0

    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri( vi)
        ti_old = iTri_old( vi,iti)
        ti_new = ti_old2ti_new( ti_old)
        mesh%iTri( vi,iti) = ti_new
      END DO
    END DO

    ! Shuffle triangle data: Tri, TriC, Tricc
    ! =======================================

    Tri_old   = mesh%Tri
    TriC_old  = mesh%TriC
    Tricc_old = mesh%Tricc

    mesh%Tri   = 0
    mesh%TriC  = 0
    mesh%Tricc = 0._dp

    DO ti_new = 1, mesh%nTri

      ! This new triangle corresponds to this old triangle
      ti_old = ti_new2ti_old( ti_new)

      ! Tri
      mesh%Tri( ti_new,:) = Tri_old( ti_old,:)

      ! TriC
      DO n = 1, 3
        tj_old = TriC_old( ti_old,n)
        IF (tj_old == 0) THEN
          tj_new = 0
        ELSE
          tj_new = ti_old2ti_new( tj_old)
        END IF
        mesh%TriC( ti_new,n) = tj_new
      END DO

      ! Tricc
      mesh%Tricc( ti_new,:) = Tricc_old( ti_old,:)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE enforce_contiguous_process_domains_triangles

! == Polygons describing specific regions

  SUBROUTINE calc_polygon_Pine_Island_Glacier( poly)
    ! Return a polygon enveloping the Pine Island Glacier catchment basin
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Pine_Island_Glacier'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 42,2))

    poly(  1,:) = [ -1.64e6_dp, -3.4e5_dp]
    poly(  2,:) = [-1.60e6_dp, -3.5e5_dp]
    poly(  3,:) = [-1.55e6_dp, -3.4e5_dp]
    poly(  4,:) = [-1.50e6_dp, -3.2e5_dp]
    poly(  5,:) = [-1.45e6_dp, -2.9e5_dp]
    poly(  6,:) = [-1.40e6_dp, -2.5e5_dp]
    poly(  7,:) = [-1.37e6_dp, -2.0e5_dp]
    poly(  8,:) = [-1.34e6_dp, -1.7e5_dp]
    poly(  9,:) = [-1.30e6_dp, -1.6e5_dp]
    poly( 10,:) = [-1.26e6_dp, -1.6e5_dp]
    poly( 11,:) = [-1.22e6_dp, -1.7e5_dp]
    poly( 12,:) = [-1.18e6_dp, -1.75e5_dp]
    poly( 13,:) = [-1.14e6_dp, -1.75e5_dp]
    poly( 14,:) = [-1.11e6_dp, -1.72e5_dp]
    poly( 15,:) = [-1.09e6_dp, -1.6e5_dp]
    poly( 16,:) = [-1.085e6_dp, -1.4e5_dp]
    poly( 17,:) = [-1.09e6_dp, -1.2e5_dp]
    poly( 18,:) = [-1.1e6_dp, -1.0e5_dp]
    poly( 19,:) = [-1.13e6_dp, -0.7e5_dp]
    poly( 20,:) = [-1.17e6_dp, -0.4e5_dp]
    poly( 21,:) = [-1.21e6_dp, -0.2e5_dp]
    poly( 22,:) = [-1.26e6_dp, -0.0e5_dp]
    poly( 23,:) = [-1.32e6_dp, 0.1e5_dp]
    poly( 24,:) = [-1.45e6_dp, 0.1e5_dp]
    poly( 25,:) = [-1.48e6_dp, 0.15e5_dp]
    poly( 26,:) = [-1.51e6_dp, 0.35e5_dp]
    poly( 27,:) = [-1.53e6_dp, 0.75e5_dp]
    poly( 28,:) = [-1.55e6_dp, 0.95e5_dp]
    poly( 29,:) = [-1.58e6_dp, 0.1e6_dp]
    poly( 30,:) = [-1.62e6_dp, 0.11e6_dp]
    poly( 31,:) = [-1.65e6_dp, 0.12e6_dp]
    poly( 32,:) = [-1.67e6_dp, 0.10e6_dp]
    poly( 33,:) = [-1.69e6_dp, 0.9e5_dp]
    poly( 34,:) = [-1.71e6_dp, 0.5e5_dp]
    poly( 35,:) = [-1.74e6_dp, 0.1e5_dp]
    poly( 36,:) = [-1.75e6_dp, -0.5e5_dp]
    poly( 37,:) = [-1.75e6_dp, -0.15e6_dp]
    poly( 38,:) = [-1.71e6_dp, -0.19e6_dp]
    poly( 39,:) = [-1.66e6_dp, -0.2e6_dp]
    poly( 40,:) = [-1.64e6_dp, -0.21e6_dp]
    poly( 41,:) = [-1.63e6_dp, -0.23e6_dp]
    poly( 42,:) = [-1.63e6_dp, -0.29e6_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Pine_Island_Glacier

  SUBROUTINE calc_polygon_Thwaites_Glacier( poly)
    ! Return a polygon enveloping the Pine Island Glacier catchment basin
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Thwaites_Glacier'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 44,2))

    poly(  1,:) = [-1.6e6_dp, -5.4e5_dp]
    poly(  2,:) = [-1.55e6_dp, -5.4e5_dp]
    poly(  3,:) = [-1.50e6_dp, -5.5e5_dp]
    poly(  4,:) = [-1.45e6_dp, -5.6e5_dp]
    poly(  5,:) = [-1.40e6_dp, -5.65e5_dp]
    poly(  6,:) = [-1.37e6_dp, -5.75e5_dp]
    poly(  7,:) = [-1.35e6_dp, -6e5_dp]
    poly(  8,:) = [-1.35e6_dp, -6.5e5_dp]
    poly(  9,:) = [-1.34e6_dp, -6.9e5_dp]
    poly( 10,:) = [-1.32e6_dp, -7.3e5_dp]
    poly( 11,:) = [-1.29e6_dp, -7.6e5_dp]
    poly( 12,:) = [-1.25e6_dp, -7.8e5_dp]
    poly( 13,:) = [-1.22e6_dp, -7.8e5_dp]
    poly( 14,:) = [-1.20e6_dp, -7.6e5_dp]
    poly( 15,:) = [-1.18e6_dp, -7.4e5_dp]
    poly( 16,:) = [-1.15e6_dp, -6.9e5_dp]
    poly( 17,:) = [-1.14e6_dp, -6.4e5_dp]
    poly( 18,:) = [-1.14e6_dp, -5.9e5_dp]
    poly( 19,:) = [-1.11e6_dp, -5.6e5_dp]
    poly( 20,:) = [-1.08e6_dp, -5.5e5_dp]
    poly( 21,:) = [-1.04e6_dp, -5.4e5_dp]
    poly( 22,:) = [-1.01e6_dp, -5.2e5_dp]
    poly( 23,:) = [-0.99e6_dp, -5.0e5_dp]
    poly( 24,:) = [-0.99e6_dp, -4.6e5_dp]
    poly( 25,:) = [-1.02e6_dp, -4.4e5_dp]
    poly( 26,:) = [-1.04e6_dp, -4.2e5_dp]
    poly( 27,:) = [-1.06e6_dp, -3.9e5_dp]
    poly( 28,:) = [-1.07e6_dp, -3.5e5_dp]
    poly( 29,:) = [-1.07e6_dp, -3.2e5_dp]
    poly( 30,:) = [-1.09e6_dp, -2.8e5_dp]
    poly( 31,:) = [-1.12e6_dp, -2.5e5_dp]
    poly( 32,:) = [-1.15e6_dp, -2.2e5_dp]
    poly( 33,:) = [-1.18e6_dp, -1.9e5_dp]
    poly( 34,:) = [-1.22e6_dp, -1.7e5_dp]
    poly( 35,:) = [-1.26e6_dp, -1.6e5_dp]
    poly( 36,:) = [-1.30e6_dp, -1.6e5_dp]
    poly( 37,:) = [-1.34e6_dp, -1.7e5_dp]
    poly( 38,:) = [-1.37e6_dp, -2.0e5_dp]
    poly( 39,:) = [-1.40e6_dp, -2.5e5_dp]
    poly( 40,:) = [-1.45e6_dp, -2.9e5_dp]
    poly( 41,:) = [-1.50e6_dp, -3.2e5_dp]
    poly( 42,:) = [-1.55e6_dp, -3.4e5_dp]
    poly( 43,:) = [-1.60e6_dp, -3.5e5_dp]
    poly( 44,:) = [-1.64e6_dp, -3.4e5_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Thwaites_Glacier

  SUBROUTINE calc_polygon_Amery_ice_shelf( poly)
    ! Return a polygon enveloping the Amery ice shelf catchment basin
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Amery_ice_shelf'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 11,2))

    poly(  1,:) = [2.2798e6_dp, 0.8624e6_dp]
    poly(  2,:) = [1.9637e6_dp, 0.9955e6_dp]
    poly(  3,:) = [1.4229e6_dp, 0.9234e6_dp]
    poly(  4,:) = [1.3480e6_dp, 0.7792e6_dp]
    poly(  5,:) = [1.2981e6_dp, 0.6711e6_dp]
    poly(  6,:) = [1.4340e6_dp, 0.4353e6_dp]
    poly(  7,:) = [1.6337e6_dp, 0.4742e6_dp]
    poly(  8,:) = [1.8056e6_dp, 0.5019e6_dp]
    poly(  9,:) = [1.8777e6_dp, 0.4215e6_dp]
    poly( 10,:) = [2.1079e6_dp, 0.4520e6_dp]
    poly( 11,:) = [2.3075e6_dp, 0.6711e6_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Amery_ice_shelf

  SUBROUTINE calc_polygon_Riiser_Larsen_ice_shelf( poly)
    ! Return a polygon enveloping the Riiser-Larsen ice shelf catchment basin
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Riiser_Larsen_ice_shelf'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 31,2))

    poly(  1,:) = [-0.6469e6_dp, 1.6448e6_dp]
    poly(  2,:) = [-0.6507e6_dp, 1.7370e6_dp]
    poly(  3,:) = [-0.6411e6_dp, 1.8005e6_dp]
    poly(  4,:) = [-0.5989e6_dp, 1.8370e6_dp]
    poly(  5,:) = [-0.5508e6_dp, 1.8639e6_dp]
    poly(  6,:) = [-0.5104e6_dp, 1.9081e6_dp]
    poly(  7,:) = [-0.4758e6_dp, 1.9331e6_dp]
    poly(  8,:) = [-0.4451e6_dp, 1.9542e6_dp]
    poly(  9,:) = [-0.4393e6_dp, 1.9946e6_dp]
    poly( 10,:) = [-0.3336e6_dp, 1.9696e6_dp]
    poly( 11,:) = [-0.3048e6_dp, 1.9292e6_dp]
    poly( 12,:) = [-0.2644e6_dp, 1.9081e6_dp]
    poly( 13,:) = [-0.2029e6_dp, 1.8927e6_dp]
    poly( 14,:) = [-0.1741e6_dp, 1.8716e6_dp]
    poly( 15,:) = [-0.1644e6_dp, 1.8351e6_dp]
    poly( 16,:) = [-0.1414e6_dp, 1.8043e6_dp]
    poly( 17,:) = [-0.1222e6_dp, 1.7659e6_dp]
    poly( 18,:) = [-0.1202e6_dp, 1.7313e6_dp]
    poly( 19,:) = [-0.1318e6_dp, 1.6928e6_dp]
    poly( 20,:) = [-0.1644e6_dp, 1.6640e6_dp]
    poly( 21,:) = [-0.2125e6_dp, 1.6275e6_dp]
    poly( 22,:) = [-0.2394e6_dp, 1.5948e6_dp]
    poly( 23,:) = [-0.2663e6_dp, 1.5833e6_dp]
    poly( 24,:) = [-0.3259e6_dp, 1.5813e6_dp]
    poly( 25,:) = [-0.3778e6_dp, 1.5717e6_dp]
    poly( 26,:) = [-0.4201e6_dp, 1.5640e6_dp]
    poly( 27,:) = [-0.4528e6_dp, 1.5640e6_dp]
    poly( 28,:) = [-0.4931e6_dp, 1.5660e6_dp]
    poly( 29,:) = [-0.5354e6_dp, 1.5698e6_dp]
    poly( 30,:) = [-0.5758e6_dp, 1.5871e6_dp]
    poly( 31,:) = [-0.6142e6_dp, 1.6102e6_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Riiser_Larsen_ice_shelf

  SUBROUTINE calc_polygon_Siple_Coast( poly)
    ! Return a polygon enveloping the Siple Coast area
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Siple_Coast'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 21,2))

    poly(  1,:) = [-0.6394e6_dp, -0.2184e6_dp]
    poly(  2,:) = [-0.6852e6_dp, -0.3498e6_dp]
    poly(  3,:) = [-0.7219e6_dp, -0.3101e6_dp]
    poly(  4,:) = [-0.8165e6_dp, -0.2979e6_dp]
    poly(  5,:) = [-0.8288e6_dp, -0.3681e6_dp]
    poly(  6,:) = [-0.7402e6_dp, -0.4567e6_dp]
    poly(  7,:) = [-1.0059e6_dp, -0.3803e6_dp]
    poly(  8,:) = [-1.0029e6_dp, -0.4689e6_dp]
    poly(  9,:) = [-0.9326e6_dp, -0.5514e6_dp]
    poly( 10,:) = [-0.8440e6_dp, -0.6125e6_dp]
    poly( 11,:) = [-1.0609e6_dp, -0.6033e6_dp]
    poly( 12,:) = [-0.8807e6_dp, -0.6980e6_dp]
    poly( 13,:) = [-1.0273e6_dp, -0.7652e6_dp]
    poly( 14,:) = [-1.0609e6_dp, -0.9210e6_dp]
    poly( 15,:) = [-0.9876e6_dp, -1.0737e6_dp]
    poly( 16,:) = [-0.7463e6_dp, -1.0004e6_dp]
    poly( 17,:) = [-0.6363e6_dp, -1.0981e6_dp]
    poly( 18,:) = [-0.5019e6_dp, -1.1287e6_dp]
    poly( 19,:) = [ 0.0051e6_dp, -0.8355e6_dp]
    poly( 20,:) = [-0.0132e6_dp, -0.2887e6_dp]
    poly( 21,:) = [-0.3034e6_dp, -0.1573e6_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Siple_Coast

  SUBROUTINE calc_polygon_Larsen_ice_shelf( poly)
    ! Return a polygon enveloping the Larsen C ice shelf area
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Larsen_ice_shelf'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 5,2))

    poly(  1,:) = [-2.3962e6_dp, 0.8370e6_dp]
    poly(  2,:) = [-1.9819e6_dp, 0.8482e6_dp]
    poly(  3,:) = [-1.8363e6_dp, 1.0721e6_dp]
    poly(  4,:) = [-2.4857e6_dp, 1.6880e6_dp]
    poly(  5,:) = [-2.6985e6_dp, 1.3968e6_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Larsen_ice_shelf

  SUBROUTINE calc_polygon_Transantarctic_Mountains( poly)
    ! Return a polygon enveloping the Transantarctic Mountains
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Transantarctic_Mountains'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 12,2))

    poly(  1,:) = [ 0.2911e6_dp, -1.3464e6_dp]
    poly(  2,:) = [ 0.5487e6_dp, -1.2233e6_dp]
    poly(  3,:) = [ 0.6158e6_dp, -1.1225e6_dp]
    poly(  4,:) = [ 0.5934e6_dp, -0.7978e6_dp]
    poly(  5,:) = [ 0.3695e6_dp, -0.5067e6_dp]
    poly(  6,:) = [ 0.1680e6_dp, -0.3387e6_dp]
    poly(  7,:) = [-0.1792e6_dp, -0.1708e6_dp]
    poly(  8,:) = [-0.4143e6_dp, -0.1484e6_dp]
    poly(  9,:) = [-0.4591e6_dp, -0.3947e6_dp]
    poly( 10,:) = [ 0.0672e6_dp, -0.7306e6_dp]
    poly( 11,:) = [ 0.2127e6_dp, -0.8762e6_dp]
    poly( 12,:) = [ 0.3359e6_dp, -1.0217e6_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Transantarctic_Mountains

  SUBROUTINE calc_polygon_Patagonia( poly)
    ! Return a polygon enveloping the region where the former
    ! Patagonian ice sheet peaked during the last glacial maximum
    !
    ! (based on manual analysis of the PATICE reconstruction,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do mesh refinement)
    !
    ! This assumes a very particular stereographic projection for
    ! the model domain, which as of now reads:
    !
    ! lambda_M_ANT_config    = 289.0
    ! phi_M_ANT_config       = -47.0
    ! beta_stereo_ANT_config = 71.0
    ! xmin_ANT_config        = -400000.0
    ! xmax_ANT_config        =  400000.0
    ! ymin_ANT_config        = -1110000.0
    ! ymax_ANT_config        =  1110000.0

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Patagonia'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 130,2))

    poly(  1,:) = [-0.5409e5_dp,  0.8469e5_dp]
    poly(  2,:) = [-0.1447e5_dp,  0.7109e5_dp]
    poly(  3,:) = [-0.4643e5_dp,  0.4652e5_dp]
    poly(  4,:) = [-0.8643e5_dp,  0.3775e5_dp]
    poly(  5,:) = [-0.6965e5_dp, -0.0131e5_dp]
    poly(  6,:) = [-0.3570e5_dp, -0.2868e5_dp]
    poly(  7,:) = [-0.5385e5_dp, -0.6522e5_dp]
    poly(  8,:) = [-0.7107e5_dp, -1.0253e5_dp]
    poly(  9,:) = [-0.8151e5_dp, -1.4277e5_dp]
    poly( 10,:) = [-1.0274e5_dp, -1.8160e5_dp]
    poly( 11,:) = [-0.8547e5_dp, -2.1939e5_dp]
    poly( 12,:) = [-1.2204e5_dp, -2.3569e5_dp]
    poly( 13,:) = [-1.0413e5_dp, -2.7274e5_dp]
    poly( 14,:) = [-0.6790e5_dp, -2.9327e5_dp]
    poly( 15,:) = [-1.0808e5_dp, -2.9403e5_dp]
    poly( 16,:) = [-1.2326e5_dp, -3.3331e5_dp]
    poly( 17,:) = [-0.8064e5_dp, -3.3984e5_dp]
    poly( 18,:) = [-0.8601e5_dp, -3.8041e5_dp]
    poly( 19,:) = [-0.6702e5_dp, -4.1567e5_dp]
    poly( 20,:) = [-1.0990e5_dp, -4.1746e5_dp]
    poly( 21,:) = [-1.0944e5_dp, -4.5829e5_dp]
    poly( 22,:) = [-1.0065e5_dp, -5.0013e5_dp]
    poly( 23,:) = [-0.7625e5_dp, -5.3305e5_dp]
    poly( 24,:) = [-0.6571e5_dp, -5.7201e5_dp]
    poly( 25,:) = [-0.2665e5_dp, -5.9098e5_dp]
    poly( 26,:) = [-0.6042e5_dp, -6.1652e5_dp]
    poly( 27,:) = [-0.1109e5_dp, -6.1846e5_dp]
    poly( 28,:) = [-0.2399e5_dp, -6.6280e5_dp]
    poly( 29,:) = [-0.3598e5_dp, -7.0219e5_dp]
    poly( 30,:) = [-0.0160e5_dp, -6.7598e5_dp]
    poly( 31,:) = [ 0.1932e5_dp, -6.3562e5_dp]
    poly( 32,:) = [ 0.3721e5_dp, -6.7754e5_dp]
    poly( 33,:) = [ 0.7725e5_dp, -6.8812e5_dp]
    poly( 34,:) = [ 0.5899e5_dp, -7.2985e5_dp]
    poly( 35,:) = [ 0.6939e5_dp, -7.6885e5_dp]
    poly( 36,:) = [ 1.0711e5_dp, -7.9010e5_dp]
    poly( 37,:) = [ 1.4723e5_dp, -8.0512e5_dp]
    poly( 38,:) = [ 1.8838e5_dp, -8.2100e5_dp]
    poly( 39,:) = [ 2.3419e5_dp, -8.1880e5_dp]
    poly( 40,:) = [ 1.9120e5_dp, -8.2741e5_dp]
    poly( 41,:) = [ 2.2133e5_dp, -8.5626e5_dp]
    poly( 42,:) = [ 1.7064e5_dp, -8.5981e5_dp]
    poly( 43,:) = [ 1.3405e5_dp, -8.7757e5_dp]
    poly( 44,:) = [ 1.5893e5_dp, -9.2305e5_dp]
    poly( 45,:) = [ 1.2073e5_dp, -8.9282e5_dp]
    poly( 46,:) = [ 0.8931e5_dp, -9.2499e5_dp]
    poly( 47,:) = [ 0.4329e5_dp, -8.9716e5_dp]
    poly( 48,:) = [-0.0783e5_dp, -8.7450e5_dp]
    poly( 49,:) = [-0.5789e5_dp, -8.3810e5_dp]
    poly( 50,:) = [-1.3746e5_dp, -7.7720e5_dp]
    poly( 51,:) = [-1.7175e5_dp, -7.2582e5_dp]
    poly( 52,:) = [-2.0364e5_dp, -6.9912e5_dp]
    poly( 53,:) = [-2.3851e5_dp, -6.5499e5_dp]
    poly( 54,:) = [-2.7530e5_dp, -5.3553e5_dp]
    poly( 55,:) = [-2.8841e5_dp, -4.9232e5_dp]
    poly( 56,:) = [-2.9286e5_dp, -4.4711e5_dp]
    poly( 57,:) = [-3.0745e5_dp, -3.3945e5_dp]
    poly( 58,:) = [-3.2386e5_dp, -2.9319e5_dp]
    poly( 59,:) = [-3.2507e5_dp, -2.5020e5_dp]
    poly( 60,:) = [-3.3118e5_dp, -2.0832e5_dp]
    poly( 61,:) = [-3.2731e5_dp, -1.4841e5_dp]
    poly( 62,:) = [-3.1580e5_dp, -0.9411e5_dp]
    poly( 63,:) = [-2.7306e5_dp, -0.6904e5_dp]
    poly( 64,:) = [-2.3183e5_dp, -0.1091e5_dp]
    poly( 65,:) = [-2.2407e5_dp,  0.3389e5_dp]
    poly( 66,:) = [-1.9999e5_dp,  0.6980e5_dp]
    poly( 67,:) = [-2.0187e5_dp,  1.1338e5_dp]
    poly( 68,:) = [-2.0304e5_dp,  1.5620e5_dp]
    poly( 69,:) = [-1.8445e5_dp,  2.0050e5_dp]
    poly( 70,:) = [-1.8952e5_dp,  2.5343e5_dp]
    poly( 71,:) = [-2.0799e5_dp,  2.1778e5_dp]
    poly( 72,:) = [-2.0690e5_dp,  1.7765e5_dp]
    poly( 73,:) = [-2.1957e5_dp,  1.3968e5_dp]
    poly( 74,:) = [-2.1956e5_dp,  0.9968e5_dp]
    poly( 75,:) = [-2.3816e5_dp,  0.5980e5_dp]
    poly( 76,:) = [-2.6377e5_dp,  0.2888e5_dp]
    poly( 77,:) = [-3.0360e5_dp,  0.3382e5_dp]
    poly( 78,:) = [-3.0505e5_dp,  0.7408e5_dp]
    poly( 79,:) = [-3.0559e5_dp,  1.1415e5_dp]
    poly( 80,:) = [-2.7731e5_dp,  1.4247e5_dp]
    poly( 81,:) = [-2.6206e5_dp,  1.7961e5_dp]
    poly( 82,:) = [-2.6578e5_dp,  2.1960e5_dp]
    poly( 83,:) = [-2.6453e5_dp,  3.0531e5_dp]
    poly( 84,:) = [-2.3148e5_dp,  3.5499e5_dp]
    poly( 85,:) = [-2.3681e5_dp,  4.0124e5_dp]
    poly( 86,:) = [-2.2282e5_dp,  4.5178e5_dp]
    poly( 87,:) = [-2.1428e5_dp,  4.9543e5_dp]
    poly( 88,:) = [-1.9766e5_dp,  5.3229e5_dp]
    poly( 89,:) = [-1.8472e5_dp,  5.7164e5_dp]
    poly( 90,:) = [-1.5207e5_dp,  5.9587e5_dp]
    poly( 91,:) = [-1.1078e5_dp,  5.7979e5_dp]
    poly( 92,:) = [-1.1706e5_dp,  6.2322e5_dp]
    poly( 93,:) = [-1.5765e5_dp,  6.1397e5_dp]
    poly( 94,:) = [-1.1696e5_dp,  6.3771e5_dp]
    poly( 95,:) = [-1.3889e5_dp,  6.7127e5_dp]
    poly( 96,:) = [-0.9800e5_dp,  6.8870e5_dp]
    poly( 97,:) = [-1.2072e5_dp,  7.2285e5_dp]
    poly( 98,:) = [-1.0288e5_dp,  7.5925e5_dp]
    poly( 99,:) = [-1.0287e5_dp,  7.9975e5_dp]
    poly(100,:) = [-0.7675e5_dp,  8.3210e5_dp]
    poly(101,:) = [-0.6032e5_dp,  8.6911e5_dp]
    poly(102,:) = [-0.6162e5_dp,  9.1558e5_dp]
    poly(103,:) = [-0.5947e5_dp,  9.5816e5_dp]
    poly(104,:) = [ 0.0682e5_dp,  9.7742e5_dp]
    poly(105,:) = [ 0.1766e5_dp,  9.3773e5_dp]
    poly(106,:) = [-0.2248e5_dp,  9.5232e5_dp]
    poly(107,:) = [-0.2538e5_dp,  9.1150e5_dp]
    poly(108,:) = [-0.0421e5_dp,  8.7645e5_dp]
    poly(109,:) = [-0.2659e5_dp,  8.4246e5_dp]
    poly(110,:) = [-0.2557e5_dp,  8.0233e5_dp]
    poly(111,:) = [-0.2657e5_dp,  7.6161e5_dp]
    poly(112,:) = [-0.3503e5_dp,  7.2132e5_dp]
    poly(113,:) = [-0.1780e5_dp,  6.8387e5_dp]
    poly(114,:) = [-0.3278e5_dp,  6.4600e5_dp]
    poly(115,:) = [-0.1889e5_dp,  6.0735e5_dp]
    poly(116,:) = [-0.3383e5_dp,  5.6766e5_dp]
    poly(117,:) = [-0.3034e5_dp,  5.2726e5_dp]
    poly(118,:) = [-0.4553e5_dp,  4.8891e5_dp]
    poly(119,:) = [-0.4680e5_dp,  4.4792e5_dp]
    poly(120,:) = [-0.1740e5_dp,  4.2013e5_dp]
    poly(121,:) = [-0.1911e5_dp,  3.7940e5_dp]
    poly(122,:) = [-0.3552e5_dp,  3.4222e5_dp]
    poly(123,:) = [-0.1063e5_dp,  3.0976e5_dp]
    poly(124,:) = [-0.4183e5_dp,  2.8087e5_dp]
    poly(125,:) = [-0.7996e5_dp,  2.6792e5_dp]
    poly(126,:) = [-0.4217e5_dp,  2.5454e5_dp]
    poly(127,:) = [-0.7106e5_dp,  2.2639e5_dp]
    poly(128,:) = [-0.7213e5_dp,  1.8527e5_dp]
    poly(129,:) = [-0.8023e5_dp,  1.4559e5_dp]
    poly(130,:) = [-0.6977e5_dp,  1.0662e5_dp]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Patagonia

  SUBROUTINE calc_polygon_Narsarsuaq( poly)
    ! Return a polygon enveloping the Narsarsuaq area in Southern Greenland
    !
    ! (based on manual analysis of the Rignot velocity data,
    ! not meant for basin-integrated SMB stuff or such, but accurate
    ! enough to do region-specific mesh refinement)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Narsarsuaq'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 22,2))

    poly(  1,:) = [ 0.0182e6, -3.0978e6]
    poly(  2,:) = [ 0.0544e6, -3.1334e6]
    poly(  3,:) = [ 0.0550e6, -3.1469e6]
    poly(  4,:) = [ 0.0495e6, -3.1579e6]
    poly(  5,:) = [ 0.0556e6, -3.1634e6]
    poly(  6,:) = [ 0.0495e6, -3.1720e6]
    poly(  7,:) = [ 0.0354e6, -3.1781e6]
    poly(  8,:) = [ 0.0434e6, -3.2008e6]
    poly(  9,:) = [ 0.0403e6, -3.2162e6]
    poly( 10,:) = [ 0.0219e6, -3.2107e6]
    poly( 11,:) = [ 0.0035e6, -3.2174e6]
    poly( 12,:) = [-0.0131e6, -3.2217e6]
    poly( 13,:) = [-0.0247e6, -3.2254e6]
    poly( 14,:) = [-0.0775e6, -3.2015e6]
    poly( 15,:) = [-0.1075e6, -3.1518e6]
    poly( 16,:) = [-0.1088e6, -3.1285e6]
    poly( 17,:) = [-0.0990e6, -3.1064e6]
    poly( 18,:) = [-0.0830e6, -3.0953e6]
    poly( 19,:) = [-0.0511e6, -3.0800e6]
    poly( 20,:) = [-0.0321e6, -3.0708e6]
    poly( 21,:) = [-0.0180e6, -3.0555e6]
    poly( 22,:) = [ 0.0059e6, -3.0555e6]

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Narsarsuaq

  SUBROUTINE calc_polygon_Tijn_test_ISMIP_HOM_A( poly)
    ! Return a polygon enveloping the central square of the
    ! ISMIP-HOM Experiment A domain

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_polygon_Tijn_test_ISMIP_HOM_A'

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( poly( 4,2))

    poly(  1,:) = [-0.5_dp, -0.5_dp] * C%refgeo_idealised_ISMIP_HOM_L
    poly(  2,:) = [ 0.5_dp, -0.5_dp] * C%refgeo_idealised_ISMIP_HOM_L
    poly(  3,:) = [ 0.5_dp,  0.5_dp] * C%refgeo_idealised_ISMIP_HOM_L
    poly(  4,:) = [-0.5_dp,  0.5_dp] * C%refgeo_idealised_ISMIP_HOM_L

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_polygon_Tijn_test_ISMIP_HOM_A

! == Some fun and useful tools

  SUBROUTINE mesh_add_smileyface( mesh, res, width)
    ! Add a smileyface to a mesh. Because we can.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: res
    REAL(dp),                   INTENT(IN)        :: width

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'mesh_add_smileyface'
    INTEGER                                       :: i,n
    REAL(dp)                                      :: alpha_min, r, theta, x0, xw, y0, yw
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: line
    REAL(dp), DIMENSION(2)                        :: p

    ! Add routine to path
    CALL init_routine( routine_name)

    alpha_min = 25._dp * pi / 180._dp

    x0 = (mesh%xmax + mesh%xmin) / 2._dp
    xw = (mesh%xmax - mesh%xmin) / 2._dp
    y0 = (mesh%ymax + mesh%ymin) / 2._dp
    yw = (mesh%ymax - mesh%ymin) / 2._dp

    ! Smileyface - circle
    n = 50
    r = 0.75_dp
    ALLOCATE( line( n,4))
    DO i = 1, n
      theta = 2._dp * pi * REAL( i-1,dp) / REAL( n-1,dp)
      line( i,1:2) = [x0 + r * xw * COS( theta), y0 + yw * r * SIN( theta)]
      theta = 2._dp * pi * REAL( i  ,dp) / REAL( n-1,dp)
      line( i,3:4) = [x0 + r * xw * COS( theta), y0 + yw * r * SIN( theta)]
    END DO
    CALL refine_mesh_line( mesh, line, res, width, alpha_min)
    DEALLOCATE( line)

    ! Smileyface - smile
    n = 30
    r = 0.4_dp
    ALLOCATE( line( n,4))
    DO i = 1, n
      theta = -2._dp * pi * REAL( i-1,dp) / REAL( 2*n-1,dp)
      line( i,1:2) = [x0 + r * xw * COS( theta), y0 + yw * r * SIN( theta)]
      theta = -2._dp * pi * REAL( i  ,dp) / REAL( 2*n-1,dp)
      line( i,3:4) = [x0 + r * xw * COS( theta), y0 + yw * r * SIN( theta)]
    END DO
    CALL refine_mesh_line( mesh, line, res, width, alpha_min)
    DEALLOCATE( line)

    ! Smileyface - eyes
    p = [x0 - 0.3_dp * xw, y0 + 0.4_dp * yw]
    CALL refine_mesh_point( mesh, p, res, alpha_min)
    p = [x0 + 0.3_dp * xw, y0 + 0.4_dp * yw]
    CALL refine_mesh_point( mesh, p, res, alpha_min)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE mesh_add_smileyface

  SUBROUTINE mesh_add_UFEMISM_letters( mesh, res, width)
    ! Add the UFEMISM letters to a mesh. Because we can.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: res
    REAL(dp),                   INTENT(IN)        :: width

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'mesh_add_UFEMISM_letters'
    REAL(dp)                                      :: alpha_min
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: line
    INTEGER                                       :: i,n
    REAL(dp)                                      :: x0, xw, y0, yw

    ! Add routine to path
    CALL init_routine( routine_name)

    alpha_min = 25._dp * pi / 180._dp

    x0 = (mesh%xmax + mesh%xmin) / 2._dp
    xw = (mesh%xmax - mesh%xmin) / 2._dp
    y0 = (mesh%ymax + mesh%ymin) / 2._dp
    yw = (mesh%ymax - mesh%ymin) / 2._dp

    ! UFEMISM letters - normalised coordinates (mesh domain [-1,1,-1,1])
    n = 24
    ALLOCATE( line( n,4), source = 0._dp)

    line(  1,:) = [-0.80_dp,  0.20_dp, -0.80_dp, -0.20_dp]
    line(  2,:) = [-0.80_dp, -0.20_dp, -0.65_dp, -0.20_dp]
    line(  3,:) = [-0.65_dp, -0.20_dp, -0.65_dp,  0.20_dp]
    line(  4,:) = [-0.55_dp,  0.20_dp, -0.55_dp, -0.20_dp]
    line(  5,:) = [-0.55_dp,  0.20_dp, -0.40_dp,  0.20_dp]
    line(  6,:) = [-0.55_dp,  0.00_dp, -0.40_dp,  0.00_dp]
    line(  7,:) = [-0.30_dp,  0.20_dp, -0.30_dp, -0.20_dp]
    line(  8,:) = [-0.30_dp,  0.20_dp, -0.15_dp,  0.20_dp]
    line(  9,:) = [-0.30_dp,  0.00_dp, -0.15_dp,  0.00_dp]
    line( 10,:) = [-0.30_dp, -0.20_dp, -0.15_dp, -0.20_dp]
    line( 11,:) = [-0.05_dp,  0.20_dp, -0.05_dp, -0.20_dp]
    line( 12,:) = [-0.05_dp,  0.20_dp,  0.05_dp,  0.00_dp]
    line( 13,:) = [ 0.05_dp,  0.00_dp,  0.15_dp,  0.20_dp]
    line( 14,:) = [ 0.15_dp,  0.20_dp,  0.15_dp, -0.20_dp]
    line( 15,:) = [ 0.25_dp,  0.20_dp,  0.25_dp, -0.20_dp]
    line( 16,:) = [ 0.35_dp,  0.20_dp,  0.50_dp,  0.20_dp]
    line( 17,:) = [ 0.35_dp,  0.20_dp,  0.35_dp,  0.00_dp]
    line( 18,:) = [ 0.35_dp,  0.00_dp,  0.50_dp,  0.00_dp]
    line( 19,:) = [ 0.50_dp,  0.00_dp,  0.50_dp, -0.20_dp]
    line( 20,:) = [ 0.35_dp, -0.20_dp,  0.50_dp, -0.20_dp]
    line( 21,:) = [ 0.60_dp,  0.20_dp,  0.60_dp, -0.20_dp]
    line( 22,:) = [ 0.60_dp,  0.20_dp,  0.70_dp,  0.00_dp]
    line( 23,:) = [ 0.70_dp,  0.00_dp,  0.80_dp,  0.20_dp]
    line( 24,:) = [ 0.80_dp,  0.20_dp,  0.80_dp, -0.20_dp]

    ! Scale to actual mesh domain
    DO i = 1, n
      line( i,1) = x0 + xw * line( i,1)
      line( i,2) = y0 + yw * line( i,2)
      line( i,3) = x0 + xw * line( i,3)
      line( i,4) = y0 + yw * line( i,4)
    END DO

    CALL refine_mesh_line( mesh, line, res, width, alpha_min)

    ! Clean up after yourself
    DEALLOCATE( line)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE mesh_add_UFEMISM_letters

END MODULE mesh_refinement
