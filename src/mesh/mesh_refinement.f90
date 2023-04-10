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
                                                                     cross2
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
      meets_resolution_criterion = longest_leg <= res_max

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
        meets_resolution_criterion = longest_leg <= res_max
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
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: p_line        ! Collection of line segments
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: width         ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_line'
    INTEGER                                       :: nl
    INTEGER                                       :: ti,li,it
    REAL(dp), DIMENSION(2)                        :: pp,qq,pp2,qq2,pp_cropped,qq_cropped,dd
    LOGICAL                                       :: is_valid_line
    INTEGER                                       :: tip, tiq, n, via, vib, vic, tip2, tiq2
    REAL(dp), DIMENSION(2)                        :: va,vb,vc,llis
    LOGICAL                                       :: do_cross, crosses
    INTEGER                                       :: li_min, li_max
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    nl = SIZE( p_line,1)
    ! Safety
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
            IF (longest_leg > res_max) meets_resolution_criterion = .FALSE.
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
    INTEGER                                       :: ti_in
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
    CALL refine_mesh_line( mesh, p_line, res_max, res_max, alpha_min)

    ! Initialise the refinement stack with all triangles lying (partly) inside the polygon
    mesh%refinement_stackN = 0
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
          meets_resolution_criterion = longest_leg <= res_max

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_polygon

  SUBROUTINE refine_mesh_split_encroaching_triangles( mesh, alpha_min)
    ! As a last step in any mesh refinement, make sure no triangles exist anymore
    ! that encroach on the domain border (i.e. have their circumcenter lying
    ! outside the domain)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_split_encroaching_triangles'
    INTEGER                                       :: ti, n, vi, n_vertices_on_border
    INTEGER                                       :: via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_geometry_criterion
    LOGICAL                                       :: meets_encroachment_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the refinement stack with all border triangles
    mesh%refinement_stackN = 0
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

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_split_encroaching_triangles

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
    CALL refine_mesh_split_encroaching_triangles( mesh, alpha_min)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE Lloyds_algorithm_single_iteration

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
    REAL(dp)                                      :: alpha_min, r, theta, x0, xw, y0, yw, x, y
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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE mesh_add_UFEMISM_letters

END MODULE mesh_refinement
