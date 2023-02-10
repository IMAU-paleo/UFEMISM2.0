MODULE mesh_creation

  ! Routines used to create a mesh.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE reallocate_mod                                         , ONLY: reallocate
  USE math_utilities                                         , ONLY: segment_intersection, is_in_triangle, longest_triangle_leg, smallest_triangle_angle, &
                                                                     circumcenter, lies_on_line_segment, crop_line_to_domain, geometric_center, is_in_polygon
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: extend_mesh_primary, crop_mesh_primary
  USE mesh_utilities                                         , ONLY: update_triangle_circumcenter, find_containing_triangle
  USE mesh_Delaunay                                          , ONLY: split_triangle, remove_triangle_from_refinement_stack

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Mesh refinement based on different criteria

  SUBROUTINE refine_mesh_point( mesh, POI, res_max, alpha_min)
    ! Refine a mesh based on a 0-D point criterion

    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: POI           ! Collection of line segments
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_line'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: refinement_map
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: refinement_stack
    INTEGER                                       :: refinement_stackN
    INTEGER                                       :: ti_in
    INTEGER                                       :: ti, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Iteratively refine the mesh ==
  ! =================================

    ALLOCATE( refinement_map(   mesh%nTri_mem), source = 0)
    ALLOCATE( refinement_stack( mesh%nTri_mem), source = 0)
    refinement_stackN = 0

    ! Initialise the refinement stack with the triangle containing the point
    ti_in = 1
    CALL find_containing_triangle( mesh, POI, ti_in)
    refinement_map( ti_in) = 1
    refinement_stackN = 1
    refinement_stack( 1) = ti_in

    ! Keep refining until all triangles match the criterion
    DO WHILE (refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
        CALL reallocate( refinement_map  , mesh%nTri_mem)
        CALL reallocate( refinement_stack, mesh%nTri_mem)
      END IF

      ! Check the first (and therefore likely the biggest) triangle in the stack
      ti = refinement_stack( 1)

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

        CALL split_triangle( mesh, ti, p_new, &
          refinement_map    = refinement_map   , &
          refinement_stack  = refinement_stack , &
          refinement_stackN = refinement_stackN)

        ! Find out again which triangle contains the point, and add it to the stack
        CALL find_containing_triangle( mesh, POI, ti_in)
        IF (refinement_map( ti_in) == 0) THEN
          refinement_map( ti_in) = 1
          refinement_stackN = refinement_stackN + 1
          refinement_stack( refinement_stackN) = ti_in
        END IF

      ELSE
        ! Remove triangle ti from the refinement stack
        CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, ti)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Clean up after yourself
    DEALLOCATE( refinement_map  )
    DEALLOCATE( refinement_stack)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_point

  SUBROUTINE refine_mesh_line( mesh, line, res_max, alpha_min)
    ! Refine a mesh based on a 1-D line criterion

    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: line          ! Collection of line segments
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_line'
    INTEGER                                       :: nl
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE       :: Tri_li
    INTEGER                                       :: ti,li,it
    REAL(dp), DIMENSION(2)                        :: pp,qq,pp_cropped,qq_cropped,dd
    LOGICAL                                       :: is_valid_line
    INTEGER                                       :: tip, tiq, n, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va,vb,vc,llis
    LOGICAL                                       :: do_cross, crosses
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: refinement_map
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: refinement_stack
    INTEGER                                       :: refinement_stackN
    INTEGER                                       :: li_min, li_max
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

    nl = SIZE( line,1)
    ! Safety
    IF (SIZE( line,2) /= 4) CALL crash('line must be an n-by-4 array!')

  ! == Initialise triangle-line overlap ranges ==
  ! =============================================

    ! The array Tri_li has nTri-by-2 elements, and describes the range of line segments
    ! crossing through each triangle. E.g., if Tri_li( 5,:) = [175,179], that means that
    ! triangle 5 is crosses by line segments 175 to 179.
    ! Since the lines (grounding line, calving front, etc.) are constructed in such a way
    ! that the vast majority of line segments are contiguous, generally a triangle will
    ! only overlap with a small range of them.

    ALLOCATE( Tri_li( mesh%nTri_mem, 2), source = 0)

    Tri_li( 1:mesh%nTri,1) = nl + 1
    Tri_li( 1:mesh%nTri,2) = 0

    tip = 1
    tiq = 1

    DO li = 1, nl

      ! Line endpoints
      pp = [line( li,1), line( li,2)]
      qq = [line( li,3), line( li,4)]

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
      Tri_li( tip,1) = MIN( Tri_li( tip,1),li)
      Tri_li( tip,2) = MAX( Tri_li( tip,2),li)
      Tri_li( tiq,1) = MIN( Tri_li( tiq,1),li)
      Tri_li( tiq,2) = MAX( Tri_li( tiq,2),li)

      ! If they both lie inside the same triangle, no need to trace
      IF (tip == tiq) CYCLE

      ! p and q lie in different triangles - perform a pseudo-trace
      DO WHILE (NORM2( pp - qq) > res_max)

        ! Move pp in the direction of qq
        dd = qq - pp
        dd = dd / NORM2( dd)
        dd = dd * res_max
        pp = pp + dd

        ! Find the triangle that now contains pp
        CALL find_containing_triangle( mesh, pp, tip)

        ! Add li to the overlap lists of tip
        Tri_li( tip,1) = MIN( Tri_li( tip,1),li)
        Tri_li( tip,2) = MAX( Tri_li( tip,2),li)

        ! If we've reached tiq, the trace is done
        IF (tip == tiq) EXIT

      END DO ! DO WHILE (tip /= tiq)

    END DO ! DO li = 1, nl

  ! == Iteratively refine the mesh ==
  ! =================================

    ! If a triangle overlaps with any of the line segments, check if it is small enough.
    ! If not, split it. The new triangles will inherit the old one's line segment range.
    ! Update the (now reduced) overlap range for the new triangles.

    ALLOCATE( refinement_map(   mesh%nTri_mem), source = 0)
    ALLOCATE( refinement_stack( mesh%nTri_mem), source = 0)
    refinement_stackN = 0

    ! Mark which triangles need to be refined right now
    DO ti = 1, mesh%nTri
      IF (Tri_li( ti,2) == 0) THEN
        ! This triangle does not overlap with any line segments, so it doesn't need refining
        CYCLE
      ELSE
        ! Mark this triangle for refinement
        refinement_map( ti) = 1
        refinement_stackN = refinement_stackN + 1
        refinement_stack( refinement_stackN) = ti
      END IF
    END DO

    ! Keep refining until all triangles match the criterion
    DO WHILE (refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
        CALL reallocate( refinement_map  , mesh%nTri_mem   )
        CALL reallocate( refinement_stack, mesh%nTri_mem   )
        CALL reallocate( Tri_li          , mesh%nTri_mem, 2)
      END IF

      ! Check the first (and therefore likely the biggest) triangle in the stack
      ti = refinement_stack( 1)

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
        li_min = Tri_li( ti,1)
        li_max = Tri_li( ti,2)

        ! If that's already zero, skip
        IF (li_min == nl+1 .OR. li_max == 0) THEN
          ! No line overlap anyway

          meets_resolution_criterion = .TRUE.

        ELSE  ! IF (li_min == nl+1 .OR. li_max == 0) THEN
          ! Recalculate triangle overlap range

          Tri_li( ti,:) = [nl+1,0]

          DO li = li_min, li_max

            ! Line endpoints
            pp = [line( li,1), line( li,2)]
            qq = [line( li,3), line( li,4)]

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
            crosses = crosses .OR. lies_on_line_segment( pp, qq, va, mesh%tol_dist)
            crosses = crosses .OR. lies_on_line_segment( pp, qq, vb, mesh%tol_dist)
            crosses = crosses .OR. lies_on_line_segment( pp, qq, vc, mesh%tol_dist)

            ! If so, add it to the overlap range
            IF (crosses) THEN
              Tri_li( ti,1) = MIN( Tri_li( ti,1),li)
              Tri_li( ti,2) = MAX( Tri_li( ti,2),li)
            END IF

          END DO ! DO li = li_min, li_max

          ! If this triangle overlaps with any line segments,
          ! check if it meets the resolution criterion

          meets_resolution_criterion = .TRUE.
          IF (Tri_li( ti,1) <= Tri_li( ti,2)) THEN
            longest_leg = longest_triangle_leg( va, vb, vc)
            IF (longest_leg > res_max) meets_resolution_criterion = .FALSE.
          END IF

        END IF ! IF (li_min == nl+1 .OR. li_max == 0) THEN

      END IF ! IF (.NOT. meets_geometry_criterion) THEN

      ! If either of the two criteria is not met, split the triangle
      IF (.NOT. meets_geometry_criterion .OR. .NOT. meets_resolution_criterion) THEN
        ! Split triangle ti at its circumcenter

        p_new = circumcenter( va, vb, vc)

        CALL split_triangle( mesh, ti, p_new, &
          refinement_map    = refinement_map   , &
          refinement_stack  = refinement_stack , &
          refinement_stackN = refinement_stackN, &
          Tri_li            = Tri_li)

      ELSE
        ! Remove triangle ti from the refinement stack
        CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, ti)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Clean up after yourself
    DEALLOCATE( refinement_map  )
    DEALLOCATE( refinement_stack)
    DEALLOCATE( Tri_li          )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_line

  SUBROUTINE refine_mesh_polygon( mesh, poly, res_max, alpha_min)
    ! Refine a mesh based on a 2-D polygon criterion

    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh          ! The mesh that should be refined
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: poly          ! Polygon
    REAL(dp),                   INTENT(IN)        :: res_max       ! Maximum allowed resolution for triangles crossed by any of these line segments
    REAL(dp),                   INTENT(IN)        :: alpha_min     ! Minimum allowed internal triangle angle

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'refine_mesh_line'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: refinement_map
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: refinement_stack
    INTEGER                                       :: refinement_stackN
    INTEGER                                       :: ti_in
    INTEGER                                       :: ti, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc, gc, vab, vbc, vca
    REAL(dp)                                      :: longest_leg, smallest_angle
    LOGICAL                                       :: meets_resolution_criterion
    LOGICAL                                       :: meets_geometry_criterion
    REAL(dp), DIMENSION(2)                        :: p_new

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Iteratively refine the mesh ==
  ! =================================

    ALLOCATE( refinement_map(   mesh%nTri_mem), source = 0)
    ALLOCATE( refinement_stack( mesh%nTri_mem), source = 0)
    refinement_stackN = 0

    ! Initialise the refinement stack with all triangles lying (partly) inside the polygon
    DO ti = 1, mesh%nTri

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      va  = mesh%V( via,:)
      vb  = mesh%V( vib,:)
      vc  = mesh%V( vic,:)

      ! Its geometric centre and edge midpoints
      gc = geometric_center( va, vb, vc)
      vab = (va + vb) / 2._dp
      vbc = (vb + vc) / 2._dp
      vca = (vc + va) / 2._dp

      IF (is_in_polygon( poly, va) .OR. &
          is_in_polygon( poly, vb) .OR. &
          is_in_polygon( poly, vc) .OR. &
          is_in_polygon( poly, gc) .OR. &
          is_in_polygon( poly, vab) .OR. &
          is_in_polygon( poly, vbc) .OR. &
          is_in_polygon( poly, vca)) THEN
        refinement_map( ti) = 1
        refinement_stackN = refinement_stackN + 1
        refinement_stack( refinement_stackN) = ti
      END IF

    END DO ! DO ti = 1, mesh%nTri

    ! Keep refining until all triangles match the criterion
    DO WHILE (refinement_stackN > 0)

      ! If needed, allocate more memory for the mesh
      IF (mesh%nV > mesh%nV_mem - 10 .OR. mesh%nTri > mesh%nTri_mem - 10) THEN
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
        CALL reallocate( refinement_map  , mesh%nTri_mem)
        CALL reallocate( refinement_stack, mesh%nTri_mem)
      END IF

      ! Check the first (and therefore likely the biggest) triangle in the stack
      ti = refinement_stack( 1)

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      va  = mesh%V( via,:)
      vb  = mesh%V( vib,:)
      vc  = mesh%V( vic,:)

      ! Its geometric centre and edge midpoints
      gc = geometric_center( va, vb, vc)
      vab = (va + vb) / 2._dp
      vbc = (vb + vc) / 2._dp
      vca = (vc + va) / 2._dp

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
            is_in_polygon( poly, vc) .OR. &
            is_in_polygon( poly, gc) .OR. &
            is_in_polygon( poly, vab) .OR. &
            is_in_polygon( poly, vbc) .OR. &
            is_in_polygon( poly, vca)) THEN
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

        CALL split_triangle( mesh, ti, p_new, &
          refinement_map    = refinement_map   , &
          refinement_stack  = refinement_stack , &
          refinement_stackN = refinement_stackN)

      ELSE
        ! Remove triangle ti from the refinement stack
        CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, ti)
      END IF

    END DO ! DO WHILE (refinement_stackN > 0)

    ! Crop surplus mesh memory
    CALL crop_mesh_primary( mesh)

    ! Clean up after yourself
    DEALLOCATE( refinement_map  )
    DEALLOCATE( refinement_stack)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE refine_mesh_polygon

! == Initialise the five-vertex dummy mesh

  SUBROUTINE initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)
    ! Initialises a 5-vertex, 4-triangle "dummy"  mesh:
    !
    !   v4 - - - - - - - - v3   V          nC     C             niTri   iTri          edge_index
    !   | \              / |    -1 -1      3      2  5  4        2      1  4            6
    !   |  \            /  |     1 -1      3      3  5  1        2      2  1            4
    !   |   \    t3    /   |     1  1      3      4  5  2        2      3  2            2
    !   |    \        /    |    -1  1      3      1  5  3        2      4  3            8
    !   |     \      /     |     0  0      4      1  2  3  4     4      1  2  3  4      0
    !   |      \    /      |
    !   |       \  /       |    Tri           TriC
    !   |  t4    v5    t2  |    1  2  5      2  4  0
    !   |       /  \       |    2  3  5      3  1  0
    !   |      /    \      |    3  4  5      4  2  0
    !   |     /      \     |    4  1  5      1  3  0
    !   |    /        \    |
    !   |   /    t1    \   |
    !   |  /            \  |
    !   | /              \ |
    !   v1 - - - - - - - - v2
    !
    ! NOTE: memory must already be allocated for the mesh before calling this routine

    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: xmin    ! Mesh domain
    REAL(dp),                   INTENT(IN)        :: xmax
    REAL(dp),                   INTENT(IN)        :: ymin
    REAL(dp),                   INTENT(IN)        :: ymax

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_dummy_mesh'
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Meta properties
    mesh%xmin                 = xmin    ! Boundaries of the square domain.
    mesh%xmax                 = xmax
    mesh%ymin                 = ymin
    mesh%ymax                 = ymax

    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist   = ((mesh%xmax - mesh%xmin) + (mesh%ymax-mesh%ymin)) * tol / 2._dp

    ! The four corners, plus one central vertex.
    mesh%nV              = 5

    mesh%V               = 0._dp
    mesh%V( 1,:)         = [xmin, ymin]
    mesh%V( 2,:)         = [xmax, ymin]
    mesh%V( 3,:)         = [xmax, ymax]
    mesh%V( 4,:)         = [xmin, ymax]
    mesh%V( 5,:)         = [(xmin+xmax)/2._dp, (ymin+ymax)/2._dp]

    mesh%edge_index      = 0
    mesh%edge_index(1:5) = [6, 4, 2, 8, 0]

    mesh%nC              = 0
    mesh%nC( 1:5)        = [3, 3, 3, 3, 4]

    mesh%C               = 0
    mesh%C( 1,1:4)       = [2, 5, 4, 0]
    mesh%C( 2,1:4)       = [3, 5, 1, 0]
    mesh%C( 3,1:4)       = [4, 5, 2, 0]
    mesh%C( 4,1:4)       = [1, 5, 3, 0]
    mesh%C( 5,1:4)       = [1, 2, 3, 4]

    mesh%niTri           = 0
    mesh%niTri( 1:5)     = [2, 2, 2, 2, 4]

    mesh%iTri            = 0
    mesh%iTri( 1,1:4)    = [1, 4, 0, 0]
    mesh%iTri( 2,1:4)    = [2, 1, 0, 0]
    mesh%iTri( 3,1:4)    = [3, 2, 0, 0]
    mesh%iTri( 4,1:4)    = [4, 3, 0, 0]
    mesh%iTri( 5,1:4)    = [1, 2, 3, 4]

    mesh%nTri            = 4

    mesh%Tri             = 0
    mesh%Tri( 1,:)       = [1, 2, 5]
    mesh%Tri( 2,:)       = [2, 3, 5]
    mesh%Tri( 3,:)       = [3, 4, 5]
    mesh%Tri( 4,:)       = [4, 1, 5]

    mesh%TriC            = 0
    mesh%TriC( 1,:)      = [2, 4, 0]
    mesh%TriC( 2,:)      = [3, 1, 0]
    mesh%TriC( 3,:)      = [4, 2, 0]
    mesh%TriC( 4,:)      = [1, 3, 0]

    mesh%TriCC = 0._dp
    CALL update_triangle_circumcenter( mesh, 1)
    CALL update_triangle_circumcenter( mesh, 2)
    CALL update_triangle_circumcenter( mesh, 3)
    CALL update_triangle_circumcenter( mesh, 4)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_dummy_mesh

END MODULE mesh_creation
