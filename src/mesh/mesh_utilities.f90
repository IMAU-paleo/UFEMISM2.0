MODULE mesh_utilities

  ! Generally useful functions used in mesh creation and updating.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh
  USE math_utilities                                         , ONLY: geometric_center, is_in_triangle, lies_on_line_segment, circumcenter, &
                                                                     line_from_points, line_line_intersection

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Finding the vertices of a vertex' Voronoi cell

  SUBROUTINE calc_Voronoi_cell_vertices(        mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a vertex's Voronoi cell

    IMPLICIT NONE

    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION( mesh%nC_mem,2), INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables
    INTEGER                                            :: vvi

    IF (mesh%edge_index( vi) == 0) THEN
      ! Free vertex

      CALL calc_Voronoi_cell_vertices_free(mesh, vi, Vor, nVor)

    ELSEIF (mesh%edge_index( vi) == 2 .OR. &
            mesh%edge_index( vi) == 4 .OR. &
            mesh%edge_index( vi) == 6 .OR. &
            mesh%edge_index( vi) == 8) THEN
      ! Corner vertex

      CALL calc_Voronoi_cell_vertices_corner(mesh, vi, Vor, nVor)

    ELSE
      ! Boundary vertex

      CALL calc_Voronoi_cell_vertices_edge(mesh, vi, Vor, nVor)

    END IF

    ! Safety: sometimes, Voronoi vertices end up just very slightly
    !         outside the mesh domain; move them to the boundary.
    DO vvi = 1, nVor

      ! Safety: if a Voronoi vertex is too far outside the mesh domain, crash.
      IF (Vor(vvi,1) < mesh%xmin - mesh%tol_dist .OR. &
          Vor(vvi,1) > mesh%xmax + mesh%tol_dist .OR. &
          Vor(vvi,2) < mesh%ymin - mesh%tol_dist .OR. &
          Vor(vvi,2) > mesh%ymax + mesh%tol_dist) THEN
        CALL crash('find_Voronoi_cell_vertices: found Voronoi cell vertex outside of mesh domain!')
      END IF

      Vor(vvi,1) = MAX( MIN( Vor(vvi,1), mesh%xmax), mesh%xmin)
      Vor(vvi,2) = MAX( MIN( Vor(vvi,2), mesh%ymax), mesh%ymin)

    END DO

  END SUBROUTINE calc_Voronoi_cell_vertices

  SUBROUTINE calc_Voronoi_cell_vertices_free(   mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a free vertex's Voronoi cell

    IMPLICIT NONE

    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION( mesh%nC_mem,2), INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables
    INTEGER                                            :: iti, iti_clock, iti_anti, ti, ti_clock, ti_anti
    REAL(dp), DIMENSION(2)                             :: cc, cc_clock, cc_anti

    Vor  = 0._dp
    nVor = 0

    DO iti = 1, mesh%niTri(vi)

      ! Find indices of current, clockwise neighbouring and anticlockwise
      ! neighbouring triangles.
      iti_clock = iti - 1
      IF (iti_clock == 0) iti_clock = mesh%niTri( vi)
      iti_anti  = iti + 1
      IF (iti_anti > mesh%niTri( vi)) iti_anti = 1

      ti       = mesh%iTri( vi,iti)
      ti_clock = mesh%iTri( vi,iti_clock)
      ti_anti  = mesh%iTri( vi,iti_anti)

      ! If necessary, crop (split) the circumcenter of the current triangle.
      cc       = mesh%Tricc( ti,:)
      CALL crop_circumcenter( mesh, ti, ti_clock, cc_clock)
      CALL crop_circumcenter( mesh, ti, ti_anti,  cc_anti)

      ! Add the resulting Voronoi vertex/vertices
      IF (cc_clock( 1) /= cc_anti( 1) .OR. cc_clock( 2) /= cc_anti( 2)) THEN
        nVor = nVor+1
        Vor( nVor,:) = cc_clock
        nVor = nVor+1
        Vor( nVor,:) = cc_anti
      ELSE
        nVor = nVor+1
        Vor( nVor,:) = cc
      END IF

    END DO ! DO t = 1, mesh%niTri(vi)

    ! Repeat the first Voronoi vertex
    nVor = nVor+1
    Vor( nVor,:) = Vor( 1,:)

  END SUBROUTINE calc_Voronoi_cell_vertices_free

  SUBROUTINE calc_Voronoi_cell_vertices_edge(   mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up an edge vertex's Voronoi cell

    IMPLICIT NONE

    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION( mesh%nC_mem,2), INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables
    INTEGER                                            :: iti
    REAL(dp), DIMENSION(2)                             :: cc, cc_clock, cc_anti, cc_cropped

    Vor  = 0._dp
    nVor = 0

    ! == Boundary cell ==
    ! If the first or last circumcenter lies outside of the grid, crop it.
    ! If not, add the point on the edge closest to that circumcenter as an additional Voronoi cell vertex.

    DO iti = 1, mesh%niTri( vi)

      cc = mesh%Tricc(mesh%iTri( vi,iti),:)

      IF (iti  ==  1) THEN
        ! Start by possibly adding the boundary projection of the vertex
        IF     ((mesh%edge_index( vi) == 1 .OR. mesh%edge_index( vi) == 2) .AND. cc( 2)<mesh%ymax) THEN
          nVor = nVor+1
          Vor( nVor,:) = [cc( 1), mesh%ymax]
        ELSEIF ((mesh%edge_index( vi) == 3 .OR. mesh%edge_index( vi) == 4) .AND. cc( 1)<mesh%xmax) THEN
          nVor = nVor+1
          Vor( nVor,:) = [mesh%xmax, cc( 2)]
        ELSEIF ((mesh%edge_index( vi) == 5 .OR. mesh%edge_index( vi) == 6) .AND. cc( 2)>mesh%ymin) THEN
          nVor = nVor+1
          Vor( nVor,:) = [cc( 1), mesh%ymin]
        ELSEIF ((mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 8) .AND. cc( 1)>mesh%xmin) THEN
          nVor = nVor+1
          Vor( nVor,:) = [mesh%xmin, cc( 2)]
        END IF

        ! Then add the (possibly cropped) vertex
        CALL crop_circumcenter(mesh, mesh%iTri( vi,1), mesh%iTri( vi,2), cc_cropped)
        nVor = nVor+1
        Vor( nVor,:) = cc_cropped
      END IF ! IF (iti  ==  1) THEN

      IF (iti > 1 .AND. iti < mesh%niTri( vi)) THEN
        ! Split the circumcenter
        CALL crop_circumcenter(mesh, mesh%iTri( vi,iti), mesh%iTri( vi,iti+1), cc_anti)
        CALL crop_circumcenter(mesh, mesh%iTri( vi,iti), mesh%iTri( vi,iti-1), cc_clock)
        IF (cc_anti( 1) /= cc_clock( 1) .OR. cc_anti( 2) /= cc_clock( 2)) THEN
          nVor = nVor+1
          Vor( nVor,:) = cc_clock
          nVor = nVor+1
          Vor( nVor,:) = cc_anti
        ELSE
          nVor = nVor+1
          Vor( nVor,:) = cc
        END IF
      END IF ! IF (iti > 1 .AND. iti < mesh%niTri( vi)) THEN

      IF (iti  ==  mesh%niTri( vi)) THEN
        ! First add the (possibly cropped) vertex
        CALL crop_circumcenter(mesh, mesh%iTri( vi,iti), mesh%iTri( vi,iti-1), cc_cropped)
        nVor = nVor+1
        Vor( nVor,:) = cc_cropped

        ! Then possibly add the boundary projection of the vertex
        IF     ((mesh%edge_index( vi) == 1 .OR. mesh%edge_index( vi) == 8) .AND. cc( 2)<mesh%ymax) THEN
          nVor = nVor+1
          Vor( nVor,:) = [cc( 1), mesh%ymax]
        ELSEIF ((mesh%edge_index( vi) == 3 .OR. mesh%edge_index( vi) == 2) .AND. cc( 1)<mesh%xmax) THEN
          nVor = nVor+1
          Vor( nVor,:) = [mesh%xmax, cc( 2)]
        ELSEIF ((mesh%edge_index( vi) == 5 .OR. mesh%edge_index( vi) == 4) .AND. cc( 2)>mesh%ymin) THEN
          nVor = nVor+1
          Vor( nVor,:) = [cc( 1), mesh%ymin]
        ELSEIF ((mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 6) .AND. cc( 1)>mesh%xmin) THEN
          nVor = nVor+1
          Vor( nVor,:) = [mesh%xmin, cc( 2)]
        END IF

      END IF ! IF (iti  ==  mesh%niTri( vi)) THEN

    END DO ! DO n = 1, mesh%niTri( vi)

  END SUBROUTINE calc_Voronoi_cell_vertices_edge

  SUBROUTINE calc_Voronoi_cell_vertices_corner( mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a corner vertex's Voronoi cell

    IMPLICIT NONE

    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION( mesh%nC_mem,2), INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables
    REAL(dp), DIMENSION(2)                             :: cc

    Vor  = 0._dp
    nVor = 0

    IF (mesh%niTri( vi) > 1) THEN
      ! This corner vertex has more than one triangle, can be handled by Edge version

      CALL calc_Voronoi_cell_vertices_edge(mesh, vi, Vor, nVor)

      IF     (mesh%edge_index( vi) == 2) THEN
        ! Northeast corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmax, mesh%ymax]
      ELSEIF (mesh%edge_index( vi) == 4) THEN
        ! Southeast corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmax, mesh%ymin]
      ELSEIF (mesh%edge_index( vi) == 6) THEN
        ! Southwest corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmin, mesh%ymin]
      ELSEIF (mesh%edge_index( vi) == 8) THEN
        ! Northwest corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmin, mesh%ymax]
      END IF

    ELSE
      ! This corner vertex has only a single triangle, best handled manually

      cc = mesh%Tricc(mesh%iTri( vi,1),:)

      IF     (mesh%edge_index( vi)==2) THEN
        ! Northeast corner
        nVor = 3
        Vor( 1,:) = [cc( 1), mesh%ymax]
        Vor( 2,:) = cc
        Vor( 3,:) = [mesh%xmax, cc( 2)]
      ELSEIF (mesh%edge_index( vi)==4) THEN
        ! Southeast corner
        nVor = 3
        Vor( 1,:) = [mesh%xmax, cc( 2)]
        Vor( 2,:) = cc
        Vor( 3,:) = [cc( 1), mesh%ymin]
      ELSEIF (mesh%edge_index( vi)==6) THEN
        ! Southwest corner
        nVor = 3
        Vor( 1,:) = [cc( 1), mesh%ymin]
        Vor( 2,:) = cc
        Vor( 3,:) = [mesh%xmin, cc( 2)]
      ELSEIF (mesh%edge_index( vi)==8) THEN
        ! Northwest corner
        nVor = 3
        Vor( 1,:) = [mesh%xmin, cc( 2)]
        Vor( 2,:) = cc
        Vor( 3,:) = [cc( 1), mesh%ymax]
      ELSE
        CALL crash('calc_Voronoi_cell_vertices_corner: a non-corner vertex has only one triangle? This cannot be!')
      END IF ! IF (mesh%edge_index( vi)==2) THEN

    END IF

  END SUBROUTINE calc_Voronoi_cell_vertices_corner

  SUBROUTINE crop_circumcenter( mesh, t1, t2, ccc)
    ! Crop the circumcenter of triangle t1 in the direction of t2

    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: t1, t2
    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: ccc

    REAL(dp)                                      :: la, lb, lc, le, lf, lg
    REAL(dp), DIMENSION(2)                        :: p, q

    ccc  = mesh%Tricc( t1,:)

    IF     (mesh%Tri_edge_index( t1) == 1 .AND. mesh%Tricc( t1,2) > mesh%ymax) THEN
      ! North boundary triangle
      CALL line_from_points([mesh%xmin, mesh%ymax], [mesh%xmax, mesh%ymax], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index( t1) == 3 .AND. mesh%Tricc( t1,1) > mesh%xmax) THEN
      ! East boundary triangle
      CALL line_from_points([mesh%xmax, mesh%ymax], [mesh%xmax, mesh%ymin], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index( t1) == 5 .AND. mesh%Tricc( t1,2) < mesh%ymin) THEN
      ! South boundary triangle
      CALL line_from_points([mesh%xmin, mesh%ymin], [mesh%xmax, mesh%ymin], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index( t1) == 7 .AND. mesh%Tricc( t1,1) < mesh%xmin) THEN
      ! West boundary triangle
      CALL line_from_points([mesh%xmin, mesh%ymax], [mesh%xmin, mesh%ymin], la, lb, lc)
    ELSE
      RETURN
    END IF

   p = mesh%Tricc( t1,:)
   q = mesh%Tricc( t2,:)
   CALL line_from_points( p, q, le, lf, lg)
   CALL line_line_intersection( la, lb, lc, le, lf, lg, ccc)

  END SUBROUTINE crop_circumcenter

!  SUBROUTINE calc_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
!    ! Return the endpoints of the shared Voronoi cell boundary represented by edge aci
!
!    IMPLICIT NONE
!
!    ! In/output variables
!    TYPE(type_mesh),          INTENT(IN)          :: mesh
!    INTEGER,                  INTENT(IN)          :: aci
!    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: cc1, cc2
!
!    ! Local variables
!    INTEGER                                       :: til,tir
!
!    til = mesh%Aci( aci,5)
!    tir = mesh%Aci( aci,6)
!
!    IF (mesh%edge_index_Ac( aci) > 0) THEN
!      ! Boundary segments have only one adjacent triangle
!
!      IF (til > 0) THEN
!        cc1 = mesh%Tricc( til,:)
!      ELSE
!        cc1 = mesh%Tricc( tir,:)
!      END IF
!      IF     (mesh%edge_index_Ac( aci) == 1) THEN
!        ! North
!        cc2 = [cc1(1), mesh%ymax]
!      ELSEIF (mesh%edge_index_Ac( aci) == 3) THEN
!        ! East
!        cc2 = [mesh%xmax, cc1(2)]
!      ELSEIF (mesh%edge_index_Ac( aci) == 5) THEN
!        ! South
!        cc2 = [cc1(1), mesh%ymin]
!      ELSEIF (mesh%edge_index_Ac( aci) == 7) THEN
!        ! West
!        cc2 = [mesh%xmin, cc1(2)]
!      END IF
!
!    ELSE ! IF (mesh%edge_index_Ac( aci) > 0) THEN
!
!      cc1 = mesh%Tricc( til,:)
!      cc2 = mesh%Tricc( tir,:)
!
!    END IF ! IF (mesh%edge_index_Ac( aci) > 0) THEN
!
!  END SUBROUTINE calc_shared_Voronoi_boundary

! == Some basic geometrical operations

  PURE FUNCTION is_boundary_segment( mesh, v1, v2) RESULT(isso)
   ! Determine whether or not the line between two vertices is an Edge segment

    IMPLICIT NONE

   TYPE(type_mesh),          INTENT(IN)          :: mesh
   INTEGER,                  INTENT(IN)          :: v1, v2
   LOGICAL                                       :: isso

   IF (mesh%edge_index( v1) == 0 .OR. mesh%edge_index( v2) == 0) THEN
     isso = .FALSE.
     RETURN
   END IF

   isso = .FALSE.

   IF (mesh%edge_index( v1) == 1) THEN
     IF (mesh%edge_index( v2) == 8 .OR. &
         mesh%edge_index( v2) == 1 .OR. &
         mesh%edge_index( v2) == 2) THEN
       isso = .TRUE.
     END IF
   ELSEIF (mesh%edge_index( v1) == 2) THEN
     IF (mesh%edge_index( v2) == 8 .OR. &
         mesh%edge_index( v2) == 1 .OR. &
         mesh%edge_index( v2) == 3 .OR. &
         mesh%edge_index( v2) == 4) THEN
       isso = .TRUE.
     END IF
   ELSEIF (mesh%edge_index( v1) == 3) THEN
     IF (mesh%edge_index( v2) == 2 .OR. &
         mesh%edge_index( v2) == 3 .OR. &
         mesh%edge_index( v2) == 4) THEN
       isso = .TRUE.
     END IF
   ELSEIF (mesh%edge_index( v1) == 4) THEN
     IF (mesh%edge_index( v2) == 2 .OR. &
         mesh%edge_index( v2) == 3 .OR. &
         mesh%edge_index( v2) == 5 .OR. &
         mesh%edge_index( v2) == 6) THEN
       isso = .TRUE.
     END IF
   ELSEIF (mesh%edge_index( v1) == 5) THEN
     IF (mesh%edge_index( v2) == 4 .OR. &
         mesh%edge_index( v2) == 5 .OR. &
         mesh%edge_index( v2) == 6) THEN
       isso = .TRUE.
     END IF
   ELSEIF (mesh%edge_index( v1) == 6) THEN
     IF (mesh%edge_index( v2) == 4 .OR. &
         mesh%edge_index( v2) == 5 .OR. &
         mesh%edge_index( v2) == 7 .OR. &
         mesh%edge_index( v2) == 8) THEN
       isso = .TRUE.
     END IF
   ELSEIF (mesh%edge_index( v1) == 7) THEN
     IF (mesh%edge_index( v2) == 6 .OR. &
         mesh%edge_index( v2) == 7 .OR. &
         mesh%edge_index( v2) == 8) THEN
       isso = .TRUE.
     END IF
   ELSEIF (mesh%edge_index( v1) == 8) THEN
     IF (mesh%edge_index( v2) == 6 .OR. &
         mesh%edge_index( v2) == 7 .OR. &
         mesh%edge_index( v2) == 1 .OR. &
         mesh%edge_index( v2) == 2) THEN
       isso = .TRUE.
     END IF
   END IF

  END FUNCTION is_boundary_segment

  PURE FUNCTION is_walltowall( mesh, ti) RESULT(isso)
   ! Determine whether or not a triangle is "wall to wall"
   ! (i.e. contains vertices lying on opposite domain boundaries)

    IMPLICIT NONE

   TYPE(type_mesh),          INTENT(IN)          :: mesh
   INTEGER,                  INTENT(IN)          :: ti
   LOGICAL                                       :: isso

   LOGICAL                                       :: has_north, has_south, has_east, has_west
   INTEGER                                       :: n, vi

   has_north = .FALSE.
   has_south = .FALSE.
   has_east  = .FALSE.
   has_west  = .FALSE.

   DO n = 1, 3
     vi = mesh%Tri( ti,n)
     IF (mesh%edge_index( vi) == 1) has_north = .TRUE.
     IF (mesh%edge_index( vi) == 3) has_east  = .TRUE.
     IF (mesh%edge_index( vi) == 5) has_south = .TRUE.
     IF (mesh%edge_index( vi) == 7) has_west  = .TRUE.
   END DO

   isso = .FALSE.
   IF (has_north .AND. has_south) isso = .TRUE.
   IF (has_west  .AND. has_east ) isso = .TRUE.

  END FUNCTION is_walltowall

  SUBROUTINE update_triangle_circumcenter( mesh, ti)
    ! Calculate the circumcenter of mesh triangle ti

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: v1, v2, v3, cc

    v1 = mesh%V( mesh%Tri( ti,1),:)
    v2 = mesh%V( mesh%Tri( ti,2),:)
    v3 = mesh%V( mesh%Tri( ti,3),:)
    cc = circumcenter( v1, v2, v3)

    ! If find_circumcenter yields infinity, it's because p and q have the
    ! same y-coordinate. Rearrange vertices in triangle matrix (maintaining
    ! counter-clockwise orientation) and try again.

    IF (cc( 1) > (mesh%xmax - mesh%xmin) * 1E5_dp .OR. cc( 2) > (mesh%ymax - mesh%ymin) * 1E5_dp) THEN
      mesh%Tri( ti,:) = [mesh%Tri( ti,2), mesh%Tri( ti,3), mesh%Tri( ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:)
      cc = circumcenter( v1, v2, v3)
    END IF

    IF (cc( 1) > (mesh%xmax - mesh%xmin) * 1E5_dp .OR. cc( 2) > (mesh%ymax - mesh%ymin) * 1E5_dp) THEN
      mesh%Tri( ti,:) = [mesh%Tri( ti,2), mesh%Tri( ti,3), mesh%Tri( ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:)
      cc = circumcenter( v1, v2, v3)
    END IF

    IF (cc(1) > (mesh%xmax - mesh%xmin) * 1E5_dp .OR. cc( 2) > (mesh%ymax - mesh%ymin) * 1E5_dp) THEN
      CALL warning('update_triangle_circumcenter: triangle  doesn yield a valid circumcenter!')
    END IF

    mesh%Tricc( ti,:) = cc

  END SUBROUTINE update_triangle_circumcenter

! == Some basic search operations on a mesh

  SUBROUTINE find_containing_triangle( mesh, p, ti_in)
    ! Find the triangle containing the point p. First do a "linear search":
    ! Start at initial guess ti_in. Check all neighbours of ti_in, find the one
    ! closest to p, select that one as the new ti_in. Repeat until all neighbours
    ! of ti_in are further away from p than ti_in itself.
    ! Then (if needed) search outward from there using a flood-fill algorithm.

    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: ti_in

    REAL(dp), DIMENSION(2)                        :: qq, rr, ss
    INTEGER                                       :: ncycle, t_prev
    REAL(dp), DIMENSION(2)                        :: gcti, gctc
    REAL(dp)                                      :: d, dc, dcmin
    INTEGER                                       :: tc, tcmin
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: map, stack1, stack2
    INTEGER                                       :: stackN1, stackN2
    LOGICAL                                       :: FoundIt
    INTEGER                                       :: n, ti, n2, tin

    ! If p lies outside the mesh domain, throw an error
    IF (p( 1) < mesh%xmin .OR. p( 1) > mesh%xmax .OR. p( 2) < mesh%ymin .OR. p( 2) > mesh%ymax) THEN
      CALL crash('find_containing_triangle - ERROR: point lies outside mesh domain!')
    END IF

    ! See if the initial guess is correct.
    qq = mesh%V( mesh%Tri( ti_in,1),:)
    rr = mesh%V( mesh%Tri( ti_in,2),:)
    ss = mesh%V( mesh%Tri( ti_in,3),:)
    IF (is_in_triangle( qq,rr,ss,p)) RETURN

    ! If not, start with a linear search.

  ! == Linear search ==
  ! ===================

    ncycle = 0
    t_prev = ti_in
    DO WHILE (ncycle < mesh%nTri)

      qq = mesh%V( mesh%Tri( ti_in,1),:)
      rr = mesh%V( mesh%Tri( ti_in,2),:)
      ss = mesh%V( mesh%Tri( ti_in,3),:)
      gcti = geometric_center( qq,rr,ss)
      d = NORM2( gcti - p)

      dcmin = d + 10._dp
      tcmin = 0
      DO n = 1, 3
        tc   = mesh%TriC( ti_in,n)
        IF (tc == 0)      CYCLE ! This triangle neighbour doesn't exist
        IF (tc == t_prev) CYCLE ! This triangle neighbour is the one we just came from
        qq = mesh%V( mesh%Tri( tc,1),:)
        rr = mesh%V( mesh%Tri( tc,2),:)
        ss = mesh%V( mesh%Tri( tc,3),:)
        gctc = geometric_center( qq,rr,ss)
        dc = NORM2( gctc - p)
        IF (dc < dcmin) THEN
          dcmin = dc
          tcmin = tc
        END IF
      END DO

      IF (dcmin < d) THEN
        t_prev = ti_in
        ti_in = tcmin
      ELSE
        EXIT
      END IF

    END DO ! DO WHILE (ncycle < mesh%nTri)

    ! Check if the result from the linear search is correct.
    qq = mesh%V(mesh%Tri( ti_in,1),:)
    rr = mesh%V(mesh%Tri( ti_in,2),:)
    ss = mesh%V(mesh%Tri( ti_in,3),:)
    IF (is_in_triangle( qq, rr, ss, p)) RETURN
    IF (lies_on_line_segment( qq, rr, p, mesh%tol_dist)) RETURN
    IF (lies_on_line_segment( rr, ss, p, mesh%tol_dist)) RETURN
    IF (lies_on_line_segment( ss, qq, p, mesh%tol_dist)) RETURN

    ! It's not. Perform a flood-fill style outward search.

  ! == Flood-fill search ==
  ! =======================

    ALLOCATE( map(    mesh%nTri), source = 0)
    ALLOCATE( stack1( mesh%nTri), source = 0)
    ALLOCATE( stack2( mesh%nTri), source = 0)

    map( ti_in)  = 1 ! We checked that one.

    ! Add ti_in's neighbours to the stack.
    stackN1 = 0
    DO n = 1, 3
      IF (mesh%TriC( ti_in,n) > 0) THEN
        stackN1 = stackN1 + 1
        stack1( stackN1) = mesh%TriC( ti_in,n)
      END IF
    END DO

    FoundIt = .FALSE.
    DO WHILE (.NOT. FoundIt)
      ! Check all triangles in the stack. If they're not it, add their
      ! non-checked neighbours to the new stack.

      stack2  = 0
      stackN2 = 0

      DO n = 1, stackN1
        ti = stack1( n)
        qq = mesh%V( mesh%Tri( ti,1),:)
        rr = mesh%V( mesh%Tri( ti,2),:)
        ss = mesh%V( mesh%Tri( ti,3),:)
        IF (is_in_triangle( qq, rr, ss, p)) THEN

          ! Found it!
          FoundIt = .TRUE.
          ti_in = ti
          EXIT

        ELSE ! if (is_in_triangle(ti,p))
          ! Did not find it. And add this triangle's non-checked neighbours to the new stack.

          DO n2 = 1, 3
            tin = mesh%TriC( ti,n2)
            IF (tin == 0)       CYCLE ! This neighbour doesn't exist.
            IF (map( tin) == 1) CYCLE ! This neighbour has already been checked or is already in the stack.
            stackN2 = stackN2 + 1
            stack2( stackN2) = tin
            map( tin) = 1
          END DO

        END IF ! IF (is_in_triangle(q, r, s, p, tol)) THEN
      END DO ! DO n = 1, mesh%triStackN1

      ! Cycle stacks.
      stack1  = stack2
      stackN1 = stackN2

      ! If no more non-checked neighbours could be found, terminate and throw an error.
      IF (stackN2 == 0 .AND. .NOT. FoundIt) THEN
        CALL crash('find_containing_triangle - ERROR: couldnt find triangle containing this point! (p = [{dp_01}, {dp_02}])', dp_01 = p(1), dp_02 = p(2))
      END IF

    END DO ! DO WHILE (.NOT. FoundIt)

    ! Clean up after yourself
    DEALLOCATE( map   )
    DEALLOCATE( stack1)
    DEALLOCATE( stack2)

  END SUBROUTINE find_containing_triangle

  SUBROUTINE find_containing_vertex( mesh, p, vi)
    ! Find the vertex whose Voronoi cell contains the point p, using a "linear search"
    ! Start at initial guess vi. Check all neighbours of vi, find the one
    ! closest to p, select that one as the new vi. Repeat until all neighbours
    ! of vi are further away from p than vi itself.

    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: vi

    INTEGER                                       :: ncycle, vi_prev, ci, vc, vcmin
    REAL(dp)                                      :: d, dc, dcmin


    ncycle = 0
    vi_prev = vi
    DO WHILE (ncycle < mesh%nV)

      d = NORM2( mesh%V( vi,:) - p)

      dcmin = d + 10._dp
      vcmin = 0
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi,ci)
        IF (vc == vi_prev) CYCLE ! This is the neighbour we just came from
        dc = NORM2( mesh%V( vc,:) - p)
        IF (dc < dcmin) THEN
          dcmin = dc
          vcmin = vc
        END IF
      END DO

      IF (dcmin < d) THEN
        vi_prev = vi
        vi = vcmin
      ELSE
        RETURN
      END IF

    END DO ! DO WHILE (ncycle < mesh%nV)

    ! If we reach this point, we didnt find the containing vertex - should not be possible, so throw an error
    CALL crash('find_containing_vertex - ERROR: couldnt find closest vertex!')

  END SUBROUTINE find_containing_vertex

  PURE FUNCTION is_in_Voronoi_cell( mesh, p, vi) RESULT( isso)
    ! If the point p lies closer to vertex vi than to any other vertex, then
    ! by definition it lies inside the Voronoi cell of vertex vi.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(IN)          :: vi

    ! Local variables:
    LOGICAL                                       :: isso
    REAL(dp)                                      :: dist_vi
    INTEGER                                       :: vvi, vj

    isso = .TRUE.

    dist_vi = NORM2( mesh%V( vi,:) - p)

    DO vvi = 1, mesh%nC( vi)
      vj = mesh%C( vi,vvi)
      IF (NORM2( mesh%V( vj,:) - p) < dist_vi) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO

  END FUNCTION is_in_Voronoi_cell

! == Diagnostic tools

  SUBROUTINE write_mesh_to_screen( mesh)
    ! Write a mesh to the screen. Best not to do this with large meshes.

    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER                                       :: vi, ti

    WRITE(0,*) '============================================================================'
    WRITE(0,*) ''
    WRITE(0,*) ' name = ', mesh%name
    WRITE(0,*) ' xmin = ', mesh%xmin
    WRITE(0,*) ' xmax = ', mesh%xmax
    WRITE(0,*) ' ymin = ', mesh%ymin
    WRITE(0,*) ' ymax = ', mesh%ymax
    WRITE(0,*) ''
    WRITE(0,*) ' vi    nC             C            niTri         iTri         edge_index         x              y'
    DO vi = 1, mesh%nV
      WRITE(0,'(A,I3,A,I3,A,6I3,A,I3,A,6I3,A,I3,A,F12.1,A,F12.1)') &
      ' ', vi, '   ', mesh%nC(vi), '    ', mesh%C(vi,1:6), '    ', mesh%niTri(vi), '    ', mesh%iTri(vi,1:6), '    ', mesh%edge_index(vi), &
      '    ', mesh%V(vi,1), '    ', mesh%V(vi,2)
    END DO

    WRITE(0,*) ' ti       Tri         TriC'
    DO ti = 1, mesh%nTri
      WRITE(0,'(A,I3,A,3I3,A,3I3)') &
      ' ', ti, '   ', mesh%Tri(ti,:), '    ', mesh%TriC(ti,:)
    END DO
    WRITE(0,*) '============================================================================'

  END SUBROUTINE write_mesh_to_screen

  SUBROUTINE write_mesh_to_text_file( mesh, filename)
    ! Write a mesh to a text file.

    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    CHARACTER(LEN=*),         INTENT(IN)          :: filename
    INTEGER                                       :: vi, ci, ti, fp

    ! Create a new text file
    OPEN(newUNIT  = fp, FILE = TRIM( filename), STATUS = 'REPLACE')

    ! Header
    WRITE(UNIT = fp, FMT = '(A)')       ' Mesh data'
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  xmin    = ', mesh%xmin
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  xmax    = ', mesh%xmax
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  ymin    = ', mesh%ymin
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  ymax    = ', mesh%ymax
    WRITE(UNIT = fp, FMT = '(A,I2)')    '  nC_mem = ', mesh%nC_mem
    WRITE(UNIT = fp, FMT = '(A,I6)')    '  nV      = ', mesh%nV
    WRITE(UNIT = fp, FMT = '(A,I6)')    '  nTri    = ', mesh%nTri
    WRITE(UNIT = fp, FMT = '(A)')       ''
    WRITE(UNIT = fp, FMT = '(A,I6,A,I3,A)')       'Vertex data: ', mesh%nV, ' rows, ', 2 + 1 + mesh%nC_mem + 1 + mesh%nC_mem + 1, ' columns'
    WRITE(UNIT = fp, FMT = '(A)')       ''

    ! Vertex data
    WRITE(UNIT = fp, FMT = '(A)')       'V  nC  C  niTri  iTri  edge_index'
    DO vi = 1, mesh%nV
      WRITE(UNIT = fp, FMT = '(2F24.14,I3)', ADVANCE = 'NO') mesh%V(vi,1), mesh%V(vi,2), mesh%nC(vi)
      DO ci = 1, mesh%nC_mem
        WRITE(UNIT = fp, FMT = '(I6)', ADVANCE = 'NO') mesh%C(vi,ci)
      END DO
      WRITE(UNIT = fp, FMT = '(I3)', ADVANCE = 'NO') mesh%niTri(vi)
      DO ci = 1, mesh%nC_mem
        WRITE(UNIT = fp, FMT = '(I6)', ADVANCE = 'NO') mesh%iTri(vi,ci)
      END DO
      WRITE(UNIT = fp, FMT = '(I3)', ADVANCE = 'NO') mesh%edge_index(vi)
      WRITE(UNIT = fp, FMT = '(A)') ''
    END DO
    WRITE(UNIT = fp, FMT = '(A)')       ''

    ! Triangle data
    WRITE(UNIT = fp, FMT = '(A)')       'Tri  TriC'
    DO ti = 1, mesh%nTri
      WRITE(UNIT = fp, FMT = '(6I6)') mesh%Tri(ti,1), mesh%Tri(ti,2), mesh%Tri(ti,3), mesh%TriC(ti,1), mesh%TriC(ti,2), mesh%TriC(ti,3)
    END DO

    ! Close the text file
    CLOSE(UNIT = fp)

  END SUBROUTINE write_mesh_to_text_file

  SUBROUTINE check_mesh( mesh)
    ! Check if the mesh data is self-consistent

    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER                                       :: vi, ci, vc, ci2, vc2, iti, iti2, ti, n, v1, v2, v3, ti2, n2
    LOGICAL                                       :: FoundIt

    ! == V
    ! =============================================================
    DO vi = 1, mesh%nV
      IF (mesh%V(vi,1) < mesh%xmin - mesh%tol_dist .OR. mesh%V(vi,1) > mesh%xmax + mesh%tol_dist .OR. &
          mesh%V(vi,2) < mesh%ymin - mesh%tol_dist .OR. mesh%V(vi,2) > mesh%ymax + mesh%tol_dist) THEN
        WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' outside mesh domain! (x = [', &
          mesh%xmin, ', ', mesh%V(vi,1), ',', mesh%xmax, '], y = [', mesh%ymin, ', ', mesh%V(vi,2), ',', mesh%ymax, ']'
      END IF
    END DO

    ! == nC
    ! =============================================================
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        IF (mesh%C(vi,ci) == 0) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has fewer connections than nC says!'
      END DO
      DO ci = mesh%nC(vi)+1, mesh%nC_mem
        IF (mesh%C(vi,ci) > 0)  WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has more connections than nC says!'
      END DO
    END DO

    ! == C
    ! =============================================================
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        vc = mesh%C(vi,ci)
        FoundIt = .FALSE.
        DO ci2 = 1, mesh%nC(vc)
          vc2 = mesh%C(vc,ci2)
          IF (vc2==vi) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is connected to ', vc, ', but not the other way round!'
      END DO
    END DO

    ! == niTri
    ! =============================================================
    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri(vi)
        IF (mesh%iTri(vi,iti) == 0) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has fewer iTriangles than niTri says!'
      END DO
      DO iti = mesh%niTri(vi)+1, mesh%nC_mem
        IF (mesh%iTri(vi,iti) > 0)  WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has more iTriangles than nC says!'
      END DO
    END DO

    ! == iTri
    ! =============================================================
    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri(vi)
        ti = mesh%iTri(vi,iti)
        FoundIt = .FALSE.
        DO n = 1, 3
          IF (mesh%Tri(ti,n)==vi) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' lists triangle ', ti, ' in iTri, but that triangle ', ti, ' doesnt contain vertex ', vi, '!'
      END DO

      IF (mesh%edge_index(vi) == 0) THEN

        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          n = 0
          DO iti = 1, mesh%niTri(vi)
            ti = mesh%iTri(vi,iti)
            DO iti2 = 1, mesh%niTri(vc)
              ti2 = mesh%iTri(vc,iti2)
              IF (ti==ti2) THEN
                n = n+1
                EXIT
              END IF
            END DO
          END DO
          IF (.NOT. (n==2)) WRITE(0,*) ' check_mesh - ERROR: non-edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'
        END DO

      ELSE ! IF (mesh%edge_index(vi) == 0) THEN

        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          IF (mesh%edge_index(vc)==0) THEN

            n = 0
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO iti2 = 1, mesh%niTri(vc)
                ti2 = mesh%iTri(vc,iti2)
                IF (ti==ti2) THEN
                  n = n+1
                  EXIT
                END IF
              END DO
            END DO
            IF (.NOT. (n==2)) WRITE(0,*) ' check_mesh - ERROR: edge vertex ', vi, ' and non-edge vertex ', vc, ' share ', n, ' triangles'

          ELSE ! IF (mesh%edge_index(vc)==0) THEN

            n = 0
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO iti2 = 1, mesh%niTri(vc)
                ti2 = mesh%iTri(vc,iti2)
                IF (ti==ti2) THEN
                  n = n+1
                  EXIT
                END IF
              END DO
            END DO
            IF (.NOT. is_boundary_segment( mesh, vi, vc)) CYCLE
            IF (.NOT. (n==1)) WRITE(0,*) ' check_mesh - ERROR: edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'

          END IF
        END DO

      END IF
    END DO

    ! == edge_index
    ! =============================================================
    DO vi = 1, mesh%nV

      IF (mesh%edge_index(vi) == 0) THEN

        IF (mesh%V(vi,1) <= mesh%xmin .OR. mesh%V(vi,1) >= mesh%xmax .OR. mesh%V(vi,2) <= mesh%ymin .OR. mesh%V(vi,2) >= mesh%ymax) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 0 but lies on or beyond the mesh domain boundary!'
        END IF

        ! First and last neighbours must be connected
        vc  = mesh%C(vi,1)
        vc2 = mesh%C(vi,mesh%nC(vi))

        FoundIt = .FALSE.
        DO ci = 1, mesh%nC(vc)
          IF (mesh%C(vc,ci)==vc2) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 0, but its first and last neighbours are not connected!'

      ELSEIF (mesh%edge_index(vi) == 1) THEN

        IF (ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but does not lie on the N boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==8 .OR. mesh%edge_index(vc)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==1 .OR. mesh%edge_index(vc)==2)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      ELSEIF (mesh%edge_index(vi) == 2) THEN

        IF (.NOT. vi==3) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as NE corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but does not lie on the NE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==8 .OR. mesh%edge_index(vc)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==3 .OR. mesh%edge_index(vc)==4)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      ELSEIF (mesh%edge_index(vi) == 3) THEN

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but does not lie on the E boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==2 .OR. mesh%edge_index(vc)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==3 .OR. mesh%edge_index(vc)==4)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      ELSEIF (mesh%edge_index(vi) == 4) THEN

        IF (.NOT. vi==2) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as SE corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but does not lie on the SE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==2 .OR. mesh%edge_index(vc)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==5 .OR. mesh%edge_index(vc)==6)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      ELSEIF (mesh%edge_index(vi) == 5) THEN

        IF (ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but does not lie on the S boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==4 .OR. mesh%edge_index(vc)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==5 .OR. mesh%edge_index(vc)==6)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      ELSEIF (mesh%edge_index(vi) == 6) THEN

        IF (.NOT. vi==1) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as SW corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but does not lie on the SW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==4 .OR. mesh%edge_index(vc)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==7 .OR. mesh%edge_index(vc)==8)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      ELSEIF (mesh%edge_index(vi) == 7) THEN

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but does not lie on the W boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==6 .OR. mesh%edge_index(vc)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==7 .OR. mesh%edge_index(vc)==8)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      ELSEIF (mesh%edge_index(vi) == 8) THEN

        IF (.NOT. vi==4) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as NW corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but does not lie on the NW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==6 .OR. mesh%edge_index(vc)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==1 .OR. mesh%edge_index(vc)==2)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its last connection doesnt have a matching edge_index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its first iTri doesnt have a matching Tri_edge_index!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its last iTri doesnt have a matching Tri_edge_index!'
!        END IF

      END IF

    END DO

    ! == Tri
    ! =============================================================

    DO ti = 1, mesh%nTri
      DO n = 1, 3
        vi = mesh%Tri(ti,n)
        FoundIt = .FALSE.
        DO iti = 1, mesh%niTri(vi)
          IF (mesh%iTri(vi,iti) == ti) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains vertex ', vi, ', but that vertex doesnt list ti as an iTri!'
      END DO

      v1 = mesh%Tri(ti,1)
      v2 = mesh%Tri(ti,2)
      v3 = mesh%Tri(ti,3)

      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v1)
        vc = mesh%C(v1,ci)
        IF (vc==v2) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v1, ' and ', v2, '!'

      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v1)
        vc = mesh%C(v1,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v1, ' and ', v3, '!'

      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v2)
        vc = mesh%C(v2,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v2, ' and ', v3, '!'
    END DO

    ! == TriC
    ! =============================================================

    DO ti = 1, mesh%nTri
      DO n = 1, 3
        ti2 = mesh%TriC(ti,n)
        IF (ti2 == 0) THEN
!          IF (mesh%Tri_edge_index(ti) == 0) WRITE(0,*) ' check_mesh - ERROR: non-edge triangle ', ti, ' misses a neighbour!'
          CYCLE
        END IF
        FoundIt = .FALSE.
        DO n2 = 1, 3
          IF (mesh%TriC(ti2,n2) == ti) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' is connected to ', ti2, ', but not the other way round!'
      END DO
    END DO


  END SUBROUTINE check_mesh

END MODULE mesh_utilities
