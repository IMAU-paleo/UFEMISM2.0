MODULE mesh_utilities

  ! Generally useful functions used in mesh creation and updating.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh
  USE math_utilities                                         , ONLY: geometric_center, is_in_triangle, lies_on_line_segment, circumcenter, &
                                                                     line_from_points, line_line_intersection, encroaches_upon

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Finding the vertices of a vertex' Voronoi cell

  SUBROUTINE calc_Voronoi_cell_vertices(        mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a vertex's Voronoi cell
    !
    ! NOTE: points are not repeated; if you want to calculated a loop integral
    !       around the entire Voronoi cell, be sure to include the section
    !       from Vor( nVor,:) to Vor( 1,:)!
    !
    ! NOTE: assumes that the mesh does not have any triangles whose circumcenter
    !       lies outside of the mesh domain!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables:
    INTEGER                                            :: vvi

    IF (mesh%VBI( vi) == 0) THEN
      ! Free vertex

      CALL calc_Voronoi_cell_vertices_free( mesh, vi, Vor, nVor)

    ELSEIF (mesh%VBI( vi) == 2 .OR. &
            mesh%VBI( vi) == 4 .OR. &
            mesh%VBI( vi) == 6 .OR. &
            mesh%VBI( vi) == 8) THEN
      ! Corner vertex

      CALL calc_Voronoi_cell_vertices_corner( mesh, vi, Vor, nVor)

    ELSE
      ! Vorder vertex

      CALL calc_Voronoi_cell_vertices_border( mesh, vi, Vor, nVor)

    END IF

    ! Safety: sometimes, Voronoi vertices end up just very slightly
    !         outside the mesh domain; move them to the border.
    DO vvi = 1, nVor

      ! Safety: if a Voronoi vertex is too far outside the mesh domain, crash.
      IF (Vor( vvi,1) < mesh%xmin - mesh%tol_dist .OR. &
          Vor( vvi,1) > mesh%xmax + mesh%tol_dist .OR. &
          Vor( vvi,2) < mesh%ymin - mesh%tol_dist .OR. &
          Vor( vvi,2) > mesh%ymax + mesh%tol_dist) THEN
        CALL warning('mesh domain = [{dp_01} - {dp_02}, {dp_03} - {dp_04}]', dp_01 = mesh%xmin, dp_02 = mesh%xmax, dp_03 = mesh%ymin, dp_04 = mesh%ymax)
        CALL warning('Vor = [{dp_01}, {dp_02}]', dp_01 = Vor( vvi,1), dp_02 = Vor( vvi,2))
!        CALL crash('find_Voronoi_cell_vertices: found Voronoi cell vertex outside of mesh domain!')
      END IF

      Vor( vvi,1) = MAX( MIN( Vor( vvi,1), mesh%xmax), mesh%xmin)
      Vor( vvi,2) = MAX( MIN( Vor( vvi,2), mesh%ymax), mesh%ymin)

    END DO

  END SUBROUTINE calc_Voronoi_cell_vertices

  SUBROUTINE calc_Voronoi_cell_vertices_free(   mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a free vertex's Voronoi cell

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables:
    INTEGER                                            :: iti, ti

    ! Initialise
    Vor  = 0._dp
    nVor = 0

    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi, iti)
      nVor = nVor + 1
      Vor( nVor,:) = mesh%Tricc( ti,:)
    END DO

  END SUBROUTINE calc_Voronoi_cell_vertices_free

  SUBROUTINE calc_Voronoi_cell_vertices_border(   mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a border vertex's Voronoi cell

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables:
    INTEGER                                            :: iti, ti

    ! Initialise
    Vor  = 0._dp
    nVor = 0

    ! Start with the projection of the first triangle's circumcenter on the domain border
    ti = mesh%iTri( vi,1)
    IF     (mesh%VBI( vi) == 1) THEN
      ! This vertex lies on the northern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymax]
    ELSEIF (mesh%VBI( vi) == 3) THEN
      ! This vertex lies on the eastern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmax, mesh%Tricc( ti,2)]
    ELSEIF (mesh%VBI( vi) == 5) THEN
      ! This vertex lies on the southern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymin]
    ELSEIF (mesh%VBI( vi) == 7) THEN
      ! This vertex lies on the western border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmin, mesh%Tricc( ti,2)]
    ELSE
      CALL crash('vertex does not lie on the domain border!')
    END IF

    ! Then add all the triangle circumcenters
    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi, iti)
      nVor = nVor + 1
      Vor( nVor,:) = mesh%Tricc( ti,:)
    END DO

    ! End with the projection of the last triangle's circumcenter on the domain border
    ti = mesh%iTri( vi, mesh%niTri( vi))
    IF     (mesh%VBI( vi) == 1) THEN
      ! This vertex lies on the northern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymax]
    ELSEIF (mesh%VBI( vi) == 3) THEN
      ! This vertex lies on the eastern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmax, mesh%Tricc( ti,2)]
    ELSEIF (mesh%VBI( vi) == 5) THEN
      ! This vertex lies on the southern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymin]
    ELSEIF (mesh%VBI( vi) == 7) THEN
      ! This vertex lies on the western border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmin, mesh%Tricc( ti,2)]
    ELSE
      CALL crash('vertex does not lie on the domain border!')
    END IF

  END SUBROUTINE calc_Voronoi_cell_vertices_border

  SUBROUTINE calc_Voronoi_cell_vertices_corner( mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a corner vertex's Voronoi cell

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: Vor
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables:
    INTEGER                                            :: iti, ti

    ! Initialise
    Vor  = 0._dp
    nVor = 0

    ! Start with the projection of the first triangle's circumcenter on the domain border
    ti = mesh%iTri( vi,1)
    IF     (mesh%VBI( vi) == 2) THEN
      ! This vertex lies on the northeast corner; project onto the northern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymax]
    ELSEIF (mesh%VBI( vi) == 4) THEN
      ! This vertex lies on the southeast corner; project onto the eastern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmax, mesh%Tricc( ti,2)]
    ELSEIF (mesh%VBI( vi) == 6) THEN
      ! This vertex lies on the southwest corner; project onto the southern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymin]
    ELSEIF (mesh%VBI( vi) == 8) THEN
      ! This vertex lies on the northwest corner; project onto the western border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmin, mesh%Tricc( ti,2)]
    ELSE
      CALL crash('vertex does not lie on the domain border!')
    END IF

    ! Then add all the triangle circumcenters
    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi, iti)
      nVor = nVor + 1
      Vor( nVor,:) = mesh%Tricc( ti,:)
    END DO

    ! End with the projection of the last triangle's circumcenter on the domain border
    ti = mesh%iTri( vi,1)
    IF     (mesh%VBI( vi) == 2) THEN
      ! This vertex lies on the northeast corner; project onto the eastern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmax, mesh%Tricc( ti,2)]
    ELSEIF (mesh%VBI( vi) == 4) THEN
      ! This vertex lies on the southeast corner; project onto the southern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymin]
    ELSEIF (mesh%VBI( vi) == 6) THEN
      ! This vertex lies on the southwest corner; project onto the western border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmin, mesh%Tricc( ti,2)]
    ELSEIF (mesh%VBI( vi) == 8) THEN
      ! This vertex lies on the northwest corner; project onto the northern border
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%Tricc( ti,1), mesh%ymax]
    ELSE
      CALL crash('vertex does not lie on the domain border!')
    END IF

    ! Finally, include the vertex itself (i.e. the domain corner)
    nVor = nVor + 1
    Vor( nVor,:) = mesh%V( vi,:)

  END SUBROUTINE calc_Voronoi_cell_vertices_corner

! == Some basic geometrical operations

  SUBROUTINE list_border_vertices_all( mesh, nvi_border, lvi_border)
    ! List all border vertices of the mesh, ordered counter-clockwise,
    ! starting from vertex 1 in the southwest corner.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it

    ! Initialise
    lvi_border = 0
    nvi_border = 0

    ! Trace the border
    vi = 1
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! Move to next vertex along the border
      vi = mesh%C( vi,1)

      ! If we've reached the southwest corner again, stop
      IF (vi == 1) EXIT

    END DO

  END SUBROUTINE list_border_vertices_all

  SUBROUTINE list_border_vertices_south( mesh, nvi_border, lvi_border)
    ! List all southern border vertices of the mesh, ordered east to west.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! South
    vi_start = 2
    vi_stop  = 1

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi, mesh%nC( vi))

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_south

  SUBROUTINE list_border_vertices_east( mesh, nvi_border, lvi_border)
    ! List all eastern border vertices of the mesh, ordered south to north.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! East
    vi_start = 2
    vi_stop  = 3

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi,1)

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_east

  SUBROUTINE list_border_vertices_north( mesh, nvi_border, lvi_border)
    ! List all northern border vertices of the mesh, ordered east to west.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! North
    vi_start = 3
    vi_stop  = 4

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi,1)

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_north

  SUBROUTINE list_border_vertices_west( mesh, nvi_border, lvi_border)
    ! List all western border vertices of the mesh, ordered south to north.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! West
    vi_start = 1
    vi_stop  = 4

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi, mesh%nC( vi))

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_west

  SUBROUTINE encroaches_upon_any( mesh, p, isso, vi_encroached, vj_encroached)
    ! Check if p encroaches upon a border edge.
    ! If so, return the indices [vi,vj] of the vertices spanning that edge.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p
    LOGICAL,                    INTENT(OUT)       :: isso
    INTEGER,                    INTENT(OUT)       :: vi_encroached, vj_encroached

    ! Local variables:
    INTEGER                                       :: nvi_border
    INTEGER,  DIMENSION(mesh%nV)                  :: lvi_border
    INTEGER                                       :: i, j, vi, vj
    REAL(dp), DIMENSION(2)                        :: pa, pb

    isso          = .FALSE.
    vi_encroached = 0
    vj_encroached = 0

    ! List border vertices
    CALL list_border_vertices_all( mesh, nvi_border, lvi_border)

    ! Check if p encroaches on any of the border edges
    DO i = 1, nvi_border

      j = i + 1
      IF (i == nvi_border) j = 1

      vi = lvi_border( i)
      vj = lvi_border( j)

      pa = mesh%V( vi,:)
      pb = mesh%V( vj,:)

      IF (encroaches_upon( pa, pb, p, mesh%tol_dist)) THEN
        isso          = .TRUE.
        vi_encroached = vi
        vj_encroached = vj
        EXIT
      END IF

    END DO ! DO i = 1, nvi_border

  END SUBROUTINE encroaches_upon_any

  PURE FUNCTION is_border_edge( mesh, vi, vj) RESULT(isso)
    ! Determine whether or not the edge between two vertices is a border edge

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: vi, vj

    ! Local variables:
    LOGICAL                                       :: isso

    IF (mesh%VBI( vi) == 0 .OR. mesh%VBI( vj) == 0) THEN
      isso = .FALSE.
      RETURN
    END IF

    isso = .FALSE.

    IF (mesh%VBI( vi) == 1) THEN
      IF (mesh%VBI( vj) == 8 .OR. &
          mesh%VBI( vj) == 1 .OR. &
          mesh%VBI( vj) == 2) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 2) THEN
      IF (mesh%VBI( vj) == 8 .OR. &
          mesh%VBI( vj) == 1 .OR. &
          mesh%VBI( vj) == 3 .OR. &
          mesh%VBI( vj) == 4) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 3) THEN
      IF (mesh%VBI( vj) == 2 .OR. &
          mesh%VBI( vj) == 3 .OR. &
          mesh%VBI( vj) == 4) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 4) THEN
      IF (mesh%VBI( vj) == 2 .OR. &
          mesh%VBI( vj) == 3 .OR. &
          mesh%VBI( vj) == 5 .OR. &
          mesh%VBI( vj) == 6) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 5) THEN
      IF (mesh%VBI( vj) == 4 .OR. &
          mesh%VBI( vj) == 5 .OR. &
          mesh%VBI( vj) == 6) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 6) THEN
      IF (mesh%VBI( vj) == 4 .OR. &
          mesh%VBI( vj) == 5 .OR. &
          mesh%VBI( vj) == 7 .OR. &
          mesh%VBI( vj) == 8) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 7) THEN
      IF (mesh%VBI( vj) == 6 .OR. &
          mesh%VBI( vj) == 7 .OR. &
          mesh%VBI( vj) == 8) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 8) THEN
      IF (mesh%VBI( vj) == 6 .OR. &
          mesh%VBI( vj) == 7 .OR. &
          mesh%VBI( vj) == 1 .OR. &
          mesh%VBI( vj) == 2) THEN
        isso = .TRUE.
      END IF
    END IF

  END FUNCTION is_border_edge

  PURE FUNCTION is_walltowall( mesh, ti) RESULT(isso)
   ! Determine whether or not a triangle is "wall to wall"
   ! (i.e. contains vertices lying on opposite domain borders)

    IMPLICIT NONE

    ! In/output variables:
   TYPE(type_mesh),          INTENT(IN)          :: mesh
   INTEGER,                  INTENT(IN)          :: ti
   LOGICAL                                       :: isso

    ! Local variables:
   LOGICAL                                       :: has_north, has_south, has_east, has_west
   INTEGER                                       :: n, vi

   has_north = .FALSE.
   has_south = .FALSE.
   has_east  = .FALSE.
   has_west  = .FALSE.

   DO n = 1, 3
     vi = mesh%Tri( ti,n)
     IF (mesh%VBI( vi) == 1) has_north = .TRUE.
     IF (mesh%VBI( vi) == 3) has_east  = .TRUE.
     IF (mesh%VBI( vi) == 5) has_south = .TRUE.
     IF (mesh%VBI( vi) == 7) has_west  = .TRUE.
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

! == Tools for handling the triangle refinement stack

  SUBROUTINE add_triangle_to_refinement_stack_first( mesh, ti)
    ! Add triangle ti to the top of the refinement stack

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Safety
    IF (ti == 0) CALL crash('received ti=0!')

    IF (mesh%refinement_map( ti) == 0) THEN
      mesh%refinement_map( ti) = 1
      mesh%refinement_stackN = mesh%refinement_stackN + 1
      mesh%refinement_stack( 2: mesh%nTri_mem) = mesh%refinement_stack( 1: mesh%nTri_mem-1)
      mesh%refinement_stack( 1) = ti
    END IF

  END SUBROUTINE add_triangle_to_refinement_stack_first

  SUBROUTINE add_triangle_to_refinement_stack_last( mesh, ti)
    ! Add triangle ti to the end of the refinement stack

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Safety
    IF (ti == 0) CALL crash('received ti=0!')

    IF (mesh%refinement_map( ti) == 0) THEN
      mesh%refinement_map( ti) = 1
      mesh%refinement_stackN = mesh%refinement_stackN + 1
      mesh%refinement_stack( mesh%refinement_stackN) = ti
    END IF

  END SUBROUTINE add_triangle_to_refinement_stack_last

  SUBROUTINE remove_triangle_from_refinement_stack( mesh, ti)
    ! Remove triangle ti from the refinement stack

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Local variables:
    INTEGER                                       :: i

    ! Safety
    IF (ti == 0) CALL crash('received ti=0!')

    mesh%refinement_map( ti) = 0

    i = 1
    DO WHILE (i <= mesh%refinement_stackN)

      IF (mesh%refinement_stack( i) == ti) THEN
        mesh%refinement_stack( 1:mesh%refinement_stackN) = [mesh%refinement_stack( 1:i-1), mesh%refinement_stack( i+1:mesh%refinement_stackN), 0]
        mesh%refinement_stackN = mesh%refinement_stackN - 1
      ELSE
        i = i+1
      END IF

    END DO

  END SUBROUTINE remove_triangle_from_refinement_stack

! == Some basic search operations on a mesh

  SUBROUTINE find_containing_triangle( mesh, p, ti_in)
    ! Find the triangle containing the point p. First do a "linear search":
    ! Start at initial guess ti_in. Check all neighbours of ti_in, find the one
    ! closest to p, select that one as the new ti_in. Repeat until all neighbours
    ! of ti_in are further away from p than ti_in itself.
    ! Then (if needed) search outward from there using a flood-fill algorithm.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: ti_in

    ! Local variables:
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
      CALL crash('p lies outside mesh domain!')
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
        CALL crash('couldnt find triangle containing p = [{dp_01}, {dp_02}]', dp_01 = p(1), dp_02 = p(2))
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

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: vi

    ! Local variables:
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
    CALL crash('couldnt find closest vertex!')

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

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh

    ! Local variables:
    INTEGER                                       :: vi, ti

    WRITE(0,*) '============================================================================'
    WRITE(0,*) ''
    WRITE(0,*) ' name = ', mesh%name
    WRITE(0,*) ' xmin = ', mesh%xmin
    WRITE(0,*) ' xmax = ', mesh%xmax
    WRITE(0,*) ' ymin = ', mesh%ymin
    WRITE(0,*) ' ymax = ', mesh%ymax
    WRITE(0,*) ''
    WRITE(0,*) ' vi    nC             C            niTri         iTri         VBI                x              y'
    DO vi = 1, mesh%nV
      WRITE(0,'(A,I3,A,I3,A,6I3,A,I3,A,6I3,A,I3,A,F12.1,A,F12.1)') &
      ' ', vi, '   ', mesh%nC(vi), '    ', mesh%C(vi,1:6), '    ', mesh%niTri(vi), '    ', mesh%iTri(vi,1:6), '    ', mesh%VBI(vi), &
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

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    CHARACTER(LEN=*),         INTENT(IN)          :: filename

    ! Local variables:
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
    WRITE(UNIT = fp, FMT = '(A)')       'V  nC  C  niTri  iTri  VBI'
    DO vi = 1, mesh%nV
      WRITE(UNIT = fp, FMT = '(2F24.14,I3)', ADVANCE = 'NO') mesh%V(vi,1), mesh%V(vi,2), mesh%nC(vi)
      DO ci = 1, mesh%nC_mem
        WRITE(UNIT = fp, FMT = '(I6)', ADVANCE = 'NO') mesh%C(vi,ci)
      END DO
      WRITE(UNIT = fp, FMT = '(I3)', ADVANCE = 'NO') mesh%niTri(vi)
      DO ci = 1, mesh%nC_mem
        WRITE(UNIT = fp, FMT = '(I6)', ADVANCE = 'NO') mesh%iTri(vi,ci)
      END DO
      WRITE(UNIT = fp, FMT = '(I3)', ADVANCE = 'NO') mesh%VBI(vi)
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

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh

    ! Local variables:
    INTEGER                                       :: vi, ci, vc, ci2, vc2, iti, iti2, ti, n, v1, v2, v3, ti2, n2
    LOGICAL                                       :: FoundIt
    INTEGER                                       :: tj, tivia, tivib, tivic, tjvia, tjvib, tjvic

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

      IF (mesh%VBI(vi) == 0) THEN

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

      ELSE ! IF (mesh%VBI(vi) == 0) THEN

        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          IF (mesh%VBI(vc)==0) THEN

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

          ELSE ! IF (mesh%VBI(vc)==0) THEN

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
            IF (.NOT. is_border_edge( mesh, vi, vc)) CYCLE
            IF (.NOT. (n==1)) WRITE(0,*) ' check_mesh - ERROR: edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'

          END IF
        END DO

      END IF
    END DO

    ! == border index
    ! =============================================================
    DO vi = 1, mesh%nV

      IF (mesh%VBI(vi) == 0) THEN

        IF (mesh%V(vi,1) <= mesh%xmin .OR. mesh%V(vi,1) >= mesh%xmax .OR. mesh%V(vi,2) <= mesh%ymin .OR. mesh%V(vi,2) >= mesh%ymax) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 0 but lies on or beyond the mesh domain border!'
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
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 0, but its first and last neighbours are not connected!'

      ELSEIF (mesh%VBI(vi) == 1) THEN

        IF (ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 1 but does not lie on the N border!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==8 .OR. mesh%VBI(vc)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 1 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==1 .OR. mesh%VBI(vc)==2)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 1 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==7 .OR. mesh%TriBI(ti)==8 .OR. mesh%TriBI(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 1 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==1 .OR. mesh%TriBI(ti)==2 .OR. mesh%TriBI(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 1 but its last iTri doesnt have a matching TriBI!'
!        END IF

      ELSEIF (mesh%VBI(vi) == 2) THEN

        IF (.NOT. vi==3) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as NE corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 2 but does not lie on the NE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==8 .OR. mesh%VBI(vc)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 2 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==3 .OR. mesh%VBI(vc)==4)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 2 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==1 .OR. mesh%TriBI(ti)==2 .OR. mesh%TriBI(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 2 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==1 .OR. mesh%TriBI(ti)==2 .OR. mesh%TriBI(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 2 but its last iTri doesnt have a matching TriBI!'
!        END IF

      ELSEIF (mesh%VBI(vi) == 3) THEN

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 3 but does not lie on the E border!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==2 .OR. mesh%VBI(vc)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 3 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==3 .OR. mesh%VBI(vc)==4)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 3 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==1 .OR. mesh%TriBI(ti)==2 .OR. mesh%TriBI(ti)==3)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 3 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==3 .OR. mesh%TriBI(ti)==4 .OR. mesh%TriBI(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 3 but its last iTri doesnt have a matching TriBI!'
!        END IF

      ELSEIF (mesh%VBI(vi) == 4) THEN

        IF (.NOT. vi==2) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as SE corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 4 but does not lie on the SE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==2 .OR. mesh%VBI(vc)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 4 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==5 .OR. mesh%VBI(vc)==6)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 4 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==3 .OR. mesh%TriBI(ti)==4 .OR. mesh%TriBI(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 4 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==3 .OR. mesh%TriBI(ti)==4 .OR. mesh%TriBI(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 4 but its last iTri doesnt have a matching TriBI!'
!        END IF

      ELSEIF (mesh%VBI(vi) == 5) THEN

        IF (ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 5 but does not lie on the S border!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==4 .OR. mesh%VBI(vc)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 5 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==5 .OR. mesh%VBI(vc)==6)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 5 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==3 .OR. mesh%TriBI(ti)==4 .OR. mesh%TriBI(ti)==5)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 5 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==5 .OR. mesh%TriBI(ti)==6 .OR. mesh%TriBI(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 5 but its last iTri doesnt have a matching TriBI!'
!        END IF

      ELSEIF (mesh%VBI(vi) == 6) THEN

        IF (.NOT. vi==1) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as SW corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 6 but does not lie on the SW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==4 .OR. mesh%VBI(vc)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 6 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==7 .OR. mesh%VBI(vc)==8)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 6 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==5 .OR. mesh%TriBI(ti)==6 .OR. mesh%TriBI(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 6 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==5 .OR. mesh%TriBI(ti)==6 .OR. mesh%TriBI(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 6 but its last iTri doesnt have a matching TriBI!'
!        END IF

      ELSEIF (mesh%VBI(vi) == 7) THEN

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 7 but does not lie on the W border!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==6 .OR. mesh%VBI(vc)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 7 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==7 .OR. mesh%VBI(vc)==8)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 7 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==5 .OR. mesh%TriBI(ti)==6 .OR. mesh%TriBI(ti)==7)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 7 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==7 .OR. mesh%TriBI(ti)==8 .OR. mesh%TriBI(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 7 but its last iTri doesnt have a matching TriBI!'
!        END IF

      ELSEIF (mesh%VBI(vi) == 8) THEN

        IF (.NOT. vi==4) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as NW corner!'

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 8 but does not lie on the NW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==6 .OR. mesh%VBI(vc)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 8 but its first connection doesnt have a matching border index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==1 .OR. mesh%VBI(vc)==2)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 8 but its last connection doesnt have a matching border index!'
        END IF
!        ti = mesh%iTri(vi,1)
!        IF (.NOT. (mesh%TriBI(ti)==7 .OR. mesh%TriBI(ti)==8 .OR. mesh%TriBI(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 8 but its first iTri doesnt have a matching TriBI!'
!        END IF
!        ti = mesh%iTri(vi,mesh%niTri(vi))
!        IF (.NOT. (mesh%TriBI(ti)==7 .OR. mesh%TriBI(ti)==8 .OR. mesh%TriBI(ti)==1)) THEN
!          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has border index 8 but its last iTri doesnt have a matching TriBI!'
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

      tivia = mesh%Tri( ti,1)
      tivib = mesh%Tri( ti,2)
      tivic = mesh%Tri( ti,3)
      DO tj = ti+1, mesh%nTri
        tjvia = mesh%Tri( tj,1)
        tjvib = mesh%Tri( tj,2)
        tjvic = mesh%Tri( tj,3)
        IF (ANY( [tjvia,tjvib,tjvic] == tivia) .AND. &
            ANY( [tjvia,tjvib,tjvic] == tivib) .AND. &
            ANY( [tjvia,tjvib,tjvic] == tivic)) WRITE(0,*) ' check_mesh - ERROR: triangles ', ti, ' and ', tj, ' are made up of the same vertices!'
      END DO
    END DO

    ! == TriC
    ! =============================================================

    DO ti = 1, mesh%nTri
      DO n = 1, 3
        ti2 = mesh%TriC(ti,n)
        IF (ti2 == 0) THEN
!          IF (mesh%TriBI(ti) == 0) WRITE(0,*) ' check_mesh - ERROR: non-edge triangle ', ti, ' misses a neighbour!'
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

  SUBROUTINE check_if_triangle_already_exists( mesh, ti)
    ! Check if any duplicate triangles exist

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: ti

    ! Local variables:
    INTEGER                                       :: via, vib, vic, tj

    via = mesh%Tri( ti,1)
    vib = mesh%Tri( ti,2)
    vic = mesh%Tri( ti,3)

    DO tj = 1, mesh%nTri
      IF (tj == ti) CYCLE
      IF (ANY( via == mesh%Tri( tj,:)) .AND. ANY( vib == mesh%Tri( tj,:)) .AND. ANY( vic == mesh%Tri( tj,:))) THEN
        CALL crash('duplicate triangles detected at ti = {int_01}, tj = {int_02}', int_01 = ti, int_02 = tj)
      END IF
    END DO

  END SUBROUTINE check_if_triangle_already_exists

END MODULE mesh_utilities
