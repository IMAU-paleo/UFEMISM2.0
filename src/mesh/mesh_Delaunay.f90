MODULE mesh_Delaunay

  ! Routines used for updating a Delaunay triangulation by splitting a triangle, a line or a segment.
  !
  ! All these routines are set up so that, if the mesh that goes in is self-consistent (i.e. there
  ! are no erroneous entries in the connectivity lists) and meets the Delaunay criterion everywhere,
  ! then so does the mesh that comes out.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh
  USE math_utilities                                         , ONLY: is_in_triangle, lies_on_line_segment, encroaches_upon
  USE mesh_utilities                                         , ONLY: update_triangle_circumcenter, find_containing_triangle, is_boundary_segment

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE split_triangle( mesh, ti_in_guess, p_new, refinement_map, refinement_stack, refinement_stackN, Tri_li)
     ! Add a new vertex at p_new, possibly by splitting triangle ti_in_guess

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: ti_in_guess
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p_new
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_map
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_stack
    INTEGER,                  OPTIONAL,  INTENT(INOUT) :: refinement_stackN
    INTEGER,  DIMENSION(:,:), OPTIONAL,  INTENT(INOUT) :: Tri_li

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'split_triangle'
    INTEGER                                            :: v1, v2, t1, n, nc, nnext, t_old
    INTEGER                                            :: p1, p2, p3, tn1, tn2, tn3, t_new1, t_new2, t_new3
    REAL(dp), DIMENSION(2)                             :: p
    LOGICAL                                            :: isencroached
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa,pb,pc
    INTEGER,  DIMENSION(:,:), ALLOCATABLE              :: Tri_flip
    INTEGER                                            :: nf
    LOGICAL                                            :: did_flip

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == If p_new lies outside of the mesh domain, split a boundary segment instead.
  ! ==============================================================================

    IF (p_new( 1) < mesh%xmin) THEN
      ! p_new lies to the west of the mesh domain

      ! Safety
      IF (p_new( 2) < mesh%ymin .OR. p_new( 2) > mesh%ymax) CALL crash('p_new lies way outside mesh domain!')

      ! Find the two vertices v1,v2 of the segment that must be split.
      v1 = 1
      v2 = mesh%C( v1, mesh%nC( v1))
      DO WHILE (mesh%V( v2,2) < p_new( 2))
        v1 = v2
        v2 = mesh%C( v1, mesh%nC( v1))
      END DO

      ! Safety
      IF (.NOT. (p_new( 2) >= mesh%V( v1,2) .AND. p_new( 2) <= mesh%V( v2,2))) CALL crash('couldnt find boundary segment to split!')

      p = (mesh%V( v1,:) + mesh%V( v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 1) < mesh%xmin) THEN

    IF (p_new( 1) > mesh%xmax) THEN
      ! p_new lies to the east of the mesh domain

      ! Safety
      IF (p_new( 2) < mesh%ymin .OR. p_new( 2) > mesh%ymax) CALL crash('p_new lies way outside mesh domain!')

      ! Find the two vertices v1,v2 of the segment that must be split.
      v1 = 2
      v2 = mesh%C( v1,1)
      DO WHILE (mesh%V( v2,2) < p_new( 2))
        v1 = v2
        v2 = mesh%C( v1, 1)
      END DO

      ! Safety
      IF (.NOT. (p_new( 2) >= mesh%V( v1,2) .AND. p_new( 2) <= mesh%V( v2,2))) CALL crash('couldnt find boundary segment to split!')

      p = (mesh%V( v1,:) + mesh%V( v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 1) < mesh%xmin) THEN

    IF (p_new( 2) < mesh%ymin) THEN
      ! p_new lies to the south of the mesh domain

      ! Safety
      IF (p_new( 1) < mesh%xmin .OR. p_new( 1) > mesh%xmax) CALL crash('p_new lies way outside mesh domain!')

      ! Find the two vertices v1,v2 of the segment that must be split.
      v1 = 1
      v2 = mesh%C( v1, 1)
      DO WHILE (mesh%V( v2,1) < p_new( 1))
        v1 = v2
        v2 = mesh%C( v1, 1)
      END DO

      ! Safety
      IF (.NOT. (p_new( 1) >= mesh%V( v1,1) .AND. p_new( 1) <= mesh%V( v2,1))) CALL crash('couldnt find segment to split!')

      p = (mesh%V( v1,:) + mesh%V( v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 1) < mesh%xmin) THEN

    IF (p_new( 2) > mesh%ymax) THEN
      ! p_new lies to the north of the mesh domain

      ! Safety
      IF (p_new( 1) < mesh%xmin .OR. p_new( 1) > mesh%xmax) CALL crash('p_new lies way outside mesh domain!')

      ! Find the two vertices v1,v2 of the segment that must be split.
      v1 = 1
      v2 = mesh%C( v1, mesh%nC( v1))
      DO WHILE (mesh%V( v2,1) < p_new( 1))
        v1 = v2
        v2 = mesh%C( v1, mesh%nC( v1))
      END DO

      ! Safety
      IF (.NOT. (p_new( 1) >= mesh%V( v1,1) .AND. p_new( 1) <= mesh%V( v2,1))) CALL crash('couldnt find segment to split!')

      p = (mesh%V( v1,:) + mesh%V( v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 1) < mesh%xmin) THEN

  ! == Find the triangle containing p_new.
  ! ======================================

    ! Find the triangle containing p_new
    t_old = ti_in_guess
    CALL find_containing_triangle( mesh, p_new, t_old)

    via = mesh%Tri( t_old,1)
    vib = mesh%Tri( t_old,2)
    vic = mesh%Tri( t_old,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! If the new vertex is (almost) colinear with two vertices of
    ! the containing triangle, split the line between them instead.

    IF     (lies_on_line_segment( pa, pb, p_new, mesh%tol_dist)) THEN
      CALL split_line( mesh, via, vib, p_new, refinement_map, refinement_stack, refinement_stackN, Tri_li)
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (lies_on_line_segment( pb, pc, p_new, mesh%tol_dist)) THEN
      CALL split_line( mesh, vib, vic, p_new, refinement_map, refinement_stack, refinement_stackN, Tri_li)
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (lies_on_line_segment( pc, pa, p_new, mesh%tol_dist)) THEN
      CALL split_line( mesh, vic, via, p_new, refinement_map, refinement_stack, refinement_stackN, Tri_li)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! If p_new encroaches upon a boundary segment, split that segment instead.

    IF     (is_boundary_segment( mesh, via, vib) .AND. encroaches_upon( pa, pb, p_new)) THEN
      p = (pa + pb) / 2._dp
      CALL split_segment( mesh, via, vib, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (is_boundary_segment( mesh, vib, vic) .AND. encroaches_upon( pa, pb, p_new)) THEN
      p = (pb + pc) / 2._dp
      CALL split_segment( mesh, vib, vic, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (is_boundary_segment( mesh, vic, via) .AND. encroaches_upon( pa, pb, p_new)) THEN
      p = (pc + pa) / 2._dp
      CALL split_segment( mesh, vic, via, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == Find that triangle's vertices and neighbours
    p1  = mesh%Tri(  t_old,1)
    p2  = mesh%Tri(  t_old,2)
    p3  = mesh%Tri(  t_old,3)
    tn1 = mesh%TriC( t_old,1)
    tn2 = mesh%TriC( t_old,2)
    tn3 = mesh%TriC( t_old,3)

    ! == Add the new vertex to the mesh
    mesh%nV = mesh%nV+1
    mesh%V( mesh%nV,:) = p_new
    mesh%edge_index( mesh%nV) = 0

    ! == Replace ti_old by three new triangles
    t_new1 = t_old
    t_new2 = mesh%nTri+1
    t_new3 = mesh%nTri+2
    mesh%Tri( t_new1,:) = [p1, p2, mesh%nV]
    mesh%Tri( t_new2,:) = [p2, p3, mesh%nV]
    mesh%Tri( t_new3,:) = [p3, p1, mesh%nV]
    mesh%nTri = mesh%nTri+2

    ! Add these to the refinement stack
    IF (PRESENT( refinement_map)) THEN
      CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, t_old)

      refinement_map( t_new1) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t_new1

      refinement_map( t_new2) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t_new2

      refinement_map( t_new2) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t_new2
    END IF ! IF (PRESENT( refinement_map)) THEN

    ! == Update triangle-line overlap ranges
    IF (PRESENT( Tri_li)) THEN
      Tri_li( t_new1,:) = Tri_li( t_old,:)
      Tri_li( t_new2,:) = Tri_li( t_old,:)
      Tri_li( t_new3,:) = Tri_li( t_old,:)
    END IF

    ! == Update vertex connectivity matrix
    ! p1
    DO n = 1, mesh%nC( p1)
      IF ( mesh%C( p1,n) == p2) THEN
        mesh%C( p1,:) = [mesh%C( p1,1:n), mesh%nV, mesh%C( p1,n+1:mesh%nC_mem-1)]
        mesh%nC( p1) = mesh%nC( p1)+1
        EXIT
      END IF
    END DO
    ! p2
    DO n = 1, mesh%nC( p2)
      IF ( mesh%C( p2,n) == p3) THEN
        mesh%C( p2,:) = [mesh%C( p2,1:n), mesh%nV, mesh%C( p2,n+1:mesh%nC_mem-1)]
        mesh%nC( p2) = mesh%nC( p2)+1
        EXIT
      END IF
    END DO
    ! p3
    DO n = 1, mesh%nC( p3)
      IF ( mesh%C( p3,n) == p1) THEN
        mesh%C( p3,:) = [mesh%C( p3,1:n), mesh%nV, mesh%C( p3,n+1:mesh%nC_mem-1)]
        mesh%nC( p3) = mesh%nC( p3)+1
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%C( mesh%nV,1:3) = [p1, p2, p3]
    mesh%nC( mesh%nV)    = 3

    ! == Update triangle connectivity matrix
    ! (existing) neighbours
    DO n = 1, 3
      IF ( tn1>0) THEN
        IF ( mesh%TriC( tn1,n) == t_old) mesh%TriC( tn1,n) = t_new2
      END IF
      IF ( tn2>0) THEN
        IF ( mesh%TriC( tn2,n) == t_old) mesh%TriC( tn2,n) = t_new3
      END IF
      IF ( tn3>0) THEN
        IF ( mesh%TriC( tn3,n) == t_old) mesh%TriC( tn3,n) = t_new1
      END IF
    END DO
    ! new triangles
    mesh%TriC( t_new1,:) = [t_new2, t_new3, tn3]
    mesh%TriC( t_new2,:) = [t_new3, t_new1, tn1]
    mesh%TriC( t_new3,:) = [t_new1, t_new2, tn2]

    ! == Update inverse triangle lists
    ! p1
    DO n = 1, mesh%niTri( p1)
      IF ( mesh%iTri( p1,n) == t_old) THEN
        mesh%iTri(  p1,:) = [mesh%iTri( p1,1:n-1), t_new1, t_new3, mesh%iTri( p1,n+1:mesh%nC_mem-1)]
        mesh%niTri( p1  ) = mesh%niTri( p1)+1
        EXIT
      END IF
    END DO
    ! p2
    DO n = 1, mesh%niTri( p2)
      IF ( mesh%iTri( p2,n) == t_old) THEN
        mesh%iTri(  p2,:) = [mesh%iTri( p2,1:n-1), t_new2, t_new1, mesh%iTri( p2,n+1:mesh%nC_mem-1)]
        mesh%niTri( p2  ) = mesh%niTri( p2)+1
        EXIT
      END IF
    END DO
    ! p3
    DO n = 1, mesh%niTri( p3)
      IF ( mesh%iTri( p3,n) == t_old) THEN
        mesh%iTri(  p3,:) = [mesh%iTri( p3,1:n-1), t_new3, t_new2, mesh%iTri( p3,n+1:mesh%nC_mem-1)]
        mesh%niTri( p3  ) = mesh%niTri( p3)+1
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%iTri(  mesh%nV,1:3) = [t_new1, t_new2, t_new3]
    mesh%niTri( mesh%nV    ) = 3

    ! == Update triangle circumcenters
    CALL update_triangle_circumcenter( mesh, t_new1)
    CALL update_triangle_circumcenter( mesh, t_new2)
    CALL update_triangle_circumcenter( mesh, t_new3)

    ! == Propagate flip operations outward
    ! Start with the newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    ALLOCATE( Tri_flip( mesh%nTri, 2), source = 0)
    nf = 0

    IF (tn1 > 0) THEN
      nf = nf + 1
      Tri_flip( nf,:) = [t_new2, tn1]
    END IF
    IF (tn2 > 0) THEN
      nf = nf + 1
      Tri_flip( nf,:) = [t_new3, tn2]
    END IF
    IF (tn3 > 0) THEN
      nf = nf + 1
      Tri_flip( nf,:) = [t_new1, tn3]
    END IF

    DO WHILE (nf > 0)
      CALL flip_triangle_pairs( mesh, Tri_flip, nf, did_flip, refinement_map, refinement_stack, refinement_stackN, Tri_li)
    END DO

    DEALLOCATE( Tri_flip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_triangle

  SUBROUTINE split_line(     mesh, v1a, v2a   , p_new, refinement_map, refinement_stack, refinement_stackN, Tri_li)
     ! Split the line between vertices v1a and v2a at point p_new

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: v1a, v2a
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p_new
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_map
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_stack
    INTEGER,                  OPTIONAL,  INTENT(INOUT) :: refinement_stackN
    INTEGER,  DIMENSION(:,:), OPTIONAL,  INTENT(INOUT) :: Tri_li

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'split_line'
    INTEGER                                            :: v1, v2, ti, t1, t2, n, nc, nnext, vo1, vo2
    INTEGER                                            :: t1new1, t1new2, t2new1, t2new2, t1nv1, t1nv2, t2nv1, t2nv2
    LOGICAL                                            :: AreConnected, SwitchThem
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER                                            :: li_min, li_max
    INTEGER,  DIMENSION(:,:), ALLOCATABLE              :: Tri_flip
    INTEGER                                            :: nf
    LOGICAL                                            :: did_flip

    ! Add routine to path
    CALL init_routine( routine_name)

    v1 = v1a
    v2 = v2a

    ! Check if they are even connected.
    AreConnected = .FALSE.
    DO n = 1, mesh%nC(v1)
      IF (mesh%C( v1,n) == v2) AreConnected = .TRUE.
    END DO
    IF (.NOT. AreConnected) THEN
      CALL crash('trying to split a non-existing line!')
    END IF

    ! If the line is actually a boundary segment, split that one
    IF (is_boundary_segment( mesh,v1,v2)) THEN
      p = (mesh%V( v1,:) + mesh%V( v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p, refinement_map, refinement_stack, refinement_stackN, Tri_li)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Find the triangles t1 and t2 that contain v1 and v2
    t1 = 0
    t2 = 0
    DO ti = 1, mesh%nTri
      nc = 0
      DO n = 1, 3
        IF (mesh%Tri( ti,n) == v1 .OR. mesh%Tri( ti,n) == v2) nc = nc+1
      END DO
      IF (nc == 2) THEN
        IF (t1 == 0) THEN
         t1 = ti
        ELSE
         t2 = ti
        END IF
      END IF
    END DO

    IF (t1 == 0 .OR. t2 == 0) THEN
      CALL crash('couldnt find two triangles containing both vertices!')
    END IF

    ! Order v1 and v2 anticlockwise in triangle t1
    SwitchThem = .FALSE.
    DO n = 1, 3
      nnext = n + 1
      IF (nnext == 4) nnext = 1
      IF (mesh%Tri( t1,n) == v1) THEN
        IF (mesh%Tri( t1,nnext) /= v2) SwitchThem = .TRUE.
      END IF
    END DO
    IF (SwitchThem) THEN
      v1 = v1 + v2
      v2 = v1 - v2
      v1 = v1 - v2
    END IF

    ! == Find the other relevant indices - non-shared vertices and neighbour triangles
    ! t1nv1: triangle next to t1, across from v1
    ! t1nv2: triangle next to t1, across from v2
    ! t2nv1: triangle next to t2, across from v1
    ! t2nv2: triangle next to t2, across from v2
    t1nv1 = 0
    t1nv2 = 0
    t2nv1 = 0
    t2nv2 = 0
    vo1   = 0
    vo2   = 0
    DO n = 1, 3
      IF ( mesh%Tri( t1,n) == v1) THEN
        t1nv1 = mesh%TriC( t1,n)
      ELSEIF ( mesh%Tri( t1,n) == v2) THEN
        t1nv2 = mesh%TriC( t1,n)
      ELSE
        vo1 = mesh%Tri( t1,n)
      END IF
      IF ( mesh%Tri( t2,n) == v1) THEN
        t2nv1 = mesh%TriC( t2,n)
      ELSEIF ( mesh%Tri( t2,n) == v2) THEN
        t2nv2 = mesh%TriC( t2,n)
      ELSE
        vo2 = mesh%Tri( t2,n)
      END IF
    END DO

    ! == Add new vertex to mesh
    mesh%nV = mesh%nV + 1
    mesh%V( mesh%nV,:) = p_new
    mesh%edge_index( mesh%nV) = 0

    ! == Create four new triangles
    t1new1 = t1
    t1new2 = mesh%nTri + 1
    t2new1 = t2
    t2new2 = mesh%nTri + 2
    mesh%Tri( t1new1,:) = [vo1, mesh%nV, v2]
    mesh%Tri( t1new2,:) = [vo1, v1, mesh%nV]
    mesh%Tri( t2new1,:) = [vo2, mesh%nV, v1]
    mesh%Tri( t2new2,:) = [vo2, v2, mesh%nV]
    mesh%nTri = mesh%nTri + 2

    ! Add these to the refinement stack
    IF (PRESENT( refinement_map)) THEN
      CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, t1)
      CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, t2)

      refinement_map( t1new1) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t1new1

      refinement_map( t1new2) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t1new2

      refinement_map( t2new1) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t2new1

      refinement_map( t2new2) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t2new2
    END IF ! IF (PRESENT( refinement_map)) THEN

    ! == Update triangle-line overlap ranges
    IF (PRESENT( Tri_li)) THEN
      li_min = MIN( Tri_li( t1,1), Tri_li( t2,1))
      li_max = MAX( Tri_li( t1,2), Tri_li( t2,2))
      Tri_li( t1new1,:) = [li_min, li_max]
      Tri_li( t1new2,:) = [li_min, li_max]
      Tri_li( t2new1,:) = [li_min, li_max]
      Tri_li( t2new2,:) = [li_min, li_max]
    END IF

    ! == Find circumcenters
    CALL update_triangle_circumcenter( mesh, t1new1)
    CALL update_triangle_circumcenter( mesh, t1new2)
    CALL update_triangle_circumcenter( mesh, t2new1)
    CALL update_triangle_circumcenter( mesh, t2new2)

    ! == Update inverse triangle lists
    ! vo1
    DO n = 1, mesh%niTri( vo1)
      IF (mesh%iTri( vo1,n) == t1) THEN
        mesh%iTri( vo1,:) = [mesh%iTri( vo1,1:n-1), t1new2, t1new1, mesh%iTri( vo1,n+1:mesh%nC_mem-1)]
        mesh%niTri( vo1) = mesh%niTri( vo1)+1
        EXIT
      END IF
    END DO
    ! v1
    DO n = 1, mesh%niTri( v1)
      IF (mesh%iTri( v1,n) == t1) mesh%iTri( v1,n) = t1new2
      IF (mesh%iTri( v1,n) == t2) mesh%iTri( v1,n) = t2new1
    END DO
    ! vo2
    DO n = 1, mesh%niTri( vo2)
      IF (mesh%iTri( vo2,n) == t2) THEN
        mesh%iTri( vo2,:) = [mesh%iTri( vo2,1:n-1), t2new2, t2new1, mesh%iTri( vo2,n+1:mesh%nC_mem-1)]
        mesh%niTri( vo2) = mesh%niTri( vo2)+1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%niTri( v2)
      IF (mesh%iTri( v2,n) == t1) mesh%iTri( v2,n) = t1new1
      IF (mesh%iTri( v2,n) == t2) mesh%iTri( v2,n) = t2new2
    END DO
    ! new vertex
    mesh%iTri(  mesh%nV,1:4) = [t1new2, t2new1, t2new2, t1new1]
    mesh%niTri( mesh%nV) = 4

    ! == Update triangle connectivity matrix
    ! t1nv2
    IF (t1nv2 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( t1nv2,n) == t1) mesh%TriC( t1nv2,n) = t1new2
      END DO
    END IF
    ! t1nv1
    IF (t1nv1 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( t1nv1,n) == t1) mesh%TriC( t1nv1,n) = t1new1
      END DO
    END IF
    ! t2nv2
    IF (t2nv2 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( t2nv2,n) == t2) mesh%TriC( t2nv2,n) = t2new1
      END DO
    END IF
    ! t2nv1
    IF (t2nv1 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( t2nv1,n) == t2) mesh%TriC( t2nv1,n) = t2new2
      END DO
    END IF

    ! The four new triangles
    mesh%TriC( t1new1,:) = [t2new2, t1nv1, t1new2]
    mesh%TriC( t1new2,:) = [t2new1, t1new1, t1nv2]
    mesh%TriC( t2new1,:) = [t1new2, t2nv2, t2new2]
    mesh%TriC( t2new2,:) = [t1new1, t2new1, t2nv1]

    ! == Update vertex connectivity lists
    ! po1
    DO n = 1, mesh%nC( vo1)
      IF (mesh%C( vo1,n) == v1) THEN
        mesh%C( vo1,:) = [mesh%C( vo1,1:n), mesh%nV, mesh%C( vo1,n+1:mesh%nC_mem-1)]
        mesh%nC( vo1) = mesh%nC( vo1)+1
        EXIT
      END IF
    END DO
    ! v1
    DO n = 1, mesh%nC( v1)
      IF (mesh%C( v1,n) == v2) THEN
        mesh%C( v1,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! vo2
    DO n = 1, mesh%nC( vo2)
      IF (mesh%C( vo2,n) == v2) THEN
        mesh%C( vo2,:) = [mesh%C( vo2,1:n), mesh%nV, mesh%C( vo2,n+1:mesh%nC_mem-1)]
        mesh%nC( vo2) = mesh%nC( vo2)+1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%nC( v2)
      IF (mesh%C( v2,n) == v1) THEN
        mesh%C( v2,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%C(  mesh%nV,1:4) = [v1, vo2, v2, vo1]
    mesh%nC( mesh%nV) = 4

    ! == Propagate flip operations outward
    ! Start with the newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    ALLOCATE( Tri_flip( mesh%nTri, 2), source = 0)
    nf = 0

    IF (t1nv1 > 0) THEN
      nf=nf + 1
      Tri_flip( nf,:) = [t1nv1, t1new1]
    END IF
    IF (t1nv2 > 0) THEN
      nf=nf + 1
      Tri_flip( nf,:) = [t1nv2, t1new2]
    END IF
    IF (t2nv1 > 0) THEN
      nf=nf + 1
      Tri_flip( nf,:) = [t2nv1, t2new2]
    END IF
    IF (t2nv2 > 0) THEN
      nf=nf + 1
      Tri_flip( nf,:) = [t2nv2, t2new1]
    END IF

    DO WHILE (nf > 0)
      CALL flip_triangle_pairs( mesh, Tri_flip, nf, did_flip, refinement_map, refinement_stack, refinement_stackN, Tri_li)
    END DO

    DEALLOCATE( Tri_flip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_line

  SUBROUTINE split_segment(  mesh, v1a, v2a   , p_new, refinement_map, refinement_stack, refinement_stackN, Tri_li)
    ! Split an Edge segment in two, adding a vertex halfway.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: v1a, v2a
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p_new
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_map
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_stack
    INTEGER,                  OPTIONAL,  INTENT(INOUT) :: refinement_stackN
    INTEGER,  DIMENSION(:,:), OPTIONAL,  INTENT(INOUT) :: Tri_li

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'split_segment'
    INTEGER                                            :: v1, v2, ti, t1, n, nc, nnext, tnv1, tnv2, vo
    INTEGER                                            :: tnew1, tnew2
    LOGICAL                                            :: SwitchThem
    INTEGER,  DIMENSION(:,:), ALLOCATABLE              :: Tri_flip
    INTEGER                                            :: nf
    LOGICAL                                            :: did_flip

    ! Add routine to path
    CALL init_routine( routine_name)

    v1 = v1a
    v2 = v2a

    ! Find the triangle t1 that contains v1 and v2
    t1 = 0
    DO ti = 1, mesh%nTri
      nc = 0
      DO n = 1, 3
        IF (mesh%Tri( ti,n) == v1 .OR. mesh%Tri( ti,n) == v2) nc = nc + 1
      END DO
      IF (nc == 2) t1 = ti
    END DO

    IF (t1 == 0) THEN
      CALL crash('couldnt find triangle containing both vertices!')
    END IF

    IF (mesh%edge_index( v1) == 0 .OR. mesh%edge_index( v2) == 0) THEN
      CALL crash('segment isnt made up of boundary vertices!')
    END IF

    ! Order v1 and v2 anticlockwise in triangle t1
    SwitchThem = .FALSE.
    DO n = 1, 3
      nnext = n + 1
      IF (nnext == 4) nnext = 1
      IF (mesh%Tri( t1,n) == v1) THEN
        IF (mesh%Tri( t1,nnext) /= v2) SwitchThem = .TRUE.
      END IF
    END DO
    IF (SwitchThem) THEN
      v1 = v1 + v2
      v2 = v1 - v2
      v1 = v1 - v2
    END IF

    ! == Find the other relevant indices - non-shared vertex and neighbour triangles
    ! tnv1: triangle next to t1, across from v1
    ! tnv2: triangle next to t1, across from v2
    tnv1 = 0
    tnv2 = 0
    vo   = 0
    DO n = 1, 3
      IF     (mesh%Tri( t1,n) == v1) THEN
        tnv1 = mesh%TriC( t1,n)
      ELSEIF (mesh%Tri( t1,n) == v2) THEN
        tnv2 = mesh%TriC( t1,n)
      ELSE
        vo = mesh%Tri( t1,n)
      END IF
    END DO

    ! == Add new vertex to mesh
    mesh%nV = mesh%nV + 1
    mesh%V( mesh%nV,:) = p_new

    ! == Determine edge index
    IF (       mesh%edge_index( v1) == 1) THEN
      mesh%edge_index(mesh%nV) = 1
    ELSEIF (   mesh%edge_index( v1) == 2) THEN
      IF (     mesh%edge_index( v2) == 1 .OR. mesh%edge_index( v2) == 8) THEN
        mesh%edge_index(mesh%nV) = 1
      ELSEIF ( mesh%edge_index( v2) == 3 .OR. mesh%edge_index( v2) == 4) THEN
        mesh%edge_index(mesh%nV) = 3
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSEIF (   mesh%edge_index( v1) == 3) THEN
      mesh%edge_index(mesh%nV) = 3
    ELSEIF (   mesh%edge_index( v1) == 4) THEN
      IF (     mesh%edge_index( v2) == 3 .OR. mesh%edge_index( v2) == 2) THEN
        mesh%edge_index(mesh%nV) = 3
      ELSEIF ( mesh%edge_index( v2) == 5 .OR. mesh%edge_index( v2) == 6) THEN
        mesh%edge_index(mesh%nV) = 5
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSEIF (   mesh%edge_index( v1) == 5) THEN
      mesh%edge_index(mesh%nV) = 5
    ELSEIF (   mesh%edge_index( v1) == 6) THEN
      IF (     mesh%edge_index( v2) == 5 .OR. mesh%edge_index( v2) == 4) THEN
        mesh%edge_index(mesh%nV) = 5
      ELSEIF ( mesh%edge_index( v2) == 7 .OR. mesh%edge_index( v2) == 8) THEN
        mesh%edge_index(mesh%nV) = 7
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSEIF (   mesh%edge_index( v1) == 7) THEN
      mesh%edge_index(mesh%nV) = 7
    ELSEIF (   mesh%edge_index( v1) == 8) THEN
      IF (     mesh%edge_index( v2) == 7 .OR. mesh%edge_index( v2) == 6) THEN
        mesh%edge_index(mesh%nV) = 7
      ELSEIF ( mesh%edge_index( v2) == 1 .OR. mesh%edge_index( v2) == 2) THEN
        mesh%edge_index(mesh%nV) = 1
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSE
      CALL crash('edge indices of v1 and v2 dont make sense!')
    END IF

    ! == Create two new triangles
    tnew1 = t1
    tnew2 = mesh%nTri + 1
    mesh%Tri( tnew1,:) = [v1, mesh%nV, vo]
    mesh%Tri( tnew2,:) = [v2, vo, mesh%nV]
    mesh%nTri = mesh%nTri + 1

    ! Add these to the refinement stack
    IF (PRESENT( refinement_map)) THEN
      CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, t1)

      refinement_map( tnew1) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = tnew1

      refinement_map( tnew2) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = tnew2
    END IF ! IF (PRESENT( refinement_map)) THEN

    ! == Update triangle-line overlap ranges
    IF (PRESENT( Tri_li)) THEN
      Tri_li( tnew1,:) = Tri_li( t1,:)
      Tri_li( tnew2,:) = Tri_li( t1,:)
    END IF

    ! == Update triangle circumcenters
    CALL update_triangle_circumcenter( mesh, tnew1)
    CALL update_triangle_circumcenter( mesh, tnew2)

    ! == Update inverse triangle lists
    ! vo
    DO n = 1, mesh%niTri( vo)
      IF (mesh%iTri( vo,n) == t1) THEN
        mesh%iTri(  vo,:) = [mesh%iTri( vo,1:n-1), tnew1, tnew2, mesh%iTri( vo,n+1:mesh%nC_mem-1)]
        mesh%niTri( vo  ) = mesh%niTri( vo) + 1
        EXIT
      END IF
    END DO
    ! v1 - nothing changes here
    ! v2
    DO n = 1, mesh%niTri(v2)
      IF (mesh%iTri( v2,n) == t1) mesh%iTri( v2,n) = tnew2
    END DO
    ! new vertex
    mesh%iTri(  mesh%nV,1:2) = [tnew2, tnew1]
    mesh%niTri( mesh%nV) = 2

    ! == Update triangle connectivity lists
    ! tnv1
    IF (tnv1 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tnv1,n) == t1) mesh%TriC( tnv1,n) = tnew2
      END DO
    END IF
    ! tnv2
    IF (tnv2 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tnv2,n) == t1) mesh%TriC( tnv2,n) = tnew1
      END DO
    END IF

    ! The two new triangles
    mesh%TriC( tnew1,:) = [tnew2, tnv2, 0]
    mesh%TriC( tnew2,:) = [tnew1, 0, tnv1]

    ! == Update vertex connectivity lists
    ! vo
    DO n = 1, mesh%nC( vo)
      IF (mesh%C( vo,n) == v1) THEN
        mesh%C(  vo,:) = [mesh%C( vo,1:n), mesh%nV, mesh%C( vo,n+1:mesh%nC_mem-1)]
        mesh%nC( vo  ) = mesh%nC( vo)+1
        EXIT
      END IF
    END DO
    ! v1
    DO n = 1, mesh%nC( v1)
      IF (mesh%C( v1,n) == v2) THEN
        mesh%C( v1,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%nC( v2)
      IF (mesh%C( v2,n) == v1) THEN
        mesh%C( v2,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%C(  mesh%nV,1:3) = [v2, vo, v1]
    mesh%nC( mesh%nV) = 3

    ! == Propagate flip operations outward
    ! Start with the newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    ALLOCATE( Tri_flip( mesh%nTri, 2), source = 0)
    nf = 0

    IF (tnv2 > 0) THEN
      nf = nf + 1
      Tri_flip( nf,:) = [tnew1, tnv2]
    END IF
    IF (tnv1 > 0) THEN
      nf = nf + 1
      Tri_flip( nf,:) = [tnew2, tnv1]
    END IF

    DO WHILE (nf > 0)
      CALL flip_triangle_pairs( mesh, Tri_flip, nf, did_flip, refinement_map, refinement_stack, refinement_stackN, Tri_li)
    END DO

    DEALLOCATE( Tri_flip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_segment

  SUBROUTINE move_vertex(    mesh, vi         , p_new, refinement_map, refinement_stack, refinement_stackN, Tri_li)
    ! Move vertex vi of the mesh to point p_new

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p_new
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_map
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_stack
    INTEGER,                  OPTIONAL,  INTENT(INOUT) :: refinement_stackN
    INTEGER,  DIMENSION(:,:), OPTIONAL,  INTENT(INOUT) :: Tri_li

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'move_vertex'
    INTEGER                                            :: iti, ti, t1, t2, n
    INTEGER,  DIMENSION(:,:), ALLOCATABLE              :: Tri_flip
    INTEGER                                            :: nf
    LOGICAL                                            :: did_flip, did_flip_pair

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( Tri_flip( mesh%nTri, 2), source = 0)
    nf = 0

    ! Move the vertex
    mesh%V( vi,:) = p_new

    ! Update surrounding triangle circumcentres
    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      CALL update_triangle_circumcenter( mesh, ti)
    END DO

    ! Update triangulation
    did_flip = .TRUE.
    DO WHILE (did_flip)

      did_flip = .FALSE.

      Tri_flip = 0
      nf       = 0

      DO iti = 1, mesh%niTri( vi)
        t1 = mesh%iTri( vi,iti)
        DO n = 1, 3
          t2 = mesh%TriC( t1,n)
          IF (t2 > 0) THEN
            nf = nf + 1
            Tri_flip( nf,:) = [t1,t2]
          END IF
        END DO
      END DO

      DO WHILE (nf > 0)
        CALL flip_triangle_pairs( mesh, Tri_flip, nf, did_flip_pair, refinement_map, refinement_stack, refinement_stackN, Tri_li)
        IF (did_flip_pair) did_flip = .TRUE.
      END DO

    END DO ! DO WHILE (did_flip)

    DEALLOCATE( Tri_flip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE move_vertex

  SUBROUTINE flip_triangle_pairs( mesh, Tri_flip, nf, did_flip, refinement_map, refinement_stack, refinement_stackN, Tri_li)
    ! Flip adjacent triangles, if possible and neccesary. Add new triangle
    ! pairs to the list.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,  DIMENSION( mesh%nTri,2),   INTENT(INOUT) :: Tri_flip
    INTEGER,                             INTENT(INOUT) :: nf
    LOGICAL,                             INTENT(OUT)   :: did_flip
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_map
    INTEGER,  DIMENSION(:  ), OPTIONAL,  INTENT(INOUT) :: refinement_stack
    INTEGER,                  OPTIONAL,  INTENT(INOUT) :: refinement_stackN
    INTEGER,  DIMENSION(:,:), OPTIONAL,  INTENT(INOUT) :: Tri_li

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_triangle_pairs'
    INTEGER                                            :: t1, t2, n, vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2
    LOGICAL                                            :: n1to2, n2to1, FlipThem
    INTEGER                                            :: li_min, li_max

    ! Add routine to path
    CALL init_routine( routine_name)

    did_flip = .FALSE.

    t1 = Tri_flip( 1,1)
    t2 = Tri_flip( 1,2)

    IF (t1 == 0 .OR. t2 == 0) THEN
      CALL crash('received t=0!')
    END IF

    ! == First, check if the two are really adjacent. If not, that's because of an earlier flip operation.
    ! The current one is now redundant, so remove it from the list.
    n1to2 = .FALSE.
    n2to1 = .FALSE.
    DO n = 1, 3
      IF (mesh%TriC( t1,n)==t2) n1to2 = .TRUE.
      IF (mesh%TriC( t2,n)==t1) n2to1 = .TRUE.
    END DO

    IF ((n1to2 .AND. .NOT. n2to1) .OR. (n2to1 .AND. .NOT. n1to2)) THEN
      CALL crash('somethings really wrong with the triangle connectivity matrix!')
    END IF
    IF (.NOT. n1to2 .AND. .NOT. n2to1) THEN
      ! The two triangles are no longer connected; remove them from the list and return.
      Tri_flip( 1:mesh%nTri-1,:) = Tri_flip( 2:mesh%nTri,:)
      nf = nf - 1
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == Check if a flip is necessary
    ! If not, remove the pair from the flip list.
    CALL need_flipping( mesh, t1, t2, FlipThem, vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2)

    IF (.NOT. FlipThem) THEN
      ! The two triangles do not need to be flipped; remove them from the list and return.
      Tri_flip( 1:mesh%nTri-1,:) = Tri_flip( 2:mesh%nTri,:)
      nf = nf - 1
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == Flip them

    did_flip = .TRUE.

    ! WRITE(0,'(A,I6,A,I6)') '     Flipping triangles ', t1, ' and ', t2

    ! == Update the triangle matrix
    mesh%Tri( t1,:) = [vo1, v1, vo2]
    mesh%Tri( t2,:) = [vo2, v2, vo1]

    ! Add the neighbouring triangles to the refinement stack
    IF (PRESENT( refinement_map)) THEN
      CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, t1)
      CALL remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, t2)

      refinement_map( t1) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t1

      refinement_map( t2) = 1
      refinement_stackN = refinement_stackN + 1
      refinement_stack( refinement_stackN) = t2
    END IF ! IF (PRESENT( refinement_map)) THEN

    ! == Update triangle-line overlap ranges
    IF (PRESENT( Tri_li)) THEN
      li_min = MIN( Tri_li( t1,1), Tri_li( t2,1))
      li_max = MAX( Tri_li( t1,2), Tri_li( t2,2))
      Tri_li( t1,:) = [li_min, li_max]
      Tri_li( t2,:) = [li_min, li_max]
    END IF

    ! == Update the triangle connectivity matrix
    ! t1nv1
    IF (t1nv1 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( t1nv1,n) == t1) mesh%TriC( t1nv1,n) = t2
      END DO
    END IF
    ! t1nv2: nothing changes
    ! t2nv1: nothing changes
    ! t2nv2
    IF (t2nv2 > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( t2nv2,n) == t2) mesh%TriC( t2nv2,n) = t1
      END DO
    END IF
    ! The two new triangles
    mesh%TriC( t1,:) = [t2nv2, t2, t1nv2]
    mesh%TriC( t2,:) = [t1nv1, t1, t2nv1]

    ! == Update inverse triangle lists
    ! v1
    DO n = 1, mesh%niTri( v1)
      IF ( mesh%iTri( v1,n) == t2) THEN
        mesh%iTri( v1,:) = [mesh%iTri( v1,1:n-1), mesh%iTri( v1,n+1:mesh%nC_mem), 0]
        mesh%niTri( v1) = mesh%niTri( v1)-1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%niTri( v2)
      IF ( mesh%iTri( v2,n) == t1) THEN
        mesh%iTri( v2,:) = [mesh%iTri( v2,1:n-1), mesh%iTri( v2,n+1:mesh%nC_mem), 0]
        mesh%niTri( v2) = mesh%niTri( v2)-1
        EXIT
      END IF
    END DO
    ! vo1
    DO n = 1, mesh%niTri( vo1)
      IF ( mesh%iTri( vo1,n) == t1) THEN
        mesh%iTri( vo1,:) = [mesh%iTri( vo1,1:n), t2, mesh%iTri( vo1,n+1:mesh%nC_mem-1)]
        mesh%niTri( vo1) = mesh%niTri( vo1) + 1
        EXIT
      END IF
    END DO
    ! vo2
    DO n = 1, mesh%niTri( vo2)
      IF ( mesh%iTri( vo2,n) == t2) THEN
        mesh%iTri( vo2,:) = [mesh%iTri( vo2,1:n), t1, mesh%iTri( vo2,n+1:mesh%nC_mem-1)]
        mesh%niTri( vo2) = mesh%niTri( vo2) + 1
        EXIT
      END IF
    END DO

    ! == Update vertex connectivity lists
    ! v1
    DO n = 1, mesh%nC( v1)
      IF ( mesh%C( v1,n) == v2) THEN
        mesh%C(  v1,:) = [mesh%C( v1,1:n-1), mesh%C( v1,n+1:mesh%nC_mem), 0]
        mesh%nC( v1  ) = mesh%nC( v1) - 1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%nC( v2)
      IF ( mesh%C( v2,n) == v1) THEN
        mesh%C(  v2,:) = [mesh%C( v2,1:n-1), mesh%C( v2,n+1:mesh%nC_mem), 0]
        mesh%nC( v2  ) = mesh%nC( v2) - 1
        EXIT
      END IF
    END DO
    ! vo1
    DO n = 1, mesh%nC( vo1)
      IF ( mesh%C( vo1,n) == v1) THEN
        mesh%C(  vo1,:) = [mesh%C( vo1,1:n), vo2, mesh%C( vo1,n+1:mesh%nC_mem-1)]
        mesh%nC( vo1  ) = mesh%nC( vo1) + 1
        EXIT
      END IF
    END DO
    ! vo2
    DO n = 1, mesh%nC( vo2)
      IF ( mesh%C( vo2,n) == v2) THEN
        mesh%C(  vo2,:) = [mesh%C( vo2,1:n), vo1, mesh%C( vo2,n+1:mesh%nC_mem-1)]
        mesh%nC( vo2  ) = mesh%nC( vo2) + 1
        EXIT
      END IF
    END DO

    ! == Update triangle circumcenters
   CALL update_triangle_circumcenter( mesh, t1)
   CALL update_triangle_circumcenter( mesh, t2)

    ! == Remove current triangle pair from flip list, add 4 new ones
    Tri_flip( 1:mesh%nTri-1,:) = Tri_flip( 2:mesh%nTri,:)
    nf = nf-1

    IF ( t1nv1 > 0) THEN
      Tri_flip( 2:mesh%nTri,:) = Tri_flip( 1:mesh%nTri-1,:)
      Tri_flip( 1,:) = [t1nv1, t2]
      nf = nf + 1
    END IF
    IF ( t1nv2 > 0) THEN
      Tri_flip( 2:mesh%nTri,:) = Tri_flip( 1:mesh%nTri-1,:)
      Tri_flip( 1,:) = [t1nv2, t1]
      nf = nf + 1
    END IF
    IF ( t2nv1 > 0) THEN
      Tri_flip( 2:mesh%nTri,:) = Tri_flip( 1:mesh%nTri-1,:)
      Tri_flip( 1,:) = [t2nv1, t2]
      nf = nf + 1
    END IF
    IF ( t2nv2 > 0) THEN
      Tri_flip( 2:mesh%nTri,:) = Tri_flip( 1:mesh%nTri-1,:)
      Tri_flip( 1,:) = [t2nv2, t1]
      nf = nf + 1
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE flip_triangle_pairs

  SUBROUTINE need_flipping( mesh, t1, t2, isso, vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2)
    ! Check if triangle pair [t1,t2] meets the local Delaunay criterion. If not, they must be flipped.
    !
    ! Also return some general info about the local mesh geometry, which will come in handy later.
    ! Shared vertices v1 and v2 (sorted clockwise in t1), non-shared
    ! vertices vo1 and vo2, neigbours to t1 t1nv1 (across from v1) and
    ! t1nv2 (across from v2) and neighbours to t2 t2nv1 and t2nv2 (idem)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: t1, t2
    LOGICAL,                             INTENT(OUT)   :: isso
    INTEGER,                             INTENT(OUT)   :: vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'need_flipping'
    REAL(dp), DIMENSION(2)                             :: p, q, r, s
    INTEGER                                            :: n, n1, n2, nnext
    LOGICAL                                            :: isint2, SwitchThem

    ! Add routine to path
    CALL init_routine( routine_name)

    v1  = 0
    v2  = 0
    vo1 = 0
    DO n1 = 1, 3
      isint2 = .FALSE.
      DO n2 = 1, 3
        IF (mesh%Tri( t1,n1) == mesh%Tri( t2,n2)) THEN
          isint2 = .TRUE.
          IF (v1 == 0) THEN
            v1 = mesh%Tri( t1,n1)
          ELSE
            v2 = mesh%Tri( t1,n1)
          END IF
        END IF
      END DO ! DO n2 = 1, 3
      IF (.NOT. isint2) vo1 = mesh%Tri( t1,n1)
    END DO ! DO n1 = 1, 3
    DO n = 1, 3
      IF (mesh%Tri( t2,n) /= v1 .AND. mesh%Tri( t2,n) /= v2) vo2 = mesh%Tri( t2,n)
    END DO ! DO n = 1, 3

    ! Order v1 and v2 anticlockwise in triangle t1
    SwitchThem = .FALSE.
    DO n = 1, 3
      nnext = n + 1
      IF (nnext == 4) nnext = 1
      IF (mesh%Tri( t1,n) == v1) THEN
        IF (mesh%Tri( t1,nnext) /= v2) SwitchThem = .TRUE.
      END IF
    END DO ! DO n = 1, 3
    IF (SwitchThem) THEN
      v1 = v1 + v2
      v2 = v1 - v2
      v1 = v1 - v2
    END IF

    ! == Find neighbour triangles
    DO n = 1, 3
      IF (mesh%Tri( t1,n) == v1) t1nv1 = mesh%TriC( t1,n)
      IF (mesh%Tri( t1,n) == v2) t1nv2 = mesh%TriC( t1,n)
      IF (mesh%Tri( t2,n) == v1) t2nv1 = mesh%TriC( t2,n)
      IF (mesh%Tri( t2,n) == v2) t2nv2 = mesh%TriC( t2,n)
    END DO

    ! == Determine if triangle pair t1,t2 requires flipping
    isso = .FALSE.
!    IF (norm2(mesh%V(vo2,:) - mesh%Tricc(t1,:)) < norm2(mesh%V(mesh%Tri(t1,1),:) - mesh%Tricc(t1,:)) - mesh%tol_dist .OR. &
!        norm2(mesh%V(vo1,:) - mesh%Tricc(t2,:)) < norm2(mesh%V(mesh%Tri(t2,1),:) - mesh%Tricc(t2,:)) - mesh%tol_dist) THEN
    IF ( NORM2( mesh%V( vo2,:) - mesh%Tricc( t1,:)) < NORM2( mesh%V( mesh%Tri( t1,1),:) - mesh%Tricc( t1,:)) .OR. &
        NORM2( mesh%V( vo1,:) - mesh%Tricc( t2,:)) < NORM2( mesh%V( mesh%Tri( t2,1),:) - mesh%Tricc( t2,:))) THEN
      isso = .TRUE.
    END IF

    ! If the outer angle at v1 or v2 is concave, don't flip.
    ! Check this by checking if v1 lies inside the triangle
    ! [vo2,vo2,v2], or the other way round.
    p = mesh%V( vo1,:)
    q = mesh%V( vo2,:)
    r = mesh%V( v1,:)
    s = mesh%V( v2,:)
    IF  (is_in_triangle( p, q, s, r) .OR. &
         is_in_triangle( p, q, r, s)) THEN
      isso = .FALSE.
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE need_flipping

  SUBROUTINE remove_triangle_from_refinement_stack( refinement_map, refinement_stack, refinement_stackN, ti)
    ! Remove triangle ti from the refinement map and stack
    ! (So that triangles in the stack can be sort-of sorted from big to small)

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,  DIMENSION(:),              INTENT(INOUT) :: refinement_map
    INTEGER,  DIMENSION(:),              INTENT(INOUT) :: refinement_stack
    INTEGER,                             INTENT(INOUT) :: refinement_stackN
    INTEGER,                             INTENT(IN)    :: ti

    ! Local variables:
    INTEGER                                            :: i
    LOGICAL                                            :: foundit

    ! Map
!    IF (refinement_map( ti) == 0) CALL crash('triangle ti was not marked for refinement on the map!')
    refinement_map( ti) = 0

    ! Stack
    foundit = .FALSE.
    DO i = 1, refinement_stackN
      IF (refinement_stack( i) == ti) THEN
        foundit = .TRUE.
        refinement_stack( i:refinement_stackN-1) = refinement_stack( i+1:refinement_stackN)
        refinement_stack(   refinement_stackN  ) = 0
        refinement_stackN = refinement_stackN - 1
        EXIT
      END IF
    END DO
!    IF (.NOT. foundit) CALL crash('triangle ti could not be found in the refinement stack!')

  END SUBROUTINE remove_triangle_from_refinement_stack

END MODULE mesh_Delaunay
