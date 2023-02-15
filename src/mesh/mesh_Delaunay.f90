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
  USE math_utilities                                         , ONLY: is_in_triangle, lies_on_line_segment, line_from_points, line_line_intersection, &
                                                                     perpENDicular_bisector_from_line, encroaches_upon
  USE mesh_utilities                                         , ONLY: update_triangle_circumcenter, find_containing_triangle, is_boundary_segment, &
                                                                     encroaches_upon_any, check_mesh, add_triangle_to_refinement_stack_first, &
                                                                     add_triangle_to_refinement_stack_last, check_if_triangle_already_exists, &
                                                                     write_mesh_to_text_file

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  LOGICAL :: do_debug = .FALSE.

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE split_triangle( mesh, ti_in_guess, p_new)
    ! Add a vertex at p_new, splitting triangle ti into three new ones
    !
    ! Provide a guess ti_in_guess for which triangle we think contains p (can be wrong,
    ! but guessing near the correct one speeds up the code).
    !
    ! If p_new coincides with a line or boundary segment, split that instead.
    !
    ! When going in, the local geometry looks like this:
    !
    !   \ /           \ /           \ /
    ! - -o----------- vic -----------o- -
    !   / \           / \           / \
    !      \   tb    /   \   ta    /
    !       \       /     \       /
    !        \     /       \     /
    !         \   /   ti    \   /
    !          \ /           \ /
    !      - - via --------- vib - -
    !          / \           / \
    !             \   tc    /
    !              \       /
    !               \     /
    !                \   /
    !                 \ /
    !               - -o- -
    !                 / \
    !
    ! When coming out, it looks like this:
    !
    !   \ /           \ /           \ /
    ! _ _o____________ vc ___________o_ _
    !   / \           /|\           / \
    !      \   tb    / | \   ta    /
    !       \       /  |  \       /
    !        \     /t1 | t2\     /
    !         \   /  _ vk_  \   /
    !          \ / /   t3  \ \ /
    !      _ _ va ___________ vb _ _
    !          / \           / \
    !             \   tc    /
    !              \       /
    !               \     /
    !                \   /
    !                 \ /
    !               _ _o_ _
    !                 / \

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti_in_guess
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'split_triangle'
    LOGICAL                                       :: isso
    REAL(dp), DIMENSION(2)                        :: p
    INTEGER                                       :: ci, iti, n, nf, t1, t2, t3, ti, tia, tib, tic
    REAL(dp), DIMENSION(2)                        :: va, vb, vc
    INTEGER                                       :: vi, via, vib, vic, vik, vj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If p_new encroaches upon a boundary segment, split that instead
    CALL encroaches_upon_any( mesh, p_new, isso, vi, vj)
    IF (isso) THEN
      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      CALL split_segment( mesh, vi, vj, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == If p_new lies outside of the mesh domain, split a boundary segment instead.
    ! ==============================================================================

    IF (p_new( 1) < mesh%xmin) THEN
      ! p_new lies to the west of the mesh domain

      ! Safety
      IF (p_new( 2) < mesh%ymin .OR. p_new( 2) > mesh%ymax) THEN
        CALL crash('p_new lies way outside mesh domain!')
      END IF

      ! Find the two vertices vi,vj of the segment that must be split.
      vi = 1
      vj = mesh%C( vi, mesh%nC( vi))
      DO WHILE (mesh%V( vj,2) < p_new( 2))
        vi = vj
        vj = mesh%C( vi, mesh%nC( vi))
      END DO

      ! Safety
      IF (.NOT. (p_new( 2) >= mesh%V( vi,2) .AND. p_new( 2) <= mesh%V( vj,2))) THEN
        CALL crash('couldnt find boundary segment to split!')
      END IF

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      CALL split_segment( mesh, vi, vj, p)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 1) < mesh%xmin)

    IF (p_new( 1) > mesh%xmax) THEN
      ! p_new lies to the east of the mesh domain

      ! Safety
      IF (p_new( 2) < mesh%ymin .OR. p_new( 2) > mesh%ymax) THEN
        CALL crash('p_new lies way outside mesh domain!')
      END IF

      ! Find the two vertices vi,vj of the segment that must be split.
      vi = 2
      vj = mesh%C( vi,1)
      DO WHILE (mesh%V( vj,2) < p_new( 2))
        vi = vj
        vj = mesh%C( vi, 1)
      END DO

      ! Safety
      IF (.NOT. (p_new( 2) >= mesh%V( vi,2) .AND. p_new( 2) <= mesh%V( vj,2))) THEN
        CALL crash('couldnt find boundary segment to split!')
      END IF

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      CALL split_segment( mesh, vi, vj, p)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 1) > mesh%xmax)

    IF (p_new( 2) < mesh%ymin) THEN
      ! p_new lies to the south of the mesh domain

      ! Safety
      IF (p_new( 1) < mesh%xmin .OR. p_new( 1) > mesh%xmax) THEN
        CALL crash('p_new lies way outside mesh domain!')
      END IF

      ! Find the two vertices vi,vj of the segment that must be split.
      vi = 1
      vj = mesh%C( vi, 1)
      DO WHILE (mesh%V( vj,1) < p_new( 1))
        vi = vj
        vj = mesh%C( vi, 1)
      END DO

      ! Safety
      IF (.NOT. (p_new( 1) >= mesh%V( vi,1) .AND. p_new( 1) <= mesh%V( vj,1))) THEN
        CALL crash('couldnt find segment to split!')
      END IF

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      CALL split_segment( mesh, vi, vj, p)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 2) < mesh%ymin)

    IF (p_new( 2) > mesh%ymax) THEN
      ! p_new lies to the north of the mesh domain

      ! Safety
      IF (p_new( 1) < mesh%xmin .OR. p_new( 1) > mesh%xmax) THEN
        CALL crash('p_new lies way outside mesh domain!')
      END IF

      ! Find the two vertices vi,vj of the segment that must be split.
      vi = 1
      vj = mesh%C( vi, mesh%nC( vi))
      DO WHILE (mesh%V( vj,1) < p_new( 1))
        vi = vj
        vj = mesh%C( vi, mesh%nC( vi))
      END DO

      ! Safety
      IF (.NOT. (p_new( 1) >= mesh%V( vi,1) .AND. p_new( 1) <= mesh%V( vj,1))) THEN
        CALL crash('couldnt find segment to split!')
      END IF

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      CALL split_segment( mesh, vi, vj, p)

      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (p_new( 2) > mesh%ymax)

    ! == Find the triangle containing p_new.
    ! ======================================

    ! Find the triangle containing p_new
    ti = ti_in_guess
    CALL find_containing_triangle( mesh, p_new, ti)

    ! The indices of the three vertices [a,b,c] spanning t_old
    via = mesh%Tri( ti,1)
    vib = mesh%Tri( ti,2)
    vic = mesh%Tri( ti,3)

    ! The coordinates of the three vertices [a,b,c] spanning t_old
    va  = mesh%V( via,:)
    vb  = mesh%V( vib,:)
    vc  = mesh%V( vic,:)

    ! Indice of the triangles neighbouring t_old
    tia = mesh%TriC( ti,1)
    tib = mesh%TriC( ti,2)
    tic = mesh%TriC( ti,3)

    ! == Check IF p_new encroaches upon a boundary segment. If so, split it.
    ! ======================================================================

    IF (is_boundary_segment( mesh, via, vib) .AND. encroaches_upon( va, vb, p_new)) THEN
      p = (mesh%V( via,:) + mesh%V( vib,:)) / 2._dp
      CALL split_segment( mesh, via, vib, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (is_boundary_segment( mesh, vib, vic) .AND. encroaches_upon( va, vb, p_new)) THEN
      p = (mesh%V( vib,:) + mesh%V( vic,:)) / 2._dp
      CALL split_segment( mesh, vib, vic, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (is_boundary_segment( mesh, vic, via) .AND. encroaches_upon( va, vb, p_new)) THEN
      p = (mesh%V( vic,:) + mesh%V( via,:)) / 2._dp
      CALL split_segment( mesh, vic, via, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == If the new vertex is (almost) colinear with two vertices of
    !    the containing triangle, split the line between them instead.
    ! ================================================================

    IF     (lies_on_line_segment( va, vb, p_new, mesh%tol_dist)) THEN
      CALL split_line( mesh, via, vib, p_new)
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (lies_on_line_segment( vb, vc, p_new, mesh%tol_dist)) THEN
      CALL split_line( mesh, vib, vic, p_new)
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (lies_on_line_segment( vc, va, p_new, mesh%tol_dist)) THEN
      CALL split_line( mesh, vic, via, p_new)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == All safety checks passes; split the triangle ti into three new ones.
    ! =======================================================================

    ! DENK DROM
    IF (do_debug) CALL warning('splitting triangle {int_01}', int_01 = ti)

    ! == V, Tri

    ! Create a new vertex vik at p_new
    mesh%nV = mesh%nV + 1
    vik = mesh%nV
    mesh%V( vik,:) = p_new

    ! Let triangle t1 be spanned by [vik,vic,via]
    t1 = ti
    mesh%Tri( t1,:) = [vik, vic, via]

    ! Let triangle t2 be spanned by [vik,vib,vic]
    mesh%nTri = mesh%nTri + 1
    t2 = mesh%nTri
    mesh%Tri( t2,:) = [vik,vib,vic]

    ! Let triangle t3 be spanned by [vik,via,vib]
    mesh%nTri = mesh%nTri + 1
    t3 = mesh%nTri
    mesh%Tri( t3,:) = [vik,via,vib]

    ! DENK DROM
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t1)
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t2)
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t3)

    ! == nC, C

    ! via: connection to vik is addded between vib and vic
    mesh%nC( via) = mesh%nC( via) + 1
    DO ci = 1, mesh%nC( via)
      IF (mesh%C( via,ci) == vib) THEN
        mesh%C( via,:) = [mesh%C( via,1:ci), vik, mesh%C( via,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vib: connection to vik is addded between vic and via
    mesh%nC( vib) = mesh%nC( vib) + 1
    DO ci = 1, mesh%nC( vib)
      IF (mesh%C( vib,ci) == vic) THEN
        mesh%C( vib,:) = [mesh%C( vib,1:ci), vik, mesh%C( vib,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vic: connection to vik is addded between via and vib
    mesh%nC( vic) = mesh%nC( vic) + 1
    DO ci = 1, mesh%nC( vic)
      IF (mesh%C( vic,ci) == via) THEN
        mesh%C( vic,:) = [mesh%C( vic,1:ci), vik, mesh%C( vic,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vik: connected to [via, vib, vic]
    mesh%nC( vik) = 3
    mesh%C( vik,1:3) = [via, vib, vic]

    ! == niTri, iTri

    ! via: ti is replaced by [t3,t1]
    mesh%niTri( via) = mesh%niTri( via) + 1
    DO iti = 1, mesh%niTri( via)
      IF (mesh%iTri( via,iti) == ti) THEN
        mesh%iTri( via,:) = [mesh%iTri( via,1:iti-1), t3, t1, mesh%iTri( via,iti+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vib: ti is replaced by [t2,t3]
    mesh%niTri( vib) = mesh%niTri( vib) + 1
    DO iti = 1, mesh%niTri( vib)
      IF (mesh%iTri( vib,iti) == ti) THEN
        mesh%iTri( vib,:) = [mesh%iTri( vib,1:iti-1), t2, t3, mesh%iTri( vib,iti+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vic: ti is replaced by [t1,t2]
    mesh%niTri( vic) = mesh%niTri( vic) + 1
    DO iti = 1, mesh%niTri( vic)
      IF (mesh%iTri( vic,iti) == ti) THEN
        mesh%iTri( vic,:) = [mesh%iTri( vic,1:iti-1), t1, t2, mesh%iTri( vic,iti+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vik: surrounded by [t1, t3, t2]
    mesh%niTri( vik) = 3
    mesh%iTri( vik,1:3) = [t1, t3, t2]

    ! == Boundary index

    mesh%VBI( vik) = 0

    ! == TriC

    ! tia: ti is replaced by t2
    IF (tia > 0) THEN
      DO n = 1,3
        IF (mesh%TriC( tia,n) == ti) THEN
          mesh%TriC( tia,n) = t2
          EXIT
        END IF
      END DO
    END IF
    ! tib: ti is replaced by t1
    IF (tib > 0) THEN
      DO n = 1,3
        IF (mesh%TriC( tib,n) == ti) THEN
          mesh%TriC( tib,n) = t1
          EXIT
        END IF
      END DO
    END IF
    ! tic: ti is replaced by t3
    IF (tic > 0) THEN
      DO n = 1,3
        IF (mesh%TriC( tic,n) == ti) THEN
          mesh%TriC( tic,n) = t3
          EXIT
        END IF
      END DO
    END IF
    ! t1: [tib, t3, t2]
    mesh%TriC( t1,:) = [tib, t3, t2]
    ! t2: [tia, t1, t3]
    mesh%TriC( t2,:) = [tia, t1, t3]
    ! t3: [tic, t2, t1]
    mesh%TriC( t3,:) = [tic, t2, t1]

    ! == Tricc

    CALL update_triangle_circumcenter( mesh, t1)
    CALL update_triangle_circumcenter( mesh, t2)
    CALL update_triangle_circumcenter( mesh, t3)

    ! == Refinement data

    ! Add the three new triangles to the refinement map and stack
    CALL add_triangle_to_refinement_stack_last( mesh, t1)
    CALL add_triangle_to_refinement_stack_last( mesh, t2)
    CALL add_triangle_to_refinement_stack_last( mesh, t3)

    ! Update triangle-line overlap ranges
    mesh%Tri_li( t1,:) = mesh%Tri_li( ti,:)
    mesh%Tri_li( t2,:) = mesh%Tri_li( ti,:)
    mesh%Tri_li( t3,:) = mesh%Tri_li( ti,:)

    ! == Finished splitting the triangle. Iteratively flip any triangle pairs
    !    in the neighbourhood that do not meet the local Delaunay criterion.
    ! ======================================================================

    ! Start with the four newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Tri_flip_list( :,:) = 0
    nf = 0

    IF (tia > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t2, tia]
    END IF
    IF (tib > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t1, tib]
    END IF
    IF (tic > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t3, tic]
    END IF

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere
    CALL flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_triangle

  SUBROUTINE split_line( mesh, vi, vj, p_new)
    ! Split the triangles t_left and t_right adjacent to line [vi,vj]
    ! into four new ones. If [vi,vj] is a boundary segment, split that instead.
    !
    ! When going in, the local geometry looks like this:
    !
    !   \ /            \ /            \ /
    ! - -o----------- vtop ------------o- -
    !   / \            / \            / \
    !      \          /   \          /
    !       \  ttopl /     \ ttopr  /
    !        \      /       \      /
    !         \    /   ttop  \    /
    !          \  /           \  /
    !       - - vi ----------- vj - -
    !          /  \           /  \
    !         /    \         /    \
    !        /      \  tbot /      \
    !       /        \     /        \
    !      /   tbotl  \   /  tbotr   \
    !   \ /            \ /            \ /
    ! - -o----------- vbot ------------o- -
    !   / \            / \            / \
    !
    ! When coming out, it looks like this:
    !
    !   \ /            \ /            \ /
    ! _ _o_____________vtop____________o_ _
    !   / \            /|\            / \
    !      \          / | \          /
    !       \  ttopl /  |  \ ttopr  /
    !        \      /   |   \      /
    !         \    / t1 | t2 \    /
    !          \  /     |     \  /
    !       _ _ vi ___ vk ____ vj _ _
    !          /  \     |     /  \
    !         /    \ t3 | t4 /    \
    !        /      \   |   /      \
    !       /        \  |  /        \
    !      /   tbotl  \ | /  tbotr   \
    !   \ /            \|/            \ /
    ! _ _o_____________vbot____________o_ _
    !   / \            / \            / \

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi, vj
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'split_line'
    LOGICAL                                       :: are_connected_ij, are_connected_ji
    INTEGER                                       :: ci, cj, iti, n, n1, n2, n3, nf
    REAL(dp), DIMENSION(2)                        :: p, pa, pb
    INTEGER                                       :: t1, t2, t3, t4, ti, tit, tib, titl, titr, tibl, tibr, via, vib, vit, vk

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Safety
    ! =========

    ! Check IF vi and vj are even connected
    are_connected_ij = .FALSE.
    are_connected_ji = .FALSE.
    DO ci = 1, mesh%nC( vi)
      IF (mesh%C( vi,ci) == vj) are_connected_ij = .TRUE.
    END DO
    DO cj = 1, mesh%nC( vj)
      IF (mesh%C( vj,cj) == vi) are_connected_ji = .TRUE.
    END DO
    IF (are_connected_ij .AND. are_connected_ji) THEN
      ! Both vi and vj list each other as neighbours; all is well
    ELSEIF (are_connected_ij) THEN
      ! vi lists vj as a neighbour, not not the other way round - mesh inconsistency!
      CALL crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vi, int_02 = vj)
    ELSEIF (are_connected_ji) THEN
      ! vj lists vi as a neighbour, not not the other way round - mesh inconsistency!
      CALL crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vj, int_02 = vi)
    ELSE
      ! Neither vi nor vj lists the other as a neighbour
      CALL crash('{int_01} and {int_02} are not connected!', int_01 = vi, int_02 = vj)
    END IF

    ! Check IF p_new actually lies on [vi,vj]
    pa = mesh%V( vi,:)
    pb = mesh%V( vj,:)
    IF (.NOT. lies_on_line_segment( pa, pb, p_new, mesh%tol_dist)) THEN
      CALL crash('p does not lie on [{int_01}-{int_02}!', int_01 = vi, int_02 = vj)
    END IF

    ! If [vi,vj] is a boundary segment, split that instead
    IF (is_boundary_segment( mesh, vi, vj)) THEN
      p = (pa + pb) / 2._dp
      CALL split_segment( mesh, vi, vj, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == All safety checks passes; split the triangles t_top and t_bot adjacent
    !    to line [vi,vj] into four new ones
    ! =======================================================================

    ! DENK DROM
    IF (do_debug) CALL warning('splitting line [{int_01}-{int_02}}]', int_01 = vi, int_02 = vj)

    ! == Find the local neighbourhood of vertices and triangles

    ! Find the triangles tit and tib that are above and below the line vi-vj,
    ! respectively (see diagram).

    tit = 0
    tib = 0

    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF     (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) THEN
          tit = ti
        ELSEIF (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi) THEN
          tib = ti
        END IF
      END DO
    END DO

    ! Safety
    IF (tit == 0 .OR. tib == 0) THEN
      CALL crash('couldnt find triangles adjacent to [{int_01}-{int_02}!', int_01 = vi, int_02 = vj)
    END IF

    ! Find the vertices vit and vib that are on the opposite corners of tit and
    ! tib, respectively (see diagram).

    vit = 0
    vib = 0

    DO n1 = 1, 3
      n2 = n1 + 1
      IF (n2 == 4) n2 = 1
      n3 = n2 + 1
      IF (n3 == 4) n3 = 1
      IF (mesh%Tri( tit,n1) == vi .AND. mesh%Tri( tit,n2) == vj) THEN
        vit = mesh%Tri( tit,n3)
      END IF
      IF (mesh%Tri( tib,n1) == vj .AND. mesh%Tri( tib,n2) == vi) THEN
        vib = mesh%Tri( tib,n3)
      END IF
    END DO

    ! Safety
    IF (vit == 0 .OR. vib == 0) THEN
      CALL crash('couldnt find vertices opposite from [{int_01}-{int_02}!', int_01 = vi, int_02 = vj)
    END IF

    ! Find the triangles titl, titr, tibl, and tibr that are adjacent to tit
    ! and tir, respectively (see diagram).

    titl = 0
    titr = 0
    tibl = 0
    tibr = 0

    DO n = 1,3
      IF (mesh%Tri( tit,n) == vj) THEN
        titl = mesh%TriC( tit,n)
      END IF
      IF (mesh%Tri( tit,n) == vi) THEN
        titr = mesh%TriC( tit,n)
      END IF
      IF (mesh%Tri( tib,n) == vj) THEN
        tibl = mesh%TriC( tib,n)
      END IF
      IF (mesh%Tri( tib,n) == vi) THEN
        tibr = mesh%TriC( tib,n)
      END IF
    END DO

    ! == V, Tri

    ! Create a new vertex vik at p_new
    mesh%nV = mesh%nV + 1
    vk = mesh%nV
    mesh%V( vk,:) = p_new

    ! Let triangle t1 be spanned by [vk,vit,vi]
    t1 = tit
    mesh%Tri( t1,:) = [vk, vit, vi]

    ! Let triangle t2 be spanned by [vk,vj,vit]
    t2 = tib
    mesh%Tri( t2,:) = [vk, vj, vit]

    ! Let triangle t3 be spanned by [vk,vi,vib]
    mesh%nTri = mesh%nTri + 1
    t3 = mesh%nTri
    mesh%Tri( t3,:) = [vk, vi, vib]

    ! Let triangle t4 be spanned by [vk,vib,vj]
    mesh%nTri = mesh%nTri + 1
    t4 = mesh%nTri
    mesh%Tri( t4,:) = [vk, vib, vj]

    ! DENK DROM
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t1)
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t2)
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t3)
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t4)

    ! == nC, C

    ! vi: connection to vj is replaced by vk
    DO ci = 1, mesh%nC( vi)
      IF (mesh%C( vi,ci) == vj) THEN
        mesh%C( vi,ci) = vk
        EXIT
      END IF
    END DO
    ! vj: connection to vi is replaced by vk
    DO ci = 1, mesh%nC( vj)
      IF (mesh%C( vj,ci) == vi) THEN
        mesh%C( vj,ci) = vk
        EXIT
      END IF
    END DO
    ! vit: connection to vk is added between those to vi and vj
    mesh%nC( vit) = mesh%nC( vit) + 1
    DO ci = 1, mesh%nC( vit)
      IF (mesh%C( vit,ci) == vi) THEN
        mesh%C( vit,:) = [mesh%C( vit,1:ci), vk, mesh%C( vit,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vib: connection to vk is added between those to vj and vi
    mesh%nC( vib) = mesh%nC( vib) + 1
    DO ci = 1, mesh%nC( vib)
      IF (mesh%C( vib,ci) == vj) THEN
        mesh%C( vib,:) = [mesh%C( vib,1:ci), vk, mesh%C( vib,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vik: [vi,vib,vj,vit]
    mesh%nC( vk    ) = 4
    mesh%C(  vk,1:4) = [vi, vib, vj, vit]

    ! == niTri, iTri

    ! vi: tit is replaced by t1, tib is replaced by t3
    DO iti = 1, mesh%niTri( vi)
      IF     (mesh%iTri( vi,iti) == tit) THEN
        mesh%iTri( vi,iti) = t1
      ELSEIF (mesh%iTri( vi,iti) == tib) THEN
        mesh%iTri( vi,iti) = t3
      END IF
    END DO
    ! vj: tit is replaced by t2, tib is replaced by t4
    DO iti = 1, mesh%niTri( vj)
      IF     (mesh%iTri( vj,iti) == tit) THEN
        mesh%iTri( vj,iti) = t2
      ELSEIF (mesh%iTri( vj,iti) == tib) THEN
        mesh%iTri( vj,iti) = t4
      END IF
    END DO
    ! vit: tit is replaced by [t1,t2]
    mesh%niTri( vit) = mesh%niTri( vit) + 1
    DO iti = 1, mesh%niTri( vit)
      IF (mesh%iTri( vit,iti) == tit) THEN
        mesh%iTri( vit,:) = [mesh%iTri( vit,1:iti-1), t1, t2, mesh%iTri( vit,iti+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vib: tib is replaced by [t4,t3]
    mesh%niTri( vib) = mesh%niTri( vib) + 1
    DO iti = 1, mesh%niTri( vib)
      IF (mesh%iTri( vib,iti) == tib) THEN
        mesh%iTri( vib,:) = [mesh%iTri( vib,1:iti-1), t4, t3, mesh%iTri( vib,iti+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vik: [t1,t3,t4,t2]
    mesh%niTri( vk    ) = 4
    mesh%iTri(  vk,1:4) = [t1, t3, t4, t2]

    ! == Boundary index

    mesh%VBI( vk) = 0

    ! == TriC

    ! t1: [titl,t3,t2]
    mesh%TriC( t1,:) = [titl, t3, t2]
    ! t2: [titr,t1,t4]
    mesh%TriC( t2,:) = [titr, t1, t4]
    ! t3: [tibl,t4,t1]
    mesh%TriC( t3,:) = [tibl, t4, t1]
    ! t4: [tibr,t2,t3]
    mesh%TriC( t4,:) = [tibr, t2, t3]
    ! titl: tit is replaced by t1
    IF (titl > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( titl,n) == tit) THEN
          mesh%TriC( titl,n) = t1
        END IF
      END DO
    END IF
    ! titr: tit is replaced by t2
    IF (titr > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( titr,n) == tit) THEN
          mesh%TriC( titr,n) = t2
        END IF
      END DO
    END IF
    ! tibl: tib is replaced by t3
    IF (tibl > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tibl,n) == tib) THEN
          mesh%TriC( tibl,n) = t3
        END IF
      END DO
    END IF
    ! tibr: tib is replaced by t4
    IF (tibr > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tibr,n) == tib) THEN
          mesh%TriC( tibr,n) = t4
        END IF
      END DO
    END IF

    ! == Tricc

    CALL update_triangle_circumcenter( mesh, t1)
    CALL update_triangle_circumcenter( mesh, t2)
    CALL update_triangle_circumcenter( mesh, t3)
    CALL update_triangle_circumcenter( mesh, t4)

    ! == Refinement data

    ! Add the four new triangles to the refinement map and stack
    CALL add_triangle_to_refinement_stack_last( mesh, t1)
    CALL add_triangle_to_refinement_stack_last( mesh, t2)
    CALL add_triangle_to_refinement_stack_last( mesh, t3)
    CALL add_triangle_to_refinement_stack_last( mesh, t4)

    ! Update triangle-line overlap ranges
    mesh%Tri_li( t1,:) = mesh%Tri_li( tit,:)
    mesh%Tri_li( t2,:) = mesh%Tri_li( tit,:)
    mesh%Tri_li( t3,:) = mesh%Tri_li( tib,:)
    mesh%Tri_li( t4,:) = mesh%Tri_li( tib,:)

    ! == Finished splitting the line. Iteratively flip any triangle pairs
    !    in the neighbourhood that do not meet the local Delaunay criterion.
    ! ======================================================================

    ! Start with the four newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Tri_flip_list( :,:) = 0
    nf = 0

    IF (titl > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t1, titl]
    END IF
    IF (titr > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t2, titr]
    END IF
    IF (tibl > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t3, tibl]
    END IF
    IF (tibr > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t4, tibr]
    END IF

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere
    CALL flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_line

  SUBROUTINE split_segment( mesh, vi, vj, p_new)
    ! Split the triangle ti adjacent to boundary segment [vi,vj] into two new ones
    !
    ! When going in, the local geometry looks like this:
    !
    !   \ /           \ /           \ /
    ! - -o----------- vic -----------o- -
    !   / \           / \           / \
    !      \         /   \         /
    !       \  tib  /     \  tia  /
    !        \     /   ti  \     /
    !         \   /         \   /
    !          \ /           \ /
    !      - - via --------- vib - -
    !
    ! (With either via=vi, vib=vj, or the other way round)
    !
    ! When coming out, it looks like this:
    !
    !   \ /           \ /           \ /
    ! - -o----------- vic -----------o- -
    !   / \           /|\           / \
    !      \         / | \         /
    !       \  tib  /  |  \  tia  /
    !        \     /   |   \     /
    !         \   / t2 | t2 \   /
    !          \ /     |     \ /
    !      - - via -- vik -- vib - -

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi, vj
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'split_segment'
    LOGICAL                                       :: are_connected_ij, are_connected_ji
    INTEGER                                       :: ci, cj, iti, n, n1, n2, n3, nf
    REAL(dp), DIMENSION(2)                        :: pa, pb
    INTEGER                                       :: t1, t2, ti, tia, tib, tic, tii, via, vib, vic, vik

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Safety
    ! =========

    ! Check IF vi and vj are even connected
    are_connected_ij = .FALSE.
    are_connected_ji = .FALSE.
    DO ci = 1, mesh%nC( vi)
      IF (mesh%C( vi,ci) == vj) are_connected_ij = .TRUE.
    END DO
    DO cj = 1, mesh%nC( vj)
      IF (mesh%C( vj,cj) == vi) are_connected_ji = .TRUE.
    END DO
    IF (are_connected_ij .AND. are_connected_ji) THEN
      ! Both vi and vj list each other as neighbours; all is well
    ELSEIF (are_connected_ij) THEN
      ! vi lists vj as a neighbour, not not the other way round - mesh inconsistency!
      CALL crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vi, int_02 = vj)
    ELSEIF (are_connected_ji) THEN
      ! vj lists vi as a neighbour, not not the other way round - mesh inconsistency!
      CALL crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vj, int_02 = vi)
    ELSE
      ! Neither vi nor vj lists the other as a neighbour
      CALL crash('{int_01} and {int_02} are not connected!', int_01 = vi, int_02 = vj)
    END IF

    ! Check IF p_new actually lies on [vi,vj]
    pa = mesh%V( vi,:)
    pb = mesh%V( vj,:)
    IF (.NOT. lies_on_line_segment( pa, pb, p_new, mesh%tol_dist)) THEN
      CALL crash('p does not lie on [{int_01}-{int_02}]!', int_01 = vi, int_02 = vj)
    END IF

    ! Check IF [vi,vj] is actually a boundary segment
    IF (.NOT. is_boundary_segment( mesh, vi, vj)) THEN
      CALL crash('[{int_01}-{int_02}] is not a boundary segment!', int_01 = vi, int_02 = vj)
    END IF

    ! == All safety checks passes; split the triangle ti adjacent to boundary
    !    segment [vi,vj] into two new ones
    ! =======================================================================

    ! DENK DROM
    IF (do_debug) CALL warning('splitting segment [{int_01}-{int_02}}]', int_01 = vi, int_02 = vj)

    ! == Find the local neighbourhood of vertices and triangles

    ! Find the triangle t1 adjacent to the boundary segment [vi,vj]

    ti = 0
    DO iti = 1, mesh%niTri( vi)
      tii = mesh%iTri( vi,iti)
      IF (ANY( mesh%Tri( tii,:) == vj)) THEN
        ti = tii
      END IF
    END DO

    IF (ti == 0) THEN
      CALL crash('couldnt find triangle containing {int_01} and {int_02}!', int_01 = vi, int_02 = vj)
    END IF

    ! Let ti be spanned by vertices [via,vib,vic], such that via and
    ! vib lie on the boundary, and vic in the interior. Note: either via = vi
    ! and vib = vj, or via = vj and vib = vi, but its easier from here on to
    ! have the three vertices ordered counter-clockwise.
    ! Let tia and tib be the triangles adjacent to ti across from via and vib,
    ! respectively. Note: it is possible that at most one of those doesnt exist
    ! (but not both!)

    via = 0
    vib = 0
    vic = 0

    tia = 0
    tib = 0
    tic = 0

    DO n1 = 1, 3
      n2 = n1 + 1
      IF (n2 == 4) n2 = 1
      n3 = n2 + 1
      IF (n3 == 4) n3 = 1
      IF ((mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) .OR. &
          (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi)) THEN
        via = mesh%Tri(  ti,n1)
        vib = mesh%Tri(  ti,n2)
        vic = mesh%Tri(  ti,n3)
        tia = mesh%TriC( ti,n1)
        tib = mesh%TriC( ti,n2)
        tic = mesh%TriC( ti,n3)
      END IF
    END DO

    ! Safety
    IF (via == 0 .OR. vib == 0 .OR. vic == 0) CALL crash('mesh%Tri doesnt make sense!')
    IF (tia == 0 .AND. tib == 0) CALL crash('triangle has only one neighbour!')
    IF (tic > 0) CALL crash('triangle doesnt appear to be a boundary triangle!')

    ! == V, Tri

    ! Create a new vertex vik at p_new
    mesh%nV = mesh%nV + 1
    vik = mesh%nV
    mesh%V( vik,:) = p_new

    ! Let triangle t1 be spanned by [via,vik,vic]
    t1 = ti
    mesh%Tri( t1,:) = [via, vik, vic]

    ! Let triangle t2 be spanned by [vik,vib,vic]
    mesh%nTri = mesh%nTri + 1
    t2 = mesh%nTri
    mesh%Tri( t2,:) = [vik,vib,vic]

    ! DENK DROM
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t1)
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t2)

    ! == nC, C

    ! via: vik replaces vib as the first connection
    mesh%C( via,1) = vik
    ! vib: vik reaplces via as the last connection
    mesh%C( vib, mesh%nC( vib)) = vik
    ! vic: vik comes in between via and vib
    mesh%nC( vic  ) = mesh%nC( vic) + 1
    DO ci = 1, mesh%nC( vic)
      IF (mesh%C( vic,ci) == via) THEN
        mesh%C(  vic,:) = [mesh%C( vic,1:ci), vik, mesh%C( vic,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO

    ! vik
    mesh%nC( vik) = 3
    mesh%C( vik,1:3) = [vib, vic, via]

    ! == niTri, iTri

    ! Inverse triangle list of via: ti is replaced by t1 as the first one
    mesh%iTri( via, 1) = t1
    ! Inverse triangle list of vib: ti is replaced by t2 as the last one
    mesh%iTri( vib, mesh%niTri( vib)) = t2
    ! Inverse triangle list of vic: replace ti by [t1,t2]
    mesh%niTri( vic  ) = mesh%niTri( vic) + 1
    DO iti = 1, mesh%niTri( vic)
      IF (mesh%iTri( vic,iti) == ti) THEN
        mesh%iTri(  vic,:) = [mesh%iTri( vic,1:iti-1), t1, t2, mesh%iTri( vic,iti+1:mesh%nC_mem-1)]
      END IF
    END DO
    ! Inverse triangle list of vik: [t2,t1]
    mesh%niTri( vik) = 2
    mesh%iTri( vik,1:2) = [t2,t1]

    ! == Boundary index

    IF     ((mesh%VBI( via) == 8 .OR. mesh%VBI( via) == 1 .OR. mesh%VBI( via) == 2) .AND. &
            (mesh%VBI( vib) == 8 .OR. mesh%VBI( vib) == 1 .OR. mesh%VBI( vib) == 2)) THEN
      ! North
      mesh%VBI( vik) = 1
    ELSEIF ((mesh%VBI( via) == 2 .OR. mesh%VBI( via) == 3 .OR. mesh%VBI( via) == 4) .AND. &
            (mesh%VBI( vib) == 2 .OR. mesh%VBI( vib) == 3 .OR. mesh%VBI( vib) == 4)) THEN
      ! East
      mesh%VBI( vik) = 3
    ELSEIF ((mesh%VBI( via) == 4 .OR. mesh%VBI( via) == 5 .OR. mesh%VBI( via) == 6) .AND. &
            (mesh%VBI( vib) == 4 .OR. mesh%VBI( vib) == 5 .OR. mesh%VBI( vib) == 6)) THEN
      ! South
      mesh%VBI( vik) = 5
    ELSEIF ((mesh%VBI( via) == 6 .OR. mesh%VBI( via) == 7 .OR. mesh%VBI( via) == 8) .AND. &
            (mesh%VBI( vib) == 6 .OR. mesh%VBI( vib) == 7 .OR. mesh%VBI( vib) == 8)) THEN
      ! West
      mesh%VBI( vik) = 7
    ELSE
      CALL crash('edge indices of via and vib dont make sense!')
    END IF

    ! == TriC

    ! Triangle t1 is adjacent to t2 and tib
    mesh%TriC( t1,:) = [t2,tib,0]
    ! Triangle t2 is adjacent to tia and t1
    mesh%TriC( t2,:) = [tia,t1,0]
    ! Triangle tia is now adjacent to t2 instead of ti
    IF (tia > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tia,n) == ti) THEN
          mesh%TriC( tia,n) = t2
        END IF
      END DO
    END IF
    ! Triangle tib is now adjacent to t1 instead of ti
    IF (tib > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tib,n) == ti) THEN
          mesh%TriC( tib,n) = t1
        END IF
      END DO
    END IF

    ! == Tricc

    CALL update_triangle_circumcenter( mesh, t1)
    CALL update_triangle_circumcenter( mesh, t2)

    ! == Refinement data

    ! Add the two new triangles to the refinement map and stack
    CALL add_triangle_to_refinement_stack_last( mesh, t1)
    CALL add_triangle_to_refinement_stack_last( mesh, t2)

    ! Update triangle-line overlap ranges
    mesh%Tri_li( t1,:) = mesh%Tri_li( ti,:)
    mesh%Tri_li( t2,:) = mesh%Tri_li( ti,:)

    ! == Finished splitting the segment. Iteratively flip any triangle pairs
    !    in the neighbourhood that do not meet the local Delaunay criterion.
    ! ======================================================================

    ! Start with the two newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Tri_flip_list( :,:) = 0
    nf = 0

    IF (tia > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t2, tia]
    END IF
    IF (tib > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t1, tib]
    END IF

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere
    CALL flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_segment

  SUBROUTINE move_vertex( mesh, vi, p)
    ! Move vertex vi of the mesh to point p

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'move_vertex'
    INTEGER                                       :: nf, iti, ti, n, tj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Move the vertex
    mesh%V( vi,:) = p

    ! Update surrounding triangle circumcentres
    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      CALL update_triangle_circumcenter( mesh, ti)
    END DO

    ! Update triangulation
    mesh%Tri_flip_list = 0
    nf = 0

    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      DO n = 1, 3
        tj = mesh%TriC( ti,n)
        IF (tj > 0) THEN
          nf = nf + 1
          mesh%Tri_flip_list( nf,:) = [ti,tj]
        END IF
      END DO
    END DO

    CALL flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE move_vertex

  SUBROUTINE flip_triangles_until_Delaunay( mesh, nf)
    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(INOUT)     :: nf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'flip_triangles_until_Delaunay'
    INTEGER                                       :: ti, tj
    LOGICAL                                       :: are_connected_ij, are_connected_ji
    INTEGER                                       :: n, vi, vj, vii

    ! Add routine to path
    CALL init_routine( routine_name)

    DO WHILE (nf > 0)

      ! Take the last triangle pair from the list
      ti = mesh%Tri_flip_list( nf,1)
      tj = mesh%Tri_flip_list( nf,2)
      nf = nf - 1

      ! Safety
      IF (ti == 0 .OR. tj == 0) THEN
        CALL crash('found ti=0 in mesh%Tri_flip_list!')
      END IF

      ! Check IF these two triangles are still connected (they might have
      ! become disconnected due to an earlier flip operation
      are_connected_ij = .FALSE.
      are_connected_ji = .FALSE.
      DO n = 1, 3
        IF (mesh%TriC( ti,n) == tj) are_connected_ij = .TRUE.
        IF (mesh%TriC( tj,n) == ti) are_connected_ji = .TRUE.
      END DO
      IF (.NOT. are_connected_ij .AND. .NOT. are_connected_ij) THEN
        ! These two triangles are no longer connected
        CYCLE
      ELSEIF ((are_connected_ij .AND. .NOT. are_connected_ji) .OR. (.NOT. are_connected_ij .AND. are_connected_ji)) THEN
        ! Safety
        CALL crash('inconsistency in TriC!')
      END IF

      ! If they do not meet the local Delaunay criterion, flip them, and add
      ! any new triangle pairs to the flip list
      IF (.NOT. are_Delaunay( mesh, ti, tj)) THEN
        ! Flip them
        CALL flip_triangle_pair( mesh, ti, tj, nf)
      END IF ! IF .NOT. are_Delaunay( mesh, ti, tj)

    END DO ! while (nf>0)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE flip_triangles_until_Delaunay

  SUBROUTINE flip_triangle_pair( mesh, ti, tj, nf)
    ! Flip the triangle pair ti-tj, supposedly because it doesn't meet the
    ! local Delaunay criterion
    !
    ! When going in, the local geometry looks like this:
    !
    !   \ /           \ /           \ /
    ! - -o----------- vic ---------- o- -
    !   / \           / \           / \
    !      \  tib    /   \   tia   /
    !       \       /     \       /
    !        \     /       \     /
    !         \   /   ti    \   /
    !          \ /           \ /
    !      - - via -------- vib - -
    !          / \           / \
    !         /   \         /   \
    !        /     \  tj   /     \
    !       /       \     /       \
    !      /   tjb   \   /   tja   \
    !   \ /           \ /           \ /
    ! - -o----------- vid ---------- o- -
    !   / \           / \           / \
    !
    ! When coming out, it looks like this:
    !
    !   \ /           \ /           \ /
    ! - -o----------- vic ---------- o- -
    !   / \           /|\           / \
    !      \  tib    / | \   tia   /
    !       \       /  |  \       /
    !        \     /   |   \     /
    !         \   / t1 | t2 \   /
    !          \ /     |     \ /
    !      - - via     |     vib - -
    !          / \     |     / \
    !         /   \    |    /   \
    !        /     \   |   /     \
    !       /       \  |  /       \
    !      /   tjb   \ | /   tja   \
    !   \ /           \|/           \ /
    ! - -o----------- vid ---------- o- -
    !   / \           / \           / \

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti,tj
    INTEGER,                    INTENT(INOUT)     :: nf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'flip_triangle_pair'
    LOGICAL                                       :: are_connected_ij, are_connected_ji
    INTEGER                                       :: n, vi, vj, vii
    LOGICAL                                       :: is_in_tj
    INTEGER                                       :: via, vib, vic, vid
    INTEGER                                       :: ci, iti
    INTEGER                                       :: li_min, li_max
    INTEGER                                       :: n1, n2, n3
    INTEGER                                       :: tia, tib, tja, tjb, t1, t2, tii
    LOGICAL                                       :: via_has_ti, via_has_tj
    LOGICAL                                       :: vib_has_ti, vib_has_tj
    LOGICAL                                       :: vic_has_ti, vic_has_tj
    LOGICAL                                       :: vid_has_ti, vid_has_tj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (ti == 0 .OR. tj == 0) THEN
      CALL crash('Found ti=0 in mesh%Tri_flip_list!')
    END IF

    ! Check IF these two triangles are connected
    are_connected_ij = .FALSE.
    are_connected_ji = .FALSE.
    DO n = 1, 3
      IF (mesh%TriC( ti,n) == tj) are_connected_ij = .TRUE.
      IF (mesh%TriC( tj,n) == ti) are_connected_ji = .TRUE.
    END DO
    IF (.NOT. are_connected_ij .AND. .NOT. are_connected_ij) THEN
      ! These two triangles are not connected
      CALL crash('{int_01} and {int_02} are not connected!', int_01 = ti, int_02 = tj)
    ELSEIF (are_connected_ij .AND. .NOT. are_connected_ji .OR. .NOT. are_connected_ij .AND. are_connected_ji) THEN
      ! One of them lists the other as a neighbour, but not vice versa
      CALL crash('inconsistency in TriC!')
    END IF

    ! Find the two vertices vi and vj that are shared by ti and tj

    vi = 0
    vj = 0

    DO n = 1, 3
      vii = mesh%Tri( ti,n)
      is_in_tj = .FALSE.
      DO n2 = 1, 3
        IF (mesh%Tri( tj,n2) == vii) THEN
          is_in_tj = .TRUE.
          EXIT
        END IF
      END DO
      IF (is_in_tj) THEN
        IF (vi == 0) THEN
          vi = vii
        ELSE
          vj = vii
        END IF
      END IF
    END DO

    ! Safety
    IF (vi == 0 .OR. vj == 0) THEN
      CALL crash('couldnt find two shared vertices!')
    END IF

    ! Find via,vib,vic,vid (see diagram)

    via = 0
    vib = 0
    vic = 0
    vid = 0

    DO n1 = 1, 3

      n2 = n1 + 1
      IF (n2 == 4) n2 = 1
      n3 = n2 + 1
      IF (n3 == 4) n3 = 1

      IF ((mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) .OR. &
          (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi)) THEN
        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)
      END IF

      IF ((mesh%Tri( tj,n1) == vi .AND. mesh%Tri( tj,n2) == vj) .OR. &
          (mesh%Tri( tj,n1) == vj .AND. mesh%Tri( tj,n2) == vi)) THEN
        vid = mesh%Tri( tj,n3)
      END IF

    END DO

    ! Safety
    IF (via == 0 .OR. vib == 0 .OR. vic == 0 .OR. vid == 0) THEN
      CALL crash('couldnt figure out local geometry!')
    END IF
    IF (via == vib .OR. via == vic .OR. via == vid .OR. &
                        vib == vic .OR. vib == vid .OR. &
                                        vic == vid) THEN
      CALL crash('found duplicate vertices!')
    END IF

    via_has_ti = .FALSE.
    via_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( via)
      IF     (mesh%iTri( via,iti) == ti) THEN
        via_has_ti = .TRUE.
      ELSEIF (mesh%iTri( via,iti) == tj) THEN
        via_has_tj = .TRUE.
      END IF
    END DO
    IF (.NOT. via_has_ti .OR. .NOT. via_has_tj) THEN
      CALL crash('inconsistent mesh geometry!')
    END IF

    vib_has_ti = .FALSE.
    vib_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( vib)
      IF     (mesh%iTri( vib,iti) == ti) THEN
        vib_has_ti = .TRUE.
      ELSEIF (mesh%iTri( vib,iti) == tj) THEN
        vib_has_tj = .TRUE.
      END IF
    END DO
    IF (.NOT. vib_has_ti .OR. .NOT. vib_has_tj) THEN
      CALL crash('inconsistent mesh geometry!')
    END IF

    vic_has_ti = .FALSE.
    vic_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( vic)
      IF     (mesh%iTri( vic,iti) == ti) THEN
        vic_has_ti = .TRUE.
      ELSEIF (mesh%iTri( vic,iti) == tj) THEN
        vic_has_tj = .TRUE.
      END IF
    END DO
    IF (.NOT. vic_has_ti .OR. vic_has_tj) THEN
      CALL crash('inconsistent mesh geometry!')
    END IF

    vid_has_ti = .FALSE.
    vid_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( vid)
      IF     (mesh%iTri( vid,iti) == ti) THEN
        vid_has_ti = .TRUE.
      ELSEIF (mesh%iTri( vid,iti) == tj) THEN
        vid_has_tj = .TRUE.
      END IF
    END DO
    IF (vid_has_ti .OR. .NOT. vid_has_tj) THEN
      CALL crash('inconsistent mesh geometry!')
    END IF

    ! Find triangles tia,tib,tja,tjb (see diagram)

    tia = 0
    DO iti = 1, mesh%niTri( vic)
      tii = mesh%iTri( vic,iti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF (mesh%Tri( tii,n1) == vic .AND. mesh%Tri( tii,n2) == vib) THEN
          tia = tii
          EXIT
        END IF
      END DO
      IF (tia > 0) EXIT
    END DO

    tib = 0
    DO iti = 1, mesh%niTri( via)
      tii = mesh%iTri( via,iti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF (mesh%Tri( tii,n1) == via .AND. mesh%Tri( tii,n2) == vic) THEN
          tib = tii
          EXIT
        END IF
      END DO
      IF (tib > 0) EXIT
    END DO

    tja = 0
    DO iti = 1, mesh%niTri( vib)
      tii = mesh%iTri( vib,iti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF (mesh%Tri( tii,n1) == vib .AND. mesh%Tri( tii,n2) == vid) THEN
          tja = tii
          EXIT
        END IF
      END DO
      IF (tja > 0) EXIT
    END DO

    tjb = 0
    DO iti = 1, mesh%niTri( vid)
      tii = mesh%iTri( vid,iti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF (mesh%Tri( tii,n1) == vid .AND. mesh%Tri( tii,n2) == via) THEN
          tjb = tii
          EXIT
        END IF
      END DO
      IF (tjb > 0) EXIT
    END DO

    ! Safety
    IF (tia > 0 .AND. tib > 0) THEN
      IF (tia == tib) THEN
        CALL write_mesh_to_text_file( mesh, 'crashmesh.txt')
        CALL crash('found duplicate triangles!')
      END IF
    END IF
    IF (tia > 0 .AND. tja > 0) THEN
      IF (tia == tja) THEN
        CALL write_mesh_to_text_file( mesh, 'crashmesh.txt')
        CALL crash('found duplicate triangles!')
      END IF
    END IF
    IF (tia > 0 .AND. tjb > 0) THEN
      IF (tia == tjb) THEN
        CALL write_mesh_to_text_file( mesh, 'crashmesh.txt')
        CALL crash('found duplicate triangles!')
      END IF
    END IF
    IF (tib > 0 .AND. tja > 0) THEN
      IF (tib == tja) THEN
        CALL write_mesh_to_text_file( mesh, 'crashmesh.txt')
        CALL crash('found duplicate triangles!')
      END IF
    END IF
    IF (tib > 0 .AND. tjb > 0) THEN
      IF (tib == tjb) THEN
        CALL write_mesh_to_text_file( mesh, 'crashmesh.txt')
        CALL crash('found duplicate triangles!')
      END IF
    END IF
    IF (tja > 0 .AND. tjb > 0) THEN
      IF (tja == tjb) THEN
        CALL write_mesh_to_text_file( mesh, 'crashmesh.txt')
        CALL crash('found duplicate triangles!')
      END IF
    END IF

    ! == All safety checks passes; flip the triangle pair ti-tj
    ! =========================================================

    ! DENK DROM
    IF (do_debug) CALL warning('flipping triangles [{int_01}-{int_02}}]', int_01 = ti, int_02 = tj)

    ! == V, Tri

    ! Let triangle t1 be spanned by [via, vid, vic]
    t1 = ti
    mesh%Tri( t1,:) = [via, vid, vic]

    ! Let triangle t2 be spanned by [vib, vic, vid]
    t2 = tj
    mesh%Tri( t2,:) = [vib, vic, vid]

    ! DENK DROM
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t1)
    IF (do_debug) CALL check_if_triangle_already_exists( mesh, t2)

    ! == nC, C

    ! via: connection to vib is removed
    DO ci = 1, mesh%nC( via)
      IF (mesh%C( via,ci) == vib) THEN
        mesh%C( via,:) = [mesh%C( via,1:ci-1), mesh%C( via,ci+1:mesh%nC_mem), 0]
        mesh%nC( via) = mesh%nC( via) - 1
        EXIT
      END IF
    END DO
    ! vib: connection to via is removed
    DO ci = 1, mesh%nC( vib)
      IF (mesh%C( vib,ci) == via) THEN
        mesh%C( vib,:) = [mesh%C( vib,1:ci-1), mesh%C( vib,ci+1:mesh%nC_mem), 0]
        mesh%nC( vib) = mesh%nC( vib) - 1
        EXIT
      END IF
    END DO
    ! vic: a connection to vid is added between via and vib
    mesh%nC( vic) = mesh%nC( vic) + 1
    DO ci = 1, mesh%nC( vic)
      IF (mesh%C( vic,ci) == via) THEN
        mesh%C( vic,:) = [mesh%C( vic,1:ci), vid, mesh%C( vic,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vid: a connection to vic is added between vib and via
    mesh%nC( vid) = mesh%nC( vid) + 1
    DO ci = 1, mesh%nC( vid)
      IF (mesh%C( vid,ci) == vib) THEN
        mesh%C( vid,:) = [mesh%C( vid,1:ci), vic, mesh%C( vid,ci+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO

    ! == niTri, iTri

    ! via: tj,ti are replaced by t1
    DO iti = 1, mesh%niTri( via)
      IF (mesh%iTri( via,iti) == tj) THEN
        mesh%iTri( via,:) = [mesh%iTri( via,1:iti-1), t1, mesh%iTri( via,iti+2:mesh%nC_mem), 0]
        mesh%niTri( via) = mesh%niTri( via) - 1
        EXIT
      END IF
    END DO
    ! vib: ti,tj are replaced by t2
    DO iti = 1, mesh%niTri( vib)
      IF (mesh%iTri( vib,iti) == ti) THEN
        mesh%iTri( vib,:) = [mesh%iTri( vib,1:iti-1), t2, mesh%iTri( vib,iti+2:mesh%nC_mem), 0]
        mesh%niTri( vib) = mesh%niTri( vib) - 1
        EXIT
      END IF
    END DO
    ! vic: ti is replaced by t1,t2
    mesh%niTri( vic) = mesh%niTri( vic) + 1
    DO iti = 1, mesh%niTri( vic)
      IF (mesh%iTri( vic,iti) == ti) THEN
        mesh%iTri( vic,:) = [mesh%iTri( vic,1:iti-1), t1, t2, mesh%iTri( vic,iti+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO
    ! vid: tj is replaced by t2,t1
    mesh%niTri( vid) = mesh%niTri( vid) + 1
    DO iti = 1, mesh%niTri( vid)
      IF (mesh%iTri( vid,iti) == tj) THEN
        mesh%iTri( vid,:) = [mesh%iTri( vid,1:iti-1), t2, t1, mesh%iTri( vid,iti+1:mesh%nC_mem-1)]
        EXIT
      END IF
    END DO

    ! == Boundary index

    ! No changes.

    ! == TriC

    ! tia: ti is replaced by t2
    IF (tia > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tia,n) == ti) THEN
          mesh%TriC( tia,n) = t2
          EXIT
        END IF
      END DO
    END IF
    ! tib: ti is replaced by t1
    IF (tib > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tib,n) == ti) THEN
          mesh%TriC( tib,n) = t1
          EXIT
        END IF
      END DO
    END IF
    ! tja: tj is replaced by t2
    IF (tja > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tja,n) == tj) THEN
          mesh%TriC( tja,n) = t2
          EXIT
        END IF
      END DO
    END IF
    ! tjb: tj is replaced by t1
    IF (tjb > 0) THEN
      DO n = 1, 3
        IF (mesh%TriC( tjb,n) == tj) THEN
          mesh%TriC( tjb,n) = t1
          EXIT
        END IF
      END DO
    END IF
    ! t1 is adjacent to [t2, tib, tjb]
    mesh%TriC( t1,:) = [t2, tib, tjb]
    ! t2 is adjacent to [t1, tja, tia]
    mesh%TriC( t2,:) = [t1, tja, tia]

    ! == Tricc

    CALL update_triangle_circumcenter( mesh, t1)
    CALL update_triangle_circumcenter( mesh, t2)

    ! == Refinement data

    ! Add the two new triangles to the refinement map and stack
    CALL add_triangle_to_refinement_stack_first( mesh, t1)
    CALL add_triangle_to_refinement_stack_first( mesh, t2)

    ! Update triangle-line overlap ranges
    li_min = MIN( mesh%Tri_li( ti,1), mesh%Tri_li( tj,1))
    li_max = MAX( mesh%Tri_li( ti,2), mesh%Tri_li( tj,2))
    mesh%Tri_li( t1,:) = [li_min, li_max]
    mesh%Tri_li( t2,:) = [li_min, li_max]

    ! Add the four new pairs to the flip list
    IF (tia > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tia, t2]
    END IF
    IF (tib > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tib, t1]
    END IF
    IF (tja > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tja, t2]
    END IF
    IF (tjb > 0) THEN
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tjb, t1]
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE flip_triangle_pair

  FUNCTION are_Delaunay( mesh, ti, tj) RESULT( isso)
    ! Check if triangle pair ti-tj meets the local Delaunay criterion
    !
    ! The local geometry looks like this:
    !
    !       vic
    !       / \
    !      /   \
    !     / ti  \
    !    /       \
    !  via ----- vib
    !    \       /
    !     \ tj  /
    !      \   /
    !       \ /
    !       vid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti,tj
    LOGICAL                                       :: isso

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'are_Delaunay'
    LOGICAL                                       :: are_connected_ij, are_connected_ji
    INTEGER                                       :: n, vi, vj, vii, n1, n2, n3, iti
    LOGICAL                                       :: is_in_tj
    INTEGER                                       :: via, vib, vic, vid
    LOGICAL                                       :: via_has_ti, via_has_tj
    LOGICAL                                       :: vib_has_ti, vib_has_tj
    LOGICAL                                       :: vic_has_ti, vic_has_tj
    LOGICAL                                       :: vid_has_ti, vid_has_tj
    REAL(dp), DIMENSION(2)                        :: va, vb, vc, vd, cci, ccj
    REAL(dp)                                      :: ccri, ccrj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (ti == 0 .OR. tj == 0) THEN
      CALL crash('Found ti=0 in mesh%Tri_flip_list!')
    END IF

    ! Check if these two triangles are connected
    are_connected_ij = .FALSE.
    are_connected_ji = .FALSE.
    DO n = 1, 3
      IF (mesh%TriC( ti,n) == tj) are_connected_ij = .TRUE.
      IF (mesh%TriC( tj,n) == ti) are_connected_ji = .TRUE.
    END DO
    IF (.NOT. are_connected_ij .AND. .NOT. are_connected_ij) THEN
      ! These two triangles are not connected
      CALL crash('{int_01} and {int_02} are not connected!', int_01 = ti, int_02 = tj)
    ELSEIF (are_connected_ij .AND. .NOT. are_connected_ji .OR. .NOT. are_connected_ij .AND. are_connected_ji) THEN
      ! One of them lists the other as a neighbour, but not vice versa
      CALL crash('inconsistency in TriC!')
    END IF

    ! Find the two vertices vi and vj that are shared by ti and tj

    vi = 0
    vj = 0

    DO n = 1, 3
      vii = mesh%Tri( ti,n)
      is_in_tj = .FALSE.
      DO n2 = 1, 3
        IF (mesh%Tri( tj,n2) == vii) THEN
          is_in_tj = .TRUE.
          EXIT
        END IF
      END DO
      IF (is_in_tj) THEN
        IF (vi == 0) THEN
          vi = vii
        ELSE
          vj = vii
        END IF
      END IF
    END DO

    ! Safety
    IF (vi == 0 .OR. vj == 0) THEN
      CALL crash('couldnt find two shared vertices!')
    END IF

    ! Find via,vib,vic,vid (see diagram)

    via = 0
    vib = 0
    vic = 0
    vid = 0

    DO n1 = 1, 3

      n2 = n1 + 1
      IF (n2 == 4) n2 = 1
      n3 = n2 + 1
      IF (n3 == 4) n3 = 1

      IF ((mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) .OR. &
          (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi)) THEN
        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)
      END IF

      IF ((mesh%Tri( tj,n1) == vi .AND. mesh%Tri( tj,n2) == vj) .OR. &
          (mesh%Tri( tj,n1) == vj .AND. mesh%Tri( tj,n2) == vi)) THEN
        vid = mesh%Tri( tj,n3)
      END IF

    END DO

    ! Safety
    IF (via == 0 .OR. vib == 0 .OR. vic == 0 .OR. vid == 0) THEN
      CALL crash('couldnt figure out local geometry!')
    END IF

    via_has_ti = .FALSE.
    via_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( via)
      IF     (mesh%iTri( via,iti) == ti) THEN
        via_has_ti = .TRUE.
      ELSEIF (mesh%iTri( via,iti) == tj) THEN
        via_has_tj = .TRUE.
      END IF
    END DO
    IF (.NOT. via_has_ti) CALL crash('inconsistent mesh geometry! (via doesnt have ti as an itriangle)')
    IF (.NOT. via_has_tj) CALL crash('inconsistent mesh geometry! (via doesnt have tj as an itriangle)')

    vib_has_ti = .FALSE.
    vib_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( vib)
      IF     (mesh%iTri( vib,iti) == ti) THEN
        vib_has_ti = .TRUE.
      ELSEIF (mesh%iTri( vib,iti) == tj) THEN
        vib_has_tj = .TRUE.
      END IF
    END DO
    IF (.NOT. vib_has_ti) CALL crash('inconsistent mesh geometry! (vib doesnt have ti as an itriangle)')
    IF (.NOT. vib_has_tj) CALL crash('inconsistent mesh geometry! (vib doesnt have tj as an itriangle)')

    vic_has_ti = .FALSE.
    vic_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( vic)
      IF     (mesh%iTri( vic,iti) == ti) THEN
        vic_has_ti = .TRUE.
      ELSEIF (mesh%iTri( vic,iti) == tj) THEN
        vic_has_tj = .TRUE.
      END IF
    END DO
    IF (.NOT. vic_has_ti) CALL crash('inconsistent mesh geometry! (vic doesnt have ti as an itriangle)')
    IF (      vic_has_tj) THEN
      CALL warning('ti = [{int_01}, {int_02}, {int_03}], tj = [{int_04}, {int_05}, {int_06}]', &
        int_01 = mesh%Tri( ti,1), int_02 = mesh%Tri( ti,2), int_03 = mesh%Tri( ti,3), &
        int_04 = mesh%Tri( tj,1), int_05 = mesh%Tri( tj,2), int_06 = mesh%Tri( tj,3))
    END IF
    IF (      vic_has_tj) CALL crash('inconsistent mesh geometry! (vic has tj as an itriangle)')

    vid_has_ti = .FALSE.
    vid_has_tj = .FALSE.
    DO iti = 1, mesh%niTri( vid)
      IF     (mesh%iTri( vid,iti) == ti) THEN
        vid_has_ti = .TRUE.
      ELSEIF (mesh%iTri( vid,iti) == tj) THEN
        vid_has_tj = .TRUE.
      END IF
    END DO
    IF (      vid_has_ti) CALL crash('inconsistent mesh geometry! (vid has ti as an itriangle)')
    IF (.NOT. vid_has_tj) CALL crash('inconsistent mesh geometry! (vid doesnt have tj as an itriangle)')

    ! Check if ti-tj meets the Delaunay criterion

    va = mesh%V( via,:)
    vb = mesh%V( vib,:)
    vc = mesh%V( vic,:)
    vd = mesh%V( vid,:)

    cci = mesh%Tricc( ti,:)
    ccj = mesh%Tricc( tj,:)

    ccri = NORM2( va - cci)
    ccrj = NORM2( va - ccj)

    isso = .TRUE.

    IF     (NORM2( vd - cci) < ccri) THEN
      ! vid lies inside the circumcircle of ti
      isso = .FALSE.
    ELSEIF (NORM2( vc - ccj) < ccrj) THEN
      ! vic lies inside the circumcircle of tj
      isso = .FALSE.
    END IF

    ! If the outer angle at via or vib is concave, don't flip.
    ! Check this by checking if via lies inside the triangle
    ! [vib,vic,vid], or the other way round.

    IF  (is_in_triangle( vb, vc, vd, va) .OR. &
         is_in_triangle( va, vd, vc, vb)) THEN
      isso = .TRUE.
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END FUNCTION are_Delaunay

END MODULE mesh_Delaunay
