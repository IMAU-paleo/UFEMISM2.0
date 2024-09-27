module mesh_Delaunay

  ! Routines used for updating a Delaunay triangulation by splitting a triangle, an edge, or a border edge.
  !
  ! All these routines are set up so that, if the mesh that goes in is self-consistent (i.e. there
  ! are no erroneous entries in the connectivity lists) and meets the Delaunay criterion everywhere,
  ! then so does the mesh that comes out.

! ===== Preamble =====
! ====================

  use mpi
  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging                  , only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use mesh_types                                             , only: type_mesh
  use math_utilities                                         , only: is_in_triangle, lies_on_line_segment, line_from_points, line_line_intersection, &
                                                                     perpendicular_bisector_from_line, encroaches_upon
  use mesh_utilities                                         , only: update_triangle_circumcenter, find_containing_triangle, is_border_edge, &
                                                                     encroaches_upon_any, check_mesh, add_triangle_to_refinement_stack_first, &
                                                                     add_triangle_to_refinement_stack_last, check_if_triangle_already_exists, &
                                                                     write_mesh_to_text_file

  implicit none

! ===== Global variables =====
! ============================

  logical :: do_debug = .false.

contains

! ===== subroutines =====
! =======================

  subroutine split_triangle( mesh, ti_in_guess, p_new)
    ! Add a vertex at p_new, splitting triangle ti into three new ones
    !
    ! Provide a guess ti_in_guess for which triangle we think contains p (can be wrong,
    ! but guessing near the correct one speeds up the code).
    !
    ! if p_new coincides with a line or border edge, split that instead.
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

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    integer,                    intent(in)        :: ti_in_guess
    real(dp), dimension(2),     intent(in)        :: p_new

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'split_triangle'
    logical                                       :: isso
    real(dp), dimension(2)                        :: p
    integer                                       :: ci, iti, n, nf, t1, t2, t3, ti, tia, tib, tic
    real(dp), dimension(2)                        :: va, vb, vc
    integer                                       :: vi, via, vib, vic, vik, vj
    integer                                       :: li_min, li_max

    ! Add routine to path
    call init_routine( routine_name)

    ! if p_new encroaches upon a border edge, split that instead
    call encroaches_upon_any( mesh, p_new, isso, vi, vj)
    if (isso) then
      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      call split_border_edge( mesh, vi, vj, p)
      call finalise_routine( routine_name)
      return
    end if

    ! == if p_new lies outside of the mesh domain, split a border edge instead.
    ! =========================================================================

    if (p_new( 1) < mesh%xmin) then
      ! p_new lies to the west of the mesh domain

      ! Safety
      if (p_new( 2) < mesh%ymin .or. p_new( 2) > mesh%ymax) then
        call crash('p_new lies way outside mesh domain!')
      end if

      ! Find the two vertices vi,vj of the border edge that must be split.
      vi = 1
      vj = mesh%C( vi, mesh%nC( vi))
      do while (mesh%V( vj,2) < p_new( 2))
        vi = vj
        vj = mesh%C( vi, mesh%nC( vi))
      end do

      ! Safety
      if (.not. (p_new( 2) >= mesh%V( vi,2) .and. p_new( 2) <= mesh%V( vj,2))) then
        call crash('couldnt find border edge to split!')
      end if

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      call split_border_edge( mesh, vi, vj, p)

      call finalise_routine( routine_name)
      return

    end if ! if (p_new( 1) < mesh%xmin)

    if (p_new( 1) > mesh%xmax) then
      ! p_new lies to the east of the mesh domain

      ! Safety
      if (p_new( 2) < mesh%ymin .or. p_new( 2) > mesh%ymax) then
        call crash('p_new lies way outside mesh domain!')
      end if

      ! Find the two vertices vi,vj of the border edge that must be split.
      vi = 2
      vj = mesh%C( vi,1)
      do while (mesh%V( vj,2) < p_new( 2))
        vi = vj
        vj = mesh%C( vi, 1)
      end do

      ! Safety
      if (.not. (p_new( 2) >= mesh%V( vi,2) .and. p_new( 2) <= mesh%V( vj,2))) then
        call crash('couldnt find border edge to split!')
      end if

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      call split_border_edge( mesh, vi, vj, p)

      call finalise_routine( routine_name)
      return

    end if ! if (p_new( 1) > mesh%xmax)

    if (p_new( 2) < mesh%ymin) then
      ! p_new lies to the south of the mesh domain

      ! Safety
      if (p_new( 1) < mesh%xmin .or. p_new( 1) > mesh%xmax) then
        call crash('p_new lies way outside mesh domain!')
      end if

      ! Find the two vertices vi,vj of the border edge that must be split.
      vi = 1
      vj = mesh%C( vi, 1)
      do while (mesh%V( vj,1) < p_new( 1))
        vi = vj
        vj = mesh%C( vi, 1)
      end do

      ! Safety
      if (.not. (p_new( 1) >= mesh%V( vi,1) .and. p_new( 1) <= mesh%V( vj,1))) then
        call crash('couldnt find border edge to split!')
      end if

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      call split_border_edge( mesh, vi, vj, p)

      call finalise_routine( routine_name)
      return

    end if ! if (p_new( 2) < mesh%ymin)

    if (p_new( 2) > mesh%ymax) then
      ! p_new lies to the north of the mesh domain

      ! Safety
      if (p_new( 1) < mesh%xmin .or. p_new( 1) > mesh%xmax) then
        call crash('p_new lies way outside mesh domain!')
      end if

      ! Find the two vertices vi,vj of the border edge that must be split.
      vi = 1
      vj = mesh%C( vi, mesh%nC( vi))
      do while (mesh%V( vj,1) < p_new( 1))
        vi = vj
        vj = mesh%C( vi, mesh%nC( vi))
      end do

      ! Safety
      if (.not. (p_new( 1) >= mesh%V( vi,1) .and. p_new( 1) <= mesh%V( vj,1))) then
        call crash('couldnt find border edge to split!')
      end if

      p = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
      call split_border_edge( mesh, vi, vj, p)

      call finalise_routine( routine_name)
      return

    end if ! if (p_new( 2) > mesh%ymax)

    ! == Find the triangle containing p_new.
    ! ======================================

    ! Find the triangle containing p_new
    ti = ti_in_guess
    call find_containing_triangle( mesh, p_new, ti)

    ! The indices of the three vertices [a,b,c] spanning ti
    via = mesh%Tri( ti,1)
    vib = mesh%Tri( ti,2)
    vic = mesh%Tri( ti,3)

    ! The coordinates of the three vertices [a,b,c] spanning ti
    va  = mesh%V( via,:)
    vb  = mesh%V( vib,:)
    vc  = mesh%V( vic,:)

    ! Indice of the triangles neighbouring ti
    tia = mesh%TriC( ti,1)
    tib = mesh%TriC( ti,2)
    tic = mesh%TriC( ti,3)

    ! == Check if p_new encroaches upon a border edge. if so, split it.
    ! ======================================================================

    if (is_border_edge( mesh, via, vib) .and. encroaches_upon( va, vb, p_new, mesh%tol_dist)) then
      p = (mesh%V( via,:) + mesh%V( vib,:)) / 2._dp
      call split_border_edge( mesh, via, vib, p)
      call finalise_routine( routine_name)
      return
    end if

    if (is_border_edge( mesh, vib, vic) .and. encroaches_upon( va, vb, p_new, mesh%tol_dist)) then
      p = (mesh%V( vib,:) + mesh%V( vic,:)) / 2._dp
      call split_border_edge( mesh, vib, vic, p)
      call finalise_routine( routine_name)
      return
    end if

    if (is_border_edge( mesh, vic, via) .and. encroaches_upon( va, vb, p_new, mesh%tol_dist)) then
      p = (mesh%V( vic,:) + mesh%V( via,:)) / 2._dp
      call split_border_edge( mesh, vic, via, p)
      call finalise_routine( routine_name)
      return
    end if

    ! == if the new vertex is (almost) colinear with two vertices of
    !    the containing triangle, split the line between them instead.
    ! ================================================================

    if     (lies_on_line_segment( va, vb, p_new, mesh%tol_dist)) then
      call split_edge( mesh, via, vib, p_new)
      call finalise_routine( routine_name)
      return
    elseif (lies_on_line_segment( vb, vc, p_new, mesh%tol_dist)) then
      call split_edge( mesh, vib, vic, p_new)
      call finalise_routine( routine_name)
      return
    elseif (lies_on_line_segment( vc, va, p_new, mesh%tol_dist)) then
      call split_edge( mesh, vic, via, p_new)
      call finalise_routine( routine_name)
      return
    end if

    ! == All safety checks passes; split the triangle ti into three new ones.
    ! =======================================================================

    ! DENK DROM
    if (do_debug) call warning('splitting triangle {int_01}', int_01 = ti)

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
    if (do_debug) call check_if_triangle_already_exists( mesh, t1)
    if (do_debug) call check_if_triangle_already_exists( mesh, t2)
    if (do_debug) call check_if_triangle_already_exists( mesh, t3)

    ! == nC, C

    ! via: connection to vik is addded between vib and vic
    mesh%nC( via) = mesh%nC( via) + 1
    do ci = 1, mesh%nC( via)
      if (mesh%C( via,ci) == vib) then
        mesh%C( via,:) = [mesh%C( via,1:ci), vik, mesh%C( via,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vib: connection to vik is addded between vic and via
    mesh%nC( vib) = mesh%nC( vib) + 1
    do ci = 1, mesh%nC( vib)
      if (mesh%C( vib,ci) == vic) then
        mesh%C( vib,:) = [mesh%C( vib,1:ci), vik, mesh%C( vib,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vic: connection to vik is addded between via and vib
    mesh%nC( vic) = mesh%nC( vic) + 1
    do ci = 1, mesh%nC( vic)
      if (mesh%C( vic,ci) == via) then
        mesh%C( vic,:) = [mesh%C( vic,1:ci), vik, mesh%C( vic,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vik: connected to [via, vib, vic]
    mesh%nC( vik) = 3
    mesh%C( vik,1:3) = [via, vib, vic]

    ! == niTri, iTri

    ! via: ti is replaced by [t3,t1]
    mesh%niTri( via) = mesh%niTri( via) + 1
    do iti = 1, mesh%niTri( via)
      if (mesh%iTri( via,iti) == ti) then
        mesh%iTri( via,:) = [mesh%iTri( via,1:iti-1), t3, t1, mesh%iTri( via,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vib: ti is replaced by [t2,t3]
    mesh%niTri( vib) = mesh%niTri( vib) + 1
    do iti = 1, mesh%niTri( vib)
      if (mesh%iTri( vib,iti) == ti) then
        mesh%iTri( vib,:) = [mesh%iTri( vib,1:iti-1), t2, t3, mesh%iTri( vib,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vic: ti is replaced by [t1,t2]
    mesh%niTri( vic) = mesh%niTri( vic) + 1
    do iti = 1, mesh%niTri( vic)
      if (mesh%iTri( vic,iti) == ti) then
        mesh%iTri( vic,:) = [mesh%iTri( vic,1:iti-1), t1, t2, mesh%iTri( vic,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vik: surrounded by [t1, t3, t2]
    mesh%niTri( vik) = 3
    mesh%iTri( vik,1:3) = [t1, t3, t2]

    ! == Border index

    mesh%VBI( vik) = 0

    ! == TriC

    ! tia: ti is replaced by t2
    if (tia > 0) then
      do n = 1,3
        if (mesh%TriC( tia,n) == ti) then
          mesh%TriC( tia,n) = t2
          exit
        end if
      end do
    end if
    ! tib: ti is replaced by t1
    if (tib > 0) then
      do n = 1,3
        if (mesh%TriC( tib,n) == ti) then
          mesh%TriC( tib,n) = t1
          exit
        end if
      end do
    end if
    ! tic: ti is replaced by t3
    if (tic > 0) then
      do n = 1,3
        if (mesh%TriC( tic,n) == ti) then
          mesh%TriC( tic,n) = t3
          exit
        end if
      end do
    end if
    ! t1: [tib, t3, t2]
    mesh%TriC( t1,:) = [tib, t3, t2]
    ! t2: [tia, t1, t3]
    mesh%TriC( t2,:) = [tia, t1, t3]
    ! t3: [tic, t2, t1]
    mesh%TriC( t3,:) = [tic, t2, t1]

    ! == Tricc

    call update_triangle_circumcenter( mesh, t1)
    call update_triangle_circumcenter( mesh, t2)
    call update_triangle_circumcenter( mesh, t3)

    ! == Refinement data

    ! Add the three new triangles to the refinement map and stack
    call add_triangle_to_refinement_stack_last( mesh, t1)
    call add_triangle_to_refinement_stack_last( mesh, t2)
    call add_triangle_to_refinement_stack_last( mesh, t3)

    ! Update triangle-line overlap ranges
    li_min = mesh%Tri_li( ti,1)
    li_max = mesh%Tri_li( ti,2)
    mesh%Tri_li( t1,:) = [li_min, li_max]
    mesh%Tri_li( t2,:) = [li_min, li_max]
    mesh%Tri_li( t3,:) = [li_min, li_max]

    ! == Finished splitting the triangle. Iteratively flip any triangle pairs
    !    in the neighbourhood that do not meet the local Delaunay criterion.
    ! ======================================================================

    ! Start with the four newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Tri_flip_list( :,:) = 0
    nf = 0

    if (tia > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t2, tia]
    end if
    if (tib > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t1, tib]
    end if
    if (tic > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t3, tic]
    end if

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere
    call flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine split_triangle

  subroutine split_edge( mesh, vi, vj, p_new)
    ! Split the triangles t_left and t_right adjacent to line [vi,vj]
    ! into four new ones. if [vi,vj] is a border edge, split that instead.
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

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    integer,                    intent(in)        :: vi, vj
    real(dp), dimension(2),     intent(in)        :: p_new

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'split_edge'
    logical                                       :: are_connected_ij, are_connected_ji
    integer                                       :: ci, cj, iti, n, n1, n2, n3, nf
    real(dp), dimension(2)                        :: p, pa, pb
    integer                                       :: t1, t2, t3, t4, ti, tit, tib, titl, titr, tibl, tibr, vib, vit, vk
    integer                                       :: li_min, li_max

    ! Add routine to path
    call init_routine( routine_name)

    ! == Safety
    ! =========

    ! Check if vi and vj are even connected
    are_connected_ij = .false.
    are_connected_ji = .false.
    do ci = 1, mesh%nC( vi)
      if (mesh%C( vi,ci) == vj) are_connected_ij = .true.
    end do
    do cj = 1, mesh%nC( vj)
      if (mesh%C( vj,cj) == vi) are_connected_ji = .true.
    end do
    if (are_connected_ij .and. are_connected_ji) then
      ! Both vi and vj list each other as neighbours; all is well
    elseif (are_connected_ij) then
      ! vi lists vj as a neighbour, not not the other way round - mesh inconsistency!
      call crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vi, int_02 = vj)
    elseif (are_connected_ji) then
      ! vj lists vi as a neighbour, not not the other way round - mesh inconsistency!
      call crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vj, int_02 = vi)
    else
      ! Neither vi nor vj lists the other as a neighbour
      call crash('{int_01} and {int_02} are not connected!', int_01 = vi, int_02 = vj)
    end if

    ! Check if p_new actually lies on [vi,vj]
    pa = mesh%V( vi,:)
    pb = mesh%V( vj,:)
    if (.not. lies_on_line_segment( pa, pb, p_new, mesh%tol_dist)) then
      call crash('p does not lie on [{int_01}-{int_02}!', int_01 = vi, int_02 = vj)
    end if

    ! if [vi,vj] is a border edge, split that instead
    if (is_border_edge( mesh, vi, vj)) then
      p = (pa + pb) / 2._dp
      call split_border_edge( mesh, vi, vj, p)
      call finalise_routine( routine_name)
      return
    end if

    ! == All safety checks passes; split the triangles t_top and t_bot adjacent
    !    to line [vi,vj] into four new ones
    ! =======================================================================

    ! DENK DROM
    if (do_debug) call warning('splitting line [{int_01}-{int_02}}]', int_01 = vi, int_02 = vj)

    ! == Find the local neighbourhood of vertices and triangles

    ! Find the triangles tit and tib that are above and below the line vi-vj,
    ! respectively (see diagram).

    tit = 0
    tib = 0

    do iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        if     (mesh%Tri( ti,n1) == vi .and. mesh%Tri( ti,n2) == vj) then
          tit = ti
        elseif (mesh%Tri( ti,n1) == vj .and. mesh%Tri( ti,n2) == vi) then
          tib = ti
        end if
      end do
    end do

    ! Safety
    if (tit == 0 .or. tib == 0) then
      call crash('couldnt find triangles adjacent to [{int_01}-{int_02}!', int_01 = vi, int_02 = vj)
    end if

    ! Find the vertices vit and vib that are on the opposite corners of tit and
    ! tib, respectively (see diagram).

    vit = 0
    vib = 0

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1
      if (mesh%Tri( tit,n1) == vi .and. mesh%Tri( tit,n2) == vj) then
        vit = mesh%Tri( tit,n3)
      end if
      if (mesh%Tri( tib,n1) == vj .and. mesh%Tri( tib,n2) == vi) then
        vib = mesh%Tri( tib,n3)
      end if
    end do

    ! Safety
    if (vit == 0 .or. vib == 0) then
      call crash('couldnt find vertices opposite from [{int_01}-{int_02}!', int_01 = vi, int_02 = vj)
    end if

    ! Find the triangles titl, titr, tibl, and tibr that are adjacent to tit
    ! and tir, respectively (see diagram).

    titl = 0
    titr = 0
    tibl = 0
    tibr = 0

    do n = 1,3
      if (mesh%Tri( tit,n) == vj) then
        titl = mesh%TriC( tit,n)
      end if
      if (mesh%Tri( tit,n) == vi) then
        titr = mesh%TriC( tit,n)
      end if
      if (mesh%Tri( tib,n) == vj) then
        tibl = mesh%TriC( tib,n)
      end if
      if (mesh%Tri( tib,n) == vi) then
        tibr = mesh%TriC( tib,n)
      end if
    end do

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
    if (do_debug) call check_if_triangle_already_exists( mesh, t1)
    if (do_debug) call check_if_triangle_already_exists( mesh, t2)
    if (do_debug) call check_if_triangle_already_exists( mesh, t3)
    if (do_debug) call check_if_triangle_already_exists( mesh, t4)

    ! == nC, C

    ! vi: connection to vj is replaced by vk
    do ci = 1, mesh%nC( vi)
      if (mesh%C( vi,ci) == vj) then
        mesh%C( vi,ci) = vk
        exit
      end if
    end do
    ! vj: connection to vi is replaced by vk
    do ci = 1, mesh%nC( vj)
      if (mesh%C( vj,ci) == vi) then
        mesh%C( vj,ci) = vk
        exit
      end if
    end do
    ! vit: connection to vk is added between those to vi and vj
    mesh%nC( vit) = mesh%nC( vit) + 1
    do ci = 1, mesh%nC( vit)
      if (mesh%C( vit,ci) == vi) then
        mesh%C( vit,:) = [mesh%C( vit,1:ci), vk, mesh%C( vit,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vib: connection to vk is added between those to vj and vi
    mesh%nC( vib) = mesh%nC( vib) + 1
    do ci = 1, mesh%nC( vib)
      if (mesh%C( vib,ci) == vj) then
        mesh%C( vib,:) = [mesh%C( vib,1:ci), vk, mesh%C( vib,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vik: [vi,vib,vj,vit]
    mesh%nC( vk    ) = 4
    mesh%C(  vk,1:4) = [vi, vib, vj, vit]

    ! == niTri, iTri

    ! vi: tit is replaced by t1, tib is replaced by t3
    do iti = 1, mesh%niTri( vi)
      if     (mesh%iTri( vi,iti) == tit) then
        mesh%iTri( vi,iti) = t1
      elseif (mesh%iTri( vi,iti) == tib) then
        mesh%iTri( vi,iti) = t3
      end if
    end do
    ! vj: tit is replaced by t2, tib is replaced by t4
    do iti = 1, mesh%niTri( vj)
      if     (mesh%iTri( vj,iti) == tit) then
        mesh%iTri( vj,iti) = t2
      elseif (mesh%iTri( vj,iti) == tib) then
        mesh%iTri( vj,iti) = t4
      end if
    end do
    ! vit: tit is replaced by [t1,t2]
    mesh%niTri( vit) = mesh%niTri( vit) + 1
    do iti = 1, mesh%niTri( vit)
      if (mesh%iTri( vit,iti) == tit) then
        mesh%iTri( vit,:) = [mesh%iTri( vit,1:iti-1), t1, t2, mesh%iTri( vit,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vib: tib is replaced by [t4,t3]
    mesh%niTri( vib) = mesh%niTri( vib) + 1
    do iti = 1, mesh%niTri( vib)
      if (mesh%iTri( vib,iti) == tib) then
        mesh%iTri( vib,:) = [mesh%iTri( vib,1:iti-1), t4, t3, mesh%iTri( vib,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vik: [t1,t3,t4,t2]
    mesh%niTri( vk    ) = 4
    mesh%iTri(  vk,1:4) = [t1, t3, t4, t2]

    ! == Border index

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
    if (titl > 0) then
      do n = 1, 3
        if (mesh%TriC( titl,n) == tit) then
          mesh%TriC( titl,n) = t1
        end if
      end do
    end if
    ! titr: tit is replaced by t2
    if (titr > 0) then
      do n = 1, 3
        if (mesh%TriC( titr,n) == tit) then
          mesh%TriC( titr,n) = t2
        end if
      end do
    end if
    ! tibl: tib is replaced by t3
    if (tibl > 0) then
      do n = 1, 3
        if (mesh%TriC( tibl,n) == tib) then
          mesh%TriC( tibl,n) = t3
        end if
      end do
    end if
    ! tibr: tib is replaced by t4
    if (tibr > 0) then
      do n = 1, 3
        if (mesh%TriC( tibr,n) == tib) then
          mesh%TriC( tibr,n) = t4
        end if
      end do
    end if

    ! == Tricc

    call update_triangle_circumcenter( mesh, t1)
    call update_triangle_circumcenter( mesh, t2)
    call update_triangle_circumcenter( mesh, t3)
    call update_triangle_circumcenter( mesh, t4)

    ! == Refinement data

    ! Add the four new triangles to the refinement map and stack
    call add_triangle_to_refinement_stack_last( mesh, t1)
    call add_triangle_to_refinement_stack_last( mesh, t2)
    call add_triangle_to_refinement_stack_last( mesh, t3)
    call add_triangle_to_refinement_stack_last( mesh, t4)

    ! Update triangle-line overlap ranges
    li_min = min( mesh%Tri_li( tit,1), mesh%Tri_li( tib,1))
    li_max = max( mesh%Tri_li( tit,2), mesh%Tri_li( tib,2))
    mesh%Tri_li( t1,:) = [li_min, li_max]
    mesh%Tri_li( t2,:) = [li_min, li_max]
    mesh%Tri_li( t3,:) = [li_min, li_max]
    mesh%Tri_li( t4,:) = [li_min, li_max]

    ! == Finished splitting the line. Iteratively flip any triangle pairs
    !    in the neighbourhood that do not meet the local Delaunay criterion.
    ! ======================================================================

    ! Start with the four newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Tri_flip_list( :,:) = 0
    nf = 0

    if (titl > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t1, titl]
    end if
    if (titr > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t2, titr]
    end if
    if (tibl > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t3, tibl]
    end if
    if (tibr > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t4, tibr]
    end if

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere
    call flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine split_edge

  subroutine split_border_edge( mesh, vi, vj, p_new)
    ! Split the triangle ti adjacent to border edge [vi,vj] into two new ones
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
    type(type_mesh),            intent(inout)     :: mesh
    integer,                    intent(in)        :: vi, vj
    real(dp), dimension(2),     intent(in)        :: p_new

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'split_border_edge'
    logical                                       :: are_connected_ij, are_connected_ji
    integer                                       :: ci, cj, iti, n, n1, n2, n3, nf
    real(dp), dimension(2)                        :: pa, pb
    integer                                       :: t1, t2, ti, tia, tib, tic, tii, via, vib, vic, vik
    integer                                       :: li_min, li_max

    ! Add routine to path
    call init_routine( routine_name)

    ! == Safety
    ! =========

    ! Check if vi and vj are even connected
    are_connected_ij = .false.
    are_connected_ji = .false.
    do ci = 1, mesh%nC( vi)
      if (mesh%C( vi,ci) == vj) are_connected_ij = .true.
    end do
    do cj = 1, mesh%nC( vj)
      if (mesh%C( vj,cj) == vi) are_connected_ji = .true.
    end do
    if (are_connected_ij .and. are_connected_ji) then
      ! Both vi and vj list each other as neighbours; all is well
    elseif (are_connected_ij) then
      ! vi lists vj as a neighbour, not not the other way round - mesh inconsistency!
      call crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vi, int_02 = vj)
    elseif (are_connected_ji) then
      ! vj lists vi as a neighbour, not not the other way round - mesh inconsistency!
      call crash('mesh inconsistency: {int_01} lists {int_02} as a neighbour, but not the other way round!', int_01 = vj, int_02 = vi)
    else
      ! Neither vi nor vj lists the other as a neighbour
      call crash('{int_01} and {int_02} are not connected!', int_01 = vi, int_02 = vj)
    end if

    ! Check if p_new actually lies on [vi,vj]
    pa = mesh%V( vi,:)
    pb = mesh%V( vj,:)
    if (.not. lies_on_line_segment( pa, pb, p_new, mesh%tol_dist)) then
      call crash('p does not lie on [{int_01}-{int_02}]!', int_01 = vi, int_02 = vj)
    end if

    ! Check if [vi,vj] is actually a border edge
    if (.not. is_border_edge( mesh, vi, vj)) then
      call crash('[{int_01}-{int_02}] is not a border edge!', int_01 = vi, int_02 = vj)
    end if

    ! == All safety checks passes; split the triangle ti adjacent to border
    !    edge [vi,vj] into two new ones
    ! =======================================================================

    ! DENK DROM
    if (do_debug) call warning('splitting border edge [{int_01}-{int_02}}]', int_01 = vi, int_02 = vj)

    ! == Find the local neighbourhood of vertices and triangles

    ! Find the triangle t1 adjacent to the border edge [vi,vj]

    ti = 0
    do iti = 1, mesh%niTri( vi)
      tii = mesh%iTri( vi,iti)
      if (ANY( mesh%Tri( tii,:) == vj)) then
        ti = tii
      end if
    end do

    if (ti == 0) then
      call crash('couldnt find triangle containing {int_01} and {int_02}!', int_01 = vi, int_02 = vj)
    end if

    ! Let ti be spanned by vertices [via,vib,vic], such that via and
    ! vib lie on the border, and vic in the interior. Note: either via = vi
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

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1
      if ((mesh%Tri( ti,n1) == vi .and. mesh%Tri( ti,n2) == vj) .or. &
          (mesh%Tri( ti,n1) == vj .and. mesh%Tri( ti,n2) == vi)) then
        via = mesh%Tri(  ti,n1)
        vib = mesh%Tri(  ti,n2)
        vic = mesh%Tri(  ti,n3)
        tia = mesh%TriC( ti,n1)
        tib = mesh%TriC( ti,n2)
        tic = mesh%TriC( ti,n3)
      end if
    end do

    ! Safety
    if (via == 0 .or. vib == 0 .or. vic == 0) call crash('mesh%Tri doesnt make sense!')
    if (tia == 0 .and. tib == 0) call crash('triangle has only one neighbour!')
    if (tic > 0) call crash('triangle doesnt appear to be a border triangle!')

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
    if (do_debug) call check_if_triangle_already_exists( mesh, t1)
    if (do_debug) call check_if_triangle_already_exists( mesh, t2)

    ! == nC, C

    ! via: vik replaces vib as the first connection
    mesh%C( via,1) = vik
    ! vib: vik reaplces via as the last connection
    mesh%C( vib, mesh%nC( vib)) = vik
    ! vic: vik comes in between via and vib
    mesh%nC( vic  ) = mesh%nC( vic) + 1
    do ci = 1, mesh%nC( vic)
      if (mesh%C( vic,ci) == via) then
        mesh%C(  vic,:) = [mesh%C( vic,1:ci), vik, mesh%C( vic,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do

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
    do iti = 1, mesh%niTri( vic)
      if (mesh%iTri( vic,iti) == ti) then
        mesh%iTri(  vic,:) = [mesh%iTri( vic,1:iti-1), t1, t2, mesh%iTri( vic,iti+1:mesh%nC_mem-1)]
      end if
    end do
    ! Inverse triangle list of vik: [t2,t1]
    mesh%niTri( vik) = 2
    mesh%iTri( vik,1:2) = [t2,t1]

    ! == Border index

    if     ((mesh%VBI( via) == 8 .or. mesh%VBI( via) == 1 .or. mesh%VBI( via) == 2) .and. &
            (mesh%VBI( vib) == 8 .or. mesh%VBI( vib) == 1 .or. mesh%VBI( vib) == 2)) then
      ! North
      mesh%VBI( vik) = 1
    elseif ((mesh%VBI( via) == 2 .or. mesh%VBI( via) == 3 .or. mesh%VBI( via) == 4) .and. &
            (mesh%VBI( vib) == 2 .or. mesh%VBI( vib) == 3 .or. mesh%VBI( vib) == 4)) then
      ! East
      mesh%VBI( vik) = 3
    elseif ((mesh%VBI( via) == 4 .or. mesh%VBI( via) == 5 .or. mesh%VBI( via) == 6) .and. &
            (mesh%VBI( vib) == 4 .or. mesh%VBI( vib) == 5 .or. mesh%VBI( vib) == 6)) then
      ! South
      mesh%VBI( vik) = 5
    elseif ((mesh%VBI( via) == 6 .or. mesh%VBI( via) == 7 .or. mesh%VBI( via) == 8) .and. &
            (mesh%VBI( vib) == 6 .or. mesh%VBI( vib) == 7 .or. mesh%VBI( vib) == 8)) then
      ! West
      mesh%VBI( vik) = 7
    else
      call crash('edge indices of via and vib dont make sense!')
    end if

    ! == TriC

    ! Triangle t1 is adjacent to t2 and tib
    mesh%TriC( t1,:) = [t2,tib,0]
    ! Triangle t2 is adjacent to tia and t1
    mesh%TriC( t2,:) = [tia,t1,0]
    ! Triangle tia is now adjacent to t2 instead of ti
    if (tia > 0) then
      do n = 1, 3
        if (mesh%TriC( tia,n) == ti) then
          mesh%TriC( tia,n) = t2
        end if
      end do
    end if
    ! Triangle tib is now adjacent to t1 instead of ti
    if (tib > 0) then
      do n = 1, 3
        if (mesh%TriC( tib,n) == ti) then
          mesh%TriC( tib,n) = t1
        end if
      end do
    end if

    ! == Tricc

    call update_triangle_circumcenter( mesh, t1)
    call update_triangle_circumcenter( mesh, t2)

    ! == Refinement data

    ! Add the two new triangles to the refinement map and stack
    call add_triangle_to_refinement_stack_last( mesh, t1)
    call add_triangle_to_refinement_stack_last( mesh, t2)

    ! Update triangle-line overlap ranges
    li_min = mesh%Tri_li( ti,1)
    li_max = mesh%Tri_li( ti,2)
    mesh%Tri_li( t1,:) = [li_min, li_max]
    mesh%Tri_li( t2,:) = [li_min, li_max]

    ! == Finished splitting the border edge. Iteratively flip any triangle pairs
    !    in the neighbourhood that do not meet the local Delaunay criterion.
    ! ======================================================================

    ! Start with the two newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Tri_flip_list( :,:) = 0
    nf = 0

    if (tia > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t2, tia]
    end if
    if (tib > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [t1, tib]
    end if

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere
    call flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine split_border_edge

  subroutine move_vertex( mesh, vi, p)
    ! Move vertex vi of the mesh to point p

    IMPLICIT NONE

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    integer,                    intent(in)        :: vi
    real(dp), dimension(2),     intent(in)        :: p

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'move_vertex'
    integer                                       :: nf, iti, ti, n, tj

    ! Add routine to path
    call init_routine( routine_name)

    ! Move the vertex
    mesh%V( vi,:) = p

    ! Update surrounding triangle circumcentres
    do iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      call update_triangle_circumcenter( mesh, ti)
    end do

    ! Update triangulation
    mesh%Tri_flip_list = 0
    nf = 0

    do iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        if (tj > 0) then
          nf = nf + 1
          mesh%Tri_flip_list( nf,:) = [ti,tj]
        end if
      end do
    end do

    call flip_triangles_until_Delaunay( mesh, nf)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine move_vertex

  subroutine flip_triangles_until_Delaunay( mesh, nf)
    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere

    IMPLICIT NONE

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    integer,                    intent(inout)     :: nf

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'flip_triangles_until_Delaunay'
    integer                                       :: ti, tj
    logical                                       :: are_connected_ij, are_connected_ji
    integer                                       :: n

    ! Add routine to path
    call init_routine( routine_name)

    do while (nf > 0)

      ! Take the last triangle pair from the list
      ti = mesh%Tri_flip_list( nf,1)
      tj = mesh%Tri_flip_list( nf,2)
      nf = nf - 1

      ! Safety
      if (ti == 0 .or. tj == 0) then
        call crash('found ti=0 in mesh%Tri_flip_list!')
      end if

      ! Check if these two triangles are still connected (they might have
      ! become disconnected due to an earlier flip operation
      are_connected_ij = .false.
      are_connected_ji = .false.
      do n = 1, 3
        if (mesh%TriC( ti,n) == tj) are_connected_ij = .true.
        if (mesh%TriC( tj,n) == ti) are_connected_ji = .true.
      end do
      if (.not. are_connected_ij .and. .not. are_connected_ij) then
        ! These two triangles are no longer connected
        CYCLE
      elseif ((are_connected_ij .and. .not. are_connected_ji) .or. (.not. are_connected_ij .and. are_connected_ji)) then
        ! Safety
        call crash('inconsistency in TriC!')
      end if

      ! if they do not meet the local Delaunay criterion, flip them, and add
      ! any new triangle pairs to the flip list
      if (.not. are_Delaunay( mesh, ti, tj)) then
        ! Flip them
        call flip_triangle_pair( mesh, ti, tj, nf)
      end if ! if .not. are_Delaunay( mesh, ti, tj)

    end do ! while (nf>0)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_triangles_until_Delaunay

  subroutine flip_triangle_pair( mesh, ti, tj, nf)
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
    type(type_mesh),            intent(inout)     :: mesh
    integer,                    intent(in)        :: ti,tj
    integer,                    intent(inout)     :: nf

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'flip_triangle_pair'
    logical                                       :: are_connected_ij, are_connected_ji
    integer                                       :: n, vi, vj, vii
    logical                                       :: is_in_tj
    integer                                       :: via, vib, vic, vid
    integer                                       :: ci, iti
    integer                                       :: n1, n2, n3
    integer                                       :: tia, tib, tja, tjb, t1, t2, tii
    logical                                       :: via_has_ti, via_has_tj
    logical                                       :: vib_has_ti, vib_has_tj
    logical                                       :: vic_has_ti, vic_has_tj
    logical                                       :: vid_has_ti, vid_has_tj
    integer                                       :: li_min, li_max
    logical                                       :: foundit

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (ti == 0 .or. tj == 0) then
      call crash('Found ti=0 in mesh%Tri_flip_list!')
    end if

    ! Check if these two triangles are connected
    are_connected_ij = .false.
    are_connected_ji = .false.
    do n = 1, 3
      if (mesh%TriC( ti,n) == tj) are_connected_ij = .true.
      if (mesh%TriC( tj,n) == ti) are_connected_ji = .true.
    end do
    if (.not. are_connected_ij .and. .not. are_connected_ij) then
      ! These two triangles are not connected
      call crash('{int_01} and {int_02} are not connected!', int_01 = ti, int_02 = tj)
    elseif (are_connected_ij .and. .not. are_connected_ji .or. .not. are_connected_ij .and. are_connected_ji) then
      ! One of them lists the other as a neighbour, but not vice versa
      call crash('inconsistency in TriC!')
    end if

    ! Find the two vertices vi and vj that are shared by ti and tj

    vi = 0
    vj = 0

    do n = 1, 3
      vii = mesh%Tri( ti,n)
      is_in_tj = .false.
      do n2 = 1, 3
        if (mesh%Tri( tj,n2) == vii) then
          is_in_tj = .true.
          exit
        end if
      end do
      if (is_in_tj) then
        if (vi == 0) then
          vi = vii
        else
          vj = vii
        end if
      end if
    end do

    ! Safety
    if (vi == 0 .or. vj == 0) then
      call crash('couldnt find two shared vertices!')
    end if

    ! Find via,vib,vic,vid (see diagram)

    via = 0
    vib = 0
    vic = 0
    vid = 0

    do n1 = 1, 3

      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      if ((mesh%Tri( ti,n1) == vi .and. mesh%Tri( ti,n2) == vj) .or. &
          (mesh%Tri( ti,n1) == vj .and. mesh%Tri( ti,n2) == vi)) then
        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)
      end if

      if ((mesh%Tri( tj,n1) == vi .and. mesh%Tri( tj,n2) == vj) .or. &
          (mesh%Tri( tj,n1) == vj .and. mesh%Tri( tj,n2) == vi)) then
        vid = mesh%Tri( tj,n3)
      end if

    end do

    ! Safety
    if (via == 0 .or. vib == 0 .or. vic == 0 .or. vid == 0) then
      call crash('couldnt figure out local geometry!')
    end if
    if (via == vib .or. via == vic .or. via == vid .or. &
                        vib == vic .or. vib == vid .or. &
                                        vic == vid) then
      call crash('found duplicate vertices!')
    end if

    via_has_ti = .false.
    via_has_tj = .false.
    do iti = 1, mesh%niTri( via)
      if     (mesh%iTri( via,iti) == ti) then
        via_has_ti = .true.
      elseif (mesh%iTri( via,iti) == tj) then
        via_has_tj = .true.
      end if
    end do
    if (.not. via_has_ti .or. .not. via_has_tj) then
      call crash('inconsistent mesh geometry!')
    end if

    vib_has_ti = .false.
    vib_has_tj = .false.
    do iti = 1, mesh%niTri( vib)
      if     (mesh%iTri( vib,iti) == ti) then
        vib_has_ti = .true.
      elseif (mesh%iTri( vib,iti) == tj) then
        vib_has_tj = .true.
      end if
    end do
    if (.not. vib_has_ti .or. .not. vib_has_tj) then
      call crash('inconsistent mesh geometry!')
    end if

    vic_has_ti = .false.
    vic_has_tj = .false.
    do iti = 1, mesh%niTri( vic)
      if     (mesh%iTri( vic,iti) == ti) then
        vic_has_ti = .true.
      elseif (mesh%iTri( vic,iti) == tj) then
        vic_has_tj = .true.
      end if
    end do
    if (.not. vic_has_ti .or. vic_has_tj) then
      call crash('inconsistent mesh geometry!')
    end if

    vid_has_ti = .false.
    vid_has_tj = .false.
    do iti = 1, mesh%niTri( vid)
      if     (mesh%iTri( vid,iti) == ti) then
        vid_has_ti = .true.
      elseif (mesh%iTri( vid,iti) == tj) then
        vid_has_tj = .true.
      end if
    end do
    if (vid_has_ti .or. .not. vid_has_tj) then
      call crash('inconsistent mesh geometry!')
    end if

    ! Find triangles tia,tib,tja,tjb (see diagram)

    tia = 0
    do iti = 1, mesh%niTri( vic)
      tii = mesh%iTri( vic,iti)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        if (mesh%Tri( tii,n1) == vic .and. mesh%Tri( tii,n2) == vib) then
          tia = tii
          exit
        end if
      end do
      if (tia > 0) exit
    end do

    tib = 0
    do iti = 1, mesh%niTri( via)
      tii = mesh%iTri( via,iti)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        if (mesh%Tri( tii,n1) == via .and. mesh%Tri( tii,n2) == vic) then
          tib = tii
          exit
        end if
      end do
      if (tib > 0) exit
    end do

    tja = 0
    do iti = 1, mesh%niTri( vib)
      tii = mesh%iTri( vib,iti)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        if (mesh%Tri( tii,n1) == vib .and. mesh%Tri( tii,n2) == vid) then
          tja = tii
          exit
        end if
      end do
      if (tja > 0) exit
    end do

    tjb = 0
    do iti = 1, mesh%niTri( vid)
      tii = mesh%iTri( vid,iti)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        if (mesh%Tri( tii,n1) == vid .and. mesh%Tri( tii,n2) == via) then
          tjb = tii
          exit
        end if
      end do
      if (tjb > 0) exit
    end do

    ! Safety
    if (tia > 0 .and. tib > 0) then
      if (tia == tib) then
        call write_mesh_to_text_file( mesh, 'crashmesh.txt')
        call crash('found duplicate triangles!')
      end if
    end if
    if (tia > 0 .and. tja > 0) then
      if (tia == tja) then
        call write_mesh_to_text_file( mesh, 'crashmesh.txt')
        call crash('found duplicate triangles!')
      end if
    end if
    if (tia > 0 .and. tjb > 0) then
      if (tia == tjb) then
        call write_mesh_to_text_file( mesh, 'crashmesh.txt')
        call crash('found duplicate triangles!')
      end if
    end if
    if (tib > 0 .and. tja > 0) then
      if (tib == tja) then
        call write_mesh_to_text_file( mesh, 'crashmesh.txt')
        call crash('found duplicate triangles!')
      end if
    end if
    if (tib > 0 .and. tjb > 0) then
      if (tib == tjb) then
        call write_mesh_to_text_file( mesh, 'crashmesh.txt')
        call crash('found duplicate triangles!')
      end if
    end if
    if (tja > 0 .and. tjb > 0) then
      if (tja == tjb) then
        call write_mesh_to_text_file( mesh, 'crashmesh.txt')
        call crash('found duplicate triangles!')
      end if
    end if

    ! == All safety checks passes; flip the triangle pair ti-tj
    ! =========================================================

    ! DENK DROM
    if (do_debug) call warning('flipping triangles [{int_01}-{int_02}}]', int_01 = ti, int_02 = tj)

    ! == V, Tri

    ! Let triangle t1 be spanned by [via, vid, vic]
    t1 = ti
    mesh%Tri( t1,:) = [via, vid, vic]

    ! Let triangle t2 be spanned by [vib, vic, vid]
    t2 = tj
    mesh%Tri( t2,:) = [vib, vic, vid]

    ! DENK DROM
    if (do_debug) call check_if_triangle_already_exists( mesh, t1)
    if (do_debug) call check_if_triangle_already_exists( mesh, t2)

    ! == nC, C

    ! via: connection to vib is removed
    do ci = 1, mesh%nC( via)
      if (mesh%C( via,ci) == vib) then
        mesh%C( via,:) = [mesh%C( via,1:ci-1), mesh%C( via,ci+1:mesh%nC_mem), 0]
        mesh%nC( via) = mesh%nC( via) - 1
        exit
      end if
    end do
    ! vib: connection to via is removed
    do ci = 1, mesh%nC( vib)
      if (mesh%C( vib,ci) == via) then
        mesh%C( vib,:) = [mesh%C( vib,1:ci-1), mesh%C( vib,ci+1:mesh%nC_mem), 0]
        mesh%nC( vib) = mesh%nC( vib) - 1
        exit
      end if
    end do
    ! vic: a connection to vid is added between via and vib
    mesh%nC( vic) = mesh%nC( vic) + 1
    do ci = 1, mesh%nC( vic)
      if (mesh%C( vic,ci) == via) then
        mesh%C( vic,:) = [mesh%C( vic,1:ci), vid, mesh%C( vic,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vid: a connection to vic is added between vib and via
    mesh%nC( vid) = mesh%nC( vid) + 1
    do ci = 1, mesh%nC( vid)
      if (mesh%C( vid,ci) == vib) then
        mesh%C( vid,:) = [mesh%C( vid,1:ci), vic, mesh%C( vid,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do

    ! == niTri, iTri

    ! via: tj,ti are replaced by t1
    ! First remove ti
    do iti = 1, mesh%niTri( via)
      if (mesh%iTri( via,iti) == ti) then
        mesh%iTri( via,:) = [mesh%iTri( via,1:iti-1), mesh%iTri( via,iti+1:mesh%nC_mem), 0]
        mesh%niTri( via) = mesh%niTri( via) - 1
        exit
      end if
    end do
    ! then replace tj by t1
    do iti = 1, mesh%niTri( via)
      if (mesh%iTri( via,iti) == tj) then
        mesh%iTri( via,iti) = t1
        exit
      end if
    end do

    ! vib: ti,tj are replaced by t2
    ! First remove tj
    do iti = 1, mesh%niTri( vib)
      if (mesh%iTri( vib,iti) == tj) then
        mesh%iTri( vib,:) = [mesh%iTri( vib,1:iti-1), mesh%iTri( vib,iti+1:mesh%nC_mem), 0]
        mesh%niTri( vib) = mesh%niTri( vib) - 1
        exit
      end if
    end do
    ! then replace ti by t2
    do iti = 1, mesh%niTri( vib)
      if (mesh%iTri( vib,iti) == ti) then
        mesh%iTri( vib,iti) = t2
        exit
      end if
    end do

    ! vic: ti is replaced by t1,t2
    mesh%niTri( vic) = mesh%niTri( vic) + 1
    do iti = 1, mesh%niTri( vic)
      if (mesh%iTri( vic,iti) == ti) then
        mesh%iTri( vic,:) = [mesh%iTri( vic,1:iti-1), t1, t2, mesh%iTri( vic,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do

    ! vid: tj is replaced by t2,t1
    mesh%niTri( vid) = mesh%niTri( vid) + 1
    do iti = 1, mesh%niTri( vid)
      if (mesh%iTri( vid,iti) == tj) then
        foundit = .true.
        mesh%iTri( vid,:) = [mesh%iTri( vid,1:iti-1), t2, t1, mesh%iTri( vid,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do

    ! == Border index

    ! No changes.

    ! == TriC

    ! tia: ti is replaced by t2
    if (tia > 0) then
      do n = 1, 3
        if (mesh%TriC( tia,n) == ti) then
          mesh%TriC( tia,n) = t2
          exit
        end if
      end do
    end if
    ! tib: ti is replaced by t1
    if (tib > 0) then
      do n = 1, 3
        if (mesh%TriC( tib,n) == ti) then
          mesh%TriC( tib,n) = t1
          exit
        end if
      end do
    end if
    ! tja: tj is replaced by t2
    if (tja > 0) then
      do n = 1, 3
        if (mesh%TriC( tja,n) == tj) then
          mesh%TriC( tja,n) = t2
          exit
        end if
      end do
    end if
    ! tjb: tj is replaced by t1
    if (tjb > 0) then
      do n = 1, 3
        if (mesh%TriC( tjb,n) == tj) then
          mesh%TriC( tjb,n) = t1
          exit
        end if
      end do
    end if
    ! t1 is adjacent to [t2, tib, tjb]
    mesh%TriC( t1,:) = [t2, tib, tjb]
    ! t2 is adjacent to [t1, tja, tia]
    mesh%TriC( t2,:) = [t1, tja, tia]

    ! == Tricc

    call update_triangle_circumcenter( mesh, t1)
    call update_triangle_circumcenter( mesh, t2)

    ! == Refinement data

    ! Add the two new triangles to the refinement map and stack
    call add_triangle_to_refinement_stack_first( mesh, t1)
    call add_triangle_to_refinement_stack_first( mesh, t2)

    ! Update triangle-line overlap ranges
    li_min = min( mesh%Tri_li( ti,1), mesh%Tri_li( tj,1))
    li_max = max( mesh%Tri_li( ti,2), mesh%Tri_li( tj,2))
    mesh%Tri_li( t1,:) = [li_min, li_max]
    mesh%Tri_li( t2,:) = [li_min, li_max]

    ! Add the four new pairs to the flip list
    if (tia > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tia, t2]
    end if
    if (tib > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tib, t1]
    end if
    if (tja > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tja, t2]
    end if
    if (tjb > 0) then
      nf = nf + 1
      mesh%Tri_flip_list( nf,:) = [tjb, t1]
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_triangle_pair

  function are_Delaunay( mesh, ti, tj) result( isso)
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

    ! In/output variables:
    type(type_mesh),            intent(in)        :: mesh
    integer,                    intent(in)        :: ti,tj
    logical                                       :: isso

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'are_Delaunay'
    logical                                       :: are_connected_ij, are_connected_ji
    integer                                       :: n, vi, vj, vii, n1, n2, n3, iti
    logical                                       :: is_in_tj
    integer                                       :: via, vib, vic, vid
    logical                                       :: via_has_ti, via_has_tj
    logical                                       :: vib_has_ti, vib_has_tj
    logical                                       :: vic_has_ti, vic_has_tj
    logical                                       :: vid_has_ti, vid_has_tj
    real(dp), dimension(2)                        :: va, vb, vc, vd, cci, ccj
    real(dp)                                      :: ccri, ccrj

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (ti == 0 .or. tj == 0) then
      call crash('Found ti=0 in mesh%Tri_flip_list!')
    end if

    ! Check if these two triangles are connected
    are_connected_ij = .false.
    are_connected_ji = .false.
    do n = 1, 3
      if (mesh%TriC( ti,n) == tj) are_connected_ij = .true.
      if (mesh%TriC( tj,n) == ti) are_connected_ji = .true.
    end do
    if (.not. are_connected_ij .and. .not. are_connected_ij) then
      ! These two triangles are not connected
      call crash('{int_01} and {int_02} are not connected!', int_01 = ti, int_02 = tj)
    elseif (are_connected_ij .and. .not. are_connected_ji .or. .not. are_connected_ij .and. are_connected_ji) then
      ! One of them lists the other as a neighbour, but not vice versa
      call crash('inconsistency in TriC!')
    end if

    ! Find the two vertices vi and vj that are shared by ti and tj

    vi = 0
    vj = 0

    do n = 1, 3
      vii = mesh%Tri( ti,n)
      is_in_tj = .false.
      do n2 = 1, 3
        if (mesh%Tri( tj,n2) == vii) then
          is_in_tj = .true.
          exit
        end if
      end do
      if (is_in_tj) then
        if (vi == 0) then
          vi = vii
        else
          vj = vii
        end if
      end if
    end do

    ! Safety
    if (vi == 0 .or. vj == 0) then
      call crash('couldnt find two shared vertices!')
    end if

    ! Find via,vib,vic,vid (see diagram)

    via = 0
    vib = 0
    vic = 0
    vid = 0

    do n1 = 1, 3

      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      if ((mesh%Tri( ti,n1) == vi .and. mesh%Tri( ti,n2) == vj) .or. &
          (mesh%Tri( ti,n1) == vj .and. mesh%Tri( ti,n2) == vi)) then
        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)
      end if

      if ((mesh%Tri( tj,n1) == vi .and. mesh%Tri( tj,n2) == vj) .or. &
          (mesh%Tri( tj,n1) == vj .and. mesh%Tri( tj,n2) == vi)) then
        vid = mesh%Tri( tj,n3)
      end if

    end do

    ! Safety
    if (via == 0 .or. vib == 0 .or. vic == 0 .or. vid == 0) then
      call crash('couldnt figure out local geometry!')
    end if

    via_has_ti = .false.
    via_has_tj = .false.
    do iti = 1, mesh%niTri( via)
      if     (mesh%iTri( via,iti) == ti) then
        via_has_ti = .true.
      elseif (mesh%iTri( via,iti) == tj) then
        via_has_tj = .true.
      end if
    end do
    if (.not. via_has_ti) call crash('inconsistent mesh geometry! (via doesnt have ti as an itriangle)')
    if (.not. via_has_tj) call crash('inconsistent mesh geometry! (via doesnt have tj as an itriangle)')

    vib_has_ti = .false.
    vib_has_tj = .false.
    do iti = 1, mesh%niTri( vib)
      if     (mesh%iTri( vib,iti) == ti) then
        vib_has_ti = .true.
      elseif (mesh%iTri( vib,iti) == tj) then
        vib_has_tj = .true.
      end if
    end do
    if (.not. vib_has_ti) call crash('inconsistent mesh geometry! (vib doesnt have ti as an itriangle)')
    if (.not. vib_has_tj) call crash('inconsistent mesh geometry! (vib doesnt have tj as an itriangle)')

    vic_has_ti = .false.
    vic_has_tj = .false.
    do iti = 1, mesh%niTri( vic)
      if     (mesh%iTri( vic,iti) == ti) then
        vic_has_ti = .true.
      elseif (mesh%iTri( vic,iti) == tj) then
        vic_has_tj = .true.
      end if
    end do
    if (.not. vic_has_ti) call crash('inconsistent mesh geometry! (vic doesnt have ti as an itriangle)')
    if (      vic_has_tj) then
      call warning('ti = [{int_01}, {int_02}, {int_03}], tj = [{int_04}, {int_05}, {int_06}]', &
        int_01 = mesh%Tri( ti,1), int_02 = mesh%Tri( ti,2), int_03 = mesh%Tri( ti,3), &
        int_04 = mesh%Tri( tj,1), int_05 = mesh%Tri( tj,2), int_06 = mesh%Tri( tj,3))
    end if
    if (      vic_has_tj) call crash('inconsistent mesh geometry! (vic has tj as an itriangle)')

    vid_has_ti = .false.
    vid_has_tj = .false.
    do iti = 1, mesh%niTri( vid)
      if     (mesh%iTri( vid,iti) == ti) then
        vid_has_ti = .true.
      elseif (mesh%iTri( vid,iti) == tj) then
        vid_has_tj = .true.
      end if
    end do
    if (      vid_has_ti) call crash('inconsistent mesh geometry! (vid has ti as an itriangle)')
    if (.not. vid_has_tj) call crash('inconsistent mesh geometry! (vid doesnt have tj as an itriangle)')

    ! Check if ti-tj meets the Delaunay criterion

    va = mesh%V( via,:)
    vb = mesh%V( vib,:)
    vc = mesh%V( vic,:)
    vd = mesh%V( vid,:)

    cci = mesh%Tricc( ti,:)
    ccj = mesh%Tricc( tj,:)

    ccri = norm2( va - cci)
    ccrj = norm2( va - ccj)

    isso = .true.

    if     (norm2( vd - cci) < ccri) then
      ! vid lies inside the circumcircle of ti
      isso = .false.
    elseif (norm2( vc - ccj) < ccrj) then
      ! vic lies inside the circumcircle of tj
      isso = .false.
    end if

    ! if the outer angle at via or vib is concave, don't flip.
    ! Check this by checking if via lies inside the triangle
    ! [vib,vic,vid], or the other way round.

    if  (is_in_triangle( vb, vc, vd, va) .or. &
         is_in_triangle( va, vd, vc, vb)) then
      isso = .true.
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function are_Delaunay

end module mesh_Delaunay
