module split_edges

  ! Split an edge of the mesh, and update the Delaunay triangulation accordingly.

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use plane_geometry, only: lies_on_line_segment
  use mesh_utilities, only: is_border_edge, update_triangle_circumcenter, add_triangle_to_refinement_stack_last
  use split_border_edges, only: split_border_edge
  use flip_triangles, only: flip_triangles_until_Delaunay

  implicit none

  logical :: do_debug = .false.

  private

  public :: split_edge

contains

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
    type(type_mesh),         intent(inout) :: mesh
    integer,                 intent(in)    :: vi, vj
    real(dp), dimension(2),  intent(in)    :: p_new

    ! Local variables:
    character(len=256), parameter :: routine_name = 'split_edge'
    integer                       :: ci, iti, n, n1, n2, n3, nf
    real(dp), dimension(2)        :: p, pa, pb
    integer                       :: t1, t2, t3, t4, ti, tit, tib, titl, titr, tibl, tibr, vib, vit, vk
    integer                       :: li_min, li_max

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Safety - assert that vi and vj are valid vertices
    call assert( test_ge_le( vi, 1, mesh%nV), 'invalid value for vi')
    call assert( test_ge_le( vj, 1, mesh%nV), 'invalid value for vj')

    ! Safety - assert that vi and vj are connected both ways
    call assert( test_mesh_vertices_are_neighbours( mesh, vi, vj), 'vi and vj are not connected')

    ! Safety - assert that p_new lies on the line [vi,vj]
    pa = mesh%V( vi,:)
    pb = mesh%V( vj,:)
    call assert( lies_on_line_segment( pa, pb, p_new, mesh%tol_dist), 'p does not lie on line [vi,vj]')
#endif

    ! If [vi,vj] is a border edge, split that instead
    if (is_border_edge( mesh, vi, vj)) then
      p = (pa + pb) / 2._dp
      call split_border_edge( mesh, vi, vj, p)
      call finalise_routine( routine_name)
      return
    end if

    ! Split the triangles t_top and t_bot adjacent to line [vi,vj] into four new ones
    ! ===============================================================================

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

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_ge_le( tit, 1, mesh%nTri), 'couldnt find valid triangle tit adjacent to [vi,vj]')
    call assert( test_ge_le( tib, 1, mesh%nTri), 'couldnt find valid triangle tib adjacent to [vi,vj]')
#endif

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

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_ge_le( vit, 1, mesh%nV), 'couldnt find valid vertex vit opposite [vi,vj]')
    call assert( test_ge_le( vib, 1, mesh%nV), 'couldnt find valid vertex vib opposite [vi,vj]')
#endif

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

#if (DO_ASSERTIONS)
    ! Safety - check if everything went alright and we didnt create any duplicate triangles
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t1), 'a triangle with the vertices of t1 already exists')
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t2), 'a triangle with the vertices of t2 already exists')
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t3), 'a triangle with the vertices of t3 already exists')
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t4), 'a triangle with the vertices of t4 already exists')
#endif

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

end module split_edges
