module split_border_edges

  ! Split a border edge of the mesh, and update the Delaunay triangulation accordingly.

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use plane_geometry, only: lies_on_line_segment
  use mesh_utilities, only: is_border_edge, update_triangle_circumcenter, add_triangle_to_refinement_stack_last
  use flip_triangles, only: flip_triangles_until_Delaunay

  implicit none

  logical :: do_debug = .false.

  private

  public :: split_border_edge

contains

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

    ! In/output variables:
    type(type_mesh),        intent(inout) :: mesh
    integer,                intent(in)    :: vi, vj
    real(dp), dimension(2), intent(in)    :: p_new

    ! Local variables:
    character(len=256), parameter :: routine_name = 'split_border_edge'
    integer                       :: ci, iti, n, n1, n2, n3, nf
    real(dp), dimension(2)        :: pa, pb
    integer                       :: t1, t2, ti, tia, tib, tic, tii, via, vib, vic, vik
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

    ! Safety - assert that [vi,vj] is a border edge
    call assert( is_border_edge( mesh, vi, vj), '[vi,vj] is not a border edge')
#endif

    ! Split the triangle ti adjacent to borderedge [vi,vj] into two new ones
    ! ======================================================================

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

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_ge_le( ti, 1, mesh%nTri), 'couldnt find valid triangle ti containing vi and vj')
#endif

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

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_ge_le( via, 1, mesh%nV), 'invalid value for via in mesh%Tri')
    call assert( test_ge_le( vib, 1, mesh%nV), 'invalid value for vib in mesh%Tri')
    call assert( test_ge_le( vic, 1, mesh%nV), 'invalid value for vic in mesh%Tri')
    call assert( (tia /= 0 .or. tib /= 0), 'triangle ti has only one neighbour')
    call assert( tic == 0, 'triangle ti does not appear to be a border triangle (since it has three neighbours)')
#endif

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

#if (DO_ASSERTIONS)
    ! Safety - check if everything went alright and we didnt create any duplicate triangles
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t1), 'a triangle with the vertices of t1 already exists')
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t2), 'a triangle with the vertices of t2 already exists')
#endif

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

end module split_border_edges
