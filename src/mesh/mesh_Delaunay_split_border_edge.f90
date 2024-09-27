module mesh_Delaunay_split_border_edge

  ! Split a border edge of the mesh, and update the Delaunay triangulation accordingly.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use math_utilities, only: lies_on_line_segment
  use mesh_utilities, only: is_border_edge, check_if_triangle_already_exists, update_triangle_circumcenter, &
    add_triangle_to_refinement_stack_last
  use mesh_Delaunay_flip_triangles, only: flip_triangles_until_Delaunay

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

end module mesh_Delaunay_split_border_edge
