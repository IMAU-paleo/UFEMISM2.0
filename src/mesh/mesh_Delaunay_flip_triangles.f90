module mesh_Delaunay_flip_triangles

  ! Iteratively flip triangle pairs until the global Delaunay criterion is satisfied

  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use mesh_Delaunay_check_Delaunay_criterion, only: are_Delaunay
  use mesh_utilities, only: write_mesh_to_text_file, check_if_triangle_already_exists, &
    update_triangle_circumcenter, add_triangle_to_refinement_stack_first

  implicit none

  logical :: do_debug = .false.

  private

  public :: flip_triangles_until_Delaunay

contains

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

end module mesh_Delaunay_flip_triangles
