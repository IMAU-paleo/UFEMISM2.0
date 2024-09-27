module mesh_Delaunay_split_triangle

  ! Split a triangle of the mesh, and update the Delaunay triangulation accordingly.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use math_utilities, only: lies_on_line_segment
  use mesh_utilities, only: encroaches_upon_any, find_containing_triangle, encroaches_upon, is_border_edge, &
    check_if_triangle_already_exists, update_triangle_circumcenter, add_triangle_to_refinement_stack_last
  use mesh_Delaunay_split_edge, only: split_edge
  use mesh_Delaunay_split_border_edge, only: split_border_edge
  use mesh_Delaunay_flip_triangles, only: flip_triangles_until_Delaunay

  implicit none

  logical :: do_debug = .false.

contains

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

end module mesh_Delaunay_split_triangle
