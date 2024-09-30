module line_tracing_triangles

  ! Line tracing algorithm through mesh triangles

  use precisions, only: dp
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use line_tracing_basic, only: add_integrals_to_single_row
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use math_utilities, only: lies_on_line_segment, segment_intersection, crop_line_to_domain, &
    line_integral_xdy, line_integral_mxydx, line_integral_xydy, is_in_triangle
  use mesh_utilities, only: find_containing_triangle

  implicit none

  private

  public :: trace_line_tri

contains

  !> Trace the line [pq] through the triangles of the mesh and calculate
  !> the three line integrals for the line segments inside the different triangles
  subroutine trace_line_tri( mesh, p, q, single_row, count_coincidences, ti_hint)

    ! In/output variables
    type(type_mesh),                        intent(in)    :: mesh
    real(dp), dimension(2),                 intent(in)    :: p,q
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    logical,                                intent(in)    :: count_coincidences
    integer,                                intent(inout) :: ti_hint

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'trace_line_tri'
    real(dp), dimension(2)         :: pp, qq
    logical                        :: is_valid_line
    logical                        :: finished
    integer                        :: n_cycles
    integer                        :: ti_in, vi_on, ei_on
    real(dp), dimension(2)         :: p_next
    integer                        :: ti_left
    logical                        :: coincides
    real(dp)                       :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    call init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the mesh domain
    call crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    if (.not. is_valid_line) then
      ! [pq] doesn't pass through the mesh domain anywhere
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    call trace_line_tri_start( mesh, pp, ti_hint, ti_in, vi_on, ei_on)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      if     (ti_in  > 0) then
        ! p lies inside triangle ti_in
        call trace_line_tri_ti( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      elseif (vi_on  > 0) then
        ! p lies on vertex vi_on
        call trace_line_tri_vi( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      elseif (ei_on > 0) then
        ! p lies on edge ei_on
        call trace_line_tri_ei( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      else
        call crash('coincidence indicators dont make sense!')
      end if

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, mesh%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, mesh%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, mesh%tol_dist)

      ! Add them to the results structure
      if (norm2( p_next - pp) > mesh%tol_dist) then
        call add_integrals_to_single_row( single_row, ti_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      end if

      ! cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      if (n_cycles > mesh%nV) then
        call crash('trace_line_tri - iterative tracer got stuck!')
      end if

      ! Update ti_hint, for more efficiency
      ti_hint = ti_left

    end do ! do while (.not. finished)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine trace_line_tri

  !> Initialise the coincidence indicators for the point p, i.e. check if p either...
  !>    - lies inside triangle ti_in, ...
  !>    - lies on vertex vi_on, or...
  !>    - lies on edge ei_on
  subroutine trace_line_tri_start( mesh, p, ti_hint, ti_in, vi_on, ei_on)

    ! In/output variables:
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p
    integer,                intent(inout) :: ti_hint
    integer,                intent(out)   :: ti_in
    integer,                intent(out)   :: vi_on
    integer,                intent(out)   :: ei_on

    ! Local variables:
    integer                :: via, vib, vic
    real(dp), dimension(2) :: pa, pb, pc
    integer                :: vvi, vj, ei

    ! Initialise
    ti_in  = 0
    vi_on  = 0
    ei_on = 0

    ! Find the triangle containing p
    call find_containing_triangle( mesh, p, ti_hint)

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_hint,1)
    vib = mesh%Tri( ti_hint,2)
    vic = mesh%Tri( ti_hint,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Check if p lies on any of the three vertices
    if     (norm2( pa - p) < mesh%tol_dist) then
      ! p lies on via
      vi_on = via
      return
    elseif (norm2( pb - p) < mesh%tol_dist) then
      ! p lies on vib
      vi_on = vib
      return
    elseif (norm2( pc - p) < mesh%tol_dist) then
      ! p lies on vic
      vi_on = vic
      return
    end if

    ! Check if p lies on any of the three edges
    if     (lies_on_line_segment( pa, pb, p, mesh%tol_dist)) then
      ! p lies on the edge connecting via and vib
      do vvi = 1, mesh%nC( via)
        vj = mesh%C(  via,vvi)
        ei = mesh%VE( via,vvi)
        if (vj == vib) then
          ei_on = ei
          return
        end if
      end do
    elseif (lies_on_line_segment( pb, pc, p, mesh%tol_dist)) then
      ! p lies on the edge connecting vib and vic
      do vvi = 1, mesh%nC( vib)
        vj = mesh%C(  vib,vvi)
        ei = mesh%VE( vib,vvi)
        if (vj == vic) then
          ei_on = ei
          return
        end if
      end do
    elseif (lies_on_line_segment( pc, pa, p, mesh%tol_dist)) then
      ! p lies on the edge connecting vic and via
      do vvi = 1, mesh%nC( vic)
        vj = mesh%C(  vic,vvi)
        ei = mesh%VE( vic,vvi)
        if (vj == via) then
          ei_on = ei
          return
        end if
      end do
    end if

    ! if p lies not on the vertices or edges of the triangle, then it must lie inside of it
    ti_in = ti_hint

  end subroutine trace_line_tri_start

  !> Given the line [pq], where p lies inside triangle ti_in,
  !> find the point p_next where [pq] crosses into the next triangle.
  subroutine trace_line_tri_ti( mesh, p, q, &
    p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)

    ! In/output variables
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p,q
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(inout) :: ti_in
    integer,                intent(inout) :: vi_on
    integer,                intent(inout) :: ei_on
    integer,                intent(out)   :: ti_left
    logical,                intent(out)   :: coincides
    logical,                intent(out)   :: finished

    ! Local variables:
    integer                :: via, vib, vic
    real(dp), dimension(2) :: pa, pb, pc
    integer                :: vvi, vj, ei
    real(dp), dimension(2) :: llis
    logical                :: do_cross

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_in,1)
    vib = mesh%Tri( ti_in,2)
    vic = mesh%Tri( ti_in,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Safety
    if (ti_in == 0 .or. vi_on > 0 .or. ei_on > 0) then
      call crash('trace_line_tri_ti - coincidence indicators dont make sense!')
    end if
    if (.not. is_in_triangle( pa, pb, pc, p)) then
      call crash('trace_line_tri_ti - p does not lie inside triangle ti_in!')
    end if

    ! Check if q lies inside the same triangle
    if (is_in_triangle( pa, pb, pc, q)) then
      ! q lies inside the same triangle
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on vertex via
    if (norm2( pa - q) < mesh%tol_dist) then
      ! q lies on vertex via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on vertex vib
    if (norm2( pb - q) < mesh%tol_dist) then
      ! q lies on vertex vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on vertex vic
    if (norm2( pc - q) < mesh%tol_dist) then
      ! q lies on vertex vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on edge via-vib
    if (lies_on_line_segment( pa, pb, q, mesh%tol_dist)) then
      ! q lies on edge via-vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on edge vib-vic
    if (lies_on_line_segment( pb, pc, q, mesh%tol_dist)) then
      ! q lies on edge vib-vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on edge vic-via
    if (lies_on_line_segment( pc, pa, q, mesh%tol_dist)) then
      ! q lies on edge vic-via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through via
    if (lies_on_line_segment( p, q, pa, mesh%tol_dist)) then
      ! [pq] passes through via
      p_next    = pa
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = via
      ei_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through vib
    if (lies_on_line_segment( p, q, pb, mesh%tol_dist)) then
      ! [pq] passes through vib
      p_next    = pb
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vib
      ei_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through vic
    if (lies_on_line_segment( p, q, pc, mesh%tol_dist)) then
      ! [pq] passes through vic
      p_next    = pc
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vic
      ei_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] crosses edge via-vib
    call segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
    if (do_cross) then
      ! [pq] crosses edge [via,vib]
      ! Find the edge connecting via and vib
      do vvi = 1, mesh%nC( via)
        vj = mesh%C(  via,vvi)
        ei = mesh%VE( via,vvi)
        if (vj == vib) then
          ei_on = ei
          exit
        end if
      end do
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] crosses edge vib-vic
    call segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
    if (do_cross) then
      ! [pq] crosses edge [vib,vic]
      ! Find the edge connecting vib and vic
      do vvi = 1, mesh%nC( vib)
        vj = mesh%C(  vib,vvi)
        ei = mesh%VE( vib,vvi)
        if (vj == vic) then
          ei_on = ei
          exit
        end if
      end do
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] crosses edge vic-via
    call segment_intersection( p, q, pc, pa, llis, do_cross, mesh%tol_dist)
    if (do_cross) then
      ! [pq] crosses edge [vic,via]
      ! Find the edge connecting vic and via
      do vvi = 1, mesh%nC( vic)
        vj = mesh%C(  vic,vvi)
        ei = mesh%VE( vic,vvi)
        if (vj == via) then
          ei_on = ei
          exit
        end if
      end do
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_tri_ti - reached the unreachable end of the subroutine!')

  end subroutine trace_line_tri_ti

  !> Given the line [pq], where p lies on vertex vi_on,
  !> find the point p_next where [pq] crosses into the next triangle.
  subroutine trace_line_tri_vi( mesh, p, q, &
    p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)

    ! In/output variables
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p,q
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(inout) :: ti_in
    integer,                intent(inout) :: vi_on
    integer,                intent(inout) :: ei_on
    integer,                intent(out)   :: ti_left
    logical,                intent(out)   :: coincides
    logical,                intent(out)   :: finished

    ! Local variables:
    integer                                            :: via, vib, vic
    real(dp), dimension(2)                             :: pa, pb, pc, pv
    integer                                            :: vvi, vj, ei, n1, n2, n3, vti, ti
    real(dp), dimension(2)                             :: llis
    logical                                            :: do_cross

    ! Safety
    if (ti_in > 0 .or. vi_on == 0 .or. ei_on > 0) then
      call crash('trace_line_tri_vi - coincidence indicators dont make sense!')
    end if
    if (norm2( p - mesh%V( vi_on,:)) > mesh%tol_dist) then
      call crash('trace_line_tri_vi - p does not lie on vertex vi_on!')
    end if

    ! Check if q lies on any of the edges originating in this vertex
    do vvi = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,vvi)
      ei = mesh%VE( vi_on,vvi)
      pv = mesh%V( vj,:)
      if (norm2( mesh%V( vj,:) - q) < mesh%tol_dist .or. &
          lies_on_line_segment( p, pv, q, mesh%tol_dist)) then
        ! q lies on edge ei, connecting vi_on and vj
        if (mesh%EV( ei,1) == vi_on) then
          ti_left = mesh%ETri( ei,1)
        else
          ti_left = mesh%ETri( ei,2)
        end if
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        coincides = .true.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies inside any of the triangles surrounding vi_on
    do vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)
      if (is_in_triangle( pa, pb, pc, q) .or. &
          lies_on_line_segment( pa, pb, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pb, pc, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pc, pa, q, mesh%tol_dist)) then
        ! q lies inside adjacent triangle ti
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = ti
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if [pq] passes through any of the neighbouring vertices
    do vvi = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,vvi)
      ei = mesh%VE( vi_on,vvi)
      pv = mesh%V( vj,:)
      if (lies_on_line_segment( p, q, pv, mesh%tol_dist)) then
        ! [pq] passes through neighbouring vertex vj, which is connected to vi_on by edge ei
        p_next    = pv
        ti_in     = 0
        vi_on     = vj
        ei_on     = 0
        if (mesh%EV( ei,1) == vi_on) then
          ti_left = mesh%ETri( ei,1)
        else
          ti_left = mesh%ETri( ei,2)
        end if
        coincides = .true.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] exits into any of the adjacent triangles
    do vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        n3 = n2 + 1
        if (n3 == 4) n3 = 1
        if (mesh%Tri( ti,n1) == vi_on) then
          vib = mesh%Tri( ti,n2)
          vic = mesh%Tri( ti,n3)
          pb  = mesh%V( vib,:)
          pc  = mesh%V( vic,:)
          ! Find the opposite triangle edge
          ei = 0
          do vvi = 1, mesh%nC( vib)
            vj = mesh%C( vib,vvi)
            if (vj == vic) then
              ei = mesh%VE( vib,vvi)
              exit
            end if
          end do
          call segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
          if (do_cross) then
            ! [pq] exits triangle ti through the opposite edge ei
            p_next    = llis
            ti_in     = 0
            vi_on     = 0
            ei_on     = ei
            ti_left   = ti
            coincides = .false.
            finished  = .false.
            return
          end if
        end if
      end do
    end do

    ! This point should not be reachable!
    call crash('trace_line_tri_vi - reached the unreachable end of the subroutine!')

  end subroutine trace_line_tri_vi

  !> Given the line [pq], where p lies on edge ei,
  !> find the point p_next where [pq] crosses into the next triangle.
  subroutine trace_line_tri_ei( mesh, p, q, &
    p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)

    ! In/output variables
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p,q
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(inout) :: ti_in
    integer,                intent(inout) :: vi_on
    integer,                intent(inout) :: ei_on
    integer,                intent(out)   :: ti_left
    logical,                intent(out)   :: coincides
    logical,                intent(out)   :: finished

    ! Local variables:
    integer                :: via, vib, vil, vir, til, tir
    real(dp), dimension(2) :: pa, pb, pl, pr
    integer                :: vvi, vj, ei
    real(dp), dimension(2) :: llis
    logical                :: do_cross

    ! Some more info about this edge
    via = mesh%EV(   ei_on,1)
    vib = mesh%EV(   ei_on,2)
    vil = mesh%EV(   ei_on,3)
    vir = mesh%EV(   ei_on,4)
    til = mesh%ETri( ei_on,1)
    tir = mesh%ETri( ei_on,2)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    if (vil > 0) pl  = mesh%V( vil,:)
    if (vir > 0) pr  = mesh%V( vir,:)

    ! Safety
    if (ti_in > 0 .or. vi_on > 0 .or. ei_on == 0) then
      call crash('trace_line_tri_ei - coincidence indicators dont make sense!')
    end if
    if (.not. lies_on_line_segment( pa, pb, p, mesh%tol_dist)) then
      call crash('trace_line_tri_ei - p does not lie on edge ei_on!')
    end if

    ! Check if q lies on the same edge in the direction of via
    if (lies_on_line_segment( p, pa, q, mesh%tol_dist)) then
      ! q lies on the same edge in the direction of via
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      ti_left   = tir
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the same edge in the direction of vib
    if (lies_on_line_segment( p, pb, q, mesh%tol_dist)) then
      ! q lies on the same edge in the direction of vib
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      ti_left   = til
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside either of the two adjacent triangles
    if (til > 0) then
      if (is_in_triangle( pa, pb, pl, q)) then
        ! q lies inside triangle til
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .true.
        return
      end if
    end if
    if (tir > 0) then
      if (is_in_triangle( pa, pr, pb, q)) then
        ! q lies inside triangle tir
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .true.
        return
      end if
    end if

    ! Check if [pq] passes through pa
    if (lies_on_line_segment( p, q, pa, mesh%tol_dist)) then
      ! [pq] passes through pa
      p_next    = pa
      ti_in     = 0
      vi_on     = via
      ei_on     = 0
      ti_left   = tir
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through pb
    if (lies_on_line_segment( p, q, pb, mesh%tol_dist)) then
      ! [pq] passes through pb
      p_next    = pb
      ti_in     = 0
      vi_on     = vib
      ei_on     = 0
      ti_left   = til
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through pl
    if (til > 0) then
      if (lies_on_line_segment( p, q, pl, mesh%tol_dist)) then
        ! [pq] passes through pl
        p_next    = pl
        ti_in     = 0
        vi_on     = vil
        ei_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] passes through pr
    if (tir > 0) then
      if (lies_on_line_segment( p, q, pr, mesh%tol_dist)) then
        ! [pq] passes through pr
        p_next    = pr
        ti_in     = 0
        vi_on     = vir
        ei_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [via,vil]
    if (til > 0) then
      call segment_intersection( p, q, pa, pl, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [via,vil]
        ! Find the edge connecting via and vil
        do vvi = 1, mesh%nC( via)
          vj = mesh%C(  via,vvi)
          ei = mesh%VE( via,vvi)
          if (vj == vil) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [vil,vib]
    if (til > 0) then
      call segment_intersection( p, q, pl, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [vil,vib]
        ! Find the edge connecting vil and vib
        do vvi = 1, mesh%nC( vil)
          vj = mesh%C(  vil,vvi)
          ei = mesh%VE( vil,vvi)
          if (vj == vib) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [via,vir]
    if (tir > 0) then
      call segment_intersection( p, q, pa, pr, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [via,vir]
        ! Find the edge connecting via and vir
        do vvi = 1, mesh%nC( via)
          vj  = mesh%C(    via,vvi)
          ei = mesh%VE( via,vvi)
          if (vj == vir) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [vir,vib]
    if (tir > 0) then
      call segment_intersection( p, q, pr, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [vir,vib]
        ! Find the edge connecting vir and vib
        do vvi = 1, mesh%nC( vir)
          vj = mesh%C(  vir,vvi)
          ei = mesh%VE( vir,vvi)
          if (vj == vib) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! This point should not be reachable!
    call crash('trace_line_tri_ei - reached the unreachable end of the subroutine!')

  end subroutine trace_line_tri_ei

end module line_tracing_triangles
