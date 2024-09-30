module line_tracing_Voronoi

  ! Line tracing algorithm through mesh Voronoi cells

  use precisions, only: dp
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use line_tracing_basic, only: add_integrals_to_single_row
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use math_utilities, only: lies_on_line_segment, segment_intersection, crop_line_to_domain, &
    line_integral_xdy, line_integral_mxydx, line_integral_xydy
  use mesh_utilities, only: find_shared_Voronoi_boundary, is_in_Voronoi_cell, calc_Voronoi_cell, &
    find_containing_vertex

  implicit none

  private

  public :: trace_line_Vor

contains

  !> Trace the line [pq] through the Voronoi cells of the mesh and calculate
  !> the three line integrals for the line segments inside the different Voronoi cells
  subroutine trace_line_Vor( mesh, p, q, single_row, count_coincidences, vi_hint)

    ! In/output variables
    type(type_mesh),                        intent(in)    :: mesh
    real(dp), dimension(2),                 intent(in)    :: p,q
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    logical,                                intent(in)    :: count_coincidences
    integer,                                intent(inout) :: vi_hint

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'trace_line_Vor'
    real(dp), dimension(2)         :: pp,qq
    logical                        :: is_valid_line
    logical                        :: finished
    integer                        :: n_cycles
    integer                        :: vi_in, ti_on, ei_on
    real(dp), dimension(2)         :: p_next
    integer                        :: vi_left
    logical                        :: coincides
    real(dp)                       :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    call init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the mesh domain
    call crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    if (.not.is_valid_line) then
      ! [pq] doesn't pass through the mesh domain anywhere
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    call trace_line_Vor_start( mesh, pp, vi_hint, vi_in, ti_on, ei_on)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not.finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      if     (vi_in  > 0) then
        ! p lies inside the Voronoi cell of vertex vi_in
        call trace_line_Vor_vi( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      elseif (ti_on  > 0) then
        ! p lies on the circumcentre of triangle ti_on
        call trace_line_Vor_ti( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      elseif (ei_on > 0) then
        ! p lies on the shared Voronoi cell boundary represented by edge ei_on
        call trace_line_Vor_ei( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      end if

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, mesh%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, mesh%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, mesh%tol_dist)

      ! Add them to the results structure
      if (norm2( p_next - pp) > mesh%tol_dist) then
        call add_integrals_to_single_row( single_row, vi_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      end if

      ! cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      if (n_cycles > mesh%nV) then
        call crash('trace_line_Vor - iterative tracer got stuck!')
      end if

      ! Update vi_hint, for more efficiency
      vi_hint = vi_left

    end do ! do while (.not.finished)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine trace_line_Vor

  !> Initialise the coincidence indicators for the point p, i.e. check if p either...
  !>    - lies inside the Voronoi cell of vertex vi_in, ...
  !>    - lies on the circumcentre of triangle ti_on, or...
  !>    - lies on the shared Voronoi cell boundary represented by edge ei_on
  subroutine trace_line_Vor_start( mesh, p, vi_hint, vi_in, ti_on, ei_on)

    ! In/output variables
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p
    integer,                intent(inout) :: vi_hint
    integer,                intent(out)   :: vi_in
    integer,                intent(out)   :: ti_on
    integer,                intent(out)   :: ei_on

    ! Local variables:
    integer                :: vti, ti, vei, ei
    real(dp), dimension(2) :: cc1, cc2

    ! Initialise
    vi_in = 0
    ti_on = 0
    ei_on = 0

    ! Find the vertex whose Voronoi cell contains p
    call find_containing_vertex( mesh, p, vi_hint)

    ! Check if p lies on any of the surrounding triangles' circumcentres
    do vti = 1, mesh%niTri( vi_hint)
      ti = mesh%iTri( vi_hint,vti)
      if (norm2( mesh%Tricc( ti,:) - p) < mesh%tol_dist) then
        ! p lies on the circumcentre of triangle ti
        ti_on = ti
        return
      end if
    end do

    ! Check if p lies on any of the shared Voronoi boundaries
    do vei = 1, mesh%nC( vi_hint)
      ei = mesh%VE( vi_hint,vei)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, p, mesh%tol_dist)) then
        ! p lies on the shared Voronoi cell boundary represented by edge ei
        ei_on = ei
        return
      end if
    end do

    ! if p lies not on the boundary of the Voronoi cell, then it must lie inside of it
    vi_in = vi_hint

  end subroutine trace_line_Vor_start

  !> Given the line [pq], where p lies inside the Voronoi cell of vertex vi_in,
  !> find the point p_next where [pq] crosses into the next Voronoi cell.
  subroutine trace_line_Vor_vi(  mesh, p, q, &
    p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)

    ! In/output variables
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p,q
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(inout) :: vi_in
    integer,                intent(inout) :: ti_on
    integer,                intent(inout) :: ei_on
    integer,                intent(out)   :: vi_left
    logical,                intent(out)   :: coincides
    logical,                intent(out)   :: finished

    ! Local variables:
    integer                             :: vti, ti, ei, vori, vorj, ci, vj
    real(dp), dimension(2)              :: r, llis, pa, pb
    real(dp)                            :: dx
    logical                             :: do_cross
    real(dp), dimension( mesh%nC_mem,2) :: Vor
    integer,  dimension( mesh%nC_mem  ) :: Vor_vi
    integer,  dimension( mesh%nC_mem  ) :: Vor_ti
    integer                             :: nVor

    ! Safety
    if (vi_in == 0 .or. ti_on > 0 .or. ei_on > 0 .or. (.not. is_in_Voronoi_cell( mesh, p, vi_in))) then
      call crash('trace_line_Vor_vi - coincidence indicators dont make sense!')
    end if

    ! Check if q lies inside the same Voronoi cell
    if (is_in_Voronoi_cell( mesh, q, vi_in)) then
      ! q lies inside the same Voronoi cell
      p_next    = q
      vi_left   = vi_in
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on the boundary of this Voronoi cell
    dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 200
    call calc_Voronoi_cell( mesh, vi_in, dx, Vor, Vor_vi, Vor_ti, nVor)
    do vori = 1, nVor
      vorj = vori + 1
      if (vorj == nVor + 1) vorj = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori,:)
      pb = Vor( vorj,:)
      if (norm2( q - pa) < mesh%tol_dist .or. lies_on_line_segment( pa, pb, q, mesh%tol_dist)) then
        ! q lies on the boundary of the same Voronoi cell
        p_next    = q
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if [pq] passes through any of the surrounding triangles' circumcentres
    do vti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,vti)
      r  = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, r, mesh%tol_dist)) then
        ! [pq] passes through this triangle's circumcentre
        p_next    = mesh%Tricc( ti,:)
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] passes through any of the shared Voronoi boundaries
    do vori = 1, nVor
      vorj = vori + 1
      if (vorj == nVor + 1) vorj = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori,:)
      pb = Vor( vorj,:)
      ! The other vertex sharing this Voronoi cell boundary
      vj = Vor_vi( vori)
      ! The edge representing this shared Voronoi cell boundary
      ei = 0
      do ci = 1, mesh%nC( vi_in)
        if (mesh%C( vi_in,ci) == vj) then
          ei = mesh%VE( vi_in,ci)
          exit
        end if
      end do
      ! Safety
      if (ei == 0) call crash('couldnt find edge between vi and vj!')

      call segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes into the Voronoi cell of vj
        p_next    = llis
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! This point should not be reachable!
    call crash('trace_line_Vor_vi - reached the unreachable end of the subroutine!')

  end subroutine trace_line_Vor_vi

  !> Given the line [pq], where p lies on the circumcentre of triangle ti_on,
  !> find the point p_next where [pq] crosses into the next Voronoi cell.
  subroutine trace_line_Vor_ti(  mesh, p, q, &
    p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)

    ! In/output variables
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p,q
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(inout) :: vi_in
    integer,                intent(inout) :: ti_on
    integer,                intent(inout) :: ei_on
    integer,                intent(out)   :: vi_left
    logical,                intent(out)   :: coincides
    logical,                intent(out)   :: finished

    ! Local variables:
    integer                :: via, vib, vic, vvi, vj, ei, acab, acbc, acca, tj
    real(dp), dimension(2) :: cc, cc1, cc2, llis
    logical                :: do_cross

    ! Safety
    if (vi_in > 0 .or. ti_on == 0 .or. ei_on > 0 .or. norm2( mesh%Tricc( ti_on,:) - p) > mesh%tol_dist) then
      call crash('trace_line_Vor_ti - coincidence indicators dont make sense!')
    end if

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_on,1)
    vib = mesh%Tri( ti_on,2)
    vic = mesh%Tri( ti_on,3)

    ! Find the three Voronoi cell boundaries that meet here
    acab = 0
    do vvi = 1, mesh%nC( via)
      vj = mesh%C(  via,vvi)
      ei = mesh%VE( via,vvi)
      if (vj == vib) then
        acab = ei
        exit
      end if
    end do
    acbc = 0
    do vvi = 1, mesh%nC( vib)
      vj = mesh%C(  vib,vvi)
      ei = mesh%VE( vib,vvi)
      if (vj == vic) then
        acbc = ei
        exit
      end if
    end do
    acca = 0
    do vvi = 1, mesh%nC( vic)
      vj = mesh%C(  vic,vvi)
      ei = mesh%VE( vic,vvi)
      if (vj == via) then
        acca = ei
        exit
      end if
    end do

    ! Check if q lies on the Voronoi cell boundary separating via from vib
    call find_shared_Voronoi_boundary( mesh, acab, cc1, cc2)
    if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .or. &
        norm2( cc1 - q) < mesh%tol_dist .or. &
        norm2( cc2 - q) < mesh%tol_dist) then
      ! q lies on the Voronoi cell boundary separating via from vib
      if (mesh%ETri( acab,1) == ti_on) then
        vi_left = mesh%EV( acab,2)
      else
        vi_left = mesh%EV( acab,1)
      end if
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the Voronoi cell boundary separating vib from vic
    call find_shared_Voronoi_boundary( mesh, acbc, cc1, cc2)
    if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .or. &
        norm2( cc1 - q) < mesh%tol_dist .or. &
        norm2( cc2 - q) < mesh%tol_dist) then
      ! q lies on the Voronoi cell boundary separating vib from vic
      if (mesh%ETri( acbc,1) == ti_on) then
        vi_left = mesh%EV( acbc,2)
      else
        vi_left = mesh%EV( acbc,1)
      end if
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the Voronoi cell boundary separating vic from via
    call find_shared_Voronoi_boundary( mesh, acca, cc1, cc2)
    if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .or. &
        norm2( cc1 - q) < mesh%tol_dist .or. &
        norm2( cc2 - q) < mesh%tol_dist) then
      ! q lies on the Voronoi cell boundary separating vic from via
      if (mesh%ETri( acca,1) == ti_on) then
        vi_left = mesh%EV( acca,2)
      else
        vi_left = mesh%EV( acca,1)
      end if
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside any of the three adjacent Voronoi cells
    if (is_in_Voronoi_cell( mesh, q, via)) then
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .false.
      finished  = .true.
      return
    end if
    if (is_in_Voronoi_cell( mesh, q, vib)) then
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .false.
      finished  = .true.
      return
    end if
    if (is_in_Voronoi_cell( mesh, q, vic)) then
      ! q lies inside the Voronoi cell of vic
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vic
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the circumcentre of any of the three neighbouring triangles
    tj = mesh%TriC( ti_on,1)
    if (tj > 0) then
      cc = mesh%Tricc( tj,:)
      if (lies_on_line_segment( p, q, cc, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = vic
        coincides = .true.
        finished  = .false.
        return
      end if
    end if

    tj = mesh%TriC( ti_on,2)
    if (tj > 0) then
      cc = mesh%Tricc( tj,:)
      if (lies_on_line_segment( p, q, cc, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = via
        coincides = .true.
        finished  = .false.
        return
      end if
    end if

    tj = mesh%TriC( ti_on,3)
    if (tj > 0) then
      cc = mesh%Tricc( tj,:)
      if (lies_on_line_segment( p, q, cc, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = vib
        coincides = .true.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses the boundary of the Voronoi cell of via
    do vvi = 1, mesh%nC( via)
      vj = mesh%C( via,vvi)
      if (vj == vib .or. vj == vic) cycle
      ei = mesh%VE( via,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = via
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] crosses the boundary of the Voronoi cell of vib
    do vvi = 1, mesh%nC( vib)
      vj = mesh%C( vib,vvi)
      if (vj == via .or. vj == vic) cycle
      ei = mesh%VE( vib,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vib
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] crosses the boundary of the Voronoi cell of vic
    do vvi = 1, mesh%nC( vic)
      vj = mesh%C( vic,vvi)
      if (vj == via .or. vj == vib) cycle
      ei = mesh%VE( vic,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vic
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! This point should not be reachable!
    call crash('trace_line_Vor_ti - reached the unreachable end of the subroutine!')

  end subroutine trace_line_Vor_ti

  !> Given the line [pq], where p lies on the shared Voronoi boundary represented by edge ei_on,
  !> find the point p_next where [pq] crosses into the next Voronoi cell.
  subroutine trace_line_Vor_ei( mesh, p, q, &
    p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)

    ! In/output variables
    type(type_mesh),        intent(in)    :: mesh
    real(dp), dimension(2), intent(in)    :: p,q
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(inout) :: vi_in
    integer,                intent(inout) :: ti_on
    integer,                intent(inout) :: ei_on
    integer,                intent(out)   :: vi_left
    logical,                intent(out)   :: coincides
    logical,                intent(out)   :: finished

    ! Local variables:
    integer                :: via, vib, vil, vir, til, tir, vvi, ei, vti, ti
    real(dp), dimension(2) :: cc1, cc2, ccl, ccr, llis
    logical                :: do_cross

    ! Find the endpoints of this shared Voronoi boundary
    call find_shared_Voronoi_boundary( mesh, ei_on, cc1, cc2)

    ! Safety
    if (vi_in > 0 .or. ti_on > 0 .or. ei_on == 0 .or. (.not. lies_on_line_segment( cc1, cc2, p, mesh%tol_dist))) then
      call crash('trace_line_Vor_ei - coincidence indicators dont make sense!')
    end if

    ! A bit more detail is needed
    via = mesh%EV(   ei_on,1)
    vib = mesh%EV(   ei_on,2)
    vil = mesh%EV(   ei_on,3)
    vir = mesh%EV(   ei_on,4)
    til = mesh%ETri( ei_on,1)
    tir = mesh%ETri( ei_on,2)

    if (til == 0) then
      ! Apparently ei lies on the domain border and has no triangle on its left-hand side
      ccr = cc1
      ccl = cc2
    elseif (tir == 0) then
      ! Apparently ei lies on the domain border and has no triangle on its right-hand side
      ccl = cc1
      ccr = cc2
    else
      ! ei lies in the interior and has triangles on both sides
      ccl = mesh%Tricc( til,:)
      ccr = mesh%Tricc( tir,:)
    end if

    ! Check if q coincides with ccl
    if (norm2( ccl - q) < mesh%tol_dist) then
      ! q coincides with ccl
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q coincides with ccr
    if (norm2( ccr - q) < mesh%tol_dist) then
      ! q coincides with ccr
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through ccl
    if (lies_on_line_segment( p, q, ccl, mesh%tol_dist)) then
      ! [pq] passes through ccl
      p_next    = ccl
      vi_in     = 0
      ti_on     = til
      ei_on     = 0
      vi_left   = via
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through ccr
    if (lies_on_line_segment( p, q, ccr, mesh%tol_dist)) then
      ! [pq] passes through ccr
      p_next    = ccr
      vi_in     = 0
      ti_on     = tir
      ei_on     = 0
      vi_left   = vib
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if q lies inside the Voronoi cell of via
    if (is_in_Voronoi_cell( mesh, q, via)) then
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the Voronoi cell of vib
    if (is_in_Voronoi_cell( mesh, q, vib)) then
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on the circumcentre of any of the triangles surrounding via
    do vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      if (norm2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) then
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = via
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies on the circumcentre of any of the triangles surrounding vib
    do vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      if (norm2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) then
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = vib
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies on boundary of the Voronoi cell of via
    do vvi = 1, mesh%nC( via)
      ei = mesh%VE( via,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) then
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = via
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies on boundary of the Voronoi cell of vib
    do vvi = 1, mesh%nC( vib)
      ei = mesh%VE( vib,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) then
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = vib
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if pq crosses the circumcentre of any of the triangles surrounding via
    do vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      cc1 = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of triangle ti
        p_next    = cc1
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        vi_left   = via
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if pq crosses the boundary of the Voronoi cell of via
    do vvi = 1, mesh%nC( via)
      ei = mesh%VE( via,vvi)
      if (ei == ei_on) cycle
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = via
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if pq crosses the circumcentre of any of the triangles surrounding vib
    do vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      cc1 = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of triangle ti
        p_next    = cc1
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        vi_left   = vib
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if pq crosses the boundary of the Voronoi cell of vib
    do vvi = 1, mesh%nC( vib)
      ei = mesh%VE( vib,vvi)
      if (ei == ei_on) cycle
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vib
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! This point should not be reachable!
    call crash('trace_line_Vor_ei - reached the unreachable end of the subroutine!')

  end subroutine trace_line_Vor_ei

end module line_tracing_Voronoi
