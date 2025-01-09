module line_tracing_triangles

  ! Line tracing algorithm through mesh triangles

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use line_tracing_basic
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use plane_geometry, only: lies_on_line_segment, segment_intersection, crop_line_to_domain, &
    is_in_triangle
  use line_integrals, only: line_integral_xdy, line_integral_mxydx, line_integral_xydy
  use mesh_utilities, only: find_containing_triangle

  implicit none

  private

  public :: trace_line_tri
  public :: trace_line_tri_start, trace_line_tri_ti, trace_line_tri_vi, trace_line_tri_ei

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
    real(dp), dimension(2)    :: pp, qq
    logical                   :: is_valid_line
    logical                   :: finished
    integer                   :: n_cycles
    type(type_coinc_ind_mesh) :: coinc_ind
    real(dp), dimension(2)    :: p_next
    integer                   :: ti_left
    logical                   :: coincides
    real(dp)                  :: LI_xdy, LI_mxydx, LI_xydy

    ! Crop the line [pq] so that it lies within the mesh domain
    call crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, &
      pp, qq, is_valid_line)

    if (.not. is_valid_line) then
      ! [pq] doesn't pass through the mesh domain anywhere
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    call trace_line_tri_start( mesh, pp, ti_hint, coinc_ind)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      select case (coinc_ind%grid)
      case (no_value)
        call crash('trace_line_tri - coincidence indicator doesnt make sense')
      case default
        call crash('trace_line_tri - coincidence indicator doesnt make sense')
      case( a_grid)
        call trace_line_tri_vi( mesh, pp, qq, p_next, coinc_ind, ti_left, coincides, finished)
      case (b_grid)
        call trace_line_tri_ti( mesh, pp, qq, p_next, coinc_ind, ti_left, coincides, finished)
      case (c_grid)
        call trace_line_tri_ei( mesh, pp, qq, p_next, coinc_ind, ti_left, coincides, finished)
      end select

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
      if (ti_left > 0) ti_hint = ti_left

    end do ! do while (.not. finished)

  end subroutine trace_line_tri

  !> Initialise the coincidence indicators for the point p, i.e. check if p either...
  !>    - lies inside triangle ti_in, ...
  !>    - lies on vertex vi_on, or...
  !>    - lies on edge ei_on
  subroutine trace_line_tri_start( mesh, p, ti_hint, coinc_ind)

    ! In/output variables:
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p
    integer,                   intent(inout) :: ti_hint
    type(type_coinc_ind_mesh), intent(out)   :: coinc_ind

    ! Local variables:
    integer                :: via, vib, vic
    real(dp), dimension(2) :: pa, pb, pc
    integer                :: vvi, vj

#if (DO_ASSERTIONS)
    ! Safety - check if p lies inside the mesh domain
    call assert( &
      test_ge_le( p(1), mesh%xmin, mesh%xmax) .and. &
      test_ge_le( p(2), mesh%ymin, mesh%ymax), &
      'trace_line_tri_start - p lies outside the mesh domain')
#endif

    ! Initialise
    coinc_ind%grid = no_value
    coinc_ind%i    = 0

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
      coinc_ind%grid = a_grid
      coinc_ind%i    = via
    elseif (norm2( pb - p) < mesh%tol_dist) then
      ! p lies on vib
      coinc_ind%grid = a_grid
      coinc_ind%i    = vib
    elseif (norm2( pc - p) < mesh%tol_dist) then
      ! p lies on vic
      coinc_ind%grid = a_grid
      coinc_ind%i    = vic

    ! Check if p lies on any of the three edges
    elseif (lies_on_line_segment( pa, pb, p, mesh%tol_dist)) then
      ! p lies on the edge connecting via and vib

      coinc_ind%grid = c_grid
      do vvi = 1, mesh%nC( via)
        vj = mesh%C(  via,vvi)
        if (vj == vib) then
          coinc_ind%i = mesh%VE( via,vvi)
        end if
      end do
#if (DO_ASSERTIONS)
      ! Safety
      call assert( coinc_ind%i > 0, &
        'trace_line_tri_start - couldnt find edge ei connecting vertices vi,vj')
#endif

    elseif (lies_on_line_segment( pb, pc, p, mesh%tol_dist)) then
      ! p lies on the edge connecting vib and vic

      coinc_ind%grid = c_grid
      do vvi = 1, mesh%nC( vib)
        vj = mesh%C(  vib,vvi)
        if (vj == vic) then
          coinc_ind%i = mesh%VE( vib,vvi)
          exit
        end if
      end do
#if (DO_ASSERTIONS)
      ! Safety
      call assert( coinc_ind%i > 0, &
        'trace_line_tri_start - couldnt find edge ei connecting vertices vi,vj')
#endif

    elseif (lies_on_line_segment( pc, pa, p, mesh%tol_dist)) then
      ! p lies on the edge connecting vic and via

      coinc_ind%grid = c_grid
      do vvi = 1, mesh%nC( vic)
        vj = mesh%C(  vic,vvi)
        if (vj == via) then
          coinc_ind%i = mesh%VE( vic,vvi)
          exit
        end if
      end do
#if (DO_ASSERTIONS)
      ! Safety
      call assert( coinc_ind%i > 0, &
        'trace_line_tri_start - couldnt find edge ei connecting vertices vi,vj')
#endif

    ! Check if p lies inside triangle ti
    elseif (is_in_triangle( pa, pb, pc, p)) then
      coinc_ind%grid = b_grid
      coinc_ind%i    = ti_hint

    else
      call crash('trace_line_tri_start - couldnt find where p is on the mesh')
    end if

  end subroutine trace_line_tri_start

  !> Given the line [pq], where p lies inside triangle ti_in,
  !> find the point p_next where [pq] crosses into the next triangle.
  subroutine trace_line_tri_ti( mesh, p, q, p_next, coinc_ind, ti_left, coincides, finished)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished

    ! Local variables:
    integer                :: ti_in, via, vib, vic
    real(dp), dimension(2) :: pa, pb, pc
    logical                :: q_in_ti, q_on_vi, q_on_ei, pq_through_vi, pq_through_ei

#if (DO_ASSERTIONS)
    ! Safety
    call assert( coinc_ind%grid == b_grid, 'trace_line_tri_ti - coincidence grid is not b_grid')
    call assert( test_ge_le( coinc_ind%i, 1, mesh%nTri), 'trace_line_tri_ti - invalid value for ti')
#endif

    ti_in = coinc_ind%i

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_in,1)
    vib = mesh%Tri( ti_in,2)
    vic = mesh%Tri( ti_in,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

#if (DO_ASSERTIONS)
    ! Safety
    call assert( is_in_triangle( pa, pb, pc, p), 'trace_line_tri_ti - p does not lie in triangle ti')
#endif

    ! Check if q lies inside the same triangle ti
    call trace_line_tri_ti_q_in_ti( q, ti_in, pa, pb, pc, &
      p_next, coinc_ind, ti_left, coincides, finished, q_in_ti)
    if (q_in_ti) return

    ! Check if q lies on any of the three vertices spanning triangle ti
    call trace_line_tri_ti_q_on_vi( mesh, q, ti_in, &
      p_next, coinc_ind, ti_left, coincides, finished, q_on_vi)
    if (q_on_vi) return

    ! Check if q lies on any of the three edges spanning triangle ti
    call trace_line_tri_ti_q_on_ei( mesh, q, ti_in, &
      p_next, coinc_ind, ti_left, coincides, finished, q_on_ei)
    if (q_on_ei) return

    ! Check if pq exits triangle ti through one of its vertices
    call trace_line_tri_ti_pq_through_vi( mesh, p, q, ti_in, &
      p_next, coinc_ind, ti_left, coincides, finished, pq_through_vi)
    if (pq_through_vi) return

    ! Check if pq exits triangle ti through one of its edges
    call trace_line_tri_ti_pq_through_ei( mesh, p, q, ti_in, &
      p_next, coinc_ind, ti_left, coincides, finished, pq_through_ei)
    if (pq_through_ei) return

    if (.not. (q_in_ti .or. q_on_vi .or. q_on_ei .or. pq_through_vi .or. pq_through_ei)) then
      call crash('trace_line_tri_ti - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_tri_ti

  subroutine trace_line_tri_ti_q_in_ti( q, ti_in, pa, pb, pc, &
    p_next, coinc_ind, ti_left, coincides, finished, q_in_ti)

    ! In/output variables
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: ti_in
    real(dp), dimension(2),    intent(in)    :: pa, pb, pc
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_in_ti

    q_in_ti = .false.

    if (is_in_triangle( pa, pb, pc, q)) then
      ! q lies inside triangle ti
      p_next         = q
      ti_left        = ti_in
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coincides      = .false.
      finished       = .true.
      q_in_ti        = .true.
    end if

  end subroutine trace_line_tri_ti_q_in_ti

  subroutine trace_line_tri_ti_q_on_vi( mesh, q, ti_in, &
    p_next, coinc_ind, ti_left, coincides, finished, q_on_vi)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: ti_in
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_vi

    ! Local variables
    integer :: n, vi

    q_on_vi = .false.

    do n = 1, 3
      vi = mesh%Tri( ti_in,n)
      if (norm2( mesh%V( vi,:) - q) < mesh%tol_dist) then
        ! q lies on vertex vi
        p_next         = q
        ti_left        = ti_in
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        coincides      = .false.
        finished       = .true.
        q_on_vi        = .true.
        exit
      end if
    end do

  end subroutine trace_line_tri_ti_q_on_vi

  subroutine trace_line_tri_ti_q_on_ei( mesh, q, ti_in, &
    p_next, coinc_ind, ti_left, coincides, finished, q_on_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: ti_in
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_ei

    ! Local variables
    integer                :: n1, n2
    real(dp), dimension(2) :: pa, pb

    q_on_ei = .false.

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      pa = mesh%V( mesh%Tri( ti_in,n1),:)
      pb = mesh%V( mesh%Tri( ti_in,n2),:)
      if (lies_on_line_segment( pa, pb, q, mesh%tol_dist)) then
        ! q lies on edge via-vib
        p_next         = q
        ti_left        = ti_in
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        coincides      = .false.
        finished       = .true.
        q_on_ei        = .true.
        exit
      end if
    end do

  end subroutine trace_line_tri_ti_q_on_ei

  subroutine trace_line_tri_ti_pq_through_vi( mesh, p, q, ti_in, &
    p_next, coinc_ind, ti_left, coincides, finished, pq_through_vi)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: ti_in
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_vi

    ! Local variables
    integer                :: n, vi
    real(dp), dimension(2) :: pvi

    pq_through_vi = .false.

    do n = 1, 3
      vi = mesh%Tri( ti_in,n)
      pvi = mesh%V( vi,:)
      if (lies_on_line_segment( p, q, pvi, mesh%tol_dist)) then
        ! [pq] passes through vi
        p_next         = pvi
        ti_left        = ti_in
        coinc_ind%grid = a_grid
        coinc_ind%i    = vi
        coincides      = .false.
        finished       = .false.
        pq_through_vi  = .true.
        exit
      end if
    end do

  end subroutine trace_line_tri_ti_pq_through_vi

  subroutine trace_line_tri_ti_pq_through_ei( mesh, p, q, ti_in, &
    p_next, coinc_ind, ti_left, coincides, finished, pq_through_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: ti_in
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ei

    ! Local variables
    integer                :: n1, n2, vi, vj, ci, vk
    real(dp), dimension(2) :: pvi, pvj, llis
    logical                :: do_cross

    pq_through_ei = .false.

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1

      vi = mesh%Tri( ti_in,n1)
      vj = mesh%Tri( ti_in,n2)

      pvi = mesh%V( vi,:)
      pvj = mesh%V( vj,:)

      call segment_intersection( p, q, pvi, pvj, llis, do_cross, mesh%tol_dist)

      if (do_cross) then
        ! [pq] crosses edge [via,vib]
        coinc_ind%grid = c_grid
        ! Find the edge connecting via and vib
        do ci = 1, mesh%nC( vi)
          vk = mesh%C(  vi,ci)
          if (vk == vj) then
            coinc_ind%i = mesh%VE( vi,ci)
            exit
          end if
        end do
        p_next    = llis
        ti_left   = ti_in
        coincides = .false.
        finished  = .false.
        pq_through_ei = .true.
        exit

      end if
    end do

  end subroutine trace_line_tri_ti_pq_through_ei

  !> Given the line [pq], where p lies on vertex vi_on,
  !> find the point p_next where [pq] crosses into the next triangle.
  subroutine trace_line_tri_vi( mesh, p, q, &
    p_next, coinc_ind, ti_left, coincides, finished)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished

    ! Local variables:
    logical :: q_on_ei, q_in_ti, pq_through_vi, pq_through_ei

#if (DO_ASSERTIONS)
    ! Safety
    call assert( coinc_ind%grid == a_grid, 'trace_line_tri_vi - coincidence grid is not a_grid')
    call assert( test_ge_le( coinc_ind%i, 1, mesh%nV), 'trace_line_tri_vi - invalid value for vi')
    call assert( test_le( norm2( p - mesh%V( coinc_ind%i,:)), mesh%tol_dist), &
      'trace_line_tri_vi - p doesnt lie on vertex vi')
#endif

    ! Check if q lies on any of the edges originating in this vertex
    call trace_line_tri_vi_q_on_ei( mesh, p, q, &
      p_next, coinc_ind, ti_left, coincides, finished, q_on_ei)
    if (q_on_ei) return

    ! Check if q lies inside any of the triangles surrounding vi_on
    call trace_line_tri_vi_q_in_ti( mesh, q, &
      p_next, coinc_ind, ti_left, coincides, finished, q_in_ti)
    if (q_in_ti) return

    ! Check if [pq] passes through any of the neighbouring vertices
    call trace_line_tri_vi_pq_through_vi( mesh, p, q, &
      p_next, coinc_ind, ti_left, coincides, finished, pq_through_vi)
    if (pq_through_vi) return

    ! Check if [pq] exits into any of the adjacent triangles
    call trace_line_tri_vi_pq_through_ei( mesh, p, q, &
      p_next, coinc_ind, ti_left, coincides, finished, pq_through_ei)
    if (pq_through_ei) return

    if (.not. (q_on_ei .or. q_in_ti.or. pq_through_vi .or. pq_through_ei)) then
      call crash('trace_line_tri_vi - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_tri_vi

  subroutine trace_line_tri_vi_q_on_ei( mesh, p, q, &
    p_next, coinc_ind, ti_left, coincides, finished, q_on_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_ei

    ! Local variables
    integer                :: vi_on, ci, vj, ei
    real(dp), dimension(2) :: pvi

    q_on_ei = .false.

    vi_on = coinc_ind%i

    do ci = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,ci)
      ei = mesh%VE( vi_on,ci)
      pvi = mesh%V( vj,:)
      if (norm2( mesh%V( vj,:) - q) < mesh%tol_dist .or. &
          lies_on_line_segment( p, pvi, q, mesh%tol_dist)) then
        ! q lies on edge ei, connecting vi_on and vj
        if (mesh%EV( ei,1) == vi_on) then
          ti_left = mesh%ETri( ei,1)
        else
          ti_left = mesh%ETri( ei,2)
        end if
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        coincides      = .true.
        finished       = .true.
        q_on_ei        = .true.
        exit
      end if
    end do

  end subroutine trace_line_tri_vi_q_on_ei

  subroutine trace_line_tri_vi_q_in_ti( mesh, q, &
    p_next, coinc_ind, ti_left, coincides, finished, q_in_ti)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_in_ti

    ! Local variables
    integer                :: vi_on, iti, ti, via, vib, vic
    real(dp), dimension(2) :: pa, pb, pc

    q_in_ti = .false.

    vi_on = coinc_ind%i

    do iti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,iti)
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
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        ti_left        = ti
        coincides      = .false.
        finished       = .true.
        q_in_ti        = .true.
        exit
      end if
    end do

  end subroutine trace_line_tri_vi_q_in_ti

  subroutine trace_line_tri_vi_pq_through_vi( mesh, p, q, &
    p_next, coinc_ind, ti_left, coincides, finished, pq_through_vi)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_vi

    ! Local variables
    integer                :: vi_on, ci, vj, ei
    real(dp), dimension(2) :: pv

    pq_through_vi = .false.

    vi_on = coinc_ind%i

    do ci = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,ci)
      ei = mesh%VE( vi_on,ci)
      pv = mesh%V( vj,:)
      if (lies_on_line_segment( p, q, pv, mesh%tol_dist)) then
        ! [pq] passes through neighbouring vertex vj, which is connected to vi_on by edge ei
        p_next         = pv
        if (mesh%EV( ei,1) == vi_on) then
          ti_left = mesh%ETri( ei,1)
        else
          ti_left = mesh%ETri( ei,2)
        end if
        coinc_ind%grid = a_grid
        coinc_ind%i    = vj
        coincides      = .true.
        finished       = .false.
        pq_through_vi  = .true.
        exit
      end if
    end do

  end subroutine trace_line_tri_vi_pq_through_vi

  subroutine trace_line_tri_vi_pq_through_ei( mesh, p, q, &
    p_next, coinc_ind, ti_left, coincides, finished, pq_through_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ei

    ! Local variables
    integer                :: vi_on, vti, ti, n1, n2, n3, vib, vic, ei, vvi, vj
    real(dp), dimension(2) :: pb, pc, llis
    logical                :: do_cross

    pq_through_ei = .false.

    vi_on = coinc_ind%i

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
            p_next         = llis
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei
            ti_left        = ti
            coincides      = .false.
            finished       = .false.
            pq_through_ei  = .true.
            return
          end if
        end if
      end do
    end do

  end subroutine trace_line_tri_vi_pq_through_ei

  !> Given the line [pq], where p lies on edge ei,
  !> find the point p_next where [pq] crosses into the next triangle.
  subroutine trace_line_tri_ei( mesh, p, q, &
    p_next, coinc_ind, ti_left, coincides, finished)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished

    ! Local variables:
    integer                :: ei_on, via, vib, vil, vir, til, tir
    real(dp), dimension(2) :: pa, pb, pl, pr
    logical                :: q_on_ei, q_in_ti, pq_through_vi, pq_through_ei

#if (DO_ASSERTIONS)
    call assert( coinc_ind%grid == c_grid, 'trace_line_tri_ei - coincidence grid is not c_grid')
    call assert( test_ge_le( coinc_ind%i, 1, mesh%nE), 'trace_line_tri_ei - invalid value for ei')
#endif

    ei_on = coinc_ind%i

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

#if (DO_ASSERTIONS)
    call assert( lies_on_line_segment( pa, pb, p, mesh%tol_dist), &
      'trace_line_tri_ei - p does not lie on edge ei')
#endif

    ! Check if q lies on the same edge as p
    call trace_line_tri_ei_q_on_ei( mesh, p, q, til, tir, pa, pb, &
      p_next, coinc_ind, ti_left, coincides, finished, q_on_ei)
    if (q_on_ei) return

    ! Check if q lies inside either of the two adjacent triangles
    call trace_line_tri_ei_q_in_ti( mesh, q, til, tir, pa, pb, pr, pl, &
      p_next, coinc_ind, ti_left, coincides, finished, q_in_ti)
    if (q_in_ti) return

    ! Check if pq passes through any of the (up to) four vertices surrounding edge ei
    call trace_line_tri_ei_pq_through_vi( mesh, p, q, via, vib, vil, vir, til, tir, &
      pa, pb, pr, pl, p_next, coinc_ind, ti_left, coincides, finished, pq_through_vi)
    if (pq_through_vi) return

    ! Check if pq passes through any of the (up to) four edges surrounding edge ei
    call trace_line_tri_ei_pq_through_ei( mesh, p, q, via, vib, vil, vir, til, tir, &
      pa, pb, pr, pl, p_next, coinc_ind, ti_left, coincides, finished, pq_through_ei)
    if (pq_through_ei) return

    if (.not. (q_on_ei .or. q_in_ti .or. pq_through_vi .or. pq_through_ei)) then
      call crash('trace_line_tri_ei - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_tri_ei

  subroutine trace_line_tri_ei_q_on_ei( mesh, p, q, til, tir, pa, pb, &
    p_next, coinc_ind, ti_left, coincides, finished, q_on_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: til, tir
    real(dp), dimension(2),    intent(in)    :: pa, pb
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_ei

    q_on_ei = .false.

    ! Check if q lies on the same edge in the direction of via
    if (lies_on_line_segment( p, pa, q, mesh%tol_dist) .or. &
      norm2( q - pa) < mesh%tol_dist) then
      ! q lies on the same edge in the direction of via
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      ti_left        = tir
      coincides      = .true.
      finished       = .true.
      q_on_ei        = .true.
      return
    end if

    ! Check if q lies on the same edge in the direction of vib
    if (lies_on_line_segment( p, pb, q, mesh%tol_dist) .or. &
      norm2( q - pb) < mesh%tol_dist) then
      ! q lies on the same edge in the direction of vib
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      ti_left        = til
      coincides      = .true.
      finished       = .true.
      q_on_ei        = .true.
      return
    end if

  end subroutine trace_line_tri_ei_q_on_ei

  subroutine trace_line_tri_ei_q_in_ti( mesh, q, til, tir, pa, pb, pr, pl, &
    p_next, coinc_ind, ti_left, coincides, finished, q_in_ti)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: til, tir
    real(dp), dimension(2),    intent(in)    :: pa, pb, pr, pl
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_in_ti

    q_in_ti = .false.

    if (til > 0) then
      if (is_in_triangle( pa, pb, pl, q) .or. &
          lies_on_line_segment( pa, pb, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pb, pl, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pl, pa, q, mesh%tol_dist)) then
        ! q lies inside triangle til
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        ti_left        = til
        coincides      = .false.
        finished       = .true.
        q_in_ti        = .true.
        return
      end if
    end if

    if (tir > 0) then
      if (is_in_triangle( pa, pr, pb, q) .or. &
          lies_on_line_segment( pa, pr, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pr, pb, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pb, pa, q, mesh%tol_dist)) then
        ! q lies inside triangle tir
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        ti_left        = tir
        coincides      = .false.
        finished       = .true.
        q_in_ti        = .true.
        return
      end if
    end if

  end subroutine trace_line_tri_ei_q_in_ti

  subroutine trace_line_tri_ei_pq_through_vi( mesh, p, q, via, vib, vil, vir, til, tir, &
    pa, pb, pr, pl, p_next, coinc_ind, ti_left, coincides, finished, pq_through_vi)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: via, vib, vil, vir, til, tir
    real(dp), dimension(2),    intent(in)    :: pa, pb, pr, pl
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_vi

    pq_through_vi = .false.

    ! Check if [pq] passes through pa
    if (lies_on_line_segment( p, q, pa, mesh%tol_dist)) then
      ! [pq] passes through pa
      p_next         = pa
      coinc_ind%grid = a_grid
      coinc_ind%i    = via
      ti_left        = tir
      coincides      = .true.
      finished       = .false.
      pq_through_vi  = .true.
      return
    end if

    ! Check if [pq] passes through pb
    if (lies_on_line_segment( p, q, pb, mesh%tol_dist)) then
      ! [pq] passes through pb
      p_next         = pb
      coinc_ind%grid = a_grid
      coinc_ind%i    = vib
      ti_left        = til
      coincides      = .true.
      finished       = .false.
      pq_through_vi  = .true.
      return
    end if

    ! Check if [pq] passes through pl
    if (til > 0) then
      if (lies_on_line_segment( p, q, pl, mesh%tol_dist)) then
        ! [pq] passes through pl
        p_next         = pl
        coinc_ind%grid = a_grid
        coinc_ind%i    = vil
        ti_left        = til
        coincides      = .false.
        finished       = .false.
        pq_through_vi  = .true.
        return
      end if
    end if

    ! Check if [pq] passes through pr
    if (tir > 0) then
      if (lies_on_line_segment( p, q, pr, mesh%tol_dist)) then
        ! [pq] passes through pr
        p_next         = pr
        coinc_ind%grid = a_grid
        coinc_ind%i    = vir
        ti_left        = tir
        coincides      = .false.
        finished       = .false.
        pq_through_vi  = .true.
        return
      end if
    end if

  end subroutine trace_line_tri_ei_pq_through_vi

  subroutine trace_line_tri_ei_pq_through_ei( mesh, p, q, via, vib, vil, vir, til, tir, &
    pa, pb, pr, pl, p_next, coinc_ind, ti_left, coincides, finished, pq_through_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: via, vib, vil, vir, til, tir
    real(dp), dimension(2),    intent(in)    :: pa, pb, pr, pl
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: ti_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ei

    ! Local variables
    integer                :: ci, vj
    logical                :: do_cross
    real(dp), dimension(2) :: llis

    pq_through_ei = .false.

    ! Check if [pq] crosses edge [via,vil]
    if (til > 0) then
      call segment_intersection( p, q, pa, pl, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [via,vil]
        coinc_ind%grid = c_grid
        ! Find the edge connecting via and vil
        do ci = 1, mesh%nC( via)
          vj = mesh%C( via,ci)
          if (vj == vil) then
            coinc_ind%i = mesh%VE( via,ci)
            exit
          end if
        end do
        p_next        = llis
        ti_left       = til
        coincides     = .false.
        finished      = .false.
        pq_through_ei = .true.
        return
      end if
    end if

    ! Check if [pq] crosses edge [vil,vib]
    if (til > 0) then
      call segment_intersection( p, q, pl, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [vil,vib]
        coinc_ind%grid = c_grid
        ! Find the edge connecting vil and vib
        do ci = 1, mesh%nC( vil)
          vj = mesh%C( vil,ci)
          if (vj == vib) then
            coinc_ind%i = mesh%VE( vil,ci)
            exit
          end if
        end do
        p_next        = llis
        ti_left       = til
        coincides     = .false.
        finished      = .false.
        pq_through_ei = .true.
        return
      end if
    end if

    ! Check if [pq] crosses edge [via,vir]
    if (tir > 0) then
      call segment_intersection( p, q, pa, pr, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [via,vir]
        coinc_ind%grid = c_grid
        ! Find the edge connecting via and vir
        do ci = 1, mesh%nC( via)
          vj  = mesh%C( via,ci)
          if (vj == vir) then
            coinc_ind%i = mesh%VE( via,ci)
            exit
          end if
        end do
        p_next        = llis
        ti_left       = tir
        coincides     = .false.
        finished      = .false.
        pq_through_ei = .true.
        return
      end if
    end if

    ! Check if [pq] crosses edge [vir,vib]
    if (tir > 0) then
      call segment_intersection( p, q, pr, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [vir,vib]
        coinc_ind%grid = c_grid
        ! Find the edge connecting vir and vib
        do ci = 1, mesh%nC( vir)
          vj = mesh%C( vir,ci)
          if (vj == vib) then
            coinc_ind%i = mesh%VE( vir,ci)
            exit
          end if
        end do
        p_next        = llis
        ti_left       = tir
        coincides     = .false.
        finished      = .false.
        pq_through_ei = .true.
        return
      end if
    end if

  end subroutine trace_line_tri_ei_pq_through_ei

end module line_tracing_triangles
