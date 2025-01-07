module line_tracing_Voronoi

  ! Line tracing algorithm through mesh Voronoi cells

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use line_tracing_basic
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use plane_geometry, only: lies_on_line_segment, segment_intersection, crop_line_to_domain
  use line_integrals, only: line_integral_xdy, line_integral_mxydx, line_integral_xydy
  use mesh_utilities, only: find_shared_Voronoi_boundary, is_in_Voronoi_cell, calc_Voronoi_cell, &
    find_containing_vertex

  implicit none

  private

  public :: trace_line_Vor, trace_line_Vor_start, trace_line_Vor_vi, trace_line_Vor_ti, trace_line_Vor_ei

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
    real(dp), dimension(2)    :: pp,qq
    logical                   :: is_valid_line
    logical                   :: finished
    integer                   :: n_cycles
    type(type_coinc_ind_mesh) :: coinc_ind
    real(dp), dimension(2)    :: p_next
    integer                   :: vi_left
    logical                   :: coincides
    real(dp)                  :: LI_xdy, LI_mxydx, LI_xydy

    ! Crop the line [pq] so that it lies within the mesh domain
    call crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    if (.not.is_valid_line) then
      ! [pq] doesn't pass through the mesh domain anywhere
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    call trace_line_Vor_start( mesh, pp, vi_hint, coinc_ind)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not.finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      select case (coinc_ind%grid)
      case default
        call crash('trace_line_Vor - coincidence indicator doesnt make sense')
      case (no_value)
        call crash('trace_line_Vor - coincidence indicator doesnt make sense')
      case (a_grid)
        call trace_line_Vor_vi( mesh, pp, qq, p_next, coinc_ind, vi_left, coincides, finished)
      case (b_grid)
        call trace_line_Vor_ti( mesh, pp, qq, p_next, coinc_ind, vi_left, coincides, finished)
      case (c_grid)
        call trace_line_Vor_ei( mesh, pp, qq, p_next, coinc_ind, vi_left, coincides, finished)
      end select

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

    end do

  end subroutine trace_line_Vor

  !> Initialise the coincidence indicators for the point p, i.e. check if p either...
  !>    - lies inside the Voronoi cell of vertex vi_in, ...
  !>    - lies on the circumcentre of triangle ti_on, or...
  !>    - lies on the shared Voronoi cell boundary represented by edge ei_on
  subroutine trace_line_Vor_start( mesh, p, vi_hint, coinc_ind)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p
    integer,                   intent(inout) :: vi_hint
    type(type_coinc_ind_mesh), intent(out)   :: coinc_ind

    ! Local variables:
    integer                :: iti, ti, vei, ei
    real(dp), dimension(2) :: cc1, cc2

    ! Initialise
    coinc_ind%grid = no_value
    coinc_ind%i    = 0

    ! Find the vertex whose Voronoi cell contains p
    call find_containing_vertex( mesh, p, vi_hint)

    ! Check if p lies on any of the surrounding triangles' circumcentres
    do iti = 1, mesh%niTri( vi_hint)
      ti = mesh%iTri( vi_hint,iti)
      if (norm2( mesh%Tricc( ti,:) - p) < mesh%tol_dist) then
        ! p lies on the circumcentre of triangle ti
        coinc_ind%grid = b_grid
        coinc_ind%i    = ti
        return
      end if
    end do

    ! Check if p lies on any of the shared Voronoi boundaries
    do vei = 1, mesh%nC( vi_hint)
      ei = mesh%VE( vi_hint,vei)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, p, mesh%tol_dist)) then
        ! p lies on the shared Voronoi cell boundary represented by edge ei
        coinc_ind%grid = c_grid
        coinc_ind%i    = ei
        return
      end if
    end do

    ! If p lies not on the boundary of the Voronoi cell, then it must lie inside of it
    coinc_ind%grid = a_grid
    coinc_ind%i    = vi_hint

  end subroutine trace_line_Vor_start

  !> Given the line [pq], where p lies inside the Voronoi cell of vertex vi_in,
  !> find the point p_next where [pq] crosses into the next Voronoi cell.
  subroutine trace_line_Vor_vi(  mesh, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished

    ! Local variables:
    integer                             :: vi_in
    real(dp)                            :: dx
    real(dp), dimension( mesh%nC_mem,2) :: Vor
    integer,  dimension( mesh%nC_mem  ) :: Vor_vi
    integer,  dimension( mesh%nC_mem  ) :: Vor_ti
    integer                             :: nVor
    logical                             :: q_in_vi, q_on_ei, pq_through_ti, pq_through_ei

#if (DO_ASSERTIONS)
    call assert( coinc_ind%grid == a_grid, 'trace_line_Vor_vi - coincidence grid is not a_grid')
    call assert( test_ge_le( coinc_ind%i, 1, mesh%nV), 'trace_line_Vor_vi - invalid value for vi')
    call assert( is_in_Voronoi_cell( mesh, p, coinc_ind%i), &
      'trace_line_Vor_vi - p does not lie in the Voronoi cell of vertex vi')
#endif

    vi_in = coinc_ind%i

    ! Check if q lies inside the same Voronoi cell
    call trace_line_Vor_vi_q_in_vi(  mesh, vi_in, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_in_vi)
    if (q_in_vi) return

    dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 200
    call calc_Voronoi_cell( mesh, vi_in, dx, Vor, Vor_vi, Vor_ti, nVor)

    ! Check if q lies on the boundary of this Voronoi cell
    call trace_line_Vor_vi_q_on_ei(  mesh, vi_in, nVor, Vor, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_on_ei)
    if (q_on_ei) return

    ! Check if [pq] passes through any of the surrounding triangles' circumcentres
    call trace_line_Vor_vi_pq_through_ti(  mesh, vi_in, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, pq_through_ti)
    if (pq_through_ti) return

    ! Check if [pq] passes through any of the shared Voronoi boundaries
    call trace_line_Vor_vi_pq_through_ei(  mesh, vi_in, nVor, Vor, Vor_vi, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, pq_through_ei)
    if (pq_through_ei) return

    if (.not. (q_in_vi .or. q_on_ei .or. pq_through_ti .or. pq_through_ei)) then
      call crash('trace_line_Vor_vi - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_Vor_vi

  subroutine trace_line_Vor_vi_q_in_vi(  mesh, vi_in, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_in_vi)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: vi_in
    real(dp), dimension(2),    intent(in)    :: q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_in_vi

    q_in_vi = .false.

    if (is_in_Voronoi_cell( mesh, q, vi_in)) then
      ! q lies inside the same Voronoi cell
      p_next         = q
      vi_left        = vi_in
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coincides      = .false.
      finished       = .true.
      q_in_vi        = .true.
    end if

  end subroutine trace_line_Vor_vi_q_in_vi

  subroutine trace_line_Vor_vi_q_on_ei(  mesh, vi_in, nVor, Vor, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_on_ei)

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    integer,                             intent(in)    :: vi_in
    integer,                             intent(in)    :: nVor
    real(dp), dimension( mesh%nC_mem,2), intent(in)    :: Vor
    real(dp), dimension(2),              intent(in)    :: q
    real(dp), dimension(2),              intent(out)   :: p_next
    type(type_coinc_ind_mesh),           intent(inout) :: coinc_ind
    integer,                             intent(out)   :: vi_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished
    logical,                             intent(out)   :: q_on_ei

    ! Local variables:
    integer                :: vori1, vori2
    real(dp), dimension(2) :: pa, pb

    q_on_ei = .false.

    do vori1 = 1, nVor
      vori2 = vori1 + 1
      if (vori2 == nVor + 1) vori2 = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori1,:)
      pb = Vor( vori2,:)
      if (norm2( q - pa) < mesh%tol_dist .or. lies_on_line_segment( pa, pb, q, mesh%tol_dist)) then
        ! q lies on the boundary of the same Voronoi cell
        p_next         = q
        vi_left        = vi_in
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        coincides      = .false.
        finished       = .true.
        q_on_ei        = .true.
        exit
      end if
    end do

  end subroutine trace_line_Vor_vi_q_on_ei

  subroutine trace_line_Vor_vi_pq_through_ti(  mesh, vi_in, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, pq_through_ti)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: vi_in
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ti

    ! Local variables
    integer                :: iti, ti
    real(dp), dimension(2) :: r

    pq_through_ti = .false.

    do iti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,iti)
      r  = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, r, mesh%tol_dist)) then
        ! [pq] passes through this triangle's circumcentre
        p_next         = mesh%Tricc( ti,:)
        vi_left        = vi_in
        coinc_ind%grid = b_grid
        coinc_ind%i    = ti
        coincides      = .false.
        finished       = .false.
        pq_through_ti  = .true.
        exit
      end if
    end do

  end subroutine trace_line_Vor_vi_pq_through_ti

  subroutine trace_line_Vor_vi_pq_through_ei(  mesh, vi_in, nVor, Vor, Vor_vi, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, pq_through_ei)

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    integer,                             intent(in)    :: vi_in
    integer,                             intent(in)    :: nVor
    real(dp), dimension( mesh%nC_mem,2), intent(in)    :: Vor
    integer , dimension( mesh%nC_mem  ), intent(in)    :: Vor_vi
    real(dp), dimension(2),              intent(in)    :: p, q
    real(dp), dimension(2),              intent(out)   :: p_next
    type(type_coinc_ind_mesh),           intent(inout) :: coinc_ind
    integer,                             intent(out)   :: vi_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished
    logical,                             intent(out)   :: pq_through_ei

    ! Local variables
    integer                :: vori1, vori2, vj, ei, ci
    real(dp), dimension(2) :: pa, pb, llis
    logical                :: do_cross

    pq_through_ei = .false.

    do vori1 = 1, nVor
      vori2 = vori1 + 1
      if (vori2 == nVor + 1) vori2 = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori1,:)
      pb = Vor( vori2,:)
      ! The other vertex sharing this Voronoi cell boundary
      vj = Vor_vi( vori1)
      ! The edge representing this shared Voronoi cell boundary
      ei = 0
      do ci = 1, mesh%nC( vi_in)
        if (mesh%C( vi_in,ci) == vj) then
          ei = mesh%VE( vi_in,ci)
          exit
        end if
      end do
      ! Safety
      if (ei == 0) call crash('trace_line_Vor_vi_pq_through_ei - couldnt find edge between vi and vj!')

      call segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes into the Voronoi cell of vj
        p_next         = llis
        vi_left        = vi_in
        coinc_ind%grid = c_grid
        coinc_ind%i    = ei
        coincides      = .false.
        finished       = .false.
        pq_through_ei  = .true.
        exit
      end if
    end do

  end subroutine trace_line_Vor_vi_pq_through_ei

  !> Given the line [pq], where p lies on the circumcentre of triangle ti_on,
  !> find the point p_next where [pq] crosses into the next Voronoi cell.
  subroutine trace_line_Vor_ti(  mesh, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished

    ! Local variables:
    integer :: ti_on
    logical :: q_on_ei, q_in_vi, pq_through_ti, pq_through_ei

#if (DO_ASSERTIONS)
    call assert( coinc_ind%grid == b_grid, 'trace_line_Vor_ti - coincidence grid is not b_grid')
    call assert( test_ge_le( coinc_ind%i, 1, mesh%nTri), 'trace_line_Vor_ti - invalid value for ti')
    call assert( norm2( mesh%Tricc( coinc_ind%i,:) - p) <= mesh%tol_dist, &
      'trace_line_Vor_ti - p does not lie on the circumcentre of ti')
#endif

    ti_on = coinc_ind%i

    ! Check if q lies on any of the thre Voronoi cell boundaries
    ! originating at the circumcentre  of triangle ti
    call trace_line_Vor_ti_q_on_ei(  mesh, ti_on, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_on_ei)
    if (q_on_ei) return

    ! Check if q lies inside any of the three adjacent Voronoi cells
    call trace_line_Vor_ti_q_in_vi( mesh, ti_on, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_in_vi)
    if (q_in_vi) return

    ! Check if [pq] passes through the circumcentre of any of the three neighbouring triangles
    call trace_line_Vor_ti_pq_through_ti( mesh, ti_on, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, pq_through_ti)
    if (pq_through_ti) return

    ! Check if [pq] crosses the boundary of one of the adjacent Voronoi cells
    call trace_line_Vor_ti_pq_through_ei( mesh, ti_on, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, pq_through_ei)
    if (pq_through_ei) return

    if (.not. (q_on_ei .or. q_in_vi .or. pq_through_ti .or. pq_through_ei)) then
      call crash('trace_line_Vor_ti - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_Vor_ti

  subroutine trace_line_Vor_ti_q_on_ei( mesh, ti_on, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_on_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: ti_on
    real(dp), dimension(2),    intent(in)    :: q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_ei

    ! Local variables
    integer                :: n1, n2, n3, via, vib, vic, ci, vj, ei
    real(dp), dimension(2) :: cc1, cc2

    q_on_ei = .false.

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      via = mesh%Tri( ti_on,n1)
      vib = mesh%Tri( ti_on,n2)
      vic = mesh%Tri( ti_on,n3)

      do ci = 1, mesh%nC( via)
        vj = mesh%C(  via,ci)
        if (vj == vib) then
          ei = mesh%VE( via,ci)
          exit
        end if
      end do

      ! Check if q lies on the Voronoi cell boundary separating via from vib
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .or. &
          norm2( cc1 - q) < mesh%tol_dist .or. &
          norm2( cc2 - q) < mesh%tol_dist) then
        ! q lies on the Voronoi cell boundary separating via from vib
        if (mesh%ETri( ei,1) == ti_on) then
          vi_left = mesh%EV( ei,2)
        else
          vi_left = mesh%EV( ei,1)
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

  end subroutine trace_line_Vor_ti_q_on_ei

  subroutine trace_line_Vor_ti_q_in_vi( mesh, ti_on, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_in_vi)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: ti_on
    real(dp), dimension(2),    intent(in)    :: q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_in_vi

    ! Local variables
    integer :: n, vi

    q_in_vi = .false.

    do n = 1, 3
      vi = mesh%Tri( ti_on,n)
      if (is_in_Voronoi_cell( mesh, q, vi)) then
        ! q lies inside the Voronoi cell of vertex vi
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        vi_left        = vi
        coincides      = .false.
        finished       = .true.
        q_in_vi        = .true.
        exit
      end if
    end do

  end subroutine trace_line_Vor_ti_q_in_vi

  subroutine trace_line_Vor_ti_pq_through_ti( mesh, ti_on, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, pq_through_ti)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: ti_on
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ti

    ! Local variables
    integer                :: n1, n2, n3, via, vib, vic, tj, iti0, iti1, iti2, tj0, tj1, tj2
    real(dp), dimension(2) :: cc, cc2

    pq_through_ti = .false.

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      via = mesh%Tri( ti_on,n1)
      vib = mesh%Tri( ti_on,n2)
      vic = mesh%Tri( ti_on,n3)

      tj = mesh%TriC( ti_on,n1)

      if (tj > 0) then
        cc = mesh%Tricc( tj,:)
        if (lies_on_line_segment( p, q, cc, mesh%tol_dist)) then
          ! [pq] passes through the circumcentre of this neighbouring triangle
          p_next         = cc
          coinc_ind%grid = b_grid
          coinc_ind%i    = tj
          vi_left        = vic
          coincides      = .true.
          finished       = .false.
          pq_through_ti  = .true.
          exit
        end if
      end if

      ! Check if [pq] passes through the Voronoi vertices spanning the
      ! boundaries of the Voronoi cells adjoining the triangle circumcenter
      ! that p lies upon

      do iti1 = 1, mesh%niTri( via)
        iti0 = iti1 - 1
        if (iti0 == 0) iti0 = mesh%niTri( via)
        iti2 = iti1 + 1
        if (iti2 == mesh%niTri( via) + 1) iti2 = 1
        tj0 = mesh%iTri( via,iti0)
        tj1 = mesh%iTri( via,iti1)
        tj2 = mesh%iTri( via,iti2)
        if (tj0 == ti_on .or. tj1 == ti_on .or. tj2 == ti_on) cycle

        cc2 = mesh%Tricc( tj1,:)
        if (lies_on_line_segment( p, q, cc2, mesh%tol_dist)) then
          ! [pq] passes through this Voronoi vertex (i.e. triangle circumcentre)
          p_next         = cc2
          coinc_ind%grid = b_grid
          coinc_ind%i    = tj1
          vi_left        = via
          coincides      = .false.
          finished       = .false.
          pq_through_ti  = .true.
          return
        end if

      end do
    end do

  end subroutine trace_line_Vor_ti_pq_through_ti

  subroutine trace_line_Vor_ti_pq_through_ei( mesh, ti_on, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, pq_through_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: ti_on
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ei

    ! Local variables
    integer                :: n1, n2, n3, via, vib, vic, ci, vj, ei
    real(dp), dimension(2) :: cc1, cc2, llis
    logical                :: do_cross

    pq_through_ei = .false.

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      via = mesh%Tri( ti_on,n1)
      vib = mesh%Tri( ti_on,n2)
      vic = mesh%Tri( ti_on,n3)

      do ci = 1, mesh%nC( via)
        vj = mesh%C( via,ci)
        if (vj == vib .or. vj == vic) cycle
        ei = mesh%VE( via,ci)
        call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
        call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
        if (do_cross) then
          ! [pq] crosses this part of the boundary of the Voronoi cell of via
          p_next         = llis
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei
          vi_left        = via
          coincides      = .false.
          finished       = .false.
          pq_through_ei  = .true.
          exit
        end if
      end do

    end do

  end subroutine trace_line_Vor_ti_pq_through_ei

  !> Given the line [pq], where p lies on the shared Voronoi boundary represented by edge ei_on,
  !> find the point p_next where [pq] crosses into the next Voronoi cell.
  subroutine trace_line_Vor_ei( mesh, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    real(dp), dimension(2),    intent(in)    :: p,q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished

    ! Local variables:
    integer                :: ei_on, via, vib, vil, vir, til, tir
    real(dp), dimension(2) :: cc1, cc2, ccl, ccr
    logical                :: q_on_ei, q_in_vi, q_on_ti, q_on_other_ei
    logical                :: pq_through_ei, pq_through_ti, pq_through_other_ei

#if (DO_ASSERTIONS)
    call assert( coinc_ind%grid == c_grid, 'trace_line_Vor_ei - coincidence grid is not c_grid')
    call assert( test_ge_le( coinc_ind%i, 1, mesh%nE), 'trace_line_Vor_ei - invalid value for ei')
#endif

    ei_on = coinc_ind%i

    ! Find the endpoints of this shared Voronoi boundary
    call find_shared_Voronoi_boundary( mesh, ei_on, cc1, cc2)

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

#if (DO_ASSERTIONS)
    call assert( lies_on_line_segment( ccl, ccr, p, mesh%tol_dist), &
      'trace_line_Vor_ei - p does not lie on the shared Voronoi boundary represented by ei')
#endif

    ! Check if q lies on the same shared Voronoi cell boundary
    call trace_line_Vor_ei_q_on_ei( mesh, via, vib, ccl, ccr, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_on_ei)
    if (q_on_ei) return

    ! Check if q lies inside either of the two adjacent Voronoi cells
    call trace_line_Vor_ei_q_in_vi( mesh, via, vib, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_in_vi)
    if (q_in_vi) return

    ! Check if q lies on the circumcentres of any of the triangles surrounding
    ! either of the two adjacent vertices
    call trace_line_Vor_ei_q_on_ti( mesh, via, vib, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_on_ti)
    if (q_on_ti) return

    ! Check if q lies on the boundary of the Voronoi cells of
    ! either of the two adjacent vertices
    call trace_line_Vor_ei_q_on_other_ei( mesh, via, vib, q, &
      p_next, coinc_ind, vi_left, coincides, finished, q_on_other_ei)
    if (q_on_other_ei) return

    ! Check if pq passes through either of the two adjacent triangle circumcentres
    call trace_line_Vor_ei_pq_through_ei( mesh, via, vib, til, tir, ccl, ccr, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, pq_through_ei)
    if (pq_through_ei) return

    ! Check if pq crosses the circumcentre of any of the triangles surrounding
    ! either of the two adjacent vertices
    call trace_line_Vor_ei_pq_through_ti( mesh, via, vib, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, pq_through_ti)
    if (pq_through_ti) return

    ! Check if pq crosses the boundary of the Voronoi cells of
    ! either of the two adjacent vertices
    call trace_line_Vor_ei_pq_through_other_ei( mesh, ei_on, via, vib, p, q, &
      p_next, coinc_ind, vi_left, coincides, finished, pq_through_other_ei)
    if (pq_through_other_ei) return

    if (.not. (q_on_ei .or. q_in_vi .or. q_on_ti .or. q_on_other_ei .or. &
      pq_through_ei .or. pq_through_ti .or. pq_through_other_ei)) then
      call crash('trace_line_Vor_ei - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_Vor_ei

  subroutine trace_line_Vor_ei_q_on_ei( mesh, via, vib, ccl, ccr, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_on_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: via, vib
    real(dp), dimension(2),    intent(in)    :: ccl, ccr
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_ei

    q_on_ei = .false.

    ! Check if q lies on the same shared Voronoi cell boundary, in the direction of ccl
    if (lies_on_line_segment( p, ccl, q, mesh%tol_dist) .or. &
      norm2( ccl - q) < mesh%tol_dist) then
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      vi_left        = via
      coincides      = .true.
      finished       = .true.
      q_on_ei        = .true.
      return
    end if

    ! Check if q lies on the same shared Voronoi cell boundary, in the direction of ccr
    if (lies_on_line_segment( p, ccr, q, mesh%tol_dist) .or. &
      norm2( ccr - q) < mesh%tol_dist) then
      ! q coincides with ccr
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      vi_left        = vib
      coincides      = .true.
      finished       = .true.
      q_on_ei        = .true.
      return
    end if

  end subroutine trace_line_Vor_ei_q_on_ei

  subroutine trace_line_Vor_ei_q_in_vi( mesh, via, vib, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_in_vi)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: via, vib
    real(dp), dimension(2),    intent(in)    :: q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_in_vi

    q_in_vi = .false.

    if (is_in_Voronoi_cell( mesh, q, via)) then
      ! q lies inside the Voronoi cell of via
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      vi_left        = via
      coincides      = .false.
      finished       = .true.
      q_in_vi        = .true.
    elseif (is_in_Voronoi_cell( mesh, q, vib)) then
      ! q lies inside the Voronoi cell of vib
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      vi_left        = vib
      coincides      = .false.
      finished       = .true.
      q_in_vi        = .true.
    end if

  end subroutine trace_line_Vor_ei_q_in_vi

  subroutine trace_line_Vor_ei_q_on_ti( mesh, via, vib, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_on_ti)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: via, vib
    real(dp), dimension(2),    intent(in)    :: q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_ti

    ! Local variables
    integer :: iti, ti

    q_on_ti = .false.

    ! Check if q lies on the circumcentre of any of the triangles surrounding via
    do iti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,iti)
      if (norm2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) then
        ! q lies on this triangle's circumcentre
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        vi_left        = via
        coincides      = .false.
        finished       = .true.
        q_on_ti        = .true.
        return
      end if
    end do

    ! Check if q lies on the circumcentre of any of the triangles surrounding vib
    do iti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,iti)
      if (norm2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) then
        ! q lies on this triangle's circumcentre
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        vi_left        = vib
        coincides      = .false.
        finished       = .true.
        q_on_ti        = .true.
        return
      end if
    end do

  end subroutine trace_line_Vor_ei_q_on_ti

  subroutine trace_line_Vor_ei_q_on_other_ei( mesh, via, vib, q, &
    p_next, coinc_ind, vi_left, coincides, finished, q_on_other_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: via, vib
    real(dp), dimension(2),    intent(in)    :: q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: q_on_other_ei

    ! Local variables
    integer                :: ci, ei
    real(dp), dimension(2) :: cc1, cc2

    q_on_other_ei = .false.

    ! Check if q lies on boundary of the Voronoi cell of via
    do ci = 1, mesh%nC( via)
      ei = mesh%VE( via,ci)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) then
        ! q lies on this shared Voronoi boundary
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        vi_left        = via
        coincides      = .false.
        finished       = .true.
        q_on_other_ei  = .true.
        return
      end if
    end do

    ! Check if q lies on boundary of the Voronoi cell of vib
    do ci = 1, mesh%nC( vib)
      ei = mesh%VE( vib,ci)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) then
        ! q lies on this shared Voronoi boundary
        p_next         = q
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        vi_left        = vib
        coincides      = .false.
        finished       = .true.
        q_on_other_ei  = .true.
        return
      end if
    end do

  end subroutine trace_line_Vor_ei_q_on_other_ei

  subroutine trace_line_Vor_ei_pq_through_ei( mesh, via, vib, til, tir, ccl, ccr, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, pq_through_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: via, vib, til, tir
    real(dp), dimension(2),    intent(in)    :: ccl, ccr
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ei

    pq_through_ei = .false.

    if (lies_on_line_segment( p, q, ccl, mesh%tol_dist)) then
      ! [pq] passes through ccl
      p_next         = ccl
      coinc_ind%grid = b_grid
      coinc_ind%i    = til
      vi_left        = via
      coincides      = .true.
      finished       = .false.
      pq_through_ei  = .true.
    elseif (lies_on_line_segment( p, q, ccr, mesh%tol_dist)) then
      ! [pq] passes through ccr
      p_next         = ccr
      coinc_ind%grid = b_grid
      coinc_ind%i    = tir
      vi_left        = vib
      coincides      = .true.
      finished       = .false.
      pq_through_ei  = .true.
    end if

  end subroutine trace_line_Vor_ei_pq_through_ei

  subroutine trace_line_Vor_ei_pq_through_ti( mesh, via, vib, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, pq_through_ti)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: via, vib
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_ti

    ! Local variables
    integer                :: vti, ti
    real(dp), dimension(2) :: cc1

    pq_through_ti = .false.

    ! Check if pq crosses the circumcentre of any of the triangles surrounding via
    do vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      cc1 = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of triangle ti
        p_next         = cc1
        coinc_ind%grid = b_grid
        coinc_ind%i    = ti
        vi_left        = via
        coincides      = .false.
        finished       = .false.
        pq_through_ti  = .true.
        return
      end if
    end do

    ! Check if pq crosses the circumcentre of any of the triangles surrounding vib
    do vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      cc1 = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of triangle ti
        p_next         = cc1
        coinc_ind%grid = b_grid
        coinc_ind%i    = ti
        vi_left        = vib
        coincides      = .false.
        finished       = .false.
        pq_through_ti  = .true.
        return
      end if
    end do

  end subroutine trace_line_Vor_ei_pq_through_ti

  subroutine trace_line_Vor_ei_pq_through_other_ei( mesh, ei_on, via, vib, p, q, &
    p_next, coinc_ind, vi_left, coincides, finished, pq_through_other_ei)

    ! In/output variables
    type(type_mesh),           intent(in)    :: mesh
    integer,                   intent(in)    :: ei_on, via, vib
    real(dp), dimension(2),    intent(in)    :: p, q
    real(dp), dimension(2),    intent(out)   :: p_next
    type(type_coinc_ind_mesh), intent(inout) :: coinc_ind
    integer,                   intent(out)   :: vi_left
    logical,                   intent(out)   :: coincides
    logical,                   intent(out)   :: finished
    logical,                   intent(out)   :: pq_through_other_ei

    ! Local variables
    integer                :: ci, ei
    real(dp), dimension(2) :: cc1, cc2, llis
    logical                :: do_cross

    pq_through_other_ei = .false.

    ! Check if pq crosses the boundary of the Voronoi cell of via
    do ci = 1, mesh%nC( via)
      ei = mesh%VE( via,ci)
      if (ei == ei_on) cycle
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next              = llis
        coinc_ind%grid      = c_grid
        coinc_ind%i         = ei
        vi_left             = via
        coincides           = .false.
        finished            = .false.
        pq_through_other_ei = .true.
        return
      end if
    end do

    ! Check if pq crosses the boundary of the Voronoi cell of vib
    do ci = 1, mesh%nC( vib)
      ei = mesh%VE( vib,ci)
      if (ei == ei_on) cycle
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next              = llis
        coinc_ind%grid      = c_grid
        coinc_ind%i         = ei
        vi_left             = vib
        coincides           = .false.
        finished            = .false.
        pq_through_other_ei = .true.
        return
      end if
    end do

  end subroutine trace_line_Vor_ei_pq_through_other_ei

end module line_tracing_Voronoi
