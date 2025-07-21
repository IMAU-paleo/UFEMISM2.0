module line_tracing_grid

  ! Line tracing algorithm through square grid cells

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use remapping_types, only: type_single_row_mapping_matrices
  use line_tracing_basic
  use grid_types, only: type_grid
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use plane_geometry, only: lies_on_line_segment, segment_intersection, crop_line_to_domain
  use line_integrals, only: line_integral_xdy, line_integral_mxydx, line_integral_xydy

  implicit none

  private

  public :: trace_line_grid
  public :: trace_line_grid_start, trace_line_grid_a, trace_line_grid_b, trace_line_grid_cx, trace_line_grid_cy

contains

  !> Trace the line [pq] through the grid and calculate the three
  !> line integrals for the line segments inside the different grid cells.
  subroutine trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! In/output variables
    type(type_grid),                        intent(in)    :: grid
    real(dp), dimension(2),                 intent(in)    :: p,q
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    logical,                                intent(in)    :: count_coincidences

    ! Local variables:
    real(dp), dimension(2)    :: pp,qq
    logical                   :: pq_passes_through_domain
    logical                   :: finished
    integer                   :: n_cycles
    type(type_coinc_ind_grid) :: coinc_ind
    real(dp), dimension(2)    :: p_next
    integer                   :: n_left
    logical                   :: coincides
    real(dp)                  :: LI_xdy, LI_mxydx, LI_xydy

    ! Crop the line [pq] so that it lies within the domain
    call crop_line_to_domain( p, q, grid%xmin, grid%xmax, grid%ymin, grid%ymax, grid%tol_dist, &
      pp, qq, pq_passes_through_domain)

    if (.not. pq_passes_through_domain) then
      ! [pq] doesn't pass through the domain anywhere
      return
    end if

    ! Initialise the coincidence indicator for the point p
    call trace_line_grid_start( grid, pp, coinc_ind)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      select case (coinc_ind%grid)
      case default
        call crash('trace_line_grid - coincidence indicator doesnt make sense')
      case (no_value)
        call crash('trace_line_grid - coincidence indicator doesnt make sense')
      case (a_grid)
        call trace_line_grid_a( grid, pp, qq, coinc_ind, p_next, n_left, coincides, finished)
      case (b_grid)
        call trace_line_grid_b( grid, pp, qq, coinc_ind, p_next, n_left, coincides, finished)
      case (cx_grid)
        call trace_line_grid_cx( grid, pp, qq, coinc_ind, p_next, n_left, coincides, finished)
      case (cy_grid)
        call trace_line_grid_cy( grid, pp, qq, coinc_ind, p_next, n_left, coincides, finished)
      end select

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, grid%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, grid%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, grid%tol_dist)

      ! Add them to the results structure
      if (norm2( p_next - pp) > grid%tol_dist) then
        call add_integrals_to_single_row( single_row, n_left, &
          LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      end if

      ! Cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      if (n_cycles > grid%n) then
        call crash('trace_line_grid - iterative tracer got stuck!')
      end if

    end do

  end subroutine trace_line_grid

  !> Initialise the coincidence indicator for the point p
  subroutine trace_line_grid_start( grid, p, coinc_ind)

    ! In/output variables
    type(type_grid),           intent(in)  :: grid
    real(dp), dimension(2),    intent(in)  :: p
    type(type_coinc_ind_grid), intent(out) :: coinc_ind

    ! Local variables:
    integer  :: i,j
    real(dp) :: xl,xu,yl,yu

#if (DO_ASSERTIONS)
    ! Safety - check if p lies inside the grid domain
    call assert( &
      test_ge_le( p(1), grid%xmin - grid%dx / 2._dp, grid%xmax + grid%dx / 2._dp) .and. &
      test_ge_le( p(2), grid%ymin - grid%dx / 2._dp, grid%ymax + grid%dx / 2._dp), &
      'trace_line_grid_start - p lies outside the grid domain')
#endif

    ! Initialise
    coinc_ind%grid = no_value
    coinc_ind%i = 0
    coinc_ind%j = 0

    ! Find the grid cell containing p
    i = max( 1, min( grid%nx, 1 + floor( (p(1) - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
    j = max( 1, min( grid%ny, 1 + floor( (p(2) - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

#if (DO_ASSERTIONS)
    call assert( test_ge_le( p(1), xl - grid%tol_dist, xu + grid%tol_dist) .and. test_ge_le( p(2), yl - grid%tol_dist, yu + grid%tol_dist), &
      'trace_line_grid_start - couldnt find grid cell containing p')
#endif

    ! Check if p lies on either of the four surrounding b-grid points
    if (test_tol( p(1), xl, grid%tol_dist) .and. test_tol( p(2), yl, grid%tol_dist)) then
      ! p lies on the southwest corner of a-grid cell [i,j] (i.e. b-grid point [i-1,j-1])
      coinc_ind%grid = b_grid
      coinc_ind%i = i-1
      coinc_ind%j = j-1
    elseif (test_tol( p(1), xl, grid%tol_dist) .and. test_tol( p(2), yu, grid%tol_dist)) then
      ! p lies on the northwest corner of a-grid cell [i,j] (i.e. b-grid point [i-1,j])
      coinc_ind%grid = b_grid
      coinc_ind%i = i-1
      coinc_ind%j = j
    elseif (test_tol( p(1), xu, grid%tol_dist) .and. test_tol( p(2), yl, grid%tol_dist)) then
      ! p lies on the southeast corner of a-grid cell [i,j] (i.e. b-grid point [i,j-1])
      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j-1
    elseif (test_tol( p(1), xu, grid%tol_dist) .and. test_tol( p(2), yu, grid%tol_dist)) then
      ! p lies on the northeast corner of a-grid cell [i,j] (i.e. b-grid point [i,j])
      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

    ! Check if p lies on any of the four borders
    elseif (test_tol( p(1), xl, grid%tol_dist)) then
      ! p coincides with the western border of a-grid cell [i,j] (i.e. cx-grid point [i-1,j])
      coinc_ind%grid = cx_grid
      coinc_ind%i = i-1
      coinc_ind%j = j
    elseif (test_tol( p(1), xu, grid%tol_dist)) then
      ! p coincides with the eastern border of a-grid cell [i,j] (i.e. cx-grid point [i,j])
      coinc_ind%grid = cx_grid
      coinc_ind%i = i
      coinc_ind%j = j
    elseif (test_tol( p(2), yl, grid%tol_dist)) then
      ! p coincides with the southern border of a-grid cell [i,j] (i.e. cy-grid point [i,j-1])
      coinc_ind%grid = cy_grid
      coinc_ind%i = i
      coinc_ind%j = j-1
    elseif (test_tol( p(2), yu, grid%tol_dist)) then
      ! p coincides with the northern border of a-grid cell [i,j] (i.e. cy-grid point [i,j])
      coinc_ind%grid = cy_grid
      coinc_ind%i = i
      coinc_ind%j = j

    ! p doesn't lie on the corners or borders, so it must lie inside the grid cell
    elseif (test_ge_le( p(1), xl, xu) .and. test_ge_le( p(2), yl, yu)) then
      coinc_ind%grid = a_grid
      coinc_ind%i = i
      coinc_ind%j = j
    else
      call crash('trace_line_grid_start - couldnt find where p is on the grid')
    end if

  end subroutine trace_line_grid_start

  !> Given the line [pq], where p lies inside grid cell aij_in,
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_a( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: xl,xu,yl,yu
    real(dp), dimension(2) :: sw,nw,se,ne
    logical                :: q_in_a, pq_through_b, pq_through_cx_cy

#if (DO_ASSERTIONS)
    ! Safety
    call assert( coinc_ind%grid == a_grid .and. &
      test_ge_le( coinc_ind%i, 1, grid%nx) .and. test_ge_le( coinc_ind%j, 1, grid%ny), &
      'trace_line_grid_a - coincidence indicator doesnt make sense')
#endif

    i = coinc_ind%i
    j = coinc_ind%j

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! This grid cell's corners
    sw = [xl,yl]
    nw = [xl,yu]
    se = [xu,yl]
    ne = [xu,yu]

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_ge_le( p(1), xl, xu) .and. test_ge_le( p(2), yl, yu), &
      'trace_line_grid_a - p doesnt lie in aij')
#endif

    ! Check if q lies inside the same a-grid cell [i,j]
    call trace_line_grid_a_q_in_a( grid, q, i, j, xl, xu, yl, yu, &
      coinc_ind, p_next, n_left, coincides, finished, q_in_a)
    if (q_in_a) return

    ! Check if pq exits a-gridl cell [i,j] through one of its corners (i.e. the b-grid points)
    call trace_line_grid_a_pq_through_b( grid, p, q, i, j, sw, nw, se, ne, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_b)
    if (pq_through_b) return

    ! Check if pq exits a-gridl cell [i,j] through one of its borders (i.e. the cx/cy-grid points)
    call trace_line_grid_a_pq_through_cx_cy( grid, p, q, i, j, sw, nw, se, ne, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_cx_cy)
    if (pq_through_cx_cy) return

    if (.not. (q_in_a .or. pq_through_b .or. pq_through_cx_cy)) then
      call crash('trace_line_grid_a - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_grid_a

  subroutine trace_line_grid_a_q_in_a( grid, q, i, j, xl, xu, yl, yu, &
    coinc_ind, p_next, n_left, coincides, finished, q_in_a)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: i,j
    real(dp),                  intent(in)    :: xl, xu, yl, yu
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, q_in_a

    q_in_a    = .false.
    p_next    = [0._dp, 0._dp]
    n_left    = 0
    coincides = .false.
    finished  = .false.

    ! Check if q lies inside the same grid cell
    if (test_ge_le( q( 1), xl - grid%tol_dist, xu + grid%tol_dist) .and. &
        test_ge_le( q( 2), yl - grid%tol_dist, yu + grid%tol_dist)) then
      ! q lies inside the same grid cell
      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .true.
    end if

  end subroutine trace_line_grid_a_q_in_a

  subroutine trace_line_grid_a_pq_through_b( grid, p, q, i, j, sw, nw, se, ne, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_b)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_b

    pq_through_b = .false.
    p_next       = [0._dp, 0._dp]
    n_left       = 0
    coincides    = .false.
    finished     = .false.

    if (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits a-grid cell [i,j] through the southwest corner (i.e. b-grid point [i-1,j-1])

      pq_through_b   = .true.
      p_next         = sw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits a-grid cell [i,j] through the northwest corner (i.e. b-grid point [i-1,j])

      pq_through_b   = .true.
      p_next         = nw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits a-grid cell [i,j] through the southeast corner (i.e. b-grid point [i,j-1])

      pq_through_b   = .true.
      p_next         = se
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits a-grid cell [i,j] through the northeast corner (i.e. b-grid point [i,j])

      pq_through_b   = .true.
      p_next         = ne
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    end if

  end subroutine trace_line_grid_a_pq_through_b

  subroutine trace_line_grid_a_pq_through_cx_cy( grid, p, q, i, j, sw, nw, se, ne, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_cx_cy)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_cx_cy

    ! Local variables
    real(dp), dimension(2) :: llis
    logical                :: do_cross

    pq_through_cx_cy = .false.
    p_next           = [0._dp, 0._dp]
    n_left           = 0
    coincides        = .false.
    finished         = .false.

    call segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits a-grid cell [i,j] through the western border (i.e. cx-grid point [i-1,j])

      pq_through_cx_cy = .true.
      p_next           = llis
      coinc_ind%grid   = cx_grid
      coinc_ind%i      = i-1
      coinc_ind%j      = j
      n_left           = grid%ij2n( i,j)
      coincides        = .false.
      finished         = .false.
      return

    end if

    call segment_intersection( p, q, se, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits a-grid cell [i,j] through the eastern border (i.e. cx-grid point [i,j])

      pq_through_cx_cy = .true.
      p_next           = llis
      coinc_ind%grid   = cx_grid
      coinc_ind%i      = i
      coinc_ind%j      = j
      n_left           = grid%ij2n( i,j)
      coincides        = .false.
      finished         = .false.
      return

    end if

    call segment_intersection( p, q, sw, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits a-grid cell [i,j] through the southern border (i.e. cy-grid point [i,j-1])

      pq_through_cx_cy = .true.
      p_next           = llis
      coinc_ind%grid   = cy_grid
      coinc_ind%i      = i
      coinc_ind%j      = j-1
      n_left           = grid%ij2n( i,j)
      coincides        = .false.
      finished         = .false.
      return

    end if

    call segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits a-grid cell [i,j] through the northern border (i.e. cy-grid point [i,j])

      pq_through_cx_cy = .true.
      p_next           = llis
      coinc_ind%grid   = cy_grid
      coinc_ind%i      = i
      coinc_ind%j      = j
      n_left           = grid%ij2n( i,j)
      coincides        = .false.
      finished         = .false.
      return

    end if

  end subroutine trace_line_grid_a_pq_through_cx_cy

  !> Given the line [pq], where p lies on b-grid point bij_on
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_b( grid, p, q, &
    coinc_ind, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: x, y, xl, xu, yl, yu
    real(dp), dimension(2) :: sw, nw, se, ne, ww, ee, ss, nn
    logical                :: q_on_cx_cy, q_in_a, pq_through_b, pq_through_a

#if (DO_ASSERTIONS)
    ! Safety
    call assert( coinc_ind%grid == b_grid .and. &
      test_ge_le( coinc_ind%i, 1, grid%nx-1) .and. test_ge_le( coinc_ind%j, 1, grid%ny-1), &
      'trace_line_grid_b - coincidence indicator doesnt make sense!')
#endif

    i = coinc_ind%i
    j = coinc_ind%j

    ! The eight surrounding b-grid points spanning the four surrounding a-grid cells
    x  = grid%x( i) + grid%dx / 2._dp
    y  = grid%y( j) + grid%dx / 2._dp
    xl = x - grid%dx
    xu = x + grid%dx
    yl = y - grid%dx
    yu = y + grid%dx

    sw = [xl,yl]
    ww = [xl,y ]
    nw = [xl,yu]
    ss = [x ,yl]
    nn = [x ,yu]
    se = [xu,yl]
    ee = [xu,y ]
    ne = [xu,yu]

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_tol( p(1), x, grid%tol_dist) .and. test_tol( p(2), y, grid%tol_dist), &
      'trace_line_grid_b - p doesnt lie on bij')
#endif

    ! Check if q lies on any of the four surrounding edges (i.e. cx/cy-grid points)
    call trace_line_grid_b_q_on_cx_cy( grid, q, i, j, &
      x, y, xl, xu, yl, yu, &
      coinc_ind, p_next, n_left, coincides, finished, q_on_cx_cy)
    if (q_on_cx_cy) return

    ! Check if q lies inside any of the four surrounding a-grid cells
    call trace_line_grid_b_q_in_a( grid, q, i, j, &
      x, y, xl, xu, yl, yu, &
      coinc_ind, p_next, n_left, coincides, finished, q_in_a)
    if (q_in_a) return

    ! Check if pq passes through any of the 8 surrounding b-grid points
    call trace_line_grid_b_pq_through_b( grid, p, q, i, j, &
      sw, nw, se, ne, ww, ee, ss, nn, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_b)
    if (pq_through_b) return

    ! Check if pq exits any of the surrounding a-grid cells through their opposite borders
    call trace_line_grid_b_pq_through_a( grid, p, q, i, j, &
      sw, nw, se, ne, ww, ee, ss, nn, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_a)
    if (pq_through_a) return

    if (.not. (q_on_cx_cy .or. q_in_a .or. pq_through_b .or. pq_through_a)) then
      call crash('trace_line_grid_b - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_grid_b

  subroutine trace_line_grid_b_q_on_cx_cy( grid, q, i, j, &
    x, y, xl, xu, yl, yu, &
    coinc_ind, p_next, n_left, coincides, finished, q_on_cx_cy)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: i,j
    real(dp),                  intent(in)    :: x, y, xl, xu, yl, yu
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, q_on_cx_cy

    q_on_cx_cy = .false.
    p_next     = [0._dp, 0._dp]
    n_left     = 0
    coincides  = .false.
    finished   = .false.

    if (test_ge_le( q(1), xl - grid%tol_dist, x + grid%tol_dist) .and. &
        test_tol( q( 2), y, grid%tol_dist)) then
      ! q lies on the edge to the west (i.e. cy-grid point [i,j])

      q_on_cx_cy     = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j)
      coincides      = .true.
      finished       = .true.

    elseif (test_ge_le( q(1), x - grid%tol_dist, xu + grid%tol_dist) .and. &
        test_tol( q( 2), y, grid%tol_dist)) then
      ! q lies on the edge to the east (i.e. cy-grid point [i+1,j])

      q_on_cx_cy     = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i+1,j+1)
      coincides      = .true.
      finished       = .true.

    elseif (test_ge_le( q(2), yl - grid%tol_dist, y + grid%tol_dist) .and. &
        test_tol( q( 1), x, grid%tol_dist)) then
      ! q lies on the edge to the east (i.e. cx-grid point [i,j])

      q_on_cx_cy     = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i+1,j)
      coincides      = .true.
      finished       = .true.

    elseif (test_ge_le( q(2), y - grid%tol_dist, yu + grid%tol_dist) .and. &
        test_tol( q( 1), x, grid%tol_dist)) then
      ! q lies on the edge to the east (i.e. cx-grid point [i,j+1])

      q_on_cx_cy     = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j+1)
      coincides      = .true.
      finished       = .true.

    end if

  end subroutine trace_line_grid_b_q_on_cx_cy

  subroutine trace_line_grid_b_q_in_a( grid, q, i, j, &
    x, y, xl, xu, yl, yu, &
    coinc_ind, p_next, n_left, coincides, finished, q_in_a)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: i,j
    real(dp),                  intent(in)    :: x, y, xl, xu, yl, yu
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, q_in_a

    q_in_a    = .false.
    p_next    = [0._dp, 0._dp]
    n_left    = 0
    coincides = .false.
    finished  = .false.

    if (test_ge_le( q(1), xl - grid%tol_dist, x + grid%tol_dist) .and. &
        test_ge_le( q(2), yl - grid%tol_dist, y + grid%tol_dist)) then
      ! q lies inside the a-grid cell to the southwest

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i ,j  )
      coincides      = .false.
      finished       = .true.

    elseif (test_ge_le( q(1), xl - grid%tol_dist, x  + grid%tol_dist) .and. &
        test_ge_le( q(2), y  - grid%tol_dist, yu + grid%tol_dist)) then
      ! q lies inside the a-grid cell to the northwest

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .true.

    elseif (test_ge_le( q(1), x  - grid%tol_dist, xu + grid%tol_dist) .and. &
            test_ge_le( q(2), yl - grid%tol_dist, y  + grid%tol_dist)) then
      ! q lies inside the a-grid cell to the southeast

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i+1,j  )
      coincides      = .false.
      finished       = .true.

    elseif (test_ge_le( q(1), x - grid%tol_dist, xu + grid%tol_dist) .and. &
            test_ge_le( q(2), y - grid%tol_dist, yu + grid%tol_dist)) then
      ! q lies inside the a-grid cell to the northeast

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i+1,j+1)
      coincides      = .false.
      finished       = .true.

    end if

  end subroutine trace_line_grid_b_q_in_a

  subroutine trace_line_grid_b_pq_through_b( grid, p, q, i, j, &
    sw, nw, se, ne, ww, ee, ss, nn, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_b)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne, ww, ee, ss, nn
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_b

    pq_through_b = .false.
    p_next     = [0._dp, 0._dp]
    n_left     = 0
    coincides  = .false.
    finished   = .false.

    if (lies_on_line_segment( p, q, ww, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the west

      pq_through_b   = .true.
      p_next         = ww
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, nn, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the north

      pq_through_b   = .true.
      p_next         = nn
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, ee, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the east

      pq_through_b   = .true.
      p_next         = ee
      coinc_ind%grid = b_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i+1,j+1)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, ss, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the south

      pq_through_b   = .true.
      p_next         = ss
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i+1,j)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the southwest

      pq_through_b   = .true.
      p_next         = sw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the northwest

      pq_through_b   = .true.
      p_next         = nw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the southeast

      pq_through_b   = .true.
      p_next         = se
      coinc_ind%grid = b_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the northeast

      pq_through_b   = .true.
      p_next         = ne
      coinc_ind%grid = b_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i+1,j+1)
      coincides      = .false.
      finished       = .false.

    end if

  end subroutine trace_line_grid_b_pq_through_b

  subroutine trace_line_grid_b_pq_through_a( grid, p, q, i, j, &
    sw, nw, se, ne, ww, ee, ss, nn, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_a)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne, ww, ee, ss, nn
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_a

    ! Local variables:
    real(dp), dimension(2) :: llis
    logical                :: do_cross

    pq_through_a = .false.
    p_next     = [0._dp, 0._dp]
    n_left     = 0
    coincides  = .false.
    finished   = .false.

    call segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i,j+1] (to the northwest) through
      ! its western border (i.e. cx-grid point [i-1,j+1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i,j+1] (to the northwest) through
      ! its northern border (i.e. cy-grid point [i,j+1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i+1,j+1] (to the northeast) through
      ! its northern border (i.e. cy-grid point [i+1,j+1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i+1,j+1)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i+1,j+1] (to the northeast) through
      ! its eastern border (i.e. cx-grid point [i+1,j+1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i+1,j+1)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i+1,j] (to the southeast) through
      ! its eastern border (i.e. cx-grid point [i+1,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i+1,j] (to the southeast) through
      ! its southern border (i.e. cy-grid point [i+1,j-1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i,j] (to the southwest) through
      ! its southern border (i.e. cy-grid point [i,j-1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! pq exits a-grid cell [i,j] (to the southwest) through
      ! its western border (i.e. cx-grid point [i-1,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

  end subroutine trace_line_grid_b_pq_through_a

  !> Given the line [pq], where p lies on cx-grid edge cxij_on
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_cx( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: x,yl,yu
    real(dp), dimension(2) :: sw,nw,se,ne,ss,nn
    logical                :: q_on_cx, q_in_a, pq_through_b, pq_through_a

#if (DO_ASSERTIONS)
    ! Safety
    call assert( coinc_ind%grid == cx_grid .and. &
      test_ge_le( coinc_ind%i, 1, grid%nx-1) .and. test_ge_le( coinc_ind%j, 1, grid%ny), &
      'trace_line_grid_cx - coincidence indicator doesnt make sense')
#endif

    i = coinc_ind%i
    j = coinc_ind%j

    ! This c-grid edge
    x  = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! The b-grid points spanning the two adjacent a-grid cells
    sw = [x - grid%dx, yl]
    nw = [x - grid%dx, yu]
    ss = [x          , yl]
    nn = [x          , yu]
    se = [x + grid%dx, yl]
    ne = [x + grid%dx, yu]

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_tol( p(1), x, grid%tol_dist) .and. &
      test_ge_le( p(2), yl, yu), &
      'trace_line_grid_cx - p doesnt lie on cxij')
#endif

    ! Check if q lies on the same cx-grid point as p
    call trace_line_grid_cx_q_on_cx( grid, p, q, i, j, x, yl, yu, &
      coinc_ind, p_next, n_left, coincides, finished, q_on_cx)
    if (q_on_cx) return

    ! Check if q lies inside either of the two neighbouring a-grid cells
    call trace_line_grid_cx_q_in_a( grid, q, i, j, x, yl, yu, &
      coinc_ind, p_next, n_left, coincides, finished, q_in_a)
    if (q_in_a) return

    ! Check if pq passes through any of the 6 surrounding b-grid points
    call trace_line_grid_cx_pq_through_b( grid, p, q, i, j, &
      sw, nw, se, ne, ss, nn, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_b)
    if (pq_through_b) return

    ! Check if pq exits either of the two neighbouring a-grid cells through their borders
    call trace_line_grid_cx_pq_through_a( grid, p, q, i, j, &
      sw, nw, se, ne, ss, nn, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_a)
    if (pq_through_a) return

    if (.not. (q_on_cx .or. q_in_a .or. pq_through_b .or. pq_through_a)) then
      call crash('trace_line_grid_cx - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_grid_cx

  subroutine trace_line_grid_cx_q_on_cx( grid, p, q, i, j, x, yl, yu, &
    coinc_ind, p_next, n_left, coincides, finished, q_on_cx)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p, q
    integer,                   intent(in)    :: i,j
    real(dp),                  intent(in)    :: x, yl, yu
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, q_on_cx

    q_on_cx   = .false.
    p_next    = [0._dp, 0._dp]
    n_left    = 0
    coincides = .false.
    finished  = .false.

    if (test_tol( q(1), x, grid%tol_dist) .and. test_ge_le( q(2), yl - grid%tol_dist, p(2))) then
      ! q lies on the same cx-grid cell in the southern direction

      q_on_cx        = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i+1,j)
      coincides      = .true.
      finished       = .true.

    elseif (test_tol( q(1), x, grid%tol_dist) .and. test_ge_le( q(2), p(2), yu + grid%tol_dist)) then
      ! q lies on the same cx-grid cell in the northern direction

      q_on_cx        = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j)
      coincides      = .true.
      finished       = .true.

    end if

  end subroutine trace_line_grid_cx_q_on_cx

  subroutine trace_line_grid_cx_q_in_a( grid, q, i, j, x, yl, yu, &
    coinc_ind, p_next, n_left, coincides, finished, q_in_a)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: i,j
    real(dp),                  intent(in)    :: x, yl, yu
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, q_in_a

    q_in_a   = .false.
    p_next    = [0._dp, 0._dp]
    n_left    = 0
    coincides = .false.
    finished  = .false.

    if (test_ge_le( q(1), x - grid%dx - grid%tol_dist, x + grid%tol_dist) .and. &
        test_ge_le( q(2), yl - grid%tol_dist, yu + grid%tol_dist)) then
      ! q lies inside the grid cell to the west (i.e. a-grid cell [i,j])

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .true.

    elseif (test_ge_le( q(1), x - grid%tol_dist, x + grid%dx + grid%tol_dist) .and. &
            test_ge_le( q(2), yl - grid%tol_dist, yu + grid%tol_dist)) then
      ! q lies inside the grid cell to the east (i.e. a-grid cell [i+1,j])

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .true.

    end if

  end subroutine trace_line_grid_cx_q_in_a

  subroutine trace_line_grid_cx_pq_through_b( grid, p, q, i, j, &
    sw, nw, se, ne, ss, nn, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_b)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p, q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne, ss, nn
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_b

    pq_through_b = .false.
    p_next       = [0._dp, 0._dp]
    n_left       = 0
    coincides    = .false.
    finished     = .false.

    if (lies_on_line_segment( p, q, ss, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the south (i.e. b-grid point [i,j-1])

      pq_through_b   = .true.
      p_next         = ss
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i+1,j)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, nn, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the north (i.e. b-grid point [i,j])

      pq_through_b   = .true.
      p_next         = nn
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits the grid cell to the west (i.e. a-grid cell [i,j])
      ! through the b-grid point to the northwest (i.e. b-grid point [i-1,j])

      pq_through_b   = .true.
      p_next         = nw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits the grid cell to the west (i.e. a-grid cell [i,j])
      ! through the b-grid point to the southwest (i.e. b-grid point [i-1,j-1])

      pq_through_b   = .true.
      p_next         = sw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits the grid cell to the east (i.e. a-grid cell [i+1,j])
      ! through the b-grid point to the northeast (i.e. b-grid point [i+1,j])

      pq_through_b   = .true.
      p_next         = ne
      coinc_ind%grid = b_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits the grid cell to the east (i.e. a-grid cell [i+1,j]
      ! through the b-grid point to the southeast (i.e. b-grid point [i+1,j-1])

      pq_through_b   = .true.
      p_next         = se
      coinc_ind%grid = b_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.

    end if

  end subroutine trace_line_grid_cx_pq_through_b

  subroutine trace_line_grid_cx_pq_through_a( grid, p, q, i, j, &
    sw, nw, se, ne, ss, nn, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_a)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p, q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne, ss, nn
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_a

    ! Local variables
    real(dp), dimension(2) :: llis
    logical                :: do_cross

    pq_through_a = .false.
    p_next       = [0._dp, 0._dp]
    n_left       = 0
    coincides    = .false.
    finished     = .false.

    call segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the west (i.e. a-grid cell [i,j])
      ! through its southern boundary (i.e. cy-grid point [i,j-1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the west (i.e. a-grid cell [i,j])
      ! through its western boundary (i.e. cx-grid point [i-1,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the west (i.e. a-grid cell [i,j])
      ! through its northern boundary (i.e. cy-grid point [i,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the east (i.e. a-grid cell [i+1,j])
      ! through its northern boundary (i.e. cy-grid point [i+1,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, ne, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the east (i.e. a-grid cell [i+1,j])
      ! through its eastern boundary (i.e. cx-grid point [i+1,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the east (i.e. a-grid cell [i+1,j])
      ! through its southern boundary (i.e. cy-grid point [i+1,j-1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i+1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i+1,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

  end subroutine trace_line_grid_cx_pq_through_a

  !> Given the line [pq], where p lies on cy-grid edge cyij_on
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_cy( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p,q
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: xl,xu,y
    real(dp), dimension(2) :: sw,nw,se,ne,ww,ee
    logical                :: q_on_cy, q_in_a, pq_through_b, pq_through_a

#if (DO_ASSERTIONS)
    ! Safety
    call assert( coinc_ind%grid == cy_grid .and. &
      test_ge_le( coinc_ind%i, 1, grid%nx) .and. test_ge_le( coinc_ind%j, 1, grid%ny-1), &
      'trace_line_grid_cy - coincidence indicator doesnt make sense')
#endif

    i = coinc_ind%i
    j = coinc_ind%j

    ! This c-grid edge
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    y  = grid%y( j) + grid%dx / 2._dp

    ! The b-grid points spanning the two adjacent a-grid cells
    sw = [xl, y - grid%dx]
    se = [xu, y - grid%dx]
    ww = [xl, y          ]
    ee = [xu, y          ]
    nw = [xl, y + grid%dx]
    ne = [xu, y + grid%dx]

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_ge_le( p(1), xl, xu) .and. &
      test_tol( p(2), y, grid%tol_dist), &
      'trace_line_grid_cy - p doesnt lie on cyij')
#endif

    ! Check if q lies on the same cy-grid point as p
    call trace_line_grid_cy_q_on_cy( grid, p, q, i, j, y, xl, xu, &
      coinc_ind, p_next, n_left, coincides, finished, q_on_cy)
    if (q_on_cy) return

    ! Check if q lies inside either of the two neighbouring a-grid cells
    call trace_line_grid_cy_q_in_a( grid, q, i, j, y, xl, xu, &
      coinc_ind, p_next, n_left, coincides, finished, q_in_a)
    if (q_in_a) return

    ! Check if pq passes through any of the 6 surrounding b-grid points
    call trace_line_grid_cy_pq_through_b( grid, p, q, i, j, sw, nw, se, ne, ww, ee, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_b)
    if (pq_through_b) return

    ! Check if pq exits either of the two neighbouring a-grid cells through their borders
    call trace_line_grid_cy_pq_through_a( grid, p, q, i, j, sw, nw, se, ne, ww, ee, &
      coinc_ind, p_next, n_left, coincides, finished, pq_through_a)
    if (pq_through_a) return

    if (.not. (q_on_cy .or. q_in_a .or. pq_through_b .or. pq_through_a)) then
      call crash('trace_line_grid_cy - couldnt find out where pq goes from here')
    end if

  end subroutine trace_line_grid_cy

  subroutine trace_line_grid_cy_q_on_cy( grid, p, q, i, j, y, xl, xu, &
    coinc_ind, p_next, n_left, coincides, finished, q_on_cy)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p, q
    integer,                   intent(in)    :: i,j
    real(dp),                  intent(in)    :: y, xl, xu
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, q_on_cy

    q_on_cy   = .false.
    p_next    = [0._dp, 0._dp]
    n_left    = 0
    coincides = .false.
    finished  = .false.

    if (test_tol( q(2), y, grid%tol_dist) .and. test_ge_le( q(1), xl - grid%tol_dist, p(1))) then
      ! q lies on the same cy-grid cell in the western direction

      q_on_cy        = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j)
      coincides      = .true.
      finished       = .true.

    elseif (test_tol( q(2), y, grid%tol_dist) .and. test_ge_le( q(1), p(1), xu + grid%tol_dist)) then
      ! q lies on the same cy-grid cell in the eastern direction

      q_on_cy        = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j+1)
      coincides      = .true.
      finished       = .true.

    end if

  end subroutine trace_line_grid_cy_q_on_cy

  subroutine trace_line_grid_cy_q_in_a( grid, q, i, j, y, xl, xu, &
    coinc_ind, p_next, n_left, coincides, finished, q_in_a)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: q
    integer,                   intent(in)    :: i,j
    real(dp),                  intent(in)    :: y, xl, xu
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, q_in_a

    q_in_a   = .false.
    p_next    = [0._dp, 0._dp]
    n_left    = 0
    coincides = .false.
    finished  = .false.

    if (test_ge_le( q(1), xl - grid%tol_dist, xu + grid%tol_dist) .and. &
        test_ge_le( q(2), y - grid%dx - grid%tol_dist, y + grid%tol_dist)) then
      ! q lies inside the grid cell to the south (i.e. a-grid cell [i,j])

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .true.

    elseif (test_ge_le( q(1), xl - grid%tol_dist, xu + grid%tol_dist) .and. &
            test_ge_le( q(2), y - grid%tol_dist, y + grid%dx + grid%tol_dist)) then
      ! q lies inside the grid cell to the north (i.e. a-grid cell [i,j+1])

      q_in_a         = .true.
      p_next         = q
      coinc_ind%grid = no_value
      coinc_ind%i    = 0
      coinc_ind%j    = 0
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .true.

    end if

  end subroutine trace_line_grid_cy_q_in_a

  subroutine trace_line_grid_cy_pq_through_b( grid, p, q, i, j, sw, nw, se, ne, ww, ee, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_b)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p, q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne, ww, ee
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_b

    pq_through_b = .false.
    p_next       = [0._dp, 0._dp]
    n_left       = 0
    coincides    = .false.
    finished     = .false.

    if (lies_on_line_segment( p, q, ww, grid%tol_dist))  then
      ! [pq] passes through the b-grid point to the west (i.e. b-grid point [i-1,j])

      pq_through_b   = .true.
      p_next         = ww
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, ee, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the east (i.e. b-grid point [i,j])

      pq_through_b   = .true.
      p_next         = ee
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j+1)
      coincides      = .true.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits the grid cell to the north (i.e. a-grid cell [i,j+1]
      ! through the b-grid point to the northwest (i.e. b-grid point [i-1,j+1])

      pq_through_b   = .true.
      p_next         = nw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits the grid cell to the north (i.e. a-grid cell [i,j+1]
      ! through the b-grid point to the northeast (i.e. b-grid point [i,j+1])

      pq_through_b   = .true.
      p_next         = ne
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits the grid cell to the south (i.e. a-grid cell [i,j]
      ! through the b-grid point to the southwest (i.e. b-grid point [i-1,j-1])

      pq_through_b   = .true.
      p_next         = sw
      coinc_ind%grid = b_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    elseif (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits the grid cell to the south (i.e. a-grid cell [i,j]
      ! through the b-grid point to the southeast (i.e. b-grid point [i,j-1])

      pq_through_b   = .true.
      p_next         = se
      coinc_ind%grid = b_grid
      coinc_ind%i    = i
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.

    end if

  end subroutine trace_line_grid_cy_pq_through_b

  subroutine trace_line_grid_cy_pq_through_a( grid, p, q, i, j, sw, nw, se, ne, ww, ee, &
    coinc_ind, p_next, n_left, coincides, finished, pq_through_a)

    ! In/output variables
    type(type_grid),           intent(in)    :: grid
    real(dp), dimension(2),    intent(in)    :: p, q
    integer,                   intent(in)    :: i,j
    real(dp), dimension(2),    intent(in)    :: sw, nw, se, ne, ww, ee
    type(type_coinc_ind_grid), intent(inout) :: coinc_ind
    real(dp), dimension(2),    intent(out)   :: p_next
    integer,                   intent(out)   :: n_left
    logical,                   intent(out)   :: coincides, finished, pq_through_a

    ! Local variables:
    real(dp), dimension(2) :: llis
    logical                :: do_cross

    pq_through_a = .false.
    p_next       = [0._dp, 0._dp]
    n_left       = 0
    coincides    = .false.
    finished     = .false.

    call segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the north (i.e. a-grid cell [i,j+1])
      ! through its western boundary (i.e. cx-grid point [i-1,j+1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the north (i.e. a-grid cell [i,j+1])
      ! through its northern boundary (i.e. cy-grid point [i,j+1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the north (i.e. a-grid cell [i,j+1])
      ! through its eastern boundary (i.e. cx-grid point [i,j+1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i
      coinc_ind%j    = j+1
      n_left         = grid%ij2n( i,j+1)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the south (i.e. a-grid cell [i,j])
      ! through its eastern boundary (i.e. cx-grid point [i,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, se, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the south (i.e. a-grid cell [i,j])
      ! through its southern boundary (i.e. cy-grid point [i,j-1])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cy_grid
      coinc_ind%i    = i
      coinc_ind%j    = j-1
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

    call segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the grid cell to the south (i.e. a-grid cell [i,j])
      ! through its western boundary (i.e. cx-grid point [i-1,j])

      pq_through_a   = .true.
      p_next         = llis
      coinc_ind%grid = cx_grid
      coinc_ind%i    = i-1
      coinc_ind%j    = j
      n_left         = grid%ij2n( i,j)
      coincides      = .false.
      finished       = .false.
      return

    end if

  end subroutine trace_line_grid_cy_pq_through_a

end module line_tracing_grid
