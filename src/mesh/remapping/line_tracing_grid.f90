module line_tracing_grid

  ! Line tracing algorithm through square grid cells

  use precisions, only: dp
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use line_tracing_basic, only: add_integrals_to_single_row
  use grid_types, only: type_grid
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use math_utilities, only: lies_on_line_segment, segment_intersection, crop_line_to_domain, &
    line_integral_xdy, line_integral_mxydx, line_integral_xydy

  implicit none

  private

  public :: trace_line_grid

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
    character(len=1024), parameter :: routine_name = 'trace_line_Vor'
    real(dp)                       :: xmin, xmax, ymin, ymax
    real(dp), dimension(2)         :: pp,qq
    logical                        :: is_valid_line
    logical                        :: finished
    integer                        :: n_cycles
    integer,  dimension(2)         :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2)         :: p_next
    integer                        :: n_left
    logical                        :: coincides
    real(dp)                       :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    call init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the domain
    xmin = grid%xmin !- grid%dx / 2._dp
    xmax = grid%xmax !+ grid%dx / 2._dp
    ymin = grid%ymin !- grid%dx / 2._dp
    ymax = grid%ymax !+ grid%dx / 2._dp
    call crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, grid%tol_dist, pp, qq, is_valid_line)

    if (.not. is_valid_line) then
      ! [pq] doesn't pass through the domain anywhere
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside grid cell aij_in, ...
    !    - lies on the b-grid point bij_on, or...
    !    - lies on the edge cij_on
    call trace_line_grid_start( grid, pp, aij_in, bij_on, cxij_on, cyij_on)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      if     (aij_in(  1) > 0 .or. aij_in(  2) > 0) then
        ! p lies inside a-grid cell aij_in
        call trace_line_grid_a(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      elseif (bij_on(  1) > 0 .or. bij_on(  2) > 0) then
        ! p lies on b-grid point bij_on
        call trace_line_grid_b(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      elseif (cxij_on( 1) > 0 .or. cxij_on( 2) > 0) then
        ! p lies on cx-grid edge cxij_on
        call trace_line_grid_cx( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      elseif (cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
        ! p lies on cy-grid edge cyij_on
        call trace_line_grid_cy( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      else
        call crash('found no coincidence indicators!')
      end if

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, grid%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, grid%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, grid%tol_dist)

      ! Add them to the results structure
      if (norm2( p_next - pp) > grid%tol_dist) then
        call add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      end if

      ! cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      if (n_cycles > grid%n) then
        call crash('trace_line_grid - iterative tracer got stuck!')
      end if

    end do ! do while (.not. finished)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine trace_line_grid

  !> Initialise the coincidence indicators for the point p, i.e. check if p either...
  !>    - lies inside grid cell aij_in, ...
  !>    - lies on the b-grid point bij_on, or...
  !>    - lies on the edge cij_on
  subroutine trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)

    ! In/output variables
    type(type_grid),        intent(in)    :: grid
    real(dp), dimension(2), intent(in)    :: p
    integer,  dimension(2), intent(out)   :: aij_in, bij_on, cxij_on, cyij_on

    ! Local variables:
    integer  :: i,j
    real(dp) :: xl,xu,yl,yu

    ! Initialise
    aij_in  = [0,0]
    bij_on  = [0,0]
    cxij_on = [0,0]
    cyij_on = [0,0]

    ! Find the grid cell containing p
    i = 1 + floor( (p(1) - grid%xmin + grid%dx / 2._dp) / grid%dx)
    j = 1 + floor( (p(2) - grid%ymin + grid%dx / 2._dp) / grid%dx)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! Check if p lies on either of the four surrounding b-grid points
    if     (i > 1       .and. j > 1       .and. abs( p(1) - xl) < grid%tol_dist .and. abs( p(2) - yl) < grid%tol_dist) then
      ! p coincides with the southwest corner
      bij_on = [i-1,j-1]
      return
    elseif (i > 1       .and. j < grid%ny .and. abs( p(1) - xl) < grid%tol_dist .and. abs( p(2) - yu) < grid%tol_dist) then
      ! p coincides with the northwest corner
      bij_on = [i-1,j  ]
      return
    elseif (i < grid%nx .and. j < 1       .and. abs( p(1) - xu) < grid%tol_dist .and. abs( p(2) - yl) < grid%tol_dist) then
      ! p coincides with the southeast corner
      bij_on = [i  ,j-1]
      return
    elseif (i < grid%nx .and. j < grid%ny .and. abs( p(1) - xu) < grid%tol_dist .and. abs( p(2) - yu) < grid%tol_dist) then
      ! p coincides with the northeast corner
      bij_on = [i  ,j  ]
      return
    end if

    ! Check if p lies on any of the four borders
    if     (i > 1       .and. abs( p(1) - xl) < grid%tol_dist) then
      ! p coincides with the western border
      cxij_on = [i-1,j  ]
      return
    elseif (i < grid%nx .and. abs( p(1) - xu) < grid%tol_dist) then
      ! p coincides with the eastern border
      cxij_on = [i  ,j  ]
      return
    elseif (j > 1       .and. abs( p(2) - yl) < grid%tol_dist) then
      ! p coincides with the southern border
      cyij_on = [i  ,j-1]
      return
    elseif (j < grid%ny .and. abs( p(2) - yu) < grid%tol_dist) then
      ! p coincides with the northern border
      cyij_on = [i  ,j  ]
      return
    end if

    ! p doesn't lie on the corners or borders, so it must lie inside the grid cell
    aij_in = [i,j]

  end subroutine trace_line_grid_start

  !> Given the line [pq], where p lies inside grid cell aij_in,
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_a( grid, p, q, &
    aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),        intent(in)    :: grid
    real(dp), dimension(2), intent(in)    :: p,q
    integer,  dimension(2), intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(out)   :: n_left
    logical,                intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: xl,xu,yl,yu
    real(dp), dimension(2) :: sw,nw,se,ne
    logical                :: do_cross
    real(dp), dimension(2) :: llis

    ! Safety
    if ((aij_in( 1) == 0 .and. aij_in( 2) == 0) .or. cxij_on( 1) > 0 .or. cxij_on( 2) > 0 .or. &
        bij_on( 1) > 0 .or. bij_on( 2) > 0 .or. cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      call crash('trace_line_grid_a - coincidence indicators dont make sense!')
    end if

    i = aij_in( 1)
    j = aij_in( 2)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! More safety
    if (p(1) < xl .or. p(1) > xu .or. p(2) < yl .or. p(2) > yu) then
      call crash('trace_line_grid_a - coincidence indicators dont make sense!')
    end if

    ! Check if q lies inside the same grid cell
    if (q(1) >= xl - grid%tol_dist .and. &
        q(1) <= xu + grid%tol_dist .and. &
        q(2) >= yl - grid%tol_dist .and. &
        q(2) <= yu + grid%tol_dist) then
      ! q lies inside the same grid cell
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if pq passes through any of the four corners
    sw = [xl,yl]
    nw = [xl,yu]
    se = [xu,yl]
    ne = [xu,yu]

    if (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits this grid cell through the southwest corner
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    if (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits this grid cell through the northwest corner
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    if (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits this grid cell through the southeast corner
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    if (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits this grid cell through the northeast corner
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through any of the four boundaries
    call segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    call segment_intersection( p, q, se, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    call segment_intersection( p, q, sw, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    call segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j  ]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_a - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_a

  !> Given the line [pq], where p lies on b-grid point bij_on
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_b( grid, p, q, &
    aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),        intent(in)    :: grid
    real(dp), dimension(2), intent(in)    :: p,q
    integer,  dimension(2), intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(out)   :: n_left
    logical,                intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: x,y,xl,xu,yl,yu
    real(dp), dimension(2) :: sw,nw,se,ne,ww,ee,ss,nn
    logical                :: do_cross
    real(dp), dimension(2) :: llis

    ! Safety
    if (aij_in( 1) > 0 .or. aij_in( 2) > 0 .or. cxij_on( 1) > 0 .or. cxij_on( 2) > 0 .or. &
    (bij_on( 1) == 0 .and. bij_on( 2) == 0) .or. cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
      call crash('trace_line_grid_b - coincidence indicators dont make sense!')
    end if

    i = bij_on( 1)
    j = bij_on( 2)

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

    ! More safety
    if (abs( p(1) - x) > grid%tol_dist .or. abs( p(2) - y) > grid%tol_dist) then
      call crash('trace_line_grid_b - coincidence indicators dont make sense!')
    end if

    ! Check if q lies on the cy-grid edge to the west
    if (q(1) < x + grid%tol_dist .and. q(1) > xl - grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the cy-grid edge to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the cy-grid edge to the east
    if (q(1) > x - grid%tol_dist .and. q(1) < xu + grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the cy-grid edge to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the cx-grid edge to the south
    if (q(2) < y + grid%tol_dist .and. q(2) > yl - grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the cx-grid edge to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the cx-grid edge to the north
    if (q(2) > y - grid%tol_dist .and. q(2) < yu + grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the cx-grid edge to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the northwest
    if (q(1) > xl - grid%tol_dist .and. q(1) < x  + grid%tol_dist .and. &
        q(2) > y  - grid%tol_dist .and. q(2) < yu + grid%tol_dist) then
      ! q lies inside the a-grid cell to the northwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the northeast
    if (q(1) > x  - grid%tol_dist .and. q(1) < xu + grid%tol_dist .and. &
        q(2) > y  - grid%tol_dist .and. q(2) < yu + grid%tol_dist) then
      ! q lies inside the a-grid cell to the northeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the southeast
    if (q(1) > x  - grid%tol_dist .and. q(1) < xu + grid%tol_dist .and. &
        q(2) > yl - grid%tol_dist .and. q(2) < y  + grid%tol_dist) then
      ! q lies inside the a-grid cell to the southeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the southwest
    if (q(1) > xl - grid%tol_dist .and. q(1) < x  + grid%tol_dist .and. &
        q(2) > yl - grid%tol_dist .and. q(2) < y  + grid%tol_dist) then
      ! q lies inside the a-grid cell to the southwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the west
    if (lies_on_line_segment( p, q, ww, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the north
    if (lies_on_line_segment( p, q, nn, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the west
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i  ,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the east
    if (lies_on_line_segment( p, q, ee, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i+1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the south
    if (lies_on_line_segment( p, q, ss, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northwest through its western boundary
    call segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northwest through its northern boundary
    call segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northwest through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northeast through its northern boundary
    call segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northeast through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j+1]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northeast through its eastern boundary
    call segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southeast through its eastern boundary
    call segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southeast through its southern boundary
    call segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southeast through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southwest through its southern boundary
    call segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southwest through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southwest through its western boundary
    call segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_b - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_b

  !> Given the line [pq], where p lies on cx-grid edge cxij_on
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_cx( grid, p, q, &
    aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),        intent(in)    :: grid
    real(dp), dimension(2), intent(in)    :: p,q
    integer,  dimension(2), intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(out)   :: n_left
    logical,                intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: x,yl,yu
    real(dp), dimension(2) :: sw,nw,se,ne,ss,nn
    logical                :: do_cross
    real(dp), dimension(2) :: llis

    ! Safety
    if (aij_in( 1) > 0 .or. aij_in( 2) > 0 .or. (cxij_on( 1) == 0 .and. cxij_on( 2) == 0) .or. &
    bij_on( 1) > 0 .or. bij_on( 2) > 0 .or. cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
      call crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    end if

    i = cxij_on( 1)
    j = cxij_on( 2)

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

    ! More safety
    if (p(2) < yl .or. p(2) > yu .or. abs( p(1) - x) > grid%tol_dist) then
      call crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    end if

    ! Check if q lies on the same cx-grid cell in the southern direction
    if (q(2) < p(2) .and. q(2) >= yl - grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the same cx-grid cell in the southern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the same cx-grid cell in the northern direction
    if (q(2) > p(2) .and. q(2) <= yu + grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the same cx-grid cell in the northern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the west
    if (q(2) >= yl - grid%tol_dist .and. q(2) <= yu + grid%tol_dist .and. &
        q(1) >= x - grid%dx - grid%tol_dist .and. q(1) <= x + grid%tol_dist) then
      ! q lies inside the grid cell to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the east
    if (q(2) >= yl - grid%tol_dist .and. q(2) <= yu + grid%tol_dist .and. &
        q(1) <= x + grid%dx + grid%tol_dist .and. q(1) >= x - grid%tol_dist) then
      ! q lies inside the grid cell to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the south
    if (lies_on_line_segment( p, q, ss, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the north
    if (lies_on_line_segment( p, q, nn, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the north
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through the b-grid point to the northwest
    if (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through the b-grid point to the southwest
    if (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through the b-grid point to the northeast
    if (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i+1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through the b-grid point to the southeast
    if (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i+1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through its southern boundary
    call segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the west through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through its western boundary
    call segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the west through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through its northern boundary
    call segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the west through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through its northern boundary
    call segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the east through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through its eastern boundary
    call segment_intersection( p, q, ne, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the east through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through its southern boundary
    call segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the east through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_cx - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_cx

  !> Given the line [pq], where p lies on cy-grid edge cyij_on
  !> find the point p_next where [pq] crosses into the next grid cell.
  subroutine trace_line_grid_cy( grid, p, q, &
    aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)

    ! In/output variables
    type(type_grid),        intent(in)    :: grid
    real(dp), dimension(2), intent(in)    :: p,q
    integer,  dimension(2), intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2), intent(out)   :: p_next
    integer,                intent(out)   :: n_left
    logical,                intent(out)   :: coincides, finished

    ! Local variables:
    integer                :: i,j
    real(dp)               :: xl,xu,y
    real(dp), dimension(2) :: sw,nw,se,ne,ww,ee
    logical                :: do_cross
    real(dp), dimension(2) :: llis

    ! Safety
    if (aij_in( 1) > 0 .or. aij_in( 2) > 0 .or. cxij_on( 1) > 0 .or. cxij_on( 2) > 0 .or. &
        bij_on( 1) > 0 .or. bij_on( 2) > 0 .or. (cyij_on( 1) == 0 .and. cyij_on( 2) == 0)) then
      call crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    end if

    i = cyij_on( 1)
    j = cyij_on( 2)

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

    ! More safety
    if (p(1) < xl .or. p(1) > xu .or. abs( p(2) - y) > grid%tol_dist) then
      call crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    end if

    ! Check if q lies on the same cy-grid cell in the western direction
    if (q(1) < p(1) .and. q(1) >= xl - grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the same cy-grid cell in the western direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the same cy-grid cell in the eastern direction
    if (q(1) > p(1) .and. q(1) <= xu + grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the same cy-grid cell in the eastern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the south
    if (q(1) >= xl - grid%tol_dist .and. q(1) <= xu + grid%tol_dist .and. &
        q(2) >= y - grid%dx - grid%tol_dist .and. q(2) <= y + grid%tol_dist) then
      ! q lies inside the grid cell to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the north
    if (q(1) >= xl - grid%tol_dist .and. q(1) <= xu + grid%tol_dist .and. &
        q(2) <= y + grid%dx + grid%tol_dist .and. q(2) >= y - grid%tol_dist) then
      ! q lies inside the grid cell to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the west
    if (lies_on_line_segment( p, q, ww, grid%tol_dist))  then
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the east
    if (lies_on_line_segment( p, q, ee, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northwest
    if (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northeast
    if (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southwest
    if (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southeast
    if (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through its western boundary
    call segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the north through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through its northern boundary
    call segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the north through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through its eastern boundary
    call segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the north through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through its eastern boundary
    call segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the south through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through its southern boundary
    call segment_intersection( p, q, se, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the south through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through its western boundary
    call segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the south through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_cy - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_cy

end module line_tracing_grid
