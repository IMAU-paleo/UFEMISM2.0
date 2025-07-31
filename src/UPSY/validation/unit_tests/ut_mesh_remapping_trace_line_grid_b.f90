module ut_mesh_remapping_trace_line_grid_b

  ! Unit tests for mesh functions - remapping - trace_line_grid_b

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use grid_types, only: type_grid
  use line_tracing_basic
  use line_tracing_grid
  use mpi_basic, only: par

  implicit none

  private

  public :: test_trace_line_grid_b

contains

  subroutine test_trace_line_grid_b( test_name_parent, grid)
    ! Test the trace_line_grid_b subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_b'
    character(len=1024), parameter :: test_name_local = 'b'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid_b_q_on_cx_cy  ( test_name, grid)
    call test_trace_line_grid_b_q_in_a      ( test_name, grid)
    call test_trace_line_grid_b_pq_through_b( test_name, grid)
    call test_trace_line_grid_b_pq_through_a( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_b

  subroutine test_trace_line_grid_b_q_on_cx_cy( test_name_parent, grid)
    ! Test the trace_line_grid_b subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_b_q_on_cx_cy'
    character(len=1024), parameter :: test_name_local = 'q_on_cx_cy'
    character(len=1024)            :: test_name
    integer                        :: i,j
    real(dp)                       :: x, y
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_on_cx_cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_cx_cy = .true.

    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      x = grid%x( i) + grid%dx / 2._dp
      y = grid%y( j) + grid%dx / 2._dp

      p = [x,y]

    ! q lies on the edge to the west (i.e. cy-grid point [i,j])
    ! =========================================================

      q = [x - grid%dx / 2._dp, y]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_on_cx_cy = verified_q_on_cx_cy .and. &
        test_tol( p_next, q, grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .true.)

    ! q lies on the edge to the east (i.e. cy-grid point [i+1,j])
    ! ===========================================================

      q = [x + grid%dx / 2._dp, y]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_on_cx_cy = verified_q_on_cx_cy .and. &
        test_tol( p_next, q, grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j+1) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .true.)

    ! q lies on the edge to the south (i.e. cx-grid point [i,j])
    ! ==========================================================

      q = [x, y - grid%dx / 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_on_cx_cy = verified_q_on_cx_cy .and. &
        test_tol( p_next, q, grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .true.)

    ! q lies on the edge to the north (i.e. cx-grid point [i,j+1])
    ! ============================================================

      q = [x, y + grid%dx / 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_on_cx_cy = verified_q_on_cx_cy .and. &
        test_tol( p_next, q, grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .true.)

    end do
    end do

    call unit_test( verified_q_on_cx_cy, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_b_q_on_cx_cy

  subroutine test_trace_line_grid_b_q_in_a( test_name_parent, grid)
    ! Test the trace_line_grid_b subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_b_q_in_a'
    character(len=1024), parameter :: test_name_local = 'q_in_a'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,iiq,jjq
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_in_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_in_a = .true.

    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      x = grid%x( i) + grid%dx / 2._dp
      y = grid%y( j) + grid%dx / 2._dp

      p = [x,y]

    ! q lies inside the a-grid cell to the southwest (i.e. a-grid cell [i,j])
    ! =======================================================================

      xmin = x - grid%dx
      xmax = xmin + grid%dx
      ymin = y - grid%dx
      ymax = ymin + grid%dx

      n_sub = 20
      do iiq = 1, n_sub
      do jjq = 1, n_sub

        q( 1) = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
        q( 2) = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)

        ! If q lies on the cx/cy edges, skip
        if (test_tol( q( 1), x, grid%tol_dist) .or. test_tol( q( 2), y, grid%tol_dist)) then
          cycle
        end if

        coinc_ind%grid = b_grid
        coinc_ind%i = i
        coinc_ind%j = j

        call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

        verified_q_in_a = verified_q_in_a .and. &
          test_tol( p_next, q, grid%tol_dist) .and. &
          n_left == grid%ij2n( i,j) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .true.)

      end do
      end do

    ! q lies inside the a-grid cell to the northwest (i.e. a-grid cell [i,j+1])
    ! =========================================================================

      xmin = x - grid%dx
      xmax = xmin + grid%dx
      ymin = y
      ymax = ymin + grid%dx

      n_sub = 20
      do iiq = 1, n_sub
      do jjq = 1, n_sub

        q( 1) = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
        q( 2) = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)

        ! If q lies on the cx/cy edges, skip
        if (test_tol( q( 1), x, grid%tol_dist) .or. test_tol( q( 2), y, grid%tol_dist)) then
          cycle
        end if

        coinc_ind%grid = b_grid
        coinc_ind%i = i
        coinc_ind%j = j

        call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

        verified_q_in_a = verified_q_in_a .and. &
          test_tol( p_next, q, grid%tol_dist) .and. &
          n_left == grid%ij2n( i,j+1) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .true.)

      end do
      end do

    ! q lies inside the a-grid cell to the southeast (i.e. a-grid cell [i+1,j])
    ! =========================================================================

      xmin = x
      xmax = xmin + grid%dx
      ymin = y - grid%dx
      ymax = ymin + grid%dx

      n_sub = 20
      do iiq = 1, n_sub
      do jjq = 1, n_sub

        q( 1) = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
        q( 2) = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)

        ! If q lies on the cx/cy edges, skip
        if (test_tol( q( 1), x, grid%tol_dist) .or. test_tol( q( 2), y, grid%tol_dist)) then
          cycle
        end if

        coinc_ind%grid = b_grid
        coinc_ind%i = i
        coinc_ind%j = j

        call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

        verified_q_in_a = verified_q_in_a .and. &
          test_tol( p_next, q, grid%tol_dist) .and. &
          n_left == grid%ij2n( i+1,j) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .true.)

      end do
      end do

    ! q lies inside the a-grid cell to the northeast (i.e. a-grid cell [i+1,j+1])
    ! ===========================================================================

      xmin = x
      xmax = xmin + grid%dx
      ymin = y
      ymax = ymin + grid%dx

      n_sub = 20
      do iiq = 1, n_sub
      do jjq = 1, n_sub

        q( 1) = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
        q( 2) = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)

        ! If q lies on the cx/cy edges, skip
        if (test_tol( q( 1), x, grid%tol_dist) .or. test_tol( q( 2), y, grid%tol_dist)) then
          cycle
        end if

        coinc_ind%grid = b_grid
        coinc_ind%i = i
        coinc_ind%j = j

        call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

        verified_q_in_a = verified_q_in_a .and. &
          test_tol( p_next, q, grid%tol_dist) .and. &
          n_left == grid%ij2n( i+1,j+1) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .true.)

      end do
      end do

    end do
    end do

    call unit_test( verified_q_in_a, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_b_q_in_a

  subroutine test_trace_line_grid_b_pq_through_b( test_name_parent, grid)
    ! Test the trace_line_grid_b subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_b_q_through_b'
    character(len=1024), parameter :: test_name_local = 'pq_through_b'
    character(len=1024)            :: test_name
    integer                        :: i,j
    real(dp)                       :: x, y
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_through_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_through_b = .true.

    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      x = grid%x( i) + grid%dx / 2._dp
      y = grid%y( j) + grid%dx / 2._dp

      p = [x,y]

    ! pq passes through the b-grid point to the west (i.e. b-grid point [i-1,j])
    ! ==========================================================================

      q = [x - grid%dx * 2._dp, y]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x - grid%dx, y], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i-1 .and. &
        coinc_ind%j == j .and. &
        n_left == grid%ij2n( i,j) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the east (i.e. b-grid point [i+1,j])
    ! ==========================================================================

      q = [x + grid%dx * 2._dp, y]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x + grid%dx, y], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j .and. &
        n_left == grid%ij2n( i+1,j+1) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the south (i.e. b-grid point [i,j-1])
    ! ===========================================================================

      q = [x, y - grid%dx * 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x, y - grid%dx], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i .and. &
        coinc_ind%j == j-1 .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the north (i.e. b-grid point [i,j+1])
    ! ===========================================================================

      q = [x, y + grid%dx * 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x, y + grid%dx], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i .and. &
        coinc_ind%j == j+1 .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the southwest (i.e. b-grid point [i-1,j-1])
    ! =================================================================================

      q = [x - grid%dx * 2._dp, y - grid%dx * 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x - grid%dx, y - grid%dx], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i-1 .and. &
        coinc_ind%j == j-1 .and. &
        n_left == grid%ij2n( i,j) .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the northwest (i.e. b-grid point [i-1,j+1])
    ! =================================================================================

      q = [x - grid%dx * 2._dp, y + grid%dx * 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x - grid%dx, y + grid%dx], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i-1 .and. &
        coinc_ind%j == j+1 .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the southeast (i.e. b-grid point [i+1,j-1])
    ! =================================================================================

      q = [x + grid%dx * 2._dp, y - grid%dx * 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x + grid%dx, y - grid%dx], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j-1 .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the southeast (i.e. b-grid point [i+1,j-1])
    ! =================================================================================

      q = [x + grid%dx * 2._dp, y + grid%dx * 2._dp]

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x + grid%dx, y + grid%dx], grid%tol_dist) .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j+1 .and. &
        n_left == grid%ij2n( i+1,j+1) .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    end do
    end do

    call unit_test( verified_q_through_b, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_b_pq_through_b

  subroutine test_trace_line_grid_b_pq_through_a( test_name_parent, grid)
    ! Test the trace_line_grid_b subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_b_pq_through_a'
    character(len=1024), parameter :: test_name_local = 'pq_through_a'
    character(len=1024)            :: test_name
    integer                        :: i,j
    real(dp)                       :: x, y
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_through_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_through_a = .true.

    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      x = grid%x( i) + grid%dx / 2._dp
      y = grid%y( j) + grid%dx / 2._dp

      p = [x,y]

    ! pq exits a-grid cell [i,j] (to the southwest) through its western border (i.e. cx-grid point [i-1,j])
    ! ====================================================================================================

      q( 1) = x - grid%dx * 2._dp
      q( 2) = y - grid%dx

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx, y - grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i-1 .and. &
        coinc_ind%j == j .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq exits a-grid cell [i,j] (to the southwest) through its southern border (i.e. cy-grid point [i,j-1])
    ! =====================================================================================================

      q( 1) = x - grid%dx
      q( 2) = y - grid%dx * 2._dp

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx / 2._dp, y - grid%dx], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j) .and. &
        coinc_ind%grid == cy_grid .and. &
        coinc_ind%i == i .and. &
        coinc_ind%j == j-1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq exits a-grid cell [i,j+1] (to the northwest) through its western border (i.e. cx-grid point [i-1,j+1])
    ! ========================================================================================================

      q( 1) = x - grid%dx * 2._dp
      q( 2) = y + grid%dx

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx, y + grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i-1 .and. &
        coinc_ind%j == j+1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq exits a-grid cell [i,j+1] (to the northwest) through its northern border (i.e. cy-grid point [i,j+1])
    ! =======================================================================================================

      q( 1) = x - grid%dx
      q( 2) = y + grid%dx * 2._dp

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx / 2._dp, y + grid%dx], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        coinc_ind%grid == cy_grid .and. &
        coinc_ind%i == i .and. &
        coinc_ind%j == j+1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq exits a-grid cell [i+1,j] (to the southeast) through its eastern border (i.e. cx-grid point [i+1,j])
    ! ======================================================================================================

      q( 1) = x + grid%dx * 2._dp
      q( 2) = y - grid%dx

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x + grid%dx, y - grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq exits a-grid cell [i+1,j] (to the southeast) through its southern border (i.e. cy-grid point [i+1,j-1])
    ! =========================================================================================================

      q( 1) = x + grid%dx
      q( 2) = y - grid%dx * 2._dp

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x + grid%dx / 2._dp, y - grid%dx], grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        coinc_ind%grid == cy_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j-1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq exits a-grid cell [i+1,j+1] (to the northeast) through its eastern border (i.e. cx-grid point [i+1,j+1])
    ! ==========================================================================================================

      q( 1) = x + grid%dx * 2._dp
      q( 2) = y + grid%dx

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x + grid%dx, y + grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j+1) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j+1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq exits a-grid cell [i+1,j+1] (to the northeast) through its northern border (i.e. cy-grid point [i+1,j+1])
    ! ===========================================================================================================

      q( 1) = x + grid%dx
      q( 2) = y + grid%dx * 2._dp

      coinc_ind%grid = b_grid
      coinc_ind%i = i
      coinc_ind%j = j

      call trace_line_grid_b( grid, p, q, coinc_ind, p_next, n_left, coincides, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x + grid%dx / 2._dp, y + grid%dx], grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j+1) .and. &
        coinc_ind%grid == cy_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j+1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    end do
    end do

    call unit_test( verified_q_through_a, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_b_pq_through_a

end module ut_mesh_remapping_trace_line_grid_b