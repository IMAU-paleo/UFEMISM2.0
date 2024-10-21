module ut_mesh_remapping_trace_line_grid_a

  ! Unit tests for mesh functions - remapping - trace_line_grid_a

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

  public :: test_trace_line_grid_a

contains

  subroutine test_trace_line_grid_a( test_name_parent, grid)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_a'
    character(len=1024), parameter :: test_name_local = 'a'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid_a_q_inside_a         ( test_name, grid)
    call test_trace_line_grid_a_pq_exits_through_b ( test_name, grid)
    call test_trace_line_grid_a_pq_exits_through_cx( test_name, grid)
    call test_trace_line_grid_a_pq_exits_through_cy( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_a

  subroutine test_trace_line_grid_a_q_inside_a( test_name_parent, grid)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_a_q_inside_a'
    character(len=1024), parameter :: test_name_local = 'q_inside_a'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,iip,jjp,iiq,jjq
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, xmintol, xmaxtol, ymintol, ymaxtol
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_inside_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_inside_a = .true.

    do i = 1, grid%nx
    do j = 1, grid%ny

      x = grid%x( i)
      y = grid%y( j)

      xmin = x - grid%dx / 2._dp
      xmax = x + grid%dx / 2._dp
      ymin = y - grid%dx / 2._dp
      ymax = y + grid%dx / 2._dp

      xmintol = xmin + grid%tol_dist * 2._dp
      xmaxtol = xmax - grid%tol_dist * 2._dp
      ymintol = ymin + grid%tol_dist * 2._dp
      ymaxtol = ymax - grid%tol_dist * 2._dp

      ! Loop over a set of points spread out within this a-grid cell
      n_sub = 20
      do iip = 1, n_sub
      do jjp = 1, n_sub

        p( 1) = xmintol + (xmaxtol - xmintol) * real( iip-1,dp) / real( n_sub-1,dp)
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)

        do iiq = 1, n_sub
        do jjq = 1, n_sub

          if (iiq == iip .and. jjq == jjp) cycle

          q( 1) = xmintol + (xmaxtol - xmintol) * real( iiq-1,dp) / real( n_sub-1,dp)
          q( 2) = ymintol + (ymaxtol - ymintol) * real( jjq-1,dp) / real( n_sub-1,dp)

          coinc_ind%grid = a_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_a( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_q_inside_a = verified_q_inside_a .and. &
            test_tol( p_next, q, grid%tol_dist) .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            coinc_ind%j == 0 .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do
        end do

      end do
      end do

    end do
    end do

    call unit_test( verified_q_inside_a, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_a_q_inside_a

  subroutine test_trace_line_grid_a_pq_exits_through_b( test_name_parent, grid)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_a_pq_exits_through_b'
    character(len=1024), parameter :: test_name_local = 'pq_exits_through_b'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,iip,jjp
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, xmintol, xmaxtol, ymintol, ymaxtol
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_exits_through_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_exits_through_b = .true.

    do i = 1, grid%nx
    do j = 1, grid%ny

      x = grid%x( i)
      y = grid%y( j)

      xmin = x - grid%dx / 2._dp
      xmax = x + grid%dx / 2._dp
      ymin = y - grid%dx / 2._dp
      ymax = y + grid%dx / 2._dp

      xmintol = xmin + grid%tol_dist * 2._dp
      xmaxtol = xmax - grid%tol_dist * 2._dp
      ymintol = ymin + grid%tol_dist * 2._dp
      ymaxtol = ymax - grid%tol_dist * 2._dp

      ! Loop over a set of points spread out within this a-grid cell
      n_sub = 20
      do iip = 1, n_sub
      do jjp = 1, n_sub

        p( 1) = xmintol + (xmaxtol - xmintol) * real( iip-1,dp) / real( n_sub-1,dp)
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)

      ! pq exits through the southwest corner (i.e. b-grid point [i-1,j-1])
      ! ===================================================================

        q(1) = 2*xmin - p(1)
        q(2) = 2*ymin - p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_b = verified_pq_exits_through_b .and. &
          test_tol( p_next, [xmin,ymin], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i-1 .and. &
          coinc_ind%j == j-1 .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq exits through the northwest corner (i.e. b-grid point [i-1,j])
      ! ===================================================================

        q(1) = 2*xmin - p(1)
        q(2) = 2*ymax - p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_b = verified_pq_exits_through_b .and. &
          test_tol( p_next, [xmin,ymax], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i-1 .and. &
          coinc_ind%j == j .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq exits through the southeast corner (i.e. b-grid point [i,j-1])
      ! =================================================================

        q(1) = 2*xmax - p(1)
        q(2) = 2*ymin - p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_b = verified_pq_exits_through_b .and. &
          test_tol( p_next, [xmax,ymin], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i .and. &
          coinc_ind%j == j-1 .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq exits through the northeast corner (i.e. b-grid point [i,j])
      ! ===============================================================

        q(1) = 2*xmax - p(1)
        q(2) = 2*ymax - p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_b = verified_pq_exits_through_b .and. &
          test_tol( p_next, [xmax,ymax], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i .and. &
          coinc_ind%j == j .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      end do
      end do

    end do
    end do

    call unit_test( verified_pq_exits_through_b, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_a_pq_exits_through_b

  subroutine test_trace_line_grid_a_pq_exits_through_cx( test_name_parent, grid)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_a_pq_exits_through_cx'
    character(len=1024), parameter :: test_name_local = 'pq_exits_through_cx'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,iip,jjp
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, xmintol, xmaxtol, ymintol, ymaxtol
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_exits_through_cx

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_exits_through_cx = .true.

    do i = 1, grid%nx
    do j = 1, grid%ny

      x = grid%x( i)
      y = grid%y( j)

      xmin = x - grid%dx / 2._dp
      xmax = x + grid%dx / 2._dp
      ymin = y - grid%dx / 2._dp
      ymax = y + grid%dx / 2._dp

      xmintol = xmin + grid%tol_dist * 2._dp
      xmaxtol = xmax - grid%tol_dist * 2._dp
      ymintol = ymin + grid%tol_dist * 2._dp
      ymaxtol = ymax - grid%tol_dist * 2._dp

      ! Loop over a set of points spread out within this a-grid cell
      n_sub = 20
      do iip = 1, n_sub
      do jjp = 1, n_sub

        p( 1) = xmintol + (xmaxtol - xmintol) * real( iip-1,dp) / real( n_sub-1,dp)
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)

      ! pq exits through the western cell border (i.e. cx-grid point [i-1,j])
      ! =====================================================================

        q(1) = 2*xmin - p(1)
        q(2) = p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_cx = verified_pq_exits_through_cx .and. &
          test_tol( p_next, [xmin,p(2)], grid%tol_dist) .and. &
          coinc_ind%grid == cx_grid .and. &
          coinc_ind%i == i-1 .and. &
          coinc_ind%j == j .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq exits through the eastern cell border (i.e. cx-grid point [i,j])
      ! ===================================================================

        q(1) = 2*xmax - p(1)
        q(2) = p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_cx = verified_pq_exits_through_cx .and. &
          test_tol( p_next, [xmax,p(2)], grid%tol_dist) .and. &
          coinc_ind%grid == cx_grid .and. &
          coinc_ind%i == i .and. &
          coinc_ind%j == j .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      end do
      end do

    end do
    end do

    call unit_test( verified_pq_exits_through_cx, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_a_pq_exits_through_cx

  subroutine test_trace_line_grid_a_pq_exits_through_cy( test_name_parent, grid)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_a_pq_exits_through_cy'
    character(len=1024), parameter :: test_name_local = 'pq_exits_through_cy'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,iip,jjp
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, xmintol, xmaxtol, ymintol, ymaxtol
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_exits_through_cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_exits_through_cy = .true.

    do i = 1, grid%nx
    do j = 1, grid%ny

      x = grid%x( i)
      y = grid%y( j)

      xmin = x - grid%dx / 2._dp
      xmax = x + grid%dx / 2._dp
      ymin = y - grid%dx / 2._dp
      ymax = y + grid%dx / 2._dp

      xmintol = xmin + grid%tol_dist * 2._dp
      xmaxtol = xmax - grid%tol_dist * 2._dp
      ymintol = ymin + grid%tol_dist * 2._dp
      ymaxtol = ymax - grid%tol_dist * 2._dp

      ! Loop over a set of points spread out within this a-grid cell
      n_sub = 20
      do iip = 1, n_sub
      do jjp = 1, n_sub

        p( 1) = xmintol + (xmaxtol - xmintol) * real( iip-1,dp) / real( n_sub-1,dp)
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)

      ! pq exits through the southern cell border (i.e. cy-grid point [i,j-1])
      ! ======================================================================

        q(1) = p(1)
        q(2) = 2*ymin - p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_cy = verified_pq_exits_through_cy .and. &
          test_tol( p_next, [p(1),ymin], grid%tol_dist) .and. &
          coinc_ind%grid == cy_grid .and. &
          coinc_ind%i == i .and. &
          coinc_ind%j == j-1 .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq exits through the northern cell border (i.e. cy-grid point [i,j])
      ! ====================================================================

        q(1) = p(1)
        q(2) = 2*ymax - p(2)

        coinc_ind%grid = a_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_a( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_exits_through_cy = verified_pq_exits_through_cy .and. &
          test_tol( p_next, [p(1),ymax], grid%tol_dist) .and. &
          coinc_ind%grid == cy_grid .and. &
          coinc_ind%i == i .and. &
          coinc_ind%j == j .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      end do
      end do

    end do
    end do

    call unit_test( verified_pq_exits_through_cy, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_a_pq_exits_through_cy

end module ut_mesh_remapping_trace_line_grid_a