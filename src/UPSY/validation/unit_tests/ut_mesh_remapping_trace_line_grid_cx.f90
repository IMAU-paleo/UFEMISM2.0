module ut_mesh_remapping_trace_line_grid_cx

  ! Unit tests for mesh functions - remapping - trace_line_grid_cx

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

  public :: test_trace_line_grid_cx

contains

  subroutine test_trace_line_grid_cx( test_name_parent, grid)
    ! Test the trace_line_grid_cx subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_cx'
    character(len=1024), parameter :: test_name_local = 'cx'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid_cx_q_on_cx     ( test_name, grid)
    call test_trace_line_grid_cx_q_in_a      ( test_name, grid)
    call test_trace_line_grid_cx_pq_through_b( test_name, grid)
    call test_trace_line_grid_cx_pq_through_a( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_cx

  subroutine test_trace_line_grid_cx_q_on_cx( test_name_parent, grid)
    ! Test the trace_line_grid_cx subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_cx_q_on_cx'
    character(len=1024), parameter :: test_name_local = 'q_on_cx'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,jjp,jjq
    real(dp)                       :: x, y, ymintol, ymaxtol, xp, yp, xq, yq
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_on_cx

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_cx = .true.

    do i = 2, grid%nx-2
    do j = 2, grid%ny-1

    ! p lies on the border between a-grid cells [i,j] and [i+1,j],
    ! which is represented here by cx-grid point [i,j]
    ! ================================================

      x = (grid%x( i) + grid%x( i+1)) / 2._dp
      y = grid%y( j)
      ymintol = y - grid%dx / 2._dp + grid%tol_dist * 2._dp
      ymaxtol = y + grid%dx / 2._dp - grid%tol_dist * 2._dp

      n_sub = 20
      do jjp = 1, n_sub

        xp = x
        yp = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

      ! q lies on the same border, to the south of p
      ! ============================================

        do jjq = 1, jjp-1

          xq = x
          yq = ymintol + (ymaxtol - ymintol) * real( jjq-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_q_on_cx = verified_q_on_cx .and. &
            test_tol( p_next, q, grid%tol_dist) .and. &
            n_left == grid%ij2n( i+1,j) .and. &
            (coincides .eqv. .true.) .and. &
            (finished .eqv. .true.)

        end do

        ! q lies on the same border, to the north of p
        ! ============================================

        do jjq = jjp+1, n_sub

          xq = x
          yq = ymintol + (ymaxtol - ymintol) * real( jjq-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_q_on_cx = verified_q_on_cx .and. &
            test_tol( p_next, q, grid%tol_dist) .and. &
            n_left == grid%ij2n( i,j) .and. &
            (coincides .eqv. .true.) .and. &
            (finished .eqv. .true.)

        end do

      end do

    end do
    end do

    call unit_test( verified_q_on_cx, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_cx_q_on_cx

  subroutine test_trace_line_grid_cx_q_in_a( test_name_parent, grid)
    ! Test the trace_line_grid_cx subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_cx_q_in_a'
    character(len=1024), parameter :: test_name_local = 'q_in_a'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,jjp,iiq,jjq
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, ymintol, ymaxtol, xp, yp, xq, yq
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

    do i = 2, grid%nx-2
    do j = 2, grid%ny-1

    ! p lies on the border between a-grid cells [i,j] and [i+1,j],
    ! which is represented here by cx-grid point [i,j]
    ! ================================================

      x = (grid%x( i) + grid%x( i+1)) / 2._dp
      y = grid%y( j)
      ymintol = y - grid%dx / 2._dp + grid%tol_dist * 2._dp
      ymaxtol = y + grid%dx / 2._dp - grid%tol_dist * 2._dp

      n_sub = 20
      do jjp = 1, n_sub

        xp = x
        yp = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

      ! q lies inside a-grid cell [i,j], to the west
      ! ============================================

        xmin = grid%x( i) - grid%dx / 2._dp
        xmax = grid%x( i) + grid%dx / 2._dp - grid%tol_dist * 2._dp
        ymin = grid%y( j) - grid%dx / 2._dp
        ymax = grid%y( j) + grid%dx / 2._dp

        do iiq = 1, n_sub
        do jjq = 1, n_sub

          xq = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
          yq = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_q_in_a = verified_q_in_a .and. &
            test_tol( p_next, q, grid%tol_dist) .and. &
            n_left == grid%ij2n( i,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do
        end do

      ! q lies inside a-grid cell [i+1,j], to the east
      ! ==============================================

        xmin = grid%x( i+1) - grid%dx / 2._dp + grid%tol_dist * 2._dp
        xmax = grid%x( i+1) + grid%dx / 2._dp
        ymin = grid%y( j) - grid%dx / 2._dp
        ymax = grid%y( j) + grid%dx / 2._dp

        do iiq = 1, n_sub
        do jjq = 1, n_sub

          xq = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
          yq = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_q_in_a = verified_q_in_a .and. &
            test_tol( p_next, q, grid%tol_dist) .and. &
            n_left == grid%ij2n( i+1,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do
        end do

      end do

    end do
    end do

    call unit_test( verified_q_in_a, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_cx_q_in_a

  subroutine test_trace_line_grid_cx_pq_through_b( test_name_parent, grid)
    ! Test the trace_line_grid_cx subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_cx_pq_through_b'
    character(len=1024), parameter :: test_name_local = 'pq_through_b'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,jjp
    real(dp)                       :: x, y, xx, yy, ymintol, ymaxtol, xp, yp, xq, yq
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_b = .true.

    do i = 2, grid%nx-2
    do j = 2, grid%ny-1

    ! p lies on the border between a-grid cells [i,j] and [i+1,j],
    ! which is represented here by cx-grid point [i,j]
    ! ================================================

      x = (grid%x( i) + grid%x( i+1)) / 2._dp
      y = grid%y( j)
      ymintol = y - grid%dx / 2._dp + grid%tol_dist * 2._dp
      ymaxtol = y + grid%dx / 2._dp - grid%tol_dist * 2._dp

      n_sub = 20
      do jjp = 1, n_sub

        xp = x
        yp = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

      ! pq passes through b-grid point [i,j], to the north
      ! ==================================================

        xq = x
        yq = grid%y( j+1)
        q = [xq,yq]

        coinc_ind%grid = cx_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_cx( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_through_b = verified_pq_through_b .and. &
          test_tol( p_next, [x, y + grid%dx / 2._dp], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i .and. &
          coinc_ind%j == j .and. &
          n_left == grid%ij2n( i,j) .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

      ! pq passes through b-grid point [i,j-1], to the south
      ! ====================================================

        xq = x
        yq = grid%y( j-1)
        q = [xq,yq]

        coinc_ind%grid = cx_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_cx( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_through_b = verified_pq_through_b .and. &
          test_tol( p_next, [x, y - grid%dx / 2._dp], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i .and. &
          coinc_ind%j == j-1 .and. &
          n_left == grid%ij2n( i+1,j) .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

      ! pq passes through b-grid point [i-1,j-1], to the southwest
      ! ==========================================================

        xx = grid%x( i-1) + grid%dx / 2._dp
        yy = grid%y( j-1) + grid%dx / 2._dp
        q( 1) = 2._dp * xx - p( 1)
        q( 2) = 2._dp * yy - p( 2)

        coinc_ind%grid = cx_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_cx( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_through_b = verified_pq_through_b .and. &
          test_tol( p_next, [xx, yy], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i-1 .and. &
          coinc_ind%j == j-1 .and. &
          n_left == grid%ij2n( i,j) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq passes through b-grid point [i-1,j], to the northwest
      ! ==========================================================

        xx = grid%x( i-1) + grid%dx / 2._dp
        yy = grid%y( j  ) + grid%dx / 2._dp
        q( 1) = 2._dp * xx - p( 1)
        q( 2) = 2._dp * yy - p( 2)

        coinc_ind%grid = cx_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_cx( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_through_b = verified_pq_through_b .and. &
          test_tol( p_next, [xx, yy], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i-1 .and. &
          coinc_ind%j == j .and. &
          n_left == grid%ij2n( i,j) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq passes through b-grid point [i+1,j-1], to the southeast
      ! ==========================================================

        xx = grid%x( i+1) + grid%dx / 2._dp
        yy = grid%y( j-1) + grid%dx / 2._dp
        q( 1) = 2._dp * xx - p( 1)
        q( 2) = 2._dp * yy - p( 2)

        coinc_ind%grid = cx_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_cx( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_through_b = verified_pq_through_b .and. &
          test_tol( p_next, [xx, yy], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i+1 .and. &
          coinc_ind%j == j-1 .and. &
          n_left == grid%ij2n( i+1,j) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      ! pq passes through b-grid point [i+1,j], to the northeast
      ! ==========================================================

        xx = grid%x( i+1) + grid%dx / 2._dp
        yy = grid%y( j  ) + grid%dx / 2._dp
        q( 1) = 2._dp * xx - p( 1)
        q( 2) = 2._dp * yy - p( 2)

        coinc_ind%grid = cx_grid
        coinc_ind%i    = i
        coinc_ind%j    = j

        call trace_line_grid_cx( grid, p, q, &
          coinc_ind, p_next, n_left, coincides, finished)

        verified_pq_through_b = verified_pq_through_b .and. &
          test_tol( p_next, [xx, yy], grid%tol_dist) .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == i+1 .and. &
          coinc_ind%j == j .and. &
          n_left == grid%ij2n( i+1,j) .and. &
          (coincides .eqv. .false.) .and. &
          (finished .eqv. .false.)

      end do

    end do
    end do

    call unit_test( verified_pq_through_b, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_cx_pq_through_b

  subroutine test_trace_line_grid_cx_pq_through_a( test_name_parent, grid)
    ! Test the trace_line_grid_cx subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_cx_pq_through_a'
    character(len=1024), parameter :: test_name_local = 'pq_through_a'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,jjp,iiq,jjq
    real(dp)                       :: x, y, xx, yy, xmintol, xmaxtol, ymintol, ymaxtol, xp, yp, xq, yq
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_a = .true.

    do i = 2, grid%nx-2
    do j = 2, grid%ny-1

    ! p lies on the border between a-grid cells [i,j] and [i+1,j],
    ! which is represented here by cx-grid point [i,j]
    ! ================================================

      x = (grid%x( i) + grid%x( i+1)) / 2._dp
      y = grid%y( j)
      ymintol = y - grid%dx / 2._dp + grid%tol_dist * 2._dp
      ymaxtol = y + grid%dx / 2._dp - grid%tol_dist * 2._dp

      n_sub = 20
      do jjp = 1, n_sub

        xp = x
        yp = ymintol + (ymaxtol - ymintol) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

      ! pq passes through the southern border of the a-grid cell to the west (i.e. cy-grid point [i,j-1])
      ! =================================================================================================

        yy = grid%y( j) - grid%dx / 2._dp - grid%tol_dist * 2._dp
        xmintol = grid%x( i) - grid%dx / 2._dp + grid%tol_dist * 3._dp
        xmaxtol = grid%x( i) + grid%dx / 2._dp - grid%tol_dist * 3._dp

        do iiq = 1, n_sub

          xq = xmintol + (xmaxtol - xmintol) * real( iiq-1,dp) / real( n_sub-1,dp)
          yq = yy
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_pq_through_a = verified_pq_through_a .and. &
            coinc_ind%grid == cy_grid .and. &
            coinc_ind%i == i .and. &
            coinc_ind%j == j-1 .and. &
            n_left == grid%ij2n( i,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      ! pq passes through the northern border of the a-grid cell to the west (i.e. cy-grid point [i,j])
      ! ===============================================================================================

        yy = grid%y( j) + grid%dx / 2._dp + grid%tol_dist * 2._dp
        xmintol = grid%x( i) - grid%dx / 2._dp + grid%tol_dist * 3._dp
        xmaxtol = grid%x( i) + grid%dx / 2._dp - grid%tol_dist * 3._dp

        do iiq = 1, n_sub

          xq = xmintol + (xmaxtol - xmintol) * real( iiq-1,dp) / real( n_sub-1,dp)
          yq = yy
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_pq_through_a = verified_pq_through_a .and. &
            coinc_ind%grid == cy_grid .and. &
            coinc_ind%i == i .and. &
            coinc_ind%j == j .and. &
            n_left == grid%ij2n( i,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      ! pq passes through the southern border of the a-grid cell to the east (i.e. cy-grid point [i+1,j-1])
      ! ===================================================================================================

        yy = grid%y( j) - grid%dx / 2._dp - grid%tol_dist * 2._dp
        xmintol = grid%x( i+1) - grid%dx / 2._dp + grid%tol_dist * 3._dp
        xmaxtol = grid%x( i+1) + grid%dx / 2._dp - grid%tol_dist * 3._dp

        do iiq = 1, n_sub

          xq = xmintol + (xmaxtol - xmintol) * real( iiq-1,dp) / real( n_sub-1,dp)
          yq = yy
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_pq_through_a = verified_pq_through_a .and. &
            coinc_ind%grid == cy_grid .and. &
            coinc_ind%i == i+1 .and. &
            coinc_ind%j == j-1 .and. &
            n_left == grid%ij2n( i+1,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      ! pq passes through the northern border of the a-grid cell to the east (i.e. cy-grid point [i+1,j])
      ! =================================================================================================

        yy = grid%y( j) + grid%dx / 2._dp + grid%tol_dist * 2._dp
        xmintol = grid%x( i+1) - grid%dx / 2._dp + grid%tol_dist * 3._dp
        xmaxtol = grid%x( i+1) + grid%dx / 2._dp - grid%tol_dist * 3._dp

        do iiq = 1, n_sub

          xq = xmintol + (xmaxtol - xmintol) * real( iiq-1,dp) / real( n_sub-1,dp)
          yq = yy
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_pq_through_a = verified_pq_through_a .and. &
            coinc_ind%grid == cy_grid .and. &
            coinc_ind%i == i+1 .and. &
            coinc_ind%j == j .and. &
            n_left == grid%ij2n( i+1,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      ! pq passes through the western border of the a-grid cell to the west (i.e. cx-grid point [i-1,j])
      ! =================================================================================================

        xx = grid%x( i) - grid%dx / 2._dp - grid%tol_dist * 2._dp
        ymintol = grid%y( j) - grid%dx / 2._dp + grid%tol_dist * 3._dp
        ymaxtol = grid%y( j) + grid%dx / 2._dp - grid%tol_dist * 3._dp

        do jjq = 1, n_sub

          xq = xx
          yq = ymintol + (ymaxtol - ymintol) * real( jjq-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_pq_through_a = verified_pq_through_a .and. &
            coinc_ind%grid == cx_grid .and. &
            coinc_ind%i == i-1 .and. &
            coinc_ind%j == j .and. &
            n_left == grid%ij2n( i,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      ! pq passes through the eastern border of the a-grid cell to the east (i.e. cx-grid point [i+1,j])
      ! =================================================================================================

        xx = grid%x( i+1) + grid%dx / 2._dp + grid%tol_dist * 2._dp
        ymintol = grid%y( j) - grid%dx / 2._dp + grid%tol_dist * 3._dp
        ymaxtol = grid%y( j) + grid%dx / 2._dp - grid%tol_dist * 3._dp

        do jjq = 1, n_sub

          xq = xx
          yq = ymintol + (ymaxtol - ymintol) * real( jjq-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          coinc_ind%grid = cx_grid
          coinc_ind%i    = i
          coinc_ind%j    = j

          call trace_line_grid_cx( grid, p, q, &
            coinc_ind, p_next, n_left, coincides, finished)

          verified_pq_through_a = verified_pq_through_a .and. &
            coinc_ind%grid == cx_grid .and. &
            coinc_ind%i == i+1 .and. &
            coinc_ind%j == j .and. &
            n_left == grid%ij2n( i+1,j) .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      end do

    end do
    end do

    call unit_test( verified_pq_through_a, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_cx_pq_through_a

end module ut_mesh_remapping_trace_line_grid_cx