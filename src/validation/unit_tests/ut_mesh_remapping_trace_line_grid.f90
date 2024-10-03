module ut_mesh_remapping_trace_line_grid

  ! Unit tests for mesh functions - remapping - trace_line_grid

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use grid_types, only: type_grid
  use grid_basic, only: setup_square_grid
  use mpi_basic, only: par

  implicit none

  private

  public :: test_trace_line_grid

  type type_coinc_ind_grid
    integer :: grid
    integer :: i,j
  end type type_coinc_ind_grid

  integer, parameter :: no_value = -1
  integer, parameter :: a_grid   = 1
  integer, parameter :: b_grid   = 2
  integer, parameter :: cx_grid  = 3
  integer, parameter :: cy_grid  = 4

contains

  subroutine test_trace_line_grid( test_name_parent)
    ! Test the trace_line_grid subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid'
    character(len=1024), parameter :: test_name_local = 'trace_line_grid'
    character(len=1024)            :: test_name
    real(dp)                       :: xmin, xmax, ymin, ymax, dx
    character(len=1024)            :: name
    type(type_grid)                :: grid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Initialise square grid
    name = 'test_grid'
    xmin = 0._dp
    xmax = 10._dp
    ymin = 0._dp
    ymax = 10._dp
    dx   = 1._dp
    call setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid)

    call test_trace_line_grid_start( test_name, grid)
    call test_trace_line_grid_a    ( test_name, grid)
    call test_trace_line_grid_b    ( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid

  subroutine test_trace_line_grid_start( test_name_parent, grid)
    ! Test the trace_line_grid_start subroutine

    use line_tracing_grid

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_start'
    character(len=1024), parameter :: test_name_local = 'trace_line_grid_start'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid_start_p_inside_a( test_name, grid)
    call test_trace_line_grid_start_p_on_b    ( test_name, grid)
    call test_trace_line_grid_start_p_on_cx   ( test_name, grid)
    call test_trace_line_grid_start_p_on_cy   ( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_start

  subroutine test_trace_line_grid_start_p_inside_a( test_name_parent, grid)
    ! Test if trace_line_grid_start is able to identify points p that lie inside a-grid cells

    use line_tracing_grid, only: trace_line_grid_start

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid), intent(in ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_start_p_inside_a'
    character(len=1024), parameter :: test_name_local = 'p_inside_a'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,ii,jj
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, xmintol, xmaxtol, ymintol, ymaxtol
    real(dp), dimension(2)         :: p
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    logical                        :: verified_p_inside_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_inside_a = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_start_p_inside_a'//&
      ' only loops over grid interior; trace_line_grid_start'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

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
      n_sub = 50
      do ii = 1, n_sub
      do jj = 1, n_sub

        p( 1) = xmintol + (xmaxtol - xmintol) * real( ii-1,dp) / real( n_sub-1,dp)
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jj-1,dp) / real( n_sub-1,dp)

        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)

        verified_p_inside_a = verified_p_inside_a .and. &
          coinc_ind%grid == a_grid .and. &
          coinc_ind%i    == i .and. &
          coinc_ind%j    == j

      end do
      end do

    end do
    end do

    call unit_test( verified_p_inside_a, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_start_p_inside_a

  subroutine test_trace_line_grid_start_p_on_b( test_name_parent, grid)
    ! Test if trace_line_grid_start is able to identify points p that lie on b-grid points

    use line_tracing_grid, only: trace_line_grid_start

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid), intent(in ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_start_p_on_b'
    character(len=1024), parameter :: test_name_local = 'p_on_b'
    character(len=1024)            :: test_name
    integer                        :: i,j
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax
    real(dp), dimension(2)         :: p
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    logical                        :: verified_p_on_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_on_b = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_start_p_on_b'//&
      ' only loops over grid interior; trace_line_grid_start'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

      x = grid%x( i)
      y = grid%y( j)

      xmin = x - grid%dx / 2._dp
      xmax = x + grid%dx / 2._dp
      ymin = y - grid%dx / 2._dp
      ymax = y + grid%dx / 2._dp

      ! Southwest (i.e. b-grid point [i-1,j-1])
      p = [xmin, ymin]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i-1 .and. &
        coinc_ind%j    == j-1

      ! Northwest (i.e. b-grid point [i-1,j])
      p = [xmin, ymax]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i-1 .and. &
        coinc_ind%j    == j

      ! Southeast (i.e. b-grid point [i,j-1])
      p = [xmax, ymin]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i .and. &
        coinc_ind%j    == j-1

      ! Northeast (i.e. b-grid point [i,j])
      p = [xmax, ymax]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i .and. &
        coinc_ind%j    == j

    end do
    end do

    call unit_test( verified_p_on_b, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_start_p_on_b

  subroutine test_trace_line_grid_start_p_on_cx( test_name_parent, grid)
    ! Test if trace_line_grid_start is able to identify points p that lie on cx-grid points

    use line_tracing_grid, only: trace_line_grid_start

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid), intent(in ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_start_p_on_cx'
    character(len=1024), parameter :: test_name_local = 'p_on_cx'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,jj
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, ymintol, ymaxtol
    real(dp), dimension(2)         :: p
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    logical                        :: verified_p_on_cx

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_on_cx = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_start_p_on_cx'//&
      ' only loops over grid interior; trace_line_grid_start'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

      x = grid%x( i)
      y = grid%y( j)

      xmin = x - grid%dx / 2._dp
      xmax = x + grid%dx / 2._dp
      ymin = y - grid%dx / 2._dp
      ymax = y + grid%dx / 2._dp

      ymintol = ymin + grid%tol_dist * 2._dp
      ymaxtol = ymax - grid%tol_dist * 2._dp

      ! West (i.e. cx-grid point [i-1,j])
      p( 1) = xmin
      do jj = 1, n_sub
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jj-1,dp) / real( n_sub-1,dp)
        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
        verified_p_on_cx = verified_p_on_cx .and. &
          coinc_ind%grid == cx_grid .and. &
          coinc_ind%i    == i-1 .and. &
          coinc_ind%j    == j
      end do

      ! East (i.e. cx-grid point [i,j])
      p( 1) = xmax
      do jj = 1, n_sub
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jj-1,dp) / real( n_sub-1,dp)
        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
        verified_p_on_cx = verified_p_on_cx .and. &
          coinc_ind%grid == cx_grid .and. &
          coinc_ind%i    == i .and. &
          coinc_ind%j    == j
      end do

    end do
    end do

    call unit_test( verified_p_on_cx, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_start_p_on_cx

  subroutine test_trace_line_grid_start_p_on_cy( test_name_parent, grid)
    ! Test if trace_line_grid_start is able to identify points p that lie on cy-grid points

    use line_tracing_grid, only: trace_line_grid_start

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid), intent(in ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_start_p_on_cy'
    character(len=1024), parameter :: test_name_local = 'p_on_cy'
    character(len=1024)            :: test_name
    integer                        :: i,j,n_sub,ii
    real(dp)                       :: x, y, xmin, xmax, ymin, ymax, xmintol, xmaxtol
    real(dp), dimension(2)         :: p
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    logical                        :: verified_p_on_cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_on_cy = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_start_p_on_cy'//&
      ' only loops over grid interior; trace_line_grid_start'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

      x = grid%x( i)
      y = grid%y( j)

      xmin = x - grid%dx / 2._dp
      xmax = x + grid%dx / 2._dp
      ymin = y - grid%dx / 2._dp
      ymax = y + grid%dx / 2._dp

      xmintol = xmin + grid%tol_dist * 2._dp
      xmaxtol = xmax - grid%tol_dist * 2._dp

      ! South (i.e. cy-grid point [i,j-1])
      p( 2) = ymin
      do ii = 1, n_sub
        p( 1) = xmintol + (xmaxtol - xmintol) * real( ii-1,dp) / real( n_sub-1,dp)
        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
        verified_p_on_cy = verified_p_on_cy .and. &
          coinc_ind%grid == cy_grid .and. &
          coinc_ind%i    == i .and. &
          coinc_ind%j    == j-1
      end do

      ! North (i.e. cy-grid point [i,j])
      p( 2) = ymax
      do ii = 1, n_sub
        p( 1) = xmintol + (xmaxtol - xmintol) * real( ii-1,dp) / real( n_sub-1,dp)
        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, .false.)
        verified_p_on_cy = verified_p_on_cy .and. &
          coinc_ind%grid == cy_grid .and. &
          coinc_ind%i    == i .and. &
          coinc_ind%j    == j
      end do

    end do
    end do

    call unit_test( verified_p_on_cy, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_start_p_on_cy

  subroutine test_trace_line_grid_a( test_name_parent, grid)
    ! Test the trace_line_grid_a subroutine

    use line_tracing_grid

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_a'
    character(len=1024), parameter :: test_name_local = 'trace_line_grid_a'
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

    use line_tracing_grid

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
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_inside_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_inside_a = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_a_q_inside_a'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

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

          aij_in  = [i,j]
          bij_on  = [0,0]
          cxij_on = [0,0]
          cyij_on = [0,0]

          call trace_line_grid_a( grid, p, q, &
            aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
          coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

    use line_tracing_grid

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
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_exits_through_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_exits_through_b = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_a_pq_exits_through_b'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

    use line_tracing_grid

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
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_exits_through_cx

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_exits_through_cx = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_a_pq_exits_through_cx'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

    use line_tracing_grid

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
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_exits_through_cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_exits_through_cy = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_a_pq_exits_through_cy'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 2, grid%nx-1
    do j = 2, grid%ny-1

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [i,j]
        bij_on  = [0,0]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_a( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

  subroutine test_trace_line_grid_b( test_name_parent, grid)
    ! Test the trace_line_grid_b subroutine

    use line_tracing_grid

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_b'
    character(len=1024), parameter :: test_name_local = 'trace_line_grid_b'
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

    use line_tracing_grid

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
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_on_cx_cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_cx_cy = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_b_q_on_cx_cy'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      x = grid%x( i) + grid%dx / 2._dp
      y = grid%y( j) + grid%dx / 2._dp

      p = [x,y]

    ! q lies on the edge to the west (i.e. cy-grid point [i,j])
    ! =========================================================

      q = [x - grid%dx / 2._dp, y]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_on_cx_cy = verified_q_on_cx_cy .and. &
        test_tol( p_next, q, grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .true.)

    ! q lies on the edge to the east (i.e. cy-grid point [i+1,j])
    ! ===========================================================

      q = [x + grid%dx / 2._dp, y]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_on_cx_cy = verified_q_on_cx_cy .and. &
        test_tol( p_next, q, grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j+1) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .true.)

    ! q lies on the edge to the south (i.e. cx-grid point [i,j])
    ! ==========================================================

      q = [x, y - grid%dx / 2._dp]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_on_cx_cy = verified_q_on_cx_cy .and. &
        test_tol( p_next, q, grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        (coincides .eqv. .true.) .and. &
        (finished .eqv. .true.)

    ! q lies on the edge to the north (i.e. cx-grid point [i,j+1])
    ! ============================================================

      q = [x, y + grid%dx / 2._dp]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

    use line_tracing_grid

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
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_in_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_in_a = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_b_q_in_a'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
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

        aij_in  = [0,0]
        bij_on  = [i,j]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_b( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [0,0]
        bij_on  = [i,j]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_b( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [0,0]
        bij_on  = [i,j]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_b( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

        aij_in  = [0,0]
        bij_on  = [i,j]
        cxij_on = [0,0]
        cyij_on = [0,0]

        call trace_line_grid_b( grid, p, q, &
          aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

    use line_tracing_grid

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_b_q_through_b'
    character(len=1024), parameter :: test_name_local = 'q_through_b'
    character(len=1024)            :: test_name
    integer                        :: i,j
    real(dp)                       :: x, y
    real(dp), dimension(2)         :: p, q, p_next
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_through_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_through_b = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_b_q_through_b'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      x = grid%x( i) + grid%dx / 2._dp
      y = grid%y( j) + grid%dx / 2._dp

      p = [x,y]

    ! pq passes through the b-grid point to the west (i.e. b-grid point [i-1,j])
    ! ==========================================================================

      q = [x - grid%dx * 2._dp, y]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ! coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x - grid%dx, y - grid%dx], grid%tol_dist) .and. &
        bij_on( 1) == i-1 .and. &
        bij_on( 2) == j-1 .and. &
        n_left == grid%ij2n( i,j) .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the northwest (i.e. b-grid point [i-1,j+1])
    ! =================================================================================

      q = [x - grid%dx * 2._dp, y + grid%dx * 2._dp]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ! coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x - grid%dx, y + grid%dx], grid%tol_dist) .and. &
        bij_on( 1) == i-1 .and. &
        bij_on( 2) == j+1 .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the southeast (i.e. b-grid point [i+1,j-1])
    ! =================================================================================

      q = [x + grid%dx * 2._dp, y - grid%dx * 2._dp]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ! coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x + grid%dx, y - grid%dx], grid%tol_dist) .and. &
        bij_on( 1) == i+1 .and. &
        bij_on( 2) == j-1 .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! pq passes through the b-grid point to the southeast (i.e. b-grid point [i+1,j-1])
    ! =================================================================================

      q = [x + grid%dx * 2._dp, y + grid%dx * 2._dp]

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ! coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_b = verified_q_through_b .and. &
        test_tol( p_next, [x + grid%dx, y + grid%dx], grid%tol_dist) .and. &
        bij_on( 1) == i+1 .and. &
        bij_on( 2) == j+1 .and. &
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

    use line_tracing_grid

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
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    integer                        :: n_left
    logical                        :: coincides, finished
    logical                        :: verified_q_through_a

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_through_a = .true.

    ! DENK DROM
    if (par%master) call warning('DENK DROM - test_trace_line_grid_b_pq_through_a'//&
      ' only loops over grid interior; trace_line_grid_a'//&
      ' cannot yet handle the border cells properly; fix this!')
    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      x = grid%x( i) + grid%dx / 2._dp
      y = grid%y( j) + grid%dx / 2._dp

      p = [x,y]

    ! q exits a-grid cell [i,j] (to the southwest) through its western border (i.e. cx-grid point [i-1,j])
    ! ====================================================================================================

      q( 1) = x - grid%dx * 2._dp
      q( 2) = y - grid%dx

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx, y - grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i-1 .and. &
        coinc_ind%j == j .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! q exits a-grid cell [i,j] (to the southwest) through its southern border (i.e. cy-grid point [i,j-1])
    ! =====================================================================================================

      q( 1) = x - grid%dx
      q( 2) = y - grid%dx * 2._dp

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx / 2._dp, y - grid%dx], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j) .and. &
        coinc_ind%grid == cy_grid .and. &
        coinc_ind%i == i .and. &
        coinc_ind%j == j-1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! q exits a-grid cell [i,j+1] (to the northwest) through its western border (i.e. cx-grid point [i-1,j+1])
    ! ========================================================================================================

      q( 1) = x - grid%dx * 2._dp
      q( 2) = y + grid%dx

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx, y + grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i-1 .and. &
        coinc_ind%j == j+1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! q exits a-grid cell [i,j+1] (to the northwest) through its northern border (i.e. cy-grid point [i,j+1])
    ! =======================================================================================================

      q( 1) = x - grid%dx
      q( 2) = y + grid%dx * 2._dp

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x - grid%dx / 2._dp, y + grid%dx], grid%tol_dist) .and. &
        n_left == grid%ij2n( i,j+1) .and. &
        coinc_ind%grid == cy_grid .and. &
        coinc_ind%i == i .and. &
        coinc_ind%j == j+1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! q exits a-grid cell [i+1,j] (to the southeast) through its eastern border (i.e. cx-grid point [i+1,j])
    ! ======================================================================================================

      q( 1) = x + grid%dx * 2._dp
      q( 2) = y - grid%dx

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x + grid%dx, y - grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! q exits a-grid cell [i+1,j] (to the southeast) through its southern border (i.e. cy-grid point [i+1,j-1])
    ! =========================================================================================================

      q( 1) = x + grid%dx
      q( 2) = y - grid%dx * 2._dp

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x + grid%dx / 2._dp, y - grid%dx], grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j) .and. &
        coinc_ind%grid == cy_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j-1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! q exits a-grid cell [i+1,j+1] (to the northeast) through its eastern border (i.e. cx-grid point [i+1,j+1])
    ! ==========================================================================================================

      q( 1) = x + grid%dx * 2._dp
      q( 2) = y + grid%dx

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

      verified_q_through_a = verified_q_through_a .and. &
        test_tol( p_next, [x + grid%dx, y + grid%dx / 2._dp], grid%tol_dist) .and. &
        n_left == grid%ij2n( i+1,j+1) .and. &
        coinc_ind%grid == cx_grid .and. &
        coinc_ind%i == i+1 .and. &
        coinc_ind%j == j+1 .and. &
        (coincides .eqv. .false.) .and. &
        (finished .eqv. .false.)

    ! q exits a-grid cell [i+1,j+1] (to the northeast) through its northern border (i.e. cy-grid point [i+1,j+1])
    ! ===========================================================================================================

      q( 1) = x + grid%dx
      q( 2) = y + grid%dx * 2._dp

      aij_in  = [0,0]
      bij_on  = [i,j]
      cxij_on = [0,0]
      cyij_on = [0,0]

      call trace_line_grid_b( grid, p, q, &
        aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished)

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

  function old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on, finished) result( coinc_ind)
    integer, dimension(2), intent(in) :: aij_in, bij_on, cxij_on, cyij_on
    logical,               intent(in) :: finished
    type(type_coinc_ind_grid) :: coinc_ind

    if (aij_in( 1) > 0 .or. aij_in( 2) > 0) then
      coinc_ind%grid = a_grid
      coinc_ind%i = aij_in( 1)
      coinc_ind%j = aij_in( 2)
    elseif (bij_on(  1) > 0 .or. bij_on(  2) > 0) then
      coinc_ind%grid = b_grid
      coinc_ind%i = bij_on( 1)
      coinc_ind%j = bij_on( 2)
    elseif (cxij_on( 1) > 0 .or. cxij_on( 2) > 0) then
      coinc_ind%grid = cx_grid
      coinc_ind%i = cxij_on( 1)
      coinc_ind%j = cxij_on( 2)
    elseif (cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
      coinc_ind%grid = cy_grid
      coinc_ind%i = cyij_on( 1)
      coinc_ind%j = cyij_on( 2)
    else
      if (finished) then
        coinc_ind%grid = no_value
        coinc_ind%i    = 0
        coinc_ind%j    = 0
      else
        call crash('old2new_coinc_ind found no coincidence indicators!')
      end if
    end if

  end function old2new_coinc_ind

end module ut_mesh_remapping_trace_line_grid