module ut_mesh_remapping

  ! Unit tests for mesh functions - remapping.

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

  public :: test_remapping

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

  subroutine test_remapping( test_name_parent)
    ! Test the remapping subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_remapping'
    character(len=1024), parameter :: test_name_local = 'remapping'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remapping

  subroutine test_trace_line_grid( test_name_parent)
    ! Test the trace_line_grid subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid'
    character(len=1024), parameter :: test_name_local = 'trace_line_grid'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid_start( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid

  subroutine test_trace_line_grid_start( test_name_parent)
    ! Test the trace_line_grid_start subroutine

    use line_tracing_grid

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_grid_start'
    character(len=1024), parameter :: test_name_local = 'trace_line_grid_start'
    character(len=1024)            :: test_name
    real(dp)                       :: xmin, xmax, ymin, ymax, dx
    character(len=1024)            :: name
    type(type_grid)                :: grid
    integer                        :: i,j,n_sub,ii,jj
    real(dp)                       :: x,y,xmintol,xmaxtol,ymintol,ymaxtol
    real(dp), dimension(2)         :: p
    integer, dimension(2)          :: aij_in, bij_on, cxij_on, cyij_on
    type(type_coinc_ind_grid)      :: coinc_ind
    logical                        :: verified_p_inside_a
    logical                        :: verified_p_on_b
    logical                        :: verified_p_on_cx
    logical                        :: verified_p_on_cy

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

    verified_p_inside_a  = .true.
    verified_p_on_b  = .true.
    verified_p_on_cx = .true.
    verified_p_on_cy = .true.

    n_sub = 50

    ! DENK DROM
    if (par%master) call warning('DENK DROM - only loops over grid interior; trace_line_grid_start'//&
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

    ! p lies inside an a-grid cell
    ! ============================

      ! Loop over a set of points spread out within this a-grid cell
      do ii = 1, n_sub
      do jj = 1, n_sub
        p( 1) = xmintol + (xmaxtol - xmintol) * real( ii-1,dp) / real( n_sub-1,dp)
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jj-1,dp) / real( n_sub-1,dp)
        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
        verified_p_inside_a = verified_p_inside_a .and. &
          coinc_ind%grid == a_grid .and. &
          coinc_ind%i    == i .and. &
          coinc_ind%j    == j
      end do
      end do

    ! p lies on a b-grid point
    ! ========================

      ! Southwest (i.e. b-grid point [i-1,j-1])
      p = [xmin, ymin]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i-1 .and. &
        coinc_ind%j    == j-1

      ! Northwest (i.e. b-grid point [i-1,j])
      p = [xmin, ymax]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i-1 .and. &
        coinc_ind%j    == j

      ! Southeast (i.e. b-grid point [i,j-1])
      p = [xmax, ymin]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i .and. &
        coinc_ind%j    == j-1

      ! Northeast (i.e. b-grid point [i,j])
      p = [xmax, ymax]
      call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
      coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
      verified_p_on_b = verified_p_on_b .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i    == i .and. &
        coinc_ind%j    == j

    ! p lies on a cx-grid point
    ! =========================

      ! West (i.e. cx-grid point [i-1,j])
      p( 1) = xmin
      do jj = 1, n_sub
        p( 2) = ymintol + (ymaxtol - ymintol) * real( jj-1,dp) / real( n_sub-1,dp)
        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
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
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
        verified_p_on_cx = verified_p_on_cx .and. &
          coinc_ind%grid == cx_grid .and. &
          coinc_ind%i    == i .and. &
          coinc_ind%j    == j
      end do

    ! p lies on a cy-grid point
    ! =========================

      ! South (i.e. cy-grid point [i,j-1])
      p( 2) = ymin
      do ii = 1, n_sub
        p( 1) = xmintol + (xmaxtol - xmintol) * real( ii-1,dp) / real( n_sub-1,dp)
        call trace_line_grid_start( grid, p, aij_in, bij_on, cxij_on, cyij_on)
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
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
        coinc_ind = old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on)
        verified_p_on_cy = verified_p_on_cy .and. &
          coinc_ind%grid == cy_grid .and. &
          coinc_ind%i    == i .and. &
          coinc_ind%j    == j
      end do

    end do
    end do

    call unit_test( verified_p_inside_a, trim(test_name)//'/p_inside_a')
    call unit_test( verified_p_on_b    , trim(test_name)//'/p_on_b')
    call unit_test( verified_p_on_cx   , trim(test_name)//'/p_on_cx')
    call unit_test( verified_p_on_cy   , trim(test_name)//'/p_on_cy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_start

  function old2new_coinc_ind( aij_in, bij_on, cxij_on, cyij_on) result( coinc_ind)
    integer, dimension(2), intent(in) :: aij_in, bij_on, cxij_on, cyij_on
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
      call crash('old2new_coinc_ind found no coincidence indicators!')
    end if

  end function old2new_coinc_ind

end module ut_mesh_remapping