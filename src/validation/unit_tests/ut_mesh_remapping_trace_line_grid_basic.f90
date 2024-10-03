module ut_mesh_remapping_trace_line_grid_basic

  ! Unit tests for mesh functions - remapping - trace_line_grid - basic stuff

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash

  implicit none

  private

  public :: type_coinc_ind_grid, no_value, a_grid, b_grid, cx_grid, cy_grid
  public :: old2new_coinc_ind

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

end module ut_mesh_remapping_trace_line_grid_basic