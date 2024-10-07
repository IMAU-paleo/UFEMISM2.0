module line_tracing_basic

  ! Some very basic functionality for the line tracing algorithms

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash
  use remapping_types, only: type_single_row_mapping_matrices

  implicit none

  private

  public :: add_integrals_to_single_row
  public :: type_coinc_ind_grid, no_value, a_grid, b_grid, cx_grid, cy_grid
  public :: old2new_coinc_ind, coinc_ind2aij_in, coinc_ind2bij_on, coinc_ind2cxij_on, coinc_ind2cyij_on

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

  !> Add the values for a single row of the three line-integral matrices
  subroutine add_integrals_to_single_row( single_row, index_left, &
    LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

    ! In/output variables
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    integer,                                intent(in)    :: index_left
    real(dp),                               intent(in)    :: LI_xdy, LI_mxydx, LI_xydy
    logical,                                intent(in)    :: coincides, count_coincidences

    ! Local variables:
    logical :: do_add_integrals, is_listed
    integer :: i, i_add

    ! Check whether we actually need to add the line integrals
    do_add_integrals = .true.
    if (coincides .and. (.not. count_coincidences)) do_add_integrals = .false.

    ! Check if an entry from this left-hand vertex is already listed
    is_listed = .false.
    i_add     = 0

    do i = 1, single_row%n
      if (single_row%index_left( i) == index_left) then
        is_listed = .true.
        i_add     = i
        exit
      end if
    end do
    if (.not. is_listed) then
      single_row%n = single_row%n + 1
      i_add = single_row%n
    end if

    ! Add data
    single_row%index_left( i_add) = index_left
    if (do_add_integrals) then
      single_row%LI_xdy(   i_add) = single_row%LI_xdy(   i_add) + LI_xdy
      single_row%LI_mxydx( i_add) = single_row%LI_mxydx( i_add) + LI_mxydx
      single_row%LI_xydy(  i_add) = single_row%LI_xydy(  i_add) + LI_xydy
    end if

    ! if necessary, extend memory
    if (single_row%n > single_row%n_max - 10) call extend_single_row_memory( single_row, 100)

  end subroutine add_integrals_to_single_row

  !> Extend memory for a single row of the three line-integral matrices
  subroutine extend_single_row_memory( single_row, n_extra)

    ! In/output variables
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    integer,                                intent(in)    :: n_extra

    ! Local variables:
    integer                             :: n
    integer,  dimension(:), allocatable :: index_left_temp
    real(dp), dimension(:), allocatable :: LI_xdy_temp, LI_mxydx_temp, LI_xydy_temp

    n = single_row%n

    ! allocate temporary memory
    allocate( index_left_temp( n))
    allocate( LI_xdy_temp(     n))
    allocate( LI_mxydx_temp(   n))
    allocate( LI_xydy_temp(    n))

    ! Copy data to temporary memory
    index_left_temp = single_row%index_left( 1:n)
    LI_xdy_temp     = single_row%LI_xdy(     1:n)
    LI_mxydx_temp   = single_row%LI_mxydx(   1:n)
    LI_xydy_temp    = single_row%LI_xydy(    1:n)

    ! deallocate memory
    deallocate( single_row%index_left)
    deallocate( single_row%LI_xdy    )
    deallocate( single_row%LI_mxydx  )
    deallocate( single_row%LI_xydy   )

    ! allocate new, extended memory
    single_row%n_max = single_row%n_max + n_extra
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! Copy data back from temporary memory
    single_row%index_left( 1:n) = index_left_temp
    single_row%LI_xdy(     1:n) = LI_xdy_temp
    single_row%LI_mxydx(   1:n) = LI_mxydx_temp
    single_row%LI_xydy(    1:n) = LI_xydy_temp

  end subroutine extend_single_row_memory

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

  function coinc_ind2aij_in( coinc_ind) result( aij_in)
    type(type_coinc_ind_grid), intent(in) :: coinc_ind
    integer, dimension(2) :: aij_in

    if (coinc_ind%grid == a_grid) then
      aij_in = [coinc_ind%i, coinc_ind%j]
    else
      aij_in = [0,0]
    end if

  end function coinc_ind2aij_in

  function coinc_ind2bij_on( coinc_ind) result( bij_on)
    type(type_coinc_ind_grid), intent(in) :: coinc_ind
    integer, dimension(2) :: bij_on

    if (coinc_ind%grid == b_grid) then
      bij_on = [coinc_ind%i, coinc_ind%j]
    else
      bij_on = [0,0]
    end if

  end function coinc_ind2bij_on

  function coinc_ind2cxij_on( coinc_ind) result( cxij_on)
    type(type_coinc_ind_grid), intent(in) :: coinc_ind
    integer, dimension(2) :: cxij_on

    if (coinc_ind%grid == cx_grid) then
      cxij_on = [coinc_ind%i, coinc_ind%j]
    else
      cxij_on = [0,0]
    end if

  end function coinc_ind2cxij_on

  function coinc_ind2cyij_on( coinc_ind) result( cyij_on)
    type(type_coinc_ind_grid), intent(in) :: coinc_ind
    integer, dimension(2) :: cyij_on

    if (coinc_ind%grid == cy_grid) then
      cyij_on = [coinc_ind%i, coinc_ind%j]
    else
      cyij_on = [0,0]
    end if

  end function coinc_ind2cyij_on

end module line_tracing_basic
