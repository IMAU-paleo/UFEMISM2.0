module line_tracing_basic

  ! Some very basic functionality for the line tracing algorithms

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash
  use remapping_types, only: type_single_row_mapping_matrices

  implicit none

  private

  public :: add_integrals_to_single_row
  public :: type_coinc_ind_grid, type_coinc_ind_mesh, no_value, a_grid, b_grid, c_grid, cx_grid, cy_grid
  public :: coinc_ind_mesh_old2new, coinc_ind_mesh_new2old

  type type_coinc_ind_grid
    integer :: grid
    integer :: i,j
  end type type_coinc_ind_grid

  type type_coinc_ind_mesh
    integer :: grid
    integer :: i
  end type type_coinc_ind_mesh

  integer, parameter :: no_value = -1
  integer, parameter :: a_grid   = 1
  integer, parameter :: b_grid   = 2
  integer, parameter :: c_grid   = 3
  integer, parameter :: cx_grid  = 4
  integer, parameter :: cy_grid  = 5

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

  subroutine coinc_ind_mesh_old2new( vi_in, ti_on, ei_on, finished, coinc_ind)
    integer, intent(in) :: vi_in, ti_on, ei_on
    logical, intent(in) :: finished
    type(type_coinc_ind_mesh), intent(out) :: coinc_ind
    if (vi_in > 0) then
      call assert( ti_on == 0, 'vi_in > 0 and ti_on > 0')
      call assert( ei_on == 0, 'vi_in > 0 and ei_on > 0')
      coinc_ind%grid = a_grid
      coinc_ind%i = vi_in
    elseif (ti_on > 0) then
      call assert( vi_in == 0, 'ti_on > 0 and vi_in > 0')
      call assert( ei_on == 0, 'ti_on > 0 and ei_on > 0')
      coinc_ind%grid = b_grid
      coinc_ind%i = ti_on
    elseif (ei_on > 0) then
      call assert( vi_in == 0, 'ei_on > 0 and vi_in > 0')
      call assert( ti_on == 0, 'ei_on > 0 and ti_on > 0')
      coinc_ind%grid = c_grid
      coinc_ind%i = ei_on
    elseif (finished) then
      coinc_ind%grid = no_value
      coinc_ind%i = 0
    else
      call crash('vi_in, ti_on, ei_on are all zero but finished = false')
    end if
  end subroutine coinc_ind_mesh_old2new
  subroutine coinc_ind_mesh_new2old( coinc_ind, vi_in, ti_on, ei_on)
    type(type_coinc_ind_mesh), intent(in) :: coinc_ind
    integer, intent(out) :: vi_in, ti_on, ei_on
    vi_in = 0
    ti_on = 0
    ei_on = 0
    select case (coinc_ind%grid)
    case default
      call crash('invalid grid in coincidence indicator')
    case (a_grid)
      vi_in = coinc_ind%i
    case (b_grid)
      ti_on = coinc_ind%i
    case (c_grid)
      ei_on = coinc_ind%i
    case (no_value)
    end select
  end subroutine coinc_ind_mesh_new2old

end module line_tracing_basic
