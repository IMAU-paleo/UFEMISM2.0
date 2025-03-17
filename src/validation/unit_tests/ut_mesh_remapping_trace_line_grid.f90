module ut_mesh_remapping_trace_line_grid

  ! Unit tests for mesh functions - remapping - trace_line_grid

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use grid_types, only: type_grid
  use line_tracing_basic
  use line_tracing_grid
  use mpi_basic, only: par
  use remapping_types, only: type_single_row_mapping_matrices

  implicit none

  private

  public :: test_trace_line_grid

contains

  subroutine test_trace_line_grid( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid'
    character(len=1024), parameter         :: test_name_local = 'full'
    character(len=1024)                    :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_grid_westeast          ( test_name, grid)
    call test_trace_line_grid_eastwest          ( test_name, grid)
    call test_trace_line_grid_southnorth        ( test_name, grid)
    call test_trace_line_grid_northsouth        ( test_name, grid)
    call test_trace_line_grid_southwestnortheast( test_name, grid)
    call test_trace_line_grid_northeastsouthwest( test_name, grid)
    call test_trace_line_grid_northwestsoutheast( test_name, grid)
    call test_trace_line_grid_southeastnorthwest( test_name, grid)
    call test_trace_line_grid_single_cell       ( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid

  subroutine test_trace_line_grid_westeast( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_westeast'
    character(len=1024), parameter         :: test_name_local = 'westeast'
    character(len=1024)                    :: test_name
    integer                                :: i,j
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_coinciding, verified_notcoinciding
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! [pq] coincides with a grid line
    ! ===============================

    verified_coinciding = .true.

    do j = 2, grid%ny-1

      q = [grid%x( grid%nx) - grid%dx / 2._dp, grid%y( j) - grid%dx / 2._dp]

      do i = 2, grid%nx-2

        p = [grid%x( i) - grid%dx / 2._dp, grid%y( j) - grid%dx / 2._dp]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_coinciding = verified_coinciding .and. &
            i_left == i-1+ii .and. &
            j_left == j
        end do

        ! Don't count coincidences
        ! ========================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .false.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results (same as before, but now the line integrals
        ! should not be included)
        do ii = 1, single_row%n
          verified_coinciding = verified_coinciding .and. &
            test_eq( single_row%LI_xdy  ( ii), 0._dp) .and. &
            test_eq( single_row%LI_xydy ( ii), 0._dp) .and. &
            test_eq( single_row%LI_mxydx( ii), 0._dp)
        end do

      end do
    end do

    call unit_test( verified_coinciding, trim(test_name) // '/pq_coincides')

    ! [pq] coincides with a grid line
    ! ===============================

    verified_notcoinciding = .true.

    do j = 2, grid%ny-1

      q = [grid%x( grid%nx), grid%y( j)]

      do i = 2, grid%nx-2

        p = [grid%x( i), grid%y( j)]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_notcoinciding = verified_notcoinciding .and. &
            i_left == i-1+ii .and. &
            j_left == j
        end do

      end do
    end do

    call unit_test( verified_notcoinciding, trim(test_name) // '/pq_coincides_not')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_westeast

  subroutine test_trace_line_grid_eastwest( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_eastwest'
    character(len=1024), parameter         :: test_name_local = 'eastwest'
    character(len=1024)                    :: test_name
    integer                                :: i,j
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_coinciding, verified_notcoinciding
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! [pq] coincides with a grid line
    ! ===============================

    verified_coinciding = .true.

    do j = grid%ny-1, 2, -1

      q = [grid%x( 1) + grid%dx / 2._dp, grid%y( j) + grid%dx / 2._dp]

      do i = grid%nx-2, 2, -1

        p = [grid%x( i) + grid%dx / 2._dp, grid%y( j) + grid%dx / 2._dp]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_coinciding = verified_coinciding .and. &
            i_left == i+1-ii .and. &
            j_left == j
        end do

        ! Don't count coincidences
        ! ========================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .false.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results (same as before, but now the line integrals
        ! should not be included)
        do ii = 1, single_row%n
          verified_coinciding = verified_coinciding .and. &
            test_eq( single_row%LI_xdy  ( ii), 0._dp) .and. &
            test_eq( single_row%LI_xydy ( ii), 0._dp) .and. &
            test_eq( single_row%LI_mxydx( ii), 0._dp)
        end do

      end do
    end do

    call unit_test( verified_coinciding, trim(test_name) // '/pq_coincides')

    ! [pq] coincides with a grid line
    ! ===============================

    verified_notcoinciding = .true.

    do j = grid%ny-1, 2, -1

      q = [grid%x( 1), grid%y( j)]

      do i = grid%nx-1, 3, -1

        p = [grid%x( i), grid%y( j)]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_notcoinciding = verified_notcoinciding .and. &
            i_left == i+1-ii .and. &
            j_left == j
        end do

      end do
    end do

    call unit_test( verified_notcoinciding, trim(test_name) // '/pq_coincides_not')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_eastwest

  subroutine test_trace_line_grid_southnorth( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_southnorth'
    character(len=1024), parameter         :: test_name_local = 'southnorth'
    character(len=1024)                    :: test_name
    integer                                :: i,j
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_coinciding, verified_notcoinciding
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! [pq] coincides with a grid line
    ! ===============================

    verified_coinciding = .true.

    do i = 2, grid%nx-1

      q = [grid%x( i) - grid%dx / 2._dp, grid%y( grid%ny) - grid%dx / 2._dp]

      do j = 2, grid%ny-2

        p = [grid%x( i) - grid%dx / 2._dp, grid%y( j) - grid%dx / 2._dp]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_coinciding = verified_coinciding .and. &
            i_left == i-1 .and. &
            j_left == j-1+ii
        end do

        ! Don't count coincidences
        ! ========================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .false.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results (same as before, but now the line integrals
        ! should not be included)
        do ii = 1, single_row%n
          verified_coinciding = verified_coinciding .and. &
            test_eq( single_row%LI_xdy  ( ii), 0._dp) .and. &
            test_eq( single_row%LI_xydy ( ii), 0._dp) .and. &
            test_eq( single_row%LI_mxydx( ii), 0._dp)
        end do

      end do
    end do

    call unit_test( verified_coinciding, trim(test_name) // '/pq_coincides')

    ! [pq] coincides with a grid line
    ! ===============================

    verified_notcoinciding = .true.

    do i = 2, grid%nx-1

      q = [grid%x( i), grid%y( grid%ny)]

      do j = 2, grid%ny-2

        p = [grid%x( i), grid%y( j)]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_notcoinciding = verified_notcoinciding .and. &
            i_left == i .and. &
            j_left == j-1+ii
        end do

      end do
    end do

    call unit_test( verified_notcoinciding, trim(test_name) // '/pq_coincides_not')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_southnorth

  subroutine test_trace_line_grid_northsouth( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_northsouth'
    character(len=1024), parameter         :: test_name_local = 'northsouth'
    character(len=1024)                    :: test_name
    integer                                :: i,j
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_coinciding, verified_notcoinciding
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! [pq] coincides with a grid line
    ! ===============================

    verified_coinciding = .true.

    do i = 2, grid%nx-1

      q = [grid%x( i) + grid%dx / 2._dp, grid%y( 1) + grid%dx / 2._dp]

      do j = 2, grid%ny-2

        p = [grid%x( i) + grid%dx / 2._dp, grid%y( j) + grid%dx / 2._dp]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_coinciding = verified_coinciding .and. &
            i_left == i+1 .and. &
            j_left == j+1-ii
        end do

        ! Don't count coincidences
        ! ========================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .false.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results (same as before, but now the line integrals
        ! should not be included)
        do ii = 1, single_row%n
          verified_coinciding = verified_coinciding .and. &
            test_eq( single_row%LI_xdy  ( ii), 0._dp) .and. &
            test_eq( single_row%LI_xydy ( ii), 0._dp) .and. &
            test_eq( single_row%LI_mxydx( ii), 0._dp)
        end do

      end do
    end do

    call unit_test( verified_coinciding, trim(test_name) // '/pq_coincides')

    ! [pq] coincides with a grid line
    ! ===============================

    verified_notcoinciding = .true.

    do i = 2, grid%nx-1

      q = [grid%x( i), grid%y( 1)]

      do j = 2, grid%ny-2

        p = [grid%x( i), grid%y( j)]

        ! Count coincidences
        ! ==================

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_grid( grid, p, q, single_row, count_coincidences)

        ! Check results
        do ii = 1, single_row%n
          n_left = single_row%index_left( ii)
          i_left = grid%n2ij( n_left,1)
          j_left = grid%n2ij( n_left,2)
          ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
          verified_notcoinciding = verified_notcoinciding .and. &
            i_left == i .and. &
            j_left == j+1-ii
        end do

      end do
    end do

    call unit_test( verified_notcoinciding, trim(test_name) // '/pq_coincides_not')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_northsouth

  subroutine test_trace_line_grid_southwestnortheast( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_southwestnortheast'
    character(len=1024), parameter         :: test_name_local = 'southwestnortheast'
    character(len=1024)                    :: test_name
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_pq_through_b, verified_pq_left_of_b, verified_pq_right_of_b
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! This only works on a square domain
    call assert( test_eq( grid%nx, grid%ny),'This test requires a square domain')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! pq passes through b-grid points
    ! ===============================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_through_b = .true.

    p = [grid%xmin, grid%ymin]
    q = [grid%xmax, grid%ymax]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
      verified_pq_through_b = verified_pq_through_b .and. &
        i_left == ii .and. &
        j_left == ii
    end do

    call unit_test( verified_pq_through_b, trim(test_name) // '/pq_through_b')

    ! pq passes left of the b-grid points
    ! ===================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_left_of_b = .true.

    p = [grid%xmin, grid%ymin + grid%dx / 4._dp]
    q = [grid%xmax - grid%dx / 4._dp, grid%ymax]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      !   floor( real( ii+1,dp) / 2._dp), ',', floor( real( ii+2,dp) / 2._dp), '])'
      verified_pq_left_of_b = verified_pq_left_of_b .and. &
        i_left == floor( real( ii+1,dp) / 2._dp) .and. &
        j_left == floor( real( ii+2,dp) / 2._dp)
    end do

    call unit_test( verified_pq_left_of_b, trim(test_name) // '/pq_left_of_b')

    ! pq passes right of the b-grid points
    ! ====================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_right_of_b = .true.

    p = [grid%xmin + grid%dx / 4._dp, grid%ymin]
    q = [grid%xmax, grid%ymax - grid%dx / 4._dp]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      !   floor( real( ii+2,dp) / 2._dp), ',', floor( real( ii+1,dp) / 2._dp), '])'
      verified_pq_right_of_b = verified_pq_right_of_b .and. &
        i_left == floor( real( ii+2,dp) / 2._dp) .and. &
        j_left == floor( real( ii+1,dp) / 2._dp)
    end do

    call unit_test( verified_pq_right_of_b, trim(test_name) // '/pq_right_of_b')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_southwestnortheast

  subroutine test_trace_line_grid_northeastsouthwest( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_northeastsouthwest'
    character(len=1024), parameter         :: test_name_local = 'northeastsouthwest'
    character(len=1024)                    :: test_name
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_pq_through_b, verified_pq_left_of_b, verified_pq_right_of_b
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! This only works on a square domain
    call assert( test_eq( grid%nx, grid%ny),'This test requires a square domain')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! pq passes through b-grid points
    ! ===============================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_through_b = .true.

    p = [grid%xmax, grid%ymax]
    q = [grid%xmin, grid%ymin]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
      verified_pq_through_b = verified_pq_through_b .and. &
        i_left == grid%nx+1-ii .and. &
        j_left == grid%ny+1-ii
    end do

    call unit_test( verified_pq_through_b, trim(test_name) // '/pq_through_b')

    ! pq passes to the left of the b-grid points
    ! ==========================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_left_of_b = .true.

    p = [grid%xmax, grid%ymax - grid%dx / 4._dp]
    q = [grid%xmin + grid%dx / 4._dp, grid%ymin]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      !   grid%nx - floor( real( ii-1,dp) / 2._dp), ',', grid%nx - floor( real( ii,dp) / 2._dp), '])'
      verified_pq_left_of_b = verified_pq_left_of_b .and. &
        i_left == grid%nx - floor( real( ii-1,dp) / 2._dp) .and. &
        j_left == grid%nx - floor( real( ii  ,dp) / 2._dp)
    end do

    call unit_test( verified_pq_left_of_b, trim(test_name) // '/pq_left_of_b')

    ! pq passes to the right of the b-grid points
    ! ===========================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_right_of_b = .true.

    p = [grid%xmax - grid%dx / 4._dp, grid%ymax]
    q = [grid%xmin, grid%ymin + grid%dx / 4._dp]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      !   grid%nx - floor( real( ii,dp) / 2._dp), ',', grid%nx - floor( real( ii-1,dp) / 2._dp), '])'
      verified_pq_right_of_b = verified_pq_right_of_b .and. &
        i_left == grid%nx - floor( real( ii  ,dp) / 2._dp) .and. &
        j_left == grid%nx - floor( real( ii-1,dp) / 2._dp)
    end do

    call unit_test( verified_pq_right_of_b, trim(test_name) // '/pq_right_of_b')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_northeastsouthwest

  subroutine test_trace_line_grid_northwestsoutheast( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_northwestsoutheast'
    character(len=1024), parameter         :: test_name_local = 'northwestsoutheast'
    character(len=1024)                    :: test_name
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_pq_through_b, verified_pq_left_of_b, verified_pq_right_of_b
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! This only works on a square domain
    call assert( test_eq( grid%nx, grid%ny),'This test requires a square domain')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! pq passes through b-grid points
    ! ===============================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_through_b = .true.

    p = [grid%xmin, grid%ymax]
    q = [grid%xmax, grid%ymin]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
      verified_pq_through_b = verified_pq_through_b .and. &
        i_left == ii .and. &
        j_left == grid%ny+1-ii
    end do

    call unit_test( verified_pq_through_b, trim(test_name) // '/pq_through_b')

    ! pq passes left of the b-grid points
    ! ===================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_left_of_b = .true.

    p = [grid%xmin + grid%dx / 4._dp, grid%ymax]
    q = [grid%xmax, grid%ymin + grid%dx / 4._dp]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      !   floor( real( ii+2,dp) / 2._dp), ',', grid%ny - floor( real( ii-1,dp) / 2._dp), '])'
      verified_pq_left_of_b = verified_pq_left_of_b .and. &
        i_left == floor( real( ii+2,dp) / 2._dp) .and. &
        j_left == grid%ny - floor( real( ii-1,dp) / 2._dp)
    end do

    call unit_test( verified_pq_left_of_b, trim(test_name) // '/pq_left_of_b')

    ! pq passes right of the b-grid points
    ! ====================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_right_of_b = .true.

    p = [grid%xmin, grid%ymax - grid%dx / 4._dp]
    q = [grid%xmax - grid%dx / 4._dp, grid%ymin]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      !   floor( real( ii+1,dp) / 2._dp), ',', grid%ny - floor( real( ii,dp) / 2._dp), '])'
      verified_pq_right_of_b = verified_pq_right_of_b .and. &
        i_left == floor( real( ii+1,dp) / 2._dp) .and. &
        j_left == grid%ny - floor( real( ii,dp) / 2._dp)
    end do

    call unit_test( verified_pq_right_of_b, trim(test_name) // '/pq_right_of_b')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_northwestsoutheast

  subroutine test_trace_line_grid_southeastnorthwest( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_southeastnorthwest'
    character(len=1024), parameter         :: test_name_local = 'southeastnorthwest'
    character(len=1024)                    :: test_name
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_pq_through_b, verified_pq_left_of_b, verified_pq_right_of_b
    logical                                :: count_coincidences
    integer                                :: ii, n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! This only works on a square domain
    call assert( test_eq( grid%nx, grid%ny),'This test requires a square domain')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! pq passes through b-grid points
    ! ===============================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_through_b = .true.

    p = [grid%xmax, grid%ymin]
    q = [grid%xmin, grid%ymax]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, ']'
      verified_pq_through_b = verified_pq_through_b .and. &
        i_left == grid%nx+1-ii .and. &
        j_left == ii
    end do

    call unit_test( verified_pq_through_b, trim(test_name) // '/pq_through_b')

    ! pq passes left of the b-grid points
    ! ===================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_left_of_b = .true.

    p = [grid%xmax, grid%ymin + grid%dx / 4._dp]
    q = [grid%xmin + grid%dx / 4._dp, grid%ymax]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      !   grid%nx - floor( real( ii-1,dp) / 2._dp), ',', floor( real( ii+2,dp) / 2._dp), '])'
      verified_pq_left_of_b = verified_pq_left_of_b .and. &
        i_left == grid%nx - floor( real( ii-1,dp) / 2._dp) .and. &
        j_left == floor( real( ii+2,dp) / 2._dp)
    end do

    call unit_test( verified_pq_left_of_b, trim(test_name) // '/pq_left_of_b')

    ! pq passes right of the b-grid points
    ! ====================================

    single_row%n          = 0
    single_row%index_left = 0
    single_row%LI_xdy     = 0._dp
    single_row%LI_mxydx   = 0._dp
    single_row%LI_xydy    = 0._dp

    verified_pq_right_of_b = .true.

    p = [grid%xmax - grid%dx / 4._dp, grid%ymin]
    q = [grid%xmin, grid%ymax - grid%dx / 4._dp]

    count_coincidences = .true.
    call trace_line_grid( grid, p, q, single_row, count_coincidences)

    ! Check results
    do ii = 1, single_row%n
      n_left = single_row%index_left( ii)
      i_left = grid%n2ij( n_left,1)
      j_left = grid%n2ij( n_left,2)
      ! if (par%primary) write(0,*) '[i,j]-left(', ii, ') = [', i_left, ',', j_left, '] (attempt [', &
      ! grid%ny - floor( real( ii,dp) / 2._dp), ',', floor( real( ii+1,dp) / 2._dp), '])'
      verified_pq_right_of_b = verified_pq_right_of_b .and. &
        i_left == grid%ny - floor( real( ii,dp) / 2._dp) .and. &
        j_left == floor( real( ii+1,dp) / 2._dp)
    end do

    call unit_test( verified_pq_right_of_b, trim(test_name) // '/pq_right_of_b')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_southeastnorthwest

  subroutine test_trace_line_grid_single_cell( test_name_parent, grid)
    ! Test the trace_line_grid subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_grid_single_cell'
    character(len=1024), parameter         :: test_name_local = 'single_cell'
    character(len=1024)                    :: test_name
    integer                                :: i,j
    real(dp), dimension(2)                 :: sw, se, ne, nw
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified
    logical                                :: count_coincidences
    integer                                :: n_left, i_left, j_left

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = grid%nx * 2
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    verified = .true.

    do i = 2, grid%nx-1
      do j = 2, grid%ny-1

        sw = [grid%x( i) - grid%dx / 2._dp, grid%y( j) - grid%dx / 2._dp]
        se = [grid%x( i) + grid%dx / 2._dp, grid%y( j) - grid%dx / 2._dp]
        ne = [grid%x( i) + grid%dx / 2._dp, grid%y( j) + grid%dx / 2._dp]
        nw = [grid%x( i) - grid%dx / 2._dp, grid%y( j) + grid%dx / 2._dp]

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        ! Trace the border of this grid cell
        count_coincidences = .true.
        call trace_line_grid( grid, sw, se, single_row, count_coincidences)
        call trace_line_grid( grid, se, ne, single_row, count_coincidences)
        call trace_line_grid( grid, ne, nw, single_row, count_coincidences)
        call trace_line_grid( grid, nw, sw, single_row, count_coincidences)

        ! Check results
        n_left = single_row%index_left( 1)
        i_left = grid%n2ij( n_left,1)
        j_left = grid%n2ij( n_left,2)
        ! if (par%primary) write(0,*) '[i,j]-left = [', i_left, ',', j_left, '] (LI = ', &
        !   single_row%LI_xdy( 1), '; ', grid%dx**2._dp, ')'
        verified = verified .and. &
          single_row%n == 1 .and. &
          i_left == i .and. &
          j_left == j .and. &
          test_tol( single_row%LI_xdy( 1), grid%dx**2._dp, grid%dx**2._dp * 1E-12_dp)

      end do
    end do

    call unit_test( verified, trim(test_name)//'/all')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_grid_single_cell

end module ut_mesh_remapping_trace_line_grid