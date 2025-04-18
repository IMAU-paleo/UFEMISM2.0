module grid_basic

  ! Functions for working with simple square x/y-grids

  use precisions, only: dp
  use grid_types, only: type_grid
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use parameters
  use petsc_basic, only: mat_CSR2petsc
  use reallocate_mod, only: reallocate
  use interpolation, only: linint_points
  use projections, only: inverse_oblique_sg_projection
  use mpi_distributed_memory, only: partition_list
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, finalise_matrix_CSR_dist, &
    add_entry_CSR_dist, deallocate_matrix_CSR_dist
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary, distribute_gridded_data_from_primary
  use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_MIN, MPI_MAX

  implicit none

contains

! == Basic square grid functionality

  subroutine setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid, lambda_M, phi_M, beta_stereo)
    !< Set up a square grid that covers the specified domain

    ! In/output variables:
    character(len=256), intent(in   ) :: name
    real(dp),           intent(in   ) :: xmin, xmax, ymin, ymax        ! [m] Domain
    real(dp),           intent(in   ) :: dx                            ! [m] Resolution
    type(type_grid),    intent(  out) :: grid
    real(dp), optional, intent(in   ) :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    character(len=256), parameter :: routine_name = 'setup_square_grid'
    real(dp)                      :: xmid, ymid
    integer                       :: nx_pos, ny_pos
    integer                       :: i,j,ii,jj

    ! Add routine to path
    call init_routine( routine_name)

    ! Name
    grid%name = name

    ! Resolution
    grid%dx = dx

    ! Determine how many grid points are needed to fully cover the domain

    xmid = (xmin + xmax) / 2._dp
    nx_pos = 0
    do while (xmid + real( nx_pos,dp) * grid%dx + grid%dx / 2._dp < xmax)
      nx_pos = nx_pos + 1
    end do
    grid%nx = 1 + 2 * nx_pos

    ymid = (ymin + ymax) / 2._dp
    ny_pos = 0
    do while (ymid + real( ny_pos,dp) * grid%dx + grid%dx / 2._dp < ymax)
      ny_pos = ny_pos + 1
    end do
    grid%ny = 1 + 2 * ny_pos

    ! Fill in x and y
    allocate( grid%x( grid%nx))
    allocate( grid%y( grid%ny))

    do i = 1, grid%nx
      ii = i - (nx_pos+1)
      grid%x( i) = xmid + real( ii,dp) * grid%dx
    end do

    do j = 1, grid%ny
      jj = j - (ny_pos+1)
      grid%y( j) = ymid + real( jj,dp) * grid%dx
    end do

    ! Calculate secondary grid geometry data
    call calc_secondary_grid_data( grid, lambda_M, phi_M, beta_stereo)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_square_grid

  subroutine deallocate_grid( grid)
    !< deallocate memory for a grid object

    ! In/output variables:
    type(type_grid), intent(inout) :: grid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'deallocate_grid'

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( grid%x   )) deallocate( grid%x   )
    if (allocated( grid%y   )) deallocate( grid%y   )
    if (allocated( grid%lon )) deallocate( grid%lon )
    if (allocated( grid%lat )) deallocate( grid%lat )
    if (allocated( grid%ij2n)) deallocate( grid%ij2n)
    if (allocated( grid%n2ij)) deallocate( grid%n2ij)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_grid

  subroutine check_if_grids_are_identical( grid1, grid2, isso)
    !< Check if two grids are identical

    ! In/output variables:
    type(type_grid), intent(in   ) :: grid1, grid2
    logical,         intent(  out) :: isso

    ! Local variables:
    character(len=256), parameter :: routine_name = 'check_if_grids_are_identical'
    real(dp), parameter           :: rtol = 1E-9_dp
    integer                       :: i,j

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    ! Size
    if (grid1%nx /= grid2%nx .or. grid1%ny /= grid2%ny) then
      isso = .false.
      goto 888
    end if

    ! Coordinates
    do i = 1, grid1%nx
      if (abs(grid1%x( i) - grid2%x( i)) > abs(rtol * max(grid1%x(i), grid2%x(i)))) then
        isso = .false.
        goto 888
      end if
    end do
    do j = 1, grid1%ny
      if (abs(grid1%y( j) - grid2%y( j)) > abs(rtol * max(grid1%y(j), grid2%y(j)))) then
        isso = .false.
        goto 888
      end if
    end do

    ! Finalise routine path
888 call finalise_routine( routine_name)

  end subroutine check_if_grids_are_identical

  subroutine calc_secondary_grid_data( grid, lambda_M, phi_M, beta_stereo)
    !< Calculate secondary geometry data for a square grid

    ! In/output variables:
    type(type_grid),    intent(inout) :: grid
    real(dp), optional, intent(in   ) :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    character(len=256), parameter :: routine_name = 'calc_secondary_grid_data'
    integer                       :: i,j
    real(dp), parameter           :: tol = 1E-9_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Resolution
    grid%dx = abs( grid%x( 2) - grid%x( 1))

#if (DO_ASSERTIONS)
    ! Safety
    do i = 1, grid%nx-1
      if (1._dp - abs(grid%x( i+1) - grid%x( i)) / grid%dx > 1E-6_dp) then
        call crash( trim( grid%name) // '" has an irregular x-dimension!')
      end if
    end do
    do j = 1, grid%ny-1
      if (1._dp - abs(grid%y( j+1) - grid%y( j)) / grid%dx > 1E-6_dp) then
        call crash( trim( grid%name) // '" has an irregular y-dimension!')
      end if
    end do
#endif

    ! Domain size
    grid%xmin = minval( grid%x)
    grid%xmax = maxval( grid%x)
    grid%ymin = minval( grid%y)
    grid%ymax = maxval( grid%y)

    ! Tolerance; points lying within this distance of each other are treated as identical
    grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

    ! Conversion tables for grid-form vs. vector-form data
    call calc_field_to_vector_form_translation_tables( grid)

    ! Parallelisation
    call setup_grid_parallelisation( grid)

    ! Lon/lat-coordinates
    if (present( lambda_M) .or. present( phi_M) .or. present( beta_stereo)) then

      ! Safety
      if (.not. (present( lambda_M) .and. present( phi_M) .and. present( beta_stereo))) then
        call crash('need lambda_M, phi_M, and beta_stereo!')
      end if

      ! allocate memory
      allocate( grid%lon( grid%nx, grid%ny))
      allocate( grid%lat( grid%nx, grid%ny))

      ! Calculate lon/lat-coordinates for each grid point using the provided oblique stereographic projection parameters
      do i = 1, grid%nx
      do j = 1, grid%ny
        call inverse_oblique_sg_projection( grid%x( i), grid%y( j), lambda_M, phi_M, beta_stereo, grid%lon( i,j), grid%lat( i,j))
      end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_secondary_grid_data

  subroutine calc_matrix_operators_grid( grid, M_ddx_CSR, M_ddy_CSR)
    !< Calculate matrix operators for partial derivatives on a regular grid
    !< (needed for conservative remapping)

    ! In/output variables:
    type(type_grid),                 intent(in   ) :: grid
    type(type_sparse_matrix_CSR_dp), intent(  out) :: M_ddx_CSR, M_ddy_CSR

    ! Local variables:
    character(len=256), parameter   :: routine_name = 'calc_matrix_operators_grid'
    integer                         :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                         :: row, i, j, col
    real(dp)                        :: valpos, valneg

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Matrix size
    ncols           = grid%n      ! from
    ncols_loc       = grid%n_loc
    nrows           = grid%n      ! to
    nrows_loc       = grid%n_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( M_ddx_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( M_ddy_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! == Fill matrix coefficients
    ! ===========================

    valpos =  1._dp / (2._dp * grid%dx)
    valneg = -1._dp / (2._dp * grid%dx)

    do row = grid%n1, grid%n2

      ! Grid indices
      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)

      ! Skip the border
      if (i == 1 .or. i == grid%nx .or. j == 1 .or. j == grid%ny) then
        M_ddx_CSR%ptr( row+1) = M_ddx_CSR%ptr( row)
        M_ddy_CSR%ptr( row+1) = M_ddy_CSR%ptr( row)
        cycle
      end if

    ! == d/dx

      ! Left
      col = grid%ij2n( i-1,j)
      M_ddx_CSR%nnz = M_ddx_CSR%nnz + 1
      M_ddx_CSR%ind( M_ddx_CSR%nnz) = col
      M_ddx_CSR%val( M_ddx_CSR%nnz) = valneg

      ! Left
      col = grid%ij2n( i+1,j)
      M_ddx_CSR%nnz = M_ddx_CSR%nnz + 1
      M_ddx_CSR%ind( M_ddx_CSR%nnz) = col
      M_ddx_CSR%val( M_ddx_CSR%nnz) = valpos

      ! Ptr
      M_ddx_CSR%ptr( row+1) = M_ddx_CSR%ptr( row) + 2

    ! == d/dy

      ! Left
      col = grid%ij2n( i,j-1)
      M_ddy_CSR%nnz = M_ddy_CSR%nnz + 1
      M_ddy_CSR%ind( M_ddy_CSR%nnz) = col
      M_ddy_CSR%val( M_ddy_CSR%nnz) = valneg

      ! Left
      col = grid%ij2n( i,j+1)
      M_ddy_CSR%nnz = M_ddy_CSR%nnz + 1
      M_ddy_CSR%ind( M_ddy_CSR%nnz) = col
      M_ddy_CSR%val( M_ddy_CSR%nnz) = valpos

      ! Ptr
      M_ddy_CSR%ptr( row+1) = M_ddy_CSR%ptr( row) + 2

    end do ! DO row = grid%n1, grid%n2

    call finalise_matrix_CSR_dist( M_ddx_CSR)
    call finalise_matrix_CSR_dist( M_ddy_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_grid

  subroutine calc_field_to_vector_form_translation_tables( grid)
    !< Calculate grid-cell-to-matrix-row translation tables

    ! In/output variables
    type(type_grid), intent(inout) :: grid

    ! Local variables:
    character(len=256), parameter :: routine_name = 'calc_field_to_vector_form_translation_tables'
    integer                       :: i,j,n

    ! Add routine to path
    call init_routine( routine_name)

    ! Total number of grid cells
    grid%n = grid%nx * grid%ny

    ! Allocate memory
    if (allocated( grid%ij2n)) deallocate( grid%ij2n)
    if (allocated( grid%n2ij)) deallocate( grid%n2ij)
    allocate( grid%ij2n( grid%nx, grid%ny), source = 0)
    allocate( grid%n2ij( grid%n , 2      ), source = 0)

    ! Fill in tables
    n = 0
    do i = 1, grid%nx
    do j = 1, grid%ny
      n = n + 1
      grid%ij2n( i,j) = n
      grid%n2ij( n,:) = [i,j]
    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_field_to_vector_form_translation_tables

  subroutine setup_grid_parallelisation( grid)

    ! In/output variables
    type(type_grid), intent(inout) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_grid_parallelisation'
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Parallelisation domains
    call partition_list( grid%n, par%i, par%n, grid%n1, grid%n2)
    grid%n_loc = grid%n2 + 1 - grid%n1

    ! Parallel array info
    grid%pai%n  = grid%n

    grid%pai%i1    = grid%n1
    grid%pai%i2    = grid%n2
    grid%pai%n_loc = grid%n_loc

    call MPI_ALLREDUCE( grid%pai%i1, grid%pai%i1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( grid%pai%i2, grid%pai%i2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    grid%pai%n_node = grid%pai%i2_node + 1 - grid%pai%i1_node

    grid%pai%i1_nih = grid%pai%i1_node
    grid%pai%i2_nih = grid%pai%i2_node
    grid%pai%n_nih  = grid%pai%n_node

    grid%pai%i1_hle = 0
    grid%pai%i2_hle = -1
    grid%pai%n_hle  = 0

    grid%pai%i1_hli = 0
    grid%pai%i2_hli = -1
    grid%pai%n_hli  = 0

    grid%pai%i1_hre = 0
    grid%pai%i2_hre = -1
    grid%pai%n_hre  = 0

    grid%pai%i1_hri = 0
    grid%pai%i2_hri = -1
    grid%pai%n_hri  = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_grid_parallelisation

! == Calculate contour lines and polygons from gridded data (for mesh generation)

  subroutine calc_grid_contour_as_line( grid, d, f, line, mask)
    !< Calculate a contour line at level f for data d on a square grid.
    !< Generate the contour line in UFEMISM line-segment format (i.e. unordered
    !< individual line segments).

    ! In/output variables
    type(type_grid),                         intent(in   ) :: grid
    real(dp), dimension(:,:  ),              intent(in   ) :: d
    real(dp),                                intent(in   ) :: f
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: line
    logical,  dimension(:,:  ), optional,    intent(in   ) :: mask

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_grid_contour_as_line'
    real(dp), parameter                     :: tol = 1E-5_dp
    real(dp), dimension(:,:  ), allocatable :: d_scaled
    integer                                 :: i,j
    logical,  dimension(:,:  ), allocatable :: mask_loc
    integer                                 :: n_max, n
    real(dp)                                :: d_sw, d_nw, d_se, d_ne
    real(dp)                                :: xw, xe, xs, xn, yw, ye, ys, yn
    logical                                 :: do_cross_w, do_cross_e, do_cross_s, do_cross_n

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d,1) /= grid%nx .or. size( d,2) /= grid%ny) call crash('d is not nx-by-ny!')

    ! Trivial case: if all values of d are greater or smaller than f,
    ! the contour line is empty
    if (minval( d) >= f .or. maxval( d) <= f) then
      allocate( line( 0,0))
      call finalise_routine( routine_name)
      return
    end if

    ! Shift d so the contour lies at d_scaled = 0
    allocate( d_scaled( grid%nx, grid%ny))
    d_scaled = d - f

    ! Set the mask to optionally skip certain grid cells
    allocate( mask_loc( grid%nx, grid%ny))
    if (present( mask)) then
      mask_loc = mask
    else
      mask_loc =  .true.
    end if

    n_max = 1000
    allocate( line( n_max, 4))

    n = 0
    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      ! Skip this grid cell if told so
      if ((.not. mask_loc( i  ,j  )) .and. &
          (.not. mask_loc( i  ,j+1)) .and. &
          (.not. mask_loc( i+1,j  )) .and. &
          (.not. mask_loc( i+1,j+1))) cycle

      ! Extend allocated memory IF needed
      if (n > n_max - 10) then
        n_max = n + 1000
        call reallocate( line, n_max, 4)
      end if

      ! The four corners of the b-grid cell
      d_sw = d_scaled( i  ,j  )
      d_nw = d_scaled( i  ,j+1)
      d_se = d_scaled( i+1,j  )
      d_ne = d_scaled( i+1,j+1)

      xw = grid%x( i  )
      xe = grid%x( i+1)
      ys = grid%y( j  )
      yn = grid%y( j+1)

      ! If all four corners are above/below the level, no line here
      if ((d_sw >= 0._dp .and. d_nw >= 0._dp .and. d_se >= 0._dp .and. d_ne >= 0._dp) .or. &
          (d_sw <= 0._dp .and. d_nw <= 0._dp .and. d_se <= 0._dp .and. d_ne <= 0._dp)) then
        cycle
      end if

      ! Add tolerances to keep line lengths finite
      if (d_sw >= 0._dp) then
        d_sw = max( d_sw, tol)
      else
        d_sw = min( d_sw, tol)
      end if
      if (d_nw >= 0._dp) then
        d_nw = max( d_nw, tol)
      else
        d_nw = min( d_nw, tol)
      end if
      if (d_se >= 0._dp) then
        d_se = max( d_se, tol)
      else
        d_se = min( d_se, tol)
      end if
      if (d_ne >= 0._dp) then
        d_ne = max( d_ne, tol)
      else
        d_ne = min( d_ne, tol)
      end if

      ! Find boundary crossings

      do_cross_w = .false.
      do_cross_e = .false.
      do_cross_s = .false.
      do_cross_n = .false.

      if (d_sw * d_nw < 0._dp) then
        ! The contour crosses the western boundary of this b-grid cell
        do_cross_w = .true.
        yw = linint_points( ys, yn, d_sw, d_nw, 0._dp)
      end if

      if (d_se * d_ne < 0._dp) then
        ! The contour crosses the eastern boundary of this b-grid cell
        do_cross_e = .true.
        ye = linint_points( ys, yn, d_se, d_ne, 0._dp)
      end if

      if (d_nw * d_ne < 0._dp) then
        ! The contour crosses the northern boundary of this b-grid cell
        do_cross_n = .true.
        xn = linint_points( xw, xe, d_nw, d_ne, 0._dp)
      end if

      if (d_sw * d_se < 0._dp) then
        ! The contour crosses the southern boundary of this b-grid cell
        do_cross_s = .true.
        xs = linint_points( xw, xe, d_sw, d_se, 0._dp)
      end if

      ! Add line segments
      n = n + 1
      if     (do_cross_w) then
        if     (do_cross_e) then
          ! From west to east
          line( n,:) = [xw,yw,xe,ye]
        elseif (do_cross_s) then
          ! From west to south
          line( n,:) = [xw,yw,xs,ys]
        elseif (do_cross_n) then
          ! From west to north
          line( n,:) = [xw,yw,xn,yn]
        else
          call crash('found only a crossing at the western boundary!')
        end if
      elseif (do_cross_e) then
        IF     (do_cross_s) then
          ! From east to south
          line( n,:) = [xe,ye,xs,ys]
        elseif (do_cross_n) then
          ! From east to north
          line( n,:) = [xe,ye,xn,yn]
        else
          call crash('found only a crossing at the eastern boundary!')
        end if
      elseif (do_cross_s) then
        IF     (do_cross_n) then
          ! From south to north
          line( n,:) = [xs,ys,xn,yn]
        else
          call crash('found only a crossing at the southern boundary!')
        end if
      elseif (do_cross_n) then
          call crash('found only a crossing at the northern boundary!')
      else
        call crash('whaa!')
      end if

    end do
    end do

    ! Crop memory
    call reallocate( line, n, 4)

    ! Clean up after yourself
    deallocate( d_scaled)
    deallocate( mask_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grid_contour_as_line

  subroutine poly2line( poly, line)
    !< Convert a multi-polygon to UFEMISM line-segment format (i.e. unordered
    !< individual line segments).

    ! In/output variables
    real(dp), dimension(:,:  ),              intent(in   ) :: poly
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: line

    ! Local variables:
    character(len=256), parameter :: routine_name = 'poly2line'
    integer                       :: n,i1,i2

    ! Add routine to path
    call init_routine( routine_name)

    n = size( poly,1)

    ! allocate memory for line segments
    allocate( line( n,4))

    ! Convert polygon to line-segment format
    do i1 = 1, n
      i2 = i1 + 1
      if (i2 == n+1) i2 = 1
      line( i1,:) = [poly( i1,1), poly( i1,2), poly( i2,1), poly( i2,2)]
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine poly2line

  subroutine calc_grid_mask_as_polygons( grid, mask, poly_mult)
    !< Calculate a set of polygon enveloping all TRUE-valued mask cells

    ! In/output variables
    type(type_grid),                         intent(in   ) :: grid
    logical,  dimension(:,:  ),              intent(in   ) :: mask
    real(dp), dimension(:,:  ), allocatable, intent(  out) :: poly_mult

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_grid_mask_as_polygons'
    logical,  dimension(:,:), allocatable :: mask_loc
    integer                               :: i,j
    real(dp), dimension(:,:), allocatable :: poly
    integer                               :: n_poly, n_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( mask,1) /= grid%nx .or. size( mask,2) /= grid%ny) call crash('incorrect data dimensions!')

    ! Trivial case for no TRUE values at all
    if (.not. any( mask)) then
      allocate( poly_mult( 0,0))
      call finalise_routine( routine_name)
      return
    end if

    ! Make a local copy of the logical mask
    allocate( mask_loc( grid%nx, grid%ny))
    mask_loc = mask

    ! Initialise poly_mult and poly
    allocate( poly_mult( grid%nx*grid%ny,2))
    allocate( poly(      grid%nx*grid%ny,2))
    n_tot = 0

    ! Calculate polygons for all TRUE regions of the mask
    do i = 1, grid%nx
    do j = 1, grid%ny

      if (mask_loc( i,j)) then
        ! Found a seed for a TRUE region

        ! Calculate a polygon enveloping this TRUE region, and
        ! remove the region from the mask
        call calc_grid_mask_as_polygon( grid, mask_loc, i, j, poly, n_poly)

        ! Add this polygon to poly_mult
        poly_mult( n_tot+1,1) = real( n_poly,dp)
        poly_mult( n_tot+1,2) = 0._dp
        poly_mult( n_tot+2:n_tot+1+n_poly,:) = poly( 1:n_poly,:)
        n_tot = n_tot + 1 + n_poly

      end if ! if (mask_loc( i,j)) then

    end do
    end do

    ! Crop memory
    call reallocate( poly_mult, n_tot, 2)

    ! Clean up after yourself
    deallocate( mask_loc)
    deallocate( poly)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grid_mask_as_polygons

  subroutine calc_grid_mask_as_polygon( grid, mask, i0, j0, poly, n_poly)
    !< Calculate a polygon enveloping the set of TRUE-valued mask cells around i0,j0,
    !< and remove that set of grid cells from the mask

    ! In/output variables
    type(type_grid),          intent(in   ) :: grid
    logical,  dimension(:,:), intent(inout) :: mask
    integer,                  intent(in   ) :: i0,j0
    real(dp), dimension(:,:), intent(  out) :: poly
    integer,                  intent(  out) :: n_poly

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_grid_mask_as_polygon'
    logical,  dimension(:,:), allocatable :: mask_ext
    integer,  dimension(:,:), allocatable :: map
    integer,  dimension(:,:), allocatable :: stack
    integer                               :: stackN
    integer                               :: i,j,ii,jj,i_sw,j_sw
    character(len=256)                    :: dir, dir_prev
    integer                               :: it

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Safety
    if (size( mask,1) /= grid%nx .or. size( mask,2) /= grid%ny) call crash('incorrect data dimensions!')
    if (.not. mask( i0,j0)) call crash('seed at i0,j0 is not TRUE!')
#endif

    ! Pad a row of FALSEs around the mask so the set of TRUEs around i0,j0 is really an island
    allocate( mask_ext( 0:grid%nx+1, 0:grid%ny+1), source = .false.)
    mask_ext( 1:grid%nx, 1:grid%ny) = mask

    ! Use a flood-fill algorithm to find the map of same-valued grid cells around mask cell i0,j0
    allocate( map(   0:grid%nx+1, 0:grid%ny+1) , source = 0)
    allocate( stack( (grid%nx+2)*(grid%ny+2),2))

    map( i0,j0) = 1
    stackN = 1
    stack( 1,:) = [i0,j0]

    i_sw = 0
    j_sw = 0

    do while (stackN > 0)

      ! Take the last element from the stack
      i = stack( stackN,1)
      j = stack( stackN,2)
      stackN = stackN - 1

      ! Mark it as mapped
      map( i,j) = 2

      ! Remove it from the input mask
      mask( i,j) = .false.

      ! Add all non-mapped, non-stacked, TRUE-valued neighbours to the stack
      DO ii = max( 0, i-1), min( grid%nx+1, i+1)
      DO jj = max( 0, j-1), min( grid%ny+1, j+1)
        if (ii /= i .and. jj /= j) cycle ! Don't include diagonal neighbours
        if (map( ii,jj) == 0 .and. mask_ext( ii,jj)) then
          ! Add this neighbour to the stack
          stackN = stackN + 1
          stack( stackN,:) = [ii,jj]
          ! Mark this neighbour on the map as stacked
          map( ii,jj) = 1
        end if
      end do
      end do

      ! If it is a southwest corner, save it as a starting point for the outline tracer
      if (.not. mask_ext( i-1,j) .and. .not. mask_ext( i,j-1)) then
        i_sw = i
        j_sw = j
      end if

    end do ! do while (stackN > 0)
    ! Safety
    if (i_sw == 0 .or. j_sw == 0) call crash('couldnt find starting SW corner!')

    ! Start at the southwest corner we found earlier
    ii = i_sw
    jj = j_sw
    n_poly = 2
    poly( 1,:) = [grid%x( ii) - grid%dx/2._dp, grid%y( jj) - grid%dx/2._dp]
    poly( 2,:) = [grid%x( ii) + grid%dx/2._dp, grid%y( jj) - grid%dx/2._dp]
    dir = 'east'

    ! Trace the outline of the mapped grid cells.
    it = 0
    do while (.true.)

      ! Safety
      it = it + 1
      if (it > grid%nx*grid%ny) call crash('outline tracer got stuck!')

      ! cycle direction
      dir_prev = dir

      ! Check which way we go next
      if     (dir_prev == 'east') then
        ! We were going east, so [ii,jj] is to the north of us, and we can
        ! continue north, east, or south

        if     (map( ii+1,jj) == 0) then
          ! Go north
          dir = 'north'
        elseif (map( ii+1,jj-1) == 0) then
          ! Go east
          dir = 'east'
          ii = ii+1
        elseif (map( ii,jj-1) == 0) then
          ! Go south
          dir = 'south'
          ii = ii+1
          jj = jj-1
        else
          call crash('outline tracer got stuck while going east!')
        end if

      elseif (dir_prev == 'north') then
        ! We were going north, so [ii,jj] is to the west of us, and we can
        ! continue west, north, or east

        if     (map( ii,jj+1) == 0) then
          ! Go west
          dir = 'west'
        elseif (map( ii+1,jj+1) == 0) then
          ! Go north
          dir = 'north'
          jj = jj+1
        elseif (map( ii+1,jj) == 0) then
          ! Go east
          dir = 'east'
          ii = ii+1
          jj = jj+1
        else
          call crash('outline tracer got stuck while going north!')
        end if

      elseif (dir_prev == 'west') then
        ! We were going west, so [ii,jj] is to the south of us, and we can
        ! continue south, west, or north

        if     (map( ii-1,jj) == 0) then
          ! Go south
          dir = 'south'
        elseif (map( ii-1,jj+1) == 0) then
          ! Go west
          dir = 'west'
          ii = ii-1
        elseif (map( ii,jj+1) == 0) then
          ! Go north
          dir = 'north'
          ii = ii-1
          jj = jj+1
        else
          call crash('outline tracer got stuck while going north!')
        end if

      elseif (dir_prev == 'south') then
        ! We were going south, so [ii,jj] is to the east of us, and we can
        ! continue east, south, or west

        if     (map( ii,jj-1) == 0) then
          ! Go east
          dir = 'east'
        elseif (map( ii-1,jj-1) == 0) then
          ! Go south
          dir = 'south'
          jj = jj-1
        elseif (map( ii-1,jj) == 0) then
          ! Go west
          dir = 'west'
          ii = ii-1
          jj = jj-1
        else
          call crash('outline tracer got stuck while going north!')
        end if

      else
        call crash('unknown dir_prev "' // trim( dir_prev) // '"!')
      end if

      ! Add new vertex to the polygon
      n_poly = n_poly+1
      if     (dir == 'east') then
        poly( n_poly,:) = poly( n_poly-1,:) + [grid%dx, 0._dp]
      elseif (dir == 'north') then
        poly( n_poly,:) = poly( n_poly-1,:) + [0._dp, grid%dx]
      elseif (dir == 'west') then
        poly( n_poly,:) = poly( n_poly-1,:) + [-grid%dx, 0._dp]
      elseif (dir == 'south') then
        poly( n_poly,:) = poly( n_poly-1,:) + [0._dp, -grid%dx]
      else
        call crash('unknown dir "' // trim( dir) // '"!')
      end if

      ! If we've reached the starting point again, we're done
      if (norm2( poly( n_poly,:) - poly( 1,:)) < grid%dx / 10._dp) then
        ! Don't double count
        n_poly = n_poly-1
        exit
      end if

    end do ! do while (.true.)

    ! Clean up after yourself
    deallocate( mask_ext)
    deallocate( map)
    deallocate( stack)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grid_mask_as_polygon

end module grid_basic
