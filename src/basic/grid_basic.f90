MODULE grid_basic

  ! Data types and subroutines for working with simple square grids

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE parameters
  USE petsc_basic                                            , ONLY: perr, mat_CSR2petsc
  USE reallocate_mod                                         , ONLY: reallocate
  USE math_utilities                                         , ONLY: linint_points, inverse_oblique_sg_projection
  USE mpi_distributed_memory                                 , ONLY: partition_list, distribute_from_master_dp_1D, gather_to_master_dp_1D
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  TYPE type_grid
    ! A square grid

    ! Basic properties
    CHARACTER(LEN=256)                      :: name                          !           A nice name tag
    INTEGER                                 :: nx                            !           Number of grid cells in the x-direction
    INTEGER                                 :: ny                            !           Number of grid cells in the x-direction
    INTEGER                                 :: n                             !           Total number of grid cells (= nx * ny)
    REAL(dp)                                :: dx                            ! [m]       Resolution (square grid, so dy = dx)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: x                             ! [m]       x-coordinates
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: y                             ! [m]       y-coordinates
    REAL(dp)                                :: xmin                          ! [m]       x and y range of the square covered by the grid
    REAL(dp)                                :: xmax                          ! [m]
    REAL(dp)                                :: ymin                          ! [m]
    REAL(dp)                                :: ymax                          ! [m]

    ! Remapping data
    REAL(dp)                                :: tol_dist                      ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: ij2n, n2ij                    !           Conversion table for grid-form vs. vector-form data

    ! Lon/lat-coordinates
    REAL(dp)                                :: lambda_M
    REAL(dp)                                :: phi_M
    REAL(dp)                                :: beta_stereo
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: lon, lat

    ! Parallelisation
    INTEGER                                 :: n1,n2,n_loc                   ! Matrix rows owned by each process

  END TYPE type_grid

CONTAINS

! ===== Subroutines ======
! ========================

! == Basic square grid functionality

  SUBROUTINE setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid, lambda_M, phi_M, beta_stereo)
    ! Set up a square grid that covers the specified domain

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: name
    REAL(dp),                            INTENT(IN)    :: xmin, xmax, ymin, ymax        ! [m] Domain
    REAL(dp),                            INTENT(IN)    :: dx                            ! [m] Resolution
    TYPE(type_grid),                     INTENT(OUT)   :: grid
    REAL(dp),                  OPTIONAL, INTENT(IN)    :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_square_grid'
    REAL(dp)                                           :: xmid, ymid
    INTEGER                                            :: nx_pos, ny_pos
    INTEGER                                            :: i,j,ii,jj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Name
    grid%name = name

    ! Resolution
    grid%dx = dx

    ! Determine how many grid points are needed to fully cover the domain

    xmid = (xmin + xmax) / 2._dp
    nx_pos = 0
    DO WHILE (xmid + REAL( nx_pos,dp) * grid%dx + grid%dx / 2._dp < xmax)
      nx_pos = nx_pos + 1
    END DO
    grid%nx = 1 + 2 * nx_pos

    ymid = (ymin + ymax) / 2._dp
    ny_pos = 0
    DO WHILE (ymid + REAL( ny_pos,dp) * grid%dx + grid%dx / 2._dp < ymax)
      ny_pos = ny_pos + 1
    END DO
    grid%ny = 1 + 2 * ny_pos

    ! Fill in x and y
    ALLOCATE( grid%x( grid%nx))
    ALLOCATE( grid%y( grid%ny))

    DO i = 1, grid%nx
      ii = i - (nx_pos+1)
      grid%x( i) = xmid + REAL( ii,dp) * grid%dx
    END DO

    DO j = 1, grid%ny
      jj = j - (ny_pos+1)
      grid%y( j) = ymid + REAL( jj,dp) * grid%dx
    END DO

    ! Calculate secondary grid geometry data
    CALL calc_secondary_grid_data( grid, lambda_M, phi_M, beta_stereo)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_square_grid

  SUBROUTINE deallocate_grid( grid)
    ! Deallocate memory for a grid object

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (ALLOCATED( grid%x   )) DEALLOCATE( grid%x   )
    IF (ALLOCATED( grid%y   )) DEALLOCATE( grid%y   )
    IF (ALLOCATED( grid%lon )) DEALLOCATE( grid%lon )
    IF (ALLOCATED( grid%lat )) DEALLOCATE( grid%lat )
    IF (ALLOCATED( grid%ij2n)) DEALLOCATE( grid%ij2n)
    IF (ALLOCATED( grid%n2ij)) DEALLOCATE( grid%n2ij)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_grid

  SUBROUTINE check_if_grids_are_identical( grid1, grid2, isso)
    ! Check if two grids are identical

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid1, grid2
    LOGICAL,                             INTENT(OUT)   :: isso

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_if_grids_are_identical'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    isso = .TRUE.

    ! Size
    IF (grid1%nx /= grid2%nx .OR. grid1%ny /= grid2%ny) THEN
      isso = .FALSE.
      RETURN
    END IF

    ! Coordinates
    DO i = 1, grid1%nx
      IF (ABS( 1._dp - (grid1%x( i) / grid2%x( i))) > tol) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO
    DO j = 1, grid1%ny
      IF (ABS( 1._dp - (grid1%y( j) / grid2%y( j))) > tol) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_if_grids_are_identical

  SUBROUTINE calc_secondary_grid_data( grid, lambda_M, phi_M, beta_stereo)
    ! Calculate secondary geometry data for a square grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp),                  OPTIONAL, INTENT(IN)    :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_secondary_grid_data'
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Resolution
    grid%dx   = ABS( grid%x( 2) - grid%x( 1))

    ! Safety
    DO i = 1, grid%nx-1
      IF (1._dp - ABS(grid%x( i+1) - grid%x( i)) / grid%dx > 1E-6_dp) CALL crash( TRIM( grid%name) // '" has an irregular x-dimension!')
    END DO
    DO j = 1, grid%ny-1
      IF (1._dp - ABS(grid%y( j+1) - grid%y( j)) / grid%dx > 1E-6_dp) CALL crash( TRIM( grid%name) // '" has an irregular y-dimension!')
    END DO

    ! Domain size
    grid%xmin = MINVAL( grid%x)
    grid%xmax = MAXVAL( grid%x)
    grid%ymin = MINVAL( grid%y)
    grid%ymax = MAXVAL( grid%y)

    ! Tolerance; points lying within this distance of each other are treated as identical
    grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

    ! Conversion tables for grid-form vs. vector-form data
    CALL calc_field_to_vector_form_translation_tables( grid)

    ! Lon/lat-coordinates
    IF (PRESENT( lambda_M) .OR. PRESENT( phi_M) .OR. PRESENT( beta_stereo)) THEN

      ! Safety
      IF (.NOT. (PRESENT( lambda_M) .AND. PRESENT( phi_M) .AND. PRESENT( beta_stereo))) CALL crash('need lambda_M, phi_M, and beta_stereo!')

      ! Allocate memory
      ALLOCATE( grid%lon( grid%nx, grid%ny))
      ALLOCATE( grid%lat( grid%nx, grid%ny))

      ! Calculate lon/lat-coordinates for each grid point using the provided oblique stereographic projection parameters
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), lambda_M, phi_M, beta_stereo, grid%lon( i,j), grid%lat( i,j))
      END DO
      END DO

    END IF ! IF (PRESENT( lambda_M) .OR. PRESENT( phi_M) .OR. PRESENT( beta_stereo)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_secondary_grid_data

  SUBROUTINE calc_matrix_operators_grid( grid, M_ddx, M_ddy)
    ! Calculate matrix operators for partial derivatives on a regular grid (needed for conservative remapping)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(tMat),                          INTENT(OUT)   :: M_ddx, M_ddy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_grid'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_ddx_CSR, M_ddy_CSR
    INTEGER                                            :: row, i, j, col
    REAL(dp)                                           :: val

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = grid%n      ! from
    ncols_loc       = grid%n_loc
    nrows           = grid%n      ! to
    nrows_loc       = grid%n_loc
    nnz_per_row_est = 3
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( M_ddx_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! == Fill matrix coefficients
  ! ===========================

    DO row = grid%n1, grid%n2

      ! Grid indices
      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)

    ! == d/dx

      IF     (i == 1) THEN
        ! Use a second-order accurate three-point one-sided differencing scheme on the border

        ! i
        col = grid%ij2n( i,j)
        val = -1.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

        ! i+1
        col = grid%ij2n( i+1,j)
        val = 2._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

        ! i+2
        col = grid%ij2n( i+2,j)
        val = -0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

      ELSEIF (i == grid%nx) THEN
        ! Use a second-order accurate three-point one-sided differencing scheme on the border

        ! i
        col = grid%ij2n( i,j)
        val = 1.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

        ! i-1
        col = grid%ij2n( i-1,j)
        val = -2._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

        ! i-2
        col = grid%ij2n( i-2,j)
        val = 0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

      ELSE
        ! Use a second-order accurate two-sided differencing scheme in the interior

        ! i-1
        col = grid%ij2n( i-1,j)
        val = -0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

        ! i+1
        col = grid%ij2n( i+1,j)
        val = 0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx_CSR, row, col, val)

      END IF

    ! == d/dy

      IF     (j == 1) THEN
        ! Use a second-order accurate three-point one-sided differencing scheme on the border

        ! j
        col = grid%ij2n( i,j)
        val = -1.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

        ! j+1
        col = grid%ij2n( i,j+1)
        val = 2._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

        ! j+2
        col = grid%ij2n( i,j+2)
        val = -0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

      ELSEIF (j == grid%ny) THEN
        ! Use a second-order accurate three-point one-sided differencing scheme on the border

        ! j
        col = grid%ij2n( i,j)
        val = 1.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

        ! j-1
        col = grid%ij2n( i,j-1)
        val = -2._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

        ! j-2
        col = grid%ij2n( i,j-2)
        val = 0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

      ELSE
        ! Use a second-order accurate two-sided differencing scheme in the interior

        ! j-1
        col = grid%ij2n( i,j-1)
        val = -0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

        ! j+1
        col = grid%ij2n( i,j+1)
        val = 0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy_CSR, row, col, val)

      END IF

    END DO ! DO row = grid%n1, grid%n2

    ! Convert to PETSc format
    CALL mat_CSR2petsc( M_ddx_CSR, M_ddx)
    CALL mat_CSR2petsc( M_ddy_CSR, M_ddy)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_ddx_CSR)
    CALL deallocate_matrix_CSR_dist( M_ddy_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_grid

  SUBROUTINE calc_field_to_vector_form_translation_tables( grid)
    ! Calculate grid-cell-to-matrix-row translation tables

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                         INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_field_to_vector_form_translation_tables'
    INTEGER                                                :: i,j,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total number of grid cells
    grid%n = grid%nx * grid%ny

    ! Allocate memory
    IF (ALLOCATED( grid%ij2n)) DEALLOCATE( grid%ij2n)
    IF (ALLOCATED( grid%n2ij)) DEALLOCATE( grid%n2ij)
    ALLOCATE( grid%ij2n( grid%nx, grid%ny), source = 0)
    ALLOCATE( grid%n2ij( grid%n , 2      ), source = 0)

    ! Fill in tables
    n = 0
    DO i = 1, grid%nx
    DO j = 1, grid%ny
      n = n + 1
      grid%ij2n( i,j) = n
      grid%n2ij( n,:) = [i,j]
    END DO
    END DO

    ! Parallelisation domains
    CALL partition_list( grid%n, par%i, par%n, grid%n1, grid%n2)
    grid%n_loc = grid%n2 + 1 - grid%n1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_field_to_vector_form_translation_tables

! == Calculate contour lines and polygons from gridded data (for mesh generation)

  SUBROUTINE calc_grid_contour_as_line( grid, d, f, line, mask)
    ! Calculate a contour line at level f for data d on a square grid.
    ! Generate the contour line in UFEMISM line-segment format (i.e. unordered
    ! individual line segments).

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d
    REAL(dp),                                INTENT(IN)    :: f
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: line
    LOGICAL,  DIMENSION(:,:  ), OPTIONAL,    INTENT(IN)    :: mask

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_grid_contour_as_line'
    REAL(dp), PARAMETER                                    :: tol = 1E-5_dp
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_scaled
    INTEGER                                                :: i,j
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE                :: mask_loc
    INTEGER                                                :: n_max, n
    REAL(dp)                                               :: d_sw, d_nw, d_se, d_ne
    REAL(dp)                                               :: xw, xe, xs, xn, yw, ye, ys, yn
    LOGICAL                                                :: do_cross_w, do_cross_e, do_cross_s, do_cross_n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= grid%nx .OR. SIZE( d,2) /= grid%ny) CALL crash('d is not nx-by-ny!')

    ! Trivial case: if all values of d are greater or smaller than f,
    ! the contour line is empty
    IF (MINVAL( d) >= f .OR. MAXVAL( d) <= f) THEN
      ALLOCATE( line( 0,0))
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Shift d so the contour lies at d_scaled = 0
    ALLOCATE( d_scaled( grid%nx, grid%ny))
    d_scaled = d - f

    ! Set the mask to optionally skip certain grid cells
    ALLOCATE( mask_loc( grid%nx, grid%ny))
    IF (PRESENT( mask)) THEN
      mask_loc = mask
    ELSE
      mask_loc =  .TRUE.
    END IF

    n_max = 1000
    ALLOCATE( line( n_max, 4))

    n = 0
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny-1

      ! Skip this grid cell if told so
      IF ((.NOT. mask_loc( i  ,j  )) .AND. &
          (.NOT. mask_loc( i  ,j+1)) .AND. &
          (.NOT. mask_loc( i+1,j  )) .AND. &
          (.NOT. mask_loc( i+1,j+1))) CYCLE

      ! Extend allocated memory IF needed
      IF (n > n_max - 10) THEN
        n_max = n + 1000
        CALL reallocate( line, n_max, 4)
      END IF

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
      IF ((d_sw >= 0._dp .AND. d_nw >= 0._dp .AND. d_se >= 0._dp .AND. d_ne >= 0._dp) .OR. &
          (d_sw <= 0._dp .AND. d_nw <= 0._dp .AND. d_se <= 0._dp .AND. d_ne <= 0._dp)) THEN
        CYCLE
      END IF

      ! Add tolerances to keep line lengths finite
      IF (d_sw >= 0._dp) THEN
        d_sw = MAX( d_sw, tol)
      ELSE
        d_sw = MIN( d_sw, tol)
      END IF
      IF (d_nw >= 0._dp) THEN
        d_nw = MAX( d_nw, tol)
      ELSE
        d_nw = MIN( d_nw, tol)
      END IF
      IF (d_se >= 0._dp) THEN
        d_se = MAX( d_se, tol)
      ELSE
        d_se = MIN( d_se, tol)
      END IF
      IF (d_ne >= 0._dp) THEN
        d_ne = MAX( d_ne, tol)
      ELSE
        d_ne = MIN( d_ne, tol)
      END IF

      ! Find boundary crossings

      do_cross_w = .FALSE.
      do_cross_e = .FALSE.
      do_cross_s = .FALSE.
      do_cross_n = .FALSE.

      IF (d_sw * d_nw < 0._dp) THEN
        ! The contour crosses the western boundary of this b-grid cell
        do_cross_w = .TRUE.
        yw = linint_points( ys, yn, d_sw, d_nw, 0._dp)
      END IF

      IF (d_se * d_ne < 0._dp) THEN
        ! The contour crosses the eastern boundary of this b-grid cell
        do_cross_e = .TRUE.
        ye = linint_points( ys, yn, d_se, d_ne, 0._dp)
      END IF

      IF (d_nw * d_ne < 0._dp) THEN
        ! The contour crosses the northern boundary of this b-grid cell
        do_cross_n = .TRUE.
        xn = linint_points( xw, xe, d_nw, d_ne, 0._dp)
      END IF

      IF (d_sw * d_se < 0._dp) THEN
        ! The contour crosses the southern boundary of this b-grid cell
        do_cross_s = .TRUE.
        xs = linint_points( xw, xe, d_sw, d_se, 0._dp)
      END IF

      ! Add line segments
      n = n + 1
      IF     (do_cross_w) THEN
        IF     (do_cross_e) THEN
          ! From west to east
          line( n,:) = [xw,yw,xe,ye]
        ELSEIF (do_cross_s) THEN
          ! From west to south
          line( n,:) = [xw,yw,xs,ys]
        ELSEIF (do_cross_n) THEN
          ! From west to north
          line( n,:) = [xw,yw,xn,yn]
        ELSE
          CALL crash('found only a crossing at the western boundary!')
        END IF
      ELSEIF (do_cross_e) THEN
        IF     (do_cross_s) THEN
          ! From east to south
          line( n,:) = [xe,ye,xs,ys]
        ELSEIF (do_cross_n) THEN
          ! From east to north
          line( n,:) = [xe,ye,xn,yn]
        ELSE
          CALL crash('found only a crossing at the eastern boundary!')
        END IF
      ELSEIF (do_cross_s) THEN
        IF     (do_cross_n) THEN
          ! From south to north
          line( n,:) = [xs,ys,xn,yn]
        ELSE
          CALL crash('found only a crossing at the southern boundary!')
        END IF
      ELSEIF (do_cross_n) THEN
          CALL crash('found only a crossing at the northern boundary!')
      ELSE
        CALL crash('whaa!')
      END IF

    END DO
    END DO

    ! Crop memory
    CALL reallocate( line, n, 4)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grid_contour_as_line

  SUBROUTINE poly2line( poly, line)
    ! Convert a multi-polygon to UFEMISM line-segment format (i.e. unordered
    ! individual line segments).

    IMPLICIT NONE

    ! In/output variables
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: poly
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: line

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'poly2line'
    INTEGER                                                :: n,i1,i2

    ! Add routine to path
    CALL init_routine( routine_name)

    n = SIZE( poly,1)

    ! Allocate memory for line segments
    ALLOCATE( line( n,4))

    ! Convert polygon to line-segment format
    DO i1 = 1, n
      i2 = i1 + 1
      IF (i2 == n+1) i2 = 1
      line( i1,:) = [poly( i1,1), poly( i1,2), poly( i2,1), poly( i2,2)]
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE poly2line

  SUBROUTINE calc_grid_mask_as_polygons( grid, mask, poly_mult)
    ! Calculate a set of polygon enveloping all TRUE-valued mask cells

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                         INTENT(IN)    :: grid
    LOGICAL,  DIMENSION(:,:  ),              INTENT(IN)    :: mask
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly_mult

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_grid_mask_as_polygons'
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE                :: mask_loc
    INTEGER                                                :: i,j
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: poly
    INTEGER                                                :: n_poly, n_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( mask,1) /= grid%nx .OR. SIZE( mask,2) /= grid%ny) CALL crash('incorrect data dimensions!')

    ! Trivial case for no TRUE values at all
    IF (.NOT. ANY( mask)) THEN
      ALLOCATE( poly_mult( 0,0))
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Make a local copy of the logical mask
    ALLOCATE( mask_loc( grid%nx, grid%ny))
    mask_loc = mask

    ! Initialise poly_mult and poly
    ALLOCATE( poly_mult( grid%nx*grid%ny,2))
    ALLOCATE( poly(      grid%nx*grid%ny,2))
    n_tot = 0

    ! Calculate polygons for all TRUE regions of the mask
    DO i = 1, grid%nx
    DO j = 1, grid%ny

      IF (mask_loc( i,j)) THEN
        ! Found a seed for a TRUE region

        ! Calculate a polygon enveloping this TRUE region, and
        ! remove the region from the mask
        CALL calc_grid_mask_as_polygon( grid, mask_loc, i, j, poly, n_poly)

        ! Add this polygon to poly_mult
        poly_mult( n_tot+1,1) = REAL( n_poly,dp)
        poly_mult( n_tot+1,2) = 0._dp
        poly_mult( n_tot+2:n_tot+1+n_poly,:) = poly( 1:n_poly,:)
        n_tot = n_tot + 1 + n_poly

      END IF ! IF (mask_loc( i,j)) THEN

    END DO
    END DO

    ! Crop memory
    CALL reallocate( poly_mult, n_tot, 2)

    ! Clean up after yourself
    DEALLOCATE( mask_loc)
    DEALLOCATE( poly)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grid_mask_as_polygons

  SUBROUTINE calc_grid_mask_as_polygon( grid, mask, i0, j0, poly, n_poly)
    ! Calculate a polygon enveloping the set of TRUE-valued mask cells around i0,j0,
    ! and remove that set of grid cells from the mask

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                         INTENT(IN)    :: grid
    LOGICAL,  DIMENSION(:,:  ),              INTENT(INOUT) :: mask
    INTEGER,                                 INTENT(IN)    :: i0,j0
    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: poly
    INTEGER,                                 INTENT(OUT)   :: n_poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_grid_mask_as_polygon'
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE                :: mask_ext
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE                :: map
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE                :: stack
    INTEGER                                                :: stackN
    INTEGER                                                :: i,j,ii,jj,i_sw,j_sw
    CHARACTER(LEN=256)                                     :: dir, dir_prev
    INTEGER                                                :: it

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( mask,1) /= grid%nx .OR. SIZE( mask,2) /= grid%ny) CALL crash('incorrect data dimensions!')
    IF (.NOT. mask( i0,j0)) CALL crash('seed at i0,j0 is not TRUE!')

    ! Pad a row of FALSEs around the mask so the set of TRUEs around i0,j0 is really an island
    ALLOCATE( mask_ext( 0:grid%nx+1, 0:grid%ny+1), source = .FALSE.)
    mask_ext( 1:grid%nx, 1:grid%ny) = mask

    ! Use a flood-fill algorithm to find the map of same-valued grid cells around mask cell i0,j0
    ALLOCATE( map(   0:grid%nx+1, 0:grid%ny+1) , source = 0)
    ALLOCATE( stack( (grid%nx+2)*(grid%ny+2),2))

    map( i0,j0) = 1
    stackN = 1
    stack( 1,:) = [i0,j0]

    i_sw = 0
    j_sw = 0

    DO WHILE (stackN > 0)

      ! Take the last element from the stack
      i = stack( stackN,1)
      j = stack( stackN,2)
      stackN = stackN - 1

      ! Mark it as mapped
      map( i,j) = 2

      ! Remove it from the input mask
      mask( i,j) = .FALSE.

      ! Add all non-mapped, non-stacked, TRUE-valued neighbours to the stack
      DO ii = MAX( 0, i-1), MIN( grid%nx+1, i+1)
      DO jj = MAX( 0, j-1), MIN( grid%ny+1, j+1)
        IF (ii /= i .AND. jj /= j) CYCLE ! Don't include diagonal neighbours
        IF (map( ii,jj) == 0 .AND. mask_ext( ii,jj)) THEN
          ! Add this neighbour to the stack
          stackN = stackN + 1
          stack( stackN,:) = [ii,jj]
          ! Mark this neighbour on the map as stacked
          map( ii,jj) = 1
        END IF
      END DO
      END DO

      ! If it is a southwest corner, save it as a starting point for the outline tracer
      IF (.NOT. mask_ext( i-1,j) .AND. .NOT. mask_ext( i,j-1)) THEN
        i_sw = i
        j_sw = j
      END IF

    END DO ! DO WHILE (stackN > 0)
    ! Safety
    IF (i_sw == 0 .OR. j_sw == 0) CALL crash('couldnt find starting SW corner!')

    ! Start at the southwest corner we found earlier
    ii = i_sw
    jj = j_sw
    n_poly = 2
    poly( 1,:) = [grid%x( ii) - grid%dx/2._dp, grid%y( jj) - grid%dx/2._dp]
    poly( 2,:) = [grid%x( ii) + grid%dx/2._dp, grid%y( jj) - grid%dx/2._dp]
    dir = 'east'

    ! Trace the outline of the mapped grid cells.
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > grid%nx*grid%ny) CALL crash('outline tracer got stuck!')

      ! Cycle direction
      dir_prev = dir

      ! Check which way we go next
      IF     (dir_prev == 'east') THEN
        ! We were going east, so [ii,jj] is to the north of us, and we can
        ! continue north, east, or south

        IF     (map( ii+1,jj) == 0) THEN
          ! Go north
          dir = 'north'
        ELSEIF (map( ii+1,jj-1) == 0) THEN
          ! Go east
          dir = 'east'
          ii = ii+1
        ELSEIF (map( ii,jj-1) == 0) THEN
          ! Go south
          dir = 'south'
          ii = ii+1
          jj = jj-1
        ELSE
          CALL crash('outline tracer got stuck while going east!')
        END IF

      ELSEIF (dir_prev == 'north') THEN
        ! We were going north, so [ii,jj] is to the west of us, and we can
        ! continue west, north, or east

        IF     (map( ii,jj+1) == 0) THEN
          ! Go west
          dir = 'west'
        ELSEIF (map( ii+1,jj+1) == 0) THEN
          ! Go north
          dir = 'north'
          jj = jj+1
        ELSEIF (map( ii+1,jj) == 0) THEN
          ! Go east
          dir = 'east'
          ii = ii+1
          jj = jj+1
        ELSE
          CALL crash('outline tracer got stuck while going north!')
        END IF

      ELSEIF (dir_prev == 'west') THEN
        ! We were going west, so [ii,jj] is to the south of us, and we can
        ! continue south, west, or north

        IF     (map( ii-1,jj) == 0) THEN
          ! Go south
          dir = 'south'
        ELSEIF (map( ii-1,jj+1) == 0) THEN
          ! Go west
          dir = 'west'
          ii = ii-1
        ELSEIF (map( ii,jj+1) == 0) THEN
          ! Go north
          dir = 'north'
          ii = ii-1
          jj = jj+1
        ELSE
          CALL crash('outline tracer got stuck while going north!')
        END IF

      ELSEIF (dir_prev == 'south') THEN
        ! We were going south, so [ii,jj] is to the east of us, and we can
        ! continue east, south, or west

        IF     (map( ii,jj-1) == 0) THEN
          ! Go east
          dir = 'east'
        ELSEIF (map( ii-1,jj-1) == 0) THEN
          ! Go south
          dir = 'south'
          jj = jj-1
        ELSEIF (map( ii-1,jj) == 0) THEN
          ! Go west
          dir = 'west'
          ii = ii-1
          jj = jj-1
        ELSE
          CALL crash('outline tracer got stuck while going north!')
        END IF

      ELSE
        CALL crash('unknown dir_prev "' // TRIM( dir_prev) // '"!')
      END IF

      ! Add new vertex to the polygon
      n_poly = n_poly+1
      IF     (dir == 'east') THEN
        poly( n_poly,:) = poly( n_poly-1,:) + [grid%dx, 0._dp]
      ELSEIF (dir == 'north') THEN
        poly( n_poly,:) = poly( n_poly-1,:) + [0._dp, grid%dx]
      ELSEIF (dir == 'west') THEN
        poly( n_poly,:) = poly( n_poly-1,:) + [-grid%dx, 0._dp]
      ELSEIF (dir == 'south') THEN
        poly( n_poly,:) = poly( n_poly-1,:) + [0._dp, -grid%dx]
      ELSE
        CALL crash('unknown dir "' // TRIM( dir) // '"!')
      END IF

      ! If we've reached the starting point again, we're done
      IF (NORM2( poly( n_poly,:) - poly( 1,:)) < grid%dx / 10._dp) THEN
        ! Don't double count
        n_poly = n_poly-1
        EXIT
      END IF

    END DO ! DO WHILE (.TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grid_mask_as_polygon

! == Subroutines for manipulating gridded data in distributed memory

  SUBROUTINE distribute_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 2-D gridded data field from the Master.
    ! Input from Master: total data field in field form
    ! Output to all: partial data in vector form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:    ),                        INTENT(IN)    :: d_grid

    ! Output variables:
    REAL(dp), DIMENSION(:      ),                        INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_gridded_data_from_master_dp_2D'
    INTEGER                                                            :: n,i,j
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_total

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Convert gridded data to vector form
    IF (par%master) THEN

      ! Allocate memory
      ALLOCATE( d_grid_vec_total( grid%n), source = 0._dp)

      ! Convert to vector form
      DO n = 1, grid%n
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        d_grid_vec_total( n) = d_grid( i,j)
      END DO

    END IF ! IF (par%master) THEN

    ! Distribute vector-form data to the processes
    CALL distribute_from_master_dp_1D( d_grid_vec_total, d_grid_vec_partial)

    ! Clean up after yourself
    IF (par%master) DEALLOCATE( d_grid_vec_total)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_gridded_data_from_master_dp_2D

  SUBROUTINE distribute_gridded_data_from_master_dp_3D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 3-D gridded data field from the Master.
    ! Input from Master: total data field in field form
    ! Output to all: partial data in vector form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_grid

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ),                        INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_gridded_data_from_master_dp_3D'
    INTEGER                                                            :: k
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: d_grid_2D
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_partial_2D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (par%master .AND. SIZE( d_grid,3) /= SIZE( d_grid_vec_partial,2)) CALL crash('vector sizes dont match!')

    ! Allocate memory
    IF (par%master) ALLOCATE( d_grid_2D( SIZE( d_grid,1), SIZE( d_grid,2)), source = 0._dp)
    ALLOCATE( d_grid_vec_partial_2D( SIZE( d_grid_vec_partial,1)), source = 0._dp)

    ! Treat each layer as a separate 2-D field
    DO k = 1, SIZE( d_grid_vec_partial,2)
      IF (par%master) d_grid_2D = d_grid( :,:,k)
      CALL distribute_gridded_data_from_master_dp_2D( grid, d_grid_2D, d_grid_vec_partial_2D)
      d_grid_vec_partial( :,k) = d_grid_vec_partial_2D
    END DO

    ! Clean up after yourself
    IF (par%master) DEALLOCATE( d_grid_2D)
    DEALLOCATE( d_grid_vec_partial_2D)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_gridded_data_from_master_dp_3D

  SUBROUTINE gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial, d_grid)
    ! Gather a 2-D gridded data field to the Master.
    ! Input from all: partial data in vector form
    ! Output to Master: total data field in field form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:      ),                        INTENT(IN)    :: d_grid_vec_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ),                        INTENT(OUT)   :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_gridded_data_to_master_dp_2D'
    INTEGER                                                            :: n,i,j
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_total

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    IF (par%master) ALLOCATE( d_grid_vec_total( grid%n), source = 0._dp)

    ! Gather data
    CALL gather_to_master_dp_1D( d_grid_vec_partial, d_grid_vec_total)

    ! Convert to grid form
    IF (par%master) THEN
      DO n = 1, grid%n
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        d_grid( i,j) = d_grid_vec_total( n)
      END DO
    END IF ! IF (par%master) THEN

    ! Clean up after yourself
    IF (par%master) DEALLOCATE( d_grid_vec_total)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_gridded_data_to_master_dp_2D

  SUBROUTINE gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)
    ! Gather a 3-D gridded data field to the Master.
    ! Input from all: partial data in vector form
    ! Output to Master: total data field in field form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:    ),                        INTENT(IN)    :: d_grid_vec_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(OUT)   :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_gridded_data_to_master_dp_3D'
    INTEGER                                                            :: k
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: d_grid_2D
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_partial_2D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (par%master .AND. SIZE( d_grid,3) /= SIZE( d_grid_vec_partial,2)) CALL crash('vector sizes dont match!')

    ! Allocate memory
    IF (par%master) ALLOCATE( d_grid_2D( grid%nx, grid%ny), source = 0._dp)
    ALLOCATE( d_grid_vec_partial_2D( grid%n_loc), source = 0._dp)

    ! Treat each layer as a separate 2-D field
    DO k = 1, SIZE( d_grid_vec_partial,2)
      d_grid_vec_partial_2D = d_grid_vec_partial( :,k)
      CALL gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial_2D, d_grid_2D)
      IF (par%master) d_grid( :,:,k) = d_grid_2D
    END DO

    ! Clean up after yourself
    IF (par%master) DEALLOCATE( d_grid_2D)
    DEALLOCATE( d_grid_vec_partial_2D)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_gridded_data_to_master_dp_3D

END MODULE grid_basic
