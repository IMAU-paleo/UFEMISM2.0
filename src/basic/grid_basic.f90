MODULE grid_basic

  ! Data types and subroutines for working with simple square grids

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE parameters
  USE petsc_basic                                            , ONLY: perr, mat_CSR2petsc
  USE reallocate_mod                                         , ONLY: reallocate
  USE math_utilities                                         , ONLY: linint_points
  USE mpi_distributed_memory                                 , ONLY: partition_list, distribute_from_master_dp_1D
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

    ! Parallelisation
    INTEGER                                 :: n1,n2,n_loc                   ! Matrix rows owned by each process

  END TYPE type_grid

  TYPE type_grid_lonlat
    ! A lon/lat grid

    ! Basic properties
    CHARACTER(LEN=256)                      :: name                          !           A nice name tag
    INTEGER                                 :: nlon                          !           Number of grid cells in the longitude direction
    INTEGER                                 :: nlat                          !           Number of grid cells in the latitude direction
    INTEGER                                 :: n                             !           Total number of grid cells (= nx * ny)
    REAL(dp)                                :: dlon                          ! [degrees] Resolution in the longitude direction
    REAL(dp)                                :: dlat                          ! [degrees] Resolution in the latitude direction
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lon                           ! [degrees east ] Longitude of each grid point
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lat                           ! [degrees north] Latitude  of each grid point
    REAL(dp)                                :: lonmin                        ! [degrees east ] Lon/lat range covered by the grid
    REAL(dp)                                :: lonmax                        ! [degrees east ]
    REAL(dp)                                :: latmin                        ! [degrees north]
    REAL(dp)                                :: latmax                        ! [degrees north]

    ! Remapping data
    REAL(dp)                                :: tol_dist                      ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: ij2n, n2ij                    !           Conversion tables for grid-form vs. vector-form data

    ! Parallelisation
    INTEGER                                 :: n1,n2,n_loc                   ! Matrix rows owned by each process

  END TYPE type_grid_lonlat

CONTAINS

! ===== Subroutines ======
! ========================

! == Calculate contour lines for mesh generation from gridded/meshed data

  SUBROUTINE calc_grid_contour_as_line( grid, d, f, line)
    ! Calculate a contour line at level f for data d on a square grid%
    ! Generate the contour line in UFEMISM line format (i.e. unordered
    ! individual line segments).

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d
    REAL(dp),                                INTENT(IN)    :: f
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: line

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'mesh_add_smileyface'
    REAL(dp), PARAMETER                                    :: tol = 1E-5_dp
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_scaled
    REAL(dp)                                               :: d_scaled_min, d_scaled_max
    INTEGER                                                :: i,j
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

    ! Scale d so that it runs from -1 to 1, with 0 representing f
    d_scaled = d - f
    d_scaled_min = MINVAL( d_scaled)
    d_scaled_max = MAXVAL( d_scaled)
    DO i = 1, grid%nx
      DO j = 1, grid%ny
        IF (d_scaled( i,j) > 0._dp) THEN
          d_scaled( i,j) = d_scaled( i,j) /           d_scaled_max
        ELSE
          d_scaled( i,j) = d_scaled( i,j) / (-1._dp * d_scaled_min)
        END IF
      END DO
    END DO

    n_max = 1000
    ALLOCATE( line( n_max, 4))

    n = 0
    DO i = 1, grid%nx-1
      DO j = 1, grid%ny-1

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

! == Calculate d/dx, d/dy operators on a square grid (used in grid-to-mesh remapping)

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

    END DO ! DO n = grid%n1, grid%n2

    ! Convert to PETSc format
    CALL mat_CSR2petsc( M_ddx_CSR, M_ddx)
    CALL mat_CSR2petsc( M_ddy_CSR, M_ddy)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_ddx_CSR)
    CALL deallocate_matrix_CSR_dist( M_ddy_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_grid

! == Set up field-form -to- vector-form translation tables

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

! == Distribute gridded data from the Master

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

END MODULE grid_basic
