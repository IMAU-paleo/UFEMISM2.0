MODULE mesh_zeta

  ! Routines used in performing the vertical coordinate transformation z -> zeta

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE reallocate_mod                                         , ONLY: reallocate
  USE mesh_types                                             , ONLY: type_mesh
  USE CSR_matrix_basics                            , ONLY: allocate_matrix_CSR_loc, add_entry_CSR_dist
  use shape_functions, only: calc_shape_functions_1D_reg_2nd_order

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  ! Initialising the scaled vertical coordinate zeta
  ! ================================================

  subroutine initialise_scaled_vertical_coordinate( mesh)
    ! Initialise the scaled vertical coordinate zeta

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_scaled_vertical_coordinate'

    ! Add routine to path
    call init_routine( routine_name)

    allocate( mesh%zeta(      mesh%nz  ))
    allocate( mesh%zeta_stag( mesh%nz-1))

    ! Calculate zeta values
    select case (mesh%choice_zeta_grid)
    case default
      call crash('unknown choice_zeta_grid "' // trim( mesh%choice_zeta_grid) // '"')
    case ('regular')
      call initialise_scaled_vertical_coordinate_regular( mesh)
    case ('irregular_log')
      call initialise_scaled_vertical_coordinate_irregular_log( mesh)
    case ('old_15_layer_zeta')
      call initialise_scaled_vertical_coordinate_old_15_layer( mesh)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_scaled_vertical_coordinate

  SUBROUTINE initialise_scaled_vertical_coordinate_regular( mesh)
    ! Initialise the scaled vertical coordinate zeta
    ! Regular zeta grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_scaled_vertical_coordinate_regular'
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Fill zeta values
    DO k = 1, mesh%nz
      mesh%zeta( k) = REAL( k-1,dp) / REAL( mesh%nz-1,dp)
    END DO

    ! Calculate zeta_stag
    mesh%zeta_stag = (mesh%zeta( 1:mesh%nz-1) + mesh%zeta( 2:mesh%nz)) / 2._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_scaled_vertical_coordinate_regular

  SUBROUTINE initialise_scaled_vertical_coordinate_irregular_log( mesh)
    ! Initialise the scaled vertical coordinate zeta
    !
    ! This scheme ensures that the ratio between subsequent grid spacings is
    ! constant, and that the ratio between the first (surface) and last (basal)
    ! layer thickness is (approximately) equal to R

    ! In/output variables:
    type(type_mesh),    intent(inout) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_scaled_vertical_coordinate_irregular_log'
    INTEGER                                            :: k
    REAL(dp)                                           :: sigma, sigma_stag

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (mesh%zeta_irregular_log_R <= 0._dp) CALL crash('R should be positive!')

    ! Exception: R = 1 implies a regular grid, but the equation becomes 0/0
    IF (mesh%zeta_irregular_log_R == 1._dp) THEN
      CALL initialise_scaled_vertical_coordinate_regular( mesh)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    DO k = 1, mesh%nz
      ! Regular grid
      sigma = REAL( k-1,dp) / REAL( mesh%nz-1,dp)
      mesh%zeta( mesh%nz + 1 - k) = 1._dp - (mesh%zeta_irregular_log_R**sigma - 1._dp) / (mesh%zeta_irregular_log_R - 1._dp)
      ! Staggered grid
      IF (k < mesh%nz) THEN
        sigma_stag = sigma + 0.5_dp / REAL( mesh%nz-1,dp)
        mesh%zeta_stag( mesh%nz - k) = 1._dp - (mesh%zeta_irregular_log_R**sigma_stag - 1._dp) / (mesh%zeta_irregular_log_R - 1._dp)
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_scaled_vertical_coordinate_irregular_log

  SUBROUTINE initialise_scaled_vertical_coordinate_old_15_layer( mesh)
    ! Set up the vertical zeta grid
    !
    ! Use the old irregularly-spaced 15 layers from ANICE

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_scaled_vertical_coordinate_old_15_layer'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (mesh%nz /= 15) CALL crash('only works when nz = 15!')

    mesh%zeta = (/ 0.00_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp, 0.50_dp, 0.60_dp, 0.70_dp, 0.80_dp, 0.90_dp, 0.925_dp, 0.95_dp, 0.975_dp, 0.99_dp, 1.00_dp /)

    ! Calculate zeta_stag
    mesh%zeta_stag = (mesh%zeta( 1:mesh%nz-1) + mesh%zeta( 2:mesh%nz)) / 2._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_scaled_vertical_coordinate_old_15_layer

  ! Some mathematical operations on a function of zeta
  ! ==================================================

  PURE FUNCTION integrate_from_zeta_is_one_to_zeta_is_zetap( zeta, f) RESULT( integral_f)
    ! This subroutine integrates f from zeta( k=nz) = 1 (which corresponds to the ice base)
    ! to the level zetap = zeta( k) for all values of k
    !
    ! NOTE: if the integrand f is positive, the integral is negative because the integration is in
    ! the opposite zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the bottom up to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f( k+1) and f( k):
    !  integral_f( k) = integral_f( k+1) + 0.5 * (f( k+1) + f( k)) * (-dzeta)
    ! with dzeta = zeta( k+1) - zeta( k). So for f > 0,  integral_f < 0.
    !
    ! Heiko Goelzer (h.goelzer@uu.nl) Jan 2016

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp), DIMENSION( SIZE( zeta,1))                :: integral_f

    ! Local variables:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    integral_f( nz) = 0._dp

    DO k = nz-1, 1, -1
      integral_f( k) = integral_f( k+1) - 0.5_dp * (f( k+1) + f( k)) * (zeta( k+1) - zeta( k))
    END DO

  END FUNCTION integrate_from_zeta_is_one_to_zeta_is_zetap

  PURE FUNCTION integrate_from_zeta_is_zero_to_zeta_is_zetap( zeta, f) RESULT( integral_f)
    ! This subroutine integrates f from zeta( k=1) = 0 (which corresponds to the ice surface)
    ! to the level zetap = zeta( k) for all values of k
    !
    ! If the integrand f is positive, the integral is positive because the integration is in
    ! the zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the top down to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f( k) and f( k-1):
    ! integral_f( k) = integral_f( k-1) + 0.5 * (f( k) + f( k-1)) * (dzeta)
    ! with dzeta = zeta( k+1) - zeta( k).
    !
    ! Heiko Goelzer (h.goelzer@uu.nl) Jan 2016

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp), DIMENSION( SIZE( zeta,1))                :: integral_f

    ! Local variables:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    integral_f( 1) = 0._dp
    DO k = 2, nz, 1
      integral_f( k) = integral_f( k-1) + 0.5_dp * (f( k) + f( k-1)) * (zeta( k) - zeta( k-1))
    END DO

  END FUNCTION integrate_from_zeta_is_zero_to_zeta_is_zetap

  PURE FUNCTION integrate_over_zeta( zeta, f) RESULT( integral_f)
    ! Integrate f over zeta from zeta = 1 (the ice base) to zeta = 0 (the ice surface)
    !
    ! NOTE: defined so that if f is positive, then integral_f is positive too.

    IMPLICIT NONE

    ! In/output variable:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp)                                           :: integral_f

    ! Local variable:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    ! Initial value is zero
    integral_f = 0._dp

    ! Intermediate values include sum of all previous values
    ! Take current value as average between points

    DO k = 2, nz
       integral_f = integral_f + 0.5_dp * (f( k) + f( k-1)) * (zeta( k) - zeta( k-1))
    END DO

  END FUNCTION integrate_over_zeta

  PURE FUNCTION vertical_average( zeta, f) RESULT( average_f)
    ! Calculate the vertical average of any given function f defined at the vertical zeta grid.
    !
    ! The integration is in the direction of the positive zeta-axis from zeta( k=1) = 0 up to zeta( k=nz) = 1.
    ! Numerically: the average between layer k and k+1 is calculated and multiplied by the distance between those
    ! layers k and k+1, which is imediately the weight factor for this contribution because de total layer distance
    ! is scaled to 1. The sum of all the weighted contribution gives average_f the vertical average of f.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp)                                           :: average_f

    ! Local variables:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    average_f = 0._dp

    DO k = 1, nz-1
      average_f = average_f + 0.5_dp * (f( k+1) + f( k)) * (zeta( k+1) - zeta( k))
    END DO

  END FUNCTION vertical_average

  ! Gradient operators in the zeta dimension
  ! ========================================

  SUBROUTINE calc_vertical_operators_reg_1D( mesh)
    ! Calculate mapping and d/dzeta operators in the 1-D vertical column
    !
    ! NOTE: since its well possible that the number of cores running the model
    !       exceeds the number of vertical layers, these matrices are stored in
    !       process-local memory

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                           INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_vertical_operators_reg_1D'
    INTEGER                                                      :: ncols, nrows, nnz_per_row_est, nnz_est
    INTEGER                                                      :: k
    REAL(dp)                                                     :: z
    INTEGER                                                      :: n_c
    INTEGER,  DIMENSION(2)                                       :: i_c
    REAL(dp), DIMENSION(2)                                       :: z_c
    INTEGER                                                      :: i,kn
    REAL(dp)                                                     :: Nfz_i, Nfzz_i
    REAL(dp), DIMENSION(2)                                       :: Nfz_c, Nfzz_c

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nz      ! from
    nrows           = mesh%nz      ! to
    nnz_per_row_est = 3
    nnz_est         = nrows * nnz_per_row_est

    CALL allocate_matrix_CSR_loc( mesh%M_ddzeta_k_k_1D  , nrows, ncols, nnz_est)
    CALL allocate_matrix_CSR_loc( mesh%M_d2dzeta2_k_k_1D, nrows, ncols, nnz_est)

    DO k = 1, mesh%nz

      ! Source: layer k
      z = mesh%zeta( k)

      ! Neighbours: layers k-1 and k+1
      IF     (k == 1) THEN
        n_c = 2
        i_c = [2, 3]
        z_c = [mesh%zeta( 2), mesh%zeta( 3)]
      ELSEIF (k == mesh%nz) THEN
        n_c = 2
        i_c = [mesh%nz-2, mesh%nz-1]
        z_c = [mesh%zeta( mesh%nz-2), mesh%zeta( mesh%nz-1)]
      ELSE
        n_c = 2
        i_c = [k-1, k+1]
        z_c = [mesh%zeta( k-1), mesh%zeta( k+1)]
      END IF

      ! Calculate shape functions
      CALL calc_shape_functions_1D_reg_2nd_order( z, n_c, n_c, z_c, Nfz_i, Nfzz_i, Nfz_c, Nfzz_c)

      ! Diagonal element: shape function for the source point
      CALL add_entry_CSR_dist( mesh%M_ddzeta_k_k_1D  , k, k, Nfz_i )
      CALL add_entry_CSR_dist( mesh%M_d2dzeta2_k_k_1D, k, k, Nfzz_i)

      ! Off-diagonal elements: shape functions for neighbours
      DO i = 1, n_c
        kn = i_c( i)
        CALL add_entry_CSR_dist( mesh%M_ddzeta_k_k_1D  , k, kn, Nfz_c(  i))
        CALL add_entry_CSR_dist( mesh%M_d2dzeta2_k_k_1D, k, kn, Nfzz_c( i))
      END DO

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_operators_reg_1D

  SUBROUTINE calc_vertical_operators_stag_1D( mesh)
    ! Calculate mapping and d/dzeta operators in the 1-D vertical column
    !
    ! NOTE: since its well possible that the number of cores running the model
    !       exceeds the number of vertical layers, these matrices are stored in
    !       process-local memory

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                           INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_vertical_operators_stag_1D'
    INTEGER                                                      :: ncols, nrows, nnz_per_row_est, nnz_est
    INTEGER                                                      :: ks
    REAL(dp)                                                     :: dzeta
    INTEGER                                                      :: k, row, col
    INTEGER                                                      ::  k_lo,  k_hi,  ks_lo,  ks_hi
    REAL(dp)                                                     :: wk_lo, wk_hi, wks_lo, wks_hi

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == k (regular) to ks (staggered)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nz      ! from
    nrows           = mesh%nz-1    ! to
    nnz_per_row_est = 2
    nnz_est         = nrows * nnz_per_row_est

    CALL allocate_matrix_CSR_loc( mesh%M_map_k_ks_1D   , nrows, ncols, nnz_est)
    CALL allocate_matrix_CSR_loc( mesh%M_ddzeta_k_ks_1D, nrows, ncols, nnz_est)

    DO ks = 1, mesh%nz-1

      ! Indices of neighbouring grid points
      k_lo = ks
      k_hi = ks + 1

      ! Local grid spacing
      dzeta = mesh%zeta( k_hi) - mesh%zeta( k_lo)

      ! Linear interpolation weights
      wk_lo = (mesh%zeta( k_hi) - mesh%zeta_stag( ks)) / dzeta
      wk_hi = 1._dp - wk_lo

      ! Left-hand neighbour
      row = ks
      col = k_lo
      CALL add_entry_CSR_dist( mesh%M_map_k_ks_1D   , row, col, wk_lo         )
      CALL add_entry_CSR_dist( mesh%M_ddzeta_k_ks_1D, row, col, -1._dp / dzeta)

      ! Right-hand neighbour
      row = ks
      col = k_hi
      CALL add_entry_CSR_dist( mesh%M_map_k_ks_1D   , row, col, wk_hi         )
      CALL add_entry_CSR_dist( mesh%M_ddzeta_k_ks_1D, row, col,  1._dp / dzeta)

    END DO

  ! == ks (staggered) to k (regular)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nz-1    ! from
    nrows           = mesh%nz      ! to
    nnz_per_row_est = 2
    nnz_est         = nrows * nnz_per_row_est

    CALL allocate_matrix_CSR_loc( mesh%M_map_ks_k_1D   , nrows, ncols, nnz_est)
    CALL allocate_matrix_CSR_loc( mesh%M_ddzeta_ks_k_1D, nrows, ncols, nnz_est)

    DO k = 1, mesh%nz

      ! Indices of neighbouring grid points
      IF     (k == 1) THEN
        ks_lo = 1
        ks_hi = 2
      ELSEIF (k == mesh%nz) THEN
        ks_lo = mesh%nz - 2
        ks_hi = mesh%nz - 1
      ELSE
        ks_lo = k - 1
        ks_hi = k
      END IF

      ! Local grid spacing
      dzeta = mesh%zeta_stag( ks_hi) - mesh%zeta_stag( ks_lo)

      ! Linear interpolation weights
      wks_lo = (mesh%zeta_stag( ks_hi) - mesh%zeta( k)) / dzeta
      wks_hi = 1._dp - wks_lo

      ! Left-hand neighbour
      row = k
      col = ks_lo
      CALL add_entry_CSR_dist( mesh%M_map_ks_k_1D   , row, col, wks_lo        )
      CALL add_entry_CSR_dist( mesh%M_ddzeta_ks_k_1D, row, col, -1._dp / dzeta)

      ! Right-hand neighbour
      row = k
      col = ks_hi
      CALL add_entry_CSR_dist( mesh%M_map_ks_k_1D   , row, col, wks_hi        )
      CALL add_entry_CSR_dist( mesh%M_ddzeta_ks_k_1D, row, col,  1._dp / dzeta)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_operators_stag_1D

  SUBROUTINE calc_zeta_operators_tridiagonal( mesh)
    ! Calculate zeta operators in tridiagonal form for efficient use in thermodynamics

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_zeta_operators_tridiagonal'
    INTEGER                                                      :: k,i,ii1,ii2,ii,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for the three diagonals of both matrices
    ALLOCATE( mesh%M_ddzeta_k_k_ldiag(   mesh%nz-1))
    ALLOCATE( mesh%M_ddzeta_k_k_diag(    mesh%nz  ))
    ALLOCATE( mesh%M_ddzeta_k_k_udiag(   mesh%nz-1))
    ALLOCATE( mesh%M_d2dzeta2_k_k_ldiag( mesh%nz-1))
    ALLOCATE( mesh%M_d2dzeta2_k_k_diag(  mesh%nz  ))
    ALLOCATE( mesh%M_d2dzeta2_k_k_udiag( mesh%nz-1))

    ! Fill in coefficients

    ! d/dzeta
    DO k = 1, mesh%nz

      ! Lower diagonal
      IF (k > 1) THEN
        ! Initialise
        mesh%M_ddzeta_k_k_ldiag( k-1) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_ddzeta_k_k_1D%ptr( i)
        ii2 = mesh%M_ddzeta_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_ddzeta_k_k_1D%ind( ii)
          IF (j == k-1) mesh%M_ddzeta_k_k_ldiag( k-1) = mesh%M_ddzeta_k_k_1D%val( ii)
        END DO
      END IF

      ! Central diagonal
      ! Initialise
      mesh%M_ddzeta_k_k_diag( k) = 0._dp
      ! Find matrix element at [k,k-1]
      i = k
      ii1 = mesh%M_ddzeta_k_k_1D%ptr( i)
      ii2 = mesh%M_ddzeta_k_k_1D%ptr( i+1)-1
      DO ii = ii1, ii2
        j = mesh%M_ddzeta_k_k_1D%ind( ii)
        IF (j == k) mesh%M_ddzeta_k_k_diag( k) = mesh%M_ddzeta_k_k_1D%val( ii)
      END DO

      ! Upper diagonal
      IF (k < mesh%nz) THEN
        ! Initialise
        mesh%M_ddzeta_k_k_udiag( k) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_ddzeta_k_k_1D%ptr( i)
        ii2 = mesh%M_ddzeta_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_ddzeta_k_k_1D%ind( ii)
          IF (j == k+1) mesh%M_ddzeta_k_k_udiag( k) = mesh%M_ddzeta_k_k_1D%val( ii)
        END DO
      END IF

    END DO

    ! d2/dzeta2
    DO k = 1, mesh%nz

      ! Lower diagonal
      IF (k > 1) THEN
        ! Initialise
        mesh%M_d2dzeta2_k_k_ldiag( k-1) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_d2dzeta2_k_k_1D%ptr( i)
        ii2 = mesh%M_d2dzeta2_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_d2dzeta2_k_k_1D%ind( ii)
          IF (j == k-1) mesh%M_d2dzeta2_k_k_ldiag( k-1) = mesh%M_d2dzeta2_k_k_1D%val( ii)
        END DO
      END IF

      ! Central diagonal
      ! Initialise
      mesh%M_d2dzeta2_k_k_diag( k) = 0._dp
      ! Find matrix element at [k,k-1]
      i = k
      ii1 = mesh%M_d2dzeta2_k_k_1D%ptr( i)
      ii2 = mesh%M_d2dzeta2_k_k_1D%ptr( i+1)-1
      DO ii = ii1, ii2
        j = mesh%M_d2dzeta2_k_k_1D%ind( ii)
        IF (j == k) mesh%M_d2dzeta2_k_k_diag( k) = mesh%M_d2dzeta2_k_k_1D%val( ii)
      END DO

      ! Upper diagonal
      IF (k < mesh%nz) THEN
        ! Initialise
        mesh%M_d2dzeta2_k_k_udiag( k) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_d2dzeta2_k_k_1D%ptr( i)
        ii2 = mesh%M_d2dzeta2_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_d2dzeta2_k_k_1D%ind( ii)
          IF (j == k+1) mesh%M_d2dzeta2_k_k_udiag( k) = mesh%M_d2dzeta2_k_k_1D%val( ii)
        END DO
      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_zeta_operators_tridiagonal

END MODULE mesh_zeta