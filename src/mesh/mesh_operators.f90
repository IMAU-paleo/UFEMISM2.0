MODULE mesh_operators

  ! Routines for calculating matrix operators on the mesh.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE math_utilities                                         , ONLY: calc_matrix_inverse_2_by_2, calc_matrix_inverse_3_by_3, calc_matrix_inverse_5_by_5
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist
  USE mesh_utilities                                         , ONLY: extend_group_single_iteration_a, extend_group_single_iteration_b, &
                                                                     extend_group_single_iteration_c
  USE petsc_basic                                            , ONLY: perr, multiply_CSR_matrix_with_vector_1D, multiply_CSR_matrix_with_vector_2D

  IMPLICIT NONE

CONTAINS

! ===== Subroutines for applying mesh operators =====
! ===================================================

  ! NOTE: not all possible combinations are written, as only a few
  !       are ever needed.

  ! Gradients between a-grid and a-grid
  SUBROUTINE ddx_a_a_2D( mesh, d_a, ddx_a)
    ! ddx a 2-D data field from the a-grid (vertices) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddx_a_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddx_a,1) /= mesh%nV_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddxping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_a, d_a, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_a_2D

  SUBROUTINE ddx_a_a_3D( mesh, d_a, ddx_a)
    ! ddx a 3-D data field from the a-grid (vertices) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddx_a_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddx_a,1) /= mesh%nV_loc .OR. SIZE( d_a,2) /= SIZE( ddx_a,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddxping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_ddx_a_a, d_a, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_a_3D

  SUBROUTINE ddy_a_a_2D( mesh, d_a, ddy_a)
    ! ddy a 2-D data field from the a-grid (vertices) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddy_a_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddy_a,1) /= mesh%nV_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddyping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_a, d_a, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_a_2D

  SUBROUTINE ddy_a_a_3D( mesh, d_a, ddy_a)
    ! ddy a 3-D data field from the a-grid (vertices) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddy_a_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddy_a,1) /= mesh%nV_loc .OR. SIZE( d_a,2) /= SIZE( ddy_a,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddyping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_ddy_a_a, d_a, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_a_3D

  ! Mapping between a-grid and b-grid
  SUBROUTINE map_a_b_2D( mesh, d_a, d_b)
    ! Map a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'map_a_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_a_b, d_a, d_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_a_b_2D

  SUBROUTINE map_a_b_3D( mesh, d_a, d_b)
    ! Map a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'map_a_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc .OR. SIZE( d_a,2) /= SIZE( d_b,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_map_a_b, d_a, d_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_a_b_3D

  SUBROUTINE map_b_a_2D( mesh, d_b, d_a)
    ! Map a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'map_b_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_map_b_a, d_b, d_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_b_a_2D

  SUBROUTINE map_b_a_3D( mesh, d_b, d_a)
    ! Map a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'map_b_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc .OR. SIZE( d_a,2) /= SIZE( d_b,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_map_b_a, d_b, d_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_b_a_3D

  ! Gradients between a-grid and b-grid
  SUBROUTINE ddx_a_b_2D( mesh, d_a, ddx_b)
    ! ddx a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddx_a_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddx_b,1) /= mesh%nTri_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddxping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_b_2D

  SUBROUTINE ddx_a_b_3D( mesh, d_a, ddx_b)
    ! ddx a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddx_a_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddx_b,1) /= mesh%nTri_loc .OR. SIZE( d_a,2) /= SIZE( ddx_b,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddxping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_b_3D

  SUBROUTINE ddx_b_a_2D( mesh, d_b, ddx_a)
    ! ddx a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddx_b_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( ddx_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddxping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_a, d_b, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_a_2D

  SUBROUTINE ddx_b_a_3D( mesh, d_b, ddx_a)
    ! ddx a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddx_b_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( ddx_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc .OR. SIZE( ddx_a,2) /= SIZE( d_b,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddxping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_ddx_b_a, d_b, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_a_3D

  SUBROUTINE ddy_a_b_2D( mesh, d_a, ddy_b)
    ! ddy a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddy_a_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddy_b,1) /= mesh%nTri_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddyping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_b_2D

  SUBROUTINE ddy_a_b_3D( mesh, d_a, ddy_b)
    ! ddy a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddy_a_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV_loc .OR. SIZE( ddy_b,1) /= mesh%nTri_loc .OR. SIZE( d_a,2) /= SIZE( ddy_b,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddyping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_b_3D

  SUBROUTINE ddy_b_a_2D( mesh, d_b, ddy_a)
    ! ddy a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)     :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)    :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddy_b_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( ddy_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddyping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_a, d_b, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_a_2D

  SUBROUTINE ddy_b_a_3D( mesh, d_b, ddy_a)
    ! ddy a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)     :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)     :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)    :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'ddy_b_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( ddy_a,1) /= mesh%nV_loc .OR. SIZE( d_b,1) /= mesh%nTri_loc .OR. SIZE( ddy_a,2) /= SIZE( d_b,2)) THEN
      CALL crash('vector and matrix sizes dont match!')
    END IF

    ! Perform the ddyping operation as a matrix multiplication
    CALL multiply_CSR_matrix_with_vector_2D( mesh%M_ddy_b_a, d_b, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_a_3D

! ===== Subroutines for calculating mesh operators =====
! ======================================================

  SUBROUTINE calc_all_matrix_operators_mesh( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between all the grids (a,b,c) on the mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_all_matrix_operators_mesh'
    INTEGER                                            :: nz

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Matrix operators
    nz = 4 ! DENK DROM - need to replace this with actual zeta thing!
    CALL calc_field_to_vector_form_translation_tables( mesh, nz)

    CALL calc_matrix_operators_mesh_a_a(               mesh)
    CALL calc_matrix_operators_mesh_a_b(               mesh)
    CALL calc_matrix_operators_mesh_a_c(               mesh)

    CALL calc_matrix_operators_mesh_b_a(               mesh)
    CALL calc_matrix_operators_mesh_b_b(               mesh)
    CALL calc_matrix_operators_mesh_b_c(               mesh)

    CALL calc_matrix_operators_mesh_c_a(               mesh)
    CALL calc_matrix_operators_mesh_c_b(               mesh)
    CALL calc_matrix_operators_mesh_c_c(               mesh)

    CALL calc_matrix_operators_mesh_b_b_2nd_order(     mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_all_matrix_operators_mesh

! == Calculate mapping and gradient operators between the a-, b-, and c-grids

  SUBROUTINE calc_matrix_operators_mesh_a_a( mesh)
    ! Calculate d/dx, and d/dy matrix operators between the a-grid (vertices) and the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_a_a'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: vi
    REAL(dp)                                           :: x, y
    INTEGER                                            :: vj
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nV)                       :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 2

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nV      ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_ddx_a_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_a_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_ddx_a_a%i1, mesh%M_ddx_a_a%i2

      ! The vertex represented by this matrix row
      vi = mesh%n2vi( row)
      x  = mesh%V( vi,1)
      y  = mesh%V( vi,2)

      ! Clean up previous map
      DO i = 1, stackN
        vj = stack( i)
        map( vj) = 0
      END DO

      ! Initialise the list of neighbours: just vi itself
      map( vi)  = 1
      stackN    = 1
      stack( 1) = vi

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_a( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        vj = stack( i)
        IF (vj == vi) CYCLE
        n_c = n_c + 1
        i_c( n_c) = vj
        x_c( n_c) = mesh%V( vj,1)
        y_c( n_c) = mesh%V( vj,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( mesh%M_ddx_a_a, row, row, Nfx_i)
      CALL add_entry_CSR_dist( mesh%M_ddy_a_a, row, row, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        vj = i_c( i)
        col = mesh%vi2n( vj)
        CALL add_entry_CSR_dist( mesh%M_ddx_a_a, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_a_a, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_a_a

  SUBROUTINE calc_matrix_operators_mesh_a_b( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_a_b'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, vi
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nV)                       :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nTri      ! to
    nrows_loc       = mesh%nTri_loc
    nnz_per_row_est = 3
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_map_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_map_a_b%i1, mesh%M_map_a_b%i2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Clean up previous map
      DO i = 1, stackN
        vi = stack( i)
        map( vi) = 0
      END DO

      ! Initialise the list of neighbours: the three vertices spanning ti
      stackN = 0
      DO n = 1, 3
        vi = mesh%Tri( ti,n)
        map( vi) = 1
        stackN = stackN + 1
        stack( stackN) = vi
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN < n_neighbours_min)
        CALL extend_group_single_iteration_a( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        vi = stack( i)
        n_c = n_c + 1
        i_c( n_c) = vi
        x_c( n_c) = mesh%V( vi,1)
        y_c( n_c) = mesh%V( vi,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        vi = i_c( i)
        col = mesh%vi2n( vi)
        CALL add_entry_CSR_dist( mesh%M_map_a_b, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( mesh%M_ddx_a_b, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_a_b, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_a_b

  SUBROUTINE calc_matrix_operators_mesh_a_c( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the c-grid (edges)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_a_c'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ei
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, vi
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nV)                       :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 4
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_map_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_map_a_c%i1, mesh%M_map_a_c%i2

      ! The edge represented by this matrix row
      ei = mesh%n2ei( row)
      x  = mesh%E( ei,1)
      y  = mesh%E( ei,2)

      ! Clean up previous map
      DO i = 1, stackN
        vi = stack( i)
        map( vi) = 0
      END DO

      ! Initialise the list of neighbours: the three or four vertices adjacent to edge ei
      stackN = 0
      DO n = 1, 4
        vi = mesh%EV( ei,n)
        IF (vi == 0) CYCLE
        map( vi) = 1
        stackN = stackN + 1
        stack( stackN) = vi
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN < n_neighbours_min)
        CALL extend_group_single_iteration_a( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        vi = stack( i)
        n_c = n_c + 1
        i_c( n_c) = vi
        x_c( n_c) = mesh%V( vi,1)
        y_c( n_c) = mesh%V( vi,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        vi = i_c( i)
        col = mesh%vi2n( vi)
        CALL add_entry_CSR_dist( mesh%M_map_a_c, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( mesh%M_ddx_a_c, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_a_c, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_a_c

  SUBROUTINE calc_matrix_operators_mesh_b_a( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_a'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: vi
    REAL(dp)                                           :: x, y
    INTEGER                                            :: iti, ti
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nTri)                     :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri    ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_map_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_map_b_a%i1, mesh%M_map_b_a%i2

      ! The vertex represented by this matrix row
      vi = mesh%n2vi( row)
      x  = mesh%V( vi,1)
      y  = mesh%V( vi,2)

      ! Clean up previous map
      DO i = 1, stackN
        ti = stack( i)
        map( ti) = 0
      END DO

      ! Initialise the list of neighbours: all the triangles surrounding vi
      stackN = 0
      DO iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        map( ti) = 1
        stackN = stackN + 1
        stack( stackN) = ti
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        ti = stack( i)
        n_c = n_c + 1
        i_c( n_c) = ti
        x_c( n_c) = mesh%TriGC( ti,1)
        y_c( n_c) = mesh%TriGC( ti,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        ti = i_c( i)
        col = mesh%ti2n( ti)
        CALL add_entry_CSR_dist( mesh%M_map_b_a, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( mesh%M_ddx_b_a, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_b_a, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_a

  SUBROUTINE calc_matrix_operators_mesh_b_b( mesh)
    ! Calculate d/dx, and d/dy matrix operators between the b-grid (triangles) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_b'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: tj
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nTri)                     :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 2

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri    ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nTri    ! to
    nrows_loc       = mesh%nTri_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_ddx_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_ddx_b_b%i1, mesh%M_ddx_b_b%i2

      ! The triangle represented by this matrix row
      ti = mesh%n2ti( row)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Clean up previous map
      DO i = 1, stackN
        tj = stack( i)
        map( tj) = 0
      END DO

      ! Initialise the list of neighbours: just ti itself
      map( ti)  = 1
      stackN    = 1
      stack( 1) = ti

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        tj = stack( i)
        IF (tj == ti) CYCLE
        n_c = n_c + 1
        i_c( n_c) = tj
        x_c( n_c) = mesh%TriGC( tj,1)
        y_c( n_c) = mesh%TriGC( tj,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( mesh%M_ddx_b_b, row, row, Nfx_i)
      CALL add_entry_CSR_dist( mesh%M_ddy_b_b, row, row, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        tj = i_c( i)
        col = mesh%ti2n( tj)
        CALL add_entry_CSR_dist( mesh%M_ddx_b_b, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_b_b, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_b

  SUBROUTINE calc_matrix_operators_mesh_b_b_2nd_order( mesh)
    ! Calculate 2nd-order accurate d/dx, d/dy, d2/dx2, d2/dxdy, and d2/dy2 matrix operators between the b-grid (triangles) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_b_2nd_order'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: tj
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nTri)                     :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 5

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri    ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nTri    ! to
    nrows_loc       = mesh%nTri_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M2_ddx_b_b   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_ddy_b_b   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dx2_b_b , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dxdy_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dy2_b_b , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))
    ALLOCATE( Nfxx_c( n_neighbours_max))
    ALLOCATE( Nfxy_c( n_neighbours_max))
    ALLOCATE( Nfyy_c( n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M2_ddx_b_b%i1, mesh%M2_ddx_b_b%i2

      ! The triangle represented by this matrix row
      ti = mesh%n2ti( row)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Clean up previous map
      DO i = 1, stackN
        tj = stack( i)
        map( tj) = 0
      END DO

      ! Initialise the list of neighbours: just ti itself
      map( ti)  = 1
      stackN    = 1
      stack( 1) = ti

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        tj = stack( i)
        IF (tj == ti) CYCLE
        n_c = n_c + 1
        i_c( n_c) = tj
        x_c( n_c) = mesh%TriGC( tj,1)
        y_c( n_c) = mesh%TriGC( tj,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_2nd_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i, Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( mesh%M2_ddx_b_b   , row, row, Nfx_i )
      CALL add_entry_CSR_dist( mesh%M2_ddy_b_b   , row, row, Nfy_i )
      CALL add_entry_CSR_dist( mesh%M2_d2dx2_b_b , row, row, Nfxx_i)
      CALL add_entry_CSR_dist( mesh%M2_d2dxdy_b_b, row, row, Nfxy_i)
      CALL add_entry_CSR_dist( mesh%M2_d2dy2_b_b , row, row, Nfyy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        tj = i_c( i)
        col = mesh%ti2n( tj)
        CALL add_entry_CSR_dist( mesh%M2_ddx_b_b   , row, col, Nfx_c(  i))
        CALL add_entry_CSR_dist( mesh%M2_ddy_b_b   , row, col, Nfy_c(  i))
        CALL add_entry_CSR_dist( mesh%M2_d2dx2_b_b , row, col, Nfxx_c( i))
        CALL add_entry_CSR_dist( mesh%M2_d2dxdy_b_b, row, col, Nfxy_c( i))
        CALL add_entry_CSR_dist( mesh%M2_d2dy2_b_b , row, col, Nfyy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    DEALLOCATE( Nfxx_c)
    DEALLOCATE( Nfxy_c)
    DEALLOCATE( Nfyy_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_b_2nd_order

  SUBROUTINE calc_matrix_operators_mesh_b_c( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the c-grid (edges)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_c'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ei
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, ti
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nTri)                     :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri      ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 6
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_map_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_map_b_c%i1, mesh%M_map_b_c%i2

      ! The edge represented by this matrix row
      ei = mesh%n2ei( row)
      x  = mesh%E( ei,1)
      y  = mesh%E( ei,2)

      ! Clean up previous map
      DO i = 1, stackN
        ti = stack( i)
        map( ti) = 0
      END DO

      ! Initialise the list of neighbours: the two triangles adjacent to edge ei
      stackN = 0
      DO n = 1, 2
        ti = mesh%ETri( ei,n)
        IF (ti == 0) CYCLE
        map( ti) = 1
        stackN = stackN + 1
        stack( stackN ) = ti
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        ti = stack( i)
        n_c = n_c + 1
        i_c( n_c) = ti
        x_c( n_c) = mesh%TriGC( ti,1)
        y_c( n_c) = mesh%TriGC( ti,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        ti = i_c( i)
        col = mesh%ti2n( ti)
        CALL add_entry_CSR_dist( mesh%M_map_b_c, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( mesh%M_ddx_b_c, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_b_c, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_c

  SUBROUTINE calc_matrix_operators_mesh_c_a( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_c_a'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: vi
    REAL(dp)                                           :: x, y
    INTEGER                                            :: ci, ei
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nE)                       :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i,j
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nE      ! from
    ncols_loc       = mesh%nE_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_map_c_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_c_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_c_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_map_c_a%i1, mesh%M_map_c_a%i2

      ! The vertex represented by this matrix row
      vi = mesh%n2vi( row)
      x  = mesh%V( vi,1)
      y  = mesh%V( vi,2)

      ! Clean up previous map
      DO i = 1, stackN
        ei = stack( i)
        map( ei) = 0
      END DO

      ! Initialise the list of neighbours: all the edges surrounding vi
      stackN = 0
      DO ci = 1, mesh%nC( vi)
        ei = mesh%VE( vi,ci)
        map( ei) = 1
        stackN = stackN + 1
        stack( stackN) = ei
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN < n_neighbours_min)
        CALL extend_group_single_iteration_c( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        ei = stack( i)
        n_c = n_c + 1
        i_c( n_c) = ei
        x_c( n_c) = mesh%E( ei,1)
        y_c( n_c) = mesh%E( ei,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        ei = i_c( i)
        col = mesh%ei2n( ei)
        CALL add_entry_CSR_dist( mesh%M_map_c_a, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( mesh%M_ddx_c_a, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_c_a, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_c_a

  SUBROUTINE calc_matrix_operators_mesh_c_b( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_c_b'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, n2, vi, vj, ci, ei
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nE)                       :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nE        ! from
    ncols_loc       = mesh%nE_loc
    nrows           = mesh%nTri      ! to
    nrows_loc       = mesh%nTri_loc
    nnz_per_row_est = 3
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_map_c_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_c_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_c_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_map_c_b%i1, mesh%M_map_c_b%i2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Clean up previous map
      DO i = 1, stackN
        ei = stack( i)
        map( ei) = 0
      END DO

      ! Initialise the list of neighbours: the three edges spanning ti
      stackN = 0
      DO n = 1, 3
        n2 = n+1
        IF (n2 == 4) n2 = 1
        vi = mesh%Tri( ti,n )
        vj = mesh%Tri( ti,n2)
        ei = 0
        DO ci = 1, mesh%nC( vi)
          IF (mesh%C( vi,ci) == vj) THEN
            ei = mesh%VE( vi,ci)
            EXIT
          END IF
        END DO
        ! Safety
        IF (ei == 0) CALL crash('couldnt find edge connecting vi and vj!')
        map( ei) = 2
        stackN = stackN + 1
        stack( stackN) = ei
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN < n_neighbours_min)
        CALL extend_group_single_iteration_c( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        ei = stack( i)
        n_c = n_c + 1
        i_c( n_c) = ei
        x_c( n_c) = mesh%E( ei,1)
        y_c( n_c) = mesh%E( ei,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        ei = i_c( i)
        col = mesh%ei2n( ei)
        CALL add_entry_CSR_dist( mesh%M_map_c_b, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( mesh%M_ddx_c_b, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_c_b, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_c_b

  SUBROUTINE calc_matrix_operators_mesh_c_c( mesh)
    ! Calculate d/dx, and d/dy matrix operators between the c-grid (edges) and the c-grid (edges)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_c_c'
    INTEGER                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ei
    REAL(dp)                                           :: x, y
    INTEGER                                            :: ej
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER,  DIMENSION(mesh%nE)                       :: map, stack
    INTEGER                                            :: stackN
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 2

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nE      ! from
    ncols_loc       = mesh%nE_loc
    nrows           = mesh%nE      ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( mesh%M_ddx_c_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_c_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    DO row = mesh%M_ddx_c_c%i1, mesh%M_ddx_c_c%i2

      ! The edge represented by this matrix row
      ei = mesh%n2ei( row)
      x  = mesh%E( ei,1)
      y  = mesh%E( ei,2)

      ! Clean up previous map
      DO i = 1, stackN
        ej = stack( i)
        map( ej) = 0
      END DO

      ! Initialise the list of neighbours: just ei itself
      map( ei)  = 1
      stackN    = 1
      stack( 1) = ei

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (stackN - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_c( mesh, map, stack, stackN)
      END DO

      ! Safety
      IF (stackN - 1 > n_neighbours_max) CALL crash('expanded local neighbourhood too far!')

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, stackN
        IF (n_c == n_neighbours_max) EXIT
        ej = stack( i)
        IF (ej == ei) CYCLE
        n_c = n_c + 1
        i_c( n_c) = ej
        x_c( n_c) = mesh%E( ej,1)
        y_c( n_c) = mesh%E( ej,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( mesh%M_ddx_c_c, row, row, Nfx_i)
      CALL add_entry_CSR_dist( mesh%M_ddy_c_c, row, row, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        ej = i_c( i)
        col = mesh%ei2n( ej)
        CALL add_entry_CSR_dist( mesh%M_ddx_c_c, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( mesh%M_ddy_c_c, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_c_c

! == Calculate field-to-vector-form translation tables

  SUBROUTINE calc_field_to_vector_form_translation_tables( mesh, nz)
    ! Calculate grid-cell-to-matrix-row translation tables

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: nz

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_field_to_vector_form_translation_tables'
    INTEGER                                            :: vi,ti,ei,k,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Grid sizes
    mesh%nna     = mesh%nV
    mesh%nnauv   = mesh%nV            * 2
    mesh%nnak    = mesh%nV   *  nz
    mesh%nnakuv  = mesh%nV   *  nz    * 2
    mesh%nnaks   = mesh%nV   * (nz-1)
    mesh%nnaksuv = mesh%nV   * (nz-1) * 2

    mesh%nnb     = mesh%nTri
    mesh%nnbuv   = mesh%nTri          * 2
    mesh%nnbk    = mesh%nTri *  nz
    mesh%nnbkuv  = mesh%nTri *  nz    * 2
    mesh%nnbks   = mesh%nTri * (nz-1)
    mesh%nnbksuv = mesh%nTri * (nz-1) * 2

    mesh%nnc     = mesh%nE
    mesh%nncuv   = mesh%nE            * 2
    mesh%nnck    = mesh%nE   *  nz
    mesh%nnckuv  = mesh%nE   *  nz    * 2
    mesh%nncks   = mesh%nE   * (nz-1)
    mesh%nncksuv = mesh%nE   * (nz-1) * 2

    ! Allocate shared memory
    ALLOCATE( mesh%n2vi(     mesh%nna       ), source = 0)
    ALLOCATE( mesh%n2viuv(   mesh%nnauv  , 2), source = 0)
    ALLOCATE( mesh%n2vik(    mesh%nnak   , 2), source = 0)
    ALLOCATE( mesh%n2vikuv(  mesh%nnakuv , 3), source = 0)
    ALLOCATE( mesh%n2viks(   mesh%nnaks  , 2), source = 0)
    ALLOCATE( mesh%n2viksuv( mesh%nnaksuv, 3), source = 0)

    ALLOCATE( mesh%n2ti(     mesh%nnb       ), source = 0)
    ALLOCATE( mesh%n2tiuv(   mesh%nnbuv  , 2), source = 0)
    ALLOCATE( mesh%n2tik(    mesh%nnbk   , 2), source = 0)
    ALLOCATE( mesh%n2tikuv(  mesh%nnbkuv , 3), source = 0)
    ALLOCATE( mesh%n2tiks(   mesh%nnbks  , 2), source = 0)
    ALLOCATE( mesh%n2tiksuv( mesh%nnbksuv, 3), source = 0)

    ALLOCATE( mesh%n2ei(     mesh%nnc       ), source = 0)
    ALLOCATE( mesh%n2eiuv(   mesh%nncuv  , 2), source = 0)
    ALLOCATE( mesh%n2eik(    mesh%nnck   , 2), source = 0)
    ALLOCATE( mesh%n2eikuv(  mesh%nnckuv , 3), source = 0)
    ALLOCATE( mesh%n2eiks(   mesh%nncks  , 2), source = 0)
    ALLOCATE( mesh%n2eiksuv( mesh%nncksuv, 3), source = 0)

    ALLOCATE( mesh%vi2n(     mesh%nV           ), source = 0)
    ALLOCATE( mesh%viuv2n(   mesh%nV        , 2), source = 0)
    ALLOCATE( mesh%vik2n(    mesh%nV  , nz     ), source = 0)
    ALLOCATE( mesh%vikuv2n(  mesh%nV  , nz  , 2), source = 0)
    ALLOCATE( mesh%viks2n(   mesh%nV  , nz-1   ), source = 0)
    ALLOCATE( mesh%viksuv2n( mesh%nV  , nz-1, 2), source = 0)

    ALLOCATE( mesh%ti2n(     mesh%nTri         ), source = 0)
    ALLOCATE( mesh%tiuv2n(   mesh%nTri      , 2), source = 0)
    ALLOCATE( mesh%tik2n(    mesh%nTri, nz     ), source = 0)
    ALLOCATE( mesh%tikuv2n(  mesh%nTri, nz  , 2), source = 0)
    ALLOCATE( mesh%tiks2n(   mesh%nTri, nz-1   ), source = 0)
    ALLOCATE( mesh%tiksuv2n( mesh%nTri, nz-1, 2), source = 0)

    ALLOCATE( mesh%ei2n(     mesh%nE           ), source = 0)
    ALLOCATE( mesh%eiuv2n(   mesh%nE        , 2), source = 0)
    ALLOCATE( mesh%eik2n(    mesh%nE  , nz     ), source = 0)
    ALLOCATE( mesh%eikuv2n(  mesh%nE  , nz  , 2), source = 0)
    ALLOCATE( mesh%eiks2n(   mesh%nE  , nz-1   ), source = 0)
    ALLOCATE( mesh%eiksuv2n( mesh%nE  , nz-1, 2), source = 0)

! == a-grid (vertices)

  ! == 2-D

    ! == scalar

      n = 0
      DO vi = 1, mesh%nV
        n = n+1
        mesh%vi2n( vi) = n
        mesh%n2vi( n ) = vi
      END DO

    ! == vector

      n = 0
      DO vi = 1, mesh%nV
        DO uv = 1, 2
          n = n+1
          mesh%viuv2n( vi,uv) = n
          mesh%n2viuv( n,1) = vi
          mesh%n2viuv( n,2) = uv
        END DO
      END DO

  ! == 3-D regular

    ! == scalar

      n = 0
      DO vi = 1, mesh%nV
        DO k = 1, nz
          n = n+1
          mesh%vik2n( vi,k) = n
          mesh%n2vik( n,1) = vi
          mesh%n2vik( n,2) = k
        END DO
      END DO

    ! == vector

      n = 0
      DO vi = 1, mesh%nV
        DO k = 1, nz
          DO uv = 1, 2
            n = n+1
            mesh%vikuv2n( vi,k,uv) = n
            mesh%n2vikuv( n,1) = vi
            mesh%n2vikuv( n,2) = k
            mesh%n2vikuv( n,3) = uv
          END DO
        END DO
      END DO

  ! == 3-D staggered

    ! == scalar

      n = 0
      DO vi = 1, mesh%nV
        DO ks = 1, nz-1
          n = n+1
          mesh%viks2n( vi,ks) = n
          mesh%n2viks( n,1) = vi
          mesh%n2viks( n,2) = ks
        END DO
      END DO

    ! == vector

      n = 0
      DO vi = 1, mesh%nV
        DO ks = 1, nz-1
          DO uv = 1, 2
            n = n+1
            mesh%viksuv2n( vi,ks,uv) = n
            mesh%n2viksuv( n,1) = vi
            mesh%n2viksuv( n,2) = ks
            mesh%n2viksuv( n,3) = uv
          END DO
        END DO
      END DO

! == b-grid (triangles)

  ! == 2-D

    ! == scalar

      n = 0
      DO ti = 1, mesh%nTri
        n = n+1
        mesh%ti2n( ti) = n
        mesh%n2ti( n ) = ti
      END DO

    ! == vector

      n = 0
      DO ti = 1, mesh%nTri
        DO uv = 1, 2
          n = n+1
          mesh%tiuv2n( ti,uv) = n
          mesh%n2tiuv( n,1) = ti
          mesh%n2tiuv( n,2) = uv
        END DO
      END DO

  ! == 3-D regular

    ! == scalar

      n = 0
      DO ti = 1, mesh%nTri
        DO k = 1, nz
          n = n+1
          mesh%tik2n( ti,k) = n
          mesh%n2tik( n,1) = ti
          mesh%n2tik( n,2) = k
        END DO
      END DO

    ! == vector

      n = 0
      DO ti = 1, mesh%nTri
        DO k = 1, nz
          DO uv = 1, 2
            n = n+1
            mesh%tikuv2n( ti,k,uv) = n
            mesh%n2tikuv( n,1) = ti
            mesh%n2tikuv( n,2) = k
            mesh%n2tikuv( n,3) = uv
          END DO
        END DO
      END DO

  ! == 3-D staggered

    ! == scalar

      n = 0
      DO ti = 1, mesh%nTri
        DO ks = 1, nz-1
          n = n+1
          mesh%tiks2n( ti,ks) = n
          mesh%n2tiks( n,1) = ti
          mesh%n2tiks( n,2) = ks
        END DO
      END DO

    ! == vector

      n = 0
      DO ti = 1, mesh%nTri
        DO ks = 1, nz-1
          DO uv = 1, 2
            n = n+1
            mesh%tiksuv2n( ti,ks,uv) = n
            mesh%n2tiksuv( n,1) = ti
            mesh%n2tiksuv( n,2) = ks
            mesh%n2tiksuv( n,3) = uv
          END DO
        END DO
      END DO

! == c-grid (edges)

  ! == 2-D

    ! == scalar

      n = 0
      DO ei = 1, mesh%nE
        n = n+1
        mesh%ei2n( ei) = n
        mesh%n2ei( n ) = ei
      END DO

    ! == vector

      n = 0
      DO ei = 1, mesh%nE
        DO uv = 1, 2
          n = n+1
          mesh%eiuv2n( ei,uv) = n
          mesh%n2eiuv( n,1) = ei
          mesh%n2eiuv( n,2) = uv
        END DO
      END DO

  ! == 3-D regular

    ! == scalar

      n = 0
      DO ei = 1, mesh%nE
        DO k = 1, nz
          n = n+1
          mesh%eik2n( ei,k) = n
          mesh%n2eik( n,1) = ei
          mesh%n2eik( n,2) = k
        END DO
      END DO

    ! == vector

      n = 0
      DO ei = 1, mesh%nE
        DO k = 1, nz
          DO uv = 1, 2
            n = n+1
            mesh%eikuv2n( ei,k,uv) = n
            mesh%n2eikuv( n,1) = ei
            mesh%n2eikuv( n,2) = k
            mesh%n2eikuv( n,3) = uv
          END DO
        END DO
      END DO

  ! == 3-D staggered

    ! == scalar

      n = 0
      DO ei = 1, mesh%nE
        DO ks = 1, nz-1
          n = n+1
          mesh%eiks2n( ei,ks) = n
          mesh%n2eiks( n,1) = ei
          mesh%n2eiks( n,2) = ks
        END DO
      END DO

    ! == vector

      n = 0
      DO ei = 1, mesh%nE
        DO ks = 1, nz-1
          DO uv = 1, 2
            n = n+1
            mesh%eiksuv2n( ei,ks,uv) = n
            mesh%n2eiksuv( n,1) = ei
            mesh%n2eiksuv( n,2) = ks
            mesh%n2eiksuv( n,3) = uv
          END DO
        END DO
      END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_field_to_vector_form_translation_tables

! == Shape functions

  SUBROUTINE calc_shape_functions_1D_reg_2nd_order( x, n_max, n_c, x_c, Nfx_i, Nfxx_i, Nfx_c, Nfxx_c)
    ! Calculate shape functions...
    ! ...in one dimension...
    ! ...on the regular grid (i.e. f is known)...
    ! ...to 2nd-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x          ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c        ! Coordinates of the surrounding points where we know f
    REAL(dp),                            INTENT(OUT)   :: Nfx_i      ! d/dx   shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfxx_i     ! d2/dx2 shape function for the point [x,y]
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfxx_c     ! d2/dx2 shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, w
    REAL(dp), DIMENSION(2,2)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 2) CALL crash('calc_shape_functions_1D_reg_2nd_order needs at least 2 neighbours!')

    ! Calculate distances relative to x
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (ABS( dx( ci))**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + w(ci)**2 *       dx( ci)    *       dx( ci)
      ATWTWA( 1,2) = ATWTWA( 1,2) + w(ci)**2 *       dx( ci)    * 1/2 * dx( ci)**2

      ATWTWA( 2,1) = ATWTWA( 2,1) + w(ci)**2 * 1/2 * dx( ci)**2 *       dx( ci)
      ATWTWA( 2,2) = ATWTWA( 2,2) + w(ci)**2 * 1/2 * dx( ci)**2 * 1/2 * dx( ci)**2
    END DO

    ! Invert ATWTWA to find M
    M = calc_matrix_inverse_2_by_2( ATWTWA)

    ! Calculate shape functions
    Nfx_c   = 0._dp
    Nfxx_c  = 0._dp
    DO ci = 1, n_c
      Nfx_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *        dx( ci)   ) + &
        (M( 1,2) * 1/2  * dx( ci)**2))
      Nfxx_c(  ci) = w( ci)**2 * ( &
        (M( 2,1) *        dx( ci)   ) + &
        (M( 2,2) * 1/2  * dx( ci)**2))
    END DO

    Nfx_i  = -SUM( Nfx_c )
    Nfxx_i = -SUM( Nfxx_c)

  END SUBROUTINE calc_shape_functions_1D_reg_2nd_order

  SUBROUTINE calc_shape_functions_1D_stag_2nd_order( x, n_max, n_c, x_c, Nf_c, Nfx_c)
    ! Calculate shape functions...
    ! ...in one dimension...
    ! ...on the staggered grid (i.e. f is not known)...
    ! ...to 2nd-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x          ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c      ! Coordinates of the surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nf_c       ! map    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, w
    REAL(dp), DIMENSION(2,2)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 2) CALL crash('calc_shape_functions_1D_stag_2nd_order needs at least 2 neighbours!')

    ! Calculate distances relative to x
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (ABS( dx( ci))**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + (w( ci)**2 * 1       * 1      )
      ATWTWA( 1,2) = ATWTWA( 1,2) + (w( ci)**2 * 1       * dx( ci))

      ATWTWA( 2,1) = ATWTWA( 2,1) + (w( ci)**2 * dx( ci) * 1      )
      ATWTWA( 2,2) = ATWTWA( 2,2) + (w( ci)**2 * dx( ci) * dx( ci))
    END DO

    ! Invert ATWTWA to find M
    M = calc_matrix_inverse_2_by_2( ATWTWA)

    ! Calculate shape functions
    Nf_c   = 0._dp
    Nfx_c  = 0._dp
    DO ci = 1, n_c
      Nf_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *        1         ) + &
        (M( 1,2) *        dx( ci)**2))
      Nfx_c(  ci) = w( ci)**2 * ( &
        (M( 2,1) *        1         ) + &
        (M( 2,2) *        dx( ci)**2))
    END DO

  END SUBROUTINE calc_shape_functions_1D_stag_2nd_order

  SUBROUTINE calc_shape_functions_2D_reg_1st_order( x, y, n_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)
    ! Calculate shape functions...
    ! ...in two dimensions...
    ! ...on the regular grid (i.e. f is known)...
    ! ...to 1st-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x, y       ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c, y_c   ! Coordinates of the surrounding points where we know f
    REAL(dp),                            INTENT(OUT)   :: Nfx_i      ! d/dx   shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfy_i      ! d/dy   shape function for the point [x,y]
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfy_c      ! d/dy   shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, dy, w
    REAL(dp), DIMENSION(2,2)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 2) CALL crash('calc_shape_functions_2D_reg_1st_order needs at least 2 neighbours!')

    ! Calculate distances relative to [x,y]
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + w(ci)**2 *       dx( ci)    *       dx( ci)
      ATWTWA( 1,2) = ATWTWA( 1,2) + w(ci)**2 *       dx( ci)    *       dy( ci)

      ATWTWA( 2,1) = ATWTWA( 2,1) + w(ci)**2 *       dy( ci)    *       dx( ci)
      ATWTWA( 2,2) = ATWTWA( 2,2) + w(ci)**2 *       dy( ci)    *       dy( ci)
    END DO

    ! Invert ATWTWA to find M
    M = calc_matrix_inverse_2_by_2( ATWTWA)

    ! Calculate shape functions
    Nfx_c   = 0._dp
    Nfy_c   = 0._dp
    DO ci = 1, n_c
      Nfx_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *        dx( ci)   ) + &
        (M( 1,2) *        dy( ci)   ))
      Nfy_c(   ci) = w( ci)**2 * ( &
        (M( 2,1) *        dx( ci)   ) + &
        (M( 2,2) *        dy( ci)   ))
    END DO

    Nfx_i  = -SUM( Nfx_c )
    Nfy_i  = -SUM( Nfy_c )

  END SUBROUTINE calc_shape_functions_2D_reg_1st_order

  SUBROUTINE calc_shape_functions_2D_reg_2nd_order( x, y, n_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i, Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c)
    ! Calculate shape functions...
    ! ...in two dimensions...
    ! ...on the regular grid (i.e. f is known)...
    ! ...to 2nd-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x, y       ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c, y_c   ! Coordinates of the surrounding points where we know f
    REAL(dp),                            INTENT(OUT)   :: Nfx_i      ! d/dx    shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfy_i      ! d/dy    shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfxx_i     ! d2/dx2  shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfxy_i     ! d2/dxdy shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfyy_i     ! d2/dxy2 shape function for the point [x,y]
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfy_c      ! d/dy    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfxx_c     ! d2/dx2  shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfxy_c     ! d2/dxdy shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfyy_c     ! d2/dy2  shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, dy, w
    REAL(dp), DIMENSION(5,5)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 5) CALL crash('calc_shape_functions_2D_reg_2nd_order needs at least 2 neighbours!')

    ! Calculate distances relative to [x,y]
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c

      ATWTWA( 1,1) = ATWTWA( 1,1) + w( ci)**2 *       dx( ci)                 *       dx( ci)
      ATWTWA( 1,2) = ATWTWA( 1,2) + w( ci)**2 *       dx( ci)                 *                    dy( ci)
      ATWTWA( 1,3) = ATWTWA( 1,3) + w( ci)**2 *       dx( ci)                 * 1/2 * dx( ci)**2
      ATWTWA( 1,4) = ATWTWA( 1,4) + w( ci)**2 *       dx( ci)                 *       dx( ci)    * dy( ci)
      ATWTWA( 1,5) = ATWTWA( 1,5) + w( ci)**2 *       dx( ci)                 * 1/2 *              dy( ci)**2

      ATWTWA( 2,1) = ATWTWA( 2,1) + w( ci)**2 *                    dy( ci)    *       dx( ci)
      ATWTWA( 2,2) = ATWTWA( 2,2) + w( ci)**2 *                    dy( ci)    *                    dy( ci)
      ATWTWA( 2,3) = ATWTWA( 2,3) + w( ci)**2 *                    dy( ci)    * 1/2 * dx( ci)**2
      ATWTWA( 2,4) = ATWTWA( 2,4) + w( ci)**2 *                    dy( ci)    *       dx( ci)    * dy( ci)
      ATWTWA( 2,5) = ATWTWA( 2,5) + w( ci)**2 *                    dy( ci)    * 1/2 *              dy( ci)**2

      ATWTWA( 3,1) = ATWTWA( 3,1) + w( ci)**2 * 1/2 * dx( ci)**2              *       dx( ci)
      ATWTWA( 3,2) = ATWTWA( 3,2) + w( ci)**2 * 1/2 * dx( ci)**2              *                    dy( ci)
      ATWTWA( 3,3) = ATWTWA( 3,3) + w( ci)**2 * 1/2 * dx( ci)**2              * 1/2 * dx( ci)**2
      ATWTWA( 3,4) = ATWTWA( 3,4) + w( ci)**2 * 1/2 * dx( ci)**2              *       dx( ci)    * dy( ci)
      ATWTWA( 3,5) = ATWTWA( 3,5) + w( ci)**2 * 1/2 * dx( ci)**2              * 1/2 *              dy( ci)**2

      ATWTWA( 4,1) = ATWTWA( 4,1) + w( ci)**2 *       dx( ci)    * dy( ci)    *       dx( ci)
      ATWTWA( 4,2) = ATWTWA( 4,2) + w( ci)**2 *       dx( ci)    * dy( ci)    *                    dy( ci)
      ATWTWA( 4,3) = ATWTWA( 4,3) + w( ci)**2 *       dx( ci)    * dy( ci)    * 1/2 * dx( ci)**2
      ATWTWA( 4,4) = ATWTWA( 4,4) + w( ci)**2 *       dx( ci)    * dy( ci)    *       dx( ci)    * dy( ci)
      ATWTWA( 4,5) = ATWTWA( 4,5) + w( ci)**2 *       dx( ci)    * dy( ci)    * 1/2 *              dy( ci)**2

      ATWTWA( 5,1) = ATWTWA( 5,1) + w( ci)**2 * 1/2 *              dy( ci)**2 *       dx( ci)
      ATWTWA( 5,2) = ATWTWA( 5,2) + w( ci)**2 * 1/2 *              dy( ci)**2 *                    dy( ci)
      ATWTWA( 5,3) = ATWTWA( 5,3) + w( ci)**2 * 1/2 *              dy( ci)**2 * 1/2 * dx( ci)**2
      ATWTWA( 5,4) = ATWTWA( 5,4) + w( ci)**2 * 1/2 *              dy( ci)**2 *       dx( ci)    * dy( ci)
      ATWTWA( 5,5) = ATWTWA( 5,5) + w( ci)**2 * 1/2 *              dy( ci)**2 * 1/2 *              dy( ci)**2

    END DO

    ! Invert ATWTWA to find M
    M = calc_matrix_inverse_5_by_5( ATWTWA)

    ! Calculate shape functions

    Nfx_c   = 0._dp
    Nfy_c   = 0._dp
    Nfxx_c  = 0._dp
    Nfxy_c  = 0._dp
    Nfyy_c  = 0._dp

    DO ci = 1, n_c

      Nfx_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *       dx( ci)                ) + &
        (M( 1,2) *                    dy( ci)   ) + &
        (M( 1,3) * 1/2 * dx( ci)**2             ) + &
        (M( 1,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 1,5) * 1/2 *              dy( ci)**2))

      Nfy_c(   ci) = w( ci)**2 * ( &
        (M( 2,1) *       dx( ci)                ) + &
        (M( 2,2) *                    dy( ci)   ) + &
        (M( 2,3) * 1/2 * dx( ci)**2             ) + &
        (M( 2,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 2,5) * 1/2 *              dy( ci)**2))

      Nfxx_c(   ci) = w( ci)**2 * ( &
        (M( 3,1) *       dx( ci)                ) + &
        (M( 3,2) *                    dy( ci)   ) + &
        (M( 3,3) * 1/2 * dx( ci)**2             ) + &
        (M( 3,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 3,5) * 1/2 *              dy( ci)**2))

      Nfxy_c(   ci) = w( ci)**2 * ( &
        (M( 4,1) *       dx( ci)                ) + &
        (M( 4,2) *                    dy( ci)   ) + &
        (M( 4,3) * 1/2 * dx( ci)**2             ) + &
        (M( 4,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 4,5) * 1/2 *              dy( ci)**2))

      Nfyy_c(   ci) = w( ci)**2 * ( &
        (M( 5,1) *       dx( ci)                ) + &
        (M( 5,2) *                    dy( ci)   ) + &
        (M( 5,3) * 1/2 * dx( ci)**2             ) + &
        (M( 5,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 5,5) * 1/2 *              dy( ci)**2))

    END DO

    Nfx_i   = -SUM( Nfx_c  )
    Nfy_i   = -SUM( Nfy_c  )
    Nfxx_i  = -SUM( Nfxx_c )
    Nfxy_i  = -SUM( Nfxy_c )
    Nfyy_i  = -SUM( Nfyy_c )

  END SUBROUTINE calc_shape_functions_2D_reg_2nd_order

  SUBROUTINE calc_shape_functions_2D_stag_1st_order( x, y, n_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)
    ! Calculate shape functions...
    ! ...in two dimensions...
    ! ...on the staggered grid (i.e. f is not known)...
    ! ...to 1st-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x, y       ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c, y_c   ! Coordinates of the surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nf_c       ! map    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfy_c      ! d/dy   shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, dy, w
    REAL(dp), DIMENSION(3,3)                           :: ATWTWA, M
    REAL(dp) :: det

    ! Safety
    IF (n_c < 3) CALL crash('calc_shape_functions_2D_stag_1st_order needs at least 3 neighbours!')

    ! Calculate distances relative to [x,y]
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + (w( ci)**2 * 1._dp   * 1._dp  )
      ATWTWA( 1,2) = ATWTWA( 1,2) + (w( ci)**2 * 1._dp   * dx( ci))
      ATWTWA( 1,3) = ATWTWA( 1,3) + (w( ci)**2 * 1._dp   * dy( ci))

      ATWTWA( 2,1) = ATWTWA( 2,1) + (w( ci)**2 * dx( ci) * 1._dp  )
      ATWTWA( 2,2) = ATWTWA( 2,2) + (w( ci)**2 * dx( ci) * dx( ci))
      ATWTWA( 2,3) = ATWTWA( 2,3) + (w( ci)**2 * dx( ci) * dy( ci))

      ATWTWA( 3,1) = ATWTWA( 3,1) + (w( ci)**2 * dy( ci) * 1._dp  )
      ATWTWA( 3,2) = ATWTWA( 3,2) + (w( ci)**2 * dy( ci) * dx( ci))
      ATWTWA( 3,3) = ATWTWA( 3,3) + (w( ci)**2 * dy( ci) * dy( ci))
    END DO

    ! Invert ATWTWA to find M
    M = calc_matrix_inverse_3_by_3( ATWTWA)

    ! Calculate shape functions
    Nf_c    = 0._dp
    Nfx_c   = 0._dp
    Nfy_c   = 0._dp
    DO ci = 1, n_c
      Nf_c(  ci) = w( ci)**2 * ( &
        (M( 1,1) * 1._dp  ) + &
        (M( 1,2) * dx( ci)) + &
        (M( 1,3) * dy( ci)))
      Nfx_c(  ci) = w( ci)**2 * ( &
        (M( 2,1) * 1._dp  ) + &
        (M( 2,2) * dx( ci)) + &
        (M( 2,3) * dy( ci)))
      Nfy_c(  ci) = w( ci)**2 * ( &
        (M( 3,1) * 1._dp  ) + &
        (M( 3,2) * dx( ci)) + &
        (M( 3,3) * dy( ci)))
    END DO

  END SUBROUTINE calc_shape_functions_2D_stag_1st_order

END MODULE mesh_operators
