MODULE mesh_operators

  ! Routines for calculating matrix operators on the mesh.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh
  USE basic_data_types                                       , ONLY: type_sparse_matrix_CSR_dp
  USE math_utilities                                         , ONLY: calc_matrix_inverse_2_by_2, calc_matrix_inverse_3_by_3, calc_matrix_inverse_5_by_5
  USE CSR_sparse_matrix_utilities                            , ONLY: allocate_matrix_CSR_dist, add_entry_CSR_dist
  USE mesh_utilities                                         , ONLY: extend_group_single_iteration_a, extend_group_single_iteration_b, &
                                                                     extend_group_single_iteration_c

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Calculate mapping and gradient operators between the a-, b-, and c-grids

  SUBROUTINE calc_matrix_operators_mesh_a_a( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_a_a'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    nrows           = mesh%nV      ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_ddx_a_a, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_a_a, nrows, ncols, nnz_est_proc)

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
      vi = mesh%n2vi( row,1)
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
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 3
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_map_a_b, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_a_b, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_a_b, nrows, ncols, nnz_est_proc)

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
      ti = mesh%n2ti( row,1)
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
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    nrows           = mesh%nE        ! to
    nnz_per_row_est = 4
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_map_a_c, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_a_c, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_a_c, nrows, ncols, nnz_est_proc)

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
      ei = mesh%n2ei( row,1)
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
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    nrows           = mesh%nV      ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_map_b_a, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_b_a, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_b_a, nrows, ncols, nnz_est_proc)

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
      vi = mesh%n2vi( row,1)
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
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_b'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, tj
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
    ncols           = mesh%nTri      ! from
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 3
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_ddx_b_b, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_b_b, nrows, ncols, nnz_est_proc)

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

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row,1)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Clean up previous map
      DO i = 1, stackN
        ti = stack( i)
        map( ti) = 0
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
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_b_2nd_order'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    INTEGER                                            :: row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, tj
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
    ncols           = mesh%nTri      ! from
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 10
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M2_ddx_b_b   , nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_ddy_b_b   , nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dx2_b_b , nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dxdy_b_b, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dy2_b_b , nrows, ncols, nnz_est_proc)

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

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row,1)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Clean up previous map
      DO i = 1, stackN
        ti = stack( i)
        map( ti) = 0
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
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    nrows           = mesh%nE        ! to
    nnz_per_row_est = 6
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_map_b_c, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_b_c, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_b_c, nrows, ncols, nnz_est_proc)

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
      ei = mesh%n2ei( row,1)
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
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    nrows           = mesh%nV      ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_map_c_a, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_c_a, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_c_a, nrows, ncols, nnz_est_proc)

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
      vi = mesh%n2vi( row,1)
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
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 3
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_map_c_b, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddx_c_b, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_c_b, nrows, ncols, nnz_est_proc)

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
      ti = mesh%n2ti( row,1)
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
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the c-grid (edges)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_c_c'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
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
    ncols           = mesh%nE     ! from
    nrows           = mesh%nE     ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( mesh%M_ddx_c_c, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( mesh%M_ddy_c_c, nrows, ncols, nnz_est_proc)

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

      ! The vertex represented by this matrix row
      ei = mesh%n2ei( row,1)
      x  = mesh%E( ei,1)
      y  = mesh%E( ei,2)

      ! Clean up previous map
      DO i = 1, stackN
        ei = stack( i)
        map( ei) = 0
      END DO

      ! Initialise the list of neighbours: just aci itself
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

  SUBROUTINE calc_field_to_vector_form_translation_tables( mesh)
    ! Calculate grid-cell-to-matrix-row translation tables

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_grid_cell_to_matrix_row_translation_tables'
    INTEGER                                            :: vi,ti,ei,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Grid sizes
    mesh%nna   = mesh%nV
    mesh%nnb   = mesh%nTri
    mesh%nnc   = mesh%nE
    mesh%nnauv = mesh%nV   * 2
    mesh%nnbuv = mesh%nTri * 2
    mesh%nncuv = mesh%nE   * 2

    ! Allocate memory
    ALLOCATE( mesh%n2vi(   mesh%nna    , 1), source = 0)
    ALLOCATE( mesh%n2ti(   mesh%nnb    , 1), source = 0)
    ALLOCATE( mesh%n2ei(   mesh%nnc    , 1), source = 0)
    ALLOCATE( mesh%n2viuv( mesh%nnauv  , 2), source = 0)
    ALLOCATE( mesh%n2tiuv( mesh%nnbuv  , 2), source = 0)
    ALLOCATE( mesh%n2eiuv( mesh%nncuv  , 2), source = 0)

    ALLOCATE( mesh%vi2n(   mesh%nV        ), source = 0)
    ALLOCATE( mesh%ti2n(   mesh%nTri      ), source = 0)
    ALLOCATE( mesh%ei2n(   mesh%nE        ), source = 0)
    ALLOCATE( mesh%viuv2n( mesh%nV     , 2), source = 0)
    ALLOCATE( mesh%tiuv2n( mesh%nTri   , 2), source = 0)
    ALLOCATE( mesh%eiuv2n( mesh%nE     , 2), source = 0)

  ! == a-grid (vertices)

    ! scalar

    n = 0
    DO vi = 1, mesh%nV
      n = n+1
      mesh%vi2n( vi) = n
      mesh%n2vi( n,1) = vi
    END DO

    ! vector

    n = 0
    DO vi = 1, mesh%nV
      DO uv = 1, 2
        n = n+1
        mesh%viuv2n( vi,uv) = n
        mesh%n2viuv( n,1) = vi
        mesh%n2viuv( n,2) = uv
      END DO
    END DO

  ! == b-grid (triangles)

    ! scalar

    n = 0
    DO ti = 1, mesh%nTri
      n = n+1
      mesh%ti2n( ti) = n
      mesh%n2ti( n,1) = ti
    END DO

    ! vector

    n = 0
    DO ti = 1, mesh%nTri
      DO uv = 1, 2
        n = n+1
        mesh%tiuv2n( ti,uv) = n
        mesh%n2tiuv( n,1) = ti
        mesh%n2tiuv( n,2) = uv
      END DO
    END DO

  ! == c-grid (edges)

    ! scalar

    n = 0
    DO ei = 1, mesh%nE
      n = n+1
      mesh%ei2n( ei) = n
      mesh%n2ei( n,1) = ei
    END DO

    ! vector

    n = 0
    DO ei = 1, mesh%nE
      DO uv = 1, 2
        n = n+1
        mesh%eiuv2n( ei,uv) = n
        mesh%n2eiuv( n,1) = ei
        mesh%n2eiuv( n,2) = uv
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
