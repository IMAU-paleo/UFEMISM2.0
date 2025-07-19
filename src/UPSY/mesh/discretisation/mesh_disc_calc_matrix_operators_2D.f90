module mesh_disc_calc_matrix_operators_2D

  ! Routines for calculating 2-D matrix operators on the mesh.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, &
    finalise_matrix_CSR_dist
  use mesh_utilities, only: extend_group_single_iteration_a, extend_group_single_iteration_b, &
    extend_group_single_iteration_c
  use shape_functions, only: calc_shape_functions_2D_reg_1st_order, &
    calc_shape_functions_2D_reg_2nd_order, calc_shape_functions_2D_stag_1st_order
  use mesh_translation_tables, only: calc_field_to_vector_form_translation_tables
  use mesh_zeta, only: calc_vertical_operators_reg_1D, calc_vertical_operators_stag_1D, &
    calc_zeta_operators_tridiagonal

  implicit none

  private

  public :: calc_all_matrix_operators_mesh, calc_matrix_operators_mesh_a_b

contains

subroutine calc_all_matrix_operators_mesh( mesh)
  ! Calculate mapping, d/dx, and d/dy matrix operators between all the grids (a,b,c) on the mesh

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_all_matrix_operators_mesh'

  ! Add routine to path
  call init_routine( routine_name)

  ! Matrix operators
  call calc_field_to_vector_form_translation_tables( mesh)

  call calc_matrix_operators_mesh_a_a( mesh)
  call calc_matrix_operators_mesh_a_b( mesh)

  call calc_matrix_operators_mesh_b_a( mesh)
  call calc_matrix_operators_mesh_b_b( mesh)

  call calc_matrix_operators_mesh_b_b_2nd_order( mesh)

  ! Calculate the 1-D zeta operators (needed for thermodynamics)
  call calc_vertical_operators_reg_1D(  mesh)
  call calc_vertical_operators_stag_1D( mesh)

  ! Zeta operators in tridiagonal form for efficient use in thermodynamics
  call calc_zeta_operators_tridiagonal( mesh)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_all_matrix_operators_mesh

subroutine calc_matrix_operators_mesh_a_a( mesh)
  ! Calculate d/dx, and d/dy matrix operators between the a-grid (vertices) and the a-grid (vertices)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_a_a'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
  integer                                 :: row
  integer                                 :: vi
  real(dp)                                :: x, y
  integer                                 :: vj
  integer                                 :: n_neighbours_min
  integer                                 :: n_neighbours_max
  integer,  dimension(mesh%nV)            :: map, stack
  integer                                 :: stackN
  integer                                 :: i
  integer                                 :: n_c
  integer,  dimension(:    ), allocatable :: i_c
  real(dp), dimension(:    ), allocatable :: x_c, y_c
  real(dp)                                :: Nfx_i, Nfy_i
  real(dp), dimension(:    ), allocatable :: Nfx_c, Nfy_c
  LOGICAL                                 :: succeeded
  integer                                 :: col

  ! Add routine to path
  call init_routine( routine_name)

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

  call allocate_matrix_CSR_dist( mesh%M_ddx_a_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_V, pai_y = mesh%pai_V)
  call allocate_matrix_CSR_dist( mesh%M_ddy_a_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_V, pai_y = mesh%pai_V)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for map, stack, and shape functions
  n_neighbours_max = mesh%nC_mem**2
  allocate( i_c(    n_neighbours_max))
  allocate( x_c(    n_neighbours_max))
  allocate( y_c(    n_neighbours_max))
  allocate( Nfx_c(  n_neighbours_max))
  allocate( Nfy_c(  n_neighbours_max))

  map    = 0
  stack  = 0
  stackN = 0

  do row = mesh%M_ddx_a_a%i1, mesh%M_ddx_a_a%i2

    ! The vertex represented by this matrix row
    vi = mesh%n2vi( row)
    x  = mesh%V( vi,1)
    y  = mesh%V( vi,2)

    ! Clean up previous map
    do i = 1, stackN
      vj = stack( i)
      map( vj) = 0
    end do

    ! Initialise the list of neighbours: just vi itself
    map( vi)  = 1
    stackN    = 1
    stack( 1) = vi

    ! Extend outward until enough neighbours are found to calculate the shape functions
    do while (stackN - 1 < n_neighbours_min)
      call extend_group_single_iteration_a( mesh, map, stack, stackN)
      ! Safety
      if (stackN - 1 > n_neighbours_max) call crash('expanded local neighbourhood too far!')
    end do

    ! Calculate shape functions; if this fails, add more neighbours until it succeeds
    succeeded = .false.
    do while (.not. succeeded)

      ! Get the coordinates of the neighbours
      n_c = 0
      do i = 1, stackN
        if (n_c == n_neighbours_max) exit
        vj = stack( i)
        if (vj == vi) cycle
        n_c = n_c + 1
        i_c( n_c) = vj
        x_c( n_c) = mesh%V( vj,1)
        y_c( n_c) = mesh%V( vj,2)
      end do

      ! Calculate shape functions
      call calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, &
        n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c, succeeded)

      ! if the shape functions couldnt be calculated, include more neighbours and try again
      if (.not. succeeded) call extend_group_single_iteration_a( mesh, map, stack, stackN)

    end do ! do while (.not. succeeded)

    ! Fill them into the matrices

    ! Diagonal elements: shape functions for the home element
    call add_entry_CSR_dist( mesh%M_ddx_a_a, row, row, Nfx_i)
    call add_entry_CSR_dist( mesh%M_ddy_a_a, row, row, Nfy_i)

    ! Off-diagonal elements: shape functions for the neighbours
    do i = 1, n_c
      vj = i_c( i)
      col = mesh%vi2n( vj)
      call add_entry_CSR_dist( mesh%M_ddx_a_a, row, col, Nfx_c( i))
      call add_entry_CSR_dist( mesh%M_ddy_a_a, row, col, Nfy_c( i))
    end do

  end do ! do row = row1, row2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_ddx_a_a)
  call finalise_matrix_CSR_dist( mesh%M_ddy_a_a)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_matrix_operators_mesh_a_a

subroutine calc_matrix_operators_mesh_a_b( mesh)
  ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the b-grid (triangles)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_a_b'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
  integer                                 :: row
  integer                                 :: ti
  real(dp)                                :: x, y
  integer                                 :: n, vi
  integer                                 :: n_neighbours_min
  integer                                 :: n_neighbours_max
  integer,  dimension(mesh%nV)            :: map, stack
  integer                                 :: stackN
  integer                                 :: i
  integer                                 :: n_c
  integer,  dimension(:    ), allocatable :: i_c
  real(dp), dimension(:    ), allocatable :: x_c, y_c
  real(dp), dimension(:    ), allocatable :: Nf_c, Nfx_c, Nfy_c
  LOGICAL                                 :: succeeded
  integer                                 :: col

  ! Add routine to path
  call init_routine( routine_name)

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

  call allocate_matrix_CSR_dist( mesh%M_map_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_V, pai_y = mesh%pai_Tri)
  call allocate_matrix_CSR_dist( mesh%M_ddx_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_V, pai_y = mesh%pai_Tri)
  call allocate_matrix_CSR_dist( mesh%M_ddy_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_V, pai_y = mesh%pai_Tri)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for map, stack, and shape functions
  n_neighbours_max = mesh%nC_mem**2
  allocate( i_c(    n_neighbours_max))
  allocate( x_c(    n_neighbours_max))
  allocate( y_c(    n_neighbours_max))
  allocate( Nf_c(   n_neighbours_max))
  allocate( Nfx_c(  n_neighbours_max))
  allocate( Nfy_c(  n_neighbours_max))

  map    = 0
  stack  = 0
  stackN = 0

  do row = mesh%M_map_a_b%i1, mesh%M_map_a_b%i2

    ! The vertex represented by this matrix row
    ti = mesh%n2ti( row)
    x  = mesh%TriGC( ti,1)
    y  = mesh%TriGC( ti,2)

    ! Clean up previous map
    do i = 1, stackN
      vi = stack( i)
      map( vi) = 0
    end do

    ! Initialise the list of neighbours: the three vertices spanning ti
    stackN = 0
    do n = 1, 3
      vi = mesh%Tri( ti,n)
      map( vi) = 1
      stackN = stackN + 1
      stack( stackN) = vi
    end do

    ! Extend outward until enough neighbours are found to calculate the shape functions
    do while (stackN < n_neighbours_min)
      call extend_group_single_iteration_a( mesh, map, stack, stackN)
      ! Safety
      if (stackN - 1 > n_neighbours_max) call crash('expanded local neighbourhood too far!')
    end do

    ! Calculate shape functions; if this fails, add more neighbours until it succeeds
    succeeded = .false.
    do while (.not. succeeded)

      ! Get the coordinates of the neighbours
      n_c = 0
      do i = 1, stackN
        if (n_c == n_neighbours_max) exit
        vi = stack( i)
        n_c = n_c + 1
        i_c( n_c) = vi
        x_c( n_c) = mesh%V( vi,1)
        y_c( n_c) = mesh%V( vi,2)
      end do

      ! Calculate shape functions
      call calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c, succeeded)

      ! if the shape functions couldnt be calculated, include more neighbours and try again
      if (.not. succeeded) call extend_group_single_iteration_a( mesh, map, stack, stackN)

    end do ! do while (.not. succeeded)

    ! Fill them into the matrices

    ! Off-diagonal elements: shape functions for the neighbours
    do i = 1, n_c
      vi = i_c( i)
      col = mesh%vi2n( vi)
      call add_entry_CSR_dist( mesh%M_map_a_b, row, col, Nf_c(  i))
      call add_entry_CSR_dist( mesh%M_ddx_a_b, row, col, Nfx_c( i))
      call add_entry_CSR_dist( mesh%M_ddy_a_b, row, col, Nfy_c( i))
    end do

  end do ! do row = row1, row2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_map_a_b)
  call finalise_matrix_CSR_dist( mesh%M_ddx_a_b)
  call finalise_matrix_CSR_dist( mesh%M_ddy_a_b)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_matrix_operators_mesh_a_b

subroutine calc_matrix_operators_mesh_b_a( mesh)
  ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the a-grid (vertices)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_b_a'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
  integer                                 :: row
  integer                                 :: vi
  real(dp)                                :: x, y
  integer                                 :: iti, ti
  integer                                 :: n_neighbours_min
  integer                                 :: n_neighbours_max
  integer,  dimension(mesh%nTri)          :: map, stack
  integer                                 :: stackN
  integer                                 :: i
  integer                                 :: n_c
  integer,  dimension(:    ), allocatable :: i_c
  real(dp), dimension(:    ), allocatable :: x_c, y_c
  real(dp), dimension(:    ), allocatable :: Nf_c, Nfx_c, Nfy_c
  LOGICAL                                 :: succeeded
  integer                                 :: col

  ! Add routine to path
  call init_routine( routine_name)

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

  call allocate_matrix_CSR_dist( mesh%M_map_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_V)
  call allocate_matrix_CSR_dist( mesh%M_ddx_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_V)
  call allocate_matrix_CSR_dist( mesh%M_ddy_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_V)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for map, stack, and shape functions
  n_neighbours_max = mesh%nC_mem**2
  allocate( i_c(    n_neighbours_max))
  allocate( x_c(    n_neighbours_max))
  allocate( y_c(    n_neighbours_max))
  allocate( Nf_c(   n_neighbours_max))
  allocate( Nfx_c(  n_neighbours_max))
  allocate( Nfy_c(  n_neighbours_max))

  map    = 0
  stack  = 0
  stackN = 0

  do row = mesh%M_map_b_a%i1, mesh%M_map_b_a%i2

    ! The vertex represented by this matrix row
    vi = mesh%n2vi( row)
    x  = mesh%V( vi,1)
    y  = mesh%V( vi,2)

    ! Clean up previous map
    do i = 1, stackN
      ti = stack( i)
      map( ti) = 0
    end do

    ! Initialise the list of neighbours: all the triangles surrounding vi
    stackN = 0
    do iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      map( ti) = 1
      stackN = stackN + 1
      stack( stackN) = ti
    end do

    ! Extend outward until enough neighbours are found to calculate the shape functions
    do while (stackN < n_neighbours_min)
      call extend_group_single_iteration_b( mesh, map, stack, stackN)
      ! Safety
      if (stackN - 1 > n_neighbours_max) call crash('expanded local neighbourhood too far!')
    end do

    ! Calculate shape functions; if this fails, add more neighbours until it succeeds
    succeeded = .false.
    do while (.not. succeeded)

      ! Get the coordinates of the neighbours
      n_c = 0
      do i = 1, stackN
        if (n_c == n_neighbours_max) exit
        ti = stack( i)
        n_c = n_c + 1
        i_c( n_c) = ti
        x_c( n_c) = mesh%TriGC( ti,1)
        y_c( n_c) = mesh%TriGC( ti,2)
      end do

      ! Calculate shape functions
      call calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c, succeeded)

      ! if the shape functions couldnt be calculated, include more neighbours and try again
      if (.not. succeeded) call extend_group_single_iteration_b( mesh, map, stack, stackN)

    end do ! do while (.not. succeeded)

    ! Fill them into the matrices

    ! Off-diagonal elements: shape functions for the neighbours
    do i = 1, n_c
      ti = i_c( i)
      col = mesh%ti2n( ti)
      call add_entry_CSR_dist( mesh%M_map_b_a, row, col, Nf_c(  i))
      call add_entry_CSR_dist( mesh%M_ddx_b_a, row, col, Nfx_c( i))
      call add_entry_CSR_dist( mesh%M_ddy_b_a, row, col, Nfy_c( i))
    end do

  end do ! do row = row1, row2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_map_b_a)
  call finalise_matrix_CSR_dist( mesh%M_ddx_b_a)
  call finalise_matrix_CSR_dist( mesh%M_ddy_b_a)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_matrix_operators_mesh_b_a

subroutine calc_matrix_operators_mesh_b_b( mesh)
  ! Calculate d/dx, and d/dy matrix operators between the b-grid (triangles) and the b-grid (triangles)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_b_b'
  integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
  integer                                 :: row
  integer                                 :: ti
  real(dp)                                :: x, y
  integer                                 :: tj
  integer                                 :: n_neighbours_min
  integer                                 :: n_neighbours_max
  integer,  dimension(mesh%nTri)          :: map, stack
  integer                                 :: stackN
  integer                                 :: i
  integer                                 :: n_c
  integer,  dimension(:    ), allocatable :: i_c
  real(dp), dimension(:    ), allocatable :: x_c, y_c
  real(dp)                                :: Nfx_i, Nfy_i
  real(dp), dimension(:    ), allocatable :: Nfx_c, Nfy_c
  LOGICAL                                 :: succeeded
  integer                                 :: col

  ! Add routine to path
  call init_routine( routine_name)

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

  call allocate_matrix_CSR_dist( mesh%M_ddx_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
  pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)
  call allocate_matrix_CSR_dist( mesh%M_ddy_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
  pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for map, stack, and shape functions
  n_neighbours_max = mesh%nC_mem**2
  allocate( i_c(    n_neighbours_max))
  allocate( x_c(    n_neighbours_max))
  allocate( y_c(    n_neighbours_max))
  allocate( Nfx_c(  n_neighbours_max))
  allocate( Nfy_c(  n_neighbours_max))

  map    = 0
  stack  = 0
  stackN = 0

  do row = mesh%M_ddx_b_b%i1, mesh%M_ddx_b_b%i2

    ! The triangle represented by this matrix row
    ti = mesh%n2ti( row)
    x  = mesh%TriGC( ti,1)
    y  = mesh%TriGC( ti,2)

    ! Clean up previous map
    do i = 1, stackN
      tj = stack( i)
      map( tj) = 0
    end do

    ! Initialise the list of neighbours: just ti itself
    map( ti)  = 1
    stackN    = 1
    stack( 1) = ti

    ! Extend outward until enough neighbours are found to calculate the shape functions
    do while (stackN - 1 < n_neighbours_min)
      call extend_group_single_iteration_b( mesh, map, stack, stackN)
      ! Safety
      if (stackN - 1 > n_neighbours_max) call crash('expanded local neighbourhood too far!')
    end do

    ! Calculate shape functions; if this fails, add more neighbours until it succeeds
    succeeded = .false.
    do while (.not. succeeded)

      ! Get the coordinates of the neighbours
      n_c = 0
      do i = 1, stackN
        if (n_c == n_neighbours_max) exit
        tj = stack( i)
        if (tj == ti) cycle
        n_c = n_c + 1
        i_c( n_c) = tj
        x_c( n_c) = mesh%TriGC( tj,1)
        y_c( n_c) = mesh%TriGC( tj,2)
      end do

      ! Calculate shape functions
      call calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c, succeeded)

      ! if the shape functions couldnt be calculated, include more neighbours and try again
      if (.not. succeeded) call extend_group_single_iteration_b( mesh, map, stack, stackN)

    end do ! do while (.not. succeeded)

    ! Fill them into the matrices

    ! Diagonal elements: shape functions for the home element
    call add_entry_CSR_dist( mesh%M_ddx_b_b, row, row, Nfx_i)
    call add_entry_CSR_dist( mesh%M_ddy_b_b, row, row, Nfy_i)

    ! Off-diagonal elements: shape functions for the neighbours
    do i = 1, n_c
      tj = i_c( i)
      col = mesh%ti2n( tj)
      call add_entry_CSR_dist( mesh%M_ddx_b_b, row, col, Nfx_c( i))
      call add_entry_CSR_dist( mesh%M_ddy_b_b, row, col, Nfy_c( i))
    end do

  end do ! do row = row1, row2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M_ddx_b_b)
  call finalise_matrix_CSR_dist( mesh%M_ddy_b_b)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_matrix_operators_mesh_b_b

subroutine calc_matrix_operators_mesh_b_b_2nd_order( mesh)
  ! Calculate 2nd-order accurate d/dx, d/dy, d2/dx2, d2/dxdy, and d2/dy2 matrix operators between the b-grid (triangles) and the b-grid (triangles)

  ! In/output variables:
  type(type_mesh),                     intent(inout) :: mesh

  ! Local variables:
  character(len=256), parameter                      :: routine_name = 'calc_matrix_operators_mesh_b_b_2nd_order'
  integer                                            :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
  integer                                            :: row
  integer                                            :: ti
  real(dp)                                           :: x, y
  integer                                            :: tj
  integer                                            :: n_neighbours_min
  integer                                            :: n_neighbours_max
  integer,  dimension(mesh%nTri)                     :: map, stack
  integer                                            :: stackN
  integer                                            :: i
  integer                                            :: n_c
  integer,  dimension(:    ), allocatable            :: i_c
  real(dp), dimension(:    ), allocatable            :: x_c, y_c
  real(dp)                                           :: Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i
  real(dp), dimension(:    ), allocatable            :: Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c
  LOGICAL                                            :: succeeded
  integer                                            :: col

  ! Add routine to path
  call init_routine( routine_name)

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

  call allocate_matrix_CSR_dist( mesh%M2_ddx_b_b   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)
  call allocate_matrix_CSR_dist( mesh%M2_ddy_b_b   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)
  call allocate_matrix_CSR_dist( mesh%M2_d2dx2_b_b , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)
  call allocate_matrix_CSR_dist( mesh%M2_d2dxdy_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)
  call allocate_matrix_CSR_dist( mesh%M2_d2dy2_b_b , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
    pai_x = mesh%pai_Tri, pai_y = mesh%pai_Tri)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

  ! allocate memory for map, stack, and shape functions
  n_neighbours_max = mesh%nC_mem**2
  allocate( i_c(    n_neighbours_max))
  allocate( x_c(    n_neighbours_max))
  allocate( y_c(    n_neighbours_max))
  allocate( Nfx_c(  n_neighbours_max))
  allocate( Nfy_c(  n_neighbours_max))
  allocate( Nfxx_c( n_neighbours_max))
  allocate( Nfxy_c( n_neighbours_max))
  allocate( Nfyy_c( n_neighbours_max))

  map    = 0
  stack  = 0
  stackN = 0

  do row = mesh%M2_ddx_b_b%i1, mesh%M2_ddx_b_b%i2

    ! The triangle represented by this matrix row
    ti = mesh%n2ti( row)
    x  = mesh%TriGC( ti,1)
    y  = mesh%TriGC( ti,2)

    ! Clean up previous map
    do i = 1, stackN
      tj = stack( i)
      map( tj) = 0
    end do

    ! Initialise the list of neighbours: just ti itself
    map( ti)  = 1
    stackN    = 1
    stack( 1) = ti

    ! Extend outward until enough neighbours are found to calculate the shape functions
    do while (stackN - 1 < n_neighbours_min)
      call extend_group_single_iteration_b( mesh, map, stack, stackN)
      ! Safety
      if (stackN - 1 > n_neighbours_max) call crash('expanded local neighbourhood too far!')
    end do

    ! Calculate shape functions; if this fails, add more neighbours until it succeeds
    succeeded = .false.
    do while (.not. succeeded)

      ! Get the coordinates of the neighbours
      n_c = 0
      do i = 1, stackN
        if (n_c == n_neighbours_max) exit
        tj = stack( i)
        if (tj == ti) cycle
        n_c = n_c + 1
        i_c( n_c) = tj
        x_c( n_c) = mesh%TriGC( tj,1)
        y_c( n_c) = mesh%TriGC( tj,2)
      end do

      ! Calculate shape functions
      call calc_shape_functions_2D_reg_2nd_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i, Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c, succeeded)

      ! if the shape functions couldnt be calculated, include more neighbours and try again
      if (.not. succeeded) call extend_group_single_iteration_b( mesh, map, stack, stackN)

    end do ! do while (.not. succeeded)

    ! Fill them into the matrices

    ! Diagonal elements: shape functions for the home element
    call add_entry_CSR_dist( mesh%M2_ddx_b_b   , row, row, Nfx_i )
    call add_entry_CSR_dist( mesh%M2_ddy_b_b   , row, row, Nfy_i )
    call add_entry_CSR_dist( mesh%M2_d2dx2_b_b , row, row, Nfxx_i)
    call add_entry_CSR_dist( mesh%M2_d2dxdy_b_b, row, row, Nfxy_i)
    call add_entry_CSR_dist( mesh%M2_d2dy2_b_b , row, row, Nfyy_i)

    ! Off-diagonal elements: shape functions for the neighbours
    do i = 1, n_c
      tj = i_c( i)
      col = mesh%ti2n( tj)
      call add_entry_CSR_dist( mesh%M2_ddx_b_b   , row, col, Nfx_c(  i))
      call add_entry_CSR_dist( mesh%M2_ddy_b_b   , row, col, Nfy_c(  i))
      call add_entry_CSR_dist( mesh%M2_d2dx2_b_b , row, col, Nfxx_c( i))
      call add_entry_CSR_dist( mesh%M2_d2dxdy_b_b, row, col, Nfxy_c( i))
      call add_entry_CSR_dist( mesh%M2_d2dy2_b_b , row, col, Nfyy_c( i))
    end do

  end do ! do row = row1, row2

  ! Crop matrix memory
  call finalise_matrix_CSR_dist( mesh%M2_ddx_b_b   )
  call finalise_matrix_CSR_dist( mesh%M2_ddy_b_b   )
  call finalise_matrix_CSR_dist( mesh%M2_d2dx2_b_b )
  call finalise_matrix_CSR_dist( mesh%M2_d2dxdy_b_b)
  call finalise_matrix_CSR_dist( mesh%M2_d2dy2_b_b )

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_matrix_operators_mesh_b_b_2nd_order

end module mesh_disc_calc_matrix_operators_2D
