module mesh_operators

  ! Routines for calculating matrix operators on the mesh.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use math_utilities, only: calc_shape_functions_2D_reg_1st_order, calc_shape_functions_2D_reg_2nd_order, &
    calc_shape_functions_2D_stag_1st_order
  use CSR_sparse_matrix_utilities, only: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, &
    read_single_row_CSR_dist, crop_matrix_CSR_dist, deallocate_matrix_CSR_dist
  use mesh_utilities, only: extend_group_single_iteration_a, extend_group_single_iteration_b, &
    extend_group_single_iteration_c
  use petsc_basic, only: multiply_CSR_matrix_with_vector_1D, multiply_CSR_matrix_with_vector_2D
  use mesh_zeta, only: calc_vertical_operators_reg_1D, calc_vertical_operators_stag_1D, calc_zeta_operators_tridiagonal
  use ice_model_types, only: type_ice_model
  use mpi_distributed_memory, only: gather_to_all_dp_2D

  implicit none

contains

! ===== subroutines for applying mesh operators =====
! ===================================================

  ! NOTE: not all possible combinations are written, as only a few
  !       are ever needed.

  ! Gradients between a-grid and a-grid
  subroutine ddx_a_a_2D( mesh, d_a, ddx_a)
    ! ddx a 2-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(  out) :: ddx_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddx_a,1) /= mesh%nV_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_a, d_a, ddx_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_a_2D

  subroutine ddx_a_a_3D( mesh, d_a, ddx_a)
    ! ddx a 3-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: ddx_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddx_a,1) /= mesh%nV_loc .or. size( d_a,2) /= size( ddx_a,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_ddx_a_a, d_a, ddx_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_a_3D

  subroutine ddy_a_a_2D( mesh, d_a, ddy_a)
    ! ddy a 2-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: ddy_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddy_a,1) /= mesh%nV_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_a, d_a, ddy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_a_2D

  subroutine ddy_a_a_3D( mesh, d_a, ddy_a)
    ! ddy a 3-D data field from the a-grid (vertices) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    real(dp), dimension(:,:  ),   intent(in   ) :: d_a
    real(dp), dimension(:,:  ),   intent(out  ) :: ddy_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddy_a,1) /= mesh%nV_loc .or. size( d_a,2) /= size( ddy_a,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_ddy_a_a, d_a, ddy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_a_3D

  ! Mapping between a-grid and b-grid
  subroutine map_a_b_2D( mesh, d_a, d_b)
    ! Map a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: d_b

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_a_b_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_a_b, d_a, d_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_a_b_2D

  subroutine map_a_b_3D( mesh, d_a, d_b)
    ! Map a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: d_b

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_a_b_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc .or. size( d_a,2) /= size( d_b,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_map_a_b, d_a, d_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_a_b_3D

  subroutine map_b_a_2D( mesh, d_b, d_a)
    ! Map a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_b
    real(dp), dimension(:    ), intent(out  ) :: d_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_b_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_map_b_a, d_b, d_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_b_a_2D

  subroutine map_b_a_3D( mesh, d_b, d_a)
    ! Map a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_b
    real(dp), dimension(:,:  ), intent(out  ) :: d_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'map_b_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc .or. size( d_a,2) /= size( d_b,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_map_b_a, d_b, d_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_b_a_3D

  ! Gradients between a-grid and b-grid
  subroutine ddx_a_b_2D( mesh, d_a, ddx_b)
    ! ddx a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: ddx_b

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_b_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddx_b,1) /= mesh%nTri_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_b_2D

  subroutine ddx_a_b_3D( mesh, d_a, ddx_b)
    ! ddx a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: ddx_b

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_a_b_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddx_b,1) /= mesh%nTri_loc .or. size( d_a,2) /= size( ddx_b,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_a_b_3D

  subroutine ddx_b_a_2D( mesh, d_b, ddx_a)
    ! ddx a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_b
    real(dp), dimension(:    ), intent(out  ) :: ddx_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_b_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( ddx_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddx_b_a, d_b, ddx_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_b_a_2D

  subroutine ddx_b_a_3D( mesh, d_b, ddx_a)
    ! ddx a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_b
    real(dp), dimension(:,:  ), intent(out  ) :: ddx_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddx_b_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( ddx_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc .or. size( ddx_a,2) /= size( d_b,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddxping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_ddx_b_a, d_b, ddx_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddx_b_a_3D

  subroutine ddy_a_b_2D( mesh, d_a, ddy_b)
    ! ddy a 2-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_a
    real(dp), dimension(:    ), intent(out  ) :: ddy_b

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_b_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddy_b,1) /= mesh%nTri_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_b_2D

  subroutine ddy_a_b_3D( mesh, d_a, ddy_b)
    ! ddy a 3-D data field from the a-grid (vertices) to the b-grid (triangles)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_a
    real(dp), dimension(:,:  ), intent(out  ) :: ddy_b

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_a_b_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_a,1) /= mesh%nV_loc .or. size( ddy_b,1) /= mesh%nTri_loc .or. size( d_a,2) /= size( ddy_b,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_a_b_3D

  subroutine ddy_b_a_2D( mesh, d_b, ddy_a)
    ! ddy a 2-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:    ), intent(in   ) :: d_b
    real(dp), dimension(:    ), intent(out  ) :: ddy_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_b_a_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( ddy_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_1D( mesh%M_ddy_b_a, d_b, ddy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_b_a_2D

  subroutine ddy_b_a_3D( mesh, d_b, ddy_a)
    ! ddy a 3-D data field from the b-grid (triangles) to the a-grid (vertices)

    ! In/output variables:
    type(type_mesh),            intent(in   ) :: mesh
    real(dp), dimension(:,:  ), intent(in   ) :: d_b
    real(dp), dimension(:,:  ), intent(out  ) :: ddy_a

    ! Local variables:
    character(len=256), parameter :: routine_name = 'ddy_b_a_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( ddy_a,1) /= mesh%nV_loc .or. size( d_b,1) /= mesh%nTri_loc .or. size( ddy_a,2) /= size( d_b,2)) then
      call crash('vector and matrix sizes dont match!')
    end if

    ! Perform the ddyping operation as a matrix multiplication
    call multiply_CSR_matrix_with_vector_2D( mesh%M_ddy_b_a, d_b, ddy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine ddy_b_a_3D

! ===== subroutines for calculating 2-D mesh operators =====
! ==========================================================

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
    call calc_matrix_operators_mesh_a_c( mesh)

    call calc_matrix_operators_mesh_b_a( mesh)
    call calc_matrix_operators_mesh_b_b( mesh)
    call calc_matrix_operators_mesh_b_c( mesh)

    call calc_matrix_operators_mesh_c_a( mesh)
    call calc_matrix_operators_mesh_c_b( mesh)
    call calc_matrix_operators_mesh_c_c( mesh)

    call calc_matrix_operators_mesh_b_b_2nd_order( mesh)

    ! Calculate the 1-D zeta operators (needed for thermodynamics)
    call calc_vertical_operators_reg_1D(  mesh)
    call calc_vertical_operators_stag_1D( mesh)

    ! Zeta operators in tridiagonal form for efficient use in thermodynamics
    call calc_zeta_operators_tridiagonal( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_all_matrix_operators_mesh

! == Calculate mapping and gradient operators between the a-, b-, and c-grids

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

    call allocate_matrix_CSR_dist( mesh%M_ddx_a_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_a_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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
    call crop_matrix_CSR_dist( mesh%M_ddx_a_a)
    call crop_matrix_CSR_dist( mesh%M_ddy_a_a)

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

    call allocate_matrix_CSR_dist( mesh%M_map_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddx_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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
    call crop_matrix_CSR_dist( mesh%M_map_a_b)
    call crop_matrix_CSR_dist( mesh%M_ddx_a_b)
    call crop_matrix_CSR_dist( mesh%M_ddy_a_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_mesh_a_b

  subroutine calc_matrix_operators_mesh_a_c( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the c-grid (edges)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_a_c'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                 :: row
    integer                                 :: ei
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
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 4
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( mesh%M_map_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddx_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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

    do row = mesh%M_map_a_c%i1, mesh%M_map_a_c%i2

      ! The edge represented by this matrix row
      ei = mesh%n2ei( row)
      x  = mesh%E( ei,1)
      y  = mesh%E( ei,2)

      ! Clean up previous map
      do i = 1, stackN
        vi = stack( i)
        map( vi) = 0
      end do

      ! Initialise the list of neighbours: the three or four vertices adjacent to edge ei
      stackN = 0
      do n = 1, 4
        vi = mesh%EV( ei,n)
        if (vi == 0) cycle
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
        call add_entry_CSR_dist( mesh%M_map_a_c, row, col, Nf_c(  i))
        call add_entry_CSR_dist( mesh%M_ddx_a_c, row, col, Nfx_c( i))
        call add_entry_CSR_dist( mesh%M_ddy_a_c, row, col, Nfy_c( i))
      end do

    end do ! do row = row1, row2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_map_a_c)
    call crop_matrix_CSR_dist( mesh%M_ddx_a_c)
    call crop_matrix_CSR_dist( mesh%M_ddy_a_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_mesh_a_c

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

    call allocate_matrix_CSR_dist( mesh%M_map_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddx_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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
    call crop_matrix_CSR_dist( mesh%M_map_b_a)
    call crop_matrix_CSR_dist( mesh%M_ddx_b_a)
    call crop_matrix_CSR_dist( mesh%M_ddy_b_a)

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

    call allocate_matrix_CSR_dist( mesh%M_ddx_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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
    call crop_matrix_CSR_dist( mesh%M_ddx_b_b)
    call crop_matrix_CSR_dist( mesh%M_ddy_b_b)

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

    call allocate_matrix_CSR_dist( mesh%M2_ddx_b_b   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_ddy_b_b   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_d2dx2_b_b , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_d2dxdy_b_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_d2dy2_b_b , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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
    call crop_matrix_CSR_dist( mesh%M2_ddx_b_b   )
    call crop_matrix_CSR_dist( mesh%M2_ddy_b_b   )
    call crop_matrix_CSR_dist( mesh%M2_d2dx2_b_b )
    call crop_matrix_CSR_dist( mesh%M2_d2dxdy_b_b)
    call crop_matrix_CSR_dist( mesh%M2_d2dy2_b_b )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_mesh_b_b_2nd_order

  subroutine calc_matrix_operators_mesh_b_c( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the c-grid (edges)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_b_c'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                 :: row
    integer                                 :: ei
    real(dp)                                :: x, y
    integer                                 :: n, ti
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
    ncols           = mesh%nTri      ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 6
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( mesh%M_map_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddx_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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

    do row = mesh%M_map_b_c%i1, mesh%M_map_b_c%i2

      ! The edge represented by this matrix row
      ei = mesh%n2ei( row)
      x  = mesh%E( ei,1)
      y  = mesh%E( ei,2)

      ! Clean up previous map
      do i = 1, stackN
        ti = stack( i)
        map( ti) = 0
      end do

      ! Initialise the list of neighbours: the two triangles adjacent to edge ei
      stackN = 0
      do n = 1, 2
        ti = mesh%ETri( ei,n)
        if (ti == 0) cycle
        map( ti) = 1
        stackN = stackN + 1
        stack( stackN ) = ti
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
        call add_entry_CSR_dist( mesh%M_map_b_c, row, col, Nf_c(  i))
        call add_entry_CSR_dist( mesh%M_ddx_b_c, row, col, Nfx_c( i))
        call add_entry_CSR_dist( mesh%M_ddy_b_c, row, col, Nfy_c( i))
      end do

    end do ! do row = row1, row2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_map_b_c)
    call crop_matrix_CSR_dist( mesh%M_ddx_b_c)
    call crop_matrix_CSR_dist( mesh%M_ddy_b_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_mesh_b_c

  subroutine calc_matrix_operators_mesh_c_a( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the a-grid (vertices)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_c_a'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                 :: row
    integer                                 :: vi
    real(dp)                                :: x, y
    integer                                 :: ci, ei
    integer                                 :: n_neighbours_min
    integer                                 :: n_neighbours_max
    integer,  dimension(mesh%nE)            :: map, stack
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
    ncols           = mesh%nE      ! from
    ncols_loc       = mesh%nE_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( mesh%M_map_c_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddx_c_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_c_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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

    do row = mesh%M_map_c_a%i1, mesh%M_map_c_a%i2

      ! The vertex represented by this matrix row
      vi = mesh%n2vi( row)
      x  = mesh%V( vi,1)
      y  = mesh%V( vi,2)

      ! Clean up previous map
      do i = 1, stackN
        ei = stack( i)
        map( ei) = 0
      end do

      ! Initialise the list of neighbours: all the edges surrounding vi
      stackN = 0
      do ci = 1, mesh%nC( vi)
        ei = mesh%VE( vi,ci)
        map( ei) = 1
        stackN = stackN + 1
        stack( stackN) = ei
      end do

      ! Extend outward until enough neighbours are found to calculate the shape functions
      do while (stackN < n_neighbours_min)
        call extend_group_single_iteration_c( mesh, map, stack, stackN)
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
          ei = stack( i)
          n_c = n_c + 1
          i_c( n_c) = ei
          x_c( n_c) = mesh%E( ei,1)
          y_c( n_c) = mesh%E( ei,2)
        end do

        ! Calculate shape functions
        call calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c, succeeded)

        ! if the shape functions couldnt be calculated, include more neighbours and try again
        if (.not. succeeded) call extend_group_single_iteration_c( mesh, map, stack, stackN)

      end do ! do while (.not. succeeded)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      do i = 1, n_c
        ei = i_c( i)
        col = mesh%ei2n( ei)
        call add_entry_CSR_dist( mesh%M_map_c_a, row, col, Nf_c(  i))
        call add_entry_CSR_dist( mesh%M_ddx_c_a, row, col, Nfx_c( i))
        call add_entry_CSR_dist( mesh%M_ddy_c_a, row, col, Nfy_c( i))
      end do

    end do ! do row = row1, row2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_map_c_a)
    call crop_matrix_CSR_dist( mesh%M_ddx_c_a)
    call crop_matrix_CSR_dist( mesh%M_ddy_c_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_mesh_c_a

  subroutine calc_matrix_operators_mesh_c_b( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the b-grid (triangles)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_c_b'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                 :: row
    integer                                 :: ti
    real(dp)                                :: x, y
    integer                                 :: n, n2, vi, vj, ci, ei
    integer                                 :: n_neighbours_min
    integer                                 :: n_neighbours_max
    integer,  dimension(mesh%nE)            :: map, stack
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
    ncols           = mesh%nE        ! from
    ncols_loc       = mesh%nE_loc
    nrows           = mesh%nTri      ! to
    nrows_loc       = mesh%nTri_loc
    nnz_per_row_est = 3
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( mesh%M_map_c_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddx_c_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_c_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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

    do row = mesh%M_map_c_b%i1, mesh%M_map_c_b%i2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Clean up previous map
      do i = 1, stackN
        ei = stack( i)
        map( ei) = 0
      end do

      ! Initialise the list of neighbours: the three edges spanning ti
      stackN = 0
      do n = 1, 3
        n2 = n+1
        if (n2 == 4) n2 = 1
        vi = mesh%Tri( ti,n )
        vj = mesh%Tri( ti,n2)
        ei = 0
        do ci = 1, mesh%nC( vi)
          if (mesh%C( vi,ci) == vj) then
            ei = mesh%VE( vi,ci)
            exit
          end if
        end do
        ! Safety
        if (ei == 0) call crash('couldnt find edge connecting vi and vj!')
        map( ei) = 2
        stackN = stackN + 1
        stack( stackN) = ei
      end do

      ! Extend outward until enough neighbours are found to calculate the shape functions
      do while (stackN < n_neighbours_min)
        call extend_group_single_iteration_c( mesh, map, stack, stackN)
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
          ei = stack( i)
          n_c = n_c + 1
          i_c( n_c) = ei
          x_c( n_c) = mesh%E( ei,1)
          y_c( n_c) = mesh%E( ei,2)
        end do

        ! Calculate shape functions
        call calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c, succeeded)

        ! if the shape functions couldnt be calculated, include more neighbours and try again
        if (.not. succeeded) call extend_group_single_iteration_c( mesh, map, stack, stackN)

      end do ! do while (.not. succeeded)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      do i = 1, n_c
        ei = i_c( i)
        col = mesh%ei2n( ei)
        call add_entry_CSR_dist( mesh%M_map_c_b, row, col, Nf_c(  i))
        call add_entry_CSR_dist( mesh%M_ddx_c_b, row, col, Nfx_c( i))
        call add_entry_CSR_dist( mesh%M_ddy_c_b, row, col, Nfy_c( i))
      end do

    end do ! do row = row1, row2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_map_c_b)
    call crop_matrix_CSR_dist( mesh%M_ddx_c_b)
    call crop_matrix_CSR_dist( mesh%M_ddy_c_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_mesh_c_b

  subroutine calc_matrix_operators_mesh_c_c( mesh)
    ! Calculate d/dx, and d/dy matrix operators between the c-grid (edges) and the c-grid (edges)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_matrix_operators_mesh_c_c'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                 :: row
    integer                                 :: ei
    real(dp)                                :: x, y
    integer                                 :: ej
    integer                                 :: n_neighbours_min
    integer                                 :: n_neighbours_max
    integer,  dimension(mesh%nE)            :: map, stack
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
    ncols           = mesh%nE      ! from
    ncols_loc       = mesh%nE_loc
    nrows           = mesh%nE      ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( mesh%M_ddx_c_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_c_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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

    do row = mesh%M_ddx_c_c%i1, mesh%M_ddx_c_c%i2

      ! The edge represented by this matrix row
      ei = mesh%n2ei( row)
      x  = mesh%E( ei,1)
      y  = mesh%E( ei,2)

      ! Clean up previous map
      do i = 1, stackN
        ej = stack( i)
        map( ej) = 0
      end do

      ! Initialise the list of neighbours: just ei itself
      map( ei)  = 1
      stackN    = 1
      stack( 1) = ei

      ! Extend outward until enough neighbours are found to calculate the shape functions
      do while (stackN - 1 < n_neighbours_min)
        call extend_group_single_iteration_c( mesh, map, stack, stackN)
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
          ej = stack( i)
          if (ej == ei) cycle
          n_c = n_c + 1
          i_c( n_c) = ej
          x_c( n_c) = mesh%E( ej,1)
          y_c( n_c) = mesh%E( ej,2)
        end do

        ! Calculate shape functions
        call calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c, succeeded)

        ! if the shape functions couldnt be calculated, include more neighbours and try again
        if (.not. succeeded) call extend_group_single_iteration_c( mesh, map, stack, stackN)

      end do ! do while (.not. succeeded)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      call add_entry_CSR_dist( mesh%M_ddx_c_c, row, row, Nfx_i)
      call add_entry_CSR_dist( mesh%M_ddy_c_c, row, row, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      do i = 1, n_c
        ej = i_c( i)
        col = mesh%ei2n( ej)
        call add_entry_CSR_dist( mesh%M_ddx_c_c, row, col, Nfx_c( i))
        call add_entry_CSR_dist( mesh%M_ddy_c_c, row, col, Nfy_c( i))
      end do

    end do ! do row = row1, row2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_ddx_c_c)
    call crop_matrix_CSR_dist( mesh%M_ddy_c_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_matrix_operators_mesh_c_c

! == Calculate field-to-vector-form translation tables

  subroutine calc_field_to_vector_form_translation_tables( mesh)
    ! Calculate grid-cell-to-matrix-row translation tables

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter :: routine_name = 'calc_field_to_vector_form_translation_tables'
    integer                       :: nz,vi,ti,ei,k,ks,uv,n

    ! Add routine to path
    call init_routine( routine_name)

    nz = mesh%nz

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

    ! allocate shared memory
    allocate( mesh%n2vi(     mesh%nna       ), source = 0)
    allocate( mesh%n2viuv(   mesh%nnauv  , 2), source = 0)
    allocate( mesh%n2vik(    mesh%nnak   , 2), source = 0)
    allocate( mesh%n2vikuv(  mesh%nnakuv , 3), source = 0)
    allocate( mesh%n2viks(   mesh%nnaks  , 2), source = 0)
    allocate( mesh%n2viksuv( mesh%nnaksuv, 3), source = 0)

    allocate( mesh%n2ti(     mesh%nnb       ), source = 0)
    allocate( mesh%n2tiuv(   mesh%nnbuv  , 2), source = 0)
    allocate( mesh%n2tik(    mesh%nnbk   , 2), source = 0)
    allocate( mesh%n2tikuv(  mesh%nnbkuv , 3), source = 0)
    allocate( mesh%n2tiks(   mesh%nnbks  , 2), source = 0)
    allocate( mesh%n2tiksuv( mesh%nnbksuv, 3), source = 0)

    allocate( mesh%n2ei(     mesh%nnc       ), source = 0)
    allocate( mesh%n2eiuv(   mesh%nncuv  , 2), source = 0)
    allocate( mesh%n2eik(    mesh%nnck   , 2), source = 0)
    allocate( mesh%n2eikuv(  mesh%nnckuv , 3), source = 0)
    allocate( mesh%n2eiks(   mesh%nncks  , 2), source = 0)
    allocate( mesh%n2eiksuv( mesh%nncksuv, 3), source = 0)

    allocate( mesh%vi2n(     mesh%nV           ), source = 0)
    allocate( mesh%viuv2n(   mesh%nV        , 2), source = 0)
    allocate( mesh%vik2n(    mesh%nV  , nz     ), source = 0)
    allocate( mesh%vikuv2n(  mesh%nV  , nz  , 2), source = 0)
    allocate( mesh%viks2n(   mesh%nV  , nz-1   ), source = 0)
    allocate( mesh%viksuv2n( mesh%nV  , nz-1, 2), source = 0)

    allocate( mesh%ti2n(     mesh%nTri         ), source = 0)
    allocate( mesh%tiuv2n(   mesh%nTri      , 2), source = 0)
    allocate( mesh%tik2n(    mesh%nTri, nz     ), source = 0)
    allocate( mesh%tikuv2n(  mesh%nTri, nz  , 2), source = 0)
    allocate( mesh%tiks2n(   mesh%nTri, nz-1   ), source = 0)
    allocate( mesh%tiksuv2n( mesh%nTri, nz-1, 2), source = 0)

    allocate( mesh%ei2n(     mesh%nE           ), source = 0)
    allocate( mesh%eiuv2n(   mesh%nE        , 2), source = 0)
    allocate( mesh%eik2n(    mesh%nE  , nz     ), source = 0)
    allocate( mesh%eikuv2n(  mesh%nE  , nz  , 2), source = 0)
    allocate( mesh%eiks2n(   mesh%nE  , nz-1   ), source = 0)
    allocate( mesh%eiksuv2n( mesh%nE  , nz-1, 2), source = 0)

    ! == a-grid (vertices)

      ! == 2-D

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        n = n+1
        mesh%vi2n( vi) = n
        mesh%n2vi( n ) = vi
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do uv = 1, 2
          n = n+1
          mesh%viuv2n( vi,uv) = n
          mesh%n2viuv( n,1) = vi
          mesh%n2viuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        do k = 1, nz
          n = n+1
          mesh%vik2n( vi,k) = n
          mesh%n2vik( n,1) = vi
          mesh%n2vik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%vikuv2n( vi,k,uv) = n
            mesh%n2vikuv( n,1) = vi
            mesh%n2vikuv( n,2) = k
            mesh%n2vikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        do ks = 1, nz-1
          n = n+1
          mesh%viks2n( vi,ks) = n
          mesh%n2viks( n,1) = vi
          mesh%n2viks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%viksuv2n( vi,ks,uv) = n
            mesh%n2viksuv( n,1) = vi
            mesh%n2viksuv( n,2) = ks
            mesh%n2viksuv( n,3) = uv
          end do
        end do
      end do

    ! == b-grid (triangles)

      ! == 2-D

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        n = n+1
        mesh%ti2n( ti) = n
        mesh%n2ti( n ) = ti
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do uv = 1, 2
          n = n+1
          mesh%tiuv2n( ti,uv) = n
          mesh%n2tiuv( n,1) = ti
          mesh%n2tiuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        do k = 1, nz
          n = n+1
          mesh%tik2n( ti,k) = n
          mesh%n2tik( n,1) = ti
          mesh%n2tik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%tikuv2n( ti,k,uv) = n
            mesh%n2tikuv( n,1) = ti
            mesh%n2tikuv( n,2) = k
            mesh%n2tikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        do ks = 1, nz-1
          n = n+1
          mesh%tiks2n( ti,ks) = n
          mesh%n2tiks( n,1) = ti
          mesh%n2tiks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%tiksuv2n( ti,ks,uv) = n
            mesh%n2tiksuv( n,1) = ti
            mesh%n2tiksuv( n,2) = ks
            mesh%n2tiksuv( n,3) = uv
          end do
        end do
      end do

    ! == c-grid (edges)

      ! == 2-D

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        n = n+1
        mesh%ei2n( ei) = n
        mesh%n2ei( n ) = ei
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do uv = 1, 2
          n = n+1
          mesh%eiuv2n( ei,uv) = n
          mesh%n2eiuv( n,1) = ei
          mesh%n2eiuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        do k = 1, nz
          n = n+1
          mesh%eik2n( ei,k) = n
          mesh%n2eik( n,1) = ei
          mesh%n2eik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%eikuv2n( ei,k,uv) = n
            mesh%n2eikuv( n,1) = ei
            mesh%n2eikuv( n,2) = k
            mesh%n2eikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        do ks = 1, nz-1
          n = n+1
          mesh%eiks2n( ei,ks) = n
          mesh%n2eiks( n,1) = ei
          mesh%n2eiks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%eiksuv2n( ei,ks,uv) = n
            mesh%n2eiksuv( n,1) = ei
            mesh%n2eiksuv( n,2) = ks
            mesh%n2eiksuv( n,3) = uv
          end do
        end do
      end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_field_to_vector_form_translation_tables

! ===== subroutines for calculating 3-D mesh operators =====
! ==========================================================

  ! == Calculate the operator matrices

  subroutine calc_3D_matrix_operators_mesh( mesh, ice)
    ! Calculate all 3-D gradient operators in Cartesian coordinates

    ! In/output variables:
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in   )    :: ice

    ! Local variables:
    character(len=256), parameter :: routine_name = 'calc_3D_matrix_operators_mesh'

    ! Add routine to path
    call init_routine( routine_name)

    ! bk to ak (for calculating the horizontal stretch/shear strain rates in the BPA)
    call calc_3D_matrix_operators_mesh_bk_ak( mesh, ice)

    ! ak to bk (for calculating the horizontal gradients of the effective viscosity in the BPA)
    call calc_3D_matrix_operators_mesh_ak_bk( mesh, ice)

    ! bk to bks (for calculating the vertical shear strain rates in the BPA)
    call calc_3D_matrix_operators_mesh_bk_bks( mesh, ice)

    ! bks to bk (for calculating the vertical gradient of the effective viscosity in the BPA)
    call calc_3D_matrix_operators_mesh_bks_bk( mesh, ice)

    ! Map between the bks-grid and the ak-grid (for calculating strain rates in the BPA)
    call calc_3D_mapping_operator_mesh_bks_ak( mesh)
    call calc_3D_mapping_operator_mesh_ak_bks( mesh)

    ! bk to bk (for constructing the BPA stiffness matrix)
    call calc_3D_matrix_operators_mesh_bk_bk( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_matrix_operators_mesh

  subroutine calc_3D_matrix_operators_mesh_bk_ak( mesh, ice)
    ! Calculate all 3-D gradient operators in Cartesian coordinates
    !
    ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
    ! are defined as:
    !
    !   xh   = x
    !   yh   = y
    !   zeta = (Hs - z) / Hi
    !
    ! Applying the chain rule to the gradient operators d/dx, d/dy yields:
    !
    !    d/dx = d/dxh + dzeta/dx d/dzeta
    !    d/dy = d/dyh + dzeta/dy d/dzeta
    !
    ! Theoretically, we could convert the basic d/dx, d/dy, d/dzeta matrix operators to PETSc format, perform
    ! matrix multiplications on them (so e.g. M_ddx = M_ddxh + D( dzeta/dx) M_ddzeta), and then convert
    ! the result back to CSR format. However, this is rather cumbersome to do (especially because
    ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
    ! so we'd need to convert them to act on the 3-D mesh first).
    !
    ! Since the resulting operators  are relatively easy to interpret, we can just calculate
    ! their coefficients directly, which is done here.

    ! In/output variables:
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bk_ak'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
    integer                                 :: vi
    integer,  dimension(:    ), allocatable :: single_row_vi_ind
    integer                                 :: single_row_vi_nnz
    real(dp), dimension(:    ), allocatable :: single_row_map_val
    real(dp), dimension(:    ), allocatable :: single_row_ddxh_val
    real(dp), dimension(:    ), allocatable :: single_row_ddyh_val
    integer                                 :: k
    integer,  dimension(:    ), allocatable :: single_row_k_ind
    integer                                 :: single_row_k_nnz
    real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
    integer                                 :: row_vik
    integer                                 :: ii,col_ti,ti,jj,kk,col_tikk
    real(dp)                                :: dzeta_dx, dzeta_dy
    real(dp)                                :: c_map, c_ddxh, c_ddyh, c_ddzeta
    real(dp)                                :: c_ddx, c_ddy

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Deallocate existing matrices if necessary
    if (allocateD( mesh%M_ddx_bk_ak%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddx_bk_ak)
    if (allocateD( mesh%M_ddy_bk_ak%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddy_bk_ak)

    ! Matrix size
    ncols           = mesh%nTri     * mesh%nz ! from
    ncols_loc       = mesh%nTri_loc * mesh%nz
    nrows           = mesh%nV       * mesh%nz ! to
    nrows_loc       = mesh%nV_loc   * mesh%nz
    nnz_est_proc    = mesh%M_map_b_a%nnz * mesh%nz * 3

    call allocate_matrix_CSR_dist( mesh%M_ddx_bk_ak, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_bk_ak, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for single matrix rows
    allocate( single_row_vi_ind(   mesh%nC_mem*2))
    allocate( single_row_map_val(  mesh%nC_mem*2))
    allocate( single_row_ddxh_val( mesh%nC_mem*2))
    allocate( single_row_ddyh_val( mesh%nC_mem*2))

    allocate( single_row_k_ind(      mesh%nz))
    allocate( single_row_ddzeta_val( mesh%nz))

    ! Loop over all vertices
    do vi = mesh%vi1, mesh%vi2

      ! Read coefficients from the 2-D gradient operators for this triangle
      call read_single_row_CSR_dist( mesh%M_map_b_a, vi, single_row_vi_ind, single_row_map_val , single_row_vi_nnz)
      call read_single_row_CSR_dist( mesh%M_ddx_b_a, vi, single_row_vi_ind, single_row_ddxh_val, single_row_vi_nnz)
      call read_single_row_CSR_dist( mesh%M_ddy_b_a, vi, single_row_vi_ind, single_row_ddyh_val, single_row_vi_nnz)

      ! Loop over all layers
      do k = 1, mesh%nz

        ! Vertex vi, layer k corresponds to this matrix row
        row_vik = mesh%vik2n( vi,k)

        ! Read coefficients from the zeta gradient operators for this layer
        call read_single_row_CSR_dist( mesh%M_ddzeta_k_k_1D, k, single_row_k_ind, single_row_ddzeta_val, single_row_k_nnz)

        ! Gradients of zeta at vertex vi, layer k
        dzeta_dx = ice%dzeta_dx_ak( vi,k)
        dzeta_dy = ice%dzeta_dy_ak( vi,k)

        ! Loop over the entire 3-D local neighbourhood, calculate
        ! coefficients for all 3-D matrix operators

        do ii = 1, single_row_vi_nnz

          col_ti = single_row_vi_ind( ii)
          ti = mesh%n2ti( col_ti)

          ! Coefficients for horizontal gradient matrix operators
          c_map  = single_row_map_val(  ii)
          c_ddxh = single_row_ddxh_val( ii)
          c_ddyh = single_row_ddyh_val( ii)

          do jj = 1, single_row_k_nnz

            kk = single_row_k_ind( jj)

            ! Triangle ti, layer kk corresponds to this matrix row
            col_tikk = mesh%tik2n( ti,kk)

            ! Coefficients for vertical gradient matrix operators
            c_ddzeta = single_row_ddzeta_val( jj)

            ! Calculate coefficients
            c_ddx = 0._dp
            c_ddy = 0._dp

            ! Horizontal-only part
            if (kk == k) then
              c_ddx = c_ddx + c_ddxh                      ! Now:  d/dx   =  d/dxh ...
              c_ddy = c_ddy + c_ddyh                      ! Now:  d/dy   =  d/dyh ...
            end if ! if (kk == k) then

            ! Mixed part
            c_ddx = c_ddx + dzeta_dx * c_map * c_ddzeta   ! Now:  d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
            c_ddy = c_ddy + dzeta_dy * c_map * c_ddzeta   ! Now:  d/dy   =  d/dyh    +  dzeta/dy   d/dzeta

            ! Add to CSR matrices
            call add_entry_CSR_dist( mesh%M_ddx_bk_ak, row_vik, col_tikk, c_ddx)
            call add_entry_CSR_dist( mesh%M_ddy_bk_ak, row_vik, col_tikk, c_ddy)

          end do ! do jj = 1, single_row_k_nnz
        end do ! do ii = 1, single_row_vi_nnz

      end do ! do k = 1, mesh%nz
    end do ! do vi = mesh%vi1, mesh%vi2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_ddx_bk_ak)
    call crop_matrix_CSR_dist( mesh%M_ddy_bk_ak)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_matrix_operators_mesh_bk_ak

  subroutine calc_3D_matrix_operators_mesh_ak_bk( mesh, ice)
    ! Calculate all 3-D gradient operators in Cartesian coordinates
    !
    ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
    ! are defined as:
    !
    !   xh   = x
    !   yh   = y
    !   zeta = (Hs - z) / Hi
    !
    ! Applying the chain rule to the gradient operators d/dx, d/dy yields:
    !
    !    d/dx = d/dxh + dzeta/dx d/dzeta
    !    d/dy = d/dyh + dzeta/dy d/dzeta
    !
    ! Theoretically, we could convert the basic d/dx, d/dy, d/dzeta matrix operators to PETSc format, perform
    ! matrix multiplications on them (so e.g. M_ddx = M_ddxh + D( dzeta/dx) M_ddzeta), and then convert
    ! the result back to CSR format. However, this is rather cumbersome to do (especially because
    ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
    ! so we'd need to convert them to act on the 3-D mesh first).
    !
    ! Since the resulting operators  are relatively easy to interpret, we can just calculate
    ! their coefficients directly, which is done here.

    ! In/output variables:
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_ak_bk'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
    integer                                 :: ti
    integer,  dimension(:    ), allocatable :: single_row_ti_ind
    integer                                 :: single_row_ti_nnz
    real(dp), dimension(:    ), allocatable :: single_row_map_val
    real(dp), dimension(:    ), allocatable :: single_row_ddxh_val
    real(dp), dimension(:    ), allocatable :: single_row_ddyh_val
    integer                                 :: k
    integer,  dimension(:    ), allocatable :: single_row_k_ind
    integer                                 :: single_row_k_nnz
    real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
    integer                                 :: row_tik
    integer                                 :: ii,col_vi,vi,jj,kk,col_vikk
    real(dp)                                :: dzeta_dx, dzeta_dy
    real(dp)                                :: c_map, c_ddxh, c_ddyh, c_ddzeta
    real(dp)                                :: c_ddx, c_ddy

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Deallocate existing matrices if necessary
    if (allocateD( mesh%M_ddx_ak_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddx_ak_bk)
    if (allocateD( mesh%M_ddy_ak_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddy_ak_bk)

    ! Matrix size
    ncols           = mesh%nV       * mesh%nz ! from
    ncols_loc       = mesh%nV_loc   * mesh%nz
    nrows           = mesh%nTri     * mesh%nz ! to
    nrows_loc       = mesh%nTri_loc * mesh%nz
    nnz_est_proc    = mesh%M_map_a_b%nnz * mesh%nz * 3

    call allocate_matrix_CSR_dist( mesh%M_ddx_ak_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddy_ak_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for single matrix rows
    allocate( single_row_ti_ind(   mesh%nC_mem*2))
    allocate( single_row_map_val(  mesh%nC_mem*2))
    allocate( single_row_ddxh_val( mesh%nC_mem*2))
    allocate( single_row_ddyh_val( mesh%nC_mem*2))

    allocate( single_row_k_ind(      mesh%nz))
    allocate( single_row_ddzeta_val( mesh%nz))

    ! Loop over all triangles
    do ti = mesh%ti1, mesh%ti2

      ! Read coefficients from the 2-D gradient operators for this triangle
      call read_single_row_CSR_dist( mesh%M_map_a_b, ti, single_row_ti_ind, single_row_map_val , single_row_ti_nnz)
      call read_single_row_CSR_dist( mesh%M_ddx_a_b, ti, single_row_ti_ind, single_row_ddxh_val, single_row_ti_nnz)
      call read_single_row_CSR_dist( mesh%M_ddy_a_b, ti, single_row_ti_ind, single_row_ddyh_val, single_row_ti_nnz)

      ! Loop over all layers
      do k = 1, mesh%nz

        ! Triangle ti, layer k corresponds to this matrix row
        row_tik = mesh%tik2n( ti,k)

        ! Read coefficients from the zeta gradient operators for this layer
        call read_single_row_CSR_dist( mesh%M_ddzeta_k_k_1D, k, single_row_k_ind, single_row_ddzeta_val, single_row_k_nnz)

        ! Gradients of zeta at triangle ti, layer k
        dzeta_dx = ice%dzeta_dx_bk( ti,k)
        dzeta_dy = ice%dzeta_dy_bk( ti,k)

        ! Loop over the entire 3-D local neighbourhood, calculate
        ! coefficients for all 3-D matrix operators

        do ii = 1, single_row_ti_nnz

          col_vi = single_row_ti_ind( ii)
          vi = mesh%n2vi( col_vi)

          ! Coefficients for horizontal gradient matrix operators
          c_map  = single_row_map_val(  ii)
          c_ddxh = single_row_ddxh_val( ii)
          c_ddyh = single_row_ddyh_val( ii)

          do jj = 1, single_row_k_nnz

            kk = single_row_k_ind( jj)

            ! Vertex vi, layer kk corresponds to this matrix row
            col_vikk = mesh%vik2n( vi,kk)

            ! Coefficients for vertical gradient matrix operators
            c_ddzeta = single_row_ddzeta_val( jj)

            ! Calculate coefficients
            c_ddx = 0._dp
            c_ddy = 0._dp

            ! Horizontal-only part
            if (kk == k) then
              c_ddx = c_ddx + c_ddxh                      ! Now:  d/dx   =  d/dxh ...
              c_ddy = c_ddy + c_ddyh                      ! Now:  d/dy   =  d/dyh ...
            end if ! if (kk == k) then

            ! Mixed part
            c_ddx = c_ddx + dzeta_dx * c_map * c_ddzeta   ! Now:  d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
            c_ddy = c_ddy + dzeta_dy * c_map * c_ddzeta   ! Now:  d/dy   =  d/dyh    +  dzeta/dy   d/dzeta

            ! Add to CSR matrices
            call add_entry_CSR_dist( mesh%M_ddx_ak_bk, row_tik, col_vikk, c_ddx)
            call add_entry_CSR_dist( mesh%M_ddy_ak_bk, row_tik, col_vikk, c_ddy)

          end do ! do jj = 1, single_row_k_nnz
        end do ! do ii = 1, single_row_ti_nnz

      end do ! do k = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_ddx_ak_bk)
    call crop_matrix_CSR_dist( mesh%M_ddy_ak_bk)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_matrix_operators_mesh_ak_bk

  subroutine calc_3D_matrix_operators_mesh_bk_bks( mesh, ice)
    ! Calculate all 3-D gradient operators in Cartesian coordinates
    !
    ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
    ! are defined as:
    !
    !   xh   = x
    !   yh   = y
    !   zeta = (Hs - z) / Hi
    !
    ! Applying the chain rule to the gradient operators d/dz, d2/dz2 yields:
    !
    !    d/dz   =  dzeta/dz     d/dzeta
    !   d2/dz2  = (dzeta/dz)^2 d2/dzeta2
    !
    ! Theoretically, we could convert all these basic matrix operators to PETSc format, perform
    ! matrix multiplications on them (so e.g. M_d2dxhdzeta = M_ddxh * M_ddzeta), and then convert
    ! the result back to CSR format. However, this is rather cumbersome to do (especially because
    ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
    ! so we'd need to convert them to act on the 3-D mesh first).
    !
    ! Since the resulting operators  are relatively easy to interpret, we can just calculate
    ! their coefficients directly, which is done here.

    ! In/output variables:
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bk_bks'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                 :: ti
    integer                                 :: ks
    integer,  dimension(:    ), allocatable :: single_row_ks_ind
    integer                                 :: single_row_ks_nnz
    real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
    integer                                 :: row_tiks
    integer                                 :: jj,k,col_tik
    real(dp)                                :: dzeta_dz
    real(dp)                                :: c_ddzeta
    real(dp)                                :: c_ddz

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Deallocate existing matrices if necessary
    if (allocateD( mesh%M_ddz_bk_bks%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddz_bk_bks)

    ! Matrix size
    ncols           = mesh%nTri     *  mesh%nz    ! from
    ncols_loc       = mesh%nTri_loc *  mesh%nz
    nrows           = mesh%nTri     * (mesh%nz-1) ! to
    nrows_loc       = mesh%nTri_loc * (mesh%nz-1)
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( mesh%M_ddz_bk_bks, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for single matrix rows
    allocate( single_row_ks_ind(     mesh%nz))
    allocate( single_row_ddzeta_val( mesh%nz))

    ! Loop over all triangles
    do ti = mesh%ti1, mesh%ti2

      ! Loop over all staggered layers
      do ks = 1, mesh%nz-1

        ! Triangle ti, staggered layer ks corresponds to this matrix row
        row_tiks = mesh%tiks2n( ti,ks)

        ! Read coefficients from the zeta gradient operators for this staggered layer
        call read_single_row_CSR_dist( mesh%M_ddzeta_k_ks_1D, ks, single_row_ks_ind, single_row_ddzeta_val, single_row_ks_nnz)

        ! Gradients of zeta at triangle ti, staggered layer ks
        dzeta_dz = ice%dzeta_dz_bks( ti,ks)

        ! Loop over the vertical local neighbourhood, calculate
        ! coefficients for all 3-D matrix operators

        do jj = 1, single_row_ks_nnz

          k = single_row_ks_ind( jj)

          ! Triangle tj, layer kk corresponds to this matrix row
          col_tik = mesh%tik2n( ti,k)

          ! Coefficients for vertical gradient matrix operators
          c_ddzeta = single_row_ddzeta_val( jj)

          ! Calculate coefficients
          c_ddz = dzeta_dz * c_ddzeta

          ! Add to CSR matrices
          call add_entry_CSR_dist( mesh%M_ddz_bk_bks, row_tiks, col_tik, c_ddz)

        end do ! do jj = 1, single_row_ks_nnz

      end do ! do ks = 1, mesh%nz-1
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    DEallocate( single_row_ks_ind)
    DEallocate( single_row_ddzeta_val)

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_ddz_bk_bks)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_matrix_operators_mesh_bk_bks

  subroutine calc_3D_matrix_operators_mesh_bks_bk( mesh, ice)
    ! Calculate all 3-D gradient operators in Cartesian coordinates
    !
    ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
    ! are defined as:
    !
    !   xh   = x
    !   yh   = y
    !   zeta = (Hs - z) / Hi
    !
    ! Applying the chain rule to the gradient operators d/dz, d2/dz2 yields:
    !
    !    d/dz   =  dzeta/dz     d/dzeta
    !   d2/dz2  = (dzeta/dz)^2 d2/dzeta2
    !
    ! Theoretically, we could convert all these basic matrix operators to PETSc format, perform
    ! matrix multiplications on them (so e.g. M_d2dxhdzeta = M_ddxh * M_ddzeta), and then convert
    ! the result back to CSR format. However, this is rather cumbersome to do (especially because
    ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
    ! so we'd need to convert them to act on the 3-D mesh first).
    !
    ! Since the resulting operators  are relatively easy to interpret, we can just calculate
    ! their coefficients directly, which is done here.

    ! In/output variables:
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bks_bk'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                 :: ti
    integer                                 :: k
    integer,  dimension(:    ), allocatable :: single_row_k_ind
    integer                                 :: single_row_k_nnz
    real(dp), dimension(:    ), allocatable :: single_row_map_val
    real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
    integer                                 :: row_tik
    integer                                 :: jj,ks,col_tiks
    real(dp)                                :: dzeta_dz
    real(dp)                                :: c_ddzeta
    real(dp)                                :: c_ddz, c_map

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Deallocate existing matrices if necessary
    if (allocateD( mesh%M_map_bks_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_map_bks_bk)
    if (allocateD( mesh%M_ddz_bks_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M_ddz_bks_bk)

    ! Matrix size
    ncols           = mesh%nTri     * (mesh%nz-1) ! from
    ncols_loc       = mesh%nTri_loc * (mesh%nz-1)
    nrows           = mesh%nTri     *  mesh%nz    ! to
    nrows_loc       = mesh%nTri_loc *  mesh%nz
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( mesh%M_map_bks_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M_ddz_bks_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for single matrix rows
    allocate( single_row_k_ind(      mesh%nz))
    allocate( single_row_map_val(    mesh%nz))
    allocate( single_row_ddzeta_val( mesh%nz))

    ! Loop over all triangles
    do ti = mesh%ti1, mesh%ti2

      ! Loop over all layers
      do k = 1, mesh%nz

        ! Triangle ti, layer k corresponds to this matrix row
        row_tik = mesh%tik2n( ti,k)

        ! Read coefficients from the zeta gradient operators for this staggered layer
        call read_single_row_CSR_dist( mesh%M_map_ks_k_1D   , k, single_row_k_ind, single_row_map_val   , single_row_k_nnz)
        call read_single_row_CSR_dist( mesh%M_ddzeta_ks_k_1D, k, single_row_k_ind, single_row_ddzeta_val, single_row_k_nnz)

        ! Gradients of zeta at triangle ti, layer k
        dzeta_dz = ice%dzeta_dz_bk( ti,k)

        ! Loop over the vertical local neighbourhood, calculate
        ! coefficients for all 3-D matrix operators

        do jj = 1, single_row_k_nnz

          ks = single_row_k_ind( jj)

          ! Triangle tj, layer kk corresponds to this matrix row
          col_tiks = mesh%tiks2n( ti,ks)

          ! Coefficients for vertical gradient matrix operators
          c_map    = single_row_map_val(    jj)
          c_ddzeta = single_row_ddzeta_val( jj)

          ! Calculate coefficients
          c_ddz = dzeta_dz * c_ddzeta

          ! Add to CSR matrices
          call add_entry_CSR_dist( mesh%M_map_bks_bk, row_tik, col_tiks, c_map)
          call add_entry_CSR_dist( mesh%M_ddz_bks_bk, row_tik, col_tiks, c_ddz)

        end do ! do jj = 1, single_row_k_nnz

      end do ! do ks = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_map_bks_bk)
    call crop_matrix_CSR_dist( mesh%M_ddz_bks_bk)

    ! Clean up after yourself
    DEallocate( single_row_k_ind)
    DEallocate( single_row_map_val)
    DEallocate( single_row_ddzeta_val)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_matrix_operators_mesh_bks_bk

  subroutine calc_3D_mapping_operator_mesh_bks_ak( mesh)
    ! Calculate mapping operator from the bks-grid to the ak-grid

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_mapping_operator_mesh_bks_ak'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
    integer                                 :: vi
    integer,  dimension(:    ), allocatable :: single_row_vi_ind
    integer                                 :: single_row_vi_nnz
    real(dp), dimension(:    ), allocatable :: single_row_map_b_a_val
    integer                                 :: k
    integer,  dimension(:    ), allocatable :: single_row_k_ind
    integer                                 :: single_row_k_nnz
    real(dp), dimension(:    ), allocatable :: single_row_map_ks_k_val
    integer                                 :: row_vik
    integer                                 :: ii,col_ti,ti,jj,ks,col_tiks
    real(dp)                                :: c_map_b_a, c_map_ks_k
    real(dp)                                :: c_map

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Deallocate existing matrices if necessary
    if (allocateD( mesh%M_map_bks_ak%ptr)) call deallocate_matrix_CSR_dist( mesh%M_map_bks_ak)

    ! Matrix size
    ncols           = mesh%nTri     * (mesh%nz-1) ! from
    ncols_loc       = mesh%nTri_loc * (mesh%nz-1)
    nrows           = mesh%nV       *  mesh%nz    ! to
    nrows_loc       = mesh%nV_loc   *  mesh%nz
    nnz_est_proc    = mesh%M_map_b_a%nnz * mesh%nz * 2

    call allocate_matrix_CSR_dist( mesh%M_map_bks_ak, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for single matrix rows
    allocate( single_row_vi_ind(      mesh%nC_mem*2))
    allocate( single_row_map_b_a_val( mesh%nC_mem*2))

    allocate( single_row_k_ind(        mesh%nz      ))
    allocate( single_row_map_ks_k_val( mesh%nz      ))

    ! Loop over all vertices
    do vi = mesh%vi1, mesh%vi2

      ! Read coefficients from the 2-D gradient operators for this vertex
      call read_single_row_CSR_dist( mesh%M_map_b_a, vi, single_row_vi_ind, single_row_map_b_a_val, single_row_vi_nnz)

      ! Loop over all layers
      do k = 1, mesh%nz

        ! Vertex vi, layer k corresponds to this matrix row
        row_vik = mesh%vik2n( vi,k)

        ! Read coefficients from the zeta gradient operators for this layer
        call read_single_row_CSR_dist( mesh%M_map_ks_k_1D, k, single_row_k_ind, single_row_map_ks_k_val, single_row_k_nnz)

        ! Loop over the entire 3-D local neighbourhood, calculate
        ! coefficients for all 3-D matrix operators

        do ii = 1, single_row_vi_nnz

          col_ti = single_row_vi_ind( ii)
          ti = mesh%n2ti( col_ti)

          ! Coefficients for horizontal gradient matrix operators
          c_map_b_a = single_row_map_b_a_val( ii)

          do jj = 1, single_row_k_nnz

            ks = single_row_k_ind( jj)

            ! Triangle ti, staggered layer ks corresponds to this matrix row
            col_tiks = mesh%tiks2n( ti,ks)

            ! Coefficients for vertical gradient matrix operators
            c_map_ks_k = single_row_map_ks_k_val( jj)

            ! Calculate coefficients
            c_map = c_map_b_a * c_map_ks_k

            ! Add to CSR matrices
            call add_entry_CSR_dist( mesh%M_map_bks_ak, row_vik, col_tiks, c_map)

          end do ! do jj = 1, single_row_k_nnz
        end do ! do ii = 1, single_row_vi_nnz

      end do ! do k = 1, mesh%nz
    end do ! do vi = mesh%vi1, mesh%vi2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_map_bks_ak)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_mapping_operator_mesh_bks_ak

  subroutine calc_3D_mapping_operator_mesh_ak_bks( mesh)
    ! Calculate mapping operator from the ak-grid to the bks-grid

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_mapping_operator_mesh_ak_bks'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
    integer                                 :: ti
    integer,  dimension(:    ), allocatable :: single_row_ti_ind
    integer                                 :: single_row_ti_nnz
    real(dp), dimension(:    ), allocatable :: single_row_map_a_b_val
    integer                                 :: ks
    integer,  dimension(:    ), allocatable :: single_row_ks_ind
    integer                                 :: single_row_ks_nnz
    real(dp), dimension(:    ), allocatable :: single_row_map_k_ks_val
    integer                                 :: row_tiks
    integer                                 :: ii,col_vi,vi,jj,k,col_vik
    real(dp)                                :: c_map_a_b, c_map_k_ks
    real(dp)                                :: c_map

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Deallocate existing matrices if necessary
    if (allocateD( mesh%M_map_ak_bks%ptr)) call deallocate_matrix_CSR_dist( mesh%M_map_ak_bks)

    ! Matrix size
    ncols           = mesh%nV       *  mesh%nz    ! from
    ncols_loc       = mesh%nV_loc   *  mesh%nz
    nrows           = mesh%nTri     * (mesh%nz-1) ! to
    nrows_loc       = mesh%nTri_loc * (mesh%nz-1)
    nnz_est_proc    = mesh%M_map_a_b%nnz * (mesh%nz-1) * 2

    call allocate_matrix_CSR_dist( mesh%M_map_ak_bks, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for single matrix rows
    allocate( single_row_ti_ind(      mesh%nC_mem*2))
    allocate( single_row_map_a_b_val( mesh%nC_mem*2))

    allocate( single_row_ks_ind(       mesh%nz      ))
    allocate( single_row_map_k_ks_val( mesh%nz      ))

    ! Loop over all triangles
    do ti = mesh%ti1, mesh%ti2

      ! Read coefficients from the 2-D gradient operators for this triangle
      call read_single_row_CSR_dist( mesh%M_map_a_b, ti, single_row_ti_ind, single_row_map_a_b_val, single_row_ti_nnz)

      ! Loop over all staggered layers
      do ks = 1, mesh%nz-1

        ! Triangle ti, staggered layer ks corresponds to this matrix row
        row_tiks = mesh%tiks2n( ti,ks)

        ! Read coefficients from the zeta gradient operators for this staggered layer
        call read_single_row_CSR_dist( mesh%M_map_k_ks_1D, ks, single_row_ks_ind, single_row_map_k_ks_val, single_row_ks_nnz)

        ! Loop over the entire 3-D local neighbourhood, calculate
        ! coefficients for all 3-D matrix operators

        do ii = 1, single_row_ti_nnz

          col_vi = single_row_ti_ind( ii)
          vi = mesh%n2vi( col_vi)

          ! Coefficients for horizontal gradient matrix operators
          c_map_a_b = single_row_map_a_b_val( ii)

          do jj = 1, single_row_ks_nnz

            k = single_row_ks_ind( jj)

            ! Vertex vi, layer k corresponds to this matrix row
            col_vik = mesh%vik2n( vi,k)

            ! Coefficients for vertical gradient matrix operators
            c_map_k_ks = single_row_map_k_ks_val( jj)

            ! Calculate coefficients
            c_map = c_map_a_b * c_map_k_ks

            ! Add to CSR matrices
            call add_entry_CSR_dist( mesh%M_map_ak_bks, row_tiks, col_vik, c_map)

          end do ! do jj = 1, single_row_ks_nnz
        end do ! do ii = 1, single_row_ti_nnz

      end do ! do ks = 1, mesh%nz-1
    end do !     do ti = mesh%ti1, mesh%ti2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M_map_ak_bks)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_mapping_operator_mesh_ak_bks

  subroutine calc_3D_matrix_operators_mesh_bk_bk( mesh, ice)
    ! Calculate all 3-D gradient operators in Cartesian coordinates
    !
    ! The basic operators are defined in transformed coordinates [xh, yh, zeta], which
    ! are defined as:
    !
    !   xh   = x
    !   yh   = y
    !   zeta = (Hs - z) / Hi
    !
    ! Applying the chain rule to the gradient operators d/dx, d/dy, d/dz, d2/dx2, d2/dxdy,
    ! d2/dy2, d2/dz2 yields:
    !
    !    d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
    !    d/dy   =  d/dyh    +  dzeta/dy   d/dzeta
    !    d/dz   =              dzeta/dz   d/dzeta
    !   d2/dx2  = d2/dxh2   + d2zeta/dx2  d/dzeta + (dzeta/dx)^2       d2/dzeta2 + 2 dzeta/dx d2/dxhdzeta
    !   d2/dxdy = d2/dxhdyh + d2zeta/dxdy d/dzeta +  dzeta/dx dzeta/dy d2/dzeta2 +   dzeta/dx d2/dyhdzeta + dzeta/dy dxhdzeta
    !   d2/dy2  = d2/dyh2   + d2zeta/dy2  d/dzeta + (dzeta/dy)^2       d2/dzeta2 + 2 dzeta/dy d2/dyhdzeta
    !   d2/dz2  =                                   (dzeta/dz)^2       d2/dzeta2
    !
    ! The d/dxh, d/dyh, d2/dxh2, d2/dxhdyh, d2/dyh2 operators all act in the horizontal plane, so
    ! calculating them requires only information from neighbouring triangles in the same horizontal
    ! layer k. The d/dzeta, d2/dzeta2 operators only act in the vertical column, so they only
    ! require information from the two adjacent horizontal layers k-1, k+1 at the same triangle ti.
    ! Only the two "mixed" gradient operators, d2/dxhdzeta, d2/dyhdzeta, require information from
    ! all neighbouring triangles in both adjacent layers.
    !
    ! Theoretically, we could convert all these basic matrix operators to PETSc format, perform
    ! matrix multiplications on them (so e.g. M_d2dxhdzeta = M_ddxh * M_ddzeta), and then convert
    ! the result back to CSR format. However, this is rather cumbersome to do (especially because
    ! the basic operators are defined in 2-D for the horizontal ones and in 1-D for the vertical,
    ! so we'd need to convert them to act on the 3-D mesh first).
    !
    ! Since the resulting operators  are relatively easy to interpret, we can just calculate
    ! their coefficients directly, which is done here.

    ! In/output variables:
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_matrix_operators_mesh_bk_bk'
    integer                                 :: ncols, nrows, ncols_loc, nrows_loc, nnz_est_proc
    integer                                 :: ti
    integer,  dimension(:    ), allocatable :: single_row_ti_ind
    integer                                 :: single_row_ti_nnz
    real(dp), dimension(:    ), allocatable :: single_row_ddxh_val
    real(dp), dimension(:    ), allocatable :: single_row_ddyh_val
    real(dp), dimension(:    ), allocatable :: single_row_d2dxh2_val
    real(dp), dimension(:    ), allocatable :: single_row_d2dxhdyh_val
    real(dp), dimension(:    ), allocatable :: single_row_d2dyh2_val
    integer                                 :: k
    integer,  dimension(:    ), allocatable :: single_row_k_ind
    integer                                 :: single_row_k_nnz
    real(dp), dimension(:    ), allocatable :: single_row_ddzeta_val
    real(dp), dimension(:    ), allocatable :: single_row_d2dzeta2_val
    integer                                 :: row_tik
    integer                                 :: ii,col_tj,tj,jj,kk,col_tjkk
    real(dp)                                :: dzeta_dx, dzeta_dy, dzeta_dz, d2zeta_dx2, d2zeta_dxdy, d2zeta_dy2
    real(dp)                                :: c_ddxh, c_ddyh, c_d2dxh2, c_d2dxhdyh, c_d2dyh2, c_ddzeta, c_d2dzeta2
    real(dp)                                :: c_ddx, c_ddy, c_ddz, c_d2dx2, c_d2dxdy, c_d2dy2, c_d2dz2

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Deallocate existing matrices if necessary
    if (allocateD( mesh%M2_ddx_bk_bk%ptr   )) call deallocate_matrix_CSR_dist( mesh%M2_ddx_bk_bk)
    if (allocateD( mesh%M2_ddy_bk_bk%ptr   )) call deallocate_matrix_CSR_dist( mesh%M2_ddy_bk_bk)
    if (allocateD( mesh%M2_ddz_bk_bk%ptr   )) call deallocate_matrix_CSR_dist( mesh%M2_ddz_bk_bk)
    if (allocateD( mesh%M2_d2dx2_bk_bk%ptr )) call deallocate_matrix_CSR_dist( mesh%M2_d2dx2_bk_bk)
    if (allocateD( mesh%M2_d2dxdy_bk_bk%ptr)) call deallocate_matrix_CSR_dist( mesh%M2_d2dxdy_bk_bk)
    if (allocateD( mesh%M2_d2dy2_bk_bk%ptr )) call deallocate_matrix_CSR_dist( mesh%M2_d2dy2_bk_bk)
    if (allocateD( mesh%M2_d2dz2_bk_bk%ptr )) call deallocate_matrix_CSR_dist( mesh%M2_d2dz2_bk_bk)

    ! Matrix size
    ncols           = mesh%nTri     * mesh%nz ! from
    ncols_loc       = mesh%nTri_loc * mesh%nz
    nrows           = mesh%nTri     * mesh%nz ! to
    nrows_loc       = mesh%nTri_loc * mesh%nz
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz * mesh%nz * 3

    call allocate_matrix_CSR_dist( mesh%M2_ddx_bk_bk   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_ddy_bk_bk   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_ddz_bk_bk   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_d2dx2_bk_bk , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_d2dxdy_bk_bk, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_d2dy2_bk_bk , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( mesh%M2_d2dz2_bk_bk , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for single matrix rows
    allocate( single_row_ti_ind(       mesh%nC_mem*2))
    allocate( single_row_ddxh_val(     mesh%nC_mem*2))
    allocate( single_row_ddyh_val(     mesh%nC_mem*2))
    allocate( single_row_d2dxh2_val(   mesh%nC_mem*2))
    allocate( single_row_d2dxhdyh_val( mesh%nC_mem*2))
    allocate( single_row_d2dyh2_val(   mesh%nC_mem*2))

    allocate( single_row_k_ind(        mesh%nz      ))
    allocate( single_row_ddzeta_val(   mesh%nz      ))
    allocate( single_row_d2dzeta2_val( mesh%nz      ))

    ! Loop over all triangles
    do ti = mesh%ti1, mesh%ti2

      ! Read coefficients from the 2-D gradient operators for this triangle
      call read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ti_ind, single_row_ddxh_val    , single_row_ti_nnz)
      call read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ti_ind, single_row_ddyh_val    , single_row_ti_nnz)
      call read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ti_ind, single_row_d2dxh2_val  , single_row_ti_nnz)
      call read_single_row_CSR_dist( mesh%M2_d2dxdy_b_b, ti, single_row_ti_ind, single_row_d2dxhdyh_val, single_row_ti_nnz)
      call read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ti_ind, single_row_d2dyh2_val  , single_row_ti_nnz)

      ! Loop over all layers
      do k = 1, mesh%nz

        ! Triangle ti, layer k corresponds to this matrix row
        row_tik = mesh%tik2n( ti,k)

        ! Read coefficients from the zeta gradient operators for this layer
        call read_single_row_CSR_dist( mesh%M_ddzeta_k_k_1D  , k, single_row_k_ind, single_row_ddzeta_val  , single_row_k_nnz)
        call read_single_row_CSR_dist( mesh%M_d2dzeta2_k_k_1D, k, single_row_k_ind, single_row_d2dzeta2_val, single_row_k_nnz)

        ! Gradients of zeta at triangle ti, layer k
        dzeta_dx    = ice%dzeta_dx_bk(    ti,k)
        dzeta_dy    = ice%dzeta_dy_bk(    ti,k)
        dzeta_dz    = ice%dzeta_dz_bk(    ti,k)
        d2zeta_dx2  = ice%d2zeta_dx2_bk(  ti,k)
        d2zeta_dxdy = ice%d2zeta_dxdy_bk( ti,k)
        d2zeta_dy2  = ice%d2zeta_dy2_bk(  ti,k)

        ! Loop over the entire 3-D local neighbourhood, calculate
        ! coefficients for all 3-D matrix operators

        do ii = 1, single_row_ti_nnz

          col_tj = single_row_ti_ind( ii)
          tj = mesh%n2ti( col_tj)

          ! Coefficients for horizontal gradient matrix operators
          c_ddxh     = single_row_ddxh_val(     ii)
          c_ddyh     = single_row_ddyh_val(     ii)
          c_d2dxh2   = single_row_d2dxh2_val(   ii)
          c_d2dxhdyh = single_row_d2dxhdyh_val( ii)
          c_d2dyh2   = single_row_d2dyh2_val(   ii)

          do jj = 1, single_row_k_nnz

            kk = single_row_k_ind( jj)

            ! Triangle tj, layer kk corresponds to this matrix row
            col_tjkk = mesh%tik2n( tj,kk)

            ! Coefficients for vertical gradient matrix operators
            c_ddzeta   = single_row_ddzeta_val(   jj)
            c_d2dzeta2 = single_row_d2dzeta2_val( jj)

            ! Calculate coefficients
            c_ddx    = 0._dp
            c_ddy    = 0._dp
            c_ddz    = 0._dp
            c_d2dx2  = 0._dp
            c_d2dxdy = 0._dp
            c_d2dy2  = 0._dp
            c_d2dz2  = 0._dp

            ! Horizontal-only part
            if (kk == k) then
              c_ddx    = c_ddx    + c_ddxh                                                            ! Now:  d/dx   =  d/dxh ...
              c_ddy    = c_ddy    + c_ddyh                                                            ! Now:  d/dy   =  d/dyh ...
              c_d2dx2  = c_d2dx2  + c_d2dxh2                                                          ! Now: d2/dx2  = d2/dxh2 ...
              c_d2dxdy = c_d2dxdy + c_d2dxhdyh                                                        ! Now: d2/dxdy = d2/dxhdyh ...
              c_d2dy2  = c_d2dy2  + c_d2dyh2                                                          ! Now: d2/dy2  = d2/dyh2 ...
            end if ! if (kk == k) then

            ! Vertical-only part
            if (tj == ti) then
              c_ddx    = c_ddx    + dzeta_dx    * c_ddzeta                                            ! Now:  d/dx   =  d/dxh    +  dzeta/dx   d/dzeta
              c_ddy    = c_ddy    + dzeta_dy    * c_ddzeta                                            ! Now:  d/dy   =  d/dyh    +  dzeta/dy   d/dzeta
              c_ddz    = c_ddz    + dzeta_dz    * c_ddzeta                                            ! Now:  d/dz   =              dzeta/dz   d/dzeta
              c_d2dx2  = c_d2dx2  + d2zeta_dx2  * c_ddzeta + dzeta_dx * dzeta_dx * c_d2dzeta2         ! Now: d2/dx2  = d2/dxh2   + d2zeta/dx2  d/dzeta + (dzeta/dx)^2       d2/dzeta2 + ...
              c_d2dxdy = c_d2dxdy + d2zeta_dxdy * c_ddzeta + dzeta_dx * dzeta_dy * c_d2dzeta2         ! Now: d2/dxdy = d2/dxhdyh + d2zeta/dxdy d/dzeta +  dzeta/dx dzeta/dy d2/dzeta2 + ...
              c_d2dy2  = c_d2dy2  + d2zeta_dy2  * c_ddzeta + dzeta_dy * dzeta_dy * c_d2dzeta2         ! Now: d2/dy2  = d2/dyh2   + d2zeta/dy2  d/dzeta + (dzeta/dy)^2       d2/dzeta2 + ...
              c_d2dz2  = c_d2dz2  + dzeta_dz**2 * c_d2dzeta2                                          ! Now: d2/dz2  =                                   (dzeta/dz)^2       d2/dzeta2
            end if ! if (tj == ti) then

            ! Mixed part
            c_d2dx2  = c_d2dx2  + 2._dp * dzeta_dx * c_ddxh * c_ddzeta                                ! Now: d2/dx2  = d2/dxh2   + d2zeta/dx2  d/dzeta + (dzeta/dx)^2       d2/dzeta2 + 2 dzeta/dx d2/dxhdzeta
            c_d2dxdy = c_d2dxdy +         dzeta_dx * c_ddyh * c_ddzeta + dzeta_dy * c_ddxh * c_ddzeta ! Now: d2/dxdy = d2/dxhdyh + d2zeta/dxdy d/dzeta +  dzeta/dx dzeta/dy d2/dzeta2 +   dzeta/dx d2/dyhdzeta + dzeta/dy dxhdzeta
            c_d2dy2  = c_d2dy2  + 2._dp * dzeta_dy * c_ddyh * c_ddzeta                                ! Now: d2/dy2  = d2/dyh2   + d2zeta/dy2  d/dzeta + (dzeta/dy)^2       d2/dzeta2 + 2 dzeta/dy d2/dyhdzeta

            ! Add to CSR matrices
            call add_entry_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, col_tjkk, c_ddx   )
            call add_entry_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, col_tjkk, c_ddy   )
            call add_entry_CSR_dist( mesh%M2_ddz_bk_bk   , row_tik, col_tjkk, c_ddz   )
            call add_entry_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, col_tjkk, c_d2dx2 )
            call add_entry_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, col_tjkk, c_d2dxdy)
            call add_entry_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, col_tjkk, c_d2dy2 )
            call add_entry_CSR_dist( mesh%M2_d2dz2_bk_bk , row_tik, col_tjkk, c_d2dz2 )

          end do ! do jj = 1, single_row_k_nnz
        end do ! do ii = 1, single_row_ti_nnz

      end do ! do k = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Crop matrix memory
    call crop_matrix_CSR_dist( mesh%M2_ddx_bk_bk)
    call crop_matrix_CSR_dist( mesh%M2_ddy_bk_bk)
    call crop_matrix_CSR_dist( mesh%M2_ddz_bk_bk)
    call crop_matrix_CSR_dist( mesh%M2_d2dx2_bk_bk)
    call crop_matrix_CSR_dist( mesh%M2_d2dxdy_bk_bk)
    call crop_matrix_CSR_dist( mesh%M2_d2dy2_bk_bk)
    call crop_matrix_CSR_dist( mesh%M2_d2dz2_bk_bk)

    ! Clean up after yourself
    DEallocate( single_row_ti_ind)
    DEallocate( single_row_ddxh_val)
    DEallocate( single_row_ddyh_val)
    DEallocate( single_row_d2dxh2_val)
    DEallocate( single_row_d2dxhdyh_val)
    DEallocate( single_row_d2dyh2_val)
    DEallocate( single_row_k_ind)
    DEallocate( single_row_ddzeta_val)
    DEallocate( single_row_d2dzeta2_val)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_matrix_operators_mesh_bk_bk

  ! Apply them to calculate actual gradients

  subroutine calc_3D_gradient_bk_ak( mesh, AA, d_bk, grad_d_ak)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: d_bk
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz), intent(out  ) :: grad_d_ak

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bk_ak'
    real(dp), dimension(:,:  ), allocatable :: d_bk_tot
    integer                                 :: vi,k,row_vik,ii1,ii2,ii,col_tikk,ti,kk
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nV_loc   * mesh%nz .or. &
        AA%m     /= mesh%nV       * mesh%nz .or. &
        AA%n_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%n     /= mesh%nTri     * mesh%nz .or. &
        size(      d_bk,1) /= mesh%nTri_loc .or. size(      d_bk,2) /= mesh%nz .or. &
        size( grad_d_ak,1) /= mesh%nV_loc   .or. size( grad_d_ak,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match!')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bk_tot( mesh%nTri, mesh%nz))

    ! Gather data
    call gather_to_all_dp_2D( d_bk, d_bk_tot)

    ! Calculate gradient
    do vi = mesh%vi1, mesh%vi2
    do k  = 1, mesh%nz

      ! Vertex vi, layer k corresponds to this matrix row
      row_vik = mesh%vik2n( vi,k)

      ! Initialise
      grad_d_ak( vi,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_vik)
      ii2 = AA%ptr( row_vik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tikk = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle ti, layer kk
        ti = mesh%n2tik( col_tikk,1)
        kk = mesh%n2tik( col_tikk,2)

        ! Add contribution
        grad_d_ak( vi,k) = grad_d_ak( vi,k) + cAA * d_bk_tot( ti,kk)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bk_ak

  subroutine calc_3D_gradient_ak_bk( mesh, AA, d_ak, grad_d_bk)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz), intent(in   ) :: d_ak
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(out  ) :: grad_d_bk

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_ak_bk'
    real(dp), dimension(:,:  ), allocatable :: d_ak_tot
    integer                                 :: ti,k,row_tik,ii1,ii2,ii,col_vikk,vi,kk
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%m     /= mesh%nTri     * mesh%nz .or. &
        AA%n_loc /= mesh%nV_loc   * mesh%nz .or. &
        AA%n     /= mesh%nV       * mesh%nz .or. &
        size(      d_ak,1) /= mesh%nV_loc   .or. size(      d_ak,2) /= mesh%nz .or. &
        size( grad_d_bk,1) /= mesh%nTri_loc .or. size( grad_d_bk,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match!')
    end if

    ! allocate memory for gathered vector x
    allocate( d_ak_tot( mesh%nV, mesh%nz))

    ! Gather data
    call gather_to_all_dp_2D( d_ak, d_ak_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Initialise
      grad_d_bk( ti,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tik)
      ii2 = AA%ptr( row_tik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_vikk = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to vertex vi, layer kk
        vi = mesh%n2vik( col_vikk,1)
        kk = mesh%n2vik( col_vikk,2)

        ! Add contribution
        grad_d_bk( ti,k) = grad_d_bk( ti,k) + cAA * d_ak_tot( vi,kk)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_ak_bk

  subroutine calc_3D_gradient_bk_bks( mesh, AA, d_bk, grad_d_bks)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                  intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                  intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz),   intent(in   ) :: d_bk
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1), intent(out  ) :: grad_d_bks

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bk_bks'
    real(dp), dimension(:,:  ), allocatable :: d_bk_tot
    integer                                 :: ti,ks,row_tiks,ii1,ii2,ii,col_tjk,tj,k
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%m     /= mesh%nTri     * (mesh%nz-1) .or. &
        AA%n_loc /= mesh%nTri_loc *  mesh%nz    .or. &
        AA%n     /= mesh%nTri     *  mesh%nz    .or. &
        size(      d_bk ,1) /= mesh%nTri_loc .or. size(      d_bk ,2) /= mesh%nz .or. &
        size( grad_d_bks,1) /= mesh%nTri_loc .or. size( grad_d_bks,2) /= mesh%nz-1) then
      call crash('matrix and vector sizes dont match!')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bk_tot( mesh%nTri, mesh%nz))

    ! Gather data
    call gather_to_all_dp_2D( d_bk, d_bk_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do ks = 1, mesh%nz-1

      ! Triangle ti, staggered layer ks corresponds to this matrix row
      row_tiks = mesh%tiks2n( ti,ks)

      ! Initialise
      grad_d_bks( ti,ks) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tiks)
      ii2 = AA%ptr( row_tiks+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tjk = AA%ind( ii)
        cAA     = AA%val( ii)

        ! This matrix column corresponds to triangle tj, layer k
        tj = mesh%n2tik( col_tjk,1)
        k  = mesh%n2tik( col_tjk,2)

        ! Add contribution
        grad_d_bks( ti,ks) = grad_d_bks( ti,ks) + cAA * d_bk_tot( tj,k)

      end do ! do ii = ii1, ii2

    end do ! do ks  = 1, mesh%nz-1
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bk_bks

  subroutine calc_3D_gradient_bks_bk( mesh, AA, d_bks, grad_d_bk)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                  intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                  intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1), intent(in   ) :: d_bks
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz),   intent(out  ) :: grad_d_bk

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bks_bk'
    real(dp), dimension(:,:  ), allocatable :: d_bks_tot
    integer                                 :: ti,k,row_tik,ii1,ii2,ii,col_tjks,tj,ks
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc *  mesh%nz    .or. &
        AA%m     /= mesh%nTri     *  mesh%nz    .or. &
        AA%n_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%n     /= mesh%nTri     * (mesh%nz-1) .or. &
        size(      d_bks,1) /= mesh%nTri_loc .or. size(      d_bks,2) /= mesh%nz-1 .or. &
        size( grad_d_bk ,1) /= mesh%nTri_loc .or. size( grad_d_bk ,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match!')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bks_tot( mesh%nTri, mesh%nz-1))

    ! Gather data
    call gather_to_all_dp_2D( d_bks, d_bks_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do k = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Initialise
      grad_d_bk( ti,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tik)
      ii2 = AA%ptr( row_tik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tjks = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle tj, staggered layer ks
        tj = mesh%n2tiks( col_tjks,1)
        ks = mesh%n2tiks( col_tjks,2)

        ! Add contribution
        grad_d_bk( ti,k) = grad_d_bk( ti,k) + cAA * d_bks_tot( tj,ks)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bks_bk

  subroutine map_bks_ak( mesh, AA, d_bks, grad_d_ak)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                     intent(in   )    :: mesh
    type(type_sparse_matrix_CSR_dp),     intent(in   )    :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1),        intent(in   )    :: d_bks
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz),          intent(out  )    :: grad_d_ak

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'map_bks_ak'
    real(dp), dimension(:,:  ), allocatable            :: d_bks_tot
    integer                                            :: vi,k,row_vik,ii1,ii2,ii,col_tiks,ti,ks
    real(dp)                                           :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nV_loc   *  mesh%nz    .or. &
        AA%m     /= mesh%nV       *  mesh%nz    .or. &
        AA%n_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%n     /= mesh%nTri     * (mesh%nz-1) .or. &
        size(      d_bks,1) /= mesh%nTri_loc .or. size(      d_bks,2) /= mesh%nz-1 .or. &
        size( grad_d_ak ,1) /= mesh%nV_loc   .or. size( grad_d_ak ,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match!')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bks_tot( mesh%nTri, mesh%nz-1))

    ! Gather data
    call gather_to_all_dp_2D( d_bks, d_bks_tot)

    ! Calculate gradient
    do vi = mesh%vi1, mesh%vi2
    do k = 1, mesh%nz

      ! Vertex vi, layer k corresponds to this matrix row
      row_vik = mesh%vik2n( vi,k)

      ! Initialise
      grad_d_ak( vi,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_vik)
      ii2 = AA%ptr( row_vik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tiks = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle ti, staggered layer ks
        ti = mesh%n2tiks( col_tiks,1)
        ks = mesh%n2tiks( col_tiks,2)

        ! Add contribution
        grad_d_ak( vi,k) = grad_d_ak( vi,k) + cAA * d_bks_tot( ti,ks)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_bks_ak

  subroutine map_ak_bks( mesh, AA, d_ak, grad_d_bks)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                  intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                  intent(in   ) :: AA
    real(dp), dimension(mesh%vi1:mesh%vi2,mesh%nz),   intent(in   ) :: d_ak
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz-1), intent(out  ) :: grad_d_bks

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'map_ak_bks'
    real(dp), dimension(:,:  ), allocatable :: d_ak_tot
    integer                                 :: ti,ks,row_tiks,ii1,ii2,ii,col_vik,vi,k
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * (mesh%nz-1) .or. &
        AA%m     /= mesh%nTri     * (mesh%nz-1) .or. &
        AA%n_loc /= mesh%nV_loc   *  mesh%nz    .or. &
        AA%n     /= mesh%nV       *  mesh%nz    .or. &
        size(      d_ak ,1) /= mesh%nV_loc   .or. size(      d_ak ,2) /= mesh%nz .or. &
        size( grad_d_bks,1) /= mesh%nTri_loc .or. size( grad_d_bks,2) /= mesh%nz-1) then
      call crash('matrix and vector sizes dont match!')
    end if

    ! allocate memory for gathered vector x
    allocate( d_ak_tot( mesh%nV, mesh%nz))

    ! Gather data
    call gather_to_all_dp_2D( d_ak, d_ak_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do ks = 1, mesh%nz-1

      ! Triangle ti, staggered layer ks corresponds to this matrix row
      row_tiks = mesh%tiks2n( ti,ks)

      ! Initialise
      grad_d_bks( ti,ks) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tiks)
      ii2 = AA%ptr( row_tiks+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_vik = AA%ind( ii)
        cAA     = AA%val( ii)

        ! This matrix column corresponds to vertex vi, layer k
        vi = mesh%n2vik( col_vik,1)
        k  = mesh%n2vik( col_vik,2)

        ! Add contribution
        grad_d_bks( ti,ks) = grad_d_bks( ti,ks) + cAA * d_ak_tot( vi,k)

      end do ! do ii = ii1, ii2

    end do ! do ks = 1, mesh%nz-1
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_ak_bks

  subroutine calc_3D_gradient_bk_bk( mesh, AA, d_bk, grad_d_bk)
    ! Apply a 3-D gradient operator to a 3-D data field

    ! In- and output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: d_bk
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(out  ) :: grad_d_bk

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_3D_gradient_bk_bk'
    real(dp), dimension(:,:  ), allocatable :: d_bk_tot
    integer                                 :: ti,k,row_tik,ii1,ii2,ii,col_tjkk,tj,kk
    real(dp)                                :: cAA

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety: check sizes
    if (AA%m_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%m     /= mesh%nTri     * mesh%nz .or. &
        AA%n_loc /= mesh%nTri_loc * mesh%nz .or. &
        AA%n     /= mesh%nTri     * mesh%nz .or. &
        size(      d_bk,1) /= mesh%nTri_loc .or. size(      d_bk,2) /= mesh%nz .or. &
        size( grad_d_bk,1) /= mesh%nTri_loc .or. size( grad_d_bk,2) /= mesh%nz) then
      call crash('matrix and vector sizes dont match!')
    end if

    ! allocate memory for gathered vector x
    allocate( d_bk_tot( mesh%nTri, mesh%nz))

    ! Gather data
    call gather_to_all_dp_2D( d_bk, d_bk_tot)

    ! Calculate gradient
    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! Triangle ti, layer k corresponds to this matrix row
      row_tik = mesh%tik2n( ti,k)

      ! Initialise
      grad_d_bk( ti,k) = 0._dp

      ! Loop over all contributing 3-D neighbours
      ii1 = AA%ptr( row_tik)
      ii2 = AA%ptr( row_tik+1) - 1

      do ii = ii1, ii2

        ! Read matrix coefficient
        col_tjkk = AA%ind( ii)
        cAA      = AA%val( ii)

        ! This matrix column corresponds to triangle tj, layer kk
        tj = mesh%n2tik( col_tjkk,1)
        kk = mesh%n2tik( col_tjkk,2)

        ! Add contribution
        grad_d_bk( ti,k) = grad_d_bk( ti,k) + cAA * d_bk_tot( tj,kk)

      end do ! do ii = ii1, ii2

    end do ! do k  = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_gradient_bk_bk

end module mesh_operators
