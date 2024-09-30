module create_maps_grid_mesh

  ! Create remapping objects between a square Cartesian grid and a mesh.

#include <petsc/finclude/petscksp.h>
  use petscksp
  use mpi_basic, only: par, sync
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_utilities, only: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, &
    deallocate_matrix_CSR_dist, add_entry_CSR_dist, add_empty_row_CSR_dist
  use remapping_types, only: type_single_row_mapping_matrices, type_map
  use math_utilities, only: is_in_triangle
  use mesh_utilities, only: find_containing_triangle, calc_Voronoi_cell, is_in_Voronoi_cell
  use petsc_basic, only: mat_CSR2petsc
  use line_tracing_grid, only: trace_line_grid
  use line_tracing_triangles, only: trace_line_tri
  use line_tracing_Voronoi, only: trace_line_Vor
  use grid_basic, only: calc_matrix_operators_grid

  implicit none

  private

  public :: create_map_from_xy_grid_to_mesh
  public :: create_map_from_xy_grid_to_mesh_triangles
  public :: create_map_from_mesh_to_xy_grid

contains

  subroutine create_map_from_xy_grid_to_mesh( grid, mesh, map)
    ! Create a new mapping object from an x/y-grid to a mesh.
    !
    ! By default uses 2nd-order conservative remapping.
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as having an area smaller
    !       than ten times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_xy_grid_to_mesh'
    type(PetscErrorCode)                               :: perr
    logical                                            :: count_coincidences
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR
    integer,  dimension(:    ), allocatable            :: mask_do_simple_average
    integer                                            :: vi
    real(dp), dimension( mesh%nC_mem,2)                :: Vor
    integer,  dimension( mesh%nC_mem  )                :: Vor_vi
    integer,  dimension( mesh%nC_mem  )                :: Vor_ti
    integer                                            :: nVor
    integer                                            :: vori1, vori2
    real(dp), dimension(2)                             :: p, q
    integer                                            :: k, i, j, kk, vj
    real(dp)                                           :: xl, xu, yl, yu
    real(dp), dimension(2)                             :: sw, se, nw, ne
    integer                                            :: vi_hint
    real(dp)                                           :: xmin, xmax, ymin, ymax
    integer                                            :: il, iu, jl, ju
    type(type_single_row_mapping_matrices)             :: single_row_Vor, single_row_grid
    type(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    type(tMat)                                         :: w0    , w1x    , w1y
    integer                                            :: row, k1, k2, col
    real(dp)                                           :: A_overlap_tot
    type(tMat)                                         :: grid_M_ddx, grid_M_ddy
    type(tMat)                                         :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = grid%name
    map%name_dst  = mesh%name
    map%method    = '2nd_order_conservative'

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix size
    nrows           = mesh%nV  ! to
    nrows_loc       = mesh%nV_loc
    ncols           = grid%n   ! from
    ncols_loc       = grid%n_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%A) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%A))) )

    call allocate_matrix_CSR_dist( A_xdy_a_g_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_a_g_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_a_g_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for single row results
    single_row_Vor%n_max = nnz_per_row_max
    single_row_Vor%n     = 0
    allocate( single_row_Vor%index_left( single_row_Vor%n_max))
    allocate( single_row_Vor%LI_xdy(     single_row_Vor%n_max))
    allocate( single_row_Vor%LI_mxydx(   single_row_Vor%n_max))
    allocate( single_row_Vor%LI_xydy(    single_row_Vor%n_max))

    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    allocate( single_row_grid%index_left( single_row_grid%n_max))
    allocate( single_row_grid%LI_xdy(     single_row_grid%n_max))
    allocate( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    allocate( single_row_grid%LI_xydy(    single_row_grid%n_max))

    allocate( mask_do_simple_average( mesh%nV), source = 0)

    ! Calculate line integrals around all Voronoi cells
    do row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      if (mesh%A( vi) < 10._dp * grid%dx**2) then
        ! This Voronoi cell is small enough to warrant a proper line integral

        mask_do_simple_average( vi) = 0

        ! Clean up single row results
        single_row_Vor%n          = 0
        single_row_Vor%index_left = 0
        single_row_Vor%LI_xdy     = 0._dp
        single_row_Vor%LI_mxydx   = 0._dp
        single_row_Vor%LI_xydy    = 0._dp

        ! Integrate around the complete Voronoi cell boundary
        call calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)
        do vori1 = 1, nVor
          vori2 = vori1 + 1
          if (vori2 > nVor) vori2 = 1
          p = Vor( vori1,:)
          q = Vor( vori2,:)
          count_coincidences = .true.
          call trace_line_grid( grid, p, q, single_row_Vor, count_coincidences)
        end do

        ! Safety
        if (single_row_Vor%n == 0) call crash('couldnt find any grid cells overlapping with this small Voronoi cell!')

        ! Next integrate around the grid cells overlapping with this Voronoi cell
        do k = 1, single_row_Vor%n

          ! Clean up single row results
          single_row_grid%n          = 0
          single_row_grid%index_left = 0
          single_row_grid%LI_xdy     = 0._dp
          single_row_grid%LI_mxydx   = 0._dp
          single_row_grid%LI_xydy    = 0._dp

          ! The grid cell
          col = single_row_Vor%index_left( k)
          i   = grid%n2ij( col,1)
          j   = grid%n2ij( col,2)

          xl = grid%x( i) - grid%dx / 2._dp
          xu = grid%x( i) + grid%dx / 2._dp
          yl = grid%y( j) - grid%dx / 2._dp
          yu = grid%y( j) + grid%dx / 2._dp

          sw = [xl,yl]
          nw = [xl,yu]
          se = [xu,yl]
          ne = [xu,yu]

          ! Integrate around the grid cell
          vi_hint = vi
          count_coincidences = .false.
          call trace_line_Vor( mesh, sw, se, single_row_grid, count_coincidences, vi_hint)
          call trace_line_Vor( mesh, se, ne, single_row_grid, count_coincidences, vi_hint)
          call trace_line_Vor( mesh, ne, nw, single_row_grid, count_coincidences, vi_hint)
          call trace_line_Vor( mesh, nw, sw, single_row_grid, count_coincidences, vi_hint)

          ! Safety
          if (single_row_grid%n == 0) call crash('couldnt find any Voronoi cells overlapping with this grid cell!')

          ! Add contribution for this particular triangle
          do kk = 1, single_row_grid%n
            vj = single_row_grid%index_left( kk)
            if (vj == vi) then
              ! Add contribution to this triangle
              single_row_Vor%LI_xdy(   k) = single_row_Vor%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
              single_row_Vor%LI_mxydx( k) = single_row_Vor%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
              single_row_Vor%LI_xydy(  k) = single_row_Vor%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
              exit
            end if
          end do ! do kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          call add_entry_CSR_dist( A_xdy_a_g_CSR  , row, col, single_row_Vor%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_a_g_CSR, row, col, single_row_Vor%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_a_g_CSR , row, col, single_row_Vor%LI_xydy(  k))

        end do ! do k = 1, single_row_Vor%n

      else ! if (mesh%A( vi) < 10._dp * grid%dx**2) then
        ! This Voronoi cell is big enough that we can just average over the grid cells it contains

        mask_do_simple_average( vi) = 1

        ! Clean up single row results
        single_row_Vor%n = 0

        ! Find the square of grid cells enveloping this Voronoi cell
        call calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)

        xmin = minval( Vor( 1:nVor,1))
        xmax = maxval( Vor( 1:nVor,1))
        ymin = minval( Vor( 1:nVor,2))
        ymax = maxval( Vor( 1:nVor,2))

        il = max( 1, min( grid%nx, 1 + floor( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        iu = max( 1, min( grid%nx, 1 + floor( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        jl = max( 1, min( grid%ny, 1 + floor( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx) ))
        ju = max( 1, min( grid%ny, 1 + floor( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

        ! Check which of the grid cells in this square lie inside the triangle
        do i = il, iu
        do j = jl, ju

          col = grid%ij2n( i,j)
          p   = [grid%x( i), grid%y( j)]

          if (is_in_Voronoi_cell( mesh, p, vi)) then
            ! This grid cell lies inside the triangle; add it to the single row
            single_row_Vor%n = single_row_Vor%n + 1
            single_row_Vor%index_left( single_row_Vor%n) = col
            single_row_Vor%LI_xdy(     single_row_Vor%n) = grid%dx**2
            single_row_Vor%LI_mxydx(   single_row_Vor%n) = grid%x( i) * grid%dx**2
            single_row_Vor%LI_xydy(    single_row_Vor%n) = grid%y( j) * grid%dx**2
          end if

        end do
        end do

        ! Safety
        if (single_row_Vor%n == 0) call crash('couldnt find any grid cells overlapping with this big Voronoi cell!')

        ! Add entries to the big matrices
        do k = 1, single_row_Vor%n
          col = single_row_Vor%index_left( k)
          call add_entry_CSR_dist( A_xdy_a_g_CSR  , vi, col, single_row_Vor%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_a_g_CSR, vi, col, single_row_Vor%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_a_g_CSR , vi, col, single_row_Vor%LI_xydy(  k))
        end do

      end if ! if (mesh%A( vi) < 4._dp * grid%dx**2) then

    end do ! do vi = mesh%vi1, mesh%vi2
    call sync

    ! Clean up after yourself
    deallocate( single_row_Vor%index_left )
    deallocate( single_row_Vor%LI_xdy     )
    deallocate( single_row_Vor%LI_mxydx   )
    deallocate( single_row_Vor%LI_xydy    )

    deallocate( single_row_grid%index_left )
    deallocate( single_row_grid%LI_xdy     )
    deallocate( single_row_grid%LI_mxydx   )
    deallocate( single_row_grid%LI_xydy    )

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    ! ==============================================================

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      k1 = A_xdy_a_g_CSR%ptr( row  )
      k2 = A_xdy_a_g_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_a_g_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_a_g_CSR%ind( k)
        call add_entry_CSR_dist( w0_CSR, row, col, A_xdy_a_g_CSR%val( k) / A_overlap_tot)
      end do

      if (mask_do_simple_average( vi) == 0) then
        ! For small vertices, include the gradient terms

        do k = k1, k2
          col = A_xdy_a_g_CSR%ind( k)
          ! Grid cell
          i = grid%n2ij( col,1)
          j = grid%n2ij( col,2)
          call add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_a_g_CSR%val( k) / A_overlap_tot) - (grid%x( i) * w0_CSR%val( k)))
          call add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_a_g_CSR%val(  k) / A_overlap_tot) - (grid%y( j) * w0_CSR%val( k)))
        end do

      else
        ! For large vertices, don't include the gradient terms

        call add_empty_row_CSR_dist( w1x_CSR, row)
        call add_empty_row_CSR_dist( w1y_CSR, row)

      end if ! if (mask_do_simple_average( vi) == 0) then

    end do ! do row = mesh%vi1, mesh%vi2
    call sync

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( A_xdy_a_g_CSR  )
    call deallocate_matrix_CSR_dist( A_mxydx_a_g_CSR)
    call deallocate_matrix_CSR_dist( A_xydy_a_g_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix

    call calc_matrix_operators_grid( grid, grid_M_ddx, grid_M_ddy)

    call MatDuplicate( w0, MAT_COPY_VALUES, map%M, perr)
    call MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatDestroy( grid_M_ddx    , perr)
    call MatDestroy( grid_M_ddy    , perr)
    call MatDestroy( w0            , perr)
    call MatDestroy( w1x           , perr)
    call MatDestroy( w1y           , perr)

    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    call MatDestroy( M1, perr)
    call MatDestroy( M2, perr)
    deallocate( mask_do_simple_average)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_xy_grid_to_mesh

  subroutine create_map_from_xy_grid_to_mesh_triangles( grid, mesh, map)
    ! Create a new mapping object from an x/y-grid to the triangles of a mesh.
    !
    ! By default uses 2nd-order conservative remapping.
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as having an area smaller
    !       than ten times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_xy_grid_to_mesh_triangles'
    type(PetscErrorCode)                               :: perr
    logical                                            :: count_coincidences
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR
    integer,  dimension(:    ), allocatable            :: mask_do_simple_average
    integer                                            :: ti
    integer                                            :: n1, n2, via, vib, vic
    real(dp), dimension(2)                             :: p, q
    integer                                            :: k, i, j, kk, tj
    real(dp)                                           :: xl, xu, yl, yu
    real(dp), dimension(2)                             :: sw, se, nw, ne
    integer                                            :: ti_hint
    real(dp)                                           :: xmin, xmax, ymin, ymax
    integer                                            :: il, iu, jl, ju
    type(type_single_row_mapping_matrices)             :: single_row_tri, single_row_grid
    real(dp), dimension(2)                             :: pa, pb, pc
    type(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    type(tMat)                                         :: w0    , w1x    , w1y
    integer                                            :: row, k1, k2, col
    real(dp)                                           :: A_overlap_tot
    type(tMat)                                         :: grid_M_ddx, grid_M_ddy
    type(tMat)                                         :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = grid%name
    map%name_dst  = trim( mesh%name) // '_triangles'
    map%method    = '2nd_order_conservative'

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix size
    nrows           = mesh%nTri  ! to
    nrows_loc       = mesh%nTri_loc
    ncols           = grid%n     ! from
    ncols_loc       = grid%n_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%TriA) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%TriA))) )

    call allocate_matrix_CSR_dist( A_xdy_b_g_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_b_g_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_b_g_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for single row results
    single_row_tri%n_max = nnz_per_row_max
    single_row_tri%n     = 0
    allocate( single_row_tri%index_left( single_row_tri%n_max))
    allocate( single_row_tri%LI_xdy(     single_row_tri%n_max))
    allocate( single_row_tri%LI_mxydx(   single_row_tri%n_max))
    allocate( single_row_tri%LI_xydy(    single_row_tri%n_max))

    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    allocate( single_row_grid%index_left( single_row_grid%n_max))
    allocate( single_row_grid%LI_xdy(     single_row_grid%n_max))
    allocate( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    allocate( single_row_grid%LI_xydy(    single_row_grid%n_max))

    allocate( mask_do_simple_average( mesh%nTri), source = 0)

    ! Calculate line integrals around all triangles
    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      if (mesh%TriA( ti) < 10._dp * grid%dx**2) then
        ! This triangle is small enough to warrant a proper line integral

        mask_do_simple_average( ti) = 0

        ! Clean up single row results
        single_row_tri%n          = 0
        single_row_tri%index_left = 0
        single_row_tri%LI_xdy     = 0._dp
        single_row_tri%LI_mxydx   = 0._dp
        single_row_tri%LI_xydy    = 0._dp

        ! Integrate around the triangle
        do n1 = 1, 3
          n2 = n1 + 1
          if (n2 == 4) n2 = 1
          via = mesh%Tri( ti,n1)
          vib = mesh%Tri( ti,n2)
          p = mesh%V( via,:)
          q = mesh%V( vib,:)
          count_coincidences = .true.
          call trace_line_grid( grid, p, q, single_row_tri, count_coincidences)
        end do

        ! Safety
        if (single_row_tri%n == 0) call crash('couldnt find any grid cells overlapping with this small triangle!')

        ! Next integrate around the grid cells overlapping with this triangle
        do k = 1, single_row_tri%n

          ! Clean up single row results
          single_row_grid%n          = 0
          single_row_grid%index_left = 0
          single_row_grid%LI_xdy     = 0._dp
          single_row_grid%LI_mxydx   = 0._dp
          single_row_grid%LI_xydy    = 0._dp

          ! The grid cell
          col = single_row_tri%index_left( k)
          i   = grid%n2ij( col,1)
          j   = grid%n2ij( col,2)

          xl = grid%x( i) - grid%dx / 2._dp
          xu = grid%x( i) + grid%dx / 2._dp
          yl = grid%y( j) - grid%dx / 2._dp
          yu = grid%y( j) + grid%dx / 2._dp

          sw = [xl,yl]
          nw = [xl,yu]
          se = [xu,yl]
          ne = [xu,yu]

          ! Integrate around the grid cell
          ti_hint = ti
          count_coincidences = .false.
          call trace_line_tri( mesh, sw, se, single_row_grid, count_coincidences, ti_hint)
          call trace_line_tri( mesh, se, ne, single_row_grid, count_coincidences, ti_hint)
          call trace_line_tri( mesh, ne, nw, single_row_grid, count_coincidences, ti_hint)
          call trace_line_tri( mesh, nw, sw, single_row_grid, count_coincidences, ti_hint)

          ! Safety
          if (single_row_grid%n == 0) call crash('couldnt find any triangles overlapping with this grid cell!')

          ! Add contribution for this particular triangle
          do kk = 1, single_row_grid%n
            tj = single_row_grid%index_left( kk)
            if (tj == ti) then
              ! Add contribution to this triangle
              single_row_tri%LI_xdy(   k) = single_row_tri%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
              single_row_tri%LI_mxydx( k) = single_row_tri%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
              single_row_tri%LI_xydy(  k) = single_row_tri%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
              exit
            end if
          end do ! do kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          call add_entry_CSR_dist( A_xdy_b_g_CSR  , row, col, single_row_tri%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_b_g_CSR, row, col, single_row_tri%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_b_g_CSR , row, col, single_row_tri%LI_xydy(  k))

        end do ! do k = 1, single_row_tri%n

      else ! if (mesh%TriA( ti) < 10._dp * grid%dx**2) then
        ! This triangle is big enough that we can just average over the grid cells it contains

        mask_do_simple_average( ti) = 1

        ! Clean up single row results
        single_row_tri%n = 0

        ! Find the square of grid cells enveloping this triangle
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        pa = mesh%V( via,:)
        pb = mesh%V( vib,:)
        pc = mesh%V( vic,:)

        xmin = min( min( pa( 1), pb( 1)), pc( 1))
        xmax = max( max( pa( 1), pb( 1)), pc( 1))
        ymin = min( min( pa( 2), pb( 2)), pc( 2))
        ymax = max( max( pa( 2), pb( 2)), pc( 2))

        il = max( 1, min( grid%nx, 1 + floor( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        iu = max( 1, min( grid%nx, 1 + floor( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        jl = max( 1, min( grid%ny, 1 + floor( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx) ))
        ju = max( 1, min( grid%ny, 1 + floor( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

        ! Check which of the grid cells in this square lie inside the triangle
        do i = il, iu
        do j = jl, ju

          col = grid%ij2n( i,j)
          p   = [grid%x( i), grid%y( j)]

          if (is_in_triangle( pa, pb, pc, p)) then
            ! This grid cell lies inside the triangle; add it to the single row
            single_row_tri%n = single_row_tri%n + 1
            single_row_tri%index_left( single_row_tri%n) = col
            single_row_tri%LI_xdy(     single_row_tri%n) = grid%dx**2
            single_row_tri%LI_mxydx(   single_row_tri%n) = grid%x( i) * grid%dx**2
            single_row_tri%LI_xydy(    single_row_tri%n) = grid%y( j) * grid%dx**2
          end if

        end do
        end do

        ! Safety
        if (single_row_tri%n == 0) call crash('couldnt find any grid cells overlapping with this big triangle!')

        ! Add entries to the big matrices
        do k = 1, single_row_tri%n
          col = single_row_tri%index_left( k)
          call add_entry_CSR_dist( A_xdy_b_g_CSR  , ti, col, single_row_tri%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_b_g_CSR, ti, col, single_row_tri%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_b_g_CSR , ti, col, single_row_tri%LI_xydy(  k))
        end do

      end if ! if (mesh%TriA( ti) < 10._dp * grid%dx**2) then

    end do ! do row = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    deallocate( single_row_tri%index_left )
    deallocate( single_row_tri%LI_xdy     )
    deallocate( single_row_tri%LI_mxydx   )
    deallocate( single_row_tri%LI_xydy    )

    deallocate( single_row_grid%index_left )
    deallocate( single_row_grid%LI_xdy     )
    deallocate( single_row_grid%LI_mxydx   )
    deallocate( single_row_grid%LI_xydy    )

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    ! ==============================================================

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      k1 = A_xdy_b_g_CSR%ptr( row  )
      k2 = A_xdy_b_g_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_b_g_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_b_g_CSR%ind( k)
        call add_entry_CSR_dist( w0_CSR, row, col, A_xdy_b_g_CSR%val( k) / A_overlap_tot)
      end do

      if (mask_do_simple_average( ti) == 0) then
        ! For small triangles, include the gradient terms

        do k = k1, k2
          col = A_xdy_b_g_CSR%ind( k)
          ! Grid cell
          i = grid%n2ij( col,1)
          j = grid%n2ij( col,2)
          call add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_b_g_CSR%val( k) / A_overlap_tot) - (grid%x( i) * w0_CSR%val( k)))
          call add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_b_g_CSR%val(  k) / A_overlap_tot) - (grid%y( j) * w0_CSR%val( k)))
        end do

      else
        ! For large triangles, don't include the gradient terms

        call add_empty_row_CSR_dist( w1x_CSR, row)
        call add_empty_row_CSR_dist( w1y_CSR, row)

      end if ! if (mask_do_simple_average( vi) == 0) then

    end do ! do row = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( A_xdy_b_g_CSR  )
    call deallocate_matrix_CSR_dist( A_mxydx_b_g_CSR)
    call deallocate_matrix_CSR_dist( A_xydy_b_g_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix

    call calc_matrix_operators_grid( grid, grid_M_ddx, grid_M_ddy)

    call MatDuplicate( w0, MAT_COPY_VALUES, map%M, perr)
    call MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatDestroy( grid_M_ddx    , perr)
    call MatDestroy( grid_M_ddy    , perr)
    call MatDestroy( w0            , perr)
    call MatDestroy( w1x           , perr)
    call MatDestroy( w1y           , perr)

    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    call MatDestroy( M1, perr)
    call MatDestroy( M2, perr)
    deallocate( mask_do_simple_average)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_xy_grid_to_mesh_triangles

  subroutine create_map_from_mesh_to_xy_grid( mesh, grid, map)
    ! Create a new mapping object from a mesh to an x/y-grid.
    !
    ! By default uses 2nd-order conservative remapping.
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as having an area smaller
    !       than ten times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is generally very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_mesh_to_xy_grid'
    integer                                            :: ierr
    type(PetscErrorCode)                               :: perr
    logical                                            :: count_coincidences
    integer,  dimension(:,:  ), allocatable            :: overlaps_with_small_triangle, containing_triangle
    integer                                            :: row, ti
    integer                                            :: via, vib, vic
    real(dp), dimension(2)                             :: pa, pb, pc
    real(dp)                                           :: xmin, xmax, ymin, ymax
    integer                                            :: il, iu, jl, ju
    integer                                            :: i, j, n_ext, ii, jj
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR
    type(type_single_row_mapping_matrices)             :: single_row_grid, single_row_Tri
    integer                                            :: ti_hint
    real(dp), dimension(2)                             :: p
    real(dp)                                           :: xl, xu, yl, yu
    real(dp), dimension(2)                             :: sw, se, nw, ne
    integer                                            :: k, kk, nn
    real(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy
    type(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    type(tMat)                                         :: w0    , w1x    , w1y
    integer                                            :: k1, k2, col
    real(dp)                                           :: A_overlap_tot
    type(tMat)                                         :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    type(tMat)                                         :: M1, M2
    integer                                            :: ncols_row
    integer,  dimension(:    ), allocatable            :: cols_row
    real(dp), dimension(:    ), allocatable            :: vals_row
    logical                                            :: has_value

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh%name
    map%name_dst  = grid%name
    map%method    = '2nd_order_conservative'

    ! == Find all grid cells that overlap with small triangles
    ! ========================================================

    allocate( overlaps_with_small_triangle( grid%nx, grid%ny), source = 0)
    allocate( containing_triangle(          grid%nx, grid%ny), source = 0)

    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)

      ! The square enveloping this triangle
      xmin = min( min( pa(1), pb(1)), pc(1))
      xmax = max( max( pa(1), pb(1)), pc(1))
      ymin = min( min( pa(2), pb(2)), pc(2))
      ymax = max( max( pa(2), pb(2)), pc(2))

      ! The square of grid cells enveloping this triangle
      il = 1 + floor( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx)
      iu = 1 + floor( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx)
      jl = 1 + floor( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx)
      ju = 1 + floor( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx)

      il = max( 1      , il - 1)
      iu = min( grid%nx, iu + 1)
      jl = max( 1      , jl - 1)
      ju = min( grid%ny, ju + 1)

      if (mesh%TriA( ti) < 10._dp * grid%dx**2) then
        ! This triangle is small; mark all grid cells it overlaps with

        ! Mark all these grid cells
        do i = il, iu
        do j = jl, ju
          overlaps_with_small_triangle( i,j) = 1
        end do
        end do

      else
        ! This triangle is large; mark all grid cells it contains

        ! Mark all these grid cells
        do i = il, iu
        do j = jl, ju
          p = [grid%x( i), grid%y( j)]
          if (is_in_triangle( pa, pb, pc, p)) then
            containing_triangle( i,j) = ti
          end if
        end do
        end do

      end if ! if (mesh%TriA( ti) < 4._dp * grid%dx**2) then

    end do

     ! Reduce results across the processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, overlaps_with_small_triangle, grid%nx * grid%ny, MPI_integer, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, containing_triangle         , grid%nx * grid%ny, MPI_integer, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Treat grid cells that possibly were not yet marked before
    do row = grid%n1, grid%n2

      ! Grid cell indices
      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)

      if (containing_triangle( i,j) == 0 .and. overlaps_with_small_triangle( i,j) == 0) then
        ! This grid cell does not overlap with a small triangle, but was not yet marked
        ! as being contained inside a large one; find the large triangle containing it.

        ! For efficiency, find the nearest grid cell that does list which large
        ! triangle contains it; use that as a hint for the triangle search
        n_ext = 0
        ti_hint = 0
        do while (ti_hint == 0)
          n_ext = n_ext+1
          ! Safety
          if (n_ext > max( grid%nx, grid%ny)) exit
          il = max( 1      , i-n_ext)
          iu = min( grid%nx, i+n_ext)
          jl = max( 1      , j-n_ext)
          ju = min( grid%ny, j+n_ext)
          do ii = il, iu
          do jj = jl, ju
            if (containing_triangle( ii,jj) > 0) then
              ti_hint = containing_triangle( ii,jj)
              exit
            end if
          end do
          if (ti_hint > 0) exit
          end do
        end do
        if (ti_hint == 0) ti_hint = 1

        ! Find the triangle containing this grid cell
        p = [max( mesh%xmin, min( mesh%xmax, grid%x( i) )), max( mesh%ymin, min( mesh%ymax, grid%y( j) ))]
        call find_containing_triangle( mesh, p, ti_hint)
        containing_triangle( i,j) = ti_hint

      end if

    end do
    call sync

    ! == Integrate around all grid cells that overlap with small triangles
    ! ====================================================================

    ! Initialise the three matrices using the native UFEMISM CSR-matrix format

    ! Matrix size
    nrows           = grid%n     ! to
    nrows_loc       = grid%n_loc
    ncols           = mesh%nTri  ! from
    ncols_loc       = mesh%nTri_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%TriA) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%TriA))) )

    call allocate_matrix_CSR_dist( A_xdy_g_b_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_g_b_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_g_b_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for single row results
    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    allocate( single_row_grid%index_left( single_row_grid%n_max))
    allocate( single_row_grid%LI_xdy(     single_row_grid%n_max))
    allocate( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    allocate( single_row_grid%LI_xydy(    single_row_grid%n_max))

    single_row_Tri%n_max = nnz_per_row_max
    single_row_Tri%n     = 0
    allocate( single_row_Tri%index_left( single_row_Tri%n_max))
    allocate( single_row_Tri%LI_xdy(     single_row_Tri%n_max))
    allocate( single_row_Tri%LI_mxydx(   single_row_Tri%n_max))
    allocate( single_row_Tri%LI_xydy(    single_row_Tri%n_max))

    ti_hint = 1

    do row = grid%n1, grid%n2

      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)
      p = [grid%x( i), grid%y( j)]

      ! The four sides of the grid cell
      xl = grid%x( i) - grid%dx / 2._dp
      xu = grid%x( i) + grid%dx / 2._dp
      yl = grid%y( j) - grid%dx / 2._dp
      yu = grid%y( j) + grid%dx / 2._dp

      ! if this grid cell lies entirely outside of the mesh domain, use
      ! nearest-neighbour extrapolation

      do while (xl <= mesh%xmin)
        i = i+1
        p( 1) = grid%x( i)
        xl = grid%x( i) - grid%dx / 2._dp
        xu = grid%x( i) + grid%dx / 2._dp
        if (i > grid%nx) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do
      do while (xu >= mesh%xmax)
        i = i-1
        p( 1) = grid%x( i)
        xl = grid%x( i) - grid%dx / 2._dp
        xu = grid%x( i) + grid%dx / 2._dp
        if (i < 1) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do
      do while (yl <= mesh%ymin)
        j = j+1
        p( 2) = grid%y( j)
        yl = grid%y( j) - grid%dx / 2._dp
        yu = grid%y( j) + grid%dx / 2._dp
        if (j > grid%ny) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do
      do while (yu >= mesh%ymax)
        j = j-1
        p( 2) = grid%y( j)
        yl = grid%y( j) - grid%dx / 2._dp
        yu = grid%y( j) + grid%dx / 2._dp
        if (j < 1) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do

      if (overlaps_with_small_triangle( i,j) == 1) then
        ! This grid cell overlaps with a small triangle; integrate around it, and around
        ! all triangles overlapping with it

        sw = [xl, yl]
        nw = [xl, yu]
        se = [xu, yl]
        ne = [xu, yu]

        ! Clear the single row results
        single_row_grid%n          = 0
        single_row_grid%index_left = 0
        single_row_grid%LI_xdy     = 0._dp
        single_row_grid%LI_mxydx   = 0._dp
        single_row_grid%LI_xydy    = 0._dp

        ! Integrate over all four sides
        count_coincidences = .true.
        call trace_line_tri( mesh, sw, se, single_row_grid, count_coincidences, ti_hint)
        call trace_line_tri( mesh, se, ne, single_row_grid, count_coincidences, ti_hint)
        call trace_line_tri( mesh, ne, nw, single_row_grid, count_coincidences, ti_hint)
        call trace_line_tri( mesh, nw, sw, single_row_grid, count_coincidences, ti_hint)

        ! Next, integrate around all the triangles overlapping with this grid cell
        do k = 1, single_row_grid%n

          ti = single_row_grid%index_left( k)

          col = mesh%ti2n( ti)

          ! The three vertices spanning this triangle
          via = mesh%Tri( ti,1)
          vib = mesh%Tri( ti,2)
          vic = mesh%Tri( ti,3)

          pa  = mesh%V( via,:)
          pb  = mesh%V( vib,:)
          pc  = mesh%V( vic,:)

          ! Clear the single row results
          single_row_Tri%n = 0
          single_row_Tri%index_left = 0
          single_row_Tri%LI_xdy     = 0._dp
          single_row_Tri%LI_mxydx   = 0._dp
          single_row_Tri%LI_xydy    = 0._dp

          ! Integrate over all three triangle sides
          count_coincidences = .false.
          call trace_line_grid( grid, pa, pb, single_row_Tri, count_coincidences)
          call trace_line_grid( grid, pb, pc, single_row_Tri, count_coincidences)
          call trace_line_grid( grid, pc, pa, single_row_Tri, count_coincidences)

          ! Add contribution for this particular grid cell
          do kk = 1, single_row_Tri%n
            nn = single_row_Tri%index_left( kk)
            if (nn == row) then
              ! Add contribution to this triangle
              single_row_grid%LI_xdy(   k) = single_row_grid%LI_xdy(   k) + single_row_Tri%LI_xdy(   kk)
              single_row_grid%LI_mxydx( k) = single_row_grid%LI_mxydx( k) + single_row_Tri%LI_mxydx( kk)
              single_row_grid%LI_xydy(  k) = single_row_grid%LI_xydy(  k) + single_row_Tri%LI_xydy(  kk)
              exit
            end if
          end do ! do kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          call add_entry_CSR_dist( A_xdy_g_b_CSR  , row, col, single_row_grid%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_g_b_CSR, row, col, single_row_grid%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_g_b_CSR , row, col, single_row_grid%LI_xydy(  k))

        end do ! do k = 1, single_row_grid%n

      else ! if (overlaps_with_small_triangle( i,j) == 1) then
        ! This grid cell does not overlap with a small triangle; use only the
        ! contribution from the nearest triangle

        ti_hint = containing_triangle( i,j)

        col = mesh%ti2n( ti_hint)

        LI_xdy   = grid%dx**2
        LI_mxydx = grid%dx**2 * grid%x( i)
        LI_xydy  = grid%dx**2 * grid%y( j)

        call add_entry_CSR_dist( A_xdy_g_b_CSR  , row, col, LI_xdy  )
        call add_entry_CSR_dist( A_mxydx_g_b_CSR, row, col, LI_mxydx)
        call add_entry_CSR_dist( A_xydy_g_b_CSR , row, col, LI_xydy )

      end if ! if (overlaps_with_small_triangle( i,j) == 1) then

    end do ! do n = grid%n1, grid%n2

    ! Clean up after yourself
    deallocate( overlaps_with_small_triangle)
    deallocate( containing_triangle         )

    deallocate( single_row_grid%index_left )
    deallocate( single_row_grid%LI_xdy     )
    deallocate( single_row_grid%LI_mxydx   )
    deallocate( single_row_grid%LI_xydy    )

    deallocate( single_row_Tri%index_left )
    deallocate( single_row_Tri%LI_xdy     )
    deallocate( single_row_Tri%LI_mxydx   )
    deallocate( single_row_Tri%LI_xydy    )

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    ! ==============================================================

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = grid%n1, grid%n2

      k1 = A_xdy_g_b_CSR%ptr( row  )
      k2 = A_xdy_g_b_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_g_b_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_g_b_CSR%ind( k)
        ti = mesh%n2ti( col)
        call add_entry_CSR_dist( w0_CSR , row, col,  A_xdy_g_b_CSR%val(   k) / A_overlap_tot)
        call add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_g_b_CSR%val( k) / A_overlap_tot) - (mesh%TriGC( ti,1) * w0_CSR%val( k)))
        call add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_g_b_CSR%val(  k) / A_overlap_tot) - (mesh%TriGC( ti,2) * w0_CSR%val( k)))
      end do

    end do ! do row = mesh%vi1, mesh%vi2
    call sync

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( A_xdy_g_b_CSR  )
    call deallocate_matrix_CSR_dist( A_mxydx_g_b_CSR)
    call deallocate_matrix_CSR_dist( A_xydy_g_b_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix
    ! ==============================

    ! Convert matrices to PETSc format
    call mat_CSR2petsc( mesh%M_map_a_b, M_map_a_b)
    call mat_CSR2petsc( mesh%M_ddx_a_b, M_ddx_a_b)
    call mat_CSR2petsc( mesh%M_ddy_a_b, M_ddy_a_b)

    ! M = (w0 * M_map_a_b) + (w1x * M_ddx_a_b) + (w1y * M_ddy_a_b)
    call MatMatMult( w0,  M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, map%M, perr)
    call MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1   , perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2   , perr)
    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    call MatDestroy( w0            , perr)
    call MatDestroy( w1x           , perr)
    call MatDestroy( w1y           , perr)
    call MatDestroy( M_map_a_b     , perr)
    call MatDestroy( M_ddx_a_b     , perr)
    call MatDestroy( M_ddy_a_b     , perr)
    call MatDestroy( M1            , perr)
    call MatDestroy( M2            , perr)

    ! Safety: check if all grid cells get values
    ! ==========================================

    allocate( cols_row( nnz_per_row_max))
    allocate( vals_row( nnz_per_row_max))

    do row = grid%n1, grid%n2

      ! w0
      call MatGetRow( map%M, row-1, ncols_row, cols_row, vals_row, perr)
      if (ncols_row == 0) call crash('ncols == 0!')
      has_value = .false.
      do k = 1, ncols_row
        if (vals_row( k) /= 0._dp) has_value = .true.
      end do
      if (.not. has_value) call crash('only zeroes!')
      call MatRestoreRow( map%M, row-1, ncols_row, cols_row, vals_row, perr)

    end do
    call sync

    ! Clean up after yourself
    deallocate( cols_row)
    deallocate( vals_row)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_xy_grid

end module create_maps_grid_mesh
