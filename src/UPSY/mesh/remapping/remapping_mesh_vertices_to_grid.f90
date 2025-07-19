module remapping_mesh_vertices_to_grid

  ! Create remapping objects between a square Cartesian grid and a mesh.

#include <petsc/finclude/petscksp.h>
  use petscksp
  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, finalise_matrix_CSR_dist, &
    deallocate_matrix_CSR_dist, add_entry_CSR_dist, add_empty_row_CSR_dist
  use remapping_types, only: type_single_row_mapping_matrices, type_map
  use plane_geometry, only: is_in_triangle
  use mesh_utilities, only: find_containing_triangle, calc_Voronoi_cell, is_in_Voronoi_cell
  use petsc_basic, only: mat_CSR2petsc
  use line_tracing_grid, only: trace_line_grid
  use line_tracing_triangles, only: trace_line_tri
  use line_tracing_Voronoi, only: trace_line_Vor
  use netcdf_output

  implicit none

  private

  public :: create_map_from_mesh_vertices_to_xy_grid, calc_approximate_overlaps, calc_A_matrices, &
    calc_w_matrices, dump_grid_and_mesh_to_netcdf, delete_grid_and_mesh_netcdf_dump_files

contains

  subroutine create_map_from_mesh_vertices_to_xy_grid( mesh, grid, output_dir, map)
    !< Create a new mapping object from a mesh to an x/y-grid.

    ! By default uses 2nd-order conservative remapping.

    ! In/output variables
    type(type_mesh),  intent(in   ) :: mesh
    type(type_grid),  intent(in   ) :: grid
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'create_map_from_mesh_vertices_to_xy_grid'
    integer, dimension(grid%nx, grid%ny)   :: overlaps_with_small_triangle, containing_triangle
    type(type_sparse_matrix_CSR_dp)        :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR
    type(tMat)                             :: w0, w1x, w1y
    character(len=1024)                    :: filename_grid, filename_mesh

    ! Add routine to path
    call init_routine( routine_name)

    call dump_grid_and_mesh_to_netcdf( grid, mesh, output_dir, filename_grid, filename_mesh)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = mesh%name
    map%name_dst  = grid%name
    map%method    = '2nd_order_conservative'

    call calc_approximate_overlaps( mesh, grid, &
      overlaps_with_small_triangle, containing_triangle)

    call calc_A_matrices( mesh, grid, &
      overlaps_with_small_triangle, containing_triangle, &
      A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR)

    call calc_w_matrices( mesh, grid, &
      A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR, w0, w1x, w1y)

    call calc_remapping_matrix( mesh, w0, w1x, w1y, map%M)

    call check_remapping_matrix_validity( mesh, grid, map%M)

    call delete_grid_and_mesh_netcdf_dump_files( filename_grid, filename_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_vertices_to_xy_grid

  subroutine calc_approximate_overlaps( mesh, grid, &
    overlaps_with_small_triangle, containing_triangle)
    !< Determine which grid cells overlap with which small, and with which large triangles

    ! In/output variables
    type(type_mesh),                      intent(in   ) :: mesh
    type(type_grid),                      intent(in   ) :: grid
    integer, dimension(grid%nx, grid%ny), intent(  out) :: overlaps_with_small_triangle
    integer, dimension(grid%nx, grid%ny), intent(  out) :: containing_triangle

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_approximate_overlaps'
    integer                                :: ierr
    integer                                :: row, ti
    integer                                :: via, vib, vic
    real(dp), dimension(2)                 :: pa, pb, pc, p
    real(dp)                               :: xmin, xmax, ymin, ymax
    integer                                :: il, iu, jl, ju
    integer                                :: i, j, n_ext, ii, jj, ti_hint

    ! Add routine to path
    call init_routine( routine_name)

    overlaps_with_small_triangle = 0
    containing_triangle          = 0

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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_approximate_overlaps

  subroutine calc_A_matrices( mesh, grid, &
    overlaps_with_small_triangle, containing_triangle, &
    A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR)
    !< Calculate the A-matrices for the mesh-to-grid remapping operator

    ! In/output variables
    type(type_mesh),                      intent(in   ) :: mesh
    type(type_grid),                      intent(in   ) :: grid
    integer, dimension(grid%nx, grid%ny), intent(in   ) :: overlaps_with_small_triangle
    integer, dimension(grid%nx, grid%ny), intent(in   ) :: containing_triangle
    type(type_sparse_matrix_CSR_dp),      intent(  out) :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_A_matrices'
    logical                                :: count_coincidences
    integer                                :: row, ti
    integer                                :: via, vib, vic
    real(dp), dimension(2)                 :: pa, pb, pc
    integer                                :: i, j
    integer                                :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_single_row_mapping_matrices) :: single_row_grid, single_row_Tri
    integer                                :: ti_hint
    real(dp), dimension(2)                 :: p
    real(dp)                               :: xl, xu, yl, yu
    real(dp), dimension(2)                 :: sw, se, nw, ne
    integer                                :: k, kk, nn, col
    real(dp)                               :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    call init_routine( routine_name)

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

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_A_matrices

  subroutine calc_w_matrices( mesh, grid, &
    A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR, w0, w1x, w1y)
    !< Calculate the w-matrices for the mesh-to-grid remapping operator

    ! In/output variables
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_grid),                 intent(in   ) :: grid
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR
    type(tMat)                     , intent(  out) :: w0, w1x, w1y

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_w_matrices'
    integer                         :: row, ti
    integer                         :: nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc
    integer                         :: k
    type(type_sparse_matrix_CSR_dp) :: w0_CSR, w1x_CSR, w1y_CSR
    integer                         :: k1, k2, col
    real(dp)                        :: A_overlap_tot

    ! Add routine to path
    call init_routine( routine_name)

    nrows        = A_xdy_g_b_CSR%m
    ncols        = A_xdy_g_b_CSR%n
    nrows_loc    = A_xdy_g_b_CSR%m_loc
    ncols_loc    = A_xdy_g_b_CSR%n_loc
    nnz_est_proc = A_xdy_g_b_CSR%nnz

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

    end do

    call finalise_matrix_CSR_dist( w0_CSR)
    call finalise_matrix_CSR_dist( w1x_CSR)
    call finalise_matrix_CSR_dist( w1y_CSR)

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_w_matrices

  subroutine calc_remapping_matrix( mesh, w0, w1x, w1y, M)
    !< Calculate the mesh-to-grid-vertices remapping matrix M

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh
    type(tMat),      intent(inout) :: w0, w1x, w1y
    type(tMat),      intent(  out) :: M

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_remapping_matrix'
    type(PetscErrorCode)           :: perr
    type(tMat)                     :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    type(tMat)                     :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert matrices to PETSc format
    call mat_CSR2petsc( mesh%M_map_a_b, M_map_a_b)
    call mat_CSR2petsc( mesh%M_ddx_a_b, M_ddx_a_b)
    call mat_CSR2petsc( mesh%M_ddy_a_b, M_ddy_a_b)

    ! M = (w0 * M_map_a_b) + (w1x * M_ddx_a_b) + (w1y * M_ddy_a_b)

    ! M = (w0 * M_map_a_b)
    call MatMatMult( w0,  M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M, perr)

    ! M = (w0 * M_map_a_b) + (w1x * M_ddx_a_b)
    call MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)
    call MatAXPY( M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)

    ! M = (w0 * M_map_a_b) + (w1x * M_ddx_a_b) + (w1y * M_ddy_a_b)
    call MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)
    call MatAXPY( M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_remapping_matrix

  subroutine check_remapping_matrix_validity( mesh, grid, M)
    !< Safety: check if all grid cells get values

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh
    type(type_grid), intent(in   ) :: grid
    type(tMat),      intent(in   ) :: M

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'check_remapping_matrix_validity'
    type(PetscErrorCode)                        :: perr
    integer                                     :: k, row
    integer                                     :: nnz_per_row_max
    integer                                     :: ncols_row
    integer,  dimension(:), allocatable, target :: cols_row
    real(dp), dimension(:), allocatable, target :: vals_row
    integer,  dimension(:), pointer             :: cols_row_
    real(dp), dimension(:), pointer             :: vals_row_
    logical                                     :: has_value

    ! Add routine to path
    call init_routine( routine_name)

    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%TriA) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%TriA))) )

    allocate( cols_row( nnz_per_row_max))
    allocate( vals_row( nnz_per_row_max))

    cols_row_ => cols_row
    vals_row_ => vals_row

    do row = grid%n1, grid%n2

      call MatGetRow( M, row-1, ncols_row, cols_row_, vals_row_, perr)

      if (ncols_row == 0) call crash('ncols == 0!')

      has_value = .false.
      do k = 1, ncols_row
        if (vals_row_( k) /= 0._dp) has_value = .true.
      end do
      if (.not. has_value) call crash('only zeroes!')

      call MatRestoreRow( M, row-1, ncols_row, cols_row_, vals_row_, perr)

    end do

    deallocate( cols_row)
    deallocate( vals_row)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_remapping_matrix_validity

  subroutine dump_grid_and_mesh_to_netcdf( grid, mesh, output_dir, filename_grid, filename_mesh)
    !< Dump grid and mesh to NetCDF files that will be deleted after the remapping
    !< operator has been successfully calculated. If the remapping crashes, the NetCDF
    !< files will still be there, so that Tijn can use them to figure out the bug.

    ! In/output variables
    type(type_grid),  intent(in   ) :: grid
    type(type_mesh),  intent(in   ) :: mesh
    character(len=*), intent(in   ) :: output_dir
    character(len=*), intent(  out) :: filename_grid, filename_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'dump_grid_and_mesh_to_netcdf'

    ! Add routine to path
    call init_routine( routine_name)

    filename_grid = trim(output_dir) // '/grid2mesh_grid_dump.nc'
    filename_mesh = trim(output_dir) // '/grid2mesh_mesh_dump.nc'

    call save_xy_grid_as_netcdf( filename_grid, grid)
    call save_mesh_as_netcdf(    filename_mesh, mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine dump_grid_and_mesh_to_netcdf

  subroutine delete_grid_and_mesh_netcdf_dump_files( filename_grid, filename_mesh)

    ! In/output variables
    character(len=*), intent(in   ) :: filename_grid, filename_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_grid_and_mesh_netcdf_dump_files'
    integer                        :: stat

    ! Add routine to path
    call init_routine( routine_name)

    ! Delete grid & mesh netcdf dumps
    if (par%primary) then
      open(unit = 1234, iostat = stat, file = filename_grid, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
      open(unit = 1234, iostat = stat, file = filename_mesh, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_grid_and_mesh_netcdf_dump_files

end module remapping_mesh_vertices_to_grid
