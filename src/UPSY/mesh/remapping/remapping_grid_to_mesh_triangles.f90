module remapping_grid_to_mesh_triangles

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
  use grid_basic, only: calc_matrix_operators_grid
  use netcdf_output

  implicit none

  private

  public :: create_map_from_xy_grid_to_mesh_triangles

contains

  subroutine create_map_from_xy_grid_to_mesh_triangles( grid, mesh, output_dir, map)
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
    type(type_grid),  intent(in   ) :: grid
    type(type_mesh),  intent(in   ) :: mesh
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'create_map_from_xy_grid_to_mesh_triangles'
    logical, dimension(mesh%ti1:mesh%ti2) :: lies_outside_grid_domain
    logical, dimension(mesh%ti1:mesh%ti2) :: is_large_triangle
    type(type_sparse_matrix_CSR_dp)       :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR
    type(tMat)                            :: w0, w1x, w1y
    type(type_sparse_matrix_CSR_dp)       :: grid_M_ddx_CSR, grid_M_ddy_CSR
    character(len=1024)                   :: filename_grid, filename_mesh

    ! Add routine to path
    call init_routine( routine_name)

    call dump_grid_and_mesh_to_netcdf( grid, mesh, output_dir, filename_grid, filename_mesh)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = grid%name
    map%name_dst  = trim( mesh%name) // '_triangles'
    map%method    = '2nd_order_conservative'

    call find_triangles_outside_grid_domain( grid, mesh, lies_outside_grid_domain)
    call find_large_triangles( grid, mesh, is_large_triangle)

    call calc_A_matrices( grid, mesh, &
      lies_outside_grid_domain, is_large_triangle, &
      A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR)

    call calc_w_matrices( grid, mesh, &
      lies_outside_grid_domain, is_large_triangle, &
      A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR, w0, w1x, w1y)

    call calc_matrix_operators_grid( grid, grid_M_ddx_CSR, grid_M_ddy_CSR)

    call calc_remapping_matrix( w0, w1x, w1y, grid_M_ddx_CSR, grid_M_ddy_CSR, map%M)

    call delete_grid_and_mesh_netcdf_dump_files( filename_grid, filename_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_xy_grid_to_mesh_triangles

  subroutine find_triangles_outside_grid_domain( grid, mesh, lies_outside_grid_domain)
    !< Identify mesh vertices whose Voronoi cell lies partially outside the grid domain
    !< (so they can be skipped, as remapping is untrustworthy there)

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%ti1:mesh%ti2), intent(  out) :: lies_outside_grid_domain

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_vertices_outside_grid_domain'
    integer                        :: ti, n, vi

    ! Add routine to path
    call init_routine( routine_name)

    lies_outside_grid_domain = .false.
    do ti = mesh%ti1, mesh%ti2
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        if (mesh%V( vi,1) <= grid%xmin - grid%dx/2._dp .or. &
            mesh%V( vi,1) >= grid%xmax + grid%dx/2._dp .or. &
            mesh%V( vi,2) <= grid%ymin - grid%dx/2._dp .or. &
            mesh%V( vi,2) >= grid%ymax + grid%dx/2._dp) then
          lies_outside_grid_domain( ti) = .true.
        end if
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_triangles_outside_grid_domain

  subroutine find_large_triangles( grid, mesh, is_large_triangle)
    !< Identify triangles that are large enough that we can simply average
    !< over the overlapping grid cells, without needing to calculate all the line integrals

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%ti1:mesh%ti2), intent(  out) :: is_large_triangle

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_large_triangles'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    is_large_triangle = .false.
    do ti = mesh%ti1, mesh%ti2
      if (mesh%TriA( ti) >= 10._dp * grid%dx**2) then
        is_large_triangle( ti) = .true.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_large_triangles

  subroutine calc_A_matrices( grid, mesh, &
    lies_outside_grid_domain, is_large_triangle, &
    A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR)
    !< Calculate the A-matrices for the grid-to-mesh remapping operator

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%ti1:mesh%ti2), intent(in   ) :: lies_outside_grid_domain
    logical, dimension(mesh%ti1:mesh%ti2), intent(in   ) :: is_large_triangle
    type(type_sparse_matrix_CSR_dp),       intent(  out) :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_A_matrices'
    integer                                :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    integer                                :: row, ti
    type(type_single_row_mapping_matrices) :: single_row_tri, single_row_grid

    ! Add routine to path
    call init_routine( routine_name)

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

    ! Calculate line integrals around all triangles
    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      ! If this triangle lies (partially) outside the domain of the grid, skip it.
      if (lies_outside_grid_domain( ti)) then
        call add_empty_row_CSR_dist( A_xdy_b_g_CSR,   row)
        call add_empty_row_CSR_dist( A_mxydx_b_g_CSR, row)
        call add_empty_row_CSR_dist( A_xydy_b_g_CSR,  row)
        cycle
      end if

      if (is_large_triangle( ti)) then
        call calc_A_matrices_large_triangle( grid, mesh, ti, &
          single_row_tri, A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR)
      else
        call calc_A_matrices_small_triangle( grid, mesh, row, ti, &
          single_row_tri, single_row_grid, A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR)
      end if

    end do

    call finalise_matrix_CSR_dist( A_xdy_b_g_CSR)
    call finalise_matrix_CSR_dist( A_mxydx_b_g_CSR)
    call finalise_matrix_CSR_dist( A_xydy_b_g_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_A_matrices

  subroutine calc_A_matrices_small_triangle( grid, mesh, row, ti, &
    single_row_tri, single_row_grid, A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR)
    !< Calculate a single row of the A-matrices, for a small triangle
    !< (i.e., calculating the line integrals rather than simply averaging over the grid cells)

    ! In/output variables
    type(type_grid),                        intent(in   ) :: grid
    type(type_mesh),                        intent(in   ) :: mesh
    integer,                                intent(in   ) :: row, ti
    type(type_single_row_mapping_matrices), intent(inout) :: single_row_tri, single_row_grid
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR

    ! Local variables:
    logical                        :: count_coincidences
    integer                        :: n1, n2, via, vib
    real(dp), dimension(2)         :: p, q
    integer                        :: k, i, j, kk, tj, col
    real(dp)                       :: xl, xu, yl, yu
    real(dp), dimension(2)         :: sw, se, nw, ne
    integer                        :: ti_hint

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
    if (single_row_tri%n == 0) then
      call crash('couldnt find any grid cells overlapping with small triangle ti = {int_01}', int_01 = ti)
    end if

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
      end do

      ! Add entries to the big matrices
      call add_entry_CSR_dist( A_xdy_b_g_CSR  , row, col, single_row_tri%LI_xdy(   k))
      call add_entry_CSR_dist( A_mxydx_b_g_CSR, row, col, single_row_tri%LI_mxydx( k))
      call add_entry_CSR_dist( A_xydy_b_g_CSR , row, col, single_row_tri%LI_xydy(  k))

    end do

  end subroutine calc_A_matrices_small_triangle

  subroutine calc_A_matrices_large_triangle( grid, mesh, ti, &
    single_row_tri, A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR)
    !< Calculate a single row of the A-matrices, for a large triangle
    !< (i.e., simply averaging over the grid cells rather than calculating the line integrals)

    ! In/output variables
    type(type_grid),                        intent(in   ) :: grid
    type(type_mesh),                        intent(in   ) :: mesh
    integer,                                intent(in   ) :: ti
    type(type_single_row_mapping_matrices), intent(inout) :: single_row_tri
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR

    ! Local variables:
    integer                :: via, vib, vic
    real(dp), dimension(2) :: p
    integer                :: k, i, j, col
    real(dp)               :: xmin, xmax, ymin, ymax
    integer                :: il, iu, jl, ju
    real(dp), dimension(2) :: pa, pb, pc

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

  end subroutine calc_A_matrices_large_triangle

  subroutine calc_w_matrices( grid, mesh, &
    lies_outside_grid_domain, is_large_triangle, &
    A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR, w0, w1x, w1y)
    !< Calculate the w-matrices for the grid-to-mesh remapping operator

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%ti1:mesh%ti2), intent(in   ) :: lies_outside_grid_domain
    logical, dimension(mesh%ti1:mesh%ti2), intent(in   ) :: is_large_triangle
    type(type_sparse_matrix_CSR_dp),       intent(in   ) :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR
    type(tMat),                            intent(  out) :: w0, w1x, w1y

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_w_matrices'
    integer                         :: nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc
    integer                         :: ti
    integer                         :: k, i, j
    type(type_sparse_matrix_CSR_dp) :: w0_CSR, w1x_CSR, w1y_CSR
    integer                         :: row, k1, k2, col
    real(dp)                        :: A_overlap_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    ! ==============================================================

    nrows        = A_xdy_b_g_CSR%m
    ncols        = A_xdy_b_g_CSR%n
    nrows_loc    = A_xdy_b_g_CSR%m_loc
    ncols_loc    = A_xdy_b_g_CSR%n_loc
    nnz_est_proc = A_xdy_b_g_CSR%nnz

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      if (lies_outside_grid_domain( ti)) then
        ! Skip these
        call add_empty_row_CSR_dist( w0_CSR,  row)
        call add_empty_row_CSR_dist( w1x_CSR, row)
        call add_empty_row_CSR_dist( w1y_CSR, row)
        cycle
      end if

      k1 = A_xdy_b_g_CSR%ptr( row  )
      k2 = A_xdy_b_g_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_b_g_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_b_g_CSR%ind( k)
        call add_entry_CSR_dist( w0_CSR, row, col, A_xdy_b_g_CSR%val( k) / A_overlap_tot)
      end do

      if (.not. is_large_triangle( ti)) then
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

      end if

    end do

    call finalise_matrix_CSR_dist( w0_CSR )
    call finalise_matrix_CSR_dist( w1x_CSR)
    call finalise_matrix_CSR_dist( w1y_CSR)

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_w_matrices

  subroutine calc_remapping_matrix( w0, w1x, w1y, grid_M_ddx_CSR, grid_M_ddy_CSR, M)
    !< Calculate the grid-to-mesh-vertices remapping matrix M

    ! In/output variables
    type(tMat),                      intent(in   ) :: w0, w1x, w1y
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: grid_M_ddx_CSR, grid_M_ddy_CSR
    type(tMat),                      intent(  out) :: M

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_remapping_matrix'
    type(PetscErrorCode)           :: perr
    type(tMat)                      :: grid_M_ddx, grid_M_ddy, M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    call mat_CSR2petsc( grid_M_ddx_CSR, grid_M_ddx)
    call mat_CSR2petsc( grid_M_ddy_CSR, grid_M_ddy)

    ! M = w0 + w1x * M_ddx + w1y * M_ddy

    ! M = w0
    call MatDuplicate( w0, MAT_COPY_VALUES, M, perr)

    ! M = w0 + w1x * M_ddx
    call MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)
    call MatAXPY( M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)

    ! M = w0 + w1x * M_ddx + w1y * M_ddy
    call MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)
    call MatAXPY( M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_remapping_matrix

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

end module remapping_grid_to_mesh_triangles
