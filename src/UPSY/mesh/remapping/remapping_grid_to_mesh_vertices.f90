module remapping_grid_to_mesh_vertices

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

  public :: create_map_from_xy_grid_to_mesh_vertices

contains

  subroutine create_map_from_xy_grid_to_mesh_vertices( grid, mesh, output_dir, map)
    !< Create a new mapping object from an x/y-grid to the vertices of a mesh.

    ! By default uses 2nd-order conservative remapping.
    !
    ! NOTE: the current implementation is a compromise. For "small" vertices (defined as having a Voronoi cell smaller
    !       than ten times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" vertices (defined as
    !       all the rest), the result is very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    ! In/output variables
    type(type_grid),  intent(in   ) :: grid
    type(type_mesh),  intent(in   ) :: mesh
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'create_map_from_xy_grid_to_mesh_vertices'
    logical, dimension(mesh%vi1:mesh%vi2) :: lies_outside_grid_domain
    logical, dimension(mesh%vi1:mesh%vi2) :: is_large_vertex
    type(type_sparse_matrix_CSR_dp)       :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR
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
    map%name_dst  = mesh%name
    map%method    = '2nd_order_conservative'

    call find_vertices_outside_grid_domain( grid, mesh, lies_outside_grid_domain)
    call find_large_vertices( grid, mesh, is_large_vertex)

    call calc_A_matrices( grid, mesh, &
      lies_outside_grid_domain, is_large_vertex, &
      A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR)

    call calc_w_matrices( grid, mesh, &
      lies_outside_grid_domain, is_large_vertex, &
      A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR, w0, w1x, w1y)

    call calc_matrix_operators_grid( grid, grid_M_ddx_CSR, grid_M_ddy_CSR)

    call calc_remapping_matrix( w0, w1x, w1y, grid_M_ddx_CSR, grid_M_ddy_CSR, map%M)

    call delete_grid_and_mesh_netcdf_dump_files( filename_grid, filename_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_xy_grid_to_mesh_vertices

  subroutine find_vertices_outside_grid_domain( grid, mesh, lies_outside_grid_domain)
    !< Identify mesh vertices whose Voronoi cell lies partially outside the grid domain
    !< (so they can be skipped, as remapping is untrustworthy there)

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(  out) :: lies_outside_grid_domain

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_vertices_outside_grid_domain'
    integer                        :: vi, vvi, vori

    ! Add routine to path
    call init_routine( routine_name)

    lies_outside_grid_domain = .false.
    do vi = mesh%vi1, mesh%vi2
      do vvi = 1, mesh%nVVor( vi)
        vori = mesh%VVor( vi,vvi)
        if (mesh%Vor( vori,1) <= grid%xmin - grid%dx/2._dp .or. &
            mesh%Vor( vori,1) >= grid%xmax + grid%dx/2._dp .or. &
            mesh%Vor( vori,2) <= grid%ymin - grid%dx/2._dp .or. &
            mesh%Vor( vori,2) >= grid%ymax + grid%dx/2._dp) then
          lies_outside_grid_domain( vi) = .true.
        end if
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_vertices_outside_grid_domain

  subroutine find_large_vertices( grid, mesh, is_large_vertex)
    !< Identify vertices whose Voronoi cell is large enough that we can simply average
    !< over the overlapping grid cells, without needing to calculate all the line integrals

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(  out) :: is_large_vertex

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_large_vertices'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    is_large_vertex = .false.
    do vi = mesh%vi1, mesh%vi2
      if (mesh%A( vi) >= 10._dp * grid%dx**2) then
        is_large_vertex( vi) = .true.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_large_vertices

  subroutine calc_A_matrices( grid, mesh, &
    lies_outside_grid_domain, is_large_vertex, &
    A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR)
    !< Calculate the A-matrices for the grid-to-mesh remapping operator

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: lies_outside_grid_domain
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: is_large_vertex
    type(type_sparse_matrix_CSR_dp),       intent(  out) :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_A_matrices'
    integer                                :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    integer                                :: vi, row
    type(type_single_row_mapping_matrices) :: single_row_Vor, single_row_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Matrix sizes
    ncols           = grid%n       ! from
    ncols_loc       = grid%n_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%A) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%A))) )

    call allocate_matrix_CSR_dist( A_xdy_a_g_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_a_g_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_a_g_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for single row results
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

    ! Calculate line integrals around all Voronoi cells
    do row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      ! If the Voronoi cell of this vertex lies (partially) outside the domain of the grid, skip it.
      if (lies_outside_grid_domain( vi)) then
        call add_empty_row_CSR_dist( A_xdy_a_g_CSR,   row)
        call add_empty_row_CSR_dist( A_mxydx_a_g_CSR, row)
        call add_empty_row_CSR_dist( A_xydy_a_g_CSR,  row)
        cycle
      end if

      ! Clean up single row results
      single_row_Vor%n          = 0
      single_row_Vor%index_left = 0
      single_row_Vor%LI_xdy     = 0._dp
      single_row_Vor%LI_mxydx   = 0._dp
      single_row_Vor%LI_xydy    = 0._dp

      if (is_large_vertex( vi)) then
        call calc_A_matrices_large_vertex( grid, mesh, row, vi, &
         single_row_Vor, A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR)
      else
        call calc_A_matrices_small_vertex( grid, mesh, row, vi, &
          single_row_Vor, single_row_grid, A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR)
      end if

    end do

    call finalise_matrix_CSR_dist( A_xdy_a_g_CSR)
    call finalise_matrix_CSR_dist( A_mxydx_a_g_CSR)
    call finalise_matrix_CSR_dist( A_xydy_a_g_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_A_matrices

  subroutine calc_A_matrices_small_vertex( grid, mesh, row, vi, &
    single_row_Vor, single_row_grid, A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR)
    !< Calculate a single row of the A-matrices, for a small Voronoi cell
    !< (i.e., calculating the line integrals rather than simply averaging over the grid cells)

    ! In/output variables
    type(type_grid),                        intent(in   ) :: grid
    type(type_mesh),                        intent(in   ) :: mesh
    integer,                                intent(in   ) :: row, vi
    type(type_single_row_mapping_matrices), intent(inout) :: single_row_Vor, single_row_grid
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR

    ! Local variables:
    logical                                :: count_coincidences
    integer                                :: col
    real(dp), dimension( mesh%nC_mem,2)    :: Vor
    integer,  dimension( mesh%nC_mem  )    :: Vor_vi
    integer,  dimension( mesh%nC_mem  )    :: Vor_ti
    integer                                :: nVor
    integer                                :: vori1, vori2
    real(dp), dimension(2)                 :: p, q
    integer                                :: k, i, j, kk, vj
    real(dp)                               :: xl, xu, yl, yu
    real(dp), dimension(2)                 :: sw, se, nw, ne
    integer                                :: vi_hint

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
    if (single_row_Vor%n == 0) then
      call crash('couldnt find any grid cells overlapping with the small ' // &
        'Voronoi cell of vertex {int_01}', int_01 = vi)
    end if

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

      ! Add contribution for this particular Voronoi cell
      do kk = 1, single_row_grid%n
        vj = single_row_grid%index_left( kk)
        if (vj == vi) then
          ! Add contribution to this Voronoi cell
          single_row_Vor%LI_xdy(   k) = single_row_Vor%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
          single_row_Vor%LI_mxydx( k) = single_row_Vor%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
          single_row_Vor%LI_xydy(  k) = single_row_Vor%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
          exit
        end if
      end do

      ! Add entries to the big matrices
      call add_entry_CSR_dist( A_xdy_a_g_CSR  , row, col, single_row_Vor%LI_xdy(   k))
      call add_entry_CSR_dist( A_mxydx_a_g_CSR, row, col, single_row_Vor%LI_mxydx( k))
      call add_entry_CSR_dist( A_xydy_a_g_CSR , row, col, single_row_Vor%LI_xydy(  k))

    end do

  end subroutine calc_A_matrices_small_vertex

  subroutine calc_A_matrices_large_vertex( grid, mesh, row, vi, &
    single_row_Vor, A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR)
    !< Calculate a single row of the A-matrices, for a large vertex
    !< (i.e., simply averaging over the grid cells rather than calculating the line integrals)

    ! In/output variables
    type(type_grid),                        intent(in   ) :: grid
    type(type_mesh),                        intent(in   ) :: mesh
    integer,                                intent(in   ) :: row, vi
    type(type_single_row_mapping_matrices), intent(inout) :: single_row_Vor
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR

    ! Local variables
    integer                                :: col
    real(dp), dimension( mesh%nC_mem,2)    :: Vor
    integer,  dimension( mesh%nC_mem  )    :: Vor_vi
    integer,  dimension( mesh%nC_mem  )    :: Vor_ti
    integer                                :: nVor
    real(dp), dimension(2)                 :: p
    integer                                :: k, i, j
    real(dp)                               :: xmin, xmax, ymin, ymax
    integer                                :: il, iu, jl, ju

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

    ! Check which of the grid cells in this square lie inside the Voronoi cell
    do i = il, iu
    do j = jl, ju

      col = grid%ij2n( i,j)
      p   = [grid%x( i), grid%y( j)]

      if (is_in_Voronoi_cell( mesh, p, vi)) then
        ! This grid cell lies inside the Voronoi cell; add it to the single row
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
      call add_entry_CSR_dist( A_xdy_a_g_CSR  , row, col, single_row_Vor%LI_xdy(   k))
      call add_entry_CSR_dist( A_mxydx_a_g_CSR, row, col, single_row_Vor%LI_mxydx( k))
      call add_entry_CSR_dist( A_xydy_a_g_CSR , row, col, single_row_Vor%LI_xydy(  k))
    end do

  end subroutine calc_A_matrices_large_vertex

  subroutine calc_w_matrices( grid, mesh, &
    lies_outside_grid_domain, is_large_vertex, &
    A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR, w0, w1x, w1y)
    !< Calculate the w-matrices for the grid-to-mesh remapping operator

    ! In/output variables
    type(type_grid),                       intent(in   ) :: grid
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: lies_outside_grid_domain
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: is_large_vertex
    type(type_sparse_matrix_CSR_dp),       intent(in   ) :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR
    type(tMat),                            intent(  out) :: w0, w1x, w1y

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_w_matrices'
    integer                         :: nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc
    integer                         :: vi
    integer                         :: k, i, j
    type(type_sparse_matrix_CSR_dp) :: w0_CSR, w1x_CSR, w1y_CSR
    integer                         :: row, k1, k2, col
    real(dp)                        :: A_overlap_tot

    ! Add routine to path
    call init_routine( routine_name)

    nrows        = A_xdy_a_g_CSR%m
    ncols        = A_xdy_a_g_CSR%n
    nrows_loc    = A_xdy_a_g_CSR%m_loc
    ncols_loc    = A_xdy_a_g_CSR%n_loc
    nnz_est_proc = A_xdy_a_g_CSR%nnz

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      if (lies_outside_grid_domain( vi)) then
        ! Skip these
        call add_empty_row_CSR_dist( w0_CSR,  row)
        call add_empty_row_CSR_dist( w1x_CSR, row)
        call add_empty_row_CSR_dist( w1y_CSR, row)
        cycle
      end if

      k1 = A_xdy_a_g_CSR%ptr( row  )
      k2 = A_xdy_a_g_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_a_g_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_a_g_CSR%ind( k)
        call add_entry_CSR_dist( w0_CSR, row, col, A_xdy_a_g_CSR%val( k) / A_overlap_tot)
      end do

      if (.not. is_large_vertex( vi)) then
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
    character(len=1024), parameter  :: routine_name = 'calc_remapping_matrix'
    type(PetscErrorCode)            :: perr
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

end module remapping_grid_to_mesh_vertices
