module remapping_mesh_to_mesh

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, finalise_matrix_CSR_dist, &
    add_empty_row_CSR_dist, add_entry_CSR_dist, deallocate_matrix_CSR_dist
  use plane_geometry, only: triangle_area
  use mesh_utilities, only: calc_Voronoi_cell, find_containing_triangle, find_containing_vertex
  use petsc_basic, only: mat_CSR2petsc, mat_petsc2CSR
  use line_tracing_triangles, only: trace_line_tri
  use line_tracing_Voronoi, only: trace_line_Vor
  use netcdf_output

  implicit none

  private

  public :: create_map_from_mesh_to_mesh_nearest_neighbour, create_map_from_mesh_to_mesh_trilin, &
    create_map_from_mesh_to_mesh_2nd_order_conservative, &
    create_map_from_mesh_tri_to_mesh_tri_2nd_order_conservative

contains

  !> Create a new mapping object from a mesh to a mesh using nearest-neighbour interpolation.
  subroutine create_map_from_mesh_to_mesh_nearest_neighbour( mesh_src, mesh_dst, output_dir, map)

    ! In/output variables
    type(type_mesh),  intent(in   ) :: mesh_src
    type(type_mesh),  intent(in   ) :: mesh_dst
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_map_from_mesh_to_mesh_nearest_neighbour'
    integer                         :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp) :: M_CSR
    integer                         :: row, vi_dst
    real(dp), dimension(2)          :: p
    integer                         :: vi_src, col
    character(len=1024)             :: filename_mesh_src, filename_mesh_dst
    integer                         :: stat

    ! Add routine to path
    call init_routine( routine_name)

    ! Dump the two meshes to NetCDF. If the remapping crashes, having these available will
    ! help Tijn to find the error. If not, then they will be deleted at the end of this routine.
    filename_mesh_src = trim(output_dir) // '/mesh2mesh_nn_mesh_src_dump.nc'
    filename_mesh_dst = trim(output_dir) // '/mesh2mesh_nn_mesh_dst_dump.nc'
    call save_mesh_as_netcdf( filename_mesh_src, mesh_src)
    call save_mesh_as_netcdf( filename_mesh_dst, mesh_dst)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = 'nearest_neighbour'

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    nrows_loc       = mesh_dst%nV_loc
    ncols           = mesh_src%nV   ! from
    ncols_loc       = mesh_src%nV_loc
    nnz_per_row_max = 1
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src vertex containing them
    vi_src = 1
    do row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      call find_containing_vertex( mesh_src, p, vi_src)

      col = mesh_src%vi2n( vi_src)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, row, col, 1._dp)

    end do

    call finalise_matrix_CSR_dist( M_CSR)

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Delete mesh netcdf dumps
    if (par%primary) then
      open(unit = 1234, iostat = stat, file = filename_mesh_src, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
      open(unit = 1234, iostat = stat, file = filename_mesh_dst, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_nearest_neighbour

  !> Create a new mapping object from a mesh to a mesh using trilinear interpolation.
  subroutine create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, output_dir, map)

    ! In/output variables
    type(type_mesh),  intent(in   ) :: mesh_src
    type(type_mesh),  intent(in   ) :: mesh_dst
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_map_from_mesh_to_mesh_trilin'
    integer                         :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp) :: M_CSR
    integer                         :: row, vi_dst
    real(dp), dimension(2)          :: p
    integer                         :: ti_src, via, vib, vic
    real(dp), dimension(2)          :: pa, pb, pc
    real(dp)                        :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc, wa, wb, wc
    integer                         :: cola, colb, colc
    character(len=1024)             :: filename_mesh_src, filename_mesh_dst
    integer                         :: stat

    ! Add routine to path
    call init_routine( routine_name)

    ! Dump the two meshes to NetCDF. If the remapping crashes, having these available will
    ! help Tijn to find the error. If not, then they will be deleted at the end of this routine.
    filename_mesh_src = trim(output_dir) // '/mesh2mesh_trilin_mesh_src_dump.nc'
    filename_mesh_dst = trim(output_dir) // '/mesh2mesh_trilin_mesh_dst_dump.nc'
    call save_mesh_as_netcdf( filename_mesh_src, mesh_src)
    call save_mesh_as_netcdf( filename_mesh_dst, mesh_dst)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = 'trilin'

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    nrows_loc       = mesh_dst%nV_loc
    ncols           = mesh_src%nV   ! from
    ncols_loc       = mesh_src%nV_loc
    nnz_per_row_max = 3
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    ti_src = 1
    do row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      call find_containing_triangle( mesh_src, p, ti_src)

      ! Calculate the trilinear interpolation weights
      via = mesh_src%Tri( ti_src,1)
      vib = mesh_src%Tri( ti_src,2)
      vic = mesh_src%Tri( ti_src,3)

      pa  = mesh_src%V( via,:)
      pb  = mesh_src%V( vib,:)
      pc  = mesh_src%V( vic,:)

      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_abc = Atri_abp + Atri_bcp + Atri_cap

      wa = Atri_bcp / Atri_abc
      wb = Atri_cap / Atri_abc
      wc = Atri_abp / Atri_abc

      ! Matrix columns corresponding to these three vertices
      cola = mesh_src%vi2n( via)
      colb = mesh_src%vi2n( vib)
      colc = mesh_src%vi2n( vic)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, row, cola, wa)
      call add_entry_CSR_dist( M_CSR, row, colb, wb)
      call add_entry_CSR_dist( M_CSR, row, colc, wc)

    end do

    call finalise_matrix_CSR_dist( M_CSR)

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Delete mesh netcdf dumps
    if (par%primary) then
      open(unit = 1234, iostat = stat, file = filename_mesh_src, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
      open(unit = 1234, iostat = stat, file = filename_mesh_dst, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_trilin

  !> Create a new mapping object from a mesh to a mesh using 2nd-order conservative interpolation.
  subroutine create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, output_dir, map)

    ! In/output variables
    type(type_mesh),  intent(in   ) :: mesh_src
    type(type_mesh),  intent(in   ) :: mesh_dst
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'create_map_from_mesh_to_mesh_2nd_order_conservative'
    type(tMat)                          :: A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b
    type(tMat)                          :: w0, w1x, w1y
    type(tMat)                          :: M_cons_1st_order
    character(len=1024)                 :: filename_mesh_src, filename_mesh_dst
    integer                             :: stat

    ! Add routine to path
    call init_routine( routine_name)

    ! Dump the two meshes to NetCDF. If the remapping crashes, having these available will
    ! help Tijn to find the error. If not, then they will be deleted at the end of this routine.
    filename_mesh_src = trim( output_dir) // '/mesh_src'
    filename_mesh_dst = trim( output_dir) // '/mesh_dst'
    call save_mesh_as_netcdf( filename_mesh_src, mesh_src)
    call save_mesh_as_netcdf( filename_mesh_dst, mesh_dst)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = '2nd_order_conservative'

    call calc_A_matrices( mesh_src, mesh_dst, A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b)

    call calc_w_matrices( mesh_src, mesh_dst, A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b, w0, w1x, w1y)

    call calc_remapping_matrix( mesh_src, w0, w1x, w1y, M_cons_1st_order, map%M)

    call correct_mesh_to_mesh_map( mesh_src, mesh_dst, output_dir, M_cons_1st_order, map%M)

    ! Delete mesh netcdf dumps
    if (par%primary) then
      open(unit = 1234, iostat = stat, file = filename_mesh_src, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
      open(unit = 1234, iostat = stat, file = filename_mesh_dst, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_2nd_order_conservative

  subroutine create_map_from_mesh_tri_to_mesh_tri_2nd_order_conservative( mesh_src, mesh_dst, output_dir, map)

    ! In/output variables
    type(type_mesh),  intent(in   ) :: mesh_src
    type(type_mesh),  intent(in   ) :: mesh_dst
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'create_map_from_mesh_tri_to_mesh_tri_2nd_order_conservative'
    type(tMat)                          :: A_xdy, A_mxydx, A_xydy
    type(tMat)                          :: w0, w1x, w1y
    type(tMat)                          :: M_cons_1st_order
    character(len=1024)                 :: filename_mesh_src, filename_mesh_dst
    integer                             :: stat

    ! Add routine to path
    call init_routine( routine_name)

    ! Dump the two meshes to NetCDF. If the remapping crashes, having these available will
    ! help Tijn to find the error. If not, then they will be deleted at the end of this routine.
    filename_mesh_src = trim(output_dir) // '/mesh2mesh_tri_cons_mesh_src_dump.nc'
    filename_mesh_dst = trim(output_dir) // '/mesh2mesh_tri_cons_mesh_dst_dump.nc'
    call save_mesh_as_netcdf( filename_mesh_src, mesh_src)
    call save_mesh_as_netcdf( filename_mesh_dst, mesh_dst)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = trim(mesh_src%name) // '_triangles'
    map%name_dst  = trim(mesh_dst%name) // '_triangles'
    map%method    = '2nd_order_conservative'

    call calc_A_matrices_tri( mesh_src, mesh_dst, A_xdy, A_mxydx, A_xydy)

    call calc_w_matrices_tri( mesh_src, mesh_dst, A_xdy, A_mxydx, A_xydy, w0, w1x, w1y)

    call calc_remapping_matrix_tri( mesh_src, w0, w1x, w1y, M_cons_1st_order, map%M)

    ! call correct_mesh_to_mesh_map( mesh_src, mesh_dst, M_cons_1st_order, map%M)

    ! Delete mesh netcdf dumps
    if (par%primary) then
      open(unit = 1234, iostat = stat, file = filename_mesh_src, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
      open(unit = 1234, iostat = stat, file = filename_mesh_dst, status = 'old')
      if (stat == 0) close(1234, status = 'delete')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_tri_to_mesh_tri_2nd_order_conservative

  subroutine calc_A_matrices( mesh_src, mesh_dst, A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b)
    !< Calculate the A-matrices for the mesh-to-mesh remapping operator

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh_src
    type(type_mesh), intent(in   ) :: mesh_dst
    type(tMat),      intent(  out) :: A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_A_matrices'
    logical                        :: count_coincidences
    type(tMat)                     :: A_xdy_b_a  , A_mxydx_b_a  , A_xydy_b_a
    type(tMat)                     :: A_xdy_b_a_T, A_mxydx_b_a_T, A_xydy_b_a_T
    type(PetscErrorCode)           :: perr

    ! Add routine to path
    call init_routine( routine_name)

    ! Integrate around the Voronoi cells of the destination mesh through the triangles of the source mesh
    count_coincidences = .true.
    call integrate_Voronoi_cells_through_triangles( mesh_dst, mesh_src, A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b, count_coincidences)

    ! Integrate around the triangles of the source mesh through the Voronoi cells of the destination mesh
    count_coincidences = .false.
    call integrate_triangles_through_Voronoi_cells( mesh_src, mesh_dst, A_xdy_b_a, A_mxydx_b_a, A_xydy_b_a, count_coincidences)

    ! Transpose line integral matrices
    call MatCreateTranspose( A_xdy_b_a  , A_xdy_b_a_T  , perr)
    call MatCreateTranspose( A_mxydx_b_a, A_mxydx_b_a_T, perr)
    call MatCreateTranspose( A_xydy_b_a , A_xydy_b_a_T , perr)

    ! Combine line integrals around areas of overlap to get surface integrals over areas of overlap
    call MatAXPY( A_xdy_a_b  , 1._dp, A_xdy_b_a_T  , UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( A_mxydx_a_b, 1._dp, A_mxydx_b_a_T, UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( A_xydy_a_b , 1._dp, A_xydy_b_a_T , UNKNOWN_NONZERO_PATTERN, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_A_matrices

  subroutine calc_A_matrices_tri( mesh_src, mesh_dst, A_xdy, A_mxydx, A_xydy)
    !< Calculate the A-matrices for the mesh-to-mesh triangles remapping operator

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh_src
    type(type_mesh), intent(in   ) :: mesh_dst
    type(tMat),      intent(  out) :: A_xdy, A_mxydx, A_xydy

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_A_matrices_tri'
    logical                        :: count_coincidences
    type(tMat)                     :: A_xdy_src_dst  , A_mxydx_src_dst  , A_xydy_src_dst
    type(tMat)                     :: A_xdy_src_dst_T, A_mxydx_src_dst_T, A_xydy_src_dst_T
    type(PetscErrorCode)           :: perr

    ! Add routine to path
    call init_routine( routine_name)

    ! Integrate around the triangles of the destination mesh through the triangles of the source mesh
    count_coincidences = .true.
    call integrate_triangles_through_triangles( mesh_dst, mesh_src, A_xdy, A_mxydx, A_xydy, count_coincidences)

    ! Integrate around the triangles of the source mesh through the triangles of the destination mesh
    count_coincidences = .false.
    call integrate_triangles_through_triangles( mesh_src, mesh_dst, A_xdy_src_dst, A_mxydx_src_dst, A_xydy_src_dst, count_coincidences)

    ! Transpose line integral matrices
    call MatCreateTranspose( A_xdy_src_dst  , A_xdy_src_dst_T  , perr)
    call MatCreateTranspose( A_mxydx_src_dst, A_mxydx_src_dst_T, perr)
    call MatCreateTranspose( A_xydy_src_dst , A_xydy_src_dst_T , perr)

    ! Combine line integrals around areas of overlap to get surface integrals over areas of overlap
    call MatAXPY( A_xdy  , 1._dp, A_xdy_src_dst_T  , UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( A_mxydx, 1._dp, A_mxydx_src_dst_T, UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( A_xydy , 1._dp, A_xydy_src_dst_T , UNKNOWN_NONZERO_PATTERN, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_A_matrices_tri

  subroutine calc_w_matrices( mesh_src, mesh_dst, A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b, w0, w1x, w1y)
    !< Calculate the w-matrices for the mesh-to-mesh remapping operator

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh_src
    type(type_mesh), intent(in   ) :: mesh_dst
    type(tMat),      intent(in   ) :: A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b
    type(tMat),      intent(  out) :: w0, w1x, w1y

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'calc_w_matrices'
    type(PetscErrorCode)                        :: perr
    integer                                     :: nnz_per_row_max
    integer                                     :: istart, iend, n, k, ti
    integer                                     :: ncols
    integer,  dimension(:), allocatable, target :: cols
    real(dp), dimension(:), allocatable, target :: vals, w0_row, w1x_row, w1y_row
    integer,  dimension(:), pointer             :: cols_
    real(dp), dimension(:), pointer             :: vals_, w0_row_, w1x_row_, w1y_row_
    real(dp)                                    :: A_overlap_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    call MatConvert( A_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w0, perr)
    call MatConvert( A_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1x, perr)
    call MatConvert( A_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1y, perr)

    ! Estimate maximum number of non-zeros per row (i.e. maximum number of grid cells overlapping with a mesh triangle)
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_src%TriA) / minval( mesh_dst%A   )), &
                                    ceiling( 2._dp * maxval( mesh_dst%A   ) / minval( mesh_src%TriA))) )

    ! allocate memory for a single matrix row
    allocate( cols(    nnz_per_row_max))
    allocate( vals(    nnz_per_row_max))
    allocate( w0_row(  nnz_per_row_max))
    allocate( w1x_row( nnz_per_row_max))
    allocate( w1y_row( nnz_per_row_max))

    cols_    => cols
    vals_    => vals
    w0_row_  => w0_row
    w1x_row_ => w1x_row
    w1y_row_ => w1y_row

    call MatGetOwnershipRange( A_xdy_a_b  , istart, iend, perr)

    do n = istart+1, iend ! +1 because PETSc indexes from 0

      ! Calculate area of overlap
      call MatGetRow( A_xdy_a_b, n-1, ncols, cols_, vals_, perr)
      A_overlap_tot = sum( vals_( 1:ncols))
      call MatRestoreRow( A_xdy_a_b, n-1, ncols, cols_, vals_, perr)

      ! Skip vertices with zero overlap (which can happen if the boundary
      ! of their Voronoi cell coincides with that of this one)
      if (A_overlap_tot <= tiny( A_overlap_tot) * 16._dp) cycle

      ! w0
      call MatGetRow( A_xdy_a_b, n-1, ncols, cols_, vals_, perr)
      do k = 1, ncols
        w0_row_( k) = vals_( k) / A_overlap_tot
        call MatSetValues( w0, 1, [n-1], 1, [cols_( k)], [w0_row_( k)], INSERT_VALUES, perr)
      end do
      call MatRestoreRow( A_xdy_a_b, n-1, ncols, cols_, vals_, perr)

      ! w1x
      call MatGetRow( A_mxydx_a_b, n-1, ncols, cols_, vals_, perr)
      do k = 1, ncols
        ti = cols_( k)+1
        w1x_row_( k) = (vals_( k) / A_overlap_tot) - (mesh_src%TriGC( ti,1) * w0_row( k))
        call MatSetValues( w1x, 1, [n-1], 1, [cols_( k)], [w1x_row_( k)], INSERT_VALUES, perr)
      end do
      call MatRestoreRow( A_mxydx_a_b, n-1, ncols, cols_, vals_, perr)

      ! w1y
      call MatGetRow( A_xydy_a_b, n-1, ncols, cols_, vals_, perr)
      do k = 1, ncols
        ti = cols_( k)+1
        w1y_row_( k) = (vals_( k) / A_overlap_tot) - (mesh_src%TriGC( ti,2) * w0_row( k))
        call MatSetValues( w1y, 1, [n-1], 1, [cols_( k)], [w1y_row_( k)], INSERT_VALUES, perr)
      end do
      call MatRestoreRow( A_xydy_a_b, n-1, ncols, cols_, vals_, perr)

    end do

    call MatAssemblyBegin( w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, perr)

    deallocate( cols)
    deallocate( vals)
    deallocate( w0_row)
    deallocate( w1x_row)
    deallocate( w1y_row)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_w_matrices

  subroutine calc_w_matrices_tri( mesh_src, mesh_dst, A_xdy, A_mxydx, A_xydy, w0, w1x, w1y)
    !< Calculate the w-matrices for the mesh-to-mesh triangles remapping operator

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh_src
    type(type_mesh), intent(in   ) :: mesh_dst
    type(tMat),      intent(in   ) :: A_xdy, A_mxydx, A_xydy
    type(tMat),      intent(  out) :: w0, w1x, w1y

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'calc_w_matrices_tri'
    type(PetscErrorCode)                        :: perr
    integer                                     :: nnz_per_row_max
    integer                                     :: istart, iend, n, k, ti
    integer                                     :: ncols
    integer,  dimension(:), allocatable, target :: cols
    real(dp), dimension(:), allocatable, target :: vals, w0_row, w1x_row, w1y_row
    integer,  dimension(:), pointer             :: cols_
    real(dp), dimension(:), pointer             :: vals_, w0_row_, w1x_row_, w1y_row_
    real(dp)                                    :: A_overlap_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    call MatConvert( A_xdy, MATAIJ, MAT_INITIAL_MATRIX, w0, perr)
    call MatConvert( A_xdy, MATAIJ, MAT_INITIAL_MATRIX, w1x, perr)
    call MatConvert( A_xdy, MATAIJ, MAT_INITIAL_MATRIX, w1y, perr)

    ! Estimate maximum number of non-zeros per row (i.e. maximum number of grid cells overlapping with a mesh triangle)
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_src%TriA) / minval( mesh_dst%A   )), &
                                    ceiling( 2._dp * maxval( mesh_dst%A   ) / minval( mesh_src%TriA))) )

    ! allocate memory for a single matrix row
    allocate( cols(    nnz_per_row_max))
    allocate( vals(    nnz_per_row_max))
    allocate( w0_row(  nnz_per_row_max))
    allocate( w1x_row( nnz_per_row_max))
    allocate( w1y_row( nnz_per_row_max))

    cols_    => cols
    vals_    => vals
    w0_row_  => w0_row
    w1x_row_ => w1x_row
    w1y_row_ => w1y_row

    call MatGetOwnershipRange( A_xdy  , istart, iend, perr)

    do n = istart+1, iend ! +1 because PETSc indexes from 0

      ! Calculate area of overlap
      call MatGetRow( A_xdy, n-1, ncols, cols_, vals_, perr)
      A_overlap_tot = sum( vals_( 1:ncols))
      call MatRestoreRow( A_xdy, n-1, ncols, cols_, vals_, perr)

      ! Skip vertices with zero overlap (which can happen if the boundary
      ! of their Voronoi cell coincides with that of this one)
      if (A_overlap_tot <= tiny( A_overlap_tot) * 16._dp) cycle

      ! w0
      call MatGetRow( A_xdy, n-1, ncols, cols_, vals_, perr)
      do k = 1, ncols
        w0_row_( k) = vals_( k) / A_overlap_tot
        call MatSetValues( w0, 1, [n-1], 1, [cols_( k)], [w0_row_( k)], INSERT_VALUES, perr)
      end do
      call MatRestoreRow( A_xdy, n-1, ncols, cols_, vals_, perr)

      ! w1x
      call MatGetRow( A_mxydx, n-1, ncols, cols_, vals_, perr)
      do k = 1, ncols
        ti = cols_( k)+1
        w1x_row_( k) = (vals_( k) / A_overlap_tot) - (mesh_src%TriGC( ti,1) * w0_row( k))
        call MatSetValues( w1x, 1, [n-1], 1, [cols_( k)], [w1x_row_( k)], INSERT_VALUES, perr)
      end do
      call MatRestoreRow( A_mxydx, n-1, ncols, cols_, vals_, perr)

      ! w1y
      call MatGetRow( A_xydy, n-1, ncols, cols_, vals_, perr)
      do k = 1, ncols
        ti = cols_( k)+1
        w1y_row_( k) = (vals_( k) / A_overlap_tot) - (mesh_src%TriGC( ti,2) * w0_row( k))
        call MatSetValues( w1y, 1, [n-1], 1, [cols_( k)], [w1y_row_( k)], INSERT_VALUES, perr)
      end do
      call MatRestoreRow( A_xydy, n-1, ncols, cols_, vals_, perr)

    end do

    call MatAssemblyBegin( w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, perr)

    deallocate( cols)
    deallocate( vals)
    deallocate( w0_row)
    deallocate( w1x_row)
    deallocate( w1y_row)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_w_matrices_tri

  subroutine calc_remapping_matrix( mesh_src, w0, w1x, w1y, M_cons_1st_order, M)
    !< Calculate the mesh-to-mesh remapping matrix M

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh_src
    type(tMat),      intent(in   ) :: w0, w1x, w1y
    type(tMat),      intent(  out) :: M_cons_1st_order, M

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_remapping_matrix'
    type(PetscErrorCode)           :: perr
    type(tMat)                     :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    type(tMat)                     :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. allocated( mesh_src%vi2n)) then
      call crash('matrix operators for mesh "' // trim( mesh_src%name) // '" have not been calculated!')
    end if

    ! Convert matrices to PETSc format
    call mat_CSR2petsc( mesh_src%M_map_a_b, M_map_a_b)
    call mat_CSR2petsc( mesh_src%M_ddx_a_b, M_ddx_a_b)
    call mat_CSR2petsc( mesh_src%M_ddy_a_b, M_ddy_a_b)

    ! 1st-order = w0 * map_a_b
    call MatMatMult( w0 , M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M_cons_1st_order, perr)

    ! 2nd-order = 1st-order + w1x * ddx_a_b + w1y * ddy_a_b
    call MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatConvert( M_cons_1st_order, MATAIJ, MAT_INITIAL_MATRIX, M, perr)
    call MatAXPY( M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_remapping_matrix

  subroutine calc_remapping_matrix_tri( mesh_src, w0, w1x, w1y, M_cons_1st_order, M)
    !< Calculate the mesh-to-mesh triangles remapping matrix M

    ! In/output variables
    type(type_mesh), intent(in   ) :: mesh_src
    type(tMat),      intent(in   ) :: w0, w1x, w1y
    type(tMat),      intent(  out) :: M_cons_1st_order, M

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_remapping_matrix_tri'
    type(PetscErrorCode)           :: perr
    type(tMat)                     :: M_ddx_b_b, M_ddy_b_b
    type(tMat)                     :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. allocated( mesh_src%vi2n)) then
      call crash('matrix operators for mesh "' // trim( mesh_src%name) // '" have not been calculated!')
    end if

    ! Convert matrices to PETSc format
    call mat_CSR2petsc( mesh_src%M_ddx_b_b, M_ddx_b_b)
    call mat_CSR2petsc( mesh_src%M_ddy_b_b, M_ddy_b_b)

    ! 1st-order = w0
    call MatConvert( w0, MATAIJ, MAT_INITIAL_MATRIX, M_cons_1st_order, perr)

    ! 2nd-order = 1st-order + w1x * ddx_b_b + w1y * ddy_b_b
    call MatMatMult( w1x, M_ddx_b_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, M_ddy_b_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatConvert( M_cons_1st_order, MATAIJ, MAT_INITIAL_MATRIX, M, perr)
    call MatAXPY( M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_remapping_matrix_tri

  subroutine correct_mesh_to_mesh_map( mesh_src, mesh_dst, output_dir, M_cons_1st_order, M_cons_2nd_order)
    !< Apply some final corrections to the 2nd-order conservative mesh-to-mesh remapping operator
    !
    ! - set remapped data to zero on the domain border
    ! - use direct copying for identical vertices

    ! In/output variables
    type(type_mesh),  intent(in)    :: mesh_src
    type(type_mesh),  intent(in)    :: mesh_dst
    character(len=*), intent(in   ) :: output_dir
    type(tMat),       intent(in)    :: M_cons_1st_order
    type(tMat),       intent(inout) :: M_cons_2nd_order

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'correct_mesh_to_mesh_map'
    integer                                 :: perr
    type(type_sparse_matrix_CSR_dp)         :: M_cons_1st_order_CSR
    type(type_sparse_matrix_CSR_dp)         :: M_cons_2nd_order_CSR
    integer                                 :: i, vi_dst, k1, k2, k
    integer                                 :: j, vi_src
    logical                                 :: do_direct_copy
    integer                                 :: vi_src_copy
    logical                                 :: Voronoi_cells_are_identical
    real(dp), dimension( mesh_src%nC_mem,2) :: Vor_src   , Vor_dst
    integer,  dimension( mesh_src%nC_mem  ) :: Vor_src_vi, Vor_dst_vi
    integer,  dimension( mesh_src%nC_mem  ) :: Vor_src_ti, Vor_dst_ti
    integer                                 :: nVor_src  , nVor_dst
    integer                                 :: vori
    type(type_map)                          :: map_trilin
    type(type_sparse_matrix_CSR_dp)         :: M_trilin_CSR
    logical,  dimension( mesh_dst%nV)       :: isgood_1st_order
    logical,  dimension( mesh_dst%nV)       :: isgood_2nd_order
    integer                                 :: kk1,kk2,kk

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert matrices to CSR format for easier handling
    call mat_petsc2CSR( M_cons_1st_order, M_cons_1st_order_CSR)
    call mat_petsc2CSR( M_cons_2nd_order, M_cons_2nd_order_CSR)

    ! == Set to zero
    !
    ! 2nd-order conservative doesn't work all that well on the
    ! domain border; just set the result to zero there.

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      if (mesh_dst%VBI( vi_dst) > 0) then
        do k = k1, k2
          M_cons_2nd_order_CSR%val( k) = 0._dp
        end do
      end if

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! == Direct copying
    !
    ! With the new mesh generation code (of UFE2.0), many vertices away from the grounding line
    ! remain unchanged after a mesh update. The vertex-to-triangle-to-vertex remapping is slightly
    ! diffusive, so instead we can just copy data directly for those vertices.

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do_direct_copy = .false.
      vi_src_copy    = 0

      ! Loop over all src vertices contributing to this dst vertex
      do k = k1, k2

        j = M_cons_2nd_order_CSR%ind( k)
        vi_src = mesh_src%n2vi( j)

        if (norm2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%R( vi_dst) / 1E2_dp) then
          ! These vertices coincide; check if their Voronoi cells are identical

          Voronoi_cells_are_identical = .true.

          call calc_Voronoi_cell( mesh_src, vi_src, 0._dp, Vor_src, Vor_src_vi, Vor_src_ti, nVor_src)
          call calc_Voronoi_cell( mesh_dst, vi_dst, 0._dp, Vor_dst, Vor_dst_vi, Vor_dst_ti, nVor_dst)

          if (nVor_src /= nVor_dst) then
            Voronoi_cells_are_identical = .false.
          else
            do vori = 1, nVor_src
              if (norm2( Vor_src( vori,:) - Vor_dst( vori,:)) > mesh_dst%R( vi_dst) / 1E2_dp) then
                Voronoi_cells_are_identical = .false.
              end if
            end do
          end if

          if (Voronoi_cells_are_identical) then
            ! These two vertices have identical Voronoi cells; use direct copying
            do_direct_copy = .true.
            vi_src_copy    = vi_src
            exit
          end if ! if (Voronoi_cells_are_identical) then

        end if ! if (norm2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%tol_dist) then

      end do ! do k = k1, k2

      ! if a source vertex with an identical Voronoi cell to this dst vertex was
      ! found, copy data from that vertex directly
      if (do_direct_copy) then
        ! Loop over all src vertices contributing to this dst vertex; set all
        ! contributions to zero except the one we're copying (which is set to 1)

        do k = k1, k2

          j = M_cons_2nd_order_CSR%ind( k)
          vi_src = mesh_src%n2vi( j)

          if (vi_src == vi_src_copy) then
            M_cons_2nd_order_CSR%val( k) = 1._dp
          else
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end if

        end do ! do k = k1, k2

      end if ! if (do_direct_copy) then

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! == Remapping errors
    !
    ! On very rare occasions, the remapping operator is just wrong, likely due to round-off
    ! errors in determining if vertices coincide or not. Usually, the 1st-order operator
    ! is fine. if that one fails too, just replace the answer with trilinear interpolation.
    !
    ! Faulty operators can be detected by negative coefficients in the remapping matrix,
    ! which automatically violate conservation of extreme values.

    ! Calculate the trilinear interpolation operator to serve as a back-up
    call create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, output_dir, map_trilin)
    call mat_petsc2CSR( map_trilin%M, M_trilin_CSR)

    ! Find faulty operators in the 1st-order conservative remapping operator
    isgood_1st_order = .true.
    do i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

      k1 = M_cons_1st_order_CSR%ptr( i)
      k2 = M_cons_1st_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do k = k1, k2
        if (M_cons_1st_order_CSR%val( k) < 0._dp) then
          isgood_1st_order( vi_dst) = .false.
        end if
      end do

    end do ! do i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

    ! Find faulty operators in the 2nd-order conservative remapping operator
    isgood_2nd_order = .true.
    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do k = k1, k2
        if (M_cons_2nd_order_CSR%val( k) < 0._dp) then
          isgood_2nd_order( vi_dst) = .false.
        end if
      end do

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Replace faulty operators in the 2nd-order conservative remapping operator

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      vi_dst = mesh_dst%n2vi( i)

      if (.not. isgood_2nd_order( vi_dst)) then
        ! Replace this faulty operator

        if (isgood_1st_order( vi_dst)) then
          ! Replace with the 1st-order conservative remapping operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          do k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end do

          ! then copy values from the 1st-order operator
          kk1 = M_cons_1st_order_CSR%ptr( i)
          kk2 = M_cons_1st_order_CSR%ptr( i+1) - 1
          do kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_cons_1st_order_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_cons_1st_order_CSR%val( kk)
          end do

        else ! if (isgood_1st_order( vi_dst)) then
          ! Replace with the trilinear interpolation operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          do k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end do

          ! then copy values from the trilinear interpolation operator
          kk1 = M_trilin_CSR%ptr( i)
          kk2 = M_trilin_CSR%ptr( i+1) - 1
          do kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_trilin_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_trilin_CSR%val( kk)
          end do

        end if ! if (isgood_1st_order( vi_dst)) then

      end if ! if (.not. isgood_2nd_order( vi_dst)) then

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Convert back to PETSc format
    call MatDestroy( M_cons_2nd_order, perr)
    call mat_CSR2petsc( M_cons_2nd_order_CSR, M_cons_2nd_order)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( M_cons_1st_order_CSR)
    call deallocate_matrix_CSR_dist( M_cons_2nd_order_CSR)
    call deallocate_matrix_CSR_dist( M_trilin_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine correct_mesh_to_mesh_map

  !> Integrate around the triangles of mesh_tri through the Voronoi cells of mesh_Vor
  subroutine integrate_triangles_through_Voronoi_cells( mesh_tri, mesh_Vor, &
    A_xdy_b_a, A_mxydx_b_a, A_xydy_b_a, count_coincidences)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_tri
    type(type_mesh), intent(in)    :: mesh_Vor
    type(tMat),      intent(out)   :: A_xdy_b_a
    type(tMat),      intent(out)   :: A_mxydx_b_a
    type(tMat),      intent(out)   :: A_xydy_b_a
    logical,         intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'integrate_triangles_through_Voronoi_cells'
    integer                                :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)        :: A_xdy_b_a_CSR, A_mxydx_b_a_CSR, A_xydy_b_a_CSR
    type(type_single_row_mapping_matrices) :: single_row
    integer                                :: via, vib, vic, ti, vi_hint, k
    real(dp), dimension(2)                 :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_tri%nTri  ! to
    nrows_loc       = mesh_tri%nTri_loc
    ncols           = mesh_Vor%nV    ! from
    ncols_loc       = mesh_Vor%nV_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_tri%TriA) / minval( mesh_Vor%A   )), &
                                    ceiling( 2._dp * maxval( mesh_Vor%A   ) / minval( mesh_tri%TriA)) ))

    call allocate_matrix_CSR_dist( A_xdy_b_a_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_b_a_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_b_a_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! == Trace all the line segments to fill the matrices
    ! ===================================================

    vi_hint = 1

    do ti = mesh_tri%ti1, mesh_tri%ti2

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! The three vertices spanning this triangle
      via = mesh_tri%Tri( ti,1)
      vib = mesh_tri%Tri( ti,2)
      vic = mesh_tri%Tri( ti,3)

      ! Integrate over the three triangle sides
      p = mesh_tri%V( via,:)
      q = mesh_tri%V( vib,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vib,:)
      q = mesh_tri%V( vic,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vic,:)
      q = mesh_tri%V( via,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      ! Add the results for this triangle to the sparse matrix
      if (single_row%n == 0) then
        call add_empty_row_CSR_dist( A_xdy_b_a_CSR  , ti)
        call add_empty_row_CSR_dist( A_mxydx_b_a_CSR, ti)
        call add_empty_row_CSR_dist( A_xydy_b_a_CSR , ti)
      else
        do k = 1, single_row%n
          call add_entry_CSR_dist( A_xdy_b_a_CSR  , ti, single_row%index_left( k), single_row%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_b_a_CSR, ti, single_row%index_left( k), single_row%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_b_a_CSR , ti, single_row%index_left( k), single_row%LI_xydy(  k))
        end do
      end if

    end do

    call finalise_matrix_CSR_dist( A_xdy_b_a_CSR  )
    call finalise_matrix_CSR_dist( A_mxydx_b_a_CSR)
    call finalise_matrix_CSR_dist( A_xydy_b_a_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( A_xdy_b_a_CSR  , A_xdy_b_a  )
    call mat_CSR2petsc( A_mxydx_b_a_CSR, A_mxydx_b_a)
    call mat_CSR2petsc( A_xydy_b_a_CSR , A_xydy_b_a )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_triangles_through_Voronoi_cells

  !> Integrate around the grid cells of the grid through the triangles of the mesh
  subroutine integrate_Voronoi_cells_through_triangles( mesh_Vor, mesh_tri, &
    A_xdy_a_b, A_mxydx_a_b, A_xydy_a_b, count_coincidences)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_Vor
    type(type_mesh), intent(in)    :: mesh_tri
    type(tMat),      intent(out)   :: A_xdy_a_b
    type(tMat),      intent(out)   :: A_mxydx_a_b
    type(tMat),      intent(out)   :: A_xydy_a_b
    logical,         intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'integrate_Voronoi_cells_through_triangles'
    integer                                 :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)         :: A_xdy_a_b_CSR, A_mxydx_a_b_CSR, A_xydy_a_b_CSR
    type(type_single_row_mapping_matrices)  :: single_row
    integer                                 :: vi, vori1, vori2, k, ti_hint
    real(dp), dimension( mesh_Vor%nC_mem,2) :: Vor
    integer,  dimension( mesh_Vor%nC_mem  ) :: Vor_vi
    integer,  dimension( mesh_Vor%nC_mem  ) :: Vor_ti
    integer                                 :: nVor
    real(dp), dimension(2)                  :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix size
    nrows           = mesh_Vor%nV    ! to
    nrows_loc       = mesh_Vor%nV_loc
    ncols           = mesh_tri%nTri  ! from
    ncols_loc       = mesh_tri%nTri_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_tri%TriA) / minval( mesh_vor%A   )), &
                                    ceiling( 2._dp * maxval( mesh_vor%A   ) / minval( mesh_tri%TriA)) ))

    call allocate_matrix_CSR_dist( A_xdy_a_b_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_a_b_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_a_b_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! == Trace all the line segments to fill the matrices
    ! ===================================================

    ti_hint = 1

    do vi = mesh_Vor%vi1, mesh_Vor%vi2

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! Integrate over the complete Voronoi cell boundary
      call calc_Voronoi_cell( mesh_Vor, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)
      do vori1 = 1, nVor
        vori2 = vori1 + 1
        if (vori2 > nVor) vori2 = 1
        p = Vor( vori1,:)
        q = Vor( vori2,:)
        call trace_line_tri( mesh_tri, p, q, single_row, count_coincidences, ti_hint)
      end do

      ! Add the results for this triangle to the sparse matrix
      if (single_row%n == 0) then
        call add_empty_row_CSR_dist( A_xdy_a_b_CSR  , vi)
        call add_empty_row_CSR_dist( A_mxydx_a_b_CSR, vi)
        call add_empty_row_CSR_dist( A_xydy_a_b_CSR , vi)
      else
        do k = 1, single_row%n
          call add_entry_CSR_dist( A_xdy_a_b_CSR  , vi, single_row%index_left( k), single_row%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_a_b_CSR, vi, single_row%index_left( k), single_row%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_a_b_CSR , vi, single_row%index_left( k), single_row%LI_xydy(  k))
        end do
      end if

    end do

    call finalise_matrix_CSR_dist( A_xdy_a_b_CSR  )
    call finalise_matrix_CSR_dist( A_mxydx_a_b_CSR)
    call finalise_matrix_CSR_dist( A_xydy_a_b_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( A_xdy_a_b_CSR  , A_xdy_a_b  )
    call mat_CSR2petsc( A_mxydx_a_b_CSR, A_mxydx_a_b)
    call mat_CSR2petsc( A_xydy_a_b_CSR , A_xydy_a_b )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_Voronoi_cells_through_triangles

  !> Integrate around the triangles of mesh_top through the triangles of mesh_bot
  subroutine integrate_triangles_through_triangles( mesh_top, mesh_bot, &
    A_xdy_top_bot, A_mxydx_top_bot, A_xydy_top_bot, count_coincidences)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh_top
    type(type_mesh), intent(in)    :: mesh_bot
    type(tMat),      intent(out)   :: A_xdy_top_bot
    type(tMat),      intent(out)   :: A_mxydx_top_bot
    type(tMat),      intent(out)   :: A_xydy_top_bot
    logical,         intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'integrate_triangles_through_triangles'
    integer                                :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)        :: A_xdy_top_bot_CSR, A_mxydx_top_bot_CSR, A_xydy_top_bot_CSR
    type(type_single_row_mapping_matrices) :: single_row
    integer                                :: via, vib, vic, ti_top, ti_bot_hint, k
    real(dp), dimension(2)                 :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_top%nTri  ! to
    nrows_loc       = mesh_top%nTri_loc
    ncols           = mesh_bot%nTri    ! from
    ncols_loc       = mesh_bot%nTri_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_top%TriA) / minval( mesh_bot%TriA)), &
                                    ceiling( 2._dp * maxval( mesh_bot%TriA) / minval( mesh_top%TriA)) ))

    call allocate_matrix_CSR_dist( A_xdy_top_bot_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_top_bot_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_top_bot_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single top-mesh triangle through the bottom mesh's triangles
    single_row%n_max = 100
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! == Trace all the line segments to fill the matrices
    ! ===================================================

    ti_bot_hint = 1

    do ti_top = mesh_top%ti1, mesh_top%ti2

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! The three vertices spanning this top-mesh triangle
      via = mesh_top%Tri( ti_top,1)
      vib = mesh_top%Tri( ti_top,2)
      vic = mesh_top%Tri( ti_top,3)

      ! Integrate over the three top-mesh triangle sides
      p = mesh_top%V( via,:)
      q = mesh_top%V( vib,:)
      call trace_line_tri( mesh_bot, p, q, single_row, count_coincidences, ti_bot_hint)

      p = mesh_top%V( vib,:)
      q = mesh_top%V( vic,:)
      call trace_line_tri( mesh_bot, p, q, single_row, count_coincidences, ti_bot_hint)

      p = mesh_top%V( vic,:)
      q = mesh_top%V( via,:)
      call trace_line_tri( mesh_bot, p, q, single_row, count_coincidences, ti_bot_hint)

      ! Add the results for this triangle to the sparse matrix
      if (single_row%n == 0) then
        call add_empty_row_CSR_dist( A_xdy_top_bot_CSR  , ti_top)
        call add_empty_row_CSR_dist( A_mxydx_top_bot_CSR, ti_top)
        call add_empty_row_CSR_dist( A_xydy_top_bot_CSR , ti_top)
      else
        do k = 1, single_row%n
          call add_entry_CSR_dist( A_xdy_top_bot_CSR  , ti_top, single_row%index_left( k), single_row%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_top_bot_CSR, ti_top, single_row%index_left( k), single_row%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_top_bot_CSR , ti_top, single_row%index_left( k), single_row%LI_xydy(  k))
        end do
      end if

    end do

    call finalise_matrix_CSR_dist( A_xdy_top_bot_CSR  )
    call finalise_matrix_CSR_dist( A_mxydx_top_bot_CSR)
    call finalise_matrix_CSR_dist( A_xydy_top_bot_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( A_xdy_top_bot_CSR  , A_xdy_top_bot  )
    call mat_CSR2petsc( A_mxydx_top_bot_CSR, A_mxydx_top_bot)
    call mat_CSR2petsc( A_xydy_top_bot_CSR , A_xydy_top_bot )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_triangles_through_triangles

end module remapping_mesh_to_mesh
