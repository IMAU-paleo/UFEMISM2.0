module remapping_mesh_triangles_to_grid

#include <petsc/finclude/petscksp.h>
  use petscksp
  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use remapping_types, only: type_map
  use petsc_basic, only: mat_CSR2petsc
  use remapping_mesh_vertices_to_grid, only: calc_approximate_overlaps, calc_A_matrices, &
    calc_w_matrices, dump_grid_and_mesh_to_netcdf, delete_grid_and_mesh_netcdf_dump_files

  implicit none

  private

  public :: create_map_from_mesh_triangles_to_xy_grid

contains

  subroutine create_map_from_mesh_triangles_to_xy_grid( mesh, grid, output_dir, map)
    !< Create a new mapping object from a mesh to an x/y-grid.

    ! By default uses 1st-order conservative remapping.

    ! In/output variables
    type(type_mesh),  intent(in   ) :: mesh
    type(type_grid),  intent(in   ) :: grid
    character(len=*), intent(in   ) :: output_dir
    type(type_map),   intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'create_map_from_mesh_triangles_to_xy_grid'
    integer, dimension(grid%nx, grid%ny)   :: overlaps_with_small_triangle, containing_triangle
    type(type_sparse_matrix_CSR_dp)        :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR
    type(tMat)                             :: w1x, w1y
    character(len=1024)                    :: filename_grid, filename_mesh

    ! Add routine to path
    call init_routine( routine_name)

    call dump_grid_and_mesh_to_netcdf( grid, mesh, output_dir, filename_grid, filename_mesh)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = trim(mesh%name) // '_triangles'
    map%name_dst  = grid%name
    map%method    = '2nd_order_conservative'

    call calc_approximate_overlaps( mesh, grid, &
      overlaps_with_small_triangle, containing_triangle)

    call calc_A_matrices( mesh, grid, &
      overlaps_with_small_triangle, containing_triangle, &
      A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR)

    call calc_w_matrices( mesh, grid, &
      A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR, map%M, w1x, w1y)

    call delete_grid_and_mesh_netcdf_dump_files( filename_grid, filename_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_triangles_to_xy_grid

end module remapping_mesh_triangles_to_grid
