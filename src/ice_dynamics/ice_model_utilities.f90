module ice_model_utilities
  !< Generally useful functions used by the ice model.

  use mpi
  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use model_configuration, only: C
  use parameters, only: ice_density, seawater_density
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use remapping_main, only: Atlas
  use SMB_model_types, only: type_SMB_model
  use BMB_model_types, only: type_BMB_model
  use LMB_model_types, only: type_LMB_model
  use AMB_model_types, only: type_AMB_model
  use netcdf_io_main
  use mesh_utilities, only: interpolate_to_point_dp_2D, extrapolate_Gaussian
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  use plane_geometry, only: is_in_polygon, triangle_area
  use mpi_distributed_memory, only: gather_to_all
  use ice_geometry_basics, only: is_floating
  use projections, only: oblique_sg_projection
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, &
    ddx_b_a_2D, ddy_b_a_2D
  use create_maps_grid_mesh, only: create_map_from_xy_grid_to_mesh, create_map_from_xy_grid_to_mesh_triangles
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_master
  use bedrock_cumulative_density_functions
  use subgrid_grounded_fractions_main
  use masks_mod
  use zeta_gradients
  use subgrid_ice_margin
  use ice_thickness_safeties
  use inversion_utilities

  implicit none

contains

! == Masks
! ========

! == Effective ice thickness
! ==========================

! == Zeta gradients
! =================

! == No-ice mask
! ==============

! == Ice thickness modification
! =============================

! == Trivia
! =========

! == Target dHi_dt initialisation
! ===============================

! == Target uabs_surf initialisation
! ==================================

end module ice_model_utilities
