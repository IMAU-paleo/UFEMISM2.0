module netcdf_input
  !< Wrapper module to reduce the number of "use" statements in other modules

  use netcdf_determine_indexing
  use netcdf_setup_grid_mesh_from_file
  use netcdf_read_field_from_mesh_file
  use netcdf_read_field_from_series_file
  use netcdf_read_field_from_lonlat_grid_file
  use netcdf_read_field_from_xy_grid_file

end module netcdf_input
