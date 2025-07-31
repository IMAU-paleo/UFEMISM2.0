module netcdf_output
  !< Wrapper module to reduce the number of "use" statements in other modules

  use netcdf_generate_numbered_filename
  use netcdf_add_write_scalar_variables
  use netcdf_add_basic_dimensions
  use netcdf_add_field_grid
  use netcdf_add_field_mesh
  use netcdf_setup_grid_mesh_in_file
  use netcdf_write_field_grid
  use netcdf_write_field_mesh
  use netcdf_resource_tracking

end module netcdf_output
