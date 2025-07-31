module netcdf_basic
  !< Wrapper module to reduce the number of "use" statements in other modules

  use netcdf_basic_wrappers
  use netcdf_field_name_options
  use netcdf_inquire_dimensions
  use netcdf_inquire_grid_mesh
  use netcdf_read_var_primary
  use netcdf_check_dimensions
  use netcdf_write_var_primary
  use netcdf_check_fields
  use netcdf_find_timeframe
  use netcdf_save_single_variables

end module netcdf_basic
