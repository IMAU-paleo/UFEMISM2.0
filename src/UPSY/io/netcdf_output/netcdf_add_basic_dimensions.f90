module netcdf_add_basic_dimensions

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use netcdf_basic
  use netcdf, only: NF90_INT, NF90_DOUBLE, NF90_UNLIMITED

  implicit none

  private

  public :: add_time_dimension_to_file, add_month_dimension_to_file, add_zeta_dimension_to_file, &
    add_depth_dimension_to_file

contains

  subroutine add_time_dimension_to_file( filename, ncid)
    !< Add a time dimension and variable to an existing NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_time_dimension_to_file'
    integer                        :: id_dim_time
    integer                        :: id_var_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Create time dimension
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_time), NF90_UNLIMITED, id_dim_time)

    ! Create time variable
    call create_variable(  filename, ncid, get_first_option_from_list( field_name_options_time), NF90_DOUBLE, (/ id_dim_time /), id_var_time)
    call add_attribute_char( filename, ncid, id_var_time, 'long_name', 'Time')
    call add_attribute_char( filename, ncid, id_var_time, 'units', 'years')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_time_dimension_to_file

  subroutine add_month_dimension_to_file( filename, ncid)
    !< Add a month dimension and variable to an existing NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_month_dimension_to_file'
    integer                        :: id_dim_month
    integer                        :: id_var_month

    ! Add routine to path
    call init_routine( routine_name)

    ! Create month dimension
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_month), 12, id_dim_month)

    ! Create month variable
    call create_variable(  filename, ncid, get_first_option_from_list( field_name_options_month), NF90_INT, (/ id_dim_month /), id_var_month)
    call add_attribute_char( filename, ncid, id_var_month, 'long_name', 'Month')
    call add_attribute_char( filename, ncid, id_var_month, 'units', '1-12')
    call add_attribute_char( filename, ncid, id_var_month, 'description', '1 = Jan, 2 = Feb, ..., 12 = Dec')

    ! Write month variable
    call write_var_primary( filename, ncid, id_var_month, (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_month_dimension_to_file

  subroutine add_zeta_dimension_to_file( filename, ncid, zeta)
    !< Add a zeta dimension and variable to an existing NetCDF file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    real(dp), dimension(:), intent(in   ) :: zeta

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_zeta_dimension_to_file'
    integer                        :: nz
    integer                        :: id_dim_zeta
    integer                        :: id_var_zeta

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( zeta,1)

    ! Create month dimension
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_zeta), nz, id_dim_zeta)

    ! Create month variable
    call create_variable(  filename, ncid, get_first_option_from_list( field_name_options_zeta), NF90_DOUBLE, (/ id_dim_zeta /), id_var_zeta)
    call add_attribute_char( filename, ncid, id_var_zeta, 'long_name', 'Scaled vertical coordinate')
    call add_attribute_char( filename, ncid, id_var_zeta, 'units', '0-1')
    call add_attribute_char( filename, ncid, id_var_zeta, 'transformation', 'zeta = (h - z) / H; zeta = 0 at the ice surface; zeta = 1 at the ice base')

    ! Write month variable
    call write_var_primary( filename, ncid, id_var_zeta, zeta)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_zeta_dimension_to_file

  subroutine add_depth_dimension_to_file( filename, ncid, depth)
    !< Add a depth dimension and variable to an existing NetCDF file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    real(dp), dimension(:), intent(in   ) :: depth

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_depth_dimension_to_file'
    integer                        :: nz
    integer                        :: id_dim_depth
    integer                        :: id_var_depth

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( depth,1)

    ! Create month dimension
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_depth), nz, id_dim_depth)

    ! Create month variable
    call create_variable(  filename, ncid, get_first_option_from_list( field_name_options_depth), NF90_DOUBLE, (/ id_dim_depth /), id_var_depth)
    call add_attribute_char( filename, ncid, id_var_depth, 'long_name', 'Depth')
    call add_attribute_char( filename, ncid, id_var_depth, 'units', 'meters')

    ! Write month variable
    call write_var_primary( filename, ncid, id_var_depth, depth)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_depth_dimension_to_file

end module netcdf_add_basic_dimensions
