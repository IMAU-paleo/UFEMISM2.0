module netcdf_read_field_from_orca_grid_file
  !< Read data fields from an ORCA grid file

  use mpi
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid_orca
  use permute_mod
  use flip_mod
  use grid_orca_basic, only: distribute_orca_gridded_data_from_master_dp_3D, deallocate_orca_grid
  use netcdf_basic
  use netcdf_setup_grid_mesh_from_file
  use netcdf_determine_indexing

  implicit none

  private

  public :: read_field_from_orca_file_3D_ocean

contains

  subroutine read_field_from_orca_file_3D_ocean( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 3-D ocean data field from a NetCDF file on an ORCA grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_orca_file_3D_ocean'
    integer                                   :: ncid
    type(type_grid_orca)                      :: grid_loc
    integer                                   :: ndepth_loc
    real(dp), dimension(:), allocatable       :: depth_loc
    integer                                   :: id_var
    character(len=1024)                       :: var_name
    character(len=1024)                       :: indexing, londir, latdir
    real(dp), dimension(:,:,:  ), allocatable :: d_grid
    real(dp), dimension(:,:,:,:), allocatable :: d_grid_with_time
    integer                                   :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! TODO EL

    ! Set up the grid from the file
    ! call setup_orca_grid_from_file( filename, ncid, grid_loc)

    ! Set up the vertical coordinate depth from the file
    ! call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)

    ! Look for the specified variable in the file
    ! call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    ! if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    ! call check_orca_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    ! if (par%master) allocate( d_grid( grid_loc%ny, grid_loc%nx, ndepth_loc))

    ! Read data from file
    ! if (.not. present( time_to_read)) then
    !   call read_var_master( filename, ncid, id_var, d_grid)
    ! else
      ! allocate memory
      ! if (par%master) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, ndepth_loc, 1))
      ! Find out which timeframe to read
      ! call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      ! call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, ndepth_loc, 1 /) )
      ! Copy to output memory
      ! if (par%master) d_grid = d_grid_with_time( :,:,:,1)
      ! Clean up after yourself
      ! if (par%master) deallocate( d_grid_with_time)
    ! end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    ! call distribute_orca_gridded_data_from_master_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    ! if (par%master) deallocate( d_grid)
    ! call deallocate_orca_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_orca_file_3D_ocean

end module netcdf_read_field_from_orca_grid_file
