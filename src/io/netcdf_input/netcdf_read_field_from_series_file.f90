module netcdf_read_field_from_series_file
  use mpi
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mpi_distributed_memory, only: distribute_from_master
  use netcdf_determine_indexing
  use netcdf_basic
  use grid_lonlat_basic
  use netcdf_setup_grid_mesh_from_file
  use grid_types, only: type_grid_lonlat, type_grid_lat
  use flip_mod

  implicit none

  private

  public :: read_field_from_series_file_monthly, read_field_from_lat_file_1D_monthly, read_time_from_file

contains

  subroutine read_field_from_series_file_monthly( filename, field_name_options, series, time_to_read)
    !< Read a 1-D data field (e.g., time series) from a NetCDF file

    ! In/output variables:
    character(len=*),          intent( in) :: filename
    character(len=*),          intent( in) :: field_name_options
    real(dp), dimension(:),    intent(out) :: series
    real(dp),                  intent( in) :: time_to_read

    ! Local variables
    character(len=1024), parameter         :: routine_name = 'read_field_from_series_file_monthly'
    integer                                :: ncid
    integer                                :: id_var
    character(len=1024)                    :: var_name
    real(dp), dimension(:), allocatable    :: monthly_cycle
    integer                                :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read data from file
    ! ======================
    
    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the file has a valid month dimension
    call check_month( filename, ncid)

    !Allocate memory for time series
    if (par%master) allocate( monthly_cycle( 12))

    ! Find out which timeframe to read
    call find_timeframe( filename, ncid, time_to_read, ti)

    ! Read data
    call read_var_master( filename, ncid, id_var, monthly_cycle, start=(/ 1, ti /), count=(/ 12, 1 /) )

    ! Copy to output memory
    series = monthly_cycle

    ! Clean up after yourself
    if (par%master) deallocate(monthly_cycle)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_series_file_monthly

  subroutine read_field_from_lat_file_1D_monthly( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 1-D monthly data field from a NetCDF file on a lat-grid, copying the values along the 'lon' dimension

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_lat_file_1D_monthly'
    integer                                   :: ncid
    type(type_grid_lat)                       :: vec_loc
    type(type_grid_lonlat)                    :: grid_loc
    integer                                   :: id_var
    character(len=1024)                       :: var_name
    character(len=1024)                       :: latdir
    real(dp), dimension(:,:  ),   allocatable :: d_vec
    real(dp), dimension(:,:,:),   allocatable :: d_vec_with_time
    real(dp), dimension(:,:,:  ), allocatable :: d_grid
    real(dp), dimension(:,:,:,:), allocatable :: d_grid_with_time
    integer                                   :: ti, i

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    call setup_lonlat_grid_from_lat_file( filename, ncid, grid_loc, vec_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the file has a valid month dimension
    call check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    call check_lat_grid_field_dp_1D_monthly( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%master) allocate( d_vec( vec_loc%nlat, 12))
    if (par%master) allocate( d_grid( grid_loc%nlon, grid_loc%nlat, 12))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_master( filename, ncid, id_var, d_vec)
      
      ! copy along the longitudes
      if (par%master) then
        do i = 1, grid_loc%nlon
            d_grid(i,:,:) = d_vec
        end do
      end if
    else
      ! allocate memory
      !if (par%master) allocate( d_grid_with_time( grid_loc%nlon, vec_loc%nlat, 12, 1)) ! probably not used
      if (par%master) allocate( d_vec_with_time( vec_loc%nlat, 12, 1))

      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)

      ! Read data
      call read_var_master( filename, ncid, id_var, d_vec_with_time, start = (/ 1, 1, ti /), count = (/ vec_loc%nlat, 12, 1 /) )
      
      ! Copy to output memory, replicating along the longitudes
      if (par%master) then
        do i = 1, grid_loc%nlon
            d_grid(i,:,:) = d_vec_with_time(:,:,1)
        end do

        ! Copy to output memory = already done in the step above
        !d_grid = d_grid_with_time( :,:,1)

        ! Clean up after yourself
        deallocate( d_vec_with_time)
      end if
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Perform necessary corrections to the gridded data
    ! ====================================================

    ! Determine indexing and dimension directions
    call determine_lat_indexing( filename, ncid, var_name, latdir)

    if     (latdir == 'normal') then
      ! No need to do anything
    elseif (latdir == 'reverse') then
      call flip( grid_loc%lat)
      if (par%master) call flip( d_grid, 2)
    else
      call crash('unknown latdir = "' // trim( latdir) // '"!')
    end if

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_lonlat_gridded_data_from_master_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_grid)
    call deallocate_lonlat_grid( grid_loc)
    call deallocate_lat_grid(vec_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_lat_file_1D_monthly

  subroutine read_time_from_file( filename, time_from_file)
    ! Find the timeframe in the file that is closest to the desired time.
    ! if the file has no time dimension or variable, throw an error.

    ! In/output variables:
    character(len=*),               intent(in   )  :: filename
    real(dp), dimension(:),         intent(out   ) :: time_from_file

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'read_time_from_file'
    integer                             :: ncid
    integer                             :: nt, id_dim_time, id_var_time
    !real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Check if the file contains a valid time dimension and variable
    call check_time( filename, ncid)

    ! inquire size of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! inquire time variable ID
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! allocate memory
    !if (par%master) allocate( time_from_file( nt)) ! TODO: do we need to allocate time?

    ! Read time from file
    call read_var_master( filename, ncid, id_var_time, time_from_file)
    call MPI_BCAST( time_from_file, nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_time_from_file


end module netcdf_read_field_from_series_file
