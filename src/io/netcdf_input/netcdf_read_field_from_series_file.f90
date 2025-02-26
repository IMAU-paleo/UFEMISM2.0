module netcdf_read_field_from_series_file
  use mpi
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mpi_distributed_memory, only: distribute_from_master
  use netcdf_basic
  use netcdf_setup_grid_mesh_from_file

  implicit none

  private

  public :: read_field_from_series_file_monthly

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
    allocate( monthly_cycle( 12, 1))

    ! Find out which timeframe to read
    call find_timeframe( filename, ncid, time_to_read, ti)

    ! Read data
    call read_var_master( filename, ncid, id_var, monthly_cycle, start=(/ 1, ti), count=(/ 12, 1 /) )

    ! Copy to output memory
    series = monthly_cycle

    ! Clean up after yourself
    if (par%master) deallocate (monthly_cycle)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_series_file_monthly
