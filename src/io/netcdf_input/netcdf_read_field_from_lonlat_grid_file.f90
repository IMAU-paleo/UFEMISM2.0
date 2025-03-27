module netcdf_read_field_from_lonlat_grid_file
  !< Read data fields from a lon/lat-grid file

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid_lonlat
  use permute_mod
  use flip_mod
  use grid_lonlat_basic, only: distribute_lonlat_gridded_data_from_primary_dp_2D, &
    distribute_lonlat_gridded_data_from_primary_dp_3D, deallocate_lonlat_grid
  use netcdf_basic
  use netcdf_setup_grid_mesh_from_file
  use netcdf_determine_indexing

  implicit none

  private

  public :: read_field_from_lonlat_file_2D, read_field_from_lonlat_file_2D_monthly, &
    read_field_from_lonlat_file_3D, read_field_from_lonlat_file_3D_ocean

contains

  subroutine read_field_from_lonlat_file_2D( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on a lon/lat-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    real(dp), dimension(:),           intent(  out) :: d_grid_vec_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_lonlat_file_2D'
    integer                                 :: ncid
    type(type_grid_lonlat)                  :: grid_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    character(len=1024)                     :: indexing, londir, latdir
    real(dp), dimension(:,:  ), allocatable :: d_grid
    real(dp), dimension(:,:,:), allocatable :: d_grid_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the grid from the file
    call setup_lonlat_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    if (indexing == 'lonlat') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlon, grid_loc%nlat))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlat, grid_loc%nlon))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Perform necessary corrections to the gridded data
    ! ====================================================

    ! Indexing
    if     (indexing == 'lonlat') then
      ! No need to do anything
    elseif (indexing == 'latlon') then
      if (par%primary) call permute( d_grid, map = [2,1])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
    if     (latdir == 'normal') then
      ! No need to do anything
    elseif (latdir == 'reverse') then
      call flip( grid_loc%lat)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown latdir = "' // trim( latdir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_lonlat_gridded_data_from_primary_dp_2D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the primary

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid)
    call deallocate_lonlat_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_lonlat_file_2D

  subroutine read_field_from_lonlat_file_2D_monthly( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 2-D monthly data field from a NetCDF file on a lon/lat-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_lonlat_file_2D_monthly'
    integer                                   :: ncid
    type(type_grid_lonlat)                    :: grid_loc
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

    ! Set up the grid from the file
    call setup_lonlat_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the file has a valid month dimension
    call check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    call check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    if (indexing == 'lonlat') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlon, grid_loc%nlat, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, 12, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlat, grid_loc%nlon, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, 12, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Perform necessary corrections to the gridded data
    ! ====================================================

    ! Indexing
    if     (indexing == 'lonlat') then
      ! No need to do anything
    elseif (indexing == 'latlon') then
      if (par%primary) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
    if     (latdir == 'normal') then
      ! No need to do anything
    elseif (latdir == 'reverse') then
      call flip( grid_loc%lat)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown latdir = "' // trim( latdir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_lonlat_gridded_data_from_primary_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid)
    call deallocate_lonlat_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_lonlat_file_2D_monthly

  subroutine read_field_from_lonlat_file_3D( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 3-D data field from a NetCDF file on a lon/lat-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_lonlat_file_3D'
    integer                                   :: ncid
    type(type_grid_lonlat)                    :: grid_loc
    integer                                   :: nzeta_loc
    real(dp), dimension(:), allocatable       :: zeta_loc
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

    ! Set up the grid from the file
    call setup_lonlat_grid_from_file( filename, ncid, grid_loc)

    ! Set up the vertical coordinate zeta from the file
    call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    if (indexing == 'lonlat') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlon, grid_loc%nlat, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, nzeta_loc, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlat, grid_loc%nlon, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, nzeta_loc, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Perform necessary corrections to the gridded data
    ! ====================================================

    ! Indexing
    if     (indexing == 'lonlat') then
      ! No need to do anything
    elseif (indexing == 'latlon') then
      if (par%primary) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
    if     (latdir == 'normal') then
      ! No need to do anything
    elseif (latdir == 'reverse') then
      call flip( grid_loc%lat)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown latdir = "' // trim( latdir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_lonlat_gridded_data_from_primary_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid)
    call deallocate_lonlat_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_lonlat_file_3D

  subroutine read_field_from_lonlat_file_3D_ocean( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 3-D ocean data field from a NetCDF file on a lon/lat-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_lonlat_file_3D_ocean'
    integer                                   :: ncid
    type(type_grid_lonlat)                    :: grid_loc
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

    ! Set up the grid from the file
    call setup_lonlat_grid_from_file( filename, ncid, grid_loc)

    ! Set up the vertical coordinate depth from the file
    call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_lonlat_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)

    if (indexing == 'lonlat') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlon, grid_loc%nlat, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, ndepth_loc, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nlat, grid_loc%nlon, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, ndepth_loc, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Perform necessary corrections to the gridded data
    ! ====================================================

    ! Indexing
    if     (indexing == 'lonlat') then
      ! No need to do anything
    elseif (indexing == 'latlon') then
      if (par%primary) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
    if     (latdir == 'normal') then
      ! No need to do anything
    elseif (latdir == 'reverse') then
      call flip( grid_loc%lat)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown latdir = "' // trim( latdir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_lonlat_gridded_data_from_primary_dp_3D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid)
    call deallocate_lonlat_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_lonlat_file_3D_ocean

end module netcdf_read_field_from_lonlat_grid_file
