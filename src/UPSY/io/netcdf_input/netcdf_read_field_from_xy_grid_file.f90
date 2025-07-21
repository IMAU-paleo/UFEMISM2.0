module netcdf_read_field_from_xy_grid_file
  !< Read data fields from a x/y-grid file

  use mpi_f08, only:
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary
  use grid_basic, only: deallocate_grid
  use permute_mod
  use flip_mod
  use netcdf_basic
  use netcdf_setup_grid_mesh_from_file
  use netcdf_determine_indexing

  implicit none

  private

  public :: read_field_from_xy_file_int_2D, read_field_from_xy_file_int_3D, &
    read_field_from_xy_file_dp_2D, read_field_from_xy_file_dp_2D_monthly, &
    read_field_from_xy_file_dp_3D, read_field_from_xy_file_dp_3D_ocean

contains

subroutine read_field_from_xy_file_int_2D( filename, field_name_options, &
  d_grid_vec_partial, time_to_read)
  !< Read a 2-D data field from a NetCDF file on an x/y-grid

  ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

  ! In/output variables:
  character(len=*),                 intent(in   ) :: filename
  character(len=*),                 intent(in   ) :: field_name_options
  integer,  dimension(:),           intent(  out) :: d_grid_vec_partial
  real(dp),               optional, intent(in   ) :: time_to_read

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_field_from_xy_file_int_2D'
  integer                                 :: ncid
  type(type_grid)                         :: grid_loc
  integer                                 :: id_var
  character(len=1024)                     :: var_name
  character(len=1024)                     :: indexing, xdir, ydir
  integer,  dimension(:,:  ), allocatable :: d_grid
  integer,  dimension(:,:,:), allocatable :: d_grid_with_time
  integer                                 :: ti

  ! Add routine to path
  call init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

  ! Open the NetCDF file
  call open_existing_netcdf_file_for_reading( filename, ncid)

  ! Set up the grid from the file
  call setup_xy_grid_from_file( filename, ncid, grid_loc)

  ! Look for the specified variable in the file
  call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
  if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

  ! Check if the variable has the required dimensions
  call check_xy_grid_field_int_2D( filename, ncid, var_name, should_have_time = present( time_to_read))

  ! Determine indexing and dimension directions
  call determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

  if (indexing == 'xy') then

    ! allocate memory
    if (par%primary) then
      allocate( d_grid( grid_loc%nx, grid_loc%ny))
    else
      allocate( d_grid( 0,0))
    end if

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_grid)
    else
      ! allocate memory
      if (par%primary) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 1 /) )
      ! Copy to output memory
      if (par%primary) d_grid = d_grid_with_time( :,:,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_grid_with_time)
    end if

  elseif (indexing == 'yx') then

    ! allocate memory
    if (par%primary) then
      allocate( d_grid( grid_loc%ny, grid_loc%nx))
    else
      allocate( d_grid( 0,0))
    end if

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_grid)
    else
      ! allocate memory
      if (par%primary) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 1 /) )
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
  if     (indexing == 'xy') then
    ! No need to do anything
  elseif (indexing == 'yx') then
    if (par%primary) call permute( d_grid, map = [2,1])
  else
    call crash('unknown indexing = "' // trim( indexing) // '"!')
  end if

  ! xdir
  if     (xdir == 'normal') then
    ! No need to do anything
  elseif (xdir == 'reverse') then
    call flip( grid_loc%x)
    if (par%primary) call flip( d_grid, 1)
  else
    call crash('unknown xdir = "' // trim( xdir) // '"!')
  end if

  ! ydir
  if     (ydir == 'normal') then
    ! No need to do anything
  elseif (ydir == 'reverse') then
    call flip( grid_loc%y)
    if (par%primary) call flip( d_grid, 2)
  else
    call crash('unknown ydir = "' // trim( ydir) // '"!')
  end if

  ! == Distribute gridded data from the primary to all processes in partial vector form
  ! ==================================================================================

  ! Distribute data
  call distribute_gridded_data_from_primary( grid_loc, d_grid, d_grid_vec_partial)

  ! Clean up after yourself
  call deallocate_grid( grid_loc)
  deallocate( d_grid)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_field_from_xy_file_int_2D

subroutine read_field_from_xy_file_int_3D( filename, field_name_options, &
  d_grid_vec_partial, time_to_read)
  !< Read a 3-D data field from a NetCDF file on an x/y-grid

  ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

  ! In/output variables:
  character(len=*),                   intent(in   ) :: filename
  character(len=*),                   intent(in   ) :: field_name_options
  integer,  dimension(:,:),           intent(  out) :: d_grid_vec_partial
  real(dp),                 optional, intent(in   ) :: time_to_read

  ! Local variables:
  character(len=1024), parameter            :: routine_name = 'read_field_from_xy_file_int_3D'
  integer                                   :: ncid
  type(type_grid)                           :: grid_loc
  integer                                   :: nzeta_loc
  real(dp), dimension(:    ), allocatable   :: zeta_loc
  integer                                   :: id_var
  character(len=1024)                       :: var_name
  character(len=1024)                       :: indexing, xdir, ydir
  integer,  dimension(:,:,:  ), allocatable :: d_grid
  integer,  dimension(:,:,:,:), allocatable :: d_grid_with_time
  integer                                   :: ti

  ! Add routine to path
  call init_routine( routine_name)

  ! == Read grid and data from file
  ! ===============================

  ! Open the NetCDF file
  call open_existing_netcdf_file_for_reading( filename, ncid)

  ! Set up the grid from the file
  call setup_xy_grid_from_file( filename, ncid, grid_loc)

  ! Set up the vertical coordinate zeta from the file
  call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

  ! Look for the specified variable in the file
  call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
  if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

  ! Check if the variable has the required dimensions
  call check_xy_grid_field_int_3D( filename, ncid, var_name, should_have_time = present( time_to_read))

  ! Determine indexing and dimension directions
  call determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

  if (indexing == 'xy') then

    ! allocate memory
    if (par%primary) allocate( d_grid( grid_loc%nx, grid_loc%ny, nzeta_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_grid)
    else
      ! allocate memory
      if (par%primary) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, nzeta_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, nzeta_loc, 1 /) )
      ! Copy to output memory
      if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_grid_with_time)
    end if

  elseif (indexing == 'yx') then

    ! allocate memory
    if (par%primary) allocate( d_grid( grid_loc%ny, grid_loc%nx, nzeta_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_grid)
    else
      ! allocate memory
      if (par%primary) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, nzeta_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, nzeta_loc, 1 /) )
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
  if     (indexing == 'xy') then
    ! No need to do anything
  elseif (indexing == 'yx') then
    if (par%primary) call permute( d_grid, map = [2,1,3])
  else
    call crash('unknown indexing = "' // trim( indexing) // '"!')
  end if

  ! xdir
  if     (xdir == 'normal') then
    ! No need to do anything
  elseif (xdir == 'reverse') then
    call flip( grid_loc%x)
    if (par%primary) call flip( d_grid, 1)
  else
    call crash('unknown xdir = "' // trim( xdir) // '"!')
  end if

  ! ydir
  if     (ydir == 'normal') then
    ! No need to do anything
  elseif (ydir == 'reverse') then
    call flip( grid_loc%y)
    if (par%primary) call flip( d_grid, 2)
  else
    call crash('unknown ydir = "' // trim( ydir) // '"!')
  end if

  ! == Distribute gridded data from the primary to all processes in partial vector form
  ! ==================================================================================

  ! Distribute data
  call distribute_gridded_data_from_primary( grid_loc, d_grid, d_grid_vec_partial)

  ! Clean up after yourself
  if (par%primary) deallocate( d_grid)
  call deallocate_grid( grid_loc)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_field_from_xy_file_int_3D

  subroutine read_field_from_xy_file_dp_2D( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on an x/y-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    real(dp), dimension(:),           intent(  out) :: d_grid_vec_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_xy_file_dp_2D'
    integer                                 :: ncid
    type(type_grid)                         :: grid_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    character(len=1024)                     :: indexing, xdir, ydir
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
    call setup_xy_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    if (indexing == 'xy') then

      ! allocate memory
      if (par%primary) then
        allocate( d_grid( grid_loc%nx, grid_loc%ny))
      else
        allocate( d_grid( 0,0))
      end if

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%primary) then
        allocate( d_grid( grid_loc%ny, grid_loc%nx))
      else
        allocate( d_grid( 0,0))
      end if

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 1 /) )
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
    if     (indexing == 'xy') then
      ! No need to do anything
    elseif (indexing == 'yx') then
      if (par%primary) call permute( d_grid, map = [2,1])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_primary( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    call deallocate_grid( grid_loc)
    deallocate( d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_dp_2D

  subroutine read_field_from_xy_file_dp_2D_monthly( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 2-D monthly data field from a NetCDF file on an x/y-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_xy_file_dp_2D_monthly'
    integer                                   :: ncid
    type(type_grid)                           :: grid_loc
    integer                                   :: id_var
    character(len=1024)                       :: var_name
    character(len=1024)                       :: indexing, xdir, ydir
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
    call setup_xy_grid_from_file( filename, ncid, grid_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the file has a valid month dimension
    call check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    call check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    if (indexing == 'xy') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nx, grid_loc%ny, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 12, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%ny, grid_loc%nx, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 12, 1 /) )
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
    if     (indexing == 'xy') then
      ! No need to do anything
    elseif (indexing == 'yx') then
      if (par%primary) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_primary( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid)
    call deallocate_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_dp_2D_monthly

  subroutine read_field_from_xy_file_dp_3D( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 3-D data field from a NetCDF file on an x/y-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_xy_file_dp_3D'
    integer                                   :: ncid
    type(type_grid)                           :: grid_loc
    integer                                   :: nzeta_loc
    real(dp), dimension(:    ), allocatable   :: zeta_loc
    integer                                   :: id_var
    character(len=1024)                       :: var_name
    character(len=1024)                       :: indexing, xdir, ydir
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
    call setup_xy_grid_from_file( filename, ncid, grid_loc)

    ! Set up the vertical coordinate zeta from the file
    call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    if (indexing == 'xy') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nx, grid_loc%ny, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, nzeta_loc, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%ny, grid_loc%nx, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, nzeta_loc, 1 /) )
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
    if     (indexing == 'xy') then
      ! No need to do anything
    elseif (indexing == 'yx') then
      if (par%primary) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_primary( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid)
    call deallocate_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_dp_3D

  subroutine read_field_from_xy_file_dp_3D_ocean( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    ! Read a 3-D data ocean field from a NetCDF file on an x/y-grid
    !
    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_xy_file_dp_3D_ocean'
    integer                                   :: ncid
    type(type_grid)                           :: grid_loc
    integer                                   :: ndepth_loc
    real(dp), dimension(:    ), allocatable   :: depth_loc
    integer                                   :: id_var
    character(len=1024)                       :: var_name
    character(len=1024)                       :: indexing, xdir, ydir
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
    call setup_xy_grid_from_file( filename, ncid, grid_loc)

    ! Set up the vertical coordinate depth from the file
    call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_xy_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! Determine indexing and dimension directions
    call determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)

    if (indexing == 'xy') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%nx, grid_loc%ny, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, ndepth_loc, 1 /) )
        ! Copy to output memory
        if (par%primary) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%primary) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%primary) allocate( d_grid( grid_loc%ny, grid_loc%nx, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_primary( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%primary) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, ndepth_loc, 1 /) )
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
    if     (indexing == 'xy') then
      ! No need to do anything
    elseif (indexing == 'yx') then
      if (par%primary) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%primary) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%primary) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_primary( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid)
    call deallocate_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_dp_3D_ocean

end module netcdf_read_field_from_xy_grid_file
