module netcdf_write_field_grid
  !< Write data to a field in a grid-based NetCDF file

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary
  use netcdf_basic

  implicit none

  private

  public :: write_to_field_multopt_grid_int_2D, write_to_field_multopt_grid_int_3D, &
    write_to_field_multopt_grid_dp_2D, write_to_field_multopt_grid_dp_2D_monthly, &
    write_to_field_multopt_grid_dp_3D, write_to_field_multopt_grid_dp_3D_ocean, &
    write_to_field_multopt_grid_int_2D_notime, write_to_field_multopt_grid_int_3D_notime, &
    write_to_field_multopt_grid_dp_2D_notime, write_to_field_multopt_grid_dp_2D_monthly_notime, &
    write_to_field_multopt_grid_dp_3D_notime, write_to_field_multopt_grid_dp_3D_ocean_notime


contains

  ! Write to fields with a time dimension

  subroutine write_to_field_multopt_grid_int_2D( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 2-D data field defined on a grid to a NetCDF file variable on the same grid

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_grid),        intent(in   ) :: grid
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    integer,  dimension(:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_grid_int_2D'
    integer                                 :: id_var, id_dim_time, ti
    character(len=1024)                     :: var_name
    integer,  dimension(:,:  ), allocatable :: d_grid
    integer,  dimension(:,:,:), allocatable :: d_grid_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_int_2D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_grid_with_time( grid%nx, grid%ny,1))
      d_grid_with_time( :,:,1) = d_grid
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_int_2D

  subroutine write_to_field_multopt_grid_int_3D( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 3-D data field defined on a grid to a NetCDF file variable on the same grid

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    integer,  dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'write_to_field_multopt_grid_int_3D'
    integer                                   :: id_var, id_dim_time, ti, nz
    character(len=1024)                       :: var_name
    integer,  dimension(:,:,:  ), allocatable :: d_grid
    integer,  dimension(:,:,:,:), allocatable :: d_grid_with_time

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( d_grid_vec_partial,2)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_int_3D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, nz))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_grid_with_time( grid%nx, grid%ny,nz,1))
      d_grid_with_time( :,:,:,1) = d_grid
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nz, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_int_3D

  subroutine write_to_field_multopt_grid_dp_2D( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 2-D data field defined on a grid to a NetCDF file variable on the same grid

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_grid),        intent(in   ) :: grid
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_grid_dp_2D'
    integer                                 :: id_var, id_dim_time, ti
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_grid
    real(dp), dimension(:,:,:), allocatable :: d_grid_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_grid_with_time( grid%nx, grid%ny,1))
      d_grid_with_time( :,:,1) = d_grid
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_2D

  subroutine write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 2-D monthly data field defined on a grid to a NetCDF file variable on the same grid

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'write_to_field_multopt_grid_dp_2D_monthly'
    integer                                   :: id_var, id_dim_time, ti
    character(len=1024)                       :: var_name
    real(dp), dimension(:,:,:  ), allocatable :: d_grid
    real(dp), dimension(:,:,:,:), allocatable :: d_grid_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, 12))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_grid_with_time( grid%nx, grid%ny,12,1))
      d_grid_with_time( :,:,:,1) = d_grid
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, 12, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_2D_monthly

  subroutine write_to_field_multopt_grid_dp_3D( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 3-D data field defined on a grid to a NetCDF file variable on the same grid

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'write_to_field_multopt_grid_dp_3D'
    integer                                   :: id_var, id_dim_time, ti, nz
    character(len=1024)                       :: var_name
    real(dp), dimension(:,:,:  ), allocatable :: d_grid
    real(dp), dimension(:,:,:,:), allocatable :: d_grid_with_time

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( d_grid_vec_partial,2)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, nz))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_grid_with_time( grid%nx, grid%ny,nz,1))
      d_grid_with_time( :,:,:,1) = d_grid
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nz, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_3D

  subroutine write_to_field_multopt_grid_dp_3D_ocean( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 3-D ocean data field defined on a grid to a NetCDF file variable on the same grid

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'write_to_field_multopt_grid_dp_3D_ocean'
    integer                                   :: id_var, id_dim_time, ti, nz
    character(len=1024)                       :: var_name
    real(dp), dimension(:,:,:  ), allocatable :: d_grid
    real(dp), dimension(:,:,:,:), allocatable :: d_grid_with_time

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( d_grid_vec_partial,2)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, nz))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_grid_with_time( grid%nx, grid%ny,nz,1))
      d_grid_with_time( :,:,:,1) = d_grid
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nz, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_3D_ocean

  ! Write to fields without a time dimension

  subroutine write_to_field_multopt_grid_int_2D_notime( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 2-D data field defined on a grid to a NetCDF file variable on the same grid

    ! In/output variables:
    type(type_grid),        intent(in   ) :: grid
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    integer,  dimension(:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_grid_int_2D_notime'
    integer                               :: id_var
    character(len=1024)                   :: var_name
    integer,  dimension(:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_int_2D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_int_2D_notime

  subroutine write_to_field_multopt_grid_int_3D_notime( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 3-D data field defined on a grid to a NetCDF file variable on the same grid

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    integer,  dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_grid_int_3D_notime'
    integer                                 :: id_var, nz
    character(len=1024)                     :: var_name
    integer,  dimension(:,:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( d_grid_vec_partial,2)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_int_3D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, nz))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_int_3D_notime

  subroutine write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 2-D data field defined on a grid to a NetCDF file variable on the same grid

    ! In/output variables:
    type(type_grid),        intent(in   ) :: grid
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_grid_dp_2D_notime'
    integer                               :: id_var
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_2D_notime

  subroutine write_to_field_multopt_grid_dp_2D_monthly_notime( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 2-D monthly data field defined on a grid to a NetCDF file variable on the same grid

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_grid_dp_2D_monthly_notime'
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, 12))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_2D_monthly_notime

  subroutine write_to_field_multopt_grid_dp_3D_notime( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 3-D data field defined on a grid to a NetCDF file variable on the same grid

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_grid_dp_3D_notime'
    integer                                 :: id_var, nz
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( d_grid_vec_partial,2)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, nz))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_3D_notime

  subroutine write_to_field_multopt_grid_dp_3D_ocean_notime( grid, filename, ncid, &
    field_name_options, d_grid_vec_partial)
    !< Write a 3-D data field defined on a grid to a NetCDF file variable on the same grid

    ! In/output variables:
    type(type_grid),          intent(in   ) :: grid
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_grid_dp_3D_ocean_notime'
    integer                                 :: id_var, nz
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    nz = size( d_grid_vec_partial,2)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_xy_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the primary
    if (par%primary) allocate( d_grid( grid%nx, grid%ny, nz))
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_grid_dp_3D_ocean_notime

end module netcdf_write_field_grid
