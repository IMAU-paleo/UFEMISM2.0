module netcdf_input

! ===== Flexible reading of input files =====
! ===========================================
!
! These routines allow for flexible reading of input files. The three top-level functions
! allow you to read 2-D, 2-D monthly, and 3-D data fields from NetCDF files. The files can
! contain the data on a regular x/y-grid, a regular lon/lat-grid, or a mesh. The routines
! will automatically detect which one of these it is, and select the appropriate
! subroutine. The data will also be automatically mapped to the provided model mesh.

  use mpi
  use precisions, only: dp
  use mpi_basic, only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use mpi_distributed_memory, only: distribute_from_master
  use grid_basic, only: type_grid, calc_secondary_grid_data, deallocate_grid
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_master
  use grid_lonlat_basic, only: type_grid_lonlat, calc_lonlat_field_to_vector_form_translation_tables, &
    distribute_lonlat_gridded_data_from_master_dp_2D, deallocate_lonlat_grid, &
    distribute_lonlat_gridded_data_from_master_dp_3D
  use permute_mod, only: permute
  use flip_mod, only: flip
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary, deallocate_mesh
  use mesh_utilities, only: check_mesh
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_parallel_creation, only: broadcast_mesh
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use remapping_main, only: map_from_xy_grid_to_mesh_2D, map_from_lonlat_grid_to_mesh_2D, map_from_mesh_to_mesh_2D, &
    map_from_xy_grid_to_mesh_3D, map_from_lonlat_grid_to_mesh_3D, map_from_mesh_to_mesh_3D, &
    map_from_vertical_to_vertical_2D_ocean

  use netcdf, only: NF90_MAX_VAR_DIMS
  use netcdf_field_name_options
  use netcdf_inquire_grid_mesh
  use netcdf_read_var_master
  use netcdf_check_dimensions
  use netcdf_basic_wrappers
  use netcdf_check_fields
  use netcdf_find_timeframe

  implicit none

contains

  ! ===== Top-level functions =====
  ! ===============================

  ! Read and map to mesh
  subroutine read_field_from_file_2D( filename, field_name_options, &
    mesh, d_partial, time_to_read)
    !< Read a data field from a NetCDF file, and map it to the model mesh.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    character(len=*),       intent(in   ) :: field_name_options
    type(type_mesh),        intent(in   ) :: mesh
    real(dp), dimension(:), intent(  out) :: d_partial
    real(dp), optional,     intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'read_field_from_file_2D'
    logical                             :: file_exists
    logical                             :: has_xy_grid, has_lonlat_grid, has_mesh
    integer                             :: ncid
    type(type_grid)                     :: grid_from_file
    type(type_grid_lonlat)              :: grid_lonlat_from_file
    type(type_mesh)                     :: mesh_from_file
    real(dp), dimension(:), allocatable :: d_grid_vec_partial_from_file
    real(dp), dimension(:), allocatable :: d_grid_lonlat_vec_partial_from_file
    real(dp), dimension(:), allocatable :: d_mesh_partial_from_file
    character(len=1024), parameter      :: method_mesh2mesh = '2nd_order_conservative'

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) then
      call crash('file "' // trim( filename) // '" not found!')
    end if

    ! Find out on what kind of grid the file is defined
    call inquire_xy_grid(     filename, has_xy_grid    )
    call inquire_lonlat_grid( filename, has_lonlat_grid)
    call inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    if (has_xy_grid     .and. has_lonlat_grid) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    if (has_xy_grid     .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a mesh!')
    if (has_lonlat_grid .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! Choose the appropriate subroutine
    if (has_xy_grid) then
      ! Data is provided on an x/y-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_xy_grid_from_file( filename, ncid, grid_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_vec_partial_from_file( grid_from_file%n1: grid_from_file%n2))

      ! Read gridded data
      call read_field_from_xy_file_2D( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_xy_grid_to_mesh_2D( grid_from_file, mesh, d_grid_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_grid( grid_from_file)
      deallocate( d_grid_vec_partial_from_file)

    elseif (has_lonlat_grid) then
      ! Data is provided on a lon/lat-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_lonlat_grid_from_file( filename, ncid, grid_lonlat_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_lonlat_vec_partial_from_file( grid_lonlat_from_file%n1: grid_lonlat_from_file%n2))

      ! Read gridded data
      call read_field_from_lonlat_file_2D( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_lonlat_grid_to_mesh_2D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_lonlat_grid( grid_lonlat_from_file)
      deallocate( d_grid_lonlat_vec_partial_from_file)

    elseif (has_mesh) then
      ! Data is provided on a mesh

      ! Set up the mesh from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_mesh_from_file( filename, ncid, mesh_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_mesh_partial_from_file( mesh_from_file%vi1: mesh_from_file%vi2))

      ! Read meshed data
      call read_field_from_mesh_file_2D( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_mesh_to_mesh_2D( mesh_from_file, mesh, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      call deallocate_mesh( mesh_from_file)
      deallocate( d_mesh_partial_from_file)

    else
      call crash('file "' // trim( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_2D

  subroutine read_field_from_file_2D_monthly( filename, field_name_options, &
    mesh, d_partial, time_to_read)
    !< Read a data field from a NetCDF file, and map it to the model mesh.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    ! In/output variables:
    character(len=*),         intent(in   ) :: filename
    character(len=*),         intent(in   ) :: field_name_options
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(  out) :: d_partial
    real(dp), optional,       intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_file_2D_monthly'
    logical                               :: file_exists
    logical                               :: has_xy_grid, has_lonlat_grid, has_mesh
    integer                               :: ncid
    type(type_grid)                       :: grid_from_file
    type(type_grid_lonlat)                :: grid_lonlat_from_file
    type(type_mesh)                       :: mesh_from_file
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_from_file
    real(dp), dimension(:,:), allocatable :: d_grid_lonlat_vec_partial_from_file
    real(dp), dimension(:,:), allocatable :: d_mesh_partial_from_file
    character(len=1024), parameter        :: method_mesh2mesh = '2nd_order_conservative'

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) then
      call crash('file "' // trim( filename) // '" not found!')
    end if

    ! Find out on what kind of grid the file is defined
    call inquire_xy_grid(     filename, has_xy_grid    )
    call inquire_lonlat_grid( filename, has_lonlat_grid)
    call inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    if (has_xy_grid     .and. has_lonlat_grid) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    if (has_xy_grid     .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a mesh!')
    if (has_lonlat_grid .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! Choose the appropriate subroutine
    if (has_xy_grid) then
      ! Data is provided on an x/y-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_xy_grid_from_file( filename, ncid, grid_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_vec_partial_from_file( grid_from_file%n1: grid_from_file%n2,12))

      ! Read gridded data
      call read_field_from_xy_file_2D_monthly( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, d_grid_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_grid( grid_from_file)
      deallocate( d_grid_vec_partial_from_file)

    elseif (has_lonlat_grid) then
      ! Data is provided on a lon/lat-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_lonlat_grid_from_file( filename, ncid, grid_lonlat_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_lonlat_vec_partial_from_file( grid_lonlat_from_file%n1: grid_lonlat_from_file%n2,12))

      ! Read gridded data
      call read_field_from_lonlat_file_2D_monthly( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_lonlat_grid_to_mesh_3D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_lonlat_grid( grid_lonlat_from_file)
      deallocate( d_grid_lonlat_vec_partial_from_file)

    elseif (has_mesh) then
      ! Data is provided on a mesh

      ! Set up the mesh from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_mesh_from_file( filename, ncid, mesh_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_mesh_partial_from_file( mesh_from_file%vi1: mesh_from_file%vi2,12))

      ! Read meshed data
      call read_field_from_mesh_file_2D_monthly( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_mesh_to_mesh_3D( mesh_from_file, mesh, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      call deallocate_mesh( mesh_from_file)
      deallocate( d_mesh_partial_from_file)

    else
      call crash('file "' // trim( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_2D_monthly

  subroutine read_field_from_file_3D( filename, field_name_options, &
    mesh, d_partial, time_to_read, nzeta, zeta)
    !< Read a data field from a NetCDF file, and map it to the model mesh.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.
    !
    ! NOTE: assumes the vertical grid of the input file is identical to that of the model!

    ! In/output variables:
    character(len=*),                              intent(in   ) :: filename
    character(len=*),                              intent(in   ) :: field_name_options
    type(type_mesh),                               intent(in   ) :: mesh
    real(dp), dimension(:,:),                      intent(  out) :: d_partial
    real(dp), optional,                            intent(in   ) :: time_to_read
    integer,                             optional, intent(  out) :: nzeta
    real(dp), dimension(:), allocatable, optional, intent(  out) :: zeta

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_file_3D'
    logical                               :: file_exists
    logical                               :: has_xy_grid, has_lonlat_grid, has_mesh
    integer                               :: ncid
    type(type_grid)                       :: grid_from_file
    type(type_grid_lonlat)                :: grid_lonlat_from_file
    type(type_mesh)                       :: mesh_from_file
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_from_file
    real(dp), dimension(:,:), allocatable :: d_grid_lonlat_vec_partial_from_file
    real(dp), dimension(:,:), allocatable :: d_mesh_partial_from_file
    character(len=1024), parameter        :: method_mesh2mesh = '2nd_order_conservative'
    integer                               :: nzeta_loc
    real(dp), dimension(:), allocatable   :: zeta_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) then
      call crash('file "' // trim( filename) // '" not found!')
    end if

    ! Find out on what kind of grid the file is defined
    call inquire_xy_grid(     filename, has_xy_grid    )
    call inquire_lonlat_grid( filename, has_lonlat_grid)
    call inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    if (has_xy_grid     .and. has_lonlat_grid) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    if (has_xy_grid     .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a mesh!')
    if (has_lonlat_grid .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! Choose the appropriate subroutine
    if (has_xy_grid) then
      ! Data is provided on an x/y-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_xy_grid_from_file( filename, ncid, grid_from_file)
      call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_vec_partial_from_file( grid_from_file%n1: grid_from_file%n2, nzeta_loc))

      ! Read gridded data
      call read_field_from_xy_file_3D( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, d_grid_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_grid( grid_from_file)
      deallocate( d_grid_vec_partial_from_file)

    elseif (has_lonlat_grid) then
      ! Data is provided on a lon/lat-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_lonlat_grid_from_file( filename, ncid, grid_lonlat_from_file)
      call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_lonlat_vec_partial_from_file( grid_lonlat_from_file%n1: grid_lonlat_from_file%n2, nzeta_loc))

      ! Read gridded data
      call read_field_from_lonlat_file_3D( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_lonlat_grid_to_mesh_3D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_lonlat_grid( grid_lonlat_from_file)
      deallocate( d_grid_lonlat_vec_partial_from_file)

    elseif (has_mesh) then
      ! Data is provided on a mesh

      ! Set up the mesh from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_mesh_from_file( filename, ncid, mesh_from_file)
      call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)
      call close_netcdf_file( ncid)

      ! allocate memory for meshed data
      allocate( d_mesh_partial_from_file( mesh_from_file%vi1: mesh_from_file%vi2, nzeta_loc))

      ! Read meshed data
      call read_field_from_mesh_file_3D( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_mesh_to_mesh_3D( mesh_from_file, mesh, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      call deallocate_mesh( mesh_from_file)
      deallocate( d_mesh_partial_from_file)

    else
      call crash('file "' // trim( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    end if

    ! if so specified, return the zeta read from file as output
    if (present( nzeta) .or. present( zeta)) then
      ! Safety
      if (.not. present( nzeta) .or. .not. present( zeta)) &
        call crash('should ask for both nzeta and zeta!')
      nzeta = nzeta_loc
      allocate( zeta( nzeta))
      zeta = zeta_loc
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_3D

  subroutine read_field_from_file_3D_ocean( filename, field_name_options, &
    mesh, d_partial, time_to_read, ndepth, depth)
    !< Read a data field from a NetCDF file, and map it to the model mesh.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    ! In/output variables:
    character(len=*),                              intent(in   ) :: filename
    character(len=*),                              intent(in   ) :: field_name_options
    type(type_mesh),                               intent(in   ) :: mesh
    real(dp), dimension(:,:),                      intent(  out) :: d_partial
    real(dp), optional,                            intent(in   ) :: time_to_read
    integer ,                            optional, intent(  out) :: ndepth
    real(dp), dimension(:), allocatable, optional, intent(  out) :: depth

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_file_3D_ocean'
    logical                               :: file_exists
    logical                               :: has_xy_grid, has_lonlat_grid, has_mesh
    integer                               :: ncid
    type(type_grid)                       :: grid_from_file
    type(type_grid_lonlat)                :: grid_lonlat_from_file
    type(type_mesh)                       :: mesh_from_file
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_from_file
    real(dp), dimension(:,:), allocatable :: d_grid_lonlat_vec_partial_from_file
    real(dp), dimension(:,:), allocatable :: d_mesh_partial_from_file
    real(dp), dimension(:,:), allocatable :: d_partial_raw_layers
    character(len=1024), parameter        :: method_mesh2mesh = '2nd_order_conservative'
    integer                               :: ndepth_loc
    real(dp), dimension(:), allocatable   :: depth_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) then
      call crash('file "' // trim( filename) // '" not found!')
    end if

    ! Find out on what kind of grid the file is defined
    call inquire_xy_grid(     filename, has_xy_grid    )
    call inquire_lonlat_grid( filename, has_lonlat_grid)
    call inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    if (has_xy_grid     .and. has_lonlat_grid) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    if (has_xy_grid     .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both an x/y-grid and a mesh!')
    if (has_lonlat_grid .and. has_mesh       ) call crash('file "' // trim( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! Choose the appropriate subroutine
    if (has_xy_grid) then
      ! Data is provided on an x/y-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_xy_grid_from_file( filename, ncid, grid_from_file)
      call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_vec_partial_from_file( grid_from_file%n1: grid_from_file%n2, ndepth_loc))

      ! Read gridded data
      call read_field_from_xy_file_3D_ocean( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! allocate memory for meshed data using all data layers
      allocate( d_partial_raw_layers( mesh%vi1:mesh%vi2, ndepth_loc))

      ! Remap data horizontally
      call map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, d_grid_vec_partial_from_file, d_partial_raw_layers)

      ! Remap data vertically
      call map_from_vertical_to_vertical_2D_ocean( mesh, depth_loc, C%z_ocean, d_partial_raw_layers, d_partial)

      ! Clean up after yourself
      call deallocate_grid( grid_from_file)
      deallocate( d_grid_vec_partial_from_file)
      deallocate( d_partial_raw_layers)

    elseif (has_lonlat_grid) then
      ! Data is provided on a lon/lat-grid

      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_lonlat_grid_from_file( filename, ncid, grid_lonlat_from_file)
      call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_lonlat_vec_partial_from_file( grid_lonlat_from_file%n1: grid_lonlat_from_file%n2, ndepth_loc))

      ! Read gridded data
      call read_field_from_lonlat_file_3D_ocean( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, time_to_read = time_to_read)

      ! allocate memory for meshed data using all data layers
      allocate( d_partial_raw_layers( mesh%vi1:mesh%vi2, ndepth_loc))

      ! Remap data horizontally
      call map_from_lonlat_grid_to_mesh_3D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial_raw_layers)

      ! Remap data vertically
      call map_from_vertical_to_vertical_2D_ocean( mesh, depth_loc, C%z_ocean, d_partial_raw_layers, d_partial)

      ! Clean up after yourself
      call deallocate_lonlat_grid( grid_lonlat_from_file)
      deallocate( d_grid_lonlat_vec_partial_from_file)
      deallocate( d_partial_raw_layers)

    elseif (has_mesh) then
      ! Data is provided on a mesh

      ! Set up the mesh from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_mesh_from_file( filename, ncid, mesh_from_file)
      call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)
      call close_netcdf_file( ncid)

      ! allocate memory for meshed data
      allocate( d_mesh_partial_from_file( mesh_from_file%vi1: mesh_from_file%vi2, ndepth_loc))

      ! Read meshed data
      call read_field_from_mesh_file_3D_ocean( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! allocate memory for meshed data using all data layers
      allocate( d_partial_raw_layers( mesh%vi1:mesh%vi2, ndepth_loc))

      ! Remap data horizontally
      call map_from_mesh_to_mesh_3D( mesh_from_file, mesh, d_mesh_partial_from_file, d_partial_raw_layers, method = method_mesh2mesh)

      ! Remap data vertically
      call map_from_vertical_to_vertical_2D_ocean( mesh, depth_loc, C%z_ocean, d_partial_raw_layers, d_partial)

      ! Clean up after yourself
      call deallocate_mesh( mesh_from_file)
      deallocate( d_mesh_partial_from_file)
      deallocate( d_partial_raw_layers)

    else
      call crash('file "' // trim( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    end if

    ! if so specified, return the depth read from file as output
    if (present( ndepth) .or. present( depth)) then
      ! Safety
      if (.not. present( ndepth) .or. .not. present( depth)) call crash('should ask for both ndepth and depth!')
      ndepth = ndepth_loc
      allocate( depth( ndepth))
      depth = depth_loc
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_3D_ocean

  ! ===== Medium-level functions =====
  ! ==================================

  ! Read data fields from an x/y-grid file
  subroutine read_field_from_xy_file_2D( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on an x/y-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    real(dp), dimension(:),           intent(  out) :: d_grid_vec_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_xy_file_2D'
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
      if (par%master) then
        allocate( d_grid( grid_loc%nx, grid_loc%ny))
      else
        allocate( d_grid( 0,0))
      end if

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%master) then
        allocate( d_grid( grid_loc%ny, grid_loc%nx))
      else
        allocate( d_grid( 0,0))
      end if

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%master) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_master( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    call deallocate_grid( grid_loc)
    deallocate( d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_2D

  subroutine read_field_from_xy_file_2D_int( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on an x/y-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    integer,  dimension(:),           intent(  out) :: d_grid_vec_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_xy_file_2D_int'
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
      if (par%master) then
        allocate( d_grid( grid_loc%nx, grid_loc%ny))
      else
        allocate( d_grid( 0,0))
      end if

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%master) then
        allocate( d_grid( grid_loc%ny, grid_loc%nx))
      else
        allocate( d_grid( 0,0))
      end if

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%master) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_master( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    call deallocate_grid( grid_loc)
    deallocate( d_grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_2D_int

  subroutine read_field_from_xy_file_2D_monthly( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 2-D monthly data field from a NetCDF file on an x/y-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_xy_file_2D_monthly'
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
      if (par%master) allocate( d_grid( grid_loc%nx, grid_loc%ny, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, 12, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%master) allocate( d_grid( grid_loc%ny, grid_loc%nx, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, 12, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%master) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_master( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_grid)
    call deallocate_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_2D_monthly

  subroutine read_field_from_xy_file_3D( filename, field_name_options, &
    d_grid_vec_partial, time_to_read)
    !< Read a 3-D data field from a NetCDF file on an x/y-grid

    ! NOTE: the grid should be read before, and memory allocated for d_grid_vec_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_grid_vec_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter            :: routine_name = 'read_field_from_xy_file_3D'
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
      if (par%master) allocate( d_grid( grid_loc%nx, grid_loc%ny, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, nzeta_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%master) allocate( d_grid( grid_loc%ny, grid_loc%nx, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, nzeta_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%master) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_master( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_grid)
    call deallocate_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_3D

  subroutine read_field_from_xy_file_3D_ocean( filename, field_name_options, &
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
    character(len=1024), parameter            :: routine_name = 'read_field_from_xy_file_3D_ocean'
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
      if (par%master) allocate( d_grid( grid_loc%nx, grid_loc%ny, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nx, grid_loc%ny, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nx, grid_loc%ny, ndepth_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'yx') then

      ! allocate memory
      if (par%master) allocate( d_grid( grid_loc%ny, grid_loc%nx, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%ny, grid_loc%nx, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%ny, grid_loc%nx, ndepth_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! xdir
    if     (xdir == 'normal') then
      ! No need to do anything
    elseif (xdir == 'reverse') then
      call flip( grid_loc%x)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown xdir = "' // trim( xdir) // '"!')
    end if

    ! ydir
    if     (ydir == 'normal') then
      ! No need to do anything
    elseif (ydir == 'reverse') then
      call flip( grid_loc%y)
      if (par%master) call flip( d_grid, 2)
    else
      call crash('unknown ydir = "' // trim( ydir) // '"!')
    end if

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_gridded_data_from_master( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_grid)
    call deallocate_grid( grid_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_xy_file_3D_ocean

  ! Read data fields from a lon/lat-grid file
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
      if (par%master) allocate( d_grid( grid_loc%nlon, grid_loc%nlat))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%master) allocate( d_grid( grid_loc%nlat, grid_loc%nlon))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
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
    call distribute_lonlat_gridded_data_from_master_dp_2D( grid_loc, d_grid, d_grid_vec_partial)

    ! Clean up gridded data on the master

    ! Clean up after yourself
    if (par%master) deallocate( d_grid)
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
      if (par%master) allocate( d_grid( grid_loc%nlon, grid_loc%nlat, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, 12, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%master) allocate( d_grid( grid_loc%nlat, grid_loc%nlon, 12))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, 12, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, 12, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
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
      if (par%master) allocate( d_grid( grid_loc%nlon, grid_loc%nlat, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, nzeta_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%master) allocate( d_grid( grid_loc%nlat, grid_loc%nlon, nzeta_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, nzeta_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, nzeta_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
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
      if (par%master) allocate( d_grid( grid_loc%nlon, grid_loc%nlat, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlon, grid_loc%nlat, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlon, grid_loc%nlat, ndepth_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
      end if

    elseif (indexing == 'latlon') then

      ! allocate memory
      if (par%master) allocate( d_grid( grid_loc%nlat, grid_loc%nlon, ndepth_loc))

      ! Read data from file
      if (.not. present( time_to_read)) then
        call read_var_master( filename, ncid, id_var, d_grid)
      else
        ! allocate memory
        if (par%master) allocate( d_grid_with_time( grid_loc%nlat, grid_loc%nlon, ndepth_loc, 1))
        ! Find out which timeframe to read
        call find_timeframe( filename, ncid, time_to_read, ti)
        ! Read data
        call read_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid_loc%nlat, grid_loc%nlon, ndepth_loc, 1 /) )
        ! Copy to output memory
        if (par%master) d_grid = d_grid_with_time( :,:,:,1)
        ! Clean up after yourself
        if (par%master) deallocate( d_grid_with_time)
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
      if (par%master) call permute( d_grid, map = [2,1,3])
    else
      call crash('unknown indexing = "' // trim( indexing) // '"!')
    end if

    ! londir
    if     (londir == 'normal') then
      ! No need to do anything
    elseif (londir == 'reverse') then
      call flip( grid_loc%lon)
      if (par%master) call flip( d_grid, 1)
    else
      call crash('unknown londir = "' // trim( londir) // '"!')
    end if

    ! latdir
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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_lonlat_file_3D_ocean

  ! Read data fields from a mesh file
  subroutine read_field_from_mesh_file_2D( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    real(dp), dimension(:),           intent(  out) :: d_mesh_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_2D'
    integer                               :: ncid
    type(type_mesh)                       :: mesh_loc
    integer                               :: id_var
    character(len=1024)                   :: var_name
    real(dp), dimension(:  ), allocatable :: d_mesh
    real(dp), dimension(:,:), allocatable :: d_mesh_with_time
    integer                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%master) allocate( d_mesh( mesh_loc%nV))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_master( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%master) allocate( d_mesh_with_time( mesh_loc%nV, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_master( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, ti /), count = (/ mesh_loc%nV, 1 /) )
      ! Copy to output memory
      if (par%master) d_mesh = d_mesh_with_time( :,1)
      ! Clean up after yourself
      if (par%master) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_2D

  subroutine read_field_from_mesh_file_2D_b( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    real(dp), dimension(:),           intent(  out) :: d_mesh_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_2D_b'
    integer                               :: ncid
    type(type_mesh)                       :: mesh_loc
    integer                               :: id_var
    character(len=1024)                   :: var_name
    real(dp), dimension(:  ), allocatable :: d_mesh
    real(dp), dimension(:,:), allocatable :: d_mesh_with_time
    integer                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%master) allocate( d_mesh( mesh_loc%nTri))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_master( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%master) allocate( d_mesh_with_time( mesh_loc%nTri, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_master( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, ti /), count = (/ mesh_loc%nTri, 1 /) )
      ! Copy to output memory
      if (par%master) d_mesh = d_mesh_with_time( :,1)
      ! Clean up after yourself
      if (par%master) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_2D_b

  subroutine read_field_from_mesh_file_2D_monthly( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D monthly data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_2D_monthly'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the file has a valid month dimension
    call check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%master) allocate( d_mesh( mesh_loc%nV, 12))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_master( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%master) allocate( d_mesh_with_time( mesh_loc%nV, 12, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_master( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, 12, 1 /) )
      ! Copy to output memory
      if (par%master) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%master) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_2D_monthly

  subroutine read_field_from_mesh_file_3D( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 3-D data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_3D'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: nzeta_loc
    real(dp), dimension(:), allocatable     :: zeta_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Set up the vertical coordinate zeta from the file
    call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%master) allocate( d_mesh( mesh_loc%nV, nzeta_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_master( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%master) allocate( d_mesh_with_time( mesh_loc%nV, nzeta_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_master( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, nzeta_loc, 1 /) )
      ! Copy to output memory
      if (par%master) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%master) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_3D

  subroutine read_field_from_mesh_file_3D_b( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 3-D data field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_3D'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: nzeta_loc
    real(dp), dimension(:), allocatable     :: zeta_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Set up the vertical coordinate zeta from the file
    call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%master) allocate( d_mesh( mesh_loc%nTri, nzeta_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_master( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%master) allocate( d_mesh_with_time( mesh_loc%nTri, nzeta_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_master( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nTri, nzeta_loc, 1 /) )
      ! Copy to output memory
      if (par%master) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%master) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_3D_b

  subroutine read_field_from_mesh_file_3D_ocean( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 3-D ocean data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_3D_ocean'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: ndepth_loc
    real(dp), dimension(:), allocatable     :: depth_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Set up the vertical coordinate depth from the file
    call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%master) allocate( d_mesh( mesh_loc%nV, ndepth_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_master( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%master) allocate( d_mesh_with_time( mesh_loc%nV, ndepth_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_master( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, ndepth_loc, 1 /) )
      ! Copy to output memory
      if (par%master) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%master) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_3D_ocean

  subroutine read_field_from_mesh_file_3D_CDF( filename, field_name_options, &
    d_mesh_partial)
    !< Read a cumulative density function field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_3D_CDF'
    integer                               :: ncid
    integer                               :: id_dim_vi, id_var, id_dim_bins, nbins_loc, nV_loc
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_mesh

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Get number of mesh vertices
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi, dim_length = nV_loc)

    ! Check that number of bins in file match the ones in the config file
    call inquire_dim_multopt( filename, ncid, 'bin', id_dim_bins, dim_length = nbins_loc)
    if (nbins_loc /= C%subgrid_bedrock_cdf_nbins) then
      call crash('number of CDF bins in external file ({int_01}) does not match subgrid_bedrock_cdf_nbins in config file!', int_01 = nbins_loc)
    end if

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! allocate memory
    if (par%master) allocate( d_mesh( nV_loc, C%subgrid_bedrock_cdf_nbins))

    ! Read data from file
    call read_var_master( filename, ncid, id_var, d_mesh)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_3D_CDF

  subroutine read_field_from_mesh_file_3D_b_CDF( filename, field_name_options, &
    d_mesh_partial)
    !< Read a cumulative density function field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_3D_b_CDF'
    integer                                 :: ncid
    integer                                 :: id_dim_ti, id_var, id_dim_bins, nbins_loc, nTri_loc
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Get number of mesh triangles
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti, dim_length = nTri_loc)

    ! Check that number of bins in file match the ones in the config file
    call inquire_dim_multopt( filename, ncid, 'bin', id_dim_bins, dim_length = nbins_loc)
    if (nbins_loc /= C%subgrid_bedrock_cdf_nbins) then
      call crash('number of CDF bins in external file ({int_01}) does not match subgrid_bedrock_cdf_nbins in config file!', int_01 = nbins_loc)
    end if

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! allocate memory
    if (par%master) allocate( d_mesh( nTri_loc, C%subgrid_bedrock_cdf_nbins))

    ! Read data from file
    call read_var_master( filename, ncid, id_var, d_mesh)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the master to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_master( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%master) deallocate( d_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_3D_b_CDF

  subroutine read_field_from_file_0D( filename, field_name_options, d, time_to_read)
    !< Read a 0-D data field from a NetCDF file

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    character(len=*),           intent(in   ) :: field_name_options
    real(dp),                   intent(  out) :: d
    real(dp),         optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'read_field_from_file_0D'
    integer                                :: ncid
    integer                                :: id_var, id_dim_time
    character(len=1024)                    :: var_name
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                :: ti
    real(dp), dimension(1)                 :: d_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    if (present( time_to_read)) then
      ! Assume the file has a time dimension, and we're reading from time_to_read

      ! Check if the file has a time dimension and variable
      call check_time( filename, ncid)

      ! Inquire variable info
      call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

      ! Inquire file time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

      ! Check if the variable has time as a dimension
      if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      if (.not. ANY( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      ! Inquire length of time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

      ! Read the data
      call read_var_master( filename, ncid, id_var, d_with_time, start = (/ ti /), count = (/ 1 /))
      if (par%master) d = d_with_time( 1)

      ! Broadcast to all processes
      call MPI_BCAST( d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    else ! if (present( time_to_read)) then
      ! Assume the file has no time dimension and we're just reading the data directly

      ! Inquire variable info
      call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

      ! Check if the variable has time as a dimension
      if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      ! Read the data
      call read_var_master( filename, ncid, id_var, d)

      ! Broadcast to all processes
      call MPI_BCAST( d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end if ! if (present( time_to_read)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_0D

  ! ===== Set up grids/mesh from a NetCDF file =====
  ! ================================================

  subroutine setup_xy_grid_from_file( filename, ncid, grid)
    !< Set up an x/y-grid from a NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    type(type_grid),  intent(  out) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_xy_grid_from_file'
    real(dp), parameter            :: tol = 1E-9_dp
    integer                        :: id_dim_x, id_dim_y
    integer                        :: id_var_x, id_var_y

    ! Add routine to path
    call init_routine( routine_name)

    ! Give the grid a nice name
    grid%name = 'xy_grid_from_file_"' // trim( filename) // '"'

    ! Check grid dimensions and variables for validity
    call check_x( filename, ncid)
    call check_y( filename, ncid)

    ! Inquire x and y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x, dim_length = grid%nx)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y, dim_length = grid%ny)

    ! allocate memory for x and y
    allocate( grid%x( grid%nx))
    allocate( grid%y( grid%ny))

    ! Inquire x and y variables
    call inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    call inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Read x and y
    call read_var_master(  filename, ncid, id_var_x, grid%x)
    call read_var_master(  filename, ncid, id_var_y, grid%y)

    ! Broadcast x and y from the master to the other processes
    call MPI_BCAST( grid%x, grid%nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( grid%y, grid%ny, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Calculate secondary grid geometry data
    call calc_secondary_grid_data( grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_xy_grid_from_file

  subroutine setup_lonlat_grid_from_file( filename, ncid, grid)
    !< Set up a lon/lat-grid from a NetCDF file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    type(type_grid_lonlat), intent(  out) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_lonlat_grid_from_file'
    integer                        :: id_dim_lon, id_dim_lat
    integer                        :: id_var_lon, id_var_lat

    ! Add routine to path
    call init_routine( routine_name)

    ! Give the grid a nice name
    grid%name = 'lonlat_grid_from_file_"' // trim( filename) // '"'

    ! Check grid dimensions and variables for validity
    call check_lon( filename, ncid)
    call check_lat( filename, ncid)

    ! Inquire lon and lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = grid%nlon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = grid%nlat)

    ! allocate memory for lon and lat
    allocate( grid%lon( grid%nlon))
    allocate( grid%lat( grid%nlat))

    ! Inquire lon and lat variables
    call inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read x and y
    call read_var_master( filename, ncid, id_var_lon, grid%lon)
    call read_var_master( filename, ncid, id_var_lat, grid%lat)

    ! Broadcast x and y from the master to the other processes
    call MPI_BCAST( grid%lon, grid%nlon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( grid%lat, grid%nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Secondary data
    call calc_lonlat_field_to_vector_form_translation_tables( grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_lonlat_grid_from_file

  subroutine setup_mesh_from_file( filename, ncid, mesh)
    !< Set up a mesh from a NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    type(type_mesh),  intent(  out) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_mesh_from_file'
    character(len=1024)            :: name
    integer                        :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    integer                        :: nV_mem, nTri_mem, nC_mem, n_two, n_three
    integer                        :: id_var_xmin, id_var_xmax, id_var_ymin, id_var_ymax, id_var_tol_dist, id_var_lambda_M, id_var_phi_M, id_var_beta_stereo
    integer                        :: id_var_V, id_var_nC, id_var_C, id_var_niTri, id_var_iTri, id_var_VBI
    integer                        :: id_var_Tri, id_var_Tricc, id_var_TriC
    real(dp), parameter            :: tol = 1E-9_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Give the mesh a nice name
    name = 'mesh_from_file_"' // trim( filename) // '"'

    ! Check mesh dimensions and variables for validity
    call check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   , dim_length = nV_mem  )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   , dim_length = nTri_mem)
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   , dim_length = nC_mem  )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  , dim_length = n_two   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three, dim_length = n_three )

    ! allocate memory for the mesh
    if (par%master) then
      call allocate_mesh_primary( mesh, name, nV_mem, nTri_mem, nC_mem)
      mesh%nV   = mesh%nV_mem
      mesh%nTri = mesh%nTri_mem
    end if

    ! == Inquire mesh variables
    ! =========================

    ! Metadata
    call inquire_var_multopt( filename, ncid, 'xmin'                           , id_var_xmin          )
    call inquire_var_multopt( filename, ncid, 'xmax'                           , id_var_xmax          )
    call inquire_var_multopt( filename, ncid, 'ymin'                           , id_var_ymin          )
    call inquire_var_multopt( filename, ncid, 'ymax'                           , id_var_ymax          )
    call inquire_var_multopt( filename, ncid, 'tol_dist'                       , id_var_tol_dist      )
    call inquire_var_multopt( filename, ncid, 'lambda_M'                       , id_var_lambda_M      )
    call inquire_var_multopt( filename, ncid, 'phi_M'                          , id_var_phi_M         )
    call inquire_var_multopt( filename, ncid, 'beta_stereo'                    , id_var_beta_stereo   )

    ! Vertex data
    call inquire_var_multopt( filename, ncid, field_name_options_V             , id_var_V             )
    call inquire_var_multopt( filename, ncid, field_name_options_nC            , id_var_nC            )
    call inquire_var_multopt( filename, ncid, field_name_options_C             , id_var_C             )
    call inquire_var_multopt( filename, ncid, field_name_options_niTri         , id_var_niTri         )
    call inquire_var_multopt( filename, ncid, field_name_options_iTri          , id_var_iTri          )
    call inquire_var_multopt( filename, ncid, field_name_options_VBI           , id_var_VBI           )

    ! Triangle data
    call inquire_var_multopt( filename, ncid, field_name_options_Tri           , id_var_Tri           )
    call inquire_var_multopt( filename, ncid, field_name_options_Tricc         , id_var_Tricc         )
    call inquire_var_multopt( filename, ncid, field_name_options_TriC          , id_var_TriC          )

    ! == Read mesh data
    ! =================

    ! Metadata
    call read_var_master(  filename, ncid, id_var_xmin          , mesh%xmin          )
    call read_var_master(  filename, ncid, id_var_xmax          , mesh%xmax          )
    call read_var_master(  filename, ncid, id_var_ymin          , mesh%ymin          )
    call read_var_master(  filename, ncid, id_var_ymax          , mesh%ymax          )
    call read_var_master(  filename, ncid, id_var_tol_dist      , mesh%tol_dist      )
    call read_var_master(  filename, ncid, id_var_lambda_M      , mesh%lambda_M      )
    call read_var_master(  filename, ncid, id_var_phi_M         , mesh%phi_M         )
    call read_var_master(  filename, ncid, id_var_beta_stereo   , mesh%beta_stereo   )

    ! Vertex data
    call read_var_master(  filename, ncid, id_var_V             , mesh%V             )
    call read_var_master( filename, ncid, id_var_nC            , mesh%nC            )
    call read_var_master( filename, ncid, id_var_C             , mesh%C             )
    call read_var_master( filename, ncid, id_var_niTri         , mesh%niTri         )
    call read_var_master( filename, ncid, id_var_iTri          , mesh%iTri          )
    call read_var_master( filename, ncid, id_var_VBI           , mesh%VBI           )

    ! Triangle data
    call read_var_master( filename, ncid, id_var_Tri           , mesh%Tri           )
    call read_var_master(  filename, ncid, id_var_Tricc         , mesh%Tricc         )
    call read_var_master( filename, ncid, id_var_TriC          , mesh%TriC          )

    ! Safety - check if the mesh data read from NetCDF makes sense
    if (par%master) call check_mesh( mesh)

    ! Broadcast read mesh from the master to the other processes
    call broadcast_mesh( mesh)

    ! Calculate secondary mesh data
    call calc_all_secondary_mesh_data( mesh, mesh%lambda_M, mesh%phi_M, mesh%beta_stereo)

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_mesh_from_file

  subroutine setup_zeta_from_file( filename, ncid, nzeta, zeta)
    !< Set up a zeta coordinate from a NetCDF file

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(  out) :: nzeta
    real(dp), dimension(:), allocatable, intent(  out) ::  zeta

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_zeta_from_file'
    integer                        :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    call init_routine( routine_name)

    ! Check zeta dimension and variable for validity
    call check_zeta( filename, ncid)

    ! Inquire zeta dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta, dim_length = nzeta)

    ! Inquire zeta variable
    call inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var_zeta)

    ! allocate memory
    allocate( zeta( nzeta))

    ! Read zeta from file
    call read_var_master( filename, ncid, id_var_zeta, zeta)

    ! Broadcast zeta from master to all other processes
    call MPI_BCAST( zeta, nzeta, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_zeta_from_file

  subroutine setup_depth_from_file( filename, ncid, ndepth, depth)
    !< Set up a depth coordinate from a NetCDF file

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(  out) :: ndepth
    real(dp), dimension(:), allocatable, intent(  out) ::  depth

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_depth_from_file'
    integer                        :: id_dim_depth, id_var_depth

    ! Add routine to path
    call init_routine( routine_name)

    ! Check depth dimension and variable for validity
    call check_depth( filename, ncid)

    ! Inquire depth dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth, dim_length = ndepth)

    ! Inquire depth variable
    call inquire_var_multopt( filename, ncid, field_name_options_depth, id_var_depth)

    ! allocate memory
    allocate( depth( ndepth))

    ! Read depth from file
    call read_var_master( filename, ncid, id_var_depth, depth)

    ! Broadcast depth from master to all other processes
    call MPI_BCAST( depth, ndepth, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_depth_from_file

  ! ===== Determine indexing and dimension directions =====
  ! =======================================================

  subroutine determine_xy_indexing( filename, ncid, var_name, indexing, xdir, ydir)
    !< Determine the indexing and dimension directions of a variable in an x/y-grid file

    ! In/output variables:
    character(len=*),    intent(in   ) :: filename
    integer,             intent(in   ) :: ncid
    character(len=*),    intent(in   ) :: var_name
    character(len=1024), intent(  out) :: indexing
    character(len=1024), intent(  out) :: xdir
    character(len=1024), intent(  out) :: ydir

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'determine_xy_indexing'
    integer                                :: id_dim_x, id_dim_y
    integer                                :: nx, ny
    real(dp), dimension(:), allocatable    :: x, y
    integer                                :: id_var_x, id_var_y, id_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if the x and y dimensions and variables of this file are valid
    call check_x( filename, ncid)
    call check_y( filename, ncid)

    ! Inquire x and y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x, dim_length = nx)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y, dim_length = ny)

    ! allocate memory for x and y
    allocate( x( nx))
    allocate( y( ny))

    ! Inquire x and y variables
    call inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    call inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Read x and y
    call read_var_master(  filename, ncid, id_var_x, x)
    call read_var_master(  filename, ncid, id_var_y, y)
    call MPI_BCAST( x, nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( y, ny, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Determine directions of x and y
    if (x( 2) > x( 1)) then
      xdir = 'normal'
    else
      xdir = 'reverse'
    end if
    if (y( 2) > y( 1)) then
      ydir = 'normal'
    else
      ydir = 'reverse'
    end if

    ! Inquire dimensions of the specified field variable
    call inquire_var_multopt( filename, ncid, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    if     (dims_of_var( 1) == id_dim_x .and. dims_of_var( 2) == id_dim_y) then
      indexing = 'xy'
    elseif (dims_of_var( 1) == id_dim_y .and. dims_of_var( 2) == id_dim_x) then
      indexing = 'yx'
    else
      call crash('x and y are not the first two dimensions of variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if

    ! Clean up after yourself
    deallocate( x)
    deallocate( y)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_xy_indexing

  subroutine determine_lonlat_indexing( filename, ncid, var_name, indexing, londir, latdir)
    !< Determine the indexing and dimension directions of a variable in a lon/lat-grid file

    ! In/output variables:
    character(len=*),    intent(in   ) :: filename
    integer,             intent(in   ) :: ncid
    character(len=*),    intent(in   ) :: var_name
    character(len=1024), intent(  out) :: indexing
    character(len=1024), intent(  out) :: londir
    character(len=1024), intent(  out) :: latdir

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'determine_lonlat_indexing'
    integer                                :: id_dim_lon, id_dim_lat
    integer                                :: nlon, nlat
    real(dp), dimension(:), allocatable    :: lon, lat
    integer                                :: id_var_lon, id_var_lat, id_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if the lon and lat dimensions and variables of this file are valid
    call check_lon( filename, ncid)
    call check_lat( filename, ncid)

    ! Inquire lon and lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = nlon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = nlat)

    ! allocate memory for lon and lat
    allocate( lon( nlon))
    allocate( lat( nlat))

    ! Inquire lon and lon variables
    call inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read lon and lat
    call read_var_master(  filename, ncid, id_var_lon, lon)
    call read_var_master(  filename, ncid, id_var_lat, lat)
    call MPI_BCAST( lon, nlon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( lat, nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Determine directions of x and y
    if (lon( 2) > lon( 1)) then
      londir = 'normal'
    else
      londir = 'reverse'
    end if
    if (lat( 2) > lat( 1)) then
      latdir = 'normal'
    else
      latdir = 'reverse'
    end if

    ! Inquire dimensions of the specified field variable
    call inquire_var_multopt( filename, ncid, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    if     (dims_of_var( 1) == id_dim_lon .and. dims_of_var( 2) == id_dim_lat) then
      indexing = 'lonlat'
    elseif (dims_of_var( 1) == id_dim_lat .and. dims_of_var( 2) == id_dim_lon) then
      indexing = 'latlon'
    else
      call crash('longitude and latitude are not the first two dimensions of variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if

    ! Clean up after yourself
    deallocate( lon)
    deallocate( lat)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_lonlat_indexing

end module netcdf_input
