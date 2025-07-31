module read_and_remap_field_from_file
!< Routines for flexibly reading 2-D, 2-D monthly, 3-D, and 3-D ocean data fields from NetCDF files.

! The files can contain the data on a regular x/y-grid, a regular lon/lat-grid, or a mesh.
! The routines will automatically detect which one of these it is, and select the
! appropriate subroutine. The data will also be automatically mapped to the provided model mesh.

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_DOUBLE_PRECISION
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, insert_val_into_string_int, warning
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid, type_grid_lonlat, type_grid_lat
  use remapping_main
  use grid_basic, only: deallocate_grid
  use grid_lonlat_basic, only: deallocate_lonlat_grid, deallocate_lat_grid
  use mesh_memory, only: deallocate_mesh
  use netcdf_basic
  use netcdf_setup_grid_mesh_from_file
  use netcdf_read_field_from_mesh_file
  use netcdf_read_field_from_series_file
  use netcdf_read_field_from_lonlat_grid_file
  use netcdf_read_field_from_xy_grid_file
  use netcdf_read_field_from_series_file
  use netcdf, only: NF90_MAX_VAR_DIMS

  implicit none

  private

  public :: read_field_from_file_2D, read_field_from_file_2D_monthly, read_field_from_file_3D, &
    read_field_from_file_3D_ocean, read_field_from_file_0D, read_field_from_file_2D_b

contains

  ! Read and map to mesh
  subroutine read_field_from_file_2D( filename, field_name_options, &
    mesh, output_dir, d_partial, time_to_read)
    !< Read a data field from a NetCDF file, and map it to the model mesh.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    character(len=*),       intent(in   ) :: field_name_options
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: output_dir
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
      call read_field_from_xy_file_dp_2D( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_xy_grid_to_mesh_2D( grid_from_file, mesh, output_dir, d_grid_vec_partial_from_file, d_partial)

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
      call read_field_from_mesh_file_dp_2D( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_mesh_to_mesh_2D( mesh_from_file, mesh, output_dir, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      call deallocate_mesh( mesh_from_file)
      deallocate( d_mesh_partial_from_file)

    else
      call crash('file "' // trim( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_2D

  subroutine read_field_from_file_2D_b( filename, field_name_options, &
    mesh, output_dir, d_partial, time_to_read)
    !< Read a data field from a NetCDF file, and map it to the model mesh triangles.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    character(len=*),       intent(in   ) :: field_name_options
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: output_dir
    real(dp), dimension(:), intent(  out) :: d_partial
    real(dp), optional,     intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'read_field_from_file_2D_b'
    logical                             :: file_exists
    logical                             :: has_xy_grid, has_lonlat_grid, has_mesh
    integer                             :: ncid
    type(type_grid)                     :: grid_from_file
    type(type_mesh)                     :: mesh_from_file
    real(dp), dimension(:), allocatable :: d_grid_vec_partial_from_file
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
      call read_field_from_xy_file_dp_2D( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_xy_grid_to_mesh_triangles_2D( grid_from_file, mesh, output_dir, d_grid_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_grid( grid_from_file)
      deallocate( d_grid_vec_partial_from_file)

    elseif (has_lonlat_grid) then
      ! Data is provided on a lon/lat-grid

      call crash('remapping from lon/lat-grid to mesh triangles not supported')

    elseif (has_mesh) then
      ! Data is provided on a mesh

      ! Set up the mesh from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_mesh_from_file( filename, ncid, mesh_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_mesh_partial_from_file( mesh_from_file%ti1: mesh_from_file%ti2))

      ! Read meshed data
      call read_field_from_mesh_file_dp_2D_b( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_mesh_tri_to_mesh_tri_2D( mesh_from_file, mesh, output_dir, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

      ! Clean up after yourself
      call deallocate_mesh( mesh_from_file)
      deallocate( d_mesh_partial_from_file)

    else
      call crash('file "' // trim( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_2D_b

  subroutine read_field_from_file_2D_monthly( filename, field_name_options, &
    mesh, output_dir, d_partial, time_to_read)
    !< Read a data field from a NetCDF file, and map it to the model mesh.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    ! In/output variables:
    character(len=*),         intent(in   ) :: filename
    character(len=*),         intent(in   ) :: field_name_options
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: output_dir
    real(dp), dimension(:,:), intent(  out) :: d_partial
    real(dp), optional,       intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_file_2D_monthly'
    logical                               :: file_exists
    logical                               :: has_xy_grid, has_lonlat_grid, has_mesh, has_lat_grid
    integer                               :: ncid
    type(type_grid)                       :: grid_from_file
    type(type_grid_lonlat)                :: grid_lonlat_from_file
    type(type_grid_lat)                   :: grid_lat_from_file
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
    call inquire_lat_grid(    filename, has_lat_grid)
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
      call read_field_from_xy_file_dp_2D_monthly( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, output_dir, d_grid_vec_partial_from_file, d_partial)

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

    elseif (has_lat_grid) then

      ! Data is provided on a lat-only grid
      ! Set up the grid from the file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_lonlat_grid_from_lat_file( filename, ncid, grid_lonlat_from_file, grid_lat_from_file)
      call close_netcdf_file( ncid)

      ! allocate memory for gridded data
      allocate( d_grid_lonlat_vec_partial_from_file( grid_lonlat_from_file%n1: grid_lonlat_from_file%n2,12))

      ! Read gridded data
      call read_field_from_lat_file_1D_monthly( filename, field_name_options, d_grid_lonlat_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_lonlat_grid_to_mesh_3D( grid_lonlat_from_file, mesh, d_grid_lonlat_vec_partial_from_file, d_partial)

      ! Clean up after yourself
      call deallocate_lonlat_grid( grid_lonlat_from_file)
      call deallocate_lat_grid( grid_lat_from_file)
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
      call read_field_from_mesh_file_dp_2D_monthly( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_mesh_to_mesh_3D( mesh_from_file, mesh, output_dir, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

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
    mesh, output_dir, d_partial, time_to_read, nzeta, zeta)
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
    character(len=*),                              intent(in   ) :: output_dir
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
      call read_field_from_xy_file_dp_3D( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, output_dir, d_grid_vec_partial_from_file, d_partial)

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
      call read_field_from_mesh_file_dp_3D( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! Remap data
      call map_from_mesh_to_mesh_3D( mesh_from_file, mesh, output_dir, d_mesh_partial_from_file, d_partial, method = method_mesh2mesh)

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
    mesh, output_dir, z_ocean, d_partial, time_to_read, ndepth, depth)
    !< Read a data field from a NetCDF file, and map it to the model mesh.

    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    ! In/output variables:
    character(len=*),                              intent(in   ) :: filename
    character(len=*),                              intent(in   ) :: field_name_options
    type(type_mesh),                               intent(in   ) :: mesh
    character(len=*),                              intent(in   ) :: output_dir
    real(dp), dimension(:),                        intent(in   ) :: z_ocean
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
      call read_field_from_xy_file_dp_3D_ocean( filename, field_name_options, d_grid_vec_partial_from_file, time_to_read = time_to_read)

      ! allocate memory for meshed data using all data layers
      allocate( d_partial_raw_layers( mesh%vi1:mesh%vi2, ndepth_loc))

      ! Remap data horizontally
      call map_from_xy_grid_to_mesh_3D( grid_from_file, mesh, output_dir, d_grid_vec_partial_from_file, d_partial_raw_layers)

      ! Remap data vertically
      call map_from_vertical_to_vertical_2D_ocean( mesh, depth_loc, z_ocean, d_partial_raw_layers, d_partial)

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
      call map_from_vertical_to_vertical_2D_ocean( mesh, depth_loc, z_ocean, d_partial_raw_layers, d_partial)

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
      call read_field_from_mesh_file_dp_3D_ocean( filename, field_name_options, d_mesh_partial_from_file, time_to_read = time_to_read)

      ! allocate memory for meshed data using all data layers
      allocate( d_partial_raw_layers( mesh%vi1:mesh%vi2, ndepth_loc))

      ! Remap data horizontally
      call map_from_mesh_to_mesh_3D( mesh_from_file, mesh, output_dir, d_mesh_partial_from_file, d_partial_raw_layers, method = method_mesh2mesh)

      ! Remap data vertically
      call map_from_vertical_to_vertical_2D_ocean( mesh, depth_loc, z_ocean, d_partial_raw_layers, d_partial)

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

  subroutine read_field_from_file_0D( filename, field_name_options, d, time_to_read)
    !< Read a 0-D data field from a NetCDF file

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    character(len=*),           intent(in   ) :: field_name_options
    real(dp),                   intent(  out) :: d
    real(dp),         optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'read_field_from_file_0D'
    integer                                :: ncid, ierr
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

      ! Find timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)

      ! Read the data
      call read_var_primary( filename, ncid, id_var, d_with_time, start = (/ ti /), count = (/ 1 /))
      if (par%primary) d = d_with_time( 1)

      ! Broadcast to all processes
      call MPI_BCAST( d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    else ! if (present( time_to_read)) then
      ! Assume the file has no time dimension and we're just reading the data directly

      ! Inquire variable info
      call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

      ! Check if the variable has time as a dimension
      if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      ! Read the data
      call read_var_primary( filename, ncid, id_var, d)

      ! Broadcast to all processes
      call MPI_BCAST( d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end if ! if (present( time_to_read)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_file_0D

end module read_and_remap_field_from_file
