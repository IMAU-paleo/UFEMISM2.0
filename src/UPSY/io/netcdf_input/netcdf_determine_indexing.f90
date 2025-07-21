module netcdf_determine_indexing
  !< Determine dimension indexing and directions in a NetCDF file

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_DOUBLE_PRECISION
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use netcdf_basic
  use netcdf, only: NF90_MAX_VAR_DIMS

  implicit none

  private

  public :: determine_xy_indexing, determine_lonlat_indexing, determine_lat_indexing

contains

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
    integer                                :: ierr
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
    call read_var_primary(  filename, ncid, id_var_x, x)
    call read_var_primary(  filename, ncid, id_var_y, y)
    call MPI_BCAST( x(:), nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( y(:), ny, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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
    integer                                :: ierr
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
    call read_var_primary(  filename, ncid, id_var_lon, lon)
    call read_var_primary(  filename, ncid, id_var_lat, lat)
    call MPI_BCAST( lon(:), nlon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( lat(:), nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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

  subroutine determine_lat_indexing( filename, ncid, var_name, latdir)

    !< Determine the indexing and dimension directions of a variable in a lat-only grid file

    ! In/output variables:
    character(len=*),    intent(in   ) :: filename
    integer,             intent(in   ) :: ncid
    character(len=*),    intent(in   ) :: var_name
    character(len=1024), intent(  out) :: latdir

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'determine_lat_indexing'
    integer                                :: id_dim_lat
    integer                                :: nlat
    real(dp), dimension(:), allocatable    :: lat
    integer                                :: ierr
    integer                                :: id_var_lat, id_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if the lat dimension and variables of this file are valid
    call check_lat( filename, ncid)

    ! Inquire lat dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = nlat)

    ! allocate memory for lon and lat
    allocate( lat( nlat))

    ! Inquire lon and lon variables
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read lon and lat
    call read_var_primary(  filename, ncid, id_var_lat, lat)
    call MPI_BCAST( lat, nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Determine direction of y
    if (lat( 2) > lat( 1)) then
      latdir = 'normal'
    else
      latdir = 'reverse'
    end if

    ! Inquire dimensions of the specified field variable
    call inquire_var_multopt( filename, ncid, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    if (dims_of_var( 1) == id_dim_lat) then
      ! do nothing, it's all good
    else
      call crash('latitude is not the first dimension of variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if

    ! Clean up after yourself
    deallocate( lat)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_lat_indexing

end module netcdf_determine_indexing
