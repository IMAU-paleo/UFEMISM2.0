module netcdf_basic

  ! Basic NetCDF routines
  ! =====================
  !
  ! These routines handle some very basic stuff; opening and closing
  ! NetCDF files in different modes, inquiring dimensions and variables,
  ! creating dimensions, variables and attributes;.
  !
  ! Also contains the routines for setting up an x/y-grid, lon/lat-grid,
  ! or a mesh from a NetCDF file.

! ===== Preamble =====
! ====================

  use assertions_basic
  use mpi
  use precisions                                             , only: dp, int8
  use mpi_basic                                              , only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging                  , only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C, git_commit_hash

  ! Import  NetCDF functionality
  use netcdf, only: NF90_NOERR, NF90_OPEN, NF90_CLOSE, NF90_NOWRITE, NF90_INQ_DIMID, NF90_inquire_dimension, &
              NF90_INQ_VARID, NF90_inquire_VARIABLE, NF90_MAX_VAR_DIMS, NF90_GET_VAR, &
              NF90_CREATE, NF90_NOCLOBBER, NF90_NETCDF4, NF90_endDEF, NF90_REDEF, NF90_DEF_DIM, NF90_DEF_VAR, &
              NF90_PUT_ATT, NF90_WRITE, NF90_INT, NF90_FLOAT, NF90_DOUBLE, NF90_PUT_VAR, NF90_UNLIMITED, &
              NF90_inquire_ATTRIBUTE, NF90_SHARE, NF90_GLOBAL, NF90_INT64

  implicit none

  ! NetCDF error code
  integer :: nerr

  ! Possible names for different dimensions and variables
  ! =====================================================

  ! Different options for the name of a dimension or variable can now be tried.
  ! They are separated by a double vertical bar ||

  ! dimensions
  character(len=1024), parameter :: field_name_options_x              = 'x||X||x1||X1||nx||NX||x-coordinate||X-coordinate||easting||Easting'
  character(len=1024), parameter :: field_name_options_y              = 'y||Y||y1||Y1||ny||NY||y-coordinate||Y-coordinate||northing||Northing'
  character(len=1024), parameter :: field_name_options_zeta           = 'zeta||Zeta'
  character(len=1024), parameter :: field_name_options_lon            = 'lon||Lon||long||Long||longitude||Longitude'
  character(len=1024), parameter :: field_name_options_lat            = 'lat||Lat||latitude||Latitude'
  character(len=1024), parameter :: field_name_options_time           = 'time||Time||t||nt'
  character(len=1024), parameter :: field_name_options_month          = 'month||Month'
  character(len=1024), parameter :: field_name_options_depth          = 'depth||Depth'

  ! Mesh
  character(len=1024), parameter :: field_name_options_dim_nV         = 'vi'
  character(len=1024), parameter :: field_name_options_dim_nTri       = 'ti'
  character(len=1024), parameter :: field_name_options_dim_nC_mem     = 'ci'
  character(len=1024), parameter :: field_name_options_dim_nE         = 'ei'
  character(len=1024), parameter :: field_name_options_dim_nVor       = 'vori'
  character(len=1024), parameter :: field_name_options_dim_two        = 'two'
  character(len=1024), parameter :: field_name_options_dim_three      = 'three'
  character(len=1024), parameter :: field_name_options_dim_four       = 'four'

  character(len=1024), parameter :: field_name_options_V              = 'V'
  character(len=1024), parameter :: field_name_options_Tri            = 'Tri'
  character(len=1024), parameter :: field_name_options_nC             = 'nC'
  character(len=1024), parameter :: field_name_options_C              = 'C'
  character(len=1024), parameter :: field_name_options_niTri          = 'niTri'
  character(len=1024), parameter :: field_name_options_iTri           = 'iTri'
  character(len=1024), parameter :: field_name_options_VBI            = 'VBI'
  character(len=1024), parameter :: field_name_options_Tricc          = 'Tricc'
  character(len=1024), parameter :: field_name_options_TriC           = 'TriC'
  character(len=1024), parameter :: field_name_options_TriBI          = 'TriBI'
  character(len=1024), parameter :: field_name_options_E              = 'E'
  character(len=1024), parameter :: field_name_options_VE             = 'VE'
  character(len=1024), parameter :: field_name_options_EV             = 'EV'
  character(len=1024), parameter :: field_name_options_ETri           = 'ETri'
  character(len=1024), parameter :: field_name_options_EBI            = 'EBI'
  character(len=1024), parameter :: field_name_options_vi2vori        = 'vi2vori'
  character(len=1024), parameter :: field_name_options_ti2vori        = 'ti2vori'
  character(len=1024), parameter :: field_name_options_ei2vori        = 'ei2vori'
  character(len=1024), parameter :: field_name_options_vori2vi        = 'vori2vi'
  character(len=1024), parameter :: field_name_options_vori2ti        = 'vori2ti'
  character(len=1024), parameter :: field_name_options_vori2ei        = 'vori2ei'
  character(len=1024), parameter :: field_name_options_Vor            = 'Vor'
  character(len=1024), parameter :: field_name_options_VornC          = 'VornC'
  character(len=1024), parameter :: field_name_options_VorC           = 'VorC'
  character(len=1024), parameter :: field_name_options_nVVor          = 'nVVor'
  character(len=1024), parameter :: field_name_options_VVor           = 'VVor'
  character(len=1024), parameter :: field_name_options_TriGC          = 'TriGC'
  character(len=1024), parameter :: field_name_options_A              = 'A'
  character(len=1024), parameter :: field_name_options_R              = 'R'

  ! Variables
  character(len=1024), parameter :: field_name_options_Hi             = 'Hi||thickness||lithk'
  character(len=1024), parameter :: field_name_options_Hb             = 'Hb||bed||topg'
  character(len=1024), parameter :: field_name_options_Hs             = 'Hs||surface||orog'
  character(len=1024), parameter :: field_name_options_SL             = 'SL'
  character(len=1024), parameter :: field_name_options_dHb            = 'dHb'
  character(len=1024), parameter :: field_name_options_Ti             = 'Ti'
  character(len=1024), parameter :: field_name_options_T_ocean        = 'T_ocean||t_ocean||t_an'
  character(len=1024), parameter :: field_name_options_S_ocean        = 'S_ocean||s_ocean||s_an'

contains

  ! Check if a file contains all the variables and dimensions describing
  ! an x/y-grid, a lon/lat-grid, or a mesh
  subroutine inquire_xy_grid( filename, has_xy_grid)
    ! inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular x/y-grid.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    logical,          intent(  out) :: has_xy_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_xy_grid'
    integer                        :: ncid
    integer                        :: id_dim_x, id_dim_y
    integer                        :: id_var_x, id_var_y

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for x and y dimensions and variables
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)
    call inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    call inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Check if everything is there
    has_xy_grid = .true.

    if (id_dim_x == -1) has_xy_grid = .false.
    if (id_dim_y == -1) has_xy_grid = .false.
    if (id_var_x == -1) has_xy_grid = .false.
    if (id_var_y == -1) has_xy_grid = .false.

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_xy_grid

  subroutine inquire_lonlat_grid( filename, has_lonlat_grid)
    ! inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular lon/lat-grid.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    logical,          intent(  out) :: has_lonlat_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_lonlat_grid'
    integer                        :: ncid
    integer                        :: id_dim_lon, id_dim_lat
    integer                        :: id_var_lon, id_var_lat

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for x and y dimensions and variables
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)
    call inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Check if everything is there
    has_lonlat_grid = .true.

    if (id_dim_lon == -1) has_lonlat_grid = .false.
    if (id_dim_lat == -1) has_lonlat_grid = .false.
    if (id_var_lon == -1) has_lonlat_grid = .false.
    if (id_var_lat == -1) has_lonlat_grid = .false.

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_lonlat_grid

  subroutine inquire_mesh( filename, has_mesh)
    ! inquire if a NetCDF file contains all the dimensions and variables
    ! describing a mesh.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    logical,          intent(  out) :: has_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_mesh'
    integer                        :: ncid
    integer                        :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    integer                        :: id_var_V, id_var_nC, id_var_C, id_var_niTri, id_var_iTri, id_var_VBI
    integer                        :: id_var_Tri, id_var_Tricc, id_var_TriC, id_var_TrIBI

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three)

    ! inquire mesh variables
    call inquire_var_multopt( filename, ncid, field_name_options_V         , id_var_V    )
    call inquire_var_multopt( filename, ncid, field_name_options_nC        , id_var_nC   )
    call inquire_var_multopt( filename, ncid, field_name_options_C         , id_var_C    )
    call inquire_var_multopt( filename, ncid, field_name_options_niTri     , id_var_niTri)
    call inquire_var_multopt( filename, ncid, field_name_options_iTri      , id_var_iTri )
    call inquire_var_multopt( filename, ncid, field_name_options_VBI       , id_var_VBI  )
    call inquire_var_multopt( filename, ncid, field_name_options_Tri       , id_var_Tri  )
    call inquire_var_multopt( filename, ncid, field_name_options_Tricc     , id_var_Tricc)
    call inquire_var_multopt( filename, ncid, field_name_options_TriC      , id_var_TriC )
    call inquire_var_multopt( filename, ncid, field_name_options_TriBI     , id_var_TriBI)

    ! Check if everything is there
    has_mesh = .true.

    if (id_dim_vi    == -1) has_mesh = .false.
    if (id_dim_ti    == -1) has_mesh = .false.
    if (id_dim_ci    == -1) has_mesh = .false.
    if (id_dim_two   == -1) has_mesh = .false.
    if (id_dim_three == -1) has_mesh = .false.

    if (id_var_V     == -1) has_mesh = .false.
    if (id_var_nC    == -1) has_mesh = .false.
    if (id_var_C     == -1) has_mesh = .false.
    if (id_var_niTri == -1) has_mesh = .false.
    if (id_var_iTri  == -1) has_mesh = .false.
    if (id_var_VBI   == -1) has_mesh = .false.
    if (id_var_Tri   == -1) has_mesh = .false.
    if (id_var_Tricc == -1) has_mesh = .false.
    if (id_var_TriC  == -1) has_mesh = .false.
    if (id_var_TriBI == -1) has_mesh = .false.

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_mesh

  ! inquire if a file contains the variable and dimension for
  ! zeta, z_ocean, time, or months
  subroutine inquire_zeta( filename, ncid, has_zeta)
    ! inquire if a NetCDF file contains a zeta dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    logical,          intent(  out) :: has_zeta

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_zeta'
    integer                        :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Look for zeta dimension and variable
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)
    call inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var_zeta)

    ! Check if everything is there
    has_zeta = .true.

    if (id_dim_zeta == -1) has_zeta = .false.
    if (id_var_zeta == -1) has_zeta = .false.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_zeta

  subroutine inquire_month( filename, ncid, has_month)
    ! inquire if a NetCDF file contains a month dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    logical,          intent(  out) :: has_month

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_month'
    integer                        :: id_dim_month, id_var_month

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Look for month dimension and variable
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)
    call inquire_var_multopt( filename, ncid, field_name_options_month, id_var_month)

    ! Check if everything is there
    has_month = .true.

    if (id_dim_month == -1) has_month = .false.
    if (id_var_month == -1) has_month = .false.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_month

  subroutine inquire_time( filename, ncid, has_time)
    ! inquire if a NetCDF file contains a time dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    logical,          intent(  out) :: has_time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_time'
    integer                        :: id_dim_time, id_var_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Look for time dimension and variable
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! Check if everything is there
    has_time = .true.

    if (id_dim_time == -1) has_time = .false.
    if (id_var_time == -1) has_time = .false.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_time

  subroutine find_timeframe( filename, ncid, time, ti)
    ! Find the timeframe in the file that is closest to the desired time.
    ! if the file has no time dimension or variable, throw an error.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    real(dp),         intent(in   ) :: time
    integer,          intent(  out) :: ti

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'find_timeframe'
    integer                             :: nt, id_dim_time, id_var_time
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: tii
    real(dp)                            :: dt_min

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file contains a valid time dimension and variable
    call check_time( filename, ncid)

    ! inquire size of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! inquire time variable ID
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! allocate memory
    allocate( time_from_file( nt))

    ! Read time from file
    call read_var_master_dp_1D( filename, ncid, id_var_time, time_from_file)
    call MPI_BCAST( time_from_file, nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Find timeframe closest to desired time
    if (time_from_file( 1) > time) then
      ! Desired time beyond lower limit
      call warning('desired timeframe at t = {dp_01} before start of file time for file "' &
        // trim( filename) // '"; reading data from t = {dp_02} instead!', &
        dp_01 = time, dp_02 = time_from_file( 1))
      ti = 1
    elseif (time_from_file( nt) < time) then
      ! Desired time beyond upper limit
      call warning('desired timeframe at t = {dp_01} after end of file time for file "' &
        // trim( filename) // '"; reading data from t = {dp_02} instead!', &
        dp_01 = time, dp_02 = time_from_file( nt))
      ti = nt
    else
      ! Desired time is within the file time
      dt_min = huge( 1._dp)
      do tii = 1, nt
        if (abs( time_from_file( tii) - time) < dt_min) then
          ti = tii
          dt_min = abs( time_from_file( tii) - time)
        end if
      end do
      if (dt_min > 0._dp) then
        call warning('desired timeframe at t = {dp_01} not present in file "' &
          // trim( filename) // '"; reading data from closest match at t = {dp_02} instead!', &
          dp_01 = time, dp_02 = time_from_file( ti))
      end if
    end if

    ! Clean up after yourself
    deallocate( time_from_file)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_timeframe

! ===== Safety checks on variables and dimensions =====
! =====================================================

  ! x/y-grid dimensions
  subroutine check_x( filename, ncid)
    ! Check if this file contains a valid x dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_x'
    integer                                 :: id_dim
    integer                                 :: n
    character(len=1024)                     :: dim_name
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    real(dp), dimension(:), allocatable     :: x
    real(dp)                                :: dx, dxp
    integer                                 :: i

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid x dimension could be found in file "' // trim( filename) // '"!')
    if (n == NF90_UNLIMITED) call crash('dimension "' // trim( dim_name) // '" in file "' // trim( filename) // '" is unlimited!')
    if (n < 1) call crash('dimension "' // trim( dim_name) // '" in file "' // trim( filename) // '" has length {int_01}!', int_01  = n)

    ! inquire variable
    call inquire_var_multopt( filename, ncid, field_name_options_x, id_var, &
      var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid x variable could be found in file "' // trim( filename) // '"!')

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) call crash('variable "' // trim( var_name) // &
      '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check variable dimension
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (dims_of_var( 1) /= id_dim) call crash('variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '" does not have ' // trim( dim_name) // ' as a dimension!')

    ! allocate memory
    allocate( x( n))

    ! Read variable
    call read_var_master_dp_1D( filename, ncid, id_var, x)

    if (par%master) call assert( (.not. any( isnan( x))), 'found NaNs in x')

    ! Check grid spacing
    if (par%master) then
      dx = x( 2) - x( 1)
      do i = 2, n
        dxp = x( i) - x( i-1)
        if (abs( 1._dp - dxp / dx) > 1E-5_dp) call crash('x coordinate in file "' // trim( filename) // '" is irregular!')
      end do
    end if ! if (par%master) then
    call sync

    ! Clean up after yourself
    deallocate( x)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_x

  subroutine check_y( filename, ncid)
    ! Check if this file contains a valid y dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_y'
    integer                                 :: id_dim
    integer                                 :: n
    character(len=1024)                     :: dim_name
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    real(dp), dimension(:), allocatable     :: y
    real(dp)                                :: dy, dyp
    integer                                 :: i

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid y dimension could be found in file "' // trim( filename) // '"!')
    if (n == NF90_UNLIMITED) call crash('dimension "' // trim( dim_name) // '" in file "' // trim( filename) // '" is unlimited!')
    if (n < 1) call crash('dimension "' // trim( dim_name) // '" in file "' // trim( filename) // '" has length {int_01}!', int_01  = n)

    ! inquire variable
    call inquire_var_multopt( filename, ncid, field_name_options_y, id_var, &
      var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid y variable could be found in file "' // trim( filename) // '"!')

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) call crash('variable "' // trim( var_name) // &
      '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check variable dimension
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (dims_of_var( 1) /= id_dim) call crash('variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '" does not have ' // trim( dim_name) // ' as a dimension!')

    ! allocate memory
    allocate( y( n))

    ! Read variable
    call read_var_master_dp_1D( filename, ncid, id_var, y)

    if (par%master) call assert( (.not. any( isnan( y))), 'found NaNs in y')

    ! Check grid spacing
    if (par%master) then
      dy = y( 2) - y( 1)
      do i = 2, n
        dyp = y( i) - y( i-1)
        if (abs( 1._dp - dyp / dy) > 1E-5_dp) call crash('y coordinate in file "' // trim( filename) // '" is irregular!')
      end do
    end if
    call sync

    ! Clean up after yourself
    deallocate( y)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_y

  ! lon/lat-grid dimensions
  subroutine check_lon( filename, ncid)
    ! Check if this file contains a valid longitude dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_lon'
    integer                                 :: id_dim
    integer                                 :: n
    character(len=1024)                     :: dim_name
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    real(dp), dimension(:), allocatable     :: lon
    real(dp)                                :: dlon, dlonp
    integer                                 :: i

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid longitude dimension could be found in file "' // trim( filename) // '"!')
    if (n == NF90_UNLIMITED) call crash('longitude dimension in file "' // trim( filename) // '" is unlimited!')
    if (n < 1) call crash('longitude dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n)

    ! inquire variable
    call inquire_var_multopt( filename, ncid, field_name_options_lon, id_var, &
      var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid longitude variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) call crash('longitude variable in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 1) call crash('longitude variable in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (dims_of_var( 1) /= id_dim) call crash('longitude variable in file "' // trim( filename) // '" does not have longitude as a dimension!')

    ! allocate memory
    allocate( lon( n))

    ! Read variable
    call read_var_master_dp_1D( filename, ncid, id_var, lon)

    if (par%master) call assert( (.not. any( isnan( lon))), 'found NaNs in lon')

    ! Check grid spacing
    if (par%master) then
      dlon = lon( 2) - lon( 1)
      do i = 2, n
        dlonp = lon( i) - lon( i-1)
        if (abs( 1._dp - dlonp / dlon) > 1E-5_dp) call crash('longitude coordinate in file "' // trim( filename) // '" is irregular!')
      end do
    end if ! if (par%master) then
    call sync

    ! Clean up after yourself
    deallocate( lon)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lon

  subroutine check_lat( filename, ncid)
    ! Check if this file contains a valid latitude dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_lat'
    integer                                 :: id_dim
    integer                                 :: n
    character(len=1024)                     :: dim_name
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    real(dp), dimension(:), allocatable     :: lat
    real(dp)                                :: dlat, dlatp
    integer                                 :: i

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid latitude dimension could be found in file "' // trim( filename) // '"!')
    if (n == NF90_UNLIMITED) call crash('latitude dimension in file "' // trim( filename) // '" is unlimited!')
    if (n < 1) call crash('latitude dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n)

    ! inquire variable
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid latitude variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) call crash('latitude variable in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 1) call crash('latitude variable in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (dims_of_var( 1) /= id_dim) call crash('latitude variable in file "' // trim( filename) // '" does not have latitude as a dimension!')

    ! allocate memory
    allocate( lat( n))

    ! Read variable
    call read_var_master_dp_1D( filename, ncid, id_var, lat)

    if (par%master) call assert( (.not. any( isnan( lat))), 'found NaNs in lat')

    ! Check grid spacing
    if (par%master) then
      dlat = lat( 2) - lat( 1)
      do i = 2, n
        dlatp = lat( i) - lat( i-1)
        if (abs( 1._dp - dlatp / dlat) > 1E-5_dp) call crash('latitude coordinate in file "' // trim( filename) // '" is irregular!')
      end do
    end if ! if (par%master) then
    call sync

    ! Clean up after yourself
    deallocate( lat)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lat

  ! Mesh dimensions
  subroutine check_mesh_dimensions( filename, ncid)
    ! Check if this file contains valid mesh dimensions and variables

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_mesh_dimensions'
    integer                                 :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    integer                                 :: nV, nTri, nC_mem, n_two, n_three
    character(len=1024)                     :: dim_name_vi, dim_name_ti, dim_name_ci, dim_name_two, dim_name_three
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! == inquire dimensions
    ! =====================

    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   , dim_length = nV     , dim_name = dim_name_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   , dim_length = nTri   , dim_name = dim_name_ti   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   , dim_length = nC_mem , dim_name = dim_name_ci   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  , dim_length = n_two  , dim_name = dim_name_two  )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three, dim_length = n_three, dim_name = dim_name_three)

    ! Safety checks on dimensions
    if (id_dim_vi    == -1) call crash('no valid vi    dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_ti    == -1) call crash('no valid ti    dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_ci    == -1) call crash('no valid ci    dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_two   == -1) call crash('no valid two   dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_three == -1) call crash('no valid three dimension could be found in file "' // trim( filename) // '"!')

    if (nV      == NF90_UNLIMITED) call crash('vi    dimension in file "' // trim( filename) // '" is unlimited!')
    if (nTri    == NF90_UNLIMITED) call crash('ti    dimension in file "' // trim( filename) // '" is unlimited!')
    if (nC_mem  == NF90_UNLIMITED) call crash('ci    dimension in file "' // trim( filename) // '" is unlimited!')
    if (n_two   == NF90_UNLIMITED) call crash('two   dimension in file "' // trim( filename) // '" is unlimited!')
    if (n_three == NF90_UNLIMITED) call crash('three dimension in file "' // trim( filename) // '" is unlimited!')

    if (nV      <  1) call crash('vi    dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = nV     )
    if (nTri    <  1) call crash('ti    dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = nTri   )
    if (nC_mem  <  1) call crash('ci    dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = nC_mem )
    if (n_two   /= 2) call crash('two   dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n_two  )
    if (n_three /= 3) call crash('three dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n_three)

    ! == inquire variables
    ! ====================

    ! Metadata
    ! ========

    ! xmin
    call inquire_var_multopt( filename, ncid, 'xmin', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid xmin variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! xmax
    call inquire_var_multopt( filename, ncid, 'xmax', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid xmax variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! ymin
    call inquire_var_multopt( filename, ncid, 'ymin', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid ymin variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! ymax
    call inquire_var_multopt( filename, ncid, 'ymax', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid ymax variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! tol_dist
    call inquire_var_multopt( filename, ncid, 'tol_dist', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid tol_dist variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! lambda_M
    call inquire_var_multopt( filename, ncid, 'lambda_M', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid lambda_M variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! phi_M
    call inquire_var_multopt( filename, ncid, 'phi_M', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid phi_M variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! beta_stereo
    call inquire_var_multopt( filename, ncid, 'beta_stereo', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    if (id_var == -1) call crash('no valid beta_stereo variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Vertex data
    ! ===========

    ! V
    call inquire_var_multopt( filename, ncid, field_name_options_V, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid V variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. (dims_of_var( 1) == id_dim_vi .and. dims_of_var( 2) == id_dim_two)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi and two as dimensions!')

    ! nC
    call inquire_var_multopt( filename, ncid, field_name_options_nC, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid nC variable could be found in file "' // trim( filename) // '"!')
    if (.not. var_type == NF90_INT) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. dims_of_var( 1) == id_dim_vi) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')

    ! C
    call inquire_var_multopt( filename, ncid, field_name_options_C, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid C variable could be found in file "' // trim( filename) // '"!')
    if (.not. var_type == NF90_INT) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. (dims_of_var( 1) == id_dim_vi .and. dims_of_var( 2) == id_dim_ci)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi and ci as dimensions!')

    ! niTri
    call inquire_var_multopt( filename, ncid, field_name_options_niTri, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid niTri variable could be found in file "' // trim( filename) // '"!')
    if (.not. var_type == NF90_INT) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. dims_of_var( 1) == id_dim_vi) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')

    ! iTri
    call inquire_var_multopt( filename, ncid, field_name_options_iTri, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid iTri variable could be found in file "' // trim( filename) // '"!')
    if (.not. var_type == NF90_INT) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. (dims_of_var( 1) == id_dim_vi .and. dims_of_var( 2) == id_dim_ci)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi and ci as dimensions!')

    ! VBI
    call inquire_var_multopt( filename, ncid, field_name_options_VBI, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid VBI variable could be found in file "' // trim( filename) // '"!')
    if (.not. var_type == NF90_INT) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. dims_of_var( 1) == id_dim_vi) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')

    ! Triangle data
    ! =============

    ! Tri
    call inquire_var_multopt( filename, ncid, field_name_options_Tri, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid Tri variable could be found in file "' // trim( filename) // '"!')
    if (.not. var_type == NF90_INT) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. (dims_of_var( 1) == id_dim_ti .and. dims_of_var( 2) == id_dim_three)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti and three as dimensions!')

    ! Tricc
    call inquire_var_multopt( filename, ncid, field_name_options_Tricc, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid Tricc variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. (dims_of_var( 1) == id_dim_ti .and. dims_of_var( 2) == id_dim_two)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti and two as dimensions!')

    ! TriC
    call inquire_var_multopt( filename, ncid, field_name_options_TriC, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid TriC variable could be found in file "' // trim( filename) // '"!')
    if (.not. var_type == NF90_INT) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. (dims_of_var( 1) == id_dim_ti .and. dims_of_var( 2) == id_dim_three)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti and three as dimensions!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_dimensions

  ! Zeta, z_ocean, month, time dimensions
  subroutine check_zeta( filename, ncid)
    ! Check if this file contains a valid zeta dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_zeta'
    integer                                 :: id_dim
    integer                                 :: n
    character(len=1024)                     :: dim_name
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    real(dp), dimension(:), allocatable     :: zeta
    integer                                 :: k

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid zeta dimension could be found in file "' // trim( filename) // '"!')
    if (n == NF90_UNLIMITED) call crash('zeta dimension in file "' // trim( filename) // '" is unlimited!')
    if (n < 1) call crash('zeta dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n)

    ! inquire variable
    call inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var, &
      var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid zeta variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) call crash('zeta variable in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 1) call crash('zeta variable in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (dims_of_var( 1) /= id_dim) call crash('zeta variable in file "' // trim( filename) // '" does not have zeta as a dimension!')

    ! allocate memory
    allocate( zeta( n))

    ! Read variable
    call read_var_master_dp_1D( filename, ncid, id_var, zeta)

    ! Check validity
    if (par%master) then
      call assert( (.not. any( isnan( zeta))), 'found NaNs in zeta')

      if (zeta( 1) /= 0._dp) call crash('zeta in file "' // trim( filename) // '" does not start at zero!')
      if (zeta( n) /= 1._dp) call crash('zeta in file "' // trim( filename) // '" does not end at one!')

      do k = 2, n
        if (zeta( k) <= zeta( k-1)) call crash('zeta in file "' // trim( filename) // '" does not increase monotonously!')
      end do
    end if
    call sync

    ! Clean up after yourself
    deallocate( zeta)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_zeta

  subroutine check_month( filename, ncid)
    ! Check if this file contains a valid month dimension (we don't really care about the variable)

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'check_month'
    integer                        :: id_dim
    integer                        :: n
    character(len=1024)            :: dim_name

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid month dimension could be found in file "' // trim( filename) // '"!')
    if (n == NF90_UNLIMITED) call crash('month dimension in file "' // trim( filename) // '" is unlimited!')
    if (n /= 12) call crash('month dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_month

  subroutine check_time( filename, ncid)
    ! Check if this file contains a valid time dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_time'
    integer                                 :: id_dim
    integer                                 :: n
    character(len=1024)                     :: dim_name
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    real(dp), dimension(:), allocatable     :: time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid time dimension could be found in file "' // trim( filename) // '"!')
    if (n < 0) call crash('time dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n)

    ! inquire variable
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var, &
      var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid time variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) call crash('time variable in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 1) call crash('time variable in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (dims_of_var( 1) /= id_dim) call crash('time variable in file "' // trim( filename) // '" does not have time as a dimension!')

    ! For new output files, time is still empty. if it's not, check if entries are valid
    if (n > 0) then

      ! allocate memory
      allocate( time( n))

      ! Read variable
      call read_var_master_dp_1D( filename, ncid, id_var, time)

      ! Check validity
      if (par%master) call assert( (.not. any( isnan( time))), 'found NaN in time')

      ! Clean up after yourself
      deallocate( time)

    end if ! if (n > 0) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_time

  subroutine check_depth( filename, ncid)
    ! Check if this file contains a valid depth dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'check_depth'
    integer                                 :: id_dim
    integer                                 :: n
    character(len=1024)                     :: dim_name
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    real(dp), dimension(:), allocatable     :: depth
    integer                                 :: k

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    if (id_dim == -1) call crash('no valid depth dimension could be found in file "' // trim( filename) // '"!')
    if (n == NF90_UNLIMITED) call crash('depth dimension in file "' // trim( filename) // '" is unlimited!')
    if (n < 1) call crash('depth dimension in file "' // trim( filename) // '" has length n = {int_01}!', int_01  = n)

    ! inquire variable
    call inquire_var_multopt( filename, ncid, field_name_options_depth, id_var, &
      var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no valid depth variable could be found in file "' // trim( filename) // '"!')
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) call crash('depth variable in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    if (ndims_of_var /= 1) call crash('depth variable in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (dims_of_var( 1) /= id_dim) call crash('depth variable in file "' // trim( filename) // '" does not have depth as a dimension!')

    ! allocate memory
    allocate( depth( n))

    ! Read variable
    call read_var_master_dp_1D( filename, ncid, id_var, depth)

    ! Check validity
    if (par%master) then
      call assert( (.not. any( isnan( depth))), 'found NaNs in depth')

      do k = 2, n
        if (depth( k) <= depth( k-1)) call crash('depth in file "' // trim( filename) // '" does not increase monotonously!')
      end do
    end if
    call sync

    ! Clean up after yourself
    deallocate( depth)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_depth

  ! x/y-grid field variables
  subroutine check_xy_grid_field_int_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_int_2D'
    integer                                :: id_dim_x, id_dim_y, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x( filename, ncid)
    call check_y( filename, ncid)

    ! inquire x,y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has x,y as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_int_2D

  subroutine check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_2D'
    integer                                :: id_dim_x, id_dim_y, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x( filename, ncid)
    call check_y( filename, ncid)

    ! inquire x,y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has x,y as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_2D

  subroutine check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_2D_monthly'
    integer                                :: id_dim_x, id_dim_y, id_dim_month, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x(     filename, ncid)
    call check_y(     filename, ncid)
    call check_month( filename, ncid)

    ! inquire x,y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x    )
    call inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y    )
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')
    if (.not. any( dims_of_var == id_dim_month)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have month as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has x,y,m as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y,m as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y,m as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_2D_monthly

  subroutine check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_3D'
    integer                                :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x(    filename, ncid)
    call check_y(    filename, ncid)
    call check_zeta( filename, ncid)

    ! inquire x,y,zeta dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    call inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has x,y,zeta as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y,zeta as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y,zeta as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_3D

  subroutine check_xy_grid_field_dp_3D_ocean(       filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_3D_ocean'
    integer                                :: id_dim_x, id_dim_y, id_dim_depth, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x(     filename, ncid)
    call check_y(     filename, ncid)
    call check_depth( filename, ncid)

    ! inquire x,y,depth dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x   )
    call inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y   )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')
    if (.not. any( dims_of_var == id_dim_depth)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have depth as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has x,y,depth as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y,depth as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y,depth as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_3D_ocean

  ! lon/lat-grid field variables
  subroutine check_lonlat_grid_field_int_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_int_2D'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon( filename, ncid)
    call check_lat( filename, ncid)

    ! inquire lon,lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has lon,lat as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_int_2D

  subroutine check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_2D'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon( filename, ncid)
    call check_lat( filename, ncid)

    ! inquire lon,lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has lon,lat as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_2D

  subroutine check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_2D_monthly'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_month, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon(   filename, ncid)
    call check_lat(   filename, ncid)
    call check_month( filename, ncid)

    ! inquire lon,lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon  )
    call inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat  )
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_month)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have month as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has lon,lat,m as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat,m as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat,m as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_2D_monthly

  subroutine check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_3D'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon(  filename, ncid)
    call check_lat(  filename, ncid)
    call check_zeta( filename, ncid)

    ! inquire lon,lat,zeta dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    call inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has lon,lat,zeta as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat,zeta as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat,zeta as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_3D

  subroutine check_lonlat_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_3D_ocean'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_depth, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon(   filename, ncid)
    call check_lat(   filename, ncid)
    call check_depth( filename, ncid)

    ! inquire lon,lat,depth dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon )
    call inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_depth)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have depth as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has lon,lat,depth as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat,depth as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat,depth as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_3D_ocean

  ! mesh field variables
  subroutine check_mesh_field_int_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_int_2D'
    integer                                :: id_dim_vi, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_int_2D

  subroutine check_mesh_field_int_2D_b( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_int_2D_b'
    integer                                :: id_dim_ti, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ti)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_int_2D_b

  subroutine check_mesh_field_int_2D_c( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_int_2D_c'
    integer                                :: id_dim_ei, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ei)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ei as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_int_2D_c

  subroutine check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D'
    integer                                :: id_dim_vi, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D

  subroutine check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D_b'
    integer                                :: id_dim_ti, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ti)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D_b

  subroutine check_mesh_field_dp_2D_c( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D_c'
    integer                                :: id_dim_ei, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ei)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ei as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D_c

  subroutine check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly mesh variable by this name
    !
    ! NOTE: this is 2-D monthly in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D_monthly'
    integer                                :: id_dim_vi, id_dim_month, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_month(           filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_month , id_dim_month)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')
    if (.not. any( dims_of_var == id_dim_month)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have month as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has vi,m as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi,m as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi,m as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D_monthly

  subroutine check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_3D'
    integer                                :: id_dim_vi, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_zeta(            filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta  , id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has vi,zeta as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi,zeta as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi,zeta as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_3D

  subroutine check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_3D'
    integer                                :: id_dim_ti, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_zeta(            filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta    , id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ti  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has ti,zeta as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have ti,zeta as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have ti,zeta as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_3D_b

  subroutine check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_3D_ocean'
    integer                                :: id_dim_vi, id_dim_depth, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_depth(           filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth , id_dim_depth)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')
    if (.not. any( dims_of_var == id_dim_depth)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have depth as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has vi,depth as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi,depth as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi,depth as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_3D_ocean

! ===== Flexible looking for dimensions and variables =====
! =========================================================

  ! Look for dimensions
  subroutine inquire_dim_multopt( filename, ncid, dim_name_options, id_dim, dim_length, dim_name)
    ! inquire if this file contains a dimension by name of dim_name.
    ! if so, return its length and identifier. if not, return -1 for both.
    !
    ! Supports providing multiple options for the dimension name, separated by two
    ! vertical bars || e.g. if we're looking for an X-dimension, we could do something like:
    !
    ! call inquire_dim_multopt( ncid, dim_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', dim_length, id_dim)
    !
    ! if more than one match is found, crash.

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    character(len=*),                    intent(in   ) :: dim_name_options
    integer,                             intent(  out) :: id_dim
    integer,             optional,       intent(  out) :: dim_length
    character(len=1024), optional,       intent(  out) :: dim_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_dim_multopt'
    character(len=1024)            :: dim_name_options_parsed
    character(len=1024)            :: dim_name_options_redux
    integer                        :: i, n_matches
    integer                        :: dim_length_try, dim_length_match
    integer                        :: id_dim_try, id_dim_match
    character(len=1024)            :: dim_name_try, dim_name_match

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Parse field name options
    call parse_field_name_options( dim_name_options, dim_name_options_parsed)

    ! Try all options provided in dim_name_options

    dim_name_options_redux = trim( dim_name_options_parsed)
    n_matches = 0

    do while (.true.)

      i = index( dim_name_options_redux, '||')

      if (i > 0) then
        ! More than one option is left over; take the last one

        dim_name_try = dim_name_options_redux( 1:i-1)
        dim_name_options_redux = dim_name_options_redux( i+2:len_trim( dim_name_options_redux))

      else
        ! Only one option is left over

        dim_name_try = dim_name_options_redux
        dim_name_options_redux( 1:len( dim_name_options_redux)) = ''

      end if

      ! Try the selected name option
      call inquire_dim( filename, ncid, dim_name_try, dim_length_try, id_dim_try)

      if (id_dim_try == -1) then
        ! No dimension by this name was found; try the next option
      else
        ! A dimension by this name was found; hurray!
        n_matches  = n_matches + 1
        dim_length_match = dim_length_try
        id_dim_match     = id_dim_try
        dim_name_match   = dim_name_try
      end if

      ! if the list of options is now empty, exit
      if (len_trim( dim_name_options_redux) == 0) EXIT

    end do

    if (n_matches == 0) then
      ! None of the optional dimension names were found in the NetCDF file
      dim_length_match = -1
      id_dim_match     = -1
    elseif (n_matches > 1) then
      ! More than one match was found
      call crash('more than one of the provided dimension names were found in file "' // trim( filename) // '"!')
    else
      ! We found exactly one match; hurray!
    end if

    ! Copy to output arguments
    id_dim = id_dim_match
    if (present( dim_name  )) dim_name   = dim_name_match
    if (present( dim_length)) dim_length = dim_length_match

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_dim_multopt

  ! Look for variables
  subroutine inquire_var_multopt( filename, ncid, var_name_options, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    ! inquire if this file contains a variable by name of var_name.
    ! if so, return its identifier. if not, return -1.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! call inquire_var_multopt( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', id_var)
    !
    ! if more than one match is found, crash.

    ! In/output variables:
    character(len=*),                                 intent(in   ) :: filename
    integer,                                          intent(in   ) :: ncid
    character(len=*),                                 intent(in   ) :: var_name_options
    integer,                                          intent(  out) :: id_var
    character(len=1024),                    optional, intent(  out) :: var_name
    integer,                                optional, intent(  out) :: var_type
    integer,                                optional, intent(  out) :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS), optional, intent(  out) :: dims_of_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_var_multopt'
    character(len=1024)            :: var_name_options_parsed
    character(len=1024)            :: var_name_options_redux
    integer                        :: i, n_matches, id_var_try
    character(len=1024)            :: var_name_try, var_name_match

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Parse field name options
    call parse_field_name_options( var_name_options, var_name_options_parsed)

    ! Try all options provided in var_name_options

    var_name_options_redux = trim( var_name_options_parsed)
    n_matches = 0

    do while (.true.)

      i = index( var_name_options_redux, '||')

      if (i > 0) then
        ! More than one option is left over; take the last one

        var_name_try = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:len_trim( var_name_options_redux))

      else
        ! Only one option is left over

        var_name_try = trim( var_name_options_redux)
        var_name_options_redux( 1:len( var_name_options_redux)) = ''

      end if

      ! Try the selected name option
      call inquire_var( filename, ncid, var_name_try, id_var_try)

      if (id_var_try == -1) then
        ! No variable by this name was found; try the next option
      else
        ! A variable by this name was found; hurray!
        n_matches      = n_matches + 1
        id_var         = id_var_try
        var_name_match = var_name_try
      end if

      ! if the list of options is now empty, exit
      if (len_trim( var_name_options_redux) == 0) EXIT

    end do

    if (n_matches == 0) then
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    elseif (n_matches > 1) then
      ! More than one match was found
      call crash('more than one of the provided variable names were found in file "' // trim( filename) // '"!')
    else
      ! We found exactly one match. inquire additional info on this variable.
      call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    end if

    ! Copy to output arguments
    if (present( var_name)) var_name = var_name_match

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_var_multopt

! ===== Parse flexible dimension/variable names =====
! ===================================================

  subroutine parse_field_name_options( field_name_options, field_name_options_parsed)
    ! Check if a default set of field name options should be used.

    ! In/output variables:
    character(len=*),    intent(in   ) :: field_name_options
    character(len=1024), intent(  out) :: field_name_options_parsed

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'parse_field_name_options'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    field_name_options_parsed = field_name_options

    if (index( field_name_options,'default_options_') > 0) then
      ! Use one of the default options

      ! dimensions
      if     (field_name_options == 'default_options_x') then
        field_name_options_parsed = field_name_options_x
      elseif (field_name_options == 'default_options_y') then
        field_name_options_parsed = field_name_options_y
      elseif (field_name_options == 'default_options_zeta') then
        field_name_options_parsed = field_name_options_zeta
      elseif (field_name_options == 'default_options_lon') then
        field_name_options_parsed = field_name_options_lon
      elseif (field_name_options == 'default_options_lat') then
        field_name_options_parsed = field_name_options_lat
      elseif (field_name_options == 'default_options_time') then
        field_name_options_parsed = field_name_options_time

      ! Variables
      elseif (field_name_options == 'default_options_Hi') then
        field_name_options_parsed = field_name_options_Hi
      elseif (field_name_options == 'default_options_Hb') then
        field_name_options_parsed = field_name_options_Hb
      elseif (field_name_options == 'default_options_Hs') then
        field_name_options_parsed = field_name_options_Hs
      elseif (field_name_options == 'default_options_SL') then
        field_name_options_parsed = field_name_options_SL

      ! Unrecognised default options
      else
        call crash('unregocnised default field name option "' // trim( field_name_options) // '"')
      end if

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine parse_field_name_options

  function get_first_option_from_list( field_name_options) result( field_name)
    ! Get the first option from a list of field name options

    ! In/output variables:
    character(len=*),  intent(in   ) :: field_name_options
    character(len=1024)              :: field_name

    ! Local variables:
    integer :: i

    field_name( 1:256) = ' '

    i = index( field_name_options,'||')

    if (i > 0) then
      field_name = field_name_options( 1:i-1)
    else
      field_name = trim( field_name_options)
    end if

  end function get_first_option_from_list

! ===== Read data from variables =====
! ====================================

  ! NOTE: only the Master actually reads data! Distributing to other processes
  !       must be done afterward

  subroutine read_var_master_int_0D(  filename, ncid, id_var, d)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    integer,          intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_var_master_int_0D'
    character(len=1024)            :: var_name
    integer                        :: var_type
    integer                        :: ndims_of_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_int_0D

  subroutine read_var_master_int_1D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                     intent(in   ) :: filename
    integer,                              intent(in   ) :: ncid
    integer,                              intent(in   ) :: id_var
    integer,  dimension(:    ), optional, intent(  out) :: d
    integer,  dimension(1    ), optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_int_1D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(1)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_int_1D

  subroutine read_var_master_int_2D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    integer,                            intent(in   ) :: ncid
    integer,                            intent(in   ) :: id_var
    integer,  dimension(:,:), optional, intent(  out) :: d
    integer,  dimension(2),   optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_int_2D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(2)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_int_2D

  subroutine read_var_master_int_3D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(in   ) :: id_var
    integer, dimension(:,:,:), optional, intent(  out) :: d
    integer, dimension(3),     optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_int_3D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(3)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_int_3D

  subroutine read_var_master_int_4D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                      intent(in   ) :: filename
    integer,                               intent(in   ) :: ncid
    integer,                               intent(in   ) :: id_var
    integer, dimension(:,:,:,:), optional, intent(  out) :: d
    integer, dimension(4),       optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_int_4D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(4)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_int_4D

  subroutine read_var_master_dp_0D( filename, ncid, id_var, d)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    real(dp),         intent(  out) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_var_master_dp_0D'
    character(len=1024)            :: var_name
    integer                        :: var_type
    integer                        :: ndims_of_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_dp_0D

  subroutine read_var_master_dp_1D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    integer,                          intent(in   ) :: ncid
    integer,                          intent(in   ) :: id_var
    real(dp), dimension(:), optional, intent(  out) :: d
    integer,  dimension(1), optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_dp_1D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(1)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_dp_1D

  subroutine read_var_master_dp_2D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(in   ) :: id_var
    real(dp), dimension(:,:  ), optional,intent(  out) :: d
    integer,  dimension(2),     optional,intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_dp_2D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(2)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_dp_2D

  subroutine read_var_master_dp_3D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                     intent(in   ) :: filename
    integer,                              intent(in   ) :: ncid
    integer,                              intent(in   ) :: id_var
    real(dp), dimension(:,:,:), optional, intent(  out) :: d
    integer,  dimension(3),     optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_dp_3D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(3)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_dp_3D

  subroutine read_var_master_dp_4D( filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    ! In/output variables:
    character(len=*),                       intent(in   ) :: filename
    integer,                                intent(in   ) :: ncid
    integer,                                intent(in   ) :: id_var
    real(dp), dimension(:,:,:,:), optional, intent(  out) :: d
    integer,  dimension(4),       optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_var_master_dp_4D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(4)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    end do

    ! Read the data
    if (par%master) then
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      if (nerr /= NF90_NOERR) call crash('NF90_GET_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_var_master_dp_4D

! ===== Write data to variables =====
! ===================================

  ! NOTE: only the Master actually writes data! Gathering from other processes
  !       must be done beforehand

  subroutine write_var_master_int_0D( filename, ncid, id_var, d)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    integer,          intent(in   ) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_var_master_int_0D'
    character(len=1024)            :: var_name
    integer                        :: var_type
    integer                        :: ndims_of_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_int_0D

  subroutine write_var_master_int_1D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(in   ) :: id_var
    integer,  dimension(:    ), optional,intent(in   ) :: d
    integer,  dimension(1    ), optional,intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'write_var_master_int_1D'
    character(len=1024)                                :: var_name
    integer                                            :: var_type
    integer                                            :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS)            :: dims_of_var
    integer                                            :: di
    character(len=1024)                                :: dim_name
    integer                                            :: dim_length
    integer,  dimension( 1)                            :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = 1
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if ( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

       ! Check if this dimension is large enough to read this amount of data
       if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_int_1D

  subroutine write_var_master_int_2D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    integer,                            intent(in   ) :: ncid
    integer,                            intent(in   ) :: id_var
    integer,  dimension(:,:), optional, intent(in   ) :: d
    integer,  dimension(2  ), optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_int_2D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(2)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
            trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_int_2D

  subroutine write_var_master_int_3D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(in   ) :: id_var
    integer, dimension(:,:,:), optional, intent(in   ) :: d
    integer, dimension(3),     optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_int_3D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(3)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
          trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_int_3D

  subroutine write_var_master_int_4D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                      intent(in   ) :: filename
    integer,                               intent(in   ) :: ncid
    integer,                               intent(in   ) :: id_var
    integer, dimension(:,:,:,:), optional, intent(in   ) :: d
    integer, dimension(4),       optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_int_4D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(4)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
          trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if
    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_int_4D

  subroutine write_var_master_int8_2D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                        intent(in   ) :: filename
    integer,                                 intent(in   ) :: ncid
    integer,                                 intent(in   ) :: id_var
    integer(int8), dimension(:,:), optional, intent(in   ) :: d
    integer,       dimension(2),   optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_int8_2D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(2)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_INT64)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (par%master .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
            trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_int8_2D

  subroutine write_var_master_dp_0D( filename, ncid, id_var, d)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    real(dp),         intent(in   ) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_var_master_dp_0D'
    character(len=1024)            :: var_name
    integer                        :: var_type
    integer                        :: ndims_of_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_dp_0D

  subroutine write_var_master_dp_1D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                intent(in   ) :: filename
    integer,                         intent(in   ) :: ncid
    integer,                         intent(in   ) :: id_var
    real(dp), dimension(:), optional,intent(in   ) :: d
    integer,  dimension(1), optional,intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_dp_1D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(1)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied =  1
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if ( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (var_name /= 'time') then
          ! Exception for time, because there the dimension is usually unlimited
          if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
            trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
        end if
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_dp_1D

  subroutine write_var_master_dp_2D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    integer,                            intent(in   ) :: ncid
    integer,                            intent(in   ) :: id_var
    real(dp), dimension(:,:), optional, intent(in   ) :: d
    integer,  dimension(2),   optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_dp_2D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(2)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%master) then
        if( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_dp_2D

  subroutine write_var_master_dp_3D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                     intent(in   ) :: filename
    integer,                              intent(in   ) :: ncid
    integer,                              intent(in   ) :: id_var
    real(dp), dimension(:,:,:), optional, intent(in   ) :: d
    integer,  dimension(3),     optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_dp_3D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(3)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_dp_3D

  subroutine write_var_master_dp_4D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    ! In/output variables:
    character(len=*),                       intent(in   ) :: filename
    integer,                                intent(in   ) :: ncid
    integer,                                intent(in   ) :: id_var
    real(dp), dimension(:,:,:,:), optional, intent(in   ) :: d
    integer,  dimension(4),       optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_master_dp_4D'
    character(len=1024)                     :: var_name
    integer                                 :: var_type
    integer                                 :: ndims_of_var
    integer,  dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    integer                                 :: di
    character(len=1024)                     :: dim_name
    integer                                 :: dim_length
    integer, dimension(4)                   :: start_applied, count_applied

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%master .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%master .and. ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1, 1 /)
    end if
    if (par%master .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%master .and. any( count_applied == 0)) call crash('count must be positive!')

    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do

    ! Write the data
    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_VAR failed for variable "' // trim( var_name) // '" in file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_master_dp_4D

! ===== Basic NetCDF wrapper functions =====
! ==========================================

  ! inquire dimensions and variables
  subroutine inquire_dim( filename, ncid, dim_name, dim_length, id_dim)
    ! inquire if this file contains a dimension by name of dim_name.
    ! if so, return its length and identifier; if not, return -1 for both.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: dim_name
    integer,          intent(  out) :: dim_length
    integer,          intent(  out) :: id_dim

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_dim'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    if (par%master) then

      ! Check if a dimension of this name exists in the file
      nerr = NF90_INQ_DIMID( ncid, dim_name, id_dim)

      if (nerr /= NF90_NOERR) then
        ! if a dimension by this name does not exist, return -1 for the length and ID
        id_dim     = -1
        dim_length = -1
      else
        ! if a dimension by this name exists, find its length
        nerr = NF90_inquire_dimension( ncid, id_dim, len = dim_length)
        if (nerr /= NF90_NOERR) call crash('NF90_inquire_dimension failed for file "' // trim( filename) // '"!')
      end if

    end if ! if (par%master) then
    call sync

    call MPI_BCAST( id_dim    , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( dim_length, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_dim

  subroutine inquire_dim_info( filename, ncid, id_dim, dim_name, dim_length)
    ! inquire some info of a dimension

    ! In/output variables:
    character(len=*),              intent(in   ) :: filename
    integer,                       intent(in   ) :: ncid
    integer,                       intent(in   ) :: id_dim
    character(len=1024), optional, intent(  out) :: dim_name
    integer,             optional, intent(  out) :: dim_length

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_dim_info'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    if (par%master) then
      ! inquire some info on this variable
      nerr = NF90_inquire_dimension( ncid, id_dim, name = dim_name, len = dim_length)
      if (nerr /= NF90_NOERR) call crash('NF90_inquire_dimension failed for file "' // trim( filename) // '"!')
    end if ! if (par%master) then

    if (present( dim_name  )) call MPI_BCAST( dim_name  , 256, MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    if (present( dim_length)) call MPI_BCAST( dim_length, 1  , MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_dim_info

  subroutine inquire_var( filename, ncid, var_name, id_var)
    ! inquire if this file contains a variable by name of var_name.
    ! if so, return its identifier. if not, return -1 for the identifier.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: var_name
    integer,          intent(  out) :: id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_var'
    character                      :: dummy1

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! To prevent "unused variable" compiler warnings
    dummy1 = filename( 1:1)

    if (par%master) then

      ! Check if a variable of this name exists in the file
      nerr = NF90_INQ_VARID( ncid, var_name, id_var)

      if (nerr /= NF90_NOERR) then
        ! if a variable by this name does not exist, return -1 for the ID
        id_var = -1
      end if

    end if ! if (par%master) then
    call sync

    call MPI_BCAST( id_var, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_var

  subroutine inquire_var_info( filename, ncid, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    ! inquire some info of a variable

    ! In/output variables:
    character(len=*),                                 intent(in   ) :: filename
    integer,                                          intent(in   ) :: ncid
    integer,                                          intent(in   ) :: id_var
    character(len=1024),                    optional, intent(  out) :: var_name
    integer,                                optional, intent(  out) :: var_type
    integer,                                optional, intent(  out) :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS), optional, intent(  out) :: dims_of_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_var_info'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    if (par%master) then
      ! inquire some info on this variable
      nerr = NF90_inquire_VARIABLE( ncid, id_var, name = var_name, xtype = var_type, ndims = ndims_of_var, dimids = dims_of_var)
      if (nerr /= NF90_NOERR) call crash('NF90_inquire_VARIABLE failed for file "' // trim( filename) // '"!')
    end if ! if (par%master) then

    if (present( var_name    )) call MPI_BCAST( var_name    , 256              , MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    if (present( var_type    )) call MPI_BCAST( var_type    , 1                , MPI_integer, 0, MPI_COMM_WORLD, ierr)
    if (present( ndims_of_var)) call MPI_BCAST( ndims_of_var, 1                , MPI_integer, 0, MPI_COMM_WORLD, ierr)
    if (present(  dims_of_var)) call MPI_BCAST( dims_of_var , NF90_MAX_VAR_DIMS, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_var_info

  ! Create new NetCDF file
  subroutine create_new_netcdf_file_for_writing( filename, ncid)
    ! Create a new NetCDF file in the specified location for writing.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(  out) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_new_netcdf_file_for_writing'
    logical                        :: file_exists

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if this file already exists
    if (par%master) then
      inquire( exist = file_exists, file = trim( filename))
      if (file_exists) then
        call crash('file "' // trim( filename) // '" already exists!')
      end if
    end if

    ! Create the NetCDF file
    if (par%master) then
      nerr = NF90_CREATE( filename, ior( NF90_NOCLOBBER, NF90_NETCDF4), ncid)
      if (nerr /= NF90_NOERR) call crash('NF90_CREATE failed for file "' // trim( filename) // '"!')
    end if ! if (par%master) then
    call MPI_BCAST( ncid, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Add some very basic info about the current simulation to the header
    call add_attribute_char( filename, ncid, NF90_GLOBAL, 'git commit hash', git_commit_hash)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_new_netcdf_file_for_writing

  ! Create dimensions, variables, and attributes
  subroutine create_dimension( filename, ncid, dim_name, dim_length, id_dim)
    ! Create a new dimension in a NetCDF file.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: dim_name
    integer,          intent(in   ) :: dim_length
    integer,          intent(  out) :: id_dim

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_dimension'
    integer                        :: dim_length_present

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Safety: check if a dimension by this name is already present in this file
    call inquire_dim( filename, ncid, dim_name, dim_length_present, id_dim)
    if (id_dim /= -1) then
      !call crash('file "' // trim( filename) // '" already contains dimension "' // trim( dim_name) // '"!')
      call finalise_routine( routine_name)
      RETURN
    end if

    ! Add the dimension
    if (par%master) then
      nerr = NF90_DEF_DIM( ncid, dim_name, dim_length, id_dim)
      if (nerr /= NF90_NOERR) call crash('NF90_DEF_DIM failed for file "' // trim( filename) // '"!')
    end if
    call sync

    call MPI_BCAST( id_dim, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_dimension

  subroutine create_variable( filename, ncid, var_name, var_type, dim_ids, id_var)
    ! Create a new variable in a NetCDF file.

    ! In/output variables:
    character(len=*),      intent(in   ) :: filename
    integer,               intent(in   ) :: ncid
    character(len=*),      intent(in   ) :: var_name
    integer,               intent(in   ) :: var_type
    integer, dimension(:), intent(in   ) :: dim_ids
    integer,               intent(  out) :: id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_variable'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Safety: check if a variable by this name is already present in this file
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var /= -1) then
      !call crash('file "' // trim( filename) // '" already contains variable "' // trim( var_name) // '"!')
      call finalise_routine( routine_name)
      return
    end if

    ! Add the variable
    if (par%master) then
      nerr = NF90_DEF_VAR( ncid, name = var_name, xtype = var_type, dimids = dim_ids, varid = id_var)
      if (nerr /= NF90_NOERR) call crash('NF90_DEF_VAR failed for file "' // trim( filename) // '"!')
    end if
    call sync

    call MPI_BCAST( id_var, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_variable

  subroutine create_scalar_variable( filename, ncid, var_name, var_type, id_var)
    ! Create a new scalar variable in a NetCDF file.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: var_name
    integer,          intent(in   ) :: var_type
    integer,          intent(  out) :: id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_scalar_variable'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Safety: check if a variable by this name is already present in this file
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var /= -1) call crash('file "' // trim( filename) // '" already contains variable "' // trim( var_name) // '"!')

    ! Add the variable
    if (par%master) then
      nerr = NF90_DEF_VAR( ncid, name = var_name, xtype = var_type, varid = id_var)
      if (nerr /= NF90_NOERR) call crash('NF90_DEF_VAR failed for file "' // trim( filename) // '"!')
    end if
    call sync

    call MPI_BCAST( id_var, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_scalar_variable

  subroutine add_attribute_int( filename, ncid, id_var, att_name, att_val)
    ! Add an integer-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    character(len=*), intent(in   ) :: att_name
    integer,          intent(in   ) :: att_val

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_attribute_int'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Add the attribute
    if (par%master) then
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_ATT failed for file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_attribute_int

  subroutine add_attribute_dp( filename, ncid, id_var, att_name, att_val)
    ! Add a double-precision-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    character(len=*), intent(in   ) :: att_name
    real(dp),         intent(in   ) :: att_val

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_attribute_dp'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Add the attribute
    if (par%master) then
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_ATT failed for file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_attribute_dp

  subroutine add_attribute_char( filename, ncid, id_var, att_name, att_val)
    ! Add a character-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    character(len=*), intent(in   ) :: att_name
    character(len=*), intent(in   ) :: att_val

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_attribute_char'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Add the attribute
    if (par%master) then
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      if (nerr /= NF90_NOERR) call crash('NF90_PUT_ATT failed for file "' // trim( filename) // '"!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_attribute_char

  ! Open and close a NetCDF file
  subroutine open_existing_netcdf_file_for_reading( filename, ncid)
    ! Open the NetCDF file in the specified location for reading only,
    ! and return its identifier.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(  out) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'open_netcdf_file_for_reading'
    logical                        :: file_exists

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) then
      call crash('file "' // trim( filename) // '" not found!')
    end if

    ! Open the NetCDF file with read-only access
    if (par%master) then
      nerr = NF90_OPEN( trim( filename), NF90_NOWRITE, ncid)
      if (nerr /= NF90_NOERR) call crash('NF90_OPEN failed for file "' // trim( filename) // '"!')
    end if ! if (par%master) then

    call MPI_BCAST( ncid, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine open_existing_netcdf_file_for_reading

  subroutine open_existing_netcdf_file_for_writing( filename, ncid)
    ! Open an existing NetCDF file in data mode
    ! In data mode, no new dimensions, variables, or attributes can be created,
    ! but data can be written to existing variables.
    ! When opening an existing NetCDF file, it is by default in data mode.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(  out) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'open_existing_netcdf_file_for_writing'
    logical                        :: file_exists

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) then
      call crash('file "' // trim( filename) // '" not found!')
    end if

    ! Open the NetCDF file with read+write access
    if (par%master) then
      nerr = NF90_OPEN( trim( filename), ior( NF90_WRITE, NF90_SHARE), ncid)
      if (nerr /= NF90_NOERR) call crash('NF90_OPEN failed for file "' // trim( filename) // '" beeperdebeep!')
    end if ! if (par%master) then

    call MPI_BCAST( ncid, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine open_existing_netcdf_file_for_writing

  subroutine close_netcdf_file( ncid)
    ! Close an opened NetCDF file

    ! In/output variables:
    integer, intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'close_netcdf_file'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Close netCDF file:
    if (par%master) then
      nerr = NF90_CLOSE( ncid)
      if (nerr /= NF90_NOERR) call crash('NF90_CLOSE failed!')
    end if ! if (par%master) then
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine close_netcdf_file

end module netcdf_basic
