MODULE netcdf_basic

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

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C, git_commit_hash
  USE math_utilities                                         , ONLY: check_for_NaN_dp_0D, check_for_NaN_int_0D, &
                                                                     check_for_NaN_dp_1D, check_for_NaN_int_1D, &
                                                                     check_for_NaN_dp_2D, check_for_NaN_int_2D, &
                                                                     check_for_NaN_dp_3D, check_for_NaN_int_3D, &
                                                                     check_for_NaN_dp_4D, check_for_NaN_int_4D

  ! Import  NetCDF functionality
  USE netcdf, ONLY: NF90_NOERR, NF90_OPEN, NF90_CLOSE, NF90_NOWRITE, NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
              NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, NF90_MAX_VAR_DIMS, NF90_GET_VAR, &
              NF90_CREATE, NF90_NOCLOBBER, NF90_NETCDF4, NF90_ENDDEF, NF90_REDEF, NF90_DEF_DIM, NF90_DEF_VAR, &
              NF90_PUT_ATT, NF90_WRITE, NF90_INT, NF90_FLOAT, NF90_DOUBLE, NF90_PUT_VAR, NF90_UNLIMITED, &
              NF90_INQUIRE_ATTRIBUTE, NF90_SHARE, NF90_GLOBAL

  IMPLICIT NONE

  ! NetCDF error code
  INTEGER :: nerr

  ! Possible names for different dimensions and variables
  ! =====================================================

  ! Different options for the name of a dimension or variable can now be tried.
  ! They are separated by a double vertical bar ||

  ! Dimensions
  CHARACTER(LEN=256), PARAMETER :: field_name_options_x              = 'x||X||x1||X1||nx||NX||x-coordinate||X-coordinate||easting||Easting'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_y              = 'y||Y||y1||Y1||ny||NY||y-coordinate||Y-coordinate||northing||Northing'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_zeta           = 'zeta||Zeta'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_lon            = 'lon||Lon||long||Long||longitude||Longitude'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_lat            = 'lat||Lat||latitude||Latitude'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_time           = 'time||Time||t||nt'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_month          = 'month||Month'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_depth          = 'depth||Depth'

  ! Mesh
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dim_nV         = 'vi'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dim_nTri       = 'ti'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dim_nC_mem     = 'ci'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dim_nE         = 'ei'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dim_two        = 'two'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dim_three      = 'three'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dim_four       = 'four'

  CHARACTER(LEN=256), PARAMETER :: field_name_options_V              = 'V'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Tri            = 'Tri'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_nC             = 'nC'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_C              = 'C'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_niTri          = 'niTri'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_iTri           = 'iTri'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_VBI            = 'VBI'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Tricc          = 'Tricc'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_TriC           = 'TriC'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_TriBI          = 'TriBI'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_E              = 'E'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_VE             = 'VE'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_EV             = 'EV'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_ETri           = 'ETri'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_EBI            = 'EBI'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_TriGC          = 'TriGC'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_A              = 'A'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_R              = 'R'

  ! Variables
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Hi             = 'Hi||thickness||lithk'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Hb             = 'Hb||bed||topg'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Hs             = 'Hs||surface||orog'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_SL             = 'SL'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dHb            = 'dHb'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Ti             = 'Ti'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_T_ocean        = 'T_ocean||t_ocean||t_an'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_S_ocean        = 'S_ocean||s_ocean||s_an'

CONTAINS

  ! Check if a file contains all the variables and dimensions describing
  ! an x/y-grid, a lon/lat-grid, or a mesh
  SUBROUTINE inquire_xy_grid( filename, has_xy_grid)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular x/y-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_xy_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_xy_grid'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var_x, id_var_y

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for x and y dimensions and variables
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)
    CALL inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    CALL inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Check if everything is there
    has_xy_grid = .TRUE.

    IF (id_dim_x              == -1) has_xy_grid = .FALSE.
    IF (id_dim_y              == -1) has_xy_grid = .FALSE.
    IF (id_var_x              == -1) has_xy_grid = .FALSE.
    IF (id_var_y              == -1) has_xy_grid = .FALSE.

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_xy_grid

  SUBROUTINE inquire_lonlat_grid( filename, has_lonlat_grid)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular lon/lat-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_lonlat_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_lonlat_grid'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: id_var_lon, id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for x and y dimensions and variables
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)
    CALL inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    CALL inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Check if everything is there
    has_lonlat_grid = .TRUE.

    IF (id_dim_lon            == -1) has_lonlat_grid = .FALSE.
    IF (id_dim_lat            == -1) has_lonlat_grid = .FALSE.
    IF (id_var_lon            == -1) has_lonlat_grid = .FALSE.
    IF (id_var_lat            == -1) has_lonlat_grid = .FALSE.

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_lonlat_grid

  SUBROUTINE inquire_mesh( filename, has_mesh)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a mesh.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_mesh'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    INTEGER                                            :: id_var_V, id_var_nC, id_var_C, id_var_niTri, id_var_iTri, id_var_VBI
    INTEGER                                            :: id_var_Tri, id_var_Tricc, id_var_TriC, id_var_TrIBI

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three)

    ! Inquire mesh variables
    CALL inquire_var_multopt( filename, ncid, field_name_options_V         , id_var_V    )
    CALL inquire_var_multopt( filename, ncid, field_name_options_nC        , id_var_nC   )
    CALL inquire_var_multopt( filename, ncid, field_name_options_C         , id_var_C    )
    CALL inquire_var_multopt( filename, ncid, field_name_options_niTri     , id_var_niTri)
    CALL inquire_var_multopt( filename, ncid, field_name_options_iTri      , id_var_iTri )
    CALL inquire_var_multopt( filename, ncid, field_name_options_VBI       , id_var_VBI  )
    CALL inquire_var_multopt( filename, ncid, field_name_options_Tri       , id_var_Tri  )
    CALL inquire_var_multopt( filename, ncid, field_name_options_Tricc     , id_var_Tricc)
    CALL inquire_var_multopt( filename, ncid, field_name_options_TriC      , id_var_TriC )
    CALL inquire_var_multopt( filename, ncid, field_name_options_TriBI     , id_var_TriBI)

    ! Check if everything is there
    has_mesh = .TRUE.

    IF (id_dim_vi    == -1) has_mesh = .FALSE.
    IF (id_dim_ti    == -1) has_mesh = .FALSE.
    IF (id_dim_ci    == -1) has_mesh = .FALSE.
    IF (id_dim_two   == -1) has_mesh = .FALSE.
    IF (id_dim_three == -1) has_mesh = .FALSE.

    IF (id_var_V     == -1) has_mesh = .FALSE.
    IF (id_var_nC    == -1) has_mesh = .FALSE.
    IF (id_var_C     == -1) has_mesh = .FALSE.
    IF (id_var_niTri == -1) has_mesh = .FALSE.
    IF (id_var_iTri  == -1) has_mesh = .FALSE.
    IF (id_var_VBI   == -1) has_mesh = .FALSE.
    IF (id_var_Tri   == -1) has_mesh = .FALSE.
    IF (id_var_Tricc == -1) has_mesh = .FALSE.
    IF (id_var_TriC  == -1) has_mesh = .FALSE.
    IF (id_var_TriBI == -1) has_mesh = .FALSE.

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_mesh

  ! Inquire if a file contains the variable and dimension for
  ! zeta, z_ocean, time, or months
  SUBROUTINE inquire_zeta( filename, ncid, has_zeta)
    ! Inquire if a NetCDF file contains a zeta dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    LOGICAL,                             INTENT(OUT)   :: has_zeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_zeta'
    INTEGER                                            :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for zeta dimension and variable
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)
    CALL inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var_zeta)

    ! Check if everything is there
    has_zeta = .TRUE.

    IF (id_dim_zeta           == -1) has_zeta = .FALSE.
    IF (id_var_zeta           == -1) has_zeta = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_zeta

  SUBROUTINE inquire_month( filename, ncid, has_month)
    ! Inquire if a NetCDF file contains a month dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    LOGICAL,                             INTENT(OUT)   :: has_month

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_month'
    INTEGER                                            :: id_dim_month, id_var_month

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for month dimension and variable
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)
    CALL inquire_var_multopt( filename, ncid, field_name_options_month, id_var_month)

    ! Check if everything is there
    has_month = .TRUE.

    IF (id_dim_month           == -1) has_month = .FALSE.
    IF (id_var_month           == -1) has_month = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_month

  SUBROUTINE inquire_time( filename, ncid, has_time)
    ! Inquire if a NetCDF file contains a time dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    LOGICAL,                             INTENT(OUT)   :: has_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_time'
    INTEGER                                            :: id_dim_time, id_var_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for time dimension and variable
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
    CALL inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! Check if everything is there
    has_time = .TRUE.

    IF (id_dim_time           == -1) has_time = .FALSE.
    IF (id_var_time           == -1) has_time = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_time

  SUBROUTINE find_timeframe( filename, ncid, time, ti)
    ! Find the timeframe in the file that is closest to the desired time.
    ! If the file has no time dimension or variable, throw an error.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    REAL(dp),                            INTENT(IN)    :: time
    INTEGER,                             INTENT(OUT)   :: ti

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_timeframe'
    INTEGER                                            :: nt, id_dim_time, id_var_time
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: time_from_file
    INTEGER                                            :: tii
    REAL(dp)                                           :: dt_min

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file contains a valid time dimension and variable
    CALL check_time( filename, ncid)

    ! Inquire size of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! Inquire time variable ID
    CALL inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! Allocate memory
    ALLOCATE( time_from_file( nt))

    ! Read time from file
    CALL read_var_master_dp_1D( filename, ncid, id_var_time, time_from_file)
    CALL MPI_BCAST( time_from_file, nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Find timeframe closest to desired time
    IF (time_from_file( 1) > time) THEN
      ! Desired time beyond lower limit
      CALL warning('desired timeframe at t = {dp_01} before start of file time for file "' // TRIM( filename) // '"; reading data from t = {dp_02} instead!', dp_01 = time, dp_02 = time_from_file( 1))
      ti = 1
    ELSEIF (time_from_file( nt) < time) THEN
      ! Desired time beyond upper limit
      CALL warning('desired timeframe at t = {dp_01} after end of file time for file "' // TRIM( filename) // '"; reading data from t = {dp_02} instead!', dp_01 = time, dp_02 = time_from_file( nt))
      ti = nt
    ELSE
      ! Desired time is within the file time
      dt_min = HUGE( 1._dp)
      DO tii = 1, nt
        IF (ABS( time_from_file( tii) - time) < dt_min) THEN
          ti = tii
          dt_min = ABS( time_from_file( tii) - time)
        END IF
      END DO
      IF (dt_min > 0._dp) THEN
        CALL warning('desired timeframe at t = {dp_01} not present in file "' // TRIM( filename) // '"; reading data from closest match at t = {dp_02} instead!', dp_01 = time, dp_02 = time_from_file( ti))
      END IF
    END IF

    ! Clean up after yourself
    DEALLOCATE( time_from_file)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_timeframe

! ===== Safety checks on variables and dimensions =====
! =====================================================

  ! x/y-grid dimensions
  SUBROUTINE check_x( filename, ncid)
    ! Check if this file contains a valid x dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_x'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x
    REAL(dp)                                           :: dx, dxp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" has length {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_x, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid x variable could be found in file "' // TRIM( filename) // '"!')

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('variable "' // TRIM( var_name) // &
      '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check variable dimension
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ' // TRIM( dim_name) // ' as a dimension!')

    ! Allocate memory
    ALLOCATE( x( n))

    ! Read variable
    CALL read_var_master_dp_1D( filename, ncid, id_var, x)

    ! Check validity
    IF (par%master) CALL check_for_NaN_dp_1D( x, 'x')

    ! Check grid spacing
    IF (par%master) THEN
      dx = x( 2) - x( 1)
      DO i = 2, n
        dxp = x( i) - x( i-1)
        IF (ABS( 1._dp - dxp / dx) > 1E-5_dp) CALL crash('x coordinate in file "' // TRIM( filename) // '" is irregular!')
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( x)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_x

  SUBROUTINE check_y( filename, ncid)
    ! Check if this file contains a valid y dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_y'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: y
    REAL(dp)                                           :: dy, dyp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" has length {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_y, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid y variable could be found in file "' // TRIM( filename) // '"!')

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('variable "' // TRIM( var_name) // &
      '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check variable dimension
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ' // TRIM( dim_name) // ' as a dimension!')

    ! Allocate memory
    ALLOCATE( y( n))

    ! Read variable
    CALL read_var_master_dp_1D( filename, ncid, id_var, y)

    ! Check validity
    IF (par%master) CALL check_for_NaN_dp_1D( y, 'y')

    ! Check grid spacing
    IF (par%master) THEN
      dy = y( 2) - y( 1)
      DO i = 2, n
        dyp = y( i) - y( i-1)
        IF (ABS( 1._dp - dyp / dy) > 1E-5_dp) CALL crash('y coordinate in file "' // TRIM( filename) // '" is irregular!')
      END DO
    END IF
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( y)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_y

  ! lon/lat-grid dimensions
  SUBROUTINE check_lon( filename, ncid)
    ! Check if this file contains a valid longitude dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lon'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: lon
    REAL(dp)                                           :: dlon, dlonp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid longitude dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('longitude dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('longitude dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_lon, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid longitude variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('longitude variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('longitude variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('longitude variable in file "' // TRIM( filename) // '" does not have longitude as a dimension!')

    ! Allocate memory
    ALLOCATE( lon( n))

    ! Read variable
    CALL read_var_master_dp_1D( filename, ncid, id_var, lon)

    ! Check validity
    IF (par%master) CALL check_for_NaN_dp_1D( lon, 'lon')

    ! Check grid spacing
    IF (par%master) THEN
      dlon = lon( 2) - lon( 1)
      DO i = 2, n
        dlonp = lon( i) - lon( i-1)
        IF (ABS( 1._dp - dlonp / dlon) > 1E-5_dp) CALL crash('longitude coordinate in file "' // TRIM( filename) // '" is irregular!')
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( lon)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lon

  SUBROUTINE check_lat( filename, ncid)
    ! Check if this file contains a valid latitude dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lat'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: lat
    REAL(dp)                                           :: dlat, dlatp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid latitude dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('latitude dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('latitude dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_lat, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid latitude variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('latitude variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('latitude variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('latitude variable in file "' // TRIM( filename) // '" does not have latitude as a dimension!')

    ! Allocate memory
    ALLOCATE( lat( n))

    ! Read variable
    CALL read_var_master_dp_1D( filename, ncid, id_var, lat)

    ! Check validity
    IF (par%master) CALL check_for_NaN_dp_1D( lat, 'lat')

    ! Check grid spacing
    IF (par%master) THEN
      dlat = lat( 2) - lat( 1)
      DO i = 2, n
        dlatp = lat( i) - lat( i-1)
        IF (ABS( 1._dp - dlatp / dlat) > 1E-5_dp) CALL crash('latitude coordinate in file "' // TRIM( filename) // '" is irregular!')
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( lat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lat

  ! Mesh dimensions
  SUBROUTINE check_mesh_dimensions( filename, ncid)
    ! Check if this file contains valid mesh dimensions and variables

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_dimensions'
    INTEGER                                            :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    INTEGER                                            :: nV, nTri, nC_mem, n_two, n_three
    CHARACTER(LEN=256)                                 :: dim_name_vi, dim_name_ti, dim_name_ci, dim_name_two, dim_name_three
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

  ! == Inquire dimensions
  ! =====================

    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   , dim_length = nV     , dim_name = dim_name_vi   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   , dim_length = nTri   , dim_name = dim_name_ti   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   , dim_length = nC_mem , dim_name = dim_name_ci   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  , dim_length = n_two  , dim_name = dim_name_two  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three, dim_length = n_three, dim_name = dim_name_three)

    ! Safety checks on dimensions
    IF (id_dim_vi    == -1) CALL crash('no valid vi    dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_ti    == -1) CALL crash('no valid ti    dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_ci    == -1) CALL crash('no valid ci    dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_two   == -1) CALL crash('no valid two   dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_three == -1) CALL crash('no valid three dimension could be found in file "' // TRIM( filename) // '"!')

    IF (nV      == NF90_UNLIMITED) CALL crash('vi    dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (nTri    == NF90_UNLIMITED) CALL crash('ti    dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (nC_mem  == NF90_UNLIMITED) CALL crash('ci    dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n_two   == NF90_UNLIMITED) CALL crash('two   dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n_three == NF90_UNLIMITED) CALL crash('three dimension in file "' // TRIM( filename) // '" is unlimited!')

    IF (nV      <  1) CALL crash('vi    dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = nV     )
    IF (nTri    <  1) CALL crash('ti    dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = nTri   )
    IF (nC_mem  <  1) CALL crash('ci    dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = nC_mem )
    IF (n_two   /= 2) CALL crash('two   dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n_two  )
    IF (n_three /= 3) CALL crash('three dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n_three)

  ! == Inquire variables
  ! ====================

    ! Metadata
    ! ========

    ! xmin
    CALL inquire_var_multopt( filename, ncid, 'xmin', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid xmin variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! xmax
    CALL inquire_var_multopt( filename, ncid, 'xmax', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid xmax variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! ymin
    CALL inquire_var_multopt( filename, ncid, 'ymin', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid ymin variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! ymax
    CALL inquire_var_multopt( filename, ncid, 'ymax', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid ymax variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! tol_dist
    CALL inquire_var_multopt( filename, ncid, 'tol_dist', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid tol_dist variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! lambda_M
    CALL inquire_var_multopt( filename, ncid, 'lambda_M', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid lambda_M variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! phi_M
    CALL inquire_var_multopt( filename, ncid, 'phi_M', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid phi_M variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! beta_stereo
    CALL inquire_var_multopt( filename, ncid, 'beta_stereo', id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var)
    IF (id_var == -1) CALL crash('no valid beta_stereo variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Vertex data
    ! ===========

    ! V
    CALL inquire_var_multopt( filename, ncid, field_name_options_V, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid V variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. (dims_of_var( 1) == id_dim_vi .AND. dims_of_var( 2) == id_dim_two)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi and two as dimensions!')

    ! nC
    CALL inquire_var_multopt( filename, ncid, field_name_options_nC, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid nC variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. var_type == NF90_INT) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. dims_of_var( 1) == id_dim_vi) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')

    ! C
    CALL inquire_var_multopt( filename, ncid, field_name_options_C, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid C variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. var_type == NF90_INT) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. (dims_of_var( 1) == id_dim_vi .AND. dims_of_var( 2) == id_dim_ci)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi and ci as dimensions!')

    ! niTri
    CALL inquire_var_multopt( filename, ncid, field_name_options_niTri, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid niTri variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. var_type == NF90_INT) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. dims_of_var( 1) == id_dim_vi) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')

    ! iTri
    CALL inquire_var_multopt( filename, ncid, field_name_options_iTri, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid iTri variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. var_type == NF90_INT) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. (dims_of_var( 1) == id_dim_vi .AND. dims_of_var( 2) == id_dim_ci)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi and ci as dimensions!')

    ! VBI
    CALL inquire_var_multopt( filename, ncid, field_name_options_VBI, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid VBI variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. var_type == NF90_INT) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. dims_of_var( 1) == id_dim_vi) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')

    ! Triangle data
    ! =============

    ! Tri
    CALL inquire_var_multopt( filename, ncid, field_name_options_Tri, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid Tri variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. var_type == NF90_INT) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. (dims_of_var( 1) == id_dim_ti .AND. dims_of_var( 2) == id_dim_three)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ti and three as dimensions!')

    ! Tricc
    CALL inquire_var_multopt( filename, ncid, field_name_options_Tricc, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid Tricc variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. (dims_of_var( 1) == id_dim_ti .AND. dims_of_var( 2) == id_dim_two)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ti and two as dimensions!')

    ! TriC
    CALL inquire_var_multopt( filename, ncid, field_name_options_TriC, id_var, var_name = var_name, &
      var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid TriC variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. var_type == NF90_INT) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. (dims_of_var( 1) == id_dim_ti .AND. dims_of_var( 2) == id_dim_three)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ti and three as dimensions!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_dimensions

  ! Zeta, z_ocean, month, time dimensions
  SUBROUTINE check_zeta( filename, ncid)
    ! Check if this file contains a valid zeta dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_zeta'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: zeta
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('zeta dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('zeta dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid zeta variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('zeta variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('zeta variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('zeta variable in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    ! Allocate memory
    ALLOCATE( zeta( n))

    ! Read variable
    CALL read_var_master_dp_1D( filename, ncid, id_var, zeta)

    ! Check validity
    IF (par%master) THEN
      CALL check_for_NaN_dp_1D( zeta, 'zeta')

      IF (zeta( 1) /= 0._dp) CALL crash('zeta in file "' // TRIM( filename) // '" does not start at zero!')
      IF (zeta( n) /= 1._dp) CALL crash('zeta in file "' // TRIM( filename) // '" does not end at one!')

      DO k = 2, n
        IF (zeta( k) <= zeta( k-1)) CALL crash('zeta in file "' // TRIM( filename) // '" does not increase monotonously!')
      END DO
    END IF
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( zeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_zeta

  SUBROUTINE check_month( filename, ncid)
    ! Check if this file contains a valid month dimension (we don't really care about the variable)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_month'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('month dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n /= 12) CALL crash('month dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_month

  SUBROUTINE check_time( filename, ncid)
    ! Check if this file contains a valid time dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_time'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid time dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n < 0) CALL crash('time dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_time, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid time variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('time variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('time variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('time variable in file "' // TRIM( filename) // '" does not have time as a dimension!')

    ! For new output files, time is still empty. If it's not, check if entries are valid
    IF (n > 0) THEN

      ! Allocate memory
      ALLOCATE( time( n))

      ! Read variable
      CALL read_var_master_dp_1D( filename, ncid, id_var, time)

      ! Check validity
      IF (par%master) CALL check_for_NaN_dp_1D( time, 'time')

      ! Clean up after yourself
      DEALLOCATE( time)

    END IF ! IF (n > 0) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_time

  SUBROUTINE check_depth( filename, ncid)
    ! Check if this file contains a valid depth dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_depth'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: depth
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid depth dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('depth dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('depth dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multopt( filename, ncid, field_name_options_depth, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid depth variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('depth variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('depth variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('depth variable in file "' // TRIM( filename) // '" does not have depth as a dimension!')

    ! Allocate memory
    ALLOCATE( depth( n))

    ! Read variable
    CALL read_var_master_dp_1D( filename, ncid, id_var, depth)

    ! Check validity
    IF (par%master) THEN
      CALL check_for_NaN_dp_1D( depth, 'depth')

      DO k = 2, n
        IF (depth( k) <= depth( k-1)) CALL crash('depth in file "' // TRIM( filename) // '" does not increase monotonously!')
      END DO
    END IF
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( depth)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_depth

  ! x/y-grid field variables
  SUBROUTINE check_xy_grid_field_int_2D(            filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_int_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Inquire x,y dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. var_type == NF90_INT) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has x,y as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_int_2D

  SUBROUTINE check_xy_grid_field_dp_2D(             filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_dp_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Inquire x,y dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has x,y as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_dp_2D

  SUBROUTINE check_xy_grid_field_dp_2D_monthly(     filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_dp_2D_monthly'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x(     filename, ncid)
    CALL check_y(     filename, ncid)
    CALL check_month( filename, ncid)

    ! Inquire x,y dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x    )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y    )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_month)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have month as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has x,y,m as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y,m as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y,m as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_dp_2D_monthly

  SUBROUTINE check_xy_grid_field_dp_3D(             filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_dp_3D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x(    filename, ncid)
    CALL check_y(    filename, ncid)
    CALL check_zeta( filename, ncid)

    ! Inquire x,y,zeta dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x   )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y   )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_zeta)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has x,y,zeta as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y,zeta as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y,zeta as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_dp_3D

  SUBROUTINE check_xy_grid_field_dp_3D_ocean(       filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_dp_3D_ocean'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_depth, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x(     filename, ncid)
    CALL check_y(     filename, ncid)
    CALL check_depth( filename, ncid)

    ! Inquire x,y,depth dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x    )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y    )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_depth)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have depth as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has x,y,depth as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y,depth as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y,depth as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_dp_3D_ocean

  ! lon/lat-grid field variables
  SUBROUTINE check_lonlat_grid_field_int_2D(        filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_int_2D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon( filename, ncid)
    CALL check_lat( filename, ncid)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. var_type == NF90_INT) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has lon,lat as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_int_2D

  SUBROUTINE check_lonlat_grid_field_dp_2D(         filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_2D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon( filename, ncid)
    CALL check_lat( filename, ncid)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has lon,lat as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_2D

  SUBROUTINE check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_2D_monthly'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_month, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon(   filename, ncid)
    CALL check_lat(   filename, ncid)
    CALL check_month( filename, ncid)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_month)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have month as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has lon,lat,m as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat,m as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat,m as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_2D_monthly

  SUBROUTINE check_lonlat_grid_field_dp_3D(         filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_3D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_zeta, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon(  filename, ncid)
    CALL check_lat(  filename, ncid)
    CALL check_zeta( filename, ncid)

    ! Inquire lon,lat,zeta dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_zeta)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has lon,lat,zeta as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat,zeta as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat,zeta as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_3D

  SUBROUTINE check_lonlat_grid_field_dp_3D_ocean(   filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_3D_ocean'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_depth, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon(   filename, ncid)
    CALL check_lat(   filename, ncid)
    CALL check_depth( filename, ncid)

    ! Inquire lon,lat,depth dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_depth)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have depth as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has lon,lat,depth as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat,depth as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat,depth as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_3D_ocean

  ! mesh field variables
  SUBROUTINE check_mesh_field_int_2D(               filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_int_2D'
    INTEGER                                            :: id_dim_vi, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. var_type == NF90_INT) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_vi)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 1) THEN
          ! The variable only has vi as a dimension
        ELSE
          IF (ndims_of_var == 2) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has two dimensions, but the second one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_int_2D

  SUBROUTINE check_mesh_field_int_2D_b(             filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_int_2D_b'
    INTEGER                                            :: id_dim_ti, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. var_type == NF90_INT) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_ti)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ti as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 1) THEN
          ! The variable only has vi as a dimension
        ELSE
          IF (ndims_of_var == 2) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has two dimensions, but the second one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_int_2D_b

  SUBROUTINE check_mesh_field_int_2D_c(             filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_int_2D_c'
    INTEGER                                            :: id_dim_ei, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. var_type == NF90_INT) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_ei)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ei as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 1) THEN
          ! The variable only has vi as a dimension
        ELSE
          IF (ndims_of_var == 2) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has two dimensions, but the second one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_int_2D_c

  SUBROUTINE check_mesh_field_dp_2D(                filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_dp_2D'
    INTEGER                                            :: id_dim_vi, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_vi)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 1) THEN
          ! The variable only has vi as a dimension
        ELSE
          IF (ndims_of_var == 2) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has two dimensions, but the second one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_dp_2D

  SUBROUTINE check_mesh_field_dp_2D_b(              filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_dp_2D_b'
    INTEGER                                            :: id_dim_ti, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_ti)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ti as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 1) THEN
          ! The variable only has vi as a dimension
        ELSE
          IF (ndims_of_var == 2) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has two dimensions, but the second one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_dp_2D_b

  SUBROUTINE check_mesh_field_dp_2D_c(              filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_dp_2D_c'
    INTEGER                                            :: id_dim_ei, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_ei)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ei as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 1) THEN
          ! The variable only has vi as a dimension
        ELSE
          IF (ndims_of_var == 2) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has two dimensions, but the second one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_dp_2D_c

  SUBROUTINE check_mesh_field_dp_2D_monthly(        filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly mesh variable by this name
    !
    ! NOTE: this is 2-D monthly in the physical sense, so a 2-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_dp_2D_monthly'
    INTEGER                                            :: id_dim_vi, id_dim_month, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)
    CALL check_month(           filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month , id_dim_month)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_vi   )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_month)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have month as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has vi,m as dimensions
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi,m as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi,m as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_dp_2D_monthly

  SUBROUTINE check_mesh_field_dp_3D(                filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_dp_3D'
    INTEGER                                            :: id_dim_vi, id_dim_zeta, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)
    CALL check_zeta(            filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta  , id_dim_zeta)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_vi  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_zeta)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has vi,zeta as dimensions
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi,zeta as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi,zeta as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_dp_3D

  SUBROUTINE check_mesh_field_dp_3D_b(              filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_dp_3D'
    INTEGER                                            :: id_dim_ti, id_dim_zeta, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)
    CALL check_zeta(            filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta    , id_dim_zeta)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_ti  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ti as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_zeta)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has ti,zeta as dimensions
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have ti,zeta as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have ti,zeta as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_dp_3D_b

  SUBROUTINE check_mesh_field_dp_3D_ocean(          filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_mesh_field_dp_3D_ocean'
    INTEGER                                            :: id_dim_vi, id_dim_depth, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid mesh dimensions and variables
    CALL check_mesh_dimensions( filename, ncid)
    CALL check_depth(           filename, ncid)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_depth , id_dim_depth)

    ! Inquire variable
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check mesh dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_vi   )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have vi as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_depth)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have depth as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has vi,depth as dimensions
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have vi,depth as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename, ncid)

        ! Inquire the time dimension
        CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have vi,depth as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_mesh_field_dp_3D_ocean

! ===== Flexible looking for dimensions and variables =====
! =========================================================

  ! Look for dimensions
  SUBROUTINE inquire_dim_multopt( filename, ncid, dim_name_options, id_dim, dim_length, dim_name)
    ! Inquire if this file contains a dimension by name of dim_name.
    ! If so, return its length and identifier. If not, return -1 for both.
    !
    ! Supports providing multiple options for the dimension name, separated by two
    ! vertical bars || e.g. if we're looking for an X-dimension, we could do something like:
    !
    ! CALL inquire_dim_multopt( ncid, dim_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', dim_length, id_dim)
    !
    ! IF more than one match is found, crash.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name_options
    INTEGER,                             INTENT(OUT)   :: id_dim

    INTEGER,                                INTENT(OUT), OPTIONAL :: dim_length
    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: dim_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim_multopt'
    CHARACTER(LEN=256)                                 :: dim_name_options_parsed
    CHARACTER(LEN=256)                                 :: dim_name_options_redux
    INTEGER                                            :: i, n_matches
    INTEGER                                            :: dim_length_try, dim_length_match
    INTEGER                                            :: id_dim_try, id_dim_match
    CHARACTER(LEN=256)                                 :: dim_name_try, dim_name_match

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Parse field name options
    CALL parse_field_name_options( dim_name_options, dim_name_options_parsed)

    ! Try all options provided in dim_name_options

    dim_name_options_redux = TRIM( dim_name_options_parsed)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( dim_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        dim_name_try = dim_name_options_redux( 1:i-1)
        dim_name_options_redux = dim_name_options_redux( i+2:LEN_TRIM( dim_name_options_redux))

      ELSE
        ! Only one option is left over

        dim_name_try = dim_name_options_redux
        dim_name_options_redux( 1:LEN( dim_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_dim( filename, ncid, dim_name_try, dim_length_try, id_dim_try)

      IF (id_dim_try == -1) THEN
        ! No dimension by this name was found; try the next option
      ELSE
        ! A dimension by this name was found; hurray!
        n_matches  = n_matches + 1
        dim_length_match = dim_length_try
        id_dim_match     = id_dim_try
        dim_name_match   = dim_name_try
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( dim_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional dimension names were found in the NetCDF file
      dim_length_match = -1
      id_dim_match     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided dimension names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Copy to output arguments
    id_dim = id_dim_match
    IF (PRESENT( dim_name  )) dim_name   = dim_name_match
    IF (PRESENT( dim_length)) dim_length = dim_length_match

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim_multopt

  ! Look for variables
  SUBROUTINE inquire_var_multopt( filename, ncid, var_name_options, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    ! Inquire if this file contains a variable by name of var_name.
    ! If so, return its identifier. If not, return -1.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL inquire_var_multopt( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', id_var)
    !
    ! IF more than one match is found, crash.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    INTEGER,                             INTENT(OUT)   :: id_var

    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: var_name
    INTEGER,                                INTENT(OUT), OPTIONAL :: var_type
    INTEGER,                                INTENT(OUT), OPTIONAL :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS), INTENT(OUT), OPTIONAL :: dims_of_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var_multopt'
    CHARACTER(LEN=256)                                 :: var_name_options_parsed
    CHARACTER(LEN=256)                                 :: var_name_options_redux
    INTEGER                                            :: i, n_matches, id_var_try
    CHARACTER(LEN=256)                                 :: var_name_try, var_name_match

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Parse field name options
    CALL parse_field_name_options( var_name_options, var_name_options_parsed)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options_parsed)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name_try = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name_try = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( filename, ncid, var_name_try, id_var_try)

      IF (id_var_try == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches      = n_matches + 1
        id_var         = id_var_try
        var_name_match = var_name_try
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match. Inquire additional info on this variable.
      CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    END IF

    ! Copy to output arguments
    IF (PRESENT( var_name)) var_name = var_name_match

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var_multopt

! ===== Parse flexible dimension/variable names =====
! ===================================================

  SUBROUTINE parse_field_name_options( field_name_options, field_name_options_parsed)
    ! Check if a default set of field name options should be used.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=256),                  INTENT(OUT)   :: field_name_options_parsed

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'parse_field_name_options'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    field_name_options_parsed = field_name_options

    IF (INDEX( field_name_options,'default_options_') > 0) THEN
      ! Use one of the default options

      ! Dimensions
      IF     (field_name_options == 'default_options_x') THEN
        field_name_options_parsed = field_name_options_x
      ELSEIF (field_name_options == 'default_options_y') THEN
        field_name_options_parsed = field_name_options_y
      ELSEIF (field_name_options == 'default_options_zeta') THEN
        field_name_options_parsed = field_name_options_zeta
      ELSEIF (field_name_options == 'default_options_lon') THEN
        field_name_options_parsed = field_name_options_lon
      ELSEIF (field_name_options == 'default_options_lat') THEN
        field_name_options_parsed = field_name_options_lat
      ELSEIF (field_name_options == 'default_options_time') THEN
        field_name_options_parsed = field_name_options_time

      ! Variables
      ELSEIF (field_name_options == 'default_options_Hi') THEN
        field_name_options_parsed = field_name_options_Hi
      ELSEIF (field_name_options == 'default_options_Hb') THEN
        field_name_options_parsed = field_name_options_Hb
      ELSEIF (field_name_options == 'default_options_Hs') THEN
        field_name_options_parsed = field_name_options_Hs
      ELSEIF (field_name_options == 'default_options_SL') THEN
        field_name_options_parsed = field_name_options_SL

      ! Unrecognised default options
      ELSE
        CALL crash('unregocnised default field name option "' // TRIM( field_name_options) // '"')
      END IF

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE parse_field_name_options

  FUNCTION get_first_option_from_list( field_name_options) RESULT( field_name)
    ! Get the first option from a list of field name options

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=256)                                 :: field_name

    ! Local variables:
    INTEGER                                            :: i

    field_name( 1:256) = ' '

    i = INDEX( field_name_options,'||')

    IF (i > 0) THEN
      field_name = field_name_options( 1:i-1)
    ELSE
      field_name = TRIM( field_name_options)
    END IF

  END FUNCTION get_first_option_from_list

! ===== Read data from variables =====
! ====================================

  ! NOTE: only the Master actually reads data! Distributing to other processes
  !       must be done afterward

  SUBROUTINE read_var_master_int_0D(  filename, ncid, id_var, d)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,                             INTENT(OUT)   :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_int_0D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_int_0D

  SUBROUTINE read_var_master_int_1D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:    ), optional,INTENT(OUT)   :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_int_1D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_int_1D

  SUBROUTINE read_var_master_int_2D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:  ), optional,INTENT(OUT)   :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_int_2D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_int_2D

  SUBROUTINE read_var_master_int_3D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:), optional,INTENT(OUT)   :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_int_3D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_int_3D

  SUBROUTINE read_var_master_int_4D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:,:),optional,INTENT(OUT)  :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_int_4D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_int_4D

  SUBROUTINE read_var_master_dp_0D(  filename, ncid, id_var, d)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp),                            INTENT(OUT)   :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_dp_0D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_dp_0D

  SUBROUTINE read_var_master_dp_1D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:    ), optional,INTENT(OUT)   :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_dp_1D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_dp_1D

  SUBROUTINE read_var_master_dp_2D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:  ), Optional,INTENT(OUT)   :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_dp_2D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_dp_2D

  SUBROUTINE read_var_master_dp_3D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:), optional,INTENT(OUT)   :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_dp_3D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_dp_3D

  SUBROUTINE read_var_master_dp_4D(  filename, ncid, id_var, d, start, count)
    ! Read data from a NetCDF file
    !
    ! NOTE: only the Master actually reads data! Distributing to other processes
    !       must be done afterward

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:,:), optional, INTENT(OUT)   :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_master_dp_4D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master .AND. count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_master_dp_4D

! ===== Write data to variables =====
! ===================================

  ! NOTE: only the Master actually writes data! Gathering from other processes
  !       must be done beforehand

  SUBROUTINE write_var_master_int_0D(  filename, ncid, id_var, d)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,                             INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_int_0D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_int_0D

  SUBROUTINE write_var_master_int_1D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:    ), optional,INTENT(IN)    :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_int_1D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = 1
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      IF (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if ( count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

       ! Check if this dimension is large enough to read this amount of data
       IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_int_1D

  SUBROUTINE write_var_master_int_2D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:  ), optional,INTENT(IN)    :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_int_2D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      IF (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        IF (par%master .AND. start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
            TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_int_2D

  SUBROUTINE write_var_master_int_3D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:), optional,INTENT(IN)    :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_int_3D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      IF (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
          TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_int_3D

  SUBROUTINE write_var_master_int_4D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:,:),optional,INTENT(IN)   :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_int_4D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      IF (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if( count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
          TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if
    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_int_4D

  SUBROUTINE write_var_master_dp_0D(  filename, ncid, id_var, d)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp),                            INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_dp_0D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_dp_0D

  SUBROUTINE write_var_master_dp_1D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:    ), optional,INTENT(IN)    :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_dp_1D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied =  1
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      IF (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if ( count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        IF (var_name /= 'time') THEN
          ! Exception for time, because there the dimension is usually unlimited
          IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
            TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
        end if
      END IF

    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_dp_1D

  SUBROUTINE write_var_master_dp_2D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:  ), optional,INTENT(IN)    :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_dp_2D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (par%master) then
        if( count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
        IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_dp_2D

  SUBROUTINE write_var_master_dp_3D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:), optional,INTENT(IN)    :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_dp_3D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      IF (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if( count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_dp_3D

  SUBROUTINE write_var_master_dp_4D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file
    !
    ! NOTE: only the Master actually writes data! Gathering from other processes
    !       must be done beforehand

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:,:),optional,INTENT(IN)   :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_master_dp_4D'
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (par%master .AND. .NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (par%master .AND. ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    if (par%master .and. .not. present(d)) call crash('d needs to be present on master')

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (par%master .AND. ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      if (par%master) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    END IF
    IF (par%master .AND. ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      IF (par%master) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    END DO

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_master_dp_4D

! ===== Basic NetCDF wrapper functions =====
! ==========================================

  ! Inquire dimensions and variables
  SUBROUTINE inquire_dim( filename, ncid, dim_name, dim_length, id_dim)
    ! Inquire if this file contains a dimension by name of dim_name.
    ! If so, return its length and identifier; if not, return -1 for both.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name
    INTEGER,                             INTENT(OUT)   :: dim_length
    INTEGER,                             INTENT(OUT)   :: id_dim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    IF (par%master) THEN

      ! Check if a dimension of this name exists in the file
      nerr = NF90_INQ_DIMID( ncid, dim_name, id_dim)

      IF (nerr /= NF90_NOERR) THEN
        ! If a dimension by this name does not exist, return -1 for the length and ID
        id_dim     = -1
        dim_length = -1
      ELSE
        ! If a dimension by this name exists, find its length
        nerr = NF90_INQUIRE_DIMENSION( ncid, id_dim, len = dim_length)
        IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_DIMENSION failed for file "' // TRIM( filename) // '"!')
      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( id_dim    , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( dim_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim

  SUBROUTINE inquire_dim_info( filename, ncid, id_dim, dim_name, dim_length)
    ! Inquire some info of a dimension

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_dim

    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: dim_name
    INTEGER,                                INTENT(OUT), OPTIONAL :: dim_length

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim_info'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    IF (par%master) THEN
      ! Inquire some info on this variable
      nerr = NF90_INQUIRE_DIMENSION( ncid, id_dim, name = dim_name, len = dim_length)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_DIMENSION failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    IF (PRESENT( dim_name  )) CALL MPI_BCAST( dim_name  , 256, MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT( dim_length)) CALL MPI_BCAST( dim_length, 1  , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim_info

  SUBROUTINE inquire_var( filename, ncid, var_name, id_var)
    ! Inquire if this file contains a variable by name of var_name.
    ! If so, return its identifier. If not, return -1 for the identifier.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    INTEGER,                             INTENT(OUT)   :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var'
    CHARACTER                                          :: dummy1

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! To prevent "unused variable" compiler warnings
    dummy1 = filename( 1:1)

    IF (par%master) THEN

      ! Check if a variable of this name exists in the file
      nerr = NF90_INQ_VARID( ncid, var_name, id_var)

      IF (nerr /= NF90_NOERR) THEN
        ! If a variable by this name does not exist, return -1 for the ID
        id_var = -1
      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( id_var, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var

  SUBROUTINE inquire_var_info( filename, ncid, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    ! Inquire some info of a variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var

    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: var_name
    INTEGER,                                INTENT(OUT), OPTIONAL :: var_type
    INTEGER,                                INTENT(OUT), OPTIONAL :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS), INTENT(OUT), OPTIONAL :: dims_of_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var_info'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    IF (par%master) THEN
      ! Inquire some info on this variable
      nerr = NF90_INQUIRE_VARIABLE( ncid, id_var, name = var_name, xtype = var_type, ndims = ndims_of_var, dimids = dims_of_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_VARIABLE failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    IF (PRESENT( var_name    )) CALL MPI_BCAST( var_name    , 256              , MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT( var_type    )) CALL MPI_BCAST( var_type    , 1                , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT( ndims_of_var)) CALL MPI_BCAST( ndims_of_var, 1                , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT(  dims_of_var)) CALL MPI_BCAST( dims_of_var , NF90_MAX_VAR_DIMS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var_info

  ! Create new NetCDF file
  SUBROUTINE create_new_netcdf_file_for_writing( filename, ncid)
    ! Create a new NetCDF file in the specified location for writing.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_new_netcdf_file_for_writing'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if this file already exists
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF
    END IF

    ! Create the NetCDF file
    IF (par%master) THEN
      nerr = NF90_CREATE( filename, IOR( NF90_NOCLOBBER, NF90_NETCDF4), ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_CREATE failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Add some very basic info about the current simulation to the header
    call add_attribute_char( filename, ncid, NF90_GLOBAL, 'git commit hash', git_commit_hash)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_new_netcdf_file_for_writing

  ! Create dimensions, variables, and attributes
  SUBROUTINE create_dimension( filename, ncid, dim_name, dim_length, id_dim)
    ! Create a new dimension in a NetCDF file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name
    INTEGER,                             INTENT(IN)    :: dim_length
    INTEGER,                             INTENT(OUT)   :: id_dim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_dimension'
    INTEGER                                            :: dim_length_present

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety: check if a dimension by this name is already present in this file
    CALL inquire_dim( filename, ncid, dim_name, dim_length_present, id_dim)
    IF (id_dim /= -1) THEN
      !CALL crash('file "' // TRIM( filename) // '" already contains dimension "' // TRIM( dim_name) // '"!')
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Add the dimension
    IF (par%master) THEN
      nerr = NF90_DEF_DIM( ncid, dim_name, dim_length, id_dim)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_DEF_DIM failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    CALL MPI_BCAST( id_dim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_dimension

  SUBROUTINE create_variable( filename, ncid, var_name, var_type, dim_ids, id_var)
    ! Create a new variable in a NetCDF file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    INTEGER,                             INTENT(IN)    :: var_type
    INTEGER, DIMENSION(:),               INTENT(IN)    :: dim_ids
    INTEGER,                             INTENT(OUT)   :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_variable'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety: check if a variable by this name is already present in this file
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var /= -1) THEN
      !CALL crash('file "' // TRIM( filename) // '" already contains variable "' // TRIM( var_name) // '"!')
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Add the variable
    IF (par%master) THEN
      nerr = NF90_DEF_VAR( ncid, name = var_name, xtype = var_type, dimids = dim_ids, varid = id_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_DEF_VAR failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    CALL MPI_BCAST( id_var, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_variable

  SUBROUTINE create_scalar_variable( filename, ncid, var_name, var_type, id_var)
    ! Create a new scalar variable in a NetCDF file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    INTEGER,                             INTENT(IN)    :: var_type
    INTEGER,                             INTENT(OUT)   :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_scalar_variable'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety: check if a variable by this name is already present in this file
    CALL inquire_var( filename, ncid, var_name, id_var)
    IF (id_var /= -1) CALL crash('file "' // TRIM( filename) // '" already contains variable "' // TRIM( var_name) // '"!')

    ! Add the variable
    IF (par%master) THEN
      nerr = NF90_DEF_VAR( ncid, name = var_name, xtype = var_type, varid = id_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_DEF_VAR failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    CALL MPI_BCAST( id_var, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_scalar_variable

  SUBROUTINE add_attribute_int( filename, ncid, id_var, att_name, att_val)
    ! Add an integer-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_name
    INTEGER,                             INTENT(IN)    :: att_val

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_attribute_int'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Add the attribute
    IF (par%master) THEN
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_ATT failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_attribute_int

  SUBROUTINE add_attribute_dp( filename, ncid, id_var, att_name, att_val)
    ! Add a double-precision-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_name
    REAL(dp),                            INTENT(IN)    :: att_val

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_attribute_dp'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Add the attribute
    IF (par%master) THEN
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_ATT failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_attribute_dp

  SUBROUTINE add_attribute_char( filename, ncid, id_var, att_name, att_val)
    ! Add a character-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    INTEGER,                             INTENT(IN)    :: id_var
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_val

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_attribute_char'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Add the attribute
    IF (par%master) THEN
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_ATT failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_attribute_char

  ! Open and close a NetCDF file
  SUBROUTINE open_existing_netcdf_file_for_reading( filename, ncid)
    ! Open the NetCDF file in the specified location for reading only,
    ! and return its identifier.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'open_netcdf_file_for_reading'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Open the NetCDF file with read-only access
    IF (par%master) THEN
      nerr = NF90_OPEN( TRIM( filename), NF90_NOWRITE, ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_OPEN failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE open_existing_netcdf_file_for_reading

  SUBROUTINE open_existing_netcdf_file_for_writing( filename, ncid)
    ! Open an existing NetCDF file in data mode
    ! In data mode, no new dimensions, variables, or attributes can be created,
    ! but data can be written to existing variables.
    ! When opening an existing NetCDF file, it is by default in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'open_existing_netcdf_file_for_writing'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Open the NetCDF file with read+write access
    IF (par%master) THEN
      nerr = NF90_OPEN( TRIM( filename), IOR( NF90_WRITE, NF90_SHARE), ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_OPEN failed for file "' // TRIM( filename) // '" beeperdebeep!')
    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE open_existing_netcdf_file_for_writing

  SUBROUTINE close_netcdf_file( ncid)
    ! Close an opened NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'close_netcdf_file'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Close netCDF file:
    IF (par%master) THEN
      nerr = NF90_CLOSE( ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_CLOSE failed!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE close_netcdf_file

END MODULE netcdf_basic
