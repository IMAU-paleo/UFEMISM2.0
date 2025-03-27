module netcdf_write_var_primary
  !< Write data to variables
  !< NOTE: only the primary actually writes data! Gathering from other processes must be done beforehand

  use precisions, only: dp, int8
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use netcdf_field_name_options
  use netcdf_basic_wrappers
  use netcdf, only: NF90_MAX_VAR_DIMS, NF90_DOUBLE, NF90_FLOAT, NF90_INT, NF90_PUT_VAR, NF90_INT64

  implicit none

  private

  public :: write_var_primary

  interface write_var_primary
    procedure write_var_primary_int_0D
    procedure write_var_primary_int_1D
    procedure write_var_primary_int_2D
    procedure write_var_primary_int_3D
    procedure write_var_primary_int_4D
    procedure write_var_primary_int8_2D
    procedure write_var_primary_dp_0D
    procedure write_var_primary_dp_1D
    procedure write_var_primary_dp_2D
    procedure write_var_primary_dp_3D
    procedure write_var_primary_dp_4D
  end interface write_var_primary

contains

  subroutine write_var_primary_int_0D( filename, ncid, id_var, d)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    integer,          intent(in   ) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_var_primary_int_0D'
    character(len=1024)            :: var_name
    integer                        :: var_type
    integer                        :: ndims_of_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_int_0D

  subroutine write_var_primary_int_1D(  filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                     intent(in   ) :: filename
    integer,                              intent(in   ) :: ncid
    integer,                              intent(in   ) :: id_var
    integer,  dimension(:    ),           intent(in   ) :: d
    integer,  dimension(1    ), optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'write_var_primary_int_1D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = 1
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if ( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_int_1D

  subroutine write_var_primary_int_2D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    integer,                            intent(in   ) :: ncid
    integer,                            intent(in   ) :: id_var
    integer,  dimension(:,:),           intent(in   ) :: d
    integer,  dimension(2  ), optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_int_2D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
            trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_int_2D

  subroutine write_var_primary_int_3D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(in   ) :: id_var
    integer, dimension(:,:,:),           intent(in   ) :: d
    integer, dimension(3),     optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_int_3D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1 /)
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
          trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_int_3D

  subroutine write_var_primary_int_4D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                      intent(in   ) :: filename
    integer,                               intent(in   ) :: ncid
    integer,                               intent(in   ) :: id_var
    integer, dimension(:,:,:,:),           intent(in   ) :: d
    integer, dimension(4),       optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_int_4D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_INT)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1, 1 /)
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
          trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if
    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_int_4D

  subroutine write_var_primary_int8_2D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                        intent(in   ) :: filename
    integer,                                 intent(in   ) :: ncid
    integer,                                 intent(in   ) :: id_var
    integer(int8), dimension(:,:),           intent(in   ) :: d
    integer,       dimension(2),   optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_int8_2D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_INT64)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
            trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_int8_2D

  subroutine write_var_primary_dp_0D( filename, ncid, id_var, d)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    real(dp),         intent(in   ) :: d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_var_primary_dp_0D'
    character(len=1024)            :: var_name
    integer                        :: var_type
    integer                        :: ndims_of_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 0) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_dp_0D

  subroutine write_var_primary_dp_1D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    integer,                          intent(in   ) :: ncid
    integer,                          intent(in   ) :: id_var
    real(dp), dimension(:),           intent(in   ) :: d
    integer,  dimension(1), optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_dp_1D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied =  1
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
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
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_dp_1D

  subroutine write_var_primary_dp_2D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    integer,                            intent(in   ) :: ncid
    integer,                            intent(in   ) :: id_var
    real(dp), dimension(:,:),           intent(in   ) :: d
    integer,  dimension(2),   optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_dp_2D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1 /)
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      if (par%primary) then
        if( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_dp_2D

  subroutine write_var_primary_dp_3D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                     intent(in   ) :: filename
    integer,                              intent(in   ) :: ncid
    integer,                              intent(in   ) :: id_var
    real(dp), dimension(:,:,:),           intent(in   ) :: d
    integer,  dimension(3),     optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_dp_3D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1 /)
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if( count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_dp_3D

  subroutine write_var_primary_dp_4D( filename, ncid, id_var, d, start, count)
    ! Write data to a NetCDF file

    ! In/output variables:
    character(len=*),                       intent(in   ) :: filename
    integer,                                intent(in   ) :: ncid
    integer,                                intent(in   ) :: id_var
    real(dp), dimension(:,:,:,:),           intent(in   ) :: d
    integer,  dimension(4),       optional, intent(in   ) :: start, count

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_var_primary_dp_4D'
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

#if (DO_ASSERTIONS)
    ! inquire some info on this variable
    call inquire_var_info( filename, ncid, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (par%primary .and. .not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) &
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    if (par%primary .and. ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
#endif

    ! Set start and count
    if (present( start)) then
      start_applied = start
    else
      start_applied = (/ 1, 1, 1, 1 /)
    end if
    if (par%primary .and. any( start_applied == 0)) call crash('start must be positive!')

    if (present( count)) then
      count_applied = count
    else
      if (par%primary) then
        count_applied = shape(d)
      else
        count_applied = 1
      end if
    end if
    if (par%primary .and. any( count_applied == 0)) call crash('count must be positive!')

#if (DO_ASSERTIONS)
    ! Check sizes of dimensions
    do di = 1, ndims_of_var

      ! Check size of this dimension in the file
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      if (par%primary) then
        ! Check if the combination of dimension size, start, and count, matches the size of d
        if(count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

        ! Check if this dimension is large enough to read this amount of data
        if (start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
        trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      end if

    end do
#endif

    ! Write the data
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied), &
        filename = filename, dimvarname = var_name)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_primary_dp_4D

end module netcdf_write_var_primary