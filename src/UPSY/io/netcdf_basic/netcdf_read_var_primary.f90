module netcdf_read_var_primary
  !< Read data from variables
  !< NOTE: only the primary actually reads data! Distributing to other processes must be done afterward

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use netcdf_field_name_options
  use netcdf_basic_wrappers
  use netcdf, only: NF90_MAX_VAR_DIMS, NF90_DOUBLE, NF90_FLOAT, NF90_INT, NF90_GET_VAR

  implicit none

  private

  public :: read_var_primary

  interface read_var_primary
    procedure read_var_primary_int_0D
    procedure read_var_primary_int_1D
    procedure read_var_primary_int_2D
    procedure read_var_primary_int_3D
    procedure read_var_primary_int_4D
    procedure read_var_primary_dp_0D
    procedure read_var_primary_dp_1D
    procedure read_var_primary_dp_2D
    procedure read_var_primary_dp_3D
    procedure read_var_primary_dp_4D
  end interface read_var_primary

contains

subroutine read_var_primary_int_0D(  filename, ncid, id_var, d)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*), intent(in   ) :: filename
  integer,          intent(in   ) :: ncid
  integer,          intent(in   ) :: id_var
  integer,          intent(  out) :: d

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'read_var_primary_int_0D'
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
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_int_0D

subroutine read_var_primary_int_1D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                 intent(in   ) :: filename
  integer,                          intent(in   ) :: ncid
  integer,                          intent(in   ) :: id_var
  integer,  dimension(:),           intent(  out) :: d
  integer,  dimension(1), optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_int_1D'
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

  ! Set start and count
  if (present( start)) then
    start_applied = start
  else
    start_applied = (/ 1 /)
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_int_1D

subroutine read_var_primary_int_2D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                   intent(in   ) :: filename
  integer,                            intent(in   ) :: ncid
  integer,                            intent(in   ) :: id_var
  integer,  dimension(:,:),           intent(  out) :: d
  integer,  dimension(2),   optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_int_2D'
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_int_2D

subroutine read_var_primary_int_3D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                    intent(in   ) :: filename
  integer,                             intent(in   ) :: ncid
  integer,                             intent(in   ) :: id_var
  integer, dimension(:,:,:),           intent(  out) :: d
  integer, dimension(3),     optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_int_3D'
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_int_3D

subroutine read_var_primary_int_4D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                      intent(in   ) :: filename
  integer,                               intent(in   ) :: ncid
  integer,                               intent(in   ) :: id_var
  integer, dimension(:,:,:,:),           intent(  out) :: d
  integer, dimension(4),       optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_int_4D'
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_int_4D

subroutine read_var_primary_dp_0D( filename, ncid, id_var, d)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*), intent(in   ) :: filename
  integer,          intent(in   ) :: ncid
  integer,          intent(in   ) :: id_var
  real(dp),         intent(  out) :: d

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'read_var_primary_dp_0D'
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
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_dp_0D

subroutine read_var_primary_dp_1D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                 intent(in   ) :: filename
  integer,                          intent(in   ) :: ncid
  integer,                          intent(in   ) :: id_var
  real(dp), dimension(:),           intent(  out) :: d
  integer,  dimension(1), optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_dp_1D'
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

  ! Set start and count
  if (present( start)) then
    start_applied = start
  else
    start_applied = (/ 1 /)
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_dp_1D

subroutine read_var_primary_dp_2D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                     intent(in   ) :: filename
  integer,                              intent(in   ) :: ncid
  integer,                              intent(in   ) :: id_var
  real(dp), dimension(:,:  ),           intent(  out) :: d
  integer,  dimension(2),     optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_dp_2D'
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_dp_2D

subroutine read_var_primary_dp_3D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                     intent(in   ) :: filename
  integer,                              intent(in   ) :: ncid
  integer,                              intent(in   ) :: id_var
  real(dp), dimension(:,:,:),           intent(  out) :: d
  integer,  dimension(3),     optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_dp_3D'
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_dp_3D

subroutine read_var_primary_dp_4D( filename, ncid, id_var, d, start, count)
  ! Read data from a NetCDF file
  !
  ! NOTE: only the primary actually reads data! Distributing to other processes
  !       must be done afterward

  ! In/output variables:
  character(len=*),                       intent(in   ) :: filename
  integer,                                intent(in   ) :: ncid
  integer,                                intent(in   ) :: id_var
  real(dp), dimension(:,:,:,:),           intent(  out) :: d
  integer,  dimension(4),       optional, intent(in   ) :: start, count

  ! Local variables:
  character(len=1024), parameter          :: routine_name = 'read_var_primary_dp_4D'
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

  ! Check sizes of dimensions
  do di = 1, ndims_of_var

    ! Check size of this dimension in the file
    call inquire_dim_info( filename, ncid, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

    ! Check if the combination of dimension size, start, and count, matches the size of d
    if (par%primary .and. count_applied( di) /= SIZE( d,di)) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // trim( filename) // &
      '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

    ! Check if this dimension is large enough to read this amount of data
    if (par%primary .and. start_applied( di) + count_applied( di) - 1 > dim_length) call crash('error for dimension "' // trim( dim_name) // '" of variable "' // trim( var_name) // '" in file "' // &
      trim( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

  end do

  ! Read the data
  if (par%primary) then
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var, d, start_applied, count_applied), &
      filename = filename, dimvarname = var_name)
  end if

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine read_var_primary_dp_4D

end module netcdf_read_var_primary
