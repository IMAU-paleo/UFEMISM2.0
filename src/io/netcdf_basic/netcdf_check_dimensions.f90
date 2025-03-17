module netcdf_check_dimensions
  !< Check that the dimension variables in a file are valid

  use assertions_basic
  use mpi_basic, only: par, sync
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use netcdf, only: NF90_MAX_VAR_DIMS, NF90_DOUBLE, NF90_FLOAT, NF90_UNLIMITED, NF90_INT
  use netcdf_field_name_options
  use netcdf_read_var_primary

  implicit none

  private

  public :: check_x, check_y, check_lon, check_lat, check_mesh_dimensions, check_zeta, &
    check_month, check_time, check_depth

contains

subroutine check_x( filename, ncid)
  !< Check if this file contains a valid x dimension and variable

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
  call read_var_primary( filename, ncid, id_var, x)

  if (par%primary) call assert( (.not. any( isnan( x))), 'found NaNs in x')

  ! Check grid spacing
  if (par%primary) then
    dx = x( 2) - x( 1)
    do i = 2, n
      dxp = x( i) - x( i-1)
      if (abs( 1._dp - dxp / dx) > 1E-5_dp) call crash('x coordinate in file "' // trim( filename) // '" is irregular!')
    end do
  end if
  call sync

  ! Clean up after yourself
  deallocate( x)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine check_x

subroutine check_y( filename, ncid)
  !< Check if this file contains a valid y dimension and variable

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
  call read_var_primary( filename, ncid, id_var, y)

  if (par%primary) call assert( (.not. any( isnan( y))), 'found NaNs in y')

  ! Check grid spacing
  if (par%primary) then
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

subroutine check_lon( filename, ncid)
  !< Check if this file contains a valid longitude dimension and variable

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
  call read_var_primary( filename, ncid, id_var, lon)

  if (par%primary) call assert( (.not. any( isnan( lon))), 'found NaNs in lon')

  ! Check grid spacing
  if (par%primary) then
    dlon = lon( 2) - lon( 1)
    do i = 2, n
      dlonp = lon( i) - lon( i-1)
      if (abs( 1._dp - dlonp / dlon) > 1E-5_dp) call crash('longitude coordinate in file "' // trim( filename) // '" is irregular!')
    end do
  end if
  call sync

  ! Clean up after yourself
  deallocate( lon)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine check_lon

subroutine check_lat( filename, ncid)
  !< Check if this file contains a valid latitude dimension and variable

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
  call read_var_primary( filename, ncid, id_var, lat)

  if (par%primary) call assert( (.not. any( isnan( lat))), 'found NaNs in lat')

  ! Check grid spacing
  if (par%primary) then
    dlat = lat( 2) - lat( 1)
    do i = 2, n
      dlatp = lat( i) - lat( i-1)
      if (abs( 1._dp - dlatp / dlat) > 1E-5_dp) call crash('latitude coordinate in file "' // trim( filename) // '" is irregular!')
    end do
  end if
  call sync

  ! Clean up after yourself
  deallocate( lat)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine check_lat

subroutine check_mesh_dimensions( filename, ncid)
  !< Check if this file contains valid mesh dimensions and variables

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

subroutine check_zeta( filename, ncid)
  !< Check if this file contains a valid zeta dimension and variable

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
  call read_var_primary( filename, ncid, id_var, zeta)

  ! Check validity
  if (par%primary) then
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
  !< Check if this file contains a valid month dimension (we don't really care about the variable)

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
  !< Check if this file contains a valid time dimension and variable

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
    call read_var_primary( filename, ncid, id_var, time)

    ! Check validity
    if (par%primary) call assert( (.not. any( isnan( time))), 'found NaN in time')

    ! Clean up after yourself
    deallocate( time)

  end if ! if (n > 0) then

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine check_time

subroutine check_depth( filename, ncid)
  !< Check if this file contains a valid depth dimension and variable

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
  call read_var_primary( filename, ncid, id_var, depth)

  ! Check validity
  if (par%primary) then
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

end module netcdf_check_dimensions