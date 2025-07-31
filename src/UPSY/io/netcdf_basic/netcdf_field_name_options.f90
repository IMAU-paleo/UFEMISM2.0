module netcdf_field_name_options
  !< Lists of options for flexibly-named dimensions/variables in NetCDF in/output files

  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use netcdf, only: NF90_MAX_VAR_DIMS
  use netcdf_basic_wrappers, only: inquire_var, inquire_var_info, inquire_dim

  implicit none

  private

  public :: field_name_options_x, field_name_options_y, field_name_options_zeta, field_name_options_lon, &
    field_name_options_lat, field_name_options_time, field_name_options_month, field_name_options_depth, &
    field_name_options_dim_nV, field_name_options_dim_nTri, field_name_options_dim_nC_mem, &
    field_name_options_dim_nE, field_name_options_dim_nVor, field_name_options_dim_two, field_name_options_dim_three, &
    field_name_options_dim_four, field_name_options_V, field_name_options_Tri, field_name_options_nC, &
    field_name_options_C, field_name_options_niTri, field_name_options_iTri, field_name_options_VBI, &
    field_name_options_Tricc, field_name_options_TriC, field_name_options_TriBI, field_name_options_E, &
    field_name_options_VE, field_name_options_EV, field_name_options_ETri, field_name_options_TriE, field_name_options_EBI, &
    field_name_options_vi2vori, field_name_options_ti2vori, field_name_options_ei2vori, field_name_options_vori2vi, &
    field_name_options_vori2ti, field_name_options_vori2ei, field_name_options_Vor, field_name_options_VornC, &
    field_name_options_VorC, field_name_options_nVVor, field_name_options_VVor, field_name_options_TriGC, field_name_options_TriA, &
    field_name_options_A, field_name_options_EA, field_name_options_R, field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, &
    field_name_options_SL, field_name_options_dHb, field_name_options_Ti, field_name_options_T_ocean, &
    field_name_options_S_ocean, field_name_options_dT_ocean, field_name_options_insolation, field_name_options_sealevel, field_name_options_GI, &
    field_name_options_CO2

  public :: inquire_dim_multopt, inquire_var_multopt, get_first_option_from_list

  ! Dimensions
  character(len=1024), parameter :: field_name_options_x              = 'x||X||x1||X1||nx||NX||x-coordinate||X-coordinate||easting||Easting'
  character(len=1024), parameter :: field_name_options_y              = 'y||Y||y1||Y1||ny||NY||y-coordinate||Y-coordinate||northing||Northing'
  character(len=1024), parameter :: field_name_options_zeta           = 'zeta||Zeta'
  character(len=1024), parameter :: field_name_options_lon            = 'lon||Lon||long||Long||longitude||Longitude'
  character(len=1024), parameter :: field_name_options_lat            = 'lat||Lat||latitude||Latitude'
  character(len=1024), parameter :: field_name_options_time           = 'time||Time||t||nt'
  character(len=1024), parameter :: field_name_options_month          = 'month||Month'
  character(len=1024), parameter :: field_name_options_depth          = 'depth||Depth'

  ! Mesh data
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
  character(len=1024), parameter :: field_name_options_TriE           = 'TriE'
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
  character(len=1024), parameter :: field_name_options_TriA           = 'TriA'
  character(len=1024), parameter :: field_name_options_A              = 'A'
  character(len=1024), parameter :: field_name_options_EA             = 'EA'
  character(len=1024), parameter :: field_name_options_R              = 'R'

  ! Ice model variables
  character(len=1024), parameter :: field_name_options_Hi             = 'Hi||thickness||lithk||ice_thickness'
  character(len=1024), parameter :: field_name_options_Hb             = 'Hb||bed||topg||bed_topography'
  character(len=1024), parameter :: field_name_options_Hs             = 'Hs||surface||orog||surface_topography'
  character(len=1024), parameter :: field_name_options_SL             = 'SL'
  character(len=1024), parameter :: field_name_options_dHb            = 'dHb'
  character(len=1024), parameter :: field_name_options_Ti             = 'Ti'
  character(len=1024), parameter :: field_name_options_T_ocean        = 'T_ocean||t_ocean||t_an||votemper'
  character(len=1024), parameter :: field_name_options_S_ocean        = 'S_ocean||s_ocean||s_an||vosaline'
  character(len=1024), parameter :: field_name_options_dT_ocean       = 'dT||dT_ocean||dTo'

  ! Global forcing variables
  character(len=1024), parameter :: field_name_options_insolation     = 'Q_TOA'
  character(len=1024), parameter :: field_name_options_sealevel       = 'SL||sea_level||sl'
  character(len=1024), parameter :: field_name_options_GI             = 'GI||gi||Glacial_Index||glacial_index||GlacialIndex'
  character(len=1024), parameter :: field_name_options_CO2            = 'CO2||co2'

contains

  subroutine inquire_dim_multopt( filename, ncid, dim_name_options, id_dim, dim_length, dim_name)
    !< Inquire if this file contains a dimension by name of the dim_name options.
    !< if so, return its length and identifier. if not, return -1 for both.

    ! Supports providing multiple options for the dimension name, separated by two
    ! vertical bars || e.g. if we're looking for an X-dimension, we could do something like:
    !
    ! call inquire_dim_multopt( ncid, dim_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', dim_length, id_dim)
    !
    ! If more than one match is found, crash.

    ! In/output variables:
    character(len=*),              intent(in   ) :: filename
    integer,                       intent(in   ) :: ncid
    character(len=*),              intent(in   ) :: dim_name_options
    integer,                       intent(  out) :: id_dim
    integer,             optional, intent(  out) :: dim_length
    character(len=1024), optional, intent(  out) :: dim_name

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
      if (len_trim( dim_name_options_redux) == 0) exit

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

  subroutine inquire_var_multopt( filename, ncid, var_name_options, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    !< Inquire if this file contains a variable by name of var_name.
    !< if so, return its identifier. if not, return -1.

    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! call inquire_var_multopt( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', id_var)
    !
    ! If more than one match is found, crash.

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
      if (len_trim( var_name_options_redux) == 0) exit

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

  subroutine parse_field_name_options( field_name_options, field_name_options_parsed)
    !< Check if a default set of field name options should be used.

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

      select case (field_name_options)

      case default
        call crash('unregocnised default field name option "' // trim( field_name_options) // '"')

      ! Dimensions
      case ('default_options_x')
        field_name_options_parsed = field_name_options_x
      case ('default_options_y')
        field_name_options_parsed = field_name_options_y
      case ('default_options_zeta')
        field_name_options_parsed = field_name_options_zeta
      case ('default_options_lon')
        field_name_options_parsed = field_name_options_lon
      case ('default_options_lat')
        field_name_options_parsed = field_name_options_lat
      case ('default_options_time')
        field_name_options_parsed = field_name_options_time

      ! Variables
      case ('default_options_Hi')
        field_name_options_parsed = field_name_options_Hi
      case ('default_options_Hb')
        field_name_options_parsed = field_name_options_Hb
      case ('default_options_Hs')
        field_name_options_parsed = field_name_options_Hs
      case ('default_options_SL')
        field_name_options_parsed = field_name_options_SL

      end select

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine parse_field_name_options

  function get_first_option_from_list( field_name_options) result( field_name)
    !< Get the first option from a list of field name options

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

end module netcdf_field_name_options
