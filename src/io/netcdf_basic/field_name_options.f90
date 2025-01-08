module field_name_options
  ! Lists of options for flexibly-named dimensions/variables in NetCDF in/output files

  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash

  implicit none

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

  ! Ice model variables
  character(len=1024), parameter :: field_name_options_Hi             = 'Hi||thickness||lithk'
  character(len=1024), parameter :: field_name_options_Hb             = 'Hb||bed||topg'
  character(len=1024), parameter :: field_name_options_Hs             = 'Hs||surface||orog'
  character(len=1024), parameter :: field_name_options_SL             = 'SL'
  character(len=1024), parameter :: field_name_options_dHb            = 'dHb'
  character(len=1024), parameter :: field_name_options_Ti             = 'Ti'
  character(len=1024), parameter :: field_name_options_T_ocean        = 'T_ocean||t_ocean||t_an'
  character(len=1024), parameter :: field_name_options_S_ocean        = 'S_ocean||s_ocean||s_an'

contains

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

end module field_name_options