module transects_main

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_ALLREDUCE, MPI_IN_PLACE, &
    MPI_INTEGER, MPI_SUM
  use mpi_basic, only: par, sync
  use mpi_distributed_memory, only: partition_list
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, crash
  use model_configuration, only: C
  use region_types, only: type_model_region
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use transect_types, only: type_transect
  use string_module, only: separate_strings_by_double_vertical_bars
  use netcdf_io_main
  use netcdf, only: NF90_UNLIMITED
  use netcdf_write_field_transect
  use remapping_transects, only: map_from_mesh_vertices_to_transect_2D, map_from_mesh_vertices_to_transect_3D, &
    map_from_mesh_triangles_to_transect_2D, map_from_mesh_triangles_to_transect_3D
  use parameters, only: ice_density
  use mesh_zeta, only: vertical_average
  use ice_geometry_basics, only: thickness_above_floatation
  use mpi_distributed_memory, only: gather_to_all
  use interpolation, only: linint_points
  use netcdf, only: NF90_DOUBLE

  implicit none

  private

  public :: initialise_transects, write_to_transect_netcdf_output_files

contains

  subroutine initialise_transects( region)

    ! In/output variables
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'initialise_transects'
    character(len=1024)                            :: transects_str
    character(len=1024), dimension(:), allocatable :: transect_strs
    integer                                        :: it

    ! Add routine to path
    call init_routine( routine_name)

    select case (region%name)
    case default
      call crash('Unknown region "' // trim(region%name) // '"')
    case ('NAM')
      transects_str = C%transects_NAM
    case ('EAS')
      transects_str = C%transects_EAS
    case ('GRL')
      transects_str = C%transects_GRL
    case ('ANT')
      transects_str = C%transects_ANT
    end select

    if (transects_str == '') then
      allocate( region%transects( 0))
      call finalise_routine( routine_name)
      return
    end if

    call separate_strings_by_double_vertical_bars( transects_str, transect_strs)

    allocate( region%transects( size(transect_strs)))

    do it = 1, size(region%transects)
      call initialise_transect( region%mesh, region%transects(it), transect_strs(it))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transects

  subroutine initialise_transect( mesh, transect, transect_str)

    ! In/output variables
    type(type_mesh),     intent(in   ) :: mesh
    type(type_transect), intent(  out) :: transect
    character(len=*),    intent(in   ) :: transect_str

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'initialise_transect'
    character(len=1024)                   :: source, name, filename
    real(dp), dimension(:,:), allocatable :: waypoints
    real(dp)                              :: dx

    ! Add routine to path
    call init_routine( routine_name)

    call parse_transect_str( transect_str, source, name, filename, dx)

    transect%name = name
    transect%dx   = dx

    if (par%primary) write(0,*) '  Initialising output transect ', &
      colour_string( trim( transect%name),'light blue'), '...'

    select case (source)
    case default
      call crash('invalid transect source "' // trim( source) // '"')
    case ('hardcoded')
      call initialise_transect_waypoints_hardcoded( mesh, name, waypoints)
    case ('read_from_file')
      call initialise_transect_waypoints_from_file( filename, waypoints)
    end select

    call calc_transect_vertices_from_waypoints( transect, waypoints, dx)

    transect%nz = mesh%nz
    allocate( transect%zeta( mesh%nz))
    transect%zeta = mesh%zeta

    call calc_velocity_weights( transect)

    call create_transect_netcdf_output_file( transect)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transect

  subroutine parse_transect_str( transect_str, source, name, filename, dx)

    ! The transect should be specified by a name (either one of the hard-coded options,
    ! or a direction to an external file), and a resolution, e.g.:
    !
    !   positiveyaxis,dx=5e3
    !
    ! This indicates using the hard-coded option "positiveyaxis" with a resolution of 5 km
    !
    !   file:transect_MISMIP+_crossshelf.cfg,dx=2e3
    !
    ! This indicates using the external file "transect_MISMIP+_crossshelf.cfg, with a resolution of 2 km

    ! In/output variables:
    character(len=*), intent(in   ) :: transect_str
    character(len=*), intent(  out) :: source, name, filename
    real(dp),         intent(  out) :: dx

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'parse_transect_str'
    integer                        :: i

    ! Add routine to path
    call init_routine( routine_name)

    i = index( transect_str,',dx=')
    if (i==0) call crash('invalid transect string "' // trim(transect_str) // '" - could not find resolution specification!')

    name = transect_str( 1:i-1)
    read( transect_str( i+4:len_trim(transect_str)),*) dx

    ! Separate source and name (and filename)
    i = index(name,'file:')
    if (i==0) then
      source = 'hardcoded'
    elseif (i==1) then
      source = 'read_from_file'
      filename = name(6:len_trim(name))
      i = index( filename,'/',back=.true.)
      name = filename(i+1:len_trim(filename)-4)
    else
      call crash('invalid transect string "' // trim(transect_str))
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine parse_transect_str

  subroutine initialise_transect_waypoints_hardcoded( mesh, name, waypoints)
    !< The native transect options hard-coded in UFEMISM

    ! In/output variables
    type(type_mesh),                       intent(in   ) :: mesh
    character(len=*),                      intent(in   ) :: name
    real(dp), dimension(:,:), allocatable, intent(  out) :: waypoints

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_transect_waypoints_hardcoded'

    ! Add routine to path
    call init_routine( routine_name)

    select case (name)
    case default
      call crash('unknown native transect option "' // trim(name) // '"')

    ! == Idealised ==
    ! ===============

    ! Outward from [x,y] = [0,0]

    case('east')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [mesh%xmax,0._dp]

    case('west')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [mesh%xmin,0._dp]

    case('south')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [0._dp,mesh%ymin]

    case('north')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [0._dp,mesh%ymax]

    case('northeast')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [mesh%xmax,mesh%ymax]

    case('southeast')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [mesh%xmax,mesh%ymin]

    case('northwest')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [mesh%xmin,mesh%ymax]

    case('southwest')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,0._dp]
      waypoints(2,:) = [mesh%xmin,mesh%ymin]

    ! Across the entire domain

    case('westeast')
      allocate(waypoints(2,2))
      waypoints(1,:) = [mesh%xmin,0._dp]
      waypoints(2,:) = [mesh%xmax,0._dp]

    case('southnorth')
      allocate(waypoints(2,2))
      waypoints(1,:) = [0._dp,mesh%ymin]
      waypoints(2,:) = [0._dp,mesh%ymax]

    case('ISMIP-HOM')
      allocate(waypoints(2,2))
      waypoints(1,:) = [mesh%xmin/2._dp,mesh%ymin/4._dp]
      waypoints(2,:) = [mesh%xmax/2._dp,mesh%ymin/4._dp]

    ! == Realistic ==
    ! ===============

    ! Antarctica

    case('PineIsland_centralflowline')

      allocate(waypoints(72,2))

      waypoints( 1,:) = [ -1.581444261355978e6_dp, -0.030311971888969e6_dp]
      waypoints( 2,:) = [ -1.582246435803775e6_dp, -0.035247204016772e6_dp]
      waypoints( 3,:) = [ -1.582303292344950e6_dp, -0.040246880739694e6_dp]
      waypoints( 4,:) = [ -1.582185249855441e6_dp, -0.045245487142550e6_dp]
      waypoints( 5,:) = [ -1.582081544495048e6_dp, -0.050244411546682e6_dp]
      waypoints( 6,:) = [ -1.582662969743382e6_dp, -0.055210490954081e6_dp]
      waypoints( 7,:) = [ -1.583493442845907e6_dp, -0.060141040053926e6_dp]
      waypoints( 8,:) = [ -1.584240740762885e6_dp, -0.065084879232468e6_dp]
      waypoints( 9,:) = [ -1.584851396852385e6_dp, -0.070047449044607e6_dp]
      waypoints(10,:) = [ -1.585308445570337e6_dp, -0.075026515871777e6_dp]
      waypoints(11,:) = [ -1.585578821263540e6_dp, -0.080019200218350e6_dp]
      waypoints(12,:) = [ -1.585586411846643e6_dp, -0.085019194456652e6_dp]
      waypoints(13,:) = [ -1.585392585848335e6_dp, -0.090015436192435e6_dp]
      waypoints(14,:) = [ -1.585000000000000e6_dp, -0.095000000000000e6_dp]
      waypoints(15,:) = [ -1.584607414151665e6_dp, -0.099984563807565e6_dp]
      waypoints(16,:) = [ -1.584178880875366e6_dp, -0.104966165882310e6_dp]
      waypoints(17,:) = [ -1.583919194172270e6_dp, -0.109959417610016e6_dp]
      waypoints(18,:) = [ -1.583954223315052e6_dp, -0.114959294904426e6_dp]
      waypoints(19,:) = [ -1.584220471737793e6_dp, -0.119952201049887e6_dp]
      waypoints(20,:) = [ -1.584610216958133e6_dp, -0.124936987771827e6_dp]
      waypoints(21,:) = [ -1.585093965814203e6_dp, -0.129913531456378e6_dp]
      waypoints(22,:) = [ -1.585561266089043e6_dp, -0.134891646607029e6_dp]
      waypoints(23,:) = [ -1.585996797933316e6_dp, -0.139872641689602e6_dp]
      waypoints(24,:) = [ -1.586344192008271e6_dp, -0.144860558825711e6_dp]
      waypoints(25,:) = [ -1.586579829858230e6_dp, -0.149855003219602e6_dp]
      waypoints(26,:) = [ -1.586688774455813e6_dp, -0.154853816186163e6_dp]
      waypoints(27,:) = [ -1.586720851753739e6_dp, -0.159853713289800e6_dp]
      waypoints(28,:) = [ -1.586818829539264e6_dp, -0.164852753232983e6_dp]
      waypoints(29,:) = [ -1.586978885063885e6_dp, -0.169850190799280e6_dp]
      waypoints(30,:) = [ -1.587219321903611e6_dp, -0.174844406466039e6_dp]
      waypoints(31,:) = [ -1.587355556033750e6_dp, -0.179842550147626e6_dp]
      waypoints(32,:) = [ -1.587304596740988e6_dp, -0.184842290455930e6_dp]
      waypoints(33,:) = [ -1.587008845897270e6_dp, -0.189833535935612e6_dp]
      waypoints(34,:) = [ -1.586313384448024e6_dp, -0.194784933048813e6_dp]
      waypoints(35,:) = [ -1.585415136891800e6_dp, -0.199703586454750e6_dp]
      waypoints(36,:) = [ -1.584571422291555e6_dp, -0.204631886940042e6_dp]
      waypoints(37,:) = [ -1.583927625182513e6_dp, -0.209590266240014e6_dp]
      waypoints(38,:) = [ -1.583556350183616e6_dp, -0.214576462673689e6_dp]
      waypoints(39,:) = [ -1.583458692634559e6_dp, -0.219575508883028e6_dp]
      waypoints(40,:) = [ -1.583673903094605e6_dp, -0.224570875181698e6_dp]
      waypoints(41,:) = [ -1.584160056601684e6_dp, -0.229547184533774e6_dp]
      waypoints(42,:) = [ -1.584950493050585e6_dp, -0.234484310238093e6_dp]
      waypoints(43,:) = [ -1.586093394698516e6_dp, -0.239351935513642e6_dp]
      waypoints(44,:) = [ -1.588274448157067e6_dp, -0.243851158314408e6_dp]
      waypoints(45,:) = [ -1.590063024692466e6_dp, -0.248520311768411e6_dp]
      waypoints(46,:) = [ -1.591964534323602e6_dp, -0.253144623730497e6_dp]
      waypoints(47,:) = [ -1.593581606297590e6_dp, -0.257875910894784e6_dp]
      waypoints(48,:) = [ -1.594960641571018e6_dp, -0.262681975991580e6_dp]
      waypoints(49,:) = [ -1.596396222333623e6_dp, -0.267471454864489e6_dp]
      waypoints(50,:) = [ -1.597877923431084e6_dp, -0.272246867083112e6_dp]
      waypoints(51,:) = [ -1.599403797726561e6_dp, -0.277008348746881e6_dp]
      waypoints(52,:) = [ -1.600961667134800e6_dp, -0.281759458396939e6_dp]
      waypoints(53,:) = [ -1.602563990685723e6_dp, -0.286495760671726e6_dp]
      waypoints(54,:) = [ -1.604227207405997e6_dp, -0.291211024200035e6_dp]
      waypoints(55,:) = [ -1.605962531725658e6_dp, -0.295900229836984e6_dp]
      waypoints(56,:) = [ -1.607790896495252e6_dp, -0.300553946874934e6_dp]
      waypoints(57,:) = [ -1.609721065547457e6_dp, -0.305166370037428e6_dp]
      waypoints(58,:) = [ -1.611801062234093e6_dp, -0.309713194619879e6_dp]
      waypoints(59,:) = [ -1.613991890870033e6_dp, -0.314207665655268e6_dp]
      waypoints(60,:) = [ -1.616297252580522e6_dp, -0.318644478405860e6_dp]
      waypoints(61,:) = [ -1.618714397099108e6_dp, -0.323021396542662e6_dp]
      waypoints(62,:) = [ -1.621028944400364e6_dp, -0.327453424384434e6_dp]
      waypoints(63,:) = [ -1.623343491701620e6_dp, -0.331885452226207e6_dp]
      waypoints(64,:) = [ -1.625658039002877e6_dp, -0.336317480067980e6_dp]
      waypoints(65,:) = [ -1.627972586304133e6_dp, -0.340749507909753e6_dp]
      waypoints(66,:) = [ -1.630287133605389e6_dp, -0.345181535751525e6_dp]
      waypoints(67,:) = [ -1.632601680906645e6_dp, -0.349613563593298e6_dp]
      waypoints(68,:) = [ -1.634916228207901e6_dp, -0.354045591435071e6_dp]
      waypoints(69,:) = [ -1.637230775509157e6_dp, -0.358477619276843e6_dp]
      waypoints(70,:) = [ -1.639545322810413e6_dp, -0.362909647118616e6_dp]
      waypoints(71,:) = [ -1.641859870111669e6_dp, -0.367341674960389e6_dp]
      waypoints(72,:) = [ -1.644174417412925e6_dp, -0.371773702802161e6_dp]

    case('PineIsland_groundingline')

      allocate(waypoints(2,2))
      waypoints(1,:) = [-1.605e6_dp, -2.45e5_dp]
      waypoints(2,:) = [-1.570e6_dp, -2.55e5_dp]

    case('Thwaites_centralflowline')

      allocate(waypoints(76,2))

      waypoints( 1,:) = [-1.262096857411107e6_dp, -0.438496332811367e6_dp]
      waypoints( 2,:) = [-1.267094568089579e6_dp, -0.438647620531574e6_dp]
      waypoints( 3,:) = [-1.272062926853136e6_dp, -0.439209235343157e6_dp]
      waypoints( 4,:) = [-1.277000000000000e6_dp, -0.440000000000000e6_dp]
      waypoints( 5,:) = [-1.281937073146864e6_dp, -0.440790764656843e6_dp]
      waypoints( 6,:) = [-1.286846151547105e6_dp, -0.441739948127141e6_dp]
      waypoints( 7,:) = [-1.291762313729901e6_dp, -0.442651731759614e6_dp]
      waypoints( 8,:) = [-1.296616077151230e6_dp, -0.443852140294145e6_dp]
      waypoints( 9,:) = [-1.301417963270147e6_dp, -0.445245657319859e6_dp]
      waypoints(10,:) = [-1.306169781465008e6_dp, -0.446801364182665e6_dp]
      waypoints(11,:) = [-1.310844962489247e6_dp, -0.448574125416052e6_dp]
      waypoints(12,:) = [-1.315510104371187e6_dp, -0.450373139374133e6_dp]
      waypoints(13,:) = [-1.320269834067099e6_dp, -0.451904469913843e6_dp]
      waypoints(14,:) = [-1.325202382474215e6_dp, -0.452722984551133e6_dp]
      waypoints(15,:) = [-1.330202374209151e6_dp, -0.452732075787285e6_dp]
      waypoints(16,:) = [-1.335194857178497e6_dp, -0.452458006816433e6_dp]
      waypoints(17,:) = [-1.340194796803991e6_dp, -0.452433435666387e6_dp]
      waypoints(18,:) = [-1.345194083700367e6_dp, -0.452517878118564e6_dp]
      waypoints(19,:) = [-1.350194047295371e6_dp, -0.452536958176968e6_dp]
      waypoints(20,:) = [-1.355193800248603e6_dp, -0.452586661362422e6_dp]
      waypoints(21,:) = [-1.360193800237820e6_dp, -0.452586332983246e6_dp]
      waypoints(22,:) = [-1.365191611720165e6_dp, -0.452734253187886e6_dp]
      waypoints(23,:) = [-1.370173725508693e6_dp, -0.453156795727950e6_dp]
      waypoints(24,:) = [-1.375094914783623e6_dp, -0.454041044629737e6_dp]
      waypoints(25,:) = [-1.380042523096452e6_dp, -0.454762967051344e6_dp]
      waypoints(26,:) = [-1.384985030458033e6_dp, -0.455519023253439e6_dp]
      waypoints(27,:) = [-1.389905381954831e6_dp, -0.456407922093552e6_dp]
      waypoints(28,:) = [-1.394758062859852e6_dp, -0.457612699265498e6_dp]
      waypoints(29,:) = [-1.399549605770857e6_dp, -0.459041375762121e6_dp]
      waypoints(30,:) = [-1.404412006691101e6_dp, -0.460206303781702e6_dp]
      waypoints(31,:) = [-1.409375125235378e6_dp, -0.460812483883206e6_dp]
      waypoints(32,:) = [-1.414363677683855e6_dp, -0.461150632473021e6_dp]
      waypoints(33,:) = [-1.419348305049395e6_dp, -0.461542410505466e6_dp]
      waypoints(34,:) = [-1.424325029013987e6_dp, -0.462024301136476e6_dp]
      waypoints(35,:) = [-1.429276934310823e6_dp, -0.462716134881800e6_dp]
      waypoints(36,:) = [-1.434182109406869e6_dp, -0.463685287747517e6_dp]
      waypoints(37,:) = [-1.439019008423905e6_dp, -0.464951939988269e6_dp]
      waypoints(38,:) = [-1.443849185719868e6_dp, -0.466243987700066e6_dp]
      waypoints(39,:) = [-1.448726931695150e6_dp, -0.467342893610787e6_dp]
      waypoints(40,:) = [-1.453607812480557e6_dp, -0.468427791190582e6_dp]
      waypoints(41,:) = [-1.458543806615878e6_dp, -0.469225263381994e6_dp]
      waypoints(42,:) = [-1.463442227989308e6_dp, -0.470227993678896e6_dp]
      waypoints(43,:) = [-1.468340985966536e6_dp, -0.471229078231043e6_dp]
      waypoints(44,:) = [-1.473300063257687e6_dp, -0.471867476555276e6_dp]
      waypoints(45,:) = [-1.478267381997524e6_dp, -0.471296736214254e6_dp]
      waypoints(46,:) = [-1.483195430158483e6_dp, -0.470451549061452e6_dp]
      waypoints(47,:) = [-1.488155780636239e6_dp, -0.469823119958193e6_dp]
      waypoints(48,:) = [-1.493126229684642e6_dp, -0.469280315889273e6_dp]
      waypoints(49,:) = [-1.498125332401152e6_dp, -0.469185595090182e6_dp]
      waypoints(50,:) = [-1.503125181982443e6_dp, -0.469224378649324e6_dp]
      waypoints(51,:) = [-1.508124245996874e6_dp, -0.469321120468730e6_dp]
      waypoints(52,:) = [-1.513110996404370e6_dp, -0.468957361794142e6_dp]
      waypoints(53,:) = [-1.518025699462808e6_dp, -0.468037745652817e6_dp]
      waypoints(54,:) = [-1.522926373830334e6_dp, -0.467046085053888e6_dp]
      waypoints(55,:) = [-1.527888330802762e6_dp, -0.466430469099107e6_dp]
      waypoints(56,:) = [-1.532887623379006e6_dp, -0.466514574610874e6_dp]
      waypoints(57,:) = [-1.537747238844029e6_dp, -0.467691068355910e6_dp]
      waypoints(58,:) = [-1.542254522273829e6_dp, -0.469855414927846e6_dp]
      waypoints(59,:) = [-1.546644235647556e6_dp, -0.472249243763942e6_dp]
      waypoints(60,:) = [-1.551060793902853e6_dp, -0.474593174901490e6_dp]
      waypoints(61,:) = [-1.555454997509384e6_dp, -0.476978751281392e6_dp]
      waypoints(62,:) = [-1.559828912059237e6_dp, -0.479401326669401e6_dp]
      waypoints(63,:) = [-1.564088885569191e6_dp, -0.482019078700150e6_dp]
      waypoints(64,:) = [-1.568262054133600e6_dp, -0.484773105594177e6_dp]
      waypoints(65,:) = [-1.572310964033627e6_dp, -0.487706760076468e6_dp]
      waypoints(66,:) = [-1.576326433354400e6_dp, -0.490686022759403e6_dp]
      waypoints(67,:) = [-1.580341902675174e6_dp, -0.493665285442338e6_dp]
      waypoints(68,:) = [-1.584357371995947e6_dp, -0.496644548125273e6_dp]
      waypoints(69,:) = [-1.588372841316720e6_dp, -0.499623810808208e6_dp]
      waypoints(70,:) = [-1.592388310637493e6_dp, -0.502603073491143e6_dp]
      waypoints(71,:) = [-1.596403779958266e6_dp, -0.505582336174078e6_dp]
      waypoints(72,:) = [-1.600419249279039e6_dp, -0.508561598857013e6_dp]
      waypoints(73,:) = [-1.604434718599813e6_dp, -0.511540861539948e6_dp]
      waypoints(74,:) = [-1.608450187920586e6_dp, -0.514520124222883e6_dp]
      waypoints(75,:) = [-1.612465657241359e6_dp, -0.517499386905818e6_dp]
      waypoints(76,:) = [-1.616481126562132e6_dp, -0.520478649588753e6_dp]

    case('Thwaites_groundingline')

      allocate(waypoints(2,2))
      waypoints(1,:) = [-1.520e6_dp, -4.00e5_dp]
      waypoints(2,:) = [-1.495e6_dp, -5.10e5_dp]

    ! Greenland

    case('Jakobshavn_centralflowline')
      call crash('transect "' // trim(name) // '" not implemented yet!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transect_waypoints_hardcoded

  subroutine initialise_transect_waypoints_from_file( filename, waypoints)

    ! In/output variables
    character(len=*),                      intent(in   ) :: filename
    real(dp), dimension(:,:), allocatable, intent(  out) :: waypoints

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_transect_waypoints_from_file'
    integer                        :: i, n_wp, ios, ierr
    real(dp), dimension(2)         :: wp

    ! Add routine to path
    call init_routine( routine_name)

    ! Let the primary read the file, then broadcast the data to the processes

    ! Determine number of waypoints
    if (par%primary) then
      n_wp = 0
      open( unit = 1337, file = filename, action = 'read')
      do while (.true.)
        read( unit = 1337, fmt = *, iostat = ios) wp(1), wp(2)
        if (ios /= 0) exit
        n_wp = n_wp + 1
      end do
      close( unit  = 1337)
      ! Safety
      if (n_wp < 2) call crash('invalid transect in file "' // trim( filename) // &
        '" - need at least two waypoints')
    end if
    call MPI_BCAST( n_wp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    allocate( waypoints( n_wp,2))

    ! Read waypoints
    if (par%primary) then
      open( unit = 1337, file = filename, action = 'read')
      do i = 1, n_wp
        read( unit = 1337, fmt = *, iostat = ios) waypoints( i,1), waypoints( i,2)
      end do
      close( unit  = 1337)
    end if
    call MPI_BCAST( waypoints(:,:), n_wp*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transect_waypoints_from_file

  subroutine calc_transect_vertices_from_waypoints( transect, waypoints, dx)
    !< Generate transect vertices spaced dx apart along waypoints

    ! In/output variables:
    type(type_transect),      intent(inout) :: transect
    real(dp), dimension(:,:), intent(in   ) :: waypoints
    real(dp),                 intent(in   ) :: dx

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_transect_vertices_from_waypoints'
    real(dp), dimension(size(waypoints,1)) :: wp_dist_along_transect
    integer                                :: i_wp, n_wp
    integer                                :: i
    real(dp)                               :: dist_along_transect, dist_between_waypoints, w

    ! Add routine to path
    call init_routine( routine_name)

    n_wp = size( waypoints,1)

    ! Calculate distance along transect for all the waypoints
    wp_dist_along_transect(1) = 0._dp
    do i_wp = 2, n_wp
      wp_dist_along_transect( i_wp) = wp_dist_along_transect( i_wp-1) + &
        norm2( waypoints( i_wp,:) - waypoints( i_wp-1,:))
    end do

    ! Calculate number of transect vertices
    transect%nV = floor( wp_dist_along_transect( n_wp) / dx)
    allocate( transect%V( transect%nV,2))

    ! Create transect vertices
    do i = 1, transect%nV

      dist_along_transect = real(i-1,dp) * dx
      i_wp = 2
      do while (wp_dist_along_transect( i_wp) < dist_along_transect)
        i_wp = i_wp + 1
      end do

      dist_between_waypoints = dist_along_transect - wp_dist_along_transect( i_wp-1)
      w = dist_between_waypoints / (wp_dist_along_transect( i_wp) - wp_dist_along_transect( i_wp-1))
      transect%V( i,:) = (1._dp - w) * waypoints( i_wp-1,:) + w * waypoints( i_wp,:)

    end do

    ! Parallelisation
    call partition_list( transect%nV, par%i, par%n, transect%vi1, transect%vi2)
    transect%nV_loc = transect%vi2 + 1 - transect%vi1

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_transect_vertices_from_waypoints

  subroutine calc_velocity_weights( transect)
    ! Calculate the weights needed to calculate parallel/orthogonal velocity components

    ! In/output variables:
    type(type_transect), intent(inout) :: transect

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_transect_vertices_from_waypoints'
    integer                        :: i
    real(dp), dimension(2)         :: d_par, d_ort

    ! Add routine to path
    call init_routine( routine_name)

    allocate( transect%wu_u_par( transect%vi1: transect%vi2), source = 0._dp)
    allocate( transect%wv_u_par( transect%vi1: transect%vi2), source = 0._dp)
    allocate( transect%wu_u_ort( transect%vi1: transect%vi2), source = 0._dp)
    allocate( transect%wv_u_ort( transect%vi1: transect%vi2), source = 0._dp)

    do i = transect%vi1, transect%vi2

      ! Calculate unit vector in the along-transect direction d
      if (i == 1) then
        d_par = transect%V( 2,:) - transect%V( 1,:)
      elseif (i == transect%nV) then
        d_par = transect%V( transect%nV,:) - transect%V( transect%nV-1,:)
      else
        d_par = transect%V( i+1,:) - transect%V( i-1,:)
      end if
      d_par = d_par / norm2( d_par)

      ! Calculate unit vector orthogonal to this (i.e. rotated 90 degrees clockwise)
      d_ort = [d_par(2), -d_par(1)]

      transect%wu_u_par( i) = d_par(1)
      transect%wv_u_par( i) = d_par(2)
      transect%wu_u_ort( i) = d_ort(1)
      transect%wv_u_ort( i) = d_ort(2)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_velocity_weights

  subroutine create_transect_netcdf_output_file( transect)

    ! In/output variables:
    type(type_transect), intent(inout) :: transect

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_transect_netcdf_output_file'
    character(len=1024)            :: filename
    integer                        :: ncid
    integer                        :: n, two, z, t

    ! Add routine to path
    call init_routine( routine_name)

    filename = trim( C%output_dir) // 'transect_' // trim( transect%name) // '.nc'
    transect%nc%filename = filename

    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Basic dimensions
    call create_dimension( filename, ncid, 'n'   , transect%nV,    transect%nc%id_dim_n)
    call create_dimension( filename, ncid, 'two' , 2,              transect%nc%id_dim_two)
    call create_dimension( filename, ncid, 'zeta', transect%nz,    transect%nc%id_dim_zeta)
    call create_dimension( filename, ncid, 'time', NF90_UNLIMITED, transect%nc%id_dim_time)

    ! Shorthand for dimension IDs for more readable code
    n   = transect%nc%id_dim_n
    two = transect%nc%id_dim_two
    z   = transect%nc%id_dim_zeta
    t   = transect%nc%id_dim_time

    ! Transect geometry
    call create_variable(    filename, ncid, 'V',    NF90_DOUBLE, [n, two], transect%nc%id_var_V)
    call add_attribute_char( filename, ncid, transect%nc%id_var_V, 'long_name', 'Transect vertex coordinates')
    call add_attribute_char( filename, ncid, transect%nc%id_var_V, 'units'    , 'm')
    call write_var_primary( filename, ncid, transect%nc%id_var_V   , transect%V)

    call create_variable(    filename, ncid, 'zeta', NF90_DOUBLE, [z], transect%nc%id_var_zeta)
    call add_attribute_char( filename, ncid, transect%nc%id_var_zeta, 'long_name', 'Scaled vertical coordinate')
    call add_attribute_char( filename, ncid, transect%nc%id_var_zeta, 'units', '0-1')
    call add_attribute_char( filename, ncid, transect%nc%id_var_zeta, 'transformation', 'zeta = (h - z) / H; zeta = 0 at the ice surface; zeta = 1 at the ice base')
    call write_var_primary( filename, ncid, transect%nc%id_var_zeta, transect%zeta)

    call create_variable(    filename, ncid, 'time', NF90_DOUBLE, [t], transect%nc%id_var_time)
    call add_attribute_char( filename, ncid, transect%nc%id_var_time, 'long_name', 'Time')
    call add_attribute_char( filename, ncid, transect%nc%id_var_time, 'units', 'years')

    ! Ice model variables
    call create_variable(  filename, ncid, 'Hi', NF90_DOUBLE, [n, t], transect%nc%id_var_Hi)
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hi, 'long_name', 'Ice thickness')
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hi, 'units', 'm')

    call create_variable(  filename, ncid, 'Hb', NF90_DOUBLE, [n, t], transect%nc%id_var_Hb)
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hb, 'long_name', 'Bedrock elevation')
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hb, 'units', 'm w.r.t. PD sea level')

    call create_variable(  filename, ncid, 'Hs', NF90_DOUBLE, [n, t], transect%nc%id_var_Hs)
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hs, 'long_name', 'Surface elevation')
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hs, 'units', 'm w.r.t. PD sea level')

    call create_variable(  filename, ncid, 'Hib', NF90_DOUBLE, [n, t], transect%nc%id_var_Hib)
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hib, 'long_name', 'Ice base elevation')
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hib, 'units', 'm w.r.t. PD sea level')

    call create_variable(  filename, ncid, 'SL', NF90_DOUBLE, [n, t], transect%nc%id_var_SL)
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hs, 'long_name', 'Geoid elevation')
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hs, 'units', 'm w.r.t. PD sea level')

    call create_variable(  filename, ncid, 'Hi_eff', NF90_DOUBLE, [n, t], transect%nc%id_var_Hi)
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hi, 'long_name', 'Effective ice thickness')
    call add_attribute_char( filename, ncid, transect%nc%id_var_Hi, 'units', 'm')

    call create_variable(  filename, ncid, 'Ti', NF90_DOUBLE, [n, z, t], transect%nc%id_var_Ti)
    call add_attribute_char( filename, ncid, transect%nc%id_var_Ti, 'long_name', 'Englacial temperature')
    call add_attribute_char( filename, ncid, transect%nc%id_var_Ti, 'units', 'K')

    call create_variable(  filename, ncid, 'u_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_u_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_3D, 'long_name', '3-D ice velocity in x-direction')
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_3D, 'units'    , 'm yr^-1')

    call create_variable(  filename, ncid, 'v_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_v_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_v_3D, 'long_name', '3-D ice velocity in y-direction')
    call add_attribute_char( filename, ncid, transect%nc%id_var_v_3D, 'units'    , 'm yr^-1')

    call create_variable(  filename, ncid, 'w_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_w_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_w_3D, 'long_name', '3-D ice velocity in z-direction')
    call add_attribute_char( filename, ncid, transect%nc%id_var_w_3D, 'units'    , 'm yr^-1')

    call create_variable(  filename, ncid, 'u_par_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_u_par_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_par_3D, 'long_name', '3-D ice velocity in along-transect direction')
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_par_3D, 'units'    , 'm yr^-1')

    call create_variable(  filename, ncid, 'u_ort_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_u_ort_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_ort_3D, 'long_name', '3-D ice velocity in across-transect direction')
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_ort_3D, 'units'    , 'm yr^-1')

    call create_variable(  filename, ncid, 'du_dx_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_du_dx_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_du_dx_3D, 'long_name', '3-D xx strain rate')
    call add_attribute_char( filename, ncid, transect%nc%id_var_du_dx_3D, 'units'    , 'yr^-1')

    call create_variable(  filename, ncid, 'dv_dy_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_dv_dy_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_dv_dy_3D, 'long_name', '3-D yy strain rate')
    call add_attribute_char( filename, ncid, transect%nc%id_var_dv_dy_3D, 'units'    , 'yr^-1')

    call create_variable(  filename, ncid, 'dw_dz_3D', NF90_DOUBLE, [n, z, t], transect%nc%id_var_dw_dz_3D)
    call add_attribute_char( filename, ncid, transect%nc%id_var_dw_dz_3D, 'long_name', '3-D zz strain rate')
    call add_attribute_char( filename, ncid, transect%nc%id_var_dw_dz_3D, 'units'    , 'yr^-1')

    ! Integrated quantities
    call create_variable(  filename, ncid, 'ice_mass_flux', NF90_DOUBLE, [t], transect%nc%id_var_ice_mass_flux)
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'long_name', 'Ice mass flux across transect')
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'units'    , 'kg yr^-1')

    call create_variable(  filename, ncid, 'grounding_line_distance_from_start', NF90_DOUBLE, [t], transect%nc%id_var_GL_dist_from_start)
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'long_name', 'Distance from start of transect to grounding line')
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'units'    , 'm')

    call create_variable(  filename, ncid, 'grounding_line_distance_from_end', NF90_DOUBLE, [t], transect%nc%id_var_GL_dist_from_end)
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'long_name', 'Distance from end of transect to grounding line')
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'units'    , 'm')

    call create_variable(  filename, ncid, 'calving_front_distance_from_start', NF90_DOUBLE, [t], transect%nc%id_var_CF_dist_from_start)
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'long_name', 'Distance from start of transect to calving front')
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'units'    , 'm')

    call create_variable(  filename, ncid, 'calving_front_distance_from_end', NF90_DOUBLE, [t], transect%nc%id_var_CF_dist_from_end)
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'long_name', 'Distance from end of transect to calving front')
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'units'    , 'm')

    call create_variable(  filename, ncid, 'Hi_eff_CF_from_start', NF90_DOUBLE, [t], transect%nc%id_var_CF_Hi_eff_from_start)
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'long_name', 'Effective ice thickness at calving front from start of transect')
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'units'    , 'm')

    call create_variable(  filename, ncid, 'Hi_eff_CF_from_end', NF90_DOUBLE, [t], transect%nc%id_var_CF_Hi_eff_from_end)
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'long_name', 'Effective ice thickness at calving front from end of transect')
    call add_attribute_char( filename, ncid, transect%nc%id_var_ice_mass_flux, 'units'    , 'm')

    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_transect_netcdf_output_file

  subroutine write_to_transect_netcdf_output_files( region)

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_transect_netcdf_output_files'
    integer                        :: it

    ! Add routine to path
    call init_routine( routine_name)

    do it = 1, size( region%transects)
      call write_to_transect_netcdf_output_file( region%mesh, region%ice, &
        region%transects( it), region%time)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_transect_netcdf_output_files

  subroutine write_to_transect_netcdf_output_file( mesh, ice, transect, time)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice
    type(type_transect),  intent(in   ) :: transect
    real(dp),             intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter                      :: routine_name = 'write_to_transect_netcdf_output_file'
    character(len=1024)                                 :: filename
    integer                                             :: ncid, k
    real(dp), dimension(transect%vi1:transect%vi2     ) :: tHi_partial
    real(dp), dimension(transect%vi1:transect%vi2     ) :: tHb_partial
    real(dp), dimension(transect%vi1:transect%vi2     ) :: tHs_partial
    real(dp), dimension(transect%vi1:transect%vi2     ) :: tHib_partial
    real(dp), dimension(transect%vi1:transect%vi2     ) :: tSL_partial
    real(dp), dimension(transect%vi1:transect%vi2     ) :: tHi_eff_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tTi_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tu_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tv_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tw_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tu_par_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tu_ort_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tdu_dx_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tdv_dy_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tdw_dz_partial
    ! real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tage_partial
    real(dp)                                            :: ice_mass_flux
    real(dp), dimension(C%nz)                           :: u_ort_prof
    integer                                             :: i, ierr
    real(dp)                                            :: GL_dist_from_start
    real(dp)                                            :: GL_dist_from_end
    real(dp)                                            :: CF_dist_from_start
    real(dp)                                            :: CF_dist_from_end
    real(dp)                                            :: Hi_eff_CF_from_start
    real(dp)                                            :: Hi_eff_CF_from_end

    ! Add routine to path
    call init_routine( routine_name)

    filename = transect%nc%filename

    if (par%primary) write(0,*) '  Writing to transect output file "', &
      colour_string( trim( filename), 'light blue'), '"...'

    ! Map ice model data to transect
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hi,       tHi_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hb,       tHb_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hs,       tHs_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hib,      tHib_partial,    'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%SL,       tSL_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hi_eff,   tHi_eff_partial, 'nearest_neighbour')
    call map_from_mesh_vertices_to_transect_3D ( mesh, transect, ice%Ti,       tTi_partial,     'trilin')
    call map_from_mesh_triangles_to_transect_3D( mesh, transect, ice%u_3D_b,   tu_partial)
    call map_from_mesh_triangles_to_transect_3D( mesh, transect, ice%v_3D_b,   tv_partial)
    call map_from_mesh_vertices_to_transect_3D ( mesh, transect, ice%w_3D,     tw_partial,      'trilin')
    call map_from_mesh_vertices_to_transect_3D(  mesh, transect, ice%du_dx_3D, tdu_dx_partial,  'trilin')
    call map_from_mesh_vertices_to_transect_3D(  mesh, transect, ice%dv_dy_3D, tdv_dy_partial,  'trilin')
    call map_from_mesh_vertices_to_transect_3D(  mesh, transect, ice%dw_dz_3D, tdw_dz_partial,  'trilin')

    ! Calculate parallel/orthogonal velocity components
    do k = 1, transect%nz
      tu_par_partial(:,k) = transect%wu_u_par * tu_partial(:,k) + transect%wv_u_par * tv_partial(:,k)
      tu_ort_partial(:,k) = transect%wu_u_ort * tu_partial(:,k) + transect%wv_u_ort * tv_partial(:,k)
    end do

    ! Integrate ice mass flux across transect
    ice_mass_flux = 0._dp
    do i = transect%vi1, transect%vi2
      u_ort_prof = tu_ort_partial( i,:)
      ice_mass_flux =  ice_mass_flux + &
        tHi_partial( i) * transect%dx * vertical_average( transect%zeta, u_ort_prof) * ice_density
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, ice_mass_flux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate distance of grounding line and calving front along the transect
    call calc_GL_CF_along_transect( transect, tHi_partial, tHb_partial, tSL_partial, tHi_eff_partial, &
      GL_dist_from_start, GL_dist_from_end, CF_dist_from_start, CF_dist_from_end, &
      Hi_eff_CF_from_start, Hi_eff_CF_from_end)

    ! Write transect data to output
    call open_existing_netcdf_file_for_writing( filename, ncid)

    call write_time_to_file( filename, ncid, time)

    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hi'      , tHi_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hb'      , tHb_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hs'      , tHs_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hib'     , tHib_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'SL'      , tSL_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'Ti'      , tTi_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hi_eff'  , tHi_eff_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'u_3D'    , tu_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'v_3D'    , tv_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'w_3D'    , tw_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'u_par_3D', tu_par_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'u_ort_3D', tu_ort_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'du_dx_3D', tdu_dx_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'dv_dy_3D', tdv_dy_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'dw_dz_3D', tdw_dz_partial)

    call write_to_field_multopt_dp_0D( filename, ncid, 'ice_mass_flux', ice_mass_flux)
    call write_to_field_multopt_dp_0D( filename, ncid, 'grounding_line_distance_from_start', GL_dist_from_start)
    call write_to_field_multopt_dp_0D( filename, ncid, 'grounding_line_distance_from_end',   GL_dist_from_end)
    call write_to_field_multopt_dp_0D( filename, ncid, 'calving_front_distance_from_start',  CF_dist_from_start)
    call write_to_field_multopt_dp_0D( filename, ncid, 'calving_front_distance_from_end',    CF_dist_from_end)
    call write_to_field_multopt_dp_0D( filename, ncid, 'Hi_eff_CF_from_start',               Hi_eff_CF_from_start)
    call write_to_field_multopt_dp_0D( filename, ncid, 'Hi_eff_CF_from_end',                 Hi_eff_CF_from_end)

    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_transect_netcdf_output_file

  subroutine calc_GL_CF_along_transect( transect, tHi_partial, tHb_partial, tSL_partial, tHi_eff_partial, &
    GL_dist_from_start, GL_dist_from_end, CF_dist_from_start, CF_dist_from_end, &
    Hi_eff_CF_from_start, Hi_eff_CF_from_end)

    ! In/output variables:
    type(type_transect),                            intent(in   ) :: transect
    real(dp), dimension(transect%vi1:transect%vi2), intent(in   ) :: tHi_partial
    real(dp), dimension(transect%vi1:transect%vi2), intent(in   ) :: tHb_partial
    real(dp), dimension(transect%vi1:transect%vi2), intent(in   ) :: tSL_partial
    real(dp), dimension(transect%vi1:transect%vi2), intent(in   ) :: tHi_eff_partial
    real(dp),                                       intent(  out) :: GL_dist_from_start
    real(dp),                                       intent(  out) :: GL_dist_from_end
    real(dp),                                       intent(  out) :: CF_dist_from_start
    real(dp),                                       intent(  out) :: CF_dist_from_end
    real(dp),                                       intent(  out) :: Hi_eff_CF_from_start
    real(dp),                                       intent(  out) :: Hi_eff_CF_from_end

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'calc_GL_CF_along_transect'
    real(dp), dimension(transect%vi1:transect%vi2) :: tTAF_partial
    real(dp), dimension(transect%nV)               :: tTAF_tot
    integer                                        :: i, il, ir
    real(dp)                                       :: TAF_left, TAF_right, dist_left, dist_right
    real(dp), dimension(transect%nV)               :: tHi_eff_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Locate the grounding line from either end of the transect
    do i = transect%vi1, transect%vi2
      tTAF_partial( i) = thickness_above_floatation( tHi_partial( i), tHb_partial( i), tSL_partial( i))
    end do
    call gather_to_all( tTAF_partial, tTAF_tot)

    GL_dist_from_start = 0._dp
    do il = 1, transect%nV-1
      ir = il+1
      if (tTAF_tot(il) * tTAF_tot( ir) < 0._dp) then
        TAF_left   = tTAF_tot( il)
        TAF_right  = tTAF_tot( ir)
        dist_left  = real(il,dp) * transect%dx
        dist_right = real(ir,dp) * transect%dx
        GL_dist_from_start = linint_points( dist_left, dist_right, TAF_left, TAF_right, 0._dp)
        exit
      end if
    end do

    GL_dist_from_end = 0._dp
    do il = transect%nV-1, 1, -1
      ir = il+1
      if (tTAF_tot(il) * tTAF_tot( ir) < 0._dp) then
        TAF_left   = tTAF_tot( il)
        TAF_right  = tTAF_tot( ir)
        dist_left  = real(transect%nV + 1 - il,dp) * transect%dx
        dist_right = real(transect%nV + 1 - ir,dp) * transect%dx
        GL_dist_from_end = linint_points( dist_left, dist_right, TAF_left, TAF_right, 0._dp)
        exit
      end if
    end do

    ! Locate the calving front from either end of the transect
    call gather_to_all( tHi_eff_partial, tHi_eff_tot)

    CF_dist_from_start   = 0._dp
    Hi_eff_CF_from_start = 0._dp
    do il = 1, transect%nV-1
      ir = il+1
      if ((tHi_eff_tot(il) > 0.1_dp .and. tHi_eff_tot(ir) < 0.1_dp) .or. &
          (tHi_eff_tot(il) < 0.1_dp .and. tHi_eff_tot(ir) > 0.1_dp)) then
        dist_left  = real(il,dp) * transect%dx
        dist_right = real(ir,dp) * transect%dx
        CF_dist_from_start = (dist_left + dist_right) / 2._dp
        Hi_eff_CF_from_start = max( tHi_eff_tot(il), tHi_eff_tot(ir))
        exit
      end if
    end do

    CF_dist_from_end   = 0._dp
    Hi_eff_CF_from_end = 0._dp
    do il = transect%nV-1, 1, -1
      ir = il+1
      if ((tHi_eff_tot(il) > 0.1_dp .and. tHi_eff_tot(ir) < 0.1_dp) .or. &
          (tHi_eff_tot(il) < 0.1_dp .and. tHi_eff_tot(ir) > 0.1_dp)) then
        dist_left  = real(transect%nV + 1 - il,dp) * transect%dx
        dist_right = real(transect%nV + 1 - ir,dp) * transect%dx
        CF_dist_from_end = (dist_left + dist_right) / 2._dp
        Hi_eff_CF_from_end = max( tHi_eff_tot(il), tHi_eff_tot(ir))
        exit
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_GL_CF_along_transect

end module transects_main
