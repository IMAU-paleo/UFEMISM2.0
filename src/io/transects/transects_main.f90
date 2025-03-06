module transects_main

  use mpi
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
  use remapping_main, only: map_from_mesh_vertices_to_transect_2D, map_from_mesh_vertices_to_transect_3D, &
    map_from_mesh_triangles_to_transect_2D, map_from_mesh_triangles_to_transect_3D
  use parameters, only: ice_density
  use mesh_zeta, only: vertical_average
  use ice_geometry_basics, only: thickness_above_floatation
  use mpi_distributed_memory, only: gather_to_all
  use interpolation, only: linint_points

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

    if (par%master) write(0,*) '  Initialising output transect ', &
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
      call crash('transect "' // trim(name) // '" not implemented yet!')

    case('Thwaites_centralflowline')
      call crash('transect "' // trim(name) // '" not implemented yet!')

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

    ! Let the master read the file, then broadcast the data to the processes

    ! Determine number of waypoints
    if (par%master) then
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
    if (par%master) then
      open( unit = 1337, file = filename, action = 'read')
      do i = 1, n_wp
        read( unit = 1337, fmt = *, iostat = ios) waypoints( i,1), waypoints( i,2)
      end do
      close( unit  = 1337)
    end if
    call MPI_BCAST( waypoints, n_wp*2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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
    call write_var_master( filename, ncid, transect%nc%id_var_V   , transect%V)

    call create_variable(    filename, ncid, 'zeta', NF90_DOUBLE, [z], transect%nc%id_var_zeta)
    call add_attribute_char( filename, ncid, transect%nc%id_var_zeta, 'long_name', 'Scaled vertical coordinate')
    call add_attribute_char( filename, ncid, transect%nc%id_var_zeta, 'units', '0-1')
    call add_attribute_char( filename, ncid, transect%nc%id_var_zeta, 'transformation', 'zeta = (h - z) / H; zeta = 0 at the ice surface; zeta = 1 at the ice base')
    call write_var_master( filename, ncid, transect%nc%id_var_zeta, transect%zeta)

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

    call create_variable(  filename, ncid, 'u_par', NF90_DOUBLE, [n, z, t], transect%nc%id_var_u_par)
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_par, 'long_name', '3-D ice velocity in along-transect direction')
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_par, 'units'    , 'm yr^-1')

    call create_variable(  filename, ncid, 'u_ort', NF90_DOUBLE, [n, z, t], transect%nc%id_var_u_ort)
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_ort, 'long_name', '3-D ice velocity in across-transect direction')
    call add_attribute_char( filename, ncid, transect%nc%id_var_u_ort, 'units'    , 'm yr^-1')

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
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tu_par_partial
    real(dp), dimension(transect%vi1:transect%vi2,C%nz) :: tu_ort_partial
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

    if (par%master) write(0,*) '  Writing to transect output file "', &
      colour_string( trim( filename), 'light blue'), '"...'

    ! Map ice model data to transect
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hi,     tHi_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hb,     tHb_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hs,     tHs_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hib,    tHib_partial,    'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%SL,     tSL_partial,     'trilin')
    call map_from_mesh_vertices_to_transect_2D ( mesh, transect, ice%Hi_eff, tHi_eff_partial, 'nearest_neighbour')
    call map_from_mesh_vertices_to_transect_3D ( mesh, transect, ice%Ti,     tTi_partial,     'trilin')
    call map_from_mesh_triangles_to_transect_3D( mesh, transect, ice%u_3D_b, tu_partial)
    call map_from_mesh_triangles_to_transect_3D( mesh, transect, ice%v_3D_b, tv_partial)

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

    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hi'    , tHi_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hb'    , tHb_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hs'    , tHs_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hib'   , tHib_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'SL'    , tSL_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'Ti'    , tTi_partial)
    call write_to_field_multopt_transect_dp_2D( transect, filename, ncid, 'Hi_eff', tHi_eff_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'u_par' , tu_par_partial)
    call write_to_field_multopt_transect_dp_3D( transect, filename, ncid, 'u_ort' , tu_ort_partial)

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