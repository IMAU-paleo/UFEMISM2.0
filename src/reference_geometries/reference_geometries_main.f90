module reference_geometries_main

  ! Contains the routines for setting up the three "reference geometries":
  ! - refgeo_init:   initial, used to initialise the simulation
  ! - refgeo_PD:     present-day, used to calculate sea-level contribution, isotope change, and more
  ! - refgeo_GIA_eq: GIA equilibrium, used for the GIA model

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use reference_geometry_types, only: type_reference_geometry
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid, setup_square_grid
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary
  use ice_geometry_basics, only: ice_surface_elevation
  use netcdf_io_main
  use remapping_main, only: map_from_xy_grid_to_mesh_2D, map_from_mesh_to_mesh_2D
  use preprocess_geometry, only: smooth_model_geometry, remove_Lake_Vostok, remove_Ellesmere, remove_tiny_islands
  use idealised_geometries, only: calc_idealised_geometry

  implicit none

  private

  public :: initialise_reference_geometries_raw, initialise_reference_geometries_on_model_mesh, &
    initialise_reference_geometry_raw_from_file, remap_reference_geometry_to_mesh, &
    reallocate_reference_geometry_on_mesh

contains

  ! Initialise reference geometries on the model mesh
  ! =================================================

  subroutine initialise_reference_geometries_on_model_mesh( region_name, mesh, refgeo_init, refgeo_PD, refgeo_GIAeq)
    !< Initialise all reference geometries on the model mesh

    ! In/output variables:
    character(len=3),              intent(in   ) :: region_name
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(inout) :: refgeo_init
    type(type_reference_geometry), intent(inout) :: refgeo_PD
    type(type_reference_geometry), intent(inout) :: refgeo_GIAeq

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_reference_geometries_on_model_mesh'
    character(len=1024)            :: choice_refgeo, choice_refgeo_idealised

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '     Mapping reference geometries to model mesh...'

    ! == Initial geometry
    ! ===================

    call reallocate_reference_geometry_on_mesh( mesh, refgeo_init)

    ! Get config choices for this model region
    select case (region_name)
      case default
        call crash('unknown region_name "' // trim( region_name))
      case ('NAM')
        choice_refgeo           = C%choice_refgeo_init_NAM
        choice_refgeo_idealised = C%choice_refgeo_init_idealised
      case ('EAS')
        choice_refgeo           = C%choice_refgeo_init_EAS
        choice_refgeo_idealised = C%choice_refgeo_init_idealised
      case ('GRL')
        choice_refgeo           = C%choice_refgeo_init_GRL
        choice_refgeo_idealised = C%choice_refgeo_init_idealised
      case ('ANT')
        choice_refgeo           = C%choice_refgeo_init_ANT
        choice_refgeo_idealised = C%choice_refgeo_init_idealised
    end select

    select case (choice_refgeo)
      case default
        call crash('unknown choice_refgeo "' // trim( choice_refgeo) // '"!')
      case ('read_from_file')
        ! For realistic geometries read from a file, remap them to the model mesh
        call remap_reference_geometry_to_mesh( mesh, refgeo_init)
      case ('idealised')
        ! For idealised geometries, calculate them directly on the model mesh
        call initialise_reference_geometry_idealised( mesh, choice_refgeo_idealised, refgeo_init)
    end select

    ! == Present-day geometry
    ! =======================

    call reallocate_reference_geometry_on_mesh( mesh, refgeo_PD)

    ! Get config choices for this model region
    select case (region_name)
      case default
        call crash('unknown region_name "' // trim( region_name) // '"!')
      case ('NAM')
        choice_refgeo           = C%choice_refgeo_PD_NAM
        choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      case ('EAS')
        choice_refgeo           = C%choice_refgeo_PD_EAS
        choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      case ('GRL')
        choice_refgeo           = C%choice_refgeo_PD_GRL
        choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      case ('ANT')
        choice_refgeo           = C%choice_refgeo_PD_ANT
        choice_refgeo_idealised = C%choice_refgeo_PD_idealised
    end select

    select case (choice_refgeo)
      case default
        call crash('unknown choice_refgeo "' // trim( choice_refgeo) // '"!')
      case ('read_from_file')
        ! For realistic geometries read from a file, remap them to the model mesh
        call remap_reference_geometry_to_mesh( mesh, refgeo_PD)
      case ('idealised')
        ! For idealised geometries, calculate them directly on the model mesh
        call initialise_reference_geometry_idealised( mesh, choice_refgeo_idealised, refgeo_PD)
    end select

    ! == GIA equilibrium geometry
    ! ===========================

    call reallocate_reference_geometry_on_mesh( mesh, refgeo_GIAeq)

    ! Get config choices for this model region
    select case (region_name)
      case default
        call crash('unknown region_name "' // trim( region_name) // '"!')
      case ('NAM')
        choice_refgeo           = C%choice_refgeo_GIAeq_NAM
        choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      case ('EAS')
        choice_refgeo           = C%choice_refgeo_GIAeq_EAS
        choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      case ('GRL')
        choice_refgeo           = C%choice_refgeo_GIAeq_GRL
        choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      case ('ANT')
        choice_refgeo           = C%choice_refgeo_GIAeq_ANT
        choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
    end select

    select case (choice_refgeo)
      case default
        call crash('unknown choice_refgeo "' // trim( choice_refgeo) // '"!')
      case ('read_from_file')
        ! For realistic geometries read from a file, remap them to the model mesh
        call remap_reference_geometry_to_mesh( mesh, refgeo_GIAeq)
      case ('idealised')
        ! For idealised geometries, calculate them directly on the model mesh
        call initialise_reference_geometry_idealised( mesh, choice_refgeo_idealised, refgeo_GIAeq)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometries_on_model_mesh

  subroutine remap_reference_geometry_to_mesh( mesh, refgeo)
    !< Remap reference geometry to the model mesh

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_reference_geometry_to_mesh'
    character(len=1024)            :: method
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_reference_geometry_on_mesh( mesh, refgeo)

    ! Determine if the initial geometry is provided gridded or meshed
    if (allocated( refgeo%grid_raw%x)) then
      ! Gridded

      ! Safety
      if (allocated( refgeo%mesh_raw%V)) call crash('found both grid and mesh in refgeo!')

      ! Remap data to the model mesh
      call map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, C%output_dir, refgeo%Hi_grid_raw, refgeo%Hi)
      call map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, C%output_dir, refgeo%Hb_grid_raw, refgeo%Hb)
      call map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, C%output_dir, refgeo%SL_grid_raw, refgeo%SL)

    elseif (allocated( refgeo%mesh_raw%V)) then
      ! Meshed

      ! Safety
      if (allocated( refgeo%grid_raw%x)) call crash('found both grid and mesh in refgeo!')

      ! Remap data to the model mesh
      method = '2nd_order_conservative'
      call map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, C%output_dir, refgeo%Hi_mesh_raw, refgeo%Hi, method)
      call map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, C%output_dir, refgeo%Hb_mesh_raw, refgeo%Hb, method)
      call map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, C%output_dir, refgeo%SL_mesh_raw, refgeo%SL, method)

    else
      call crash('no grid or mesh is found in refgeo!')
    end if

    ! don't remap Hs, but recalculate it after remapping Hi,Hb,SL
    do vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = ice_surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), refgeo%SL( vi))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_reference_geometry_to_mesh

  ! Initialise reference geometries on their raw input grid/mesh
  ! ============================================================

  subroutine initialise_reference_geometries_raw( region_name, refgeo_init, refgeo_PD, refgeo_GIAeq)
    !< Initialise all reference geometries on the raw grid/mesh

    ! In/output variables:
    character(len=3)             , intent(in   ) :: region_name
    type(type_reference_geometry), intent(  out) :: refgeo_init
    type(type_reference_geometry), intent(  out) :: refgeo_PD
    type(type_reference_geometry), intent(  out) :: refgeo_GIAeq

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_reference_geometries_raw'
    character(len=1024)            :: choice_refgeo
    character(len=1024)            :: choice_refgeo_idealised
    real(dp)                       :: dx_refgeo_idealised
    character(len=1024)            :: filename_refgeo
    real(dp)                       :: timeframe_refgeo

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initial geometry
    ! ===================

    ! Get config choices for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name))
    case ('NAM')
      choice_refgeo           = C%choice_refgeo_init_NAM
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_NAM
      timeframe_refgeo        = C%timeframe_refgeo_init_NAM
    case ('EAS')
      choice_refgeo           = C%choice_refgeo_init_EAS
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_EAS
      timeframe_refgeo        = C%timeframe_refgeo_init_EAS
    case ('GRL')
      choice_refgeo           = C%choice_refgeo_init_GRL
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_GRL
      timeframe_refgeo        = C%timeframe_refgeo_init_GRL
    case ('ANT')
      choice_refgeo           = C%choice_refgeo_init_ANT
      choice_refgeo_idealised = C%choice_refgeo_init_idealised
      dx_refgeo_idealised     = C%dx_refgeo_init_idealised
      filename_refgeo         = C%filename_refgeo_init_ANT
      timeframe_refgeo        = C%timeframe_refgeo_init_ANT
    end select

    ! Initialise reference geometry
    call initialise_reference_geometry_raw( region_name, 'initial', refgeo_init, &
      choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)

    ! == Present-day geometry
    ! =======================

    ! Get config choices for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name))
    case ('NAM')
      choice_refgeo           = C%choice_refgeo_PD_NAM
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_NAM
      timeframe_refgeo        = C%timeframe_refgeo_PD_NAM
    case ('EAS')
      choice_refgeo           = C%choice_refgeo_PD_EAS
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_EAS
      timeframe_refgeo        = C%timeframe_refgeo_PD_EAS
    case ('GRL')
      choice_refgeo           = C%choice_refgeo_PD_GRL
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_GRL
      timeframe_refgeo        = C%timeframe_refgeo_PD_GRL
    case ('ANT')
      choice_refgeo           = C%choice_refgeo_PD_ANT
      choice_refgeo_idealised = C%choice_refgeo_PD_idealised
      dx_refgeo_idealised     = C%dx_refgeo_PD_idealised
      filename_refgeo         = C%filename_refgeo_PD_ANT
      timeframe_refgeo        = C%timeframe_refgeo_PD_ANT
    end select

    ! Initialise reference geometry
    call initialise_reference_geometry_raw( region_name, 'present-day', refgeo_PD, &
      choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)

    ! == GIA equilibrium geometry
    ! ===========================

    ! Get config choices for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name))
    case ('NAM')
      choice_refgeo           = C%choice_refgeo_GIAeq_NAM
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_NAM
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_NAM
    case ('EAS')
      choice_refgeo           = C%choice_refgeo_GIAeq_EAS
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_EAS
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_EAS
    case ('GRL')
      choice_refgeo           = C%choice_refgeo_GIAeq_GRL
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_GRL
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_GRL
    case ('ANT')
      choice_refgeo           = C%choice_refgeo_GIAeq_ANT
      choice_refgeo_idealised = C%choice_refgeo_GIAeq_idealised
      dx_refgeo_idealised     = C%dx_refgeo_GIAeq_idealised
      filename_refgeo         = C%filename_refgeo_GIAeq_ANT
      timeframe_refgeo        = C%timeframe_refgeo_GIAeq_ANT
    end select

    ! Initialise reference geometry
    call initialise_reference_geometry_raw( region_name, 'GIA equilibrium', refgeo_GIAeq, &
      choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometries_raw

  subroutine initialise_reference_geometry_raw( region_name, refgeo_name, refgeo, &
    choice_refgeo, choice_refgeo_idealised, dx_refgeo_idealised, filename_refgeo, timeframe_refgeo)
    !< Initialise a reference geometry on the raw grid/mesh

    ! In/output variables:
    character(len=3)             , intent(in   ) :: region_name
    character(len=*)             , intent(in   ) :: refgeo_name
    type(type_reference_geometry), intent(  out) :: refgeo
    character(len=1024)          , intent(in   ) :: choice_refgeo
    character(len=1024)          , intent(in   ) :: choice_refgeo_idealised
    real(dp)                     , intent(in   ) :: dx_refgeo_idealised
    character(len=1024)          , intent(in   ) :: filename_refgeo
    real(dp)                     , intent(in   ) :: timeframe_refgeo

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_reference_geometry_raw'

    ! Add routine to path
    call init_routine( routine_name)

    ! Clean up memory if necessary
    if (allocated( refgeo%Hi_grid_raw)) deallocate( refgeo%Hi_grid_raw)
    if (allocated( refgeo%Hb_grid_raw)) deallocate( refgeo%Hb_grid_raw)
    if (allocated( refgeo%Hs_grid_raw)) deallocate( refgeo%Hs_grid_raw)
    if (allocated( refgeo%SL_grid_raw)) deallocate( refgeo%SL_grid_raw)

    select case (choice_refgeo)
    case default
      call crash('unknown choice_refgeo "' // trim( choice_refgeo))
    case ('idealised')
      call initialise_reference_geometry_raw_idealised( region_name, refgeo_name, refgeo, &
      choice_refgeo_idealised, dx_refgeo_idealised)
    case ('read_from_file')
      call initialise_reference_geometry_raw_from_file( region_name, refgeo_name, refgeo, &
      filename_refgeo, timeframe_refgeo)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_raw

  subroutine initialise_reference_geometry_raw_idealised( region_name, refgeo_name, refgeo, &
    choice_refgeo_idealised, dx_refgeo_idealised)
    !< Initialise an idealised reference geometry on the raw grid

    ! In/output variables:
    character(len=3)             , intent(in   ) :: region_name
    character(len=*)             , intent(in   ) :: refgeo_name
    type(type_reference_geometry), intent(  out) :: refgeo
    character(len=1024)          , intent(in   ) :: choice_refgeo_idealised
    real(dp)                     , intent(in   ) :: dx_refgeo_idealised

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'initialise_reference_geometry_raw_idealised'
    real(dp)                              :: xmin, xmax, ymin, ymax
    character(len=1024)                   :: name
    real(dp), dimension(:,:), allocatable :: Hi, Hb, Hs, SL
    integer                               :: i,j

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to screen
    if (par%primary) write(0,'(A)') '  Initialising ' // trim( refgeo_name) // ' geometry for model region ' // &
      colour_string( region_name,'light blue') // ' from idealised case "' // colour_string( trim( choice_refgeo_idealised),'light blue') // '"...'

    ! Get domain size for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name))
    case ('NAM')
      xmin = C%xmin_NAM
      xmax = C%xmax_NAM
      ymin = C%ymin_NAM
      ymax = C%ymax_NAM
    case ('EAS')
      xmin = C%xmin_EAS
      xmax = C%xmax_EAS
      ymin = C%ymin_EAS
      ymax = C%ymax_EAS
    case ('GRL')
      xmin = C%xmin_GRL
      xmax = C%xmax_GRL
      ymin = C%ymin_GRL
      ymax = C%ymax_GRL
    case ('ANT')
      xmin = C%xmin_ANT
      xmax = C%xmax_ANT
      ymin = C%ymin_ANT
      ymax = C%ymax_ANT
    end select

    ! Set up a square grid to generate the idealised geometry on
    name = 'temp_grid_for_idealised_geometry_lines'
    call setup_square_grid( name, xmin, xmax, ymin, ymax, dx_refgeo_idealised, refgeo%grid_raw)

    ! allocate memory for partial grid data
    allocate( refgeo%Hi_grid_raw( refgeo%grid_raw%n1: refgeo%grid_raw%n2), source = 0._dp)
    allocate( refgeo%Hb_grid_raw( refgeo%grid_raw%n1: refgeo%grid_raw%n2), source = 0._dp)
    allocate( refgeo%Hs_grid_raw( refgeo%grid_raw%n1: refgeo%grid_raw%n2), source = 0._dp)
    allocate( refgeo%SL_grid_raw( refgeo%grid_raw%n1: refgeo%grid_raw%n2), source = 0._dp)

    ! allocate memory for the full grid data on the primary
    if (par%primary) then
      allocate( Hi( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      allocate( Hb( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      allocate( Hs( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      allocate( SL( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
    end if

    ! Calculate the idealised geometry on the grid (primary only)
    if (par%primary) then
      do i = 1, refgeo%grid_raw%nx
      do j = 1, refgeo%grid_raw%ny
        call calc_idealised_geometry( refgeo%grid_raw%x( i), refgeo%grid_raw%y( j), Hi( i,j), Hb( i,j), Hs( i,j), SL( i,j), choice_refgeo_idealised)
      end do
      end do
    end if

    ! Distribute the data over the processes in vector form
    call distribute_gridded_data_from_primary( refgeo%grid_raw, Hi, refgeo%Hi_grid_raw)
    call distribute_gridded_data_from_primary( refgeo%grid_raw, Hb, refgeo%Hb_grid_raw)
    call distribute_gridded_data_from_primary( refgeo%grid_raw, Hs, refgeo%Hs_grid_raw)
    call distribute_gridded_data_from_primary( refgeo%grid_raw, SL, refgeo%SL_grid_raw)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_raw_idealised

  subroutine initialise_reference_geometry_raw_from_file( region_name, refgeo_name, refgeo, filename_refgeo, timeframe_refgeo)
    !< Initialise a raw reference geometry from a file

    ! In/output variables:
    character(len=3)             , intent(in   ) :: region_name
    character(len=*)             , intent(in   ) :: refgeo_name
    type(type_reference_geometry), intent(  out) :: refgeo
    character(len=*)             , intent(in   ) :: filename_refgeo
    real(dp)                     , intent(in   ) :: timeframe_refgeo

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_reference_geometry_raw_from_file'
    character(len=1024)           :: filename_refgeo_applied
    real(dp)                      :: timeframe_refgeo_applied
    logical                       :: has_xy_grid, has_lonlat_grid, has_mesh

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    filename_refgeo_applied  = filename_refgeo
    timeframe_refgeo_applied =  timeframe_refgeo
    if (index( filename_refgeo,'_LAST.nc') > 1) then
      call find_last_output_file( filename_refgeo_applied)
      call find_last_timeframe(   filename_refgeo_applied, timeframe_refgeo_applied)
    end if

    ! Print to screen
    if (par%primary) write(0,'(A)') '   Initialising ' // trim( refgeo_name) // ' geometry for model region ' // &
      colour_string( region_name,'light blue') // ' from file "' // colour_string( trim( filename_refgeo_applied),'light blue') // '"...'

    ! Find out on what kind of grid the file is defined
    call inquire_xy_grid(     filename_refgeo_applied, has_xy_grid    )
    call inquire_lonlat_grid( filename_refgeo_applied, has_lonlat_grid)
    call inquire_mesh(        filename_refgeo_applied, has_mesh       )

    ! Files with more than one grid are not recognised
    if (has_xy_grid     .and. has_lonlat_grid) call crash('file "' // trim( filename_refgeo_applied) // '" contains both an x/y-grid and a lon/lat-grid!')
    if (has_xy_grid     .and. has_mesh       ) call crash('file "' // trim( filename_refgeo_applied) // '" contains both an x/y-grid and a mesh!')
    if (has_lonlat_grid .and. has_mesh       ) call crash('file "' // trim( filename_refgeo_applied) // '" contains both a lon/lat-grid and a mesh!')

    ! Read the grid/mesh and data from the file
    if (has_xy_grid) then
      call initialise_reference_geometry_raw_from_file_grid( region_name, refgeo, filename_refgeo_applied, timeframe_refgeo_applied)
    elseif (has_mesh) then
      call initialise_reference_geometry_raw_from_file_mesh(              refgeo, filename_refgeo_applied, timeframe_refgeo_applied)
    else
      call crash('can only read reference geometry from gridded or meshed data files!')
    endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_raw_from_file

  subroutine initialise_reference_geometry_raw_from_file_grid( region_name, refgeo, filename_refgeo, timeframe_refgeo)
    !< Initialise a reference geometry on the raw grid from a NetCDF file

    ! In/output variables:
    character(len=3)             , intent(in   ) :: region_name
    type(type_reference_geometry), intent(  out) :: refgeo
    character(len=1024)          , intent(in   ) :: filename_refgeo
    real(dp)                     , intent(in   ) :: timeframe_refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_reference_geometry_raw_from_file_grid'
    integer                        :: ncid, id_var
    logical                        :: has_SL
    integer                        :: n

    ! Add routine to path
    call init_routine( routine_name)

    ! Set up the grid from the file
    call open_existing_netcdf_file_for_reading( filename_refgeo, ncid)
    call setup_xy_grid_from_file( filename_refgeo, ncid, refgeo%grid_raw)
    call close_netcdf_file( ncid)

    ! allocate memory for the raw gridded data
    allocate( refgeo%Hi_grid_raw( refgeo%grid_raw%n1:refgeo%grid_raw%n2))
    allocate( refgeo%Hb_grid_raw( refgeo%grid_raw%n1:refgeo%grid_raw%n2))
    allocate( refgeo%Hs_grid_raw( refgeo%grid_raw%n1:refgeo%grid_raw%n2))
    allocate( refgeo%SL_grid_raw( refgeo%grid_raw%n1:refgeo%grid_raw%n2))

    ! Check if a sea level variable exists in the file
    call open_existing_netcdf_file_for_reading( filename_refgeo, ncid)
    call inquire_var_multopt( filename_refgeo, ncid, 'default_options_SL', id_var)
    has_SL = id_var /= -1
    call close_netcdf_file( ncid)

    if (timeframe_refgeo /= 1E9_dp) then
      ! We need to read a specific time frame

      call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_grid_raw, time_to_read = timeframe_refgeo)
      call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_grid_raw, time_to_read = timeframe_refgeo)
      call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_grid_raw, time_to_read = timeframe_refgeo)

      ! if the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      if (has_SL) then
        call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_SL', refgeo%SL_grid_raw, time_to_read = timeframe_refgeo)
      else
        refgeo%SL_grid_raw = 0._dp
      end if

    else !  if (timeframe_refgeo /= 1E9_dp) then
      ! We need to read data from a time-less NetCDF file

      call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_grid_raw)
      call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_grid_raw)
      call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_grid_raw)

      ! if the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      if (has_SL) then
        call read_field_from_xy_file_dp_2D( filename_refgeo, 'default_options_SL', refgeo%SL_grid_raw)
      else
        refgeo%SL_grid_raw = 0._dp
      end if

    end if !  if (timeframe_refgeo /= 1E9_dp) then

    ! == Input data clean-up
    ! ======================

    ! if so desired, smooth the bedrock geometry
    call smooth_model_geometry( refgeo)

    ! if so specified, remove Lake Vostok from the Antarctic geometry
    if (region_name == 'ANT' .AND. C%remove_Lake_Vostok) then
      call remove_Lake_Vostok( refgeo)
    end if

    ! if so specified, remove Ellesmere Island from the Greenland geometry
    if (region_name == 'GRL' .AND. C%choice_mask_noice == 'remove_Ellesmere') then
      call remove_Ellesmere( refgeo)
    end if

    ! Remove a few islands from the Antarctic geometry
    if (region_name == 'ANT' .AND. C%choice_refgeo_PD_ANT /= 'idealised') then
      call remove_tiny_islands( refgeo)
    end if

    ! Remove extremely thin ice (especially a problem in BedMachine Greenland)
    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      if (refgeo%Hi_grid_raw( n) < C%refgeo_Hi_min) then
        refgeo%Hi_grid_raw( n) = 0._dp
      end if
    end do

    ! Assume ice thickness is now correct everywhere; recalculate surface elevation from that
    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      refgeo%Hs_grid_raw( n) = ice_surface_elevation( refgeo%Hi_grid_raw( n), refgeo%Hb_grid_raw( n), refgeo%SL_grid_raw( n))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_raw_from_file_grid

  subroutine initialise_reference_geometry_raw_from_file_mesh( refgeo, filename_refgeo, timeframe_refgeo)
    !< Initialise a reference geometry on the raw mesh from a NetCDF file

    ! In/output variables:
    type(type_reference_geometry), intent(  out) :: refgeo
    character(len=*)             , intent(in   ) :: filename_refgeo
    real(dp)                     , intent(in   ) :: timeframe_refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_reference_geometry_raw_from_file_mesh'
    integer                        :: ncid, id_var, vi
    logical                        :: has_SL

    ! Add routine to path
    call init_routine( routine_name)

    ! Set up the mesh from the file
    call open_existing_netcdf_file_for_reading( filename_refgeo, ncid)
    call setup_mesh_from_file( filename_refgeo, ncid, refgeo%mesh_raw)
    call close_netcdf_file( ncid)

    ! allocate memory for the raw meshed data
    allocate( refgeo%Hi_mesh_raw( refgeo%mesh_raw%vi1:refgeo%mesh_raw%vi2))
    allocate( refgeo%Hb_mesh_raw( refgeo%mesh_raw%vi1:refgeo%mesh_raw%vi2))
    allocate( refgeo%Hs_mesh_raw( refgeo%mesh_raw%vi1:refgeo%mesh_raw%vi2))
    allocate( refgeo%SL_mesh_raw( refgeo%mesh_raw%vi1:refgeo%mesh_raw%vi2))

    ! Check if a sea level variable exists in the file
    call open_existing_netcdf_file_for_reading( filename_refgeo, ncid)
    call inquire_var_multopt( filename_refgeo, ncid, 'default_options_SL', id_var)
    has_SL = id_var /= -1
    call close_netcdf_file( ncid)

    if (timeframe_refgeo /= 1E9_dp) then
      ! We need to read a specific time frame

      call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_mesh_raw, time_to_read = timeframe_refgeo)
      call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_mesh_raw, time_to_read = timeframe_refgeo)
      call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_mesh_raw, time_to_read = timeframe_refgeo)

      ! if the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      if (has_SL) then
        call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_SL', refgeo%SL_mesh_raw, time_to_read = timeframe_refgeo)
      else
        refgeo%SL_mesh_raw = 0._dp
      end if

    else !  if (timeframe_refgeo /= 1E9_dp) then
      ! We need to read data from a time-less NetCDF file

      call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_Hi', refgeo%Hi_mesh_raw)
      call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_Hb', refgeo%Hb_mesh_raw)
      call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_Hs', refgeo%Hs_mesh_raw)

      ! if the file has a sea-level field, read that; if not, assume present-day (i.e. zero)
      if (has_SL) then
        call read_field_from_mesh_file_dp_2D( filename_refgeo, 'default_options_SL', refgeo%SL_mesh_raw)
      else
        refgeo%SL_mesh_raw = 0._dp
      end if

    end if !  if (timeframe_refgeo /= 1E9_dp) then

    ! Remove extremely thin ice (especially a problem in BedMachine Greenland)
    do vi = refgeo%mesh_raw%vi1, refgeo%mesh_raw%vi2
      if (refgeo%Hi_mesh_raw( vi) < C%refgeo_Hi_min) then
        refgeo%Hi_mesh_raw( vi) = 0._dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_raw_from_file_mesh

  ! Initialise idealised geometries on the model mesh
  ! =================================================

  subroutine initialise_reference_geometry_idealised( mesh, choice_refgeo_idealised, refgeo)
    !< Initialise an idealised reference geometry on the model mesh

    ! In/output variables
    type(type_mesh)              , intent(in   ) :: mesh
    character(len=1024)          , intent(in   ) :: choice_refgeo_idealised
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_reference_geometry_idealised'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      call calc_idealised_geometry( mesh%V( vi,1), mesh%V( vi,2), &
        refgeo%Hi( vi), refgeo%Hb( vi), refgeo%Hs( vi), refgeo%SL( vi), &
        choice_refgeo_idealised)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_idealised

  ! Utilities
  ! =========

  subroutine reallocate_reference_geometry_on_mesh( mesh, refgeo)

    ! In/output variables
    type(type_mesh)              , intent(in   ) :: mesh
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variabler:
    character(len=1024), parameter :: routine_name = 'deallocate_reference_geometry_on_mesh'

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( refgeo%Hi)) call deallocate_reference_geometry_on_mesh( refgeo)

    allocate( refgeo%Hi( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( refgeo%Hb( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( refgeo%Hs( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( refgeo%SL( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reallocate_reference_geometry_on_mesh

  subroutine deallocate_reference_geometry_on_mesh( refgeo)

    ! In/output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_reference_geometry_on_mesh'

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( refgeo%Hi)) deallocate( refgeo%Hi)
    if (allocated( refgeo%Hb)) deallocate( refgeo%Hb)
    if (allocated( refgeo%Hs)) deallocate( refgeo%Hs)
    if (allocated( refgeo%SL)) deallocate( refgeo%SL)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_reference_geometry_on_mesh

end module reference_geometries_main
