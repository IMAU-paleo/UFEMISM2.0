module reference_geometries_main

  ! Contains the routines for setting up the three "reference geometries":
  ! - refgeo_init:   initial, used to initialise the simulation
  ! - refgeo_PD:     present-day, used to calculate sea-level contribution, isotope change, and more
  ! - refgeo_GIA_eq: GIA equilibrium, used for the GIA model

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use reference_geometry_types, only: type_reference_geometry
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid, setup_square_grid
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_master
  use smooth_gridded_data, only: smooth_Gaussian_grid
  use ice_geometry_basics, only: ice_surface_elevation, is_floating
  use projections, only: oblique_sg_projection
  use analytical_solutions, only: Halfar_dome, Bueler_dome
  use netcdf_io_main
  use remapping_main, only: map_from_xy_grid_to_mesh_2D, map_from_mesh_to_mesh_2D

  implicit none

  private

  public :: initialise_reference_geometries_raw, initialise_reference_geometries_on_model_mesh

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
    if (par%master) write(0,'(A)') '     Mapping reference geometries to model mesh...'

    ! == Initial geometry
    ! ===================

    ! deallocate existing memory if needed
    if (allocated( refgeo_init%Hi)) deallocate( refgeo_init%Hi)
    if (allocated( refgeo_init%Hb)) deallocate( refgeo_init%Hb)
    if (allocated( refgeo_init%Hs)) deallocate( refgeo_init%Hs)
    if (allocated( refgeo_init%SL)) deallocate( refgeo_init%SL)

    ! allocate memory for reference ice geometry on the model mesh
    allocate( refgeo_init%Hi( mesh%vi1:mesh%vi2))
    allocate( refgeo_init%Hb( mesh%vi1:mesh%vi2))
    allocate( refgeo_init%Hs( mesh%vi1:mesh%vi2))
    allocate( refgeo_init%SL( mesh%vi1:mesh%vi2))

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

    ! deallocate existing memory if needed
    if (allocated( refgeo_PD%Hi)) deallocate( refgeo_PD%Hi)
    if (allocated( refgeo_PD%Hb)) deallocate( refgeo_PD%Hb)
    if (allocated( refgeo_PD%Hs)) deallocate( refgeo_PD%Hs)
    if (allocated( refgeo_PD%SL)) deallocate( refgeo_PD%SL)

    ! allocate memory for reference ice geometry on the model mesh
    allocate( refgeo_PD%Hi( mesh%vi1:mesh%vi2))
    allocate( refgeo_PD%Hb( mesh%vi1:mesh%vi2))
    allocate( refgeo_PD%Hs( mesh%vi1:mesh%vi2))
    allocate( refgeo_PD%SL( mesh%vi1:mesh%vi2))

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

    ! deallocate existing memory if needed
    if (allocated( refgeo_GIAeq%Hi)) deallocate( refgeo_GIAeq%Hi)
    if (allocated( refgeo_GIAeq%Hb)) deallocate( refgeo_GIAeq%Hb)
    if (allocated( refgeo_GIAeq%Hs)) deallocate( refgeo_GIAeq%Hs)
    if (allocated( refgeo_GIAeq%SL)) deallocate( refgeo_GIAeq%SL)

    ! allocate memory for reference ice geometry on the model mesh
    allocate( refgeo_GIAeq%Hi( mesh%vi1:mesh%vi2))
    allocate( refgeo_GIAeq%Hb( mesh%vi1:mesh%vi2))
    allocate( refgeo_GIAeq%Hs( mesh%vi1:mesh%vi2))
    allocate( refgeo_GIAeq%SL( mesh%vi1:mesh%vi2))

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

    ! deallocate existing memory if needed
    if (allocated( refgeo%Hi)) deallocate( refgeo%Hi)
    if (allocated( refgeo%Hb)) deallocate( refgeo%Hb)
    if (allocated( refgeo%Hs)) deallocate( refgeo%Hs)
    if (allocated( refgeo%SL)) deallocate( refgeo%SL)

    ! allocate memory for reference ice geometry on the model mesh
    allocate( refgeo%Hi( mesh%vi1:mesh%vi2))
    allocate( refgeo%Hb( mesh%vi1:mesh%vi2))
    allocate( refgeo%Hs( mesh%vi1:mesh%vi2))
    allocate( refgeo%SL( mesh%vi1:mesh%vi2))

    ! Determine if the initial geometry is provided gridded or meshed
    if (allocated( refgeo%grid_raw%x)) then
      ! Gridded

      ! Safety
      if (allocated( refgeo%mesh_raw%V)) call crash('found both grid and mesh in refgeo!')

      ! Remap data to the model mesh
      call map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, refgeo%Hi_grid_raw, refgeo%Hi)
      call map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, refgeo%Hb_grid_raw, refgeo%Hb)
      call map_from_xy_grid_to_mesh_2D( refgeo%grid_raw, mesh, refgeo%SL_grid_raw, refgeo%SL)

    elseif (allocated( refgeo%mesh_raw%V)) then
      ! Meshed

      ! Safety
      if (allocated( refgeo%grid_raw%x)) call crash('found both grid and mesh in refgeo!')

      ! Remap data to the model mesh
      method = '2nd_order_conservative'
      call map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, refgeo%Hi_mesh_raw, refgeo%Hi, method)
      call map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, refgeo%Hb_mesh_raw, refgeo%Hb, method)
      call map_from_mesh_to_mesh_2D( refgeo%mesh_raw, mesh, refgeo%SL_mesh_raw, refgeo%SL, method)

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
    if (par%master) write(0,'(A)') '  Initialising ' // trim( refgeo_name) // ' geometry for model region ' // &
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

    ! allocate memory for the full grid data on the master
    if (par%master) then
      allocate( Hi( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      allocate( Hb( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      allocate( Hs( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
      allocate( SL( refgeo%grid_raw%nx, refgeo%grid_raw%ny), source = 0._dp)
    end if

    ! Calculate the idealised geometry on the grid (master only)
    if (par%master) then
      do i = 1, refgeo%grid_raw%nx
      do j = 1, refgeo%grid_raw%ny
        call calc_idealised_geometry( refgeo%grid_raw%x( i), refgeo%grid_raw%y( j), Hi( i,j), Hb( i,j), Hs( i,j), SL( i,j), choice_refgeo_idealised)
      end do
      end do
    end if

    ! Distribute the data over the processes in vector form
    call distribute_gridded_data_from_master( refgeo%grid_raw, Hi, refgeo%Hi_grid_raw)
    call distribute_gridded_data_from_master( refgeo%grid_raw, Hb, refgeo%Hb_grid_raw)
    call distribute_gridded_data_from_master( refgeo%grid_raw, Hs, refgeo%Hs_grid_raw)
    call distribute_gridded_data_from_master( refgeo%grid_raw, SL, refgeo%SL_grid_raw)

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
    if (par%master) write(0,'(A)') '   Initialising ' // trim( refgeo_name) // ' geometry for model region ' // &
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
    integer                        :: ncid, id_var
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
      call calc_idealised_geometry( mesh%V( vi,1), mesh%V( vi,2), refgeo%Hi( vi), refgeo%Hb( vi), refgeo%Hs( vi), refgeo%SL( vi), choice_refgeo_idealised)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_idealised

  ! Calculate various idealsed geometries
  ! =====================================

  subroutine calc_idealised_geometry( x, y, Hi, Hb, Hs, SL, choice_refgeo_idealised)
    !< Calculate an idealised geometry

    ! In/output variables:
    real(dp),            intent(in   ) :: x,y             ! [m] Coordinates
    real(dp),            intent(  out) :: Hi              ! [m] Ice thickness
    real(dp),            intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp),            intent(  out) :: Hs              ! [m] Surface elevation
    real(dp),            intent(  out) :: SL              ! [m] Sea level
    character(len=1024), intent(in   ) :: choice_refgeo_idealised

    ! Calculated the specified idealised geometry
    select case (choice_refgeo_idealised)
    case default
      call crash('unknown choice_refgeo_idealised "' // trim( choice_refgeo_idealised))
    case ('flatearth')
      call calc_idealised_geometry_flatearth( Hi, Hb, Hs, SL)
    case ('slabonaslope')
      call calc_idealised_geometry_slabonaslope( x, Hi, Hb, Hs, SL)
    case ('Halfar')
      call calc_idealised_geometry_Halfar( x, y, Hi, Hb, Hs, SL)
    case ('Bueler')
      call calc_idealised_geometry_Bueler( x, y, Hi, Hb, Hs, SL)
    case ('SSA_icestream')
      call calc_idealised_geometry_SSA_icestream( x, Hi, Hb, Hs, SL)
    case ('MISMIP_mod')
      call calc_idealised_geometry_MISMIP_mod( x, y, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_A')
      call calc_idealised_geometry_ISMIP_HOM_A( x, y, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_B')
      call calc_idealised_geometry_ISMIP_HOM_B( x, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_C', 'ISMIP-HOM_D')
      call calc_idealised_geometry_ISMIP_HOM_CD( x, Hi, Hb, Hs, SL)
    case ('ISMIP-HOM_E')
      call crash('ISMIP-HOM E is not implemented in UFEMISM!')
    case ('ISMIP-HOM_F')
      call calc_idealised_geometry_ISMIP_HOM_F( x, y, Hi, Hb, Hs, SL)
    case ('MISMIP+', 'MISMIPplus')
      call calc_idealised_geometry_MISMIPplus( x, y, Hi, Hb, Hs, SL)
    case ('calvmip_circular')
      call calc_idealised_geometry_CalvMIP_circular( x, y, Hi, Hb, Hs, SL)
    case ('calvmip_Thule')
      call calc_idealised_geometry_CalvMIP_Thule( x, y, Hi, Hb, Hs, SL)
    end select

  end subroutine calc_idealised_geometry

  subroutine calc_idealised_geometry_flatearth( Hi, Hb, Hs, SL)
    !< Calculate a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments

    ! In/output variables:
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    Hi = 0._dp
    Hb = 0._dp
    Hs = 0._dp
    SL = -10000._dp

  end subroutine calc_idealised_geometry_flatearth

  subroutine calc_idealised_geometry_slabonaslope( x, Hi, Hb, Hs, SL)
    !< Calculate a 2,000 m thick slab of ice on a flat, inclined plane

    ! In/output variables:
    real(dp), intent(in   ) :: x               ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_slabonaslope_Hi < 100._dp .or. &
         C%refgeo_idealised_slabonaslope_Hi > 10000._dp) then
      call crash('refgeo_idealised_slabonaslope_Hi has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_slabonaslope_Hi)
    end if
    if ( ABS( C%refgeo_idealised_slabonaslope_dhdx) > 0.1_dp) then
      call crash('refgeo_idealised_slabonaslope_dhdx has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_slabonaslope_dhdx)
    end if
#endif

    Hi = C%refgeo_idealised_slabonaslope_Hi
    Hb = C%refgeo_idealised_slabonaslope_dhdx * x
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_slabonaslope

  subroutine calc_idealised_geometry_Halfar( x, y, Hi, Hb, Hs, SL)
    !< Calculate the Halfar dome solution at t = 0

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%uniform_Glens_flow_factor < 1E-18_dp .OR. &
         C%uniform_Glens_flow_factor > 1E-15_dp) then
      call crash('uniform_flow_factor has unrealistic value of {dp_01}!', dp_01 = C%uniform_Glens_flow_factor)
    end if
    if ( C%Glens_flow_law_exponent < 1._dp .OR. C%Glens_flow_law_exponent > 5._dp) then
      call crash('Glens_flow_law_exponent has unrealistic value of {dp_01}!', dp_01 = C%Glens_flow_law_exponent)
    end if
    if ( C%refgeo_idealised_Halfar_H0 < 100._dp .OR. C%refgeo_idealised_Halfar_H0 > 10000._dp) then
      call crash('refgeo_idealised_Halfar_H0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Halfar_H0)
    end if
    if ( C%refgeo_idealised_Halfar_R0 < 100E3_dp .OR. C%refgeo_idealised_Halfar_R0 > 5000E3_dp) then
      call crash('refgeo_idealised_Halfar_R0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Halfar_R0)
    end if
#endif

    call Halfar_dome( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_Halfar_H0, C%refgeo_idealised_Halfar_R0, &
      x, y, 0._dp, Hi)
    Hb = 0._dp
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_Halfar

  subroutine calc_idealised_geometry_Bueler( x, y, Hi, Hb, Hs, SL)
    !< Calculate the Bueler dome solution at t = 0

    ! In/output variables:
    real(dp),                       intent(in   ) :: x,y             ! [m] Coordinates
    real(dp),                       intent(  out) :: Hi              ! [m] Ice thickness
    real(dp),                       intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp),                       intent(  out) :: Hs              ! [m] Surface elevation
    real(dp),                       intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp) :: M

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%uniform_Glens_flow_factor < 1E-18_dp .OR. &
         C%uniform_Glens_flow_factor > 1E-15_dp) then
      call crash('uniform_flow_factor has unrealistic value of {dp_01}!', dp_01 = C%uniform_Glens_flow_factor)
    end if
    if ( C%Glens_flow_law_exponent < 1._dp .OR. C%Glens_flow_law_exponent > 5._dp) then
      call crash('Glens_flow_law_exponent has unrealistic value of {dp_01}!', dp_01 = C%Glens_flow_law_exponent)
    end if
    if ( C%refgeo_idealised_Bueler_H0 < 100._dp .OR. C%refgeo_idealised_Bueler_H0 > 10000._dp) then
      call crash('refgeo_idealised_Bueler_H0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_H0)
    end if
    if ( C%refgeo_idealised_Bueler_R0 < 100E3_dp .OR. C%refgeo_idealised_Bueler_R0 > 5000E3_dp) then
      call crash('refgeo_idealised_Bueler_R0 has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_R0)
    end if
    if ( C%refgeo_idealised_Bueler_lambda < -10._dp .OR. C%refgeo_idealised_Bueler_lambda > 10._dp) then
      call crash('refgeo_idealised_Bueler_lambda has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_Bueler_lambda)
    end if
#endif

    call Bueler_dome( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_Bueler_H0, C%refgeo_idealised_Bueler_R0, &
      C%refgeo_idealised_Bueler_lambda, x, y, 0._dp, Hi, M)
    Hb = 0._dp
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_Bueler

  subroutine calc_idealised_geometry_SSA_icestream( x, Hi, Hb, Hs, SL)
    !< Calculate a thick slab of ice on a flat, inclined plane

    ! In/output variables:
    real(dp), intent(in   ) :: x               ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_SSA_icestream_Hi < 100._dp .OR. &
         C%refgeo_idealised_SSA_icestream_Hi > 10000._dp) then
      call crash('refgeo_idealised_SSA_icestream_Hi has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_SSA_icestream_Hi)
    end if
    if ( abs( C%refgeo_idealised_SSA_icestream_dhdx) > 0.1_dp) then
      call crash('refgeo_idealised_SSA_icestream_dhdx has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_SSA_icestream_dhdx)
    end if
#endif

    Hi = C%refgeo_idealised_SSA_icestream_Hi
    Hb = C%refgeo_idealised_SSA_icestream_dhdx * x
    SL = -10000._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_SSA_icestream

  subroutine calc_idealised_geometry_MISMIP_mod( x, y, Hi, Hb, Hs, SL)
    !< Calculate the MISMIP_mod cone-shaped island

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_MISMIP_mod_Hi_init < 0._dp .OR. &
         C%refgeo_idealised_MISMIP_mod_Hi_init > 10000._dp) then
      call crash('refgeo_idealised_MISMIP_mod_Hi_init has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_MISMIP_mod_Hi_init)
    end if
#endif

    if (sqrt( x**2 + y**2) > 900E3_dp) then
      Hi = 0._dp
    else
      Hi = C%refgeo_idealised_MISMIP_mod_Hi_init
    end if
    Hb = 150._dp - 400._dp * sqrt( x**2 + y**2)/ 750000._dp
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_MISMIP_mod

  subroutine calc_idealised_geometry_ISMIP_HOM_A( x, y, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment A (slab on a bumpy slope in both directions)

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 2000._dp - x * tan( 0.5_dp * pi / 180._dp)
    Hb = Hs - 1000._dp + 500._dp * sin( x * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L) &
                                 * sin( y * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L)
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_A

  subroutine calc_idealised_geometry_ISMIP_HOM_B( x, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment B (slab on a bumpy slope in only the x-directions)

    ! In/output variables:
    real(dp), intent(in   ) :: x               ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 2000._dp - x * TAN( 0.5_dp * pi / 180._dp)
    Hb = Hs - 1000._dp + 500._dp * SIN( x * 2._dp * pi / C%refgeo_idealised_ISMIP_HOM_L)
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_B

  subroutine calc_idealised_geometry_ISMIP_HOM_CD( x, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment C/D (slab on a slope)

    ! In/output variables:
    real(dp),                       intent(in   ) :: x               ! [m] Coordinates
    real(dp),                       intent(  out) :: Hi              ! [m] Ice thickness
    real(dp),                       intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp),                       intent(  out) :: Hs              ! [m] Surface elevation
    real(dp),                       intent(  out) :: SL              ! [m] Sea level

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 2000._dp - x * TAN( 0.1_dp * pi / 180._dp)
    Hb = Hs - 1000._dp
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_CD

  subroutine calc_idealised_geometry_ISMIP_HOM_F( x, y, Hi, Hb, Hs, SL)
    !< Calculate the geometry for ISMIP-HOM Experiment F (slab on another bumpy slope)

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp), parameter :: H0    = 1000._dp
    real(dp), parameter :: a0    = 100._dp
    real(dp), parameter :: sigma = 10000._dp

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_ISMIP_HOM_L < 2500._dp .OR. &
         C%refgeo_idealised_ISMIP_HOM_L > 320000._dp) then
      call crash('refgeo_idealised_ISMIP_HOM_L has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_ISMIP_HOM_L)
    end if
#endif

    Hs = 5000._dp - x * tan( 3._dp * pi / 180._dp)
    Hb = Hs - H0 + a0 * exp( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y - 0._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2) &
                 + a0 * exp( -((x + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2 + (y + 1._dp * C%refgeo_idealised_ISMIP_HOM_L)**2) / sigma**2)
    Hi = Hs - Hb
    SL = -10000._dp

  end subroutine calc_idealised_geometry_ISMIP_HOM_F

  subroutine calc_idealised_geometry_MISMIPplus( x, y, Hi, Hb, Hs, SL)
    !< Calculate the MISMIP+ geometry

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp)            :: xtilde,Bx,By
    real(dp), parameter :: B0     = -150._dp
    real(dp), parameter :: B2     = -728.8_dp
    real(dp), parameter :: B4     = 343.91_dp
    real(dp), parameter :: B6     = -50.57_dp
    real(dp), parameter :: xbar   = 300000._dp
    real(dp), parameter :: fc     = 4000._dp
    real(dp), parameter :: dc     = 500._dp
    real(dp), parameter :: wc     = 24000._dp
    real(dp), parameter :: zbdeep = -720._dp

#if (DO_ASSERTIONS)
    ! Safety
    if ( C%refgeo_idealised_MISMIPplus_Hi_init < 0._dp .OR. &
         C%refgeo_idealised_MISMIPplus_Hi_init > 10000._dp) then
      call crash('refgeo_idealised_MISMIPplus_Hi_init has unrealistic value of {dp_01}!', dp_01 = C%refgeo_idealised_MISMIPplus_Hi_init)
    end if
#endif

    xtilde = x / xbar
    Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
    By = (dc / (1 + exp(-2._dp*(y - wc)/fc))) + &
         (dc / (1 + exp( 2._dp*(y + wc)/fc)))

    if (x > 640E3_dp) then
      Hi = 0._dp
    else
      Hi = C%refgeo_idealised_MISMIPplus_Hi_init
    end if
    Hb = max( Bx + By, zbdeep)
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_MISMIPplus

  subroutine calc_idealised_geometry_CalvMIP_circular( x, y, Hi, Hb, Hs, SL)
    ! Calculate the geometry for the CalvingMIP circular domain

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp), parameter :: R  = 800E3_dp
    real(dp), parameter :: Bc = 900._dp
    real(dp), parameter :: Bl = -2000._dp
    real(dp), parameter :: Ba = 1100._dp
    real(dp), parameter :: rc = 0._dp
    real(dp)            :: radius, theta

    radius = sqrt(x**2 + y**2)
    theta  = atan2(y,x)

    Hi = 0._dp
    Hb = Bc - (Bc-Bl)*(radius-rc)**2 / (R-rc)**2
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_CalvMIP_circular

  subroutine calc_idealised_geometry_CalvMIP_Thule( x, y, Hi, Hb, Hs, SL)
    ! Calculate the geometry for the CalvingMIP Thule domain

    ! In/output variables:
    real(dp), intent(in   ) :: x,y             ! [m] Coordinates
    real(dp), intent(  out) :: Hi              ! [m] Ice thickness
    real(dp), intent(  out) :: Hb              ! [m] Bedrock elevation
    real(dp), intent(  out) :: Hs              ! [m] Surface elevation
    real(dp), intent(  out) :: SL              ! [m] Sea level

    ! Local variables:
    real(dp), parameter :: R  = 800E3_dp
    real(dp), parameter :: Bc = 900._dp
    real(dp), parameter :: Bl = -2000._dp
    real(dp), parameter :: Ba = 1100._dp
    real(dp), parameter :: rc = 0._dp
    real(dp)            :: radius, theta, l, a

    radius = sqrt(x**2 + y**2)
    theta  = atan2(y,x)
    l = R - cos( 2._dp * theta) * R / 2._dp
    a = Bc - (Bc-Bl)*(radius-rc)**2 / (R-rc)**2

    Hi = 0._dp
    Hb = Ba * cos( 3._dp * pi * radius / l) + a
    SL = 0._dp
    Hs = ice_surface_elevation( Hi, Hb, SL)

  end subroutine calc_idealised_geometry_CalvMIP_Thule

  ! Utilities
  ! =========

  subroutine smooth_model_geometry( refgeo)
    !< Apply some light smoothing to the initial geometry to improve numerical stability

    ! Input variables:
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter                             :: routine_name = 'smooth_model_geometry'
    integer                                                    :: n
    real(dp), dimension(refgeo%grid_raw%n1:refgeo%grid_raw%n2) :: Hb_old
    real(dp)                                                   :: r_smooth, dHb
    logical                                                    :: is_grounded

    ! Add routine to path
    call init_routine( routine_name)

    if ((.not. C%do_smooth_geometry) .or. C%r_smooth_geometry == 0._dp) then
      call finalise_routine( routine_name)
      return
    end if

    ! Store the unsmoothed bed topography so we can determine the smoothing anomaly later
    Hb_old = refgeo%Hb_grid_raw

    ! Geometry smoothing radius (in m)
    r_smooth = refgeo%grid_raw%dx * C%r_smooth_geometry

    ! Apply smoothing to the bed topography
    call smooth_Gaussian_grid( refgeo%grid_raw, refgeo%Hb_grid_raw, r_smooth)

    ! Correct smoothed geometry if necessary
    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2

      ! Calculate the smoothing anomaly
      dHb = refgeo%Hb_grid_raw( n) - Hb_old( n)

      is_grounded = .not. is_floating( refgeo%Hi_grid_raw( n), refgeo%Hb_grid_raw( n), refgeo%SL_grid_raw( n))

      if (is_grounded .and. refgeo%Hi_grid_raw( n) > 0._dp) then

        ! Correct the ice thickness so the ice surface remains unchanged (only relevant for grounded ice)
        refgeo%Hi_grid_raw( n) = refgeo%Hi_grid_raw( n) - dHb

        ! don't allow negative ice thickness
        refgeo%Hi_grid_raw( n) = MAX( 0._dp, refgeo%Hi_grid_raw( n))

      end if

      ! Correct the surface elevation if necessary
      refgeo%Hs_grid_raw( n) = refgeo%Hi_grid_raw( n) + MAX( refgeo%SL_grid_raw( n) - ice_density / seawater_density * refgeo%Hi_grid_raw( n), refgeo%Hb_grid_raw( n))

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_model_geometry

  subroutine remove_Lake_Vostok( refgeo)
    !< Remove Lake Vostok from Antarctic input geometry data
    !< by manually increasing ice thickness so that Hi = Hs - Hb

    ! NOTE: since UFEMISM doesn't consider subglacial lakes, Vostok simply shows
    !       up as a "dip" in the initial geometry. The model will run fine, the dip
    !       fills up in a few centuries, but it slows down the model for a while and
    !       it looks ugly, so we just remove it right away.

    ! In/output variables:
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_Lake_Vostok'
    real(dp), parameter            :: lake_Vostok_xmin = 1164250.0
    real(dp), parameter            :: lake_Vostok_xmax = 1514250.0
    real(dp), parameter            :: lake_Vostok_ymin = -470750.0
    real(dp), parameter            :: lake_Vostok_ymax = -220750.0
    integer                        :: il,iu,jl,ju
    integer                        :: n,i,j

    ! Add routine to path
    call init_routine( routine_name)

    il = 1
    do while (refgeo%grid_raw%x( il) < lake_Vostok_xmin)
      il = il+1
    end do
    iu = refgeo%grid_raw%nx
    do while (refgeo%grid_raw%x( iu) > lake_Vostok_xmax)
      iu = iu-1
    end do
    jl = 1
    do while (refgeo%grid_raw%y( jl) < lake_Vostok_ymin)
      jl = jl+1
    end do
    ju = refgeo%grid_raw%ny
    do while (refgeo%grid_raw%y( ju) > lake_Vostok_ymax)
      ju = ju-1
    end do

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)
      if (i >= il .and. i <= iu .and. j >= jl .and. j <= ju) then
        ! if we assume there's no subglacial water, then the entire column
        ! between bed and surface should be ice
        refgeo%Hi_grid_raw( n) = refgeo%Hs_grid_raw( n) - refgeo%Hb_grid_raw( n)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_Lake_Vostok

  subroutine remove_Ellesmere( refgeo)
    !< Remove ice from Ellesmere Island, which shows up in the Greenland domain

    ! NOTE: this routine assumes that Greenland input data use the ISMIP-standard projection.

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_Ellesmere'
    integer                        :: i,j,n
    real(dp), dimension(2)         :: pa_latlon, pb_latlon
    real(dp)                       :: xa,ya,xb,yb
    real(dp), dimension(2)         :: pa, pb
    real(dp)                       :: yl_ab

    ! Add routine to path
    call init_routine( routine_name)

    ! The two endpoints in lat,lon
    pa_latlon = [76.74_dp, -74.79_dp]
    pb_latlon = [82.19_dp, -60.00_dp]

    ! The two endpoints in x,y
    call oblique_sg_projection( pa_latlon(2), pa_latlon(1), C%lambda_M_GRL, C%phi_M_GRL, C%beta_stereo_GRL, xa, ya)
    call oblique_sg_projection( pb_latlon(2), pb_latlon(1), C%lambda_M_GRL, C%phi_M_GRL, C%beta_stereo_GRL, xb, yb)

    pa = [xa,ya]
    pb = [xb,yb]

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)

      ! Draw a line that separates Ellesmere from Greenland
      yl_ab = pa(2) + (refgeo%grid_raw%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))

      ! if grid cell is above the line, remove ice from it and
      ! actually sink the damn thing so it does not show up in
      ! the mesh when using a high-res based on coastlines.
      if (refgeo%grid_raw%y(j) > pa(2) .and. refgeo%grid_raw%y(j) > yl_ab .and. refgeo%grid_raw%x(i) < pb(1)) then
        refgeo%Hi_grid_raw( n) = 0._dp
        refgeo%Hb_grid_raw( n) = min( refgeo%Hb_grid_raw( n), -.1_dp)
        refgeo%Hs_grid_raw( n) = 0._dp
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_Ellesmere

  subroutine remove_tiny_islands( refgeo)
    !< Remove tiny islands from the Antarctic domain, so they do not
    !< cause unnecessary vertices there during mesh creation.

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_tiny_islands'
    integer                        :: i,j,n
    real(dp)                       :: x1, x2, y1, y2

    ! Add routine to path
    call init_routine( routine_name)

    ! Near tip of the peninsula
    ! =========================

    x1 = -2.5819e6_dp
    x2 = -2.0156e6_dp
    y1 = +2.1377e6_dp
    y2 = +2.5708e6_dp

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)

      ! if grid cell is within the lines, remove ice from it and
      ! actually sink the damn thing so it does not show up in
      ! the mesh when using high-res.
      if (refgeo%grid_raw%x(i) >=  min( x1,x2) .and. refgeo%grid_raw%x(i) <=  max( x1,x2) .and. &
          refgeo%grid_raw%y(j) >=  min( y1,y2) .and. refgeo%grid_raw%y(j) <=  max( y1,y2)) then
        refgeo%Hi_grid_raw( n) = 0._dp
        refgeo%Hb_grid_raw( n) = min( refgeo%Hb_grid_raw( n), -.1_dp)
        refgeo%Hs_grid_raw( n) = 0._dp
      end if

    end do

    ! The other ones
    ! ==============

    x1 = +0.4942e6_dp
    x2 = +0.9384e6_dp
    y1 = -2.6485e6_dp
    y2 = -2.2932e6_dp

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)

      ! if grid cell is within the lines, remove ice from it and
      ! actually sink the damn thing so it does not show up in
      ! the mesh when using high-res.
      if (refgeo%grid_raw%x(i) >=  min( x1,x2) .and. refgeo%grid_raw%x(i) <=  max( x1,x2) .and. &
          refgeo%grid_raw%y(j) >=  min( y1,y2) .and. refgeo%grid_raw%y(j) <=  max( y1,y2)) then
        refgeo%Hi_grid_raw( n) = 0._dp
        refgeo%Hb_grid_raw( n) = min( refgeo%Hb_grid_raw( n), -.1_dp)
        refgeo%Hs_grid_raw( n) = 0._dp
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_tiny_islands

end module reference_geometries_main
