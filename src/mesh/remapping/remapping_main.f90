module remapping_main

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use remapping_types, only: type_map
  use grid_types, only: type_grid, type_grid_lonlat
  use mesh_types, only: type_mesh
  use math_utilities, only: remap_cons_2nd_order_1D
  use mesh_utilities, only: extrapolate_Gaussian, check_if_meshes_are_identical
  use create_maps_grid_mesh
  use create_maps_gridlonlat_mesh
  use create_maps_mesh_mesh
  use apply_maps

  implicit none

  private

  public :: Atlas
  public :: map_from_xy_grid_to_mesh_2D, map_from_xy_grid_to_mesh_3D
  public :: map_from_xy_grid_to_mesh_triangles_2D, map_from_xy_grid_to_mesh_triangles_3D
  public :: map_from_lonlat_grid_to_mesh_2D, map_from_lonlat_grid_to_mesh_3D
  public :: map_from_mesh_to_xy_grid_2D, map_from_mesh_to_xy_grid_3D, map_from_mesh_to_xy_grid_2D_minval
  public :: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  public :: map_from_mesh_to_mesh_2D, map_from_mesh_to_mesh_3D
  public :: map_from_vertical_to_vertical_2D_ocean

contains

  ! From an x/y-grid to a mesh
  subroutine map_from_xy_grid_to_mesh_2D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(:    ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:    ),          intent(out)   :: d_mesh_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_xy_grid_to_mesh_2D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == grid%name .and. Atlas( mi)%name_dst == mesh%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_xy_grid_to_mesh_2D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_xy_grid_to_mesh_2D

  subroutine map_from_xy_grid_to_mesh_3D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 3-D data field from an x/y-grid to a mesh.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(:,:  ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_mesh_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_xy_grid_to_mesh_3D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == grid%name .and. Atlas( mi)%name_dst == mesh%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_xy_grid_to_mesh_3D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_xy_grid_to_mesh_3D

  ! From an x/y-grid to a mesh triangles
  subroutine map_from_xy_grid_to_mesh_triangles_2D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh triangles.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(:    ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:    ),          intent(out)   :: d_mesh_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_xy_grid_to_mesh_triangles_2D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == grid%name .and. Atlas( mi)%name_dst == (trim( mesh%name) // '_triangeles')) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh_triangles( grid, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_xy_grid_to_mesh_triangles_2D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_xy_grid_to_mesh_triangles_2D

  subroutine map_from_xy_grid_to_mesh_triangles_3D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 3-D data field from an x/y-grid to a mesh triangles.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(:,:  ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_mesh_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_xy_grid_to_mesh_triangles_3D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == grid%name .and. Atlas( mi)%name_dst == (trim( mesh%name) // '_triangles')) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh_triangles( grid, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_xy_grid_to_mesh_triangles_3D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_xy_grid_to_mesh_triangles_3D

  ! From a lon/lat-grid to a mesh
  subroutine map_from_lonlat_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 2-D data field from a lon/lat-grid to a mesh.

    ! In/output variables
    type(type_grid_lonlat),              intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(:    ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:    ),          intent(out)   :: d_mesh_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_lonlat_grid_to_mesh_2D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == grid%name .and. Atlas( mi)%name_dst == mesh%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_lonlat_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_lonlat_grid_to_mesh_2D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_lonlat_grid_to_mesh_2D

  subroutine map_from_lonlat_grid_to_mesh_3D( grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 3-D data field from a lon/lat-grid to a mesh.

    ! In/output variables
    type(type_grid_lonlat),              intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(:,:  ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_mesh_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_lonlat_grid_to_mesh_3D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == grid%name .and. Atlas( mi)%name_dst == mesh%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_lonlat_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_lonlat_grid_to_mesh_3D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_lonlat_grid_to_mesh_3D

  ! From a mesh to an x/y-grid
  subroutine map_from_mesh_to_xy_grid_2D( mesh, grid, d_mesh_partial, d_grid_vec_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(:    ),          intent(in)    :: d_mesh_partial
    real(dp), dimension(:    ),          intent(out)   :: d_grid_vec_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_mesh_to_xy_grid_2D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == mesh%name .and. Atlas( mi)%name_dst == grid%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_mesh_to_xy_grid( mesh, grid,Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_to_xy_grid_2D( mesh, grid, Atlas( mi), d_mesh_partial, d_grid_vec_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_to_xy_grid_2D

  subroutine map_from_mesh_to_xy_grid_3D( mesh, grid, d_mesh_partial, d_grid_vec_partial, method)
    ! Map a 3-D data field from an x/y-grid to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(:,:  ),          intent(in)    :: d_mesh_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_grid_vec_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_mesh_to_xy_grid_3D'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == mesh%name .and. Atlas( mi)%name_dst == grid%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_mesh_to_xy_grid( mesh, grid,Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_to_xy_grid_3D( mesh, grid, Atlas( mi), d_mesh_partial, d_grid_vec_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_to_xy_grid_3D

  subroutine map_from_mesh_to_xy_grid_2D_minval( mesh, grid, d_mesh_partial, d_grid_vec_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh.
    !
    ! For each grid cell, get the minimum value of all overlapping mesh vertices

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(:    ),          intent(in)    :: d_mesh_partial
    real(dp), dimension(:    ),          intent(out)   :: d_grid_vec_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_mesh_to_xy_grid_2D_minval'
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == mesh%name .and. Atlas( mi)%name_dst == grid%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_mesh_to_xy_grid( mesh, grid,Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_to_xy_grid_2D_minval( mesh, grid, Atlas( mi), d_mesh_partial, d_grid_vec_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_to_xy_grid_2D_minval

  ! From a mesh to a mesh
  subroutine map_from_mesh_to_mesh_with_reallocation_2D( mesh_src, mesh_dst, d_partial, method)
    ! Map a 2-D data field from a mesh to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    real(dp), dimension(:    ), allocatable, intent(inout) :: d_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_mesh_to_mesh_with_reallocation_2D'
    real(dp), dimension(:    ), allocatable            :: d_partial_new

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory for the remapped data field
    allocate( d_partial_new( mesh_dst%vi1: mesh_dst%vi2))

    ! Remap the data
    call map_from_mesh_to_mesh_2D( mesh_src, mesh_dst, d_partial, d_partial_new, method)

    ! Move allocation (and automatically also deallocate old memory, nice little bonus!)
    call MOVE_ALLOC( d_partial_new, d_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_to_mesh_with_reallocation_2D

  subroutine map_from_mesh_to_mesh_with_reallocation_3D( mesh_src, mesh_dst, d_partial, method)
    ! Map a 2-D data field from a mesh to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    real(dp), dimension(:,:  ), allocatable, intent(inout) :: d_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_mesh_to_mesh_with_reallocation_3D'
    real(dp), dimension(:,:  ), allocatable            :: d_partial_new

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory for the remapped data field
    allocate( d_partial_new( mesh_dst%vi1: mesh_dst%vi2, size( d_partial,2)))

    ! Remap the data
    call map_from_mesh_to_mesh_3D( mesh_src, mesh_dst, d_partial, d_partial_new, method)

    ! Move allocation (and automatically also deallocate old memory, nice little bonus!)
    call MOVE_ALLOC( d_partial_new, d_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_to_mesh_with_reallocation_3D

  subroutine map_from_mesh_to_mesh_2D( mesh_src, mesh_dst, d_src_partial, d_dst_partial, method)
    ! Map a 2-D data field from a mesh to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    real(dp), dimension(:    ),          intent(in)    :: d_src_partial
    real(dp), dimension(:    ),          intent(out)   :: d_dst_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_mesh_to_mesh_2D'
    logical                                            :: are_identical
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! if the two meshes are identical, the remapping operation is trivial
    call check_if_meshes_are_identical( mesh_src, mesh_dst, are_identical)
    if (are_identical) then
      d_dst_partial = d_src_partial
      call finalise_routine( routine_name)
      return
    end if

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == mesh_src%name .and. Atlas( mi)%name_dst == mesh_dst%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          if (present( method)) then
            SELECT CASE (method)
              CASE ('nearest_neighbour')
                call create_map_from_mesh_to_mesh_nearest_neighbour(      mesh_src, mesh_dst, Atlas( mi))
              CASE('trilin')
                call create_map_from_mesh_to_mesh_trilin(                 mesh_src, mesh_dst, Atlas( mi))
              CASE('2nd_order_conservative')
                call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
              CASE DEFAULT
                call crash('unknown remapping method "' // trim( method) // '"')
            end SELECT
          else
              call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
          end if
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_to_mesh_2D( mesh_src, mesh_dst, Atlas( mi), d_src_partial, d_dst_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_to_mesh_2D

  subroutine map_from_mesh_to_mesh_3D( mesh_src, mesh_dst, d_src_partial, d_dst_partial, method)
    ! Map a 3-D data field from a mesh to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    real(dp), dimension(:,:  ),          intent(in)    :: d_src_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_dst_partial
    character(len=*), optional,          intent(in)    :: method

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_mesh_to_mesh_3D'
    logical                                            :: are_identical
    integer                                            :: mi, mi_valid
    logical                                            :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! if the two meshes are identical, the remapping operation is trivial
    call check_if_meshes_are_identical( mesh_src, mesh_dst, are_identical)
    if (are_identical) then
      d_dst_partial = d_src_partial
      call finalise_routine( routine_name)
      return
    end if

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == mesh_src%name .and. Atlas( mi)%name_dst == mesh_dst%name) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          if (present( method)) then
            SELECT CASE (method)
              CASE ('nearest_neighbour')
                call create_map_from_mesh_to_mesh_nearest_neighbour(      mesh_src, mesh_dst, Atlas( mi))
              CASE ('trilin')
                call create_map_from_mesh_to_mesh_trilin(                 mesh_src, mesh_dst, Atlas( mi))
              CASE ('2nd_order_conservative')
                call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
              CASE DEFAULT
                call crash('unknown remapping method "' // trim( method) // '"')
            end SELECT
          else
              call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
          end if
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_to_mesh_3D( mesh_src, mesh_dst, Atlas( mi), d_src_partial, d_dst_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_to_mesh_3D

  ! From a vertical grid to another within the same mesh
  subroutine map_from_vertical_to_vertical_2D_ocean( mesh, vert_src, vert_dst, d_src_partial, d_dst_partial)
    ! Map mesh data between vertical grids

    ! Input variables:
    type(type_mesh),                                       intent(in)  :: mesh
    real(dp), dimension(:),                                intent(in)  :: vert_src
    real(dp), dimension(:),                                intent(in)  :: vert_dst
    real(dp), dimension(mesh%vi1:mesh%vi2,size(vert_src)), intent(in)  :: d_src_partial
    real(dp), dimension(mesh%vi1:mesh%vi2,size(vert_dst)), intent(out) :: d_dst_partial

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'map_from_vertical_to_vertical_2D_ocean'
    integer                            :: vi,k
    integer, dimension(:), allocatable :: z_mask_old, z_mask_new, mask_fill
    real(dp)                           :: z_floor
    real(dp)                           :: NaN
    real(dp), parameter                :: sigma = 4e4

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    NaN = -0.1234_dp

    ! allocate mask for valid points in a data column
    allocate( z_mask_old( size( vert_src)))
    allocate( z_mask_new( size( vert_dst)))

    do vi = mesh%vi1, mesh%vi2

      ! Determine local depth of the ocean floor, fill in both data masks
      if (d_src_partial( vi, size( vert_src)) == &
          d_src_partial( vi, size( vert_src))) then

        ! Ocean floor lies below the vertical limit of the provided data
        z_mask_old = 1
        z_floor = vert_src( size( vert_src)) + (vert_src( 2) - vert_src( 1))

      elseif (d_src_partial( vi,1) /= d_src_partial( vi,1)) then

        ! This grid cell isn't ocean at all
        z_mask_old = 0
        z_floor    = 0._dp
        NaN = d_src_partial( vi,1)

      else

        z_mask_old = 1
        k = size( vert_src)

        do while (d_src_partial( vi,k) /= d_src_partial( vi,k))
          z_mask_old( k) = 0
          z_floor = vert_src( k)
          k = k - 1
          NaN = d_src_partial( vi,k)
        end do

      end if

      z_mask_new = 0

      do k = 1, size(vert_dst)
        if (vert_dst( k) < z_floor) then
          z_mask_new = 1
        end if
      end do

      ! Regrid vertical column
      call remap_cons_2nd_order_1D( vert_src, z_mask_old, d_src_partial( vi,:), &
                                    vert_dst, z_mask_new, d_dst_partial( vi,:))

      ! Fill masked values with NaN
      do k = 1, size( vert_dst)
        if (z_mask_new( k) == 0) then
          d_dst_partial( vi,k) = NaN
        end if
      end do

    end do

    ! allocate mask for extrapolation
    allocate( mask_fill( mesh%vi1:mesh%vi2))

    ! Extrapolate into NaN areas independently for each layer
    do k = 1, size(vert_dst)
      ! Initialise assuming there's valid data everywhere
      mask_fill = 2
      ! Check this mesh layer for NaNs
      do vi = mesh%vi1, mesh%vi2
        if (d_dst_partial( vi,k) /= d_dst_partial( vi,k)) then
          ! if NaN, allow extrapolation here
          mask_fill( vi) = 1
        end if
      end do
      ! Fill NaN vertices within this layer
      call extrapolate_Gaussian( mesh, mask_fill, d_dst_partial(:,k), sigma)
    end do

    ! Clean up after yourself
    deallocate( z_mask_old)
    deallocate( z_mask_new)
    deallocate( mask_fill)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_vertical_to_vertical_2D_ocean

end module remapping_main
