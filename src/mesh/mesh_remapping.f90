module mesh_remapping

  ! Routines used in calculating and applying remapping operators between
  ! meshes, x/y-grids, and lon/lat-grids

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  use petscksp
  use mpi
  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, cerr, ierr, recv_status, sync
  use mpi_distributed_memory                                 , only: partition_list, gather_to_all_dp_1D
  use control_resources_and_error_messaging                  , only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use CSR_sparse_matrix_utilities                            , only: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, deallocate_matrix_CSR_dist, &
                                                                     add_empty_row_CSR_dist
  use petsc_basic                                            , only: mat_CSR2petsc, multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D, &
                                                                     MatDestroy, MatConvert, mat_petsc2CSR
  use grid_basic                                             , only: type_grid, calc_matrix_operators_grid, gather_gridded_data_to_master_dp_2D, &
                                                                     distribute_gridded_data_from_master_dp_2D, gather_gridded_data_to_master_dp_3D, &
                                                                     distribute_gridded_data_from_master_dp_3D, smooth_Gaussian_2D_grid, smooth_Gaussian_3D_grid
  use grid_lonlat_basic                                      , only: type_grid_lonlat
  use mesh_types                                             , only: type_mesh
  use math_utilities                                         , only: is_in_triangle, lies_on_line_segment, line_integral_xdy, line_integral_mxydx, &
                                                                     line_integral_xydy, crop_line_to_domain, segment_intersection, triangle_area
  use mesh_utilities                                         , only: is_in_Voronoi_cell, calc_Voronoi_cell, find_containing_vertex, find_containing_triangle, &
                                                                     find_shared_Voronoi_boundary, check_if_meshes_are_identical, set_border_vertices_to_interior_mean_dp_2D, &
                                                                     set_border_vertices_to_interior_mean_dp_3D, extrapolate_Gaussian
  use mesh_operators                                         , only: calc_all_matrix_operators_mesh

  implicit none

! ===== Derived types =====
! =========================

  type type_map
    ! A mapping object

    logical                                 :: is_in_use = .false.           ! Flag that indicates whether this map is in use
    character(len=256)                      :: name_src  = ''                ! Name of the source grid
    character(len=256)                      :: name_dst  = ''                ! Name of the destination grid
    character(len=256)                      :: method    = ''                ! Remapping method (nearest-neighbour, bilinear, 2-nd order conservative, etc.)
    type(tMat)                              :: M                             ! The actual operator matrix

  end type type_map

  type type_single_row_mapping_matrices
    ! Results from integrating around the border of a single grid cell

    integer                                 :: n_max
    integer                                 :: n
    integer,  dimension(:    ), allocatable :: index_left
    real(dp), dimension(:    ), allocatable :: LI_xdy, LI_mxydx, LI_xydy

  end type type_single_row_mapping_matrices

  ! The Atlas: the complete collection of all mapping objects.
  ! ==========================================================

  type(type_map), dimension(1000) :: Atlas

contains

! == High-level functions
! =======================

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

  ! Smoothing operations on the mesh
  subroutine smooth_Gaussian_2D( mesh, grid, d_mesh_partial, r)
    ! Use 2nd-order conservative remapping to map the 2-D data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(:    ),          intent(inout) :: d_mesh_partial
    real(dp),                            intent(in)    :: r

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'smooth_Gaussian_2D'
    real(dp), dimension(:    ), allocatable            :: d_grid_vec_partial

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    allocate( d_grid_vec_partial( grid%n_loc))

    ! Map data to the grid
    call map_from_mesh_to_xy_grid_2D( mesh, grid, d_mesh_partial, d_grid_vec_partial)

    ! Apply smoothing on the gridded data
    call smooth_Gaussian_2D_grid( grid, d_grid_vec_partial, r)

    ! Map data back to the mesh
    call map_from_xy_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Clean up after yourself
    deallocate( d_grid_vec_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_Gaussian_2D

  subroutine smooth_Gaussian_3D( mesh, grid, d_mesh_partial, r)
    ! Use 2nd-order conservative remapping to map the 3-D data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(:,:  ),          intent(inout) :: d_mesh_partial
    real(dp),                            intent(in)    :: r

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'smooth_Gaussian_3D'
    real(dp), dimension(:,:  ), allocatable            :: d_grid_vec_partial

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    allocate( d_grid_vec_partial( grid%n_loc, size( d_mesh_partial,2)))

    ! Map data to the grid
    call map_from_mesh_to_xy_grid_3D( mesh, grid, d_mesh_partial, d_grid_vec_partial)

    ! Apply smoothing on the gridded data
    call smooth_Gaussian_3D_grid( grid, d_grid_vec_partial, r)

    ! Map data back to the mesh
    call map_from_xy_grid_to_mesh_3D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Clean up after yourself
    deallocate( d_grid_vec_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_Gaussian_3D

  ! Clean up the Atlas after a mesh update
  subroutine clear_all_maps_involving_this_mesh( mesh)
    ! Clear all mapping objects involving a mesh by this name from the Atlas

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'map_from_xy_grid_to_mesh_2D'
    integer                                            :: mi, perr

    ! Add routine to path
    call init_routine( routine_name)

    do mi = 1, size( Atlas,1)
      if (Atlas( mi)%is_in_use .and. &
        (Atlas( mi)%name_src == mesh%name .or. &
         Atlas( mi)%name_dst == mesh%name .or. &
         Atlas( mi)%name_src == (trim( mesh%name) // '_triangles') .or. &
         Atlas( mi)%name_dst == (trim( mesh%name) // '_triangles'))) then
        ! This map involves the current mesh
        Atlas( mi)%is_in_use = .false.
        Atlas( mi)%name_src  = ''
        Atlas( mi)%name_dst  = ''
        call MatDestroy( Atlas( mi)%M, perr)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine clear_all_maps_involving_this_mesh

! == Apply existing mapping objects to remap data between grids
! =============================================================

  ! From an x/y-grid to a mesh
  subroutine apply_map_xy_grid_to_mesh_2D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 2-D data field from an x/y-grid to a mesh.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:    ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:    ),          intent(out)   :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_xy_grid_to_mesh_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nV_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_2D

  subroutine apply_map_xy_grid_to_mesh_3D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 3-D data field from an x/y-grid to a mesh.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:,:  ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_xy_grid_to_mesh_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nV_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc .or. &
      size( d_grid_vec_partial,2) /= size( d_mesh_partial,2)) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_3D

  ! From an x/y-grid to a mesh
  subroutine apply_map_xy_grid_to_mesh_triangles_2D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 2-D data field from an x/y-grid to a mesh triangles.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:    ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:    ),          intent(out)   :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_xy_grid_to_mesh_triangles_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nTri_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_triangles_2D

  subroutine apply_map_xy_grid_to_mesh_triangles_3D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 3-D data field from an x/y-grid to a mesh triangles.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:,:  ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_xy_grid_to_mesh_triangles_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nTri_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc .or. &
      size( d_grid_vec_partial,2) /= size( d_mesh_partial,2)) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_triangles_3D

  ! From a lon/lat-grid to a mesh
  subroutine apply_map_lonlat_grid_to_mesh_2D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 2-D data field from a lon/lat-grid to a mesh.

    ! In/output variables
    type(type_grid_lonlat),              intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:    ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:    ),          intent(out)   :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_lonlat_grid_to_mesh_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nV_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_lonlat_grid_to_mesh_2D

  subroutine apply_map_lonlat_grid_to_mesh_3D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 3-D data field from a lon/lat-grid to a mesh.

    ! In/output variables
    type(type_grid_lonlat),              intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:,:  ),          intent(in)    :: d_grid_vec_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_lonlat_grid_to_mesh_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nV_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc .or. &
      size( d_grid_vec_partial,2) /= size( d_mesh_partial,2)) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_lonlat_grid_to_mesh_3D

  ! From a mesh to an x/y-grid
  subroutine apply_map_mesh_to_xy_grid_2D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    ! Map a 2-D data field from a mesh to an x/y-grid.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:    ),          intent(in)    :: d_mesh_partial
    real(dp), dimension(:    ),          intent(out)   :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_mesh_to_xy_grid_2D'
    real(dp), dimension(:,:  ), allocatable            :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nV_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_mesh_partial, d_grid_vec_partial)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    if (par%master) then
      ! allocate memory for complete gridded data
      allocate( d_grid( grid%nx, grid%ny))
      ! Gather complete gridded data
      call gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:) = d_grid( 2        ,:)
      d_grid( grid%nx,:) = d_grid( grid%nx-1,:)
      d_grid( :,1      ) = d_grid( :,2        )
      d_grid( :,grid%ny) = d_grid( :,grid%ny-1)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      deallocate( d_grid)
    else ! if (par%master) then
      ! allocate zero memory for complete gridded data (only the master needs this)
      allocate( d_grid( 0,0))
      ! Gather complete gridded data
      call gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      deallocate( d_grid)
    end if ! if (par%master) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_to_xy_grid_2D

  subroutine apply_map_mesh_to_xy_grid_3D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    ! Map a 3-D data field from a mesh to an x/y-grid.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:,:  ),          intent(in)    :: d_mesh_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_mesh_to_xy_grid_3D'
    real(dp), dimension(:,:,:), allocatable            :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nV_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc .or. &
      size( d_mesh_partial,2) /= size( d_grid_vec_partial,2)) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_mesh_partial, d_grid_vec_partial)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    if (par%master) then
      ! allocate memory for complete gridded data
      allocate( d_grid( grid%nx, grid%ny, size( d_mesh_partial,2)))
      ! Gather complete gridded data
      call gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:,:) = d_grid( 2        ,:,:)
      d_grid( grid%nx,:,:) = d_grid( grid%nx-1,:,:)
      d_grid( :,1      ,:) = d_grid( :,2        ,:)
      d_grid( :,grid%ny,:) = d_grid( :,grid%ny-1,:)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_master_dp_3D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      deallocate( d_grid)
    else ! if (par%master) then
      ! allocate zero memory for complete gridded data (only the master needs this)
      allocate( d_grid( 0,0,0))
      ! Gather complete gridded data
      call gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_master_dp_3D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      deallocate( d_grid)
    end if ! if (par%master) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_to_xy_grid_3D

  subroutine apply_map_mesh_to_xy_grid_2D_minval( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    ! Map a 2-D data field from a mesh to an x/y-grid.
    !
    ! For each grid cell, get the minimum value of all overlapping mesh vertices

    ! In/output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_grid),                        intent(in)    :: grid
    type(type_map),                         intent(in)    :: map
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)    :: d_mesh_partial
    real(dp), dimension(grid%n1 :grid%n2 ), intent(out)   :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter                        :: routine_name = 'apply_map_mesh_to_xy_grid_2D_minval'
    real(dp), dimension(mesh%nV)                          :: d_mesh_tot
    type(type_sparse_matrix_CSR_dp)                       :: M_CSR
    integer                                               :: n,k1,k2,k,col,vi
    real(dp)                                              :: d_max, d_min

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_mesh_partial,1) /= mesh%nV_loc .or. size( d_grid_vec_partial,1) /= grid%n_loc) then
      call crash('data fields are the wrong size!')
    end if

    ! Gather global mesh data
    call gather_to_all_dp_1D( d_mesh_partial, d_mesh_tot)

    ! Convert mapping matrix from PETSc format to UFEMISM CSR format
    call mat_petsc2CSR( map%M, M_CSR)

    ! Find global maximum value of d
    d_max = maxval( d_mesh_partial)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Map data
    do n = grid%n1, grid%n2

      ! Initialise minimum as maximum
      d_min = d_max

      ! Loop over all mesh vertices that this grid cell overlaps with
      k1 = M_CSR%ptr( n)
      k2 = M_CSR%ptr( n+1)-1
      do k = k1, k2
        col = M_CSR%ind( k)
        ! This matrix row corresponds to this mesh vertex
        vi = mesh%n2vi( col)
        ! Update minimum value
        d_min = min( d_min, d_mesh_tot( vi))
      end do

      ! Fill into array
      d_grid_vec_partial( n) = d_min

    end do ! do n = grid%n1, grid%n2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_to_xy_grid_2D_minval

  ! From a mesh to a mesh
  subroutine apply_map_mesh_to_mesh_2D( mesh_src, mesh_dst, map, d_src_partial, d_dst_partial)
    ! Map a 2-D data field from a mesh to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:    ),          intent(in)    :: d_src_partial
    real(dp), dimension(:    ),          intent(out)   :: d_dst_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_mesh_to_mesh_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_src_partial,1) /= mesh_src%nV_loc .or. size( d_dst_partial,1) /= mesh_dst%nV_loc) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_src_partial, d_dst_partial)

    ! Set values of border vertices to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    call set_border_vertices_to_interior_mean_dp_2D( mesh_dst, d_dst_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_to_mesh_2D

  subroutine apply_map_mesh_to_mesh_3D( mesh_src, mesh_dst, map, d_src_partial, d_dst_partial)
    ! Map a 3-D data field from a mesh to a mesh.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    type(type_map),                      intent(in)    :: map
    real(dp), dimension(:,:  ),          intent(in)    :: d_src_partial
    real(dp), dimension(:,:  ),          intent(out)   :: d_dst_partial

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'apply_map_mesh_to_mesh_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d_src_partial,1) /= mesh_src%nV_loc .or. size( d_dst_partial,1) /= mesh_dst%nV_loc .or. &
      size( d_src_partial,2) /= size( d_dst_partial,2)) then
      call crash('data fields are the wrong size!')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_src_partial, d_dst_partial)

    ! Set values of border vertices to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    call set_border_vertices_to_interior_mean_dp_3D( mesh_dst, d_dst_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_to_mesh_3D

! == Create remapping objects
! ===========================

  subroutine create_map_from_xy_grid_to_mesh( grid, mesh, map)
    ! Create a new mapping object from an x/y-grid to a mesh.
    !
    ! By default uses 2nd-order conservative remapping.
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as having an area smaller
    !       than ten times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_xy_grid_to_mesh'
    type(PetscErrorCode)                               :: perr
    logical                                            :: count_coincidences
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR
    integer,  dimension(:    ), allocatable            :: mask_do_simple_average
    integer                                            :: vi
    real(dp), dimension( mesh%nC_mem,2)                :: Vor
    integer,  dimension( mesh%nC_mem  )                :: Vor_vi
    integer,  dimension( mesh%nC_mem  )                :: Vor_ti
    integer                                            :: nVor
    integer                                            :: vori1, vori2
    real(dp), dimension(2)                             :: p, q
    integer                                            :: k, i, j, kk, vj
    real(dp)                                           :: xl, xu, yl, yu
    real(dp), dimension(2)                             :: sw, se, nw, ne
    integer                                            :: vi_hint
    real(dp)                                           :: xmin, xmax, ymin, ymax
    integer                                            :: il, iu, jl, ju
    type(type_single_row_mapping_matrices)             :: single_row_Vor, single_row_grid
    type(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    type(tMat)                                         :: w0    , w1x    , w1y
    integer                                            :: row, k1, k2, col
    real(dp)                                           :: A_overlap_tot
    type(tMat)                                         :: grid_M_ddx, grid_M_ddy
    type(tMat)                                         :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = grid%name
    map%name_dst  = mesh%name
    map%method    = '2nd_order_conservative'

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix size
    nrows           = mesh%nV  ! to
    nrows_loc       = mesh%nV_loc
    ncols           = grid%n   ! from
    ncols_loc       = grid%n_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%A) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%A))) )

    call allocate_matrix_CSR_dist( A_xdy_a_g_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_a_g_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_a_g_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for single row results
    single_row_Vor%n_max = nnz_per_row_max
    single_row_Vor%n     = 0
    allocate( single_row_Vor%index_left( single_row_Vor%n_max))
    allocate( single_row_Vor%LI_xdy(     single_row_Vor%n_max))
    allocate( single_row_Vor%LI_mxydx(   single_row_Vor%n_max))
    allocate( single_row_Vor%LI_xydy(    single_row_Vor%n_max))

    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    allocate( single_row_grid%index_left( single_row_grid%n_max))
    allocate( single_row_grid%LI_xdy(     single_row_grid%n_max))
    allocate( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    allocate( single_row_grid%LI_xydy(    single_row_grid%n_max))

    allocate( mask_do_simple_average( mesh%nV), source = 0)

    ! Calculate line integrals around all Voronoi cells
    do row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      if (mesh%A( vi) < 10._dp * grid%dx**2) then
        ! This Voronoi cell is small enough to warrant a proper line integral

        mask_do_simple_average( vi) = 0

        ! Clean up single row results
        single_row_Vor%n          = 0
        single_row_Vor%index_left = 0
        single_row_Vor%LI_xdy     = 0._dp
        single_row_Vor%LI_mxydx   = 0._dp
        single_row_Vor%LI_xydy    = 0._dp

        ! Integrate around the complete Voronoi cell boundary
        call calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)
        do vori1 = 1, nVor
          vori2 = vori1 + 1
          if (vori2 > nVor) vori2 = 1
          p = Vor( vori1,:)
          q = Vor( vori2,:)
          count_coincidences = .true.
          call trace_line_grid( grid, p, q, single_row_Vor, count_coincidences)
        end do

        ! Safety
        if (single_row_Vor%n == 0) call crash('couldnt find any grid cells overlapping with this small Voronoi cell!')

        ! Next integrate around the grid cells overlapping with this Voronoi cell
        do k = 1, single_row_Vor%n

          ! Clean up single row results
          single_row_grid%n          = 0
          single_row_grid%index_left = 0
          single_row_grid%LI_xdy     = 0._dp
          single_row_grid%LI_mxydx   = 0._dp
          single_row_grid%LI_xydy    = 0._dp

          ! The grid cell
          col = single_row_Vor%index_left( k)
          i   = grid%n2ij( col,1)
          j   = grid%n2ij( col,2)

          xl = grid%x( i) - grid%dx / 2._dp
          xu = grid%x( i) + grid%dx / 2._dp
          yl = grid%y( j) - grid%dx / 2._dp
          yu = grid%y( j) + grid%dx / 2._dp

          sw = [xl,yl]
          nw = [xl,yu]
          se = [xu,yl]
          ne = [xu,yu]

          ! Integrate around the grid cell
          vi_hint = vi
          count_coincidences = .false.
          call trace_line_Vor( mesh, sw, se, single_row_grid, count_coincidences, vi_hint)
          call trace_line_Vor( mesh, se, ne, single_row_grid, count_coincidences, vi_hint)
          call trace_line_Vor( mesh, ne, nw, single_row_grid, count_coincidences, vi_hint)
          call trace_line_Vor( mesh, nw, sw, single_row_grid, count_coincidences, vi_hint)

          ! Safety
          if (single_row_grid%n == 0) call crash('couldnt find any Voronoi cells overlapping with this grid cell!')

          ! Add contribution for this particular triangle
          do kk = 1, single_row_grid%n
            vj = single_row_grid%index_left( kk)
            if (vj == vi) then
              ! Add contribution to this triangle
              single_row_Vor%LI_xdy(   k) = single_row_Vor%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
              single_row_Vor%LI_mxydx( k) = single_row_Vor%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
              single_row_Vor%LI_xydy(  k) = single_row_Vor%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
              exit
            end if
          end do ! do kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          call add_entry_CSR_dist( A_xdy_a_g_CSR  , row, col, single_row_Vor%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_a_g_CSR, row, col, single_row_Vor%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_a_g_CSR , row, col, single_row_Vor%LI_xydy(  k))

        end do ! do k = 1, single_row_Vor%n

      else ! if (mesh%A( vi) < 10._dp * grid%dx**2) then
        ! This Voronoi cell is big enough that we can just average over the grid cells it contains

        mask_do_simple_average( vi) = 1

        ! Clean up single row results
        single_row_Vor%n = 0

        ! Find the square of grid cells enveloping this Voronoi cell
        call calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)

        xmin = minval( Vor( 1:nVor,1))
        xmax = maxval( Vor( 1:nVor,1))
        ymin = minval( Vor( 1:nVor,2))
        ymax = maxval( Vor( 1:nVor,2))

        il = max( 1, min( grid%nx, 1 + floor( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        iu = max( 1, min( grid%nx, 1 + floor( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        jl = max( 1, min( grid%ny, 1 + floor( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx) ))
        ju = max( 1, min( grid%ny, 1 + floor( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

        ! Check which of the grid cells in this square lie inside the triangle
        do i = il, iu
        do j = jl, ju

          col = grid%ij2n( i,j)
          p   = [grid%x( i), grid%y( j)]

          if (is_in_Voronoi_cell( mesh, p, vi)) then
            ! This grid cell lies inside the triangle; add it to the single row
            single_row_Vor%n = single_row_Vor%n + 1
            single_row_Vor%index_left( single_row_Vor%n) = col
            single_row_Vor%LI_xdy(     single_row_Vor%n) = grid%dx**2
            single_row_Vor%LI_mxydx(   single_row_Vor%n) = grid%x( i) * grid%dx**2
            single_row_Vor%LI_xydy(    single_row_Vor%n) = grid%y( j) * grid%dx**2
          end if

        end do
        end do

        ! Safety
        if (single_row_Vor%n == 0) call crash('couldnt find any grid cells overlapping with this big Voronoi cell!')

        ! Add entries to the big matrices
        do k = 1, single_row_Vor%n
          col = single_row_Vor%index_left( k)
          call add_entry_CSR_dist( A_xdy_a_g_CSR  , vi, col, single_row_Vor%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_a_g_CSR, vi, col, single_row_Vor%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_a_g_CSR , vi, col, single_row_Vor%LI_xydy(  k))
        end do

      end if ! if (mesh%A( vi) < 4._dp * grid%dx**2) then

    end do ! do vi = mesh%vi1, mesh%vi2
    call sync

    ! Clean up after yourself
    deallocate( single_row_Vor%index_left )
    deallocate( single_row_Vor%LI_xdy     )
    deallocate( single_row_Vor%LI_mxydx   )
    deallocate( single_row_Vor%LI_xydy    )

    deallocate( single_row_grid%index_left )
    deallocate( single_row_grid%LI_xdy     )
    deallocate( single_row_grid%LI_mxydx   )
    deallocate( single_row_grid%LI_xydy    )

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    ! ==============================================================

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      k1 = A_xdy_a_g_CSR%ptr( row  )
      k2 = A_xdy_a_g_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_a_g_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_a_g_CSR%ind( k)
        call add_entry_CSR_dist( w0_CSR, row, col, A_xdy_a_g_CSR%val( k) / A_overlap_tot)
      end do

      if (mask_do_simple_average( vi) == 0) then
        ! For small vertices, include the gradient terms

        do k = k1, k2
          col = A_xdy_a_g_CSR%ind( k)
          ! Grid cell
          i = grid%n2ij( col,1)
          j = grid%n2ij( col,2)
          call add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_a_g_CSR%val( k) / A_overlap_tot) - (grid%x( i) * w0_CSR%val( k)))
          call add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_a_g_CSR%val(  k) / A_overlap_tot) - (grid%y( j) * w0_CSR%val( k)))
        end do

      else
        ! For large vertices, don't include the gradient terms

        call add_empty_row_CSR_dist( w1x_CSR, row)
        call add_empty_row_CSR_dist( w1y_CSR, row)

      end if ! if (mask_do_simple_average( vi) == 0) then

    end do ! do row = mesh%vi1, mesh%vi2
    call sync

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( A_xdy_a_g_CSR  )
    call deallocate_matrix_CSR_dist( A_mxydx_a_g_CSR)
    call deallocate_matrix_CSR_dist( A_xydy_a_g_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix

    call calc_matrix_operators_grid( grid, grid_M_ddx, grid_M_ddy)

    call MatDuplicate( w0, MAT_COPY_VALUES, map%M, perr)
    call MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatDestroy( grid_M_ddx    , perr)
    call MatDestroy( grid_M_ddy    , perr)
    call MatDestroy( w0            , perr)
    call MatDestroy( w1x           , perr)
    call MatDestroy( w1y           , perr)

    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    call MatDestroy( M1, perr)
    call MatDestroy( M2, perr)
    deallocate( mask_do_simple_average)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_xy_grid_to_mesh

  subroutine create_map_from_xy_grid_to_mesh_triangles( grid, mesh, map)
    ! Create a new mapping object from an x/y-grid to the triangles of a mesh.
    !
    ! By default uses 2nd-order conservative remapping.
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as having an area smaller
    !       than ten times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_xy_grid_to_mesh_triangles'
    type(PetscErrorCode)                               :: perr
    logical                                            :: count_coincidences
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR
    integer,  dimension(:    ), allocatable            :: mask_do_simple_average
    integer                                            :: ti
    integer                                            :: n1, n2, via, vib, vic
    real(dp), dimension(2)                             :: p, q
    integer                                            :: k, i, j, kk, tj
    real(dp)                                           :: xl, xu, yl, yu
    real(dp), dimension(2)                             :: sw, se, nw, ne
    integer                                            :: ti_hint
    real(dp)                                           :: xmin, xmax, ymin, ymax
    integer                                            :: il, iu, jl, ju
    type(type_single_row_mapping_matrices)             :: single_row_tri, single_row_grid
    real(dp), dimension(2)                             :: pa, pb, pc
    type(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    type(tMat)                                         :: w0    , w1x    , w1y
    integer                                            :: row, k1, k2, col
    real(dp)                                           :: A_overlap_tot
    type(tMat)                                         :: grid_M_ddx, grid_M_ddy
    type(tMat)                                         :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = grid%name
    map%name_dst  = trim( mesh%name) // '_triangles'
    map%method    = '2nd_order_conservative'

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix size
    nrows           = mesh%nTri  ! to
    nrows_loc       = mesh%nTri_loc
    ncols           = grid%n     ! from
    ncols_loc       = grid%n_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%TriA) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%TriA))) )

    call allocate_matrix_CSR_dist( A_xdy_b_g_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_b_g_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_b_g_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for single row results
    single_row_tri%n_max = nnz_per_row_max
    single_row_tri%n     = 0
    allocate( single_row_tri%index_left( single_row_tri%n_max))
    allocate( single_row_tri%LI_xdy(     single_row_tri%n_max))
    allocate( single_row_tri%LI_mxydx(   single_row_tri%n_max))
    allocate( single_row_tri%LI_xydy(    single_row_tri%n_max))

    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    allocate( single_row_grid%index_left( single_row_grid%n_max))
    allocate( single_row_grid%LI_xdy(     single_row_grid%n_max))
    allocate( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    allocate( single_row_grid%LI_xydy(    single_row_grid%n_max))

    allocate( mask_do_simple_average( mesh%nTri), source = 0)

    ! Calculate line integrals around all triangles
    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      if (mesh%TriA( ti) < 10._dp * grid%dx**2) then
        ! This triangle is small enough to warrant a proper line integral

        mask_do_simple_average( ti) = 0

        ! Clean up single row results
        single_row_tri%n          = 0
        single_row_tri%index_left = 0
        single_row_tri%LI_xdy     = 0._dp
        single_row_tri%LI_mxydx   = 0._dp
        single_row_tri%LI_xydy    = 0._dp

        ! Integrate around the triangle
        do n1 = 1, 3
          n2 = n1 + 1
          if (n2 == 4) n2 = 1
          via = mesh%Tri( ti,n1)
          vib = mesh%Tri( ti,n2)
          p = mesh%V( via,:)
          q = mesh%V( vib,:)
          count_coincidences = .true.
          call trace_line_grid( grid, p, q, single_row_tri, count_coincidences)
        end do

        ! Safety
        if (single_row_tri%n == 0) call crash('couldnt find any grid cells overlapping with this small triangle!')

        ! Next integrate around the grid cells overlapping with this triangle
        do k = 1, single_row_tri%n

          ! Clean up single row results
          single_row_grid%n          = 0
          single_row_grid%index_left = 0
          single_row_grid%LI_xdy     = 0._dp
          single_row_grid%LI_mxydx   = 0._dp
          single_row_grid%LI_xydy    = 0._dp

          ! The grid cell
          col = single_row_tri%index_left( k)
          i   = grid%n2ij( col,1)
          j   = grid%n2ij( col,2)

          xl = grid%x( i) - grid%dx / 2._dp
          xu = grid%x( i) + grid%dx / 2._dp
          yl = grid%y( j) - grid%dx / 2._dp
          yu = grid%y( j) + grid%dx / 2._dp

          sw = [xl,yl]
          nw = [xl,yu]
          se = [xu,yl]
          ne = [xu,yu]

          ! Integrate around the grid cell
          ti_hint = ti
          count_coincidences = .false.
          call trace_line_tri( mesh, sw, se, single_row_grid, count_coincidences, ti_hint)
          call trace_line_tri( mesh, se, ne, single_row_grid, count_coincidences, ti_hint)
          call trace_line_tri( mesh, ne, nw, single_row_grid, count_coincidences, ti_hint)
          call trace_line_tri( mesh, nw, sw, single_row_grid, count_coincidences, ti_hint)

          ! Safety
          if (single_row_grid%n == 0) call crash('couldnt find any triangles overlapping with this grid cell!')

          ! Add contribution for this particular triangle
          do kk = 1, single_row_grid%n
            tj = single_row_grid%index_left( kk)
            if (tj == ti) then
              ! Add contribution to this triangle
              single_row_tri%LI_xdy(   k) = single_row_tri%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
              single_row_tri%LI_mxydx( k) = single_row_tri%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
              single_row_tri%LI_xydy(  k) = single_row_tri%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
              exit
            end if
          end do ! do kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          call add_entry_CSR_dist( A_xdy_b_g_CSR  , row, col, single_row_tri%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_b_g_CSR, row, col, single_row_tri%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_b_g_CSR , row, col, single_row_tri%LI_xydy(  k))

        end do ! do k = 1, single_row_tri%n

      else ! if (mesh%TriA( ti) < 10._dp * grid%dx**2) then
        ! This triangle is big enough that we can just average over the grid cells it contains

        mask_do_simple_average( ti) = 1

        ! Clean up single row results
        single_row_tri%n = 0

        ! Find the square of grid cells enveloping this triangle
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        pa = mesh%V( via,:)
        pb = mesh%V( vib,:)
        pc = mesh%V( vic,:)

        xmin = min( min( pa( 1), pb( 1)), pc( 1))
        xmax = max( max( pa( 1), pb( 1)), pc( 1))
        ymin = min( min( pa( 2), pb( 2)), pc( 2))
        ymax = max( max( pa( 2), pb( 2)), pc( 2))

        il = max( 1, min( grid%nx, 1 + floor( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        iu = max( 1, min( grid%nx, 1 + floor( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        jl = max( 1, min( grid%ny, 1 + floor( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx) ))
        ju = max( 1, min( grid%ny, 1 + floor( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

        ! Check which of the grid cells in this square lie inside the triangle
        do i = il, iu
        do j = jl, ju

          col = grid%ij2n( i,j)
          p   = [grid%x( i), grid%y( j)]

          if (is_in_triangle( pa, pb, pc, p)) then
            ! This grid cell lies inside the triangle; add it to the single row
            single_row_tri%n = single_row_tri%n + 1
            single_row_tri%index_left( single_row_tri%n) = col
            single_row_tri%LI_xdy(     single_row_tri%n) = grid%dx**2
            single_row_tri%LI_mxydx(   single_row_tri%n) = grid%x( i) * grid%dx**2
            single_row_tri%LI_xydy(    single_row_tri%n) = grid%y( j) * grid%dx**2
          end if

        end do
        end do

        ! Safety
        if (single_row_tri%n == 0) call crash('couldnt find any grid cells overlapping with this big triangle!')

        ! Add entries to the big matrices
        do k = 1, single_row_tri%n
          col = single_row_tri%index_left( k)
          call add_entry_CSR_dist( A_xdy_b_g_CSR  , ti, col, single_row_tri%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_b_g_CSR, ti, col, single_row_tri%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_b_g_CSR , ti, col, single_row_tri%LI_xydy(  k))
        end do

      end if ! if (mesh%TriA( ti) < 10._dp * grid%dx**2) then

    end do ! do row = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    deallocate( single_row_tri%index_left )
    deallocate( single_row_tri%LI_xdy     )
    deallocate( single_row_tri%LI_mxydx   )
    deallocate( single_row_tri%LI_xydy    )

    deallocate( single_row_grid%index_left )
    deallocate( single_row_grid%LI_xdy     )
    deallocate( single_row_grid%LI_mxydx   )
    deallocate( single_row_grid%LI_xydy    )

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    ! ==============================================================

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      k1 = A_xdy_b_g_CSR%ptr( row  )
      k2 = A_xdy_b_g_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_b_g_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_b_g_CSR%ind( k)
        call add_entry_CSR_dist( w0_CSR, row, col, A_xdy_b_g_CSR%val( k) / A_overlap_tot)
      end do

      if (mask_do_simple_average( ti) == 0) then
        ! For small triangles, include the gradient terms

        do k = k1, k2
          col = A_xdy_b_g_CSR%ind( k)
          ! Grid cell
          i = grid%n2ij( col,1)
          j = grid%n2ij( col,2)
          call add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_b_g_CSR%val( k) / A_overlap_tot) - (grid%x( i) * w0_CSR%val( k)))
          call add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_b_g_CSR%val(  k) / A_overlap_tot) - (grid%y( j) * w0_CSR%val( k)))
        end do

      else
        ! For large triangles, don't include the gradient terms

        call add_empty_row_CSR_dist( w1x_CSR, row)
        call add_empty_row_CSR_dist( w1y_CSR, row)

      end if ! if (mask_do_simple_average( vi) == 0) then

    end do ! do row = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( A_xdy_b_g_CSR  )
    call deallocate_matrix_CSR_dist( A_mxydx_b_g_CSR)
    call deallocate_matrix_CSR_dist( A_xydy_b_g_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix

    call calc_matrix_operators_grid( grid, grid_M_ddx, grid_M_ddy)

    call MatDuplicate( w0, MAT_COPY_VALUES, map%M, perr)
    call MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatDestroy( grid_M_ddx    , perr)
    call MatDestroy( grid_M_ddy    , perr)
    call MatDestroy( w0            , perr)
    call MatDestroy( w1x           , perr)
    call MatDestroy( w1y           , perr)

    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    call MatDestroy( M1, perr)
    call MatDestroy( M2, perr)
    deallocate( mask_do_simple_average)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_xy_grid_to_mesh_triangles

  subroutine create_map_from_mesh_to_xy_grid( mesh, grid, map)
    ! Create a new mapping object from a mesh to an x/y-grid.
    !
    ! By default uses 2nd-order conservative remapping.
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as having an area smaller
    !       than ten times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is generally very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_grid),                     intent(in)    :: grid
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_mesh_to_xy_grid'
    type(PetscErrorCode)                               :: perr
    logical                                            :: count_coincidences
    integer,  dimension(:,:  ), allocatable            :: overlaps_with_small_triangle, containing_triangle
    integer                                            :: row, ti
    integer                                            :: via, vib, vic
    real(dp), dimension(2)                             :: pa, pb, pc
    real(dp)                                           :: xmin, xmax, ymin, ymax
    integer                                            :: il, iu, jl, ju
    integer                                            :: i, j, n_ext, ii, jj
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR
    type(type_single_row_mapping_matrices)             :: single_row_grid, single_row_Tri
    integer                                            :: ti_hint
    real(dp), dimension(2)                             :: p
    real(dp)                                           :: xl, xu, yl, yu
    real(dp), dimension(2)                             :: sw, se, nw, ne
    integer                                            :: k, kk, nn
    real(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy
    type(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    type(tMat)                                         :: w0    , w1x    , w1y
    integer                                            :: k1, k2, col
    real(dp)                                           :: A_overlap_tot
    type(tMat)                                         :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    type(tMat)                                         :: M1, M2
    integer                                            :: ncols_row
    integer,  dimension(:    ), allocatable            :: cols_row
    real(dp), dimension(:    ), allocatable            :: vals_row
    logical                                            :: has_value

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh%name
    map%name_dst  = grid%name
    map%method    = '2nd_order_conservative'

    ! == Find all grid cells that overlap with small triangles
    ! ========================================================

    allocate( overlaps_with_small_triangle( grid%nx, grid%ny), source = 0)
    allocate( containing_triangle(          grid%nx, grid%ny), source = 0)

    do row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)

      ! The square enveloping this triangle
      xmin = min( min( pa(1), pb(1)), pc(1))
      xmax = max( max( pa(1), pb(1)), pc(1))
      ymin = min( min( pa(2), pb(2)), pc(2))
      ymax = max( max( pa(2), pb(2)), pc(2))

      ! The square of grid cells enveloping this triangle
      il = 1 + floor( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx)
      iu = 1 + floor( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx)
      jl = 1 + floor( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx)
      ju = 1 + floor( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx)

      il = max( 1      , il - 1)
      iu = min( grid%nx, iu + 1)
      jl = max( 1      , jl - 1)
      ju = min( grid%ny, ju + 1)

      if (mesh%TriA( ti) < 10._dp * grid%dx**2) then
        ! This triangle is small; mark all grid cells it overlaps with

        ! Mark all these grid cells
        do i = il, iu
        do j = jl, ju
          overlaps_with_small_triangle( i,j) = 1
        end do
        end do

      else
        ! This triangle is large; mark all grid cells it contains

        ! Mark all these grid cells
        do i = il, iu
        do j = jl, ju
          p = [grid%x( i), grid%y( j)]
          if (is_in_triangle( pa, pb, pc, p)) then
            containing_triangle( i,j) = ti
          end if
        end do
        end do

      end if ! if (mesh%TriA( ti) < 4._dp * grid%dx**2) then

    end do

     ! Reduce results across the processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, overlaps_with_small_triangle, grid%nx * grid%ny, MPI_integer, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, containing_triangle         , grid%nx * grid%ny, MPI_integer, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Treat grid cells that possibly were not yet marked before
    do row = grid%n1, grid%n2

      ! Grid cell indices
      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)

      if (containing_triangle( i,j) == 0 .and. overlaps_with_small_triangle( i,j) == 0) then
        ! This grid cell does not overlap with a small triangle, but was not yet marked
        ! as being contained inside a large one; find the large triangle containing it.

        ! For efficiency, find the nearest grid cell that does list which large
        ! triangle contains it; use that as a hint for the triangle search
        n_ext = 0
        ti_hint = 0
        do while (ti_hint == 0)
          n_ext = n_ext+1
          ! Safety
          if (n_ext > max( grid%nx, grid%ny)) exit
          il = max( 1      , i-n_ext)
          iu = min( grid%nx, i+n_ext)
          jl = max( 1      , j-n_ext)
          ju = min( grid%ny, j+n_ext)
          do ii = il, iu
          do jj = jl, ju
            if (containing_triangle( ii,jj) > 0) then
              ti_hint = containing_triangle( ii,jj)
              exit
            end if
          end do
          if (ti_hint > 0) exit
          end do
        end do
        if (ti_hint == 0) ti_hint = 1

        ! Find the triangle containing this grid cell
        p = [max( mesh%xmin, min( mesh%xmax, grid%x( i) )), max( mesh%ymin, min( mesh%ymax, grid%y( j) ))]
        call find_containing_triangle( mesh, p, ti_hint)
        containing_triangle( i,j) = ti_hint

      end if

    end do
    call sync

    ! == Integrate around all grid cells that overlap with small triangles
    ! ====================================================================

    ! Initialise the three matrices using the native UFEMISM CSR-matrix format

    ! Matrix size
    nrows           = grid%n     ! to
    nrows_loc       = grid%n_loc
    ncols           = mesh%nTri  ! from
    ncols_loc       = mesh%nTri_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh%TriA) / (grid%dx**2)), &
                                    ceiling( 2._dp * (grid%dx**2) / minval( mesh%TriA))) )

    call allocate_matrix_CSR_dist( A_xdy_g_b_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_mxydx_g_b_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( A_xydy_g_b_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for single row results
    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    allocate( single_row_grid%index_left( single_row_grid%n_max))
    allocate( single_row_grid%LI_xdy(     single_row_grid%n_max))
    allocate( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    allocate( single_row_grid%LI_xydy(    single_row_grid%n_max))

    single_row_Tri%n_max = nnz_per_row_max
    single_row_Tri%n     = 0
    allocate( single_row_Tri%index_left( single_row_Tri%n_max))
    allocate( single_row_Tri%LI_xdy(     single_row_Tri%n_max))
    allocate( single_row_Tri%LI_mxydx(   single_row_Tri%n_max))
    allocate( single_row_Tri%LI_xydy(    single_row_Tri%n_max))

    ti_hint = 1

    do row = grid%n1, grid%n2

      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)
      p = [grid%x( i), grid%y( j)]

      ! The four sides of the grid cell
      xl = grid%x( i) - grid%dx / 2._dp
      xu = grid%x( i) + grid%dx / 2._dp
      yl = grid%y( j) - grid%dx / 2._dp
      yu = grid%y( j) + grid%dx / 2._dp

      ! if this grid cell lies entirely outside of the mesh domain, use
      ! nearest-neighbour extrapolation

      do while (xl <= mesh%xmin)
        i = i+1
        p( 1) = grid%x( i)
        xl = grid%x( i) - grid%dx / 2._dp
        xu = grid%x( i) + grid%dx / 2._dp
        if (i > grid%nx) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do
      do while (xu >= mesh%xmax)
        i = i-1
        p( 1) = grid%x( i)
        xl = grid%x( i) - grid%dx / 2._dp
        xu = grid%x( i) + grid%dx / 2._dp
        if (i < 1) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do
      do while (yl <= mesh%ymin)
        j = j+1
        p( 2) = grid%y( j)
        yl = grid%y( j) - grid%dx / 2._dp
        yu = grid%y( j) + grid%dx / 2._dp
        if (j > grid%ny) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do
      do while (yu >= mesh%ymax)
        j = j-1
        p( 2) = grid%y( j)
        yl = grid%y( j) - grid%dx / 2._dp
        yu = grid%y( j) + grid%dx / 2._dp
        if (j < 1) call crash('grid domain doesnt overlap with mesh domain at all!')
      end do

      if (overlaps_with_small_triangle( i,j) == 1) then
        ! This grid cell overlaps with a small triangle; integrate around it, and around
        ! all triangles overlapping with it

        sw = [xl, yl]
        nw = [xl, yu]
        se = [xu, yl]
        ne = [xu, yu]

        ! Clear the single row results
        single_row_grid%n          = 0
        single_row_grid%index_left = 0
        single_row_grid%LI_xdy     = 0._dp
        single_row_grid%LI_mxydx   = 0._dp
        single_row_grid%LI_xydy    = 0._dp

        ! Integrate over all four sides
        count_coincidences = .true.
        call trace_line_tri( mesh, sw, se, single_row_grid, count_coincidences, ti_hint)
        call trace_line_tri( mesh, se, ne, single_row_grid, count_coincidences, ti_hint)
        call trace_line_tri( mesh, ne, nw, single_row_grid, count_coincidences, ti_hint)
        call trace_line_tri( mesh, nw, sw, single_row_grid, count_coincidences, ti_hint)

        ! Next, integrate around all the triangles overlapping with this grid cell
        do k = 1, single_row_grid%n

          ti = single_row_grid%index_left( k)

          col = mesh%ti2n( ti)

          ! The three vertices spanning this triangle
          via = mesh%Tri( ti,1)
          vib = mesh%Tri( ti,2)
          vic = mesh%Tri( ti,3)

          pa  = mesh%V( via,:)
          pb  = mesh%V( vib,:)
          pc  = mesh%V( vic,:)

          ! Clear the single row results
          single_row_Tri%n = 0
          single_row_Tri%index_left = 0
          single_row_Tri%LI_xdy     = 0._dp
          single_row_Tri%LI_mxydx   = 0._dp
          single_row_Tri%LI_xydy    = 0._dp

          ! Integrate over all three triangle sides
          count_coincidences = .false.
          call trace_line_grid( grid, pa, pb, single_row_Tri, count_coincidences)
          call trace_line_grid( grid, pb, pc, single_row_Tri, count_coincidences)
          call trace_line_grid( grid, pc, pa, single_row_Tri, count_coincidences)

          ! Add contribution for this particular grid cell
          do kk = 1, single_row_Tri%n
            nn = single_row_Tri%index_left( kk)
            if (nn == row) then
              ! Add contribution to this triangle
              single_row_grid%LI_xdy(   k) = single_row_grid%LI_xdy(   k) + single_row_Tri%LI_xdy(   kk)
              single_row_grid%LI_mxydx( k) = single_row_grid%LI_mxydx( k) + single_row_Tri%LI_mxydx( kk)
              single_row_grid%LI_xydy(  k) = single_row_grid%LI_xydy(  k) + single_row_Tri%LI_xydy(  kk)
              exit
            end if
          end do ! do kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          call add_entry_CSR_dist( A_xdy_g_b_CSR  , row, col, single_row_grid%LI_xdy(   k))
          call add_entry_CSR_dist( A_mxydx_g_b_CSR, row, col, single_row_grid%LI_mxydx( k))
          call add_entry_CSR_dist( A_xydy_g_b_CSR , row, col, single_row_grid%LI_xydy(  k))

        end do ! do k = 1, single_row_grid%n

      else ! if (overlaps_with_small_triangle( i,j) == 1) then
        ! This grid cell does not overlap with a small triangle; use only the
        ! contribution from the nearest triangle

        ti_hint = containing_triangle( i,j)

        col = mesh%ti2n( ti_hint)

        LI_xdy   = grid%dx**2
        LI_mxydx = grid%dx**2 * grid%x( i)
        LI_xydy  = grid%dx**2 * grid%y( j)

        call add_entry_CSR_dist( A_xdy_g_b_CSR  , row, col, LI_xdy  )
        call add_entry_CSR_dist( A_mxydx_g_b_CSR, row, col, LI_mxydx)
        call add_entry_CSR_dist( A_xydy_g_b_CSR , row, col, LI_xydy )

      end if ! if (overlaps_with_small_triangle( i,j) == 1) then

    end do ! do n = grid%n1, grid%n2

    ! Clean up after yourself
    deallocate( overlaps_with_small_triangle)
    deallocate( containing_triangle         )

    deallocate( single_row_grid%index_left )
    deallocate( single_row_grid%LI_xdy     )
    deallocate( single_row_grid%LI_mxydx   )
    deallocate( single_row_grid%LI_xydy    )

    deallocate( single_row_Tri%index_left )
    deallocate( single_row_Tri%LI_xdy     )
    deallocate( single_row_Tri%LI_mxydx   )
    deallocate( single_row_Tri%LI_xydy    )

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    ! ==============================================================

    call allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    do row = grid%n1, grid%n2

      k1 = A_xdy_g_b_CSR%ptr( row  )
      k2 = A_xdy_g_b_CSR%ptr( row+1) - 1

      A_overlap_tot = sum( A_xdy_g_b_CSR%val( k1:k2))

      do k = k1, k2
        col = A_xdy_g_b_CSR%ind( k)
        ti = mesh%n2ti( col)
        call add_entry_CSR_dist( w0_CSR , row, col,  A_xdy_g_b_CSR%val(   k) / A_overlap_tot)
        call add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_g_b_CSR%val( k) / A_overlap_tot) - (mesh%TriGC( ti,1) * w0_CSR%val( k)))
        call add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_g_b_CSR%val(  k) / A_overlap_tot) - (mesh%TriGC( ti,2) * w0_CSR%val( k)))
      end do

    end do ! do row = mesh%vi1, mesh%vi2
    call sync

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( A_xdy_g_b_CSR  )
    call deallocate_matrix_CSR_dist( A_mxydx_g_b_CSR)
    call deallocate_matrix_CSR_dist( A_xydy_g_b_CSR )

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( w0_CSR , w0 )
    call mat_CSR2petsc( w1x_CSR, w1x)
    call mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix
    ! ==============================

    ! Convert matrices to PETSc format
    call mat_CSR2petsc( mesh%M_map_a_b, M_map_a_b)
    call mat_CSR2petsc( mesh%M_ddx_a_b, M_ddx_a_b)
    call mat_CSR2petsc( mesh%M_ddy_a_b, M_ddy_a_b)

    ! M = (w0 * M_map_a_b) + (w1x * M_ddx_a_b) + (w1y * M_ddy_a_b)
    call MatMatMult( w0,  M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, map%M, perr)
    call MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1   , perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2   , perr)
    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    call MatDestroy( w0            , perr)
    call MatDestroy( w1x           , perr)
    call MatDestroy( w1y           , perr)
    call MatDestroy( M_map_a_b     , perr)
    call MatDestroy( M_ddx_a_b     , perr)
    call MatDestroy( M_ddy_a_b     , perr)
    call MatDestroy( M1            , perr)
    call MatDestroy( M2            , perr)

    ! Safety: check if all grid cells get values
    ! ==========================================

    allocate( cols_row( nnz_per_row_max))
    allocate( vals_row( nnz_per_row_max))

    do row = grid%n1, grid%n2

      ! w0
      call MatGetRow( map%M, row-1, ncols_row, cols_row, vals_row, perr)
      if (ncols_row == 0) call crash('ncols == 0!')
      has_value = .false.
      do k = 1, ncols_row
        if (vals_row( k) /= 0._dp) has_value = .true.
      end do
      if (.not. has_value) call crash('only zeroes!')
      call MatRestoreRow( map%M, row-1, ncols_row, cols_row, vals_row, perr)

    end do
    call sync

    ! Clean up after yourself
    deallocate( cols_row)
    deallocate( vals_row)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_xy_grid

  subroutine create_map_from_lonlat_grid_to_mesh( grid, mesh, map)
    ! Create a new mapping object from a lon/lat-grid to a mesh.
    !
    ! By default uses bilinear interpolation.

    ! In/output variables:
    type(type_grid_lonlat),              intent(in)    :: grid
    type(type_mesh),                     intent(in)    :: mesh
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_lonlat_grid_to_mesh'
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: M_CSR
    integer                                            :: vi
    integer                                            :: il,iu,jl,ju
    real(dp)                                           :: wil,wiu,wjl,wju

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = grid%name
    map%name_dst  = mesh%name
    map%method    = 'bilin'

    ! == Initialise the mapping matrix using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix size
    nrows           = mesh%nV  ! to
    nrows_loc       = mesh%nV_loc
    ncols           = grid%n   ! from
    ncols_loc       = grid%n_loc
    nnz_per_row_max = 4
    nnz_est         = nnz_per_row_max * nrows
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Fill in the CSR matrix
    do vi = mesh%vi1, mesh%vi2

      ! Find enveloping lat-lon indices
      il  = max(1, min( grid%nlon-1, 1 + floor((mesh%lon( vi) - minval(grid%lon)) / (grid%lon(2) - grid%lon(1)))))
      iu  = il + 1
      wil = (grid%lon(iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
      wiu = 1._dp - wil

      ! Exception for pixels near the zero meridian
      if (mesh%lon( vi) < minval(grid%lon)) then
        il  = grid%nlon
        iu  = 1
        wil = (grid%lon( iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
        wiu = 1._dp - wil
      elseif (mesh%lon( vi) > maxval(grid%lon)) then
        il  = grid%nlon
        iu  = 1
        wiu = (mesh%lon( vi) - grid%lon( il)) / (grid%lon(2) - grid%lon(1))
        wil = 1._dp - wiu
      end if

      jl  = max(1, min( grid%nlat-1, 1 + floor((mesh%lat( vi) - minval(grid%lat)) / (grid%lat(2) - grid%lat(1)))))
      ju  = jl + 1
      wjl = (grid%lat( ju) - mesh%lat( vi)) / (grid%lat(2) - grid%lat(1))
      wju = 1 - wjl

      ! Add values to the CSR matrix
      call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( il,jl), wil * wjl)
      call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( il,ju), wil * wju)
      call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( iu,jl), wiu * wjl)
      call add_entry_CSR_dist( M_CSR, vi, grid%ij2n( iu,ju), wiu * wju)

    end do ! do vi = mesh%vi1, mesh%vi2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Clean up the Fortran versions
    call deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_lonlat_grid_to_mesh

  subroutine create_map_from_mesh_to_mesh_nearest_neighbour( mesh_src, mesh_dst, map)
    ! Create a new mapping object from a mesh to a mesh.
    !
    ! Uses nearest-neighbour interpolation.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_mesh_to_mesh_nearest_neighbour'
    integer                                            :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)                    :: M_CSR
    integer                                            :: row, vi_dst
    real(dp), dimension(2)                             :: p
    integer                                            :: vi_src, col

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = 'nearest_neighbour'

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    nrows_loc       = mesh_dst%nV_loc
    ncols           = mesh_src%nV   ! from
    ncols_loc       = mesh_src%nV_loc
    nnz_per_row_max = 1
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    vi_src = 1
    do row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      call find_containing_vertex( mesh_src, p, vi_src)

      col = mesh_src%vi2n( vi_src)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, row, col, 1._dp)

    end do ! do row = mesh_dst%vi1, mesh_dst%vi2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_nearest_neighbour

  subroutine create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, map)
    ! Create a new mapping object from a mesh to a mesh.
    !
    ! Uses trilinear interpolation.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_mesh_to_mesh_trilin'
    integer                                            :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)                    :: M_CSR
    integer                                            :: row, vi_dst
    real(dp), dimension(2)                             :: p
    integer                                            :: ti_src, via, vib, vic
    real(dp), dimension(2)                             :: pa, pb, pc
    real(dp)                                           :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc, wa, wb, wc
    integer                                            :: cola, colb, colc

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = 'trilin'

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    nrows_loc       = mesh_dst%nV_loc
    ncols           = mesh_src%nV   ! from
    ncols_loc       = mesh_src%nV_loc
    nnz_per_row_max = 3
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    ti_src = 1
    do row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      call find_containing_triangle( mesh_src, p, ti_src)

      ! Calculate the trilinear interpolation weights
      via = mesh_src%Tri( ti_src,1)
      vib = mesh_src%Tri( ti_src,2)
      vic = mesh_src%Tri( ti_src,3)

      pa  = mesh_src%V( via,:)
      pb  = mesh_src%V( vib,:)
      pc  = mesh_src%V( vic,:)

      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_abc = Atri_abp + Atri_bcp + Atri_cap

      wa = Atri_bcp / Atri_abc
      wb = Atri_cap / Atri_abc
      wc = Atri_abp / Atri_abc

      ! Matrix columns corresponding to these three vertices
      cola = mesh_src%vi2n( via)
      colb = mesh_src%vi2n( vib)
      colc = mesh_src%vi2n( vic)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, row, cola, wa)
      call add_entry_CSR_dist( M_CSR, row, colb, wb)
      call add_entry_CSR_dist( M_CSR, row, colc, wc)

    end do ! do row = mesh_dst%vi1, mesh_dst%vi2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_trilin

  subroutine create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, map)
    ! Create a new mapping object from a mesh to a mesh.
    !
    ! Uses 2nd-order conservative interpolation.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    type(type_map),                      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'create_map_from_mesh_to_mesh_2nd_order_conservative'
    type(PetscErrorCode)                               :: perr
    logical                                            :: count_coincidences
    integer                                            :: nnz_per_row_max
    type(tMat)                                         :: B_xdy_b_a  , B_mxydx_b_a  , B_xydy_b_a
    type(tMat)                                         :: B_xdy_a_b  , B_mxydx_a_b  , B_xydy_a_b
    type(tMat)                                         :: B_xdy_b_a_T, B_mxydx_b_a_T, B_xydy_b_a_T
    type(tMat)                                         :: w0, w1x, w1y
    integer                                            :: istart, iend, n, k, ti
    integer                                            :: ncols
    integer,  dimension(:    ), allocatable            :: cols
    real(dp), dimension(:    ), allocatable            :: vals, w0_row, w1x_row, w1y_row
    real(dp)                                           :: A_overlap_tot
    type(tMat)                                         :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    type(tMat)                                         :: M1, M2, M_cons_1st_order

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (map%is_in_use) call crash('this map is already in use!')

    ! == Initialise map metadata
    ! ==========================

    map%is_in_use = .true.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = '2nd_order_conservative'

    ! Integrate around the Voronoi cells of the destination mesh through the triangles of the source mesh
    count_coincidences = .true.
    call integrate_Voronoi_cells_through_triangles( mesh_dst, mesh_src, B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)

    ! Integrate around the triangles of the source mesh through the Voronoi cells of the destination mesh
    count_coincidences = .false.
    call integrate_triangles_through_Voronoi_cells( mesh_src, mesh_dst, B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)

    ! Transpose line integral matrices
    !if (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - transposing line integral matrices...'
    call MatCreateTranspose( B_xdy_b_a  , B_xdy_b_a_T  , perr)
    call MatCreateTranspose( B_mxydx_b_a, B_mxydx_b_a_T, perr)
    call MatCreateTranspose( B_xydy_b_a , B_xydy_b_a_T , perr)

    ! Combine line integrals around areas of overlap to get surface integrals over areas of overlap
    call MatAXPY( B_xdy_a_b  , 1._dp, B_xdy_b_a_T  , UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( B_mxydx_a_b, 1._dp, B_mxydx_b_a_T, UNKNOWN_NONZERO_PATTERN, perr)
    call MatAXPY( B_xydy_a_b , 1._dp, B_xydy_b_a_T , UNKNOWN_NONZERO_PATTERN, perr)

    call MatDestroy( B_xdy_b_a_T  , perr)
    call MatDestroy( B_mxydx_b_a_T, perr)
    call MatDestroy( B_xydy_b_a_T , perr)

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    call MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w0, perr)
    call MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1x, perr)
    call MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1y, perr)

    ! Estimate maximum number of non-zeros per row (i.e. maximum number of grid cells overlapping with a mesh triangle)
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_src%TriA) / minval( mesh_dst%A   )), &
                                    ceiling( 2._dp * maxval( mesh_dst%A   ) / minval( mesh_src%TriA))) )

    ! allocate memory for a single matrix row
    allocate( cols(    nnz_per_row_max))
    allocate( vals(    nnz_per_row_max))
    allocate( w0_row(  nnz_per_row_max))
    allocate( w1x_row( nnz_per_row_max))
    allocate( w1y_row( nnz_per_row_max))

    call MatGetOwnershipRange( B_xdy_a_b  , istart, iend, perr)

    do n = istart+1, iend ! +1 because PETSc indexes from 0

      ! w0
      call MatGetRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)
      A_overlap_tot = sum( vals( 1:ncols))
      do k = 1, ncols
        w0_row( k) = vals( k) / A_overlap_tot
        call MatSetValues( w0, 1, n-1, 1, cols( k), w0_row( k), INSERT_VALUES, perr)
      end do
      call MatRestoreRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)

      ! w1x
      call MatGetRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)
      do k = 1, ncols
        ti = cols( k)+1
        w1x_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,1) * w0_row( k))
        call MatSetValues( w1x, 1, n-1, 1, cols( k), w1x_row( k), INSERT_VALUES, perr)
      end do
      call MatRestoreRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)

      ! w1y
      call MatGetRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)
      do k = 1, ncols
        ti = cols( k)+1
        w1y_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,2) * w0_row( k))
        call MatSetValues( w1y, 1, n-1, 1, cols( k), w1y_row( k), INSERT_VALUES, perr)
      end do
      call MatRestoreRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)

    end do
    call MatAssemblyBegin( w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w0, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, perr)
    call MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, perr)
    call sync

    call MatDestroy( B_xdy_a_b  , perr)
    call MatDestroy( B_mxydx_a_b, perr)
    call MatDestroy( B_xydy_a_b , perr)

    ! == Calculate the remapping matrices
    ! ===================================

    ! Safety
    if (.not. allocateD( mesh_src%vi2n)) then
      call crash('matrix operators for mesh "' // trim( mesh_src%name) // '" have not been calculated!')
    end if

    ! Convert matrices to PETSc format
    call mat_CSR2petsc( mesh_src%M_map_a_b, M_map_a_b)
    call mat_CSR2petsc( mesh_src%M_ddx_a_b, M_ddx_a_b)
    call mat_CSR2petsc( mesh_src%M_ddy_a_b, M_ddy_a_b)

    ! 1st-order = w0 * map_a_b
    call MatMatMult( w0 , M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M_cons_1st_order, perr)

    ! 2nd-order = 1st-order + w1x * ddx_a_b + w1y * ddy_a_b
    call MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    call MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_real, M2, perr)

    call MatConvert( M_cons_1st_order, MATAIJ, MAT_INITIAL_MATRIX, map%M, perr)
    call MatAXPY( map%M, 1._dp, M1, DifFERENT_NONZERO_PATTERN, perr)
    call MatAXPY( map%M, 1._dp, M2, DifFERENT_NONZERO_PATTERN, perr)

    call MatDestroy( w0       , perr)
    call MatDestroy( w1x      , perr)
    call MatDestroy( w1y      , perr)
    call MatDestroy( M_map_a_b, perr)
    call MatDestroy( M_ddx_a_b, perr)
    call MatDestroy( M_ddy_a_b, perr)

    call MatDestroy( M1              , perr)
    call MatDestroy( M2              , perr)

    ! == Apply some final corrections
    ! ===============================

    call correct_mesh_to_mesh_map( mesh_src, mesh_dst, M_cons_1st_order, map%M)

    call MatDestroy( M_cons_1st_order, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_to_mesh_2nd_order_conservative

  subroutine correct_mesh_to_mesh_map( mesh_src, mesh_dst, M_cons_1st_order, M_cons_2nd_order)
    ! Apply some final corrections to the 2nd-order conservative mesh-to-mesh remapping operator:
    ! - set remapped data to zero on the domain border
    ! - use direct copying for identical vertices

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_src
    type(type_mesh),                     intent(in)    :: mesh_dst
    type(tMat),                          intent(in)    :: M_cons_1st_order
    type(tMat),                          intent(inout) :: M_cons_2nd_order

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'correct_mesh_to_mesh_map'
    integer                                            :: perr
    type(type_sparse_matrix_CSR_dp)                    :: M_cons_1st_order_CSR
    type(type_sparse_matrix_CSR_dp)                    :: M_cons_2nd_order_CSR
    integer                                            :: i, vi_dst, k1, k2, k
    integer                                            :: j, vi_src
    logical                                            :: do_direct_copy
    integer                                            :: vi_src_copy
    logical                                            :: Voronoi_cells_are_identical
    real(dp), dimension( mesh_src%nC_mem,2)            :: Vor_src   , Vor_dst
    integer,  dimension( mesh_src%nC_mem  )            :: Vor_src_vi, Vor_dst_vi
    integer,  dimension( mesh_src%nC_mem  )            :: Vor_src_ti, Vor_dst_ti
    integer                                            :: nVor_src  , nVor_dst
    integer                                            :: vori
    type(type_map)                                     :: map_trilin
    type(type_sparse_matrix_CSR_dp)                    :: M_trilin_CSR
    logical,  dimension( mesh_dst%nV)                  :: isgood_1st_order
    logical,  dimension( mesh_dst%nV)                  :: isgood_2nd_order
    integer                                            :: kk1,kk2,kk

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert matrices to CSR format for easier handling
    call mat_petsc2CSR( M_cons_1st_order, M_cons_1st_order_CSR)
    call mat_petsc2CSR( M_cons_2nd_order, M_cons_2nd_order_CSR)

    ! == Set to zero
    !
    ! 2nd-order conservative doesn't work all that well on the
    ! domain border; just set the result to zero there.

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      if (mesh_dst%VBI( vi_dst) > 0) then
        do k = k1, k2
          M_cons_2nd_order_CSR%val( k) = 0._dp
        end do
      end if

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! == Direct copying
    !
    ! With the new mesh generation code (of UFE2.0), many vertices away from the grounding line
    ! remain unchanged after a mesh update. The vertex-to-triangle-to-vertex remapping is slightly
    ! diffusive, so instead we can just copy data directly for those vertices.

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do_direct_copy = .false.
      vi_src_copy    = 0

      ! Loop over all src vertices contributing to this dst vertex
      do k = k1, k2

        j = M_cons_2nd_order_CSR%ind( k)
        vi_src = mesh_src%n2vi( j)

        if (norm2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%R( vi_dst) / 1E2_dp) then
          ! These vertices coincide; check if their Voronoi cells are identical

          Voronoi_cells_are_identical = .true.

          call calc_Voronoi_cell( mesh_src, vi_src, 0._dp, Vor_src, Vor_src_vi, Vor_src_ti, nVor_src)
          call calc_Voronoi_cell( mesh_dst, vi_dst, 0._dp, Vor_dst, Vor_dst_vi, Vor_dst_ti, nVor_dst)

          if (nVor_src /= nVor_dst) then
            Voronoi_cells_are_identical = .false.
          else
            do vori = 1, nVor_src
              if (norm2( Vor_src( vori,:) - Vor_dst( vori,:)) > mesh_dst%R( vi_dst) / 1E2_dp) then
                Voronoi_cells_are_identical = .false.
              end if
            end do
          end if

          if (Voronoi_cells_are_identical) then
            ! These two vertices have identical Voronoi cells; use direct copying
            do_direct_copy = .true.
            vi_src_copy    = vi_src
            exit
          end if ! if (Voronoi_cells_are_identical) then

        end if ! if (norm2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%tol_dist) then

      end do ! do k = k1, k2

      ! if a source vertex with an identical Voronoi cell to this dst vertex was
      ! found, copy data from that vertex directly
      if (do_direct_copy) then
        ! Loop over all src vertices contributing to this dst vertex; set all
        ! contributions to zero except the one we're copying (which is set to 1)

        do k = k1, k2

          j = M_cons_2nd_order_CSR%ind( k)
          vi_src = mesh_src%n2vi( j)

          if (vi_src == vi_src_copy) then
            M_cons_2nd_order_CSR%val( k) = 1._dp
          else
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end if

        end do ! do k = k1, k2

      end if ! if (do_direct_copy) then

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! == Remapping errors
    !
    ! On very rare occasions, the remapping operator is just wrong, likely due to round-off
    ! errors in determining if vertices coincide or not. Usually, the 1st-order operator
    ! is fine. if that one fails too, just replace the answer with trilinear interpolation.
    !
    ! Faulty operators can be detected by negative coefficients in the remapping matrix,
    ! which automatically violate conservation of extreme values.

    ! Calculate the trilinear interpolation operator to serve as a back-up
    call create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, map_trilin)
    call mat_petsc2CSR( map_trilin%M, M_trilin_CSR)

    ! Find faulty operators in the 1st-order conservative remapping operator
    isgood_1st_order = .true.
    do i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

      k1 = M_cons_1st_order_CSR%ptr( i)
      k2 = M_cons_1st_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do k = k1, k2
        if (M_cons_1st_order_CSR%val( k) < 0._dp) then
          isgood_1st_order( vi_dst) = .false.
        end if
      end do

    end do ! do i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

    ! Find faulty operators in the 2nd-order conservative remapping operator
    isgood_2nd_order = .true.
    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do k = k1, k2
        if (M_cons_2nd_order_CSR%val( k) < 0._dp) then
          isgood_2nd_order( vi_dst) = .false.
        end if
      end do

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Replace faulty operators in the 2nd-order conservative remapping operator

    do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      vi_dst = mesh_dst%n2vi( i)

      if (.not. isgood_2nd_order( vi_dst)) then
        ! Replace this faulty operator

        if (isgood_1st_order( vi_dst)) then
          ! Replace with the 1st-order conservative remapping operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          do k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end do

          ! then copy values from the 1st-order operator
          kk1 = M_cons_1st_order_CSR%ptr( i)
          kk2 = M_cons_1st_order_CSR%ptr( i+1) - 1
          do kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_cons_1st_order_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_cons_1st_order_CSR%val( kk)
          end do

        else ! if (isgood_1st_order( vi_dst)) then
          ! Replace with the trilinear interpolation operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          do k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          end do

          ! then copy values from the trilinear interpolation operator
          kk1 = M_trilin_CSR%ptr( i)
          kk2 = M_trilin_CSR%ptr( i+1) - 1
          do kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_trilin_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_trilin_CSR%val( kk)
          end do

        end if ! if (isgood_1st_order( vi_dst)) then

      end if ! if (.not. isgood_2nd_order( vi_dst)) then

    end do ! do i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Convert back to PETSc format
    call MatDestroy( M_cons_2nd_order, perr)
    call mat_CSR2petsc( M_cons_2nd_order_CSR, M_cons_2nd_order)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( M_cons_1st_order_CSR)
    call deallocate_matrix_CSR_dist( M_cons_2nd_order_CSR)
    call deallocate_matrix_CSR_dist( M_trilin_CSR)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine correct_mesh_to_mesh_map

! == Routines used in creating remapping matrices
! ===============================================

  ! Integrate around triangles/Voronoi cells through triangles/Voronoi cells
  subroutine integrate_triangles_through_Voronoi_cells( mesh_tri, mesh_Vor, B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)
    ! Integrate around the triangles of mesh_tri through the Voronoi cells of mesh_Vor

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_tri
    type(type_mesh),                     intent(in)    :: mesh_Vor
    type(tMat),                          intent(out)   :: B_xdy_b_a
    type(tMat),                          intent(out)   :: B_mxydx_b_a
    type(tMat),                          intent(out)   :: B_xydy_b_a
    logical,                             intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'integrate_triangles_through_Voronoi_cells'
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: B_xdy_b_a_CSR, B_mxydx_b_a_CSR, B_xydy_b_a_CSR
    type(type_single_row_mapping_matrices)             :: single_row
    integer                                            :: via, vib, vic, ti, vi_hint, k
    real(dp), dimension(2)                             :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_tri%nTri  ! to
    nrows_loc       = mesh_tri%nTri_loc
    ncols           = mesh_Vor%nV    ! from
    ncols_loc       = mesh_Vor%nV_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_tri%TriA) / minval( mesh_Vor%A   )), &
                                    ceiling( 2._dp * maxval( mesh_Vor%A   ) / minval( mesh_tri%TriA)) ))

    call allocate_matrix_CSR_dist( B_xdy_b_a_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_mxydx_b_a_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_xydy_b_a_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! == Trace all the line segments to fill the matrices
    ! ===================================================

    vi_hint = 1

    do ti = mesh_tri%ti1, mesh_tri%ti2

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! The three vertices spanning this triangle
      via = mesh_tri%Tri( ti,1)
      vib = mesh_tri%Tri( ti,2)
      vic = mesh_tri%Tri( ti,3)

      ! Integrate over the three triangle sides
      p = mesh_tri%V( via,:)
      q = mesh_tri%V( vib,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vib,:)
      q = mesh_tri%V( vic,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vic,:)
      q = mesh_tri%V( via,:)
      call trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      ! Add the results for this triangle to the sparse matrix
      if (single_row%n == 0) then
        call add_empty_row_CSR_dist( B_xdy_b_a_CSR  , ti)
        call add_empty_row_CSR_dist( B_mxydx_b_a_CSR, ti)
        call add_empty_row_CSR_dist( B_xydy_b_a_CSR , ti)
      else
        do k = 1, single_row%n
          call add_entry_CSR_dist( B_xdy_b_a_CSR  , ti, single_row%index_left( k), single_row%LI_xdy(   k))
          call add_entry_CSR_dist( B_mxydx_b_a_CSR, ti, single_row%index_left( k), single_row%LI_mxydx( k))
          call add_entry_CSR_dist( B_xydy_b_a_CSR , ti, single_row%index_left( k), single_row%LI_xydy(  k))
        end do
      end if

    end do ! do ti = mesh_tri%ti1, mesh_tri%ti2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( B_xdy_b_a_CSR  , B_xdy_b_a  )
    call mat_CSR2petsc( B_mxydx_b_a_CSR, B_mxydx_b_a)
    call mat_CSR2petsc( B_xydy_b_a_CSR , B_xydy_b_a )

    ! Clean up the Fortran versions
    call deallocate_matrix_CSR_dist( B_xdy_b_a_CSR  )
    call deallocate_matrix_CSR_dist( B_mxydx_b_a_CSR)
    call deallocate_matrix_CSR_dist( B_xydy_b_a_CSR )

    ! Clean up after yourself
    deallocate( single_row%index_left )
    deallocate( single_row%LI_xdy     )
    deallocate( single_row%LI_mxydx   )
    deallocate( single_row%LI_xydy    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_triangles_through_Voronoi_cells

  subroutine integrate_Voronoi_cells_through_triangles( mesh_Vor, mesh_tri, B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)
    ! Integrate around the grid cells of the grid through the triangles of the mesh

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh_Vor
    type(type_mesh),                     intent(in)    :: mesh_tri
    type(tMat),                          intent(out)   :: B_xdy_a_b
    type(tMat),                          intent(out)   :: B_mxydx_a_b
    type(tMat),                          intent(out)   :: B_xydy_a_b
    logical,                             intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'integrate_Voronoi_cells_through_triangles'
    integer                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    type(type_sparse_matrix_CSR_dp)                    :: B_xdy_a_b_CSR, B_mxydx_a_b_CSR, B_xydy_a_b_CSR
    type(type_single_row_mapping_matrices)             :: single_row
    integer                                            :: vi, vori1, vori2, k, ti_hint
    real(dp), dimension( mesh_Vor%nC_mem,2)            :: Vor
    integer,  dimension( mesh_Vor%nC_mem  )            :: Vor_vi
    integer,  dimension( mesh_Vor%nC_mem  )            :: Vor_ti
    integer                                            :: nVor
    real(dp), dimension(2)                             :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
    ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_Vor%nV    ! to
    nrows_loc       = mesh_Vor%nV_loc
    ncols           = mesh_tri%nTri  ! from
    ncols_loc       = mesh_tri%nTri_loc
    nnz_est         = 4 * max( nrows, ncols)
    nnz_est_proc    = ceiling( real( nnz_est, dp) / real( par%n, dp))
    nnz_per_row_max = max( 32, max( ceiling( 2._dp * maxval( mesh_tri%TriA) / minval( mesh_vor%A   )), &
                                    ceiling( 2._dp * maxval( mesh_vor%A   ) / minval( mesh_tri%TriA)) ))

    call allocate_matrix_CSR_dist( B_xdy_a_b_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_mxydx_a_b_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call allocate_matrix_CSR_dist( B_xydy_a_b_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! == Trace all the line segments to fill the matrices
    ! ===================================================

    ti_hint = 1

    do vi = mesh_Vor%vi1, mesh_Vor%vi2 ! +1 because PETSc indexes from 0

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! Integrate over the complete Voronoi cell boundary
      call calc_Voronoi_cell( mesh_Vor, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)
      do vori1 = 1, nVor
        vori2 = vori1 + 1
        if (vori2 > nVor) vori2 = 1
        p = Vor( vori1,:)
        q = Vor( vori2,:)
        call trace_line_tri( mesh_tri, p, q, single_row, count_coincidences, ti_hint)
      end do

      ! Add the results for this triangle to the sparse matrix
      if (single_row%n == 0) then
        call add_empty_row_CSR_dist( B_xdy_a_b_CSR  , vi)
        call add_empty_row_CSR_dist( B_mxydx_a_b_CSR, vi)
        call add_empty_row_CSR_dist( B_xydy_a_b_CSR , vi)
      else
        do k = 1, single_row%n
          call add_entry_CSR_dist( B_xdy_a_b_CSR  , vi, single_row%index_left( k), single_row%LI_xdy(   k))
          call add_entry_CSR_dist( B_mxydx_a_b_CSR, vi, single_row%index_left( k), single_row%LI_mxydx( k))
          call add_entry_CSR_dist( B_xydy_a_b_CSR , vi, single_row%index_left( k), single_row%LI_xydy(  k))
        end do
      end if

    end do ! do vi = mesh_Vor%vi1, mesh_Vor%vi2
    call sync

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( B_xdy_a_b_CSR  , B_xdy_a_b  )
    call mat_CSR2petsc( B_mxydx_a_b_CSR, B_mxydx_a_b)
    call mat_CSR2petsc( B_xydy_a_b_CSR , B_xydy_a_b )

    ! Clean up the Fortran versions
    call deallocate_matrix_CSR_dist( B_xdy_a_b_CSR  )
    call deallocate_matrix_CSR_dist( B_mxydx_a_b_CSR)
    call deallocate_matrix_CSR_dist( B_xydy_a_b_CSR )

    ! Clean up after yourself
    deallocate( single_row%index_left )
    deallocate( single_row%LI_xdy     )
    deallocate( single_row%LI_mxydx   )
    deallocate( single_row%LI_xydy    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_Voronoi_cells_through_triangles

  ! Line tracing algorithm through mesh triangles
  subroutine trace_line_tri( mesh, p, q, single_row, count_coincidences, ti_hint)
    ! Trace the line [pq] through the triangles of the mesh and calculate
    ! the three line integrals for the line segments inside the different triangles

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    logical,                             intent(in)    :: count_coincidences
    integer,                             intent(inout) :: ti_hint

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'trace_line_tri'
    real(dp), dimension(2)                             :: pp, qq
    logical                                            :: is_valid_line
    logical                                            :: finished
    integer                                            :: n_cycles
    integer                                            :: ti_in, vi_on, ei_on
    real(dp), dimension(2)                             :: p_next
    integer                                            :: ti_left
    logical                                            :: coincides
    real(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    call init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the mesh domain
    call crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    if (.not. is_valid_line) then
      ! [pq] doesn't pass through the mesh domain anywhere
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    call trace_line_tri_start( mesh, pp, ti_hint, ti_in, vi_on, ei_on)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      if     (ti_in  > 0) then
        ! p lies inside triangle ti_in
        call trace_line_tri_ti( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      elseif (vi_on  > 0) then
        ! p lies on vertex vi_on
        call trace_line_tri_vi( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      elseif (ei_on > 0) then
        ! p lies on edge ei_on
        call trace_line_tri_ei( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      else
        call crash('coincidence indicators dont make sense!')
      end if

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, mesh%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, mesh%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, mesh%tol_dist)

      ! Add them to the results structure
      if (norm2( p_next - pp) > mesh%tol_dist) then
        call add_integrals_to_single_row( single_row, ti_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      end if

      ! cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      if (n_cycles > mesh%nV) then
        call crash('trace_line_tri - iterative tracer got stuck!')
      end if

      ! Update ti_hint, for more efficiency
      ti_hint = ti_left

    end do ! do while (.not. finished)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine trace_line_tri

  subroutine trace_line_tri_start( mesh, p, ti_hint, ti_in, vi_on, ei_on)
    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside triangle ti_in, ...
    !    - lies on vertex vi_on, or...
    !    - lies on edge ei_on

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p
    integer,                             intent(inout) :: ti_hint
    integer,                             intent(out)   :: ti_in
    integer,                             intent(out)   :: vi_on
    integer,                             intent(out)   :: ei_on

    ! Local variables:
    integer                                            :: via, vib, vic
    real(dp), dimension(2)                             :: pa, pb, pc
    integer                                            :: vvi, vj, ei

    ! Initialise
    ti_in  = 0
    vi_on  = 0
    ei_on = 0

    ! Find the triangle containing p
    call find_containing_triangle( mesh, p, ti_hint)

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_hint,1)
    vib = mesh%Tri( ti_hint,2)
    vic = mesh%Tri( ti_hint,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Check if p lies on any of the three vertices
    if     (norm2( pa - p) < mesh%tol_dist) then
      ! p lies on via
      vi_on = via
      return
    elseif (norm2( pb - p) < mesh%tol_dist) then
      ! p lies on vib
      vi_on = vib
      return
    elseif (norm2( pc - p) < mesh%tol_dist) then
      ! p lies on vic
      vi_on = vic
      return
    end if

    ! Check if p lies on any of the three edges
    if     (lies_on_line_segment( pa, pb, p, mesh%tol_dist)) then
      ! p lies on the edge connecting via and vib
      do vvi = 1, mesh%nC( via)
        vj = mesh%C(  via,vvi)
        ei = mesh%VE( via,vvi)
        if (vj == vib) then
          ei_on = ei
          return
        end if
      end do
    elseif (lies_on_line_segment( pb, pc, p, mesh%tol_dist)) then
      ! p lies on the edge connecting vib and vic
      do vvi = 1, mesh%nC( vib)
        vj = mesh%C(  vib,vvi)
        ei = mesh%VE( vib,vvi)
        if (vj == vic) then
          ei_on = ei
          return
        end if
      end do
    elseif (lies_on_line_segment( pc, pa, p, mesh%tol_dist)) then
      ! p lies on the edge connecting vic and via
      do vvi = 1, mesh%nC( vic)
        vj = mesh%C(  vic,vvi)
        ei = mesh%VE( vic,vvi)
        if (vj == via) then
          ei_on = ei
          return
        end if
      end do
    end if

    ! if p lies not on the vertices or edges of the triangle, then it must lie inside of it
    ti_in = ti_hint

  end subroutine trace_line_tri_start

  subroutine trace_line_tri_ti(  mesh, p, q, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies inside triangle ti_in,
    ! find the point p_next where [pq] crosses into the next triangle.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(inout) :: ti_in
    integer,                             intent(inout) :: vi_on
    integer,                             intent(inout) :: ei_on
    integer,                             intent(out)   :: ti_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished

    ! Local variables:
    integer                                            :: via, vib, vic
    real(dp), dimension(2)                             :: pa, pb, pc
    integer                                            :: vvi, vj, ei
    real(dp), dimension(2)                             :: llis
    logical                                            :: do_cross

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_in,1)
    vib = mesh%Tri( ti_in,2)
    vic = mesh%Tri( ti_in,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Safety
    if (ti_in == 0 .or. vi_on > 0 .or. ei_on > 0) then
      call crash('trace_line_tri_ti - coincidence indicators dont make sense!')
    end if
    if (.not. is_in_triangle( pa, pb, pc, p)) then
      call crash('trace_line_tri_ti - p does not lie inside triangle ti_in!')
    end if

    ! Check if q lies inside the same triangle
    if (is_in_triangle( pa, pb, pc, q)) then
      ! q lies inside the same triangle
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on vertex via
    if (norm2( pa - q) < mesh%tol_dist) then
      ! q lies on vertex via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on vertex vib
    if (norm2( pb - q) < mesh%tol_dist) then
      ! q lies on vertex vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on vertex vic
    if (norm2( pc - q) < mesh%tol_dist) then
      ! q lies on vertex vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on edge via-vib
    if (lies_on_line_segment( pa, pb, q, mesh%tol_dist)) then
      ! q lies on edge via-vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on edge vib-vic
    if (lies_on_line_segment( pb, pc, q, mesh%tol_dist)) then
      ! q lies on edge vib-vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on edge vic-via
    if (lies_on_line_segment( pc, pa, q, mesh%tol_dist)) then
      ! q lies on edge vic-via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through via
    if (lies_on_line_segment( p, q, pa, mesh%tol_dist)) then
      ! [pq] passes through via
      p_next    = pa
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = via
      ei_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through vib
    if (lies_on_line_segment( p, q, pb, mesh%tol_dist)) then
      ! [pq] passes through vib
      p_next    = pb
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vib
      ei_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through vic
    if (lies_on_line_segment( p, q, pc, mesh%tol_dist)) then
      ! [pq] passes through vic
      p_next    = pc
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vic
      ei_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] crosses edge via-vib
    call segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
    if (do_cross) then
      ! [pq] crosses edge [via,vib]
      ! Find the edge connecting via and vib
      do vvi = 1, mesh%nC( via)
        vj = mesh%C(  via,vvi)
        ei = mesh%VE( via,vvi)
        if (vj == vib) then
          ei_on = ei
          exit
        end if
      end do
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] crosses edge vib-vic
    call segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
    if (do_cross) then
      ! [pq] crosses edge [vib,vic]
      ! Find the edge connecting vib and vic
      do vvi = 1, mesh%nC( vib)
        vj = mesh%C(  vib,vvi)
        ei = mesh%VE( vib,vvi)
        if (vj == vic) then
          ei_on = ei
          exit
        end if
      end do
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] crosses edge vic-via
    call segment_intersection( p, q, pc, pa, llis, do_cross, mesh%tol_dist)
    if (do_cross) then
      ! [pq] crosses edge [vic,via]
      ! Find the edge connecting vic and via
      do vvi = 1, mesh%nC( vic)
        vj = mesh%C(  vic,vvi)
        ei = mesh%VE( vic,vvi)
        if (vj == via) then
          ei_on = ei
          exit
        end if
      end do
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_tri_ti - reached the unreachable end of the subroutine!')

  end subroutine trace_line_tri_ti

  subroutine trace_line_tri_vi(  mesh, p, q, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies on vertex vi_on,
    ! find the point p_next where [pq] crosses into the next triangle.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(inout) :: ti_in
    integer,                             intent(inout) :: vi_on
    integer,                             intent(inout) :: ei_on
    integer,                             intent(out)   :: ti_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished

    ! Local variables:
    integer                                            :: via, vib, vic
    real(dp), dimension(2)                             :: pa, pb, pc, pv
    integer                                            :: vvi, vj, ei, n1, n2, n3, vti, ti
    real(dp), dimension(2)                             :: llis
    logical                                            :: do_cross

    ! Safety
    if (ti_in > 0 .or. vi_on == 0 .or. ei_on > 0) then
      call crash('trace_line_tri_vi - coincidence indicators dont make sense!')
    end if
    if (norm2( p - mesh%V( vi_on,:)) > mesh%tol_dist) then
      call crash('trace_line_tri_vi - p does not lie on vertex vi_on!')
    end if

    ! Check if q lies on any of the edges originating in this vertex
    do vvi = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,vvi)
      ei = mesh%VE( vi_on,vvi)
      pv = mesh%V( vj,:)
      if (norm2( mesh%V( vj,:) - q) < mesh%tol_dist .or. &
          lies_on_line_segment( p, pv, q, mesh%tol_dist)) then
        ! q lies on edge ei, connecting vi_on and vj
        if (mesh%EV( ei,1) == vi_on) then
          ti_left = mesh%ETri( ei,1)
        else
          ti_left = mesh%ETri( ei,2)
        end if
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        coincides = .true.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies inside any of the triangles surrounding vi_on
    do vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)
      if (is_in_triangle( pa, pb, pc, q) .or. &
          lies_on_line_segment( pa, pb, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pb, pc, q, mesh%tol_dist) .or. &
          lies_on_line_segment( pc, pa, q, mesh%tol_dist)) then
        ! q lies inside adjacent triangle ti
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = ti
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if [pq] passes through any of the neighbouring vertices
    do vvi = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,vvi)
      ei = mesh%VE( vi_on,vvi)
      pv = mesh%V( vj,:)
      if (lies_on_line_segment( p, q, pv, mesh%tol_dist)) then
        ! [pq] passes through neighbouring vertex vj, which is connected to vi_on by edge ei
        p_next    = pv
        ti_in     = 0
        vi_on     = vj
        ei_on     = 0
        if (mesh%EV( ei,1) == vi_on) then
          ti_left = mesh%ETri( ei,1)
        else
          ti_left = mesh%ETri( ei,2)
        end if
        coincides = .true.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] exits into any of the adjacent triangles
    do vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        n3 = n2 + 1
        if (n3 == 4) n3 = 1
        if (mesh%Tri( ti,n1) == vi_on) then
          vib = mesh%Tri( ti,n2)
          vic = mesh%Tri( ti,n3)
          pb  = mesh%V( vib,:)
          pc  = mesh%V( vic,:)
          ! Find the opposite triangle edge
          ei = 0
          do vvi = 1, mesh%nC( vib)
            vj = mesh%C( vib,vvi)
            if (vj == vic) then
              ei = mesh%VE( vib,vvi)
              exit
            end if
          end do
          call segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
          if (do_cross) then
            ! [pq] exits triangle ti through the opposite edge ei
            p_next    = llis
            ti_in     = 0
            vi_on     = 0
            ei_on     = ei
            ti_left   = ti
            coincides = .false.
            finished  = .false.
            return
          end if
        end if
      end do
    end do

    ! This point should not be reachable!
    call crash('trace_line_tri_vi - reached the unreachable end of the subroutine!')

  end subroutine trace_line_tri_vi

  subroutine trace_line_tri_ei( mesh, p, q, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies on edge ei,
    ! find the point p_next where [pq] crosses into the next triangle.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(inout) :: ti_in
    integer,                             intent(inout) :: vi_on
    integer,                             intent(inout) :: ei_on
    integer,                             intent(out)   :: ti_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished

    ! Local variables:
    integer                                            :: via, vib, vil, vir, til, tir
    real(dp), dimension(2)                             :: pa, pb, pl, pr
    integer                                            :: vvi, vj, ei
    real(dp), dimension(2)                             :: llis
    logical                                            :: do_cross

    ! Some more info about this edge
    via = mesh%EV(   ei_on,1)
    vib = mesh%EV(   ei_on,2)
    vil = mesh%EV(   ei_on,3)
    vir = mesh%EV(   ei_on,4)
    til = mesh%ETri( ei_on,1)
    tir = mesh%ETri( ei_on,2)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    if (vil > 0) pl  = mesh%V( vil,:)
    if (vir > 0) pr  = mesh%V( vir,:)

    ! Safety
    if (ti_in > 0 .or. vi_on > 0 .or. ei_on == 0) then
      call crash('trace_line_tri_ei - coincidence indicators dont make sense!')
    end if
    if (.not. lies_on_line_segment( pa, pb, p, mesh%tol_dist)) then
      call crash('trace_line_tri_ei - p does not lie on edge ei_on!')
    end if

    ! Check if q lies on the same edge in the direction of via
    if (lies_on_line_segment( p, pa, q, mesh%tol_dist)) then
      ! q lies on the same edge in the direction of via
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      ti_left   = tir
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the same edge in the direction of vib
    if (lies_on_line_segment( p, pb, q, mesh%tol_dist)) then
      ! q lies on the same edge in the direction of vib
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      ti_left   = til
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside either of the two adjacent triangles
    if (til > 0) then
      if (is_in_triangle( pa, pb, pl, q)) then
        ! q lies inside triangle til
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .true.
        return
      end if
    end if
    if (tir > 0) then
      if (is_in_triangle( pa, pr, pb, q)) then
        ! q lies inside triangle tir
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .true.
        return
      end if
    end if

    ! Check if [pq] passes through pa
    if (lies_on_line_segment( p, q, pa, mesh%tol_dist)) then
      ! [pq] passes through pa
      p_next    = pa
      ti_in     = 0
      vi_on     = via
      ei_on     = 0
      ti_left   = tir
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through pb
    if (lies_on_line_segment( p, q, pb, mesh%tol_dist)) then
      ! [pq] passes through pb
      p_next    = pb
      ti_in     = 0
      vi_on     = vib
      ei_on     = 0
      ti_left   = til
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through pl
    if (til > 0) then
      if (lies_on_line_segment( p, q, pl, mesh%tol_dist)) then
        ! [pq] passes through pl
        p_next    = pl
        ti_in     = 0
        vi_on     = vil
        ei_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] passes through pr
    if (tir > 0) then
      if (lies_on_line_segment( p, q, pr, mesh%tol_dist)) then
        ! [pq] passes through pr
        p_next    = pr
        ti_in     = 0
        vi_on     = vir
        ei_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [via,vil]
    if (til > 0) then
      call segment_intersection( p, q, pa, pl, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [via,vil]
        ! Find the edge connecting via and vil
        do vvi = 1, mesh%nC( via)
          vj = mesh%C(  via,vvi)
          ei = mesh%VE( via,vvi)
          if (vj == vil) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [vil,vib]
    if (til > 0) then
      call segment_intersection( p, q, pl, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [vil,vib]
        ! Find the edge connecting vil and vib
        do vvi = 1, mesh%nC( vil)
          vj = mesh%C(  vil,vvi)
          ei = mesh%VE( vil,vvi)
          if (vj == vib) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [via,vir]
    if (tir > 0) then
      call segment_intersection( p, q, pa, pr, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [via,vir]
        ! Find the edge connecting via and vir
        do vvi = 1, mesh%nC( via)
          vj  = mesh%C(    via,vvi)
          ei = mesh%VE( via,vvi)
          if (vj == vir) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses edge [vir,vib]
    if (tir > 0) then
      call segment_intersection( p, q, pr, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses edge [vir,vib]
        ! Find the edge connecting vir and vib
        do vvi = 1, mesh%nC( vir)
          vj = mesh%C(  vir,vvi)
          ei = mesh%VE( vir,vvi)
          if (vj == vib) then
            ei_on = ei
            exit
          end if
        end do
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .false.
        finished  = .false.
        return
      end if
    end if

    ! This point should not be reachable!
    call crash('trace_line_tri_ei - reached the unreachable end of the subroutine!')

  end subroutine trace_line_tri_ei

  ! Line tracing algorithm through mesh Voronoi cells
  subroutine trace_line_Vor( mesh, p, q, single_row, count_coincidences, vi_hint)
    ! Trace the line [pq] through the Voronoi cells of the mesh and calculate
    ! the three line integrals for the line segments inside the different Voronoi cells

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    logical,                             intent(in)    :: count_coincidences
    integer,                             intent(inout) :: vi_hint

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'trace_line_Vor'
    real(dp), dimension(2)                             :: pp,qq
    logical                                            :: is_valid_line
    logical                                            :: finished
    integer                                            :: n_cycles
    integer                                            :: vi_in, ti_on, ei_on
    real(dp), dimension(2)                             :: p_next
    integer                                            :: vi_left
    logical                                            :: coincides
    real(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    call init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the mesh domain
    call crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    if (.not.is_valid_line) then
      ! [pq] doesn't pass through the mesh domain anywhere
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    call trace_line_Vor_start( mesh, pp, vi_hint, vi_in, ti_on, ei_on)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not.finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      if     (vi_in  > 0) then
        ! p lies inside the Voronoi cell of vertex vi_in
        call trace_line_Vor_vi( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      elseif (ti_on  > 0) then
        ! p lies on the circumcentre of triangle ti_on
        call trace_line_Vor_ti( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      elseif (ei_on > 0) then
        ! p lies on the shared Voronoi cell boundary represented by edge ei_on
        call trace_line_Vor_ei( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      end if

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, mesh%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, mesh%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, mesh%tol_dist)

      ! Add them to the results structure
      if (norm2( p_next - pp) > mesh%tol_dist) then
        call add_integrals_to_single_row( single_row, vi_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      end if

      ! cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      if (n_cycles > mesh%nV) then
        call crash('trace_line_Vor - iterative tracer got stuck!')
      end if

      ! Update vi_hint, for more efficiency
      vi_hint = vi_left

    end do ! do while (.not.finished)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine trace_line_Vor

  subroutine trace_line_Vor_start( mesh, p, vi_hint, vi_in, ti_on, ei_on)
    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p
    integer,                             intent(inout) :: vi_hint
    integer,                             intent(out)   :: vi_in
    integer,                             intent(out)   :: ti_on
    integer,                             intent(out)   :: ei_on

    ! Local variables:
    integer                                            :: vti, ti, vei, ei
    real(dp), dimension(2)                             :: cc1, cc2

    ! Initialise
    vi_in = 0
    ti_on = 0
    ei_on = 0

    ! Find the vertex whose Voronoi cell contains p
    call find_containing_vertex( mesh, p, vi_hint)

    ! Check if p lies on any of the surrounding triangles' circumcentres
    do vti = 1, mesh%niTri( vi_hint)
      ti = mesh%iTri( vi_hint,vti)
      if (norm2( mesh%Tricc( ti,:) - p) < mesh%tol_dist) then
        ! p lies on the circumcentre of triangle ti
        ti_on = ti
        return
      end if
    end do

    ! Check if p lies on any of the shared Voronoi boundaries
    do vei = 1, mesh%nC( vi_hint)
      ei = mesh%VE( vi_hint,vei)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, p, mesh%tol_dist)) then
        ! p lies on the shared Voronoi cell boundary represented by edge ei
        ei_on = ei
        return
      end if
    end do

    ! if p lies not on the boundary of the Voronoi cell, then it must lie inside of it
    vi_in = vi_hint

  end subroutine trace_line_Vor_start

  subroutine trace_line_Vor_vi(  mesh, p, q, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies inside the Voronoi cell of vertex vi_in,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(inout) :: vi_in
    integer,                             intent(inout) :: ti_on
    integer,                             intent(inout) :: ei_on
    integer,                             intent(out)   :: vi_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished

    ! Local variables:
    integer                                            :: vti, ti, ei, vori, vorj, ci, vj
    real(dp), dimension(2)                             :: r, llis, pa, pb
    real(dp)                                           :: dx
    logical                                            :: do_cross
    real(dp), dimension( mesh%nC_mem,2)                :: Vor
    integer,  dimension( mesh%nC_mem  )                :: Vor_vi
    integer,  dimension( mesh%nC_mem  )                :: Vor_ti
    integer                                            :: nVor

    ! Safety
    if (vi_in == 0 .or. ti_on > 0 .or. ei_on > 0 .or. (.not. is_in_Voronoi_cell( mesh, p, vi_in))) then
      call crash('trace_line_Vor_vi - coincidence indicators dont make sense!')
    end if

    ! Check if q lies inside the same Voronoi cell
    if (is_in_Voronoi_cell( mesh, q, vi_in)) then
      ! q lies inside the same Voronoi cell
      p_next    = q
      vi_left   = vi_in
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on the boundary of this Voronoi cell
    dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 200
    call calc_Voronoi_cell( mesh, vi_in, dx, Vor, Vor_vi, Vor_ti, nVor)
    do vori = 1, nVor
      vorj = vori + 1
      if (vorj == nVor + 1) vorj = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori,:)
      pb = Vor( vorj,:)
      if (norm2( q - pa) < mesh%tol_dist .or. lies_on_line_segment( pa, pb, q, mesh%tol_dist)) then
        ! q lies on the boundary of the same Voronoi cell
        p_next    = q
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if [pq] passes through any of the surrounding triangles' circumcentres
    do vti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,vti)
      r  = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, r, mesh%tol_dist)) then
        ! [pq] passes through this triangle's circumcentre
        p_next    = mesh%Tricc( ti,:)
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] passes through any of the shared Voronoi boundaries
    do vori = 1, nVor
      vorj = vori + 1
      if (vorj == nVor + 1) vorj = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori,:)
      pb = Vor( vorj,:)
      ! The other vertex sharing this Voronoi cell boundary
      vj = Vor_vi( vori)
      ! The edge representing this shared Voronoi cell boundary
      ei = 0
      do ci = 1, mesh%nC( vi_in)
        if (mesh%C( vi_in,ci) == vj) then
          ei = mesh%VE( vi_in,ci)
          exit
        end if
      end do
      ! Safety
      if (ei == 0) call crash('couldnt find edge between vi and vj!')

      call segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes into the Voronoi cell of vj
        p_next    = llis
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! This point should not be reachable!
    call crash('trace_line_Vor_vi - reached the unreachable end of the subroutine!')

  end subroutine trace_line_Vor_vi

  subroutine trace_line_Vor_ti(  mesh, p, q, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies on the circumcentre of triangle ti_on,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(inout) :: vi_in
    integer,                             intent(inout) :: ti_on
    integer,                             intent(inout) :: ei_on
    integer,                             intent(out)   :: vi_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished

    ! Local variables:
    integer                                            :: via, vib, vic, vvi, vj, ei, acab, acbc, acca, tj
    real(dp), dimension(2)                             :: cc, cc1, cc2, llis
    logical                                            :: do_cross

    ! Safety
    if (vi_in > 0 .or. ti_on == 0 .or. ei_on > 0 .or. norm2( mesh%Tricc( ti_on,:) - p) > mesh%tol_dist) then
      call crash('trace_line_Vor_ti - coincidence indicators dont make sense!')
    end if

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_on,1)
    vib = mesh%Tri( ti_on,2)
    vic = mesh%Tri( ti_on,3)

    ! Find the three Voronoi cell boundaries that meet here
    acab = 0
    do vvi = 1, mesh%nC( via)
      vj = mesh%C(  via,vvi)
      ei = mesh%VE( via,vvi)
      if (vj == vib) then
        acab = ei
        exit
      end if
    end do
    acbc = 0
    do vvi = 1, mesh%nC( vib)
      vj = mesh%C(  vib,vvi)
      ei = mesh%VE( vib,vvi)
      if (vj == vic) then
        acbc = ei
        exit
      end if
    end do
    acca = 0
    do vvi = 1, mesh%nC( vic)
      vj = mesh%C(  vic,vvi)
      ei = mesh%VE( vic,vvi)
      if (vj == via) then
        acca = ei
        exit
      end if
    end do

    ! Check if q lies on the Voronoi cell boundary separating via from vib
    call find_shared_Voronoi_boundary( mesh, acab, cc1, cc2)
    if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .or. &
        norm2( cc1 - q) < mesh%tol_dist .or. &
        norm2( cc2 - q) < mesh%tol_dist) then
      ! q lies on the Voronoi cell boundary separating via from vib
      if (mesh%ETri( acab,1) == ti_on) then
        vi_left = mesh%EV( acab,2)
      else
        vi_left = mesh%EV( acab,1)
      end if
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the Voronoi cell boundary separating vib from vic
    call find_shared_Voronoi_boundary( mesh, acbc, cc1, cc2)
    if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .or. &
        norm2( cc1 - q) < mesh%tol_dist .or. &
        norm2( cc2 - q) < mesh%tol_dist) then
      ! q lies on the Voronoi cell boundary separating vib from vic
      if (mesh%ETri( acbc,1) == ti_on) then
        vi_left = mesh%EV( acbc,2)
      else
        vi_left = mesh%EV( acbc,1)
      end if
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the Voronoi cell boundary separating vic from via
    call find_shared_Voronoi_boundary( mesh, acca, cc1, cc2)
    if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .or. &
        norm2( cc1 - q) < mesh%tol_dist .or. &
        norm2( cc2 - q) < mesh%tol_dist) then
      ! q lies on the Voronoi cell boundary separating vic from via
      if (mesh%ETri( acca,1) == ti_on) then
        vi_left = mesh%EV( acca,2)
      else
        vi_left = mesh%EV( acca,1)
      end if
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside any of the three adjacent Voronoi cells
    if (is_in_Voronoi_cell( mesh, q, via)) then
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .false.
      finished  = .true.
      return
    end if
    if (is_in_Voronoi_cell( mesh, q, vib)) then
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .false.
      finished  = .true.
      return
    end if
    if (is_in_Voronoi_cell( mesh, q, vic)) then
      ! q lies inside the Voronoi cell of vic
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vic
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the circumcentre of any of the three neighbouring triangles
    tj = mesh%TriC( ti_on,1)
    if (tj > 0) then
      cc = mesh%Tricc( tj,:)
      if (lies_on_line_segment( p, q, cc, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = vic
        coincides = .true.
        finished  = .false.
        return
      end if
    end if

    tj = mesh%TriC( ti_on,2)
    if (tj > 0) then
      cc = mesh%Tricc( tj,:)
      if (lies_on_line_segment( p, q, cc, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = via
        coincides = .true.
        finished  = .false.
        return
      end if
    end if

    tj = mesh%TriC( ti_on,3)
    if (tj > 0) then
      cc = mesh%Tricc( tj,:)
      if (lies_on_line_segment( p, q, cc, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = vib
        coincides = .true.
        finished  = .false.
        return
      end if
    end if

    ! Check if [pq] crosses the boundary of the Voronoi cell of via
    do vvi = 1, mesh%nC( via)
      vj = mesh%C( via,vvi)
      if (vj == vib .or. vj == vic) cycle
      ei = mesh%VE( via,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = via
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] crosses the boundary of the Voronoi cell of vib
    do vvi = 1, mesh%nC( vib)
      vj = mesh%C( vib,vvi)
      if (vj == via .or. vj == vic) cycle
      ei = mesh%VE( vib,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vib
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if [pq] crosses the boundary of the Voronoi cell of vic
    do vvi = 1, mesh%nC( vic)
      vj = mesh%C( vic,vvi)
      if (vj == via .or. vj == vib) cycle
      ei = mesh%VE( vic,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vic
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! This point should not be reachable!
    call crash('trace_line_Vor_ti - reached the unreachable end of the subroutine!')

  end subroutine trace_line_Vor_ti

  subroutine trace_line_Vor_ei( mesh, p, q, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies on the shared Voronoi boundary represented by edge ei_on,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    real(dp), dimension(2),              intent(in)    :: p,q
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(inout) :: vi_in
    integer,                             intent(inout) :: ti_on
    integer,                             intent(inout) :: ei_on
    integer,                             intent(out)   :: vi_left
    logical,                             intent(out)   :: coincides
    logical,                             intent(out)   :: finished

    ! Local variables:
    integer                                            :: via, vib, vil, vir, til, tir, vvi, ei, vti, ti
    real(dp), dimension(2)                             :: cc1, cc2, ccl, ccr, llis
    logical                                            :: do_cross

    ! Find the endpoints of this shared Voronoi boundary
    call find_shared_Voronoi_boundary( mesh, ei_on, cc1, cc2)

    ! Safety
    if (vi_in > 0 .or. ti_on > 0 .or. ei_on == 0 .or. (.not. lies_on_line_segment( cc1, cc2, p, mesh%tol_dist))) then
      call crash('trace_line_Vor_ei - coincidence indicators dont make sense!')
    end if

    ! A bit more detail is needed
    via = mesh%EV(   ei_on,1)
    vib = mesh%EV(   ei_on,2)
    vil = mesh%EV(   ei_on,3)
    vir = mesh%EV(   ei_on,4)
    til = mesh%ETri( ei_on,1)
    tir = mesh%ETri( ei_on,2)

    if (til == 0) then
      ! Apparently ei lies on the domain border and has no triangle on its left-hand side
      ccr = cc1
      ccl = cc2
    elseif (tir == 0) then
      ! Apparently ei lies on the domain border and has no triangle on its right-hand side
      ccl = cc1
      ccr = cc2
    else
      ! ei lies in the interior and has triangles on both sides
      ccl = mesh%Tricc( til,:)
      ccr = mesh%Tricc( tir,:)
    end if

    ! Check if q coincides with ccl
    if (norm2( ccl - q) < mesh%tol_dist) then
      ! q coincides with ccl
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q coincides with ccr
    if (norm2( ccr - q) < mesh%tol_dist) then
      ! q coincides with ccr
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through ccl
    if (lies_on_line_segment( p, q, ccl, mesh%tol_dist)) then
      ! [pq] passes through ccl
      p_next    = ccl
      vi_in     = 0
      ti_on     = til
      ei_on     = 0
      vi_left   = via
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through ccr
    if (lies_on_line_segment( p, q, ccr, mesh%tol_dist)) then
      ! [pq] passes through ccr
      p_next    = ccr
      vi_in     = 0
      ti_on     = tir
      ei_on     = 0
      vi_left   = vib
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if q lies inside the Voronoi cell of via
    if (is_in_Voronoi_cell( mesh, q, via)) then
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the Voronoi cell of vib
    if (is_in_Voronoi_cell( mesh, q, vib)) then
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies on the circumcentre of any of the triangles surrounding via
    do vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      if (norm2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) then
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = via
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies on the circumcentre of any of the triangles surrounding vib
    do vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      if (norm2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) then
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = vib
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies on boundary of the Voronoi cell of via
    do vvi = 1, mesh%nC( via)
      ei = mesh%VE( via,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) then
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = via
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if q lies on boundary of the Voronoi cell of vib
    do vvi = 1, mesh%nC( vib)
      ei = mesh%VE( vib,vvi)
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      if (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) then
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = vib
        coincides = .false.
        finished  = .true.
        return
      end if
    end do

    ! Check if pq crosses the circumcentre of any of the triangles surrounding via
    do vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      cc1 = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of triangle ti
        p_next    = cc1
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        vi_left   = via
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if pq crosses the boundary of the Voronoi cell of via
    do vvi = 1, mesh%nC( via)
      ei = mesh%VE( via,vvi)
      if (ei == ei_on) cycle
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = via
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if pq crosses the circumcentre of any of the triangles surrounding vib
    do vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      cc1 = mesh%Tricc( ti,:)
      if (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) then
        ! [pq] passes through the circumcentre of triangle ti
        p_next    = cc1
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        vi_left   = vib
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! Check if pq crosses the boundary of the Voronoi cell of vib
    do vvi = 1, mesh%nC( vib)
      ei = mesh%VE( vib,vvi)
      if (ei == ei_on) cycle
      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      call segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      if (do_cross) then
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vib
        coincides = .false.
        finished  = .false.
        return
      end if
    end do

    ! This point should not be reachable!
    call crash('trace_line_Vor_ei - reached the unreachable end of the subroutine!')

  end subroutine trace_line_Vor_ei

  ! Line tracing algorithm through square grid cells
  subroutine trace_line_grid( grid, p, q, single_row, count_coincidences)
    ! Trace the line [pq] through the grid and calculate the three
    ! line integrals for the line segments inside the different grid cells.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(2),              intent(in)    :: p,q
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    logical,                             intent(in)    :: count_coincidences

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'trace_line_Vor'
    real(dp)                                           :: xmin, xmax, ymin, ymax
    real(dp), dimension(2)                             :: pp,qq
    logical                                            :: is_valid_line
    logical                                            :: finished
    integer                                            :: n_cycles
    integer,  dimension(2)                             :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2)                             :: p_next
    integer                                            :: n_left
    logical                                            :: coincides
    real(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    call init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the domain
    xmin = grid%xmin !- grid%dx / 2._dp
    xmax = grid%xmax !+ grid%dx / 2._dp
    ymin = grid%ymin !- grid%dx / 2._dp
    ymax = grid%ymax !+ grid%dx / 2._dp
    call crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, grid%tol_dist, pp, qq, is_valid_line)

    if (.not. is_valid_line) then
      ! [pq] doesn't pass through the domain anywhere
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside grid cell aij_in, ...
    !    - lies on the b-grid point bij_on, or...
    !    - lies on the edge cij_on
    call trace_line_grid_start( grid, pp, aij_in, bij_on, cxij_on, cyij_on)

    ! Iteratively trace the line through the mesh
    finished = .false.
    n_cycles = 0
    do while (.not. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      if     (aij_in(  1) > 0 .or. aij_in(  2) > 0) then
        ! p lies inside a-grid cell aij_in
        call trace_line_grid_a(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      elseif (bij_on(  1) > 0 .or. bij_on(  2) > 0) then
        ! p lies on b-grid point bij_on
        call trace_line_grid_b(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      elseif (cxij_on( 1) > 0 .or. cxij_on( 2) > 0) then
        ! p lies on cx-grid edge cxij_on
        call trace_line_grid_cx( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      elseif (cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
        ! p lies on cy-grid edge cyij_on
        call trace_line_grid_cy( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      else
        call crash('found no coincidence indicators!')
      end if

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, grid%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, grid%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, grid%tol_dist)

      ! Add them to the results structure
      if (norm2( p_next - pp) > grid%tol_dist) then
        call add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      end if

      ! cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      if (n_cycles > grid%n) then
        call crash('trace_line_grid - iterative tracer got stuck!')
      end if

    end do ! do while (.not. finished)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine trace_line_grid

  subroutine trace_line_grid_start( grid, p,    aij_in, bij_on, cxij_on, cyij_on)
    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside grid cell aij_in, ...
    !    - lies on the b-grid point bij_on, or...
    !    - lies on the edge cij_on

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(2),              intent(in)    :: p
    integer,  dimension(2),              intent(out)   :: aij_in, bij_on, cxij_on, cyij_on

    ! Local variables:
    integer                                            :: i,j
    real(dp)                                           :: xl,xu,yl,yu

    ! Initialise
    aij_in  = [0,0]
    bij_on  = [0,0]
    cxij_on = [0,0]
    cyij_on = [0,0]

    ! Find the grid cell containing p
    i = 1 + floor( (p(1) - grid%xmin + grid%dx / 2._dp) / grid%dx)
    j = 1 + floor( (p(2) - grid%ymin + grid%dx / 2._dp) / grid%dx)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! Check if p lies on either of the four surrounding b-grid points
    if     (i > 1       .and. j > 1       .and. abs( p(1) - xl) < grid%tol_dist .and. abs( p(2) - yl) < grid%tol_dist) then
      ! p coincides with the southwest corner
      bij_on = [i-1,j-1]
      return
    elseif (i > 1       .and. j < grid%ny .and. abs( p(1) - xl) < grid%tol_dist .and. abs( p(2) - yu) < grid%tol_dist) then
      ! p coincides with the northwest corner
      bij_on = [i-1,j  ]
      return
    elseif (i < grid%nx .and. j < 1       .and. abs( p(1) - xu) < grid%tol_dist .and. abs( p(2) - yl) < grid%tol_dist) then
      ! p coincides with the southeast corner
      bij_on = [i  ,j-1]
      return
    elseif (i < grid%nx .and. j < grid%ny .and. abs( p(1) - xu) < grid%tol_dist .and. abs( p(2) - yu) < grid%tol_dist) then
      ! p coincides with the northeast corner
      bij_on = [i  ,j  ]
      return
    end if

    ! Check if p lies on any of the four borders
    if     (i > 1       .and. abs( p(1) - xl) < grid%tol_dist) then
      ! p coincides with the western border
      cxij_on = [i-1,j  ]
      return
    elseif (i < grid%nx .and. abs( p(1) - xu) < grid%tol_dist) then
      ! p coincides with the eastern border
      cxij_on = [i  ,j  ]
      return
    elseif (j > 1       .and. abs( p(2) - yl) < grid%tol_dist) then
      ! p coincides with the southern border
      cyij_on = [i  ,j-1]
      return
    elseif (j < grid%ny .and. abs( p(2) - yu) < grid%tol_dist) then
      ! p coincides with the northern border
      cyij_on = [i  ,j  ]
      return
    end if

    ! p doesn't lie on the corners or borders, so it must lie inside the grid cell
    aij_in = [i,j]

  end subroutine trace_line_grid_start

  subroutine trace_line_grid_a(     grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies inside grid cell aij_in,
    ! find the point p_next where [pq] crosses into the next grid cell.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(2),              intent(in)    :: p,q
    integer,  dimension(2),              intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(out)   :: n_left
    logical,                             intent(out)   :: coincides, finished

    ! Local variables:
    integer                                            :: i,j
    real(dp)                                           :: xl,xu,yl,yu
    real(dp), dimension(2)                             :: sw,nw,se,ne
    logical                                            :: do_cross
    real(dp), dimension(2)                             :: llis

    ! Safety
    if ((aij_in( 1) == 0 .and. aij_in( 2) == 0) .or. cxij_on( 1) > 0 .or. cxij_on( 2) > 0 .or. &
        bij_on( 1) > 0 .or. bij_on( 2) > 0 .or. cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      call crash('trace_line_grid_a - coincidence indicators dont make sense!')
    end if

    i = aij_in( 1)
    j = aij_in( 2)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! More safety
    if (p(1) < xl .or. p(1) > xu .or. p(2) < yl .or. p(2) > yu) then
      call crash('trace_line_grid_a - coincidence indicators dont make sense!')
    end if

    ! Check if q lies inside the same grid cell
    if (q(1) >= xl - grid%tol_dist .and. &
        q(1) <= xu + grid%tol_dist .and. &
        q(2) >= yl - grid%tol_dist .and. &
        q(2) <= yu + grid%tol_dist) then
      ! q lies inside the same grid cell
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if pq passes through any of the four corners
    sw = [xl,yl]
    nw = [xl,yu]
    se = [xu,yl]
    ne = [xu,yu]

    if (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits this grid cell through the southwest corner
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    if (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits this grid cell through the northwest corner
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    if (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits this grid cell through the southeast corner
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    if (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits this grid cell through the northeast corner
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through any of the four boundaries
    call segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    call segment_intersection( p, q, se, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    call segment_intersection( p, q, sw, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    call segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits this grid cell through the northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j  ]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_a - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_a

  subroutine trace_line_grid_b(     grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on b-grid point bij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(2),              intent(in)    :: p,q
    integer,  dimension(2),              intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(out)   :: n_left
    logical,                             intent(out)   :: coincides, finished

    ! Local variables:
    integer                                            :: i,j
    real(dp)                                           :: x,y,xl,xu,yl,yu
    real(dp), dimension(2)                             :: sw,nw,se,ne,ww,ee,ss,nn
    logical                                            :: do_cross
    real(dp), dimension(2)                             :: llis

    ! Safety
    if (aij_in( 1) > 0 .or. aij_in( 2) > 0 .or. cxij_on( 1) > 0 .or. cxij_on( 2) > 0 .or. &
        (bij_on( 1) == 0 .and. bij_on( 2) == 0) .or. cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      call crash('trace_line_grid_b - coincidence indicators dont make sense!')
    end if

    i = bij_on( 1)
    j = bij_on( 2)

    ! The eight surrounding b-grid points spanning the four surrounding a-grid cells
    x  = grid%x( i) + grid%dx / 2._dp
    y  = grid%y( j) + grid%dx / 2._dp
    xl = x - grid%dx
    xu = x + grid%dx
    yl = y - grid%dx
    yu = y + grid%dx

    sw = [xl,yl]
    ww = [xl,y ]
    nw = [xl,yu]
    ss = [x ,yl]
    nn = [x ,yu]
    se = [xu,yl]
    ee = [xu,y ]
    ne = [xu,yu]

    ! More safety
    if (abs( p(1) - x) > grid%tol_dist .or. abs( p(2) - y) > grid%tol_dist) then
      call crash('trace_line_grid_b - coincidence indicators dont make sense!')
    end if

    ! Check if q lies on the cy-grid edge to the west
    if (q(1) < x + grid%tol_dist .and. q(1) > xl - grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the cy-grid edge to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the cy-grid edge to the east
    if (q(1) > x - grid%tol_dist .and. q(1) < xu + grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the cy-grid edge to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the cx-grid edge to the south
    if (q(2) < y + grid%tol_dist .and. q(2) > yl - grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the cx-grid edge to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the cx-grid edge to the north
    if (q(2) > y - grid%tol_dist .and. q(2) < yu + grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the cx-grid edge to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the northwest
    if (q(1) > xl - grid%tol_dist .and. q(1) < x  + grid%tol_dist .and. &
        q(2) > y  - grid%tol_dist .and. q(2) < yu + grid%tol_dist) then
      ! q lies inside the a-grid cell to the northwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the northeast
    if (q(1) > x  - grid%tol_dist .and. q(1) < xu + grid%tol_dist .and. &
        q(2) > y  - grid%tol_dist .and. q(2) < yu + grid%tol_dist) then
      ! q lies inside the a-grid cell to the northeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the southeast
    if (q(1) > x  - grid%tol_dist .and. q(1) < xu + grid%tol_dist .and. &
        q(2) > yl - grid%tol_dist .and. q(2) < y  + grid%tol_dist) then
      ! q lies inside the a-grid cell to the southeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the a-grid cell to the southwest
    if (q(1) > xl - grid%tol_dist .and. q(1) < x  + grid%tol_dist .and. &
        q(2) > yl - grid%tol_dist .and. q(2) < y  + grid%tol_dist) then
      ! q lies inside the a-grid cell to the southwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the west
    if (lies_on_line_segment( p, q, ww, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the north
    if (lies_on_line_segment( p, q, nn, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the west
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i  ,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the east
    if (lies_on_line_segment( p, q, ee, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i+1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the south
    if (lies_on_line_segment( p, q, ss, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northwest through its western boundary
    call segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northwest through its northern boundary
    call segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northwest through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northeast through its northern boundary
    call segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northeast through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j+1]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the northeast through its eastern boundary
    call segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the northeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southeast through its eastern boundary
    call segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southeast through its southern boundary
    call segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southeast through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southwest through its southern boundary
    call segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southwest through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the southwest through its western boundary
    call segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the southwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_b - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_b

  subroutine trace_line_grid_cx(    grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on cx-grid edge cxij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(2),              intent(in)    :: p,q
    integer,  dimension(2),              intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(out)   :: n_left
    logical,                             intent(out)   :: coincides, finished

    ! Local variables:
    integer                                            :: i,j
    real(dp)                                           :: x,yl,yu
    real(dp), dimension(2)                             :: sw,nw,se,ne,ss,nn
    logical                                            :: do_cross
    real(dp), dimension(2)                             :: llis

    ! Safety
    if (aij_in( 1) > 0 .or. aij_in( 2) > 0 .or. (cxij_on( 1) == 0 .and. cxij_on( 2) == 0) .or. &
        bij_on( 1) > 0 .or. bij_on( 2) > 0 .or. cyij_on( 1) > 0 .or. cyij_on( 2) > 0) then
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      call crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    end if

    i = cxij_on( 1)
    j = cxij_on( 2)

    ! This c-grid edge
    x  = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! The b-grid points spanning the two adjacent a-grid cells
    sw = [x - grid%dx, yl]
    nw = [x - grid%dx, yu]
    ss = [x          , yl]
    nn = [x          , yu]
    se = [x + grid%dx, yl]
    ne = [x + grid%dx, yu]

    ! More safety
    if (p(2) < yl .or. p(2) > yu .or. abs( p(1) - x) > grid%tol_dist) then
      call crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    end if

    ! Check if q lies on the same cx-grid cell in the southern direction
    if (q(2) < p(2) .and. q(2) >= yl - grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the same cx-grid cell in the southern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the same cx-grid cell in the northern direction
    if (q(2) > p(2) .and. q(2) <= yu + grid%tol_dist .and. abs( q(1) - x) < grid%tol_dist) then
      ! q lies on the same cx-grid cell in the northern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the west
    if (q(2) >= yl - grid%tol_dist .and. q(2) <= yu + grid%tol_dist .and. &
        q(1) >= x - grid%dx - grid%tol_dist .and. q(1) <= x + grid%tol_dist) then
      ! q lies inside the grid cell to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the east
    if (q(2) >= yl - grid%tol_dist .and. q(2) <= yu + grid%tol_dist .and. &
        q(1) <= x + grid%dx + grid%tol_dist .and. q(1) >= x - grid%tol_dist) then
      ! q lies inside the grid cell to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the south
    if (lies_on_line_segment( p, q, ss, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the north
    if (lies_on_line_segment( p, q, nn, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the north
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through the b-grid point to the northwest
    if (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through the b-grid point to the southwest
    if (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through the b-grid point to the northeast
    if (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i+1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through the b-grid point to the southeast
    if (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i+1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through its southern boundary
    call segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the west through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through its western boundary
    call segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the west through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the west through its northern boundary
    call segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the west through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through its northern boundary
    call segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the east through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through its eastern boundary
    call segment_intersection( p, q, ne, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the east through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the east through its southern boundary
    call segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the east through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_cx - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_cx

  subroutine trace_line_grid_cy(    grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on cy-grid edge cyij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    ! In/output variables
    type(type_grid),                     intent(in)    :: grid
    real(dp), dimension(2),              intent(in)    :: p,q
    integer,  dimension(2),              intent(inout) :: aij_in, bij_on, cxij_on, cyij_on
    real(dp), dimension(2),              intent(out)   :: p_next
    integer,                             intent(out)   :: n_left
    logical,                             intent(out)   :: coincides, finished

    ! Local variables:
    integer                                            :: i,j
    real(dp)                                           :: xl,xu,y
    real(dp), dimension(2)                             :: sw,nw,se,ne,ww,ee
    logical                                            :: do_cross
    real(dp), dimension(2)                             :: llis

    ! Safety
    if (aij_in( 1) > 0 .or. aij_in( 2) > 0 .or. cxij_on( 1) > 0 .or. cxij_on( 2) > 0 .or. &
        bij_on( 1) > 0 .or. bij_on( 2) > 0 .or. (cyij_on( 1) == 0 .and. cyij_on( 2) == 0)) then
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      call crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    end if

    i = cyij_on( 1)
    j = cyij_on( 2)

    ! This c-grid edge
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    y  = grid%y( j) + grid%dx / 2._dp

    ! The b-grid points spanning the two adjacent a-grid cells
    sw = [xl, y - grid%dx]
    se = [xu, y - grid%dx]
    ww = [xl, y          ]
    ee = [xu, y          ]
    nw = [xl, y + grid%dx]
    ne = [xu, y + grid%dx]

    ! More safety
    if (p(1) < xl .or. p(1) > xu .or. abs( p(2) - y) > grid%tol_dist) then
      call crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    end if

    ! Check if q lies on the same cy-grid cell in the western direction
    if (q(1) < p(1) .and. q(1) >= xl - grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the same cy-grid cell in the western direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies on the same cy-grid cell in the eastern direction
    if (q(1) > p(1) .and. q(1) <= xu + grid%tol_dist .and. abs( q(2) - y) < grid%tol_dist) then
      ! q lies on the same cy-grid cell in the eastern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .true.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the south
    if (q(1) >= xl - grid%tol_dist .and. q(1) <= xu + grid%tol_dist .and. &
        q(2) >= y - grid%dx - grid%tol_dist .and. q(2) <= y + grid%tol_dist) then
      ! q lies inside the grid cell to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if q lies inside the grid cell to the north
    if (q(1) >= xl - grid%tol_dist .and. q(1) <= xu + grid%tol_dist .and. &
        q(2) <= y + grid%dx + grid%tol_dist .and. q(2) >= y - grid%tol_dist) then
      ! q lies inside the grid cell to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .true.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the west
    if (lies_on_line_segment( p, q, ww, grid%tol_dist))  then
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] passes through the b-grid point to the east
    if (lies_on_line_segment( p, q, ee, grid%tol_dist)) then
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .true.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northwest
    if (lies_on_line_segment( p, q, nw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northeast
    if (lies_on_line_segment( p, q, ne, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southwest
    if (lies_on_line_segment( p, q, sw, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southeast
    if (lies_on_line_segment( p, q, se, grid%tol_dist)) then
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through its western boundary
    call segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the north through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through its northern boundary
    call segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the north through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the north through its eastern boundary
    call segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the north through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through its eastern boundary
    call segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the south through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through its southern boundary
    call segment_intersection( p, q, se, sw, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the south through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! Check if [pq] exits the a-grid cell to the south through its western boundary
    call segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    if (do_cross) then
      ! [pq] exits the a-grid cell to the south through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .false.
      finished  = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('trace_line_grid_cy - reached the unreachable end of the subroutine!')

  end subroutine trace_line_grid_cy

  ! Add the values for a single row of the three line-integral matrices
  subroutine add_integrals_to_single_row(  single_row, index_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
    ! Add the values for a single row of the three line-integral matrices

    ! In/output variables
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    integer,                                intent(in)    :: index_left
    real(dp),                               intent(in)    :: LI_xdy, LI_mxydx, LI_xydy
    logical,                                intent(in)    :: coincides, count_coincidences

    ! Local variables:
    logical                                            :: do_add_integrals, is_listed
    integer                                            :: i, i_add

    ! Check whether we actually need to add the line integrals
    do_add_integrals = .true.
    if (coincides .and. (.not. count_coincidences)) do_add_integrals = .false.

    ! Check if an entry from this left-hand vertex is already listed
    is_listed = .false.
    i_add     = 0

    do i = 1, single_row%n
      if (single_row%index_left( i) == index_left) then
        is_listed = .true.
        i_add     = i
        exit
      end if
    end do
    if (.not. is_listed) then
      single_row%n = single_row%n + 1
      i_add = single_row%n
    end if

    ! Add data
    single_row%index_left( i_add) = index_left
    if (do_add_integrals) then
      single_row%LI_xdy(   i_add) = single_row%LI_xdy(   i_add) + LI_xdy
      single_row%LI_mxydx( i_add) = single_row%LI_mxydx( i_add) + LI_mxydx
      single_row%LI_xydy(  i_add) = single_row%LI_xydy(  i_add) + LI_xydy
    end if

    ! if necessary, extend memory
    if (single_row%n > single_row%n_max - 10) call extend_single_row_memory( single_row, 100)

  end subroutine add_integrals_to_single_row

  subroutine extend_single_row_memory( single_row, n_extra)
    ! Extend memory for a single row of the three line-integral matrices

    ! In/output variables
    type(type_single_row_mapping_matrices), intent(inout) :: single_row
    integer,                                intent(in)    :: n_extra

    ! Local variables:
    integer                                            :: n
    integer,  dimension(:    ), allocatable            :: index_left_temp
    real(dp), dimension(:    ), allocatable            :: LI_xdy_temp, LI_mxydx_temp, LI_xydy_temp

    n = single_row%n

    ! allocate temporary memory
    allocate( index_left_temp( n))
    allocate( LI_xdy_temp(     n))
    allocate( LI_mxydx_temp(   n))
    allocate( LI_xydy_temp(    n))

    ! Copy data to temporary memory
    index_left_temp = single_row%index_left( 1:n)
    LI_xdy_temp     = single_row%LI_xdy(     1:n)
    LI_mxydx_temp   = single_row%LI_mxydx(   1:n)
    LI_xydy_temp    = single_row%LI_xydy(    1:n)

    ! deallocate memory
    deallocate( single_row%index_left)
    deallocate( single_row%LI_xdy    )
    deallocate( single_row%LI_mxydx  )
    deallocate( single_row%LI_xydy   )

    ! allocate new, extended memory
    single_row%n_max = single_row%n_max + n_extra
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    ! Copy data back from temporary memory
    single_row%index_left( 1:n) = index_left_temp
    single_row%LI_xdy(     1:n) = LI_xdy_temp
    single_row%LI_mxydx(   1:n) = LI_mxydx_temp
    single_row%LI_xydy(    1:n) = LI_xydy_temp

    ! deallocate temporary memory
    deallocate( index_left_temp)
    deallocate( LI_xdy_temp    )
    deallocate( LI_mxydx_temp  )
    deallocate( LI_xydy_temp   )

  end subroutine extend_single_row_memory

  ! Remapping of a 1-D variable (2nd-order conservative)
  subroutine remap_cons_2nd_order_1D( z_src, mask_src, d_src, z_dst, mask_dst, d_dst)
    ! 2nd-order conservative remapping of a 1-D variable
    !
    ! Used to remap ocean data from the provided vertical grid to the UFEMISM ocean vertical grid
    !
    ! Both z_src and z_dst can be irregular.
    !
    ! Both the src and dst data have a mask, with 0 indicating grid points where no data is defined.
    !
    ! This subroutine is serial, as it will be applied to single grid cells when remapping 3-D data fields,
    !   with the parallelisation being done by distributing the 2-D grid cells over the processes.

    ! In/output variables:
    real(dp), dimension(:    ),          intent(in)    :: z_src
    integer,  dimension(:    ),          intent(in)    :: mask_src
    real(dp), dimension(:    ),          intent(in)    :: d_src
    real(dp), dimension(:    ),          intent(in)    :: z_dst
    integer,  dimension(:    ),          intent(in)    :: mask_dst
    real(dp), dimension(:    ),          intent(out)   :: d_dst

    ! Local variables:
    logical                                            :: all_are_masked
    integer                                            :: nz_src, nz_dst
    integer                                            :: k
    real(dp), dimension(:    ), allocatable            :: ddz_src
    integer                                            :: k_src, k_dst
    real(dp)                                           :: zl_src, zu_src, zl_dst, zu_dst, z_lo, z_hi, z, d
    real(dp)                                           :: dz_overlap, dz_overlap_tot, d_int, d_int_tot
    real(dp)                                           :: dist_to_dst, dist_to_dst_min, max_dist
    integer                                            :: k_src_nearest_to_dst

    ! Initialise
    d_dst = 0._dp

    ! sizes
    nz_src = size( z_src,1)
    nz_dst = size( z_dst,1)

    ! Maximum distance on combined grids
    max_dist = maxval([ abs( z_src( nz_src) - z_src( 1)), &
                        abs( z_dst( nz_dst) - z_dst( 1)), &
                        abs( z_src( nz_src) - z_dst( 1)), &
                        abs( z_dst( nz_dst) - z_src( 1))])

    ! Exception for when the entire src field is masked
    all_are_masked = .true.
    do k = 1, nz_src
      if (mask_src( k) == 1) all_are_masked = .false.
    end do
    if (all_are_masked) return

    ! Exception for when the entire dst field is masked
    all_are_masked = .true.
    do k = 1, nz_dst
      if (mask_dst( k) == 1) all_are_masked = .false.
    end do
    if (all_are_masked) return

    ! Calculate derivative d_src/dz (one-sided differencing at the boundary, central differencing everywhere else)
    allocate( ddz_src( nz_src))
    do k = 2, nz_src-1
      ddz_src( k    ) = (d_src( k+1   ) - d_src( k-1     )) / (z_src( k+1   ) - z_src( k-1     ))
    end do
    ddz_src(  1     ) = (d_src( 2     ) - d_src( 1       )) / (z_src( 2     ) - z_src( 1       ))
    ddz_src(  nz_src) = (d_src( nz_src) - d_src( nz_src-1)) / (z_src( nz_src) - z_src( nz_src-1))

    ! Perform conservative remapping by finding regions of overlap
    ! between source and destination grid cells

    do k_dst = 1, nz_dst

      ! Skip masked grid cells
      if (mask_dst( k_dst) == 0) then
        d_dst( k_dst) = 0._dp
        cycle
      end if

      ! Find z range covered by this dst grid cell
      if (k_dst > 1) then
        zl_dst = 0.5_dp * (z_dst( k_dst - 1) + z_dst( k_dst))
      else
        zl_dst = z_dst( 1) - 0.5_dp * (z_dst( 2) - z_dst( 1))
      end if
      if (k_dst < nz_dst) then
        zu_dst = 0.5_dp * (z_dst( k_dst + 1) + z_dst( k_dst))
      else
        zu_dst = z_dst( nz_dst) + 0.5_dp * (z_dst( nz_dst) - z_dst( nz_dst-1))
      end if

      ! Find all overlapping src grid cells
      d_int_tot      = 0._dp
      dz_overlap_tot = 0._dp
      do k_src = 1, nz_src

        ! Skip masked grid cells
        if (mask_src( k_src) == 0) cycle

        ! Find z range covered by this src grid cell
        if (k_src > 1) then
          zl_src = 0.5_dp * (z_src( k_src - 1) + z_src( k_src))
        else
          zl_src = z_src( 1) - 0.5_dp * (z_src( 2) - z_src( 1))
        end if
        if (k_src < nz_src) then
          zu_src = 0.5_dp * (z_src( k_src + 1) + z_src( k_src))
        else
          zu_src = z_src( nz_src) + 0.5_dp * (z_src( nz_src) - z_src( nz_src-1))
        end if

        ! Find region of overlap
        z_lo = max( zl_src, zl_dst)
        z_hi = min( zu_src, zu_dst)
        dz_overlap = max( 0._dp, z_hi - z_lo)

        ! Calculate integral over region of overlap and add to sum
        if (dz_overlap > 0._dp) then
          z = 0.5_dp * (z_lo + z_hi)
          d = d_src( k_src) + ddz_src( k_src) * (z - z_src( k_src))
          d_int = d * dz_overlap

          d_int_tot      = d_int_tot      + d_int
          dz_overlap_tot = dz_overlap_tot + dz_overlap
        end if

      end do ! do k_src = 1, nz_src

      if (dz_overlap_tot > 0._dp) then
        ! Calculate dst value
        d_dst( k_dst) = d_int_tot / dz_overlap_tot
      else
        ! Exception for when no overlapping src grid cells were found; use nearest-neighbour extrapolation

        k_src_nearest_to_dst = 0._dp
        dist_to_dst_min      = max_dist
        do k_src = 1, nz_src
          if (mask_src( k_src) == 1) then
            dist_to_dst = abs( z_src( k_src) - z_dst( k_dst))
            if (dist_to_dst < dist_to_dst_min) then
              dist_to_dst_min      = dist_to_dst
              k_src_nearest_to_dst = k_src
            end if
          end if
        end do

        ! Safety
        if (k_src_nearest_to_dst == 0) then
          WRITE(0,*) '  remap_cons_2nd_order_1D - ERROR: couldnt find nearest neighbour on source grid!'
          call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        end if

        d_dst( k_dst) = d_src( k_src_nearest_to_dst)

      end if ! if (dz_overlap_tot > 0._dp) then

    end do ! do k_dst = 1, nz_dst

    ! Clean up after yourself
    deallocate( ddz_src)

  end subroutine remap_cons_2nd_order_1D

end module mesh_remapping
