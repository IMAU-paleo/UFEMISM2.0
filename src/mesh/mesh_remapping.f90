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
  use remapping_types, only: type_map, type_single_row_mapping_matrices
  use mesh_remapping_trace_line_basic, only: add_integrals_to_single_row
  use mesh_remapping_trace_line_triangles, only: trace_line_tri
  use mesh_remapping_trace_line_Voronoi_cells, only: trace_line_Vor
  use math_utilities, only: remap_cons_2nd_order_1D
  use mesh_remapping_create_map_grid_mesh, only: create_map_from_xy_grid_to_mesh, create_map_from_xy_grid_to_mesh_triangles, &
    create_map_from_mesh_to_xy_grid
  use mesh_remapping_create_map_gridlonlat_mesh, only: create_map_from_lonlat_grid_to_mesh
  use mesh_remapping_create_map_mesh_mesh, only: create_map_from_mesh_to_mesh_nearest_neighbour, &
    create_map_from_mesh_to_mesh_trilin, create_map_from_mesh_to_mesh_2nd_order_conservative

  implicit none

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

end module mesh_remapping
