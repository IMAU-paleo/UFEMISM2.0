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
  use mesh_remapping_trace_line_grid, only: trace_line_grid
  use mesh_remapping_trace_line_triangles, only: trace_line_tri
  use mesh_remapping_trace_line_Voronoi_cells, only: trace_line_Vor
  use math_utilities, only: remap_cons_2nd_order_1D

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

end module mesh_remapping
