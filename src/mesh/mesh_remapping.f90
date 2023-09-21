MODULE mesh_remapping

  ! Routines used in calculating and applying remapping operators between
  ! meshes, x/y-grids, and lon/lat-grids

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE mpi_distributed_memory                                 , ONLY: partition_list, gather_to_all_dp_1D
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, deallocate_matrix_CSR_dist, &
                                                                     add_empty_row_CSR_dist
  USE petsc_basic                                            , ONLY: mat_CSR2petsc, multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D, &
                                                                     MatDestroy, MatConvert, mat_petsc2CSR
  USE grid_basic                                             , ONLY: type_grid, calc_matrix_operators_grid, gather_gridded_data_to_master_dp_2D, &
                                                                     distribute_gridded_data_from_master_dp_2D, gather_gridded_data_to_master_dp_3D, &
                                                                     distribute_gridded_data_from_master_dp_3D, smooth_Gaussian_2D_grid, smooth_Gaussian_3D_grid
  USE grid_lonlat_basic                                      , ONLY: type_grid_lonlat
  USE mesh_types                                             , ONLY: type_mesh
  USE math_utilities                                         , ONLY: is_in_triangle, lies_on_line_segment, line_integral_xdy, line_integral_mxydx, &
                                                                     line_integral_xydy, crop_line_to_domain, segment_intersection, triangle_area
  USE mesh_utilities                                         , ONLY: is_in_Voronoi_cell, calc_Voronoi_cell, find_containing_vertex, find_containing_triangle, &
                                                                     find_shared_Voronoi_boundary, check_if_meshes_are_identical, set_border_vertices_to_interior_mean_dp_2D, &
                                                                     set_border_vertices_to_interior_mean_dp_3D, extrapolate_Gaussian
  USE mesh_operators                                         , ONLY: calc_all_matrix_operators_mesh

  IMPLICIT NONE

! ===== Derived types =====
! =========================

  TYPE type_map
    ! A mapping object

    LOGICAL                                 :: is_in_use = .FALSE.           ! Flag that indicates whether this map is in use
    CHARACTER(LEN=256)                      :: name_src  = ''                ! Name of the source grid
    CHARACTER(LEN=256)                      :: name_dst  = ''                ! Name of the destination grid
    CHARACTER(LEN=256)                      :: method    = ''                ! Remapping method (nearest-neighbour, bilinear, 2-nd order conservative, etc.)
    TYPE(tMat)                              :: M                             ! The actual operator matrix

  END TYPE type_map

  TYPE type_single_row_mapping_matrices
    ! Results from integrating around the border of a single grid cell

    INTEGER                                 :: n_max
    INTEGER                                 :: n
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: index_left
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: LI_xdy, LI_mxydx, LI_xydy

  END TYPE type_single_row_mapping_matrices

  ! The Atlas: the complete collection of all mapping objects.
  ! ==========================================================

  TYPE(type_map), DIMENSION(1000) :: Atlas

CONTAINS

! == High-level functions
! =======================

  ! From an x/y-grid to a mesh
  SUBROUTINE map_from_xy_grid_to_mesh_2D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_xy_grid_to_mesh_2D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == grid%name .AND. Atlas( mi)%name_dst == mesh%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_xy_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_xy_grid_to_mesh_2D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_xy_grid_to_mesh_2D

  SUBROUTINE map_from_xy_grid_to_mesh_3D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 3-D data field from an x/y-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_xy_grid_to_mesh_3D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == grid%name .AND. Atlas( mi)%name_dst == mesh%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_xy_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_xy_grid_to_mesh_3D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_xy_grid_to_mesh_3D

  ! From an x/y-grid to a mesh triangles
  SUBROUTINE map_from_xy_grid_to_mesh_triangles_2D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh triangles.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_xy_grid_to_mesh_triangles_2D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == grid%name .AND. Atlas( mi)%name_dst == (TRIM( mesh%name) // '_triangeles')) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_xy_grid_to_mesh_triangles( grid, mesh, Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_xy_grid_to_mesh_triangles_2D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_xy_grid_to_mesh_triangles_2D

  SUBROUTINE map_from_xy_grid_to_mesh_triangles_3D(     grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 3-D data field from an x/y-grid to a mesh triangles.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_xy_grid_to_mesh_triangles_3D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == grid%name .AND. Atlas( mi)%name_dst == (TRIM( mesh%name) // '_triangles')) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_xy_grid_to_mesh_triangles( grid, mesh, Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_xy_grid_to_mesh_triangles_3D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_xy_grid_to_mesh_triangles_3D

  ! From a lon/lat-grid to a mesh
  SUBROUTINE map_from_lonlat_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 2-D data field from a lon/lat-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_lonlat_grid_to_mesh_2D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == grid%name .AND. Atlas( mi)%name_dst == mesh%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_lonlat_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_lonlat_grid_to_mesh_2D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_lonlat_grid_to_mesh_2D

  SUBROUTINE map_from_lonlat_grid_to_mesh_3D( grid, mesh, d_grid_vec_partial, d_mesh_partial, method)
    ! Map a 3-D data field from a lon/lat-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_lonlat_grid_to_mesh_3D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == grid%name .AND. Atlas( mi)%name_dst == mesh%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_lonlat_grid_to_mesh( grid, mesh, Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_lonlat_grid_to_mesh_3D( grid, mesh, Atlas( mi), d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_lonlat_grid_to_mesh_3D

  ! From a mesh to an x/y-grid
  SUBROUTINE map_from_mesh_to_xy_grid_2D( mesh, grid, d_mesh_partial, d_grid_vec_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_mesh_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_grid_vec_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_mesh_to_xy_grid_2D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == mesh%name .AND. Atlas( mi)%name_dst == grid%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_mesh_to_xy_grid( mesh, grid,Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_mesh_to_xy_grid_2D( mesh, grid, Atlas( mi), d_mesh_partial, d_grid_vec_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_mesh_to_xy_grid_2D

  SUBROUTINE map_from_mesh_to_xy_grid_3D( mesh, grid, d_mesh_partial, d_grid_vec_partial, method)
    ! Map a 3-D data field from an x/y-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_mesh_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_grid_vec_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_mesh_to_xy_grid_3D'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == mesh%name .AND. Atlas( mi)%name_dst == grid%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_mesh_to_xy_grid( mesh, grid,Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_mesh_to_xy_grid_3D( mesh, grid, Atlas( mi), d_mesh_partial, d_grid_vec_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_mesh_to_xy_grid_3D

  SUBROUTINE map_from_mesh_to_xy_grid_2D_minval( mesh, grid, d_mesh_partial, d_grid_vec_partial, method)
    ! Map a 2-D data field from an x/y-grid to a mesh.
    !
    ! For each grid cell, get the minimum value of all overlapping mesh vertices

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_mesh_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_grid_vec_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_mesh_to_xy_grid_2D_minval'
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == mesh%name .AND. Atlas( mi)%name_dst == grid%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_mesh_to_xy_grid( mesh, grid,Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_mesh_to_xy_grid_2D_minval( mesh, grid, Atlas( mi), d_mesh_partial, d_grid_vec_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_mesh_to_xy_grid_2D_minval

  ! From a mesh to a mesh
  SUBROUTINE map_from_mesh_to_mesh_with_reallocation_2D( mesh_src, mesh_dst, d_partial, method)
    ! Map a 2-D data field from a mesh to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, INTENT(INOUT) :: d_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_mesh_to_mesh_with_reallocation_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_partial_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for the remapped data field
    ALLOCATE( d_partial_new( mesh_dst%vi1: mesh_dst%vi2))

    ! Remap the data
    CALL map_from_mesh_to_mesh_2D( mesh_src, mesh_dst, d_partial, d_partial_new, method)

    ! Move allocation (and automatically also deallocate old memory, nice little bonus!)
    CALL MOVE_ALLOC( d_partial_new, d_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_mesh_to_mesh_with_reallocation_2D

  SUBROUTINE map_from_mesh_to_mesh_with_reallocation_3D( mesh_src, mesh_dst, d_partial, method)
    ! Map a 2-D data field from a mesh to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(INOUT) :: d_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_mesh_to_mesh_with_reallocation_3D'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_partial_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for the remapped data field
    ALLOCATE( d_partial_new( mesh_dst%vi1: mesh_dst%vi2, SIZE( d_partial,2)))

    ! Remap the data
    CALL map_from_mesh_to_mesh_3D( mesh_src, mesh_dst, d_partial, d_partial_new, method)

    ! Move allocation (and automatically also deallocate old memory, nice little bonus!)
    CALL MOVE_ALLOC( d_partial_new, d_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_mesh_to_mesh_with_reallocation_3D

  SUBROUTINE map_from_mesh_to_mesh_2D( mesh_src, mesh_dst, d_src_partial, d_dst_partial, method)
    ! Map a 2-D data field from a mesh to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_src_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_dst_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_mesh_to_mesh_2D'
    LOGICAL                                            :: are_identical
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the two meshes are identical, the remapping operation is trivial
    CALL check_if_meshes_are_identical( mesh_src, mesh_dst, are_identical)
    IF (are_identical) THEN
      d_dst_partial = d_src_partial
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == mesh_src%name .AND. Atlas( mi)%name_dst == mesh_dst%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          IF (PRESENT( method)) THEN
            SELECT CASE (method)
              CASE ('nearest_neighbour')
                CALL create_map_from_mesh_to_mesh_nearest_neighbour(      mesh_src, mesh_dst, Atlas( mi))
              CASE('trilin')
                CALL create_map_from_mesh_to_mesh_trilin(                 mesh_src, mesh_dst, Atlas( mi))
              CASE('2nd_order_conservative')
                CALL create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
              CASE DEFAULT
                CALL crash('unknown remapping method "' // TRIM( method) // '"')
            END SELECT
          ELSE
              CALL create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
          END IF
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_mesh_to_mesh_2D( mesh_src, mesh_dst, Atlas( mi), d_src_partial, d_dst_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_mesh_to_mesh_2D

  SUBROUTINE map_from_mesh_to_mesh_3D( mesh_src, mesh_dst, d_src_partial, d_dst_partial, method)
    ! Map a 3-D data field from a mesh to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_src_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_dst_partial
    CHARACTER(LEN=*), OPTIONAL,          INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_mesh_to_mesh_3D'
    LOGICAL                                            :: are_identical
    INTEGER                                            :: mi, mi_valid
    LOGICAL                                            :: found_map, found_empty_page

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the two meshes are identical, the remapping operation is trivial
    CALL check_if_meshes_are_identical( mesh_src, mesh_dst, are_identical)
    IF (are_identical) THEN
      d_dst_partial = d_src_partial
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == mesh_src%name .AND. Atlas( mi)%name_dst == mesh_dst%name) THEN
        ! If so specified, look for a mapping object with the correct method
        IF (PRESENT( method)) THEN
          IF (Atlas( mi)%method /= method) CYCLE
        END IF
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          IF (PRESENT( method)) THEN
            SELECT CASE (method)
              CASE ('nearest_neighbour')
                CALL create_map_from_mesh_to_mesh_nearest_neighbour(      mesh_src, mesh_dst, Atlas( mi))
              CASE ('trilin')
                CALL create_map_from_mesh_to_mesh_trilin(                 mesh_src, mesh_dst, Atlas( mi))
              CASE ('2nd_order_conservative')
                CALL create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
              CASE DEFAULT
                CALL crash('unknown remapping method "' // TRIM( method) // '"')
            END SELECT
          ELSE
              CALL create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, Atlas( mi))
          END IF
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Apply the appropriate mapping object
    CALL apply_map_mesh_to_mesh_3D( mesh_src, mesh_dst, Atlas( mi), d_src_partial, d_dst_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_mesh_to_mesh_3D

  ! From a vertical grid to another within the same mesh
  SUBROUTINE map_from_vertical_to_vertical_2D_ocean( mesh, vert_src, vert_dst, d_src_partial, d_dst_partial)
    ! Map mesh data between vertical grids

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                                       INTENT(IN)  :: mesh
    REAL(dp), DIMENSION(:),                                INTENT(IN)  :: vert_src
    REAL(dp), DIMENSION(:),                                INTENT(IN)  :: vert_dst
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,SIZE(vert_src)), INTENT(IN)  :: d_src_partial
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,SIZE(vert_dst)), INTENT(OUT) :: d_dst_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER      :: routine_name = 'map_from_vertical_to_vertical_2D_ocean'
    INTEGER                            :: vi,k
    INTEGER, DIMENSION(:), ALLOCATABLE :: z_mask_old, z_mask_new, mask_fill
    REAL(dp)                           :: z_floor
    REAL(dp)                           :: NaN
    REAL(dp), PARAMETER                :: sigma = 4e4

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    NaN = -0.1234_dp

    ! Allocate mask for valid points in a data column
    ALLOCATE( z_mask_old( SIZE( vert_src)))
    ALLOCATE( z_mask_new( SIZE( vert_dst)))

    DO vi = mesh%vi1, mesh%vi2

      ! Determine local depth of the ocean floor, fill in both data masks
      IF (d_src_partial( vi, SIZE( vert_src)) == &
          d_src_partial( vi, SIZE( vert_src))) THEN

        ! Ocean floor lies below the vertical limit of the provided data
        z_mask_old = 1
        z_floor = vert_src( SIZE( vert_src)) + (vert_src( 2) - vert_src( 1))

      ELSEIF (d_src_partial( vi,1) /= d_src_partial( vi,1)) THEN

        ! This grid cell isn't ocean at all
        z_mask_old = 0
        z_floor    = 0._dp
        NaN = d_src_partial( vi,1)

      ELSE

        z_mask_old = 1
        k = SIZE( vert_src)

        DO WHILE (d_src_partial( vi,k) /= d_src_partial( vi,k))
          z_mask_old( k) = 0
          z_floor = vert_src( k)
          k = k - 1
          NaN = d_src_partial( vi,k)
        END DO

      END IF

      z_mask_new = 0

      DO k = 1, SIZE(vert_dst)
        IF (vert_dst( k) < z_floor) THEN
          z_mask_new = 1
        END IF
      END DO

      ! Regrid vertical column
      CALL remap_cons_2nd_order_1D( vert_src, z_mask_old, d_src_partial( vi,:), &
                                    vert_dst, z_mask_new, d_dst_partial( vi,:))

      ! Fill masked values with NaN
      DO k = 1, SIZE( vert_dst)
        IF (z_mask_new( k) == 0) THEN
          d_dst_partial( vi,k) = NaN
        END IF
      END DO

    END DO

    ! Allocate mask for extrapolation
    ALLOCATE( mask_fill( mesh%vi1:mesh%vi2))

    ! Extrapolate into NaN areas independently for each layer
    DO k = 1, SIZE(vert_dst)
      ! Initialise assuming there's valid data everywhere
      mask_fill = 2
      ! Check this mesh layer for NaNs
      DO vi = mesh%vi1, mesh%vi2
        IF (d_dst_partial( vi,k) /= d_dst_partial( vi,k)) THEN
          ! If NaN, allow extrapolation here
          mask_fill( vi) = 1
        END IF
      END DO
      ! Fill NaN vertices within this layer
      CALL extrapolate_Gaussian( mesh, mask_fill, d_dst_partial(:,k), sigma)
    END DO

    ! Clean up after yourself
    DEALLOCATE( z_mask_old)
    DEALLOCATE( z_mask_new)
    DEALLOCATE( mask_fill)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_from_vertical_to_vertical_2D_ocean

  ! Smoothing operations on the mesh
  SUBROUTINE smooth_Gaussian_2D( mesh, grid, d_mesh_partial, r)
    ! Use 2nd-order conservative remapping to map the 2-D data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_mesh_partial
    REAL(dp),                            INTENT(IN)    :: r

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_grid_vec_partial

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid%n_loc))

    ! Map data to the grid
    CALL map_from_mesh_to_xy_grid_2D( mesh, grid, d_mesh_partial, d_grid_vec_partial)

    ! Apply smoothing on the gridded data
    CALL smooth_Gaussian_2D_grid( grid, d_grid_vec_partial, r)

    ! Map data back to the mesh
    CALL map_from_xy_grid_to_mesh_2D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Clean up after yourself
    DEALLOCATE( d_grid_vec_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE smooth_Gaussian_2D

  SUBROUTINE smooth_Gaussian_3D( mesh, grid, d_mesh_partial, r)
    ! Use 2nd-order conservative remapping to map the 3-D data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_mesh_partial
    REAL(dp),                            INTENT(IN)    :: r

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_3D'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid_vec_partial

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( d_grid_vec_partial( grid%n_loc, SIZE( d_mesh_partial,2)))

    ! Map data to the grid
    CALL map_from_mesh_to_xy_grid_3D( mesh, grid, d_mesh_partial, d_grid_vec_partial)

    ! Apply smoothing on the gridded data
    CALL smooth_Gaussian_3D_grid( grid, d_grid_vec_partial, r)

    ! Map data back to the mesh
    CALL map_from_xy_grid_to_mesh_3D( grid, mesh, d_grid_vec_partial, d_mesh_partial)

    ! Clean up after yourself
    DEALLOCATE( d_grid_vec_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE smooth_Gaussian_3D

  ! Clean up the Atlas after a mesh update
  SUBROUTINE clear_all_maps_involving_this_mesh( mesh)
    ! Clear all mapping objects involving a mesh by this name from the Atlas

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_from_xy_grid_to_mesh_2D'
    INTEGER                                            :: mi, perr

    ! Add routine to path
    CALL init_routine( routine_name)

    DO mi = 1, SIZE( Atlas,1)
      IF (Atlas( mi)%is_in_use .AND. &
        (Atlas( mi)%name_src == mesh%name .OR. &
         Atlas( mi)%name_dst == mesh%name .OR. &
         Atlas( mi)%name_src == (TRIM( mesh%name) // '_triangles') .OR. &
         Atlas( mi)%name_dst == (TRIM( mesh%name) // '_triangles'))) THEN
        ! This map involves the current mesh
        Atlas( mi)%is_in_use = .FALSE.
        Atlas( mi)%name_src  = ''
        Atlas( mi)%name_dst  = ''
        CALL MatDestroy( Atlas( mi)%M, perr)
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE clear_all_maps_involving_this_mesh

! == Apply existing mapping objects to remap data between grids
! =============================================================

  ! From an x/y-grid to a mesh
  SUBROUTINE apply_map_xy_grid_to_mesh_2D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 2-D data field from an x/y-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_xy_grid_to_mesh_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nV_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_xy_grid_to_mesh_2D

  SUBROUTINE apply_map_xy_grid_to_mesh_3D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 3-D data field from an x/y-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_xy_grid_to_mesh_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nV_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc .OR. &
      SIZE( d_grid_vec_partial,2) /= SIZE( d_mesh_partial,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_xy_grid_to_mesh_3D

  ! From an x/y-grid to a mesh
  SUBROUTINE apply_map_xy_grid_to_mesh_triangles_2D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 2-D data field from an x/y-grid to a mesh triangles.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_xy_grid_to_mesh_triangles_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nTri_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_xy_grid_to_mesh_triangles_2D

  SUBROUTINE apply_map_xy_grid_to_mesh_triangles_3D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 3-D data field from an x/y-grid to a mesh triangles.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_xy_grid_to_mesh_triangles_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nTri_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc .OR. &
      SIZE( d_grid_vec_partial,2) /= SIZE( d_mesh_partial,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_xy_grid_to_mesh_triangles_3D

  ! From a lon/lat-grid to a mesh
  SUBROUTINE apply_map_lonlat_grid_to_mesh_2D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 2-D data field from a lon/lat-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_lonlat_grid_to_mesh_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nV_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_lonlat_grid_to_mesh_2D

  SUBROUTINE apply_map_lonlat_grid_to_mesh_3D( grid, mesh, map, d_grid_vec_partial, d_mesh_partial)
    ! Map a 3-D data field from a lon/lat-grid to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_lonlat_grid_to_mesh_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nV_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc .OR. &
      SIZE( d_grid_vec_partial,2) /= SIZE( d_mesh_partial,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_lonlat_grid_to_mesh_3D

  ! From a mesh to an x/y-grid
  SUBROUTINE apply_map_mesh_to_xy_grid_2D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    ! Map a 2-D data field from a mesh to an x/y-grid.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_mesh_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_mesh_to_xy_grid_2D'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nV_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( map%M, d_mesh_partial, d_grid_vec_partial)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    IF (par%master) THEN
      ! Allocate memory for complete gridded data
      ALLOCATE( d_grid( grid%nx, grid%ny))
      ! Gather complete gridded data
      CALL gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:) = d_grid( 2        ,:)
      d_grid( grid%nx,:) = d_grid( grid%nx-1,:)
      d_grid( :,1      ) = d_grid( :,2        )
      d_grid( :,grid%ny) = d_grid( :,grid%ny-1)
      ! Distribute complete gridded data back over the processes
      CALL distribute_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      DEALLOCATE( d_grid)
    ELSE ! IF (par%master) THEN
      ! Allocate zero memory for complete gridded data (only the master needs this)
      ALLOCATE( d_grid( 0,0))
      ! Gather complete gridded data
      CALL gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial)
      ! Distribute complete gridded data back over the processes
      CALL distribute_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      DEALLOCATE( d_grid)
    END IF ! IF (par%master) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_mesh_to_xy_grid_2D

  SUBROUTINE apply_map_mesh_to_xy_grid_3D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    ! Map a 3-D data field from a mesh to an x/y-grid.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_mesh_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_mesh_to_xy_grid_3D'
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nV_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc .OR. &
      SIZE( d_mesh_partial,2) /= SIZE( d_grid_vec_partial,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( map%M, d_mesh_partial, d_grid_vec_partial)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    IF (par%master) THEN
      ! Allocate memory for complete gridded data
      ALLOCATE( d_grid( grid%nx, grid%ny, SIZE( d_mesh_partial,2)))
      ! Gather complete gridded data
      CALL gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:,:) = d_grid( 2        ,:,:)
      d_grid( grid%nx,:,:) = d_grid( grid%nx-1,:,:)
      d_grid( :,1      ,:) = d_grid( :,2        ,:)
      d_grid( :,grid%ny,:) = d_grid( :,grid%ny-1,:)
      ! Distribute complete gridded data back over the processes
      CALL distribute_gridded_data_from_master_dp_3D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      DEALLOCATE( d_grid)
    ELSE ! IF (par%master) THEN
      ! Allocate zero memory for complete gridded data (only the master needs this)
      ALLOCATE( d_grid( 0,0,0))
      ! Gather complete gridded data
      CALL gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)
      ! Distribute complete gridded data back over the processes
      CALL distribute_gridded_data_from_master_dp_3D( grid, d_grid, d_grid_vec_partial)
      ! Clean up after yourself
      DEALLOCATE( d_grid)
    END IF ! IF (par%master) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_mesh_to_xy_grid_3D

  SUBROUTINE apply_map_mesh_to_xy_grid_2D_minval( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    ! Map a 2-D data field from a mesh to an x/y-grid.
    !
    ! For each grid cell, get the minimum value of all overlapping mesh vertices

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_grid),                        INTENT(IN)    :: grid
    TYPE(type_map),                         INTENT(IN)    :: map
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: d_mesh_partial
    REAL(dp), DIMENSION(grid%n1 :grid%n2 ), INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_map_mesh_to_xy_grid_2D_minval'
    REAL(dp), DIMENSION(mesh%nV)                          :: d_mesh_tot
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_CSR
    INTEGER                                               :: n,k1,k2,k,col,vi
    REAL(dp)                                              :: d_max, d_min

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh_partial,1) /= mesh%nV_loc .OR. SIZE( d_grid_vec_partial,1) /= grid%n_loc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Gather global mesh data
    CALL gather_to_all_dp_1D( d_mesh_partial, d_mesh_tot)

    ! Convert mapping matrix from PETSc format to UFEMISM CSR format
    CALL mat_petsc2CSR( map%M, M_CSR)

    ! Find global maximum value of d
    d_max = MAXVAL( d_mesh_partial)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, d_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Map data
    DO n = grid%n1, grid%n2

      ! Initialise minimum as maximum
      d_min = d_max

      ! Loop over all mesh vertices that this grid cell overlaps with
      k1 = M_CSR%ptr( n)
      k2 = M_CSR%ptr( n+1)-1
      DO k = k1, k2
        col = M_CSR%ind( k)
        ! This matrix row corresponds to this mesh vertex
        vi = mesh%n2vi( col)
        ! Update minimum value
        d_min = MIN( d_min, d_mesh_tot( vi))
      END DO

      ! Fill into array
      d_grid_vec_partial( n) = d_min

    END DO ! DO n = grid%n1, grid%n2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_mesh_to_xy_grid_2D_minval

  ! From a mesh to a mesh
  SUBROUTINE apply_map_mesh_to_mesh_2D( mesh_src, mesh_dst, map, d_src_partial, d_dst_partial)
    ! Map a 2-D data field from a mesh to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_src_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_dst_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_mesh_to_mesh_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_src_partial,1) /= mesh_src%nV_loc .OR. SIZE( d_dst_partial,1) /= mesh_dst%nV_loc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( map%M, d_src_partial, d_dst_partial)

    ! Set values of border vertices to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    CALL set_border_vertices_to_interior_mean_dp_2D( mesh_dst, d_dst_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_mesh_to_mesh_2D

  SUBROUTINE apply_map_mesh_to_mesh_3D( mesh_src, mesh_dst, map, d_src_partial, d_dst_partial)
    ! Map a 3-D data field from a mesh to a mesh.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(type_map),                      INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_src_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_dst_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_map_mesh_to_mesh_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_src_partial,1) /= mesh_src%nV_loc .OR. SIZE( d_dst_partial,1) /= mesh_dst%nV_loc .OR. &
      SIZE( d_src_partial,2) /= SIZE( d_dst_partial,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( map%M, d_src_partial, d_dst_partial)

    ! Set values of border vertices to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    CALL set_border_vertices_to_interior_mean_dp_3D( mesh_dst, d_dst_partial)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_map_mesh_to_mesh_3D

! == Create remapping objects
! ===========================

  SUBROUTINE create_map_from_xy_grid_to_mesh( grid, mesh, map)
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

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_map_from_xy_grid_to_mesh'
    TYPE(PetscErrorCode)                               :: perr
    LOGICAL                                            :: count_coincidences
    INTEGER                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: mask_do_simple_average
    INTEGER                                            :: vi
    REAL(dp), DIMENSION( mesh%nC_mem,2)                :: Vor
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_vi
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_ti
    INTEGER                                            :: nVor
    INTEGER                                            :: vori1, vori2
    REAL(dp), DIMENSION(2)                             :: p, q
    INTEGER                                            :: k, i, j, kk, vj
    REAL(dp)                                           :: xl, xu, yl, yu
    REAL(dp), DIMENSION(2)                             :: sw, se, nw, ne
    INTEGER                                            :: vi_hint
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    INTEGER                                            :: il, iu, jl, ju
    TYPE(type_single_row_mapping_matrices)             :: single_row_Vor, single_row_grid
    TYPE(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    TYPE(tMat)                                         :: w0    , w1x    , w1y
    INTEGER                                            :: row, k1, k2, col
    REAL(dp)                                           :: A_overlap_tot
    TYPE(tMat)                                         :: grid_M_ddx, grid_M_ddy
    TYPE(tMat)                                         :: M1, M2

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (map%is_in_use) CALL crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

    map%is_in_use = .TRUE.
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
    nnz_est         = 4 * MAX( nrows, ncols)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh%A) / (grid%dx**2)), &
                                    CEILING( 2._dp * (grid%dx**2) / MINVAL( mesh%A))) )

    CALL allocate_matrix_CSR_dist( A_xdy_a_g_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_mxydx_a_g_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_xydy_a_g_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for single row results
    single_row_Vor%n_max = nnz_per_row_max
    single_row_Vor%n     = 0
    ALLOCATE( single_row_Vor%index_left( single_row_Vor%n_max))
    ALLOCATE( single_row_Vor%LI_xdy(     single_row_Vor%n_max))
    ALLOCATE( single_row_Vor%LI_mxydx(   single_row_Vor%n_max))
    ALLOCATE( single_row_Vor%LI_xydy(    single_row_Vor%n_max))

    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    ALLOCATE( single_row_grid%index_left( single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xdy(     single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xydy(    single_row_grid%n_max))

    ALLOCATE( mask_do_simple_average( mesh%nV), source = 0)

    ! Calculate line integrals around all Voronoi cells
    DO row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      IF (mesh%A( vi) < 10._dp * grid%dx**2) THEN
        ! This Voronoi cell is small enough to warrant a proper line integral

        mask_do_simple_average( vi) = 0

        ! Clean up single row results
        single_row_Vor%n          = 0
        single_row_Vor%index_left = 0
        single_row_Vor%LI_xdy     = 0._dp
        single_row_Vor%LI_mxydx   = 0._dp
        single_row_Vor%LI_xydy    = 0._dp

        ! Integrate around the complete Voronoi cell boundary
        CALL calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)
        DO vori1 = 1, nVor
          vori2 = vori1 + 1
          IF (vori2 > nVor) vori2 = 1
          p = Vor( vori1,:)
          q = Vor( vori2,:)
          count_coincidences = .TRUE.
          CALL trace_line_grid( grid, p, q, single_row_Vor, count_coincidences)
        END DO

        ! Safety
        IF (single_row_Vor%n == 0) CALL crash('couldnt find any grid cells overlapping with this small Voronoi cell!')

        ! Next integrate around the grid cells overlapping with this Voronoi cell
        DO k = 1, single_row_Vor%n

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
          count_coincidences = .FALSE.
          CALL trace_line_Vor( mesh, sw, se, single_row_grid, count_coincidences, vi_hint)
          CALL trace_line_Vor( mesh, se, ne, single_row_grid, count_coincidences, vi_hint)
          CALL trace_line_Vor( mesh, ne, nw, single_row_grid, count_coincidences, vi_hint)
          CALL trace_line_Vor( mesh, nw, sw, single_row_grid, count_coincidences, vi_hint)

          ! Safety
          IF (single_row_grid%n == 0) CALL crash('couldnt find any Voronoi cells overlapping with this grid cell!')

          ! Add contribution for this particular triangle
          DO kk = 1, single_row_grid%n
            vj = single_row_grid%index_left( kk)
            IF (vj == vi) THEN
              ! Add contribution to this triangle
              single_row_Vor%LI_xdy(   k) = single_row_Vor%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
              single_row_Vor%LI_mxydx( k) = single_row_Vor%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
              single_row_Vor%LI_xydy(  k) = single_row_Vor%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
              EXIT
            END IF
          END DO ! DO kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          CALL add_entry_CSR_dist( A_xdy_a_g_CSR  , row, col, single_row_Vor%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_a_g_CSR, row, col, single_row_Vor%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_a_g_CSR , row, col, single_row_Vor%LI_xydy(  k))

        END DO ! DO k = 1, single_row_Vor%n

      ELSE ! IF (mesh%A( vi) < 10._dp * grid%dx**2) THEN
        ! This Voronoi cell is big enough that we can just average over the grid cells it contains

        mask_do_simple_average( vi) = 1

        ! Clean up single row results
        single_row_Vor%n = 0

        ! Find the square of grid cells enveloping this Voronoi cell
        CALL calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)

        xmin = MINVAL( Vor( 1:nVor,1))
        xmax = MAXVAL( Vor( 1:nVor,1))
        ymin = MINVAL( Vor( 1:nVor,2))
        ymax = MAXVAL( Vor( 1:nVor,2))

        il = MAX( 1, MIN( grid%nx, 1 + FLOOR( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        iu = MAX( 1, MIN( grid%nx, 1 + FLOOR( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        jl = MAX( 1, MIN( grid%ny, 1 + FLOOR( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx) ))
        ju = MAX( 1, MIN( grid%ny, 1 + FLOOR( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

        ! Check which of the grid cells in this square lie inside the triangle
        DO i = il, iu
        DO j = jl, ju

          col = grid%ij2n( i,j)
          p   = [grid%x( i), grid%y( j)]

          IF (is_in_Voronoi_cell( mesh, p, vi)) THEN
            ! This grid cell lies inside the triangle; add it to the single row
            single_row_Vor%n = single_row_Vor%n + 1
            single_row_Vor%index_left( single_row_Vor%n) = col
            single_row_Vor%LI_xdy(     single_row_Vor%n) = grid%dx**2
            single_row_Vor%LI_mxydx(   single_row_Vor%n) = grid%x( i) * grid%dx**2
            single_row_Vor%LI_xydy(    single_row_Vor%n) = grid%y( j) * grid%dx**2
          END IF

        END DO
        END DO

        ! Safety
        IF (single_row_Vor%n == 0) CALL crash('couldnt find any grid cells overlapping with this big Voronoi cell!')

        ! Add entries to the big matrices
        DO k = 1, single_row_Vor%n
          col = single_row_Vor%index_left( k)
          CALL add_entry_CSR_dist( A_xdy_a_g_CSR  , vi, col, single_row_Vor%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_a_g_CSR, vi, col, single_row_Vor%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_a_g_CSR , vi, col, single_row_Vor%LI_xydy(  k))
        END DO

      END IF ! IF (mesh%A( vi) < 4._dp * grid%dx**2) THEN

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( single_row_Vor%index_left )
    DEALLOCATE( single_row_Vor%LI_xdy     )
    DEALLOCATE( single_row_Vor%LI_mxydx   )
    DEALLOCATE( single_row_Vor%LI_xydy    )

    DEALLOCATE( single_row_grid%index_left )
    DEALLOCATE( single_row_grid%LI_xdy     )
    DEALLOCATE( single_row_grid%LI_mxydx   )
    DEALLOCATE( single_row_grid%LI_xydy    )

  ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
  ! ==============================================================

    CALL allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    DO row = mesh%vi1, mesh%vi2

      vi = mesh%n2vi( row)

      k1 = A_xdy_a_g_CSR%ptr( row  )
      k2 = A_xdy_a_g_CSR%ptr( row+1) - 1

      A_overlap_tot = SUM( A_xdy_a_g_CSR%val( k1:k2))

      DO k = k1, k2
        col = A_xdy_a_g_CSR%ind( k)
        CALL add_entry_CSR_dist( w0_CSR, row, col, A_xdy_a_g_CSR%val( k) / A_overlap_tot)
      END DO

      IF (mask_do_simple_average( vi) == 0) THEN
        ! For small vertices, include the gradient terms

        DO k = k1, k2
          col = A_xdy_a_g_CSR%ind( k)
          ! Grid cell
          i = grid%n2ij( col,1)
          j = grid%n2ij( col,2)
          CALL add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_a_g_CSR%val( k) / A_overlap_tot) - (grid%x( i) * w0_CSR%val( k)))
          CALL add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_a_g_CSR%val(  k) / A_overlap_tot) - (grid%y( j) * w0_CSR%val( k)))
        END DO

      ELSE
        ! For large vertices, don't include the gradient terms

        CALL add_empty_row_CSR_dist( w1x_CSR, row)
        CALL add_empty_row_CSR_dist( w1y_CSR, row)

      END IF ! IF (mask_do_simple_average( vi) == 0) THEN

    END DO ! DO row = mesh%vi1, mesh%vi2
    CALL sync

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( A_xdy_a_g_CSR  )
    CALL deallocate_matrix_CSR_dist( A_mxydx_a_g_CSR)
    CALL deallocate_matrix_CSR_dist( A_xydy_a_g_CSR )

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( w0_CSR , w0 )
    CALL mat_CSR2petsc( w1x_CSR, w1x)
    CALL mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix

    CALL calc_matrix_operators_grid( grid, grid_M_ddx, grid_M_ddy)

    CALL MatDuplicate( w0, MAT_COPY_VALUES, map%M, perr)
    CALL MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    CALL MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2, perr)

    CALL MatDestroy( grid_M_ddx    , perr)
    CALL MatDestroy( grid_M_ddy    , perr)
    CALL MatDestroy( w0            , perr)
    CALL MatDestroy( w1x           , perr)
    CALL MatDestroy( w1y           , perr)

    CALL MatAXPY( map%M, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( map%M, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    CALL MatDestroy( M2, perr)
    DEALLOCATE( mask_do_simple_average)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_map_from_xy_grid_to_mesh

  SUBROUTINE create_map_from_xy_grid_to_mesh_triangles( grid, mesh, map)
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

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_map_from_xy_grid_to_mesh_triangles'
    TYPE(PetscErrorCode)                               :: perr
    LOGICAL                                            :: count_coincidences
    INTEGER                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_xdy_b_g_CSR, A_mxydx_b_g_CSR, A_xydy_b_g_CSR
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: mask_do_simple_average
    INTEGER                                            :: ti
    INTEGER                                            :: n1, n2, via, vib, vic
    REAL(dp), DIMENSION(2)                             :: p, q
    INTEGER                                            :: k, i, j, kk, tj
    REAL(dp)                                           :: xl, xu, yl, yu
    REAL(dp), DIMENSION(2)                             :: sw, se, nw, ne
    INTEGER                                            :: ti_hint
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    INTEGER                                            :: il, iu, jl, ju
    TYPE(type_single_row_mapping_matrices)             :: single_row_tri, single_row_grid
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    TYPE(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    TYPE(tMat)                                         :: w0    , w1x    , w1y
    INTEGER                                            :: row, k1, k2, col
    REAL(dp)                                           :: A_overlap_tot
    TYPE(tMat)                                         :: grid_M_ddx, grid_M_ddy
    TYPE(tMat)                                         :: M1, M2

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (map%is_in_use) CALL crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

    map%is_in_use = .TRUE.
    map%name_src  = grid%name
    map%name_dst  = TRIM( mesh%name) // '_triangles'
    map%method    = '2nd_order_conservative'

  ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
  ! ===========================================================================

    ! Matrix size
    nrows           = mesh%nTri  ! to
    nrows_loc       = mesh%nTri_loc
    ncols           = grid%n     ! from
    ncols_loc       = grid%n_loc
    nnz_est         = 4 * MAX( nrows, ncols)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh%TriA) / (grid%dx**2)), &
                                    CEILING( 2._dp * (grid%dx**2) / MINVAL( mesh%TriA))) )

    CALL allocate_matrix_CSR_dist( A_xdy_b_g_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_mxydx_b_g_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_xydy_b_g_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for single row results
    single_row_tri%n_max = nnz_per_row_max
    single_row_tri%n     = 0
    ALLOCATE( single_row_tri%index_left( single_row_tri%n_max))
    ALLOCATE( single_row_tri%LI_xdy(     single_row_tri%n_max))
    ALLOCATE( single_row_tri%LI_mxydx(   single_row_tri%n_max))
    ALLOCATE( single_row_tri%LI_xydy(    single_row_tri%n_max))

    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    ALLOCATE( single_row_grid%index_left( single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xdy(     single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xydy(    single_row_grid%n_max))

    ALLOCATE( mask_do_simple_average( mesh%nTri), source = 0)

    ! Calculate line integrals around all triangles
    DO row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      IF (mesh%TriA( ti) < 10._dp * grid%dx**2) THEN
        ! This triangle is small enough to warrant a proper line integral

        mask_do_simple_average( ti) = 0

        ! Clean up single row results
        single_row_tri%n          = 0
        single_row_tri%index_left = 0
        single_row_tri%LI_xdy     = 0._dp
        single_row_tri%LI_mxydx   = 0._dp
        single_row_tri%LI_xydy    = 0._dp

        ! Integrate around the triangle
        DO n1 = 1, 3
          n2 = n1 + 1
          IF (n2 == 4) n2 = 1
          via = mesh%Tri( ti,n1)
          vib = mesh%Tri( ti,n2)
          p = mesh%V( via,:)
          q = mesh%V( vib,:)
          count_coincidences = .TRUE.
          CALL trace_line_grid( grid, p, q, single_row_tri, count_coincidences)
        END DO

        ! Safety
        IF (single_row_tri%n == 0) CALL crash('couldnt find any grid cells overlapping with this small triangle!')

        ! Next integrate around the grid cells overlapping with this triangle
        DO k = 1, single_row_tri%n

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
          count_coincidences = .FALSE.
          CALL trace_line_tri( mesh, sw, se, single_row_grid, count_coincidences, ti_hint)
          CALL trace_line_tri( mesh, se, ne, single_row_grid, count_coincidences, ti_hint)
          CALL trace_line_tri( mesh, ne, nw, single_row_grid, count_coincidences, ti_hint)
          CALL trace_line_tri( mesh, nw, sw, single_row_grid, count_coincidences, ti_hint)

          ! Safety
          IF (single_row_grid%n == 0) CALL crash('couldnt find any triangles overlapping with this grid cell!')

          ! Add contribution for this particular triangle
          DO kk = 1, single_row_grid%n
            tj = single_row_grid%index_left( kk)
            IF (tj == ti) THEN
              ! Add contribution to this triangle
              single_row_tri%LI_xdy(   k) = single_row_tri%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
              single_row_tri%LI_mxydx( k) = single_row_tri%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
              single_row_tri%LI_xydy(  k) = single_row_tri%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
              EXIT
            END IF
          END DO ! DO kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          CALL add_entry_CSR_dist( A_xdy_b_g_CSR  , row, col, single_row_tri%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_b_g_CSR, row, col, single_row_tri%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_b_g_CSR , row, col, single_row_tri%LI_xydy(  k))

        END DO ! DO k = 1, single_row_tri%n

      ELSE ! IF (mesh%TriA( ti) < 10._dp * grid%dx**2) THEN
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

        xmin = MIN( MIN( pa( 1), pb( 1)), pc( 1))
        xmax = MAX( MAX( pa( 1), pb( 1)), pc( 1))
        ymin = MIN( MIN( pa( 2), pb( 2)), pc( 2))
        ymax = MAX( MAX( pa( 2), pb( 2)), pc( 2))

        il = MAX( 1, MIN( grid%nx, 1 + FLOOR( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        iu = MAX( 1, MIN( grid%nx, 1 + FLOOR( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        jl = MAX( 1, MIN( grid%ny, 1 + FLOOR( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx) ))
        ju = MAX( 1, MIN( grid%ny, 1 + FLOOR( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

        ! Check which of the grid cells in this square lie inside the triangle
        DO i = il, iu
        DO j = jl, ju

          col = grid%ij2n( i,j)
          p   = [grid%x( i), grid%y( j)]

          IF (is_in_triangle( pa, pb, pc, p)) THEN
            ! This grid cell lies inside the triangle; add it to the single row
            single_row_tri%n = single_row_tri%n + 1
            single_row_tri%index_left( single_row_tri%n) = col
            single_row_tri%LI_xdy(     single_row_tri%n) = grid%dx**2
            single_row_tri%LI_mxydx(   single_row_tri%n) = grid%x( i) * grid%dx**2
            single_row_tri%LI_xydy(    single_row_tri%n) = grid%y( j) * grid%dx**2
          END IF

        END DO
        END DO

        ! Safety
        IF (single_row_tri%n == 0) CALL crash('couldnt find any grid cells overlapping with this big triangle!')

        ! Add entries to the big matrices
        DO k = 1, single_row_tri%n
          col = single_row_tri%index_left( k)
          CALL add_entry_CSR_dist( A_xdy_b_g_CSR  , ti, col, single_row_tri%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_b_g_CSR, ti, col, single_row_tri%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_b_g_CSR , ti, col, single_row_tri%LI_xydy(  k))
        END DO

      END IF ! IF (mesh%TriA( ti) < 10._dp * grid%dx**2) THEN

    END DO ! DO row = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    DEALLOCATE( single_row_tri%index_left )
    DEALLOCATE( single_row_tri%LI_xdy     )
    DEALLOCATE( single_row_tri%LI_mxydx   )
    DEALLOCATE( single_row_tri%LI_xydy    )

    DEALLOCATE( single_row_grid%index_left )
    DEALLOCATE( single_row_grid%LI_xdy     )
    DEALLOCATE( single_row_grid%LI_mxydx   )
    DEALLOCATE( single_row_grid%LI_xydy    )

  ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
  ! ==============================================================

    CALL allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    DO row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      k1 = A_xdy_b_g_CSR%ptr( row  )
      k2 = A_xdy_b_g_CSR%ptr( row+1) - 1

      A_overlap_tot = SUM( A_xdy_b_g_CSR%val( k1:k2))

      DO k = k1, k2
        col = A_xdy_b_g_CSR%ind( k)
        CALL add_entry_CSR_dist( w0_CSR, row, col, A_xdy_b_g_CSR%val( k) / A_overlap_tot)
      END DO

      IF (mask_do_simple_average( ti) == 0) THEN
        ! For small triangles, include the gradient terms

        DO k = k1, k2
          col = A_xdy_b_g_CSR%ind( k)
          ! Grid cell
          i = grid%n2ij( col,1)
          j = grid%n2ij( col,2)
          CALL add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_b_g_CSR%val( k) / A_overlap_tot) - (grid%x( i) * w0_CSR%val( k)))
          CALL add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_b_g_CSR%val(  k) / A_overlap_tot) - (grid%y( j) * w0_CSR%val( k)))
        END DO

      ELSE
        ! For large triangles, don't include the gradient terms

        CALL add_empty_row_CSR_dist( w1x_CSR, row)
        CALL add_empty_row_CSR_dist( w1y_CSR, row)

      END IF ! IF (mask_do_simple_average( vi) == 0) THEN

    END DO ! DO row = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( A_xdy_b_g_CSR  )
    CALL deallocate_matrix_CSR_dist( A_mxydx_b_g_CSR)
    CALL deallocate_matrix_CSR_dist( A_xydy_b_g_CSR )

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( w0_CSR , w0 )
    CALL mat_CSR2petsc( w1x_CSR, w1x)
    CALL mat_CSR2petsc( w1y_CSR, w1y)

    ! Calculate the remapping matrix

    CALL calc_matrix_operators_grid( grid, grid_M_ddx, grid_M_ddy)

    CALL MatDuplicate( w0, MAT_COPY_VALUES, map%M, perr)
    CALL MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    CALL MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2, perr)

    CALL MatDestroy( grid_M_ddx    , perr)
    CALL MatDestroy( grid_M_ddy    , perr)
    CALL MatDestroy( w0            , perr)
    CALL MatDestroy( w1x           , perr)
    CALL MatDestroy( w1y           , perr)

    CALL MatAXPY( map%M, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( map%M, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    CALL MatDestroy( M2, perr)
    DEALLOCATE( mask_do_simple_average)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_map_from_xy_grid_to_mesh_triangles

  SUBROUTINE create_map_from_mesh_to_xy_grid( mesh, grid, map)
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

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_map),                      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_map_from_mesh_to_xy_grid'
    TYPE(PetscErrorCode)                               :: perr
    LOGICAL                                            :: count_coincidences
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: overlaps_with_small_triangle, containing_triangle
    INTEGER                                            :: row, ti
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    INTEGER                                            :: il, iu, jl, ju
    INTEGER                                            :: i, j, n_ext, ii, jj
    INTEGER                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR
    TYPE(type_single_row_mapping_matrices)             :: single_row_grid, single_row_Tri
    INTEGER                                            :: ti_hint
    REAL(dp), DIMENSION(2)                             :: p
    REAL(dp)                                           :: xl, xu, yl, yu
    REAL(dp), DIMENSION(2)                             :: sw, se, nw, ne
    INTEGER                                            :: k, kk, nn
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy
    TYPE(type_sparse_matrix_CSR_dp)                    :: w0_CSR, w1x_CSR, w1y_CSR
    TYPE(tMat)                                         :: w0    , w1x    , w1y
    INTEGER                                            :: k1, k2, col
    REAL(dp)                                           :: A_overlap_tot
    TYPE(tMat)                                         :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    TYPE(tMat)                                         :: M1, M2
    INTEGER                                            :: ncols_row
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: cols_row
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: vals_row
    LOGICAL                                            :: has_value

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (map%is_in_use) CALL crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

    map%is_in_use = .TRUE.
    map%name_src  = mesh%name
    map%name_dst  = grid%name
    map%method    = '2nd_order_conservative'

  ! == Find all grid cells that overlap with small triangles
  ! ========================================================

    ALLOCATE( overlaps_with_small_triangle( grid%nx, grid%ny), source = 0)
    ALLOCATE( containing_triangle(          grid%nx, grid%ny), source = 0)

    DO row = mesh%ti1, mesh%ti2

      ti = mesh%n2ti( row)

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)

      ! The square enveloping this triangle
      xmin = MIN( MIN( pa(1), pb(1)), pc(1))
      xmax = MAX( MAX( pa(1), pb(1)), pc(1))
      ymin = MIN( MIN( pa(2), pb(2)), pc(2))
      ymax = MAX( MAX( pa(2), pb(2)), pc(2))

      ! The square of grid cells enveloping this triangle
      il = 1 + FLOOR( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx)
      iu = 1 + FLOOR( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx)
      jl = 1 + FLOOR( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx)
      ju = 1 + FLOOR( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx)

      il = MAX( 1      , il - 1)
      iu = MIN( grid%nx, iu + 1)
      jl = MAX( 1      , jl - 1)
      ju = MIN( grid%ny, ju + 1)

      IF (mesh%TriA( ti) < 10._dp * grid%dx**2) THEN
        ! This triangle is small; mark all grid cells it overlaps with

        ! Mark all these grid cells
        DO i = il, iu
        DO j = jl, ju
          overlaps_with_small_triangle( i,j) = 1
        END DO
        END DO

      ELSE
        ! This triangle is large; mark all grid cells it contains

        ! Mark all these grid cells
        DO i = il, iu
        DO j = jl, ju
          p = [grid%x( i), grid%y( j)]
          IF (is_in_triangle( pa, pb, pc, p)) THEN
            containing_triangle( i,j) = ti
          END IF
        END DO
        END DO

      END IF ! IF (mesh%TriA( ti) < 4._dp * grid%dx**2) THEN

    END DO

     ! Reduce results across the processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, overlaps_with_small_triangle, grid%nx * grid%ny, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, containing_triangle         , grid%nx * grid%ny, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Treat grid cells that possibly were not yet marked before
    DO row = grid%n1, grid%n2

      ! Grid cell indices
      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)

      IF (containing_triangle( i,j) == 0 .AND. overlaps_with_small_triangle( i,j) == 0) THEN
        ! This grid cell does not overlap with a small triangle, but was not yet marked
        ! as being contained inside a large one; find the large triangle containing it.

        ! For efficiency, find the nearest grid cell that does list which large
        ! triangle contains it; use that as a hint for the triangle search
        n_ext = 0
        ti_hint = 0
        DO WHILE (ti_hint == 0)
          n_ext = n_ext+1
          ! Safety
          IF (n_ext > MAX( grid%nx, grid%ny)) EXIT
          il = MAX( 1      , i-n_ext)
          iu = MIN( grid%nx, i+n_ext)
          jl = MAX( 1      , j-n_ext)
          ju = MIN( grid%ny, j+n_ext)
          DO ii = il, iu
          DO jj = jl, ju
            IF (containing_triangle( ii,jj) > 0) THEN
              ti_hint = containing_triangle( ii,jj)
              EXIT
            END IF
          END DO
          IF (ti_hint > 0) EXIT
          END DO
        END DO
        IF (ti_hint == 0) ti_hint = 1

        ! Find the triangle containing this grid cell
        p = [MAX( mesh%xmin, MIN( mesh%xmax, grid%x( i) )), MAX( mesh%ymin, MIN( mesh%ymax, grid%y( j) ))]
        CALL find_containing_triangle( mesh, p, ti_hint)
        containing_triangle( i,j) = ti_hint

      END IF

    END DO
    CALL sync

  ! == Integrate around all grid cells that overlap with small triangles
  ! ====================================================================

    ! Initialise the three matrices using the native UFEMISM CSR-matrix format

    ! Matrix size
    nrows           = grid%n     ! to
    nrows_loc       = grid%n_loc
    ncols           = mesh%nTri  ! from
    ncols_loc       = mesh%nTri_loc
    nnz_est         = 4 * MAX( nrows, ncols)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh%TriA) / (grid%dx**2)), &
                                    CEILING( 2._dp * (grid%dx**2) / MINVAL( mesh%TriA))) )

    CALL allocate_matrix_CSR_dist( A_xdy_g_b_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_mxydx_g_b_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_xydy_g_b_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for single row results
    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    ALLOCATE( single_row_grid%index_left( single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xdy(     single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xydy(    single_row_grid%n_max))

    single_row_Tri%n_max = nnz_per_row_max
    single_row_Tri%n     = 0
    ALLOCATE( single_row_Tri%index_left( single_row_Tri%n_max))
    ALLOCATE( single_row_Tri%LI_xdy(     single_row_Tri%n_max))
    ALLOCATE( single_row_Tri%LI_mxydx(   single_row_Tri%n_max))
    ALLOCATE( single_row_Tri%LI_xydy(    single_row_Tri%n_max))

    ti_hint = 1

    DO row = grid%n1, grid%n2

      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)
      p = [grid%x( i), grid%y( j)]

      ! The four sides of the grid cell
      xl = grid%x( i) - grid%dx / 2._dp
      xu = grid%x( i) + grid%dx / 2._dp
      yl = grid%y( j) - grid%dx / 2._dp
      yu = grid%y( j) + grid%dx / 2._dp

      ! If this grid cell lies entirely outside of the mesh domain, use
      ! nearest-neighbour extrapolation

      DO WHILE (xl <= mesh%xmin)
        i = i+1
        p( 1) = grid%x( i)
        xl = grid%x( i) - grid%dx / 2._dp
        xu = grid%x( i) + grid%dx / 2._dp
        IF (i > grid%nx) CALL crash('grid domain doesnt overlap with mesh domain at all!')
      END DO
      DO WHILE (xu >= mesh%xmax)
        i = i-1
        p( 1) = grid%x( i)
        xl = grid%x( i) - grid%dx / 2._dp
        xu = grid%x( i) + grid%dx / 2._dp
        IF (i < 1) CALL crash('grid domain doesnt overlap with mesh domain at all!')
      END DO
      DO WHILE (yl <= mesh%ymin)
        j = j+1
        p( 2) = grid%y( j)
        yl = grid%y( j) - grid%dx / 2._dp
        yu = grid%y( j) + grid%dx / 2._dp
        IF (j > grid%ny) CALL crash('grid domain doesnt overlap with mesh domain at all!')
      END DO
      DO WHILE (yu >= mesh%ymax)
        j = j-1
        p( 2) = grid%y( j)
        yl = grid%y( j) - grid%dx / 2._dp
        yu = grid%y( j) + grid%dx / 2._dp
        IF (j < 1) CALL crash('grid domain doesnt overlap with mesh domain at all!')
      END DO

      IF (overlaps_with_small_triangle( i,j) == 1) THEN
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
        count_coincidences = .TRUE.
        CALL trace_line_tri( mesh, sw, se, single_row_grid, count_coincidences, ti_hint)
        CALL trace_line_tri( mesh, se, ne, single_row_grid, count_coincidences, ti_hint)
        CALL trace_line_tri( mesh, ne, nw, single_row_grid, count_coincidences, ti_hint)
        CALL trace_line_tri( mesh, nw, sw, single_row_grid, count_coincidences, ti_hint)

        ! Next, integrate around all the triangles overlapping with this grid cell
        DO k = 1, single_row_grid%n

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
          count_coincidences = .FALSE.
          CALL trace_line_grid( grid, pa, pb, single_row_Tri, count_coincidences)
          CALL trace_line_grid( grid, pb, pc, single_row_Tri, count_coincidences)
          CALL trace_line_grid( grid, pc, pa, single_row_Tri, count_coincidences)

          ! Add contribution for this particular grid cell
          DO kk = 1, single_row_Tri%n
            nn = single_row_Tri%index_left( kk)
            IF (nn == row) THEN
              ! Add contribution to this triangle
              single_row_grid%LI_xdy(   k) = single_row_grid%LI_xdy(   k) + single_row_Tri%LI_xdy(   kk)
              single_row_grid%LI_mxydx( k) = single_row_grid%LI_mxydx( k) + single_row_Tri%LI_mxydx( kk)
              single_row_grid%LI_xydy(  k) = single_row_grid%LI_xydy(  k) + single_row_Tri%LI_xydy(  kk)
              EXIT
            END IF
          END DO ! DO kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          CALL add_entry_CSR_dist( A_xdy_g_b_CSR  , row, col, single_row_grid%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_g_b_CSR, row, col, single_row_grid%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_g_b_CSR , row, col, single_row_grid%LI_xydy(  k))

        END DO ! DO k = 1, single_row_grid%n

      ELSE ! IF (overlaps_with_small_triangle( i,j) == 1) THEN
        ! This grid cell does not overlap with a small triangle; use only the
        ! contribution from the nearest triangle

        ti_hint = containing_triangle( i,j)

        col = mesh%ti2n( ti_hint)

        LI_xdy   = grid%dx**2
        LI_mxydx = grid%dx**2 * grid%x( i)
        LI_xydy  = grid%dx**2 * grid%y( j)

        CALL add_entry_CSR_dist( A_xdy_g_b_CSR  , row, col, LI_xdy  )
        CALL add_entry_CSR_dist( A_mxydx_g_b_CSR, row, col, LI_mxydx)
        CALL add_entry_CSR_dist( A_xydy_g_b_CSR , row, col, LI_xydy )

      END IF ! IF (overlaps_with_small_triangle( i,j) == 1) THEN

    END DO ! DO n = grid%n1, grid%n2

    ! Clean up after yourself
    DEALLOCATE( overlaps_with_small_triangle)
    DEALLOCATE( containing_triangle         )

    DEALLOCATE( single_row_grid%index_left )
    DEALLOCATE( single_row_grid%LI_xdy     )
    DEALLOCATE( single_row_grid%LI_mxydx   )
    DEALLOCATE( single_row_grid%LI_xydy    )

    DEALLOCATE( single_row_Tri%index_left )
    DEALLOCATE( single_row_Tri%LI_xdy     )
    DEALLOCATE( single_row_Tri%LI_mxydx   )
    DEALLOCATE( single_row_Tri%LI_xydy    )

  ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
  ! ==============================================================

    CALL allocate_matrix_CSR_dist( w0_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( w1x_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( w1y_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    DO row = grid%n1, grid%n2

      k1 = A_xdy_g_b_CSR%ptr( row  )
      k2 = A_xdy_g_b_CSR%ptr( row+1) - 1

      A_overlap_tot = SUM( A_xdy_g_b_CSR%val( k1:k2))

      DO k = k1, k2
        col = A_xdy_g_b_CSR%ind( k)
        ti = mesh%n2ti( col)
        CALL add_entry_CSR_dist( w0_CSR , row, col,  A_xdy_g_b_CSR%val(   k) / A_overlap_tot)
        CALL add_entry_CSR_dist( w1x_CSR, row, col, (A_mxydx_g_b_CSR%val( k) / A_overlap_tot) - (mesh%TriGC( ti,1) * w0_CSR%val( k)))
        CALL add_entry_CSR_dist( w1y_CSR, row, col, (A_xydy_g_b_CSR%val(  k) / A_overlap_tot) - (mesh%TriGC( ti,2) * w0_CSR%val( k)))
      END DO

    END DO ! DO row = mesh%vi1, mesh%vi2
    CALL sync

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( A_xdy_g_b_CSR  )
    CALL deallocate_matrix_CSR_dist( A_mxydx_g_b_CSR)
    CALL deallocate_matrix_CSR_dist( A_xydy_g_b_CSR )

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( w0_CSR , w0 )
    CALL mat_CSR2petsc( w1x_CSR, w1x)
    CALL mat_CSR2petsc( w1y_CSR, w1y)

  ! Calculate the remapping matrix
  ! ==============================

    ! Convert matrices to PETSc format
    CALL mat_CSR2petsc( mesh%M_map_a_b, M_map_a_b)
    CALL mat_CSR2petsc( mesh%M_ddx_a_b, M_ddx_a_b)
    CALL mat_CSR2petsc( mesh%M_ddy_a_b, M_ddy_a_b)

    ! M = (w0 * M_map_a_b) + (w1x * M_ddx_a_b) + (w1y * M_ddy_a_b)
    CALL MatMatMult( w0,  M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, map%M, perr)
    CALL MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1   , perr)  ! This can be done more efficiently now that the non-zero structure is known...
    CALL MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2   , perr)
    CALL MatAXPY( map%M, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( map%M, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( w0            , perr)
    CALL MatDestroy( w1x           , perr)
    CALL MatDestroy( w1y           , perr)
    CALL MatDestroy( M_map_a_b     , perr)
    CALL MatDestroy( M_ddx_a_b     , perr)
    CALL MatDestroy( M_ddy_a_b     , perr)
    CALL MatDestroy( M1            , perr)
    CALL MatDestroy( M2            , perr)

  ! Safety: check if all grid cells get values
  ! ==========================================

    ALLOCATE( cols_row( nnz_per_row_max))
    ALLOCATE( vals_row( nnz_per_row_max))

    DO row = grid%n1, grid%n2

      ! w0
      CALL MatGetRow( map%M, row-1, ncols_row, cols_row, vals_row, perr)
      IF (ncols_row == 0) CALL crash('ncols == 0!')
      has_value = .FALSE.
      DO k = 1, ncols_row
        IF (vals_row( k) /= 0._dp) has_value = .TRUE.
      END DO
      IF (.NOT. has_value) CALL crash('only zeroes!')
      CALL MatRestoreRow( map%M, row-1, ncols_row, cols_row, vals_row, perr)

    END DO
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( cols_row)
    DEALLOCATE( vals_row)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_map_from_mesh_to_xy_grid

  SUBROUTINE create_map_from_lonlat_grid_to_mesh( grid, mesh, map)
    ! Create a new mapping object from a lon/lat-grid to a mesh.
    !
    ! By default uses bilinear interpolation.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_map),                      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_map_from_lonlat_grid_to_mesh'
    INTEGER                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_CSR
    INTEGER                                            :: vi
    INTEGER                                            :: il,iu,jl,ju
    REAL(dp)                                           :: wil,wiu,wjl,wju

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (map%is_in_use) CALL crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

    map%is_in_use = .TRUE.
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
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Fill in the CSR matrix
    DO vi = mesh%vi1, mesh%vi2

      ! Find enveloping lat-lon indices
      il  = MAX(1, MIN( grid%nlon-1, 1 + FLOOR((mesh%lon( vi) - MINVAL(grid%lon)) / (grid%lon(2) - grid%lon(1)))))
      iu  = il + 1
      wil = (grid%lon(iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
      wiu = 1._dp - wil

      ! Exception for pixels near the zero meridian
      IF (mesh%lon( vi) < MINVAL(grid%lon)) THEN
        il  = grid%nlon
        iu  = 1
        wil = (grid%lon( iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
        wiu = 1._dp - wil
      ELSEIF (mesh%lon( vi) > MAXVAL(grid%lon)) THEN
        il  = grid%nlon
        iu  = 1
        wiu = (mesh%lon( vi) - grid%lon( il)) / (grid%lon(2) - grid%lon(1))
        wil = 1._dp - wiu
      END IF

      jl  = MAX(1, MIN( grid%nlat-1, 1 + FLOOR((mesh%lat( vi) - MINVAL(grid%lat)) / (grid%lat(2) - grid%lat(1)))))
      ju  = jl + 1
      wjl = (grid%lat( ju) - mesh%lat( vi)) / (grid%lat(2) - grid%lat(1))
      wju = 1 - wjl

      ! Add values to the CSR matrix
      CALL add_entry_CSR_dist( M_CSR, vi, grid%ij2n( il,jl), wil * wjl)
      CALL add_entry_CSR_dist( M_CSR, vi, grid%ij2n( il,ju), wil * wju)
      CALL add_entry_CSR_dist( M_CSR, vi, grid%ij2n( iu,jl), wiu * wjl)
      CALL add_entry_CSR_dist( M_CSR, vi, grid%ij2n( iu,ju), wiu * wju)

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_CSR, map%M)

    ! Clean up the Fortran versions
    CALL deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_map_from_lonlat_grid_to_mesh

  SUBROUTINE create_map_from_mesh_to_mesh_nearest_neighbour( mesh_src, mesh_dst, map)
    ! Create a new mapping object from a mesh to a mesh.
    !
    ! Uses nearest-neighbour interpolation.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(type_map),                      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_map_from_mesh_to_mesh_nearest_neighbour'
    INTEGER                                            :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_CSR
    INTEGER                                            :: row, vi_dst
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER                                            :: vi_src, col

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (map%is_in_use) CALL crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

    map%is_in_use = .TRUE.
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

    CALL allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    vi_src = 1
    DO row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      CALL find_containing_vertex( mesh_src, p, vi_src)

      col = mesh_src%vi2n( vi_src)

      ! Add to the matrix
      CALL add_entry_CSR_dist( M_CSR, row, col, 1._dp)

    END DO ! DO row = mesh_dst%vi1, mesh_dst%vi2
    CALL sync

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_CSR, map%M)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_map_from_mesh_to_mesh_nearest_neighbour

  SUBROUTINE create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, map)
    ! Create a new mapping object from a mesh to a mesh.
    !
    ! Uses trilinear interpolation.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(type_map),                      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_map_from_mesh_to_mesh_trilin'
    INTEGER                                            :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_CSR
    INTEGER                                            :: row, vi_dst
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER                                            :: ti_src, via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    REAL(dp)                                           :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc, wa, wb, wc
    INTEGER                                            :: cola, colb, colc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (map%is_in_use) CALL crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

    map%is_in_use = .TRUE.
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

    CALL allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    ti_src = 1
    DO row = mesh_dst%vi1, mesh_dst%vi2

      vi_dst = mesh_dst%n2vi( row)

      p = mesh_dst%V( vi_dst,:)
      CALL find_containing_triangle( mesh_src, p, ti_src)

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
      CALL add_entry_CSR_dist( M_CSR, row, cola, wa)
      CALL add_entry_CSR_dist( M_CSR, row, colb, wb)
      CALL add_entry_CSR_dist( M_CSR, row, colc, wc)

    END DO ! DO row = mesh_dst%vi1, mesh_dst%vi2
    CALL sync

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_CSR, map%M)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_map_from_mesh_to_mesh_trilin

  SUBROUTINE create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_src, mesh_dst, map)
    ! Create a new mapping object from a mesh to a mesh.
    !
    ! Uses 2nd-order conservative interpolation.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(type_map),                      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_map_from_mesh_to_mesh_2nd_order_conservative'
    TYPE(PetscErrorCode)                               :: perr
    LOGICAL                                            :: count_coincidences
    INTEGER                                            :: nnz_per_row_max
    TYPE(tMat)                                         :: B_xdy_b_a  , B_mxydx_b_a  , B_xydy_b_a
    TYPE(tMat)                                         :: B_xdy_a_b  , B_mxydx_a_b  , B_xydy_a_b
    TYPE(tMat)                                         :: B_xdy_b_a_T, B_mxydx_b_a_T, B_xydy_b_a_T
    TYPE(tMat)                                         :: w0, w1x, w1y
    INTEGER                                            :: istart, iend, n, k, ti
    INTEGER                                            :: ncols
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: cols
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: vals, w0_row, w1x_row, w1y_row
    REAL(dp)                                           :: A_overlap_tot
    TYPE(tMat)                                         :: M_map_a_b, M_ddx_a_b, M_ddy_a_b
    TYPE(tMat)                                         :: M1, M2, M_cons_1st_order

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (map%is_in_use) CALL crash('this map is already in use!')

  ! == Initialise map metadata
  ! ==========================

    map%is_in_use = .TRUE.
    map%name_src  = mesh_src%name
    map%name_dst  = mesh_dst%name
    map%method    = '2nd_order_conservative'

    ! Integrate around the Voronoi cells of the destination mesh through the triangles of the source mesh
    count_coincidences = .TRUE.
    CALL integrate_Voronoi_cells_through_triangles( mesh_dst, mesh_src, B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)

    ! Integrate around the triangles of the source mesh through the Voronoi cells of the destination mesh
    count_coincidences = .FALSE.
    CALL integrate_triangles_through_Voronoi_cells( mesh_src, mesh_dst, B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)

    ! Transpose line integral matrices
    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - transposing line integral matrices...'
    CALL MatCreateTranspose( B_xdy_b_a  , B_xdy_b_a_T  , perr)
    CALL MatCreateTranspose( B_mxydx_b_a, B_mxydx_b_a_T, perr)
    CALL MatCreateTranspose( B_xydy_b_a , B_xydy_b_a_T , perr)

    ! Combine line integrals around areas of overlap to get surface integrals over areas of overlap
    CALL MatAXPY( B_xdy_a_b  , 1._dp, B_xdy_b_a_T  , UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( B_mxydx_a_b, 1._dp, B_mxydx_b_a_T, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( B_xydy_a_b , 1._dp, B_xydy_b_a_T , UNKNOWN_NONZERO_PATTERN, perr)

    CALL MatDestroy( B_xdy_b_a_T  , perr)
    CALL MatDestroy( B_mxydx_b_a_T, perr)
    CALL MatDestroy( B_xydy_b_a_T , perr)

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    CALL MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w0, perr)
    CALL MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1x, perr)
    CALL MatConvert( B_xdy_a_b, MATAIJ, MAT_INITIAL_MATRIX, w1y, perr)

    ! Estimate maximum number of non-zeros per row (i.e. maximum number of grid cells overlapping with a mesh triangle)
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh_src%TriA) / MINVAL( mesh_dst%A   )), &
                                    CEILING( 2._dp * MAXVAL( mesh_dst%A   ) / MINVAL( mesh_src%TriA))) )

    ! Allocate memory for a single matrix row
    ALLOCATE( cols(    nnz_per_row_max))
    ALLOCATE( vals(    nnz_per_row_max))
    ALLOCATE( w0_row(  nnz_per_row_max))
    ALLOCATE( w1x_row( nnz_per_row_max))
    ALLOCATE( w1y_row( nnz_per_row_max))

    CALL MatGetOwnershipRange( B_xdy_a_b  , istart, iend, perr)

    DO n = istart+1, iend ! +1 because PETSc indexes from 0

      ! w0
      CALL MatGetRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)
      A_overlap_tot = SUM( vals( 1:ncols))
      DO k = 1, ncols
        w0_row( k) = vals( k) / A_overlap_tot
        CALL MatSetValues( w0, 1, n-1, 1, cols( k), w0_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)

      ! w1x
      CALL MatGetRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        ti = cols( k)+1
        w1x_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,1) * w0_row( k))
        CALL MatSetValues( w1x, 1, n-1, 1, cols( k), w1x_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)

      ! w1y
      CALL MatGetRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        ti = cols( k)+1
        w1y_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,2) * w0_row( k))
        CALL MatSetValues( w1y, 1, n-1, 1, cols( k), w1y_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)

    END DO
    CALL MatAssemblyBegin( w0, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   w0, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, perr)
    CALL sync

    CALL MatDestroy( B_xdy_a_b  , perr)
    CALL MatDestroy( B_mxydx_a_b, perr)
    CALL MatDestroy( B_xydy_a_b , perr)

  ! == Calculate the remapping matrices
  ! ===================================

    ! Safety
    IF (.NOT. ALLOCATED( mesh_src%vi2n)) THEN
      CALL crash('matrix operators for mesh "' // TRIM( mesh_src%name) // '" have not been calculated!')
    END IF

    ! Convert matrices to PETSc format
    CALL mat_CSR2petsc( mesh_src%M_map_a_b, M_map_a_b)
    CALL mat_CSR2petsc( mesh_src%M_ddx_a_b, M_ddx_a_b)
    CALL mat_CSR2petsc( mesh_src%M_ddy_a_b, M_ddy_a_b)

    ! 1st-order = w0 * map_a_b
    CALL MatMatMult( w0 , M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_cons_1st_order, perr)

    ! 2nd-order = 1st-order + w1x * ddx_a_b + w1y * ddy_a_b
    CALL MatMatMult( w1x, M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    CALL MatMatMult( w1y, M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2, perr)

    CALL MatConvert( M_cons_1st_order, MATAIJ, MAT_INITIAL_MATRIX, map%M, perr)
    CALL MatAXPY( map%M, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( map%M, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, perr)

    CALL MatDestroy( w0       , perr)
    CALL MatDestroy( w1x      , perr)
    CALL MatDestroy( w1y      , perr)
    CALL MatDestroy( M_map_a_b, perr)
    CALL MatDestroy( M_ddx_a_b, perr)
    CALL MatDestroy( M_ddy_a_b, perr)

    CALL MatDestroy( M1              , perr)
    CALL MatDestroy( M2              , perr)

  ! == Apply some final corrections
  ! ===============================

    CALL correct_mesh_to_mesh_map( mesh_src, mesh_dst, M_cons_1st_order, map%M)

    CALL MatDestroy( M_cons_1st_order, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_map_from_mesh_to_mesh_2nd_order_conservative

  SUBROUTINE correct_mesh_to_mesh_map( mesh_src, mesh_dst, M_cons_1st_order, M_cons_2nd_order)
    ! Apply some final corrections to the 2nd-order conservative mesh-to-mesh remapping operator:
    ! - set remapped data to zero on the domain border
    ! - use direct copying for identical vertices

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(tMat),                          INTENT(IN)    :: M_cons_1st_order
    TYPE(tMat),                          INTENT(INOUT) :: M_cons_2nd_order

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'correct_mesh_to_mesh_map'
    INTEGER                                            :: perr
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_cons_1st_order_CSR
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_cons_2nd_order_CSR
    INTEGER                                            :: i, vi_dst, k1, k2, k
    INTEGER                                            :: j, vi_src
    LOGICAL                                            :: do_direct_copy
    INTEGER                                            :: vi_src_copy
    LOGICAL                                            :: Voronoi_cells_are_identical
    REAL(dp), DIMENSION( mesh_src%nC_mem,2)            :: Vor_src   , Vor_dst
    INTEGER,  DIMENSION( mesh_src%nC_mem  )            :: Vor_src_vi, Vor_dst_vi
    INTEGER,  DIMENSION( mesh_src%nC_mem  )            :: Vor_src_ti, Vor_dst_ti
    INTEGER                                            :: nVor_src  , nVor_dst
    INTEGER                                            :: vori
    TYPE(type_map)                                     :: map_trilin
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_trilin_CSR
    LOGICAL,  DIMENSION( mesh_dst%nV)                  :: isgood_1st_order
    LOGICAL,  DIMENSION( mesh_dst%nV)                  :: isgood_2nd_order
    INTEGER                                            :: kk1,kk2,kk

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Convert matrices to CSR format for easier handling
    CALL mat_petsc2CSR( M_cons_1st_order, M_cons_1st_order_CSR)
    CALL mat_petsc2CSR( M_cons_2nd_order, M_cons_2nd_order_CSR)

  ! == Set to zero
  !
  ! 2nd-order conservative doesn't work all that well on the
  ! domain border; just set the result to zero there.

    DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      IF (mesh_dst%VBI( vi_dst) > 0) THEN
        DO k = k1, k2
          M_cons_2nd_order_CSR%val( k) = 0._dp
        END DO
      END IF

    END DO ! DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

  ! == Direct copying
  !
  ! With the new mesh generation code (of UFE2.0), many vertices away from the grounding line
  ! remain unchanged after a mesh update. The vertex-to-triangle-to-vertex remapping is slightly
  ! diffusive, so instead we can just copy data directly for those vertices.

    DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      do_direct_copy = .FALSE.
      vi_src_copy    = 0

      ! Loop over all src vertices contributing to this dst vertex
      DO k = k1, k2

        j = M_cons_2nd_order_CSR%ind( k)
        vi_src = mesh_src%n2vi( j)

        IF (NORM2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%R( vi_dst) / 1E2_dp) THEN
          ! These vertices coincide; check if their Voronoi cells are identical

          Voronoi_cells_are_identical = .TRUE.

          CALL calc_Voronoi_cell( mesh_src, vi_src, 0._dp, Vor_src, Vor_src_vi, Vor_src_ti, nVor_src)
          CALL calc_Voronoi_cell( mesh_dst, vi_dst, 0._dp, Vor_dst, Vor_dst_vi, Vor_dst_ti, nVor_dst)

          IF (nVor_src /= nVor_dst) THEN
            Voronoi_cells_are_identical = .FALSE.
          ELSE
            DO vori = 1, nVor_src
              IF (NORM2( Vor_src( vori,:) - Vor_dst( vori,:)) > mesh_dst%R( vi_dst) / 1E2_dp) THEN
                Voronoi_cells_are_identical = .FALSE.
              END IF
            END DO
          END IF

          IF (Voronoi_cells_are_identical) THEN
            ! These two vertices have identical Voronoi cells; use direct copying
            do_direct_copy = .TRUE.
            vi_src_copy    = vi_src
            EXIT
          END IF ! IF (Voronoi_cells_are_identical) THEN

        END IF ! IF (NORM2( mesh_dst%V( vi_dst,:) - mesh_src%V( vi_src,:)) < mesh_dst%tol_dist) THEN

      END DO ! DO k = k1, k2

      ! If a source vertex with an identical Voronoi cell to this dst vertex was
      ! found, copy data from that vertex directly
      IF (do_direct_copy) THEN
        ! Loop over all src vertices contributing to this dst vertex; set all
        ! contributions to zero except the one we're copying (which is set to 1)

        DO k = k1, k2

          j = M_cons_2nd_order_CSR%ind( k)
          vi_src = mesh_src%n2vi( j)

          IF (vi_src == vi_src_copy) THEN
            M_cons_2nd_order_CSR%val( k) = 1._dp
          ELSE
            M_cons_2nd_order_CSR%val( k) = 0._dp
          END IF

        END DO ! DO k = k1, k2

      END IF ! IF (do_direct_copy) THEN

    END DO ! DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

  ! == Remapping errors
  !
  ! On very rare occasions, the remapping operator is just wrong, likely due to round-off
  ! errors in determining if vertices coincide or not. Usually, the 1st-order operator
  ! is fine. If that one fails too, just replace the answer with trilinear interpolation.
  !
  ! Faulty operators can be detected by negative coefficients in the remapping matrix,
  ! which automatically violate conservation of extreme values.

    ! Calculate the trilinear interpolation operator to serve as a back-up
    CALL create_map_from_mesh_to_mesh_trilin( mesh_src, mesh_dst, map_trilin)
    CALL mat_petsc2CSR( map_trilin%M, M_trilin_CSR)

    ! Find faulty operators in the 1st-order conservative remapping operator
    isgood_1st_order = .TRUE.
    DO i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

      k1 = M_cons_1st_order_CSR%ptr( i)
      k2 = M_cons_1st_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      DO k = k1, k2
        IF (M_cons_1st_order_CSR%val( k) < 0._dp) THEN
          isgood_1st_order( vi_dst) = .FALSE.
        END IF
      END DO

    END DO ! DO i = M_cons_1st_order_CSR%i1, M_cons_1st_order_CSR%i2

    ! Find faulty operators in the 2nd-order conservative remapping operator
    isgood_2nd_order = .TRUE.
    DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      k1 = M_cons_2nd_order_CSR%ptr( i)
      k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1

      vi_dst = mesh_dst%n2vi( i)

      DO k = k1, k2
        IF (M_cons_2nd_order_CSR%val( k) < 0._dp) THEN
          isgood_2nd_order( vi_dst) = .FALSE.
        END IF
      END DO

    END DO ! DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Replace faulty operators in the 2nd-order conservative remapping operator

    DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

      vi_dst = mesh_dst%n2vi( i)

      IF (.NOT. isgood_2nd_order( vi_dst)) THEN
        ! Replace this faulty operator

        IF (isgood_1st_order( vi_dst)) THEN
          ! Replace with the 1st-order conservative remapping operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          DO k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          END DO

          ! Then copy values from the 1st-order operator
          kk1 = M_cons_1st_order_CSR%ptr( i)
          kk2 = M_cons_1st_order_CSR%ptr( i+1) - 1
          DO kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_cons_1st_order_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_cons_1st_order_CSR%val( kk)
          END DO

        ELSE ! IF (isgood_1st_order( vi_dst)) THEN
          ! Replace with the trilinear interpolation operator

          ! First set all coefficients of the 2nd-order operator to zero
          k1 = M_cons_2nd_order_CSR%ptr( i)
          k2 = M_cons_2nd_order_CSR%ptr( i+1) - 1
          DO k = k1, k2
            M_cons_2nd_order_CSR%val( k) = 0._dp
          END DO

          ! Then copy values from the trilinear interpolation operator
          kk1 = M_trilin_CSR%ptr( i)
          kk2 = M_trilin_CSR%ptr( i+1) - 1
          DO kk = kk1, kk2
            k = kk - kk1 + k1
            M_cons_2nd_order_CSR%ind( k) = M_trilin_CSR%ind( kk)
            M_cons_2nd_order_CSR%val( k) = M_trilin_CSR%val( kk)
          END DO

        END IF ! IF (isgood_1st_order( vi_dst)) THEN

      END IF ! IF (.NOT. isgood_2nd_order( vi_dst)) THEN

    END DO ! DO i = M_cons_2nd_order_CSR%i1, M_cons_2nd_order_CSR%i2

    ! Convert back to PETSc format
    CALL MatDestroy( M_cons_2nd_order, perr)
    CALL mat_CSR2petsc( M_cons_2nd_order_CSR, M_cons_2nd_order)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_cons_1st_order_CSR)
    CALL deallocate_matrix_CSR_dist( M_cons_2nd_order_CSR)
    CALL deallocate_matrix_CSR_dist( M_trilin_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_mesh_to_mesh_map

! == Routines used in creating remapping matrices
! ===============================================

  ! Integrate around triangles/Voronoi cells through triangles/Voronoi cells
  SUBROUTINE integrate_triangles_through_Voronoi_cells( mesh_tri, mesh_Vor, B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)
    ! Integrate around the triangles of mesh_tri through the Voronoi cells of mesh_Vor

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_tri
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_Vor
    TYPE(tMat),                          INTENT(OUT)   :: B_xdy_b_a
    TYPE(tMat),                          INTENT(OUT)   :: B_mxydx_b_a
    TYPE(tMat),                          INTENT(OUT)   :: B_xydy_b_a
    LOGICAL,                             INTENT(IN)    :: count_coincidences

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'integrate_triangles_through_Voronoi_cells'
    INTEGER                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: B_xdy_b_a_CSR, B_mxydx_b_a_CSR, B_xydy_b_a_CSR
    TYPE(type_single_row_mapping_matrices)             :: single_row
    INTEGER                                            :: via, vib, vic, ti, vi_hint, k
    REAL(dp), DIMENSION(2)                             :: p, q

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
  ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_tri%nTri  ! to
    nrows_loc       = mesh_tri%nTri_loc
    ncols           = mesh_Vor%nV    ! from
    ncols_loc       = mesh_Vor%nV_loc
    nnz_est         = 4 * MAX( nrows, ncols)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh_tri%TriA) / MINVAL( mesh_Vor%A   )), &
                                    CEILING( 2._dp * MAXVAL( mesh_Vor%A   ) / MINVAL( mesh_tri%TriA)) ))

    CALL allocate_matrix_CSR_dist( B_xdy_b_a_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_mxydx_b_a_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_xydy_b_a_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    ALLOCATE( single_row%index_left( single_row%n_max))
    ALLOCATE( single_row%LI_xdy(     single_row%n_max))
    ALLOCATE( single_row%LI_mxydx(   single_row%n_max))
    ALLOCATE( single_row%LI_xydy(    single_row%n_max))

  ! == Trace all the line segments to fill the matrices
  ! ===================================================

    vi_hint = 1

    DO ti = mesh_tri%ti1, mesh_tri%ti2

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
      CALL trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vib,:)
      q = mesh_tri%V( vic,:)
      CALL trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vic,:)
      q = mesh_tri%V( via,:)
      CALL trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      ! Add the results for this triangle to the sparse matrix
      IF (single_row%n == 0) THEN
        CALL add_empty_row_CSR_dist( B_xdy_b_a_CSR  , ti)
        CALL add_empty_row_CSR_dist( B_mxydx_b_a_CSR, ti)
        CALL add_empty_row_CSR_dist( B_xydy_b_a_CSR , ti)
      ELSE
        DO k = 1, single_row%n
          CALL add_entry_CSR_dist( B_xdy_b_a_CSR  , ti, single_row%index_left( k), single_row%LI_xdy(   k))
          CALL add_entry_CSR_dist( B_mxydx_b_a_CSR, ti, single_row%index_left( k), single_row%LI_mxydx( k))
          CALL add_entry_CSR_dist( B_xydy_b_a_CSR , ti, single_row%index_left( k), single_row%LI_xydy(  k))
        END DO
      END IF

    END DO ! DO ti = mesh_tri%ti1, mesh_tri%ti2
    CALL sync

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( B_xdy_b_a_CSR  , B_xdy_b_a  )
    CALL mat_CSR2petsc( B_mxydx_b_a_CSR, B_mxydx_b_a)
    CALL mat_CSR2petsc( B_xydy_b_a_CSR , B_xydy_b_a )

    ! Clean up the Fortran versions
    CALL deallocate_matrix_CSR_dist( B_xdy_b_a_CSR  )
    CALL deallocate_matrix_CSR_dist( B_mxydx_b_a_CSR)
    CALL deallocate_matrix_CSR_dist( B_xydy_b_a_CSR )

    ! Clean up after yourself
    DEALLOCATE( single_row%index_left )
    DEALLOCATE( single_row%LI_xdy     )
    DEALLOCATE( single_row%LI_mxydx   )
    DEALLOCATE( single_row%LI_xydy    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_triangles_through_Voronoi_cells

  SUBROUTINE integrate_Voronoi_cells_through_triangles( mesh_Vor, mesh_tri, B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)
    ! Integrate around the grid cells of the grid through the triangles of the mesh

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_Vor
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_tri
    TYPE(tMat),                          INTENT(OUT)   :: B_xdy_a_b
    TYPE(tMat),                          INTENT(OUT)   :: B_mxydx_a_b
    TYPE(tMat),                          INTENT(OUT)   :: B_xydy_a_b
    LOGICAL,                             INTENT(IN)    :: count_coincidences

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'integrate_Voronoi_cells_through_triangles'
    INTEGER                                            :: nrows, ncols, nrows_loc, ncols_loc, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: B_xdy_a_b_CSR, B_mxydx_a_b_CSR, B_xydy_a_b_CSR
    TYPE(type_single_row_mapping_matrices)             :: single_row
    INTEGER                                            :: vi, vori1, vori2, k, ti_hint
    REAL(dp), DIMENSION( mesh_Vor%nC_mem,2)            :: Vor
    INTEGER,  DIMENSION( mesh_Vor%nC_mem  )            :: Vor_vi
    INTEGER,  DIMENSION( mesh_Vor%nC_mem  )            :: Vor_ti
    INTEGER                                            :: nVor
    REAL(dp), DIMENSION(2)                             :: p, q


    ! Add routine to path
    CALL init_routine( routine_name)
  ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
  ! ===========================================================================

    ! Matrix sise
    nrows           = mesh_Vor%nV    ! to
    nrows_loc       = mesh_Vor%nV_loc
    ncols           = mesh_tri%nTri  ! from
    ncols_loc       = mesh_tri%nTri_loc
    nnz_est         = 4 * MAX( nrows, ncols)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh_tri%TriA) / MINVAL( mesh_vor%A   )), &
                                    CEILING( 2._dp * MAXVAL( mesh_vor%A   ) / MINVAL( mesh_tri%TriA)) ))

    CALL allocate_matrix_CSR_dist( B_xdy_a_b_CSR  , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_mxydx_a_b_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_xydy_a_b_CSR , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    ALLOCATE( single_row%index_left( single_row%n_max))
    ALLOCATE( single_row%LI_xdy(     single_row%n_max))
    ALLOCATE( single_row%LI_mxydx(   single_row%n_max))
    ALLOCATE( single_row%LI_xydy(    single_row%n_max))

  ! == Trace all the line segments to fill the matrices
  ! ===================================================

    ti_hint = 1

    DO vi = mesh_Vor%vi1, mesh_Vor%vi2 ! +1 because PETSc indexes from 0

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! Integrate over the complete Voronoi cell boundary
      CALL calc_Voronoi_cell( mesh_Vor, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)
      DO vori1 = 1, nVor
        vori2 = vori1 + 1
        IF (vori2 > nVor) vori2 = 1
        p = Vor( vori1,:)
        q = Vor( vori2,:)
        CALL trace_line_tri( mesh_tri, p, q, single_row, count_coincidences, ti_hint)
      END DO

      ! Add the results for this triangle to the sparse matrix
      IF (single_row%n == 0) THEN
        CALL add_empty_row_CSR_dist( B_xdy_a_b_CSR  , vi)
        CALL add_empty_row_CSR_dist( B_mxydx_a_b_CSR, vi)
        CALL add_empty_row_CSR_dist( B_xydy_a_b_CSR , vi)
      ELSE
        DO k = 1, single_row%n
          CALL add_entry_CSR_dist( B_xdy_a_b_CSR  , vi, single_row%index_left( k), single_row%LI_xdy(   k))
          CALL add_entry_CSR_dist( B_mxydx_a_b_CSR, vi, single_row%index_left( k), single_row%LI_mxydx( k))
          CALL add_entry_CSR_dist( B_xydy_a_b_CSR , vi, single_row%index_left( k), single_row%LI_xydy(  k))
        END DO
      END IF

    END DO ! DO vi = mesh_Vor%vi1, mesh_Vor%vi2
    CALL sync

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( B_xdy_a_b_CSR  , B_xdy_a_b  )
    CALL mat_CSR2petsc( B_mxydx_a_b_CSR, B_mxydx_a_b)
    CALL mat_CSR2petsc( B_xydy_a_b_CSR , B_xydy_a_b )

    ! Clean up the Fortran versions
    CALL deallocate_matrix_CSR_dist( B_xdy_a_b_CSR  )
    CALL deallocate_matrix_CSR_dist( B_mxydx_a_b_CSR)
    CALL deallocate_matrix_CSR_dist( B_xydy_a_b_CSR )

    ! Clean up after yourself
    DEALLOCATE( single_row%index_left )
    DEALLOCATE( single_row%LI_xdy     )
    DEALLOCATE( single_row%LI_mxydx   )
    DEALLOCATE( single_row%LI_xydy    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_Voronoi_cells_through_triangles

  ! Line tracing algorithm through mesh triangles
  SUBROUTINE trace_line_tri( mesh, p, q, single_row, count_coincidences, ti_hint)
    ! Trace the line [pq] through the triangles of the mesh and calculate
    ! the three line integrals for the line segments inside the different triangles

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    LOGICAL,                             INTENT(IN)    :: count_coincidences
    INTEGER,                             INTENT(INOUT) :: ti_hint

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'trace_line_tri'
    REAL(dp), DIMENSION(2)                             :: pp, qq
    LOGICAL                                            :: is_valid_line
    LOGICAL                                            :: finished
    INTEGER                                            :: n_cycles
    INTEGER                                            :: ti_in, vi_on, ei_on
    REAL(dp), DIMENSION(2)                             :: p_next
    INTEGER                                            :: ti_left
    LOGICAL                                            :: coincides
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the mesh domain
    CALL crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    IF (.NOT. is_valid_line) THEN
      ! [pq] doesn't pass through the mesh domain anywhere
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    CALL trace_line_tri_start( mesh, pp, ti_hint, ti_in, vi_on, ei_on)

    ! Iteratively trace the line through the mesh
    finished = .FALSE.
    n_cycles = 0
    DO WHILE (.NOT. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      IF     (ti_in  > 0) THEN
        ! p lies inside triangle ti_in
        CALL trace_line_tri_ti( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      ELSEIF (vi_on  > 0) THEN
        ! p lies on vertex vi_on
        CALL trace_line_tri_vi( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      ELSEIF (ei_on > 0) THEN
        ! p lies on edge ei_on
        CALL trace_line_tri_ei( mesh, pp, qq, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
      ELSE
        CALL crash('coincidence indicators dont make sense!')
      END IF

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, mesh%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, mesh%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, mesh%tol_dist)

      ! Add them to the results structure
      IF (NORM2( p_next - pp) > mesh%tol_dist) THEN
        CALL add_integrals_to_single_row( single_row, ti_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      END IF

      ! Cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      IF (n_cycles > mesh%nV) THEN
        CALL crash('trace_line_tri - iterative tracer got stuck!')
      END IF

      ! Update ti_hint, for more efficiency
      ti_hint = ti_left

    END DO ! DO WHILE (.NOT. finished)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE trace_line_tri

  SUBROUTINE trace_line_tri_start( mesh, p, ti_hint, ti_in, vi_on, ei_on)
    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside triangle ti_in, ...
    !    - lies on vertex vi_on, or...
    !    - lies on edge ei_on

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    INTEGER,                             INTENT(INOUT) :: ti_hint
    INTEGER,                             INTENT(OUT)   :: ti_in
    INTEGER,                             INTENT(OUT)   :: vi_on
    INTEGER,                             INTENT(OUT)   :: ei_on

    ! Local variables:
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    INTEGER                                            :: vvi, vj, ei

    ! Initialise
    ti_in  = 0
    vi_on  = 0
    ei_on = 0

    ! Find the triangle containing p
    CALL find_containing_triangle( mesh, p, ti_hint)

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_hint,1)
    vib = mesh%Tri( ti_hint,2)
    vic = mesh%Tri( ti_hint,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Check if p lies on any of the three vertices
    IF     (NORM2( pa - p) < mesh%tol_dist) THEN
      ! p lies on via
      vi_on = via
      RETURN
    ELSEIF (NORM2( pb - p) < mesh%tol_dist) THEN
      ! p lies on vib
      vi_on = vib
      RETURN
    ELSEIF (NORM2( pc - p) < mesh%tol_dist) THEN
      ! p lies on vic
      vi_on = vic
      RETURN
    END IF

    ! Check if p lies on any of the three edges
    IF     (lies_on_line_segment( pa, pb, p, mesh%tol_dist)) THEN
      ! p lies on the edge connecting via and vib
      DO vvi = 1, mesh%nC( via)
        vj = mesh%C(  via,vvi)
        ei = mesh%VE( via,vvi)
        IF (vj == vib) THEN
          ei_on = ei
          RETURN
        END IF
      END DO
    ELSEIF (lies_on_line_segment( pb, pc, p, mesh%tol_dist)) THEN
      ! p lies on the edge connecting vib and vic
      DO vvi = 1, mesh%nC( vib)
        vj = mesh%C(  vib,vvi)
        ei = mesh%VE( vib,vvi)
        IF (vj == vic) THEN
          ei_on = ei
          RETURN
        END IF
      END DO
    ELSEIF (lies_on_line_segment( pc, pa, p, mesh%tol_dist)) THEN
      ! p lies on the edge connecting vic and via
      DO vvi = 1, mesh%nC( vic)
        vj = mesh%C(  vic,vvi)
        ei = mesh%VE( vic,vvi)
        IF (vj == via) THEN
          ei_on = ei
          RETURN
        END IF
      END DO
    END IF

    ! If p lies not on the vertices or edges of the triangle, then it must lie inside of it
    ti_in = ti_hint

  END SUBROUTINE trace_line_tri_start

  SUBROUTINE trace_line_tri_ti(  mesh, p, q, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies inside triangle ti_in,
    ! find the point p_next where [pq] crosses into the next triangle.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: ti_in
    INTEGER,                             INTENT(INOUT) :: vi_on
    INTEGER,                             INTENT(INOUT) :: ei_on
    INTEGER,                             INTENT(OUT)   :: ti_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    INTEGER                                            :: vvi, vj, ei
    REAL(dp), DIMENSION(2)                             :: llis
    LOGICAL                                            :: do_cross

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_in,1)
    vib = mesh%Tri( ti_in,2)
    vic = mesh%Tri( ti_in,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Safety
    IF (ti_in == 0 .OR. vi_on > 0 .OR. ei_on > 0) THEN
      CALL crash('trace_line_tri_ti - coincidence indicators dont make sense!')
    END IF
    IF (.NOT. is_in_triangle( pa, pb, pc, p)) THEN
      CALL crash('trace_line_tri_ti - p does not lie inside triangle ti_in!')
    END IF

    ! Check if q lies inside the same triangle
    IF (is_in_triangle( pa, pb, pc, q)) THEN
      ! q lies inside the same triangle
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on vertex via
    IF (NORM2( pa - q) < mesh%tol_dist) THEN
      ! q lies on vertex via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on vertex vib
    IF (NORM2( pb - q) < mesh%tol_dist) THEN
      ! q lies on vertex vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on vertex vic
    IF (NORM2( pc - q) < mesh%tol_dist) THEN
      ! q lies on vertex vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on edge via-vib
    IF (lies_on_line_segment( pa, pb, q, mesh%tol_dist)) THEN
      ! q lies on edge via-vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on edge vib-vic
    IF (lies_on_line_segment( pb, pc, q, mesh%tol_dist)) THEN
      ! q lies on edge vib-vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on edge vic-via
    IF (lies_on_line_segment( pc, pa, q, mesh%tol_dist)) THEN
      ! q lies on edge vic-via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through via
    IF (lies_on_line_segment( p, q, pa, mesh%tol_dist)) THEN
      ! [pq] passes through via
      p_next    = pa
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = via
      ei_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through vib
    IF (lies_on_line_segment( p, q, pb, mesh%tol_dist)) THEN
      ! [pq] passes through vib
      p_next    = pb
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vib
      ei_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through vic
    IF (lies_on_line_segment( p, q, pc, mesh%tol_dist)) THEN
      ! [pq] passes through vic
      p_next    = pc
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vic
      ei_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] crosses edge via-vib
    CALL segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
    IF (do_cross) THEN
      ! [pq] crosses edge [via,vib]
      ! Find the edge connecting via and vib
      DO vvi = 1, mesh%nC( via)
        vj = mesh%C(  via,vvi)
        ei = mesh%VE( via,vvi)
        IF (vj == vib) THEN
          ei_on = ei
          EXIT
        END IF
      END DO
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] crosses edge vib-vic
    CALL segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
    IF (do_cross) THEN
      ! [pq] crosses edge [vib,vic]
      ! Find the edge connecting vib and vic
      DO vvi = 1, mesh%nC( vib)
        vj = mesh%C(  vib,vvi)
        ei = mesh%VE( vib,vvi)
        IF (vj == vic) THEN
          ei_on = ei
          EXIT
        END IF
      END DO
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] crosses edge vic-via
    CALL segment_intersection( p, q, pc, pa, llis, do_cross, mesh%tol_dist)
    IF (do_cross) THEN
      ! [pq] crosses edge [vic,via]
      ! Find the edge connecting vic and via
      DO vvi = 1, mesh%nC( vic)
        vj = mesh%C(  vic,vvi)
        ei = mesh%VE( vic,vvi)
        IF (vj == via) THEN
          ei_on = ei
          EXIT
        END IF
      END DO
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_tri_ti - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_tri_ti

  SUBROUTINE trace_line_tri_vi(  mesh, p, q, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies on vertex vi_on,
    ! find the point p_next where [pq] crosses into the next triangle.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: ti_in
    INTEGER,                             INTENT(INOUT) :: vi_on
    INTEGER,                             INTENT(INOUT) :: ei_on
    INTEGER,                             INTENT(OUT)   :: ti_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc, pv
    INTEGER                                            :: vvi, vj, ei, n1, n2, n3, vti, ti
    REAL(dp), DIMENSION(2)                             :: llis
    LOGICAL                                            :: do_cross

    ! Safety
    IF (ti_in > 0 .OR. vi_on == 0 .OR. ei_on > 0) THEN
      CALL crash('trace_line_tri_vi - coincidence indicators dont make sense!')
    END IF
    IF (NORM2( p - mesh%V( vi_on,:)) > mesh%tol_dist) THEN
      CALL crash('trace_line_tri_vi - p does not lie on vertex vi_on!')
    END IF

    ! Check if q lies on any of the edges originating in this vertex
    DO vvi = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,vvi)
      ei = mesh%VE( vi_on,vvi)
      pv = mesh%V( vj,:)
      IF (NORM2( mesh%V( vj,:) - q) < mesh%tol_dist .OR. &
          lies_on_line_segment( p, pv, q, mesh%tol_dist)) THEN
        ! q lies on edge ei, connecting vi_on and vj
        IF (mesh%EV( ei,1) == vi_on) THEN
          ti_left = mesh%ETri( ei,1)
        ELSE
          ti_left = mesh%ETri( ei,2)
        END IF
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        coincides = .TRUE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check if q lies inside any of the triangles surrounding vi_on
    DO vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)
      IF (is_in_triangle( pa, pb, pc, q) .OR. &
          lies_on_line_segment( pa, pb, q, mesh%tol_dist) .OR. &
          lies_on_line_segment( pb, pc, q, mesh%tol_dist) .OR. &
          lies_on_line_segment( pc, pa, q, mesh%tol_dist)) THEN
        ! q lies inside adjacent triangle ti
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = ti
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check if [pq] passes through any of the neighbouring vertices
    DO vvi = 1, mesh%nC( vi_on)
      vj = mesh%C(  vi_on,vvi)
      ei = mesh%VE( vi_on,vvi)
      pv = mesh%V( vj,:)
      IF (lies_on_line_segment( p, q, pv, mesh%tol_dist)) THEN
        ! [pq] passes through neighbouring vertex vj, which is connected to vi_on by edge ei
        p_next    = pv
        ti_in     = 0
        vi_on     = vj
        ei_on     = 0
        IF (mesh%EV( ei,1) == vi_on) THEN
          ti_left = mesh%ETri( ei,1)
        ELSE
          ti_left = mesh%ETri( ei,2)
        END IF
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check if [pq] exits into any of the adjacent triangles
    DO vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        n3 = n2 + 1
        IF (n3 == 4) n3 = 1
        IF (mesh%Tri( ti,n1) == vi_on) THEN
          vib = mesh%Tri( ti,n2)
          vic = mesh%Tri( ti,n3)
          pb  = mesh%V( vib,:)
          pc  = mesh%V( vic,:)
          ! Find the opposite triangle edge
          ei = 0
          DO vvi = 1, mesh%nC( vib)
            vj = mesh%C( vib,vvi)
            IF (vj == vic) THEN
              ei = mesh%VE( vib,vvi)
              EXIT
            END IF
          END DO
          CALL segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
          IF (do_cross) THEN
            ! [pq] exits triangle ti through the opposite edge ei
            p_next    = llis
            ti_in     = 0
            vi_on     = 0
            ei_on     = ei
            ti_left   = ti
            coincides = .FALSE.
            finished  = .FALSE.
            RETURN
          END IF
        END IF
      END DO
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_tri_vi - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_tri_vi

  SUBROUTINE trace_line_tri_ei( mesh, p, q, p_next, ti_in, vi_on, ei_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies on edge ei,
    ! find the point p_next where [pq] crosses into the next triangle.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: ti_in
    INTEGER,                             INTENT(INOUT) :: vi_on
    INTEGER,                             INTENT(INOUT) :: ei_on
    INTEGER,                             INTENT(OUT)   :: ti_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vil, vir, til, tir
    REAL(dp), DIMENSION(2)                             :: pa, pb, pl, pr
    INTEGER                                            :: vvi, vj, ei
    REAL(dp), DIMENSION(2)                             :: llis
    LOGICAL                                            :: do_cross

    ! Some more info about this edge
    via = mesh%EV(   ei_on,1)
    vib = mesh%EV(   ei_on,2)
    vil = mesh%EV(   ei_on,3)
    vir = mesh%EV(   ei_on,4)
    til = mesh%ETri( ei_on,1)
    tir = mesh%ETri( ei_on,2)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    IF (vil > 0) pl  = mesh%V( vil,:)
    IF (vir > 0) pr  = mesh%V( vir,:)

    ! Safety
    IF (ti_in > 0 .OR. vi_on > 0 .OR. ei_on == 0) THEN
      CALL crash('trace_line_tri_ei - coincidence indicators dont make sense!')
    END IF
    IF (.NOT. lies_on_line_segment( pa, pb, p, mesh%tol_dist)) THEN
      CALL crash('trace_line_tri_ei - p does not lie on edge ei_on!')
    END IF

    ! Check if q lies on the same edge in the direction of via
    IF (lies_on_line_segment( p, pa, q, mesh%tol_dist)) THEN
      ! q lies on the same edge in the direction of via
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      ti_left   = tir
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the same edge in the direction of vib
    IF (lies_on_line_segment( p, pb, q, mesh%tol_dist)) THEN
      ! q lies on the same edge in the direction of vib
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      ei_on     = 0
      ti_left   = til
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside either of the two adjacent triangles
    IF (til > 0) THEN
      IF (is_in_triangle( pa, pb, pl, q)) THEN
        ! q lies inside triangle til
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END IF
    IF (tir > 0) THEN
      IF (is_in_triangle( pa, pr, pb, q)) THEN
        ! q lies inside triangle tir
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        ei_on     = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END IF

    ! Check if [pq] passes through pa
    IF (lies_on_line_segment( p, q, pa, mesh%tol_dist)) THEN
      ! [pq] passes through pa
      p_next    = pa
      ti_in     = 0
      vi_on     = via
      ei_on     = 0
      ti_left   = tir
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through pb
    IF (lies_on_line_segment( p, q, pb, mesh%tol_dist)) THEN
      ! [pq] passes through pb
      p_next    = pb
      ti_in     = 0
      vi_on     = vib
      ei_on     = 0
      ti_left   = til
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through pl
    IF (til > 0) THEN
      IF (lies_on_line_segment( p, q, pl, mesh%tol_dist)) THEN
        ! [pq] passes through pl
        p_next    = pl
        ti_in     = 0
        vi_on     = vil
        ei_on     = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check if [pq] passes through pr
    IF (tir > 0) THEN
      IF (lies_on_line_segment( p, q, pr, mesh%tol_dist)) THEN
        ! [pq] passes through pr
        p_next    = pr
        ti_in     = 0
        vi_on     = vir
        ei_on     = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check if [pq] crosses edge [via,vil]
    IF (til > 0) THEN
      CALL segment_intersection( p, q, pa, pl, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [via,vil]
        ! Find the edge connecting via and vil
        DO vvi = 1, mesh%nC( via)
          vj = mesh%C(  via,vvi)
          ei = mesh%VE( via,vvi)
          IF (vj == vil) THEN
            ei_on = ei
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check if [pq] crosses edge [vil,vib]
    IF (til > 0) THEN
      CALL segment_intersection( p, q, pl, pb, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [vil,vib]
        ! Find the edge connecting vil and vib
        DO vvi = 1, mesh%nC( vil)
          vj = mesh%C(  vil,vvi)
          ei = mesh%VE( vil,vvi)
          IF (vj == vib) THEN
            ei_on = ei
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check if [pq] crosses edge [via,vir]
    IF (tir > 0) THEN
      CALL segment_intersection( p, q, pa, pr, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [via,vir]
        ! Find the edge connecting via and vir
        DO vvi = 1, mesh%nC( via)
          vj  = mesh%C(    via,vvi)
          ei = mesh%VE( via,vvi)
          IF (vj == vir) THEN
            ei_on = ei
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check if [pq] crosses edge [vir,vib]
    IF (tir > 0) THEN
      CALL segment_intersection( p, q, pr, pb, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [vir,vib]
        ! Find the edge connecting vir and vib
        DO vvi = 1, mesh%nC( vir)
          vj = mesh%C(  vir,vvi)
          ei = mesh%VE( vir,vvi)
          IF (vj == vib) THEN
            ei_on = ei
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_tri_ei - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_tri_ei

  ! Line tracing algorithm through mesh Voronoi cells
  SUBROUTINE trace_line_Vor( mesh, p, q, single_row, count_coincidences, vi_hint)
  ! Trace the line [pq] through the Voronoi cells of the mesh and calculate
  ! the three line integrals for the line segments inside the different Voronoi cells

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    LOGICAL,                             INTENT(IN)    :: count_coincidences
    INTEGER,                             INTENT(INOUT) :: vi_hint

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'trace_line_Vor'
    REAL(dp), DIMENSION(2)                             :: pp,qq
    LOGICAL                                            :: is_valid_line
    LOGICAL                                            :: finished
    INTEGER                                            :: n_cycles
    INTEGER                                            :: vi_in, ti_on, ei_on
    REAL(dp), DIMENSION(2)                             :: p_next
    INTEGER                                            :: vi_left
    LOGICAL                                            :: coincides
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the mesh domain
    CALL crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    IF (.NOT.is_valid_line) THEN
      ! [pq] doesn't pass through the mesh domain anywhere
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on
    CALL trace_line_Vor_start( mesh, pp, vi_hint, vi_in, ti_on, ei_on)

    ! Iteratively trace the line through the mesh
    finished = .FALSE.
    n_cycles = 0
    DO WHILE (.NOT.finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      IF     (vi_in  > 0) THEN
        ! p lies inside the Voronoi cell of vertex vi_in
        CALL trace_line_Vor_vi( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      ELSEIF (ti_on  > 0) THEN
        ! p lies on the circumcentre of triangle ti_on
        CALL trace_line_Vor_ti( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      ELSEIF (ei_on > 0) THEN
        ! p lies on the shared Voronoi cell boundary represented by edge ei_on
        CALL trace_line_Vor_ei( mesh, pp, qq, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
      END IF

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, mesh%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, mesh%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, mesh%tol_dist)

      ! Add them to the results structure
      IF (NORM2( p_next - pp) > mesh%tol_dist) THEN
        CALL add_integrals_to_single_row( single_row, vi_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      END IF

      ! Cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      IF (n_cycles > mesh%nV) THEN
        CALL crash('trace_line_Vor - iterative tracer got stuck!')
      END IF

      ! Update vi_hint, for more efficiency
      vi_hint = vi_left

    END DO ! DO WHILE (.NOT.finished)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE trace_line_Vor

  SUBROUTINE trace_line_Vor_start( mesh, p, vi_hint, vi_in, ti_on, ei_on)
    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge ei_on

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    INTEGER,                             INTENT(INOUT) :: vi_hint
    INTEGER,                             INTENT(OUT)   :: vi_in
    INTEGER,                             INTENT(OUT)   :: ti_on
    INTEGER,                             INTENT(OUT)   :: ei_on

    ! Local variables:
    INTEGER                                            :: vti, ti, vei, ei
    REAL(dp), DIMENSION(2)                             :: cc1, cc2

    ! Initialise
    vi_in = 0
    ti_on = 0
    ei_on = 0

    ! Find the vertex whose Voronoi cell contains p
    CALL find_containing_vertex( mesh, p, vi_hint)

    ! Check if p lies on any of the surrounding triangles' circumcentres
    DO vti = 1, mesh%niTri( vi_hint)
      ti = mesh%iTri( vi_hint,vti)
      IF (NORM2( mesh%Tricc( ti,:) - p) < mesh%tol_dist) THEN
        ! p lies on the circumcentre of triangle ti
        ti_on = ti
        RETURN
      END IF
    END DO

    ! Check if p lies on any of the shared Voronoi boundaries
    DO vei = 1, mesh%nC( vi_hint)
      ei = mesh%VE( vi_hint,vei)
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      IF (lies_on_line_segment( cc1, cc2, p, mesh%tol_dist)) THEN
        ! p lies on the shared Voronoi cell boundary represented by edge ei
        ei_on = ei
        RETURN
      END IF
    END DO

    ! If p lies not on the boundary of the Voronoi cell, then it must lie inside of it
    vi_in = vi_hint

  END SUBROUTINE trace_line_Vor_start

  SUBROUTINE trace_line_Vor_vi(  mesh, p, q, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies inside the Voronoi cell of vertex vi_in,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: vi_in
    INTEGER,                             INTENT(INOUT) :: ti_on
    INTEGER,                             INTENT(INOUT) :: ei_on
    INTEGER,                             INTENT(OUT)   :: vi_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: vti, ti, ei, vori, vorj, ci, vj
    REAL(dp), DIMENSION(2)                             :: r, llis, pa, pb
    REAL(dp)                                           :: dx
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION( mesh%nC_mem,2)                :: Vor
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_vi
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_ti
    INTEGER                                            :: nVor

    ! Safety
    IF (vi_in == 0 .OR. ti_on > 0 .OR. ei_on > 0 .OR. (.NOT. is_in_Voronoi_cell( mesh, p, vi_in))) THEN
      CALL crash('trace_line_Vor_vi - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies inside the same Voronoi cell
    IF (is_in_Voronoi_cell( mesh, q, vi_in)) THEN
      ! q lies inside the same Voronoi cell
      p_next    = q
      vi_left   = vi_in
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the boundary of this Voronoi cell
    dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 200
    CALL calc_Voronoi_cell( mesh, vi_in, dx, Vor, Vor_vi, Vor_ti, nVor)
    DO vori = 1, nVor
      vorj = vori + 1
      IF (vorj == nVor + 1) vorj = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori,:)
      pb = Vor( vorj,:)
      IF (NORM2( q - pa) < mesh%tol_dist .OR. lies_on_line_segment( pa, pb, q, mesh%tol_dist)) THEN
        ! q lies on the boundary of the same Voronoi cell
        p_next    = q
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check if [pq] passes through any of the surrounding triangles' circumcentres
    DO vti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,vti)
      r  = mesh%Tricc( ti,:)
      IF (lies_on_line_segment( p, q, r, mesh%tol_dist)) THEN
        ! [pq] passes through this triangle's circumcentre
        p_next    = mesh%Tricc( ti,:)
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check if [pq] passes through any of the shared Voronoi boundaries
    DO vori = 1, nVor
      vorj = vori + 1
      IF (vorj == nVor + 1) vorj = 1
      ! The two endpoints of this section of the Voronoi cell boundary
      pa = Vor( vori,:)
      pb = Vor( vorj,:)
      ! The other vertex sharing this Voronoi cell boundary
      vj = Vor_vi( vori)
      ! The edge representing this shared Voronoi cell boundary
      ei = 0
      DO ci = 1, mesh%nC( vi_in)
        IF (mesh%C( vi_in,ci) == vj) THEN
          ei = mesh%VE( vi_in,ci)
          EXIT
        END IF
      END DO
      ! Safety
      IF (ei == 0) CALL crash('couldnt find edge between vi and vj!')

      CALL segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] passes into the Voronoi cell of vj
        p_next    = llis
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_Vor_vi - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_Vor_vi

  SUBROUTINE trace_line_Vor_ti(  mesh, p, q, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies on the circumcentre of triangle ti_on,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: vi_in
    INTEGER,                             INTENT(INOUT) :: ti_on
    INTEGER,                             INTENT(INOUT) :: ei_on
    INTEGER,                             INTENT(OUT)   :: vi_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vic, vvi, vj, ei, acab, acbc, acca, tj
    REAL(dp), DIMENSION(2)                             :: cc, cc1, cc2, llis
    LOGICAL                                            :: do_cross

    ! Safety
    IF (vi_in > 0 .OR. ti_on == 0 .OR. ei_on > 0 .OR. NORM2( mesh%Tricc( ti_on,:) - p) > mesh%tol_dist) THEN
      CALL crash('trace_line_Vor_ti - coincidence indicators dont make sense!')
    END IF

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_on,1)
    vib = mesh%Tri( ti_on,2)
    vic = mesh%Tri( ti_on,3)

    ! Find the three Voronoi cell boundaries that meet here
    acab = 0
    DO vvi = 1, mesh%nC( via)
      vj = mesh%C(  via,vvi)
      ei = mesh%VE( via,vvi)
      IF (vj == vib) THEN
        acab = ei
        EXIT
      END IF
    END DO
    acbc = 0
    DO vvi = 1, mesh%nC( vib)
      vj = mesh%C(  vib,vvi)
      ei = mesh%VE( vib,vvi)
      IF (vj == vic) THEN
        acbc = ei
        EXIT
      END IF
    END DO
    acca = 0
    DO vvi = 1, mesh%nC( vic)
      vj = mesh%C(  vic,vvi)
      ei = mesh%VE( vic,vvi)
      IF (vj == via) THEN
        acca = ei
        EXIT
      END IF
    END DO

    ! Check if q lies on the Voronoi cell boundary separating via from vib
    CALL find_shared_Voronoi_boundary( mesh, acab, cc1, cc2)
    IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .OR. &
        NORM2( cc1 - q) < mesh%tol_dist .OR. &
        NORM2( cc2 - q) < mesh%tol_dist) THEN
      ! q lies on the Voronoi cell boundary separating via from vib
      IF (mesh%ETri( acab,1) == ti_on) THEN
        vi_left = mesh%EV( acab,2)
      ELSE
        vi_left = mesh%EV( acab,1)
      END IF
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the Voronoi cell boundary separating vib from vic
    CALL find_shared_Voronoi_boundary( mesh, acbc, cc1, cc2)
    IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .OR. &
        NORM2( cc1 - q) < mesh%tol_dist .OR. &
        NORM2( cc2 - q) < mesh%tol_dist) THEN
      ! q lies on the Voronoi cell boundary separating vib from vic
      IF (mesh%ETri( acbc,1) == ti_on) THEN
        vi_left = mesh%EV( acbc,2)
      ELSE
        vi_left = mesh%EV( acbc,1)
      END IF
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the Voronoi cell boundary separating vic from via
    CALL find_shared_Voronoi_boundary( mesh, acca, cc1, cc2)
    IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .OR. &
        NORM2( cc1 - q) < mesh%tol_dist .OR. &
        NORM2( cc2 - q) < mesh%tol_dist) THEN
      ! q lies on the Voronoi cell boundary separating vic from via
      IF (mesh%ETri( acca,1) == ti_on) THEN
        vi_left = mesh%EV( acca,2)
      ELSE
        vi_left = mesh%EV( acca,1)
      END IF
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside any of the three adjacent Voronoi cells
    IF (is_in_Voronoi_cell( mesh, q, via)) THEN
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF
    IF (is_in_Voronoi_cell( mesh, q, vib)) THEN
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF
    IF (is_in_Voronoi_cell( mesh, q, vic)) THEN
      ! q lies inside the Voronoi cell of vic
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vic
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through the circumcentre of any of the three neighbouring triangles
    tj = mesh%TriC( ti_on,1)
    IF (tj > 0) THEN
      cc = mesh%Tricc( tj,:)
      IF (lies_on_line_segment( p, q, cc, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = vic
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    tj = mesh%TriC( ti_on,2)
    IF (tj > 0) THEN
      cc = mesh%Tricc( tj,:)
      IF (lies_on_line_segment( p, q, cc, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = via
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    tj = mesh%TriC( ti_on,3)
    IF (tj > 0) THEN
      cc = mesh%Tricc( tj,:)
      IF (lies_on_line_segment( p, q, cc, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        ei_on     = 0
        vi_left   = vib
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check if [pq] crosses the boundary of the Voronoi cell of via
    DO vvi = 1, mesh%nC( via)
      vj = mesh%C( via,vvi)
      IF (vj == vib .OR. vj == vic) CYCLE
      ei = mesh%VE( via,vvi)
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = via
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check if [pq] crosses the boundary of the Voronoi cell of vib
    DO vvi = 1, mesh%nC( vib)
      vj = mesh%C( vib,vvi)
      IF (vj == via .OR. vj == vic) CYCLE
      ei = mesh%VE( vib,vvi)
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vib
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check if [pq] crosses the boundary of the Voronoi cell of vic
    DO vvi = 1, mesh%nC( vic)
      vj = mesh%C( vic,vvi)
      IF (vj == via .OR. vj == vib) CYCLE
      ei = mesh%VE( vic,vvi)
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vic
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_Vor_ti - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_Vor_ti

  SUBROUTINE trace_line_Vor_ei( mesh, p, q, p_next, vi_in, ti_on, ei_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies on the shared Voronoi boundary represented by edge ei_on,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: vi_in
    INTEGER,                             INTENT(INOUT) :: ti_on
    INTEGER,                             INTENT(INOUT) :: ei_on
    INTEGER,                             INTENT(OUT)   :: vi_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vil, vir, til, tir, vvi, ei, vti, ti
    REAL(dp), DIMENSION(2)                             :: cc1, cc2, ccl, ccr, llis
    LOGICAL                                            :: do_cross

    ! Find the endpoints of this shared Voronoi boundary
    CALL find_shared_Voronoi_boundary( mesh, ei_on, cc1, cc2)

    ! Safety
    IF (vi_in > 0 .OR. ti_on > 0 .OR. ei_on == 0 .OR. (.NOT. lies_on_line_segment( cc1, cc2, p, mesh%tol_dist))) THEN
      CALL crash('trace_line_Vor_ei - coincidence indicators dont make sense!')
    END IF

    ! A bit more detail is needed
    via = mesh%EV(   ei_on,1)
    vib = mesh%EV(   ei_on,2)
    vil = mesh%EV(   ei_on,3)
    vir = mesh%EV(   ei_on,4)
    til = mesh%ETri( ei_on,1)
    tir = mesh%ETri( ei_on,2)

    IF (til == 0) THEN
      ! Apparently ei lies on the domain border and has no triangle on its left-hand side
      ccr = cc1
      ccl = cc2
    ELSEIF (tir == 0) THEN
      ! Apparently ei lies on the domain border and has no triangle on its right-hand side
      ccl = cc1
      ccr = cc2
    ELSE
      ! ei lies in the interior and has triangles on both sides
      ccl = mesh%Tricc( til,:)
      ccr = mesh%Tricc( tir,:)
    END IF

    ! Check if q coincides with ccl
    IF (NORM2( ccl - q) < mesh%tol_dist) THEN
      ! q coincides with ccl
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q coincides with ccr
    IF (NORM2( ccr - q) < mesh%tol_dist) THEN
      ! q coincides with ccr
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through ccl
    IF (lies_on_line_segment( p, q, ccl, mesh%tol_dist)) THEN
      ! [pq] passes through ccl
      p_next    = ccl
      vi_in     = 0
      ti_on     = til
      ei_on     = 0
      vi_left   = via
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through ccr
    IF (lies_on_line_segment( p, q, ccr, mesh%tol_dist)) THEN
      ! [pq] passes through ccr
      p_next    = ccr
      vi_in     = 0
      ti_on     = tir
      ei_on     = 0
      vi_left   = vib
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if q lies inside the Voronoi cell of via
    IF (is_in_Voronoi_cell( mesh, q, via)) THEN
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = via
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the Voronoi cell of vib
    IF (is_in_Voronoi_cell( mesh, q, vib)) THEN
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      ei_on     = 0
      vi_left   = vib
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the circumcentre of any of the triangles surrounding via
    DO vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      IF (NORM2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) THEN
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = via
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check if q lies on the circumcentre of any of the triangles surrounding vib
    DO vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      IF (NORM2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) THEN
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = vib
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check if q lies on boundary of the Voronoi cell of via
    DO vvi = 1, mesh%nC( via)
      ei = mesh%VE( via,vvi)
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) THEN
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = via
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check if q lies on boundary of the Voronoi cell of vib
    DO vvi = 1, mesh%nC( vib)
      ei = mesh%VE( vib,vvi)
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) THEN
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        ei_on     = 0
        vi_left   = vib
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check if pq crosses the circumcentre of any of the triangles surrounding via
    DO vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      cc1 = mesh%Tricc( ti,:)
      IF (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of triangle ti
        p_next    = cc1
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        vi_left   = via
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check if pq crosses the boundary of the Voronoi cell of via
    DO vvi = 1, mesh%nC( via)
      ei = mesh%VE( via,vvi)
      IF (ei == ei_on) CYCLE
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = via
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check if pq crosses the circumcentre of any of the triangles surrounding vib
    DO vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      cc1 = mesh%Tricc( ti,:)
      IF (lies_on_line_segment( p, q, cc1, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of triangle ti
        p_next    = cc1
        vi_in     = 0
        ti_on     = ti
        ei_on     = 0
        vi_left   = vib
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check if pq crosses the boundary of the Voronoi cell of vib
    DO vvi = 1, mesh%nC( vib)
      ei = mesh%VE( vib,vvi)
      IF (ei == ei_on) CYCLE
      CALL find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        ei_on     = ei
        vi_left   = vib
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_Vor_ei - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_Vor_ei

  ! Line tracing algorithm through square grid cells
  SUBROUTINE trace_line_grid( grid, p, q, single_row, count_coincidences)
    ! Trace the line [pq] through the grid and calculate the three
    ! line integrals for the line segments inside the different grid cells.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    LOGICAL,                             INTENT(IN)    :: count_coincidences

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'trace_line_Vor'
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    REAL(dp), DIMENSION(2)                             :: pp,qq
    LOGICAL                                            :: is_valid_line
    LOGICAL                                            :: finished
    INTEGER                                            :: n_cycles
    INTEGER,  DIMENSION(2)                             :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2)                             :: p_next
    INTEGER                                            :: n_left
    LOGICAL                                            :: coincides
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Crop the line [pq] so that it lies within the domain
    xmin = grid%xmin !- grid%dx / 2._dp
    xmax = grid%xmax !+ grid%dx / 2._dp
    ymin = grid%ymin !- grid%dx / 2._dp
    ymax = grid%ymax !+ grid%dx / 2._dp
    CALL crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, grid%tol_dist, pp, qq, is_valid_line)

    IF (.NOT. is_valid_line) THEN
      ! [pq] doesn't pass through the domain anywhere
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
    !    - lies inside grid cell aij_in, ...
    !    - lies on the b-grid point bij_on, or...
    !    - lies on the edge cij_on
    CALL trace_line_grid_start( grid, pp, aij_in, bij_on, cxij_on, cyij_on)

    ! Iteratively trace the line through the mesh
    finished = .FALSE.
    n_cycles = 0
    DO WHILE (.NOT. finished)

      ! Find the point p_next where [pq] crosses into the next Voronoi cell
      IF     (aij_in(  1) > 0 .OR. aij_in(  2) > 0) THEN
        ! p lies inside a-grid cell aij_in
        CALL trace_line_grid_a(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ELSEIF (bij_on(  1) > 0 .OR. bij_on(  2) > 0) THEN
        ! p lies on b-grid point bij_on
        CALL trace_line_grid_b(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ELSEIF (cxij_on( 1) > 0 .OR. cxij_on( 2) > 0) THEN
        ! p lies on cx-grid edge cxij_on
        CALL trace_line_grid_cx( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ELSEIF (cyij_on( 1) > 0 .OR. cyij_on( 2) > 0) THEN
        ! p lies on cy-grid edge cyij_on
        CALL trace_line_grid_cy( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
      ELSE
        CALL crash('found no coincidence indicators!')
      END IF

      ! Calculate the three line integrals
      LI_xdy   = line_integral_xdy(   pp, p_next, grid%tol_dist)
      LI_mxydx = line_integral_mxydx( pp, p_next, grid%tol_dist)
      LI_xydy  = line_integral_xydy(  pp, p_next, grid%tol_dist)

      ! Add them to the results structure
      IF (NORM2( p_next - pp) > grid%tol_dist) THEN
        CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
      END IF

      ! Cycle the pointer
      pp = p_next

      ! Safety
      n_cycles = n_cycles + 1
      IF (n_cycles > grid%n) THEN
        CALL crash('trace_line_grid - iterative tracer got stuck!')
      END IF

    END DO ! DO WHILE (.NOT. finished)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE trace_line_grid

  SUBROUTINE trace_line_grid_start( grid, p,    aij_in, bij_on, cxij_on, cyij_on)
    ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
    !    - lies inside grid cell aij_in, ...
    !    - lies on the b-grid point bij_on, or...
    !    - lies on the edge cij_on

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    INTEGER,  DIMENSION(2),              INTENT(OUT)   :: aij_in, bij_on, cxij_on, cyij_on

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: xl,xu,yl,yu

    ! Initialise
    aij_in  = [0,0]
    bij_on  = [0,0]
    cxij_on = [0,0]
    cyij_on = [0,0]

    ! Find the grid cell containing p
    i = 1 + FLOOR( (p(1) - grid%xmin + grid%dx / 2._dp) / grid%dx)
    j = 1 + FLOOR( (p(2) - grid%ymin + grid%dx / 2._dp) / grid%dx)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! Check IF p lies on either of the four surrounding b-grid points
    IF     (i > 1       .AND. j > 1       .AND. abs( p(1) - xl) < grid%tol_dist .AND. abs( p(2) - yl) < grid%tol_dist) THEN
      ! p coincides with the southwest corner
      bij_on = [i-1,j-1]
      RETURN
    ELSEIF (i > 1       .AND. j < grid%ny .AND. abs( p(1) - xl) < grid%tol_dist .AND. abs( p(2) - yu) < grid%tol_dist) THEN
      ! p coincides with the northwest corner
      bij_on = [i-1,j  ]
      RETURN
    ELSEIF (i < grid%nx .AND. j < 1       .AND. abs( p(1) - xu) < grid%tol_dist .AND. abs( p(2) - yl) < grid%tol_dist) THEN
      ! p coincides with the southeast corner
      bij_on = [i  ,j-1]
      RETURN
    ELSEIF (i < grid%nx .AND. j < grid%ny .AND. abs( p(1) - xu) < grid%tol_dist .AND. abs( p(2) - yu) < grid%tol_dist) THEN
      ! p coincides with the northeast corner
      bij_on = [i  ,j  ]
      RETURN
    END IF

    ! Check IF p lies on any of the four borders
    IF     (i > 1       .AND. abs( p(1) - xl) < grid%tol_dist) THEN
      ! p coincides with the western border
      cxij_on = [i-1,j  ]
      RETURN
    ELSEIF (i < grid%nx .AND. abs( p(1) - xu) < grid%tol_dist) THEN
      ! p coincides with the eastern border
      cxij_on = [i  ,j  ]
      RETURN
    ELSEIF (j > 1       .AND. abs( p(2) - yl) < grid%tol_dist) THEN
      ! p coincides with the southern border
      cyij_on = [i  ,j-1]
      RETURN
    ELSEIF (j < grid%ny .AND. abs( p(2) - yu) < grid%tol_dist) THEN
      ! p coincides with the northern border
      cyij_on = [i  ,j  ]
      RETURN
    END IF

    ! p doesn't lie on the corners or borders, so it must lie inside the grid cell
    aij_in = [i,j]

  END SUBROUTINE trace_line_grid_start

  SUBROUTINE trace_line_grid_a(     grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies inside grid cell aij_in,
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: xl,xu,yl,yu
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF ((aij_in( 1) == 0 .AND. aij_in( 2) == 0) .OR. cxij_on( 1) > 0 .OR. cxij_on( 2) > 0 .OR. &
        bij_on( 1) > 0 .OR. bij_on( 2) > 0 .OR. cyij_on( 1) > 0 .OR. cyij_on( 2) > 0) THEN
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      CALL crash('trace_line_grid_a - coincidence indicators dont make sense!')
    END IF

    i = aij_in( 1)
    j = aij_in( 2)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! More safety
    IF (p(1) < xl .OR. p(1) > xu .OR. p(2) < yl .OR. p(2) > yu) THEN
      CALL crash('trace_line_grid_a - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies inside the same grid cell
    IF (q(1) >= xl - grid%tol_dist .AND. &
        q(1) <= xu + grid%tol_dist .AND. &
        q(2) >= yl - grid%tol_dist .AND. &
        q(2) <= yu + grid%tol_dist) THEN
      ! q lies inside the same grid cell
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if pq passes through any of the four corners
    sw = [xl,yl]
    nw = [xl,yu]
    se = [xu,yl]
    ne = [xu,yu]

    IF (lies_on_line_segment( p, q, sw, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the southwest corner
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    IF (lies_on_line_segment( p, q, nw, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the northwest corner
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    IF (lies_on_line_segment( p, q, se, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the southeast corner
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    IF (lies_on_line_segment( p, q, ne, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the northeast corner
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through any of the four boundaries
    CALL segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    CALL segment_intersection( p, q, se, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    CALL segment_intersection( p, q, sw, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    CALL segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j  ]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_a - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_a

  SUBROUTINE trace_line_grid_b(     grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on b-grid point bij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: x,y,xl,xu,yl,yu
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne,ww,ee,ss,nn
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF (aij_in( 1) > 0 .OR. aij_in( 2) > 0 .OR. cxij_on( 1) > 0 .OR. cxij_on( 2) > 0 .OR. &
        (bij_on( 1) == 0 .AND. bij_on( 2) == 0) .OR. cyij_on( 1) > 0 .OR. cyij_on( 2) > 0) THEN
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      CALL crash('trace_line_grid_b - coincidence indicators dont make sense!')
    END IF

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
    IF (abs( p(1) - x) > grid%tol_dist .OR. abs( p(2) - y) > grid%tol_dist) THEN
      CALL crash('trace_line_grid_b - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies on the cy-grid edge to the west
    IF (q(1) < x + grid%tol_dist .AND. q(1) > xl - grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the cy-grid edge to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the cy-grid edge to the east
    IF (q(1) > x - grid%tol_dist .AND. q(1) < xu + grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the cy-grid edge to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the cx-grid edge to the south
    IF (q(2) < y + grid%tol_dist .AND. q(2) > yl - grid%tol_dist .AND. abs( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the cx-grid edge to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the cx-grid edge to the north
    IF (q(2) > y - grid%tol_dist .AND. q(2) < yu + grid%tol_dist .AND. abs( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the cx-grid edge to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the northwest
    IF (q(1) > xl - grid%tol_dist .AND. q(1) < x  + grid%tol_dist .AND. &
        q(2) > y  - grid%tol_dist .AND. q(2) < yu + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the northwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the northeast
    IF (q(1) > x  - grid%tol_dist .AND. q(1) < xu + grid%tol_dist .AND. &
        q(2) > y  - grid%tol_dist .AND. q(2) < yu + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the northeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the southeast
    IF (q(1) > x  - grid%tol_dist .AND. q(1) < xu + grid%tol_dist .AND. &
        q(2) > yl - grid%tol_dist .AND. q(2) < y  + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the southeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the southwest
    IF (q(1) > xl - grid%tol_dist .AND. q(1) < x  + grid%tol_dist .AND. &
        q(2) > yl - grid%tol_dist .AND. q(2) < y  + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the southwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the west
    IF (lies_on_line_segment( p, q, ww, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the north
    IF (lies_on_line_segment( p, q, nn, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the west
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i  ,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the east
    IF (lies_on_line_segment( p, q, ee, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i+1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the south
    IF (lies_on_line_segment( p, q, ss, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northwest through its western boundary
    CALL segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northwest through its northern boundary
    CALL segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northwest through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northeast through its northern boundary
    CALL segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northeast through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j+1]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northeast through its eastern boundary
    CALL segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southeast through its eastern boundary
    CALL segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southeast through its southern boundary
    CALL segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southeast through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southwest through its southern boundary
    CALL segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southwest through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southwest through its western boundary
    CALL segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_b - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_b

  SUBROUTINE trace_line_grid_cx(    grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on cx-grid edge cxij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: x,yl,yu
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne,ss,nn
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF (aij_in( 1) > 0 .OR. aij_in( 2) > 0 .OR. (cxij_on( 1) == 0 .AND. cxij_on( 2) == 0) .OR. &
        bij_on( 1) > 0 .OR. bij_on( 2) > 0 .OR. cyij_on( 1) > 0 .OR. cyij_on( 2) > 0) THEN
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      CALL crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    END IF

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
    IF (p(2) < yl .OR. p(2) > yu .OR. ABS( p(1) - x) > grid%tol_dist) THEN
      CALL crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    END IF

    ! Check IF q lies on the same cx-grid cell in the southern direction
    IF (q(2) < p(2) .AND. q(2) >= yl - grid%tol_dist .AND. ABS( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the same cx-grid cell in the southern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies on the same cx-grid cell in the northern direction
    IF (q(2) > p(2) .AND. q(2) <= yu + grid%tol_dist .AND. ABS( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the same cx-grid cell in the northern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside the grid cell to the west
    IF (q(2) >= yl - grid%tol_dist .AND. q(2) <= yu + grid%tol_dist .AND. &
        q(1) >= x - grid%dx - grid%tol_dist .AND. q(1) <= x + grid%tol_dist) THEN
      ! q lies inside the grid cell to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside the grid cell to the east
    IF (q(2) >= yl - grid%tol_dist .AND. q(2) <= yu + grid%tol_dist .AND. &
        q(1) <= x + grid%dx + grid%tol_dist .AND. q(1) >= x - grid%tol_dist) THEN
      ! q lies inside the grid cell to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF [pq] passes through the b-grid point to the south
    IF (lies_on_line_segment( p, q, ss, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] passes through the b-grid point to the north
    IF (lies_on_line_segment( p, q, nn, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the north
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through the b-grid point to the northwest
    IF (lies_on_line_segment( p, q, nw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through the b-grid point to the southwest
    IF (lies_on_line_segment( p, q, sw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through the b-grid point to the northeast
    IF (lies_on_line_segment( p, q, ne, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i+1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through the b-grid point to the southeast
    IF (lies_on_line_segment( p, q, se, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i+1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through its southern boundary
    CALL segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the west through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through its western boundary
    CALL segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the west through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through its northern boundary
    CALL segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the west through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through its northern boundary
    CALL segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the east through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through its eastern boundary
    CALL segment_intersection( p, q, ne, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the east through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through its southern boundary
    CALL segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the east through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_cx - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_cx

  SUBROUTINE trace_line_grid_cy(    grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on cy-grid edge cyij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: xl,xu,y
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne,ww,ee
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF (aij_in( 1) > 0 .OR. aij_in( 2) > 0 .OR. cxij_on( 1) > 0 .OR. cxij_on( 2) > 0 .OR. &
        bij_on( 1) > 0 .OR. bij_on( 2) > 0 .OR. (cyij_on( 1) == 0 .AND. cyij_on( 2) == 0)) THEN
      WRITE(0,*) 'aij_in = ', aij_in, ', bij_on = ', bij_on, ', cxij_on = ', cxij_on, ', cyij_on = ', cyij_on
      CALL crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    END IF

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
    IF (p(1) < xl .OR. p(1) > xu .OR. abs( p(2) - y) > grid%tol_dist) THEN
      CALL crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies on the same cy-grid cell in the western direction
    IF (q(1) < p(1) .AND. q(1) >= xl - grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the same cy-grid cell in the western direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the same cy-grid cell in the eastern direction
    IF (q(1) > p(1) .AND. q(1) <= xu + grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the same cy-grid cell in the eastern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the grid cell to the south
    IF (q(1) >= xl - grid%tol_dist .AND. q(1) <= xu + grid%tol_dist .AND. &
        q(2) >= y - grid%dx - grid%tol_dist .AND. q(2) <= y + grid%tol_dist) THEN
      ! q lies inside the grid cell to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the grid cell to the north
    IF (q(1) >= xl - grid%tol_dist .AND. q(1) <= xu + grid%tol_dist .AND. &
        q(2) <= y + grid%dx + grid%tol_dist .AND. q(2) >= y - grid%tol_dist) THEN
      ! q lies inside the grid cell to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the west
    IF (lies_on_line_segment( p, q, ww, grid%tol_dist))  THEN
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the east
    IF (lies_on_line_segment( p, q, ee, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northwest
    IF (lies_on_line_segment( p, q, nw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northeast
    IF (lies_on_line_segment( p, q, ne, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southwest
    IF (lies_on_line_segment( p, q, sw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southeast
    IF (lies_on_line_segment( p, q, se, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through its western boundary
    CALL segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the north through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through its northern boundary
    CALL segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the north through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through its eastern boundary
    CALL segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the north through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through its eastern boundary
    CALL segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the south through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through its southern boundary
    CALL segment_intersection( p, q, se, sw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the south through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through its western boundary
    CALL segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the south through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_cy - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_cy

  ! Add the values for a single row of the three line-integral matrices
  SUBROUTINE add_integrals_to_single_row(  single_row, index_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
    ! Add the values for a single row of the three line-integral matrices

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    INTEGER,                                INTENT(IN)    :: index_left
    REAL(dp),                               INTENT(IN)    :: LI_xdy, LI_mxydx, LI_xydy
    LOGICAL,                                INTENT(IN)    :: coincides, count_coincidences

    ! Local variables:
    LOGICAL                                            :: do_add_integrals, is_listed
    INTEGER                                            :: i, i_add

    ! Check whether we actually need to add the line integrals
    do_add_integrals = .TRUE.
    IF (coincides .AND. (.NOT. count_coincidences)) do_add_integrals = .FALSE.

    ! Check if an entry from this left-hand vertex is already listed
    is_listed = .FALSE.
    i_add     = 0

    DO i = 1, single_row%n
      IF (single_row%index_left( i) == index_left) THEN
        is_listed = .TRUE.
        i_add     = i
        EXIT
      END IF
    END DO
    IF (.NOT. is_listed) THEN
      single_row%n = single_row%n + 1
      i_add = single_row%n
    END IF

    ! Add data
    single_row%index_left( i_add) = index_left
    IF (do_add_integrals) THEN
      single_row%LI_xdy(   i_add) = single_row%LI_xdy(   i_add) + LI_xdy
      single_row%LI_mxydx( i_add) = single_row%LI_mxydx( i_add) + LI_mxydx
      single_row%LI_xydy(  i_add) = single_row%LI_xydy(  i_add) + LI_xydy
    END IF

    ! If necessary, extend memory
    IF (single_row%n > single_row%n_max - 10) CALL extend_single_row_memory( single_row, 100)

  END SUBROUTINE add_integrals_to_single_row

  SUBROUTINE extend_single_row_memory( single_row, n_extra)
    ! Extend memory for a single row of the three line-integral matrices

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    INTEGER,                                INTENT(IN)    :: n_extra

    ! Local variables:
    INTEGER                                            :: n
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: index_left_temp
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: LI_xdy_temp, LI_mxydx_temp, LI_xydy_temp

    n = single_row%n

    ! Allocate temporary memory
    ALLOCATE( index_left_temp( n))
    ALLOCATE( LI_xdy_temp(     n))
    ALLOCATE( LI_mxydx_temp(   n))
    ALLOCATE( LI_xydy_temp(    n))

    ! Copy data to temporary memory
    index_left_temp = single_row%index_left( 1:n)
    LI_xdy_temp     = single_row%LI_xdy(     1:n)
    LI_mxydx_temp   = single_row%LI_mxydx(   1:n)
    LI_xydy_temp    = single_row%LI_xydy(    1:n)

    ! Deallocate memory
    DEALLOCATE( single_row%index_left)
    DEALLOCATE( single_row%LI_xdy    )
    DEALLOCATE( single_row%LI_mxydx  )
    DEALLOCATE( single_row%LI_xydy   )

    ! Allocate new, extended memory
    single_row%n_max = single_row%n_max + n_extra
    ALLOCATE( single_row%index_left( single_row%n_max))
    ALLOCATE( single_row%LI_xdy(     single_row%n_max))
    ALLOCATE( single_row%LI_mxydx(   single_row%n_max))
    ALLOCATE( single_row%LI_xydy(    single_row%n_max))

    ! Copy data back from temporary memory
    single_row%index_left( 1:n) = index_left_temp
    single_row%LI_xdy(     1:n) = LI_xdy_temp
    single_row%LI_mxydx(   1:n) = LI_mxydx_temp
    single_row%LI_xydy(    1:n) = LI_xydy_temp

    ! Deallocate temporary memory
    DEALLOCATE( index_left_temp)
    DEALLOCATE( LI_xdy_temp    )
    DEALLOCATE( LI_mxydx_temp  )
    DEALLOCATE( LI_xydy_temp   )

  END SUBROUTINE extend_single_row_memory

  ! Remapping of a 1-D variable (2nd-order conservative)
  SUBROUTINE remap_cons_2nd_order_1D( z_src, mask_src, d_src, z_dst, mask_dst, d_dst)
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

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_src
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_src
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_dst
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_dst
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_dst

    ! Local variables:
    LOGICAL                                            :: all_are_masked
    INTEGER                                            :: nz_src, nz_dst
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddz_src
    INTEGER                                            :: k_src, k_dst
    REAL(dp)                                           :: zl_src, zu_src, zl_dst, zu_dst, z_lo, z_hi, z, d
    REAL(dp)                                           :: dz_overlap, dz_overlap_tot, d_int, d_int_tot
    REAL(dp)                                           :: dist_to_dst, dist_to_dst_min, max_dist
    INTEGER                                            :: k_src_nearest_to_dst

    ! Initialise
    d_dst = 0._dp

    ! Sizes
    nz_src = SIZE( z_src,1)
    nz_dst = SIZE( z_dst,1)

    ! Maximum distance on combined grids
    max_dist = MAXVAL([ ABS( z_src( nz_src) - z_src( 1)), &
                        ABS( z_dst( nz_dst) - z_dst( 1)), &
                        ABS( z_src( nz_src) - z_dst( 1)), &
                        ABS( z_dst( nz_dst) - z_src( 1))])

    ! Exception for when the entire src field is masked
    all_are_masked = .TRUE.
    DO k = 1, nz_src
      IF (mask_src( k) == 1) all_are_masked = .FALSE.
    END DO
    IF (all_are_masked) RETURN

    ! Exception for when the entire dst field is masked
    all_are_masked = .TRUE.
    DO k = 1, nz_dst
      IF (mask_dst( k) == 1) all_are_masked = .FALSE.
    END DO
    IF (all_are_masked) RETURN

    ! Calculate derivative d_src/dz (one-sided differencing at the boundary, central differencing everywhere else)
    ALLOCATE( ddz_src( nz_src))
    DO k = 2, nz_src-1
      ddz_src( k    ) = (d_src( k+1   ) - d_src( k-1     )) / (z_src( k+1   ) - z_src( k-1     ))
    END DO
    ddz_src(  1     ) = (d_src( 2     ) - d_src( 1       )) / (z_src( 2     ) - z_src( 1       ))
    ddz_src(  nz_src) = (d_src( nz_src) - d_src( nz_src-1)) / (z_src( nz_src) - z_src( nz_src-1))

    ! Perform conservative remapping by finding regions of overlap
    ! between source and destination grid cells

    DO k_dst = 1, nz_dst

      ! Skip masked grid cells
      IF (mask_dst( k_dst) == 0) THEN
        d_dst( k_dst) = 0._dp
        CYCLE
      END IF

      ! Find z range covered by this dst grid cell
      IF (k_dst > 1) THEN
        zl_dst = 0.5_dp * (z_dst( k_dst - 1) + z_dst( k_dst))
      ELSE
        zl_dst = z_dst( 1) - 0.5_dp * (z_dst( 2) - z_dst( 1))
      END IF
      IF (k_dst < nz_dst) THEN
        zu_dst = 0.5_dp * (z_dst( k_dst + 1) + z_dst( k_dst))
      ELSE
        zu_dst = z_dst( nz_dst) + 0.5_dp * (z_dst( nz_dst) - z_dst( nz_dst-1))
      END IF

      ! Find all overlapping src grid cells
      d_int_tot      = 0._dp
      dz_overlap_tot = 0._dp
      DO k_src = 1, nz_src

        ! Skip masked grid cells
        IF (mask_src( k_src) == 0) CYCLE

        ! Find z range covered by this src grid cell
        IF (k_src > 1) THEN
          zl_src = 0.5_dp * (z_src( k_src - 1) + z_src( k_src))
        ELSE
          zl_src = z_src( 1) - 0.5_dp * (z_src( 2) - z_src( 1))
        END IF
        IF (k_src < nz_src) THEN
          zu_src = 0.5_dp * (z_src( k_src + 1) + z_src( k_src))
        ELSE
          zu_src = z_src( nz_src) + 0.5_dp * (z_src( nz_src) - z_src( nz_src-1))
        END IF

        ! Find region of overlap
        z_lo = MAX( zl_src, zl_dst)
        z_hi = MIN( zu_src, zu_dst)
        dz_overlap = MAX( 0._dp, z_hi - z_lo)

        ! Calculate integral over region of overlap and add to sum
        IF (dz_overlap > 0._dp) THEN
          z = 0.5_dp * (z_lo + z_hi)
          d = d_src( k_src) + ddz_src( k_src) * (z - z_src( k_src))
          d_int = d * dz_overlap

          d_int_tot      = d_int_tot      + d_int
          dz_overlap_tot = dz_overlap_tot + dz_overlap
        END IF

      END DO ! DO k_src = 1, nz_src

      IF (dz_overlap_tot > 0._dp) THEN
        ! Calculate dst value
        d_dst( k_dst) = d_int_tot / dz_overlap_tot
      ELSE
        ! Exception for when no overlapping src grid cells were found; use nearest-neighbour extrapolation

        k_src_nearest_to_dst = 0._dp
        dist_to_dst_min      = max_dist
        DO k_src = 1, nz_src
          IF (mask_src( k_src) == 1) THEN
            dist_to_dst = ABS( z_src( k_src) - z_dst( k_dst))
            IF (dist_to_dst < dist_to_dst_min) THEN
              dist_to_dst_min      = dist_to_dst
              k_src_nearest_to_dst = k_src
            END IF
          END IF
        END DO

        ! Safety
        IF (k_src_nearest_to_dst == 0) THEN
          WRITE(0,*) '  remap_cons_2nd_order_1D - ERROR: couldnt find nearest neighbour on source grid!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF

        d_dst( k_dst) = d_src( k_src_nearest_to_dst)

      END IF ! IF (dz_overlap_tot > 0._dp) THEN

    END DO ! DO k_dst = 1, nz_dst

    ! Clean up after yourself
    DEALLOCATE( ddz_src)

  END SUBROUTINE remap_cons_2nd_order_1D

END MODULE mesh_remapping
