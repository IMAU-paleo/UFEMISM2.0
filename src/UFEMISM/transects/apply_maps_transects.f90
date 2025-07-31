module apply_maps_transects

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use transect_types, only: type_transect
  use remapping_types, only: type_map
  use petsc_basic, only: multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D

  implicit none

  private

  public :: &
    apply_map_mesh_vertices_to_transect_2D, apply_map_mesh_vertices_to_transect_3D, &
    apply_map_mesh_triangles_to_transect_2D, apply_map_mesh_triangles_to_transect_3D

contains

  !> Map a 2-D data field from the vertices of a mesh to a transect
  subroutine apply_map_mesh_vertices_to_transect_2D( mesh, transect, map, d_mesh_partial, d_transect_partial)

    ! In/output variables
    type(type_mesh),        intent(in)  :: mesh
    type(type_transect),    intent(in)  :: transect
    type(type_map),         intent(in)  :: map
    real(dp), dimension(:), intent(in)  :: d_mesh_partial
    real(dp), dimension(:), intent(out) :: d_transect_partial

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_map_mesh_vertices_to_transect_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_vertices_to_transect_2D

  !> Map a 3-D data field from the vertices of a mesh to a transect
  subroutine apply_map_mesh_vertices_to_transect_3D( mesh, transect, map, d_mesh_partial, d_transect_partial)

    ! In/output variables
    type(type_mesh),          intent(in)  :: mesh
    type(type_transect),      intent(in)  :: transect
    type(type_map),           intent(in)  :: map
    real(dp), dimension(:,:), intent(in)  :: d_mesh_partial
    real(dp), dimension(:,:), intent(out) :: d_transect_partial

    ! Local variables:
    character(len=1024), parameter                  :: routine_name = 'apply_map_mesh_vertices_to_transect_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_vertices_to_transect_3D

  !> Map a 2-D data field from the vertices of a mesh to a transect
  subroutine apply_map_mesh_triangles_to_transect_2D( mesh, transect, map, d_mesh_partial, d_transect_partial)

    ! In/output variables
    type(type_mesh),        intent(in)  :: mesh
    type(type_transect),    intent(in)  :: transect
    type(type_map),         intent(in)  :: map
    real(dp), dimension(:), intent(in)  :: d_mesh_partial
    real(dp), dimension(:), intent(out) :: d_transect_partial

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_map_mesh_triangles_to_transect_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_triangles_to_transect_2D

  !> Map a 2-D data field from the vertices of a mesh to a transect
  subroutine apply_map_mesh_triangles_to_transect_3D( mesh, transect, map, d_mesh_partial, d_transect_partial)

    ! In/output variables
    type(type_mesh),          intent(in)  :: mesh
    type(type_transect),      intent(in)  :: transect
    type(type_map),           intent(in)  :: map
    real(dp), dimension(:,:), intent(in)  :: d_mesh_partial
    real(dp), dimension(:,:), intent(out) :: d_transect_partial

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_map_mesh_triangles_to_transect_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_triangles_to_transect_3D

end module apply_maps_transects
