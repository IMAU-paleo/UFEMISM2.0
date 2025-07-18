module mesh_data_smoothing

  !> Use 2nd-order conservative remapping to map the data from the mesh
  !> to the square grid. Apply the smoothing on the gridded data, then map
  !> it back to the mesh. The numerical diffusion arising from the two mapping
  !> operations is not a problem since we're smoothing the data anyway.

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, map_from_mesh_vertices_to_xy_grid_3D, &
    map_from_xy_grid_to_mesh_2D, map_from_xy_grid_to_mesh_3D
  use smooth_gridded_data, only: smooth_Gaussian_grid

  implicit none

  private

  public :: smooth_Gaussian

  interface smooth_Gaussian
    procedure smooth_Gaussian_2D
    procedure smooth_Gaussian_3D
  end interface smooth_Gaussian

contains

  subroutine smooth_Gaussian_2D( mesh, grid, output_dir, d_mesh_partial, r)

    ! In/output variables:
    type(type_mesh),        intent(in)    :: mesh
    type(type_grid),        intent(in)    :: grid
    character(len=*),       intent(in   ) :: output_dir
    real(dp), dimension(:), intent(inout) :: d_mesh_partial
    real(dp),               intent(in)    :: r

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'smooth_Gaussian_2D'
    real(dp), dimension(:), allocatable :: d_grid_vec_partial

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    allocate( d_grid_vec_partial( grid%n_loc))

    ! Map data to the grid
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, output_dir, d_mesh_partial, d_grid_vec_partial)

    ! Apply smoothing on the gridded data
    call smooth_Gaussian_grid( grid, d_grid_vec_partial, r)

    ! Map data back to the mesh
    call map_from_xy_grid_to_mesh_2D( grid, mesh, output_dir, d_grid_vec_partial, d_mesh_partial)

    ! Clean up after yourself
    deallocate( d_grid_vec_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_Gaussian_2D

  subroutine smooth_Gaussian_3D( mesh, grid, output_dir, d_mesh_partial, r)

    ! In/output variables:
    type(type_mesh),          intent(in)    :: mesh
    type(type_grid),          intent(in)    :: grid
    character(len=*),         intent(in   ) :: output_dir
    real(dp), dimension(:,:), intent(inout) :: d_mesh_partial
    real(dp),                 intent(in)    :: r

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'smooth_Gaussian_3D'
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    allocate( d_grid_vec_partial( grid%n_loc, size( d_mesh_partial,2)))

    ! Map data to the grid
    call map_from_mesh_vertices_to_xy_grid_3D( mesh, grid, output_dir, d_mesh_partial, d_grid_vec_partial)

    ! Apply smoothing on the gridded data
    call smooth_Gaussian_grid( grid, d_grid_vec_partial, r)

    ! Map data back to the mesh
    call map_from_xy_grid_to_mesh_3D( grid, mesh, output_dir, d_grid_vec_partial, d_mesh_partial)

    ! Clean up after yourself
    deallocate( d_grid_vec_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_Gaussian_3D

end module mesh_data_smoothing
