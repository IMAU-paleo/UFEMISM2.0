module mesh_creation_main

  ! Routines used to create a mesh.

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use reduce_ice_geometry, only: reduce_gridded_ice_geometry, reduce_meshed_ice_geometry
  use mesh_creation_from_reduced_geometry, only: create_mesh_from_reduced_geometry

  implicit none

  private

  public :: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success

contains

  subroutine create_mesh_from_gridded_geometry( region_name, name, &
    grid, Hi, Hb, Hs, SL, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    !, Create a mesh from an ice geometry defined on a grid

    ! In/output variables:
    character(len=3),           intent(in   ) :: region_name
    character(len=256),         intent(in   ) :: name
    type(type_grid),            intent(in   ) :: grid
    real(dp), dimension(:    ), intent(in   ) :: Hi, Hb, Hs, SL
    real(dp),                   intent(in   ) :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in   ) :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(  out) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'create_mesh_from_gridded_geometry'
    real(dp), dimension(:,:  ), allocatable :: poly_mult_sheet
    real(dp), dimension(:,:  ), allocatable :: poly_mult_shelf
    real(dp), dimension(:,:  ), allocatable :: p_line_grounding_line
    real(dp), dimension(:,:  ), allocatable :: p_line_calving_front
    real(dp), dimension(:,:  ), allocatable :: p_line_ice_front
    real(dp), dimension(:,:  ), allocatable :: p_line_coastline

    ! Add routine to path
    call init_routine( routine_name)

    ! Reduce the gridded ice geometry to lines and polygons
    call reduce_gridded_ice_geometry( grid, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! Create a mesh from the reduced ice geometry
    call create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_gridded_geometry

  subroutine create_mesh_from_meshed_geometry( region_name, name, mesh_src, Hi, Hb, Hs, SL, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    !< Create a mesh from an ice geometry defined on a mesh

    ! In/output variables:
    character(len=3),           intent(in   ) :: region_name
    character(len=256),         intent(in   ) :: name
    type(type_mesh),            intent(in   ) :: mesh_src
    real(dp), dimension(:    ), intent(in   ) :: Hi, Hb, Hs, SL
    real(dp),                   intent(in   ) :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in   ) :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(  out) :: mesh

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'create_mesh_from_meshed_geometry'
    real(dp), dimension(:,:  ), allocatable :: poly_mult_sheet
    real(dp), dimension(:,:  ), allocatable :: poly_mult_shelf
    real(dp), dimension(:,:  ), allocatable :: p_line_grounding_line
    real(dp), dimension(:,:  ), allocatable :: p_line_calving_front
    real(dp), dimension(:,:  ), allocatable :: p_line_ice_front
    real(dp), dimension(:,:  ), allocatable :: p_line_coastline

    ! Add routine to path
    call init_routine( routine_name)

    ! Reduce the meshed ice geometry to lines and polygons
    call reduce_meshed_ice_geometry( mesh_src, Hi, Hb, Hs, SL, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline)

    ! Create a mesh from the reduced ice geometry
    call create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
      p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_meshed_geometry

  subroutine write_mesh_success( mesh)
    !< Write the mesh creation success message to the terminal

    use control_resources_and_error_messaging, only: insert_val_into_string_int, insert_val_into_string_dp

    ! In/output variables:
    type(type_mesh), intent(in) :: mesh

    ! Local variables:
    character(len=256), parameter :: routine_name = 'write_mesh_success'
    character(len=256)            :: str

    ! Add routine to path
    call init_routine( routine_name)

    str = '     Set up ' // colour_string( TRIM( mesh%name),'light blue') // ' with {int_01} vertices and {int_02} triangles' // &
      ', with a resolution of {dp_01} to {dp_02} m'
    call insert_val_into_string_int( str, '{int_01}', mesh%nV)
    call insert_val_into_string_int( str, '{int_02}', mesh%nTri)
    call insert_val_into_string_dp(  str, '{dp_01}', MINVAL( mesh%R))
    call insert_val_into_string_dp(  str, '{dp_02}', MAXVAL( mesh%R))

    if (par%primary) WRITE(0,'(A)') trim( str)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_mesh_success

end module mesh_creation_main
