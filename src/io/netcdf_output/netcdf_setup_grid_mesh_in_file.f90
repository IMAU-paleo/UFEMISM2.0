module netcdf_setup_grid_mesh_in_file
  !< Set up a grid or a mesh in a NetCDF file

  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: gather_CSR_dist_to_primary
  use netcdf_basic
  use netcdf_add_field_mesh
  use netcdf, only: NF90_DOUBLE, NF90_INT, NF90_DEF_GRP

  implicit none

  private

  public :: setup_xy_grid_in_netcdf_file, setup_mesh_in_netcdf_file, write_matrix_operators_to_netcdf_file
  public :: save_xy_grid_as_netcdf, save_mesh_as_netcdf

contains

  subroutine save_xy_grid_as_netcdf( filename, grid)

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    type(type_grid),  intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'save_xy_grid_as_netcdf'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_xy_grid_as_netcdf

  subroutine save_mesh_as_netcdf( filename, mesh)

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    type(type_mesh),  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'save_mesh_as_netcdf'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_mesh_as_netcdf

  subroutine setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    !< Set up a regular x/y-grid in an existing NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    type(type_grid),  intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_xy_grid_in_netcdf_file'
    integer                        :: id_dim_x
    integer                        :: id_dim_y
    integer                        :: id_var_x
    integer                        :: id_var_y
    integer                        :: id_var_lon
    integer                        :: id_var_lat

    ! Add routine to path
    call init_routine( routine_name)

    ! Create x/y dimensions
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_x), grid%nx, id_dim_x)
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_y), grid%ny, id_dim_y)

    ! Create and write x/y variables

    ! x
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_x), NF90_DOUBLE, (/ id_dim_x /), id_var_x)
    call add_attribute_char( filename, ncid, id_var_x, 'long_name', 'x-coordinate')
    call add_attribute_char( filename, ncid, id_var_x, 'units'    , 'm'           )
    call write_var_primary( filename, ncid, id_var_x, grid%x)

    ! y
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_y), NF90_DOUBLE, (/ id_dim_y /), id_var_y)
    call add_attribute_char( filename, ncid, id_var_y, 'long_name', 'y-coordinate')
    call add_attribute_char( filename, ncid, id_var_y, 'units'    , 'm'           )
    call write_var_primary( filename, ncid, id_var_y, grid%y)

    ! lon/lat-coordinates
    if (allocated( grid%lon) .or. allocated( grid%lat)) then

      ! Safety
      if (.not. allocated( grid%lon)) call crash('grid has lat but no lon coordinates!')
      if (.not. allocated( grid%lat)) call crash('grid has lon but no lat coordinates!')

      ! lon
      call create_variable( filename, ncid, get_first_option_from_list( field_name_options_lon), NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var_lon)
      call add_attribute_char( filename, ncid, id_var_lon, 'long_name', 'Longitude')
      call add_attribute_char( filename, ncid, id_var_lon, 'units'    , 'degrees east')
      call write_var_primary( filename, ncid, id_var_lon, grid%lon)

      ! lat
      call create_variable( filename, ncid, get_first_option_from_list( field_name_options_lat), NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var_lat)
      call add_attribute_char( filename, ncid, id_var_lat, 'long_name', 'Latitude')
      call add_attribute_char( filename, ncid, id_var_lat, 'units'    , 'degrees north')
      call write_var_primary( filename, ncid, id_var_lat, grid%lat)

    end if ! if (allocated( grid%lon)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_xy_grid_in_netcdf_file

  subroutine setup_mesh_in_netcdf_file( filename, ncid, mesh)
    !< Set up a mesh in an existing NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    type(type_mesh),  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_mesh_in_netcdf_file'

    integer :: id_dim_vi
    integer :: id_dim_ti
    integer :: id_dim_ci
    integer :: id_dim_ei
    integer :: id_dim_vori
    integer :: id_dim_two
    integer :: id_dim_three
    integer :: id_dim_four

    integer :: id_var_xmin
    integer :: id_var_xmax
    integer :: id_var_ymin
    integer :: id_var_ymax
    integer :: id_var_tol_dist
    integer :: id_var_lambda_M
    integer :: id_var_phi_M
    integer :: id_var_beta_stereo

    integer :: id_var_V
    integer :: id_var_nC
    integer :: id_var_C
    integer :: id_var_niTri
    integer :: id_var_iTri
    integer :: id_var_VBI

    integer :: id_var_Tri
    integer :: id_var_Tricc
    integer :: id_var_TriC
    integer :: id_var_TriBI

    integer :: id_var_E
    integer :: id_var_VE
    integer :: id_var_EV
    integer :: id_var_ETri
    integer :: id_var_TriE
    integer :: id_var_EBI
    integer :: id_var_EA

    integer :: id_var_vi2vori
    integer :: id_var_ti2vori
    integer :: id_var_ei2vori
    integer :: id_var_vori2vi
    integer :: id_var_vori2ti
    integer :: id_var_vori2ei
    integer :: id_var_Vor
    integer :: id_var_VornC
    integer :: id_var_VorC
    integer :: id_var_nVVor
    integer :: id_var_VVor

    integer :: id_var_TriGC
    integer :: id_var_TriA
    integer :: id_var_R
    integer :: id_var_A
    integer :: id_var_lon
    integer :: id_var_lat

    ! Add routine to path
    call init_routine( routine_name)

    ! Create mesh dimensions
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nV    ), mesh%nV    , id_dim_vi   )
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nTri  ), mesh%nTri  , id_dim_ti   )
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nC_mem), mesh%nC_mem, id_dim_ci   )
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nE    ), mesh%nE    , id_dim_ei   )
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nVor  ), mesh%nVor  , id_dim_vori )
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_two   ), 2          , id_dim_two  )
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_three ), 3          , id_dim_three)
    call create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_four  ), 4          , id_dim_four )

    ! == Create mesh variables - metadata
    ! ===================================

    ! xmin
    call create_scalar_variable( filename, ncid, 'xmin', NF90_DOUBLE, id_var_xmin)
    call add_attribute_char( filename, ncid, id_var_xmin, 'long_name'  , 'Location of western domain border')
    call add_attribute_char( filename, ncid, id_var_xmin, 'units', 'm')

    ! xmax
    call create_scalar_variable( filename, ncid, 'xmax', NF90_DOUBLE, id_var_xmax)
    call add_attribute_char( filename, ncid, id_var_xmax, 'long_name'  , 'Location of eastern domain border')
    call add_attribute_char( filename, ncid, id_var_xmax, 'units', 'm')

    ! ymin
    call create_scalar_variable( filename, ncid, 'ymin', NF90_DOUBLE, id_var_ymin)
    call add_attribute_char( filename, ncid, id_var_ymin, 'long_name'  , 'Location of southern domain border')
    call add_attribute_char( filename, ncid, id_var_ymin, 'units', 'm')

    ! ymax
    call create_scalar_variable( filename, ncid, 'ymax', NF90_DOUBLE, id_var_ymax)
    call add_attribute_char( filename, ncid, id_var_ymax, 'long_name'  , 'Location of northern domain border')
    call add_attribute_char( filename, ncid, id_var_ymax, 'units', 'm')

    ! tol_dist
    call create_scalar_variable( filename, ncid, 'tol_dist', NF90_DOUBLE, id_var_tol_dist)
    call add_attribute_char( filename, ncid, id_var_tol_dist, 'long_name'  , 'Spatial tolerance (points within this distance are assumed identical)')
    call add_attribute_char( filename, ncid, id_var_tol_dist, 'units', 'm')

    ! lambda_M
    call create_scalar_variable( filename, ncid, 'lambda_M', NF90_DOUBLE, id_var_lambda_M)
    call add_attribute_char( filename, ncid, id_var_lambda_M, 'long_name'  , 'Longitude of the pole of the oblique stereographic projection')
    call add_attribute_char( filename, ncid, id_var_lambda_M, 'units', 'degrees east')

    ! phi_M
    call create_scalar_variable( filename, ncid, 'phi_M', NF90_DOUBLE, id_var_phi_M)
    call add_attribute_char( filename, ncid, id_var_phi_M, 'long_name'  , 'Latitude of the pole of the oblique stereographic projection')
    call add_attribute_char( filename, ncid, id_var_phi_M, 'units', 'degrees north')

    ! beta_stereo
    call create_scalar_variable( filename, ncid, 'beta_stereo', NF90_DOUBLE, id_var_beta_stereo)
    call add_attribute_char( filename, ncid, id_var_beta_stereo, 'long_name'  , 'Standard parallel of the oblique stereographic projection')
    call add_attribute_char( filename, ncid, id_var_beta_stereo, 'units', 'degrees')

    ! == Create mesh variables - vertex data
    ! ======================================

    ! V
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_V             ), NF90_DOUBLE, (/ id_dim_vi, id_dim_two   /), id_var_V             )
    call add_attribute_char( filename, ncid, id_var_V, 'long_name'  , 'Vertex coordinates'         )
    call add_attribute_char( filename, ncid, id_var_V, 'units'      , 'm'                          )
    ! nC
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_nC            ), NF90_INT   , (/ id_dim_vi               /), id_var_nC            )
    call add_attribute_char( filename, ncid, id_var_nC, 'long_name'  , 'Number of vertex-vertex connections')
    ! C
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_C             ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_C             )
    call add_attribute_char( filename, ncid, id_var_C, 'long_name'  , 'Vertex-vertex connections')
    call add_attribute_char( filename, ncid, id_var_C, 'orientation', 'counter-clockwise'          )
    ! niTri
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_niTri         ), NF90_INT   , (/ id_dim_vi               /), id_var_niTri         )
    call add_attribute_char( filename, ncid, id_var_niTri, 'long_name'  , 'Number of vertex-triangle connections')
    ! iTri
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_iTri          ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_iTri          )
    call add_attribute_char( filename, ncid, id_var_iTri, 'long_name'  , 'Vertex-triangle connections')
    call add_attribute_char( filename, ncid, id_var_iTri, 'orientation', 'counter-clockwise'          )
    ! VBI
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VBI           ), NF90_INT   , (/ id_dim_vi               /), id_var_VBI           )
    call add_attribute_char( filename, ncid, id_var_VBI, 'long_name'  , 'Vertex boundary index')
    call add_attribute_char( filename, ncid, id_var_VBI, 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')

    ! == Create mesh variables - triangle data
    ! ========================================

    ! Tri
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_Tri           ), NF90_INT   , (/ id_dim_ti, id_dim_three /), id_var_Tri           )
    call add_attribute_char( filename, ncid, id_var_Tri, 'long_name'  , 'Vertex indices per triangle')
    call add_attribute_char( filename, ncid, id_var_Tri, 'orientation', 'counter-clockwise'          )
    ! Tricc
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_Tricc         ), NF90_DOUBLE, (/ id_dim_ti, id_dim_two   /), id_var_Tricc         )
    call add_attribute_char( filename, ncid, id_var_Tricc, 'long_name'  , 'Triangle circumcentre coordinates')
    call add_attribute_char( filename, ncid, id_var_Tricc, 'units'      , 'm'                          )
    ! TriC
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_TriC          ), NF90_INT   , (/ id_dim_ti, id_dim_three /), id_var_TriC          )
    call add_attribute_char( filename, ncid, id_var_TriC, 'long_name'  , 'Triangle-triangle connections')
    call add_attribute_char( filename, ncid, id_var_TriC, 'orientation', 'counter-clockwise, opposite from constituent vertices (i.e. first entry is opposite from first vertex)')
    ! TriBI
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_TriBI         ), NF90_INT   , (/ id_dim_ti               /), id_var_TriBI)
    call add_attribute_char( filename, ncid, id_var_TriBI, 'long_name'  , 'Triangle boundary index')
    call add_attribute_char( filename, ncid, id_var_TriBI, 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')

    ! == Create mesh variables - edge data
    ! ====================================

    ! E
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_E             ), NF90_DOUBLE, (/ id_dim_ei, id_dim_two   /), id_var_E             )
    call add_attribute_char( filename, ncid, id_var_E, 'long_name'  , 'Edge midpoint coordinates')
    call add_attribute_char( filename, ncid, id_var_E, 'units'      , 'm'                           )
    ! VE
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VE            ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_VE            )
    call add_attribute_char( filename, ncid, id_var_VE, 'long_name'  , 'Vertex-to-edge connectivity list')
    call add_attribute_char( filename, ncid, id_var_VE, 'orientation', 'same as vertex-vertex connectivity list')
    ! EV
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_EV            ), NF90_INT   , (/ id_dim_ei,  id_dim_four /), id_var_EV            )
    call add_attribute_char( filename, ncid, id_var_EV, 'long_name'  , 'Edge-to-vertex connectivity list')
    call add_attribute_char( filename, ncid, id_var_EV, 'orientation', 'vi,vj,vl,vr (start,end,left,right)')
    ! ETri
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_ETri          ), NF90_INT   , (/ id_dim_ei,  id_dim_two  /), id_var_ETri          )
    call add_attribute_char( filename, ncid, id_var_ETri, 'long_name'  , 'Edge-to-triangle connectivity list')
    call add_attribute_char( filename, ncid, id_var_ETri, 'orientation', 'tl,tr (left,right)')
    ! TriE
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_TriE          ), NF90_INT   , (/ id_dim_ti, id_dim_three /), id_var_TriE          )
    call add_attribute_char( filename, ncid, id_var_TriE, 'long_name'  , 'Triangle-to-edge connectivity list')
    call add_attribute_char( filename, ncid, id_var_TriE, 'orientation', 'same as triangle-triangle connectivity list')
    ! EBI
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_EBI           ), NF90_INT   , (/ id_dim_ei               /), id_var_EBI           )
    call add_attribute_char( filename, ncid, id_var_EBI, 'long_name'  , 'Edge boundary index')
    call add_attribute_char( filename, ncid, id_var_EBI, 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')
    ! EA
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_EA), NF90_DOUBLE, (/ id_dim_ei /), id_var_EA)
    call add_attribute_char( filename, ncid, id_var_EA, 'long_name'  , 'Edge area')
    call add_attribute_char( filename, ncid, id_var_EA, 'units', 'm^2')

    ! == Create mesh variables - Voronoi mesh data
    ! ============================================

    ! vi2vori
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vi2vori), NF90_INT, (/ id_dim_vi /), id_var_vi2vori)
    call add_attribute_char( filename, ncid, id_var_vi2vori, 'long_name' , 'Translation table from regular vertices to Voronoi vertices')
    ! ti2vori
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_ti2vori), NF90_INT, (/ id_dim_ti /), id_var_ti2vori)
    call add_attribute_char( filename, ncid, id_var_ti2vori, 'long_name' , 'Translation table from triangles to Voronoi vertices')
    ! ei2vori
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_ei2vori), NF90_INT, (/ id_dim_ei /), id_var_ei2vori)
    call add_attribute_char( filename, ncid, id_var_ei2vori, 'long_name' , 'Translation table from edges to Voronoi vertices')
    ! vori2vi
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vori2vi), NF90_INT, (/ id_dim_vori /), id_var_vori2vi)
    call add_attribute_char( filename, ncid, id_var_vori2vi, 'long_name' , 'Translation table from Voronoi vertices to regular vertices')
    ! vori2ti
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vori2ti), NF90_INT, (/ id_dim_vori /), id_var_vori2ti)
    call add_attribute_char( filename, ncid, id_var_vori2ti, 'long_name' , 'Translation table from Voronoi vertices to triangles')
    ! vori2ei
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vori2ei), NF90_INT, (/ id_dim_vori /), id_var_vori2ei)
    call add_attribute_char( filename, ncid, id_var_vori2ei, 'long_name' , 'Translation table from Voronoi vertices to edges')
    ! Vor
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_Vor), NF90_DOUBLE, (/ id_dim_vori, id_dim_two /), id_var_Vor)
    call add_attribute_char( filename, ncid, id_var_Vor, 'long_name' , 'Voronoi vertex coordinates')
    call add_attribute_char( filename, ncid, id_var_Vor, 'units', 'm')
    ! VornC
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VornC), NF90_INT, (/ id_dim_vori /), id_var_VornC)
    call add_attribute_char( filename, ncid, id_var_VornC, 'long_name' , 'Number of Voronoi vertex connections')
    ! VorC
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VorC), NF90_INT, (/ id_dim_vori, id_dim_three /), id_var_VorC)
    call add_attribute_char( filename, ncid, id_var_VorC, 'long_name' , 'Indices of Voronoi vertex connections')
    ! nVVor
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_nVVor), NF90_INT, (/ id_dim_vi /), id_var_nVVor)
    call add_attribute_char( filename, ncid, id_var_nVVor, 'long_name' , 'Number of Voronoi vertices spanning each Voronoi cell')
    ! VVor
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VVor), NF90_INT, (/ id_dim_vi, id_dim_ci /), id_var_VVor)
    call add_attribute_char( filename, ncid, id_var_VVor, 'long_name' , 'Indices of Voronoi vertices spanning each Voronoi cell')

    ! == Create mesh variables - secondary geometry data
    ! ==================================================

    ! TriGC
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_TriGC         ), NF90_DOUBLE, (/ id_dim_ti, id_dim_two   /), id_var_TriGC         )
    call add_attribute_char( filename, ncid, id_var_TriGC, 'long_name'  , 'Triangle geometric centre coordinates')
    call add_attribute_char( filename, ncid, id_var_TriGC, 'units'      , 'm')
    ! TriA
    call add_field_mesh_dp_2D_b_notime( filename, ncid, get_first_option_from_list( field_name_options_TriA), long_name = 'Triangle area', units = 'm^2')
    call inquire_var(                   filename, ncid, get_first_option_from_list( field_name_options_TriA), id_var_TriA)
    ! R
    call add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_R), long_name = 'Resolution', units = 'm')
    call inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_R), id_var_R)
    ! A
    call add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_A), long_name = 'Voronoi cell area', units = 'm^2')
    call inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_A), id_var_A)
    ! lon
    call add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_lon), long_name = 'Longitude', units = 'degrees east')
    call inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_lon), id_var_lon)
    ! lat
    call add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_lat), long_name = 'Latitude' , units = 'degrees north')
    call inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_lat), id_var_lat)

    ! == Write mesh data to file
    ! ==========================

    ! Metadata
    call write_var_primary(  filename, ncid, id_var_xmin       , mesh%xmin       )
    call write_var_primary(  filename, ncid, id_var_xmax       , mesh%xmax       )
    call write_var_primary(  filename, ncid, id_var_ymin       , mesh%ymin       )
    call write_var_primary(  filename, ncid, id_var_ymax       , mesh%ymax       )
    call write_var_primary(  filename, ncid, id_var_tol_dist   , mesh%tol_dist   )
    call write_var_primary(  filename, ncid, id_var_lambda_M   , mesh%lambda_M   )
    call write_var_primary(  filename, ncid, id_var_phi_M      , mesh%phi_M      )
    call write_var_primary(  filename, ncid, id_var_beta_stereo, mesh%beta_stereo)

    ! Vertex data
    call write_var_primary( filename, ncid, id_var_V    , mesh%V    )
    call write_var_primary( filename, ncid, id_var_nC   , mesh%nC   )
    call write_var_primary( filename, ncid, id_var_C    , mesh%C    )
    call write_var_primary( filename, ncid, id_var_niTri, mesh%niTri)
    call write_var_primary( filename, ncid, id_var_iTri , mesh%iTri )
    call write_var_primary( filename, ncid, id_var_VBI  , mesh%VBI  )

    ! Triangle data
    call write_var_primary( filename, ncid, id_var_Tri  , mesh%Tri  )
    call write_var_primary( filename, ncid, id_var_Tricc, mesh%Tricc)
    call write_var_primary( filename, ncid, id_var_TriC , mesh%TriC )
    call write_var_primary( filename, ncid, id_var_TriBI, mesh%TriBI)

    ! Edge data
    call write_var_primary( filename, ncid, id_var_E   , mesh%E   )
    call write_var_primary( filename, ncid, id_var_VE  , mesh%VE  )
    call write_var_primary( filename, ncid, id_var_EV  , mesh%EV  )
    call write_var_primary( filename, ncid, id_var_ETri, mesh%ETri)
    call write_var_primary( filename, ncid, id_var_TriE, mesh%TriE)
    call write_var_primary( filename, ncid, id_var_EBI , mesh%EBI )
    call write_var_primary( filename, ncid, id_var_EA  , mesh%EA  )

    ! Voronoi mesh data
    call write_var_primary( filename, ncid, id_var_vi2vori, mesh%vi2vori)
    call write_var_primary( filename, ncid, id_var_ti2vori, mesh%ti2vori)
    call write_var_primary( filename, ncid, id_var_ei2vori, mesh%ei2vori)
    call write_var_primary( filename, ncid, id_var_vori2vi, mesh%vori2vi)
    call write_var_primary( filename, ncid, id_var_vori2ti, mesh%vori2ti)
    call write_var_primary( filename, ncid, id_var_vori2ei, mesh%vori2ei)
    call write_var_primary( filename, ncid, id_var_Vor    , mesh%Vor    )
    call write_var_primary( filename, ncid, id_var_VornC  , mesh%VornC  )
    call write_var_primary( filename, ncid, id_var_VorC   , mesh%VorC   )
    call write_var_primary( filename, ncid, id_var_nVVor  , mesh%nVVor  )
    call write_var_primary( filename, ncid, id_var_VVor   , mesh%VVor   )

    ! Secondary geometry data
    call write_var_primary( filename, ncid, id_var_TriGC, mesh%TriGC)
    call write_var_primary( filename, ncid, id_var_TriA , mesh%TriA )
    call write_var_primary( filename, ncid, id_var_R    , mesh%R    )
    call write_var_primary( filename, ncid, id_var_A    , mesh%A    )
    call write_var_primary( filename, ncid, id_var_lon  , mesh%lon  )
    call write_var_primary( filename, ncid, id_var_lat  , mesh%lat  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_mesh_in_netcdf_file

  subroutine write_matrix_operators_to_netcdf_file( filename, ncid, mesh)
    !< Write all the matrix operators to the netcdf output file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(inout) :: ncid
    type(type_mesh),  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_matrix_operators_to_netcdf_file'

    call init_routine( routine_name)

    call write_mesh_translation_tables_to_netcdf_file( filename, ncid, mesh)

    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_a_a, 'M_ddx_a_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_a_b, 'M_ddx_a_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_b_a, 'M_ddx_b_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_b_b, 'M_ddx_b_b')

    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_a_a, 'M_ddy_a_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_a_b, 'M_ddy_a_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_b_a, 'M_ddy_b_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_b_b, 'M_ddy_b_b')

    ! call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_a_a, 'M_map_a_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_a_b, 'M_map_a_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_b_a, 'M_map_b_a')
    ! call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_b_b, 'M_map_b_b')

    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_ddx_b_b   , 'M2_ddx_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_ddy_b_b   , 'M2_ddy_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_d2dx2_b_b , 'M2_d2dx2_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_d2dxdy_b_b, 'M2_d2dxdy_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_d2dy2_b_b , 'M2_d2dy2_b_b')

    call finalise_routine( routine_name)

  end subroutine write_matrix_operators_to_netcdf_file

  subroutine write_mesh_translation_tables_to_netcdf_file( filename, ncid, mesh)
    !< Write the mesh translation tables to the netcdf output file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(inout) :: ncid
    type(type_mesh),  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_mesh_translation_tables_to_netcdf_file'
    integer :: id_dim_nz, id_dim_nzp1, id_dim_vi, id_dim_ti, id_dim_ei, id_dim_two, id_dim_three
    integer :: ierr
    integer :: grp_ncid
    integer :: id_dim_nna, id_dim_nnauv, id_dim_nnak, id_dim_nnaks, id_dim_nnakuv, id_dim_nnaksuv
    integer :: id_var_n2vi, id_var_n2viuv, id_var_n2vik, id_var_n2vikuv, id_var_n2viks, id_var_n2viksuv
    integer :: id_var_vi2n, id_var_viuv2n, id_var_vik2n, id_var_vikuv2n, id_var_viks2n, id_var_viksuv2n
    integer :: id_dim_nnb, id_dim_nnbuv, id_dim_nnbk, id_dim_nnbks, id_dim_nnbkuv, id_dim_nnbksuv
    integer :: id_var_n2ti, id_var_n2tiuv, id_var_n2tik, id_var_n2tikuv, id_var_n2tiks, id_var_n2tiksuv
    integer :: id_var_ti2n, id_var_tiuv2n, id_var_tik2n, id_var_tikuv2n, id_var_tiks2n, id_var_tiksuv2n
    integer :: id_dim_nnc, id_dim_nncuv, id_dim_nnck, id_dim_nncks, id_dim_nnckuv, id_dim_nncksuv
    integer :: id_var_n2ei, id_var_n2eiuv, id_var_n2eik, id_var_n2eikuv, id_var_n2eiks, id_var_n2eiksuv
    integer :: id_var_ei2n, id_var_eiuv2n, id_var_eik2n, id_var_eikuv2n, id_var_eiks2n, id_var_eiksuv2n

    call init_routine( routine_name)

    ! Create a group for them
    ierr = NF90_DEF_GRP( ncid, 'mesh_translation_tables', grp_ncid)

    call create_dimension( filename, grp_ncid, 'nz'  , mesh%nz  , id_dim_nz)
    call create_dimension( filename, grp_ncid, 'nzp1', mesh%nz+1, id_dim_nzp1)

    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_nV    ), id_dim_vi)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_nTri  ), id_dim_ti)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_nE    ), id_dim_ei)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_two   ), id_dim_two)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_three ), id_dim_three)

    ! a-grid (vertices)
    ! =================

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'nna'    , mesh%nna    , id_dim_nna)
    call create_dimension( filename, grp_ncid, 'nna'    , mesh%nna    , id_dim_nna)
    call create_dimension( filename, grp_ncid, 'nnauv'  , mesh%nnauv  , id_dim_nnauv)
    call create_dimension( filename, grp_ncid, 'nnak'   , mesh%nnak   , id_dim_nnak)
    call create_dimension( filename, grp_ncid, 'nnaks'  , mesh%nnaks  , id_dim_nnaks)
    call create_dimension( filename, grp_ncid, 'nnakuv' , mesh%nnakuv , id_dim_nnakuv)
    call create_dimension( filename, grp_ncid, 'nnaksuv', mesh%nnaksuv, id_dim_nnaksuv)

    ! Create variables
    call create_variable( filename, grp_ncid, 'n2vi'    , NF90_INT, [id_dim_nna                             ], id_var_n2vi)
    call create_variable( filename, grp_ncid, 'n2viuv'  , NF90_INT, [id_dim_nnauv  , id_dim_two             ], id_var_n2viuv)
    call create_variable( filename, grp_ncid, 'n2vik'   , NF90_INT, [id_dim_nnak   , id_dim_two             ], id_var_n2vik)
    call create_variable( filename, grp_ncid, 'n2vikuv' , NF90_INT, [id_dim_nnakuv , id_dim_three           ], id_var_n2vikuv)
    call create_variable( filename, grp_ncid, 'n2viks'  , NF90_INT, [id_dim_nnaks  , id_dim_two             ], id_var_n2viks)
    call create_variable( filename, grp_ncid, 'n2viksuv', NF90_INT, [id_dim_nnaksuv, id_dim_three           ], id_var_n2viksuv)
    call create_variable( filename, grp_ncid, 'vi2n'    , NF90_INT, [id_dim_vi                              ], id_var_vi2n)
    call create_variable( filename, grp_ncid, 'viuv2n'  , NF90_INT, [id_dim_vi                  , id_dim_two], id_var_viuv2n)
    call create_variable( filename, grp_ncid, 'vik2n'   , NF90_INT, [id_dim_vi     , id_dim_nz              ], id_var_vik2n)
    call create_variable( filename, grp_ncid, 'vikuv2n' , NF90_INT, [id_dim_vi     , id_dim_nz,   id_dim_two], id_var_vikuv2n)
    call create_variable( filename, grp_ncid, 'viks2n'  , NF90_INT, [id_dim_vi     , id_dim_nzp1            ], id_var_viks2n)
    call create_variable( filename, grp_ncid, 'viksuv2n', NF90_INT, [id_dim_vi     , id_dim_nzp1, id_dim_two], id_var_viksuv2n)

    ! Write variables
    call write_var_primary( filename, grp_ncid, id_var_n2vi    , mesh%n2vi)
    call write_var_primary( filename, grp_ncid, id_var_n2viuv  , mesh%n2viuv)
    call write_var_primary( filename, grp_ncid, id_var_n2vik   , mesh%n2vik)
    call write_var_primary( filename, grp_ncid, id_var_n2vikuv , mesh%n2vikuv)
    call write_var_primary( filename, grp_ncid, id_var_n2viks  , mesh%n2viks)
    call write_var_primary( filename, grp_ncid, id_var_n2viksuv, mesh%n2viksuv)
    call write_var_primary( filename, grp_ncid, id_var_vi2n    , mesh%vi2n)
    call write_var_primary( filename, grp_ncid, id_var_viuv2n  , mesh%viuv2n)
    call write_var_primary( filename, grp_ncid, id_var_vik2n   , mesh%vik2n)
    call write_var_primary( filename, grp_ncid, id_var_vikuv2n , mesh%vikuv2n)
    call write_var_primary( filename, grp_ncid, id_var_viks2n  , mesh%viks2n)
    call write_var_primary( filename, grp_ncid, id_var_viksuv2n, mesh%viksuv2n)

    ! b-grid (triangles)
    ! ==================

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'nnb'    , mesh%nnb    , id_dim_nnb)
    call create_dimension( filename, grp_ncid, 'nnb'    , mesh%nnb    , id_dim_nnb)
    call create_dimension( filename, grp_ncid, 'nnbuv'  , mesh%nnbuv  , id_dim_nnbuv)
    call create_dimension( filename, grp_ncid, 'nnbk'   , mesh%nnbk   , id_dim_nnbk)
    call create_dimension( filename, grp_ncid, 'nnbks'  , mesh%nnbks  , id_dim_nnbks)
    call create_dimension( filename, grp_ncid, 'nnbkuv' , mesh%nnbkuv , id_dim_nnbkuv)
    call create_dimension( filename, grp_ncid, 'nnbksuv', mesh%nnbksuv, id_dim_nnbksuv)

    ! Create variables
    call create_variable( filename, grp_ncid, 'n2ti'    , NF90_INT, [id_dim_nnb                             ], id_var_n2ti)
    call create_variable( filename, grp_ncid, 'n2tiuv'  , NF90_INT, [id_dim_nnbuv  , id_dim_two             ], id_var_n2tiuv)
    call create_variable( filename, grp_ncid, 'n2tik'   , NF90_INT, [id_dim_nnbk   , id_dim_two             ], id_var_n2tik)
    call create_variable( filename, grp_ncid, 'n2tikuv' , NF90_INT, [id_dim_nnbkuv , id_dim_three           ], id_var_n2tikuv)
    call create_variable( filename, grp_ncid, 'n2tiks'  , NF90_INT, [id_dim_nnbks  , id_dim_two             ], id_var_n2tiks)
    call create_variable( filename, grp_ncid, 'n2tiksuv', NF90_INT, [id_dim_nnbksuv, id_dim_three           ], id_var_n2tiksuv)
    call create_variable( filename, grp_ncid, 'ti2n'    , NF90_INT, [id_dim_ti                              ], id_var_ti2n)
    call create_variable( filename, grp_ncid, 'tiuv2n'  , NF90_INT, [id_dim_ti                  , id_dim_two], id_var_tiuv2n)
    call create_variable( filename, grp_ncid, 'tik2n'   , NF90_INT, [id_dim_ti     , id_dim_nz              ], id_var_tik2n)
    call create_variable( filename, grp_ncid, 'tikuv2n' , NF90_INT, [id_dim_ti     , id_dim_nz,   id_dim_two], id_var_tikuv2n)
    call create_variable( filename, grp_ncid, 'tiks2n'  , NF90_INT, [id_dim_ti     , id_dim_nzp1            ], id_var_tiks2n)
    call create_variable( filename, grp_ncid, 'tiksuv2n', NF90_INT, [id_dim_ti     , id_dim_nzp1, id_dim_two], id_var_tiksuv2n)

    ! Write variables
    call write_var_primary( filename, grp_ncid, id_var_n2ti    , mesh%n2ti)
    call write_var_primary( filename, grp_ncid, id_var_n2tiuv  , mesh%n2tiuv)
    call write_var_primary( filename, grp_ncid, id_var_n2tik   , mesh%n2tik)
    call write_var_primary( filename, grp_ncid, id_var_n2tikuv , mesh%n2tikuv)
    call write_var_primary( filename, grp_ncid, id_var_n2tiks  , mesh%n2tiks)
    call write_var_primary( filename, grp_ncid, id_var_n2tiksuv, mesh%n2tiksuv)
    call write_var_primary( filename, grp_ncid, id_var_ti2n    , mesh%ti2n)
    call write_var_primary( filename, grp_ncid, id_var_tiuv2n  , mesh%tiuv2n)
    call write_var_primary( filename, grp_ncid, id_var_tik2n   , mesh%tik2n)
    call write_var_primary( filename, grp_ncid, id_var_tikuv2n , mesh%tikuv2n)
    call write_var_primary( filename, grp_ncid, id_var_tiks2n  , mesh%tiks2n)
    call write_var_primary( filename, grp_ncid, id_var_tiksuv2n, mesh%tiksuv2n)

    ! c-grid (edges)
    ! ==============

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'nnc'    , mesh%nnc    , id_dim_nnc)
    call create_dimension( filename, grp_ncid, 'nnc'    , mesh%nnc    , id_dim_nnc)
    call create_dimension( filename, grp_ncid, 'nncuv'  , mesh%nncuv  , id_dim_nncuv)
    call create_dimension( filename, grp_ncid, 'nnck'   , mesh%nnck   , id_dim_nnck)
    call create_dimension( filename, grp_ncid, 'nncks'  , mesh%nncks  , id_dim_nncks)
    call create_dimension( filename, grp_ncid, 'nnckuv' , mesh%nnckuv , id_dim_nnckuv)
    call create_dimension( filename, grp_ncid, 'nncksuv', mesh%nncksuv, id_dim_nncksuv)

    ! Create variables
    call create_variable( filename, grp_ncid, 'n2ei'    , NF90_INT, [id_dim_nnc                             ], id_var_n2ei)
    call create_variable( filename, grp_ncid, 'n2eiuv'  , NF90_INT, [id_dim_nncuv  , id_dim_two             ], id_var_n2eiuv)
    call create_variable( filename, grp_ncid, 'n2eik'   , NF90_INT, [id_dim_nnck   , id_dim_two             ], id_var_n2eik)
    call create_variable( filename, grp_ncid, 'n2eikuv' , NF90_INT, [id_dim_nnckuv , id_dim_three           ], id_var_n2eikuv)
    call create_variable( filename, grp_ncid, 'n2eiks'  , NF90_INT, [id_dim_nncks  , id_dim_two             ], id_var_n2eiks)
    call create_variable( filename, grp_ncid, 'n2eiksuv', NF90_INT, [id_dim_nncksuv, id_dim_three           ], id_var_n2eiksuv)
    call create_variable( filename, grp_ncid, 'ei2n'    , NF90_INT, [id_dim_ei                              ], id_var_ei2n)
    call create_variable( filename, grp_ncid, 'eiuv2n'  , NF90_INT, [id_dim_ei                  , id_dim_two], id_var_eiuv2n)
    call create_variable( filename, grp_ncid, 'eik2n'   , NF90_INT, [id_dim_ei     , id_dim_nz              ], id_var_eik2n)
    call create_variable( filename, grp_ncid, 'eikuv2n' , NF90_INT, [id_dim_ei     , id_dim_nz,   id_dim_two], id_var_eikuv2n)
    call create_variable( filename, grp_ncid, 'eiks2n'  , NF90_INT, [id_dim_ei     , id_dim_nzp1            ], id_var_eiks2n)
    call create_variable( filename, grp_ncid, 'eiksuv2n', NF90_INT, [id_dim_ei     , id_dim_nzp1, id_dim_two], id_var_eiksuv2n)

    ! Write variables
    call write_var_primary( filename, grp_ncid, id_var_n2ei    , mesh%n2ei)
    call write_var_primary( filename, grp_ncid, id_var_n2eiuv  , mesh%n2eiuv)
    call write_var_primary( filename, grp_ncid, id_var_n2eik   , mesh%n2eik)
    call write_var_primary( filename, grp_ncid, id_var_n2eikuv , mesh%n2eikuv)
    call write_var_primary( filename, grp_ncid, id_var_n2eiks  , mesh%n2eiks)
    call write_var_primary( filename, grp_ncid, id_var_n2eiksuv, mesh%n2eiksuv)
    call write_var_primary( filename, grp_ncid, id_var_ei2n    , mesh%ei2n)
    call write_var_primary( filename, grp_ncid, id_var_eiuv2n  , mesh%eiuv2n)
    call write_var_primary( filename, grp_ncid, id_var_eik2n   , mesh%eik2n)
    call write_var_primary( filename, grp_ncid, id_var_eikuv2n , mesh%eikuv2n)
    call write_var_primary( filename, grp_ncid, id_var_eiks2n  , mesh%eiks2n)
    call write_var_primary( filename, grp_ncid, id_var_eiksuv2n, mesh%eiksuv2n)

    call finalise_routine( routine_name)

  end subroutine write_mesh_translation_tables_to_netcdf_file

  subroutine write_matrix_operator_to_netcdf_file( filename, ncid, A, name)
    !< Write a single matrix operator to the netcdf output file

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    integer,                          intent(inout) :: ncid
    type(type_sparse_matrix_CSR_dp),  intent(in   ) :: A
    character(len=*),                 intent(in   ) :: name

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'write_matrix_operator_to_netcdf_file'
    type(type_sparse_matrix_CSR_dp) :: A_tot
    integer                         :: ierr
    integer                         :: grp_ncid, id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz
    integer                         :: id_var_ptr, id_var_ind, id_var_val

    call init_routine( routine_name)

    ! Gather distributed matrix to the primary
    call gather_CSR_dist_to_primary( A, A_tot)

    ! Create a new NetCDF group for this matrix operator
    ierr = NF90_DEF_GRP( ncid, name, grp_ncid)

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'm'     , A_tot%m  , id_dim_m  )
    call create_dimension( filename, grp_ncid, 'mplus1', A_tot%m+1, id_dim_mp1)
    call create_dimension( filename, grp_ncid, 'n'     , A_tot%n  , id_dim_n  )
    call create_dimension( filename, grp_ncid, 'nnz'   , A_tot%nnz, id_dim_nnz)

    ! Create variables
    call create_variable( filename, grp_ncid, 'ptr', NF90_INT   , [id_dim_mp1], id_var_ptr)
    call create_variable( filename, grp_ncid, 'ind', NF90_INT   , [id_dim_nnz], id_var_ind)
    call create_variable( filename, grp_ncid, 'val', NF90_DOUBLE, [id_dim_nnz], id_var_val)

    ! Write to NetCDF
    call write_var_primary( filename, grp_ncid, id_var_ptr, A_tot%ptr              )
    call write_var_primary( filename, grp_ncid, id_var_ind, A_tot%ind( 1:A_tot%nnz))
    call write_var_primary(  filename, grp_ncid, id_var_val, A_tot%val( 1:A_tot%nnz))

    call finalise_routine( routine_name)

  end subroutine write_matrix_operator_to_netcdf_file

end module netcdf_setup_grid_mesh_in_file
