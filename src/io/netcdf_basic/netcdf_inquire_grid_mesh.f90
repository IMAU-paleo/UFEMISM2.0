module netcdf_inquire_grid_mesh
  !< Inquire if an existing NetCDF file contains dimensions and variables
  !< for an x/y-grid, a lon/lat-grid, or a mesh

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use netcdf_field_name_options
  use netcdf_basic_wrappers

  implicit none

  private

  public :: inquire_xy_grid, inquire_lonlat_grid, inquire_lat_grid, inquire_mesh

contains

  subroutine inquire_xy_grid( filename, has_xy_grid)
    !< Inquire if a NetCDF file contains all the dimensions and variables
    !< describing a regular x/y-grid.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    logical,          intent(  out) :: has_xy_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_xy_grid'
    integer                        :: ncid
    integer                        :: id_dim_x, id_dim_y
    integer                        :: id_var_x, id_var_y

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for x and y dimensions and variables
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)
    call inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    call inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Check if everything is there
    has_xy_grid = (&
      id_dim_x /= -1 .and. &
      id_dim_y /= -1 .and. &
      id_var_x /= -1 .and. &
      id_var_y /= -1)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_xy_grid

  subroutine inquire_lonlat_grid( filename, has_lonlat_grid)
    !< Inquire if a NetCDF file contains all the dimensions and variables
    !< describing a regular lon/lat-grid.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    logical,          intent(  out) :: has_lonlat_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_lonlat_grid'
    integer                        :: ncid
    integer                        :: id_dim_lon, id_dim_lat
    integer                        :: id_var_lon, id_var_lat

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for x and y dimensions and variables
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)
    call inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Check if everything is there
    has_lonlat_grid = (&
      id_dim_lon /= -1 .and. &
      id_dim_lat /= -1 .and. &
      id_var_lon /= -1 .and. &
      id_var_lat /= -1)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_lonlat_grid

  subroutine inquire_lat_grid( filename, has_lat_grid)
    !< Inquire if a NetCDF file contains the dimensions and variables
    !< describing a regular lat-only grid (like the insolation file).

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    logical,          intent(  out) :: has_lat_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_lat_grid'
    integer                        :: ncid
    integer                        :: id_dim_lat
    integer                        :: id_var_lat

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Look for x and y dimensions and variables
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Check if everything is there
    has_lat_grid = (&
      id_dim_lat /= -1 .and. &
      id_var_lat /= -1)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_lat_grid

  subroutine inquire_mesh( filename, has_mesh)
    !< Inquire if a NetCDF file contains all the dimensions and variables
    !< describing a mesh.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    logical,          intent(  out) :: has_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_mesh'
    integer                        :: ncid
    integer                        :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    integer                        :: id_var_V, id_var_nC, id_var_C, id_var_niTri, id_var_iTri, id_var_VBI
    integer                        :: id_var_Tri, id_var_Tricc, id_var_TriC, id_var_TrIBI

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nC_mem, id_dim_ci   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three)

    ! inquire mesh variables
    call inquire_var_multopt( filename, ncid, field_name_options_V         , id_var_V    )
    call inquire_var_multopt( filename, ncid, field_name_options_nC        , id_var_nC   )
    call inquire_var_multopt( filename, ncid, field_name_options_C         , id_var_C    )
    call inquire_var_multopt( filename, ncid, field_name_options_niTri     , id_var_niTri)
    call inquire_var_multopt( filename, ncid, field_name_options_iTri      , id_var_iTri )
    call inquire_var_multopt( filename, ncid, field_name_options_VBI       , id_var_VBI  )
    call inquire_var_multopt( filename, ncid, field_name_options_Tri       , id_var_Tri  )
    call inquire_var_multopt( filename, ncid, field_name_options_Tricc     , id_var_Tricc)
    call inquire_var_multopt( filename, ncid, field_name_options_TriC      , id_var_TriC )

    ! Check if everything is there
    has_mesh = (&
      id_dim_vi         /= -1 .and. &
      id_dim_ti         /= -1 .and. &
      id_dim_ci         /= -1 .and. &
      id_dim_two        /= -1 .and. &
      id_dim_three      /= -1 .and. &
      id_var_V          /= -1 .and. &
      id_var_nC         /= -1 .and. &
      id_var_C          /= -1 .and. &
      id_var_niTri      /= -1 .and. &
      id_var_iTri       /= -1 .and. &
      id_var_VBI        /= -1 .and. &
      id_var_Tri        /= -1 .and. &
      id_var_Tricc      /= -1 .and. &
      id_var_TriC       /= -1)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_mesh

end module netcdf_inquire_grid_mesh
