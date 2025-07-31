module netcdf_setup_grid_mesh_from_file
  !< Set up grids/mesh from a NetCDF file

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_DOUBLE_PRECISION
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid, type_grid_lonlat, type_grid_lat
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_parallel_creation, only: broadcast_mesh
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use grid_lonlat_basic, only: calc_lonlat_field_to_vector_form_translation_tables
  use grid_basic, only: calc_secondary_grid_data
  use netcdf_basic
  use tests_main

  implicit none

  private

  public :: setup_xy_grid_from_file, setup_lonlat_grid_from_file, setup_mesh_from_file, &
    setup_zeta_from_file, setup_depth_from_file, setup_lonlat_grid_from_lat_file

contains

  subroutine setup_xy_grid_from_file( filename, ncid, grid)
    !< Set up an x/y-grid from a NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    type(type_grid),  intent(  out) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_xy_grid_from_file'
    real(dp), parameter            :: tol = 1E-9_dp
    integer                        :: id_dim_x, id_dim_y
    integer                        :: id_var_x, id_var_y
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Give the grid a nice name
    grid%name = 'xy_grid_from_file_"' // trim( filename) // '"'

    ! Check grid dimensions and variables for validity
    call check_x( filename, ncid)
    call check_y( filename, ncid)

    ! Inquire x and y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x, dim_length = grid%nx)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y, dim_length = grid%ny)

    ! allocate memory for x and y
    allocate( grid%x( grid%nx))
    allocate( grid%y( grid%ny))

    ! Inquire x and y variables
    call inquire_var_multopt( filename, ncid, field_name_options_x, id_var_x)
    call inquire_var_multopt( filename, ncid, field_name_options_y, id_var_y)

    ! Read x and y
    call read_var_primary(  filename, ncid, id_var_x, grid%x)
    call read_var_primary(  filename, ncid, id_var_y, grid%y)

    ! Broadcast x and y from the primary to the other processes
    call MPI_BCAST( grid%x(:), grid%nx, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( grid%y(:), grid%ny, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Calculate secondary grid geometry data
    call calc_secondary_grid_data( grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_xy_grid_from_file

  subroutine setup_lonlat_grid_from_file( filename, ncid, grid)
    !< Set up a lon/lat-grid from a NetCDF file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    type(type_grid_lonlat), intent(  out) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_lonlat_grid_from_file'
    integer                        :: id_dim_lon, id_dim_lat
    integer                        :: id_var_lon, id_var_lat
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Give the grid a nice name
    grid%name = 'lonlat_grid_from_file_"' // trim( filename) // '"'

    ! Check grid dimensions and variables for validity
    call check_lon( filename, ncid)
    call check_lat( filename, ncid)

    ! Inquire lon and lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon, dim_length = grid%nlon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = grid%nlat)

    ! allocate memory for lon and lat
    allocate( grid%lon( grid%nlon))
    allocate( grid%lat( grid%nlat))

    ! Inquire lon and lat variables
    call inquire_var_multopt( filename, ncid, field_name_options_lon, id_var_lon)
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read x and y
    call read_var_primary( filename, ncid, id_var_lon, grid%lon)
    call read_var_primary( filename, ncid, id_var_lat, grid%lat)

    ! Broadcast x and y from the primary to the other processes
    call MPI_BCAST( grid%lon(:), grid%nlon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( grid%lat(:), grid%nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Secondary data
    call calc_lonlat_field_to_vector_form_translation_tables( grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_lonlat_grid_from_file

  subroutine setup_mesh_from_file( filename, ncid, mesh)
    !< Set up a mesh from a NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    type(type_mesh),  intent(inout) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_mesh_from_file'
    character(len=1024)            :: name
    integer                        :: id_dim_vi, id_dim_ti, id_dim_ci, id_dim_two, id_dim_three
    integer                        :: nV_mem, nTri_mem, n_two, n_three
    integer                        :: id_var_xmin, id_var_xmax, id_var_ymin, id_var_ymax, id_var_tol_dist, id_var_lambda_M, id_var_phi_M, id_var_beta_stereo
    integer                        :: id_var_V, id_var_nC, id_var_C, id_var_niTri, id_var_iTri, id_var_VBI
    integer                        :: id_var_Tri, id_var_Tricc, id_var_TriC

    ! Add routine to path
    call init_routine( routine_name)

    ! Give the mesh a nice name
    name = 'mesh_from_file_"' // trim( filename) // '"'

    ! Check mesh dimensions and variables for validity
    call check_mesh_dimensions( filename, ncid)

    ! Inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV    , id_dim_vi   , dim_length = nV_mem  )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri  , id_dim_ti   , dim_length = nTri_mem)
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_two   , id_dim_two  , dim_length = n_two   )
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_three , id_dim_three, dim_length = n_three )

    ! allocate memory for the mesh
    if (par%primary) then
      call allocate_mesh_primary( mesh, name, nV_mem, nTri_mem)
      mesh%nV   = mesh%nV_mem
      mesh%nTri = mesh%nTri_mem
    end if

    ! == Inquire mesh variables
    ! =========================

    ! Metadata
    call inquire_var_multopt( filename, ncid, 'xmin'                           , id_var_xmin          )
    call inquire_var_multopt( filename, ncid, 'xmax'                           , id_var_xmax          )
    call inquire_var_multopt( filename, ncid, 'ymin'                           , id_var_ymin          )
    call inquire_var_multopt( filename, ncid, 'ymax'                           , id_var_ymax          )
    call inquire_var_multopt( filename, ncid, 'tol_dist'                       , id_var_tol_dist      )
    call inquire_var_multopt( filename, ncid, 'lambda_M'                       , id_var_lambda_M      )
    call inquire_var_multopt( filename, ncid, 'phi_M'                          , id_var_phi_M         )
    call inquire_var_multopt( filename, ncid, 'beta_stereo'                    , id_var_beta_stereo   )

    ! Vertex data
    call inquire_var_multopt( filename, ncid, field_name_options_V             , id_var_V             )
    call inquire_var_multopt( filename, ncid, field_name_options_nC            , id_var_nC            )
    call inquire_var_multopt( filename, ncid, field_name_options_C             , id_var_C             )
    call inquire_var_multopt( filename, ncid, field_name_options_niTri         , id_var_niTri         )
    call inquire_var_multopt( filename, ncid, field_name_options_iTri          , id_var_iTri          )
    call inquire_var_multopt( filename, ncid, field_name_options_VBI           , id_var_VBI           )

    ! Triangle data
    call inquire_var_multopt( filename, ncid, field_name_options_Tri           , id_var_Tri           )
    call inquire_var_multopt( filename, ncid, field_name_options_Tricc         , id_var_Tricc         )
    call inquire_var_multopt( filename, ncid, field_name_options_TriC          , id_var_TriC          )

    ! == Read mesh data
    ! =================

    ! Metadata
    call read_var_primary(  filename, ncid, id_var_xmin          , mesh%xmin          )
    call read_var_primary(  filename, ncid, id_var_xmax          , mesh%xmax          )
    call read_var_primary(  filename, ncid, id_var_ymin          , mesh%ymin          )
    call read_var_primary(  filename, ncid, id_var_ymax          , mesh%ymax          )
    call read_var_primary(  filename, ncid, id_var_tol_dist      , mesh%tol_dist      )
    call read_var_primary(  filename, ncid, id_var_lambda_M      , mesh%lambda_M      )
    call read_var_primary(  filename, ncid, id_var_phi_M         , mesh%phi_M         )
    call read_var_primary(  filename, ncid, id_var_beta_stereo   , mesh%beta_stereo   )

    ! Vertex data
    call read_var_primary( filename, ncid, id_var_V             , mesh%V             )
    call read_var_primary( filename, ncid, id_var_nC            , mesh%nC            )
    call read_var_primary( filename, ncid, id_var_C             , mesh%C             )
    call read_var_primary( filename, ncid, id_var_niTri         , mesh%niTri         )
    call read_var_primary( filename, ncid, id_var_iTri          , mesh%iTri          )
    call read_var_primary( filename, ncid, id_var_VBI           , mesh%VBI           )

    ! Triangle data
    call read_var_primary( filename, ncid, id_var_Tri           , mesh%Tri           )
    call read_var_primary( filename, ncid, id_var_Tricc         , mesh%Tricc         )
    call read_var_primary( filename, ncid, id_var_TriC          , mesh%TriC          )

    ! Safety - check if the mesh data read from NetCDF makes sense
    if (par%primary) then
      if (.not. test_mesh_is_self_consistent( mesh)) call crash('invalid mesh in file ' // trim( filename))
    end if

    ! Broadcast read mesh from the primary to the other processes
    call broadcast_mesh( mesh)

    ! Calculate secondary mesh data
    call calc_all_secondary_mesh_data( mesh, mesh%lambda_M, mesh%phi_M, mesh%beta_stereo)

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_mesh_from_file

  subroutine setup_zeta_from_file( filename, ncid, nzeta, zeta)
    !< Set up a zeta coordinate from a NetCDF file

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(  out) :: nzeta
    real(dp), dimension(:), allocatable, intent(  out) ::  zeta

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_zeta_from_file'
    integer                        :: id_dim_zeta, id_var_zeta
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Check zeta dimension and variable for validity
    call check_zeta( filename, ncid)

    ! Inquire zeta dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta, dim_length = nzeta)

    ! Inquire zeta variable
    call inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var_zeta)

    ! allocate memory
    allocate( zeta( nzeta))

    ! Read zeta from file
    call read_var_primary( filename, ncid, id_var_zeta, zeta)

    ! Broadcast zeta from primary to all other processes
    call MPI_BCAST( zeta(:), nzeta, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_zeta_from_file

  subroutine setup_depth_from_file( filename, ncid, ndepth, depth)
    !< Set up a depth coordinate from a NetCDF file

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    integer,                             intent(in   ) :: ncid
    integer,                             intent(  out) :: ndepth
    real(dp), dimension(:), allocatable, intent(  out) ::  depth

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_depth_from_file'
    integer                        :: id_dim_depth, id_var_depth
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Check depth dimension and variable for validity
    call check_depth( filename, ncid)

    ! Inquire depth dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth, dim_length = ndepth)

    ! Inquire depth variable
    call inquire_var_multopt( filename, ncid, field_name_options_depth, id_var_depth)

    ! allocate memory
    allocate( depth( ndepth))

    ! Read depth from file
    call read_var_primary( filename, ncid, id_var_depth, depth)

    ! Broadcast depth from primary to all other processes
    call MPI_BCAST( depth(:), ndepth, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_depth_from_file

  subroutine setup_lonlat_grid_from_lat_file( filename, ncid, grid, vec)

    !< Set up a lat-only and a lon/lat-grid from a NetCDF file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    type(type_grid_lonlat), intent(  out) :: grid
    type(type_grid_lat),    intent(  out) :: vec

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'setup_lonlat_grid_from_lat_file'
    integer                              :: id_dim_lat
    integer                              :: id_var_lat
    integer                              :: ierr, i
    integer, parameter                   :: nlon = 360
    real(dp), parameter                  :: dlon = 1.0
    real(dp), allocatable, dimension(:)  :: lon
    CHARACTER(LEN=256)                   :: str

    ! Add routine to path
    call init_routine( routine_name)

    ! Give the grid a nice name
    grid%name = 'lonlat_grid_from_file_"' // trim( filename) // '"'

    ! Generate the vector of longitudes - hardcoded to be at 1 degree resolution
    grid%nlon = nlon
    allocate( grid%lon( grid%nlon))
    do i=1, nlon
        !lon(i) = -180 + 360 * (i - 1) / (nlon - 1) ! longitude going from -180 to + 180
        grid%lon(i) = i*dlon                             ! longitude going from 0 to 360
    end do

    ! Check latitude grid dimension and variables for validity
    call check_lat( filename, ncid)

    ! Inquire lat dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat, dim_length = grid%nlat)
    vec%nlat  = grid%nlat

    ! allocate memory for lat
    allocate( grid%lat( grid%nlat))
    allocate(  vec%lat(  vec%nlat))

    ! Inquire lat variable
    call inquire_var_multopt( filename, ncid, field_name_options_lat, id_var_lat)

    ! Read and assign y
    call read_var_primary( filename, ncid, id_var_lat, vec%lat)
    grid%lat  = vec%lat

    ! Broadcast x and y from the master to the other processes
    call MPI_BCAST( grid%lon(:), grid%nlon, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( grid%lat(:), grid%nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(  vec%lat(:),  vec%nlat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Secondary data
    call calc_lonlat_field_to_vector_form_translation_tables( grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_lonlat_grid_from_lat_file

end module netcdf_setup_grid_mesh_from_file
