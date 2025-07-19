module ct_create_test_meshes

  ! Create the suite of test meshes for the component tests.

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use tests_main
  use assertions_basic
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid
  use grid_basic, only: setup_square_grid
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform, refine_mesh_polygon
  use mesh_refinement_fun, only: mesh_add_smileyface, mesh_add_UFEMISM_letters
  use mesh_Lloyds_algorithm, only: Lloyds_algorithm_single_iteration
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use netcdf_io_main

  implicit none

  private

  public :: create_all_test_meshes_and_grids

contains

  !> Create the suite of test meshes for the component tests.
  subroutine create_all_test_meshes_and_grids( output_dir, test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=*),                            intent(in   ) :: output_dir
    character(len=*), dimension(:), allocatable, intent(  out) :: test_mesh_filenames
    character(len=*), dimension(:), allocatable, intent(  out) :: test_grid_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_all_test_meshes_and_grids'
    character(len=1024)            :: foldername_test_meshes

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '  Creating suite of test meshes and grids for component tests...'
    if (par%primary) write(0,*) ''

    ! Create output folder within the component tests folder to store the meshes
    call create_component_tests_mesh_output_folder( output_dir, foldername_test_meshes)

    ! Create all the test meshes and grids for the Antarctica domain
    call create_all_test_meshes_and_grids_Antarctica( foldername_test_meshes, &
      test_mesh_filenames, test_grid_filenames)

    if (par%primary) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_all_test_meshes_and_grids

  !> Create the output folder for the mesh component tests
  subroutine create_component_tests_mesh_output_folder( output_dir, foldername_test_meshes)

    ! In/output variables:
    character(len=*), intent(in   ) :: output_dir
    character(len=*), intent(  out) :: foldername_test_meshes

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_component_tests_mesh_output_folder'
    integer                        :: ierr
    logical                        :: ex

    ! Add routine to path
    call init_routine( routine_name)

    foldername_test_meshes = trim(output_dir) // '/test_meshes_and_grids'

    if (par%primary) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_test_meshes) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_test_meshes))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_test_meshes))

    end if
    call MPI_BCAST( foldername_test_meshes, len(foldername_test_meshes), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_component_tests_mesh_output_folder

  !> Create all the test meshes and grids for the Antarctica domain
  subroutine create_all_test_meshes_and_grids_Antarctica( foldername_test_meshes, &
    test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=*),                            intent(in)  :: foldername_test_meshes
    character(len=*), dimension(:), allocatable, intent(out) :: test_mesh_filenames
    character(len=*), dimension(:), allocatable, intent(out) :: test_grid_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_all_test_meshes_Antarctica'
    character(len=1024), parameter :: domain_name = 'Ant'
    real(dp), parameter            :: xmin                  = -3040E3_dp
    real(dp), parameter            :: xmax                  =  3040E3_dp
    real(dp), parameter            :: ymin                  = -3040E3_dp
    real(dp), parameter            :: ymax                  =  3040E3_dp
    real(dp), parameter            :: lambda_M              = 0._dp
    real(dp), parameter            :: phi_M                 = -90._dp
    real(dp), parameter            :: beta_stereo           = 71._dp

    ! Add routine to call stack
    call init_routine( routine_name)

    call create_all_test_meshes( foldername_test_meshes, test_mesh_filenames, &
      domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo)

    call create_all_test_grids( foldername_test_meshes, test_grid_filenames, &
      domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_all_test_meshes_and_grids_Antarctica

  !> Create all the test meshes for any particular domain
  subroutine create_all_test_meshes( foldername_test_meshes, test_mesh_filenames, &
    domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo)

    ! In/output variables:
    character(len=*),                            intent(in)  :: foldername_test_meshes
    character(len=*), dimension(:), allocatable, intent(out) :: test_mesh_filenames
    character(len=*),                            intent(in)  :: domain_name
    real(dp),                                    intent(in)  :: xmin, xmax, ymin, ymax
    real(dp),                                    intent(in)  :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'create_all_test_meshes'
    real(dp), parameter               :: alpha_min             = 0.4363_dp
    integer, parameter                :: nit_Lloyds_algorithm  = 2
    real(dp)                          :: res_min               = 400e3_dp
    real(dp)                          :: res_max               = 75e3_dp
    integer, dimension(4), parameter  :: nits_Lloyds_algorithm = [4,6,8,10]
    real(dp), dimension(6), parameter :: uniform_resolutions   = [400e3_dp, 300e3_dp, 200e3_dp, 150e3_dp, 100e3_dp, 75e3_dp]
    real(dp), parameter               :: uniform_resolution    = 150e3_dp
    integer                           :: i
    character(len=1)                  :: orientation

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a set of test meshes with a uniform resolution
    do i = 1, size( uniform_resolutions)
      call create_test_mesh_uniform( foldername_test_meshes, test_mesh_filenames, &
        domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        uniform_resolutions( i), alpha_min, nit_Lloyds_algorithm)
    end do

    ! Create a set of test meshes with a resolution gradient
    orientation          = 'x'
    call create_test_mesh_gradient( foldername_test_meshes, test_mesh_filenames, &
      domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
      res_min, res_max, alpha_min, nit_Lloyds_algorithm, orientation)
    orientation          = 'y'
    call create_test_mesh_gradient( foldername_test_meshes, test_mesh_filenames, &
      domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
      res_min, res_max, alpha_min, nit_Lloyds_algorithm, orientation)

    ! Create a set of otherwise identical test meshes with increasing smoothness
    do i = 1, size( nits_Lloyds_algorithm)
      call create_test_mesh_uniform( foldername_test_meshes, test_mesh_filenames, &
        domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
        uniform_resolution, alpha_min, nits_Lloyds_algorithm( i))
    end do

    ! Create a test mesh with some fun elements
    call create_test_mesh_fun( foldername_test_meshes, test_mesh_filenames, &
      domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
      uniform_resolution / 5._dp, alpha_min, nit_Lloyds_algorithm)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_all_test_meshes

  !> Create all the test grids for any particular domain
  subroutine create_all_test_grids( foldername_test_meshes, test_grid_filenames, &
    domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo)

    ! In/output variables:
    character(len=*),                            intent(in)    :: foldername_test_meshes
    character(len=*), dimension(:), allocatable, intent(inout) :: test_grid_filenames
    character(len=*),                            intent(in)    :: domain_name
    real(dp),                                    intent(in)    :: xmin, xmax, ymin, ymax
    real(dp),                                    intent(in)    :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'create_all_test_grids'
    real(dp), dimension(5), parameter :: resolutions = [128e3_dp, 64e3_dp, 32e3_dp, 16e3_dp, 8e3_dp]
    integer                           :: i
    real(dp)                          :: resolution
    character(len=1024)               :: resolution_str
    character(len=1024)               :: grid_name
    type(type_grid)                   :: grid
    character(len=1024)               :: filename
    integer                           :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a set of test grids covering the unit square with a uniform resolution
    do i = 1, size( resolutions)

      resolution = resolutions( i)
      write( resolution_str,'(es14.4)') resolution

      grid_name = 'grid_'//trim(domain_name)//'_'//trim(adjustl(resolution_str))//'_m'
      call setup_square_grid( grid_name, xmin, xmax, ymin, ymax, resolution, &
        grid, lambda_M = lambda_M, phi_M = phi_M, beta_stereo = beta_stereo)

      ! Write to NetCDF file
      filename = trim( foldername_test_meshes)//'/'//trim( grid_name)//'.nc'
      call create_new_netcdf_file_for_writing( filename, ncid)
      call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
      call close_netcdf_file( ncid)

      ! Add it to the list of test grids
      call add_filename_to_list_of_filenames( test_grid_filenames, filename)

      if (par%primary) write(0,*) '    Created test grid ', colour_string( trim( grid_name), 'light blue')

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_all_test_grids

  !> Create a test mesh with a uniform resolution
  subroutine create_test_mesh_uniform( foldername_test_meshes, test_mesh_filenames, &
    domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
    res_max, alpha_min, nit_Lloyds_algorithm)

    ! In/output variables:
    character(len=*),                            intent(in)    :: foldername_test_meshes
    character(len=*), dimension(:), allocatable, intent(inout) :: test_mesh_filenames
    character(len=*),                            intent(in)    :: domain_name
    real(dp),                                    intent(in)    :: xmin, xmax, ymin, ymax
    real(dp),                                    intent(in)    :: lambda_M, phi_M, beta_stereo
    real(dp),                                    intent(in)    :: res_max
    real(dp),                                    intent(in)    :: alpha_min
    integer,                                     intent(in)    :: nit_Lloyds_algorithm

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_test_mesh_uniform'
    type(type_mesh)                :: mesh
    character(len=14)              :: res_max_str
    character(len=3)               :: nit_Lloyds_algorithm_str
    character(len=1024)            :: mesh_name
    integer                        :: it_Lloyds_algorithm
    integer                        :: ncid
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    write( res_max_str             ,'(es14.4)') res_max
    write( nit_Lloyds_algorithm_str,'(i3)')     nit_Lloyds_algorithm
    mesh_name = 'mesh_'//trim(domain_name)//'_uniform_'//trim(adjustl(res_max_str))//'_m'//&
      '_nit_Lloyd_'//trim(adjustl(nit_Lloyds_algorithm_str))

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(mesh_name), 1000, 2000)

    ! Initialise the dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent after refine_mesh_uniform')

    ! Smooth the mesh
    do it_Lloyds_algorithm = 1, nit_Lloyds_algorithm
      call Lloyds_algorithm_single_iteration( mesh, alpha_min)
      call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent after Lloyds_algorithm_single_iteration')
    end do

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! Write to NetCDF file
    filename = trim( foldername_test_meshes)//'/'//trim( mesh_name)//'.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call write_matrix_operators_to_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Add it to the list of test meshes
    call add_filename_to_list_of_filenames( test_mesh_filenames, filename)

    if (par%primary) write(0,*) '    Created test mesh ', colour_string( trim( mesh_name), 'light blue')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_test_mesh_uniform

  !> Create a test mesh with a resolution gradient
  subroutine create_test_mesh_gradient( foldername_test_meshes, test_mesh_filenames, &
    domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
    res_min, res_max, alpha_min, nit_Lloyds_algorithm, orientation)

    ! In/output variables:
    character(len=*),                            intent(in)    :: foldername_test_meshes
    character(len=*), dimension(:), allocatable, intent(inout) :: test_mesh_filenames
    character(len=*),                            intent(in)    :: domain_name
    real(dp),                                    intent(in)    :: xmin, xmax, ymin, ymax
    real(dp),                                    intent(in)    :: lambda_M, phi_M, beta_stereo
    real(dp),                                    intent(in)    :: res_min
    real(dp),                                    intent(in)    :: res_max
    real(dp),                                    intent(in)    :: alpha_min
    integer,                                     intent(in)    :: nit_Lloyds_algorithm
    character(len=*),                            intent(in)    :: orientation

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_test_mesh_unit_square_gradient'
    type(type_mesh)                :: mesh
    character(len=14)              :: res_min_str, res_max_str
    character(len=1024)            :: mesh_name
    integer                        :: n,i
    real(dp)                       :: log_res_min, log_res_max, log_res, resolution
    real(dp)                       :: x_lo, y_lo, x_hi, y_hi
    real(dp), dimension(4,2)       :: poly
    integer                        :: it_Lloyds_algorithm
    integer                        :: ncid
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    write( res_min_str,'(es14.4)') res_min
    write( res_max_str,'(es14.4)') res_max
    mesh_name = 'mesh_'//trim(domain_name)//'_gradient_'//trim(adjustl(res_min_str))//&
      '-'//trim(adjustl(res_max_str))//'_m_'//orientation

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(mesh_name), 1000, 2000)

    ! Initialise the dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh
    n = 15
    do i = 1, n

      ! Smoothly increase resolution
      log_res_min = log( res_min)
      log_res_max = log( res_max)
      log_res = log_res_min + (log_res_max - log_res_min) * real(i-1,dp) / real(n-1,dp)
      resolution = exp( log_res)

      ! Define polygon
      select case (orientation)
      case default
        call crash('orientation should be x or y')
      case ('x')
        x_lo = mesh%xmin + (mesh%xmax - mesh%xmin) * real(i-1,dp) / real(n,dp)
        x_hi = mesh%xmax
        y_lo = mesh%ymin
        y_hi = mesh%ymax
      case ('y')
        x_lo = mesh%xmin
        x_hi = mesh%xmax
        y_lo = mesh%ymin + (mesh%ymax - mesh%ymin) * real(i-1,dp) / real(n,dp)
        y_hi = mesh%ymax
      end select
      poly( 1,:) = [x_lo, y_lo]
      poly( 2,:) = [x_hi, y_lo]
      poly( 3,:) = [x_hi, y_hi]
      poly( 4,:) = [x_lo, y_hi]

      ! Refine mesh over polygon
      call refine_mesh_polygon( mesh, poly, resolution, alpha_min)
      call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent after refine_mesh_polygon')

    end do

    ! Smooth the mesh
    do it_Lloyds_algorithm = 1, nit_Lloyds_algorithm
      call Lloyds_algorithm_single_iteration( mesh, alpha_min)
      call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent after Lloyds_algorithm_single_iteration')
    end do

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! Write to NetCDF file
    filename = trim( foldername_test_meshes)//'/'//trim( mesh_name)//'.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call write_matrix_operators_to_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Add it to the list of test meshes
    call add_filename_to_list_of_filenames( test_mesh_filenames, filename)

    if (par%primary) write(0,*) '    Created test mesh ', colour_string( trim( mesh_name), 'light blue')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_test_mesh_gradient

  !> Create a test mesh with a smileyface and the UFEMISM lettesr
  subroutine create_test_mesh_fun( foldername_test_meshes, test_mesh_filenames, &
    domain_name, xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
    res_max, alpha_min, nit_Lloyds_algorithm)

    ! In/output variables:
    character(len=*),                            intent(in)    :: foldername_test_meshes
    character(len=*), dimension(:), allocatable, intent(inout) :: test_mesh_filenames
    character(len=*),                            intent(in)    :: domain_name
    real(dp),                                    intent(in)    :: xmin, xmax, ymin, ymax
    real(dp),                                    intent(in)    :: lambda_M, phi_M, beta_stereo
    real(dp),                                    intent(in)    :: res_max
    real(dp),                                    intent(in)    :: alpha_min
    integer,                                     intent(in)    :: nit_Lloyds_algorithm

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_test_mesh_fun'
    type(type_mesh)                :: mesh
    character(len=14)              :: res_max_str
    character(len=3)               :: nit_Lloyds_algorithm_str
    character(len=1024)            :: mesh_name
    integer                        :: it_Lloyds_algorithm
    integer                        :: ncid
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    write( res_max_str             ,'(es14.4)') res_max
    write( nit_Lloyds_algorithm_str,'(i3)')     nit_Lloyds_algorithm
    mesh_name = 'mesh_'//trim(domain_name)//'_fun_'//trim(adjustl(res_max_str))//'_m'//&
      '_nit_Lloyd_'//trim(adjustl(nit_Lloyds_algorithm_str))

    ! Allocate memory
    call allocate_mesh_primary( mesh, trim(mesh_name), 1000, 2000)

    ! Initialise the dummy mesh
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the mesh around the smileyface
    call mesh_add_smileyface( mesh, res_max, res_max)
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent after mesh_add_smileyface')

    ! Refine the mesh around the UFEMISM letters
    call mesh_add_UFEMISM_letters( mesh, res_max, res_max)
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent after mesh_add_UFEMISM_letters')

    ! Smooth the mesh
    do it_Lloyds_algorithm = 1, nit_Lloyds_algorithm
      call Lloyds_algorithm_single_iteration( mesh, alpha_min)
      call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent after Lloyds_algorithm_single_iteration')
    end do

    ! Calculate secondary geometry data (needed in order to be able to write to NetCDF)
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! Write to NetCDF file
    filename = trim( foldername_test_meshes)//'/'//trim( mesh_name)//'.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call write_matrix_operators_to_netcdf_file( filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Add it to the list of test meshes
    call add_filename_to_list_of_filenames( test_mesh_filenames, filename)

    if (par%primary) write(0,*) '    Created test mesh ', colour_string( trim( mesh_name), 'light blue')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_test_mesh_fun

  !> Add a filename to an allocatable list of filenames
  subroutine add_filename_to_list_of_filenames( list_of_filenames, filename)

    ! In/output variables:
    character(len=1024), dimension(:), allocatable, intent(inout) :: list_of_filenames
    character(len=1024),                            intent(in)    :: filename

    ! Local variables:
    character(len=1024), dimension(:), allocatable :: list_of_filenames_temp

    if (.not. allocated( list_of_filenames)) then
      allocate( list_of_filenames( 1))
      list_of_filenames( 1) = trim( adjustl( filename))
    else
      ! allocate( list_of_filenames_temp( size( list_of_filenames)))
      list_of_filenames_temp = list_of_filenames
      deallocate( list_of_filenames)
      allocate( list_of_filenames( size( list_of_filenames_temp)+1))
      list_of_filenames( 1:size( list_of_filenames_temp)  ) = list_of_filenames_temp
      list_of_filenames(   size( list_of_filenames_temp)+1) = trim( adjustl( filename))
    end if

  end subroutine add_filename_to_list_of_filenames

end module ct_create_test_meshes
