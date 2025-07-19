module ut_netcdf_xy_grid

  ! Unit tests for the netcdf i/o - x/y-grid routines

  use tests_main
  use assertions_basic
  use ut_basic
  use netcdf_io_main
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, strrep
  use grid_types, only: type_grid
  use grid_basic, only: setup_square_grid
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary

  implicit none

  private

  public :: unit_tests_netcdf_xy_grid

contains

  subroutine unit_tests_netcdf_xy_grid( test_name_parent)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_netcdf_xy_grid'
    character(len=1024), parameter :: test_name_local = 'xy_grid'
    character(len=1024)            :: test_name
    character(len=1024), parameter :: grid_name = 'test_grid'
    real(dp), parameter            :: xmin = 0._dp
    real(dp), parameter            :: xmax = 10._dp
    real(dp), parameter            :: ymin = 0._dp
    real(dp), parameter            :: ymax = 10._dp
    real(dp), parameter            :: dx   = 1._dp
    type(type_grid)                :: grid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Set up a simple test grid
    call setup_square_grid( grid_name, xmin, xmax, ymin, ymax, dx, grid)

    call unit_tests_netcdf_xy_grid_int_2D       ( test_name, grid)
    call unit_tests_netcdf_xy_grid_int_3D       ( test_name, grid)
    call unit_tests_netcdf_xy_grid_dp_2D        ( test_name, grid)
    call unit_tests_netcdf_xy_grid_dp_2D_monthly( test_name, grid)
    call unit_tests_netcdf_xy_grid_dp_3D        ( test_name, grid)
    call unit_tests_netcdf_xy_grid_dp_3D_ocean  ( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_xy_grid

  subroutine unit_tests_netcdf_xy_grid_int_2D( test_name_parent, grid)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'unit_tests_netcdf_xy_grid_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:,:), allocatable :: d
    integer, dimension(:  ), allocatable :: d_grid_vec_partial
    integer                              :: i,j
    character(len=1024)                  :: filename
    integer                              :: ncid
    type(type_grid)                      :: grid_from_file
    integer, dimension(:  ), allocatable :: d_grid_vec_partial_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( grid%nx, grid%ny))
    do i = 1, grid%nx
      do j = 1, grid%ny
        d( i,j) = i*grid%ny + j
      end do
    end do

    ! Distribute data over the processes
    allocate( d_grid_vec_partial( grid%n_loc))
    call distribute_gridded_data_from_primary( grid, d, d_grid_vec_partial)

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call add_field_grid_int_2D_notime( filename, ncid, 'd')
    call write_to_field_multopt_grid_int_2D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, grid_from_file)
    call close_netcdf_file( ncid)

    allocate( d_grid_vec_partial_from_file( grid_from_file%n_loc))
    call read_field_from_xy_file_int_2D( filename, 'd', d_grid_vec_partial_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_eq( d_grid_vec_partial, d_grid_vec_partial_from_file), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_xy_grid_int_2D

  subroutine unit_tests_netcdf_xy_grid_int_3D( test_name_parent, grid)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'unit_tests_netcdf_xy_grid_int_3D'
    character(len=1024), parameter          :: test_name_local = 'int_3D'
    character(len=1024)                     :: test_name
    integer                                 :: nz
    real(dp), dimension(:    ), allocatable :: zeta
    integer,  dimension(:,:,:), allocatable :: d
    integer,  dimension(:,:  ), allocatable :: d_grid_vec_partial
    integer                                 :: i,j,k
    character(len=1024)                     :: filename
    integer                                 :: ncid
    type(type_grid)                         :: grid_from_file
    integer,  dimension(:,:  ), allocatable :: d_grid_vec_partial_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Vertical coordinate
    nz = 10
    allocate( zeta( nz))
    do k = 1, nz
      zeta( k) = real(k-1,dp) / real(nz-1,dp)
    end do

    ! Generate test data
    allocate( d( grid%nx, grid%ny, nz))
    do i = 1, grid%nx
      do j = 1, grid%ny
        do k = 1, nz
          d( i,j,k) = i*grid%ny*nz + j*nz + k
        end do
      end do
    end do

    ! Distribute data over the processes
    allocate( d_grid_vec_partial( grid%n_loc, nz))
    call distribute_gridded_data_from_primary( grid, d, d_grid_vec_partial)

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call add_zeta_dimension_to_file( filename, ncid, zeta)
    call add_field_grid_int_3D_notime( filename, ncid, 'd')
    call write_to_field_multopt_grid_int_3D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, grid_from_file)
    call close_netcdf_file( ncid)

    allocate( d_grid_vec_partial_from_file( grid_from_file%n_loc, nz))
    call read_field_from_xy_file_int_3D( filename, 'd', d_grid_vec_partial_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_eq( d_grid_vec_partial, d_grid_vec_partial_from_file), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_xy_grid_int_3D

  subroutine unit_tests_netcdf_xy_grid_dp_2D( test_name_parent, grid)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'unit_tests_netcdf_xy_grid_dp_2D'
    character(len=1024), parameter        :: test_name_local = 'dp_2D'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), allocatable :: d
    real(dp), dimension(:  ), allocatable :: d_grid_vec_partial
    integer                               :: i,j
    character(len=1024)                   :: filename
    integer                               :: ncid
    type(type_grid)                       :: grid_from_file
    real(dp), dimension(:  ), allocatable :: d_grid_vec_partial_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( grid%nx, grid%ny))
    do i = 1, grid%nx
      do j = 1, grid%ny
        d( i,j) = real(i*grid%ny,dp) + real(j,dp)
      end do
    end do

    ! Distribute data over the processes
    allocate( d_grid_vec_partial( grid%n_loc))
    call distribute_gridded_data_from_primary( grid, d, d_grid_vec_partial)

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call add_field_grid_dp_2D_notime( filename, ncid, 'd')
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, grid_from_file)
    call close_netcdf_file( ncid)

    allocate( d_grid_vec_partial_from_file( grid_from_file%n_loc))
    call read_field_from_xy_file_dp_2D( filename, 'd', d_grid_vec_partial_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d_grid_vec_partial, d_grid_vec_partial_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_xy_grid_dp_2D

  subroutine unit_tests_netcdf_xy_grid_dp_2D_monthly( test_name_parent, grid)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'unit_tests_netcdf_xy_grid_dp_2D_monthly'
    character(len=1024), parameter          :: test_name_local = 'dp_2D_monthly'
    character(len=1024)                     :: test_name
    real(dp), dimension(:,:,:), allocatable :: d
    real(dp), dimension(:,:  ), allocatable :: d_grid_vec_partial
    integer                                 :: i,j,m
    character(len=1024)                     :: filename
    integer                                 :: ncid
    type(type_grid)                         :: grid_from_file
    real(dp), dimension(:,:  ), allocatable :: d_grid_vec_partial_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( grid%nx, grid%ny, 12))
    do i = 1, grid%nx
      do j = 1, grid%ny
        do m = 1, 12
          d( i,j,m) = real(i*grid%ny*12,dp) + real(j*12,dp) + real(m,dp)
        end do
      end do
    end do

    ! Distribute data over the processes
    allocate( d_grid_vec_partial( grid%n_loc, 12))
    call distribute_gridded_data_from_primary( grid, d, d_grid_vec_partial)

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call add_month_dimension_to_file( filename, ncid)
    call add_field_grid_dp_2D_monthly_notime( filename, ncid, 'd')
    call write_to_field_multopt_grid_dp_2D_monthly_notime( grid, filename, ncid, 'd', d_grid_vec_partial)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, grid_from_file)
    call close_netcdf_file( ncid)

    allocate( d_grid_vec_partial_from_file( grid_from_file%n_loc, 12))
    call read_field_from_xy_file_dp_2D_monthly( filename, 'd', d_grid_vec_partial_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d_grid_vec_partial, d_grid_vec_partial_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_xy_grid_dp_2D_monthly

  subroutine unit_tests_netcdf_xy_grid_dp_3D( test_name_parent, grid)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'unit_tests_netcdf_xy_grid_dp_3D'
    character(len=1024), parameter          :: test_name_local = 'dp_3D'
    character(len=1024)                     :: test_name
    integer                                 :: nz
    real(dp), dimension(:    ), allocatable :: zeta
    real(dp), dimension(:,:,:), allocatable :: d
    real(dp), dimension(:,:  ), allocatable :: d_grid_vec_partial
    integer                                 :: i,j,k
    character(len=1024)                     :: filename
    integer                                 :: ncid
    type(type_grid)                         :: grid_from_file
    real(dp), dimension(:,:  ), allocatable :: d_grid_vec_partial_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Vertical coordinate
    nz = 10
    allocate( zeta( nz))
    do k = 1, nz
      zeta( k) = real(k-1,dp) / real(nz-1,dp)
    end do

    ! Generate test data
    allocate( d( grid%nx, grid%ny, nz))
    do i = 1, grid%nx
      do j = 1, grid%ny
        do k = 1, nz
          d( i,j,k) = real(i*grid%ny*nz,dp) + real(j*nz,dp) + real(k,dp)
        end do
      end do
    end do

    ! Distribute data over the processes
    allocate( d_grid_vec_partial( grid%n_loc, nz))
    call distribute_gridded_data_from_primary( grid, d, d_grid_vec_partial)

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call add_zeta_dimension_to_file( filename, ncid, zeta)
    call add_field_grid_dp_3D_notime( filename, ncid, 'd')
    call write_to_field_multopt_grid_dp_3D_notime( grid, filename, ncid, 'd', d_grid_vec_partial)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, grid_from_file)
    call close_netcdf_file( ncid)

    allocate( d_grid_vec_partial_from_file( grid_from_file%n_loc, nz))
    call read_field_from_xy_file_dp_3D( filename, 'd', d_grid_vec_partial_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d_grid_vec_partial, d_grid_vec_partial_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_xy_grid_dp_3D

  subroutine unit_tests_netcdf_xy_grid_dp_3D_ocean( test_name_parent, grid)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_grid),  intent(in) :: grid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'unit_tests_netcdf_xy_grid_dp_3D_ocean'
    character(len=1024), parameter          :: test_name_local = 'dp_3D_ocean'
    character(len=1024)                     :: test_name
    integer                                 :: nz_ocean
    real(dp), dimension(:),     allocatable :: z_ocean
    real(dp), dimension(:,:,:), allocatable :: d
    real(dp), dimension(:,:  ), allocatable :: d_grid_vec_partial
    integer                                 :: i,j,k
    character(len=1024)                     :: filename
    integer                                 :: ncid
    type(type_grid)                         :: grid_from_file
    real(dp), dimension(:,:  ), allocatable :: d_grid_vec_partial_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Vertical coordinate
    nz_ocean = 10
    if (allocated( z_ocean)) deallocate( z_ocean)
    allocate( z_ocean( nz_ocean))
    do k = 1, nz_ocean
      z_ocean( k) = real(k-1,dp) * 500._dp
    end do

    ! Generate test data
    allocate( d( grid%nx, grid%ny, nz_ocean))
    do i = 1, grid%nx
      do j = 1, grid%ny
        do k = 1, nz_ocean
          d( i,j,k) = real(i*grid%ny*nz_ocean,dp) + real(j*nz_ocean,dp) + real(k,dp)
        end do
      end do
    end do

    ! Distribute data over the processes
    allocate( d_grid_vec_partial( grid%n_loc, nz_ocean))
    call distribute_gridded_data_from_primary( grid, d, d_grid_vec_partial)

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call add_depth_dimension_to_file( filename, ncid, z_ocean)
    call add_field_grid_dp_3D_ocean_notime( filename, ncid, 'd')
    call write_to_field_multopt_grid_dp_3D_ocean_notime( grid, filename, ncid, 'd', d_grid_vec_partial)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, grid_from_file)
    call close_netcdf_file( ncid)

    allocate( d_grid_vec_partial_from_file( grid_from_file%n_loc, nz_ocean))
    call read_field_from_xy_file_dp_3D_ocean( filename, 'd', d_grid_vec_partial_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d_grid_vec_partial, d_grid_vec_partial_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_xy_grid_dp_3D_ocean

end module ut_netcdf_xy_grid