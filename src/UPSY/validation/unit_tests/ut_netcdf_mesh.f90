module ut_netcdf_mesh

  ! Unit tests for the netcdf i/o - mesh routines

  use tests_main
  use assertions_basic
  use ut_basic
  use netcdf_io_main
  use precisions, only: dp
  use parameters, only: pi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, strrep
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data

  implicit none

  private

  public :: unit_tests_netcdf_mesh

contains

  subroutine unit_tests_netcdf_mesh( test_name_parent)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_netcdf_mesh'
    character(len=1024), parameter :: test_name_local = 'mesh'
    character(len=1024)            :: test_name
    character(len=1024)            :: name
    real(dp)                       :: xmin, xmax, ymin, ymax, alpha_min, res_max
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    name = 'test_mesh'
    xmin = -500e3_dp
    xmax =  500e3_dp
    ymin = -500e3_dp
    ymax =  500e3_dp
    alpha_min = 25._dp * pi / 180._dp
    res_max = 50e3_dp

    call allocate_mesh_primary( mesh, name, 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    call unit_tests_netcdf_mesh_int_2D       ( test_name, mesh)
    call unit_tests_netcdf_mesh_int_2D_b     ( test_name, mesh)
    call unit_tests_netcdf_mesh_dp_2D        ( test_name, mesh)
    call unit_tests_netcdf_mesh_dp_2D_b      ( test_name, mesh)
    call unit_tests_netcdf_mesh_dp_2D_monthly( test_name, mesh)
    call unit_tests_netcdf_mesh_dp_3D        ( test_name, mesh)
    call unit_tests_netcdf_mesh_dp_3D_b      ( test_name, mesh)
    call unit_tests_netcdf_mesh_dp_3D_ocean  ( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh

  subroutine unit_tests_netcdf_mesh_int_2D( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'unit_tests_netcdf_mesh_int_2D'
    character(len=1024), parameter       :: test_name_local = 'int_2D'
    character(len=1024)                  :: test_name
    integer, dimension(:  ), allocatable :: d
    integer                              :: vi
    character(len=1024)                  :: filename
    integer                              :: ncid
    type(type_mesh)                      :: mesh_from_file
    integer, dimension(:  ), allocatable :: d_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( mesh%vi1:mesh%vi2))
    do vi = mesh%vi1, mesh%vi2
      d( vi) = vi
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_field_mesh_int_2D_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_int_2D_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%vi1:mesh_from_file%vi2))
    call read_field_from_mesh_file_int_2D( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_eq( d, d_from_file), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_int_2D

  subroutine unit_tests_netcdf_mesh_int_2D_b( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'unit_tests_netcdf_mesh_int_2D_b'
    character(len=1024), parameter       :: test_name_local = 'int_2D_b'
    character(len=1024)                  :: test_name
    integer, dimension(:  ), allocatable :: d
    integer                              :: ti
    character(len=1024)                  :: filename
    integer                              :: ncid
    type(type_mesh)                      :: mesh_from_file
    integer, dimension(:  ), allocatable :: d_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( mesh%ti1:mesh%ti2))
    do ti = mesh%ti1, mesh%ti2
      d( ti) = ti
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_field_mesh_int_2D_b_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_int_2D_b_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%ti1:mesh_from_file%ti2))
    call read_field_from_mesh_file_int_2D_b( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_eq( d, d_from_file), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_int_2D_b

  subroutine unit_tests_netcdf_mesh_dp_2D( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'unit_tests_netcdf_mesh_dp_2D'
    character(len=1024), parameter      :: test_name_local = 'dp_2D'
    character(len=1024)                 :: test_name
    real(dp), dimension(:), allocatable :: d
    integer                             :: vi
    character(len=1024)                 :: filename
    integer                             :: ncid
    type(type_mesh)                     :: mesh_from_file
    real(dp), dimension(:), allocatable :: d_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( mesh%vi1:mesh%vi2))
    do vi = mesh%vi1, mesh%vi2
      d( vi) = real(vi,dp)
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%vi1:mesh_from_file%vi2))
    call read_field_from_mesh_file_dp_2D( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d, d_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_dp_2D

  subroutine unit_tests_netcdf_mesh_dp_2D_b( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'unit_tests_netcdf_mesh_dp_2D_b'
    character(len=1024), parameter      :: test_name_local = 'dp_2D_b'
    character(len=1024)                 :: test_name
    real(dp), dimension(:), allocatable :: d
    integer                             :: ti
    character(len=1024)                 :: filename
    integer                             :: ncid
    type(type_mesh)                     :: mesh_from_file
    real(dp), dimension(:), allocatable :: d_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( mesh%ti1:mesh%ti2))
    do ti = mesh%ti1, mesh%ti2
      d( ti) = real(ti,dp)
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%ti1:mesh_from_file%ti2))
    call read_field_from_mesh_file_dp_2D_b( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_eq( d, d_from_file), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_dp_2D_b

  subroutine unit_tests_netcdf_mesh_dp_2D_monthly( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'unit_tests_netcdf_mesh_dp_2D_monthly'
    character(len=1024), parameter        :: test_name_local = 'dp_2D_monthly'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), allocatable :: d
    integer                               :: vi,m
    character(len=1024)                   :: filename
    integer                               :: ncid
    type(type_mesh)                       :: mesh_from_file
    real(dp), dimension(:,:), allocatable :: d_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( mesh%vi1:mesh%vi2,12))
    do vi = mesh%vi1, mesh%vi2
      do m = 1, 12
        d( vi,m) = real(vi*12,dp) + real(m,dp)
      end do
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_month_dimension_to_file( filename, ncid)
    call add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_dp_2D_monthly_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%vi1:mesh_from_file%vi2,12))
    call read_field_from_mesh_file_dp_2D_monthly( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d, d_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_dp_2D_monthly

  subroutine unit_tests_netcdf_mesh_dp_3D( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'unit_tests_netcdf_mesh_dp_3D'
    character(len=1024), parameter        :: test_name_local = 'dp_3D'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), allocatable :: d
    integer                               :: vi,k
    character(len=1024)                   :: filename
    integer                               :: ncid
    type(type_mesh)                       :: mesh_from_file
    real(dp), dimension(:,:), allocatable :: d_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( mesh%vi1:mesh%vi2,mesh%nz))
    do vi = mesh%vi1, mesh%vi2
      do k = 1, mesh%nz
        d( vi,k) = real(vi*mesh%nz,dp) + real(k,dp)
      end do
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_zeta_dimension_to_file( filename, ncid, mesh%zeta)
    call add_field_mesh_dp_3D_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_dp_3D_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%vi1:mesh_from_file%vi2,mesh_from_file%nz))
    call read_field_from_mesh_file_dp_3D( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d, d_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_dp_3D

  subroutine unit_tests_netcdf_mesh_dp_3D_b( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'unit_tests_netcdf_mesh_dp_3D_b'
    character(len=1024), parameter        :: test_name_local = 'dp_3D_b'
    character(len=1024)                   :: test_name
    real(dp), dimension(:,:), allocatable :: d
    integer                               :: ti,k
    character(len=1024)                   :: filename
    integer                               :: ncid
    type(type_mesh)                       :: mesh_from_file
    real(dp), dimension(:,:), allocatable :: d_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Generate test data
    allocate( d( mesh%ti1:mesh%ti2,mesh%nz))
    do ti = mesh%ti1, mesh%ti2
      do k = 1, mesh%nz
        d( ti,k) = real(ti*mesh%nz,dp) + real(k,dp)
      end do
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_zeta_dimension_to_file( filename, ncid, mesh%zeta)
    call add_field_mesh_dp_3D_b_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%ti1:mesh_from_file%ti2,mesh_from_file%nz))
    call read_field_from_mesh_file_dp_3D_b( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d, d_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_dp_3D_b

  subroutine unit_tests_netcdf_mesh_dp_3D_ocean( test_name_parent, mesh)
    ! Test the netcdf i/o subroutines for meshed data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'unit_tests_netcdf_mesh_dp_3D_ocean'
    character(len=1024), parameter        :: test_name_local = 'dp_3D_ocean'
    character(len=1024)                   :: test_name
    integer                               :: nz_ocean
    real(dp), dimension(:),   allocatable :: z_ocean
    real(dp), dimension(:,:), allocatable :: d
    integer                               :: vi,k
    character(len=1024)                   :: filename
    integer                               :: ncid
    type(type_mesh)                       :: mesh_from_file
    real(dp), dimension(:,:), allocatable :: d_from_file

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
    allocate( d( mesh%vi1:mesh%vi2,nz_ocean))
    do vi = mesh%vi1, mesh%vi2
      do k = 1, nz_ocean
        d( vi,k) = real(vi*nz_ocean,dp) + real(k,dp)
      end do
    end do

    ! Write to output file
    filename = trim(foldername_unit_tests_output) // '/' // trim( strrep( test_name,'/','_')) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_depth_dimension_to_file( filename, ncid, z_ocean)
    call add_field_mesh_dp_3D_ocean_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_dp_3D_ocean_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Read from file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_from_file)
    call close_netcdf_file( ncid)

    allocate( d_from_file( mesh_from_file%vi1:mesh_from_file%vi2,nz_ocean))
    call read_field_from_mesh_file_dp_3D_ocean( filename, 'd', d_from_file)

    ! Test if what we read is the same as what we wrote
    call unit_test( test_tol( d, d_from_file, 1e-5_dp), trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_mesh_dp_3D_ocean

end module ut_netcdf_mesh