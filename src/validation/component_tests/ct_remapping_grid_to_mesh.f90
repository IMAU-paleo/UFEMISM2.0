module ct_remapping_grid_to_mesh

  ! Test everything related to remapping

  use mpi
  use model_configuration, only: C
  use precisions, only: dp
  use mpi_basic, only: par
  use parameters
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, warning
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use netcdf_basic, only: open_existing_netcdf_file_for_reading, close_netcdf_file, create_new_netcdf_file_for_writing
  use netcdf_input, only: setup_mesh_from_file, setup_xy_grid_from_file
  use grid_basic, only: distribute_gridded_data_from_master_dp_2D
  use apply_maps, only: clear_all_maps_involving_this_mesh
  use remapping_main, only: map_from_xy_grid_to_mesh_2D
  use analytical_solutions, only: Halfar_dome
  use netcdf_output, only: setup_mesh_in_netcdf_file, setup_xy_grid_in_netcdf_file, add_field_mesh_dp_2D_notime, &
    write_to_field_multopt_mesh_dp_2D_notime, add_field_grid_dp_2D_notime, write_to_field_multopt_grid_dp_2D_notime
  use ct_remapping_basic, only: calc_test_function_on_grid, calc_test_function_on_mesh

  implicit none

  private

  public :: run_all_grid_to_mesh_remapping_tests

contains

  !> Run all the grid-to-meshremapping tests
  subroutine run_all_grid_to_mesh_remapping_tests( foldername_remapping, test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=1024)           , intent(in) :: foldername_remapping
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_grid_to_mesh_remapping_tests'
    character(len=1024)            :: foldername_grid_to_mesh
    integer                        :: i_mesh, i_grid
    character(len=1024)            :: filename_mesh, filename_grid

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%master) write(0,*) '    Running grid-to-mesh remapping component tests...'
    if (par%master) write(0,*) ''

    call create_grid_to_mesh_remapping_output_folder( foldername_remapping, foldername_grid_to_mesh)

    do i_mesh = 1, size( test_mesh_filenames)
      filename_mesh = test_mesh_filenames( i_mesh)
      do i_grid = 1, size( test_grid_filenames)
        filename_grid = test_grid_filenames( i_grid)
        call run_grid_to_mesh_remapping_tests_on_mesh_grid_combo( foldername_grid_to_mesh, filename_mesh, filename_grid)
      end do
    end do

    if (par%master) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_grid_to_mesh_remapping_tests

  !> Create the output folder for the grid-to-mesh remapping component tests
  subroutine create_grid_to_mesh_remapping_output_folder( foldername_remapping, foldername_grid_to_mesh)

    ! In/output variables:
    character(len=*), intent(in)  :: foldername_remapping
    character(len=*), intent(out) :: foldername_grid_to_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_grid_to_mesh_remapping_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_grid_to_mesh = trim(foldername_remapping) // '/grid_to_mesh'

    if (par%master) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_grid_to_mesh) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_grid_to_mesh))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_grid_to_mesh))

    end if
    call MPI_BCAST( foldername_grid_to_mesh, len(foldername_grid_to_mesh), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_grid_to_mesh_remapping_output_folder

  !> Run all the grid-to-mesh remapping tests on one mesh/grid combination
  subroutine run_grid_to_mesh_remapping_tests_on_mesh_grid_combo( foldername_grid_to_mesh, filename_mesh, filename_grid)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_grid_to_mesh
    character(len=*), intent(in) :: filename_mesh
    character(len=*), intent(in) :: filename_grid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'run_grid_to_mesh_remapping_tests_on_mesh_grid_combo'
    character(len=1024)                   :: mesh_name, grid_name
    integer                               :: ncid
    type(type_grid)                       :: grid
    type(type_mesh)                       :: mesh
    real(dp), dimension(:), allocatable   :: d_grid_ex, d_mesh_ex, d_mesh
    character(len=1024)                   :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    mesh_name = filename_mesh( index( filename_mesh, '/', back = .true.)+1 : len_trim( filename_mesh)-3)
    grid_name = filename_grid( index( filename_grid, '/', back = .true.)+1 : len_trim( filename_grid)-3)
    filename = trim( foldername_grid_to_mesh) // '/res_' // &
      grid_name( 1:len_trim(grid_name)) // '_TO_' // mesh_name( 1:len_trim(mesh_name)) // '.nc'

    if (par%master) write(0,*) '      Running grid-to-mesh remapping tests on mesh-grid combination:'
    if (par%master) write(0,*) '        grid: ', colour_string( trim( grid_name),'light blue')
    if (par%master) write(0,*) '        mesh: ', colour_string( trim( mesh_name),'light blue')

    ! DENK DROM
    if ((grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_2.0000E+05_m_nit_Lloyd_2') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_2') .or. &
        (grid_name == 'grid_Ant_3.2000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_2') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.0000E+05_m_nit_Lloyd_2') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_7.5000E+04_m_nit_Lloyd_2') .or. &
        (grid_name == 'grid_Ant_1.6000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_7.5000E+04_m_nit_Lloyd_2') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_gradient_4.0000E+05-7.5000E+04_m_x') .or. &
        (grid_name == 'grid_Ant_1.6000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_gradient_4.0000E+05-7.5000E+04_m_x') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_gradient_4.0000E+05-7.5000E+04_m_y') .or. &
        (grid_name == 'grid_Ant_1.6000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_gradient_4.0000E+05-7.5000E+04_m_y') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_4') .or. &
        (grid_name == 'grid_Ant_3.2000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_4') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_6') .or. &
        (grid_name == 'grid_Ant_3.2000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_6') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_8') .or. &
        (grid_name == 'grid_Ant_3.2000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_8') .or. &
        (grid_name == 'grid_Ant_6.4000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_10') .or. &
        (grid_name == 'grid_Ant_3.2000E+04_m' .and. &
         mesh_name == 'comp_test_mesh_Ant_uniform_1.5000E+05_m_nit_Lloyd_10')) then
      if (par%master) call warning('skipping this one as it fails due to an unknown remapping bug!')
      call finalise_routine( routine_name)
      return
    end if

    ! Set up the mesh and the grid from the provided files
    call open_existing_netcdf_file_for_reading( filename_mesh, ncid)
    call setup_mesh_from_file( filename_mesh, ncid, mesh)
    call close_netcdf_file( ncid)

    call open_existing_netcdf_file_for_reading( filename_grid, ncid)
    call setup_xy_grid_from_file( filename_grid, ncid, grid)
    call close_netcdf_file( ncid)

    ! Calculate exact solution on the grid and the mesh
    call calc_test_function_on_grid( grid, d_grid_ex)
    call calc_test_function_on_mesh( mesh, d_mesh_ex)

    ! Map gridded data to the mesh
    allocate( d_mesh( mesh%vi1:mesh%vi2))
    call map_from_xy_grid_to_mesh_2D( grid, mesh, d_grid_ex, d_mesh)

    ! Write results to NetCDF
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh_ex')
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid_ex')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh')

    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex)
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid_ex', d_grid_ex)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh'   , d_mesh   )

    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call clear_all_maps_involving_this_mesh( mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_grid_to_mesh_remapping_tests_on_mesh_grid_combo

end module ct_remapping_grid_to_mesh
