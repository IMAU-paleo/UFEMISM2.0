module ct_remapping_mesh_to_mesh

  ! Test everything related to remapping

  use mpi
  use model_configuration, only: C
  use precisions, only: dp
  use mpi_basic, only: par
  use parameters
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, warning
  use mesh_types, only: type_mesh
  use netcdf_basic, only: open_existing_netcdf_file_for_reading, close_netcdf_file, create_new_netcdf_file_for_writing, &
    create_dimension, create_variable
  use netcdf_input, only: setup_mesh_from_file
  use ct_remapping_basic, only: calc_test_function_on_mesh
  use remapping_main, only: map_from_mesh_to_mesh_2D
  use apply_maps, only: clear_all_maps_involving_this_mesh
  use netcdf_output, only: setup_mesh_in_netcdf_file, add_field_mesh_dp_2D_notime, write_to_field_multopt_mesh_dp_2D_notime
  use netcdf, only: NF90_DOUBLE, NF90_PUT_VAR
  use mpi_distributed_memory, only: gather_to_master_dp_1D

  implicit none

  private

  public :: run_all_mesh_to_mesh_remapping_tests

contains

  !> Run all the mesh-to-mesh remapping tests
  subroutine run_all_mesh_to_mesh_remapping_tests( foldername_remapping, test_mesh_filenames)

    ! In/output variables:
    character(len=1024)           , intent(in) :: foldername_remapping
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_mesh_to_mesh_remapping_tests'
    character(len=1024)            :: foldername_mesh_to_mesh
    integer                        :: i1, i2
    character(len=1024)            :: filename_mesh1, filename_mesh2

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%master) write(0,*) '    Running mesh-to-mesh remapping component tests...'
    if (par%master) write(0,*) ''

    call create_mesh_to_mesh_remapping_output_folder( foldername_remapping, foldername_mesh_to_mesh)

    do i1 = 1, size( test_mesh_filenames)
      filename_mesh1 = test_mesh_filenames( i1)
      do i2 = 1, size( test_mesh_filenames)
        filename_mesh2 = test_mesh_filenames( i2)
        call run_mesh_to_mesh_remapping_tests_on_mesh_mesh_combo( foldername_mesh_to_mesh, filename_mesh1, filename_mesh2)
      end do
    end do

    if (par%master) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_mesh_to_mesh_remapping_tests

  !> Create the output folder for the mesh-to-mesh remapping component tests
  subroutine create_mesh_to_mesh_remapping_output_folder( foldername_remapping, foldername_mesh_to_mesh)

    ! In/output variables:
    character(len=*), intent(in)  :: foldername_remapping
    character(len=*), intent(out) :: foldername_mesh_to_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_mesh_to_mesh_remapping_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_mesh_to_mesh = trim(foldername_remapping) // '/mesh_to_mesh'

    if (par%master) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_mesh_to_mesh) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_mesh_to_mesh))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_mesh_to_mesh))

    end if
    call MPI_BCAST( foldername_mesh_to_mesh, len(foldername_mesh_to_mesh), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_to_mesh_remapping_output_folder

  !> Run all the mesh-to-mesh remapping tests on one mesh/mesh combination
  subroutine run_mesh_to_mesh_remapping_tests_on_mesh_mesh_combo( foldername_mesh_to_mesh, filename_mesh1, filename_mesh2)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_mesh_to_mesh
    character(len=*), intent(in) :: filename_mesh1
    character(len=*), intent(in) :: filename_mesh2

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'run_mesh_to_mesh_remapping_tests_on_mesh_mesh_combo'
    character(len=1024)                   :: mesh_name1, mesh_name2
    integer                               :: ncid
    type(type_mesh)                       :: mesh1, mesh2
    real(dp), dimension(:), allocatable   :: d_mesh1_ex, d_mesh2_ex, d_mesh2
    real(dp), dimension(:), allocatable   :: d_mesh1_ex_tot
    character(len=1024)                   :: filename
    integer                               :: id_dim_mesh1_nV, id_var_mesh1_A, id_var_d_mesh1_ex
    integer                               :: nerr

    ! Add routine to call stack
    call init_routine( routine_name)

    mesh_name1 = filename_mesh1( index( filename_mesh1, '/', back = .true.)+1 : len_trim( filename_mesh1)-3)
    mesh_name2 = filename_mesh2( index( filename_mesh2, '/', back = .true.)+1 : len_trim( filename_mesh2)-3)
    filename = trim( foldername_mesh_to_mesh) // '/res_' // &
      mesh_name1( 1:len_trim(mesh_name1)) // '_TO_' // mesh_name2( 1:len_trim(mesh_name2)) // '.nc'

    if (par%master) write(0,*) '      Running mesh-to-mesh remapping tests on mesh-mesh combination:'
    if (par%master) write(0,*) '        from mesh: ', colour_string( trim( mesh_name1),'light blue')
    if (par%master) write(0,*) '          to mesh: ', colour_string( trim( mesh_name2),'light blue')

    ! Set up the mesh and the grid from the provided files
    call open_existing_netcdf_file_for_reading( filename_mesh1, ncid)
    call setup_mesh_from_file( filename_mesh1, ncid, mesh1)
    call close_netcdf_file( ncid)

    call open_existing_netcdf_file_for_reading( filename_mesh2, ncid)
    call setup_mesh_from_file( filename_mesh2, ncid, mesh2)
    call close_netcdf_file( ncid)

    ! Calculate exact solution on the grid and the mesh
    call calc_test_function_on_mesh( mesh1, d_mesh1_ex)
    call calc_test_function_on_mesh( mesh2, d_mesh2_ex)

    ! Map data to the new mesh
    allocate( d_mesh2( mesh2%nV_loc))
    call map_from_mesh_to_mesh_2D( mesh1, mesh2, d_mesh1_ex, d_mesh2)

    ! Write results to NetCDF
    ! (NOTE: deviates slightly from the style of the other tests, as there is no easy support
    !        for setting up two different meshes in the same NetCDF file. Luckily, we only
    !        need mesh%A from the old mesh.)

    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh2)

    call create_dimension( filename, ncid, 'mesh1_nV', mesh1%nV, id_dim_mesh1_nV)
    call create_variable( filename, ncid, 'mesh1_A'   , NF90_DOUBLE, (/ id_dim_mesh1_nV /), id_var_mesh1_A)
    call create_variable( filename, ncid, 'd_mesh1_ex', NF90_DOUBLE, (/ id_dim_mesh1_nV /), id_var_d_mesh1_ex)

    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh2_ex')
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh2')

    if (par%master) allocate( d_mesh1_ex_tot( mesh1%nV))
    call gather_to_master_dp_1D( d_mesh1_ex, d_mesh1_ex_tot)

    if (par%master) then
      nerr = NF90_PUT_VAR( ncid, id_var_mesh1_A   , mesh1%A)
      nerr = NF90_PUT_VAR( ncid, id_var_d_mesh1_ex, d_mesh1_ex_tot)
    end if

    call write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd_mesh2_ex', d_mesh2_ex)
    call write_to_field_multopt_mesh_dp_2D_notime( mesh2, filename, ncid, 'd_mesh2'   , d_mesh2)

    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call clear_all_maps_involving_this_mesh( mesh2)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_mesh_to_mesh_remapping_tests_on_mesh_mesh_combo

end module ct_remapping_mesh_to_mesh
