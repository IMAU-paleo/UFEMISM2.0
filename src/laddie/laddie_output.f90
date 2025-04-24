module laddie_output

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use laddie_model_types, only: type_laddie_model
  use netcdf_io_main

  implicit none

  private

  public :: create_laddie_output_file, write_to_laddie_output_file

contains

  subroutine write_to_laddie_output_file( mesh, laddie, region_name, time)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(inout) :: laddie
    character(len=3),        intent(in   ) :: region_name
    real(dp),                intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_laddie_output_file'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! If the mesh has been updated, create a new output file
    if (.not. laddie%output_file_matches_current_mesh) then
      call create_laddie_output_file( mesh, laddie, region_name)
    end if

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( laddie%output_filename, ncid)

    ! write the time to the file
    call write_time_to_file( laddie%output_filename, ncid, time)

    ! write the default data fields to the file
    call write_to_field_multopt_mesh_dp_2D(   mesh, laddie%output_filename, ncid, 'H_lad', laddie%now%H, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_filename, ncid, 'U_lad', laddie%now%U, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_filename, ncid, 'V_lad', laddie%now%V, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D(   mesh, laddie%output_filename, ncid, 'T_lad', laddie%now%T, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D(   mesh, laddie%output_filename, ncid, 'S_lad', laddie%now%S, d_is_hybrid = .true.)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_laddie_output_file

  subroutine create_laddie_output_file( mesh, laddie, region_name)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(inout) :: laddie
    character(len=3),        intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_laddie_output_file'
    character(len=1024)            :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Set filename
    filename_base = trim( C%output_dir) // 'laddie_output_' // region_name
    call generate_filename_XXXXXdotnc( filename_base, laddie%output_filename)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating laddie output file "' // colour_string( trim( laddie%output_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( laddie%output_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( laddie%output_filename, ncid, mesh)

    ! Add time dimension+variable to the file
    call add_time_dimension_to_file(  laddie%output_filename, ncid)

    ! Add the default data fields to the file
    call add_field_mesh_dp_2D(   laddie%output_filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
    call add_field_mesh_dp_2D_b( laddie%output_filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
    call add_field_mesh_dp_2D_b( laddie%output_filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
    call add_field_mesh_dp_2D(   laddie%output_filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
    call add_field_mesh_dp_2D(   laddie%output_filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

    ! Confirm that the current output file match the current model mesh
    ! (set to false whenever a new mesh is created,
    ! and set to true whenever a new output file is created)
    laddie%output_file_matches_current_mesh = .true.

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_laddie_output_file

end module laddie_output


