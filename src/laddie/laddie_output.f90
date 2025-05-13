module laddie_output

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use laddie_model_types, only: type_laddie_model
  use netcdf_io_main
  use mesh_integrate_over_domain, only: integrate_over_domain, average_over_domain
  use reallocate_mod
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_SUM, MPI_COMM_WORLD

  implicit none

  private

  public :: create_laddie_output_fields_file, create_laddie_output_scalar_file, &
            write_to_laddie_output_fields_file, write_to_laddie_output_scalar_file, &
            buffer_laddie_scalars            

  interface write_buffer_to_scalar_file_single_variable
    procedure :: write_buffer_to_scalar_file_single_variable_int
    procedure :: write_buffer_to_scalar_file_single_variable_dp
  end interface write_buffer_to_scalar_file_single_variable

contains

  subroutine write_to_laddie_output_fields_file( mesh, laddie, region_name, time)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(inout) :: laddie
    character(len=3),        intent(in   ) :: region_name
    real(dp),                intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_laddie_output_fields_file'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! If the mesh has been updated, create a new output file
    if (.not. laddie%output_fields_file_matches_current_mesh) then
      call create_laddie_output_fields_file( mesh, laddie, region_name)
    end if

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( laddie%output_fields_filename, ncid)

    ! write the time to the file
    call write_time_to_file( laddie%output_fields_filename, ncid, time)

    ! write the default data fields to the file
    call write_to_field_multopt_mesh_dp_2D(   mesh, laddie%output_fields_filename, ncid, 'H_lad', laddie%now%H, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_fields_filename, ncid, 'U_lad', laddie%now%U, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_fields_filename, ncid, 'V_lad', laddie%now%V, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D(   mesh, laddie%output_fields_filename, ncid, 'T_lad', laddie%now%T, d_is_hybrid = .true.)
    call write_to_field_multopt_mesh_dp_2D(   mesh, laddie%output_fields_filename, ncid, 'S_lad', laddie%now%S, d_is_hybrid = .true.)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_laddie_output_fields_file

  subroutine create_laddie_output_fields_file( mesh, laddie, region_name)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(inout) :: laddie
    character(len=3),        intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_laddie_output_fields_file'
    character(len=1024)            :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Set filename
    filename_base = trim( C%output_dir) // 'laddie_output_fields_' // region_name
    call generate_filename_XXXXXdotnc( filename_base, laddie%output_fields_filename)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating laddie output file "' // colour_string( trim( laddie%output_fields_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( laddie%output_fields_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( laddie%output_fields_filename, ncid, mesh)

    ! Add time dimension+variable to the file
    call add_time_dimension_to_file(  laddie%output_fields_filename, ncid)

    ! Add the default data fields to the file
    call add_field_mesh_dp_2D(   laddie%output_fields_filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
    call add_field_mesh_dp_2D_b( laddie%output_fields_filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
    call add_field_mesh_dp_2D_b( laddie%output_fields_filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
    call add_field_mesh_dp_2D(   laddie%output_fields_filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
    call add_field_mesh_dp_2D(   laddie%output_fields_filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

    ! Confirm that the current output file match the current model mesh
    ! (set to false whenever a new mesh is created,
    ! and set to true whenever a new output file is created)
    laddie%output_fields_file_matches_current_mesh = .true.

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_laddie_output_fields_file

  subroutine write_to_laddie_output_scalar_file( laddie)

    ! In/output variables
    type(type_laddie_model), intent(inout) :: laddie

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_laddie_output_scalar_file'
    character(len=1024)            :: filename
    integer                        :: ncid, n, id_dim_time, ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Shorthand for variable names
    filename = laddie%output_scalar_filename
    n        = laddie%buffer%n

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to laddie scalar output file "' // colour_string( trim( laddie%output_scalar_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( filename, ncid)

    ! Inquire number of timeframes already present in the file
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write the time to the file
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'time',              laddie%buffer%time,              n, ti+1)

    ! Write bulk scalars
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'layer_volume',      laddie%buffer%layer_volume,      n, ti+1)

    ! Reset buffer
    laddie%buffer%n = 0

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_laddie_output_scalar_file

  subroutine create_laddie_output_scalar_file( laddie, region_name)
    !< Create the scalar regional output NetCDF file

    ! In/output variables:
    type(type_laddie_model), intent(inout) :: laddie
    character(len=3),        intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_laddie_output_scalar_file'
    character(len=1024)            :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = trim( C%output_dir) // 'laddie_output_scalar_' // region_name
    call generate_filename_XXXXXdotnc( filename_base, laddie%output_scalar_filename)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating laddie scalar output file "' // colour_string( trim( laddie%output_scalar_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( laddie%output_scalar_filename, ncid)

    ! Add time dimensions+variables to the file
    call add_time_dimension_to_file( laddie%output_scalar_filename, ncid)

    ! Integrated ice geometry
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'layer_volume', long_name = 'Total mixed layer volume', units = 'm^3')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Allocate buffer
    call allocate_laddie_buffer( laddie)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_laddie_output_scalar_file

  subroutine allocate_laddie_buffer( laddie)
    !< Allocate memory to buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_laddie_model), intent(inout) :: laddie

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_laddie_buffer'
    integer                        :: n_mem

    ! Add routine to path
    call init_routine( routine_name)

    laddie%buffer%n_mem = 0
    laddie%buffer%n     = 0

    ! Only allocate memory for this on the primary
    if (par%primary) then

      n_mem = 1000
      laddie%buffer%n_mem = n_mem
      laddie%buffer%n     = 0

      allocate( laddie%buffer%time             ( n_mem), source = 0._dp)

      allocate( laddie%buffer%layer_volume     ( n_mem), source = 0._dp)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_laddie_buffer

  subroutine buffer_laddie_scalars( mesh, laddie, time)
    !< Buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(inout) :: laddie
    real(dp),                intent(in)    :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'buffer_laddie_scalars'
    integer                        :: n, vi, ierr
    real(dp)                       :: H_int

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate values
    call integrate_over_domain( mesh, laddie%now%H, H_int)

    ! Only the primary does this
    if (par%primary) then

      ! Increase timeframe count
      laddie%buffer%n = laddie%buffer%n + 1
      n = laddie%buffer%n

      ! Extend buffer memory if necessary
      if (n > laddie%buffer%n_mem - 10) call extend_laddie_buffer( laddie)

      ! Store new timeframe in buffer
      laddie%buffer%time             ( n) = time

      laddie%buffer%layer_volume     ( n) = H_int
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine buffer_laddie_scalars

  subroutine extend_laddie_buffer( laddie)
    !< Extend memory to buffer the scalar output data between output writing intervals
    !
    ! NOTE: should only be called by the primary!

    ! In/output variables:
    type(type_laddie_model), intent(inout) :: laddie

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'extend_laddie_buffer'
    integer                        :: n_mem

    ! Add routine to path
    call init_routine( routine_name)

    n_mem = laddie%buffer%n_mem * 2
    laddie%buffer%n_mem = n_mem

    call reallocate( laddie%buffer%time             , n_mem, source = 0._dp)

    call reallocate( laddie%buffer%layer_volume     , n_mem, source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extend_laddie_buffer

  subroutine write_buffer_to_scalar_file_single_variable_int( filename, ncid, var_name, d, n, ti)
    !< Write buffered scalar data of a single variable to the scalar output file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: var_name
    integer,  dimension(:), intent(in   ) :: d
    integer,                intent(in   ) :: n
    integer,                intent(in   ) :: ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_buffer_to_scalar_file_single_variable_int'
    integer                        :: id_var
    integer, dimension(1)          :: start, count
    integer,  dimension(n)         :: d_to_write

    ! Add routine to path
    call init_routine( routine_name)

    call inquire_var( filename, ncid, var_name, id_var)

    start = ti
    count = n
    d_to_write = d(1:n)

    call write_var_primary(  filename, ncid, id_var, d_to_write, start = start, count = count)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_buffer_to_scalar_file_single_variable_int

  subroutine write_buffer_to_scalar_file_single_variable_dp( filename, ncid, var_name, d, n, ti)
    !< Write buffered scalar data of a single variable to the scalar output file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: var_name
    real(dp), dimension(:), intent(in   ) :: d
    integer,                intent(in   ) :: n
    integer,                intent(in   ) :: ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_buffer_to_scalar_file_single_variable_dp'
    integer                        :: id_var
    integer, dimension(1)          :: start, count
    real(dp), dimension(n)         :: d_to_write

    ! Add routine to path
    call init_routine( routine_name)

    call inquire_var( filename, ncid, var_name, id_var)

    start = ti
    count = n
    d_to_write = d(1:n)

    call write_var_primary(  filename, ncid, id_var, d_to_write, start = start, count = count)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_buffer_to_scalar_file_single_variable_dp

end module laddie_output


