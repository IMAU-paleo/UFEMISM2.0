module laddie_output

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use parameters
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use laddie_model_types, only: type_laddie_model
  use netcdf_io_main
  use mesh_integrate_over_domain, only: integrate_over_domain, average_over_domain
  use reallocate_mod
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_MIN, MPI_SUM, MPI_COMM_WORLD
  use scalar_output_files, only: write_buffer_to_scalar_file_single_variable

  implicit none

  private

  public :: create_laddie_output_fields_file, create_laddie_output_scalar_file, &
            write_to_laddie_output_fields_file, write_to_laddie_output_scalar_file, &
            buffer_laddie_scalars

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

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( filename, ncid)

    ! Inquire number of timeframes already present in the file
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write the time to the file
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'time',              laddie%buffer%time,              n, ti+1)

    ! Write bulk scalars
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'layer_volume',      laddie%buffer%layer_volume,      n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'area_a',            laddie%buffer%area_a,            n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'area_b',            laddie%buffer%area_b,            n, ti+1)

    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'thickness_mean',    laddie%buffer%thickness_mean,    n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'thickness_max',     laddie%buffer%thickness_max,     n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'thickness_min',     laddie%buffer%thickness_min,     n, ti+1)

    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'melt_mean',         laddie%buffer%melt_mean,         n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'melt_max',          laddie%buffer%melt_max,          n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'melt_min',          laddie%buffer%melt_min,          n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'melt_tot',          laddie%buffer%melt_tot,          n, ti+1)

    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'uabs_max',          laddie%buffer%uabs_max,          n, ti+1)

    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'T_mean',            laddie%buffer%T_mean,            n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'T_max',             laddie%buffer%T_max,             n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'T_min',             laddie%buffer%T_min,             n, ti+1)

    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'S_mean',            laddie%buffer%S_mean,            n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'S_max',             laddie%buffer%S_max,             n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'S_min',             laddie%buffer%S_min,             n, ti+1)

    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'entr_tot',          laddie%buffer%entr_tot,          n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'entr_dmin_tot',     laddie%buffer%entr_dmin_tot,     n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'detr_tot',          laddie%buffer%detr_tot,          n, ti+1)
    call write_buffer_to_scalar_file_single_variable( filename, ncid, 'divQH_sum',         laddie%buffer%divQH_sum,         n, ti+1)

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
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'layer_volume',   long_name = 'Total mixed layer volume',        units = 'm^3')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'area_a',         long_name = 'Integrated floating area a-grid', units = 'm^2')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'area_b',         long_name = 'Integrated floating area b-grid', units = 'm^2')

    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'thickness_mean', long_name = 'Mean layer thickness',    units = 'm')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'thickness_max',  long_name = 'Maximum layer thickness', units = 'm')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'thickness_min',  long_name = 'Minimum layer thickness', units = 'm')

    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'melt_mean',    long_name = 'Mean melt rate',           units = 'm/yr')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'melt_max',     long_name = 'Maximum melt rate',        units = 'm/yr')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'melt_min',     long_name = 'Minimum melt rate',        units = 'm/yr')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'melt_tot',     long_name = 'Total melt rate  = BMB',   units = 'Gt/yr')

    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'uabs_max',     long_name = 'Maximum speed',            units = 'm/s')

    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'T_mean',       long_name = 'Mean temperature',         units = 'degC')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'T_max',        long_name = 'Maximum temperature',      units = 'degC')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'T_min',        long_name = 'Minimum temperature',      units = 'degC')

    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'S_mean',       long_name = 'Mean salinity',            units = 'psu')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'S_max',        long_name = 'Maximum salinity',         units = 'psu')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'S_min',        long_name = 'Minimum salinity',         units = 'psu')

    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'entr_tot',     long_name = 'Integrated entrainment',   units = 'Sv')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'entr_dmin_tot',long_name = 'Integrated entrainment for Dmin', units = 'Sv')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'detr_tot',     long_name = 'Integrated detrainment',   units = 'Sv')
    call add_field_dp_0D( laddie%output_scalar_filename, ncid, 'divQH_sum',    long_name = 'Integrated volume divergence', units = 'Sv')

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
      allocate( laddie%buffer%area_a           ( n_mem), source = 0._dp)
      allocate( laddie%buffer%area_b           ( n_mem), source = 0._dp)

      allocate( laddie%buffer%thickness_mean   ( n_mem), source = 0._dp)
      allocate( laddie%buffer%thickness_max    ( n_mem), source = 0._dp)
      allocate( laddie%buffer%thickness_min    ( n_mem), source = 0._dp)

      allocate( laddie%buffer%melt_mean        ( n_mem), source = 0._dp)
      allocate( laddie%buffer%melt_max         ( n_mem), source = 0._dp)
      allocate( laddie%buffer%melt_min         ( n_mem), source = 0._dp)
      allocate( laddie%buffer%melt_tot         ( n_mem), source = 0._dp)

      allocate( laddie%buffer%uabs_max         ( n_mem), source = 0._dp)

      allocate( laddie%buffer%T_mean           ( n_mem), source = 0._dp)
      allocate( laddie%buffer%T_max            ( n_mem), source = 0._dp)
      allocate( laddie%buffer%T_min            ( n_mem), source = 0._dp)

      allocate( laddie%buffer%S_mean           ( n_mem), source = 0._dp)
      allocate( laddie%buffer%S_max            ( n_mem), source = 0._dp)
      allocate( laddie%buffer%S_min            ( n_mem), source = 0._dp)

      allocate( laddie%buffer%entr_tot         ( n_mem), source = 0._dp)
      allocate( laddie%buffer%entr_dmin_tot    ( n_mem), source = 0._dp)
      allocate( laddie%buffer%detr_tot         ( n_mem), source = 0._dp)
      allocate( laddie%buffer%divQH_sum        ( n_mem), source = 0._dp)
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
    real(dp)                       :: H_int, H_mean, H_max, H_min
    real(dp)                       :: melt_mean, melt_max, melt_min, melt_int, melt_tot
    real(dp)                       :: Uabs_max
    real(dp)                       :: T_mean, T_max, T_min, T_int
    real(dp)                       :: S_mean, S_max, S_min, S_int
    real(dp)                       :: entr_tot, entr_dmin_tot, detr_tot, divQH_sum

    ! Add routine to path
    call init_routine( routine_name)

    ! == Calculate values ==

    ! Thickness scalars
    call integrate_over_domain( mesh, laddie%now%H, H_int)
    H_mean = H_int / laddie%area_a
    H_max = maxval( laddie%now%H, laddie%mask_a)
    call MPI_ALLREDUCE( MPI_IN_PLACE, H_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    H_min = minval( laddie%now%H, laddie%mask_a)
    call MPI_ALLREDUCE( MPI_IN_PLACE, H_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Melt scalars
    call integrate_over_domain( mesh, laddie%melt, melt_int, max_d = melt_max, min_d = melt_min)
    melt_mean = melt_int * sec_per_year / laddie%area_a
    melt_max  = melt_max * sec_per_year
    melt_min  = melt_min * sec_per_year
    melt_tot  = melt_int * sec_per_year * 1.0E-09_dp ! [Gt/yr]

    Uabs_max = maxval( sqrt( laddie%now%U**2 + laddie%now%V**2), laddie%mask_b)
    call MPI_ALLREDUCE( MPI_IN_PLACE, Uabs_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Temperature scalars
    call integrate_over_domain( mesh, laddie%now%T, T_int)
    T_mean = T_int / laddie%area_a
    T_max = maxval( laddie%now%T, laddie%mask_a)
    call MPI_ALLREDUCE( MPI_IN_PLACE, T_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    T_min = minval( laddie%now%T, laddie%mask_a)
    call MPI_ALLREDUCE( MPI_IN_PLACE, T_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Salinity scalars
    call integrate_over_domain( mesh, laddie%now%S, S_int, max_d = S_max, min_d = S_min)
    S_mean = S_int / laddie%area_a

    ! Volume fluxes
    call integrate_over_domain( mesh, laddie%entr, entr_tot)
    call integrate_over_domain( mesh, laddie%entr_dmin, entr_dmin_tot)
    call integrate_over_domain( mesh, laddie%detr, detr_tot)
    divQH_sum = sum( laddie%divQH)
    call MPI_ALLREDUCE( MPI_IN_PLACE, divQH_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! =====================

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
      laddie%buffer%area_a           ( n) = laddie%area_a
      laddie%buffer%area_b           ( n) = laddie%area_b

      laddie%buffer%thickness_mean   ( n) = H_mean
      laddie%buffer%thickness_max    ( n) = H_max
      laddie%buffer%thickness_min    ( n) = H_min

      laddie%buffer%melt_mean        ( n) = melt_mean
      laddie%buffer%melt_max         ( n) = melt_max
      laddie%buffer%melt_min         ( n) = melt_min
      laddie%buffer%melt_tot         ( n) = melt_tot

      laddie%buffer%uabs_max         ( n) = Uabs_max

      laddie%buffer%T_mean           ( n) = T_mean
      laddie%buffer%T_max            ( n) = T_max
      laddie%buffer%T_min            ( n) = T_min

      laddie%buffer%S_mean           ( n) = S_mean
      laddie%buffer%S_max            ( n) = S_max
      laddie%buffer%S_min            ( n) = S_min

      laddie%buffer%entr_tot         ( n) = entr_tot * 1.0E-6_dp
      laddie%buffer%entr_dmin_tot    ( n) = entr_dmin_tot * 1.0E-6_dp
      laddie%buffer%detr_tot         ( n) = detr_tot * 1.0E-6_dp
      laddie%buffer%divQH_sum        ( n) = divQH_sum * 1.0E-6_dp
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
    call reallocate( laddie%buffer%area_a           , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%area_b           , n_mem, source = 0._dp)

    call reallocate( laddie%buffer%thickness_mean   , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%thickness_max    , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%thickness_min    , n_mem, source = 0._dp)

    call reallocate( laddie%buffer%melt_mean        , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%melt_max         , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%melt_min         , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%melt_tot         , n_mem, source = 0._dp)

    call reallocate( laddie%buffer%uabs_max         , n_mem, source = 0._dp)

    call reallocate( laddie%buffer%T_mean           , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%T_max            , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%T_min            , n_mem, source = 0._dp)

    call reallocate( laddie%buffer%S_mean           , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%S_max            , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%S_min            , n_mem, source = 0._dp)

    call reallocate( laddie%buffer%entr_tot         , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%entr_dmin_tot    , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%detr_tot         , n_mem, source = 0._dp)
    call reallocate( laddie%buffer%divQH_sum        , n_mem, source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extend_laddie_buffer

end module laddie_output


