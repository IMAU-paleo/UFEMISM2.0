module scalar_output_files

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use region_types, only: type_model_region
  use netcdf_io_main
  use reallocate_mod

  implicit none

  private

  public :: create_scalar_regional_output_file, write_to_scalar_regional_output_file

contains

  subroutine write_to_scalar_regional_output_file( region)
    !< Write to the scalar regional output NetCDF file

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_scalar_regional_output_file'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( region%output_filename_scalar, ncid)

    ! write the time to the file
    call write_time_to_file( region%output_filename_scalar, ncid, region%time)

    ! write the default data fields to the file
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_area',          region%scalars%ice_area)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume',        region%scalars%ice_volume)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume_af',     region%scalars%ice_volume_af)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_area_PD',       region%scalars%ice_area_PD)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume_PD',     region%scalars%ice_volume_PD)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume_af_PD',  region%scalars%ice_volume_af_PD)

    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_total',         region%scalars%SMB_total)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_gr',            region%scalars%SMB_gr)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_fl',            region%scalars%SMB_fl)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_land',          region%scalars%SMB_land)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_ocean',         region%scalars%SMB_ocean)

    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_total',         region%scalars%BMB_total)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_gr',            region%scalars%BMB_gr)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_fl',            region%scalars%BMB_fl)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_land',          region%scalars%BMB_land)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_ocean',         region%scalars%BMB_ocean)

    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'LMB_total',         region%scalars%LMB_total)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'LMB_gr',            region%scalars%LMB_gr)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'LMB_fl',            region%scalars%LMB_fl)

    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_total',         region%scalars%AMB_total)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_gr',            region%scalars%AMB_gr)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_fl',            region%scalars%AMB_fl)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_land',          region%scalars%AMB_land)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_ocean',         region%scalars%AMB_ocean)

    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'gl_flux',           region%scalars%gl_flux)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'cf_gr_flux',        region%scalars%cf_gr_flux)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'cf_fl_flux',        region%scalars%cf_fl_flux)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'margin_land_flux',  region%scalars%margin_land_flux)
    call write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'margin_ocean_flux', region%scalars%margin_ocean_flux)

    ! Numerical stability info
    call write_to_field_multopt_dp_0D(  region%output_filename_scalar, ncid, 'dt_ice',           region%ice%dt_ice)
    call write_to_field_multopt_int_0D( region%output_filename_scalar, ncid, 'n_visc_its',       region%ice%n_visc_its)
    call write_to_field_multopt_int_0D( region%output_filename_scalar, ncid, 'n_Axb_its',        region%ice%n_Axb_its)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_scalar_regional_output_file

  subroutine create_scalar_regional_output_file( region)
    !< Create the scalar regional output NetCDF file

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_scalar_regional_output_file'
    character(len=1024)            :: filename_base, filename
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = trim( C%output_dir) // 'scalar_output_' // region%name
    call generate_filename_XXXXXdotnc( filename_base, filename)
    region%output_filename_scalar = filename

    ! Print to terminal
    if (par%master) write(0,'(A)') '   Creating scalar output file "' // colour_string( trim( filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Add time, zeta, and month dimensions+variables to the file
    call add_time_dimension_to_file( filename, ncid)

    ! Add the default data fields to the file

    ! Integrated ice geometry
    call add_field_dp_0D( filename, ncid, 'ice_area',          long_name = 'Total ice area', units = 'm^2')
    call add_field_dp_0D( filename, ncid, 'ice_volume',        long_name = 'Total ice volume', units = 'm s.l.e.')
    call add_field_dp_0D( filename, ncid, 'ice_volume_af',     long_name = 'Total ice volume above floatation', units = 'm s.l.e.')

    call add_field_dp_0D( filename, ncid, 'ice_area_PD',       long_name = 'Total ice area for present-day', units = 'm^2')
    call add_field_dp_0D( filename, ncid, 'ice_volume_PD',     long_name = 'Total ice volume for present-day', units = 'm s.l.e.')
    call add_field_dp_0D( filename, ncid, 'ice_volume_af_PD',  long_name = 'Total ice volume above floatation for present-day', units = 'm s.l.e.')

    ! Integrated mass fluxes
    call add_field_dp_0D( filename, ncid, 'SMB_total',         long_name = 'Area-integrated total SMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'SMB_gr',            long_name = 'Area-integrated ice sheet SMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'SMB_fl',            long_name = 'Area-integrated ice shelf SMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'SMB_land',          long_name = 'Area-integrated ice-free land SMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'SMB_ocean',         long_name = 'Area-integrated ice-free ocean SMB', units = 'Gt yr^-1')

    call add_field_dp_0D( filename, ncid, 'BMB_total',         long_name = 'Area-integrated total BMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'BMB_gr',            long_name = 'Area-integrated ice sheet BMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'BMB_fl',            long_name = 'Area-integrated ice shelf BMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'BMB_land',          long_name = 'Area-integrated ice-free land BMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'BMB_ocean',         long_name = 'Area-integrated ice-free ocean BMB', units = 'Gt yr^-1')

    call add_field_dp_0D( filename, ncid, 'LMB_total',         long_name = 'Area-integrated total LMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'LMB_gr',            long_name = 'Area-integrated ice sheet LMB', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'LMB_fl',            long_name = 'Area-integrated ice shelf LMB', units = 'Gt yr^-1')

    call add_field_dp_0D( filename, ncid, 'AMB_total',         long_name = 'Area-integrated total additional MB from other sources', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'AMB_gr',            long_name = 'Area-integrated ice sheet additional MB from other sources', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'AMB_fl',            long_name = 'Area-integrated ice shelf additional MB from other sources', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'AMB_land',          long_name = 'Area-integrated ice-free land additional MB from other sources', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'AMB_ocean',         long_name = 'Area-integrated ice-free ocean additional MB from other sources', units = 'Gt yr^-1')

    call add_field_dp_0D( filename, ncid, 'gl_flux',           long_name = 'Total lateral grounding line flux', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'cf_gr_flux',        long_name = 'Total lateral grounded calving front flux', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'cf_fl_flux',        long_name = 'Total lateral floating calving front flux', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'margin_land_flux',  long_name = 'Total lateral flux exiting the ice margin into ground', units = 'Gt yr^-1')
    call add_field_dp_0D( filename, ncid, 'margin_ocean_flux', long_name = 'Total lateral flux exiting the ice margin into water', units = 'Gt yr^-1')

    ! Numerical stability info
    call add_field_dp_0D(  filename, ncid, 'dt_ice',           long_name = 'Ice-dynamical time step', units = 'yr')
    call add_field_int_0D( filename, ncid, 'n_visc_its',       long_name = 'Number of non-linear viscosity iterations')
    call add_field_int_0D( filename, ncid, 'n_Axb_its',        long_name = 'Number of iterations in iterative solver for linearised momentum balance')

    ! Allocate memory to buffer scalar output data between output writing intervals
    call allocate_scalar_output_buffer( region)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_scalar_regional_output_file

  subroutine allocate_scalar_output_buffer( region)
    !< Allocate memory to buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_scalar_output_buffer'
    integer                        :: n_mem

    ! Add routine to path
    call init_routine( routine_name)

    ! Only allocate memory for this on the master
    if (par%master) then

      n_mem = 1000
      region%scalars%buffer%n_mem = n_mem
      region%scalars%buffer%n     = 0

      allocate( region%scalars%buffer%time             ( n_mem), source = 0._dp)

      allocate( region%scalars%buffer%ice_area         ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%ice_volume       ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%ice_volume_af    ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%ice_area_PD      ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%ice_volume_PD    ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%ice_volume_af_PD ( n_mem), source = 0._dp)

      allocate( region%scalars%buffer%SMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%SMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%SMB_fl           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%SMB_land         ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%SMB_ocean        ( n_mem), source = 0._dp)

      allocate( region%scalars%buffer%BMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%BMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%BMB_fl           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%BMB_land         ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%BMB_ocean        ( n_mem), source = 0._dp)

      allocate( region%scalars%buffer%LMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%LMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%LMB_fl           ( n_mem), source = 0._dp)

      allocate( region%scalars%buffer%AMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%AMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%AMB_fl           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%AMB_land         ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%AMB_ocean        ( n_mem), source = 0._dp)

      allocate( region%scalars%buffer%gl_flux          ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%cf_gr_flux       ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%cf_fl_flux       ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%margin_land_flux ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%margin_ocean_flux( n_mem), source = 0._dp)

      allocate( region%scalars%buffer%dt_ice           ( n_mem), source = 0._dp)
      allocate( region%scalars%buffer%n_visc_its       ( n_mem), source = 0)
      allocate( region%scalars%buffer%n_Axb_its        ( n_mem), source = 0)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_scalar_output_buffer

  subroutine buffer_scalar_output( region)
    !< Buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'buffer_scalar_output'
    integer                        :: n

    ! Add routine to path
    call init_routine( routine_name)

    ! Only the master does this
    if (par%master) then

      ! Increase timeframe count
      region%scalars%buffer%n = region%scalars%buffer%n + 1
      n = region%scalars%buffer%n

      ! Extend buffer memory if necessary
      if (n > region%scalars%buffer%n_mem - 10) call extend_scalar_output_buffer( region)

      ! Store new timeframe in buffer
      region%scalars%buffer%time             ( n) = region%time

      region%scalars%buffer%ice_area         ( n) = region%scalars%ice_area
      region%scalars%buffer%ice_volume       ( n) = region%scalars%ice_volume
      region%scalars%buffer%ice_volume_af    ( n) = region%scalars%ice_volume_af
      region%scalars%buffer%ice_area_PD      ( n) = region%scalars%ice_area_PD
      region%scalars%buffer%ice_volume_PD    ( n) = region%scalars%ice_volume_PD
      region%scalars%buffer%ice_volume_af_PD ( n) = region%scalars%ice_volume_af_PD

      region%scalars%buffer%SMB_total        ( n) = region%scalars%SMB_total
      region%scalars%buffer%SMB_gr           ( n) = region%scalars%SMB_gr
      region%scalars%buffer%SMB_fl           ( n) = region%scalars%SMB_fl
      region%scalars%buffer%SMB_land         ( n) = region%scalars%SMB_land
      region%scalars%buffer%SMB_ocean        ( n) = region%scalars%SMB_ocean

      region%scalars%buffer%BMB_total        ( n) = region%scalars%BMB_total
      region%scalars%buffer%BMB_gr           ( n) = region%scalars%BMB_gr
      region%scalars%buffer%BMB_fl           ( n) = region%scalars%BMB_fl
      region%scalars%buffer%BMB_land         ( n) = region%scalars%BMB_land
      region%scalars%buffer%BMB_ocean        ( n) = region%scalars%BMB_ocean

      region%scalars%buffer%LMB_total        ( n) = region%scalars%LMB_total
      region%scalars%buffer%LMB_gr           ( n) = region%scalars%LMB_gr
      region%scalars%buffer%LMB_fl           ( n) = region%scalars%LMB_fl

      region%scalars%buffer%AMB_total        ( n) = region%scalars%AMB_total
      region%scalars%buffer%AMB_gr           ( n) = region%scalars%AMB_gr
      region%scalars%buffer%AMB_fl           ( n) = region%scalars%AMB_fl
      region%scalars%buffer%AMB_land         ( n) = region%scalars%AMB_land
      region%scalars%buffer%AMB_ocean        ( n) = region%scalars%AMB_ocean

      region%scalars%buffer%gl_flux          ( n) = region%scalars%gl_flux
      region%scalars%buffer%cf_gr_flux       ( n) = region%scalars%cf_gr_flux
      region%scalars%buffer%cf_fl_flux       ( n) = region%scalars%cf_fl_flux
      region%scalars%buffer%margin_land_flux ( n) = region%scalars%margin_land_flux
      region%scalars%buffer%margin_ocean_flux( n) = region%scalars%margin_ocean_flux

      region%scalars%buffer%dt_ice           ( n) = region%ice%dt_ice
      region%scalars%buffer%n_visc_its       ( n) = region%ice%n_visc_its
      region%scalars%buffer%n_Axb_its        ( n) = region%ice%n_Axb_its
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine buffer_scalar_output

  subroutine extend_scalar_output_buffer( region)
    !< Extend memory to buffer the scalar output data between output writing intervals
    !
    ! NOTE: should only be called by the master!

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'extend_scalar_output_buffer'
    integer                        :: n_mem

    ! Add routine to path
    call init_routine( routine_name)

    n_mem = region%scalars%buffer%n_mem * 2
    region%scalars%buffer%n_mem = n_mem

    call reallocate( region%scalars%buffer%time             , n_mem, source = 0._dp)

    call reallocate( region%scalars%buffer%ice_area         , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%ice_volume       , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%ice_volume_af    , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%ice_area_PD      , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%ice_volume_PD    , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%ice_volume_af_PD , n_mem, source = 0._dp)

    call reallocate( region%scalars%buffer%SMB_total        , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%SMB_gr           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%SMB_fl           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%SMB_land         , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%SMB_ocean        , n_mem, source = 0._dp)

    call reallocate( region%scalars%buffer%BMB_total        , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%BMB_gr           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%BMB_fl           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%BMB_land         , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%BMB_ocean        , n_mem, source = 0._dp)

    call reallocate( region%scalars%buffer%LMB_total        , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%LMB_gr           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%LMB_fl           , n_mem, source = 0._dp)

    call reallocate( region%scalars%buffer%AMB_total        , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%AMB_gr           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%AMB_fl           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%AMB_land         , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%AMB_ocean        , n_mem, source = 0._dp)

    call reallocate( region%scalars%buffer%gl_flux          , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%cf_gr_flux       , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%cf_fl_flux       , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%margin_land_flux , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%margin_ocean_flux, n_mem, source = 0._dp)

    call reallocate( region%scalars%buffer%dt_ice           , n_mem, source = 0._dp)
    call reallocate( region%scalars%buffer%n_visc_its       , n_mem, source = 0)
    call reallocate( region%scalars%buffer%n_Axb_its        , n_mem, source = 0)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extend_scalar_output_buffer

end module scalar_output_files
