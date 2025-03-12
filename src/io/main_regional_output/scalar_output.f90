module scalar_output

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use region_types, only: type_model_region
  use netcdf_io_main

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
    call write_to_field_multopt_dp_0D(  region%output_filename_scalar, ncid, 'dt_ice',     region%ice%dt_ice)
    call write_to_field_multopt_int_0D( region%output_filename_scalar, ncid, 'n_visc_its', region%ice%n_visc_its)
    call write_to_field_multopt_int_0D( region%output_filename_scalar, ncid, 'n_Axb_its',  region%ice%n_Axb_its)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_scalar_regional_output_file

  subroutine create_scalar_regional_output_file( region)
    !< Create the scalar regional output NetCDF file

    ! In/output variables:
    type(type_model_region)                            , intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_scalar_regional_output_file'
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
    filename_base = trim( C%output_dir) // 'scalar_output_' // region%name
    call generate_filename_XXXXXdotnc( filename_base, region%output_filename_scalar)

    ! Print to terminal
    if (par%master) write(0,'(A)') '   Creating scalar output file "' // colour_string( trim( region%output_filename_scalar), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( region%output_filename_scalar, ncid)

    ! Add time, zeta, and month dimensions+variables to the file
    call add_time_dimension_to_file(  region%output_filename_scalar, ncid)

    ! Add the default data fields to the file

    ! Integrated ice geometry
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_area')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume_af')

    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_area_PD')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume_PD')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume_af_PD')

    ! Integrated mass fluxes
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_total')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_gr')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_fl')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_land')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_ocean')

    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_total')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_gr')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_fl')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_land')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_ocean')

    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'LMB_total')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'LMB_gr')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'LMB_fl')

    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_total')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_gr')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_fl')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_land')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_ocean')

    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'gl_flux')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'cf_gr_flux')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'cf_fl_flux')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'margin_land_flux')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'margin_ocean_flux')

    ! Numerical stability info
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'dt_ice')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'n_visc_its')
    call create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'n_Axb_its')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_scalar_regional_output_file

  subroutine create_scalar_regional_output_file_field( filename, ncid, choice_output_field)
    !< Add a single field to the scalar regional output NetCDF file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter                                     :: routine_name = 'create_scalar_regional_output_file_field'

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Add the specified data field to the file
    select case (choice_output_field)
      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')
      case ('none')
        ! Do nothing

      ! Total ice sheet area
      case ('ice_area')
        call add_field_dp_0D( filename, ncid, 'ice_area', long_name = 'Total ice area', units = 'm^2')

        ! Total ice sheet volume in metres of sea level equivalent
      case ('ice_volume')
        call add_field_dp_0D( filename, ncid, 'ice_volume', long_name = 'Total ice volume', units = 'm s.l.e.')

      ! Total ice sheet volume above floatation in metres of sea level equivalent
      case ('ice_volume_af')
        call add_field_dp_0D( filename, ncid, 'ice_volume_af', long_name = 'Total ice volume above floatation', units = 'm s.l.e.')

      ! Total ice sheet area for present-day
      case ('ice_area_PD')
        call add_field_dp_0D( filename, ncid, 'ice_area_PD', long_name = 'Total ice area for present-day', units = 'm^2')

        ! Total ice sheet volume in metres of sea level equivalent for present-day
      case ('ice_volume_PD')
        call add_field_dp_0D( filename, ncid, 'ice_volume_PD', long_name = 'Total ice volume for present-day', units = 'm s.l.e.')

      ! Total ice sheet volume above floatation in metres of sea level equivalent for present-day
      case ('ice_volume_af_PD')
        call add_field_dp_0D( filename, ncid, 'ice_volume_af_PD', long_name = 'Total ice volume above floatation for present-day', units = 'm s.l.e.')

      ! Total SMB integrated over the entire domain
      case ('SMB_total')
        call add_field_dp_0D( filename, ncid, 'SMB_total', long_name = 'Area-integrated total SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice sheet area
      case ('SMB_gr')
        call add_field_dp_0D( filename, ncid, 'SMB_gr', long_name = 'Area-integrated ice sheet SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice shelf area
      case ('SMB_fl')
        call add_field_dp_0D( filename, ncid, 'SMB_fl', long_name = 'Area-integrated ice shelf SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice-free land area
      case ('SMB_land')
        call add_field_dp_0D( filename, ncid, 'SMB_land', long_name = 'Area-integrated ice-free land SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice-free ocean area
      case ('SMB_ocean')
        call add_field_dp_0D( filename, ncid, 'SMB_ocean', long_name = 'Area-integrated ice-free ocean SMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire domain
      case ('BMB_total')
        call add_field_dp_0D( filename, ncid, 'BMB_total', long_name = 'Area-integrated total BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice sheet area
      case ('BMB_gr')
        call add_field_dp_0D( filename, ncid, 'BMB_gr', long_name = 'Area-integrated ice sheet BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice shelf area
      case ('BMB_fl')
        call add_field_dp_0D( filename, ncid, 'BMB_fl', long_name = 'Area-integrated ice shelf BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice-free land area
      case ('BMB_land')
        call add_field_dp_0D( filename, ncid, 'BMB_land', long_name = 'Area-integrated ice-free land BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice-free ocean area
      case ('BMB_ocean')
        call add_field_dp_0D( filename, ncid, 'BMB_ocean', long_name = 'Area-integrated ice-free ocean BMB', units = 'Gt yr^-1')

      ! Total LMB integrated over the entire domain
      case ('LMB_total')
        call add_field_dp_0D( filename, ncid, 'LMB_total', long_name = 'Area-integrated total LMB', units = 'Gt yr^-1')

      ! Total LMB integrated over the entire ice sheet area
      case ('LMB_gr')
        call add_field_dp_0D( filename, ncid, 'LMB_gr', long_name = 'Area-integrated ice sheet LMB', units = 'Gt yr^-1')

      ! Total LMB integrated over the entire ice shelf area
      case ('LMB_fl')
        call add_field_dp_0D( filename, ncid, 'LMB_fl', long_name = 'Area-integrated ice shelf LMB', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire domain
      case ('AMB_total')
        call add_field_dp_0D( filename, ncid, 'AMB_total', long_name = 'Area-integrated total additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice sheet area
      case ('AMB_gr')
        call add_field_dp_0D( filename, ncid, 'AMB_gr', long_name = 'Area-integrated ice sheet additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice shelf area
      case ('AMB_fl')
        call add_field_dp_0D( filename, ncid, 'AMB_fl', long_name = 'Area-integrated ice shelf additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice-free land area
      case ('AMB_land')
        call add_field_dp_0D( filename, ncid, 'AMB_land', long_name = 'Area-integrated ice-free land additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice-free ocean area
      case ('AMB_ocean')
        call add_field_dp_0D( filename, ncid, 'AMB_ocean', long_name = 'Area-integrated ice-free ocean additional MB from other sources', units = 'Gt yr^-1')

      ! Total flux through the grounding line
      case ('gl_flux')
        call add_field_dp_0D( filename, ncid, 'gl_flux', long_name = 'Total lateral grounding line flux', units = 'Gt yr^-1')

      ! Total flux through grounded calving fronts
      case ('cf_gr_flux')
        call add_field_dp_0D( filename, ncid, 'cf_gr_flux', long_name = 'Total lateral grounded calving front flux', units = 'Gt yr^-1')

      ! Total flux through floating calving fronts
      case ('cf_fl_flux')
        call add_field_dp_0D( filename, ncid, 'cf_fl_flux', long_name = 'Total lateral floating calving front flux', units = 'Gt yr^-1')

      ! Total flux exiting ice margins into grounded areas
      case ('margin_land_flux')
        call add_field_dp_0D( filename, ncid, 'margin_land_flux', long_name = 'Total lateral flux exiting the ice margin into ground', units = 'Gt yr^-1')

      ! Total flux exiting ice margins into marine areas
      case ('margin_ocean_flux')
        call add_field_dp_0D( filename, ncid, 'margin_ocean_flux', long_name = 'Total lateral flux exiting the ice margin into water', units = 'Gt yr^-1')

      ! Ice-dynamical time step
      case ('dt_ice')
        call add_field_dp_0D( filename, ncid, 'dt_ice', long_name = 'Ice-dynamical time step', units = 'yr')

      ! Number of non-linear viscosity iterations
      case ('n_visc_its')
        call add_field_int_0D( filename, ncid, 'n_visc_its', long_name = 'Number of non-linear viscosity iterations')

      ! Number of iterations in iterative solver for linearised momentum balance
      case ('n_Axb_its')
        call add_field_int_0D( filename, ncid, 'n_Axb_its', long_name = 'Number of iterations in iterative solver for linearised momentum balance')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_scalar_regional_output_file_field

end module scalar_output
