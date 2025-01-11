module netcdf_inquire_dimensions
  !< Inquire if an existing NetCDF file contains dimensions and variables
  !< for certain standard dimensions (zeta, time, month)

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use netcdf_field_name_options

  implicit none

  private

  public :: inquire_zeta, inquire_month, inquire_time

contains

  subroutine inquire_zeta( filename, ncid, has_zeta)
    !< Inquire if a NetCDF file contains a zeta dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    logical,          intent(  out) :: has_zeta

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_zeta'
    integer                        :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Look for zeta dimension and variable
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)
    call inquire_var_multopt( filename, ncid, field_name_options_zeta, id_var_zeta)

    ! Check if everything is there
    has_zeta = ( &
      id_dim_zeta /= -1 .and. &
      id_var_zeta /= -1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_zeta

  subroutine inquire_month( filename, ncid, has_month)
    !< Inquire if a NetCDF file contains a month dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    logical,          intent(  out) :: has_month

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_month'
    integer                        :: id_dim_month, id_var_month

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Look for month dimension and variable
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)
    call inquire_var_multopt( filename, ncid, field_name_options_month, id_var_month)

    ! Check if everything is there
    has_month = ( &
      id_dim_month /= -1 .and. &
      id_var_month /= -1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_month

  subroutine inquire_time( filename, ncid, has_time)
    !< Inquire if a NetCDF file contains a time dimension and variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    logical,          intent(  out) :: has_time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_time'
    integer                        :: id_dim_time, id_var_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Look for time dimension and variable
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! Check if everything is there
    has_time = ( &
      id_dim_time /= -1 .and. &
      id_var_time /= -1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_time

end module netcdf_inquire_dimensions