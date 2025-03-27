module netcdf_write_field_transect
  !< Write data to a field in a transect NetCDF file

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use transect_types, only: type_transect
  use mpi_distributed_memory, only: gather_to_primary
  use netcdf_basic

  implicit none

  private

  public :: write_to_field_multopt_transect_dp_2D, write_to_field_multopt_transect_dp_3D


contains

  ! Write to fields with a time dimension

  subroutine write_to_field_multopt_transect_dp_2D( transect, filename, ncid, field_name_options, d_transect_partial)
    !< Write a 2-D data field defined on a transect to a NetCDF file variable

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_transect),    intent(in   ) :: transect
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_transect_partial

    ! ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_transect_dp_2D'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    real(dp), dimension(:  ), allocatable :: d_transect
    real(dp), dimension(:,:), allocatable :: d_transect_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

    ! Gather data to the primary
    if (par%primary) allocate( d_transect( transect%nV))
    call gather_to_primary( d_transect_partial, d_transect)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_transect_with_time( transect%nV,1))
      d_transect_with_time( :,1) = d_transect
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_transect_with_time, start = (/ 1, ti /), count = (/ transect%nV, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_transect_dp_2D

  subroutine write_to_field_multopt_transect_dp_3D( transect, filename, ncid, field_name_options, d_transect_partial)
    !< Write a 3-D data field defined on a transect to a NetCDF file variable

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_transect),      intent(in   ) :: transect
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_transect_partial

    ! ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_transect_dp_3D'
    integer                                 :: id_var, id_dim_time, ti
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_transect
    real(dp), dimension(:,:,:), allocatable :: d_transect_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

    ! Gather data to the primary
    if (par%primary) allocate( d_transect( transect%nV, size( d_transect_partial,2)))
    call gather_to_primary( d_transect_partial, d_transect)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_transect_with_time( transect%nV, size( d_transect_partial,2), 1))
      d_transect_with_time( :,:,1) = d_transect
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_transect_with_time, start = (/ 1, 1, ti /), count = (/ transect%nV, size( d_transect_partial,2), 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_transect_dp_3D

end module netcdf_write_field_transect
