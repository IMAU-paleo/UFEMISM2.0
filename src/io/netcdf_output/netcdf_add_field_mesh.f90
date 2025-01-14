module netcdf_add_field_mesh
  !< Add a data field to a mesh-based NetCDF file

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use netcdf_basic
  use netcdf, only: NF90_INT, NF90_DOUBLE

  implicit none

  private

  public :: add_field_mesh_int_2D, add_field_mesh_dp_2D, add_field_mesh_dp_2D_b, &
  add_field_mesh_dp_2D_monthly, add_field_mesh_dp_3D, add_field_mesh_dp_3D_b, &
  add_field_mesh_dp_3D_ocean, add_field_mesh_int_2D_notime, add_field_mesh_int_2D_b_notime, &
  add_field_mesh_int_2D_c_notime, add_field_mesh_dp_2D_notime, add_field_mesh_dp_2D_b_notime, &
  add_field_mesh_dp_2D_c_notime, add_field_mesh_dp_2D_monthly_notime, add_field_mesh_dp_3D_notime, &
  add_field_mesh_dp_3D_b_notime, add_field_mesh_dp_3D_ocean_notime

contains

  subroutine add_field_mesh_int_2D( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_int_2D'
    integer                        :: id_dim_vi, id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    call inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    if (id_dim_vi   == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_time == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_vi, id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_int_2D

  subroutine add_field_mesh_dp_2D( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_2D'
    integer                        :: id_dim_vi, id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    call inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    if (id_dim_vi   == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_time == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_2D

  subroutine add_field_mesh_dp_2D_b( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_2D_b'
    integer                        :: id_dim_ti, id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    call inquire_dim_multopt( filename, ncid, field_name_options_time    , id_dim_time)

    ! Safety
    if (id_dim_ti   == -1) call crash('no ti dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_time == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti, id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_2D_b

  subroutine add_field_mesh_dp_2D_monthly( filename, ncid, var_name, long_name, units)
    !< Add a 2-D monthly variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_2D_monthly'
    integer                        :: id_dim_vi, id_dim_month, id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_month(           filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_month , id_dim_month)
    call inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time )

    ! Safety
    if (id_dim_vi    == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_month == -1) call crash('no month dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_time  == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_month, id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_2D_monthly

  subroutine add_field_mesh_dp_3D( filename, ncid, var_name, long_name, units)
    !< Add a 3-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_3D'
    integer                        :: id_dim_vi, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_zeta(            filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta  , id_dim_zeta)
    call inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    if (id_dim_vi   == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_zeta == -1) call crash('no zeta dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_time == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_3D

  subroutine add_field_mesh_dp_3D_b( filename, ncid, var_name, long_name, units)
    !< Add a 3-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_3D_b'
    integer                        :: id_dim_ti, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_zeta(            filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta    , id_dim_zeta)
    call inquire_dim_multopt( filename, ncid, field_name_options_time    , id_dim_time)

    ! Safety
    if (id_dim_ti   == -1) call crash('no ti dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_zeta == -1) call crash('no zeta dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_time == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_3D_b

  subroutine add_field_mesh_dp_3D_ocean( filename, ncid, var_name, long_name, units)
    !< Add a 3-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_3D_ocean'
    integer                        :: id_dim_vi, id_dim_depth, id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_depth(           filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth , id_dim_depth)
    call inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time )

    ! Safety
    if (id_dim_vi    == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_depth == -1) call crash('no depth dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_time  == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_depth, id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_3D_ocean

  subroutine add_field_mesh_int_2D_notime( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_int_2D_notime'
    integer                        :: id_dim_vi, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )

    ! Safety
    if (id_dim_vi   == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_vi /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_int_2D_notime

  subroutine add_field_mesh_int_2D_b_notime( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_int_2D_b_notime'
    integer                        :: id_dim_ti, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )

    ! Safety
    if (id_dim_ti   == -1) call crash('no ti dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_ti /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_int_2D_b( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_int_2D_b_notime

  subroutine add_field_mesh_int_2D_c_notime( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_int_2D_c_notime'
    integer                        :: id_dim_ei, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei  )

    ! Safety
    if (id_dim_ei   == -1) call crash('no ei dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_ei /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_int_2D_c( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_int_2D_c_notime

  subroutine add_field_mesh_dp_2D_notime( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_2D_notime'
    integer                        :: id_dim_vi, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )

    ! Safety
    if (id_dim_vi   == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_2D_notime

  subroutine add_field_mesh_dp_2D_b_notime( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_2D_b_notime'
    integer                        :: id_dim_ti, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )

    ! Safety
    if (id_dim_ti   == -1) call crash('no ti dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_2D_b_notime

  subroutine add_field_mesh_dp_2D_c_notime( filename, ncid, var_name, long_name, units)
    !< Add a 2-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_2D_c_notime'
    integer                        :: id_dim_ei, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei  )

    ! Safety
    if (id_dim_ei   == -1) call crash('no ei dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ei /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_2D_c( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_2D_c_notime

  subroutine add_field_mesh_dp_2D_monthly_notime( filename, ncid, var_name, long_name, units)
    !< Add a 2-D monthly variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_2D_monthly_notime'
    integer                        :: id_dim_vi, id_dim_month, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_month(           filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_month , id_dim_month)

    ! Safety
    if (id_dim_vi    == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_month == -1) call crash('no month dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_month /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_2D_monthly_notime

  subroutine add_field_mesh_dp_3D_notime( filename, ncid, var_name, long_name, units)
    !< Add a 3-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_3D_notime'
    integer                        :: id_dim_vi, id_dim_zeta, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_zeta(            filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta  , id_dim_zeta)

    ! Safety
    if (id_dim_vi   == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_zeta == -1) call crash('no zeta dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_zeta /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_3D_notime

  subroutine add_field_mesh_dp_3D_b_notime( filename, ncid, var_name, long_name, units)
    !< Add a 3-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_3D_b_notime'
    integer                        :: id_dim_ti, id_dim_zeta, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_zeta(            filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta    , id_dim_zeta)

    ! Safety
    if (id_dim_ti   == -1) call crash('no ti dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_zeta == -1) call crash('no zeta dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti, id_dim_zeta /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_3D_b_notime

  subroutine add_field_mesh_dp_3D_ocean_notime( filename, ncid, var_name, long_name, units)
    !< Add a 3-D variable to an existing NetCDF file with a mesh

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_mesh_dp_3D_ocean_notime'
    integer                        :: id_dim_vi, id_dim_depth, id_var

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check if all mesh dimensions and variables are there
    call check_depth(           filename, ncid)
    call check_mesh_dimensions( filename, ncid)
#endif

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth , id_dim_depth)

    ! Safety
    if (id_dim_vi    == -1) call crash('no vi dimension could be found in file "' // trim( filename) // '"!')
    if (id_dim_depth == -1) call crash('no depth dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_depth /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

#if (DO_ASSERTIONS)
    ! Final safety check
    call check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_mesh_dp_3D_ocean_notime

end module netcdf_add_field_mesh
