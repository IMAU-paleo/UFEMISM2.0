module netcdf_check_fields

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use netcdf_field_name_options
  use netcdf_check_dimensions
  use netcdf_basic_wrappers
  use netcdf, only: NF90_MAX_VAR_DIMS, NF90_DOUBLE, NF90_FLOAT, NF90_INT

  implicit none

  private

  public :: check_xy_grid_field_int_2D, check_xy_grid_field_int_3D, check_xy_grid_field_dp_2D, &
    check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, check_xy_grid_field_dp_3D_ocean, &
    check_lat_grid_field_dp_1D_monthly, check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, & 
    check_lonlat_grid_field_dp_2D_monthly, &
    check_lonlat_grid_field_dp_3D, check_lonlat_grid_field_dp_3D_ocean, check_mesh_field_int_2D, &
    check_mesh_field_int_2D_b, check_mesh_field_int_2D_c, check_mesh_field_dp_2D, &
    check_mesh_field_dp_2D_b, check_mesh_field_dp_2D_c, check_mesh_field_dp_2D_monthly, &
    check_mesh_field_dp_3D, check_mesh_field_dp_3D_b, check_mesh_field_dp_3D_ocean

contains

  ! x/y-grid field variables
  subroutine check_xy_grid_field_int_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_int_2D'
    integer                                :: id_dim_x, id_dim_y, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x( filename, ncid)
    call check_y( filename, ncid)

    ! inquire x,y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has x,y as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_int_2D

  subroutine check_xy_grid_field_int_3D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_int_3D'
    integer                                :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x(    filename, ncid)
    call check_y(    filename, ncid)
    call check_zeta( filename, ncid)

    ! inquire x,y,zeta dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    call inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has x,y,zeta as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y,zeta as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y,zeta as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_int_3D

  subroutine check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_2D'
    integer                                :: id_dim_x, id_dim_y, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x( filename, ncid)
    call check_y( filename, ncid)

    ! inquire x,y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has x,y as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_2D

  subroutine check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_2D_monthly'
    integer                                :: id_dim_x, id_dim_y, id_dim_month, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x(     filename, ncid)
    call check_y(     filename, ncid)
    call check_month( filename, ncid)

    ! inquire x,y dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x    )
    call inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y    )
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')
    if (.not. any( dims_of_var == id_dim_month)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have month as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has x,y,m as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y,m as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y,m as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_2D_monthly

  subroutine check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_3D'
    integer                                :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x(    filename, ncid)
    call check_y(    filename, ncid)
    call check_zeta( filename, ncid)

    ! inquire x,y,zeta dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    call inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has x,y,zeta as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y,zeta as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y,zeta as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_3D

  subroutine check_xy_grid_field_dp_3D_ocean(       filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_xy_grid_field_dp_3D_ocean'
    integer                                :: id_dim_x, id_dim_y, id_dim_depth, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid x and y dimensions and variables
    call check_x(     filename, ncid)
    call check_y(     filename, ncid)
    call check_depth( filename, ncid)

    ! inquire x,y,depth dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x   )
    call inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y   )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check x,y dimensions
    if (.not. any( dims_of_var == id_dim_x    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have x as a dimension!')
    if (.not. any( dims_of_var == id_dim_y    )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have y as a dimension!')
    if (.not. any( dims_of_var == id_dim_depth)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have depth as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has x,y,depth as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have x,y,depth as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have x,y,depth as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_xy_grid_field_dp_3D_ocean

  ! lon/lat-grid field variables
  subroutine check_lat_grid_field_dp_1D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lat_grid_field_dp_1D_monthly'
    integer                                :: id_dim_lat, id_dim_month, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lat dimension and variables
    call check_lat(   filename, ncid)
    call check_month( filename, ncid)

    ! inquire lon,lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat  )
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lat  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_month)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have month as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has lat,m as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('no-longitude variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat,m as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat,m as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lat_grid_field_dp_1D_monthly

  subroutine check_lonlat_grid_field_int_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_int_2D'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon( filename, ncid)
    call check_lat( filename, ncid)

    ! inquire lon,lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has lon,lat as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_int_2D

  subroutine check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_2D'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon( filename, ncid)
    call check_lat( filename, ncid)

    ! inquire lon,lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon, id_dim_lon)
    call inquire_dim_multopt( filename, ncid, field_name_options_lat, id_dim_lat)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has lon,lat as dimensions.
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_2D

  subroutine check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_2D_monthly'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_month, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon(   filename, ncid)
    call check_lat(   filename, ncid)
    call check_month( filename, ncid)

    ! inquire lon,lat dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon  )
    call inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat  )
    call inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_month)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have month as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has lon,lat,m as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat,m as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat,m as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_2D_monthly

  subroutine check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_3D'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon(  filename, ncid)
    call check_lat(  filename, ncid)
    call check_zeta( filename, ncid)

    ! inquire lon,lat,zeta dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    call inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has lon,lat,zeta as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat,zeta as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat,zeta as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_3D

  subroutine check_lonlat_grid_field_dp_3D_ocean( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_lonlat_grid_field_dp_3D_ocean'
    integer                                :: id_dim_lon, id_dim_lat, id_dim_depth, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid lon and lat dimensions and variables
    call check_lon(   filename, ncid)
    call check_lat(   filename, ncid)
    call check_depth( filename, ncid)

    ! inquire lon,lat,depth dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon )
    call inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth, id_dim_depth)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check lon,lat dimensions
    if (.not. any( dims_of_var == id_dim_lon  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have longitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_lat  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have latitude as a dimension!')
    if (.not. any( dims_of_var == id_dim_depth)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have depth as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 3) then
          ! The variable only has lon,lat,depth as dimensions.
        else
          if (ndims_of_var == 4) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has four dimensions, but the fourth one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have lon,lat,depth as dimensions
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 4) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have lon,lat,depth as dimensions

        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_lonlat_grid_field_dp_3D_ocean

  ! mesh field variables
  subroutine check_mesh_field_int_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_int_2D'
    integer                                :: id_dim_vi, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_int_2D

  subroutine check_mesh_field_int_2D_b( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_int_2D_b'
    integer                                :: id_dim_ti, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ti)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_int_2D_b

  subroutine check_mesh_field_int_2D_c( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_int_2D_c'
    integer                                :: id_dim_ei, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. var_type == NF90_INT) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_INT!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ei)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ei as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_int_2D_c

  subroutine check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D'
    integer                                :: id_dim_vi, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D

  subroutine check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D_b'
    integer                                :: id_dim_ti, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ti)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D_b

  subroutine check_mesh_field_dp_2D_c( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D mesh variable by this name
    !
    ! NOTE: this is 2-D in the physical sense, so a 1-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D_c'
    integer                                :: id_dim_ei, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ei)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ei as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 1) then
          ! The variable only has vi as a dimension
        else
          if (ndims_of_var == 2) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has two dimensions, but the second one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi as a dimension
        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi as a dimension

        if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D_c

  subroutine check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly mesh variable by this name
    !
    ! NOTE: this is 2-D monthly in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_2D_monthly'
    integer                                :: id_dim_vi, id_dim_month, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_month(           filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_month , id_dim_month)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')
    if (.not. any( dims_of_var == id_dim_month)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have month as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has vi,m as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi,m as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi,m as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_2D_monthly

  subroutine check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_3D'
    integer                                :: id_dim_vi, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_zeta(            filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta  , id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has vi,zeta as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi,zeta as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi,zeta as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_3D

  subroutine check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_3D'
    integer                                :: id_dim_ti, id_dim_zeta, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_zeta(            filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    call inquire_dim_multopt( filename, ncid, field_name_options_zeta    , id_dim_zeta)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_ti  )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have ti as a dimension!')
    if (.not. any( dims_of_var == id_dim_zeta)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have zeta as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has ti,zeta as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have ti,zeta as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have ti,zeta as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_3D_b

  subroutine check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time)
    ! Check if this file contains a 3-D mesh variable by this name
    !
    ! NOTE: this is 3-D in the physical sense, so a 2-D array!

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    logical,          optional, intent(in   ) :: should_have_time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'check_mesh_field_dp_3D_ocean'
    integer                                :: id_dim_vi, id_dim_depth, id_dim_time, id_var
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var
    logical                                :: file_has_time

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file has valid mesh dimensions and variables
    call check_mesh_dimensions( filename, ncid)
    call check_depth(           filename, ncid)

    ! inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    call inquire_dim_multopt( filename, ncid, field_name_options_depth , id_dim_depth)

    ! inquire variable
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var == -1) call crash('variable "' // trim( var_name) // '" could not be found in file "' // trim( filename) // '"!')

    ! inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    if (.not. (var_type == NF90_FLOAT .or. var_type == NF90_DOUBLE)) then
      call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    end if

    ! Check mesh dimensions
    if (.not. any( dims_of_var == id_dim_vi   )) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have vi as a dimension!')
    if (.not. any( dims_of_var == id_dim_depth)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have depth as a dimension!')

    if (.not. present( should_have_time)) then
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
      if (id_dim_time == -1) then
        file_has_time = .false.
      else
        file_has_time = .true.
      end if

      if (file_has_time) then
        ! Check if the variable has time as a dimension
        if (ndims_of_var == 2) then
          ! The variable only has vi,depth as dimensions
        else
          if (ndims_of_var == 3) then
            if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' &
              // trim( filename) // '" has three dimensions, but the third one is not time!')
          else
            call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          end if
        end if
      else ! if (file_has_time) then
        ! The file does not have a time dimension; the variable should only have vi,depth as dimensions
        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      end if ! if (file_has_time) then

    else ! if (.not. present( should_have_time)) then
      if (should_have_time) then
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        call check_time( filename, ncid)

        ! inquire the time dimension
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        if (ndims_of_var /= 3) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

      else ! if (should_have_time) then
        ! This variable should not have a time dimension; the variable should only have vi,depth as dimensions

        if (ndims_of_var /= 2) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      end if ! if (should_have_time) then
    end if ! if (.not. present( should_have_time)) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_mesh_field_dp_3D_ocean

end module netcdf_check_fields
