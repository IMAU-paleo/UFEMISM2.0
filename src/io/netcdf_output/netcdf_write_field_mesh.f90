module netcdf_write_field_mesh
  !< Write data to a field in a mesh-based NetCDF file

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use mpi_distributed_memory, only: gather_to_master
  use netcdf_basic

  implicit none

  private

  public :: write_to_field_multopt_mesh_int_2D, write_to_field_multopt_mesh_dp_2D, &
    write_to_field_multopt_mesh_dp_2D_b, write_to_field_multopt_mesh_dp_2D_monthly, &
    write_to_field_multopt_mesh_dp_3D, write_to_field_multopt_mesh_dp_3D_b, &
    write_to_field_multopt_mesh_dp_3D_ocean, write_to_field_multopt_mesh_int_2D_notime, &
    write_to_field_multopt_mesh_int_2D_b_notime, write_to_field_multopt_mesh_int_2D_c_notime, &
    write_to_field_multopt_mesh_dp_2D_notime, write_to_field_multopt_mesh_dp_2D_b_notime, &
    write_to_field_multopt_mesh_dp_2D_c_notime, write_to_field_multopt_mesh_dp_2D_monthly_notime, &
    write_to_field_multopt_mesh_dp_3D_notime, write_to_field_multopt_mesh_dp_3D_b_notime, &
    write_to_field_multopt_mesh_dp_3D_ocean_notime

contains

  subroutine write_to_field_multopt_mesh_int_2D( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    integer,  dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_mesh_int_2D'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    integer,  dimension(:  ), allocatable :: d_tot
    integer,  dimension(:,:), allocatable :: d_tot_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV))
    call gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_tot_with_time( mesh%nV,1))
      d_tot_with_time( :,1) = d_tot
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_int_2D

  subroutine write_to_field_multopt_mesh_dp_2D( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_mesh_dp_2D'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    real(dp), dimension(:  ), allocatable :: d_tot
    real(dp), dimension(:,:), allocatable :: d_tot_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV))
    call gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_tot_with_time( mesh%nV,1))
      d_tot_with_time( :,1) = d_tot
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_2D

  subroutine write_to_field_multopt_mesh_dp_2D_b( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_mesh_dp_2D_b'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    real(dp), dimension(:  ), allocatable :: d_tot
    real(dp), dimension(:,:), allocatable :: d_tot_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nTri))
    call gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_tot_with_time( mesh%nTri,1))
      d_tot_with_time( :,1) = d_tot
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nTri, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_2D_b

  subroutine write_to_field_multopt_mesh_dp_2D_monthly( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D monthly in the physical sense, so a 2-D array!)

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_mesh_dp_2D_monthly'
    integer                                 :: id_var, id_dim_time, ti
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_tot
    real(dp), dimension(:,:,:), allocatable :: d_tot_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV, 12))
    call gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_tot_with_time( mesh%nV,12,1))
      d_tot_with_time( :,:,1) = d_tot
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, 12, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_2D_monthly

  subroutine write_to_field_multopt_mesh_dp_3D( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 3-D in the physical sense, so a 2-D array!)

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_mesh_dp_3D'
    integer                                 :: id_var, id_dim_time, ti
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_tot
    real(dp), dimension(:,:,:), allocatable :: d_tot_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV, mesh%nz))
    call gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_tot_with_time( mesh%nV,mesh%nz,1))
      d_tot_with_time( :,:,1) = d_tot
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, mesh%nz, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_3D

  subroutine write_to_field_multopt_mesh_dp_3D_b( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 3-D in the physical sense, so a 2-D array!)

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_mesh_dp_3D_b'
    integer                                 :: id_var, id_dim_time, ti
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_tot
    real(dp), dimension(:,:,:), allocatable :: d_tot_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nTri, mesh%nz))
    call gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_tot_with_time( mesh%nTri,mesh%nz,1))
      d_tot_with_time( :,:,1) = d_tot
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nTri, mesh%nz, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_3D_b

  subroutine write_to_field_multopt_mesh_dp_3D_ocean( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 3-D ocean data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 3-D in the physical sense, so a 2-D array!)

    ! Write to the last time frame of the variable

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_to_field_multopt_mesh_dp_3D_ocean'
    integer                                 :: id_var, id_dim_time, ti
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_tot
    real(dp), dimension(:,:,:), allocatable :: d_tot_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .true.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV, C%nz_ocean))
    call gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( d_tot_with_time( mesh%nV,C%nz_ocean,1))
      d_tot_with_time( :,:,1) = d_tot
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, C%nz_ocean, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_3D_ocean

  subroutine write_to_field_multopt_mesh_int_2D_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    character(len=*),      intent(in   ) :: filename
    integer,               intent(in   ) :: ncid
    character(len=*),      intent(in   ) :: field_name_options
    integer, dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_field_multopt_mesh_int_2D_notime'
    integer                            :: id_var, id_dim_time, ti
    character(len=1024)                :: var_name
    integer, dimension(:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_int_2D_notime

  subroutine write_to_field_multopt_mesh_int_2D_b_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    character(len=*),      intent(in   ) :: filename
    integer,               intent(in   ) :: ncid
    character(len=*),      intent(in   ) :: field_name_options
    integer, dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_field_multopt_mesh_int_2D_b_notime'
    integer                            :: id_var, id_dim_time, ti
    character(len=1024)                :: var_name
    integer, dimension(:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_int_2D_b( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nTri))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_int_2D_b_notime

  subroutine write_to_field_multopt_mesh_int_2D_c_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    character(len=*),      intent(in   ) :: filename
    integer,               intent(in   ) :: ncid
    character(len=*),      intent(in   ) :: field_name_options
    integer, dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_field_multopt_mesh_int_2D_c_notime'
    integer                            :: id_var, id_dim_time, ti
    character(len=1024)                :: var_name
    integer, dimension(:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_int_2D_c( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nE))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_int_2D_c_notime

  subroutine write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_notime'
    integer                             :: id_var, id_dim_time, ti
    character(len=1024)                 :: var_name
    real(dp), dimension(:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_2D_notime

  subroutine write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_b_notime'
    integer                             :: id_var, id_dim_time, ti
    character(len=1024)                 :: var_name
    real(dp), dimension(:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nTri))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_2D_b_notime

  subroutine write_to_field_multopt_mesh_dp_2D_c_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D in the physical sense, so a 1-D array!)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: field_name_options
    real(dp), dimension(:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_c_notime'
    integer                             :: id_var, id_dim_time, ti
    character(len=1024)                 :: var_name
    real(dp), dimension(:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_2D_c( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nE))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_2D_c_notime

  subroutine write_to_field_multopt_mesh_dp_2D_monthly_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 2-D monthly in the physical sense, so a 2-D array!)

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_mesh_dp_2D_monthly_notime'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV, 12))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_2D_monthly_notime

  subroutine write_to_field_multopt_mesh_dp_3D_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 3-D in the physical sense, so a 2-D array!)

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_mesh_dp_3D_notime'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV, mesh%nz))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_3D_notime

  subroutine write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 3-D in the physical sense, so a 2-D array!)

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_mesh_dp_3D_b_notime'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nTri, mesh%nz))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_3D_b_notime

  subroutine write_to_field_multopt_mesh_dp_3D_ocean_notime( mesh, filename, ncid, &
    field_name_options, d_partial)
    !< Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    !< (Mind you, that's 3-D in the physical sense, so a 2-D array!)

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: field_name_options
    real(dp), dimension(:,:), intent(in   ) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_field_multopt_mesh_dp_3D_ocean_notime'
    integer                               :: id_var, id_dim_time, ti
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if this variable has the correct type and dimensions
    call check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .false.)
#endif

    ! Gather data to the master
    if (par%primary) allocate( d_tot( mesh%nV, C%nz_ocean))
    call gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_mesh_dp_3D_ocean_notime

end module netcdf_write_field_mesh
