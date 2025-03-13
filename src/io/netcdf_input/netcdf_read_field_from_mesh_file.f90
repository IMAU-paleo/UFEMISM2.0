module netcdf_read_field_from_mesh_file
  !< Read data fields from a mesh file

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mpi_distributed_memory, only: distribute_from_primary
  use mesh_types, only: type_mesh
  use mesh_memory, only: deallocate_mesh
  use netcdf_basic
  use netcdf_setup_grid_mesh_from_file

  implicit none

  private

  public :: read_field_from_mesh_file_int_2D, read_field_from_mesh_file_int_2D_b, &
    read_field_from_mesh_file_dp_2D, read_field_from_mesh_file_dp_2D_b, &
    read_field_from_mesh_file_dp_2D_monthly, read_field_from_mesh_file_dp_3D, &
    read_field_from_mesh_file_dp_3D_b, read_field_from_mesh_file_dp_3D_ocean, &
    read_field_from_mesh_file_CDF, read_field_from_mesh_file_CDF_b

contains

  subroutine read_field_from_mesh_file_int_2D( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    integer,  dimension(:),           intent(  out) :: d_mesh_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_int_2D'
    integer                               :: ncid
    type(type_mesh)                       :: mesh_loc
    integer                               :: id_var
    character(len=1024)                   :: var_name
    integer , dimension(:  ), allocatable :: d_mesh
    integer , dimension(:,:), allocatable :: d_mesh_with_time
    integer                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nV))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nV, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, ti /), count = (/ mesh_loc%nV, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_int_2D

  subroutine read_field_from_mesh_file_int_2D_b( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    integer , dimension(:),           intent(  out) :: d_mesh_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_int_2D_b'
    integer                               :: ncid
    type(type_mesh)                       :: mesh_loc
    integer                               :: id_var
    character(len=1024)                   :: var_name
    integer , dimension(:  ), allocatable :: d_mesh
    integer , dimension(:,:), allocatable :: d_mesh_with_time
    integer                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_int_2D_b( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nTri))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nTri, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, ti /), count = (/ mesh_loc%nTri, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_int_2D_b

  subroutine read_field_from_mesh_file_dp_2D( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    real(dp), dimension(:),           intent(  out) :: d_mesh_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_dp_2D'
    integer                               :: ncid
    type(type_mesh)                       :: mesh_loc
    integer                               :: id_var
    character(len=1024)                   :: var_name
    real(dp), dimension(:  ), allocatable :: d_mesh
    real(dp), dimension(:,:), allocatable :: d_mesh_with_time
    integer                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nV))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nV, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, ti /), count = (/ mesh_loc%nV, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_dp_2D

  subroutine read_field_from_mesh_file_dp_2D_b( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D data field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    character(len=*),                 intent(in   ) :: field_name_options
    real(dp), dimension(:),           intent(  out) :: d_mesh_partial
    real(dp),               optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_dp_2D_b'
    integer                               :: ncid
    type(type_mesh)                       :: mesh_loc
    integer                               :: id_var
    character(len=1024)                   :: var_name
    real(dp), dimension(:  ), allocatable :: d_mesh
    real(dp), dimension(:,:), allocatable :: d_mesh_with_time
    integer                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nTri))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nTri, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, ti /), count = (/ mesh_loc%nTri, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_dp_2D_b

  subroutine read_field_from_mesh_file_dp_2D_monthly( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 2-D monthly data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_dp_2D_monthly'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the file has a valid month dimension
    call check_month( filename, ncid)

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nV, 12))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nV, 12, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, 12, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_dp_2D_monthly

  subroutine read_field_from_mesh_file_dp_3D( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 3-D data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_dp_3D'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: nzeta_loc
    real(dp), dimension(:), allocatable     :: zeta_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Set up the vertical coordinate zeta from the file
    call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nV, nzeta_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nV, nzeta_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, nzeta_loc, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_dp_3D

  subroutine read_field_from_mesh_file_dp_3D_b( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 3-D data field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_dp_3D_b'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: nzeta_loc
    real(dp), dimension(:), allocatable     :: zeta_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Set up the vertical coordinate zeta from the file
    call setup_zeta_from_file( filename, ncid, nzeta_loc, zeta_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nTri, nzeta_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nTri, nzeta_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nTri, nzeta_loc, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_dp_3D_b

  subroutine read_field_from_mesh_file_dp_3D_ocean( filename, field_name_options, &
    d_mesh_partial, time_to_read)
    !< Read a 3-D ocean data field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial
    real(dp),                 optional, intent(in   ) :: time_to_read

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_dp_3D_ocean'
    integer                                 :: ncid
    type(type_mesh)                         :: mesh_loc
    integer                                 :: ndepth_loc
    real(dp), dimension(:), allocatable     :: depth_loc
    integer                                 :: id_var
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh
    real(dp), dimension(:,:,:), allocatable :: d_mesh_with_time
    integer                                 :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the file
    call setup_mesh_from_file( filename, ncid, mesh_loc)

    ! Set up the vertical coordinate depth from the file
    call setup_depth_from_file( filename, ncid, ndepth_loc, depth_loc)

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! Check if the variable has the required dimensions
    call check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = present( time_to_read))

    ! allocate memory
    if (par%primary) allocate( d_mesh( mesh_loc%nV, ndepth_loc))

    ! Read data from file
    if (.not. present( time_to_read)) then
      call read_var_primary( filename, ncid, id_var, d_mesh)
    else
      ! allocate memory
      if (par%primary) allocate( d_mesh_with_time( mesh_loc%nV, ndepth_loc, 1))
      ! Find out which timeframe to read
      call find_timeframe( filename, ncid, time_to_read, ti)
      ! Read data
      call read_var_primary( filename, ncid, id_var, d_mesh_with_time, start = (/ 1, 1, ti /), count = (/ mesh_loc%nV, ndepth_loc, 1 /) )
      ! Copy to output memory
      if (par%primary) d_mesh = d_mesh_with_time( :,:,1)
      ! Clean up after yourself
      if (par%primary) deallocate( d_mesh_with_time)
    end if

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)
    call deallocate_mesh( mesh_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_dp_3D_ocean

  subroutine read_field_from_mesh_file_CDF( filename, field_name_options, &
    d_mesh_partial)
    !< Read a cumulative density function field from a NetCDF file on a mesh

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_field_from_mesh_file_CDF'
    integer                               :: ncid
    integer                               :: id_dim_vi, id_var, id_dim_bins, nbins_loc, nV_loc
    character(len=1024)                   :: var_name
    real(dp), dimension(:,:), allocatable :: d_mesh

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Get number of mesh vertices
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi, dim_length = nV_loc)

    ! Check that number of bins in file match the ones in the config file
    call inquire_dim_multopt( filename, ncid, 'bin', id_dim_bins, dim_length = nbins_loc)
    if (nbins_loc /= C%subgrid_bedrock_cdf_nbins) then
      call crash('number of CDF bins in external file ({int_01}) does not match subgrid_bedrock_cdf_nbins in config file!', int_01 = nbins_loc)
    end if

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! allocate memory
    if (par%primary) allocate( d_mesh( nV_loc, C%subgrid_bedrock_cdf_nbins))

    ! Read data from file
    call read_var_primary( filename, ncid, id_var, d_mesh)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_CDF

  subroutine read_field_from_mesh_file_CDF_b( filename, field_name_options, &
    d_mesh_partial)
    !< Read a cumulative density function field from a NetCDF file on a mesh b-grid

    ! NOTE: the mesh should be read before, and memory allocated for d_mesh_partial!

    ! In/output variables:
    character(len=*),                   intent(in   ) :: filename
    character(len=*),                   intent(in   ) :: field_name_options
    real(dp), dimension(:,:),           intent(  out) :: d_mesh_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_field_from_mesh_file_CDF_b'
    integer                                 :: ncid
    integer                                 :: id_dim_ti, id_var, id_dim_bins, nbins_loc, nTri_loc
    character(len=1024)                     :: var_name
    real(dp), dimension(:,:  ), allocatable :: d_mesh

    ! Add routine to path
    call init_routine( routine_name)

    ! == Read grid and data from file
    ! ===============================

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Get number of mesh triangles
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti, dim_length = nTri_loc)

    ! Check that number of bins in file match the ones in the config file
    call inquire_dim_multopt( filename, ncid, 'bin', id_dim_bins, dim_length = nbins_loc)
    if (nbins_loc /= C%subgrid_bedrock_cdf_nbins) then
      call crash('number of CDF bins in external file ({int_01}) does not match subgrid_bedrock_cdf_nbins in config file!', int_01 = nbins_loc)
    end if

    ! Look for the specified variable in the file
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('couldnt find any of the options "' // trim( field_name_options) // '" in file "' // trim( filename)  // '"!')

    ! allocate memory
    if (par%primary) allocate( d_mesh( nTri_loc, C%subgrid_bedrock_cdf_nbins))

    ! Read data from file
    call read_var_primary( filename, ncid, id_var, d_mesh)

    ! Close the NetCDF file
    call close_netcdf_file( ncid)

    ! == Distribute gridded data from the primary to all processes in partial vector form
    ! ==================================================================================

    ! Distribute data
    call distribute_from_primary( d_mesh, d_mesh_partial)

    ! Clean up after yourself
    if (par%primary) deallocate( d_mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_field_from_mesh_file_CDF_b

end module netcdf_read_field_from_mesh_file
