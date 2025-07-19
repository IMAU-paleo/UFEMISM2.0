module netcdf_basic_wrappers
  !< UFEMISM wrappers for the most basic NetCDF functionality
  !< (needed to deal with parallelisation and extended error messaging)

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR, MPI_INTEGER
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mpi_basic, only: par
  use basic_model_utilities, only: git_commit_hash
  use netcdf, only: NF90_NOERR, NF90_STRERROR, NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
    NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, NF90_CREATE, NF90_DEF_DIM, NF90_DEF_VAR, &
    NF90_MAX_VAR_DIMS, NF90_PUT_ATT, NF90_OPEN, NF90_NOWRITE, NF90_WRITE, NF90_SHARE, &
    NF90_CLOSE, NF90_GLOBAL, NF90_NETCDF4, NF90_NOCLOBBER

  implicit none

  private

  public :: inquire_dim, inquire_dim_info, inquire_var, inquire_var_info, &
    create_new_netcdf_file_for_writing, create_dimension, create_variable, &
    create_scalar_variable, add_attribute_int, add_attribute_dp, add_attribute_char, &
    open_existing_netcdf_file_for_reading, open_existing_netcdf_file_for_writing, &
    close_netcdf_file, handle_netcdf_error

contains

  ! Inquire dimensions and variables
  subroutine inquire_dim( filename, ncid, dim_name, dim_length, id_dim)
    !< Inquire if this file contains a dimension by name of dim_name.
    !< If so, return its length and identifier; if not, return -1 for both.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: dim_name
    integer,          intent(  out) :: dim_length
    integer,          intent(  out) :: id_dim

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_dim'
    integer                        :: nerr, ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    if (par%primary) then

      ! Check if a dimension of this name exists in the file
      nerr = NF90_INQ_DIMID( ncid, dim_name, id_dim)

      if (nerr /= NF90_NOERR) then
        ! If a dimension by this name does not exist, return -1 for the length and ID
        id_dim     = -1
        dim_length = -1
      else
        ! If a dimension by this name exists, find its length
        call handle_netcdf_error( NF90_inquire_dimension( ncid, id_dim, len = dim_length), &
          filename = filename, dimvarname = dim_name)
      end if

    end if

    call MPI_BCAST( id_dim    , 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( dim_length, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_dim

  subroutine inquire_dim_info( filename, ncid, id_dim, dim_name, dim_length)
    !< Inquire some info of a dimension

    ! In/output variables:
    character(len=*),              intent(in   ) :: filename
    integer,                       intent(in   ) :: ncid
    integer,                       intent(in   ) :: id_dim
    character(len=1024), optional, intent(  out) :: dim_name
    integer,             optional, intent(  out) :: dim_length

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_dim_info'
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    if (par%primary) then
      call handle_netcdf_error( NF90_INQUIRE_DIMENSION( ncid, id_dim, name = dim_name, len = dim_length), &
        filename = filename)
    end if

    if (present( dim_name  )) call MPI_BCAST( dim_name  , len( dim_name), MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    if (present( dim_length)) call MPI_BCAST( dim_length, 1             , MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_dim_info

  subroutine inquire_var( filename, ncid, var_name, id_var)
    !< Inquire if this file contains a variable by name of var_name.
    !< If so, return its identifier. if not, return -1 for the identifier.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: var_name
    integer,          intent(  out) :: id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_var'
    character                      :: dummy1
    integer                        :: nerr, ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! To prevent "unused variable" compiler warnings
    dummy1 = filename( 1:1)

    if (par%primary) then

      ! Check if a variable of this name exists in the file
      nerr = NF90_INQ_VARID( ncid, var_name, id_var)

      if (nerr /= NF90_NOERR) then
        ! If a variable by this name does not exist, return -1 for the ID
        id_var = -1
      end if

    end if

    call MPI_BCAST( id_var, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_var

  subroutine inquire_var_info( filename, ncid, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    !< Inquire some info of a variable

    ! In/output variables:
    character(len=*),                                 intent(in   ) :: filename
    integer,                                          intent(in   ) :: ncid
    integer,                                          intent(in   ) :: id_var
    character(len=1024),                    optional, intent(  out) :: var_name
    integer,                                optional, intent(  out) :: var_type
    integer,                                optional, intent(  out) :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS), optional, intent(  out) :: dims_of_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'inquire_var_info'
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    if (par%primary) then
      ! Inquire some info on this variable
      call handle_netcdf_error( NF90_INQUIRE_VARIABLE( ncid, id_var, name = var_name, &
        xtype = var_type, ndims = ndims_of_var, dimids = dims_of_var), &
        filename = filename)
    end if

    if (present( var_name    )) call MPI_BCAST( var_name    , len( var_name)   , MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    if (present( var_type    )) call MPI_BCAST( var_type    , 1                , MPI_integer, 0, MPI_COMM_WORLD, ierr)
    if (present( ndims_of_var)) call MPI_BCAST( ndims_of_var, 1                , MPI_integer, 0, MPI_COMM_WORLD, ierr)
    if (present(  dims_of_var)) call MPI_BCAST( dims_of_var , NF90_MAX_VAR_DIMS, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine inquire_var_info

  ! Create new NetCDF file
  subroutine create_new_netcdf_file_for_writing( filename, ncid)
    !< Create a new NetCDF file in the specified location for writing.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(  out) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_new_netcdf_file_for_writing'
    logical                        :: file_exists
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if this file already exists
    if (par%primary) then
      inquire( exist = file_exists, file = trim( filename))
      if (file_exists) call crash('file "' // trim( filename) // '" already exists!')
    end if

    ! Create the NetCDF file
    if (par%primary) then
      call handle_netcdf_error( NF90_CREATE( filename, ior( NF90_NOCLOBBER, NF90_NETCDF4), ncid), &
        filename = filename)
    end if
    call MPI_BCAST( ncid, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Add some very basic info about the current simulation to the header
    call add_attribute_char( filename, ncid, NF90_GLOBAL, 'git commit hash', git_commit_hash)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_new_netcdf_file_for_writing

  ! Create dimensions, variables, and attributes
  subroutine create_dimension( filename, ncid, dim_name, dim_length, id_dim)
    !< Create a new dimension in a NetCDF file.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: dim_name
    integer,          intent(in   ) :: dim_length
    integer,          intent(  out) :: id_dim

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_dimension'
    integer                        :: dim_length_present, ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Safety: check if a dimension by this name is already present in this file
    call inquire_dim( filename, ncid, dim_name, dim_length_present, id_dim)
    if (id_dim /= -1) then
      !call crash('file "' // trim( filename) // '" already contains dimension "' // trim( dim_name) // '"!')
      call finalise_routine( routine_name)
      return
    end if

    ! Add the dimension
    if (par%primary) then
      call handle_netcdf_error( NF90_DEF_DIM( ncid, dim_name, dim_length, id_dim), &
        filename = filename, dimvarname = dim_name)
    end if

    call MPI_BCAST( id_dim, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_dimension

  subroutine create_variable( filename, ncid, var_name, var_type, dim_ids, id_var)
    !< Create a new variable in a NetCDF file.

    ! In/output variables:
    character(len=*),      intent(in   ) :: filename
    integer,               intent(in   ) :: ncid
    character(len=*),      intent(in   ) :: var_name
    integer,               intent(in   ) :: var_type
    integer, dimension(:), intent(in   ) :: dim_ids
    integer,               intent(  out) :: id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_variable'
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Safety: check if a variable by this name is already present in this file
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var /= -1) then
      !call crash('file "' // trim( filename) // '" already contains variable "' // trim( var_name) // '"!')
      call finalise_routine( routine_name)
      return
    end if

    ! Add the variable
    if (par%primary) then
      call handle_netcdf_error( NF90_DEF_VAR( ncid, name = var_name, xtype = var_type, &
        dimids = dim_ids, varid = id_var), filename = filename, dimvarname = var_name)
    end if

    call MPI_BCAST( id_var, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_variable

  subroutine create_scalar_variable( filename, ncid, var_name, var_type, id_var)
    !< Create a new scalar variable in a NetCDF file.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: var_name
    integer,          intent(in   ) :: var_type
    integer,          intent(  out) :: id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_scalar_variable'
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Safety: check if a variable by this name is already present in this file
    call inquire_var( filename, ncid, var_name, id_var)
    if (id_var /= -1) call crash('file "' // trim( filename) // '" already contains variable "' // trim( var_name) // '"!')

    ! Add the variable
    if (par%primary) then
      call handle_netcdf_error( NF90_DEF_VAR( ncid, name = var_name, xtype = var_type, varid = id_var), &
        filename = filename, dimvarname = var_name)
    end if

    call MPI_BCAST( id_var, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_scalar_variable

  subroutine add_attribute_int( filename, ncid, id_var, att_name, att_val)
    !< Add an integer-valued attributes to a variable.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    character(len=*), intent(in   ) :: att_name
    integer,          intent(in   ) :: att_val

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_attribute_int'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Add the attribute
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_ATT( ncid, id_var, att_name, att_val), &
        filename = filename)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_attribute_int

  subroutine add_attribute_dp( filename, ncid, id_var, att_name, att_val)
    !< Add a double-precision-valued attributes to a variable.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    character(len=*), intent(in   ) :: att_name
    real(dp),         intent(in   ) :: att_val

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_attribute_dp'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Add the attribute
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_ATT( ncid, id_var, att_name, att_val), &
        filename = filename)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_attribute_dp

  subroutine add_attribute_char( filename, ncid, id_var, att_name, att_val)
    !< Add a character-valued attributes to a variable.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    integer,          intent(in   ) :: id_var
    character(len=*), intent(in   ) :: att_name
    character(len=*), intent(in   ) :: att_val

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_attribute_char'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Add the attribute
    if (par%primary) then
      call handle_netcdf_error( NF90_PUT_ATT( ncid, id_var, att_name, att_val), &
        filename = filename)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_attribute_char

  ! Open and close a NetCDF file
  subroutine open_existing_netcdf_file_for_reading( filename, ncid)
    !< Open the NetCDF file in the specified location for reading only, and return its identifier.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(  out) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'open_netcdf_file_for_reading'
    logical                        :: file_exists
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) call crash('file "' // trim( filename) // '" not found!')

    ! Open the NetCDF file with read-only access
    if (par%primary) then
      call handle_netcdf_error( NF90_OPEN( trim( filename), NF90_NOWRITE, ncid), &
        filename = filename)
    end if

    call MPI_BCAST( ncid, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine open_existing_netcdf_file_for_reading

  subroutine open_existing_netcdf_file_for_writing( filename, ncid)
    !< Open an existing NetCDF file in data mode
    ! In data mode, no new dimensions, variables, or attributes can be created,
    ! but data can be written to existing variables.
    ! When opening an existing NetCDF file, it is by default in data mode.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(  out) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'open_existing_netcdf_file_for_writing'
    logical                        :: file_exists
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) call crash('file "' // trim( filename) // '" not found!')

    ! Open the NetCDF file with read+write access
    if (par%primary) then
      call handle_netcdf_error( NF90_OPEN( trim( filename), ior( NF90_WRITE, NF90_SHARE), ncid), &
        filename = filename)
    end if

    call MPI_BCAST( ncid, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine open_existing_netcdf_file_for_writing

  subroutine close_netcdf_file( ncid)
    !< Close an open NetCDF file

    ! In/output variables:
    integer, intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'close_netcdf_file'

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Close netCDF file:
    if (par%primary) then
      call handle_netcdf_error( NF90_CLOSE( ncid))
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine close_netcdf_file

  ! Parse and print NetCDF error message
  subroutine handle_netcdf_error( nerr, filename, dimvarname)
    !< Print the NetCDF error message and crash the model

    ! In/output variables:
    integer,                    intent(in) :: nerr
    character(len=*), optional, intent(in) :: filename, dimvarname

    ! Local variables:
    character(2048) :: str

    if (nerr == NF90_NOERR) then
      return
    else
      str = trim('NetCDF error "' // trim( NF90_STRERROR( nerr)) // '"')
      if (present( filename  )) str = trim( str) // ' (file "' // trim( filename) // '")'
      if (present( dimvarname)) str = trim( str) // ' (dimension/variable "' // trim( dimvarname) // '")'
      call crash(str)
    end if

  end subroutine handle_netcdf_error

end module netcdf_basic_wrappers