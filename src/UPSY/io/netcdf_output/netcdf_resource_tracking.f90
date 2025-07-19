module netcdf_resource_tracking

! This module contains routines for creating and writing to the
! resource tracking file, which stores the amount of computation
! time used by all the different subroutines.
!
! (See src/basic/control_resources_and_error_messaging.f90 for
! details on how the resource tracking system works)

  use precisions, only: dp
  use control_resources_and_error_messaging, only: resource_tracker, init_routine, finalise_routine, crash
  use netcdf_basic
  use netcdf_add_write_scalar_variables
  use netcdf, only: NF90_DOUBLE, NF90_INT, NF90_UNLIMITED

  implicit none

  private

  public :: create_resource_tracking_file, write_to_resource_tracking_file

  character(len=256) :: filename_resource_tracker

contains

  subroutine write_to_resource_tracking_file( time)

    ! Input variables:
    real(dp), intent(in   ) :: time

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'write_to_resource_tracking_file'
    integer                               :: n_routines, length_routine_name, i
    integer,  dimension(:  ), allocatable :: routine_name_encoded
    integer,  dimension(:,:), allocatable :: routine_names_encoded
    real(dp), dimension(:,:), allocatable :: tcomp
    integer                               :: ncid
    integer                               :: ti
    integer                               :: id_var_names, id_var_tcomp

    ! Add routine to path
    call init_routine( routine_name)

    ! Number of routines in the resource tracker
    n_routines          = size( resource_tracker)
    length_routine_name = len( resource_tracker( 1)%routine_path)

    ! allocate memory
    allocate( routine_name_encoded(             length_routine_name))
    allocate( routine_names_encoded( n_routines,length_routine_name))
    allocate( tcomp(                 n_routines,1                  ))

    ! Encode subroutine names and gather computation times
    do i = 1, n_routines
      call encode_subroutine_path_as_integer( resource_tracker( i)%routine_path, routine_name_encoded)
      routine_names_encoded( i,:) = routine_name_encoded
      tcomp( i,1) = resource_tracker( i)%tcomp
    end do

    ! Open the file for writing
    call open_existing_netcdf_file_for_writing( filename_resource_tracker, ncid)

    ! Append new timeframe to time variable
    call write_time_to_file( filename_resource_tracker, ncid, time)

    ! Find new timeframe
    call find_timeframe( filename_resource_tracker, ncid, time, ti)

    ! Find var_ids of variables to write to
    call inquire_var_multopt( filename_resource_tracker, ncid, 'routine_names_encoded', id_var_names)
    call inquire_var_multopt( filename_resource_tracker, ncid, 'tcomp'                , id_var_tcomp)

    ! Write data
    call write_var_primary( filename_resource_tracker, ncid, id_var_names, routine_names_encoded)

    ! Computation time
    call write_var_primary(  filename_resource_tracker, ncid, id_var_tcomp, tcomp, start = [1,ti], count = [n_routines,1])

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

#endif

  end subroutine write_to_resource_tracking_file

  subroutine create_resource_tracking_file( output_dir)

    ! In/output variables:
    character(len=*), intent(in) :: output_dir

    ! Local variables:
    character(len=256), parameter :: routine_name = 'create_resource_tracking_file'
    integer                       :: n_routines, length_routine_name
    integer                       :: ncid
    integer                       :: id_dim_n_routines, id_dim_time, id_dim_name_length
    integer                       :: id_var_time, id_var_names, id_var_tcomp

#if (DO_RESOURCE_TRACKING)

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Number of routines in the resource tracker
    n_routines          = size( resource_tracker)
    length_routine_name = len( resource_tracker( 1)%routine_path)

    ! Determine the file name
    filename_resource_tracker = TRIM( output_dir) // '/resource_tracking.nc'

    ! Create a new NetCDF file
    call create_new_netcdf_file_for_writing( filename_resource_tracker, ncid)

    ! Define the dimensions
    call create_dimension( filename_resource_tracker, ncid, 'n_routines' , n_routines         , id_dim_n_routines )
    call create_dimension( filename_resource_tracker, ncid, 'name_length', length_routine_name, id_dim_name_length)
    call create_dimension( filename_resource_tracker, ncid, 'time'       , NF90_UNLIMITED     , id_dim_time       )

    ! Define variables
    ! ================

    ! subroutine names
    call create_variable(    filename_resource_tracker, ncid, 'routine_names_encoded', NF90_INT, [id_dim_n_routines, id_dim_name_length], id_var_names)
    call add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'long_name', 'Encoded subroutine names')
    call add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'encoding' , 'See UFEMISM code (src/netcdf/netcdf_resource_tracking/encode_subroutine_path_as_integer)')

    ! Time
    call create_variable(    filename_resource_tracker, ncid, 'time', NF90_DOUBLE, [id_dim_time], id_var_time)
    call add_attribute_char( filename_resource_tracker, ncid, id_var_time, 'long_name', 'time')
    call add_attribute_char( filename_resource_tracker, ncid, id_var_time, 'units'    , 'years')

    ! Computation times
    call create_variable(    filename_resource_tracker, ncid, 'tcomp', NF90_DOUBLE, [id_dim_n_routines, id_dim_time], id_var_tcomp)
    call add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'long_name', 'Computation times per subroutine')
    call add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'units'    , 's')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

#endif

  end subroutine create_resource_tracking_file

  subroutine encode_subroutine_path_as_integer( subroutine_path, path_int_enc)
    !< Encode the current subroutine path as an integer array so it can be saved as a NetCDF variable

    ! Use the simplest possible encoding:
    !
    !  ' ' = -1 (empty character)
    !
    !    0 = 0
    !    1 = 1
    !    ...
    !    9 = 9
    !
    !    a = 10
    !    b = 11
    !    c = 12
    !    ...
    !    z = 36
    !
    !    A = 37
    !    B = 38
    !    C = 39
    !    ...
    !    Z = 62
    !
    !    _ = 63 (underscore)
    !    / = 64 (forward slash)
    !    ( = 65 (left  bracket)
    !    ) = 66 (right bracket)

    ! In/output variables:
    character(len=*),                         intent(in   ) :: subroutine_path
    integer, dimension(len(subroutine_path)), intent(  out) :: path_int_enc

    ! Local variables:
    integer :: i

    path_int_enc = 0

    do i = 1, len(subroutine_path)

      select case ( subroutine_path( i:i))
      case default
        call crash('unknown character in routine_path "' // trim( subroutine_path) // '"!')
      case (' ')
        path_int_enc( i) = -1
      case ('0')
        path_int_enc( i) = 0
      case ('1')
        path_int_enc( i) = 1
      case ('2')
        path_int_enc( i) = 2
      case ('3')
        path_int_enc( i) = 3
      case ('4')
        path_int_enc( i) = 4
      case ('5')
        path_int_enc( i) = 5
      case ('6')
        path_int_enc( i) = 6
      case ('7')
        path_int_enc( i) = 7
      case ('8')
        path_int_enc( i) = 8
      case ('9')
        path_int_enc( i) = 9
      case ('a')
        path_int_enc( i) = 11
      case ('b')
        path_int_enc( i) = 12
      case ('c')
        path_int_enc( i) = 13
      case ('d')
        path_int_enc( i) = 14
      case ('e')
        path_int_enc( i) = 15
      case ('f')
        path_int_enc( i) = 16
      case ('g')
        path_int_enc( i) = 17
      case ('h')
        path_int_enc( i) = 18
      case ('i')
        path_int_enc( i) = 19
      case ('j')
        path_int_enc( i) = 20
      case ('k')
        path_int_enc( i) = 21
      case ('l')
        path_int_enc( i) = 22
      case ('m')
        path_int_enc( i) = 23
      case ('n')
        path_int_enc( i) = 24
      case ('o')
        path_int_enc( i) = 25
      case ('p')
        path_int_enc( i) = 26
      case ('q')
        path_int_enc( i) = 27
      case ('r')
        path_int_enc( i) = 28
      case ('s')
        path_int_enc( i) = 29
      case ('t')
        path_int_enc( i) = 30
      case ('u')
        path_int_enc( i) = 31
      case ('v')
        path_int_enc( i) = 32
      case ('w')
        path_int_enc( i) = 33
      case ('x')
        path_int_enc( i) = 34
      case ('y')
        path_int_enc( i) = 35
      case ('z')
        path_int_enc( i) = 36
      case ('A')
        path_int_enc( i) = 37
      case ('B')
        path_int_enc( i) = 38
      case ('C')
        path_int_enc( i) = 39
      case ('D')
        path_int_enc( i) = 40
      case ('E')
        path_int_enc( i) = 41
      case ('F')
        path_int_enc( i) = 42
      case ('G')
        path_int_enc( i) = 43
      case ('H')
        path_int_enc( i) = 44
      case ('I')
        path_int_enc( i) = 45
      case ('J')
        path_int_enc( i) = 46
      case ('K')
        path_int_enc( i) = 47
      case ('L')
        path_int_enc( i) = 48
      case ('M')
        path_int_enc( i) = 49
      case ('N')
        path_int_enc( i) = 50
      case ('O')
        path_int_enc( i) = 51
      case ('P')
        path_int_enc( i) = 52
      case ('Q')
        path_int_enc( i) = 53
      case ('R')
        path_int_enc( i) = 54
      case ('S')
        path_int_enc( i) = 55
      case ('T')
        path_int_enc( i) = 56
      case ('U')
        path_int_enc( i) = 57
      case ('V')
        path_int_enc( i) = 58
      case ('W')
        path_int_enc( i) = 59
      case ('X')
        path_int_enc( i) = 60
      case ('Y')
        path_int_enc( i) = 61
      case ('Z')
        path_int_enc( i) = 62
      case ('_')
        path_int_enc( i) = 63
      case ('/')
        path_int_enc( i) = 64
      case ('(')
        path_int_enc( i) = 65
      case (')')
        path_int_enc( i) = 66
      end select

    end do

  end subroutine encode_subroutine_path_as_integer

end module netcdf_resource_tracking