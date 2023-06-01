MODULE netcdf_resource_tracking

! This module contains routines for creating and writing to the
! resource tracking file, which stores the amount of computation
! time used by all the different subroutines.
!
! (See src/basic/control_resources_and_error_messaging.f90 for
! details on how the resource tracking system works)

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string, resource_tracker
  USE model_configuration                                    , ONLY: C

  USE netcdf,        ONLY: NF90_UNLIMITED, NF90_DOUBLE, NF90_INT
  USE netcdf_basic,  ONLY: nerr, create_new_netcdf_file_for_writing, create_dimension, create_variable, add_attribute_char, close_netcdf_file, &
                           open_existing_netcdf_file_for_writing, find_timeframe, write_var_master_int_2D, write_var_master_dp_2D, inquire_var_multopt
  USE netcdf_output, ONLY: write_time_to_file

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  CHARACTER(LEN=256) :: filename_resource_tracker

CONTAINS

  SUBROUTINE write_to_resource_tracking_file( time)
    ! Write to the resource tracking output file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp),                           INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'write_to_resource_tracking_file'
    INTEGER                                           :: n_routines, length_routine_name, i
    INTEGER,  DIMENSION(:    ), ALLOCATABLE           :: routine_name_encoded
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE           :: routine_names_encoded
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE           :: tcomp
    INTEGER                                           :: ncid
    INTEGER                                           :: ti
    INTEGER                                           :: id_var_names, id_var_tcomp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Number of routines in the resource tracker
    n_routines          = SIZE( resource_tracker)
    length_routine_name = LEN( resource_tracker( 1)%routine_path)

    ! Allocate memory
    ALLOCATE( routine_name_encoded(             length_routine_name))
    ALLOCATE( routine_names_encoded( n_routines,length_routine_name))
    ALLOCATE( tcomp(                 n_routines,1                  ))

    ! Encode subroutine names and gather computation times
    DO i = 1, n_routines
      CALL encode_subroutine_path_as_integer( resource_tracker( i)%routine_path, routine_name_encoded)
      routine_names_encoded( i,:) = routine_name_encoded
      tcomp( i,1) = resource_tracker( i)%tcomp
    END DO

    ! Open the file for writing
    CALL open_existing_netcdf_file_for_writing( filename_resource_tracker, ncid)

    ! Append new timeframe to time variable
    CALL write_time_to_file( filename_resource_tracker, ncid, time)

    ! Find new timeframe
    CALL find_timeframe( filename_resource_tracker, ncid, time, ti)

    ! Find var_ids of variables to write to
    CALL inquire_var_multopt( filename_resource_tracker, ncid, 'routine_names_encoded', id_var_names)
    CALL inquire_var_multopt( filename_resource_tracker, ncid, 'tcomp'                , id_var_tcomp)

    ! Write data
    CALL write_var_master_int_2D( filename_resource_tracker, ncid, id_var_names, routine_names_encoded)

    ! Computation time
    CALL write_var_master_dp_2D(  filename_resource_tracker, ncid, id_var_tcomp, tcomp, start = [1,ti], count = [n_routines,1])

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    DEALLOCATE( routine_name_encoded )
    DEALLOCATE( routine_names_encoded)
    DEALLOCATE( tcomp                )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_resource_tracking_file

  SUBROUTINE create_resource_tracking_file
    ! Create the resource tracking output file

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'create_resource_tracking_file'
    INTEGER                                           :: n_routines, length_routine_name
    INTEGER                                           :: ncid
    INTEGER                                           :: id_dim_n_routines, id_dim_time, id_dim_name_length
    INTEGER                                           :: id_var_time, id_var_names, id_var_tcomp

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Number of routines in the resource tracker
    n_routines          = SIZE( resource_tracker)
    length_routine_name = LEN( resource_tracker( 1)%routine_path)

    ! Determine the file name
    filename_resource_tracker = TRIM( C%output_dir) // '/resource_tracking.nc'

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename_resource_tracker, ncid)

    ! Define the dimensions
    CALL create_dimension( filename_resource_tracker, ncid, 'n_routines' , n_routines         , id_dim_n_routines )
    CALL create_dimension( filename_resource_tracker, ncid, 'name_length', length_routine_name, id_dim_name_length)
    CALL create_dimension( filename_resource_tracker, ncid, 'time'       , NF90_UNLIMITED     , id_dim_time       )

    ! Define variables
    ! ================

    ! Subroutine names
    CALL create_variable(    filename_resource_tracker, ncid, 'routine_names_encoded', NF90_INT, [id_dim_n_routines, id_dim_name_length], id_var_names)
    CALL add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'long_name', 'Encoded subroutine names')
    CALL add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'encoding' , 'See UFEMISM code (src/netcdf/netcdf_resource_tracking/encode_subroutine_path_as_integer)')

    ! Time
    CALL create_variable(    filename_resource_tracker, ncid, 'time', NF90_DOUBLE, [id_dim_time], id_var_time)
    CALL add_attribute_char( filename_resource_tracker, ncid, id_var_time, 'long_name', 'time')
    CALL add_attribute_char( filename_resource_tracker, ncid, id_var_time, 'units'    , 'years')

    ! Computation times
    CALL create_variable(    filename_resource_tracker, ncid, 'tcomp', NF90_DOUBLE, [id_dim_n_routines, id_dim_time], id_var_tcomp)
    CALL add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'long_name', 'Computation times per subroutine')
    CALL add_attribute_char( filename_resource_tracker, ncid, id_var_names, 'units'    , 's')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_resource_tracking_file

  SUBROUTINE encode_subroutine_path_as_integer( subroutine_path, path_int_enc)
    ! Encode the current subroutine path as an integer array so it can be saved as a NetCDF variable
    !
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

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                         INTENT(IN)    :: subroutine_path
    INTEGER, DIMENSION(LEN(subroutine_path)), INTENT(OUT)   :: path_int_enc

    ! Local variables:
    INTEGER                                           :: i

    path_int_enc = 0

    DO i = 1, LEN(subroutine_path)

      SELECT CASE ( subroutine_path( i:i))
      CASE( ' ')
        path_int_enc( i) = -1
      CASE( '0')
        path_int_enc( i) = 0
      CASE( '1')
        path_int_enc( i) = 1
      CASE( '2')
        path_int_enc( i) = 2
      CASE( '3')
        path_int_enc( i) = 3
      CASE( '4')
        path_int_enc( i) = 4
      CASE( '5')
        path_int_enc( i) = 5
      CASE( '6')
        path_int_enc( i) = 6
      CASE( '7')
        path_int_enc( i) = 7
      CASE( '8')
        path_int_enc( i) = 8
      CASE( '9')
        path_int_enc( i) = 9
      CASE( 'a')
        path_int_enc( i) = 11
      CASE( 'b')
        path_int_enc( i) = 12
      CASE( 'c')
        path_int_enc( i) = 13
      CASE( 'd')
        path_int_enc( i) = 14
      CASE( 'e')
        path_int_enc( i) = 15
      CASE( 'f')
        path_int_enc( i) = 16
      CASE( 'g')
        path_int_enc( i) = 17
      CASE( 'h')
        path_int_enc( i) = 18
      CASE( 'i')
        path_int_enc( i) = 19
      CASE( 'j')
        path_int_enc( i) = 20
      CASE( 'k')
        path_int_enc( i) = 21
      CASE( 'l')
        path_int_enc( i) = 22
      CASE( 'm')
        path_int_enc( i) = 23
      CASE( 'n')
        path_int_enc( i) = 24
      CASE( 'o')
        path_int_enc( i) = 25
      CASE( 'p')
        path_int_enc( i) = 26
      CASE( 'q')
        path_int_enc( i) = 27
      CASE( 'r')
        path_int_enc( i) = 28
      CASE( 's')
        path_int_enc( i) = 29
      CASE( 't')
        path_int_enc( i) = 30
      CASE( 'u')
        path_int_enc( i) = 31
      CASE( 'v')
        path_int_enc( i) = 32
      CASE( 'w')
        path_int_enc( i) = 33
      CASE( 'x')
        path_int_enc( i) = 34
      CASE( 'y')
        path_int_enc( i) = 35
      CASE( 'z')
        path_int_enc( i) = 36
      CASE( 'A')
        path_int_enc( i) = 37
      CASE( 'B')
        path_int_enc( i) = 38
      CASE( 'C')
        path_int_enc( i) = 39
      CASE( 'D')
        path_int_enc( i) = 40
      CASE( 'E')
        path_int_enc( i) = 41
      CASE( 'F')
        path_int_enc( i) = 42
      CASE( 'G')
        path_int_enc( i) = 43
      CASE( 'H')
        path_int_enc( i) = 44
      CASE( 'I')
        path_int_enc( i) = 45
      CASE( 'J')
        path_int_enc( i) = 46
      CASE( 'K')
        path_int_enc( i) = 47
      CASE( 'L')
        path_int_enc( i) = 48
      CASE( 'M')
        path_int_enc( i) = 49
      CASE( 'N')
        path_int_enc( i) = 50
      CASE( 'O')
        path_int_enc( i) = 51
      CASE( 'P')
        path_int_enc( i) = 52
      CASE( 'Q')
        path_int_enc( i) = 53
      CASE( 'R')
        path_int_enc( i) = 54
      CASE( 'S')
        path_int_enc( i) = 55
      CASE( 'T')
        path_int_enc( i) = 56
      CASE( 'U')
        path_int_enc( i) = 57
      CASE( 'V')
        path_int_enc( i) = 58
      CASE( 'W')
        path_int_enc( i) = 59
      CASE( 'X')
        path_int_enc( i) = 60
      CASE( 'Y')
        path_int_enc( i) = 61
      CASE( 'Z')
        path_int_enc( i) = 62
      CASE( '_')
        path_int_enc( i) = 63
      CASE( '/')
        path_int_enc( i) = 64
      CASE( '(')
        path_int_enc( i) = 65
      CASE( ')')
        path_int_enc( i) = 66
      CASE DEFAULT
        CALL crash('unknown character in routine_path "' // TRIM( subroutine_path) // '"!')
      END SELECT

    END DO

  END SUBROUTINE encode_subroutine_path_as_integer

END MODULE netcdf_resource_tracking