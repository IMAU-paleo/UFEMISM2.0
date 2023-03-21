MODULE config_file_utilities

  ! Some subroutines to facilitate working with config files

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine, &
                                                                     remove_leading_spaces, capitalise_string

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines ======
! ========================

  SUBROUTINE check_config_file_validity( config_filename, namelist_filename)
    ! Check if the config file "config_filename" is valid
    !
    ! Do this by reading one line at a time of the config file, determining the name of the variable
    ! declared in that line, and checking if that variable also exists in the namelist file.
    !
    ! A namelist file is created earlier by writing the namelist to a text file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename
    CHARACTER(LEN=*),                INTENT(IN)        :: namelist_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_config_file_validity'
    INTEGER, PARAMETER                                 :: config_unit   = 1337
    INTEGER, PARAMETER                                 :: namelist_unit = 1338
    INTEGER                                            :: ios
    LOGICAL                                            :: found_end_of_file_config, found_end_of_file_namelist
    CHARACTER(256)                                     :: single_line_config      , single_line_namelist
    INTEGER                                            :: line_counter_config     , line_counter_namelist
    LOGICAL                                            :: found_match, found_mismatch

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the config and namelist files
    OPEN( UNIT = config_unit, FILE = config_filename, IOSTAT = ios)
    IF (ios /= 0) CALL crash('couldnt open config file "' // TRIM( config_filename) // '"!')

    ! Read one line at a time of the config file, determine the name of the variable
    ! declared in that line, and check if that variable also exists in the namelist file

    found_end_of_file_config = .FALSE.
    line_counter_config      = 0
    found_mismatch           = .FALSE.

    DO WHILE (.NOT. found_end_of_file_config)

      line_counter_config = line_counter_config + 1

      ! Read a single line from the config file
      READ( UNIT = config_unit, FMT = '(A)', IOSTAT = ios) single_line_config

      ! If we've reached the end of the file before finding the terminating forward slash, this config file is not valid.
      IF (ios < 0) CALL crash('config file "' // TRIM( config_filename) // '" is not terminated with a forward slash!')

      ! Remove all leading spaces
      CALL remove_leading_spaces( single_line_config)

      ! The variable name is the part of the string left of the first (, =, or space.
      single_line_config = single_line_config( 1:SCAN( single_line_config, '( =')-1)

      ! Get config variable in all caps for case-insensitive comparison
      CALL capitalise_string( single_line_config)

      ! The forward slash at the end terminates the config file
      IF (single_line_config == '/') THEN
        found_end_of_file_config = .TRUE.
      END IF

      ! Disregard empty lines, commented lines, and the header line
      IF (single_line_config == '' .OR. single_line_config( 1:1) == '&' .OR. single_line_config( 1:1) == '!') THEN
        CYCLE
      END IF

      ! Open the namelist file
      OPEN( UNIT = namelist_unit, FILE = namelist_filename)
      IF (ios /= 0) CALL crash('couldnt open namelist file "' // TRIM( namelist_filename) // '"!')

      ! Read all variables from the namelist file and check if any of them match the current config variable

      found_end_of_file_namelist = .FALSE.
      line_counter_namelist      = 0
      found_match                = .FALSE.

      DO WHILE ((.NOT. found_end_of_file_namelist) .AND. (.NOT. found_match))

        line_counter_namelist = line_counter_namelist + 1

        ! Read a single line from the namelist file
        READ( UNIT = namelist_unit, FMT = '(A)', IOSTAT = ios) single_line_namelist

        ! If we've reached the end of the file before finding the terminating forward slash, this namelist file is not valid.
        IF (ios < 0) CALL crash('namelist file "' // TRIM( namelist_filename) // '" is not terminated with a forward slash!')

        ! Remove all leading spaces
        CALL remove_leading_spaces( single_line_namelist)

        ! The variable name is the part of the string left of the first (, =, or space.
        single_line_namelist = single_line_namelist( 1:SCAN( single_line_namelist, '( =')-1)

        ! Get namelist variable in all caps for case-insensitive comparison
        CALL capitalise_string( single_line_namelist)

        ! The forward slash at the end terminates the config file
        IF (single_line_namelist == '/') THEN
          found_end_of_file_namelist = .TRUE.
        END IF

        ! Disregard empty lines, commented lines, and the header line
        IF (single_line_namelist == '' .OR. single_line_namelist( 1:1) == '&' .OR. single_line_namelist( 1:1) == '!') THEN
          CYCLE
        END IF

        ! Check if this namelist variable matches the config variable
        IF (single_line_namelist == single_line_config) THEN
          found_match = .TRUE.
        END IF

      END DO ! DO WHILE ((.NOT. found_end_of_file_namelist) .AND. (.NOT. found_match))

      ! If no matching variable was found in the namelist file, print an error
      IF (.NOT. found_match) CALL warning('invalid config variable "' // TRIM( single_line_config) // '" in file "' // TRIM( config_filename) // '", line {int_01}')

      ! Close the namelist file
      CLOSE( UNIT = namelist_unit)

    END DO ! DO WHILE (.NOT. found_end_of_file_config)

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! If an invalid config variable was found, crash.
    IF (found_mismatch) CALL crash('found invalid config variables!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_config_file_validity

END MODULE config_file_utilities
