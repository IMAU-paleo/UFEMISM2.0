MODULE main_configuration

  ! The different parameters that control the main model, and paths to config files
  ! controlling all model components for each model region.
  !
  ! Each config variable has two versions: one with the "_config" extension, which is
  ! an actual variable in this module only, and one without the extension, which is
  ! a field in the "C" type. The "_config" variables are used to create a NAMELIST,
  ! which makes reading an external config file really easy - anything in the file that
  ! matches a variable in the namelist overwrites the default value. After that's done,
  ! the fields in the "C" type are replaced with the values of the "_config" variables,
  ! which now have either the default values, or those specified in the external config
  ! file.
  !
  ! While this is certainly very convenient when running the model, it does make adding
  ! new config parameters a bit tedious - you have to add the "_config" variable, add it
  ! as a field in the "C" type, add it to the namelist, and let the "C" type field be
  ! overwritten in the end.
  !
  ! From UFEMISM2.0 on, config files should have ALL config variables. If not, the model
  ! will crash and tell you which ones are missing. This will promote reproducibility of
  ! results, and encourage people to be thorough when setting up their config files.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE config_file_utilities                                  , ONLY: check_config_file_validity
  USE mesh_configuration                                     , ONLY: type_CFG_mesh, initialise_mesh_config_from_file
  USE ice_configuration                                      , ONLY: type_CFG_ice,  initialise_ice_config_from_file

  IMPLICIT NONE

! ===== Configuration variables =====
! ===================================

  ! The "_config" variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!

  ! General model instructions
  ! ==========================

    ! Output directory
    LOGICAL             :: create_procedural_output_dir_config         = .TRUE.                           ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)
    CHARACTER(LEN=256)  :: fixed_output_dir_config                     = 'results_UFEMISM'                ! If not, create a directory with this name instead (stops the program if this directory already exists)

  ! == Time steps and range
  ! =======================

    REAL(dp)            :: start_time_of_run_config                    = 0.0_dp                           ! [yr] Start time of the simulations
    REAL(dp)            :: end_time_of_run_config                      = 0.0_dp                           ! [yr] End   time of the simulations
    REAL(dp)            :: dt_coupling_config                          = 100._dp                          ! [yr] Interval of coupling between the different model regions

  ! == Which model regions to simulate
  ! ==================================

    LOGICAL             :: do_NAM_config                               = .FALSE.
    LOGICAL             :: do_EAS_config                               = .FALSE.
    LOGICAL             :: do_GRL_config                               = .FALSE.
    LOGICAL             :: do_ANT_config                               = .FALSE.

  ! Paths to config files for the model components of each model region
  ! ===================================================================

    ! Mesh generation
    CHARACTER(LEN=256)  :: filename_config_mesh_NAM_config             = ''
    CHARACTER(LEN=256)  :: filename_config_mesh_EAS_config             = ''
    CHARACTER(LEN=256)  :: filename_config_mesh_GRL_config             = ''
    CHARACTER(LEN=256)  :: filename_config_mesh_ANT_config             = ''

    ! The ice-dynamical model
    CHARACTER(LEN=256)  :: filename_config_ice_NAM_config              = ''
    CHARACTER(LEN=256)  :: filename_config_ice_EAS_config              = ''
    CHARACTER(LEN=256)  :: filename_config_ice_GRL_config              = ''
    CHARACTER(LEN=256)  :: filename_config_ice_ANT_config              = ''

! ===== The CFG type =====
! ========================

  ! The "C" type, which contains all the config parameters as fields.
  ! These will all be overwritten with the values of the "_config" variables,
  ! which are either the default values specified above, are the values
  ! specified from the external config file.

  TYPE type_CFG_main
    ! The different parameters that control a UFEMISM simulation

  ! General model instructions
  ! ==========================

    ! Output directory
    LOGICAL             :: create_procedural_output_dir
    CHARACTER(LEN=256)  :: fixed_output_dir
    CHARACTER(LEN=256)  :: output_dir

  ! == Time steps and range
  ! =======================

    REAL(dp)            :: start_time_of_run
    REAL(dp)            :: end_time_of_run
    REAL(dp)            :: dt_coupling

  ! == Which model regions to simulate
  ! ==================================

    LOGICAL             :: do_NAM
    LOGICAL             :: do_EAS
    LOGICAL             :: do_GRL
    LOGICAL             :: do_ANT

  ! Paths to config files for the model components of each model region
  ! ===================================================================

    ! Mesh generation
    CHARACTER(LEN=256)  :: filename_config_mesh_NAM
    CHARACTER(LEN=256)  :: filename_config_mesh_EAS
    CHARACTER(LEN=256)  :: filename_config_mesh_GRL
    CHARACTER(LEN=256)  :: filename_config_mesh_ANT

    ! The ice-dynamical model
    CHARACTER(LEN=256)  :: filename_config_ice_NAM
    CHARACTER(LEN=256)  :: filename_config_ice_EAS
    CHARACTER(LEN=256)  :: filename_config_ice_GRL
    CHARACTER(LEN=256)  :: filename_config_ice_ANT

    ! Config structures for the different model components
    TYPE(type_CFG_mesh)           :: mesh
    TYPE(type_CFG_ice)            :: ice

  END TYPE type_CFG_main

  ! The main config structure
  TYPE(type_CFG_main)   :: C

  ! Copies containing the configuration parameters for each individual model region
  TYPE(type_CFG_main)   :: C_NAM
  TYPE(type_CFG_main)   :: C_EAS
  TYPE(type_CFG_main)   :: C_GRL
  TYPE(type_CFG_main)   :: C_ANT

CONTAINS

! ===== Subroutines ======
! ========================

  SUBROUTINE initialise_model_configuration
    ! Initialise the model configuration

    ! In/output variables:

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_model_configuration'
    CHARACTER(LEN=256)                                 :: config_filename
    CHARACTER(LEN=256)                                 :: output_dir_procedural
    LOGICAL                                            :: ex

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise main config parameters
  ! ====================================

    ! The name of the config file is provided as an input argument when calling the UFEMISM_program
    ! executable. After calling MPI_INIT, only the master process "sees" this argument, so is must be
    ! broadcast to the other processes.

    IF (par%master) THEN
      IF     (iargc() == 1) THEN
        CALL getarg( 1, config_filename)
      ELSE
        CALL crash('run UFEMISM with the path the config file as an argument, e.g. "mpi_exec  -n 2  UFEMISM_program  config-files/config_test"')
      END IF
    END IF ! IF (master) THEN
    CALL MPI_BCAST( config_filename,    256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Running UFEMISM with settings from configuration file: ', colour_string( TRIM( config_filename), 'light blue')

    ! Initialise the main config structure from the config file
    CALL initialise_main_config_from_file( config_filename, C)

  ! == Set up the output directory
  ! ==============================

    ! First get the name of the output directory (either procedural, or provided in the config file)
    C%output_dir = ' '

    IF (C%create_procedural_output_dir) THEN
      ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)

      IF (par%master) THEN
        CALL generate_procedural_output_dir_name( output_dir_procedural)
        C%output_dir( 1:LEN_TRIM( output_dir_procedural)+1) = TRIM( output_dir_procedural) // '/'
      END IF
      CALL MPI_BCAST( C%output_dir, 256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ELSE
      ! Use the provided name (return an error if this directory already exists)

      C%output_dir = TRIM( C%fixed_output_dir) // '/'

    END IF

    ! Create the directory
    IF (par%master) THEN

      ! Safety
      INQUIRE( FILE = TRIM( C%output_dir) // '/.', EXIST = ex)
      IF (ex) THEN
        CALL crash('output directory ' // TRIM( C%output_dir) // ' already exists!')
      END IF

      ! Create output directory
      CALL system('mkdir ' // TRIM( C%output_dir))

      ! Tell the user where it is
      WRITE(0,*) ''
      WRITE(0,*) ' Output directory: ', colour_string( TRIM( C%output_dir), 'light blue')

    END IF ! IF (par%master) THEN
    CALL sync

    ! Copy the config file to the output directory
    IF (par%master) THEN
      CALL system('cp ' // config_filename    // ' ' // TRIM( C%output_dir))
    END IF ! IF (master) THEN
    CALL sync

  ! == Initialise config parameters for model components per model region
  ! =====================================================================

    ! Mesh
    IF (par%master) WRITE(0,*) ''
    IF (C%do_NAM) CALL initialise_mesh_config_from_file( C%filename_config_mesh_NAM, 'NAM', C%output_dir, C_NAM%mesh)
    IF (C%do_EAS) CALL initialise_mesh_config_from_file( C%filename_config_mesh_EAS, 'EAS', C%output_dir, C_EAS%mesh)
    IF (C%do_GRL) CALL initialise_mesh_config_from_file( C%filename_config_mesh_GRL, 'GRL', C%output_dir, C_GRL%mesh)
    IF (C%do_ANT) CALL initialise_mesh_config_from_file( C%filename_config_mesh_ANT, 'ANT', C%output_dir, C_ANT%mesh)

    ! Ice
    IF (par%master) WRITE(0,*) ''
    IF (C%do_NAM) CALL initialise_ice_config_from_file(  C%filename_config_ice_NAM , 'NAM', C%output_dir, C_NAM%ice )
    IF (C%do_EAS) CALL initialise_ice_config_from_file(  C%filename_config_ice_EAS , 'EAS', C%output_dir, C_EAS%ice )
    IF (C%do_GRL) CALL initialise_ice_config_from_file(  C%filename_config_ice_GRL , 'GRL', C%output_dir, C_GRL%ice )
    IF (C%do_ANT) CALL initialise_ice_config_from_file(  C%filename_config_ice_ANT , 'ANT', C%output_dir, C_ANT%ice )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_model_configuration

  SUBROUTINE initialise_main_config_from_file( config_filename, CFG)
    ! Initialise a config structure from an external config text file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename
    TYPE(type_CFG_main),             INTENT(OUT)       :: CFG

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_main_config_from_file'
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Let each of the processors read the config file in turns so there's no access conflicts
    DO i = 0, par%n-1

      IF (i == par%i) THEN

        ! Read the external file, use a Fortran NAMELIST to overwrite the default
        ! values of the XXX_config variables
        CALL read_main_config_file( config_filename)

        ! Copy values from the XXX_config variables to the config structure
        CALL copy_main_config_variables_to_struct( CFG)

      END IF ! IF (i == par%i) THEN

      ! Make sure only one process at a time reads from / writes to disk
      CALL sync

    END DO ! DO i = 0, par%n-1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_main_config_from_file

  SUBROUTINE read_main_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.
    !
    ! Note: make sure that only one process at a time calls this subroutine!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_main_config_file'
    CHARACTER(LEN=256), PARAMETER                      :: namelist_filename = 'config_namelist_temp.txt'
    INTEGER, PARAMETER                                 :: config_unit    = 1337
    INTEGER, PARAMETER                                 :: namelist_unit  = 1338
    INTEGER                                            :: ios

    ! The NAMELIST that's used to read the external config file.
    NAMELIST /CONFIG_MAIN/&
      create_procedural_output_dir_config             , &
      fixed_output_dir_config                         , &
      start_time_of_run_config                        , &
      end_time_of_run_config                          , &
      dt_coupling_config                              , &
      do_NAM_config                                   , &
      do_EAS_config                                   , &
      do_GRL_config                                   , &
      do_ANT_config                                   , &
      filename_config_mesh_NAM_config                 , &
      filename_config_mesh_EAS_config                 , &
      filename_config_mesh_GRL_config                 , &
      filename_config_mesh_ANT_config                 , &
      filename_config_ice_NAM_config                  , &
      filename_config_ice_EAS_config                  , &
      filename_config_ice_GRL_config                  , &
      filename_config_ice_ANT_config

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write the CONFIG namelist to a temporary file
    OPEN(  UNIT = namelist_unit, FILE = TRIM( namelist_filename))
    WRITE( UNIT = namelist_unit, NML  = CONFIG_MAIN)
    CLOSE( UNIT = namelist_unit)

    ! Check the config file for validity
    CALL check_config_file_validity( config_filename, namelist_filename)

    ! Delete the temporary namelist file
    CALL system('rm -f ' // TRIM( namelist_filename))

    ! Open the config file
    OPEN(  UNIT = config_unit, FILE = TRIM( config_filename), STATUS = 'OLD', ACTION = 'READ', IOSTAT = ios)
    IF (ios /= 0) CALL crash('couldnt open config file "' // TRIM( config_filename) // '"!')

    ! Read the config file using the CONFIG namelist
    READ(  UNIT = config_unit, NML = CONFIG_MAIN, IOSTAT = ios)
    IF (ios /= 0) CALL crash('error while reading config file "' // TRIM( config_filename) // '"!')

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_main_config_file

  SUBROUTINE copy_main_config_variables_to_struct( CFG)
    ! Overwrite the values in the fields of the config structure with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values that were read from the config file.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_CFG_main),             INTENT(OUT)       :: CFG

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'copy_main_config_variables_to_struct'

    ! Add routine to path
    CALL init_routine( routine_name)

  ! General model instructions
  ! ==========================

    ! Output directory
    CFG%create_procedural_output_dir             = create_procedural_output_dir_config
    CFG%fixed_output_dir                         = fixed_output_dir_config

  ! == Time steps and range
  ! =======================

    CFG%start_time_of_run                        = start_time_of_run_config
    CFG%end_time_of_run                          = end_time_of_run_config
    CFG%dt_coupling                              = dt_coupling_config

  ! == Which model regions to simulate
  ! ==================================

    C%do_NAM                                     = do_NAM_config
    C%do_EAS                                     = do_EAS_config
    C%do_GRL                                     = do_GRL_config
    C%do_ANT                                     = do_ANT_config

  ! Paths to config files for the model components of each model region
  ! ===================================================================

    ! Mesh generation
    CFG%filename_config_mesh_NAM                 = filename_config_mesh_NAM_config
    CFG%filename_config_mesh_EAS                 = filename_config_mesh_EAS_config
    CFG%filename_config_mesh_GRL                 = filename_config_mesh_GRL_config
    CFG%filename_config_mesh_ANT                 = filename_config_mesh_ANT_config

    ! The ice-dynamical model
    CFG%filename_config_ice_NAM                  = filename_config_ice_NAM_config
    CFG%filename_config_ice_EAS                  = filename_config_ice_EAS_config
    CFG%filename_config_ice_GRL                  = filename_config_ice_GRL_config
    CFG%filename_config_ice_ANT                  = filename_config_ice_ANT_config

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE copy_main_config_variables_to_struct

  SUBROUTINE generate_procedural_output_dir_name( output_dir)
    ! Generate a procedural output directory for the current date (e.g. results_20210721_001)
    ! Keep increasing the counter at the end until a directory is available.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(256),                      INTENT(OUT)   :: output_dir

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'generate_procedural_output_dir_name'
    INTEGER,  DIMENSION(8)                             :: values
    LOGICAL                                            :: ex

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    output_dir = ' '

    ! Get current date and time
    CALL date_and_time( VALUES = values)

    ! Get proper year (assume we're still in the 21st century...)
    output_dir(1:10) = 'results_20'
    SELECT CASE( FLOOR(REAL(values(1))/10._dp)-200)
    CASE(0)
      output_dir(11:11) = '0'
    CASE(1)
      output_dir(11:11) = '1'
    CASE(2)
      output_dir(11:11) = '2'
    CASE(3)
      output_dir(11:11) = '3'
    CASE(4)
      output_dir(11:11) = '4'
    CASE(5)
      output_dir(11:11) = '5'
    CASE(6)
      output_dir(11:11) = '6'
    CASE(7)
      output_dir(11:11) = '7'
    CASE(8)
      output_dir(11:11) = '8'
    CASE(9)
      output_dir(11:11) = '9'
    CASE DEFAULT
      CALL crash('error retrieving date and time!')
    END SELECT

    SELECT CASE( MOD(values(1),10))
    CASE(0)
      output_dir(12:12) = '0'
    CASE(1)
      output_dir(12:12) = '1'
    CASE(2)
      output_dir(12:12) = '2'
    CASE(3)
      output_dir(12:12) = '3'
    CASE(4)
      output_dir(12:12) = '4'
    CASE(5)
      output_dir(12:12) = '5'
    CASE(6)
      output_dir(12:12) = '6'
    CASE(7)
      output_dir(12:12) = '7'
    CASE(8)
      output_dir(12:12) = '8'
    CASE(9)
      output_dir(12:12) = '9'
    CASE DEFAULT
      CALL crash('error retrieving date and time!')
    END SELECT

    SELECT CASE( values(2))
    CASE(1)
      output_dir(13:14) = '01'
    CASE(2)
      output_dir(13:14) = '02'
    CASE(3)
      output_dir(13:14) = '03'
    CASE(4)
      output_dir(13:14) = '04'
    CASE(5)
      output_dir(13:14) = '05'
    CASE(6)
      output_dir(13:14) = '06'
    CASE(7)
      output_dir(13:14) = '07'
    CASE(8)
      output_dir(13:14) = '08'
    CASE(9)
      output_dir(13:14) = '09'
    CASE(10)
      output_dir(13:14) = '10'
    CASE(11)
      output_dir(13:14) = '11'
    CASE(12)
      output_dir(13:14) = '12'
    CASE DEFAULT
      CALL crash('error retrieving date and time!')
    END SELECT

    SELECT CASE( FLOOR(REAL(values(3))/10._dp))
    CASE(0)
      output_dir(15:15) = '0'
    CASE(1)
      output_dir(15:15) = '1'
    CASE(2)
      output_dir(15:15) = '2'
    CASE(3)
      output_dir(15:15) = '3'
    CASE DEFAULT
      CALL crash('error retrieving date and time!')
    END SELECT

    SELECT CASE( MOD(values(3),10))
    CASE(0)
      output_dir(16:16) = '0'
    CASE(1)
      output_dir(16:16) = '1'
    CASE(2)
      output_dir(16:16) = '2'
    CASE(3)
      output_dir(16:16) = '3'
    CASE(4)
      output_dir(16:16) = '4'
    CASE(5)
      output_dir(16:16) = '5'
    CASE(6)
      output_dir(16:16) = '6'
    CASE(7)
      output_dir(16:16) = '7'
    CASE(8)
      output_dir(16:16) = '8'
    CASE(9)
      output_dir(16:16) = '9'
    CASE DEFAULT
      CALL crash('error retrieving date and time!')
    END SELECT

    output_dir(17:20) = '_001'

    INQUIRE( FILE = TRIM( output_dir) // '/.', EXIST = ex)

    DO WHILE (ex)

     IF      (output_dir(20:20) == '0') THEN
       output_dir(20:20) = '1'
     ELSE IF (output_dir(20:20) == '1') THEN
       output_dir(20:20) = '2'
     ELSE IF (output_dir(20:20) == '2') THEN
       output_dir(20:20) = '3'
     ELSE IF (output_dir(20:20) == '3') THEN
       output_dir(20:20) = '4'
     ELSE IF (output_dir(20:20) == '4') THEN
       output_dir(20:20) = '5'
     ELSE IF (output_dir(20:20) == '5') THEN
       output_dir(20:20) = '6'
     ELSE IF (output_dir(20:20) == '6') THEN
       output_dir(20:20) = '7'
     ELSE IF (output_dir(20:20) == '7') THEN
       output_dir(20:20) = '8'
     ELSE IF (output_dir(20:20) == '8') THEN
       output_dir(20:20) = '9'
     ELSE IF (output_dir(20:20) == '9') THEN
       output_dir(20:20) = '0'

       IF      (output_dir(19:19) == '0') THEN
         output_dir(19:19) = '1'
       ELSE IF (output_dir(19:19) == '1') THEN
         output_dir(19:19) = '2'
       ELSE IF (output_dir(19:19) == '2') THEN
         output_dir(19:19) = '3'
       ELSE IF (output_dir(19:19) == '3') THEN
         output_dir(19:19) = '4'
       ELSE IF (output_dir(19:19) == '4') THEN
         output_dir(19:19) = '5'
       ELSE IF (output_dir(19:19) == '5') THEN
         output_dir(19:19) = '6'
       ELSE IF (output_dir(19:19) == '6') THEN
         output_dir(19:19) = '7'
       ELSE IF (output_dir(19:19) == '7') THEN
         output_dir(19:19) = '8'
       ELSE IF (output_dir(19:19) == '8') THEN
         output_dir(19:19) = '9'
       ELSE IF (output_dir(19:19) == '9') THEN
         output_dir(19:19) = '0'

         IF      (output_dir(18:18) == '0') THEN
           output_dir(18:18) = '1'
         ELSE IF (output_dir(18:18) == '1') THEN
           output_dir(18:18) = '2'
         ELSE IF (output_dir(18:18) == '2') THEN
           output_dir(18:18) = '3'
         ELSE IF (output_dir(18:18) == '3') THEN
           output_dir(18:18) = '4'
         ELSE IF (output_dir(18:18) == '4') THEN
           output_dir(18:18) = '5'
         ELSE IF (output_dir(18:18) == '5') THEN
           output_dir(18:18) = '6'
         ELSE IF (output_dir(18:18) == '6') THEN
           output_dir(18:18) = '7'
         ELSE IF (output_dir(18:18) == '7') THEN
           output_dir(18:18) = '8'
         ELSE IF (output_dir(18:18) == '8') THEN
           output_dir(18:18) = '9'
         ELSE IF (output_dir(18:18) == '9') THEN
           output_dir(18:18) = '0'
         END IF

       END IF

     END IF

     INQUIRE( FILE = TRIM( output_dir) // '/.', EXIST = ex)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE generate_procedural_output_dir_name

END MODULE main_configuration
