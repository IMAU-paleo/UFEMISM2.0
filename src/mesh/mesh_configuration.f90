MODULE mesh_configuration

  ! The different parameters that control mesh generation
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

  IMPLICIT NONE

! ===== Configuration variables =====
! ===================================

  ! The "_config  variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!
  ! The values listed here are the default values that are used when the variable
  ! is not listed in the config file.

    ! Resolutions for different regions
    REAL(dp)            :: maximum_resolution_uniform_config           = 500._dp                          ! [m]          Maximum resolution for the entire domain
    REAL(dp)            :: maximum_resolution_grounded_ice_config      = 200._dp                          ! [m]          Maximum resolution for grounded ice
    REAL(dp)            :: maximum_resolution_floating_ice_config      = 100._dp                          ! [m]          Maximum resolution for floating ice
    REAL(dp)            :: maximum_resolution_grounding_line_config    = 30._dp                           ! [m]          Maximum resolution for the grounding line
    REAL(dp)            :: grounding_line_width_config                 = 50._dp                           ! [m]          Width of the band around the grounding line that should get this resolution
    REAL(dp)            :: maximum_resolution_calving_front_config     = 50._dp                           ! [m]          Maximum resolution for the calving front
    REAL(dp)            :: calving_front_width_config                  = 100._dp                          ! [m]          Width of the band around the calving front that should get this resolution

    ! Advanced geometry parameters
    REAL(dp)            :: alpha_min_config                            = 0.4363_dp                        ! [radians]    Smallest allowed internal triangle angle
    INTEGER             :: nit_Lloyds_algorithm_config                 = 3                                ! [-]          Number of iterations of Lloyds algorithm to be applied after refinement

    ! Memory
    INTEGER             :: nC_mem_config                               = 32                               ! [-]          How many columns of memory should be allocated for connectivity lists

    ! The scaled vertical coordinate zeta
    CHARACTER(LEN=256)  :: choice_zeta_grid_config                     = 'regular'                        ! The type of vertical grid to use; can be "regular", "irregular_log", "old_15_layer_zeta"
    INTEGER             :: nz_config                                   = 12                               ! The number of vertical layers to use
    REAL(dp)            :: zeta_irregular_log_R_config                 = 10._dp                           ! Ratio between surface and base layer spacings

! ===== The CFG type =====
! ========================

  ! The "C" type, which contains all the config parameters as fields.
  ! These will all be overwritten with the values of the "_config" variables,
  ! which are either the default values specified above, are the values
  ! specified from the external config file.

  TYPE type_CFG_mesh
    ! The different parameters that control mesh generation

    ! Resolutions for different regions
    REAL(dp)            :: maximum_resolution_uniform
    REAL(dp)            :: maximum_resolution_grounded_ice
    REAL(dp)            :: maximum_resolution_floating_ice
    REAL(dp)            :: maximum_resolution_grounding_line
    REAL(dp)            :: grounding_line_width
    REAL(dp)            :: maximum_resolution_calving_front
    REAL(dp)            :: calving_front_width

    ! Advanced geometry parameters
    REAL(dp)            :: alpha_min
    INTEGER             :: nit_Lloyds_algorithm

    ! Memory
    INTEGER             :: nC_mem

    ! The scaled vertical coordinate zeta
    CHARACTER(LEN=256)  :: choice_zeta_grid
    INTEGER             :: nz
    REAL(dp)            :: zeta_irregular_log_R

  END TYPE type_CFG_mesh

CONTAINS

! ===== Subroutines ======
! ========================

  SUBROUTINE initialise_mesh_config_from_file( config_filename, region_name, output_dir, CFG)
    ! Initialise a config structure from an external config text file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename
    CHARACTER(LEN=3),                INTENT(IN)        :: region_name
    CHARACTER(LEN=*),                INTENT(IN)        :: output_dir
    TYPE(type_CFG_mesh),             INTENT(OUT)       :: CFG

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_mesh_config_from_file'
    CHARACTER(LEN=256)                                 :: config_filename_reg
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ' Initialising mesh settings for region ', colour_string( region_name, 'light blue'), &
      ' from configuration file: ', colour_string( TRIM( config_filename), 'light blue')

    ! Let each of the processors read the config file in turns so there's no access conflicts
    DO i = 0, par%n-1

      IF (i == par%i) THEN

        ! Read the external file, use a Fortran NAMELIST to overwrite the default
        ! values of the XXX_config variables
        CALL read_mesh_config_file( config_filename)

        ! Copy values from the XXX_config variables to the config structure
        CALL copy_mesh_config_variables_to_struct( CFG)

      END IF ! IF (i == par%i) THEN

      ! Make sure only one process at a time reads from / writes to disk
      CALL sync

    END DO ! DO i = 0, par%n-1

    ! Copy the config file to the output directory
    IF (par%master) THEN
      CALL system('cp ' // TRIM( config_filename) // ' ' // TRIM( output_dir) // '/config_mesh_' // region_name // '.cfg')
    END IF ! IF (master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mesh_config_from_file

  SUBROUTINE read_mesh_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.
    !
    ! Note: make sure that only one process at a time calls this subroutine!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_mesh_config_file'
    CHARACTER(LEN=256), PARAMETER                      :: namelist_filename = 'config_namelist_temp.txt'
    INTEGER, PARAMETER                                 :: config_unit    = 1337
    INTEGER, PARAMETER                                 :: namelist_unit  = 1338
    INTEGER                                            :: ios

    ! The NAMELIST that's used to read the external config file.
    NAMELIST /CONFIG_MESH/&
      maximum_resolution_uniform_config               , &
      maximum_resolution_grounded_ice_config          , &
      maximum_resolution_floating_ice_config          , &
      maximum_resolution_grounding_line_config        , &
      grounding_line_width_config                     , &
      maximum_resolution_calving_front_config         , &
      calving_front_width_config                      , &
      alpha_min_config                                , &
      nit_Lloyds_algorithm_config                     , &
      nC_mem_config                                   , &
      choice_zeta_grid_config                         , &
      nz_config                                       , &
      zeta_irregular_log_R_config

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write the CONFIG namelist to a temporary file
    OPEN(  UNIT = namelist_unit, FILE = TRIM( namelist_filename))
    WRITE( UNIT = namelist_unit, NML  = CONFIG_MESH)
    CLOSE( UNIT = namelist_unit)

    ! Check the config file for validity
    CALL check_config_file_validity( config_filename, namelist_filename)

    ! Delete the temporary namelist file
    CALL system('rm -f ' // TRIM( namelist_filename))

    ! Open the config file
    OPEN(  UNIT = config_unit, FILE = TRIM( config_filename), STATUS = 'OLD', ACTION = 'READ', IOSTAT = ios)
    IF (ios /= 0) CALL crash('couldnt open config file "' // TRIM( config_filename) // '"!')

    ! Read the config file using the CONFIG namelist
    READ(  UNIT = config_unit, NML = CONFIG_MESH, IOSTAT = ios)
    IF (ios /= 0) CALL crash('error while reading config file "' // TRIM( config_filename) // '"!')

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_mesh_config_file

  SUBROUTINE copy_mesh_config_variables_to_struct( CFG)
    ! Overwrite the values in the fields of the config structure with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values that were read from the config file.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_CFG_mesh),             INTENT(OUT)       :: CFG

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'copy_mesh_config_variables_to_struct'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Resolutions for different regions
    CFG%maximum_resolution_uniform               = maximum_resolution_uniform_config
    CFG%maximum_resolution_grounded_ice          = maximum_resolution_grounded_ice_config
    CFG%maximum_resolution_floating_ice          = maximum_resolution_floating_ice_config
    CFG%maximum_resolution_grounding_line        = maximum_resolution_grounding_line_config
    CFG%grounding_line_width                     = grounding_line_width_config
    CFG%maximum_resolution_calving_front         = maximum_resolution_calving_front_config
    CFG%calving_front_width                      = calving_front_width_config

    ! Advanced geometry parameters
    CFG%alpha_min                                = alpha_min_config
    CFG%nit_Lloyds_algorithm                     = nit_Lloyds_algorithm_config

    ! Memory
    CFG%nC_mem                                   = nC_mem_config

    ! The scaled vertical coordinate zeta
    CFG%choice_zeta_grid                         = choice_zeta_grid_config
    CFG%nz                                       = nz_config
    CFG%zeta_irregular_log_R                     = zeta_irregular_log_R_config

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE copy_mesh_config_variables_to_struct

END MODULE mesh_configuration
