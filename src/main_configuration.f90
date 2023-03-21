MODULE mesh_configuration

  ! The main model configuration parameters

  ! The way it's done right now:
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
  ! Some day I'll figure out a more elegant solution for this...

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Configuration variables =====
! ===================================

  ! The "_config  variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!

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

  END TYPE type_CFG_mesh

CONTAINS

! ===== Subroutines ======
! ========================

END MODULE mesh_configuration
