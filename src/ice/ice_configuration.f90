MODULE ice_configuration

  ! The different parameters that control the ice-dynamical model
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

  ! == Ice dynamics - velocity
  ! ==========================

    CHARACTER(LEN=256)  :: choice_stress_balance_approximation_config  = 'SIA'                            ! Choice of stress balance approximation: "none" (= no flow, though geometry can still change due to mass balance), "SIA", "SSA", "SIA/SSA", "DIVA", "BPA"
    REAL(dp)            :: n_flow_config                               = 3.0_dp                           ! Exponent in Glen's flow law
    REAL(dp)            :: m_enh_sheet_config                          = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
    REAL(dp)            :: m_enh_shelf_config                          = 1.0_dp                           ! Ice flow enhancement factor for floating ice
    CHARACTER(LEN=256)  :: choice_hybrid_SIASSA_scheme_config          = 'add'                            ! Choice of scheme for combining SIA and SSA velocities in the hybrid approach
    LOGICAL             :: do_GL_subgrid_friction_config               = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
    REAL(dp)            :: subgrid_friction_exponent_config            = 2._dp                            ! Exponent to which f_grnd should be raised before being used to scale beta

    ! Some parameters for numerically solving the stress balance
    REAL(dp)            :: SIA_maximum_diffusivity_config              = 1E5_dp                           ! Limit the diffusivity in the SIA to this value
    REAL(dp)            :: visc_it_norm_dUV_tol_config                 = 1E-2_dp                          ! Stop criterion for the viscosity iteration: the L2-norm of successive velocity solutions should be smaller than this number
    INTEGER             :: visc_it_nit_config                          = 50                               ! Maximum number of effective viscosity iterations
    REAL(dp)            :: visc_it_relax_config                        = 0.4_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    REAL(dp)            :: epsilon_sq_0_config                         = 1E-15_dp                         ! Normalisation term so that zero velocity gives non-zero viscosity
    REAL(dp)            :: visc_eff_min_config                         = 1E3_dp                           ! Minimum value for effective viscosity
    REAL(dp)            :: beta_max_config                             = 1E20_dp                          ! Maximum value for basal friction coefficient
    REAL(dp)            :: vel_max_config                              = 5000._dp                         ! Velocities are limited to this value
    REAL(dp)            :: stressbalance_PETSc_rtol_config             = 1E-2_dp                          ! PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    REAL(dp)            :: stressbalance_PETSc_abstol_config           = 1.0_dp                           ! PETSc solver - stop criterion, absolute difference

  ! == Ice dynamics - time integration
  ! ==================================

    CHARACTER(LEN=256)  :: choice_timestepping_config                  = 'pc'                             ! Choice of timestepping method: "direct", "pc" (NOTE: 'direct' does not work with DIVA ice dynamcis!)
    CHARACTER(LEN=256)  :: choice_ice_integration_method_config        = 'explicit'                       ! Choice of ice thickness integration scheme: "none" (i.e. unchanging geometry), "explicit", "semi-implicit"
    REAL(dp)            :: dHi_PETSc_rtol_config                       = 0.001_dp                         ! dHi PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    REAL(dp)            :: dHi_PETSc_abstol_config                     = 0.001_dp                         ! dHi PETSc solver - stop criterion, absolute difference

    ! Predictor-corrector ice-thickness update
    REAL(dp)            :: pc_epsilon_config                           = 3._dp                            ! Target truncation error in dHi_dt [m/yr] (epsilon in Robinson et al., 2020, Eq. 33)
    REAL(dp)            :: pc_k_I_config                               = 0.2_dp                           ! Exponent k_I in  Robinson et al., 2020, Eq. 33
    REAL(dp)            :: pc_k_p_config                               = 0.2_dp                           ! Exponent k_p in  Robinson et al., 2020, Eq. 33
    REAL(dp)            :: pc_eta_min_config                           = 1E-8_dp                          ! Normalisation term in estimation of the truncation error (Robinson et al., Eq. 32)

  ! == Ice dynamics - boundary conditions
  ! ====================================

    CHARACTER(LEN=256)  :: BC_u_west_config                            = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain border
    CHARACTER(LEN=256)  :: BC_u_east_config                            = 'infinite'                       ! Allowed choices: "infinite", "zero", "periodic_ISMIP-HOM"
    CHARACTER(LEN=256)  :: BC_u_south_config                           = 'infinite'
    CHARACTER(LEN=256)  :: BC_u_north_config                           = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_west_config                            = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_east_config                            = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_south_config                           = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_north_config                           = 'infinite'
    CHARACTER(LEN=256)  :: BC_H_west_config                            = 'zero'                           ! Boundary conditions for ice thickness at the domain boundary
    CHARACTER(LEN=256)  :: BC_H_east_config                            = 'zero'                           ! Allowed choices:  "infinite", "zero", "ISMIP_HOM_F"
    CHARACTER(LEN=256)  :: BC_H_south_config                           = 'zero'
    CHARACTER(LEN=256)  :: BC_H_north_config                           = 'zero'

  ! == Ice dynamics - stabilisation
  ! ===============================

    CHARACTER(LEN=256)  :: choice_mask_noice_config                    = 'none'                           ! Choice of mask_noice configuration

    ! Partially fixed geometry, useful for initialisation and inversion runs
    LOGICAL             :: fixed_shelf_geometry_config                 = .FALSE.                          ! Keep geometry of floating ice fixed
    LOGICAL             :: fixed_sheet_geometry_config                 = .FALSE.                          ! Keep geometry of grounded ice fixed
    LOGICAL             :: fixed_grounding_line_config                 = .FALSE.                          ! Keep ice thickness at the grounding line fixed

! ===== The CFG type =====
! ========================

  ! The "C" type, which contains all the config parameters as fields.
  ! These will all be overwritten with the values of the "_config" variables,
  ! which are either the default values specified above, are the values
  ! specified from the external config file.

  TYPE type_CFG_ice
    ! The different parameters that control the ice-dynamical model

  ! == Ice dynamics - velocity
  ! ==========================

    CHARACTER(LEN=256)  :: choice_stress_balance_approximation
    REAL(dp)            :: n_flow
    REAL(dp)            :: m_enh_sheet
    REAL(dp)            :: m_enh_shelf
    CHARACTER(LEN=256)  :: choice_hybrid_SIASSA_scheme
    LOGICAL             :: do_GL_subgrid_friction
    REAL(dp)            :: subgrid_friction_exponent

    ! Some parameters for numerically solving the stress balance
    REAL(dp)            :: SIA_maximum_diffusivity
    REAL(dp)            :: visc_it_norm_dUV_tol
    INTEGER             :: visc_it_nit
    REAL(dp)            :: visc_it_relax
    REAL(dp)            :: epsilon_sq_0
    REAL(dp)            :: visc_eff_min
    REAL(dp)            :: beta_max
    REAL(dp)            :: vel_max
    REAL(dp)            :: stressbalance_PETSc_rtol
    REAL(dp)            :: stressbalance_PETSc_abstol

  ! == Ice dynamics - time integration
  ! ==================================

    CHARACTER(LEN=256)  :: choice_timestepping
    CHARACTER(LEN=256)  :: choice_ice_integration_method
    REAL(dp)            :: dHi_PETSc_rtol
    REAL(dp)            :: dHi_PETSc_abstol

    ! Predictor-corrector ice-thickness update
    REAL(dp)            :: pc_epsilon
    REAL(dp)            :: pc_k_I
    REAL(dp)            :: pc_k_p
    REAL(dp)            :: pc_eta_min

  ! == Ice dynamics - boundary conditions
  ! ====================================

    CHARACTER(LEN=256)  :: BC_u_west
    CHARACTER(LEN=256)  :: BC_u_east
    CHARACTER(LEN=256)  :: BC_u_south
    CHARACTER(LEN=256)  :: BC_u_north
    CHARACTER(LEN=256)  :: BC_v_west
    CHARACTER(LEN=256)  :: BC_v_east
    CHARACTER(LEN=256)  :: BC_v_south
    CHARACTER(LEN=256)  :: BC_v_north
    CHARACTER(LEN=256)  :: BC_H_west
    CHARACTER(LEN=256)  :: BC_H_east
    CHARACTER(LEN=256)  :: BC_H_south
    CHARACTER(LEN=256)  :: BC_H_north

  ! == Ice dynamics - stabilisation
  ! ===============================

    CHARACTER(LEN=256)  :: choice_mask_noice

    ! Partially fixed geometry, useful for initialisation and inversion runs
    LOGICAL             :: fixed_shelf_geometry
    LOGICAL             :: fixed_sheet_geometry
    LOGICAL             :: fixed_grounding_line

  END TYPE type_CFG_ice

CONTAINS

! ===== Subroutines ======
! ========================

  SUBROUTINE initialise_ice_config_from_file( config_filename, region_name, output_dir, CFG)
    ! Initialise a config structure from an external config text file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename
    CHARACTER(LEN=3),                INTENT(IN)        :: region_name
    CHARACTER(LEN=*),                INTENT(IN)        :: output_dir
    TYPE(type_CFG_ice),              INTENT(OUT)       :: CFG

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_config_from_file'
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ' Initialising ice settings for region ', colour_string( region_name, 'light blue'), &
      ' from configuration file: ', colour_string( TRIM( config_filename), 'light blue')

    ! Let each of the processors read the config file in turns so there's no access conflicts
    DO i = 0, par%n-1

      IF (i == par%i) THEN

        ! Read the external file, use a Fortran NAMELIST to overwrite the default
        ! values of the XXX_config variables
        CALL read_ice_config_file( config_filename)

        ! Copy values from the XXX_config variables to the config structure
        CALL copy_ice_config_variables_to_struct( CFG)

      END IF ! IF (i == par%i) THEN

      ! Make sure only one process at a time reads from / writes to disk
      CALL sync

    END DO ! DO i = 0, par%n-1

    ! Copy the config file to the output directory
    IF (par%master) THEN
      CALL system('cp ' // TRIM( config_filename) // ' ' // TRIM( output_dir) // '/config_ice_' // region_name // '.cfg')
    END IF ! IF (master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_config_from_file

  SUBROUTINE read_ice_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.
    !
    ! Note: make sure that only one process at a time calls this subroutine!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_ice_config_file'
    CHARACTER(LEN=256), PARAMETER                      :: namelist_filename = 'config_namelist_temp.txt'
    INTEGER, PARAMETER                                 :: config_unit    = 1337
    INTEGER, PARAMETER                                 :: namelist_unit  = 1338
    INTEGER                                            :: ios

    ! The NAMELIST that's used to read the external config file.
    NAMELIST /CONFIG_ICE/&
      choice_stress_balance_approximation_config      , &
      n_flow_config                                   , &
      m_enh_sheet_config                              , &
      m_enh_shelf_config                              , &
      choice_hybrid_SIASSA_scheme_config              , &
      do_GL_subgrid_friction_config                   , &
      subgrid_friction_exponent_config                , &
      SIA_maximum_diffusivity_config                  , &
      visc_it_norm_dUV_tol_config                     , &
      visc_it_nit_config                              , &
      visc_it_relax_config                            , &
      epsilon_sq_0_config                             , &
      visc_eff_min_config                             , &
      beta_max_config                                 , &
      vel_max_config                                  , &
      stressbalance_PETSc_rtol_config                 , &
      stressbalance_PETSc_abstol_config               , &
      choice_timestepping_config                      , &
      choice_ice_integration_method_config            , &
      dHi_PETSc_rtol_config                           , &
      dHi_PETSc_abstol_config                         , &
      pc_epsilon_config                               , &
      pc_k_I_config                                   , &
      pc_k_p_config                                   , &
      pc_eta_min_config                               , &
      BC_u_west_config                                , &
      BC_u_east_config                                , &
      BC_u_south_config                               , &
      BC_u_north_config                               , &
      BC_v_west_config                                , &
      BC_v_east_config                                , &
      BC_v_south_config                               , &
      BC_v_north_config                               , &
      BC_H_west_config                                , &
      BC_H_east_config                                , &
      BC_H_south_config                               , &
      BC_H_north_config                               , &
      choice_mask_noice_config                        , &
      fixed_shelf_geometry_config                     , &
      fixed_sheet_geometry_config                     , &
      fixed_grounding_line_config

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write the CONFIG namelist to a temporary file
    OPEN(  UNIT = namelist_unit, FILE = TRIM( namelist_filename))
    WRITE( UNIT = namelist_unit, NML  = CONFIG_ICE)
    CLOSE( UNIT = namelist_unit)

    ! Check the config file for validity
    CALL check_config_file_validity( config_filename, namelist_filename)

    ! Delete the temporary namelist file
    CALL system('rm -f ' // TRIM( namelist_filename))

    ! Open the config file
    OPEN(  UNIT = config_unit, FILE = TRIM( config_filename), STATUS = 'OLD', ACTION = 'READ', IOSTAT = ios)
    IF (ios /= 0) CALL crash('couldnt open config file "' // TRIM( config_filename) // '"!')

    ! Read the config file using the CONFIG namelist
    READ(  UNIT = config_unit, NML = CONFIG_ICE, IOSTAT = ios)
    IF (ios /= 0) CALL crash('error while reading config file "' // TRIM( config_filename) // '"!')

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ice_config_file

  SUBROUTINE copy_ice_config_variables_to_struct( CFG)
    ! Overwrite the values in the fields of the config structure with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values that were read from the config file.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_CFG_ice),              INTENT(OUT)       :: CFG

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'copy_ice_config_variables_to_struct'

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Ice dynamics - velocity
  ! ==========================

    CFG%choice_stress_balance_approximation      = choice_stress_balance_approximation_config
    CFG%n_flow                                   = n_flow_config
    CFG%m_enh_sheet                              = m_enh_sheet_config
    CFG%m_enh_shelf                              = m_enh_shelf_config
    CFG%choice_hybrid_SIASSA_scheme              = choice_hybrid_SIASSA_scheme_config
    CFG%do_GL_subgrid_friction                   = do_GL_subgrid_friction_config
    CFG%subgrid_friction_exponent                = subgrid_friction_exponent_config

    ! Some parameters for numerically solving the stress balance
    CFG%SIA_maximum_diffusivity                  = SIA_maximum_diffusivity_config
    CFG%visc_it_norm_dUV_tol                     = visc_it_norm_dUV_tol_config
    CFG%visc_it_nit                              = visc_it_nit_config
    CFG%visc_it_relax                            = visc_it_relax_config
    CFG%epsilon_sq_0                             = epsilon_sq_0_config
    CFG%visc_eff_min                             = visc_eff_min_config
    CFG%beta_max                                 = beta_max_config
    CFG%vel_max                                  = vel_max_config
    CFG%stressbalance_PETSc_rtol                 = stressbalance_PETSc_rtol_config
    CFG%stressbalance_PETSc_abstol               = stressbalance_PETSc_abstol_config

  ! == Ice dynamics - time integration
  ! ==================================

    CFG%choice_timestepping                      = choice_timestepping_config
    CFG%choice_ice_integration_method            = choice_ice_integration_method_config
    CFG%dHi_PETSc_rtol                           = dHi_PETSc_rtol_config
    CFG%dHi_PETSc_abstol                         = dHi_PETSc_abstol_config

    ! Predictor-corrector ice-thickness update
    CFG%pc_epsilon                               = pc_epsilon_config
    CFG%pc_k_I                                   = pc_k_I_config
    CFG%pc_k_p                                   = pc_k_p_config
    CFG%pc_eta_min                               = pc_eta_min_config

  ! == Ice dynamics - boundary conditions
  ! ====================================

    CFG%BC_u_west                                = BC_u_west_config
    CFG%BC_u_east                                = BC_u_east_config
    CFG%BC_u_south                               = BC_u_south_config
    CFG%BC_u_north                               = BC_u_north_config
    CFG%BC_v_west                                = BC_v_west_config
    CFG%BC_v_east                                = BC_v_east_config
    CFG%BC_v_south                               = BC_v_south_config
    CFG%BC_v_north                               = BC_v_north_config
    CFG%BC_H_west                                = BC_H_west_config
    CFG%BC_H_east                                = BC_H_east_config
    CFG%BC_H_south                               = BC_H_south_config
    CFG%BC_H_north                               = BC_H_north_config

  ! == Ice dynamics - stabilisation
  ! ===============================

    CFG%choice_mask_noice                        = choice_mask_noice_config

    ! Partially fixed geometry, useful for initialisation and inversion runs
    CFG%fixed_shelf_geometry                     = fixed_shelf_geometry_config
    CFG%fixed_sheet_geometry                     = fixed_sheet_geometry_config
    CFG%fixed_grounding_line                     = fixed_grounding_line_config

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE copy_ice_config_variables_to_struct

END MODULE ice_configuration
