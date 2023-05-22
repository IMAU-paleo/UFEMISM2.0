MODULE model_configuration

  ! The different parameters that control the simulation.
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
  ! NOTE: since UFEMISM 2.0, config files should list ALL config variables. This means the
  !       default values in this module are now only for illustration, and are not used
  !       anymore, so that the config file completely determines the model behaviour.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string, &
                                                                     capitalise_string, remove_leading_spaces

  IMPLICIT NONE

! ===== Configuration variables =====
! ===================================

  ! The "_config" variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!

  ! General model instructions
  ! ==========================

    ! Output directory
    LOGICAL             :: create_procedural_output_dir_config          = .TRUE.                           ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)
    CHARACTER(LEN=256)  :: fixed_output_dir_config                      = 'results_UFEMISM'                ! If not, create a directory with this name instead (stops the program if this directory already exists)

    ! Debugging
    LOGICAL             :: do_write_debug_data_config                   = .FALSE.                          ! Whether or not the debug NetCDF file should be created and written to
    LOGICAL             :: do_check_for_NaN_config                      = .FALSE.                          ! Whether or not fields should be checked for NaN values
    LOGICAL             :: do_time_display_config                       = .FALSE.                          ! Print current model time to screen

  ! == Time steps and range
  ! =======================

    REAL(dp)            :: start_time_of_run_config                     = 0.0_dp                           ! [yr] Start time of the simulations
    REAL(dp)            :: end_time_of_run_config                       = 0.0_dp                           ! [yr] End   time of the simulations
    REAL(dp)            :: dt_coupling_config                           = 100._dp                          ! [yr] Interval of coupling between the different model regions

  ! == Which model regions to simulate
  ! ==================================

    LOGICAL             :: do_NAM_config                                = .FALSE.
    LOGICAL             :: do_EAS_config                                = .FALSE.
    LOGICAL             :: do_GRL_config                                = .FALSE.
    LOGICAL             :: do_ANT_config                                = .FALSE.

  ! Do only unit tests
  ! ==================
    
    logical             :: do_unit_tests_config                         = .false.

  ! == The four model regions
  ! =========================

    ! North America
    REAL(dp)            :: lambda_M_NAM_config                          = 265._dp                          ! Longitude of the pole of the stereographic projection for the North America domain [degrees east]
    REAL(dp)            :: phi_M_NAM_config                             = 62._dp                           ! Latitude  of the pole of the stereographic projection for the North America domain [degrees north]
    REAL(dp)            :: beta_stereo_NAM_config                       = 71._dp                           ! Standard parallel     of the stereographic projection for the North America domain [degrees]
    REAL(dp)            :: xmin_NAM_config                              = -3600000._dp                     ! Western  boundary     of the North America domain [m]
    REAL(dp)            :: xmax_NAM_config                              =  3600000._dp                     ! Eastern  boundary     of the North America domain [m]
    REAL(dp)            :: ymin_NAM_config                              = -2400000._dp                     ! Southern boundary     of the North America domain [m]
    REAL(dp)            :: ymax_NAM_config                              =  2400000._dp                     ! Northern boundary     of the North America domain [m]

    ! Eurasia
    REAL(dp)            :: lambda_M_EAS_config                          = 40._dp                           ! Longitude of the pole of the stereographic projection for the Eurasia domain [degrees east]
    REAL(dp)            :: phi_M_EAS_config                             = 70._dp                           ! Latitude  of the pole of the stereographic projection for the Eurasia domain [degrees north]
    REAL(dp)            :: beta_stereo_EAS_config                       = 71._dp                           ! Standard parallel     of the stereographic projection for the Eurasia domain [degrees]
    REAL(dp)            :: xmin_EAS_config                              = -3400000._dp                     ! Western  boundary     of the Eurasia domain [m]
    REAL(dp)            :: xmax_EAS_config                              =  3400000._dp                     ! Eastern  boundary     of the Eurasia domain [m]
    REAL(dp)            :: ymin_EAS_config                              = -2080000._dp                     ! Southern boundary     of the Eurasia domain [m]
    REAL(dp)            :: ymax_EAS_config                              =  2080000._dp                     ! Northern boundary     of the Eurasia domain [m]

    ! Greenland
    REAL(dp)            :: lambda_M_GRL_config                          = -45._dp                          ! Longitude of the pole of the stereographic projection for the Greenland domain [degrees east]
    REAL(dp)            :: phi_M_GRL_config                             = 90._dp                           ! Latitude  of the pole of the stereographic projection for the Greenland domain [degrees north]
    REAL(dp)            :: beta_stereo_GRL_config                       = 70._dp                           ! Standard parallel     of the stereographic projection for the Greenland domain [degrees]
    REAL(dp)            :: xmin_GRL_config                              =  -720000._dp                     ! Western  boundary     of the Greenland domain [m]
    REAL(dp)            :: xmax_GRL_config                              =   960000._dp                     ! Eastern  boundary     of the Greenland domain [m]
    REAL(dp)            :: ymin_GRL_config                              = -3450000._dp                     ! Southern boundary     of the Greenland domain [m]
    REAL(dp)            :: ymax_GRL_config                              =  -570000._dp                     ! Northern boundary     of the Greenland domain [m]

    ! Antarctica
    REAL(dp)            :: lambda_M_ANT_config                          = 0._dp                            ! Longitude of the pole of the stereographic projection for the Antarctica domain [degrees east]
    REAL(dp)            :: phi_M_ANT_config                             = -90._dp                          ! Latitude  of the pole of the stereographic projection for the Antarctica domain [degrees north]
    REAL(dp)            :: beta_stereo_ANT_config                       = 71._dp                           ! Standard parallel     of the stereographic projection for the Antarctica domain [degrees]
    REAL(dp)            :: xmin_ANT_config                              = -3040000._dp                     ! Western  boundary     of the Antarctica domain [m]
    REAL(dp)            :: xmax_ANT_config                              =  3040000._dp                     ! Eastern  boundary     of the Antarctica domain [m]
    REAL(dp)            :: ymin_ANT_config                              = -3040000._dp                     ! Southern boundary     of the Antarctica domain [m]
    REAL(dp)            :: ymax_ANT_config                              =  3040000._dp                     ! Northern boundary     of the Antarctica domain [m]

  ! == Reference geometries (initial, present-day, and GIA equilibrium)
  ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    REAL(dp)            :: refgeo_Hi_min_config                         = 2.0_dp                           ! Remove ice thinner than this value in the reference ice geometry. Particularly useful for BedMachine Greenland, which somehow covers the entire tundra with half a meter of ice...
    LOGICAL             :: remove_Lake_Vostok_config                    = .TRUE.

    ! == Initial geometry
    ! ===================

    CHARACTER(LEN=256)  :: choice_refgeo_init_NAM_config                = 'read_from_file'                 ! Choice of present-day geometry for North America; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_init_EAS_config                = 'read_from_file'                 ! Choice of present-day geometry for Eurasia      ; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_init_GRL_config                = 'read_from_file'                 ! Choice of present-day geometry for Greenland    ; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_init_ANT_config                = 'read_from_file'                 ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    CHARACTER(LEN=256)  :: choice_refgeo_init_idealised_config          = 'flatearth'                      ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    REAL(dp)            :: dx_refgeo_init_idealised_config              = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry
    ! Path to file containing present-day geometry when choice_refgeo_init == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_refgeo_init_NAM_config              = 'data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_init_EAS_config              = 'data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_init_GRL_config              = 'data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_init_ANT_config              = 'data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_refgeo_init_NAM_config             = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    REAL(dp)            :: timeframe_refgeo_init_EAS_config             = 1E9_dp
    REAL(dp)            :: timeframe_refgeo_init_GRL_config             = 1E9_dp
    REAL(dp)            :: timeframe_refgeo_init_ANT_config             = 1E9_dp

    ! == Present-day geometry
    ! =======================

    CHARACTER(LEN=256)  :: choice_refgeo_PD_NAM_config                  = 'read_from_file'                 ! Choice of present-day geometry for North America; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_PD_EAS_config                  = 'read_from_file'                 ! Choice of present-day geometry for Eurasia      ; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_PD_GRL_config                  = 'read_from_file'                 ! Choice of present-day geometry for Greenland    ; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_PD_ANT_config                  = 'read_from_file'                 ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    CHARACTER(LEN=256)  :: choice_refgeo_PD_idealised_config            = 'flatearth'                      ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    REAL(dp)            :: dx_refgeo_PD_idealised_config                = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry
    ! Path to file containing present-day geometry when choice_refgeo_PD == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_NAM_config                = 'data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_EAS_config                = 'data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_GRL_config                = 'data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_ANT_config                = 'data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_refgeo_PD_NAM_config               = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    REAL(dp)            :: timeframe_refgeo_PD_EAS_config               = 1E9_dp
    REAL(dp)            :: timeframe_refgeo_PD_GRL_config               = 1E9_dp
    REAL(dp)            :: timeframe_refgeo_PD_ANT_config               = 1E9_dp

    ! == GIA equilibrium geometry
    ! ===========================

    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_NAM_config               = 'read_from_file'                 ! Choice of present-day geometry for North America; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_EAS_config               = 'read_from_file'                 ! Choice of present-day geometry for Eurasia      ; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_GRL_config               = 'read_from_file'                 ! Choice of present-day geometry for Greenland    ; can be "idealised", or "read_from_file"
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_ANT_config               = 'read_from_file'                 ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_idealised_config         = 'flatearth'                      ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    REAL(dp)            :: dx_refgeo_GIAeq_idealised_config             = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry
    ! Path to file containing present-day geometry when choice_refgeo_GIAeq == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_NAM_config             = 'data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_EAS_config             = 'data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_GRL_config             = 'data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_ANT_config             = 'data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_refgeo_GIAeq_NAM_config            = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    REAL(dp)            :: timeframe_refgeo_GIAeq_EAS_config            = 1E9_dp
    REAL(dp)            :: timeframe_refgeo_GIAeq_GRL_config            = 1E9_dp
    REAL(dp)            :: timeframe_refgeo_GIAeq_ANT_config            = 1E9_dp

    ! == Parameters for idealised geometries
    ! ======================================

    REAL(dp)            :: refgeo_idealised_slabonaslope_Hi_config      = 0._dp                            ! Suggested value: 2000 m
    REAL(dp)            :: refgeo_idealised_slabonaslope_dhdx_config    = 1._dp                            ! Suggested value: -0.001
    REAL(dp)            :: refgeo_idealised_Halfar_H0_config            = 0._dp                            ! Suggested value: 3000 m
    REAL(dp)            :: refgeo_idealised_Halfar_R0_config            = 0._dp                            ! Suggested value: 500E3 m
    REAL(dp)            :: refgeo_idealised_Bueler_H0_config            = 0._dp                            ! Suggested value: 3000 m
    REAL(dp)            :: refgeo_idealised_Bueler_R0_config            = 0._dp                            ! Suggested value: 500E3 m
    REAL(dp)            :: refgeo_idealised_Bueler_lambda_config        = 0._dp                            ! Suggested value: 5.0
    REAL(dp)            :: refgeo_idealised_SSA_icestream_Hi_config     = 0._dp                            ! Suggested value: 2000 m
    REAL(dp)            :: refgeo_idealised_SSA_icestream_dhdx_config   = 1._dp                            ! Suggested value: -0.001
    REAL(dp)            :: refgeo_idealised_MISMIP_mod_Hi_init_config   = -1._dp                           ! Suggested value: 100 m
    REAL(dp)            :: refgeo_idealised_ISMIP_HOM_L_config          = 0._dp                            ! Suggested value: 5E3 - 160E3 m
    REAL(dp)            :: refgeo_idealised_MISMIPplus_Hi_init_config   = -1._dp                           ! Suggested value: 100 m

  ! == Mesh generation
  ! ==================

    ! How to set up the initial mesh
    CHARACTER(LEN=256)  :: choice_initial_mesh_NAM_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'
    CHARACTER(LEN=256)  :: choice_initial_mesh_EAS_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'
    CHARACTER(LEN=256)  :: choice_initial_mesh_GRL_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'
    CHARACTER(LEN=256)  :: choice_initial_mesh_ANT_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'

    ! Paths to files containing initial meshes, if choice_initial_mesh == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_initial_mesh_NAM_config             = ''
    CHARACTER(LEN=256)  :: filename_initial_mesh_EAS_config             = ''
    CHARACTER(LEN=256)  :: filename_initial_mesh_GRL_config             = ''
    CHARACTER(LEN=256)  :: filename_initial_mesh_ANT_config             = ''

    ! Resolutions for different parts of the ice sheet
    REAL(dp)            :: maximum_resolution_uniform_config            = 800e3_dp                         ! [m]          Maximum resolution for the entire domain
    REAL(dp)            :: maximum_resolution_grounded_ice_config       = 400e3_dp                         ! [m]          Maximum resolution for grounded ice
    REAL(dp)            :: maximum_resolution_floating_ice_config       = 200e3_dp                         ! [m]          Maximum resolution for floating ice
    REAL(dp)            :: maximum_resolution_grounding_line_config     = 100e3_dp                         ! [m]          Maximum resolution for the grounding line
    REAL(dp)            :: grounding_line_width_config                  = 200e3_dp                         ! [m]          Width of the band around the grounding line that should get this resolution
    REAL(dp)            :: maximum_resolution_calving_front_config      = 100e3_dp                         ! [m]          Maximum resolution for the calving front
    REAL(dp)            :: calving_front_width_config                   = 100e3_dp                         ! [m]          Width of the band around the calving front that should get this resolution
    REAL(dp)            :: maximum_resolution_ice_front_config          = 100e3_dp                         ! [m]          Maximum resolution for the ice front
    REAL(dp)            :: ice_front_width_config                       = 200e3_dp                         ! [m]          Width of the band around the ice front that should get this resolution
    REAL(dp)            :: maximum_resolution_coastline_config          = 100e3_dp                         ! [m]          Maximum resolution for the coastline
    REAL(dp)            :: coastline_width_config                       = 200e3_dp                         ! [m]          Width of the band around the coastline that should get this resolution

    ! Regions of interest
    CHARACTER(LEN=256)  :: choice_regions_of_interest_config            = ''                               ! Regions of interest where other (higher) resolutions apply. Separated by double vertical bars "||", e.g. "PineIsland||Thwaites"
    REAL(dp)            :: ROI_maximum_resolution_uniform_config        = 100e3_dp                         ! [m]          Maximum resolution for the entire domain
    REAL(dp)            :: ROI_maximum_resolution_grounded_ice_config   = 50e3_dp                          ! [m]          Maximum resolution for grounded ice
    REAL(dp)            :: ROI_maximum_resolution_floating_ice_config   = 20e3_dp                          ! [m]          Maximum resolution for floating ice
    REAL(dp)            :: ROI_maximum_resolution_grounding_line_config = 5e3_dp                           ! [m]          Maximum resolution for the grounding line
    REAL(dp)            :: ROI_grounding_line_width_config              = 5e3_dp                           ! [m]          Width of the band around the grounding line that should get this resolution
    REAL(dp)            :: ROI_maximum_resolution_calving_front_config  = 10e3_dp                          ! [m]          Maximum resolution for the calving front
    REAL(dp)            :: ROI_calving_front_width_config               = 10e3_dp                          ! [m]          Width of the band around the calving front that should get this resolution
    REAL(dp)            :: ROI_maximum_resolution_ice_front_config      = 20e3_dp                          ! [m]          Maximum resolution for the ice front
    REAL(dp)            :: ROI_ice_front_width_config                   = 20e3_dp                          ! [m]          Width of the band around the ice front that should get this resolution
    REAL(dp)            :: ROI_maximum_resolution_coastline_config      = 50e3_dp                          ! [m]          Maximum resolution for the coastline
    REAL(dp)            :: ROI_coastline_width_config                   = 50e3_dp                          ! [m]          Width of the band around the coastline that should get this resolution

    ! Advanced geometry parameters
    LOGICAL             :: do_singlecore_mesh_creation_config           = .TRUE.                           !              Whether or not to use only a single core for mesh generation (for better reproducibility)
    REAL(dp)            :: alpha_min_config                             = 0.4363_dp                        ! [radians]    Smallest allowed internal triangle angle
    INTEGER             :: nit_Lloyds_algorithm_config                  = 3                                ! [-]          Number of iterations of Lloyds algorithm to be applied after refinement

    ! Memory
    INTEGER             :: nC_mem_config                                = 32                               ! [-]          How many columns of memory should be allocated for connectivity lists

  ! == The scaled vertical coordinate zeta
  ! ======================================

    CHARACTER(LEN=256)  :: choice_zeta_grid_config                      = 'regular'                        ! The type of vertical grid to use; can be "regular", "irregular_log", "old_15_layer_zeta"
    INTEGER             :: nz_config                                    = 12                               ! The number of vertical layers to use
    REAL(dp)            :: zeta_irregular_log_R_config                  = 10._dp                           ! Ratio between surface and base layer spacings

  ! == Ice dynamics - velocity
  ! ==========================

    ! General
    CHARACTER(LEN=256)  :: choice_stress_balance_approximation_config   = 'SIA/SSA'                        ! Choice of stress balance approximation: "none" (= no flow, though geometry can still change due to mass balance), "SIA", "SSA", "SIA/SSA", "DIVA", "BPA"
    REAL(dp)            :: n_flow_config                                = 3.0_dp                           ! Exponent in Glen's flow law
    REAL(dp)            :: m_enh_sheet_config                           = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
    REAL(dp)            :: m_enh_shelf_config                           = 1.0_dp                           ! Ice flow enhancement factor for floating ice
    CHARACTER(LEN=256)  :: choice_hybrid_SIASSA_scheme_config           = 'add'                            ! Choice of scheme for combining SIA and SSA velocities in the hybrid approach
    LOGICAL             :: do_GL_subgrid_friction_config                = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
    REAL(dp)            :: subgrid_friction_exponent_config             = 2._dp                            ! Exponent to which f_grnd should be raised before being used to scale beta

    ! Some parameters for numerically solving the stress balance
    REAL(dp)            :: SIA_maximum_diffusivity_config               = 1E5_dp                           ! Limit the diffusivity in the SIA to this value
    REAL(dp)            :: visc_it_norm_dUV_tol_config                  = 1E-2_dp                          ! Stop criterion for the viscosity iteration: the L2-norm of successive velocity solutions should be smaller than this number
    INTEGER             :: visc_it_nit_config                           = 50                               ! Maximum number of effective viscosity iterations
    REAL(dp)            :: visc_it_relax_config                         = 0.4_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    REAL(dp)            :: epsilon_sq_0_config                          = 1E-15_dp                         ! Normalisation term so that zero velocity gives non-zero viscosity
    REAL(dp)            :: visc_eff_min_config                          = 1E3_dp                           ! Minimum value for effective viscosity
    REAL(dp)            :: beta_max_config                              = 1E20_dp                          ! Maximum value for basal friction coefficient
    REAL(dp)            :: vel_max_config                               = 5000._dp                         ! Velocities are limited to this value
    REAL(dp)            :: stress_balance_PETSc_rtol_config             = 1E-2_dp                          ! PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    REAL(dp)            :: stress_balance_PETSc_abstol_config           = 1.0_dp                           ! PETSc solver - stop criterion, absolute difference

  ! == Ice dynamics - sliding law
  ! =============================

    ! Sliding laws
    CHARACTER(LEN=256)  :: choice_sliding_law_config                    = 'Zoet-Iverson'                   ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Budd", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"
    CHARACTER(LEN=256)  :: choice_idealised_sliding_law_config          = ''                               ! "ISMIP_HOM_C", "ISMIP_HOM_D", "ISMIP_HOM_E", "ISMIP_HOM_F"
    REAL(dp)            :: slid_delta_v_config                          = 1.0E-3_dp                        ! Normalisation parameter to prevent errors when velocity is zero
    REAL(dp)            :: slid_Weertman_m_config                       = 3._dp                            ! Exponent in Weertman sliding law
    REAL(dp)            :: slid_Budd_q_plastic_config                   = 0.3_dp                           ! Scaling exponent   in Budd sliding law
    REAL(dp)            :: slid_Budd_u_threshold_config                 = 100._dp                          ! Threshold velocity in Budd sliding law
    REAL(dp)            :: slid_ZI_ut_config                            = 200._dp                          ! (uniform) transition velocity used in the Zoet-Iverson sliding law [m/yr]
    REAL(dp)            :: slid_ZI_p_config                             = 5._dp                            ! Velocity exponent             used in the Zoet-Iverson sliding law

  ! == Ice dynamics - boundary conditions
  ! =====================================

    CHARACTER(LEN=256)  :: BC_u_west_config                             = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain border
    CHARACTER(LEN=256)  :: BC_u_east_config                             = 'infinite'                       ! Allowed choices: "infinite", "zero", "periodic_ISMIP-HOM"
    CHARACTER(LEN=256)  :: BC_u_south_config                            = 'infinite'
    CHARACTER(LEN=256)  :: BC_u_north_config                            = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_west_config                             = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_east_config                             = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_south_config                            = 'infinite'
    CHARACTER(LEN=256)  :: BC_v_north_config                            = 'infinite'
    CHARACTER(LEN=256)  :: BC_H_west_config                             = 'zero'                           ! Boundary conditions for ice thickness at the domain boundary
    CHARACTER(LEN=256)  :: BC_H_east_config                             = 'zero'                           ! Allowed choices:  "infinite", "zero", "ISMIP_HOM_F"
    CHARACTER(LEN=256)  :: BC_H_south_config                            = 'zero'
    CHARACTER(LEN=256)  :: BC_H_north_config                            = 'zero'

  ! == Ice dynamics - time integration
  ! ==================================

    CHARACTER(LEN=256)  :: choice_timestepping_config                   = 'pc'                             ! Choice of timestepping method: "direct", "pc" (NOTE: 'direct' does not work with DIVA ice dynamcis!)
    CHARACTER(LEN=256)  :: choice_ice_integration_method_config         = 'explicit'                       ! Choice of ice thickness integration scheme: "none" (i.e. unchanging geometry), "explicit", "semi-implicit"
    REAL(dp)            :: dHi_PETSc_rtol_config                        = 0.001_dp                         ! dHi PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    REAL(dp)            :: dHi_PETSc_abstol_config                      = 0.001_dp                         ! dHi PETSc solver - stop criterion, absolute difference

    ! Predictor-corrector ice-thickness update
    REAL(dp)            :: pc_epsilon_config                            = 3._dp                            ! Target truncation error in dHi_dt [m/yr] (epsilon in Robinson et al., 2020, Eq. 33)
    REAL(dp)            :: pc_k_I_config                                = 0.2_dp                           ! Exponent k_I in  Robinson et al., 2020, Eq. 33
    REAL(dp)            :: pc_k_p_config                                = 0.2_dp                           ! Exponent k_p in  Robinson et al., 2020, Eq. 33
    REAL(dp)            :: pc_eta_min_config                            = 1E-8_dp                          ! Normalisation term in estimation of the truncation error (Robinson et al., Eq. 32)

  ! == Ice dynamics - calving
  ! =========================

    CHARACTER(LEN=256)  :: choice_calving_law_config                    = 'threshold_thickness'            ! Choice of calving law: "none", "threshold_thickness"
    REAL(dp)            :: calving_threshold_thickness_shelf_config     = 100._dp                          ! Threshold ice thickness for ice shelf calving front in the "threshold_thickness" calving law
    REAL(dp)            :: calving_threshold_thickness_sheet_config     = 0._dp                            ! Threshold ice thickness for ice sheet calving front in the "threshold_thickness" calving law
    INTEGER             :: max_calving_rounds_config                    = 20                               ! Maximum number of calving loops during chain reaction
    LOGICAL             :: do_remove_shelves_config                     = .FALSE.                          ! If set to TRUE, all floating ice is always instantly removed (used in the ABUMIP-ABUK experiment)
    LOGICAL             :: remove_shelves_larger_than_PD_config         = .FALSE.                          ! If set to TRUE, all floating ice beyond the present-day calving front is removed (used for some Antarctic spin-ups)
    LOGICAL             :: continental_shelf_calving_config             = .FALSE.                          ! If set to TRUE, all ice beyond the continental shelf edge (set by a maximum depth) is removed
    REAL(dp)            :: continental_shelf_min_height_config          = -2000._dp                        ! Maximum depth of the continental shelf

  ! == Ice dynamics - stabilisation
  ! ===============================

    CHARACTER(LEN=256)  :: choice_mask_noice_config                     = 'none'                           ! Choice of mask_noice configuration

    ! Partially fixed geometry, useful for initialisation and inversion runs
    LOGICAL             :: fixed_shelf_geometry_config                  = .FALSE.                          ! Keep geometry of floating ice fixed
    LOGICAL             :: fixed_sheet_geometry_config                  = .FALSE.                          ! Keep geometry of grounded ice fixed
    LOGICAL             :: fixed_grounding_line_config                  = .FALSE.                          ! Keep ice thickness at the grounding line fixed

  ! == Basal hydrology
  ! ==================

    ! Basal hydrology
    CHARACTER(LEN=256)  :: choice_basal_hydrology_config                = 'Martin2011'                     ! Choice of basal hydrology model: "saturated", "Martin2011"
    REAL(dp)            :: Martin2011_hydro_Hb_min_config               = 0._dp                            ! Martin et al. (2011) basal hydrology model: low-end  Hb  value of bedrock-dependent pore-water pressure
    REAL(dp)            :: Martin2011_hydro_Hb_max_config               = 1000._dp                         ! Martin et al. (2011) basal hydrology model: high-end Hb  value of bedrock-dependent pore-water pressure


  ! == Bed roughness
  ! ==================

    CHARACTER(LEN=256)  :: choice_bed_roughness_config                  = 'uniform'                        ! "uniform", "parameterised", "read_from_file"
    CHARACTER(LEN=256)  :: choice_bed_roughness_parameterised_config    = 'Martin2011'                     ! "Martin2011", "SSA_icestream", "MISMIP+", "BIVMIP_A", "BIVMIP_B", "BIVMIP_C"
    ! Paths to files containing bed roughness fields for the chosen sliding law
    CHARACTER(LEN=256)  :: filename_bed_roughness_NAM_config            = ''
    CHARACTER(LEN=256)  :: filename_bed_roughness_EAS_config            = ''
    CHARACTER(LEN=256)  :: filename_bed_roughness_GRL_config            = ''
    CHARACTER(LEN=256)  :: filename_bed_roughness_ANT_config            = ''
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_bed_roughness_NAM_config           = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    REAL(dp)            :: timeframe_bed_roughness_EAS_config           = 1E9_dp
    REAL(dp)            :: timeframe_bed_roughness_GRL_config           = 1E9_dp
    REAL(dp)            :: timeframe_bed_roughness_ANT_config           = 1E9_dp
    ! Values for uniform bed roughness
    REAL(dp)            :: slid_Weertman_beta_sq_uniform_config         = 1.0E4_dp                         ! Uniform value for beta_sq  in Weertman sliding law
    REAL(dp)            :: slid_Coulomb_phi_fric_uniform_config         = 15._dp                           ! Uniform value for phi_fric in (regularised) Coulomb sliding law
    REAL(dp)            :: slid_Tsai2015_alpha_sq_uniform_config        = 0.5_dp                           ! Uniform value for alpha_sq in the Tsai2015 sliding law
    REAL(dp)            :: slid_Tsai2015_beta_sq_uniform_config         = 1.0E4_dp                         ! Uniform value for beta_sq  in the Tsai2015 sliding law
    REAL(dp)            :: slid_Schoof2005_alpha_sq_uniform_config      = 0.5_dp                           ! Uniform value for alpha_sq in the Schoof2005 sliding law
    REAL(dp)            :: slid_Schoof2005_beta_sq_uniform_config       = 1.0E4_dp                         ! Uniform value for beta_sq  in the Schoof2005 sliding law
    REAL(dp)            :: slid_ZI_phi_fric_uniform_config              = 15._dp                           ! Uniform value for phi_fric in the Zoet-Iverson sliding law
    ! Parameters for bed roughness parameterisations
    REAL(dp)            :: Martin2011till_phi_Hb_min_config             = -1000._dp                        ! Martin et al. (2011) bed roughness parameterisation: low-end  Hb  value of bedrock-dependent till friction angle
    REAL(dp)            :: Martin2011till_phi_Hb_max_config             = 0._dp                            ! Martin et al. (2011) bed roughness parameterisation: high-end Hb  value of bedrock-dependent till friction angle
    REAL(dp)            :: Martin2011till_phi_min_config                = 5._dp                            ! Martin et al. (2011) bed roughness parameterisation: low-end  phi value of bedrock-dependent till friction angle
    REAL(dp)            :: Martin2011till_phi_max_config                = 20._dp                           ! Martin et al. (2011) bed roughness parameterisation: high-end phi value of bedrock-dependent till friction angle

  ! == Bed roughness inversion by nudging
  ! =====================================

    ! General
    LOGICAL             :: do_bed_roughness_nudging_config              = .FALSE.                          ! Whether or not to budge the basal roughness
    REAL(dp)            :: bed_roughness_nudging_t_start_config         = -9.9E9_dp                        ! Earliest model time when nudging is allowed
    REAL(dp)            :: bed_roughness_nudging_t_end_config           = +9.9E9_dp                        ! Latest   model time when nudging is allowed
    CHARACTER(LEN=256)  :: choice_bed_roughness_nudging_method_config   = ''                               ! Choice of bed roughness nudging method
    REAL(dp)            :: bed_roughness_nudging_dt_config              = 5._dp                            ! Time step for bed roughness updates
    CHARACTER(LEN=256)  :: filename_inverted_bed_roughness_config       = 'bed_roughness_inv.nc'           ! NetCDF file where the final inverted bed roughness will be saved

    ! Parameters
    LOGICAL             :: BIVgeo_Bernales_do_smooth_config             = .FALSE.                          ! If set to TRUE, inverted basal roughness is smoothed
    REAL(dp)            :: BIVgeo_Bernales_scale_config                 = 10000._dp                        ! Scaling constant for inversion procedure [m]
    REAL(dp)            :: BIVgeo_Bernales_rsmooth_config               = 500._dp                          ! Smoothing radius for inversion procedure [m]
    REAL(dp)            :: BIVgeo_Bernales_wsmooth_config               = .01_dp                           ! Weight given to the smoothed roughness (1  = full smoothing applied)
    REAL(dp)            :: BIVgeo_Bernales_phi_min_config               = 2._dp                            ! Minimum value of phi_fric allowed during inversion
    REAL(dp)            :: BIVgeo_Bernales_phi_max_config               = 30._dp                           ! Maximum value of phi_fric allowed during inversion
    REAL(dp)            :: BIVgeo_Bernales_tol_diff_config              = 100._dp                          ! Minimum ice thickness difference [m] that triggers inversion (.OR. &)
    REAL(dp)            :: BIVgeo_Bernales_tol_frac_config              = 1.0_dp                           ! Minimum ratio between ice thickness difference and reference value that triggers inversion
    REAL(dp)            :: BIVgeo_Berends2022_tauc_config               = 10._dp                           ! Timescale       in the Berends2022 geometry/velocity-based basal inversion method [yr]
    REAL(dp)            :: BIVgeo_Berends2022_H0_config                 = 100._dp                          ! First  thickness scale in the Berends2022 geometry/velocity-based basal inversion method [m]
    REAL(dp)            :: BIVgeo_Berends2022_u0_config                 = 250._dp                          ! First  velocity  scale in the Berends2022 geometry/velocity-based basal inversion method [m/yr]
    REAL(dp)            :: BIVgeo_Berends2022_Hi_scale_config           = 300._dp                          ! Second thickness scale in the Berends2022 geometry/velocity-based basal inversion method [m]
    REAL(dp)            :: BIVgeo_Berends2022_u_scale_config            = 3000._dp                         ! Second velocity  scale in the Berends2022 geometry/velocity-based basal inversion method [m/yr]
    REAL(dp)            :: BIVgeo_Berends2022_phimin_config             = 0.1_dp                           ! Smallest allowed value for the inverted till friction angle phi
    REAL(dp)            :: BIVgeo_Berends2022_phimax_config             = 30._dp                           ! Largest  allowed value for the inverted till friction angle phi
    CHARACTER(LEN=256)  :: BIVgeo_target_velocity_filename_config       = ''                               ! NetCDF file where the target velocities are read in the CISM+ and Berends2022 geometry/velocity-based basal inversion methods

  ! == Thermodynamics and rheology
  ! ==============================

    ! Initial temperature profile
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_NAM_config    = 'Robin'                          ! Choice of initial ice temperature profile: "uniform", "linear", "Robin", "read_from_file"
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_EAS_config    = 'Robin'
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_GRL_config    = 'Robin'
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_ANT_config    = 'Robin'
    ! Uniform initial ice temperature, if choice_initial_ice_temperature == 'uniform'
    REAL(dp)            :: uniform_initial_ice_temperature_NAM_config   = 270._dp                          ! Uniform ice temperature (applied when choice_initial_ice_temperature_config = "uniform")
    REAL(dp)            :: uniform_initial_ice_temperature_EAS_config   = 270._dp
    REAL(dp)            :: uniform_initial_ice_temperature_GRL_config   = 270._dp
    REAL(dp)            :: uniform_initial_ice_temperature_ANT_config   = 270._dp
    ! Paths to files containing initial temperature fields, if
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_NAM_config  = ''
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_EAS_config  = ''
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_GRL_config  = ''
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_ANT_config  = ''
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_initial_ice_temperature_NAM_config = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    REAL(dp)            :: timeframe_initial_ice_temperature_EAS_config = 1E9_dp
    REAL(dp)            :: timeframe_initial_ice_temperature_GRL_config = 1E9_dp
    REAL(dp)            :: timeframe_initial_ice_temperature_ANT_config = 1E9_dp
    ! Thermodynamical model
    CHARACTER(LEN=256)  :: choice_thermo_model_config                   = '3D_heat_equation'               ! Choice of thermodynamical model: "none", "3D_heat_equation"
    CHARACTER(LEN=256)  :: choice_ice_heat_capacity_config              = 'Pounder1965'                    ! Choice of ice heat capacity model: "uniform", "Pounder1965"
    REAL(dp)            :: uniform_ice_heat_capacity_config             = 2009._dp                         ! Uniform ice heat capacity (applied when choice_ice_heat_capacity_config = "uniform")
    CHARACTER(LEN=256)  :: choice_ice_thermal_conductivity_config       = 'Ritz1987'                       ! Choice of ice heat capacity model: "uniform", "Ritz1987"
    REAL(dp)            :: uniform_ice_thermal_conductivity_config      = 6.626958E7_dp                    ! Uniform ice thermal conductivity (applied when choice_ice_thermal_conductivity_config = "uniform")

    ! Rheological model (relating Glen's flow parameter to ice temperature)
    CHARACTER(LEN=256)  :: choice_ice_rheology_config                   = 'Huybrechts1992'                 ! Choice of ice rheology model: "uniform", "Huybrechts1992", "MISMIP_mod"
    REAL(dp)            :: uniform_flow_factor_config                   = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")

  ! == SELEN
  ! ========

    LOGICAL             :: SELEN_run_at_t_start_config                  = .FALSE.                          ! Whether or not to run SELEN in the first coupling loop (needed for some benchmark experiments)
    INTEGER             :: SELEN_n_TDOF_iterations_config               = 1                                ! Number of Time-Dependent Ocean Function iterations
    INTEGER             :: SELEN_n_recursion_iterations_config          = 1                                ! Number of recursion iterations
    LOGICAL             :: SELEN_use_rotational_feedback_config         = .FALSE.                          ! If TRUE, rotational feedback is included
    INTEGER             :: SELEN_n_harmonics_config                     = 128                              ! Maximum number of harmonic degrees
    LOGICAL             :: SELEN_display_progress_config                = .FALSE.                          ! Whether or not to display the progress of the big loops to the screen (doesn't work on Cartesius!)

    CHARACTER(LEN=256)  :: SELEN_dir_config                             = 'data/SELEN'                     ! Directory where SELEN initial files and spherical harmonics are stored
    CHARACTER(LEN=256)  :: SELEN_global_topo_filename_config            = 'SELEN_global_topography.nc'     ! Filename for the SELEN global topography file (located in SELEN_dir)
    CHARACTER(LEN=256)  :: SELEN_TABOO_init_filename_config             = 'SELEN_TABOO_initial_file.dat'   ! Filename for the TABOO initial file           (idem                )
    CHARACTER(LEN=256)  :: SELEN_LMJ_VALUES_filename_config             = 'SELEN_lmj_values.bin'           ! Filename for the LJ and MJ values file        (idem                )

    INTEGER                  :: SELEN_irreg_time_n_config               = 15                               ! Number of entries in the irregular moving time window
    REAL(dp), DIMENSION(50)  :: SELEN_irreg_time_window_config          = &                                ! Values of entries in the irregular moving time window
   (/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)

    REAL(dp)            :: SELEN_lith_thickness_config                  = 100._dp                         ! Thickness of the elastic lithosphere [km]
    INTEGER             :: SELEN_visc_n_config                          = 3                               ! Number      of viscous asthenosphere layers
    REAL(dp), DIMENSION(3) :: SELEN_visc_prof_config                    = (/ 3._dp, 0.6_dp, 0.3_dp /)     ! Viscosities of viscous asthenosphere layers [?]

    ! Settings for the TABOO Earth deformation model
    INTEGER             :: SELEN_TABOO_CDE_config                       = 0                               ! code of the model (see taboo for explanation)
    INTEGER             :: SELEN_TABOO_TLOVE_config                     = 1                               ! Tidal love numbers yes/no
    INTEGER             :: SELEN_TABOO_DEG1_config                      = 1                               ! Tidal love numbers degree
    REAL(dp)            :: SELEN_TABOO_RCMB_config                      = 3480._dp                        ! Radius of CMB (km)


! ===== The config type =====
! ===========================

  ! The "C" type, which contains all the config parameters as fields.
  ! These will all be overwritten with the values of the "_config" variables,
  ! which are either the default values specified above, are the values
  ! specified from the external config file.

  TYPE type_config
    ! The different parameters that control a UFEMISM simulation

  ! General model instructions
  ! ==========================

    LOGICAL             :: create_procedural_output_dir
    CHARACTER(LEN=256)  :: fixed_output_dir

    ! Debugging
    LOGICAL             :: do_write_debug_data
    LOGICAL             :: do_check_for_NaN
    LOGICAL             :: do_time_display

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

  ! == Switch to do unit tests
  ! ==========================
    
    LOGICAL             :: do_unit_tests
  
  ! == The four model regions
  ! =========================

    ! North America
    REAL(dp)            :: lambda_M_NAM
    REAL(dp)            :: phi_M_NAM
    REAL(dp)            :: beta_stereo_NAM
    REAL(dp)            :: xmin_NAM
    REAL(dp)            :: xmax_NAM
    REAL(dp)            :: ymin_NAM
    REAL(dp)            :: ymax_NAM

    ! Eurasia
    REAL(dp)            :: lambda_M_EAS
    REAL(dp)            :: phi_M_EAS
    REAL(dp)            :: beta_stereo_EAS
    REAL(dp)            :: xmin_EAS
    REAL(dp)            :: xmax_EAS
    REAL(dp)            :: ymin_EAS
    REAL(dp)            :: ymax_EAS

    ! Greenland
    REAL(dp)            :: lambda_M_GRL
    REAL(dp)            :: phi_M_GRL
    REAL(dp)            :: beta_stereo_GRL
    REAL(dp)            :: xmin_GRL
    REAL(dp)            :: xmax_GRL
    REAL(dp)            :: ymin_GRL
    REAL(dp)            :: ymax_GRL

    ! Antarctica
    REAL(dp)            :: lambda_M_ANT
    REAL(dp)            :: phi_M_ANT
    REAL(dp)            :: beta_stereo_ANT
    REAL(dp)            :: xmin_ANT
    REAL(dp)            :: xmax_ANT
    REAL(dp)            :: ymin_ANT
    REAL(dp)            :: ymax_ANT

  ! == Reference geometries (initial, present-day, and GIA equilibrium)
  ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    REAL(dp)            :: refgeo_Hi_min
    LOGICAL             :: remove_Lake_Vostok

    ! == Initial geometry
    ! ===================

    CHARACTER(LEN=256)  :: choice_refgeo_init_NAM
    CHARACTER(LEN=256)  :: choice_refgeo_init_EAS
    CHARACTER(LEN=256)  :: choice_refgeo_init_GRL
    CHARACTER(LEN=256)  :: choice_refgeo_init_ANT
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    CHARACTER(LEN=256)  :: choice_refgeo_init_idealised
    REAL(dp)            :: dx_refgeo_init_idealised
    ! Path to file containing present-day geometry when choice_refgeo_init == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_refgeo_init_NAM
    CHARACTER(LEN=256)  :: filename_refgeo_init_EAS
    CHARACTER(LEN=256)  :: filename_refgeo_init_GRL
    CHARACTER(LEN=256)  :: filename_refgeo_init_ANT
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_refgeo_init_NAM
    REAL(dp)            :: timeframe_refgeo_init_EAS
    REAL(dp)            :: timeframe_refgeo_init_GRL
    REAL(dp)            :: timeframe_refgeo_init_ANT

    ! == Present-day geometry
    ! =======================

    CHARACTER(LEN=256)  :: choice_refgeo_PD_NAM
    CHARACTER(LEN=256)  :: choice_refgeo_PD_EAS
    CHARACTER(LEN=256)  :: choice_refgeo_PD_GRL
    CHARACTER(LEN=256)  :: choice_refgeo_PD_ANT
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    CHARACTER(LEN=256)  :: choice_refgeo_PD_idealised
    REAL(dp)            :: dx_refgeo_PD_idealised
    ! Path to file containing present-day geometry when choice_refgeo_PD == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_NAM
    CHARACTER(LEN=256)  :: filename_refgeo_PD_EAS
    CHARACTER(LEN=256)  :: filename_refgeo_PD_GRL
    CHARACTER(LEN=256)  :: filename_refgeo_PD_ANT
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_refgeo_PD_NAM
    REAL(dp)            :: timeframe_refgeo_PD_EAS
    REAL(dp)            :: timeframe_refgeo_PD_GRL
    REAL(dp)            :: timeframe_refgeo_PD_ANT

    ! == GIA equilibrium geometry
    ! ===========================

    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_NAM
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_EAS
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_GRL
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_ANT
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_idealised
    REAL(dp)            :: dx_refgeo_GIAeq_idealised
    ! Path to file containing present-day geometry when choice_refgeo_GIAeq == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_NAM
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_EAS
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_GRL
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_ANT
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_refgeo_GIAeq_NAM
    REAL(dp)            :: timeframe_refgeo_GIAeq_EAS
    REAL(dp)            :: timeframe_refgeo_GIAeq_GRL
    REAL(dp)            :: timeframe_refgeo_GIAeq_ANT

    ! == Parameters for idealised geometries
    ! ======================================

    REAL(dp)            :: refgeo_idealised_slabonaslope_Hi
    REAL(dp)            :: refgeo_idealised_slabonaslope_dhdx
    REAL(dp)            :: refgeo_idealised_Halfar_H0
    REAL(dp)            :: refgeo_idealised_Halfar_R0
    REAL(dp)            :: refgeo_idealised_Bueler_H0
    REAL(dp)            :: refgeo_idealised_Bueler_R0
    REAL(dp)            :: refgeo_idealised_Bueler_lambda
    REAL(dp)            :: refgeo_idealised_SSA_icestream_Hi
    REAL(dp)            :: refgeo_idealised_SSA_icestream_dhdx
    REAL(dp)            :: refgeo_idealised_MISMIP_mod_Hi_init
    REAL(dp)            :: refgeo_idealised_ISMIP_HOM_L
    REAL(dp)            :: refgeo_idealised_MISMIPplus_Hi_init

  ! == Mesh generation
  ! ==================

    ! How to set up the initial mesh
    CHARACTER(LEN=256)  :: choice_initial_mesh_NAM
    CHARACTER(LEN=256)  :: choice_initial_mesh_EAS
    CHARACTER(LEN=256)  :: choice_initial_mesh_GRL
    CHARACTER(LEN=256)  :: choice_initial_mesh_ANT

    ! Paths to files containing initial meshes, if choice_initial_mesh == 'read_from_file'
    CHARACTER(LEN=256)  :: filename_initial_mesh_NAM
    CHARACTER(LEN=256)  :: filename_initial_mesh_EAS
    CHARACTER(LEN=256)  :: filename_initial_mesh_GRL
    CHARACTER(LEN=256)  :: filename_initial_mesh_ANT

    ! Resolutions for different parts of the ice sheet
    REAL(dp)            :: maximum_resolution_uniform
    REAL(dp)            :: maximum_resolution_grounded_ice
    REAL(dp)            :: maximum_resolution_floating_ice
    REAL(dp)            :: maximum_resolution_grounding_line
    REAL(dp)            :: grounding_line_width
    REAL(dp)            :: maximum_resolution_calving_front
    REAL(dp)            :: calving_front_width
    REAL(dp)            :: maximum_resolution_ice_front
    REAL(dp)            :: ice_front_width
    REAL(dp)            :: maximum_resolution_coastline
    REAL(dp)            :: coastline_width

    ! Regions of interest
    CHARACTER(LEN=256)  :: choice_regions_of_interest
    REAL(dp)            :: ROI_maximum_resolution_uniform
    REAL(dp)            :: ROI_maximum_resolution_grounded_ice
    REAL(dp)            :: ROI_maximum_resolution_floating_ice
    REAL(dp)            :: ROI_maximum_resolution_grounding_line
    REAL(dp)            :: ROI_grounding_line_width
    REAL(dp)            :: ROI_maximum_resolution_calving_front
    REAL(dp)            :: ROI_calving_front_width
    REAL(dp)            :: ROI_maximum_resolution_ice_front
    REAL(dp)            :: ROI_ice_front_width
    REAL(dp)            :: ROI_maximum_resolution_coastline
    REAL(dp)            :: ROI_coastline_width

    ! Advanced geometry parameters
    LOGICAL             :: do_singlecore_mesh_creation
    REAL(dp)            :: alpha_min
    INTEGER             :: nit_Lloyds_algorithm

    ! Memory
    INTEGER             :: nC_mem

  ! == The scaled vertical coordinate zeta
  ! ======================================

    CHARACTER(LEN=256)  :: choice_zeta_grid
    INTEGER             :: nz
    REAL(dp)            :: zeta_irregular_log_R

  ! == Ice dynamics - velocity
  ! ==========================

    ! General
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
    REAL(dp)            :: stress_balance_PETSc_rtol
    REAL(dp)            :: stress_balance_PETSc_abstol

  ! == Ice dynamics - sliding law
  ! =============================

    ! Sliding laws
    CHARACTER(LEN=256)  :: choice_sliding_law
    CHARACTER(LEN=256)  :: choice_idealised_sliding_law
    REAL(dp)            :: slid_delta_v
    REAL(dp)            :: slid_Weertman_m
    REAL(dp)            :: slid_Budd_q_plastic
    REAL(dp)            :: slid_Budd_u_threshold
    REAL(dp)            :: slid_ZI_ut
    REAL(dp)            :: slid_ZI_p

  ! == Ice dynamics - boundary conditions
  ! =====================================

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

  ! == Ice dynamics - calving
  ! =========================

    CHARACTER(LEN=256)  :: choice_calving_law
    REAL(dp)            :: calving_threshold_thickness_shelf
    REAL(dp)            :: calving_threshold_thickness_sheet
    INTEGER             :: max_calving_rounds
    LOGICAL             :: do_remove_shelves
    LOGICAL             :: remove_shelves_larger_than_PD
    LOGICAL             :: continental_shelf_calving
    REAL(dp)            :: continental_shelf_min_height

  ! == Ice dynamics - stabilisation
  ! ===============================

    CHARACTER(LEN=256)  :: choice_mask_noice

    ! Partially fixed geometry, useful for initialisation and inversion runs
    LOGICAL             :: fixed_shelf_geometry
    LOGICAL             :: fixed_sheet_geometry
    LOGICAL             :: fixed_grounding_line

  ! == Basal hydrology
  ! ==================

    ! Basal hydrology
    CHARACTER(LEN=256)  :: choice_basal_hydrology
    REAL(dp)            :: Martin2011_hydro_Hb_min
    REAL(dp)            :: Martin2011_hydro_Hb_max

  ! == Bed roughness
  ! ==================

    CHARACTER(LEN=256)  :: choice_bed_roughness
    CHARACTER(LEN=256)  :: choice_bed_roughness_parameterised
    ! Paths to files containing bed roughness fields for the chosen sliding law
    CHARACTER(LEN=256)  :: filename_bed_roughness_NAM
    CHARACTER(LEN=256)  :: filename_bed_roughness_EAS
    CHARACTER(LEN=256)  :: filename_bed_roughness_GRL
    CHARACTER(LEN=256)  :: filename_bed_roughness_ANT
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_bed_roughness_NAM
    REAL(dp)            :: timeframe_bed_roughness_EAS
    REAL(dp)            :: timeframe_bed_roughness_GRL
    REAL(dp)            :: timeframe_bed_roughness_ANT
    ! Values for uniform bed roughness
    REAL(dp)            :: slid_Weertman_beta_sq_uniform
    REAL(dp)            :: slid_Coulomb_phi_fric_uniform
    REAL(dp)            :: slid_Tsai2015_alpha_sq_uniform
    REAL(dp)            :: slid_Tsai2015_beta_sq_uniform
    REAL(dp)            :: slid_Schoof2005_alpha_sq_uniform
    REAL(dp)            :: slid_Schoof2005_beta_sq_uniform
    REAL(dp)            :: slid_ZI_phi_fric_uniform
    ! Parameters for bed roughness parameterisations
    REAL(dp)            :: Martin2011till_phi_Hb_min
    REAL(dp)            :: Martin2011till_phi_Hb_max
    REAL(dp)            :: Martin2011till_phi_min
    REAL(dp)            :: Martin2011till_phi_max

  ! == Bed roughness inversion by nudging
  ! =====================================

    ! General
    LOGICAL             :: do_bed_roughness_nudging
    REAL(dp)            :: bed_roughness_nudging_t_start
    REAL(dp)            :: bed_roughness_nudging_t_end
    CHARACTER(LEN=256)  :: choice_bed_roughness_nudging_method
    REAL(dp)            :: bed_roughness_nudging_dt
    CHARACTER(LEN=256)  :: filename_inverted_bed_roughness

    ! Parameters
    LOGICAL             :: BIVgeo_Bernales_do_smooth
    REAL(dp)            :: BIVgeo_Bernales_scale
    REAL(dp)            :: BIVgeo_Bernales_rsmooth
    REAL(dp)            :: BIVgeo_Bernales_wsmooth
    REAL(dp)            :: BIVgeo_Bernales_phi_min
    REAL(dp)            :: BIVgeo_Bernales_phi_max
    REAL(dp)            :: BIVgeo_Bernales_tol_diff
    REAL(dp)            :: BIVgeo_Bernales_tol_frac
    REAL(dp)            :: BIVgeo_Berends2022_tauc
    REAL(dp)            :: BIVgeo_Berends2022_H0
    REAL(dp)            :: BIVgeo_Berends2022_u0
    REAL(dp)            :: BIVgeo_Berends2022_Hi_scale
    REAL(dp)            :: BIVgeo_Berends2022_u_scale
    REAL(dp)            :: BIVgeo_Berends2022_phimin
    REAL(dp)            :: BIVgeo_Berends2022_phimax
    CHARACTER(LEN=256)  :: BIVgeo_target_velocity_filename

  ! == Thermodynamics and rheology
  ! ==============================

    ! Initial temperature profile
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_NAM
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_EAS
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_GRL
    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_ANT
    ! Uniform initial ice temperature, if choice_initial_ice_temperature == 'uniform'
    REAL(dp)            :: uniform_initial_ice_temperature_NAM
    REAL(dp)            :: uniform_initial_ice_temperature_EAS
    REAL(dp)            :: uniform_initial_ice_temperature_GRL
    REAL(dp)            :: uniform_initial_ice_temperature_ANT
    ! Paths to files containing initial temperature fields, if
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_NAM
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_EAS
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_GRL
    CHARACTER(LEN=256)  :: filename_initial_ice_temperature_ANT
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    REAL(dp)            :: timeframe_initial_ice_temperature_NAM
    REAL(dp)            :: timeframe_initial_ice_temperature_EAS
    REAL(dp)            :: timeframe_initial_ice_temperature_GRL
    REAL(dp)            :: timeframe_initial_ice_temperature_ANT
    ! Thermodynamical model

    CHARACTER(LEN=256)  :: choice_thermo_model
    CHARACTER(LEN=256)  :: choice_ice_heat_capacity
    REAL(dp)            :: uniform_ice_heat_capacity
    CHARACTER(LEN=256)  :: choice_ice_thermal_conductivity
    REAL(dp)            :: uniform_ice_thermal_conductivity

    ! Rheological model (relating Glen's flow parameter to ice temperature)
    CHARACTER(LEN=256)  :: choice_ice_rheology
    REAL(dp)            :: uniform_flow_factor

  ! == SELEN
  ! ========

    LOGICAL             :: SELEN_run_at_t_start
    INTEGER             :: SELEN_n_TDOF_iterations
    INTEGER             :: SELEN_n_recursion_iterations
    LOGICAL             :: SELEN_use_rotational_feedback
    INTEGER             :: SELEN_n_harmonics
    LOGICAL             :: SELEN_display_progress

    CHARACTER(LEN=256)  :: SELEN_dir
    CHARACTER(LEN=256)  :: SELEN_global_topo_filename
    CHARACTER(LEN=256)  :: SELEN_TABOO_init_filename
    CHARACTER(LEN=256)  :: SELEN_LMJ_VALUES_filename

    INTEGER                  :: SELEN_irreg_time_n
    REAL(dp), DIMENSION(50)  :: SELEN_irreg_time_window

    REAL(dp)            :: SELEN_lith_thickness
    INTEGER             :: SELEN_visc_n
    REAL(dp), DIMENSION(3) :: SELEN_visc_prof

    ! Settings for the TABOO Earth deformation model
    INTEGER             :: SELEN_TABOO_CDE
    INTEGER             :: SELEN_TABOO_TLOVE
    INTEGER             :: SELEN_TABOO_DEG1
    REAL(dp)            :: SELEN_TABOO_RCMB

  ! == Non-configurable variables
  ! =============================

    CHARACTER(LEN=256)  :: output_dir

  END TYPE type_config

  ! The main config structure
  TYPE(type_config)   :: C

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

    IF (par%master) WRITE(0,'(A)') ''
    IF (par%master) WRITE(0,'(A)') ' Running UFEMISM with settings from configuration file: ' // colour_string( TRIM( config_filename), 'light blue')

    ! Initialise the main config structure from the config file
    CALL initialise_config_from_file( config_filename)

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
      WRITE(0,'(A)') ''
      WRITE(0,'(A)') ' Output directory: ' // colour_string( TRIM( C%output_dir), 'light blue')

    END IF ! IF (par%master) THEN
    CALL sync

    ! Copy the config file to the output directory
    IF (par%master) THEN
      CALL system('cp ' // config_filename    // ' ' // TRIM( C%output_dir))
    END IF ! IF (master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_model_configuration

  SUBROUTINE initialise_config_from_file( config_filename)
    ! Initialise a config structure from an external config text file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_config_from_file'
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Let each of the processors read the config file in turns so there's no access conflicts
    DO i = 0, par%n-1

      IF (i == par%i) THEN

        ! Read the external file, use a Fortran NAMELIST to overwrite the default
        ! values of the XXX_config variables
        CALL read_config_file( config_filename)

        ! Copy values from the XXX_config variables to the config structure
        CALL copy_config_variables_to_struct

      END IF ! IF (i == par%i) THEN

      ! Make sure only one process at a time reads from / writes to disk
      CALL sync

    END DO ! DO i = 0, par%n-1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_config_from_file

  SUBROUTINE read_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.
    !
    ! Note: make sure that only one process at a time calls this subroutine!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_config_file'
    CHARACTER(LEN=256), PARAMETER                      :: namelist_filename = 'config_namelist_temp.txt'
    INTEGER, PARAMETER                                 :: config_unit    = 1337
    INTEGER, PARAMETER                                 :: namelist_unit  = 1338
    INTEGER                                            :: ios

    ! The NAMELIST that's used to read the external config file.
    NAMELIST /CONFIG/&
      create_procedural_output_dir_config                   , &
      fixed_output_dir_config                               , &
      do_write_debug_data_config                            , &
      do_check_for_NaN_config                               , &
      do_time_display_config                                , &
      start_time_of_run_config                              , &
      end_time_of_run_config                                , &
      dt_coupling_config                                    , &
      do_NAM_config                                         , &
      do_EAS_config                                         , &
      do_GRL_config                                         , &
      do_ANT_config                                         , &
      do_unit_tests_config                                  , &
      lambda_M_NAM_config                                   , &
      phi_M_NAM_config                                      , &
      beta_stereo_NAM_config                                , &
      xmin_NAM_config                                       , &
      xmax_NAM_config                                       , &
      ymin_NAM_config                                       , &
      ymax_NAM_config                                       , &
      lambda_M_EAS_config                                   , &
      phi_M_EAS_config                                      , &
      beta_stereo_EAS_config                                , &
      xmin_EAS_config                                       , &
      xmax_EAS_config                                       , &
      ymin_EAS_config                                       , &
      ymax_EAS_config                                       , &
      lambda_M_GRL_config                                   , &
      phi_M_GRL_config                                      , &
      beta_stereo_GRL_config                                , &
      xmin_GRL_config                                       , &
      xmax_GRL_config                                       , &
      ymin_GRL_config                                       , &
      ymax_GRL_config                                       , &
      lambda_M_ANT_config                                   , &
      phi_M_ANT_config                                      , &
      beta_stereo_ANT_config                                , &
      xmin_ANT_config                                       , &
      xmax_ANT_config                                       , &
      ymin_ANT_config                                       , &
      ymax_ANT_config                                       , &
      refgeo_Hi_min_config                                  , &
      remove_Lake_Vostok_config                             , &
      choice_refgeo_init_NAM_config                         , &
      choice_refgeo_init_EAS_config                         , &
      choice_refgeo_init_GRL_config                         , &
      choice_refgeo_init_ANT_config                         , &
      choice_refgeo_init_idealised_config                   , &
      dx_refgeo_init_idealised_config                       , &
      filename_refgeo_init_NAM_config                       , &
      filename_refgeo_init_EAS_config                       , &
      filename_refgeo_init_GRL_config                       , &
      filename_refgeo_init_ANT_config                       , &
      timeframe_refgeo_init_NAM_config                      , &
      timeframe_refgeo_init_EAS_config                      , &
      timeframe_refgeo_init_GRL_config                      , &
      timeframe_refgeo_init_ANT_config                      , &
      choice_refgeo_PD_NAM_config                           , &
      choice_refgeo_PD_EAS_config                           , &
      choice_refgeo_PD_GRL_config                           , &
      choice_refgeo_PD_ANT_config                           , &
      choice_refgeo_PD_idealised_config                     , &
      dx_refgeo_PD_idealised_config                         , &
      filename_refgeo_PD_NAM_config                         , &
      filename_refgeo_PD_EAS_config                         , &
      filename_refgeo_PD_GRL_config                         , &
      filename_refgeo_PD_ANT_config                         , &
      timeframe_refgeo_PD_NAM_config                        , &
      timeframe_refgeo_PD_EAS_config                        , &
      timeframe_refgeo_PD_GRL_config                        , &
      timeframe_refgeo_PD_ANT_config                        , &
      choice_refgeo_GIAeq_NAM_config                        , &
      choice_refgeo_GIAeq_EAS_config                        , &
      choice_refgeo_GIAeq_GRL_config                        , &
      choice_refgeo_GIAeq_ANT_config                        , &
      choice_refgeo_GIAeq_idealised_config                  , &
      dx_refgeo_GIAeq_idealised_config                      , &
      filename_refgeo_GIAeq_NAM_config                      , &
      filename_refgeo_GIAeq_EAS_config                      , &
      filename_refgeo_GIAeq_GRL_config                      , &
      filename_refgeo_GIAeq_ANT_config                      , &
      timeframe_refgeo_GIAeq_NAM_config                     , &
      timeframe_refgeo_GIAeq_EAS_config                     , &
      timeframe_refgeo_GIAeq_GRL_config                     , &
      timeframe_refgeo_GIAeq_ANT_config                     , &
      refgeo_idealised_slabonaslope_Hi_config               , &
      refgeo_idealised_slabonaslope_dhdx_config             , &
      refgeo_idealised_Halfar_H0_config                     , &
      refgeo_idealised_Halfar_R0_config                     , &
      refgeo_idealised_Bueler_H0_config                     , &
      refgeo_idealised_Bueler_R0_config                     , &
      refgeo_idealised_Bueler_lambda_config                 , &
      refgeo_idealised_SSA_icestream_Hi_config              , &
      refgeo_idealised_SSA_icestream_dhdx_config            , &
      refgeo_idealised_MISMIP_mod_Hi_init_config            , &
      refgeo_idealised_ISMIP_HOM_L_config                   , &
      refgeo_idealised_MISMIPplus_Hi_init_config            , &
      choice_initial_mesh_NAM_config                        , &
      choice_initial_mesh_EAS_config                        , &
      choice_initial_mesh_GRL_config                        , &
      choice_initial_mesh_ANT_config                        , &
      filename_initial_mesh_NAM_config                      , &
      filename_initial_mesh_EAS_config                      , &
      filename_initial_mesh_GRL_config                      , &
      filename_initial_mesh_ANT_config                      , &
      maximum_resolution_uniform_config                     , &
      maximum_resolution_grounded_ice_config                , &
      maximum_resolution_floating_ice_config                , &
      maximum_resolution_grounding_line_config              , &
      grounding_line_width_config                           , &
      maximum_resolution_calving_front_config               , &
      calving_front_width_config                            , &
      maximum_resolution_ice_front_config                   , &
      ice_front_width_config                                , &
      maximum_resolution_coastline_config                   , &
      coastline_width_config                                , &
      choice_regions_of_interest_config                     , &
      ROI_maximum_resolution_uniform_config                 , &
      ROI_maximum_resolution_grounded_ice_config            , &
      ROI_maximum_resolution_floating_ice_config            , &
      ROI_maximum_resolution_grounding_line_config          , &
      ROI_grounding_line_width_config                       , &
      ROI_maximum_resolution_calving_front_config           , &
      ROI_calving_front_width_config                        , &
      ROI_maximum_resolution_ice_front_config               , &
      ROI_ice_front_width_config                            , &
      ROI_maximum_resolution_coastline_config               , &
      ROI_coastline_width_config                            , &
      do_singlecore_mesh_creation_config                    , &
      alpha_min_config                                      , &
      nit_Lloyds_algorithm_config                           , &
      nC_mem_config                                         , &
      choice_zeta_grid_config                               , &
      nz_config                                             , &
      zeta_irregular_log_R_config                           , &
      choice_stress_balance_approximation_config            , &
      n_flow_config                                         , &
      m_enh_sheet_config                                    , &
      m_enh_shelf_config                                    , &
      choice_hybrid_SIASSA_scheme_config                    , &
      do_GL_subgrid_friction_config                         , &
      subgrid_friction_exponent_config                      , &
      SIA_maximum_diffusivity_config                        , &
      visc_it_norm_dUV_tol_config                           , &
      visc_it_nit_config                                    , &
      visc_it_relax_config                                  , &
      epsilon_sq_0_config                                   , &
      visc_eff_min_config                                   , &
      beta_max_config                                       , &
      vel_max_config                                        , &
      stress_balance_PETSc_rtol_config                      , &
      stress_balance_PETSc_abstol_config                    , &
      choice_sliding_law_config                             , &
      choice_idealised_sliding_law_config                   , &
      slid_delta_v_config                                   , &
      slid_Weertman_m_config                                , &
      slid_Budd_q_plastic_config                            , &
      slid_Budd_u_threshold_config                          , &
      slid_ZI_ut_config                                     , &
      slid_ZI_p_config                                      , &
      BC_u_west_config                                      , &
      BC_u_east_config                                      , &
      BC_u_south_config                                     , &
      BC_u_north_config                                     , &
      BC_v_west_config                                      , &
      BC_v_east_config                                      , &
      BC_v_south_config                                     , &
      BC_v_north_config                                     , &
      BC_H_west_config                                      , &
      BC_H_east_config                                      , &
      BC_H_south_config                                     , &
      BC_H_north_config                                     , &
      choice_timestepping_config                            , &
      choice_ice_integration_method_config                  , &
      dHi_PETSc_rtol_config                                 , &
      dHi_PETSc_abstol_config                               , &
      pc_epsilon_config                                     , &
      pc_k_I_config                                         , &
      pc_k_p_config                                         , &
      pc_eta_min_config                                     , &
      choice_calving_law_config                             , &
      calving_threshold_thickness_shelf_config              , &
      calving_threshold_thickness_sheet_config              , &
      max_calving_rounds_config                             , &
      do_remove_shelves_config                              , &
      remove_shelves_larger_than_PD_config                  , &
      continental_shelf_calving_config                      , &
      continental_shelf_min_height_config                   , &
      choice_mask_noice_config                              , &
      fixed_shelf_geometry_config                           , &
      fixed_sheet_geometry_config                           , &
      fixed_grounding_line_config                           , &
      choice_basal_hydrology_config                         , &
      Martin2011_hydro_Hb_min_config                        , &
      Martin2011_hydro_Hb_max_config                        , &
      choice_bed_roughness_config                           , &
      choice_bed_roughness_parameterised_config             , &
      filename_bed_roughness_NAM_config                     , &
      filename_bed_roughness_EAS_config                     , &
      filename_bed_roughness_GRL_config                     , &
      filename_bed_roughness_ANT_config                     , &
      timeframe_bed_roughness_NAM_config                    , &
      timeframe_bed_roughness_EAS_config                    , &
      timeframe_bed_roughness_GRL_config                    , &
      timeframe_bed_roughness_ANT_config                    , &
      slid_Weertman_beta_sq_uniform_config                  , &
      slid_Coulomb_phi_fric_uniform_config                  , &
      slid_Tsai2015_alpha_sq_uniform_config                 , &
      slid_Tsai2015_beta_sq_uniform_config                  , &
      slid_Schoof2005_alpha_sq_uniform_config               , &
      slid_Schoof2005_beta_sq_uniform_config                , &
      slid_ZI_phi_fric_uniform_config                       , &
      Martin2011till_phi_Hb_min_config                      , &
      Martin2011till_phi_Hb_max_config                      , &
      Martin2011till_phi_min_config                         , &
      Martin2011till_phi_max_config                         , &
      do_bed_roughness_nudging_config                       , &
      bed_roughness_nudging_t_start_config                  , &
      bed_roughness_nudging_t_end_config                    , &
      choice_bed_roughness_nudging_method_config            , &
      bed_roughness_nudging_dt_config                       , &
      filename_inverted_bed_roughness_config                , &
      BIVgeo_Bernales_do_smooth_config                      , &
      BIVgeo_Bernales_scale_config                          , &
      BIVgeo_Bernales_rsmooth_config                        , &
      BIVgeo_Bernales_wsmooth_config                        , &
      BIVgeo_Bernales_phi_min_config                        , &
      BIVgeo_Bernales_phi_max_config                        , &
      BIVgeo_Bernales_tol_diff_config                       , &
      BIVgeo_Bernales_tol_frac_config                       , &
      BIVgeo_Berends2022_tauc_config                        , &
      BIVgeo_Berends2022_H0_config                          , &
      BIVgeo_Berends2022_u0_config                          , &
      BIVgeo_Berends2022_Hi_scale_config                    , &
      BIVgeo_Berends2022_u_scale_config                     , &
      BIVgeo_Berends2022_phimin_config                      , &
      BIVgeo_Berends2022_phimax_config                      , &
      BIVgeo_target_velocity_filename_config                , &
      choice_initial_ice_temperature_NAM_config             , &
      choice_initial_ice_temperature_EAS_config             , &
      choice_initial_ice_temperature_GRL_config             , &
      choice_initial_ice_temperature_ANT_config             , &
      uniform_initial_ice_temperature_NAM_config            , &
      uniform_initial_ice_temperature_EAS_config            , &
      uniform_initial_ice_temperature_GRL_config            , &
      uniform_initial_ice_temperature_ANT_config            , &
      filename_initial_ice_temperature_NAM_config           , &
      filename_initial_ice_temperature_EAS_config           , &
      filename_initial_ice_temperature_GRL_config           , &
      filename_initial_ice_temperature_ANT_config           , &
      timeframe_initial_ice_temperature_NAM_config          , &
      timeframe_initial_ice_temperature_EAS_config          , &
      timeframe_initial_ice_temperature_GRL_config          , &
      timeframe_initial_ice_temperature_ANT_config          , &
      choice_thermo_model_config                            , &
      choice_ice_heat_capacity_config                       , &
      uniform_ice_heat_capacity_config                      , &
      choice_ice_thermal_conductivity_config                , &
      uniform_ice_thermal_conductivity_config               , &
      choice_ice_rheology_config                            , &
      uniform_flow_factor_config                            , &
      SELEN_run_at_t_start_config                           , &
      SELEN_n_TDOF_iterations_config                        , &
      SELEN_n_recursion_iterations_config                   , &
      SELEN_use_rotational_feedback_config                  , &
      SELEN_n_harmonics_config                              , &
      SELEN_display_progress_config                         , &
      SELEN_dir_config                                      , &
      SELEN_global_topo_filename_config                     , &
      SELEN_TABOO_init_filename_config                      , &
      SELEN_LMJ_VALUES_filename_config                      , &
      SELEN_irreg_time_n_config                             , &
      SELEN_irreg_time_window_config                        , &
      SELEN_lith_thickness_config                           , &
      SELEN_visc_n_config                                   , &
      SELEN_visc_prof_config                                , &
      SELEN_TABOO_CDE_config                                , &
      SELEN_TABOO_TLOVE_config                              , &
      SELEN_TABOO_DEG1_config                               , &
      SELEN_TABOO_RCMB_config

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write the CONFIG namelist to a temporary file
    OPEN(  UNIT = namelist_unit, FILE = TRIM( namelist_filename))
    WRITE( UNIT = namelist_unit, NML  = CONFIG)
    CLOSE( UNIT = namelist_unit)

    ! Check the config file for validity
    CALL check_config_file_validity( config_filename, namelist_filename)

    ! Delete the temporary namelist file
    CALL system('rm -f ' // TRIM( namelist_filename))

    ! Open the config file
    OPEN(  UNIT = config_unit, FILE = TRIM( config_filename), STATUS = 'OLD', ACTION = 'READ', IOSTAT = ios)
    IF (ios /= 0) CALL crash('couldnt open config file "' // TRIM( config_filename) // '"!')

    ! Read the config file using the CONFIG namelist
    READ(  UNIT = config_unit, NML = CONFIG, IOSTAT = ios)
    IF (ios /= 0) CALL crash('error while reading config file "' // TRIM( config_filename) // '"!')

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_config_file

  SUBROUTINE copy_config_variables_to_struct
    ! Overwrite the values in the fields of the config structure with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values that were read from the config file.

    IMPLICIT NONE

    ! In/output variables:

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'copy_config_variables_to_struct'

    ! Add routine to path
    CALL init_routine( routine_name)

  ! General model instructions
  ! ==========================

    ! Output directory
    C%create_procedural_output_dir             = create_procedural_output_dir_config
    C%fixed_output_dir                         = fixed_output_dir_config

    ! Debugging
    C%do_write_debug_data                      = do_write_debug_data_config
    C%do_check_for_NaN                         = do_check_for_NaN_config
    C%do_time_display                          = do_time_display_config

  ! == Time steps and range
  ! =======================

    C%start_time_of_run                        = start_time_of_run_config
    C%end_time_of_run                          = end_time_of_run_config
    C%dt_coupling                              = dt_coupling_config

  ! == Which model regions to simulate
  ! ==================================

    C%do_NAM                                   = do_NAM_config
    C%do_EAS                                   = do_EAS_config
    C%do_GRL                                   = do_GRL_config
    C%do_ANT                                   = do_ANT_config

  ! == Switch to do unit tests
  ! ==========================
    
    C%do_unit_tests                            = do_unit_tests_config

  ! == The four model regions
  ! =========================

    ! North America
    C%lambda_M_NAM                             = lambda_M_NAM_config
    C%phi_M_NAM                                = phi_M_NAM_config
    C%beta_stereo_NAM                          = beta_stereo_NAM_config
    C%xmin_NAM                                 = xmin_NAM_config
    C%xmax_NAM                                 = xmax_NAM_config
    C%ymin_NAM                                 = ymin_NAM_config
    C%ymax_NAM                                 = ymax_NAM_config

    ! Eurasia
    C%lambda_M_EAS                             = lambda_M_EAS_config
    C%phi_M_EAS                                = phi_M_EAS_config
    C%beta_stereo_EAS                          = beta_stereo_EAS_config
    C%xmin_EAS                                 = xmin_EAS_config
    C%xmax_EAS                                 = xmax_EAS_config
    C%ymin_EAS                                 = ymin_EAS_config
    C%ymax_EAS                                 = ymax_EAS_config

    ! Greenland
    C%lambda_M_GRL                             = lambda_M_GRL_config
    C%phi_M_GRL                                = phi_M_GRL_config
    C%beta_stereo_GRL                          = beta_stereo_GRL_config
    C%xmin_GRL                                 = xmin_GRL_config
    C%xmax_GRL                                 = xmax_GRL_config
    C%ymin_GRL                                 = ymin_GRL_config
    C%ymax_GRL                                 = ymax_GRL_config

    ! Antarctica
    C%lambda_M_ANT                             = lambda_M_ANT_config
    C%phi_M_ANT                                = phi_M_ANT_config
    C%beta_stereo_ANT                          = beta_stereo_ANT_config
    C%xmin_ANT                                 = xmin_ANT_config
    C%xmax_ANT                                 = xmax_ANT_config
    C%ymin_ANT                                 = ymin_ANT_config
    C%ymax_ANT                                 = ymax_ANT_config

  ! == Reference geometries (initial, present-day, and GIA equilibrium)
  ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    C%refgeo_Hi_min                            = refgeo_Hi_min_config
    C%remove_Lake_Vostok                       = remove_Lake_Vostok_config

    ! == Initial geometry
    ! ===================

    C%choice_refgeo_init_NAM                   = choice_refgeo_init_NAM_config
    C%choice_refgeo_init_EAS                   = choice_refgeo_init_EAS_config
    C%choice_refgeo_init_GRL                   = choice_refgeo_init_GRL_config
    C%choice_refgeo_init_ANT                   = choice_refgeo_init_ANT_config
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    C%choice_refgeo_init_idealised             = choice_refgeo_init_idealised_config
    C%dx_refgeo_init_idealised                 = dx_refgeo_init_idealised_config
    ! Path to file containing present-day geometry when choice_refgeo_init == 'read_from_file'
    C%filename_refgeo_init_NAM                 = filename_refgeo_init_NAM_config
    C%filename_refgeo_init_EAS                 = filename_refgeo_init_EAS_config
    C%filename_refgeo_init_GRL                 = filename_refgeo_init_GRL_config
    C%filename_refgeo_init_ANT                 = filename_refgeo_init_ANT_config
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_refgeo_init_NAM                = timeframe_refgeo_init_NAM_config
    C%timeframe_refgeo_init_EAS                = timeframe_refgeo_init_EAS_config
    C%timeframe_refgeo_init_GRL                = timeframe_refgeo_init_GRL_config
    C%timeframe_refgeo_init_ANT                = timeframe_refgeo_init_ANT_config

    ! == Present-day geometry
    ! =======================

    C%choice_refgeo_PD_NAM                     = choice_refgeo_PD_NAM_config
    C%choice_refgeo_PD_EAS                     = choice_refgeo_PD_EAS_config
    C%choice_refgeo_PD_GRL                     = choice_refgeo_PD_GRL_config
    C%choice_refgeo_PD_ANT                     = choice_refgeo_PD_ANT_config
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    C%choice_refgeo_PD_idealised               = choice_refgeo_PD_idealised_config
    C%dx_refgeo_PD_idealised                   = dx_refgeo_PD_idealised_config
    ! Path to file containing present-day geometry when choice_refgeo_PD == 'read_from_file'
    C%filename_refgeo_PD_NAM                   = filename_refgeo_PD_NAM_config
    C%filename_refgeo_PD_EAS                   = filename_refgeo_PD_EAS_config
    C%filename_refgeo_PD_GRL                   = filename_refgeo_PD_GRL_config
    C%filename_refgeo_PD_ANT                   = filename_refgeo_PD_ANT_config
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_refgeo_PD_NAM                  = timeframe_refgeo_PD_NAM_config
    C%timeframe_refgeo_PD_EAS                  = timeframe_refgeo_PD_EAS_config
    C%timeframe_refgeo_PD_GRL                  = timeframe_refgeo_PD_GRL_config
    C%timeframe_refgeo_PD_ANT                  = timeframe_refgeo_PD_ANT_config

    ! == GIA equilibrium geometry
    ! ===========================

    C%choice_refgeo_GIAeq_NAM                  = choice_refgeo_GIAeq_NAM_config
    C%choice_refgeo_GIAeq_EAS                  = choice_refgeo_GIAeq_EAS_config
    C%choice_refgeo_GIAeq_GRL                  = choice_refgeo_GIAeq_GRL_config
    C%choice_refgeo_GIAeq_ANT                  = choice_refgeo_GIAeq_ANT_config
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    C%choice_refgeo_GIAeq_idealised            = choice_refgeo_GIAeq_idealised_config
    C%dx_refgeo_GIAeq_idealised                = dx_refgeo_GIAeq_idealised_config
    ! Path to file containing present-day geometry when choice_refgeo_GIAeq == 'read_from_file'
    C%filename_refgeo_GIAeq_NAM                = filename_refgeo_GIAeq_NAM_config
    C%filename_refgeo_GIAeq_EAS                = filename_refgeo_GIAeq_EAS_config
    C%filename_refgeo_GIAeq_GRL                = filename_refgeo_GIAeq_GRL_config
    C%filename_refgeo_GIAeq_ANT                = filename_refgeo_GIAeq_ANT_config
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_refgeo_GIAeq_NAM               = timeframe_refgeo_GIAeq_NAM_config
    C%timeframe_refgeo_GIAeq_EAS               = timeframe_refgeo_GIAeq_EAS_config
    C%timeframe_refgeo_GIAeq_GRL               = timeframe_refgeo_GIAeq_GRL_config
    C%timeframe_refgeo_GIAeq_ANT               = timeframe_refgeo_GIAeq_ANT_config

    ! == Parameters for idealised geometries
    ! ======================================

    C%refgeo_idealised_slabonaslope_Hi         = refgeo_idealised_slabonaslope_Hi_config
    C%refgeo_idealised_slabonaslope_dhdx       = refgeo_idealised_slabonaslope_dhdx_config
    C%refgeo_idealised_Halfar_H0               = refgeo_idealised_Halfar_H0_config
    C%refgeo_idealised_Halfar_R0               = refgeo_idealised_Halfar_R0_config
    C%refgeo_idealised_Bueler_H0               = refgeo_idealised_Bueler_H0_config
    C%refgeo_idealised_Bueler_R0               = refgeo_idealised_Bueler_R0_config
    C%refgeo_idealised_Bueler_lambda           = refgeo_idealised_Bueler_lambda_config
    C%refgeo_idealised_SSA_icestream_Hi        = refgeo_idealised_SSA_icestream_Hi_config
    C%refgeo_idealised_SSA_icestream_dhdx      = refgeo_idealised_SSA_icestream_dhdx_config
    C%refgeo_idealised_MISMIP_mod_Hi_init      = refgeo_idealised_MISMIP_mod_Hi_init_config
    C%refgeo_idealised_ISMIP_HOM_L             = refgeo_idealised_ISMIP_HOM_L_config
    C%refgeo_idealised_MISMIPplus_Hi_init      = refgeo_idealised_MISMIPplus_Hi_init_config

  ! == Mesh generation
  ! ==================

    ! How to set up the initial mesh
    C%choice_initial_mesh_NAM                  = choice_initial_mesh_NAM_config
    C%choice_initial_mesh_EAS                  = choice_initial_mesh_EAS_config
    C%choice_initial_mesh_GRL                  = choice_initial_mesh_GRL_config
    C%choice_initial_mesh_ANT                  = choice_initial_mesh_ANT_config

    ! Paths to files containing initial meshes, if choice_initial_mesh == 'read_from_file'
    C%filename_initial_mesh_NAM                = filename_initial_mesh_NAM_config
    C%filename_initial_mesh_EAS                = filename_initial_mesh_EAS_config
    C%filename_initial_mesh_GRL                = filename_initial_mesh_GRL_config
    C%filename_initial_mesh_ANT                = filename_initial_mesh_ANT_config

    ! Resolutions for different parts of the ice sheet
    C%maximum_resolution_uniform               = maximum_resolution_uniform_config
    C%maximum_resolution_grounded_ice          = maximum_resolution_grounded_ice_config
    C%maximum_resolution_floating_ice          = maximum_resolution_floating_ice_config
    C%maximum_resolution_grounding_line        = maximum_resolution_grounding_line_config
    C%grounding_line_width                     = grounding_line_width_config
    C%maximum_resolution_calving_front         = maximum_resolution_calving_front_config
    C%calving_front_width                      = calving_front_width_config
    C%maximum_resolution_ice_front             = maximum_resolution_ice_front_config
    C%ice_front_width                          = ice_front_width_config
    C%maximum_resolution_coastline             = maximum_resolution_coastline_config
    C%coastline_width                          = coastline_width_config

    ! Regions of interest
    C%choice_regions_of_interest               = choice_regions_of_interest_config
    C%ROI_maximum_resolution_uniform           = ROI_maximum_resolution_uniform_config
    C%ROI_maximum_resolution_grounded_ice      = ROI_maximum_resolution_grounded_ice_config
    C%ROI_maximum_resolution_floating_ice      = ROI_maximum_resolution_floating_ice_config
    C%ROI_maximum_resolution_grounding_line    = ROI_maximum_resolution_grounding_line_config
    C%ROI_grounding_line_width                 = ROI_grounding_line_width_config
    C%ROI_maximum_resolution_calving_front     = ROI_maximum_resolution_calving_front_config
    C%ROI_calving_front_width                  = ROI_calving_front_width_config
    C%ROI_maximum_resolution_ice_front         = ROI_maximum_resolution_ice_front_config
    C%ROI_ice_front_width                      = ROI_ice_front_width_config
    C%ROI_maximum_resolution_coastline         = ROI_maximum_resolution_coastline_config
    C%ROI_coastline_width                      = ROI_coastline_width_config

    ! Advanced geometry parameters
    C%do_singlecore_mesh_creation              = do_singlecore_mesh_creation_config
    C%alpha_min                                = alpha_min_config
    C%nit_Lloyds_algorithm                     = nit_Lloyds_algorithm_config

    ! Memory
    C%nC_mem                                   = nC_mem_config

  ! == The scaled vertical coordinate zeta
  ! ======================================

    C%choice_zeta_grid                         = choice_zeta_grid_config
    C%nz                                       = nz_config
    C%zeta_irregular_log_R                     = zeta_irregular_log_R_config

  ! == Ice dynamics - velocity
  ! ==========================

    ! General
    C%choice_stress_balance_approximation      = choice_stress_balance_approximation_config
    C%n_flow                                   = n_flow_config
    C%m_enh_sheet                              = m_enh_sheet_config
    C%m_enh_shelf                              = m_enh_shelf_config
    C%choice_hybrid_SIASSA_scheme              = choice_hybrid_SIASSA_scheme_config
    C%do_GL_subgrid_friction                   = do_GL_subgrid_friction_config
    C%subgrid_friction_exponent                = subgrid_friction_exponent_config

    ! Some parameters for numerically solving the stress balance
    C%SIA_maximum_diffusivity                  = SIA_maximum_diffusivity_config
    C%visc_it_norm_dUV_tol                     = visc_it_norm_dUV_tol_config
    C%visc_it_nit                              = visc_it_nit_config
    C%visc_it_relax                            = visc_it_relax_config
    C%epsilon_sq_0                             = epsilon_sq_0_config
    C%visc_eff_min                             = visc_eff_min_config
    C%beta_max                                 = beta_max_config
    C%vel_max                                  = vel_max_config
    C%stress_balance_PETSc_rtol                = stress_balance_PETSc_rtol_config
    C%stress_balance_PETSc_abstol              = stress_balance_PETSc_abstol_config

  ! == Ice dynamics - sliding law
  ! =============================

    ! Sliding laws
    C%choice_sliding_law                       = choice_sliding_law_config
    C%choice_idealised_sliding_law             = choice_idealised_sliding_law_config
    C%slid_delta_v                             = slid_delta_v_config
    C%slid_Weertman_m                          = slid_Weertman_m_config
    C%slid_Budd_q_plastic                      = slid_Budd_q_plastic_config
    C%slid_Budd_u_threshold                    = slid_Budd_u_threshold_config
    C%slid_ZI_ut                               = slid_ZI_ut_config
    C%slid_ZI_p                                = slid_ZI_p_config

  ! == Ice dynamics - boundary conditions
  ! =====================================

    C%BC_u_west                                = BC_u_west_config
    C%BC_u_east                                = BC_u_east_config
    C%BC_u_south                               = BC_u_south_config
    C%BC_u_north                               = BC_u_north_config
    C%BC_v_west                                = BC_v_west_config
    C%BC_v_east                                = BC_v_east_config
    C%BC_v_south                               = BC_v_south_config
    C%BC_v_north                               = BC_v_north_config
    C%BC_H_west                                = BC_H_west_config
    C%BC_H_east                                = BC_H_east_config
    C%BC_H_south                               = BC_H_south_config
    C%BC_H_north                               = BC_H_north_config

  ! == Ice dynamics - time integration
  ! ==================================

    C%choice_timestepping                      = choice_timestepping_config
    C%choice_ice_integration_method            = choice_ice_integration_method_config
    C%dHi_PETSc_rtol                           = dHi_PETSc_rtol_config
    C%dHi_PETSc_abstol                         = dHi_PETSc_rtol_config

    ! Predictor-corrector ice-thickness update
    C%pc_epsilon                               = pc_epsilon_config
    C%pc_k_I                                   = pc_k_I_config
    C%pc_k_p                                   = pc_k_p_config
    C%pc_eta_min                               = pc_eta_min_config

  ! == Ice dynamics - calving
  ! =========================

    C%choice_calving_law                       = choice_calving_law_config
    C%calving_threshold_thickness_shelf        = calving_threshold_thickness_shelf_config
    C%calving_threshold_thickness_sheet        = calving_threshold_thickness_sheet_config
    C%max_calving_rounds                       = max_calving_rounds_config
    C%do_remove_shelves                        = do_remove_shelves_config
    C%remove_shelves_larger_than_PD            = remove_shelves_larger_than_PD_config
    C%continental_shelf_calving                = continental_shelf_calving_config
    C%continental_shelf_min_height             = continental_shelf_min_height_config

  ! == Ice dynamics - stabilisation
  ! ===============================

    C%choice_mask_noice                        = choice_mask_noice_config

    ! Partially fixed geometry, useful for initialisation and inversion runs
    C%fixed_shelf_geometry                     = fixed_shelf_geometry_config
    C%fixed_sheet_geometry                     = fixed_sheet_geometry_config
    C%fixed_grounding_line                     = fixed_grounding_line_config

  ! == Basal hydrology
  ! ==================

    ! Basal hydrology
    C%choice_basal_hydrology                   = choice_basal_hydrology_config
    C%Martin2011_hydro_Hb_min                  = Martin2011_hydro_Hb_min_config
    C%Martin2011_hydro_Hb_max                  = Martin2011_hydro_Hb_max_config

  ! == Bed roughness
  ! ==================

    C%choice_bed_roughness                     = choice_bed_roughness_config
    C%choice_bed_roughness_parameterised       = choice_bed_roughness_parameterised_config
    ! Paths to files containing bed roughness fields for the chosen sliding law
    C%filename_bed_roughness_NAM               = filename_bed_roughness_NAM_config
    C%filename_bed_roughness_EAS               = filename_bed_roughness_EAS_config
    C%filename_bed_roughness_GRL               = filename_bed_roughness_GRL_config
    C%filename_bed_roughness_ANT               = filename_bed_roughness_ANT_config
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_bed_roughness_NAM              = timeframe_bed_roughness_NAM_config
    C%timeframe_bed_roughness_EAS              = timeframe_bed_roughness_EAS_config
    C%timeframe_bed_roughness_GRL              = timeframe_bed_roughness_GRL_config
    C%timeframe_bed_roughness_ANT              = timeframe_bed_roughness_ANT_config
    ! Values for uniform bed roughness
    C%slid_Weertman_beta_sq_uniform            = slid_Weertman_beta_sq_uniform_config
    C%slid_Coulomb_phi_fric_uniform            = slid_Coulomb_phi_fric_uniform_config
    C%slid_Tsai2015_alpha_sq_uniform           = slid_Tsai2015_alpha_sq_uniform_config
    C%slid_Tsai2015_beta_sq_uniform            = slid_Tsai2015_beta_sq_uniform_config
    C%slid_Schoof2005_alpha_sq_uniform         = slid_Schoof2005_alpha_sq_uniform_config
    C%slid_Schoof2005_beta_sq_uniform          = slid_Schoof2005_beta_sq_uniform_config
    C%slid_ZI_phi_fric_uniform                 = slid_ZI_phi_fric_uniform_config
    ! Parameters for bed roughness parameterisations
    C%Martin2011till_phi_Hb_min                = Martin2011till_phi_Hb_min_config
    C%Martin2011till_phi_Hb_max                = Martin2011till_phi_Hb_max_config
    C%Martin2011till_phi_min                   = Martin2011till_phi_min_config
    C%Martin2011till_phi_max                   = Martin2011till_phi_max_config

  ! == Bed roughness inversion by nudging
  ! =====================================

    ! General
    C%do_bed_roughness_nudging                 = do_bed_roughness_nudging_config
    C%bed_roughness_nudging_t_start            = bed_roughness_nudging_t_start_config
    C%bed_roughness_nudging_t_end              = bed_roughness_nudging_t_end_config
    C%choice_bed_roughness_nudging_method      = choice_bed_roughness_nudging_method_config
    C%bed_roughness_nudging_dt                 = bed_roughness_nudging_dt_config
    C%filename_inverted_bed_roughness          = filename_inverted_bed_roughness_config

    ! Parameters
    C%BIVgeo_Bernales_do_smooth                = BIVgeo_Bernales_do_smooth_config
    C%BIVgeo_Bernales_scale                    = BIVgeo_Bernales_scale_config
    C%BIVgeo_Bernales_rsmooth                  = BIVgeo_Bernales_rsmooth_config
    C%BIVgeo_Bernales_wsmooth                  = BIVgeo_Bernales_wsmooth_config
    C%BIVgeo_Bernales_phi_min                  = BIVgeo_Bernales_phi_min_config
    C%BIVgeo_Bernales_phi_max                  = BIVgeo_Bernales_phi_max_config
    C%BIVgeo_Bernales_tol_diff                 = BIVgeo_Bernales_tol_diff_config
    C%BIVgeo_Bernales_tol_frac                 = BIVgeo_Bernales_tol_frac_config
    C%BIVgeo_Berends2022_tauc                  = BIVgeo_Berends2022_tauc_config
    C%BIVgeo_Berends2022_H0                    = BIVgeo_Berends2022_H0_config
    C%BIVgeo_Berends2022_u0                    = BIVgeo_Berends2022_u0_config
    C%BIVgeo_Berends2022_Hi_scale              = BIVgeo_Berends2022_Hi_scale_config
    C%BIVgeo_Berends2022_u_scale               = BIVgeo_Berends2022_u_scale_config
    C%BIVgeo_Berends2022_phimin                = BIVgeo_Berends2022_phimin_config
    C%BIVgeo_Berends2022_phimax                = BIVgeo_Berends2022_phimax_config
    C%BIVgeo_target_velocity_filename          = BIVgeo_target_velocity_filename_config

  ! == Thermodynamics and rheology
  ! ==============================

    ! Initial temperature profile
    C%choice_initial_ice_temperature_NAM       = choice_initial_ice_temperature_NAM_config
    C%choice_initial_ice_temperature_EAS       = choice_initial_ice_temperature_EAS_config
    C%choice_initial_ice_temperature_GRL       = choice_initial_ice_temperature_GRL_config
    C%choice_initial_ice_temperature_ANT       = choice_initial_ice_temperature_ANT_config
    ! Uniform initial ice temperature, if choice_initial_ice_temperature == 'uniform'
    C%uniform_initial_ice_temperature_NAM      = uniform_initial_ice_temperature_NAM_config
    C%uniform_initial_ice_temperature_EAS      = uniform_initial_ice_temperature_EAS_config
    C%uniform_initial_ice_temperature_GRL      = uniform_initial_ice_temperature_GRL_config
    C%uniform_initial_ice_temperature_ANT      = uniform_initial_ice_temperature_ANT_config
    ! Paths to files containing initial temperature fields, if
    C%filename_initial_ice_temperature_NAM     = filename_initial_ice_temperature_NAM_config
    C%filename_initial_ice_temperature_EAS     = filename_initial_ice_temperature_EAS_config
    C%filename_initial_ice_temperature_GRL     = filename_initial_ice_temperature_GRL_config
    C%filename_initial_ice_temperature_ANT     = filename_initial_ice_temperature_ANT_config
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_initial_ice_temperature_NAM    = timeframe_initial_ice_temperature_NAM_config
    C%timeframe_initial_ice_temperature_EAS    = timeframe_initial_ice_temperature_EAS_config
    C%timeframe_initial_ice_temperature_GRL    = timeframe_initial_ice_temperature_GRL_config
    C%timeframe_initial_ice_temperature_ANT    = timeframe_initial_ice_temperature_ANT_config
    ! Thermodynamical model

    C%choice_thermo_model                      = choice_thermo_model_config
    C%choice_ice_heat_capacity                 = choice_ice_heat_capacity_config
    C%uniform_ice_heat_capacity                = uniform_ice_heat_capacity_config
    C%choice_ice_thermal_conductivity          = choice_ice_thermal_conductivity_config
    C%uniform_ice_thermal_conductivity         = uniform_ice_thermal_conductivity_config

    ! Rheological model (relating Glen's flow parameter to ice temperature)
    C%choice_ice_rheology                      = choice_ice_rheology_config
    C%uniform_flow_factor                      = uniform_flow_factor_config

  ! == SELEN
  ! ========

    C%SELEN_run_at_t_start                     = SELEN_run_at_t_start_config
    C%SELEN_n_TDOF_iterations                  = SELEN_n_TDOF_iterations_config
    C%SELEN_n_recursion_iterations             = SELEN_n_recursion_iterations_config
    C%SELEN_use_rotational_feedback            = SELEN_use_rotational_feedback_config
    C%SELEN_n_harmonics                        = SELEN_n_harmonics_config
    C%SELEN_display_progress                   = SELEN_display_progress_config

    C%SELEN_dir                                = SELEN_dir_config
    C%SELEN_global_topo_filename               = SELEN_global_topo_filename_config
    C%SELEN_TABOO_init_filename                = SELEN_TABOO_init_filename_config
    C%SELEN_LMJ_VALUES_filename                = SELEN_LMJ_VALUES_filename_config

    C%SELEN_irreg_time_n                       = SELEN_irreg_time_n_config
    C%SELEN_irreg_time_window                  = SELEN_irreg_time_window_config

    C%SELEN_lith_thickness                     = SELEN_lith_thickness_config
    C%SELEN_visc_n                             = SELEN_visc_n_config
    C%SELEN_visc_prof                          = SELEN_visc_prof_config

    ! Settings for the TABOO Earth deformation model
    C%SELEN_TABOO_CDE                          = SELEN_TABOO_CDE_config
    C%SELEN_TABOO_TLOVE                        = SELEN_TABOO_TLOVE_config
    C%SELEN_TABOO_DEG1                         = SELEN_TABOO_DEG1_config
    C%SELEN_TABOO_RCMB                         = SELEN_TABOO_RCMB_config

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE copy_config_variables_to_struct

  SUBROUTINE check_config_file_validity( config_filename, namelist_filename)
    ! Check if the config file "config_filename" is valid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename
    CHARACTER(LEN=*),                INTENT(IN)        :: namelist_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_config_file_validity'
    LOGICAL                                            :: ex, all_are_valid, all_are_present

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    INQUIRE( FILE = config_filename, EXIST = ex)
    IF (.NOT. ex) CALL crash('config file ' // TRIM( config_filename) // ' could not be found!')
    INQUIRE( FILE = namelist_filename, EXIST = ex)
    IF (.NOT. ex) CALL crash('namelist file ' // TRIM( namelist_filename) // ' could not be found!')

    ! Check if all the variables appearing in the config file "config_filename" are valid
    CALL check_if_all_config_variables_are_valid( config_filename, namelist_filename, all_are_valid)

    ! Check if all the expected config variables appear in the config file "config_filename"
    CALL check_if_all_expected_config_variables_are_present( config_filename, namelist_filename, all_are_present)

    ! If not all is well, crash
    IF (.NOT. (all_are_valid .AND. all_are_present)) CALL crash('config file "' // TRIM( config_filename) // '" is invalid!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_config_file_validity

  SUBROUTINE check_if_all_config_variables_are_valid( config_filename, namelist_filename, all_are_valid)
    ! Check if all the variables appearing in the config file "config_filename" are valid
    !
    ! Do this by reading one line at a time of the config file, determining the name of the variable
    ! declared in that line, and checking if that variable also exists in the namelist file.
    !
    ! The namelist file is created earlier by writing the namelist to a text file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename
    CHARACTER(LEN=*),                INTENT(IN)        :: namelist_filename
    LOGICAL,                         INTENT(OUT)       :: all_are_valid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_if_all_config_variables_are_valid'
    INTEGER, PARAMETER                                 :: config_unit   = 1337
    INTEGER, PARAMETER                                 :: namelist_unit = 1338
    INTEGER                                            :: ios
    LOGICAL                                            :: found_end_of_file_config, found_end_of_file_namelist
    CHARACTER(256)                                     :: single_line_config      , single_line_namelist
    INTEGER                                            :: line_counter_config     , line_counter_namelist
    LOGICAL                                            :: found_match

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the config file
    OPEN( UNIT = config_unit, FILE = config_filename, IOSTAT = ios)
    IF (ios /= 0) CALL crash('couldnt open config file "' // TRIM( config_filename) // '"!')

    ! Read one line at a time of the config file, determine the name of the variable
    ! declared in that line, and check if that variable also exists in the namelist file

    found_end_of_file_config = .FALSE.
    line_counter_config      = 0
    all_are_valid            = .TRUE.

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
      IF (.NOT. found_match) THEN
        all_are_valid = .FALSE.
        CALL warning('invalid config variable "' // TRIM( single_line_config) // '" in file "' // TRIM( config_filename) // '", line {int_01}')
      END IF

      ! Close the namelist file
      CLOSE( UNIT = namelist_unit)

    END DO ! DO WHILE (.NOT. found_end_of_file_config)

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_if_all_config_variables_are_valid

  SUBROUTINE check_if_all_expected_config_variables_are_present( config_filename, namelist_filename, all_are_present)
    ! Check if all the expected config variables appear in the config file "config_filename"
    !
    ! Do this by reading one line at a time of the namelist file, determining the name of the variable
    ! declared in that line, and checking if that variable also exists in the config file.
    !
    ! The namelist file is created earlier by writing the namelist to a text file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                INTENT(IN)        :: config_filename
    CHARACTER(LEN=*),                INTENT(IN)        :: namelist_filename
    LOGICAL,                         INTENT(OUT)       :: all_are_present

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_if_all_expected_config_variables_are_present'
    INTEGER, PARAMETER                                 :: config_unit   = 1337
    INTEGER, PARAMETER                                 :: namelist_unit = 1338
    INTEGER                                            :: ios
    LOGICAL                                            :: found_end_of_file_config, found_end_of_file_namelist
    CHARACTER(256)                                     :: single_line_config      , single_line_namelist
    INTEGER                                            :: line_counter_config     , line_counter_namelist
    LOGICAL                                            :: found_match

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the namelist file
    OPEN( UNIT = namelist_unit, FILE = namelist_filename, IOSTAT = ios)
    IF (ios /= 0) CALL crash('couldnt open namelist file "' // TRIM( namelist_filename) // '"!')

    ! Read one line at a time of the nanmelist file, determine the name of the variable
    ! declared in that line, and check if that variable also exists in the config file

    found_end_of_file_namelist = .FALSE.
    line_counter_namelist      = 0
    all_are_present            = .TRUE.

    DO WHILE (.NOT. found_end_of_file_namelist)

      line_counter_namelist = line_counter_namelist + 1

      ! Read a single line from the config file
      READ( UNIT = namelist_unit, FMT = '(A)', IOSTAT = ios) single_line_namelist

      ! If we've reached the end of the file before finding the terminating forward slash, this namelist file is not valid.
      IF (ios < 0) CALL crash('namelist file "' // TRIM( namelist_filename) // '" is not terminated with a forward slash!')

      ! Remove all leading spaces
      CALL remove_leading_spaces( single_line_namelist)

      ! The variable name is the part of the string left of the first (, =, or space.
      single_line_namelist = single_line_namelist( 1:SCAN( single_line_namelist, '( =')-1)

      ! Get namelist variable in all caps for case-insensitive comparison
      CALL capitalise_string( single_line_namelist)

      ! The forward slash at the end terminates the namelist file
      IF (single_line_namelist == '/') THEN
        found_end_of_file_namelist = .TRUE.
      END IF

      ! Disregard empty lines, commented lines, the header line, and the final line
      IF (single_line_namelist == '' .OR. single_line_namelist( 1:1) == '&' .OR. single_line_namelist( 1:1) == '!') THEN
        CYCLE
      END IF

      ! Open the config file
      OPEN( UNIT = config_unit, FILE = config_filename)
      IF (ios /= 0) CALL crash('couldnt open config file "' // TRIM( config_filename) // '"!')

      ! Read all variables from the config file and check if any of them match the current namelist variable

      found_end_of_file_config = .FALSE.
      line_counter_config      = 0
      found_match              = .FALSE.

      DO WHILE ((.NOT. found_end_of_file_config) .AND. (.NOT. found_match))

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

        ! Check if this namelist variable matches the config variable
        IF (single_line_config == single_line_namelist) THEN
          found_match = .TRUE.
        END IF

      END DO ! DO WHILE ((.NOT. found_end_of_file_config) .AND. (.NOT. found_match))

      ! If no matching variable was found in the config file, print an error
      IF (.NOT. found_match) THEN
        all_are_present = .FALSE.
        CALL warning('couldnt find config variable "' // TRIM( single_line_namelist) // '" in file "' // TRIM( config_filename) // '"')
      END IF

      ! Close the namelist file
      CLOSE( UNIT = config_unit)

    END DO ! DO WHILE (.NOT. found_end_of_file_namelist)

    ! Close the config file
    CLOSE( UNIT = namelist_unit)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_if_all_expected_config_variables_are_present

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

END MODULE model_configuration
