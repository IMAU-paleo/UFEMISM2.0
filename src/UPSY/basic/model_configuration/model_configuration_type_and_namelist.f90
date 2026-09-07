module model_configuration_type_and_namelist

  ! The different parameters that control the simulation.
  !
  ! Each config variable has two versions: one with the "_config" extension, which is
  ! an actual variable in this module only, and one without the extension, which is
  ! a field in the "C" type. The "_config" variables are used to create a namelist,
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

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use model_configuration_utilities, only: check_config_file_validity

  implicit none

  private

  public :: type_config, copy_config_variables_to_struct, read_config_file

  ! ===== Module variables =====
  ! ============================

  ! The module-local "XXX_config" variables, which will be collected into a namelist, and replaced
  ! by the values in the external config file. Remember the "_config" extension!

  ! General model instructions
  ! ==========================

    ! Output directory
    logical             :: create_procedural_output_dir_config          = .true.                           ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)
    character(len=1024) :: fixed_output_dir_config                      = ''                               ! If not, create a directory with this name instead (stops the program if this directory already exists)

    ! Debugging
    logical             :: do_unit_tests_config                         = .false.                          ! Whether or not to (only) perform the unit tests in the main_validation module
    logical             :: do_benchmarks_config                         = .false.                          ! Whether or not to (only) perform the benchmarks
    logical             :: do_check_for_NaN_config                      = .false.                          ! Whether or not fields should be checked for NaN values
    logical             :: do_time_display_config                       = .true.                           ! Print current model time to screen
    logical             :: do_write_matrix_operators_config             = .false.                          ! Whether or not to write the operator matrices to output
    logical             :: do_write_checksum_log_config                 = .false.                          ! Whether or not checksums should be calcuialted and written to the checksum log file

  ! == Time of simulation
  ! =====================

    real(dp)            :: start_time_of_run_config                     = 0._dp                            ! [yr] Start time of the simulation
    real(dp)            :: end_time_of_run_config                       = 0._dp                            ! [yr] End   time of the simulation

  ! == Which model regions to simulate
  ! ==================================

    real(dp)            :: dt_coupling_config                           = 0._dp                            ! [yr] Interval of coupling between the different model regions
    logical             :: do_NAM_config                                = .false.
    logical             :: do_EAS_config                                = .false.
    logical             :: do_GRL_config                                = .false.
    logical             :: do_ANT_config                                = .false.

  ! == The four model regions
  ! =========================

    ! North America
    real(dp)            :: lambda_M_NAM_config                          = 265._dp                          ! [degrees east]  [default: 265.0]   Longitude of the pole of the stereographic projection for the North America domain
    real(dp)            :: phi_M_NAM_config                             =  62._dp                          ! [degrees north] [default:  62.0]   Latitude  of the pole of the stereographic projection for the North America domain
    real(dp)            :: beta_stereo_NAM_config                       =  71._dp                          ! [degrees]       [default:  71.0]   Standard parallel     of the stereographic projection for the North America domain
    real(dp)            :: xmin_NAM_config                              = -3600000._dp                     ! [m]             [default: -3600E3] Western  boundary     of the North America domain
    real(dp)            :: xmax_NAM_config                              =  3600000._dp                     ! [m]             [default:  3600E3] Eastern  boundary     of the North America domain
    real(dp)            :: ymin_NAM_config                              = -2400000._dp                     ! [m]             [default: -2400E3] Southern boundary     of the North America domain
    real(dp)            :: ymax_NAM_config                              =  2400000._dp                     ! [m]             [default:  2400E3] Northern boundary     of the North America domain

    ! Eurasia
    real(dp)            :: lambda_M_EAS_config                          = 40._dp                           ! [degrees east]  [default:  40.0]   Longitude of the pole of the stereographic projection for the Eurasia domain
    real(dp)            :: phi_M_EAS_config                             = 70._dp                           ! [degrees north] [default:  70.0]   Latitude  of the pole of the stereographic projection for the Eurasia domain
    real(dp)            :: beta_stereo_EAS_config                       = 71._dp                           ! [degrees]       [default:  71.0]   Standard parallel     of the stereographic projection for the Eurasia domain
    real(dp)            :: xmin_EAS_config                              = -3400000._dp                     ! [m]             [default: -3400E3] Western  boundary     of the Eurasia domain
    real(dp)            :: xmax_EAS_config                              =  3400000._dp                     ! [m]             [default:  3400E3] Eastern  boundary     of the Eurasia domain
    real(dp)            :: ymin_EAS_config                              = -2080000._dp                     ! [m]             [default: -2080E3] Southern boundary     of the Eurasia domain
    real(dp)            :: ymax_EAS_config                              =  2080000._dp                     ! [m]             [default:  2080E3] Northern boundary     of the Eurasia domain

    ! Greenland
    real(dp)            :: lambda_M_GRL_config                          = -45._dp                          ! [degrees east]  [default: -45.0]   Longitude of the pole of the stereographic projection for the Greenland domain
    real(dp)            :: phi_M_GRL_config                             = 90._dp                           ! [degrees north] [default:  90.0]   Latitude  of the pole of the stereographic projection for the Greenland domain
    real(dp)            :: beta_stereo_GRL_config                       = 70._dp                           ! [degrees]       [default:  70.0]   Standard parallel     of the stereographic projection for the Greenland domain
    real(dp)            :: xmin_GRL_config                              =  -720000._dp                     ! [m]             [default:  -720E3] Western  boundary     of the Greenland domain [m]
    real(dp)            :: xmax_GRL_config                              =   960000._dp                     ! [m]             [default:   960E3] Eastern  boundary     of the Greenland domain [m]
    real(dp)            :: ymin_GRL_config                              = -3450000._dp                     ! [m]             [default: -3450E3] Southern boundary     of the Greenland domain [m]
    real(dp)            :: ymax_GRL_config                              =  -570000._dp                     ! [m]             [default:  -570E3] Northern boundary     of the Greenland domain [m]

    ! Antarctica
    real(dp)            :: lambda_M_ANT_config                          = 0._dp                            ! [degrees east]  [default:   0.0]   Longitude of the pole of the stereographic projection for the Antarctica domain
    real(dp)            :: phi_M_ANT_config                             = -90._dp                          ! [degrees north] [default: -90.0]   Latitude  of the pole of the stereographic projection for the Antarctica domain
    real(dp)            :: beta_stereo_ANT_config                       = 71._dp                           ! [degrees]       [default:  71.0]   Standard parallel     of the stereographic projection for the Antarctica domain
    real(dp)            :: xmin_ANT_config                              = -3040000._dp                     ! [m]             [default: -3040E3] Western  boundary     of the Antarctica domain [m]
    real(dp)            :: xmax_ANT_config                              =  3040000._dp                     ! [m]             [default:  3040E3] Eastern  boundary     of the Antarctica domain [m]
    real(dp)            :: ymin_ANT_config                              = -3040000._dp                     ! [m]             [default: -3040E3] Southern boundary     of the Antarctica domain [m]
    real(dp)            :: ymax_ANT_config                              =  3040000._dp                     ! [m]             [default:  3040E3] Northern boundary     of the Antarctica domain [m]

  ! == Reference geometries (initial, present-day, and GIA equilibrium)
  ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    real(dp)            :: refgeo_Hi_min_config                         = 2.0_dp                           ! [m]             [default: 2.0]     Remove ice thinner than this value in the reference ice geometry. Particularly useful for BedMachine Greenland, which somehow covers the entire tundra with half a meter of ice...
    logical             :: do_smooth_geometry_config                    = .false.                          ! Whether or not to smooth the reference bedrock
    real(dp)            :: r_smooth_geometry_config                     = 0.5_dp                           ! [m]             Geometry smoothing radius
    logical             :: remove_Lake_Vostok_config                    = .true.                           ! Whether or not to replace subglacial Lake Vostok in Antarctica with ice (recommended to set to TRUE, otherwise it will really slow down your model for the first few hundred years...)

    ! Remapping method for gridded geometry input
    character(len=1024) :: choice_refgeo_remapping_method_config        = '2nd_order_conservative'         ! Remapping method to apply to gridded geometry input. Default is 2nd_order_conservative, though for high resolution, 1st_order_conservative may be better.

    ! == Initial geometry
    ! ===================

    character(len=1024) :: choice_refgeo_init_NAM_config                = 'read_from_file'                 ! Choice of initial geometry for North America; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_init_EAS_config                = 'read_from_file'                 ! Choice of initial geometry for Eurasia      ; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_init_GRL_config                = 'read_from_file'                 ! Choice of initial geometry for Greenland    ; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_init_ANT_config                = 'read_from_file'                 ! Choice of initial geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    character(len=1024) :: choice_refgeo_init_idealised_config          = ''                               ! Choice of idealised initial geometry; see reference_geometries/calc_idealised_geometry for options
    real(dp)            :: dx_refgeo_init_idealised_config              = 5000._dp                         ! Resolution of square grid used for idealised initial geometry
    ! Path to file containing initial geometry when choice_refgeo_init == 'read_from_file'
    character(len=1024) :: filename_refgeo_init_NAM_config              = 'external/data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    character(len=1024) :: filename_refgeo_init_EAS_config              = 'external/data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    character(len=1024) :: filename_refgeo_init_GRL_config              = 'external/data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    character(len=1024) :: filename_refgeo_init_ANT_config              = 'external/data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_refgeo_init_NAM_config             = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    real(dp)            :: timeframe_refgeo_init_EAS_config             = 1E9_dp
    real(dp)            :: timeframe_refgeo_init_GRL_config             = 1E9_dp
    real(dp)            :: timeframe_refgeo_init_ANT_config             = 1E9_dp

    ! == Present-day geometry
    ! =======================

    character(len=1024) :: choice_refgeo_PD_NAM_config                  = 'read_from_file'                 ! Choice of present-day geometry for North America; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_PD_EAS_config                  = 'read_from_file'                 ! Choice of present-day geometry for Eurasia      ; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_PD_GRL_config                  = 'read_from_file'                 ! Choice of present-day geometry for Greenland    ; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_PD_ANT_config                  = 'read_from_file'                 ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    character(len=1024) :: choice_refgeo_PD_idealised_config            = ''                               ! Choice of idealised present-day geometry; see reference_geometries/calc_idealised_geometry for options
    real(dp)            :: dx_refgeo_PD_idealised_config                = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry
    ! Path to file containing present-day geometry when choice_refgeo_PD == 'read_from_file'
    character(len=1024) :: filename_refgeo_PD_NAM_config                = 'external/data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    character(len=1024) :: filename_refgeo_PD_EAS_config                = 'external/data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    character(len=1024) :: filename_refgeo_PD_GRL_config                = 'external/data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    character(len=1024) :: filename_refgeo_PD_ANT_config                = 'external/data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_refgeo_PD_NAM_config               = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    real(dp)            :: timeframe_refgeo_PD_EAS_config               = 1E9_dp
    real(dp)            :: timeframe_refgeo_PD_GRL_config               = 1E9_dp
    real(dp)            :: timeframe_refgeo_PD_ANT_config               = 1E9_dp

    ! == GIA equilibrium geometry
    ! ===========================

    character(len=1024) :: choice_refgeo_GIAeq_NAM_config               = 'read_from_file'                 ! Choice of GIA equilibrium reference geometry for North America; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_GIAeq_EAS_config               = 'read_from_file'                 ! Choice of GIA equilibrium reference geometry for Eurasia      ; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_GIAeq_GRL_config               = 'read_from_file'                 ! Choice of GIA equilibrium reference geometry for Greenland    ; can be "idealised", or "read_from_file"
    character(len=1024) :: choice_refgeo_GIAeq_ANT_config               = 'read_from_file'                 ! Choice of GIA equilibrium reference geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    character(len=1024) :: choice_refgeo_GIAeq_idealised_config         = ''                               ! Choice of idealised GIA equilibrium reference geometry; see reference_geometries/calc_idealised_geometry for options
    real(dp)            :: dx_refgeo_GIAeq_idealised_config             = 5000._dp                         ! Resolution of square grid used for idealised GIA equilibrium reference geometry
    ! Path to file containing GIA equilibrium reference geometry when choice_refgeo_GIAeq == 'read_from_file'
    character(len=1024) :: filename_refgeo_GIAeq_NAM_config             = 'external/data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    character(len=1024) :: filename_refgeo_GIAeq_EAS_config             = 'external/data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    character(len=1024) :: filename_refgeo_GIAeq_GRL_config             = 'external/data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    character(len=1024) :: filename_refgeo_GIAeq_ANT_config             = 'external/data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_refgeo_GIAeq_NAM_config            = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    real(dp)            :: timeframe_refgeo_GIAeq_EAS_config            = 1E9_dp
    real(dp)            :: timeframe_refgeo_GIAeq_GRL_config            = 1E9_dp
    real(dp)            :: timeframe_refgeo_GIAeq_ANT_config            = 1E9_dp

    ! == Parameters for idealised geometries
    ! ======================================

    real(dp)            :: refgeo_idealised_slabonaslope_Hi_config      = 0._dp                            ! Suggested value: 2000 m
    real(dp)            :: refgeo_idealised_slabonaslope_dhdx_config    = 1._dp                            ! Suggested value: -0.001
    real(dp)            :: refgeo_idealised_Halfar_H0_config            = 0._dp                            ! Suggested value: 3000 m
    real(dp)            :: refgeo_idealised_Halfar_R0_config            = 0._dp                            ! Suggested value: 500E3 m
    real(dp)            :: refgeo_idealised_Bueler_H0_config            = 0._dp                            ! Suggested value: 3000 m
    real(dp)            :: refgeo_idealised_Bueler_R0_config            = 0._dp                            ! Suggested value: 500E3 m
    real(dp)            :: refgeo_idealised_Bueler_lambda_config        = 0._dp                            ! Suggested value: 5.0
    real(dp)            :: refgeo_idealised_SSA_icestream_Hi_config     = -1._dp                           ! Suggested value: 2000 m
    real(dp)            :: refgeo_idealised_SSA_icestream_dhdx_config   = 1._dp                            ! Suggested value: -0.0003
    real(dp)            :: refgeo_idealised_SSA_icestream_L_config      = 0._dp                            ! Suggested value: 150 km
    real(dp)            :: refgeo_idealised_SSA_icestream_m_config      = 0._dp                            ! Suggested value: 1.0
    real(dp)            :: refgeo_idealised_MISMIP_mod_Hi_init_config   = -1._dp                           ! Suggested value: 100 m
    real(dp)            :: refgeo_idealised_ISMIP_HOM_L_config          = 0._dp                            ! Suggested value: 5E3 - 160E3 m
    real(dp)            :: refgeo_idealised_MISMIPplus_Hi_init_config   = -1._dp                           ! Suggested value: 100 m
    logical             :: refgeo_idealised_MISMIPplus_tune_A_config    = .false.                          ! If so, the uniform flow factor A is tuned to achieve a steady-state mid-stream grounding-line position at x = 450 km

  ! == Mesh generation
  ! ==================

    ! How to set up the initial mesh
    character(len=1024) :: choice_initial_mesh_NAM_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'
    character(len=1024) :: choice_initial_mesh_EAS_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'
    character(len=1024) :: choice_initial_mesh_GRL_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'
    character(len=1024) :: choice_initial_mesh_ANT_config               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'

    ! Paths to files containing initial meshes, if choice_initial_mesh == 'read_from_file'
    character(len=1024) :: filename_initial_mesh_NAM_config             = ''
    character(len=1024) :: filename_initial_mesh_EAS_config             = ''
    character(len=1024) :: filename_initial_mesh_GRL_config             = ''
    character(len=1024) :: filename_initial_mesh_ANT_config             = ''

    ! Resolutions for different parts of the ice sheet
    real(dp)            :: maximum_resolution_uniform_config            = 800e3_dp                         ! [m]          Maximum resolution for the entire domain
    real(dp)            :: maximum_resolution_grounded_ice_config       = 400e3_dp                         ! [m]          Maximum resolution for grounded ice
    real(dp)            :: maximum_resolution_floating_ice_config       = 200e3_dp                         ! [m]          Maximum resolution for floating ice
    real(dp)            :: maximum_resolution_grounding_line_config     = 100e3_dp                         ! [m]          Maximum resolution for the grounding line
    real(dp)            :: grounding_line_width_config                  = 200e3_dp                         ! [m]          Width of the band around the grounding line that should get this resolution
    real(dp)            :: maximum_resolution_calving_front_config      = 100e3_dp                         ! [m]          Maximum resolution for the calving front
    real(dp)            :: calving_front_width_config                   = 100e3_dp                         ! [m]          Width of the band around the calving front that should get this resolution
    real(dp)            :: maximum_resolution_ice_front_config          = 100e3_dp                         ! [m]          Maximum resolution for the ice front
    real(dp)            :: ice_front_width_config                       = 200e3_dp                         ! [m]          Width of the band around the ice front that should get this resolution
    real(dp)            :: maximum_resolution_coastline_config          = 100e3_dp                         ! [m]          Maximum resolution for the coastline
    real(dp)            :: coastline_width_config                       = 200e3_dp                         ! [m]          Width of the band around the coastline that should get this resolution

    ! Regions of interest
    character(len=1024) :: choice_regions_of_interest_config            = ''                               ! Regions of interest where other (higher) resolutions apply. Separated by double vertical bars "||", e.g. "PineIsland||Thwaites"
    real(dp)            :: ROI_maximum_resolution_uniform_config        = 100e3_dp                         ! [m]          Maximum resolution for the entire domain
    real(dp)            :: ROI_maximum_resolution_grounded_ice_config   = 50e3_dp                          ! [m]          Maximum resolution for grounded ice
    real(dp)            :: ROI_maximum_resolution_floating_ice_config   = 20e3_dp                          ! [m]          Maximum resolution for floating ice
    real(dp)            :: ROI_maximum_resolution_grounding_line_config = 5e3_dp                           ! [m]          Maximum resolution for the grounding line
    real(dp)            :: ROI_grounding_line_width_config              = 2e3_dp                           ! [m]          Width of the band around the grounding line that should get this resolution
    real(dp)            :: ROI_maximum_resolution_calving_front_config  = 10e3_dp                          ! [m]          Maximum resolution for the calving front
    real(dp)            :: ROI_calving_front_width_config               = 10e3_dp                          ! [m]          Width of the band around the calving front that should get this resolution
    real(dp)            :: ROI_maximum_resolution_ice_front_config      = 20e3_dp                          ! [m]          Maximum resolution for the ice front
    real(dp)            :: ROI_ice_front_width_config                   = 20e3_dp                          ! [m]          Width of the band around the ice front that should get this resolution
    real(dp)            :: ROI_maximum_resolution_coastline_config      = 50e3_dp                          ! [m]          Maximum resolution for the coastline
    real(dp)            :: ROI_coastline_width_config                   = 50e3_dp                          ! [m]          Width of the band around the coastline that should get this resolution

    ! Miscellaneous refinement options
    logical             :: do_refine_TransAntMounts_glaciers_config     = .false.                          !              Whether or not to refine the mesh over the Transantarctic Mountains glaciers (resolving those really helps in getting a stable Ross ice shelf)
    real(dp)            :: max_res_TransAntMounts_glaciers_config       = 5e3_dp                           ! [m]          Maximum resolution for the Transantarctic Mountains glaciers

    ! Mesh update settings
    logical             :: allow_mesh_updates_config                    = .true.                           ! [-]          Whether or not mesh updates are allowed
    real(dp)            :: dt_mesh_update_min_config                    = 50._dp                           ! [yr]         Minimum amount of time between mesh updates
    real(dp)            :: minimum_mesh_fitness_coefficient_config      = 0.95_dp                          ! [-]          If the mesh fitness coefficient drops below this threshold, trigger a mesh update
    logical             :: do_out_of_time_calving_front_relax_config    = .true.                           !              Whether or not to step out of time and relax the calving front for a few iterations after a mesh update

    ! Advanced geometry parameters
    logical             :: do_singlecore_mesh_creation_config           = .true.                           !              Whether or not to use only a single core for mesh generation (for better reproducibility)
    real(dp)            :: alpha_min_config                             = 0.4363_dp                        ! [radians]    Smallest allowed internal triangle angle (recommended value: 25 degrees = 0.4363)
    integer             :: nit_Lloyds_algorithm_config                  = 3                                ! [-]          Number of iterations of Lloyds algorithm to be applied after refinement
    real(dp)            :: mesh_resolution_tolerance_config             = 1.25_dp                          ! [-]          Factors the target resolution for trangle-size requirement. 1=strict, use >1 to avoid unnecesarily high resolution

    ! Square grid used for smoothing
    real(dp)            :: dx_square_grid_smooth_NAM_config             = 5000._dp
    real(dp)            :: dx_square_grid_smooth_EAS_config             = 5000._dp
    real(dp)            :: dx_square_grid_smooth_GRL_config             = 5000._dp
    real(dp)            :: dx_square_grid_smooth_ANT_config             = 5000._dp

  ! == The scaled vertical coordinate zeta
  ! ======================================

    character(len=1024) :: choice_zeta_grid_config                      = 'regular'                        ! The type of vertical grid to use; can be "regular", "irregular_log", "old_15_layer_zeta"
    integer             :: nz_config                                    = 12                               ! The number of vertical layers to use
    real(dp)            :: zeta_irregular_log_R_config                  = 10._dp                           ! Ratio between surface and base layer spacings

  ! == Ice dynamics - velocity
  ! ==========================

    ! General
    character(len=1024) :: choice_stress_balance_approximation_config   = 'DIVA'                           ! Choice of stress balance approximation: "none" (= no flow, though geometry can still change due to mass balance), "SIA", "SSA", "SIA/SSA", "DIVA", "BPA", "SSA_FEM_PETSc" (finite-element SSA solved with PetscFE/PetscSNES)
    character(len=1024) :: choice_hybrid_SIASSA_scheme_config           = 'add'                            ! Choice of scheme for combining SIA and SSA velocities in the hybrid approach
    logical             :: do_include_SSADIVA_crossterms_config         = .true.                           ! Whether or not to include the gradients of the effective viscosity (the "cross-terms") in the solution of the SSA/DIVA

    ! Hybrid DIVA/BPA
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_NAM_config       = ''                               ! How to determine where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for North America
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_EAS_config       = ''                               ! How to determine where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for Eurasia
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_GRL_config       = ''                               ! How to determine where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for Greenland
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_ANT_config       = ''                               ! How to determine where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for Antarctica

    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_NAM_config     = ''                               ! Path to a file containing a mask that describes where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for North America
    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_EAS_config     = ''                               ! Path to a file containing a mask that describes where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for Eurasia
    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_GRL_config     = ''                               ! Path to a file containing a mask that describes where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for Greenland
    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_ANT_config     = ''                               ! Path to a file containing a mask that describes where to solve the DIVA and where the BPA in the hybrid DIVA/BPA for Antarctica

    ! Initialisation
    character(len=1024) :: choice_initial_velocity_NAM_config           = 'zero'                           ! Can be 'zero', 'read_from_file'
    character(len=1024) :: choice_initial_velocity_EAS_config           = 'zero'
    character(len=1024) :: choice_initial_velocity_GRL_config           = 'zero'
    character(len=1024) :: choice_initial_velocity_ANT_config           = 'zero'
    ! Paths to files containing initial velocity fields
    character(len=1024) :: filename_initial_velocity_NAM_config         = ''
    character(len=1024) :: filename_initial_velocity_EAS_config         = ''
    character(len=1024) :: filename_initial_velocity_GRL_config         = ''
    character(len=1024) :: filename_initial_velocity_ANT_config         = ''
    ! Timeframes to read from the initial velocity files (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_initial_velocity_NAM_config        = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    real(dp)            :: timeframe_initial_velocity_EAS_config        = 1E9_dp
    real(dp)            :: timeframe_initial_velocity_GRL_config        = 1E9_dp
    real(dp)            :: timeframe_initial_velocity_ANT_config        = 1E9_dp

    ! Some parameters for numerically solving the stress balance
    real(dp)            :: SIA_maximum_diffusivity_config               = 1E5_dp                           ! Limit the diffusivity in the SIA to this value
    real(dp)            :: visc_it_norm_dUV_tol_config                  = 5E-5_dp                          ! Stop criterion for the viscosity iteration: the L2-norm of successive velocity solutions should be smaller than this number
    integer             :: visc_it_nit_config                           = 50                               ! Maximum number of effective viscosity iterations
    real(dp)            :: visc_it_relax_config                         = 0.2_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    real(dp)            :: visc_eff_min_config                          = 1E4_dp                           ! Minimum value for effective viscosity
    real(dp)            :: vel_max_config                               = 5000._dp                         ! Velocities are limited to this value
    character(len=1024) :: stress_balance_PETSc_KSPtype_config          = 'gmres'                          ! Which Krylov solver PETSc should use to solve the momentum balance
    character(len=1024) :: stress_balance_PETSc_PCtype_config           = 'bjacobi'                        ! Which preconditioner PETSc should use
    real(dp)            :: stress_balance_PETSc_rtol_config             = 1E-7_dp                          ! PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    real(dp)            :: stress_balance_PETSc_abstol_config           = 1E-5_dp                          ! PETSc solver - stop criterion, absolute difference

    ! Boundary conditions
    character(len=1024) :: BC_ice_front_config                          = 'infinite_slab'                  ! Boundary conditions to the momentum balance at the ice front: "infinite_slab", "ocean_pressure"
    character(len=1024) :: BC_u_west_config                             = 'infinite'                       ! Boundary conditions to the x-component of the momentum balance at the domain border: "infinite", "zero", "periodic_ISMIP-HOM"
    character(len=1024) :: BC_u_east_config                             = 'infinite'
    character(len=1024) :: BC_u_south_config                            = 'infinite'
    character(len=1024) :: BC_u_north_config                            = 'infinite'
    character(len=1024) :: BC_v_west_config                             = 'infinite'
    character(len=1024) :: BC_v_east_config                             = 'infinite'
    character(len=1024) :: BC_v_south_config                            = 'infinite'
    character(len=1024) :: BC_v_north_config                            = 'infinite'

  ! == Ice dynamics - sliding
  ! =========================

    ! General
    character(len=1024) :: choice_sliding_law_config                    = 'Zoet-Iverson'                   ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Budd", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"
    character(len=1024) :: choice_idealised_sliding_law_config          = ''                               ! "ISMIP_HOM_C", "ISMIP_HOM_D", "ISMIP_HOM_E", "ISMIP_HOM_F"

    ! Parameters for different sliding laws
    real(dp)            :: slid_Weertman_m_config                       = 3._dp                            ! Exponent in Weertman sliding law
    real(dp)            :: slid_Budd_q_plastic_config                   = 0.3_dp                           ! Scaling exponent in Budd sliding law
    real(dp)            :: slid_Budd_u_threshold_config                 = 100._dp                          ! Threshold velocity in Budd sliding law
    real(dp)            :: slid_ZI_p_config                             = 5._dp                            ! Velocity exponent used in the Zoet-Iverson sliding law
    real(dp)            :: slid_ZI_ut_config                            = 200._dp                          ! (uniform) transition velocity used in the Zoet-Iverson sliding law [m/yr]

    ! Sub-grid scaling of basal friction
    logical             :: do_GL_subgrid_friction_config                = .true.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
    character(len=1024) :: choice_subgrid_grounded_fraction_config      = 'bilin_interp_TAF+bedrock_CDF'   ! Choice of scheme to calculate the sub-grid grounded fractions: 'bilin_interp_TAF', 'bedrock_CDF', 'bilin_interp_TAF+bedrock_CDF'
    logical             :: do_read_bedrock_cdf_from_file_config         = .false.                          ! Whether or not to read the bedrock CDF from the initial mesh file. Requires choice_initial_mesh_XXX_config =  'read_from_file'!
    integer             :: subgrid_bedrock_cdf_nbins_config             = 11                               ! Number of bins to be used for sub-grid bedrock cumulative density functions
    logical             :: do_subgrid_friction_on_A_grid_config         = .false.                          ! Whether or not to apply a preliminary scaling to basal hydrology and bed roughness on the A grid, before the final scaling of beta on the b grid. The exponent of this scaling is computed based on ice thickness.
    real(dp)            :: subgrid_friction_exponent_on_B_grid_config   = 2._dp                            ! Exponent to which f_grnd should be raised before being used to scale the final value of beta on the B grid

    ! Stability
    real(dp)            :: slid_beta_max_config                         = 1E20_dp                          ! Maximum value for basal friction coefficient
    real(dp)            :: slid_delta_v_config                          = 1.0E-3_dp                        ! Normalisation parameter to prevent errors when velocity is zero

  ! == Ice dynamics - ice thickness calculation
  ! ===========================================

    ! Calculation of dH/dt
    character(len=1024) :: choice_ice_integration_method_config         = 'semi-implicit'                  ! Choice of ice thickness integration scheme: "none" (i.e. unchanging geometry), "explicit", "semi-implicit"
    real(dp)            :: dHi_semiimplicit_fs_config                   = 1.5_dp                           ! Factor for the semi-implicit ice thickness solver (0 = explicit, 0<f<1 = semi-implicit, 1 = implicit, >1 = over-implicit)
    character(len=1024) :: dHi_PETSc_KSPtype_config                     = 'gmres'                          ! Which Krylov solver PETSc should use to solve the mass continuity equation
    character(len=1024) :: dHi_PETSc_PCtype_config                      = 'bjacobi'                        ! Which preconditioner PETSc should use
    real(dp)            :: dHi_PETSc_rtol_config                        = 1E-8_dp                          ! dHi PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    real(dp)            :: dHi_PETSc_abstol_config                      = 1E-6_dp                          ! dHi PETSc solver - stop criterion, absolute difference

    ! Boundary conditions
    character(len=1024) :: BC_H_west_config                             = 'zero'                           ! Boundary conditions for ice thickness at the domain boundary
    character(len=1024) :: BC_H_east_config                             = 'zero'                           ! Allowed choices:  "infinite", "zero", "ISMIP_HOM_F"
    character(len=1024) :: BC_H_south_config                            = 'zero'
    character(len=1024) :: BC_H_north_config                            = 'zero'

  ! == Ice dynamics - target quantities
  ! ===================================

    ! Target dHi_dt
    logical             :: do_target_dHi_dt_config                      = .false.                          ! Whether or not to use a target dHi_dt field from an external file to alter the modelled one
    logical             :: do_limit_target_dHi_dt_to_SMB_config         = .false.                          ! Whether or not to limit a target dHi_dt field to available SMB in that area, so it does not melt all ice there
    real(dp)            :: target_dHi_dt_t_end_config                   = -9.9E9_dp                        ! End time of target dHi_dt application

    ! Files containing a target dHi_dt for inversions
    character(len=1024) :: filename_dHi_dt_target_NAM_config            = ''
    character(len=1024) :: filename_dHi_dt_target_EAS_config            = ''
    character(len=1024) :: filename_dHi_dt_target_GRL_config            = ''
    character(len=1024) :: filename_dHi_dt_target_ANT_config            = ''

    ! Timeframes for reading target dHi_dt from file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_dHi_dt_target_NAM_config           = 1E9_dp
    real(dp)            :: timeframe_dHi_dt_target_EAS_config           = 1E9_dp
    real(dp)            :: timeframe_dHi_dt_target_GRL_config           = 1E9_dp
    real(dp)            :: timeframe_dHi_dt_target_ANT_config           = 1E9_dp

  ! == Ice dynamics - time stepping
  ! ===============================

    ! Time stepping
    character(len=1024) :: choice_timestepping_config                   = 'pc'                             ! Choice of timestepping method: "direct", "pc" (NOTE: 'direct' does not work with DIVA ice dynamcis!)
    real(dp)            :: dt_ice_max_config                            = 10.0_dp                          ! [yr] Maximum time step of the ice dynamics model
    real(dp)            :: dt_ice_min_config                            = 0.1_dp                           ! [yr] Minimum time step of the ice dynamics model
    real(dp)            :: dt_ice_startup_phase_config                  = 0._dp                            ! [yr] Length of time window after start_time and before end_time when dt = dt_min, to ensure smooth restarts
    logical             :: do_grounded_only_adv_dt_config               = .false.                          ! If .true., only grounded ice is used in the computation of the advective time step limit (CFL condition)

    ! Predictor-corrector ice-thickness update
    real(dp)            :: pc_epsilon_config                            = 0.005_dp                         ! [m/yr] Target truncation error in dHi_dt (epsilon in Robinson et al., 2020, Eq. 33)
    real(dp)            :: pc_k_I_config                                = 0.2_dp                           ! Exponent k_I in  Robinson et al., 2020, Eq. 33
    real(dp)            :: pc_k_p_config                                = 0.2_dp                           ! Exponent k_p in  Robinson et al., 2020, Eq. 33
    real(dp)            :: pc_eta_min_config                            = 1E-8_dp                          ! Normalisation term in estimation of the truncation error (Robinson et al., Eq. 32)
    real(dp)            :: pc_max_time_step_increase_config             = 1.1_dp                           ! Each new time step is only allowed to be this much larger than the previous one
    integer             :: pc_nit_max_config                            = 50                               ! Maximum number of times a PC timestep can be repeated with reduced dt
    real(dp)            :: pc_guilty_max_config                         = 0.00_dp                          ! [%] Maximum percentage of grounded vertices allowed to have a truncation error larger than the tolerance pc_epsilon

    ! Initialisation of the predictor-corrector ice-thickness update
    character(len=1024) :: pc_choice_initialise_NAM_config              = 'zero'                           ! How to initialise the p/c scheme: 'zero', 'read_from_file'
    character(len=1024) :: pc_choice_initialise_EAS_config              = 'zero'
    character(len=1024) :: pc_choice_initialise_GRL_config              = 'zero'
    character(len=1024) :: pc_choice_initialise_ANT_config              = 'zero'
    ! Paths to files containing initial fields & values for the p/c scheme
    character(len=1024) :: filename_pc_initialise_NAM_config            = ''
    character(len=1024) :: filename_pc_initialise_EAS_config            = ''
    character(len=1024) :: filename_pc_initialise_GRL_config            = ''
    character(len=1024) :: filename_pc_initialise_ANT_config            = ''
    ! Timeframes to read from the p/c scheme initial file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_pc_initialise_NAM_config           = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    real(dp)            :: timeframe_pc_initialise_EAS_config           = 1E9_dp
    real(dp)            :: timeframe_pc_initialise_GRL_config           = 1E9_dp
    real(dp)            :: timeframe_pc_initialise_ANT_config           = 1E9_dp

  ! == Ice dynamics - calving
  ! =========================

    character(len=1024) :: choice_calving_law_config                    = 'threshold_thickness'            ! Choice of calving law: "none", "threshold_thickness"
    real(dp)            :: calving_threshold_thickness_shelf_config     = 100._dp                          ! Threshold ice thickness for ice shelf calving front in the "threshold_thickness" calving law
    real(dp)            :: calving_threshold_thickness_sheet_config     = 0._dp                            ! Threshold ice thickness for ice sheet calving front in the "threshold_thickness" calving law
    integer             :: max_calving_rounds_config                    = 20                               ! Maximum number of calving loops during chain reaction
    logical             :: do_remove_shelves_config                     = .false.                          ! If set to TRUE, all floating ice is always instantly removed (used in the ABUMIP-ABUK experiment)
    logical             :: remove_shelves_larger_than_PD_config         = .false.                          ! If set to TRUE, all floating ice beyond the present-day calving front is removed (used for some Antarctic spin-ups)
    logical             :: continental_shelf_calving_config             = .false.                          ! If set to TRUE, all ice beyond the continental shelf edge (set by a maximum depth) is removed
    real(dp)            :: continental_shelf_min_height_config          = -2000._dp                        ! Maximum depth of the continental shelf

    logical             :: do_apply_ISMIP7_fracture_mask_config         = .false.                          ! Whether or not to apply the ISMIP7 mask-based hydrofracturing forcing
    character(len=1024) :: filename_ISMIP7_fracture_mask_config         = ''                               ! The full path to the NetCDF file containing the ISMIP7 hydrofracturing mask

  ! == Ice dynamics - stabilisation
  ! ===============================

    character(len=1024) :: choice_mask_noice_config                     = 'none'                           ! Choice of mask_noice configuration
    real(dp)            :: Hi_min_config                                = 0._dp                            ! [m] Minimum ice thickness: thinner ice gets temporarily added to the no-ice mask and eventually removed
    real(dp)            :: Hi_thin_config                               = 0._dp                            ! [m] Thin-ice thickness threshold: thinner ice is considered dynamically unreliable (e.g. over steep slopes)
    logical             :: remove_ice_absent_at_PD_config               = .false.                          ! If set to TRUE, all ice not present in PD data is always instantly removed

    ! Geometry relaxation
    real(dp)            :: geometry_relaxation_t_years_config           = 0.0                              ! [yr] Amount of years (out-of-time) during which the geometry will relaxed during initialisation: no mass balance, thickness alterations, inversions, or target thinning rates will be applied during this time

    ! Mask conservation
    logical             :: do_protect_grounded_mask_config              = .false.                          ! If set to TRUE, grounded ice will not be allowed to cross the floatation threshold and will stay minimally grounded
    real(dp)            :: protect_grounded_mask_t_end_config           = -9.9E9_dp                        ! End time of grounded mask protection. After this time, it will be allowed to thin and become afloat


    ! Fix/delay ice thickness evolution
    logical             :: do_fixiness_before_start_config              = .false.                          ! Whether or not to apply fixiness values before fixiness_t_start
    real(dp)            :: fixiness_t_start_config                      = +9.9E9_dp                        ! Start time of linear transition between fixed/delayed and free evolution
    real(dp)            :: fixiness_t_end_config                        = +9.9E9_dp                        ! End   time of linear transition between fixed/delayed and free evolution
    real(dp)            :: fixiness_H_gl_gr_config                      = 0._dp                            ! Fix (1), release (0), or delay grunding line (grounded side) ice geometry evolution
    real(dp)            :: fixiness_H_gl_fl_config                      = 0._dp                            ! Fix (1), release (0), or delay grunding line (floating side) ice geometry evolution
    real(dp)            :: fixiness_H_grounded_config                   = 0._dp                            ! Fix (1), release (0), or delay grounded ice geometry evolution
    real(dp)            :: fixiness_H_floating_config                   = 0._dp                            ! Fix (1), release (0), or delay floating ice geometry evolution
    logical             :: fixiness_H_freeland_config                   = .false.                          ! Fix (.true.) ice-free vertices between fixiness_t_start and fixiness_t_end
    logical             :: fixiness_H_freeocean_config                  = .false.                          ! Fix (.true.) ice-free vertices between fixiness_t_start and fixiness_t_end

    ! Limit ice thickness evolution
    logical             :: do_limitness_before_start_config             = .false.                          ! Whether or not to apply limitness values before limitness_t_start
    real(dp)            :: limitness_t_start_config                     = +9.9E9_dp                        ! Start time of linear transition between limited and free evolution
    real(dp)            :: limitness_t_end_config                       = +9.9E9_dp                        ! End   time of linear transition between limited and free evolution
    real(dp)            :: limitness_H_gl_gr_config                     = +9.9E9_dp                        ! Maximum departure from PD ice thickness allowed at peak limitness
    real(dp)            :: limitness_H_gl_fl_config                     = +9.9E9_dp                        ! Maximum departure from PD ice thickness allowed at peak limitness
    real(dp)            :: limitness_H_grounded_config                  = +9.9E9_dp                        ! Maximum departure from PD ice thickness allowed at peak limitness
    real(dp)            :: limitness_H_floating_config                  = +9.9E9_dp                        ! Maximum departure from PD ice thickness allowed at peak limitness
    character(len=1024) :: modiness_H_style_config                      = 'none'                           ! Dynamic "modiness" limitness-like term: "none", "Ti_hom", "Ti_hom_up", "Ti_hom_down"
    real(dp)            :: modiness_T_hom_ref_config                    = 20._dp                           ! Reference Ti_hom for the modiness cold exponential style

  ! == Basal hydrology
  ! ==================

    ! Time step
    logical             :: do_asynchronous_basal_hydro_config           = .true.                           ! Whether or not the basal hydrology should be calculated asynchronously from the rest of the model; if so, use dt_basal_hydro; if not, calculate it in every time step
    real(dp)            :: dt_basal_hydro_config                        = 10._dp                           ! [yr] Time step for calculating basal hydrology

    ! Basal hydrology
    character(len=1024) :: choice_basal_hydrology_model_config          = 'Martin2011'                     ! Choice of basal hydrology model: "none", "Martin2011", "inversion", "read_from_file"
    real(dp)            :: Martin2011_hydro_Hb_min_config               = 0._dp                            ! Martin et al. (2011) basal hydrology model: low-end  Hb  value of bedrock-dependent pore-water pressure
    real(dp)            :: Martin2011_hydro_Hb_max_config               = 1000._dp                         ! Martin et al. (2011) basal hydrology model: high-end Hb  value of bedrock-dependent pore-water pressure
    real(dp)            :: basal_hydro_equil_time_config                = 0.1_dp                          ! [yr] time scale for basal hydrology to get to equilibrium
    real(dp)            :: error_function_max_effective_pressure_config = 5E6_dp                           ! Maximum effective pressure inland for the error-function model
    real(dp)            :: Leguy2014_hydro_connect_exponent_config      = 1._dp                            ! Leguy et al. (2014) hydrological connectivity of the subglacial hydrology drainage system
    real(dp)            :: Salle2025_hydro_meltwater_config             = 2.186523556E-7_dp                ! [kg m^-2 s^-1] of meltwater added to basal hydrology system
    real(dp)            :: Salle2025_W_til_max_config                   = 2.0_dp                           ! [m] Maximum allowed till water layer thickness
    real(dp)            :: Salle2025_initial_W_til_config               = 2.0_dp                           ! [m] Initial till water layer thickness
    real(dp)            :: Salle2025_initial_W_config                   = 0.01_dp                          ! [m] Initial subglacial water layer thickness
    real(dp)            :: Salle2025_Cd_config                          = 3.168874718E-11_dp               ! [m s^-1] Water leaking back from the till to the water layer above
    real(dp)            :: Salle2025_phi_config                         = 26.565_dp                        ! [degrees] till yield stress angle

    character(len=1024) :: choice_BMB_grounded_config                   = 'none'

  ! == Bed roughness
  ! ==================

    character(len=1024) :: choice_bed_roughness_config                  = 'uniform'                        ! Choice of source for friction coefficients: "uniform", "parameterised", "read_from_file"
    character(len=1024) :: choice_bed_roughness_parameterised_config    = ''                               ! "Martin2011", "SSA_icestream", "MISMIP+", "BIVMIP_A", "BIVMIP_B", "BIVMIP_C"
    ! Paths to files containing bed roughness fields for the chosen sliding law
    character(len=1024) :: filename_bed_roughness_NAM_config            = ''
    character(len=1024) :: filename_bed_roughness_EAS_config            = ''
    character(len=1024) :: filename_bed_roughness_GRL_config            = ''
    character(len=1024) :: filename_bed_roughness_ANT_config            = ''
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_bed_roughness_NAM_config           = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    real(dp)            :: timeframe_bed_roughness_EAS_config           = 1E9_dp
    real(dp)            :: timeframe_bed_roughness_GRL_config           = 1E9_dp
    real(dp)            :: timeframe_bed_roughness_ANT_config           = 1E9_dp
    ! Values for uniform bed roughness
    real(dp)            :: slid_Weertman_beta_sq_uniform_config         = 1.0E4_dp                         ! Uniform value for beta_sq  in Weertman sliding law
    real(dp)            :: slid_Coulomb_phi_fric_uniform_config         = 15._dp                           ! Uniform value for phi_fric in Coulomb sliding law
    real(dp)            :: slid_Budd_phi_fric_uniform_config            = 15._dp                           ! Uniform value for phi_fric in Budd sliding law
    real(dp)            :: slid_Tsai2015_alpha_sq_uniform_config        = 0.5_dp                           ! Uniform value for alpha_sq in the Tsai2015 sliding law
    real(dp)            :: slid_Tsai2015_beta_sq_uniform_config         = 1.0E4_dp                         ! Uniform value for beta_sq  in the Tsai2015 sliding law
    real(dp)            :: slid_Schoof2005_alpha_sq_uniform_config      = 0.5_dp                           ! Uniform value for alpha_sq in the Schoof2005 sliding law
    real(dp)            :: slid_Schoof2005_beta_sq_uniform_config       = 1.0E4_dp                         ! Uniform value for beta_sq  in the Schoof2005 sliding law
    real(dp)            :: slid_ZI_phi_fric_uniform_config              = 15._dp                           ! Uniform value for phi_fric in the Zoet-Iverson sliding law
    ! Parameters for bed roughness parameterisations
    real(dp)            :: Martin2011till_phi_Hb_min_config             = -1000._dp                        ! Martin et al. (2011) bed roughness parameterisation: low-end  Hb  value of bedrock-dependent till friction angle
    real(dp)            :: Martin2011till_phi_Hb_max_config             = 0._dp                            ! Martin et al. (2011) bed roughness parameterisation: high-end Hb  value of bedrock-dependent till friction angle
    real(dp)            :: Martin2011till_phi_min_config                = 5._dp                            ! Martin et al. (2011) bed roughness parameterisation: low-end  phi value of bedrock-dependent till friction angle
    real(dp)            :: Martin2011till_phi_max_config                = 20._dp                           ! Martin et al. (2011) bed roughness parameterisation: high-end phi value of bedrock-dependent till friction angle

  ! == Bed roughness inversion by nudging
  ! =====================================

    ! General
    logical             :: do_bed_roughness_nudging_config              = .false.                          !           Whether or not to nudge the basal roughness
    character(len=1024) :: choice_bed_roughness_nudging_method_config   = ''                               !           Choice of bed roughness nudging method
    character(len=1024) :: choice_inversion_target_geometry_config      = ''                               !           Which reference geometry to use as the target geometry (can be "init" or "PD")
    real(dp)            :: bed_roughness_nudging_dt_config              = 5._dp                            ! [yr]      Time step for bed roughness updates
    real(dp)            :: bed_roughness_nudging_t_start_config         = -9.9E9_dp                        ! [yr]      Earliest model time when nudging is allowed
    real(dp)            :: bed_roughness_nudging_t_end_config           = +9.9E9_dp                        ! [yr]      Latest   model time when nudging is allowed
    real(dp)            :: generic_bed_roughness_min_config             = 0.1_dp                           ! [?]       Smallest allowed value for the first  inverted bed roughness field
    real(dp)            :: generic_bed_roughness_max_config             = 30._dp                           ! [?]       Largest  allowed value for the first  inverted bed roughness field

    ! Bed roughness nudging model based on flowline-averaged values of H and dH/dt
    real(dp)            :: bednudge_H_dHdt_flowline_t_scale_config      = 100._dp                          ! [yr]      Timescale
    real(dp)            :: bednudge_H_dHdt_flowline_dH0_config          = 100._dp                          ! [m]       Ice thickness error scale
    real(dp)            :: bednudge_H_dHdt_flowline_dHdt0_config        = 0.6_dp                           ! [m yr^-1] Thinning rate scale
    real(dp)            :: bednudge_H_dHdt_flowline_Hi_scale_config     = 300._dp                          ! [m]       Ice thickness weight scale
    real(dp)            :: bednudge_H_dHdt_flowline_u_scale_config      = 3000._dp                         ! [m yr^-1] Ice velocity  weight scale
    real(dp)            :: bednudge_H_dHdt_flowline_r_smooth_config     = 2500._dp                         ! [m]       Radius for Gaussian filter used to smooth dC/dt as regularisation
    real(dp)            :: bednudge_H_dHdt_flowline_w_smooth_config     = 0.5_dp                           ! [-]       Relative contribution of smoothed dC/dt in regularisation

    ! Bed roughness nudging model based on local values of H and dH/dt (i.e. CISM method)
    real(dp)            :: bednudge_H_dHdt_local_H0_config              = 200._dp                          ! [m]       Ice thickness error scale
    real(dp)            :: bednudge_H_dHdt_local_tau_config             = 50._dp                           ! [yr]      Time scale
    real(dp)            :: bednudge_H_dHdt_local_L_config               = 8000._dp                         ! [m]       Length scale for the Laplacian regularisation term

    ! Bed roughness nudging model based on flowline-averaged values of H and u
    CHARACTER(len=1024) :: bednudge_H_u_flowline_file_u_target_config   = ''                               !           File containing target ice velocities
    real(dp)            :: bednudge_H_u_flowline_t_scale_config         = 10._dp                           ! [yr]      Timescale
    real(dp)            :: bednudge_H_u_flowline_H0_config              = 100._dp                          ! [m]       Ice thickness error scale
    real(dp)            :: bednudge_H_u_flowline_u0_config              = 250_dp                           ! [m yr^-1] Ice velocity  error scale
    real(dp)            :: bednudge_H_u_flowline_Hi_scale_config        = 300._dp                          ! [m]       Ice thickness weight scale
    real(dp)            :: bednudge_H_u_flowline_u_scale_config         = 3000._dp                         ! [m yr^-1] Ice velocity  weight scale
    real(dp)            :: bednudge_H_u_flowline_tau_config             = 50._dp                           ! [yr]      Time scale
    real(dp)            :: bednudge_H_u_flowline_L_config               = 2000._dp                         ! [m]       Length scale for the Laplacian regularisation term

  ! == Geothermal heat flux
  ! =======================

    character(len=1024) :: choice_geothermal_heat_flux_config           = 'uniform'                         ! Choice of geothermal heat flux; can be 'uniform' or 'read_from_file'
    real(dp)            :: uniform_geothermal_heat_flux_config          = 1.72E06_dp                        ! [J m^-2 yr^-1] Value when choice_geothermal_heat_flux == 'uniform' (1.72E06 J m^-2 yr^-1 according to Sclater et al. (1980))
    character(len=1024) :: filename_geothermal_heat_flux_config         = ''                                ! Expected units: J m^-2 s^-1 (conversion done internally). External GHF data file.

  ! == Thermodynamics
  ! =================

    ! Initial temperature profile
    character(len=1024) :: choice_initial_ice_temperature_NAM_config    = 'Robin'                          ! Choice of initial ice temperature profile: "uniform", "linear", "Robin", "read_from_file"
    character(len=1024) :: choice_initial_ice_temperature_EAS_config    = 'Robin'
    character(len=1024) :: choice_initial_ice_temperature_GRL_config    = 'Robin'
    character(len=1024) :: choice_initial_ice_temperature_ANT_config    = 'Robin'
    ! Uniform initial ice temperature, if choice_initial_ice_temperature == 'uniform'
    real(dp)            :: uniform_initial_ice_temperature_NAM_config   = 270._dp                          ! Uniform ice temperature (applied when choice_initial_ice_temperature_config = "uniform")
    real(dp)            :: uniform_initial_ice_temperature_EAS_config   = 270._dp
    real(dp)            :: uniform_initial_ice_temperature_GRL_config   = 270._dp
    real(dp)            :: uniform_initial_ice_temperature_ANT_config   = 270._dp
    ! Paths to files containing initial temperature fields, if
    character(len=1024) :: filename_initial_ice_temperature_NAM_config  = ''
    character(len=1024) :: filename_initial_ice_temperature_EAS_config  = ''
    character(len=1024) :: filename_initial_ice_temperature_GRL_config  = ''
    character(len=1024) :: filename_initial_ice_temperature_ANT_config  = ''
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_initial_ice_temperature_NAM_config = 1E9_dp                           ! Can be different from C%start_time_of_run, be careful though!
    real(dp)            :: timeframe_initial_ice_temperature_EAS_config = 1E9_dp
    real(dp)            :: timeframe_initial_ice_temperature_GRL_config = 1E9_dp
    real(dp)            :: timeframe_initial_ice_temperature_ANT_config = 1E9_dp
    ! Thermodynamical model
    character(len=1024) :: choice_thermo_model_config                   = '3D_heat_equation'               ! Choice of thermodynamical model: "none", "3D_heat_equation"
    real(dp)            :: dt_thermodynamics_config                     = 1._dp                            ! [yr] Time step for the thermodynamical model
    real(dp)            :: Hi_min_thermo_config                         = 10._dp                           ! [m]  Ice thinner than this is assumed to have a temperature equal to the annual mean surface temperature throughout the vertical column
    character(len=1024) :: choice_GL_temperature_BC_config              = 'grounded'                       ! Choice of basal boundary condition at grounding lines: 'grounded', 'subgrid', 'pmp'
    character(len=1024) :: choice_ice_heat_capacity_config              = 'Pounder1965'                    ! Choice of ice heat capacity model: "uniform", "Pounder1965"
    real(dp)            :: uniform_ice_heat_capacity_config             = 2009._dp                         ! Uniform ice heat capacity (applied when choice_ice_heat_capacity_config = "uniform")
    character(len=1024) :: choice_ice_thermal_conductivity_config       = 'Ritz1987'                       ! Choice of ice heat capacity model: "uniform", "Ritz1987"
    real(dp)            :: uniform_ice_thermal_conductivity_config      = 6.626958E7_dp                    ! Uniform ice thermal conductivity (applied when choice_ice_thermal_conductivity_config = "uniform")

  ! == Rheology and flow law
  ! =========================

    ! Flow law
    character(len=1024) :: choice_flow_law_config                       = 'Glen'                           ! Choice of flow law, relating effective viscosity to effective strain rate
    real(dp)            :: Glens_flow_law_exponent_config               = 3.0_dp                           ! Exponent in Glen's flow law
    real(dp)            :: Glens_flow_law_epsilon_sq_0_config           = 1E-8_dp                          ! Normalisation term so that zero strain rates produce a high but finite viscosity

    ! Rheology
    character(len=1024) :: choice_ice_rheology_Glen_config              = 'Huybrechts1992'                 ! Choice of ice rheology model for Glen's flow law: "uniform", "Huybrechts1992", "MISMIP_mod"
    real(dp)            :: uniform_Glens_flow_factor_config             = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")

    ! Enhancement factors
    character(len=1024) :: choice_enhancement_factor_transition_config  = 'separate'                       ! Choice of treatment of enhancement factors at transition zones: "separate", "interp"
    real(dp)            :: m_enh_sheet_config                           = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
    real(dp)            :: m_enh_shelf_config                           = 1.0_dp                           ! Ice flow enhancement factor for floating ice

  ! == Climate
  ! ==========

    ! Time step
    logical             :: do_asynchronous_climate_config               = .true.                           ! Whether or not the climate should be calculated asynchronously from the rest of the model; if so, use dt_climate; if not, calculate it in every time step
    real(dp)            :: dt_climate_config                            = 10._dp                           ! [yr] Time step for calculating climate

    ! Choice of climate model
    character(len=1024) :: choice_climate_model_NAM_config              = 'none'
    character(len=1024) :: choice_climate_model_EAS_config              = 'none'
    character(len=1024) :: choice_climate_model_GRL_config              = 'none'
    character(len=1024) :: choice_climate_model_ANT_config              = 'none'

    ! Choice of idealised climate model
    character(len=1024) :: choice_climate_model_idealised_config        = ''

    ! Choice of realistic climate model
    character(len=1024) :: choice_climate_model_realistic_config        = ''

    ! Paths to files containing fields for realistic climates
    character(len=1024) :: filename_climate_snapshot_NAM_config         = ''
    character(len=1024) :: filename_climate_snapshot_EAS_config         = ''
    character(len=1024) :: filename_climate_snapshot_GRL_config         = ''
    character(len=1024) :: filename_climate_snapshot_ANT_config         = ''

    ! Climate - global forcings
    character(len=1024) :: filename_CO2_record_config                   = ''                            ! CO2 record (netCDF file)
    character(len=1024) :: filename_d18O_record_config                  = ''                            ! d18O record (netCDF file)

    ! == Climate - fixed regional lapse rates
    logical             :: do_lapse_rate_corrections_NAM_config         = .false.                          ! whether or not to apply the lapse rates below - they seem to produce much higher SMB at the ice sheet fringes
    logical             :: do_lapse_rate_corrections_EAS_config         = .false.
    logical             :: do_lapse_rate_corrections_GRL_config         = .false.
    logical             :: do_lapse_rate_corrections_ANT_config         = .false.
    real(dp)            :: lapse_rate_temp_NAM_config                   = 7.9E-3_dp                          ! Elevation lapse rate effect on temperature [K m^-1]
    real(dp)            :: lapse_rate_temp_EAS_config                   = 7.9E-3_dp                          !
    real(dp)            :: lapse_rate_temp_GRL_config                   = 7.9E-3_dp                          !
    real(dp)            :: lapse_rate_temp_ANT_config                   = 7.9E-3_dp                          !
    real(dp)            :: lapse_rate_precip_NAM_config                 = 0.07_dp                            ! Elevation-desertification lapse rate [K^-1]
    real(dp)            :: lapse_rate_precip_EAS_config                 = 0.07_dp                            !
    real(dp)            :: lapse_rate_precip_GRL_config                 = 0.07_dp                            !
    real(dp)            :: lapse_rate_precip_ANT_config                 = 0.07_dp                            !

    ! == Climate - snapshot plus a uniform deltaT
    character(len=1024) :: filename_climate_snapshot_unif_dT_NAM_config = ''
    character(len=1024) :: filename_climate_snapshot_unif_dT_EAS_config = ''
    character(len=1024) :: filename_climate_snapshot_unif_dT_GRL_config = ''
    character(len=1024) :: filename_climate_snapshot_unif_dT_ANT_config = ''
    real(dp)            :: uniform_deltaT_NAM_config                    = 0.0_dp                          ! Temperature anomaly to be added to the T2m field [K]
    real(dp)            :: uniform_deltaT_EAS_config                    = 0.0_dp                          !
    real(dp)            :: uniform_deltaT_GRL_config                    = 0.0_dp                          !
    real(dp)            :: uniform_deltaT_ANT_config                    = 0.0_dp                          !
    real(dp)            :: precip_CC_correction_NAM_config              = 1.068_dp                        ! Precipitation correction due to increased temperatures [1.0 = no correction; 1.068 = 6.8%]
    real(dp)            :: precip_CC_correction_EAS_config              = 1.068_dp                        !
    real(dp)            :: precip_CC_correction_GRL_config              = 1.068_dp                        !
    real(dp)            :: precip_CC_correction_ANT_config              = 1.068_dp                        !

    ! == Climate - snapshot plus a transient deltaT
    character(len=1024) :: filename_climate_snapshot_trans_dT_NAM_config = ''
    character(len=1024) :: filename_climate_snapshot_trans_dT_EAS_config = ''
    character(len=1024) :: filename_climate_snapshot_trans_dT_GRL_config = ''
    character(len=1024) :: filename_climate_snapshot_trans_dT_ANT_config = ''
    character(len=1024) :: filename_atmosphere_dT_NAM_config             = ''
    character(len=1024) :: filename_atmosphere_dT_EAS_config             = ''
    character(len=1024) :: filename_atmosphere_dT_GRL_config             = ''
    character(len=1024) :: filename_atmosphere_dT_ANT_config             = ''

    ! == Climate snapshot_plus_anomalies
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_NAM_config   = ''
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_EAS_config   = ''
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_GRL_config   = ''
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_ANT_config   = ''
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_EAS_config  = ''
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_NAM_config  = ''
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_GRL_config  = ''
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_ANT_config  = ''

    ! == Climate - Insolation
    character(len=1024) :: choice_insolation_forcing_config             = 'none'                           ! 'none', 'static' or 'realistic'
    character(len=1024) :: filename_insolation_config                   = ''                               ! File with the insolation solution (Laskar 2004)
    real(dp)            :: static_insolation_time_config                = 0._dp                            ! [ka?] time to use for a static insolation

    ! == Climate matrix
    character(len=1024) :: choice_matrix_forcing_config                 = 'none'                           ! 'none', 'CO2_direct' or 'd18O_inverse_CO2'
    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    character(len=1024) :: climate_matrix_filename_PD_obs_climate_config        = ''

    ! GCM snapshots in the matrix_warm_cold option
    character(len=1024) :: climate_matrix_filename_climate_snapshot_PI_config   = ''
    character(len=1024) :: climate_matrix_filename_climate_snapshot_warm_config = ''
    character(len=1024) :: climate_matrix_filename_climate_snapshot_cold_config = ''

    real(dp)            :: climate_matrix_constant_lapserate_config    = 0.008_dp                         ! Constant atmospheric lapse rate [K m^-1]

    ! Scaling factor for CO2 vs ice weights
    real(dp)            :: climate_matrix_CO2vsice_NAM_config          = 0.5_dp                           ! Weight factor for the influence of CO2 vs ice cover on temperature
    real(dp)            :: climate_matrix_CO2vsice_EAS_config          = 0.5_dp                           ! Can be set separately for different regions
    real(dp)            :: climate_matrix_CO2vsice_GRL_config          = 0.75_dp                          ! Default values are from Berends et al, 2018
    real(dp)            :: climate_matrix_CO2vsice_ANT_config          = 0.75_dp                          ! 1.0_dp equals glacial index method

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    real(dp)            :: climate_matrix_high_CO2_level_config        = 280._dp                          ! CO2 level  pertaining to the warm climate (PI  level default)
    real(dp)            :: climate_matrix_low_CO2_level_config         = 190._dp                          ! CO2 level  pertaining to the cold climate (LGM level default)
    real(dp)            :: climate_matrix_warm_orbit_time_config       = 0._dp                            ! Orbit time pertaining to the warm climate (PI default)
    real(dp)            :: climate_matrix_cold_orbit_time_config       = -21000._dp                       ! Orbit time pertaining to the cold climate (LGM default)

    ! Whether or not to apply a bias correction to the GCM snapshots
    logical             :: climate_matrix_biascorrect_warm_config      = .true.                           ! Whether or not to apply a bias correction (modelled vs observed PI climate) to the "warm" GCM snapshot
    logical             :: climate_matrix_biascorrect_cold_config      = .true.                           ! Whether or not to apply a bias correction (modelled vs observed PI climate) to the "cold" GCM snapshot

    logical             :: climate_matrix_switch_glacial_index_precip_config = .false.                    ! If a glacial index is used for the precipitation forcing, it will only depend on CO2

    ! Settings for the ISMIP7 climate model
    character(len=1024) :: climate_ISMIP7_choice_baseline_config         = ''                               ! How to define the baseline climate for the anomalies: 'yearly' (i.e. use the provided yearly fields) or 'fixed' (i.e. use a separate, time-independent climate - probably the same present-day climate that was used for the initialisation)
    character(len=1024) :: climate_ISMIP7_filename_baseline_config       = ''                               ! Path to the separate, time-independent climate - probably the same present-day climate that was used for the initialisation
    character(len=1024) :: climate_ISMIP7_choice_refgeo_config           = ''                               ! Which reference geometry to use as the baseline for calculating delta_ts = dts/dz * delta_s: 'init', 'PD'
    character(len=1024) :: climate_ISMIP7_forcing_foldername_config          = ''                               ! Path to the directory containing the different variables directories (e.g. /path/to/base/folder, so that the climate files are located in /path/to/base/folder/acabf/version)
    character(len=1024) :: climate_ISMIP7_forcing_version_config             = ''                               ! Which version of the forcing files to use (since they often provide more than one), e.g. 'v2' means the climate files are located in /path/to/base/folder/acabf/v2. Leaving this variable empty implies that they are located in /path/to/base/folder/acabf


  ! == Ocean
  ! ========

    ! Time step
    logical             :: do_asynchronous_ocean_config                 = .true.                           ! Whether or not the ocean should be calculated asynchronously from the rest of the model; if so, use dt_climate; if not, calculate it in every time step
    real(dp)            :: dt_ocean_config                              = 10._dp                           ! [yr] Time step for calculating ocean

    ! Vertical grid
    real(dp)            :: ocean_vertical_grid_max_depth_config         = 1500._dp                         ! Maximum depth of the ocean submodel
    real(dp)            :: ocean_vertical_grid_dz_config                = 100._dp                          ! Vertical distance between ocean layers

    ! Choice of ocean model
    character(len=1024) :: choice_ocean_model_NAM_config                = 'none'
    character(len=1024) :: choice_ocean_model_EAS_config                = 'none'
    character(len=1024) :: choice_ocean_model_GRL_config                = 'none'
    character(len=1024) :: choice_ocean_model_ANT_config                = 'none'

    ! Choice of idealised ocean model
    character(len=1024) :: choice_ocean_model_idealised_config          = ''                               ! Choice of idealised ocean forcing: 'ISOMIP', 'TANH', 'LINEAR', 'uniform'
    character(len=1024) :: choice_ocean_isomip_scenario_config          = ''                               ! Scenario when using 'ISOMIP' forcing: 'WARM' or 'COLD'
    real(dp)            :: ocean_tanh_deep_temperature_config           = 1.0_dp                           ! [degC] Deep ocean temperature when using 'TANH' forcing
    real(dp)            :: ocean_tanh_thermocline_depth_config          = 100.0_dp                         ! [m]    Depth of thermocline when using 'TANH' forcing
    real(dp)            :: ocean_tanh_thermocline_scale_depth_config    = 100.0_dp                         ! [m]    Scale depth for thermocline sharpness when using 'TANH' forcing
    real(dp)            :: ocean_linear_deep_temperature_config         = -2.3_dp                          ! [degC] Deep ocean temperature when using 'LINEAR' forcing
    real(dp)            :: ocean_linear_deep_salinity_config            = 34.8_dp                          ! [psu] Deep ocean salinity when using 'LINEAR' forcing
    real(dp)            :: ocean_linear_reference_depth_config          = 2000.0_dp                        ! [m] Depth where deep values are prescribed
    real(dp)            :: ocean_lin_therm_surf_salinity_config         = 34.0_dp                          ! [psu] Surface salinity when using 'LINEAR_THERMOCLINE'
    real(dp)            :: ocean_lin_therm_deep_salinity_config         = 34.7_dp                          ! [psu] Deep salinity when using 'LINEAR_THERMOCLINE'
    real(dp)            :: ocean_lin_therm_surf_temperature_config      = -1.0_dp                          ! [degC] Surface temperature when using 'LINEAR_THERMOCLINE'
    real(dp)            :: ocean_lin_therm_deep_temperature_config      = 1.2_dp                           ! [degC] Deep temperature when using 'LINEAR_THERMOCLINE'
    real(dp)            :: ocean_lin_therm_thermocline_top_config       = 200.0_dp                         ! [m] Top of thermocline depth when using 'LINEAR_THERMOCLINE'
    real(dp)            :: ocean_lin_therm_thermocline_bottom_config    = 600.0_dp                         ! [m] Bottom of thermocline depth when using 'LINEAR_THERMOCLINE'
    real(dp)            :: ocean_uniform_T_config                       = 0.0_dp                           ! [degC] Uniform ocean temperature when using 'uniform' forcing
    real(dp)            :: ocean_uniform_S_config                       = 33.8_dp                          ! [psu] Uniform ocean salinity when using 'uniform' forcing

    ! Choice of realistic ocean model
    character(len=1024) :: choice_ocean_model_realistic_config          = ''

    ! Paths to files containing fields for realistic ocean
    character(len=1024) :: filename_ocean_snapshot_NAM_config           = ''
    character(len=1024) :: filename_ocean_snapshot_EAS_config           = ''
    character(len=1024) :: filename_ocean_snapshot_GRL_config           = ''
    character(len=1024) :: filename_ocean_snapshot_ANT_config           = ''

    character(len=1024) :: filename_ocean_warm_snapshot_NAM_config      = ''
    character(len=1024) :: filename_ocean_warm_snapshot_EAS_config      = ''
    character(len=1024) :: filename_ocean_warm_snapshot_GRL_config      = ''
    character(len=1024) :: filename_ocean_warm_snapshot_ANT_config      = ''

    character(len=1024) :: filename_ocean_cold_snapshot_NAM_config      = ''
    character(len=1024) :: filename_ocean_cold_snapshot_EAS_config      = ''
    character(len=1024) :: filename_ocean_cold_snapshot_GRL_config      = ''
    character(len=1024) :: filename_ocean_cold_snapshot_ANT_config      = ''

    ! Choice of ocean deltaT
    real(dp)            :: ocean_uniform_deltaT_NAM_config              = 0.0_dp
    real(dp)            :: ocean_uniform_deltaT_EAS_config              = 0.0_dp
    real(dp)            :: ocean_uniform_deltaT_GRL_config              = 0.0_dp
    real(dp)            :: ocean_uniform_deltaT_ANT_config              = 0.0_dp

    ! Choice of extrapolation method
    character(len=1024) :: choice_ocean_extrapolation_method_config     = 'initialisation'                 ! Method to extrapolate ocean forcing into cavities: 'initialisation'

    ! Choice of transient ocean model
    character(len=1024) :: choice_ocean_model_transient_config          = ''                                ! so far only 'deltaT' implemented

    ! Paths to files containing the deltaT record for the transient deltaT ocean model
    character(len=1024) :: filename_ocean_dT_NAM_config                 = ''
    character(len=1024) :: filename_ocean_dT_EAS_config                 = ''
    character(len=1024) :: filename_ocean_dT_GRL_config                 = ''
    character(len=1024) :: filename_ocean_dT_ANT_config                 = ''

    ! Paths to files containing the GI record for the transient GlacialIndex ocean model
    character(len=1024) :: filename_ocean_GI_NAM_config                 = ''
    character(len=1024) :: filename_ocean_GI_EAS_config                 = ''
    character(len=1024) :: filename_ocean_GI_GRL_config                 = ''
    character(len=1024) :: filename_ocean_GI_ANT_config                 = ''

    ! Settings for the snapshot_plus_anomalies ocean model
    character(len=1024) :: ocean_snp_p_anml_filename_snapshot_config    = ''                               ! File containing the ocean snapshot (e.g. World Ocean Atlas)
    character(len=1024) :: ocean_snp_p_anml_filename_anomalies_config   = ''                               ! File containing the ocean anomalies (e.g. from a GCM projection)

    ! Settings for the ISMIP7 ocean model
    character(len=1024) :: ocean_ISMIP7_forcing_foldername_config       = ''                               ! Path to the directory containing the different variables directories (e.g. /path/to/base/folder, so that the ocean files are located in /path/to/base/folder/thetao/version)
    character(len=1024) :: ocean_ISMIP7_forcing_version_config          = ''

  ! == Surface mass balance
  ! =======================

    ! Time step
    logical             :: do_asynchronous_SMB_config                   = .true.                           ! Whether or not the SMB should be calculated asynchronously from the rest of the model; if so, use dt_climate; if not, calculate it in every time step
    real(dp)            :: dt_SMB_config                                = 10._dp                           ! [yr] Time step for calculating SMB

    ! Choice of SMB model
    character(len=1024) :: choice_SMB_model_NAM_config                  = 'uniform'
    character(len=1024) :: choice_SMB_model_EAS_config                  = 'uniform'
    character(len=1024) :: choice_SMB_model_GRL_config                  = 'uniform'
    character(len=1024) :: choice_SMB_model_ANT_config                  = 'uniform'

    ! Value to be used for uniform SMB (no regional variants, only used for idealised-geometry experiments)
    real(dp)            :: uniform_SMB_config                           = 0._dp                             ! [m.i.e./yr]

    ! Choice of idealised SMB model
    character(len=1024) :: choice_SMB_model_idealised_config            = ''

    ! Prescribed SMB forcing
    character(len=1024) :: choice_SMB_prescribed_NAM_config             = ''
    character(len=1024) :: choice_SMB_prescribed_EAS_config             = ''
    character(len=1024) :: choice_SMB_prescribed_GRL_config             = ''
    character(len=1024) :: choice_SMB_prescribed_ANT_config             = ''

    ! Files containing prescribed SMB forcing
    character(len=1024) :: filename_SMB_prescribed_NAM_config           = ''
    character(len=1024) :: filename_SMB_prescribed_EAS_config           = ''
    character(len=1024) :: filename_SMB_prescribed_GRL_config           = ''
    character(len=1024) :: filename_SMB_prescribed_ANT_config           = ''

    ! Timeframes for reading prescribed SMB forcing from file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_SMB_prescribed_NAM_config          = 1E9_dp
    real(dp)            :: timeframe_SMB_prescribed_EAS_config          = 1E9_dp
    real(dp)            :: timeframe_SMB_prescribed_GRL_config          = 1E9_dp
    real(dp)            :: timeframe_SMB_prescribed_ANT_config          = 1E9_dp

    character(len=1024) :: choice_SMB_IMAUITM_init_firn_NAM_config     = 'uniform'                        ! How to initialise the firn layer in the IMAU-ITM SMB model: "uniform", "restart"
    character(len=1024) :: choice_SMB_IMAUITM_init_firn_EAS_config     = 'uniform'
    character(len=1024) :: choice_SMB_IMAUITM_init_firn_GRL_config     = 'uniform'
    character(len=1024) :: choice_SMB_IMAUITM_init_firn_ANT_config     = 'uniform'

  ! Files containing the firn model (yearly firn depth and melt)
    character(len=1024) :: filename_firn_IMAUITM_NAM_config            = ''
    character(len=1024) :: filename_firn_IMAUITM_EAS_config            = ''
    character(len=1024) :: filename_firn_IMAUITM_GRL_config            = ''
    character(len=1024) :: filename_firn_IMAUITM_ANT_config            = ''

    ! timeframe for restarting from the firn model
    real(dp)            :: timeframe_restart_firn_IMAUITM_NAM_config           = 1E9_dp
    real(dp)            :: timeframe_restart_firn_IMAUITM_EAS_config           = 1E9_dp
    real(dp)            :: timeframe_restart_firn_IMAUITM_GRL_config           = 1E9_dp
    real(dp)            :: timeframe_restart_firn_IMAUITM_ANT_config           = 1E9_dp

    ! Tuning parameters for the IMAU-ITM SMB model
    real(dp)            :: SMB_IMAUITM_initial_firn_thickness_config   = 1._dp                            ! Initial firn thickness of the IMAU-ITEM SMB model [m] (used when SMB_IMAUITM_choice_init_firn = "uniform")
    real(dp)            :: SMB_IMAUITM_C_abl_constant_NAM_config       = -49._dp                          ! 34._dp    (commented values are old ANICE defaults, but since refreezing was not calculated right
    real(dp)            :: SMB_IMAUITM_C_abl_constant_EAS_config       = -49._dp                          !            and this has since been fixed, these values will still not give the same results as
    real(dp)            :: SMB_IMAUITM_C_abl_constant_GRL_config       = -49._dp                          !            they used to in ANICE.)
    real(dp)            :: SMB_IMAUITM_C_abl_constant_ANT_config       = -49._dp
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_NAM_config             = 10._dp                           ! 10._dp
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_EAS_config             = 10._dp
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_GRL_config             = 10._dp
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_ANT_config             = 10._dp
    real(dp)            :: SMB_IMAUITM_C_abl_Q_NAM_config              = 0.0227_dp                        ! 0.513_dp
    real(dp)            :: SMB_IMAUITM_C_abl_Q_EAS_config              = 0.0227_dp
    real(dp)            :: SMB_IMAUITM_C_abl_Q_GRL_config              = 0.0227_dp
    real(dp)            :: SMB_IMAUITM_C_abl_Q_ANT_config              = 0.0227_dp
    real(dp)            :: SMB_IMAUITM_C_refr_NAM_config               = 0.051_dp                         ! 0.012_dp
    real(dp)            :: SMB_IMAUITM_C_refr_EAS_config               = 0.051_dp
    real(dp)            :: SMB_IMAUITM_C_refr_GRL_config               = 0.051_dp
    real(dp)            :: SMB_IMAUITM_C_refr_ANT_config               = 0.051_dp
    real(dp)            :: SMB_IMAUITM_albedo_water_config              = 0.1_dp
    real(dp)            :: SMB_IMAUITM_albedo_soil_config               = 0.2_dp
    real(dp)            :: SMB_IMAUITM_albedo_ice_config                = 0.5_dp
    real(dp)            :: SMB_IMAUITM_albedo_snow_config               = 0.85_dp

    ! Settings for the snapshot_plus_anomalies SMB model
    character(len=1024) :: SMB_snp_p_anml_filename_snapshot_T2m_config  = ''                               ! File containing the T2m snapshot (e.g. from a RACMO historical simulation)
    character(len=1024) :: SMB_snp_p_anml_filename_snapshot_SMB_config  = ''                               ! File containing the SMB snapshot (e.g. from a RACMO historical simulation)
    character(len=1024) :: SMB_snp_p_anml_filename_anomalies_config     = ''                               ! File containing the SMB+T2m anomalies (e.g. from a GCM projection)

    ! Settings for the ISMIP7 SMB model
    character(len=1024) :: SMB_ISMIP7_choice_SMB_baseline_config         = ''                               ! How to define the baseline SMB for the anomalies: 'yearly' (i.e. use the provided yearly acabf fields) or 'fixed' (i.e. use a separate, time-independent SMB - probably the same present-day SMB that was used for the initialisation)
    character(len=1024) :: SMB_ISMIP7_filename_SMB_baseline_fixed_config = ''                               ! Path to the separate, time-independent SMB - probably the same present-day SMB that was used for the initialisation
    logical :: SMB_ISMIP7_choice_SMB_offset_config         = .false.                                        ! Whether or not to read and apply an offset to the SMB baseline
    character(len=1024) :: SMB_ISMIP7_filename_SMB_offset_config         = ''                               ! Path to the SMB offset which is subtracted from the SMB_basline_fixed to shift the baseline period to the correct 1960-1989 during which SMB anomalies are defined as zero
    character(len=1024) :: SMB_ISMIP7_choice_refgeo_config               = ''                               ! Which reference geometry to use as the baseline for calculating delta_SMB = dSMB/dz * delta_s: 'init', 'PD'
    character(len=1024) :: SMB_ISMIP7_forcing_foldername_config          = ''                               ! Path to the directory containing the different variables directories (e.g. /path/to/base/folder, so that the SMB files are located in /path/to/base/folder/acabf/version)
    character(len=1024) :: SMB_ISMIP7_forcing_version_config             = ''                               ! Which version of the forcing files to use (since they often provide more than one), e.g. 'v2' means the SMB files are located in /path/to/base/folder/acabf/v2. Leaving this variable empty implies that they are located in /path/to/base/folder/acabf

  ! == Basal mass balance
  ! =====================

    ! Time step
    logical             :: do_asynchronous_BMB_config                   = .true.                           ! Whether or not the BMB should be calculated asynchronously from the rest of the model; if so, use dt_climate; if not, calculate it in every time step
    real(dp)            :: dt_BMB_config                                = 10._dp                           ! [yr] Time step for calculating BMB

    real(dp)            :: dt_BMB_reinit_config                         = 1000._dp                         ! [yr] Time interval when BMB should be reinitialised

    ! Hard limits on melt/refreezing rates
    real(dp)            :: BMB_maximum_allowed_melt_rate_config         = 100._dp                         ! [m/yr] Maximum allowed melt       rate   (note: positive value means melt!)
    real(dp)            :: BMB_maximum_allowed_refreezing_rate_config   = 10._dp                         ! [m/yr] Maximum allowed refreezing rate   (note: positive value means refreezing!)

    ! BMB transition phase
    logical             :: do_BMB_transition_phase_config               = .false.                          ! Whether or not the model should slowly transition from inverted BMB to modelled BMB over a specified time window (only applied when do_BMB_transition_phase_config = .true.)
    real(dp)            :: BMB_transition_phase_t_start_config          = +9.8E9_dp                        ! [yr] Start time for BMB transition phase
    real(dp)            :: BMB_transition_phase_t_end_config            = +9.9E9_dp                        ! [yr] End   time for BMB transition phase

    ! Grounding line treatment
    character(len=1024) :: choice_BMB_subgrid_config                    = 'NMP'                            ! Choice of sub-grid BMB scheme: "NMP", "FCMP", "PMP" (following Leguy et al., 2021)

    ! Choice of BMB model
    character(len=1024) :: choice_BMB_model_NAM_config                  = 'uniform'
    character(len=1024) :: choice_BMB_model_EAS_config                  = 'uniform'
    character(len=1024) :: choice_BMB_model_GRL_config                  = 'uniform'
    character(len=1024) :: choice_BMB_model_ANT_config                  = 'uniform'

    ! Choice of BMB model in ROI
    character(len=1024) :: choice_BMB_model_NAM_ROI_config              = 'identical_to_choice_BMB_model'  ! Choose BMB model in ROI, options: 'identical_to_choice_BMB_model', 'uniform', 'laddie_py'
    character(len=1024) :: choice_BMB_model_EAS_ROI_config              = 'identical_to_choice_BMB_model'  ! Choose BMB model in ROI, options: 'identical_to_choice_BMB_model', 'uniform', 'laddie_py'
    character(len=1024) :: choice_BMB_model_GRL_ROI_config              = 'identical_to_choice_BMB_model'  ! Choose BMB model in ROI, options: 'identical_to_choice_BMB_model', 'uniform', 'laddie_py'
    character(len=1024) :: choice_BMB_model_ANT_ROI_config              = 'identical_to_choice_BMB_model'  ! Choose BMB model in ROI, options: 'identical_to_choice_BMB_model', 'uniform', 'laddie_py'

    ! Prescribed BMB forcing
    character(len=1024) :: choice_BMB_prescribed_NAM_config             = ''
    character(len=1024) :: choice_BMB_prescribed_EAS_config             = ''
    character(len=1024) :: choice_BMB_prescribed_GRL_config             = ''
    character(len=1024) :: choice_BMB_prescribed_ANT_config             = ''

    ! Files containing prescribed BMB forcing
    character(len=1024) :: filename_BMB_prescribed_NAM_config           = ''
    character(len=1024) :: filename_BMB_prescribed_EAS_config           = ''
    character(len=1024) :: filename_BMB_prescribed_GRL_config           = ''
    character(len=1024) :: filename_BMB_prescribed_ANT_config           = ''

    ! Timeframes for reading prescribed BMB forcing from file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_BMB_prescribed_NAM_config          = 1E9_dp
    real(dp)            :: timeframe_BMB_prescribed_EAS_config          = 1E9_dp
    real(dp)            :: timeframe_BMB_prescribed_GRL_config          = 1E9_dp
    real(dp)            :: timeframe_BMB_prescribed_ANT_config          = 1E9_dp

    ! Choice of idealised BMB model
    character(len=1024) :: choice_BMB_model_idealised_config            = ''

    ! Choice of parameterised BMB model
    character(len=1024) :: choice_BMB_model_parameterised_config        = ''

    ! "uniform"
    real(dp)            :: uniform_BMB_config                           = 0._dp
    real(dp)            :: uniform_BMB_ROI_config                       = 0._dp
    real(dp)            :: uniform_BMB_t_start_config                   = 0._dp

    ! "parameterised"
    real(dp)            :: BMB_Favier2019_gamma_config                  = 99.32E-5
    real(dp)            :: BMB_Holland_Cmelt_config                     = 34.8_dp                          ! equivalent to 8.19e-5 if it was an exchange velocity gamma

    ! "laddie_py"
    character(len=1024) :: choice_BMB_laddie_system_config              = ''                               ! System on which the model is running: 'local_mac' or 'slurm_HPC'
    character(len=1024) :: filename_BMB_laddie_configname_config        = ''                               ! File name of basal melt provided by LADDIE
    character(len=1024) :: filename_BMB_laddie_initial_restart_config   = ''                               ! File name containing restart for laddie from laddie spinup
    character(len=1024) :: filename_BMB_laddie_initial_output_config    = ''                               ! File name containing output from laddie spinup
    character(len=1024) :: dir_BMB_laddie_model_config                  = ''                               ! Directory where laddie code is located
    character(len=1024) :: conda_activate_prompt_config                 = 'conda activate laddie'          ! Prompt to activate conda environment used for running laddie

    ! "inverted"
    real(dp)            :: BMB_inversion_t_start_config                 = -9.9E9                           ! [yr] Only nudge melt rates when the model time lies between t_start and t_end
    real(dp)            :: BMB_inversion_t_end_config                   =  9.9E9                           ! [yr]

  ! == LADDIE model
  ! ===============

    ! Parallellisation
    logical             :: do_repartition_laddie_config                 = .false.                          ! Whether or not to repartition laddie to its own mesh

    ! Remapping
    character(len=1024) :: choice_laddie_remapping_option_config        = 'no_vel'                         ! Remapping option to apply. Options: 'full', 'no_vel'

    ! Output
    logical             :: do_write_laddie_output_fields_config         = .false.                          ! Whether or not to write output fields on laddie time
    logical             :: do_write_laddie_output_scalar_config         = .false.                          ! Whether or not to write output scalars on laddie time
    real(dp)            :: time_interval_scalar_output_config           = 1._dp                            ! [days] Time interval at which to write out buffered scalars

    ! Time step
    real(dp)            :: dt_laddie_config                             = 360._dp                          ! [s] Time step for integration of laddie model
    real(dp)            :: time_duration_laddie_config                  = 6._dp                            ! [days] Duration of each run cycle
    real(dp)            :: time_duration_laddie_init_config             = 30._dp                           ! [days] Duration of initial run cycle

    ! Integration
    character(len=1024) :: choice_laddie_integration_scheme_config      = ''                               ! Choose integration scheme. Options: 'euler', 'fbrk3', 'lfra'
    real(dp)            :: laddie_fbrk3_beta1_config                    = 0.0_dp                           ! [] beta1 factor in FBRK3 integration. Must be between 0 and 1
    real(dp)            :: laddie_fbrk3_beta2_config                    = 0.0_dp                           ! [] beta2 factor in FBRK3 integration. Must be between 0 and 1
    real(dp)            :: laddie_fbrk3_beta3_config                    = 0.0_dp                           ! [] beta3 factor in FBRK3 integration. Must be between 0 and 1
    real(dp)            :: laddie_lfra_nu_config                        = 0.1_dp                           ! [] nu factor in LFRA integration. Must be between 0 and 1

    ! Momentum advection
    character(len=1024) :: choice_laddie_momentum_advection_config      = ''                               ! Choose momentum advection scheme. Options: 'none', 'upstream'

    ! Initialisation
    character(len=1024) :: choice_laddie_model_initialisation_config    = 'uniform'                        ! Initialisation method. Options: 'uniform', 'read_from_file'

    ! Initialisation 'read_from_file'
    character(len=1024) :: filename_laddie_restart_config               = ''                               ! File name containing main LADDIE variables
    real(dp)            :: timeframe_laddie_restart_config              = 1E9_dp

    ! Initialisation 'uniform'
    real(dp)            :: laddie_initial_thickness_config              = 10._dp                           ! [m] Initial value of thickness H
    real(dp)            :: laddie_initial_T_offset_config               = 0.0_dp                           ! [degC] Initial offset of T relative to ambient
    real(dp)            :: laddie_initial_S_offset_config               = -0.1_dp                          ! [PSU] Initial offset of S relative to ambient. Must be negative for stable buoyancy!

    ! Equation of state
    character(len=1024) :: choice_laddie_equation_of_state_config       = 'linear'                         ! Choose equation of state. Options: 'linear'
    real(dp)            :: uniform_laddie_eos_linear_alpha_config       = 3.733E-5_dp                      ! [K ^-1] 'linear' eos: thermal expansion coefficient
    real(dp)            :: uniform_laddie_eos_linear_beta_config        = 7.843E-4_dp                      ! [PSU ^-1] 'linear' eos: haline contraction coefficient

    ! Coriolis
    character(len=1024) :: choice_laddie_coriolis_config                = 'uniform'                        ! Choose option Coriolis parameter. Options: 'uniform', 'realistic'
    real(dp)            :: uniform_laddie_coriolis_parameter_config     = -1.37E-4_dp                      ! [s ^-1] 'linear' eos: thermal expansion coefficient

    ! Turbulent heat exchange
    character(len=1024) :: choice_laddie_gamma_config                   = 'uniform'                        ! Choose option turbulent heat exchange. Options: 'uniform', 'Jenkins1991'
    real(dp)            :: uniform_laddie_gamma_T_config                = 1.8E-4_dp                        ! [] 'uniform': gamma_T parameter. gamma_S = gamma_T/35

    ! Drag coefficients
    real(dp)            :: laddie_drag_coefficient_top_config           = 1.1E-3_dp                        ! [] Drag coefficient Cd_top in friction velocity
    real(dp)            :: laddie_drag_coefficient_mom_config           = 2.5E-3_dp                        ! [] Drag coefficient Cd_mom in momentum term

    ! Viscosity and diffusivity
    real(dp)            :: laddie_viscosity_config                      = 1.0E3_dp                         ! [m^2 s^-1] Viscosity parameter Ah
    real(dp)            :: laddie_diffusivity_config                    = 1.0E3_dp                         ! [m^2 s^-1] Diffusivity parameter Kh

    ! Entrainment
    character(len=1024) :: choice_laddie_entrainment_config             = 'Gaspar1988'                     ! Choose option entrainment parameterisation. Options: 'Holland2006', 'Gaspar1988'
    real(dp)            :: laddie_Holland2006_cl_config                 = 1.775E-2_dp                      ! [] Scaling parameter c_l in Holland2006 entrainment
    real(dp)            :: laddie_Gaspar1988_mu_config                  = 2.5_dp                           ! [] Scaling parameter mu in Gaspar1988 entrainment

    ! Stability
    real(dp)            :: laddie_thickness_minimum_config              = 2.0_dp                           ! [m] Minimum layer thickness allowed
    real(dp)            :: laddie_thickness_maximum_config              = 1500.0_dp                        ! [m] Maximum layer thickness allowed
    real(dp)            :: laddie_velocity_maximum_config               = 1.414_dp                         ! [m s^-1] Maximum velocity allowed
    real(dp)            :: laddie_buoyancy_minimum_config               = 5.0E-3_dp                        ! [kg m^-3] Minimum density difference allowed

    ! Tides
    character(len=1024) :: choice_laddie_tides_config                   = 'uniform'                        ! Choose option for tidal velocity. Options: 'uniform'
    real(dp)            :: uniform_laddie_tidal_velocity_config         = 0.1_dp                           ! [m s^-1] Uniform tidal velocity

    ! Subglacial discharge (SGD)
    character(len=1024) :: choice_laddie_SGD_config                     = 'none'                           ! Choose option for subglacial discharge. Options: 'none', 'idealised', 'read_from_file'
    character(len=1024) :: choice_laddie_SGD_idealised_config           = 'MISMIPplus_PC'                  ! Choose option for idealised SGD. Options: 'MISMIPplus_PC', 'MISMIPplus_PW', 'MISMIPplus_PE'
    real(dp)            :: laddie_SGD_flux_config                       = 50._dp                           ! [m^3 s^-1] Total subglacial discharge flux
    character(len=1024) :: filename_laddie_mask_SGD_config              = ''                               ! Gridded file containing the subglacial discharge mask
    real(dp)            :: start_time_of_applying_SGD_config            = -9E9_dp                          ! [yr] Start time of applying SGD when choice_laddie_SGD_config is 'idealised' or 'read_from_file'
    character(len=1024) :: transects_SGD_config                         = ''                               ! List of transects to use for applying SGD. Format: [file:file_path1,F=10 || file:file_path2,F=20]
    character(len=1024) :: distribute_SGD_config                        = 'single_cell'                    ! How to apply SGD from transect; single cell ('single_cell'), or distribute over 2 neighbours ('distribute_2neighbours')

  ! == Lateral mass balance
  ! =======================

    ! Time step
    real(dp)            :: dt_LMB_config                                = 10._dp                           ! [yr] Time step for calculating LMB [not implemented yet]

    ! Choice of LMB model
    character(len=1024) :: choice_LMB_model_NAM_config                  = 'uniform'
    character(len=1024) :: choice_LMB_model_EAS_config                  = 'uniform'
    character(len=1024) :: choice_LMB_model_GRL_config                  = 'uniform'
    character(len=1024) :: choice_LMB_model_ANT_config                  = 'uniform'

    ! "uniform"
    real(dp)            :: uniform_LMB_config                           = 0._dp

    ! "GlacialIndex"
    character(len=1024) :: filename_LMB_GI_NAM_config                   = ''                              ! time series file containing the glacial index record
    character(len=1024) :: filename_LMB_GI_EAS_config                   = ''                              ! time series file containing the glacial index record
    character(len=1024) :: filename_LMB_GI_GRL_config                   = ''                              ! time series file containing the glacial index record
    character(len=1024) :: filename_LMB_GI_ANT_config                   = ''                              ! time series file containing the glacial index record

    real(dp)            :: warm_LMB_NAM_config                          = -1.0_dp                         ! constant LMB value for a "warm" (GI=0.0) period
    real(dp)            :: warm_LMB_EAS_config                          = -1.0_dp                         !
    real(dp)            :: warm_LMB_GRL_config                          = -1.0_dp                         !
    real(dp)            :: warm_LMB_ANT_config                          = -1.0_dp                         !

    real(dp)            :: cold_LMB_NAM_config                          = 0.0_dp                          ! constant LMB value for a "cold" (GI=1.0) period
    real(dp)            :: cold_LMB_EAS_config                          = 0.0_dp                          !
    real(dp)            :: cold_LMB_GRL_config                          = 0.0_dp                          !
    real(dp)            :: cold_LMB_ANT_config                          = 0.0_dp                          !


  ! == Glacial isostatic adjustment
  ! ===============================

    ! General settings
    character(len=1024) :: choice_GIA_model_config                      = 'none'
    real(dp)            :: dt_GIA_config                                = 100._dp                         ! [yr] GIA model time step
    real(dp)            :: dx_GIA_config                                = 50E3_dp                         ! [m]  GIA model square grid resolution
    real(dp)            :: ELRA_lithosphere_flex_rigidity_config        = 1.0E+25                         ! [kg m^2 s^-2] Lithospheric flexural rigidity
    character(len=1024) :: choice_GIA_ELRA_relaxation_time_config       = 'uniform'                       ! 'uniform', 'read_from_file'
    real(dp)            :: ELRA_bedrock_relaxation_time_config          = 3000.0                          ! [yr] Relaxation time for bedrock adjustment
    character(len=1024) :: filename_GIA_ELRA_relaxation_time_config     = ''                              ! Path to a NetCDF file containing a laterally varying bedrock relaxation time field (variable name 'tau')
    real(dp)            :: ELRA_mantle_density_config                   = 3300.0                          ! [kg m^-3] Mantle density

  ! == Sea level
  ! ============

    character(len=1024) :: choice_sealevel_model_config                 = 'fixed'                         !     Can be "fixed", "prescribed", "eustatic", or "SELEN"
    real(dp)            :: fixed_sealevel_config                        = 0._dp                           ! [m] Fixed sea level value for the "fixed" choice
    character(len=1024) :: filename_prescribed_sealevel_config          = ''                              ! time series file containing the sea level record

  ! == SELEN
  ! ========

    logical             :: SELEN_run_at_t_start_config                  = .false.                         ! Whether or not to run SELEN in the first coupling loop (needed for some benchmark experiments)
    integer             :: SELEN_n_TDOF_iterations_config               = 1                               ! Number of Time-Dependent Ocean Function iterations
    integer             :: SELEN_n_recursion_iterations_config          = 1                               ! Number of recursion iterations
    logical             :: SELEN_use_rotational_feedback_config         = .false.                         ! If TRUE, rotational feedback is included
    integer             :: SELEN_n_harmonics_config                     = 128                             ! Maximum number of harmonic degrees
    logical             :: SELEN_display_progress_config                = .false.                         ! Whether or not to display the progress of the big loops to the screen (doesn't work on Cartesius!)

    character(len=1024) :: SELEN_dir_config                             = 'external/data/SELEN'                    ! Directory where SELEN initial files and spherical harmonics are stored
    character(len=1024) :: SELEN_global_topo_filename_config            = 'SELEN_global_topography.nc'    ! Filename for the SELEN global topography file (located in SELEN_dir)
    character(len=1024) :: SELEN_TABOO_init_filename_config             = 'SELEN_TABOO_initial_file.dat'  ! Filename for the TABOO initial file           (idem                )
    character(len=1024) :: SELEN_LMJ_VALUES_filename_config             = 'SELEN_lmj_values.bin'          ! Filename for the LJ and MJ values file        (idem                )

    integer                  :: SELEN_irreg_time_n_config               = 15                              ! Number of entries in the irregular moving time window
    real(dp), DIMENSION(50)  :: SELEN_irreg_time_window_config          = &                               ! Values of entries in the irregular moving time window
   (/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)

    real(dp)            :: SELEN_lith_thickness_config                  = 100._dp                         ! Thickness of the elastic lithosphere [km]
    integer             :: SELEN_visc_n_config                          = 3                               ! Number      of viscous asthenosphere layers
    real(dp), DIMENSION(3) :: SELEN_visc_prof_config                    = (/ 3._dp, 0.6_dp, 0.3_dp /)     ! Viscosities of viscous asthenosphere layers [?]

    ! Settings for the TABOO Earth deformation model
    integer             :: SELEN_TABOO_CDE_config                       = 0                               !     code of the model (see taboo for explanation)
    integer             :: SELEN_TABOO_TLOVE_config                     = 1                               !     Tidal love numbers yes/no
    integer             :: SELEN_TABOO_DEG1_config                      = 1                               !     Tidal love numbers degree
    real(dp)            :: SELEN_TABOO_RCMB_config                      = 3480._dp                        ! [m] Radius of CMB

  ! == Tracer tracking
  ! ==================

    character(len=1024) :: choice_tracer_tracking_model_config          = 'none'                          ! Which tracer-tracking model to use (current options: 'none', 'particles')

    ! Settings for the particle-based tracer-tracking model
    real(dp)            :: tractrackpart_dt_coupling_config             = 10._dp                          ! [yr] Coupling interval for the particle-based tracer-tracking model
    real(dp)            :: tractrackpart_dx_particle_config             = 1e3_dp                          ! [m]  Distance that tracer-tracking particles should move in a single particle time step
    real(dp)            :: tractrackpart_dt_particle_min_config         = 1._dp                           ! [yr] Minimum allowed time step for tracer-tracking particles
    real(dp)            :: tractrackpart_dt_particle_max_config         = 1000._dp                        ! [yr] Maximum allowed time step for tracer-tracking particles
    integer             :: tractrackpart_n_max_particles_config         = 1000000                         !      Maximum number of particles (i.e. for how many particles is memory allocated)
    real(dp)            :: tractrackpart_dt_new_particles_config        = 100._dp                         ! [yr] How often new batches of particles should be added
    real(dp)            :: tractrackpart_dx_new_particles_config        = 50e3_dp                         ! [m]  How far new particles should be spaced apart (square grid)
    integer             :: tractrackpart_remap_n_nearest_config         = 4                               !      Between how many nearest particles should be interpolated when mapping tracers to the mesh
    logical             :: tractrackpart_write_raw_output_config        = .false.                         !      Whether or not to write the raw particle data to a NetCDF output file
    real(dp)            :: tractrackpart_dt_raw_output_config           = 100._dp                         !      Time step for writing raw particle data to NetCDF

  ! == Output
  ! =========

    ! Basic settings
    logical             :: do_create_netcdf_output_config               = .true.                          !     Whether or not NetCDF output files should be created at all
    CHARACTER(LEN=1024) :: output_precision_config                      = 'double'                        !     Precision of floating-point output fields ('single' [32-bit], 'double' [64-bit])
    logical             :: do_compress_output_config                    = .false.                         !     Whether or not to use the NetCDF 'shuffle' and 'deflate' options to (losslessly) compress output
    real(dp)            :: dt_output_config                             = 1000._dp                        !     Time step for writing output
    real(dp)            :: dt_output_restart_config                     = 1000._dp                        !     Time step for writing restart output
    real(dp)            :: dt_output_grid_config                        = 1000._dp                        !     Time step for writing gridded output
    real(dp)            :: dx_output_grid_NAM_config                    = 40E3_dp                         ! [m] Horizontal resolution for the square grid used for output for North America
    real(dp)            :: dx_output_grid_EAS_config                    = 40E3_dp                         ! [m] Horizontal resolution for the square grid used for output for Eurasia
    real(dp)            :: dx_output_grid_GRL_config                    = 20E3_dp                         ! [m] Horizontal resolution for the square grid used for output for Greenland
    real(dp)            :: dx_output_grid_ANT_config                    = 40E3_dp                         ! [m] Horizontal resolution for the square grid used for output for Antarctica
    real(dp)            :: dx_output_grid_ROI_NAM_config                = 5E3_dp                          ! [m] Horizontal resolution for the square grid used for output for the region of interest for North America
    real(dp)            :: dx_output_grid_ROI_EAS_config                = 5E3_dp                          ! [m] Horizontal resolution for the square grid used for output for the region of interest for Eurasia
    real(dp)            :: dx_output_grid_ROI_GRL_config                = 5E3_dp                          ! [m] Horizontal resolution for the square grid used for output for the region of interest for Greenland
    real(dp)            :: dx_output_grid_ROI_ANT_config                = 5E3_dp                          ! [m] Horizontal resolution for the square grid used for output for the region of interest for Antarctica

    ! ISMIP
    logical             :: do_create_ismip_output_config                = .false.                         ! Whether or not ISMIP-specific output files should be created at all
    character(len=1024) :: ismip_scenario_name_config                   = 'ctrl'                          ! Scenario name for ISMIP output
    character(len=1024) :: ismip_group_name_config                      = 'IMAU'                          ! Group name included in ISMIP output files
    character(len=1024) :: ismip_model_name_config                      = 'UFEMISM1'                      ! Model name to which numbers can be added for variants
    character(len=4)    :: ismip_member_id_config                       = 'm001'                          ! ISM member id
    character(len=4)    :: ismip_forcing_member_id_config               = 'f001'                          ! Forcing member id
    character(len=4)    :: ismip_counter_config                         = 'C001'                          ! Counter
    character(len=1024) :: ismip_esm_name_config                        = 'CESM2-WACCM'                   ! ESM forcing
    character(len=1024) :: ismip_contact_name_config                    = ''                              ! Contact name for metadata
    character(len=1024) :: ismip_contact_email_config                   = ''                              ! Contact email for metadata
    character(len=1024) :: ismip_conventions_config                     = 'CF-1.10'                        ! CF conventions for metadata

    real(dp)            :: dt_output_ismip_config                       = 1._dp                           ! Timestep for writing ISMIP output

    ! Transects
    character(len=1024) :: transects_NAM_config                         = ''                              ! List of transects to use for North America
    character(len=1024) :: transects_EAS_config                         = ''                              ! List of transects to use for Eurasia
    character(len=1024) :: transects_GRL_config                         = ''                              ! List of transects to use for Greenland
    character(len=1024) :: transects_ANT_config                         = ''                              ! List of transects to use for Antarctica

    ! Which data fields we want to write to the main NetCDF output files
    character(len=1024) :: choice_output_field_01_config                = 'none'
    character(len=1024) :: choice_output_field_02_config                = 'none'
    character(len=1024) :: choice_output_field_03_config                = 'none'
    character(len=1024) :: choice_output_field_04_config                = 'none'
    character(len=1024) :: choice_output_field_05_config                = 'none'
    character(len=1024) :: choice_output_field_06_config                = 'none'
    character(len=1024) :: choice_output_field_07_config                = 'none'
    character(len=1024) :: choice_output_field_08_config                = 'none'
    character(len=1024) :: choice_output_field_09_config                = 'none'
    character(len=1024) :: choice_output_field_10_config                = 'none'
    character(len=1024) :: choice_output_field_11_config                = 'none'
    character(len=1024) :: choice_output_field_12_config                = 'none'
    character(len=1024) :: choice_output_field_13_config                = 'none'
    character(len=1024) :: choice_output_field_14_config                = 'none'
    character(len=1024) :: choice_output_field_15_config                = 'none'
    character(len=1024) :: choice_output_field_16_config                = 'none'
    character(len=1024) :: choice_output_field_17_config                = 'none'
    character(len=1024) :: choice_output_field_18_config                = 'none'
    character(len=1024) :: choice_output_field_19_config                = 'none'
    character(len=1024) :: choice_output_field_20_config                = 'none'
    character(len=1024) :: choice_output_field_21_config                = 'none'
    character(len=1024) :: choice_output_field_22_config                = 'none'
    character(len=1024) :: choice_output_field_23_config                = 'none'
    character(len=1024) :: choice_output_field_24_config                = 'none'
    character(len=1024) :: choice_output_field_25_config                = 'none'
    character(len=1024) :: choice_output_field_26_config                = 'none'
    character(len=1024) :: choice_output_field_27_config                = 'none'
    character(len=1024) :: choice_output_field_28_config                = 'none'
    character(len=1024) :: choice_output_field_29_config                = 'none'
    character(len=1024) :: choice_output_field_30_config                = 'none'
    character(len=1024) :: choice_output_field_31_config                = 'none'
    character(len=1024) :: choice_output_field_32_config                = 'none'
    character(len=1024) :: choice_output_field_33_config                = 'none'
    character(len=1024) :: choice_output_field_34_config                = 'none'
    character(len=1024) :: choice_output_field_35_config                = 'none'
    character(len=1024) :: choice_output_field_36_config                = 'none'
    character(len=1024) :: choice_output_field_37_config                = 'none'
    character(len=1024) :: choice_output_field_38_config                = 'none'
    character(len=1024) :: choice_output_field_39_config                = 'none'
    character(len=1024) :: choice_output_field_40_config                = 'none'
    character(len=1024) :: choice_output_field_41_config                = 'none'
    character(len=1024) :: choice_output_field_42_config                = 'none'
    character(len=1024) :: choice_output_field_43_config                = 'none'
    character(len=1024) :: choice_output_field_44_config                = 'none'
    character(len=1024) :: choice_output_field_45_config                = 'none'
    character(len=1024) :: choice_output_field_46_config                = 'none'
    character(len=1024) :: choice_output_field_47_config                = 'none'
    character(len=1024) :: choice_output_field_48_config                = 'none'
    character(len=1024) :: choice_output_field_49_config                = 'none'
    character(len=1024) :: choice_output_field_50_config                = 'none'

! ===== Configuration variables - end =====
! =========================================





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

    ! Output directory
    logical             :: create_procedural_output_dir
    character(len=1024) :: fixed_output_dir

    ! Debugging
    logical             :: do_unit_tests
    logical             :: do_benchmarks
    logical             :: do_check_for_NaN
    logical             :: do_time_display
    logical             :: do_write_matrix_operators
    logical             :: do_write_checksum_log

  ! == Time of simulation
  ! =====================

    real(dp)            :: start_time_of_run
    real(dp)            :: end_time_of_run

  ! == Which model regions to simulate
  ! ==================================

    real(dp)            :: dt_coupling
    logical             :: do_NAM
    logical             :: do_EAS
    logical             :: do_GRL
    logical             :: do_ANT

  ! == The four model regions
  ! =========================

    ! North America
    real(dp)            :: lambda_M_NAM
    real(dp)            :: phi_M_NAM
    real(dp)            :: beta_stereo_NAM
    real(dp)            :: xmin_NAM
    real(dp)            :: xmax_NAM
    real(dp)            :: ymin_NAM
    real(dp)            :: ymax_NAM

    ! Eurasia
    real(dp)            :: lambda_M_EAS
    real(dp)            :: phi_M_EAS
    real(dp)            :: beta_stereo_EAS
    real(dp)            :: xmin_EAS
    real(dp)            :: xmax_EAS
    real(dp)            :: ymin_EAS
    real(dp)            :: ymax_EAS

    ! Greenland
    real(dp)            :: lambda_M_GRL
    real(dp)            :: phi_M_GRL
    real(dp)            :: beta_stereo_GRL
    real(dp)            :: xmin_GRL
    real(dp)            :: xmax_GRL
    real(dp)            :: ymin_GRL
    real(dp)            :: ymax_GRL

    ! Antarctica
    real(dp)            :: lambda_M_ANT
    real(dp)            :: phi_M_ANT
    real(dp)            :: beta_stereo_ANT
    real(dp)            :: xmin_ANT
    real(dp)            :: xmax_ANT
    real(dp)            :: ymin_ANT
    real(dp)            :: ymax_ANT

  ! == Reference geometries (initial, present-day, and GIA equilibrium)
  ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    real(dp)            :: refgeo_Hi_min
    logical             :: do_smooth_geometry
    real(dp)            :: r_smooth_geometry
    logical             :: remove_Lake_Vostok

    ! Remapping method for gridded geometry input
    character(len=1024) :: choice_refgeo_remapping_method

    ! == Initial geometry
    ! ===================

    character(len=1024) :: choice_refgeo_init_NAM
    character(len=1024) :: choice_refgeo_init_EAS
    character(len=1024) :: choice_refgeo_init_GRL
    character(len=1024) :: choice_refgeo_init_ANT
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    character(len=1024) :: choice_refgeo_init_idealised
    real(dp)            :: dx_refgeo_init_idealised
    ! Path to file containing initial geometry when choice_refgeo_init == 'read_from_file'
    character(len=1024) :: filename_refgeo_init_NAM
    character(len=1024) :: filename_refgeo_init_EAS
    character(len=1024) :: filename_refgeo_init_GRL
    character(len=1024) :: filename_refgeo_init_ANT
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_refgeo_init_NAM
    real(dp)            :: timeframe_refgeo_init_EAS
    real(dp)            :: timeframe_refgeo_init_GRL
    real(dp)            :: timeframe_refgeo_init_ANT

    ! == Present-day geometry
    ! =======================

    character(len=1024) :: choice_refgeo_PD_NAM
    character(len=1024) :: choice_refgeo_PD_EAS
    character(len=1024) :: choice_refgeo_PD_GRL
    character(len=1024) :: choice_refgeo_PD_ANT
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    character(len=1024) :: choice_refgeo_PD_idealised
    real(dp)            :: dx_refgeo_PD_idealised
    ! Path to file containing present-day geometry when choice_refgeo_PD == 'read_from_file'
    character(len=1024) :: filename_refgeo_PD_NAM
    character(len=1024) :: filename_refgeo_PD_EAS
    character(len=1024) :: filename_refgeo_PD_GRL
    character(len=1024) :: filename_refgeo_PD_ANT
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_refgeo_PD_NAM
    real(dp)            :: timeframe_refgeo_PD_EAS
    real(dp)            :: timeframe_refgeo_PD_GRL
    real(dp)            :: timeframe_refgeo_PD_ANT

    ! == GIA equilibrium geometry
    ! ===========================

    character(len=1024) :: choice_refgeo_GIAeq_NAM
    character(len=1024) :: choice_refgeo_GIAeq_EAS
    character(len=1024) :: choice_refgeo_GIAeq_GRL
    character(len=1024) :: choice_refgeo_GIAeq_ANT
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    character(len=1024) :: choice_refgeo_GIAeq_idealised
    real(dp)            :: dx_refgeo_GIAeq_idealised
    ! Path to file containing GIA equilibrium reference geometry when choice_refgeo_GIAeq == 'read_from_file'
    character(len=1024) :: filename_refgeo_GIAeq_NAM
    character(len=1024) :: filename_refgeo_GIAeq_EAS
    character(len=1024) :: filename_refgeo_GIAeq_GRL
    character(len=1024) :: filename_refgeo_GIAeq_ANT
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_refgeo_GIAeq_NAM
    real(dp)            :: timeframe_refgeo_GIAeq_EAS
    real(dp)            :: timeframe_refgeo_GIAeq_GRL
    real(dp)            :: timeframe_refgeo_GIAeq_ANT

    ! == Parameters for idealised geometries
    ! ======================================

    real(dp)            :: refgeo_idealised_slabonaslope_Hi
    real(dp)            :: refgeo_idealised_slabonaslope_dhdx
    real(dp)            :: refgeo_idealised_Halfar_H0
    real(dp)            :: refgeo_idealised_Halfar_R0
    real(dp)            :: refgeo_idealised_Bueler_H0
    real(dp)            :: refgeo_idealised_Bueler_R0
    real(dp)            :: refgeo_idealised_Bueler_lambda
    real(dp)            :: refgeo_idealised_SSA_icestream_Hi
    real(dp)            :: refgeo_idealised_SSA_icestream_dhdx
    real(dp)            :: refgeo_idealised_SSA_icestream_L
    real(dp)            :: refgeo_idealised_SSA_icestream_m
    real(dp)            :: refgeo_idealised_MISMIP_mod_Hi_init
    real(dp)            :: refgeo_idealised_ISMIP_HOM_L
    real(dp)            :: refgeo_idealised_MISMIPplus_Hi_init
    logical             :: refgeo_idealised_MISMIPplus_tune_A

  ! == Mesh generation
  ! ==================

    ! How to set up the initial mesh
    character(len=1024) :: choice_initial_mesh_NAM
    character(len=1024) :: choice_initial_mesh_EAS
    character(len=1024) :: choice_initial_mesh_GRL
    character(len=1024) :: choice_initial_mesh_ANT

    ! Paths to files containing initial meshes, if choice_initial_mesh == 'read_from_file'
    character(len=1024) :: filename_initial_mesh_NAM
    character(len=1024) :: filename_initial_mesh_EAS
    character(len=1024) :: filename_initial_mesh_GRL
    character(len=1024) :: filename_initial_mesh_ANT

    ! Resolutions for different parts of the ice sheet
    real(dp)            :: maximum_resolution_uniform
    real(dp)            :: maximum_resolution_grounded_ice
    real(dp)            :: maximum_resolution_floating_ice
    real(dp)            :: maximum_resolution_grounding_line
    real(dp)            :: grounding_line_width
    real(dp)            :: maximum_resolution_calving_front
    real(dp)            :: calving_front_width
    real(dp)            :: maximum_resolution_ice_front
    real(dp)            :: ice_front_width
    real(dp)            :: maximum_resolution_coastline
    real(dp)            :: coastline_width

    ! Regions of interest
    character(len=1024) :: choice_regions_of_interest
    real(dp)            :: ROI_maximum_resolution_uniform
    real(dp)            :: ROI_maximum_resolution_grounded_ice
    real(dp)            :: ROI_maximum_resolution_floating_ice
    real(dp)            :: ROI_maximum_resolution_grounding_line
    real(dp)            :: ROI_grounding_line_width
    real(dp)            :: ROI_maximum_resolution_calving_front
    real(dp)            :: ROI_calving_front_width
    real(dp)            :: ROI_maximum_resolution_ice_front
    real(dp)            :: ROI_ice_front_width
    real(dp)            :: ROI_maximum_resolution_coastline
    real(dp)            :: ROI_coastline_width

    ! Miscellaneous refinement options
    logical             :: do_refine_TransAntMounts_glaciers
    real(dp)            :: max_res_TransAntMounts_glaciers

    ! Mesh update settings
    logical             :: allow_mesh_updates
    real(dp)            :: dt_mesh_update_min
    real(dp)            :: minimum_mesh_fitness_coefficient
    logical             :: do_out_of_time_calving_front_relax

    ! Advanced geometry parameters
    logical             :: do_singlecore_mesh_creation
    real(dp)            :: alpha_min
    integer             :: nit_Lloyds_algorithm
    real(dp)            :: mesh_resolution_tolerance

    ! Square grid used for smoothing
    real(dp)            :: dx_square_grid_smooth_NAM
    real(dp)            :: dx_square_grid_smooth_EAS
    real(dp)            :: dx_square_grid_smooth_GRL
    real(dp)            :: dx_square_grid_smooth_ANT

  ! == The scaled vertical coordinate zeta
  ! ======================================

    character(len=1024) :: choice_zeta_grid
    integer             :: nz
    real(dp)            :: zeta_irregular_log_R

  ! == Ice dynamics - velocity
  ! ==========================

    ! General
    character(len=1024) :: choice_stress_balance_approximation
    character(len=1024) :: choice_hybrid_SIASSA_scheme
    logical             :: do_include_SSADIVA_crossterms

    ! Hybrid DIVA/BPA
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_NAM
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_EAS
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_GRL
    character(len=1024) :: choice_hybrid_DIVA_BPA_mask_ANT

    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_NAM
    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_EAS
    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_GRL
    character(len=1024) :: filename_hybrid_DIVA_BPA_mask_ANT

    ! Initialisation
    character(len=1024) :: choice_initial_velocity_NAM
    character(len=1024) :: choice_initial_velocity_EAS
    character(len=1024) :: choice_initial_velocity_GRL
    character(len=1024) :: choice_initial_velocity_ANT
    ! Paths to files containing initial velocity fields
    character(len=1024) :: filename_initial_velocity_NAM
    character(len=1024) :: filename_initial_velocity_EAS
    character(len=1024) :: filename_initial_velocity_GRL
    character(len=1024) :: filename_initial_velocity_ANT
    ! Timeframes to read from the initial velocity files (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_initial_velocity_NAM
    real(dp)            :: timeframe_initial_velocity_EAS
    real(dp)            :: timeframe_initial_velocity_GRL
    real(dp)            :: timeframe_initial_velocity_ANT

    ! Some parameters for numerically solving the stress balance
    real(dp)            :: SIA_maximum_diffusivity
    real(dp)            :: visc_it_norm_dUV_tol
    integer             :: visc_it_nit
    real(dp)            :: visc_it_relax
    real(dp)            :: visc_eff_min
    real(dp)            :: vel_max
    character(len=1024) :: stress_balance_PETSc_KSPtype
    character(len=1024) :: stress_balance_PETSc_PCtype
    real(dp)            :: stress_balance_PETSc_rtol
    real(dp)            :: stress_balance_PETSc_abstol

    ! Boundary conditions
    character(len=1024) :: BC_ice_front
    character(len=1024) :: BC_u_west
    character(len=1024) :: BC_u_east
    character(len=1024) :: BC_u_south
    character(len=1024) :: BC_u_north
    character(len=1024) :: BC_v_west
    character(len=1024) :: BC_v_east
    character(len=1024) :: BC_v_south
    character(len=1024) :: BC_v_north

  ! == Ice dynamics - sliding
  ! =========================

    ! General
    character(len=1024) :: choice_sliding_law
    character(len=1024) :: choice_idealised_sliding_law

    ! Parameters for different sliding laws
    real(dp)            :: slid_Weertman_m
    real(dp)            :: slid_Budd_q_plastic
    real(dp)            :: slid_Budd_u_threshold
    real(dp)            :: slid_ZI_p
    real(dp)            :: slid_ZI_ut

    ! Sub-grid scaling of basal friction
    logical             :: do_GL_subgrid_friction
    character(len=1024) :: choice_subgrid_grounded_fraction
    logical             :: do_read_bedrock_cdf_from_file
    integer             :: subgrid_bedrock_cdf_nbins
    logical             :: do_subgrid_friction_on_A_grid
    real(dp)            :: subgrid_friction_exponent_on_B_grid

    ! Stability
    real(dp)            :: slid_beta_max
    real(dp)            :: slid_delta_v

  ! == Ice dynamics - ice thickness calculation
  ! ===========================================

    ! Calculation of dH/dt
    character(len=1024) :: choice_ice_integration_method
    real(dp)            :: dHi_semiimplicit_fs
    character(len=1024) :: dHi_PETSc_KSPtype
    character(len=1024) :: dHi_PETSc_PCtype
    real(dp)            :: dHi_PETSc_rtol
    real(dp)            :: dHi_PETSc_abstol

    ! Boundary conditions
    character(len=1024) :: BC_H_west
    character(len=1024) :: BC_H_east
    character(len=1024) :: BC_H_south
    character(len=1024) :: BC_H_north

  ! == Ice dynamics - target quantities
  ! ===================================

    ! Target dHi_dt
    logical             :: do_target_dHi_dt
    logical             :: do_limit_target_dHi_dt_to_SMB
    real(dp)            :: target_dHi_dt_t_end

    ! Files containing a target dHi_dt for inversions
    character(len=1024) :: filename_dHi_dt_target_NAM
    character(len=1024) :: filename_dHi_dt_target_EAS
    character(len=1024) :: filename_dHi_dt_target_GRL
    character(len=1024) :: filename_dHi_dt_target_ANT

    ! Timeframes for reading target dHi_dt from file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_dHi_dt_target_NAM
    real(dp)            :: timeframe_dHi_dt_target_EAS
    real(dp)            :: timeframe_dHi_dt_target_GRL
    real(dp)            :: timeframe_dHi_dt_target_ANT

  ! == Ice dynamics - time stepping
  ! ===============================

    ! Time stepping
    character(len=1024) :: choice_timestepping
    real(dp)            :: dt_ice_max
    real(dp)            :: dt_ice_min
    real(dp)            :: dt_ice_startup_phase
    logical             :: do_grounded_only_adv_dt

    ! Predictor-corrector ice-thickness update
    real(dp)            :: pc_epsilon
    real(dp)            :: pc_k_I
    real(dp)            :: pc_k_p
    real(dp)            :: pc_eta_min
    real(dp)            :: pc_max_time_step_increase
    integer             :: pc_nit_max
    real(dp)            :: pc_guilty_max

    ! Initialisation of the predictor-corrector ice-thickness update
    character(len=1024) :: pc_choice_initialise_NAM
    character(len=1024) :: pc_choice_initialise_EAS
    character(len=1024) :: pc_choice_initialise_GRL
    character(len=1024) :: pc_choice_initialise_ANT
    ! Paths to files containing initial fields & values for the p/c scheme
    character(len=1024) :: filename_pc_initialise_NAM
    character(len=1024) :: filename_pc_initialise_EAS
    character(len=1024) :: filename_pc_initialise_GRL
    character(len=1024) :: filename_pc_initialise_ANT
    ! Timeframes to read from the p/c scheme initial file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_pc_initialise_NAM
    real(dp)            :: timeframe_pc_initialise_EAS
    real(dp)            :: timeframe_pc_initialise_GRL
    real(dp)            :: timeframe_pc_initialise_ANT

  ! == Ice dynamics - calving
  ! =========================

    character(len=1024) :: choice_calving_law
    real(dp)            :: calving_threshold_thickness_shelf
    real(dp)            :: calving_threshold_thickness_sheet
    integer             :: max_calving_rounds
    logical             :: do_remove_shelves
    logical             :: remove_shelves_larger_than_PD
    logical             :: continental_shelf_calving
    real(dp)            :: continental_shelf_min_height

    logical             :: do_apply_ISMIP7_fracture_mask
    character(len=1024) :: filename_ISMIP7_fracture_mask

  ! == Ice dynamics - stabilisation
  ! ===============================

    character(len=1024) :: choice_mask_noice
    real(dp)            :: Hi_min
    real(dp)            :: Hi_thin
    logical             :: remove_ice_absent_at_PD

    ! Geometry relaxation
    real(dp)            :: geometry_relaxation_t_years

    ! Mask conservation
    logical             :: do_protect_grounded_mask
    real(dp)            :: protect_grounded_mask_t_end


    ! Fix/delay ice thickness evolution
    logical             :: do_fixiness_before_start
    real(dp)            :: fixiness_t_start
    real(dp)            :: fixiness_t_end
    real(dp)            :: fixiness_H_gl_gr
    real(dp)            :: fixiness_H_gl_fl
    real(dp)            :: fixiness_H_grounded
    real(dp)            :: fixiness_H_floating
    logical             :: fixiness_H_freeland
    logical             :: fixiness_H_freeocean

    ! Limit ice thickness evolution
    logical             :: do_limitness_before_start
    real(dp)            :: limitness_t_start
    real(dp)            :: limitness_t_end
    real(dp)            :: limitness_H_gl_gr
    real(dp)            :: limitness_H_gl_fl
    real(dp)            :: limitness_H_grounded
    real(dp)            :: limitness_H_floating
    character(len=1024) :: modiness_H_style
    real(dp)            :: modiness_T_hom_ref

  ! == Basal hydrology
  ! ==================

    ! Time step
    logical             :: do_asynchronous_basal_hydro
    real(dp)            :: dt_basal_hydro

    ! Basal hydrology
    character(len=1024) :: choice_basal_hydrology_model
    real(dp)            :: Martin2011_hydro_Hb_min
    real(dp)            :: Martin2011_hydro_Hb_max
    real(dp)            :: basal_hydro_equil_time
    real(dp)            :: error_function_max_effective_pressure
    real(dp)            :: Leguy2014_hydro_connect_exponent
    real(dp)            :: Salle2025_hydro_meltwater
    real(dp)            :: Salle2025_W_til_max
    real(dp)            :: Salle2025_initial_W_til
    real(dp)            :: Salle2025_initial_W
    real(dp)            :: Salle2025_Cd
    real(dp)            :: Salle2025_phi

    character(len=1024) :: choice_BMB_grounded

  ! == Bed roughness
  ! ==================

    character(len=1024) :: choice_bed_roughness
    character(len=1024) :: choice_bed_roughness_parameterised
    ! Paths to files containing bed roughness fields for the chosen sliding law
    character(len=1024) :: filename_bed_roughness_NAM
    character(len=1024) :: filename_bed_roughness_EAS
    character(len=1024) :: filename_bed_roughness_GRL
    character(len=1024) :: filename_bed_roughness_ANT
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_bed_roughness_NAM
    real(dp)            :: timeframe_bed_roughness_EAS
    real(dp)            :: timeframe_bed_roughness_GRL
    real(dp)            :: timeframe_bed_roughness_ANT
    ! Values for uniform bed roughness
    real(dp)            :: slid_Weertman_beta_sq_uniform
    real(dp)            :: slid_Coulomb_phi_fric_uniform
    real(dp)            :: slid_Budd_phi_fric_uniform
    real(dp)            :: slid_Tsai2015_alpha_sq_uniform
    real(dp)            :: slid_Tsai2015_beta_sq_uniform
    real(dp)            :: slid_Schoof2005_alpha_sq_uniform
    real(dp)            :: slid_Schoof2005_beta_sq_uniform
    real(dp)            :: slid_ZI_phi_fric_uniform
    ! Parameters for bed roughness parameterisations
    real(dp)            :: Martin2011till_phi_Hb_min
    real(dp)            :: Martin2011till_phi_Hb_max
    real(dp)            :: Martin2011till_phi_min
    real(dp)            :: Martin2011till_phi_max

  ! == Bed roughness inversion by nudging
  ! =====================================

    ! General
    logical             :: do_bed_roughness_nudging
    character(len=1024) :: choice_bed_roughness_nudging_method
    character(len=1024) :: choice_inversion_target_geometry
    real(dp)            :: bed_roughness_nudging_dt
    real(dp)            :: bed_roughness_nudging_t_start
    real(dp)            :: bed_roughness_nudging_t_end
    real(dp)            :: generic_bed_roughness_min
    real(dp)            :: generic_bed_roughness_max

    ! Bed roughness nudging model based on flowline-averaged values of H and dH/dt
    real(dp)            :: bednudge_H_dHdt_flowline_t_scale
    real(dp)            :: bednudge_H_dHdt_flowline_dH0
    real(dp)            :: bednudge_H_dHdt_flowline_dHdt0
    real(dp)            :: bednudge_H_dHdt_flowline_Hi_scale
    real(dp)            :: bednudge_H_dHdt_flowline_u_scale
    real(dp)            :: bednudge_H_dHdt_flowline_r_smooth
    real(dp)            :: bednudge_H_dHdt_flowline_w_smooth

    ! Bed roughness nudging model based on local values of H and dH/dt (i.e. CISM method)
    real(dp)            :: bednudge_H_dHdt_local_H0
    real(dp)            :: bednudge_H_dHdt_local_tau
    real(dp)            :: bednudge_H_dHdt_local_L

    ! Bed roughness nudging model based on flowline-averaged values of H and u
    CHARACTER(LEN=1024) :: bednudge_H_u_flowline_file_u_target
    real(dp)            :: bednudge_H_u_flowline_t_scale
    real(dp)            :: bednudge_H_u_flowline_H0
    real(dp)            :: bednudge_H_u_flowline_u0
    real(dp)            :: bednudge_H_u_flowline_Hi_scale
    real(dp)            :: bednudge_H_u_flowline_u_scale
    real(dp)            :: bednudge_H_u_flowline_tau
    real(dp)            :: bednudge_H_u_flowline_L

  ! == Geothermal heat flux
  ! =======================

    character(len=1024) :: choice_geothermal_heat_flux
    real(dp)            :: uniform_geothermal_heat_flux
    character(len=1024) :: filename_geothermal_heat_flux

  ! == Thermodynamics
  ! =================

    ! Initial temperature profile
    character(len=1024) :: choice_initial_ice_temperature_NAM
    character(len=1024) :: choice_initial_ice_temperature_EAS
    character(len=1024) :: choice_initial_ice_temperature_GRL
    character(len=1024) :: choice_initial_ice_temperature_ANT
    ! Uniform initial ice temperature, if choice_initial_ice_temperature == 'uniform'
    real(dp)            :: uniform_initial_ice_temperature_NAM
    real(dp)            :: uniform_initial_ice_temperature_EAS
    real(dp)            :: uniform_initial_ice_temperature_GRL
    real(dp)            :: uniform_initial_ice_temperature_ANT
    ! Paths to files containing initial temperature fields, if
    character(len=1024) :: filename_initial_ice_temperature_NAM
    character(len=1024) :: filename_initial_ice_temperature_EAS
    character(len=1024) :: filename_initial_ice_temperature_GRL
    character(len=1024) :: filename_initial_ice_temperature_ANT
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_initial_ice_temperature_NAM
    real(dp)            :: timeframe_initial_ice_temperature_EAS
    real(dp)            :: timeframe_initial_ice_temperature_GRL
    real(dp)            :: timeframe_initial_ice_temperature_ANT
    ! Thermodynamical model
    character(len=1024) :: choice_thermo_model
    real(dp)            :: dt_thermodynamics
    real(dp)            :: Hi_min_thermo
    character(len=1024) :: choice_GL_temperature_BC
    character(len=1024) :: choice_ice_heat_capacity
    real(dp)            :: uniform_ice_heat_capacity
    character(len=1024) :: choice_ice_thermal_conductivity
    real(dp)            :: uniform_ice_thermal_conductivity

  ! == Rheology and flow law
  ! =========================

    ! Flow law
    character(len=1024) :: choice_flow_law
    real(dp)            :: Glens_flow_law_exponent
    real(dp)            :: Glens_flow_law_epsilon_sq_0

    ! Rheology
    character(len=1024) :: choice_ice_rheology_Glen
    real(dp)            :: uniform_Glens_flow_factor

    ! Enhancement factors
    character(len=1024) :: choice_enhancement_factor_transition
    real(dp)            :: m_enh_sheet
    real(dp)            :: m_enh_shelf

  ! == Climate
  ! ==========

    ! Time step
    logical             :: do_asynchronous_climate
    real(dp)            :: dt_climate

    ! Choice of climate model
    character(len=1024) :: choice_climate_model_NAM
    character(len=1024) :: choice_climate_model_EAS
    character(len=1024) :: choice_climate_model_GRL
    character(len=1024) :: choice_climate_model_ANT

    ! Choice of idealised climate model
    character(len=1024) :: choice_climate_model_idealised

    ! Choice of realistic climate model
    character(len=1024) :: choice_climate_model_realistic

    ! Paths to files containing fields for realistic climates
    character(len=1024) :: filename_climate_snapshot_NAM
    character(len=1024) :: filename_climate_snapshot_EAS
    character(len=1024) :: filename_climate_snapshot_GRL
    character(len=1024) :: filename_climate_snapshot_ANT

    ! Climate - global forcings
    character(len=1024) :: filename_CO2_record
    character(len=1024) :: filename_d18O_record

    ! == Climate - fixed regional lapse rates
    logical             :: do_lapse_rate_corrections_NAM
    logical             :: do_lapse_rate_corrections_EAS
    logical             :: do_lapse_rate_corrections_GRL
    logical             :: do_lapse_rate_corrections_ANT
    real(dp)            :: lapse_rate_temp_NAM
    real(dp)            :: lapse_rate_temp_EAS
    real(dp)            :: lapse_rate_temp_GRL
    real(dp)            :: lapse_rate_temp_ANT
    real(dp)            :: lapse_rate_precip_NAM
    real(dp)            :: lapse_rate_precip_EAS
    real(dp)            :: lapse_rate_precip_GRL
    real(dp)            :: lapse_rate_precip_ANT

    ! == Climate - snapshots plus uniform deltaT
    character(len=1024) :: filename_climate_snapshot_unif_dT_NAM
    character(len=1024) :: filename_climate_snapshot_unif_dT_EAS
    character(len=1024) :: filename_climate_snapshot_unif_dT_GRL
    character(len=1024) :: filename_climate_snapshot_unif_dT_ANT
    real(dp)            :: uniform_deltaT_NAM
    real(dp)            :: uniform_deltaT_EAS
    real(dp)            :: uniform_deltaT_GRL
    real(dp)            :: uniform_deltaT_ANT
    real(dp)            :: precip_CC_correction_NAM
    real(dp)            :: precip_CC_correction_EAS
    real(dp)            :: precip_CC_correction_GRL
    real(dp)            :: precip_CC_correction_ANT

    ! == Climate - snapshot plus a transient deltaT
    character(len=1024) :: filename_climate_snapshot_trans_dT_NAM
    character(len=1024) :: filename_climate_snapshot_trans_dT_EAS
    character(len=1024) :: filename_climate_snapshot_trans_dT_GRL
    character(len=1024) :: filename_climate_snapshot_trans_dT_ANT
    character(len=1024) :: filename_atmosphere_dT_NAM
    character(len=1024) :: filename_atmosphere_dT_EAS
    character(len=1024) :: filename_atmosphere_dT_GRL
    character(len=1024) :: filename_atmosphere_dT_ANT

    ! == Climate snapshot_plus_anomalies
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_NAM
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_EAS
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_GRL
    character(len=1024) :: climate_snp_p_anml_filename_snapshot_ANT
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_NAM
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_EAS
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_GRL
    character(len=1024) :: climate_snp_p_anml_filename_anomalies_ANT

    ! == Climate - Insolation
    character(len=1024) :: choice_insolation_forcing
    character(len=1024) :: filename_insolation
    real(dp)            :: static_insolation_time

    ! == Climate matrix
    character(len=1024) :: choice_matrix_forcing
    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    character(len=1024) :: climate_matrix_filename_PD_obs_climate

    ! GCM snapshots in the matrix_warm_cold option
    character(len=1024) :: climate_matrix_filename_climate_snapshot_PI
    character(len=1024) :: climate_matrix_filename_climate_snapshot_warm
    character(len=1024) :: climate_matrix_filename_climate_snapshot_cold

    real(dp)            :: climate_matrix_constant_lapserate

    ! Scaling factor for CO2 vs ice weights
    real(dp)            :: climate_matrix_CO2vsice_NAM
    real(dp)            :: climate_matrix_CO2vsice_EAS
    real(dp)            :: climate_matrix_CO2vsice_GRL
    real(dp)            :: climate_matrix_CO2vsice_ANT

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    real(dp)            :: climate_matrix_high_CO2_level
    real(dp)            :: climate_matrix_low_CO2_level
    real(dp)            :: climate_matrix_warm_orbit_time
    real(dp)            :: climate_matrix_cold_orbit_time

    ! Whether or not to apply a bias correction to the GCM snapshots
    logical             :: climate_matrix_biascorrect_warm
    logical             :: climate_matrix_biascorrect_cold

    logical             :: climate_matrix_switch_glacial_index_precip

    ! Settings for the ISMIP7 climate model
    character(len=1024) :: climate_ISMIP7_choice_baseline
    character(len=1024) :: climate_ISMIP7_filename_baseline
    character(len=1024) :: climate_ISMIP7_choice_refgeo
    character(len=1024) :: climate_ISMIP7_forcing_foldername
    character(len=1024) :: climate_ISMIP7_forcing_version

  ! == Ocean
  ! ========

    ! Time step
    logical             :: do_asynchronous_ocean
    real(dp)            :: dt_ocean

    ! Vertical grid
    real(dp)            :: ocean_vertical_grid_max_depth
    real(dp)            :: ocean_vertical_grid_dz

    ! Choice of ocean model
    character(len=1024) :: choice_ocean_model_NAM
    character(len=1024) :: choice_ocean_model_EAS
    character(len=1024) :: choice_ocean_model_GRL
    character(len=1024) :: choice_ocean_model_ANT

    ! Choice of idealised ocean model
    character(len=1024) :: choice_ocean_model_idealised
    character(len=1024) :: choice_ocean_isomip_scenario
    real(dp)            :: ocean_tanh_deep_temperature
    real(dp)            :: ocean_tanh_thermocline_depth
    real(dp)            :: ocean_tanh_thermocline_scale_depth
    real(dp)            :: ocean_linear_deep_temperature
    real(dp)            :: ocean_linear_deep_salinity
    real(dp)            :: ocean_linear_reference_depth
    real(dp)            :: ocean_lin_therm_surf_salinity
    real(dp)            :: ocean_lin_therm_deep_salinity
    real(dp)            :: ocean_lin_therm_surf_temperature
    real(dp)            :: ocean_lin_therm_deep_temperature
    real(dp)            :: ocean_lin_therm_thermocline_top
    real(dp)            :: ocean_lin_therm_thermocline_bottom
    real(dp)            :: ocean_uniform_T
    real(dp)            :: ocean_uniform_S

    ! Choice of realistic ocean model
    character(len=1024) :: choice_ocean_model_realistic

    ! Paths to files containing fields for realistic ocean
    character(len=1024) :: filename_ocean_snapshot_NAM
    character(len=1024) :: filename_ocean_snapshot_EAS
    character(len=1024) :: filename_ocean_snapshot_GRL
    character(len=1024) :: filename_ocean_snapshot_ANT

    character(len=1024) :: filename_ocean_warm_snapshot_NAM
    character(len=1024) :: filename_ocean_warm_snapshot_EAS
    character(len=1024) :: filename_ocean_warm_snapshot_GRL
    character(len=1024) :: filename_ocean_warm_snapshot_ANT

    character(len=1024) :: filename_ocean_cold_snapshot_NAM
    character(len=1024) :: filename_ocean_cold_snapshot_EAS
    character(len=1024) :: filename_ocean_cold_snapshot_GRL
    character(len=1024) :: filename_ocean_cold_snapshot_ANT

    ! Choice of ocean deltaT
    real(dp)            :: ocean_uniform_deltaT_NAM
    real(dp)            :: ocean_uniform_deltaT_EAS
    real(dp)            :: ocean_uniform_deltaT_GRL
    real(dp)            :: ocean_uniform_deltaT_ANT

    ! Choice of extrapolation method
    character(len=1024) :: choice_ocean_extrapolation_method

  ! Choice of transient ocean model
    character(len=1024) :: choice_ocean_model_transient

    ! Paths to files containing the deltaT record for the transient deltaT ocean model
    character(len=1024) :: filename_ocean_dT_NAM
    character(len=1024) :: filename_ocean_dT_EAS
    character(len=1024) :: filename_ocean_dT_GRL
    character(len=1024) :: filename_ocean_dT_ANT

    ! Paths to files containing the GI record for the transient GlacialIndex ocean model
    character(len=1024) :: filename_ocean_GI_NAM
    character(len=1024) :: filename_ocean_GI_EAS
    character(len=1024) :: filename_ocean_GI_GRL
    character(len=1024) :: filename_ocean_GI_ANT

    ! Settings for the snapshot_plus_anomalies ocean model
    character(len=1024) :: ocean_snp_p_anml_filename_snapshot
    character(len=1024) :: ocean_snp_p_anml_filename_anomalies

    ! Settings for the ISMIP7 ocean model
    character(len=1024) :: ocean_ISMIP7_forcing_foldername
    character(len=1024) :: ocean_ISMIP7_forcing_version

  ! == Surface mass balance
  ! =======================

    ! Time step
    logical             :: do_asynchronous_SMB
    real(dp)            :: dt_SMB

    ! Choice of SMB model
    character(len=1024) :: choice_SMB_model_NAM
    character(len=1024) :: choice_SMB_model_EAS
    character(len=1024) :: choice_SMB_model_GRL
    character(len=1024) :: choice_SMB_model_ANT

    ! Value to be used for uniform SMB (no regional variants, only used for idealised-geometry experiments)
    real(dp)            :: uniform_SMB

    ! Choice of idealised SMB model
    character(len=1024) :: choice_SMB_model_idealised

    ! Prescribed SMB forcing
    character(len=1024) :: choice_SMB_prescribed_NAM
    character(len=1024) :: choice_SMB_prescribed_EAS
    character(len=1024) :: choice_SMB_prescribed_GRL
    character(len=1024) :: choice_SMB_prescribed_ANT

    ! Files containing prescribed SMB forcing
    character(len=1024) :: filename_SMB_prescribed_NAM
    character(len=1024) :: filename_SMB_prescribed_EAS
    character(len=1024) :: filename_SMB_prescribed_GRL
    character(len=1024) :: filename_SMB_prescribed_ANT

    ! Timeframes for reading prescribed SMB forcing from file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_SMB_prescribed_NAM
    real(dp)            :: timeframe_SMB_prescribed_EAS
    real(dp)            :: timeframe_SMB_prescribed_GRL
    real(dp)            :: timeframe_SMB_prescribed_ANT

    ! IMAU-ITM SMB model
    character(len=1024) :: choice_SMB_IMAUITM_init_firn_NAM
    character(len=1024) :: choice_SMB_IMAUITM_init_firn_EAS
    character(len=1024) :: choice_SMB_IMAUITM_init_firn_GRL
    character(len=1024) :: choice_SMB_IMAUITM_init_firn_ANT

    ! Files containing the firn model (yearly firn depth and melt)
    character(len=1024) :: filename_firn_IMAUITM_NAM
    character(len=1024) :: filename_firn_IMAUITM_EAS
    character(len=1024) :: filename_firn_IMAUITM_GRL
    character(len=1024) :: filename_firn_IMAUITM_ANT

    ! timeframe for restarting from the firn model
    real(dp)            :: timeframe_restart_firn_IMAUITM_NAM
    real(dp)            :: timeframe_restart_firn_IMAUITM_EAS
    real(dp)            :: timeframe_restart_firn_IMAUITM_GRL
    real(dp)            :: timeframe_restart_firn_IMAUITM_ANT

    ! Tuning parameters for the IMAU-ITM SMB model
    real(dp)            :: SMB_IMAUITM_initial_firn_thickness
    real(dp)            :: SMB_IMAUITM_C_abl_constant_NAM
    real(dp)            :: SMB_IMAUITM_C_abl_constant_EAS
    real(dp)            :: SMB_IMAUITM_C_abl_constant_GRL
    real(dp)            :: SMB_IMAUITM_C_abl_constant_ANT
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_NAM
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_EAS
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_GRL
    real(dp)            :: SMB_IMAUITM_C_abl_Ts_ANT
    real(dp)            :: SMB_IMAUITM_C_abl_Q_NAM
    real(dp)            :: SMB_IMAUITM_C_abl_Q_EAS
    real(dp)            :: SMB_IMAUITM_C_abl_Q_GRL
    real(dp)            :: SMB_IMAUITM_C_abl_Q_ANT
    real(dp)            :: SMB_IMAUITM_C_refr_NAM
    real(dp)            :: SMB_IMAUITM_C_refr_EAS
    real(dp)            :: SMB_IMAUITM_C_refr_GRL
    real(dp)            :: SMB_IMAUITM_C_refr_ANT
    real(dp)            :: SMB_IMAUITM_albedo_water
    real(dp)            :: SMB_IMAUITM_albedo_soil
    real(dp)            :: SMB_IMAUITM_albedo_ice
    real(dp)            :: SMB_IMAUITM_albedo_snow

    ! Settings for the snapshot_plus_anomalies SMB model
    character(len=1024) :: SMB_snp_p_anml_filename_snapshot_T2m
    character(len=1024) :: SMB_snp_p_anml_filename_snapshot_SMB
    character(len=1024) :: SMB_snp_p_anml_filename_anomalies

    ! Settings for the ISMIP7 SMB model
    character(len=1024) :: SMB_ISMIP7_choice_SMB_baseline
    character(len=1024) :: SMB_ISMIP7_filename_SMB_baseline_fixed
    logical             :: SMB_ISMIP7_choice_SMB_offset
    character(len=1024) :: SMB_ISMIP7_filename_SMB_offset
    character(len=1024) :: SMB_ISMIP7_choice_refgeo
    character(len=1024) :: SMB_ISMIP7_forcing_foldername
    character(len=1024) :: SMB_ISMIP7_forcing_version

  ! == Basal mass balance
  ! =====================

    ! Time step
    logical             :: do_asynchronous_BMB
    real(dp)            :: dt_BMB
    real(dp)            :: dt_BMB_reinit

    ! Hard limits on melt/refreezing rates
    real(dp)            :: BMB_maximum_allowed_melt_rate
    real(dp)            :: BMB_maximum_allowed_refreezing_rate

    ! BMB transition phase
    logical             :: do_BMB_transition_phase
    real(dp)            :: BMB_transition_phase_t_start
    real(dp)            :: BMB_transition_phase_t_end

    ! Grounding line treatment
    character(len=1024) :: choice_BMB_subgrid

    ! Choice of BMB model
    character(len=1024) :: choice_BMB_model_NAM
    character(len=1024) :: choice_BMB_model_EAS
    character(len=1024) :: choice_BMB_model_GRL
    character(len=1024) :: choice_BMB_model_ANT

    ! Choice of BMB model in ROI
    character(len=1024) :: choice_BMB_model_NAM_ROI
    character(len=1024) :: choice_BMB_model_EAS_ROI
    character(len=1024) :: choice_BMB_model_GRL_ROI
    character(len=1024) :: choice_BMB_model_ANT_ROI


    ! Prescribed BMB forcing
    character(len=1024) :: choice_BMB_prescribed_NAM
    character(len=1024) :: choice_BMB_prescribed_EAS
    character(len=1024) :: choice_BMB_prescribed_GRL
    character(len=1024) :: choice_BMB_prescribed_ANT

    ! Files containing prescribed BMB forcing
    character(len=1024) :: filename_BMB_prescribed_NAM
    character(len=1024) :: filename_BMB_prescribed_EAS
    character(len=1024) :: filename_BMB_prescribed_GRL
    character(len=1024) :: filename_BMB_prescribed_ANT

    ! Timeframes for reading prescribed BMB forcing from file (set to 1E9_dp if the file has no time dimension)
    real(dp)            :: timeframe_BMB_prescribed_NAM
    real(dp)            :: timeframe_BMB_prescribed_EAS
    real(dp)            :: timeframe_BMB_prescribed_GRL
    real(dp)            :: timeframe_BMB_prescribed_ANT

    ! Choice of idealised BMB model
    character(len=1024) :: choice_BMB_model_idealised

    ! Choice of parameterised BMB model
    character(len=1024) :: choice_BMB_model_parameterised

    ! "uniform"
    real(dp)            :: uniform_BMB
    real(dp)            :: uniform_BMB_ROI
    real(dp)            :: uniform_BMB_t_start

    ! "parameterised"
    real(dp)            :: BMB_Favier2019_gamma
    real(dp)            :: BMB_Holland_Cmelt

    ! "laddie"
    character(len=1024) :: choice_BMB_laddie_system
    character(len=1024) :: filename_BMB_laddie_configname
    character(len=1024) :: filename_BMB_laddie_initial_restart
    character(len=1024) :: filename_BMB_laddie_initial_output
    character(len=1024) :: dir_BMB_laddie_model
    character(len=1024) :: conda_activate_prompt

    ! "inverted"
    real(dp)            :: BMB_inversion_t_start
    real(dp)            :: BMB_inversion_t_end

  ! == LADDIE model
  ! ===============

    ! Parallellisation
    logical             :: do_repartition_laddie

    ! Remapping
    character(len=1024) :: choice_laddie_remapping_option

    ! Output
    logical             :: do_write_laddie_output_fields
    logical             :: do_write_laddie_output_scalar
    real(dp)            :: time_interval_scalar_output

    ! Time step
    real(dp)            :: dt_laddie
    real(dp)            :: time_duration_laddie
    real(dp)            :: time_duration_laddie_init

    ! Integration
    character(len=1024) :: choice_laddie_integration_scheme
    real(dp)            :: laddie_fbrk3_beta1
    real(dp)            :: laddie_fbrk3_beta2
    real(dp)            :: laddie_fbrk3_beta3
    real(dp)            :: laddie_lfra_nu

    ! Momentum advection
    character(len=1024) :: choice_laddie_momentum_advection

    ! Initialisation
    character(len=1024) :: choice_laddie_model_initialisation

    ! Initialisation 'read_from_file'
    character(len=1024) :: filename_laddie_restart
    real(dp)            :: timeframe_laddie_restart

    ! Initialisation 'uniform'
    real(dp)            :: laddie_initial_thickness
    real(dp)            :: laddie_initial_T_offset
    real(dp)            :: laddie_initial_S_offset

    ! Equation of state
    character(len=1024) :: choice_laddie_equation_of_state
    real(dp)            :: uniform_laddie_eos_linear_alpha
    real(dp)            :: uniform_laddie_eos_linear_beta

    ! Coriolis
    character(len=1024) :: choice_laddie_coriolis
    real(dp)            :: uniform_laddie_coriolis_parameter

    ! Turbulent heat exchange
    character(len=1024) :: choice_laddie_gamma
    real(dp)            :: uniform_laddie_gamma_T

    ! Drag coefficients
    real(dp)            :: laddie_drag_coefficient_top
    real(dp)            :: laddie_drag_coefficient_mom

    ! Viscosity and diffusivity
    real(dp)            :: laddie_viscosity
    real(dp)            :: laddie_diffusivity

    ! Entrainment
    character(len=1024) :: choice_laddie_entrainment
    real(dp)            :: laddie_Holland2006_cl
    real(dp)            :: laddie_Gaspar1988_mu

    ! Stability
    real(dp)            :: laddie_thickness_minimum
    real(dp)            :: laddie_thickness_maximum
    real(dp)            :: laddie_velocity_maximum
    real(dp)            :: laddie_buoyancy_minimum

    ! Tides
    character(len=1024) :: choice_laddie_tides
    real(dp)            :: uniform_laddie_tidal_velocity

    ! Subglacial discharge (SGD)
    character(len=1024) :: choice_laddie_SGD
    character(len=1024) :: choice_laddie_SGD_idealised
    real(dp)            :: laddie_SGD_flux
    character(len=1024) :: filename_laddie_mask_SGD
    real(dp)            :: start_time_of_applying_SGD
    character(len=1024) :: transects_SGD
    character(len=1024) :: distribute_SGD

  ! == Lateral mass balance
  ! =======================

    ! Time step
    real(dp)            :: dt_LMB

    ! Choice of LMB model
    character(len=1024) :: choice_LMB_model_NAM
    character(len=1024) :: choice_LMB_model_EAS
    character(len=1024) :: choice_LMB_model_GRL
    character(len=1024) :: choice_LMB_model_ANT

    ! "uniform"
    real(dp)            :: uniform_LMB

    ! "GlacialIndex"
    character(len=1024) :: filename_LMB_GI_NAM
    character(len=1024) :: filename_LMB_GI_EAS
    character(len=1024) :: filename_LMB_GI_GRL
    character(len=1024) :: filename_LMB_GI_ANT

    real(dp)            :: warm_LMB_NAM
    real(dp)            :: warm_LMB_EAS
    real(dp)            :: warm_LMB_GRL
    real(dp)            :: warm_LMB_ANT

    real(dp)            :: cold_LMB_NAM
    real(dp)            :: cold_LMB_EAS
    real(dp)            :: cold_LMB_GRL
    real(dp)            :: cold_LMB_ANT

  ! == Glacial isostatic adjustment
  ! ===============================

    ! General settings
    character(len=1024) :: choice_GIA_model
    real(dp)            :: dt_GIA
    real(dp)            :: dx_GIA
    real(dp)            :: ELRA_lithosphere_flex_rigidity
    character(len=1024) :: choice_GIA_ELRA_relaxation_time
    real(dp)            :: ELRA_bedrock_relaxation_time
    character(len=1024) :: filename_GIA_ELRA_relaxation_time
    real(dp)            :: ELRA_mantle_density

  ! == Sea level
  ! ============

    character(len=1024) :: choice_sealevel_model
    real(dp)            :: fixed_sealevel
    character(len=1024) :: filename_prescribed_sealevel

  ! == SELEN
  ! ========

    logical             :: SELEN_run_at_t_start
    integer             :: SELEN_n_TDOF_iterations
    integer             :: SELEN_n_recursion_iterations
    logical             :: SELEN_use_rotational_feedback
    integer             :: SELEN_n_harmonics
    logical             :: SELEN_display_progress

    character(len=1024) :: SELEN_dir
    character(len=1024) :: SELEN_global_topo_filename
    character(len=1024) :: SELEN_TABOO_init_filename
    character(len=1024) :: SELEN_LMJ_VALUES_filename

    integer                  :: SELEN_irreg_time_n
    real(dp), dimension(50)  :: SELEN_irreg_time_window

    real(dp)            :: SELEN_lith_thickness
    integer             :: SELEN_visc_n
    real(dp), dimension(3) :: SELEN_visc_prof

    ! Settings for the TABOO Earth deformation model
    integer             :: SELEN_TABOO_CDE
    integer             :: SELEN_TABOO_TLOVE
    integer             :: SELEN_TABOO_DEG1
    real(dp)            :: SELEN_TABOO_RCMB

  ! == Tracer tracking
  ! ==================

    character(len=1024) :: choice_tracer_tracking_model

    ! Settings for the particle-based tracer-tracking model
    real(dp)            :: tractrackpart_dt_coupling
    real(dp)            :: tractrackpart_dx_particle
    real(dp)            :: tractrackpart_dt_particle_min
    real(dp)            :: tractrackpart_dt_particle_max
    integer             :: tractrackpart_n_max_particles
    real(dp)            :: tractrackpart_dt_new_particles
    real(dp)            :: tractrackpart_dx_new_particles
    integer             :: tractrackpart_remap_n_nearest
    logical             :: tractrackpart_write_raw_output
    real(dp)            :: tractrackpart_dt_raw_output

  ! == Output
  ! =========

    ! Basic settings
    logical             :: do_create_netcdf_output
    CHARACTER(LEN=1024) :: output_precision
    logical             :: do_compress_output
    real(dp)            :: dt_output
    real(dp)            :: dt_output_restart
    real(dp)            :: dt_output_grid
    real(dp)            :: dx_output_grid_NAM
    real(dp)            :: dx_output_grid_EAS
    real(dp)            :: dx_output_grid_GRL
    real(dp)            :: dx_output_grid_ANT
    real(dp)            :: dx_output_grid_ROI_NAM
    real(dp)            :: dx_output_grid_ROI_EAS
    real(dp)            :: dx_output_grid_ROI_GRL
    real(dp)            :: dx_output_grid_ROI_ANT

    ! ISMIP
    logical             :: do_create_ismip_output
    character(len=1024) :: ismip_scenario_name
    character(len=1024) :: ismip_group_name
    character(len=1024) :: ismip_model_name
    character(len=4)    :: ismip_member_id
    character(len=4)    :: ismip_forcing_member_id
    character(len=4)    :: ismip_counter
    character(len=1024) :: ismip_esm_name
    character(len=1024) :: ismip_contact_name
    character(len=1024) :: ismip_contact_email
    character(len=1024) :: ismip_conventions
    real(dp)            :: dt_output_ismip

    ! Transects
    character(len=1024) :: transects_NAM
    character(len=1024) :: transects_EAS
    character(len=1024) :: transects_GRL
    character(len=1024) :: transects_ANT

    ! Which data fields we want to write to the main NetCDF output files
    character(len=1024) :: choice_output_field_01
    character(len=1024) :: choice_output_field_02
    character(len=1024) :: choice_output_field_03
    character(len=1024) :: choice_output_field_04
    character(len=1024) :: choice_output_field_05
    character(len=1024) :: choice_output_field_06
    character(len=1024) :: choice_output_field_07
    character(len=1024) :: choice_output_field_08
    character(len=1024) :: choice_output_field_09
    character(len=1024) :: choice_output_field_10
    character(len=1024) :: choice_output_field_11
    character(len=1024) :: choice_output_field_12
    character(len=1024) :: choice_output_field_13
    character(len=1024) :: choice_output_field_14
    character(len=1024) :: choice_output_field_15
    character(len=1024) :: choice_output_field_16
    character(len=1024) :: choice_output_field_17
    character(len=1024) :: choice_output_field_18
    character(len=1024) :: choice_output_field_19
    character(len=1024) :: choice_output_field_20
    character(len=1024) :: choice_output_field_21
    character(len=1024) :: choice_output_field_22
    character(len=1024) :: choice_output_field_23
    character(len=1024) :: choice_output_field_24
    character(len=1024) :: choice_output_field_25
    character(len=1024) :: choice_output_field_26
    character(len=1024) :: choice_output_field_27
    character(len=1024) :: choice_output_field_28
    character(len=1024) :: choice_output_field_29
    character(len=1024) :: choice_output_field_30
    character(len=1024) :: choice_output_field_31
    character(len=1024) :: choice_output_field_32
    character(len=1024) :: choice_output_field_33
    character(len=1024) :: choice_output_field_34
    character(len=1024) :: choice_output_field_35
    character(len=1024) :: choice_output_field_36
    character(len=1024) :: choice_output_field_37
    character(len=1024) :: choice_output_field_38
    character(len=1024) :: choice_output_field_39
    character(len=1024) :: choice_output_field_40
    character(len=1024) :: choice_output_field_41
    character(len=1024) :: choice_output_field_42
    character(len=1024) :: choice_output_field_43
    character(len=1024) :: choice_output_field_44
    character(len=1024) :: choice_output_field_45
    character(len=1024) :: choice_output_field_46
    character(len=1024) :: choice_output_field_47
    character(len=1024) :: choice_output_field_48
    character(len=1024) :: choice_output_field_49
    character(len=1024) :: choice_output_field_50

  ! == Non-configurable variables
  ! =============================

    character(len=1024) :: output_dir

  ! Total mask values (used only for diagnostic output)
  ! ===================================================

    integer             :: type_icefree_land
    integer             :: type_icefree_ocean
    integer             :: type_grounded_ice
    integer             :: type_floating_ice
    integer             :: type_groundingline_gr
    integer             :: type_groundingline_fl
    integer             :: type_calvingfront_gr
    integer             :: type_calvingfront_fl
    integer             :: type_margin
    integer             :: type_coastline

  ! Ocean layers are not user-defined, but are generated by the model
  ! =================================================================

    integer                             :: nz_ocean ! Number of ocean layers
    real(dp), dimension(:), allocatable :: z_ocean  ! Depths of ocean layers

  end type type_config

contains

  subroutine read_config_file( config_filename)
    ! Use a namelist containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.
    !
    ! Note: make sure that only one process at a time calls this subroutine!

    ! In/output variables:
    character(len=*), intent(in) :: config_filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_config_file'
    character(len=1024)            :: namelist_filename
    integer                        :: config_unit, namelist_unit
    integer                        :: ios
    integer                        :: i

    ! The namelist that's used to read the external config file.
    namelist /config/&
      create_procedural_output_dir_config                         , &
      fixed_output_dir_config                                     , &
      do_unit_tests_config                                        , &
      do_benchmarks_config                                        , &
      do_check_for_NaN_config                                     , &
      do_time_display_config                                      , &
      do_write_matrix_operators_config                            , &
      do_write_checksum_log_config                                , &
      start_time_of_run_config                                    , &
      end_time_of_run_config                                      , &
      dt_coupling_config                                          , &
      do_NAM_config                                               , &
      do_EAS_config                                               , &
      do_GRL_config                                               , &
      do_ANT_config                                               , &
      lambda_M_NAM_config                                         , &
      phi_M_NAM_config                                            , &
      beta_stereo_NAM_config                                      , &
      xmin_NAM_config                                             , &
      xmax_NAM_config                                             , &
      ymin_NAM_config                                             , &
      ymax_NAM_config                                             , &
      lambda_M_EAS_config                                         , &
      phi_M_EAS_config                                            , &
      beta_stereo_EAS_config                                      , &
      xmin_EAS_config                                             , &
      xmax_EAS_config                                             , &
      ymin_EAS_config                                             , &
      ymax_EAS_config                                             , &
      lambda_M_GRL_config                                         , &
      phi_M_GRL_config                                            , &
      beta_stereo_GRL_config                                      , &
      xmin_GRL_config                                             , &
      xmax_GRL_config                                             , &
      ymin_GRL_config                                             , &
      ymax_GRL_config                                             , &
      lambda_M_ANT_config                                         , &
      phi_M_ANT_config                                            , &
      beta_stereo_ANT_config                                      , &
      xmin_ANT_config                                             , &
      xmax_ANT_config                                             , &
      ymin_ANT_config                                             , &
      ymax_ANT_config                                             , &
      refgeo_Hi_min_config                                        , &
      do_smooth_geometry_config                                   , &
      r_smooth_geometry_config                                    , &
      remove_Lake_Vostok_config                                   , &
      choice_refgeo_remapping_method_config                       , &
      choice_refgeo_init_NAM_config                               , &
      choice_refgeo_init_EAS_config                               , &
      choice_refgeo_init_GRL_config                               , &
      choice_refgeo_init_ANT_config                               , &
      choice_refgeo_init_idealised_config                         , &
      dx_refgeo_init_idealised_config                             , &
      filename_refgeo_init_NAM_config                             , &
      filename_refgeo_init_EAS_config                             , &
      filename_refgeo_init_GRL_config                             , &
      filename_refgeo_init_ANT_config                             , &
      timeframe_refgeo_init_NAM_config                            , &
      timeframe_refgeo_init_EAS_config                            , &
      timeframe_refgeo_init_GRL_config                            , &
      timeframe_refgeo_init_ANT_config                            , &
      choice_refgeo_PD_NAM_config                                 , &
      choice_refgeo_PD_EAS_config                                 , &
      choice_refgeo_PD_GRL_config                                 , &
      choice_refgeo_PD_ANT_config                                 , &
      choice_refgeo_PD_idealised_config                           , &
      dx_refgeo_PD_idealised_config                               , &
      filename_refgeo_PD_NAM_config                               , &
      filename_refgeo_PD_EAS_config                               , &
      filename_refgeo_PD_GRL_config                               , &
      filename_refgeo_PD_ANT_config                               , &
      timeframe_refgeo_PD_NAM_config                              , &
      timeframe_refgeo_PD_EAS_config                              , &
      timeframe_refgeo_PD_GRL_config                              , &
      timeframe_refgeo_PD_ANT_config                              , &
      choice_refgeo_GIAeq_NAM_config                              , &
      choice_refgeo_GIAeq_EAS_config                              , &
      choice_refgeo_GIAeq_GRL_config                              , &
      choice_refgeo_GIAeq_ANT_config                              , &
      choice_refgeo_GIAeq_idealised_config                        , &
      dx_refgeo_GIAeq_idealised_config                            , &
      filename_refgeo_GIAeq_NAM_config                            , &
      filename_refgeo_GIAeq_EAS_config                            , &
      filename_refgeo_GIAeq_GRL_config                            , &
      filename_refgeo_GIAeq_ANT_config                            , &
      timeframe_refgeo_GIAeq_NAM_config                           , &
      timeframe_refgeo_GIAeq_EAS_config                           , &
      timeframe_refgeo_GIAeq_GRL_config                           , &
      timeframe_refgeo_GIAeq_ANT_config                           , &
      refgeo_idealised_slabonaslope_Hi_config                     , &
      refgeo_idealised_slabonaslope_dhdx_config                   , &
      refgeo_idealised_Halfar_H0_config                           , &
      refgeo_idealised_Halfar_R0_config                           , &
      refgeo_idealised_Bueler_H0_config                           , &
      refgeo_idealised_Bueler_R0_config                           , &
      refgeo_idealised_Bueler_lambda_config                       , &
      refgeo_idealised_SSA_icestream_Hi_config                    , &
      refgeo_idealised_SSA_icestream_dhdx_config                  , &
      refgeo_idealised_SSA_icestream_L_config                     , &
      refgeo_idealised_SSA_icestream_m_config                     , &
      refgeo_idealised_MISMIP_mod_Hi_init_config                  , &
      refgeo_idealised_ISMIP_HOM_L_config                         , &
      refgeo_idealised_MISMIPplus_Hi_init_config                  , &
      refgeo_idealised_MISMIPplus_tune_A_config                   , &
      choice_initial_mesh_NAM_config                              , &
      choice_initial_mesh_EAS_config                              , &
      choice_initial_mesh_GRL_config                              , &
      choice_initial_mesh_ANT_config                              , &
      filename_initial_mesh_NAM_config                            , &
      filename_initial_mesh_EAS_config                            , &
      filename_initial_mesh_GRL_config                            , &
      filename_initial_mesh_ANT_config                            , &
      maximum_resolution_uniform_config                           , &
      maximum_resolution_grounded_ice_config                      , &
      maximum_resolution_floating_ice_config                      , &
      maximum_resolution_grounding_line_config                    , &
      grounding_line_width_config                                 , &
      maximum_resolution_calving_front_config                     , &
      calving_front_width_config                                  , &
      maximum_resolution_ice_front_config                         , &
      ice_front_width_config                                      , &
      maximum_resolution_coastline_config                         , &
      coastline_width_config                                      , &
      choice_regions_of_interest_config                           , &
      ROI_maximum_resolution_uniform_config                       , &
      ROI_maximum_resolution_grounded_ice_config                  , &
      ROI_maximum_resolution_floating_ice_config                  , &
      ROI_maximum_resolution_grounding_line_config                , &
      ROI_grounding_line_width_config                             , &
      ROI_maximum_resolution_calving_front_config                 , &
      ROI_calving_front_width_config                              , &
      ROI_maximum_resolution_ice_front_config                     , &
      ROI_ice_front_width_config                                  , &
      ROI_maximum_resolution_coastline_config                     , &
      ROI_coastline_width_config                                  , &
      do_refine_TransAntMounts_glaciers_config                    , &
      max_res_TransAntMounts_glaciers_config                      , &
      allow_mesh_updates_config                                   , &
      dt_mesh_update_min_config                                   , &
      minimum_mesh_fitness_coefficient_config                     , &
      do_out_of_time_calving_front_relax_config                   , &
      do_singlecore_mesh_creation_config                          , &
      alpha_min_config                                            , &
      nit_Lloyds_algorithm_config                                 , &
      mesh_resolution_tolerance_config                            , &
      dx_square_grid_smooth_NAM_config                            , &
      dx_square_grid_smooth_EAS_config                            , &
      dx_square_grid_smooth_GRL_config                            , &
      dx_square_grid_smooth_ANT_config                            , &
      choice_zeta_grid_config                                     , &
      nz_config                                                   , &
      zeta_irregular_log_R_config                                 , &
      choice_stress_balance_approximation_config                  , &
      choice_hybrid_SIASSA_scheme_config                          , &
      do_include_SSADIVA_crossterms_config                        , &
      choice_hybrid_DIVA_BPA_mask_NAM_config                      , &
      choice_hybrid_DIVA_BPA_mask_EAS_config                      , &
      choice_hybrid_DIVA_BPA_mask_GRL_config                      , &
      choice_hybrid_DIVA_BPA_mask_ANT_config                      , &
      filename_hybrid_DIVA_BPA_mask_NAM_config                    , &
      filename_hybrid_DIVA_BPA_mask_EAS_config                    , &
      filename_hybrid_DIVA_BPA_mask_GRL_config                    , &
      filename_hybrid_DIVA_BPA_mask_ANT_config                    , &
      choice_initial_velocity_NAM_config                          , &
      choice_initial_velocity_EAS_config                          , &
      choice_initial_velocity_GRL_config                          , &
      choice_initial_velocity_ANT_config                          , &
      filename_initial_velocity_NAM_config                        , &
      filename_initial_velocity_EAS_config                        , &
      filename_initial_velocity_GRL_config                        , &
      filename_initial_velocity_ANT_config                        , &
      timeframe_initial_velocity_NAM_config                       , &
      timeframe_initial_velocity_EAS_config                       , &
      timeframe_initial_velocity_GRL_config                       , &
      timeframe_initial_velocity_ANT_config                       , &
      SIA_maximum_diffusivity_config                              , &
      visc_it_norm_dUV_tol_config                                 , &
      visc_it_nit_config                                          , &
      visc_it_relax_config                                        , &
      visc_eff_min_config                                         , &
      vel_max_config                                              , &
      stress_balance_PETSc_KSPtype_config                         , &
      stress_balance_PETSc_PCtype_config                          , &
      stress_balance_PETSc_rtol_config                            , &
      stress_balance_PETSc_abstol_config                          , &
      BC_ice_front_config                                         , &
      BC_u_west_config                                            , &
      BC_u_east_config                                            , &
      BC_u_south_config                                           , &
      BC_u_north_config                                           , &
      BC_v_west_config                                            , &
      BC_v_east_config                                            , &
      BC_v_south_config                                           , &
      BC_v_north_config                                           , &
      choice_sliding_law_config                                   , &
      choice_idealised_sliding_law_config                         , &
      slid_Weertman_m_config                                      , &
      slid_Budd_q_plastic_config                                  , &
      slid_Budd_u_threshold_config                                , &
      slid_ZI_p_config                                            , &
      slid_ZI_ut_config                                           , &
      do_GL_subgrid_friction_config                               , &
      choice_subgrid_grounded_fraction_config                     , &
      do_read_bedrock_cdf_from_file_config                        , &
      subgrid_bedrock_cdf_nbins_config                            , &
      do_subgrid_friction_on_A_grid_config                        , &
      subgrid_friction_exponent_on_B_grid_config                  , &
      slid_beta_max_config                                        , &
      slid_delta_v_config                                         , &
      choice_ice_integration_method_config                        , &
      dHi_semiimplicit_fs_config                                  , &
      dHi_PETSc_KSPtype_config                                    , &
      dHi_PETSc_PCtype_config                                     , &
      dHi_PETSc_rtol_config                                       , &
      dHi_PETSc_abstol_config                                     , &
      BC_H_west_config                                            , &
      BC_H_east_config                                            , &
      BC_H_south_config                                           , &
      BC_H_north_config                                           , &
      do_target_dHi_dt_config                                     , &
      do_limit_target_dHi_dt_to_SMB_config                        , &
      target_dHi_dt_t_end_config                                  , &
      filename_dHi_dt_target_NAM_config                           , &
      filename_dHi_dt_target_EAS_config                           , &
      filename_dHi_dt_target_GRL_config                           , &
      filename_dHi_dt_target_ANT_config                           , &
      timeframe_dHi_dt_target_NAM_config                          , &
      timeframe_dHi_dt_target_EAS_config                          , &
      timeframe_dHi_dt_target_GRL_config                          , &
      timeframe_dHi_dt_target_ANT_config                          , &
      choice_timestepping_config                                  , &
      dt_ice_max_config                                           , &
      dt_ice_min_config                                           , &
      dt_ice_startup_phase_config                                 , &
      do_grounded_only_adv_dt_config                              , &
      pc_epsilon_config                                           , &
      pc_k_I_config                                               , &
      pc_k_p_config                                               , &
      pc_eta_min_config                                           , &
      pc_max_time_step_increase_config                            , &
      pc_nit_max_config                                           , &
      pc_guilty_max_config                                        , &
      pc_choice_initialise_NAM_config                             , &
      pc_choice_initialise_EAS_config                             , &
      pc_choice_initialise_GRL_config                             , &
      pc_choice_initialise_ANT_config                             , &
      filename_pc_initialise_NAM_config                           , &
      filename_pc_initialise_EAS_config                           , &
      filename_pc_initialise_GRL_config                           , &
      filename_pc_initialise_ANT_config                           , &
      timeframe_pc_initialise_NAM_config                          , &
      timeframe_pc_initialise_EAS_config                          , &
      timeframe_pc_initialise_GRL_config                          , &
      timeframe_pc_initialise_ANT_config                          , &
      choice_calving_law_config                                   , &
      calving_threshold_thickness_shelf_config                    , &
      calving_threshold_thickness_sheet_config                    , &
      max_calving_rounds_config                                   , &
      do_remove_shelves_config                                    , &
      remove_shelves_larger_than_PD_config                        , &
      continental_shelf_calving_config                            , &
      continental_shelf_min_height_config                         , &
      do_apply_ISMIP7_fracture_mask_config                        , &
      filename_ISMIP7_fracture_mask_config                        , &
      choice_mask_noice_config                                    , &
      Hi_min_config                                               , &
      Hi_thin_config                                              , &
      remove_ice_absent_at_PD_config                              , &
      geometry_relaxation_t_years_config                          , &
      do_protect_grounded_mask_config                             , &
      protect_grounded_mask_t_end_config                          , &
      do_fixiness_before_start_config                             , &
      fixiness_t_start_config                                     , &
      fixiness_t_end_config                                       , &
      fixiness_H_gl_gr_config                                     , &
      fixiness_H_gl_fl_config                                     , &
      fixiness_H_grounded_config                                  , &
      fixiness_H_floating_config                                  , &
      fixiness_H_freeland_config                                  , &
      fixiness_H_freeocean_config                                 , &
      do_limitness_before_start_config                            , &
      limitness_t_start_config                                    , &
      limitness_t_end_config                                      , &
      limitness_H_gl_gr_config                                    , &
      limitness_H_gl_fl_config                                    , &
      limitness_H_grounded_config                                 , &
      limitness_H_floating_config                                 , &
      modiness_H_style_config                                     , &
      modiness_T_hom_ref_config                                   , &
      do_asynchronous_basal_hydro_config                          , &
      dt_basal_hydro_config                                       , &
      choice_basal_hydrology_model_config                         , &
      Martin2011_hydro_Hb_min_config                              , &
      Martin2011_hydro_Hb_max_config                              , &
      basal_hydro_equil_time_config                               , &
      error_function_max_effective_pressure_config                , &
      Leguy2014_hydro_connect_exponent_config                     , &
      Salle2025_hydro_meltwater_config                            , &
      Salle2025_W_til_max_config                                  , &
      Salle2025_initial_W_til_config                              , &
      Salle2025_initial_W_config                                  , &
      Salle2025_Cd_config                                         , &
      Salle2025_phi_config                                        , &
      choice_BMB_grounded_config                                  , &
      choice_bed_roughness_config                                 , &
      choice_bed_roughness_parameterised_config                   , &
      filename_bed_roughness_NAM_config                           , &
      filename_bed_roughness_EAS_config                           , &
      filename_bed_roughness_GRL_config                           , &
      filename_bed_roughness_ANT_config                           , &
      timeframe_bed_roughness_NAM_config                          , &
      timeframe_bed_roughness_EAS_config                          , &
      timeframe_bed_roughness_GRL_config                          , &
      timeframe_bed_roughness_ANT_config                          , &
      slid_Weertman_beta_sq_uniform_config                        , &
      slid_Coulomb_phi_fric_uniform_config                        , &
      slid_Budd_phi_fric_uniform_config                           , &
      slid_Tsai2015_alpha_sq_uniform_config                       , &
      slid_Tsai2015_beta_sq_uniform_config                        , &
      slid_Schoof2005_alpha_sq_uniform_config                     , &
      slid_Schoof2005_beta_sq_uniform_config                      , &
      slid_ZI_phi_fric_uniform_config                             , &
      Martin2011till_phi_Hb_min_config                            , &
      Martin2011till_phi_Hb_max_config                            , &
      Martin2011till_phi_min_config                               , &
      Martin2011till_phi_max_config                               , &
      do_bed_roughness_nudging_config                             , &
      choice_bed_roughness_nudging_method_config                  , &
      choice_inversion_target_geometry_config                     , &
      bed_roughness_nudging_dt_config                             , &
      bed_roughness_nudging_t_start_config                        , &
      bed_roughness_nudging_t_end_config                          , &
      generic_bed_roughness_min_config                            , &
      generic_bed_roughness_max_config                            , &
      bednudge_H_dHdt_flowline_t_scale_config                     , &
      bednudge_H_dHdt_flowline_dH0_config                         , &
      bednudge_H_dHdt_flowline_dHdt0_config                       , &
      bednudge_H_dHdt_flowline_Hi_scale_config                    , &
      bednudge_H_dHdt_flowline_u_scale_config                     , &
      bednudge_H_dHdt_flowline_r_smooth_config                    , &
      bednudge_H_dHdt_flowline_w_smooth_config                    , &
      bednudge_H_dHdt_local_H0_config                             , &
      bednudge_H_dHdt_local_tau_config                            , &
      bednudge_H_dHdt_local_L_config                              , &
      bednudge_H_u_flowline_file_u_target_config                  , &
      bednudge_H_u_flowline_t_scale_config                        , &
      bednudge_H_u_flowline_H0_config                             , &
      bednudge_H_u_flowline_u0_config                             , &
      bednudge_H_u_flowline_Hi_scale_config                       , &
      bednudge_H_u_flowline_u_scale_config                        , &
      bednudge_H_u_flowline_tau_config                            , &
      bednudge_H_u_flowline_L_config                              , &
      choice_geothermal_heat_flux_config                          , &
      uniform_geothermal_heat_flux_config                         , &
      filename_geothermal_heat_flux_config                        , &
      choice_initial_ice_temperature_NAM_config                   , &
      choice_initial_ice_temperature_EAS_config                   , &
      choice_initial_ice_temperature_GRL_config                   , &
      choice_initial_ice_temperature_ANT_config                   , &
      uniform_initial_ice_temperature_NAM_config                  , &
      uniform_initial_ice_temperature_EAS_config                  , &
      uniform_initial_ice_temperature_GRL_config                  , &
      uniform_initial_ice_temperature_ANT_config                  , &
      filename_initial_ice_temperature_NAM_config                 , &
      filename_initial_ice_temperature_EAS_config                 , &
      filename_initial_ice_temperature_GRL_config                 , &
      filename_initial_ice_temperature_ANT_config                 , &
      timeframe_initial_ice_temperature_NAM_config                , &
      timeframe_initial_ice_temperature_EAS_config                , &
      timeframe_initial_ice_temperature_GRL_config                , &
      timeframe_initial_ice_temperature_ANT_config                , &
      choice_thermo_model_config                                  , &
      dt_thermodynamics_config                                    , &
      Hi_min_thermo_config                                        , &
      choice_GL_temperature_BC_config                             , &
      choice_ice_heat_capacity_config                             , &
      uniform_ice_heat_capacity_config                            , &
      choice_ice_thermal_conductivity_config                      , &
      uniform_ice_thermal_conductivity_config                     , &
      choice_flow_law_config                                      , &
      Glens_flow_law_exponent_config                              , &
      Glens_flow_law_epsilon_sq_0_config                          , &
      choice_ice_rheology_Glen_config                             , &
      uniform_Glens_flow_factor_config                            , &
      choice_enhancement_factor_transition_config                 , &
      m_enh_sheet_config                                          , &
      m_enh_shelf_config                                          , &
      do_asynchronous_climate_config                              , &
      dt_climate_config                                           , &
      choice_climate_model_NAM_config                             , &
      choice_climate_model_EAS_config                             , &
      choice_climate_model_GRL_config                             , &
      choice_climate_model_ANT_config                             , &
      choice_climate_model_idealised_config                       , &
      choice_climate_model_realistic_config                       , &
      filename_climate_snapshot_NAM_config                        , &
      filename_climate_snapshot_EAS_config                        , &
      filename_climate_snapshot_GRL_config                        , &
      filename_climate_snapshot_ANT_config                        , &
      filename_CO2_record_config                                  , &
      filename_d18O_record_config                                 , &
      do_lapse_rate_corrections_NAM_config                        , &
      do_lapse_rate_corrections_EAS_config                        , &
      do_lapse_rate_corrections_GRL_config                        , &
      do_lapse_rate_corrections_ANT_config                        , &
      lapse_rate_temp_NAM_config                                  , &
      lapse_rate_temp_EAS_config                                  , &
      lapse_rate_temp_GRL_config                                  , &
      lapse_rate_temp_ANT_config                                  , &
      lapse_rate_precip_NAM_config                                , &
      lapse_rate_precip_EAS_config                                , &
      lapse_rate_precip_GRL_config                                , &
      lapse_rate_precip_ANT_config                                , &
      filename_climate_snapshot_unif_dT_NAM_config                , &
      filename_climate_snapshot_unif_dT_EAS_config                , &
      filename_climate_snapshot_unif_dT_GRL_config                , &
      filename_climate_snapshot_unif_dT_ANT_config                , &
      uniform_deltaT_NAM_config                                   , &
      uniform_deltaT_EAS_config                                   , &
      uniform_deltaT_GRL_config                                   , &
      uniform_deltaT_ANT_config                                   , &
      precip_CC_correction_NAM_config                             , &
      precip_CC_correction_EAS_config                             , &
      precip_CC_correction_GRL_config                             , &
      precip_CC_correction_ANT_config                             , &
      filename_climate_snapshot_trans_dT_NAM_config               , &
      filename_climate_snapshot_trans_dT_EAS_config               , &
      filename_climate_snapshot_trans_dT_GRL_config               , &
      filename_climate_snapshot_trans_dT_ANT_config               , &
      filename_atmosphere_dT_NAM_config                           , &
      filename_atmosphere_dT_EAS_config                           , &
      filename_atmosphere_dT_GRL_config                           , &
      filename_atmosphere_dT_ANT_config                           , &
      climate_snp_p_anml_filename_snapshot_NAM_config             , &
      climate_snp_p_anml_filename_snapshot_EAS_config             , &
      climate_snp_p_anml_filename_snapshot_GRL_config             , &
      climate_snp_p_anml_filename_snapshot_ANT_config             , &
      climate_snp_p_anml_filename_anomalies_EAS_config            , &
      climate_snp_p_anml_filename_anomalies_NAM_config            , &
      climate_snp_p_anml_filename_anomalies_GRL_config            , &
      climate_snp_p_anml_filename_anomalies_ANT_config            , &
      choice_insolation_forcing_config                            , &
      filename_insolation_config                                  , &
      static_insolation_time_config                               , &
      choice_matrix_forcing_config                                , &
      climate_matrix_filename_PD_obs_climate_config               , &
      climate_matrix_filename_climate_snapshot_PI_config          , &
      climate_matrix_filename_climate_snapshot_warm_config        , &
      climate_matrix_filename_climate_snapshot_cold_config        , &
      climate_matrix_constant_lapserate_config                    , &
      climate_matrix_CO2vsice_NAM_config                          , &
      climate_matrix_CO2vsice_EAS_config                          , &
      climate_matrix_CO2vsice_GRL_config                          , &
      climate_matrix_CO2vsice_ANT_config                          , &
      climate_matrix_high_CO2_level_config                        , &
      climate_matrix_low_CO2_level_config                         , &
      climate_matrix_warm_orbit_time_config                       , &
      climate_matrix_cold_orbit_time_config                       , &
      climate_matrix_biascorrect_warm_config                      , &
      climate_matrix_biascorrect_cold_config                      , &
      climate_matrix_switch_glacial_index_precip_config           , &
      climate_ISMIP7_choice_baseline_config                       , &
      climate_ISMIP7_filename_baseline_config                     , &
      climate_ISMIP7_choice_refgeo_config                         , &
      climate_ISMIP7_forcing_foldername_config                    , &
      climate_ISMIP7_forcing_version_config                       , &
      do_asynchronous_ocean_config                                , &
      dt_ocean_config                                             , &
      ocean_vertical_grid_max_depth_config                        , &
      ocean_vertical_grid_dz_config                               , &
      choice_ocean_model_NAM_config                               , &
      choice_ocean_model_EAS_config                               , &
      choice_ocean_model_GRL_config                               , &
      choice_ocean_model_ANT_config                               , &
      choice_ocean_model_idealised_config                         , &
      choice_ocean_isomip_scenario_config                         , &
      ocean_tanh_deep_temperature_config                          , &
      ocean_tanh_thermocline_depth_config                         , &
      ocean_tanh_thermocline_scale_depth_config                   , &
      ocean_linear_deep_temperature_config                        , &
      ocean_linear_deep_salinity_config                           , &
      ocean_linear_reference_depth_config                         , &
      ocean_lin_therm_surf_salinity_config                        , &
      ocean_lin_therm_deep_salinity_config                        , &
      ocean_lin_therm_surf_temperature_config                     , &
      ocean_lin_therm_deep_temperature_config                     , &
      ocean_lin_therm_thermocline_top_config                      , &
      ocean_lin_therm_thermocline_bottom_config                   , &
      ocean_uniform_T_config                                      , &
      ocean_uniform_S_config                                      , &
      choice_ocean_model_realistic_config                         , &
      filename_ocean_snapshot_NAM_config                          , &
      filename_ocean_snapshot_EAS_config                          , &
      filename_ocean_snapshot_GRL_config                          , &
      filename_ocean_snapshot_ANT_config                          , &
      filename_ocean_warm_snapshot_NAM_config                     , &
      filename_ocean_warm_snapshot_EAS_config                     , &
      filename_ocean_warm_snapshot_GRL_config                     , &
      filename_ocean_warm_snapshot_ANT_config                     , &
      filename_ocean_cold_snapshot_NAM_config                     , &
      filename_ocean_cold_snapshot_EAS_config                     , &
      filename_ocean_cold_snapshot_GRL_config                     , &
      filename_ocean_cold_snapshot_ANT_config                     , &
      ocean_uniform_deltaT_NAM_config                             , &
      ocean_uniform_deltaT_EAS_config                             , &
      ocean_uniform_deltaT_GRL_config                             , &
      ocean_uniform_deltaT_ANT_config                             , &
      choice_ocean_extrapolation_method_config                    , &
      choice_ocean_model_transient_config                         , &
      filename_ocean_dT_NAM_config                                , &
      filename_ocean_dT_EAS_config                                , &
      filename_ocean_dT_GRL_config                                , &
      filename_ocean_dT_ANT_config                                , &
      filename_ocean_GI_NAM_config                                , &
      filename_ocean_GI_EAS_config                                , &
      filename_ocean_GI_GRL_config                                , &
      filename_ocean_GI_ANT_config                                , &
      ocean_snp_p_anml_filename_snapshot_config                   , &
      ocean_snp_p_anml_filename_anomalies_config                  , &
      ocean_ISMIP7_forcing_foldername_config                      , &
      ocean_ISMIP7_forcing_version_config                         , &
      do_asynchronous_SMB_config                                  , &
      dt_SMB_config                                               , &
      choice_SMB_model_NAM_config                                 , &
      choice_SMB_model_EAS_config                                 , &
      choice_SMB_model_GRL_config                                 , &
      choice_SMB_model_ANT_config                                 , &
      uniform_SMB_config                                          , &
      choice_SMB_model_idealised_config                           , &
      choice_SMB_prescribed_NAM_config                            , &
      choice_SMB_prescribed_EAS_config                            , &
      choice_SMB_prescribed_GRL_config                            , &
      choice_SMB_prescribed_ANT_config                            , &
      filename_SMB_prescribed_NAM_config                          , &
      filename_SMB_prescribed_EAS_config                          , &
      filename_SMB_prescribed_GRL_config                          , &
      filename_SMB_prescribed_ANT_config                          , &
      timeframe_SMB_prescribed_NAM_config                         , &
      timeframe_SMB_prescribed_EAS_config                         , &
      timeframe_SMB_prescribed_GRL_config                         , &
      timeframe_SMB_prescribed_ANT_config                         , &
      choice_SMB_IMAUITM_init_firn_NAM_config                     , &
      choice_SMB_IMAUITM_init_firn_EAS_config                     , &
      choice_SMB_IMAUITM_init_firn_GRL_config                     , &
      choice_SMB_IMAUITM_init_firn_ANT_config                     , &
      filename_firn_IMAUITM_NAM_config                            , &
      filename_firn_IMAUITM_EAS_config                            , &
      filename_firn_IMAUITM_GRL_config                            , &
      filename_firn_IMAUITM_ANT_config                            , &
      timeframe_restart_firn_IMAUITM_NAM_config                   , &
      timeframe_restart_firn_IMAUITM_EAS_config                   , &
      timeframe_restart_firn_IMAUITM_GRL_config                   , &
      timeframe_restart_firn_IMAUITM_ANT_config                   , &
      SMB_IMAUITM_initial_firn_thickness_config                   , &
      SMB_IMAUITM_C_abl_constant_NAM_config                       , &
      SMB_IMAUITM_C_abl_constant_EAS_config                       , &
      SMB_IMAUITM_C_abl_constant_GRL_config                       , &
      SMB_IMAUITM_C_abl_constant_ANT_config                       , &
      SMB_IMAUITM_C_abl_Ts_NAM_config                             , &
      SMB_IMAUITM_C_abl_Ts_EAS_config                             , &
      SMB_IMAUITM_C_abl_Ts_GRL_config                             , &
      SMB_IMAUITM_C_abl_Ts_ANT_config                             , &
      SMB_IMAUITM_C_abl_Q_NAM_config                              , &
      SMB_IMAUITM_C_abl_Q_EAS_config                              , &
      SMB_IMAUITM_C_abl_Q_GRL_config                              , &
      SMB_IMAUITM_C_abl_Q_ANT_config                              , &
      SMB_IMAUITM_C_refr_NAM_config                               , &
      SMB_IMAUITM_C_refr_EAS_config                               , &
      SMB_IMAUITM_C_refr_GRL_config                               , &
      SMB_IMAUITM_C_refr_ANT_config                               , &
      SMB_IMAUITM_albedo_water_config                             , &
      SMB_IMAUITM_albedo_soil_config                              , &
      SMB_IMAUITM_albedo_ice_config                               , &
      SMB_IMAUITM_albedo_snow_config                              , &
      SMB_snp_p_anml_filename_snapshot_T2m_config                 , &
      SMB_snp_p_anml_filename_snapshot_SMB_config                 , &
      SMB_snp_p_anml_filename_anomalies_config                    , &
      SMB_ISMIP7_choice_SMB_baseline_config                       , &
      SMB_ISMIP7_forcing_foldername_config                        , &
      SMB_ISMIP7_forcing_version_config                           , &
      SMB_ISMIP7_filename_SMB_baseline_fixed_config               , &
      SMB_ISMIP7_choice_SMB_offset_config                         , &
      SMB_ISMIP7_filename_SMB_offset_config                       , &
      SMB_ISMIP7_choice_refgeo_config                             , &
      do_asynchronous_BMB_config                                  , &
      dt_BMB_config                                               , &
      dt_BMB_reinit_config                                        , &
      BMB_maximum_allowed_melt_rate_config                        , &
      BMB_maximum_allowed_refreezing_rate_config                  , &
      do_BMB_transition_phase_config                              , &
      BMB_transition_phase_t_start_config                         , &
      BMB_transition_phase_t_end_config                           , &
      choice_BMB_subgrid_config                                   , &
      choice_BMB_model_NAM_config                                 , &
      choice_BMB_model_EAS_config                                 , &
      choice_BMB_model_GRL_config                                 , &
      choice_BMB_model_ANT_config                                 , &
      choice_BMB_model_NAM_ROI_config                             , &
      choice_BMB_model_EAS_ROI_config                             , &
      choice_BMB_model_GRL_ROI_config                             , &
      choice_BMB_model_ANT_ROI_config                             , &
      choice_BMB_prescribed_NAM_config                            , &
      choice_BMB_prescribed_EAS_config                            , &
      choice_BMB_prescribed_GRL_config                            , &
      choice_BMB_prescribed_ANT_config                            , &
      filename_BMB_prescribed_NAM_config                          , &
      filename_BMB_prescribed_EAS_config                          , &
      filename_BMB_prescribed_GRL_config                          , &
      filename_BMB_prescribed_ANT_config                          , &
      timeframe_BMB_prescribed_NAM_config                         , &
      timeframe_BMB_prescribed_EAS_config                         , &
      timeframe_BMB_prescribed_GRL_config                         , &
      timeframe_BMB_prescribed_ANT_config                         , &
      choice_BMB_model_idealised_config                           , &
      choice_BMB_model_parameterised_config                       , &
      uniform_BMB_config                                          , &
      uniform_BMB_ROI_config                                      , &
      uniform_BMB_t_start_config                                  , &
      BMB_Favier2019_gamma_config                                 , &
      BMB_Holland_Cmelt_config                                    , &
      choice_BMB_laddie_system_config                             , &
      filename_BMB_laddie_configname_config                       , &
      filename_BMB_laddie_initial_restart_config                  , &
      filename_BMB_laddie_initial_output_config                   , &
      dir_BMB_laddie_model_config                                 , &
      conda_activate_prompt_config                                , &
      BMB_inversion_t_start_config                                , &
      BMB_inversion_t_end_config                                  , &
      do_repartition_laddie_config                                , &
      choice_laddie_remapping_option_config                       , &
      do_write_laddie_output_fields_config                        , &
      do_write_laddie_output_scalar_config                        , &
      time_interval_scalar_output_config                          , &
      dt_laddie_config                                            , &
      time_duration_laddie_config                                 , &
      time_duration_laddie_init_config                            , &
      choice_laddie_integration_scheme_config                     , &
      laddie_fbrk3_beta1_config                                   , &
      laddie_fbrk3_beta2_config                                   , &
      laddie_fbrk3_beta3_config                                   , &
      laddie_lfra_nu_config                                       , &
      choice_laddie_momentum_advection_config                     , &
      choice_laddie_model_initialisation_config                   , &
      filename_laddie_restart_config                              , &
      timeframe_laddie_restart_config                             , &
      laddie_initial_thickness_config                             , &
      laddie_initial_T_offset_config                              , &
      laddie_initial_S_offset_config                              , &
      choice_laddie_equation_of_state_config                      , &
      uniform_laddie_eos_linear_alpha_config                      , &
      uniform_laddie_eos_linear_beta_config                       , &
      choice_laddie_coriolis_config                               , &
      uniform_laddie_coriolis_parameter_config                    , &
      choice_laddie_gamma_config                                  , &
      uniform_laddie_gamma_T_config                               , &
      laddie_drag_coefficient_top_config                          , &
      laddie_drag_coefficient_mom_config                          , &
      laddie_viscosity_config                                     , &
      laddie_diffusivity_config                                   , &
      choice_laddie_entrainment_config                            , &
      laddie_Holland2006_cl_config                                , &
      laddie_Gaspar1988_mu_config                                 , &
      laddie_thickness_minimum_config                             , &
      laddie_thickness_maximum_config                             , &
      laddie_velocity_maximum_config                              , &
      laddie_buoyancy_minimum_config                              , &
      choice_laddie_SGD_config                                    , &
      choice_laddie_SGD_idealised_config                          , &
      laddie_SGD_flux_config                                      , &
      filename_laddie_mask_SGD_config                             , &
      start_time_of_applying_SGD_config                           , &
      transects_SGD_config                                        , &
      distribute_SGD_config                                       , &
      choice_laddie_tides_config                                  , &
      uniform_laddie_tidal_velocity_config                        , &
      dt_LMB_config                                               , &
      choice_LMB_model_NAM_config                                 , &
      choice_LMB_model_EAS_config                                 , &
      choice_LMB_model_GRL_config                                 , &
      choice_LMB_model_ANT_config                                 , &
      uniform_LMB_config                                          , &
      filename_LMB_GI_NAM_config                                  , &
      filename_LMB_GI_EAS_config                                  , &
      filename_LMB_GI_GRL_config                                  , &
      filename_LMB_GI_ANT_config                                  , &
      warm_LMB_NAM_config                                         , &
      warm_LMB_EAS_config                                         , &
      warm_LMB_GRL_config                                         , &
      warm_LMB_ANT_config                                         , &
      cold_LMB_NAM_config                                         , &
      cold_LMB_EAS_config                                         , &
      cold_LMB_GRL_config                                         , &
      cold_LMB_ANT_config                                         , &
      choice_GIA_model_config                                     , &
      dt_GIA_config                                               , &
      dx_GIA_config                                               , &
      ELRA_lithosphere_flex_rigidity_config                       , &
      choice_GIA_ELRA_relaxation_time_config                      , &
      ELRA_bedrock_relaxation_time_config                         , &
      filename_GIA_ELRA_relaxation_time_config                    , &
      ELRA_mantle_density_config                                  , &
      choice_sealevel_model_config                                , &
      fixed_sealevel_config                                       , &
      filename_prescribed_sealevel_config                         , &
      SELEN_run_at_t_start_config                                 , &
      SELEN_n_TDOF_iterations_config                              , &
      SELEN_n_recursion_iterations_config                         , &
      SELEN_use_rotational_feedback_config                        , &
      SELEN_n_harmonics_config                                    , &
      SELEN_display_progress_config                               , &
      SELEN_dir_config                                            , &
      SELEN_global_topo_filename_config                           , &
      SELEN_TABOO_init_filename_config                            , &
      SELEN_LMJ_VALUES_filename_config                            , &
      SELEN_irreg_time_n_config                                   , &
      SELEN_irreg_time_window_config                              , &
      SELEN_lith_thickness_config                                 , &
      SELEN_visc_n_config                                         , &
      SELEN_visc_prof_config                                      , &
      SELEN_TABOO_CDE_config                                      , &
      SELEN_TABOO_TLOVE_config                                    , &
      SELEN_TABOO_DEG1_config                                     , &
      SELEN_TABOO_RCMB_config                                     , &
      choice_tracer_tracking_model_config                         , &
      tractrackpart_dt_coupling_config                            , &
      tractrackpart_dx_particle_config                            , &
      tractrackpart_dt_particle_min_config                        , &
      tractrackpart_dt_particle_max_config                        , &
      tractrackpart_n_max_particles_config                        , &
      tractrackpart_dt_new_particles_config                       , &
      tractrackpart_dx_new_particles_config                       , &
      tractrackpart_remap_n_nearest_config                        , &
      tractrackpart_write_raw_output_config                       , &
      tractrackpart_dt_raw_output_config                          , &
      do_create_netcdf_output_config                              , &
      output_precision_config                                     , &
      do_compress_output_config                                   , &
      dt_output_config                                            , &
      dt_output_restart_config                                    , &
      dt_output_grid_config                                       , &
      dx_output_grid_NAM_config                                   , &
      dx_output_grid_EAS_config                                   , &
      dx_output_grid_GRL_config                                   , &
      dx_output_grid_ANT_config                                   , &
      dx_output_grid_ROI_NAM_config                               , &
      dx_output_grid_ROI_EAS_config                               , &
      dx_output_grid_ROI_GRL_config                               , &
      dx_output_grid_ROI_ANT_config                               , &
      do_create_ismip_output_config                               , &
      ismip_scenario_name_config                                  , &
      ismip_group_name_config                                     , &
      ismip_model_name_config                                     , &
      ismip_member_id_config                                      , &
      ismip_forcing_member_id_config                              , &
      ismip_counter_config                                        , &
      ismip_esm_name_config                                       , &
      ismip_contact_name_config                                   , &
      ismip_contact_email_config                                  , &
      ismip_conventions_config                                    , &
      dt_output_ismip_config                                      , &
      transects_NAM_config                                        , &
      transects_EAS_config                                        , &
      transects_GRL_config                                        , &
      transects_ANT_config                                        , &
      choice_output_field_01_config                               , &
      choice_output_field_02_config                               , &
      choice_output_field_03_config                               , &
      choice_output_field_04_config                               , &
      choice_output_field_05_config                               , &
      choice_output_field_06_config                               , &
      choice_output_field_07_config                               , &
      choice_output_field_08_config                               , &
      choice_output_field_09_config                               , &
      choice_output_field_10_config                               , &
      choice_output_field_11_config                               , &
      choice_output_field_12_config                               , &
      choice_output_field_13_config                               , &
      choice_output_field_14_config                               , &
      choice_output_field_15_config                               , &
      choice_output_field_16_config                               , &
      choice_output_field_17_config                               , &
      choice_output_field_18_config                               , &
      choice_output_field_19_config                               , &
      choice_output_field_20_config                               , &
      choice_output_field_21_config                               , &
      choice_output_field_22_config                               , &
      choice_output_field_23_config                               , &
      choice_output_field_24_config                               , &
      choice_output_field_25_config                               , &
      choice_output_field_26_config                               , &
      choice_output_field_27_config                               , &
      choice_output_field_28_config                               , &
      choice_output_field_29_config                               , &
      choice_output_field_30_config                               , &
      choice_output_field_31_config                               , &
      choice_output_field_32_config                               , &
      choice_output_field_33_config                               , &
      choice_output_field_34_config                               , &
      choice_output_field_35_config                               , &
      choice_output_field_36_config                               , &
      choice_output_field_37_config                               , &
      choice_output_field_38_config                               , &
      choice_output_field_39_config                               , &
      choice_output_field_40_config                               , &
      choice_output_field_41_config                               , &
      choice_output_field_42_config                               , &
      choice_output_field_43_config                               , &
      choice_output_field_44_config                               , &
      choice_output_field_45_config                               , &
      choice_output_field_46_config                               , &
      choice_output_field_47_config                               , &
      choice_output_field_48_config                               , &
      choice_output_field_49_config                               , &
      choice_output_field_50_config
    ! End of the config namelist

    ! Add routine to path
    call init_routine( routine_name)

    ! Generate CONFIG namelist filename
    namelist_filename = config_filename
    i = index( namelist_filename,'/')
    do while (i > 0)
      namelist_filename = namelist_filename( i+1: len_trim( namelist_filename))
      i = index( namelist_filename,'/')
    end do
    namelist_filename = 'namelist_' // trim( namelist_filename)

    ! Write the config namelist to a temporary file
    open( newunit = namelist_unit, file = trim( namelist_filename))
    write( unit = namelist_unit, nml = config)
    close( unit = namelist_unit)

    ! Check the config file for validity
    call check_config_file_validity( config_filename, namelist_filename)

    ! Delete the temporary namelist file
    call system('rm -f ' // trim( namelist_filename))

    ! Open the config file
    open( newunit = config_unit, file = trim( config_filename), &
      status = 'old', action = 'read', iostat = ios)
    if (ios /= 0) call crash('couldnt open config file "' // trim( config_filename) // '"!')

    ! Read the config file using the CONFIG namelist
    read( unit = config_unit, nml = config, iostat = ios)
    if (ios /= 0) call crash('error while reading config file "' // trim( config_filename) // '"!')

    ! Close the config file
    close( unit = config_unit)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_config_file

  subroutine copy_config_variables_to_struct( C)
    ! Overwrite the values in the fields of the config structure with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values that were read from the config file.

    ! In/output variables:
    type(type_config), intent(out) :: C

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'copy_config_variables_to_struct'

    ! Add routine to path
    call init_routine( routine_name)

    ! Copy the values of the _config variables to the C structure

    ! General model instructions
    ! ==========================

    ! Output directory
    C%create_procedural_output_dir                           = create_procedural_output_dir_config
    C%fixed_output_dir                                       = fixed_output_dir_config

    ! Debugging
    C%do_unit_tests                                          = do_unit_tests_config
    C%do_benchmarks                                          = do_benchmarks_config
    C%do_check_for_NaN                                       = do_check_for_NaN_config
    C%do_time_display                                        = do_time_display_config
    C%do_write_matrix_operators                              = do_write_matrix_operators_config
    C%do_write_checksum_log                                  = do_write_checksum_log_config

    ! == Time of simulation
    ! =====================

    C%start_time_of_run                                      = start_time_of_run_config
    C%end_time_of_run                                        = end_time_of_run_config

    ! == Which model regions to simulate
    ! ==================================

    C%dt_coupling                                            = dt_coupling_config
    C%do_NAM                                                 = do_NAM_config
    C%do_EAS                                                 = do_EAS_config
    C%do_GRL                                                 = do_GRL_config
    C%do_ANT                                                 = do_ANT_config

    ! == The four model regions
    ! =========================

    ! North America
    C%lambda_M_NAM                                           = lambda_M_NAM_config
    C%phi_M_NAM                                              = phi_M_NAM_config
    C%beta_stereo_NAM                                        = beta_stereo_NAM_config
    C%xmin_NAM                                               = xmin_NAM_config
    C%xmax_NAM                                               = xmax_NAM_config
    C%ymin_NAM                                               = ymin_NAM_config
    C%ymax_NAM                                               = ymax_NAM_config

    ! Eurasia
    C%lambda_M_EAS                                           = lambda_M_EAS_config
    C%phi_M_EAS                                              = phi_M_EAS_config
    C%beta_stereo_EAS                                        = beta_stereo_EAS_config
    C%xmin_EAS                                               = xmin_EAS_config
    C%xmax_EAS                                               = xmax_EAS_config
    C%ymin_EAS                                               = ymin_EAS_config
    C%ymax_EAS                                               = ymax_EAS_config

    ! Greenland
    C%lambda_M_GRL                                           = lambda_M_GRL_config
    C%phi_M_GRL                                              = phi_M_GRL_config
    C%beta_stereo_GRL                                        = beta_stereo_GRL_config
    C%xmin_GRL                                               = xmin_GRL_config
    C%xmax_GRL                                               = xmax_GRL_config
    C%ymin_GRL                                               = ymin_GRL_config
    C%ymax_GRL                                               = ymax_GRL_config

    ! Antarctica
    C%lambda_M_ANT                                           = lambda_M_ANT_config
    C%phi_M_ANT                                              = phi_M_ANT_config
    C%beta_stereo_ANT                                        = beta_stereo_ANT_config
    C%xmin_ANT                                               = xmin_ANT_config
    C%xmax_ANT                                               = xmax_ANT_config
    C%ymin_ANT                                               = ymin_ANT_config
    C%ymax_ANT                                               = ymax_ANT_config

    ! == Reference geometries (initial, present-day, and GIA equilibrium)
    ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    C%refgeo_Hi_min                                          = refgeo_Hi_min_config
    C%do_smooth_geometry                                     = do_smooth_geometry_config
    C%r_smooth_geometry                                      = r_smooth_geometry_config
    C%remove_Lake_Vostok                                     = remove_Lake_Vostok_config

    ! Remapping method for gridded geometry input
    C%choice_refgeo_remapping_method                         = choice_refgeo_remapping_method_config

    ! == Initial geometry
    ! ===================

    C%choice_refgeo_init_NAM                                 = choice_refgeo_init_NAM_config
    C%choice_refgeo_init_EAS                                 = choice_refgeo_init_EAS_config
    C%choice_refgeo_init_GRL                                 = choice_refgeo_init_GRL_config
    C%choice_refgeo_init_ANT                                 = choice_refgeo_init_ANT_config
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    C%choice_refgeo_init_idealised                           = choice_refgeo_init_idealised_config
    C%dx_refgeo_init_idealised                               = dx_refgeo_init_idealised_config
    ! Path to file containing initial geometry when choice_refgeo_init == 'read_from_file'
    C%filename_refgeo_init_NAM                               = filename_refgeo_init_NAM_config
    C%filename_refgeo_init_EAS                               = filename_refgeo_init_EAS_config
    C%filename_refgeo_init_GRL                               = filename_refgeo_init_GRL_config
    C%filename_refgeo_init_ANT                               = filename_refgeo_init_ANT_config
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_refgeo_init_NAM                              = timeframe_refgeo_init_NAM_config
    C%timeframe_refgeo_init_EAS                              = timeframe_refgeo_init_EAS_config
    C%timeframe_refgeo_init_GRL                              = timeframe_refgeo_init_GRL_config
    C%timeframe_refgeo_init_ANT                              = timeframe_refgeo_init_ANT_config

    ! == Present-day geometry
    ! =======================

    C%choice_refgeo_PD_NAM                                   = choice_refgeo_PD_NAM_config
    C%choice_refgeo_PD_EAS                                   = choice_refgeo_PD_EAS_config
    C%choice_refgeo_PD_GRL                                   = choice_refgeo_PD_GRL_config
    C%choice_refgeo_PD_ANT                                   = choice_refgeo_PD_ANT_config
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    C%choice_refgeo_PD_idealised                             = choice_refgeo_PD_idealised_config
    C%dx_refgeo_PD_idealised                                 = dx_refgeo_PD_idealised_config
    ! Path to file containing present-day geometry when choice_refgeo_PD == 'read_from_file'
    C%filename_refgeo_PD_NAM                                 = filename_refgeo_PD_NAM_config
    C%filename_refgeo_PD_EAS                                 = filename_refgeo_PD_EAS_config
    C%filename_refgeo_PD_GRL                                 = filename_refgeo_PD_GRL_config
    C%filename_refgeo_PD_ANT                                 = filename_refgeo_PD_ANT_config
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_refgeo_PD_NAM                                = timeframe_refgeo_PD_NAM_config
    C%timeframe_refgeo_PD_EAS                                = timeframe_refgeo_PD_EAS_config
    C%timeframe_refgeo_PD_GRL                                = timeframe_refgeo_PD_GRL_config
    C%timeframe_refgeo_PD_ANT                                = timeframe_refgeo_PD_ANT_config

    ! == GIA equilibrium geometry
    ! ===========================

    C%choice_refgeo_GIAeq_NAM                                = choice_refgeo_GIAeq_NAM_config
    C%choice_refgeo_GIAeq_EAS                                = choice_refgeo_GIAeq_EAS_config
    C%choice_refgeo_GIAeq_GRL                                = choice_refgeo_GIAeq_GRL_config
    C%choice_refgeo_GIAeq_ANT                                = choice_refgeo_GIAeq_ANT_config
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    C%choice_refgeo_GIAeq_idealised                          = choice_refgeo_GIAeq_idealised_config
    C%dx_refgeo_GIAeq_idealised                              = dx_refgeo_GIAeq_idealised_config
    ! Path to file containing GIA equilibrium reference geometry when choice_refgeo_GIAeq == 'read_from_file'
    C%filename_refgeo_GIAeq_NAM                              = filename_refgeo_GIAeq_NAM_config
    C%filename_refgeo_GIAeq_EAS                              = filename_refgeo_GIAeq_EAS_config
    C%filename_refgeo_GIAeq_GRL                              = filename_refgeo_GIAeq_GRL_config
    C%filename_refgeo_GIAeq_ANT                              = filename_refgeo_GIAeq_ANT_config
    ! Timeframe to read from the geometry file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_refgeo_GIAeq_NAM                             = timeframe_refgeo_GIAeq_NAM_config
    C%timeframe_refgeo_GIAeq_EAS                             = timeframe_refgeo_GIAeq_EAS_config
    C%timeframe_refgeo_GIAeq_GRL                             = timeframe_refgeo_GIAeq_GRL_config
    C%timeframe_refgeo_GIAeq_ANT                             = timeframe_refgeo_GIAeq_ANT_config

    ! == Parameters for idealised geometries
    ! ======================================

    C%refgeo_idealised_slabonaslope_Hi                       = refgeo_idealised_slabonaslope_Hi_config
    C%refgeo_idealised_slabonaslope_dhdx                     = refgeo_idealised_slabonaslope_dhdx_config
    C%refgeo_idealised_Halfar_H0                             = refgeo_idealised_Halfar_H0_config
    C%refgeo_idealised_Halfar_R0                             = refgeo_idealised_Halfar_R0_config
    C%refgeo_idealised_Bueler_H0                             = refgeo_idealised_Bueler_H0_config
    C%refgeo_idealised_Bueler_R0                             = refgeo_idealised_Bueler_R0_config
    C%refgeo_idealised_Bueler_lambda                         = refgeo_idealised_Bueler_lambda_config
    C%refgeo_idealised_SSA_icestream_Hi                      = refgeo_idealised_SSA_icestream_Hi_config
    C%refgeo_idealised_SSA_icestream_dhdx                    = refgeo_idealised_SSA_icestream_dhdx_config
    C%refgeo_idealised_SSA_icestream_L                       = refgeo_idealised_SSA_icestream_L_config
    C%refgeo_idealised_SSA_icestream_m                       = refgeo_idealised_SSA_icestream_m_config
    C%refgeo_idealised_MISMIP_mod_Hi_init                    = refgeo_idealised_MISMIP_mod_Hi_init_config
    C%refgeo_idealised_ISMIP_HOM_L                           = refgeo_idealised_ISMIP_HOM_L_config
    C%refgeo_idealised_MISMIPplus_Hi_init                    = refgeo_idealised_MISMIPplus_Hi_init_config
    C%refgeo_idealised_MISMIPplus_tune_A                     = refgeo_idealised_MISMIPplus_tune_A_config

    ! == Mesh generation
    ! ==================

    ! How to set up the initial mesh
    C%choice_initial_mesh_NAM                                = choice_initial_mesh_NAM_config
    C%choice_initial_mesh_EAS                                = choice_initial_mesh_EAS_config
    C%choice_initial_mesh_GRL                                = choice_initial_mesh_GRL_config
    C%choice_initial_mesh_ANT                                = choice_initial_mesh_ANT_config

    ! Paths to files containing initial meshes, if choice_initial_mesh == 'read_from_file'
    C%filename_initial_mesh_NAM                              = filename_initial_mesh_NAM_config
    C%filename_initial_mesh_EAS                              = filename_initial_mesh_EAS_config
    C%filename_initial_mesh_GRL                              = filename_initial_mesh_GRL_config
    C%filename_initial_mesh_ANT                              = filename_initial_mesh_ANT_config

    ! Resolutions for different parts of the ice sheet
    C%maximum_resolution_uniform                             = maximum_resolution_uniform_config
    C%maximum_resolution_grounded_ice                        = maximum_resolution_grounded_ice_config
    C%maximum_resolution_floating_ice                        = maximum_resolution_floating_ice_config
    C%maximum_resolution_grounding_line                      = maximum_resolution_grounding_line_config
    C%grounding_line_width                                   = grounding_line_width_config
    C%maximum_resolution_calving_front                       = maximum_resolution_calving_front_config
    C%calving_front_width                                    = calving_front_width_config
    C%maximum_resolution_ice_front                           = maximum_resolution_ice_front_config
    C%ice_front_width                                        = ice_front_width_config
    C%maximum_resolution_coastline                           = maximum_resolution_coastline_config
    C%coastline_width                                        = coastline_width_config

    ! Regions of interest
    C%choice_regions_of_interest                             = choice_regions_of_interest_config
    C%ROI_maximum_resolution_uniform                         = ROI_maximum_resolution_uniform_config
    C%ROI_maximum_resolution_grounded_ice                    = ROI_maximum_resolution_grounded_ice_config
    C%ROI_maximum_resolution_floating_ice                    = ROI_maximum_resolution_floating_ice_config
    C%ROI_maximum_resolution_grounding_line                  = ROI_maximum_resolution_grounding_line_config
    C%ROI_grounding_line_width                               = ROI_grounding_line_width_config
    C%ROI_maximum_resolution_calving_front                   = ROI_maximum_resolution_calving_front_config
    C%ROI_calving_front_width                                = ROI_calving_front_width_config
    C%ROI_maximum_resolution_ice_front                       = ROI_maximum_resolution_ice_front_config
    C%ROI_ice_front_width                                    = ROI_ice_front_width_config
    C%ROI_maximum_resolution_coastline                       = ROI_maximum_resolution_coastline_config
    C%ROI_coastline_width                                    = ROI_coastline_width_config

    ! Miscellaneous refinement options
    C%do_refine_TransAntMounts_glaciers                      = do_refine_TransAntMounts_glaciers_config
    C%max_res_TransAntMounts_glaciers                        = max_res_TransAntMounts_glaciers_config

    ! Mesh update settings
    C%allow_mesh_updates                                     = allow_mesh_updates_config
    C%dt_mesh_update_min                                     = dt_mesh_update_min_config
    C%minimum_mesh_fitness_coefficient                       = minimum_mesh_fitness_coefficient_config
    C%do_out_of_time_calving_front_relax                     = do_out_of_time_calving_front_relax_config

    ! Advanced geometry parameters
    C%do_singlecore_mesh_creation                            = do_singlecore_mesh_creation_config
    C%alpha_min                                              = alpha_min_config
    C%nit_Lloyds_algorithm                                   = nit_Lloyds_algorithm_config
    C%mesh_resolution_tolerance                              = mesh_resolution_tolerance_config

    ! Square grid used for smoothing
    C%dx_square_grid_smooth_NAM                              = dx_square_grid_smooth_NAM_config
    C%dx_square_grid_smooth_EAS                              = dx_square_grid_smooth_EAS_config
    C%dx_square_grid_smooth_GRL                              = dx_square_grid_smooth_GRL_config
    C%dx_square_grid_smooth_ANT                              = dx_square_grid_smooth_ANT_config

    ! == The scaled vertical coordinate zeta
    ! ======================================

    C%choice_zeta_grid                                       = choice_zeta_grid_config
    C%nz                                                     = nz_config
    C%zeta_irregular_log_R                                   = zeta_irregular_log_R_config

    ! == Ice dynamics - velocity
    ! ==========================

    ! General
    C%choice_stress_balance_approximation                    = choice_stress_balance_approximation_config
    C%choice_hybrid_SIASSA_scheme                            = choice_hybrid_SIASSA_scheme_config
    C%do_include_SSADIVA_crossterms                          = do_include_SSADIVA_crossterms_config

    ! Hybrid DIVA/BPA
    C%choice_hybrid_DIVA_BPA_mask_NAM                        = choice_hybrid_DIVA_BPA_mask_NAM_config
    C%choice_hybrid_DIVA_BPA_mask_EAS                        = choice_hybrid_DIVA_BPA_mask_EAS_config
    C%choice_hybrid_DIVA_BPA_mask_GRL                        = choice_hybrid_DIVA_BPA_mask_GRL_config
    C%choice_hybrid_DIVA_BPA_mask_ANT                        = choice_hybrid_DIVA_BPA_mask_ANT_config

    C%filename_hybrid_DIVA_BPA_mask_NAM                      = filename_hybrid_DIVA_BPA_mask_NAM_config
    C%filename_hybrid_DIVA_BPA_mask_EAS                      = filename_hybrid_DIVA_BPA_mask_EAS_config
    C%filename_hybrid_DIVA_BPA_mask_GRL                      = filename_hybrid_DIVA_BPA_mask_GRL_config
    C%filename_hybrid_DIVA_BPA_mask_ANT                      = filename_hybrid_DIVA_BPA_mask_ANT_config

    ! Initialisation
    C%choice_initial_velocity_NAM                            = choice_initial_velocity_NAM_config
    C%choice_initial_velocity_EAS                            = choice_initial_velocity_EAS_config
    C%choice_initial_velocity_GRL                            = choice_initial_velocity_GRL_config
    C%choice_initial_velocity_ANT                            = choice_initial_velocity_ANT_config
    ! Paths to files containing initial velocity fields
    C%filename_initial_velocity_NAM                          = filename_initial_velocity_NAM_config
    C%filename_initial_velocity_EAS                          = filename_initial_velocity_EAS_config
    C%filename_initial_velocity_GRL                          = filename_initial_velocity_GRL_config
    C%filename_initial_velocity_ANT                          = filename_initial_velocity_ANT_config
    ! Timeframes to read from the initial velocity files (set to 1E9_dp if the file has no time dimension)
    C%timeframe_initial_velocity_NAM                         = timeframe_initial_velocity_NAM_config
    C%timeframe_initial_velocity_EAS                         = timeframe_initial_velocity_EAS_config
    C%timeframe_initial_velocity_GRL                         = timeframe_initial_velocity_GRL_config
    C%timeframe_initial_velocity_ANT                         = timeframe_initial_velocity_ANT_config

    ! Some parameters for numerically solving the stress balance
    C%SIA_maximum_diffusivity                                = SIA_maximum_diffusivity_config
    C%visc_it_norm_dUV_tol                                   = visc_it_norm_dUV_tol_config
    C%visc_it_nit                                            = visc_it_nit_config
    C%visc_it_relax                                          = visc_it_relax_config
    C%visc_eff_min                                           = visc_eff_min_config
    C%vel_max                                                = vel_max_config
    C%stress_balance_PETSc_KSPtype                           = stress_balance_PETSc_KSPtype_config
    C%stress_balance_PETSc_PCtype                            = stress_balance_PETSc_PCtype_config
    C%stress_balance_PETSc_rtol                              = stress_balance_PETSc_rtol_config
    C%stress_balance_PETSc_abstol                            = stress_balance_PETSc_abstol_config

    ! Boundary conditions
    C%BC_ice_front                                           = BC_ice_front_config
    C%BC_u_west                                              = BC_u_west_config
    C%BC_u_east                                              = BC_u_east_config
    C%BC_u_south                                             = BC_u_south_config
    C%BC_u_north                                             = BC_u_north_config
    C%BC_v_west                                              = BC_v_west_config
    C%BC_v_east                                              = BC_v_east_config
    C%BC_v_south                                             = BC_v_south_config
    C%BC_v_north                                             = BC_v_north_config

    ! == Ice dynamics - sliding
    ! =========================

    ! General
    C%choice_sliding_law                                     = choice_sliding_law_config
    C%choice_idealised_sliding_law                           = choice_idealised_sliding_law_config

    ! Parameters for different sliding laws
    C%slid_Weertman_m                                        = slid_Weertman_m_config
    C%slid_Budd_q_plastic                                    = slid_Budd_q_plastic_config
    C%slid_Budd_u_threshold                                  = slid_Budd_u_threshold_config
    C%slid_ZI_p                                              = slid_ZI_p_config
    C%slid_ZI_ut                                             = slid_ZI_ut_config

    ! Sub-grid scaling of basal friction
    C%do_GL_subgrid_friction                                 = do_GL_subgrid_friction_config
    C%choice_subgrid_grounded_fraction                       = choice_subgrid_grounded_fraction_config
    C%do_read_bedrock_cdf_from_file                          = do_read_bedrock_cdf_from_file_config
    C%subgrid_bedrock_cdf_nbins                              = subgrid_bedrock_cdf_nbins_config
    C%do_subgrid_friction_on_A_grid                          = do_subgrid_friction_on_A_grid_config
    C%subgrid_friction_exponent_on_B_grid                    = subgrid_friction_exponent_on_B_grid_config

    ! Stability
    C%slid_beta_max                                          = slid_beta_max_config
    C%slid_delta_v                                           = slid_delta_v_config

    ! == Ice dynamics - ice thickness calculation
    ! ===========================================

    ! Calculation of dH/dt
    C%choice_ice_integration_method                          = choice_ice_integration_method_config
    C%dHi_semiimplicit_fs                                    = dHi_semiimplicit_fs_config
    C%dHi_PETSc_KSPtype                                      = dHi_PETSc_KSPtype_config
    C%dHi_PETSc_PCtype                                       = dHi_PETSc_PCtype_config
    C%dHi_PETSc_rtol                                         = dHi_PETSc_rtol_config
    C%dHi_PETSc_abstol                                       = dHi_PETSc_abstol_config

    ! Boundary conditions
    C%BC_H_west                                              = BC_H_west_config
    C%BC_H_east                                              = BC_H_east_config
    C%BC_H_south                                             = BC_H_south_config
    C%BC_H_north                                             = BC_H_north_config

    ! == Ice dynamics - target quantities
    ! ===================================

    ! Target dHi_dt
    C%do_target_dHi_dt                                       = do_target_dHi_dt_config
    C%do_limit_target_dHi_dt_to_SMB                          = do_limit_target_dHi_dt_to_SMB_config
    C%target_dHi_dt_t_end                                    = target_dHi_dt_t_end_config

    ! Files containing a target dHi_dt for inversions
    C%filename_dHi_dt_target_NAM                             = filename_dHi_dt_target_NAM_config
    C%filename_dHi_dt_target_EAS                             = filename_dHi_dt_target_EAS_config
    C%filename_dHi_dt_target_GRL                             = filename_dHi_dt_target_GRL_config
    C%filename_dHi_dt_target_ANT                             = filename_dHi_dt_target_ANT_config

    ! Timeframes for reading target dHi_dt from file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_dHi_dt_target_NAM                            = timeframe_dHi_dt_target_NAM_config
    C%timeframe_dHi_dt_target_EAS                            = timeframe_dHi_dt_target_EAS_config
    C%timeframe_dHi_dt_target_GRL                            = timeframe_dHi_dt_target_GRL_config
    C%timeframe_dHi_dt_target_ANT                            = timeframe_dHi_dt_target_ANT_config

    ! == Ice dynamics - time stepping
    ! ===============================

    ! Time stepping
    C%choice_timestepping                                    = choice_timestepping_config
    C%dt_ice_max                                             = dt_ice_max_config
    C%dt_ice_min                                             = dt_ice_min_config
    C%dt_ice_startup_phase                                   = dt_ice_startup_phase_config
    C%do_grounded_only_adv_dt                                = do_grounded_only_adv_dt_config

    ! Predictor-corrector ice-thickness update
    C%pc_epsilon                                             = pc_epsilon_config
    C%pc_k_I                                                 = pc_k_I_config
    C%pc_k_p                                                 = pc_k_p_config
    C%pc_eta_min                                             = pc_eta_min_config
    C%pc_max_time_step_increase                              = pc_max_time_step_increase_config
    C%pc_nit_max                                             = pc_nit_max_config
    C%pc_guilty_max                                          = pc_guilty_max_config

    ! Initialisation of the predictor-corrector ice-thickness update
    C%pc_choice_initialise_NAM                               = pc_choice_initialise_NAM_config
    C%pc_choice_initialise_EAS                               = pc_choice_initialise_EAS_config
    C%pc_choice_initialise_GRL                               = pc_choice_initialise_GRL_config
    C%pc_choice_initialise_ANT                               = pc_choice_initialise_ANT_config
    ! Paths to files containing initial fields & values for the p/c scheme
    C%filename_pc_initialise_NAM                             = filename_pc_initialise_NAM_config
    C%filename_pc_initialise_EAS                             = filename_pc_initialise_EAS_config
    C%filename_pc_initialise_GRL                             = filename_pc_initialise_GRL_config
    C%filename_pc_initialise_ANT                             = filename_pc_initialise_ANT_config
    ! Timeframes to read from the p/c scheme initial file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_pc_initialise_NAM                            = timeframe_pc_initialise_NAM_config
    C%timeframe_pc_initialise_EAS                            = timeframe_pc_initialise_EAS_config
    C%timeframe_pc_initialise_GRL                            = timeframe_pc_initialise_GRL_config
    C%timeframe_pc_initialise_ANT                            = timeframe_pc_initialise_ANT_config

    ! == Ice dynamics - calving
    ! =========================

    C%choice_calving_law                                     = choice_calving_law_config
    C%calving_threshold_thickness_shelf                      = calving_threshold_thickness_shelf_config
    C%calving_threshold_thickness_sheet                      = calving_threshold_thickness_sheet_config
    C%max_calving_rounds                                     = max_calving_rounds_config
    C%do_remove_shelves                                      = do_remove_shelves_config
    C%remove_shelves_larger_than_PD                          = remove_shelves_larger_than_PD_config
    C%continental_shelf_calving                              = continental_shelf_calving_config
    C%continental_shelf_min_height                           = continental_shelf_min_height_config

    C%do_apply_ISMIP7_fracture_mask                          = do_apply_ISMIP7_fracture_mask_config
    C%filename_ISMIP7_fracture_mask                          = filename_ISMIP7_fracture_mask_config

    ! == Ice dynamics - stabilisation
    ! ===============================

    C%choice_mask_noice                                      = choice_mask_noice_config
    C%Hi_min                                                 = Hi_min_config
    C%Hi_thin                                                = Hi_thin_config
    C%remove_ice_absent_at_PD                                = remove_ice_absent_at_PD_config

    ! Geometry relaxation
    C%geometry_relaxation_t_years                            = geometry_relaxation_t_years_config

    ! Mask conservation
    C%do_protect_grounded_mask                               = do_protect_grounded_mask_config
    C%protect_grounded_mask_t_end                            = protect_grounded_mask_t_end_config

    ! Fix/delay ice thickness evolution
    C%do_fixiness_before_start                               = do_fixiness_before_start_config
    C%fixiness_t_start                                       = fixiness_t_start_config
    C%fixiness_t_end                                         = fixiness_t_end_config
    C%fixiness_H_gl_gr                                       = fixiness_H_gl_gr_config
    C%fixiness_H_gl_fl                                       = fixiness_H_gl_fl_config
    C%fixiness_H_grounded                                    = fixiness_H_grounded_config
    C%fixiness_H_floating                                    = fixiness_H_floating_config
    C%fixiness_H_freeland                                    = fixiness_H_freeland_config
    C%fixiness_H_freeocean                                   = fixiness_H_freeocean_config

    ! Limit ice thickness evolution
    C%do_limitness_before_start                              = do_limitness_before_start_config
    C%limitness_t_start                                      = limitness_t_start_config
    C%limitness_t_end                                        = limitness_t_end_config
    C%limitness_H_gl_gr                                      = limitness_H_gl_gr_config
    C%limitness_H_gl_fl                                      = limitness_H_gl_fl_config
    C%limitness_H_grounded                                   = limitness_H_grounded_config
    C%limitness_H_floating                                   = limitness_H_floating_config
    C%modiness_H_style                                       = modiness_H_style_config
    C%modiness_T_hom_ref                                     = modiness_T_hom_ref_config

    ! == Basal hydrology
    ! ==================

    ! Time step
    C%do_asynchronous_basal_hydro                            = do_asynchronous_basal_hydro_config
    C%dt_basal_hydro                                         = dt_basal_hydro_config

    ! Basal hydrology
    C%choice_basal_hydrology_model                           = choice_basal_hydrology_model_config
    C%Martin2011_hydro_Hb_min                                = Martin2011_hydro_Hb_min_config
    C%Martin2011_hydro_Hb_max                                = Martin2011_hydro_Hb_max_config
    C%basal_hydro_equil_time                                 = basal_hydro_equil_time_config
    C%error_function_max_effective_pressure                  = error_function_max_effective_pressure_config
    C%Leguy2014_hydro_connect_exponent                       = Leguy2014_hydro_connect_exponent_config
    C%Salle2025_hydro_meltwater                              = Salle2025_hydro_meltwater_config
    C%Salle2025_W_til_max                                    = Salle2025_W_til_max_config
    C%Salle2025_initial_W_til                                = Salle2025_initial_W_til_config
    C%Salle2025_initial_W                                    = Salle2025_initial_W_config
    C%Salle2025_Cd                                           = Salle2025_Cd_config
    C%Salle2025_phi                                          = Salle2025_phi_config
    C%choice_BMB_grounded                                    = choice_BMB_grounded_config

    ! == Bed roughness
    ! ==================

    C%choice_bed_roughness                                   = choice_bed_roughness_config
    C%choice_bed_roughness_parameterised                     = choice_bed_roughness_parameterised_config
    ! Paths to files containing bed roughness fields for the chosen sliding law
    C%filename_bed_roughness_NAM                             = filename_bed_roughness_NAM_config
    C%filename_bed_roughness_EAS                             = filename_bed_roughness_EAS_config
    C%filename_bed_roughness_GRL                             = filename_bed_roughness_GRL_config
    C%filename_bed_roughness_ANT                             = filename_bed_roughness_ANT_config
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_bed_roughness_NAM                            = timeframe_bed_roughness_NAM_config
    C%timeframe_bed_roughness_EAS                            = timeframe_bed_roughness_EAS_config
    C%timeframe_bed_roughness_GRL                            = timeframe_bed_roughness_GRL_config
    C%timeframe_bed_roughness_ANT                            = timeframe_bed_roughness_ANT_config
    ! Values for uniform bed roughness
    C%slid_Weertman_beta_sq_uniform                          = slid_Weertman_beta_sq_uniform_config
    C%slid_Coulomb_phi_fric_uniform                          = slid_Coulomb_phi_fric_uniform_config
    C%slid_Budd_phi_fric_uniform                             = slid_Budd_phi_fric_uniform_config
    C%slid_Tsai2015_alpha_sq_uniform                         = slid_Tsai2015_alpha_sq_uniform_config
    C%slid_Tsai2015_beta_sq_uniform                          = slid_Tsai2015_beta_sq_uniform_config
    C%slid_Schoof2005_alpha_sq_uniform                       = slid_Schoof2005_alpha_sq_uniform_config
    C%slid_Schoof2005_beta_sq_uniform                        = slid_Schoof2005_beta_sq_uniform_config
    C%slid_ZI_phi_fric_uniform                               = slid_ZI_phi_fric_uniform_config
    ! Parameters for bed roughness parameterisations
    C%Martin2011till_phi_Hb_min                              = Martin2011till_phi_Hb_min_config
    C%Martin2011till_phi_Hb_max                              = Martin2011till_phi_Hb_max_config
    C%Martin2011till_phi_min                                 = Martin2011till_phi_min_config
    C%Martin2011till_phi_max                                 = Martin2011till_phi_max_config

    ! == Bed roughness inversion by nudging
    ! =====================================

    ! General
    C%do_bed_roughness_nudging                               = do_bed_roughness_nudging_config
    C%choice_bed_roughness_nudging_method                    = choice_bed_roughness_nudging_method_config
    C%choice_inversion_target_geometry                       = choice_inversion_target_geometry_config
    C%bed_roughness_nudging_dt                               = bed_roughness_nudging_dt_config
    C%bed_roughness_nudging_t_start                          = bed_roughness_nudging_t_start_config
    C%bed_roughness_nudging_t_end                            = bed_roughness_nudging_t_end_config
    C%generic_bed_roughness_min                              = generic_bed_roughness_min_config
    C%generic_bed_roughness_max                              = generic_bed_roughness_max_config

    ! Bed roughness nudging model based on flowline-averaged values of H and dH/dt
    C%bednudge_H_dHdt_flowline_t_scale                       = bednudge_H_dHdt_flowline_t_scale_config
    C%bednudge_H_dHdt_flowline_dH0                           = bednudge_H_dHdt_flowline_dH0_config
    C%bednudge_H_dHdt_flowline_dHdt0                         = bednudge_H_dHdt_flowline_dHdt0_config
    C%bednudge_H_dHdt_flowline_Hi_scale                      = bednudge_H_dHdt_flowline_Hi_scale_config
    C%bednudge_H_dHdt_flowline_u_scale                       = bednudge_H_dHdt_flowline_u_scale_config
    C%bednudge_H_dHdt_flowline_r_smooth                      = bednudge_H_dHdt_flowline_r_smooth_config
    C%bednudge_H_dHdt_flowline_w_smooth                      = bednudge_H_dHdt_flowline_w_smooth_config

    ! Bed roughness nudging model based on local values of H and dH/dt (i.e. CISM method)
    C%bednudge_H_dHdt_local_H0                               = bednudge_H_dHdt_local_H0_config
    C%bednudge_H_dHdt_local_tau                              = bednudge_H_dHdt_local_tau_config
    C%bednudge_H_dHdt_local_L                                = bednudge_H_dHdt_local_L_config

    ! Bed roughness nudging model based on flowline-averaged values of H and u
    C%bednudge_H_u_flowline_file_u_target                    = bednudge_H_u_flowline_file_u_target_config
    C%bednudge_H_u_flowline_t_scale                          = bednudge_H_u_flowline_t_scale_config
    C%bednudge_H_u_flowline_H0                               = bednudge_H_u_flowline_H0_config
    C%bednudge_H_u_flowline_u0                               = bednudge_H_u_flowline_u0_config
    C%bednudge_H_u_flowline_Hi_scale                         = bednudge_H_u_flowline_Hi_scale_config
    C%bednudge_H_u_flowline_u_scale                          = bednudge_H_u_flowline_u_scale_config
    C%bednudge_H_u_flowline_tau                              = bednudge_H_u_flowline_tau_config
    C%bednudge_H_u_flowline_L                                = bednudge_H_u_flowline_L_config

    ! == Geothermal heat flux
    ! =======================

    C%choice_geothermal_heat_flux                            = choice_geothermal_heat_flux_config
    C%uniform_geothermal_heat_flux                           = uniform_geothermal_heat_flux_config
    C%filename_geothermal_heat_flux                          = filename_geothermal_heat_flux_config

    ! == Thermodynamics
    ! =================

    ! Initial temperature profile
    C%choice_initial_ice_temperature_NAM                     = choice_initial_ice_temperature_NAM_config
    C%choice_initial_ice_temperature_EAS                     = choice_initial_ice_temperature_EAS_config
    C%choice_initial_ice_temperature_GRL                     = choice_initial_ice_temperature_GRL_config
    C%choice_initial_ice_temperature_ANT                     = choice_initial_ice_temperature_ANT_config
    ! Uniform initial ice temperature, if choice_initial_ice_temperature == 'uniform'
    C%uniform_initial_ice_temperature_NAM                    = uniform_initial_ice_temperature_NAM_config
    C%uniform_initial_ice_temperature_EAS                    = uniform_initial_ice_temperature_EAS_config
    C%uniform_initial_ice_temperature_GRL                    = uniform_initial_ice_temperature_GRL_config
    C%uniform_initial_ice_temperature_ANT                    = uniform_initial_ice_temperature_ANT_config
    ! Paths to files containing initial temperature fields, if
    C%filename_initial_ice_temperature_NAM                   = filename_initial_ice_temperature_NAM_config
    C%filename_initial_ice_temperature_EAS                   = filename_initial_ice_temperature_EAS_config
    C%filename_initial_ice_temperature_GRL                   = filename_initial_ice_temperature_GRL_config
    C%filename_initial_ice_temperature_ANT                   = filename_initial_ice_temperature_ANT_config
    ! Timeframes to read from the bed roughness file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_initial_ice_temperature_NAM                  = timeframe_initial_ice_temperature_NAM_config
    C%timeframe_initial_ice_temperature_EAS                  = timeframe_initial_ice_temperature_EAS_config
    C%timeframe_initial_ice_temperature_GRL                  = timeframe_initial_ice_temperature_GRL_config
    C%timeframe_initial_ice_temperature_ANT                  = timeframe_initial_ice_temperature_ANT_config
    ! Thermodynamical model
    C%choice_thermo_model                                    = choice_thermo_model_config
    C%dt_thermodynamics                                      = dt_thermodynamics_config
    C%Hi_min_thermo                                          = Hi_min_thermo_config
    C%choice_GL_temperature_BC                               = choice_GL_temperature_BC_config
    C%choice_ice_heat_capacity                               = choice_ice_heat_capacity_config
    C%uniform_ice_heat_capacity                              = uniform_ice_heat_capacity_config
    C%choice_ice_thermal_conductivity                        = choice_ice_thermal_conductivity_config
    C%uniform_ice_thermal_conductivity                       = uniform_ice_thermal_conductivity_config

    ! == Rheology and flow law
    ! =========================

    ! Flow law
    C%choice_flow_law                                        = choice_flow_law_config
    C%Glens_flow_law_exponent                                = Glens_flow_law_exponent_config
    C%Glens_flow_law_epsilon_sq_0                            = Glens_flow_law_epsilon_sq_0_config

    ! Rheology
    C%choice_ice_rheology_Glen                               = choice_ice_rheology_Glen_config
    C%uniform_Glens_flow_factor                              = uniform_Glens_flow_factor_config

    ! Enhancement factors
    C%choice_enhancement_factor_transition                   = choice_enhancement_factor_transition_config
    C%m_enh_sheet                                            = m_enh_sheet_config
    C%m_enh_shelf                                            = m_enh_shelf_config

    ! == Climate
    ! ==========

    ! Time step
    C%do_asynchronous_climate                                = do_asynchronous_climate_config
    C%dt_climate                                             = dt_climate_config

    ! Choice of climate model
    C%choice_climate_model_NAM                               = choice_climate_model_NAM_config
    C%choice_climate_model_EAS                               = choice_climate_model_EAS_config
    C%choice_climate_model_GRL                               = choice_climate_model_GRL_config
    C%choice_climate_model_ANT                               = choice_climate_model_ANT_config

    ! Choice of idealised climate model
    C%choice_climate_model_idealised                         = choice_climate_model_idealised_config

    ! Choice of realistic climate model
    C%choice_climate_model_realistic                         = choice_climate_model_realistic_config

    ! Paths to files containing fields for realistic climates
    C%filename_climate_snapshot_NAM                          = filename_climate_snapshot_NAM_config
    C%filename_climate_snapshot_EAS                          = filename_climate_snapshot_EAS_config
    C%filename_climate_snapshot_GRL                          = filename_climate_snapshot_GRL_config
    C%filename_climate_snapshot_ANT                          = filename_climate_snapshot_ANT_config

    ! Climate - global forcings
    C%filename_CO2_record                                    = filename_CO2_record_config
    C%filename_d18O_record                                   = filename_d18O_record_config

    ! Lapse rates
    C%do_lapse_rate_corrections_NAM                          = do_lapse_rate_corrections_NAM_config
    C%do_lapse_rate_corrections_EAS                          = do_lapse_rate_corrections_EAS_config
    C%do_lapse_rate_corrections_GRL                          = do_lapse_rate_corrections_GRL_config
    C%do_lapse_rate_corrections_ANT                          = do_lapse_rate_corrections_ANT_config
    C%lapse_rate_temp_NAM                                    = lapse_rate_temp_NAM_config
    C%lapse_rate_temp_EAS                                    = lapse_rate_temp_EAS_config
    C%lapse_rate_temp_GRL                                    = lapse_rate_temp_GRL_config
    C%lapse_rate_temp_ANT                                    = lapse_rate_temp_ANT_config
    C%lapse_rate_precip_NAM                                  = lapse_rate_precip_NAM_config
    C%lapse_rate_precip_EAS                                  = lapse_rate_precip_EAS_config
    C%lapse_rate_precip_GRL                                  = lapse_rate_precip_GRL_config
    C%lapse_rate_precip_ANT                                  = lapse_rate_precip_ANT_config

    ! Climate - snapshot plus uniform deltaT
    C%filename_climate_snapshot_unif_dT_NAM                  = filename_climate_snapshot_unif_dT_NAM_config
    C%filename_climate_snapshot_unif_dT_EAS                  = filename_climate_snapshot_unif_dT_EAS_config
    C%filename_climate_snapshot_unif_dT_GRL                  = filename_climate_snapshot_unif_dT_GRL_config
    C%filename_climate_snapshot_unif_dT_ANT                  = filename_climate_snapshot_unif_dT_ANT_config
    C%uniform_deltaT_NAM                                     = uniform_deltaT_NAM_config
    C%uniform_deltaT_EAS                                     = uniform_deltaT_EAS_config
    C%uniform_deltaT_GRL                                     = uniform_deltaT_GRL_config
    C%uniform_deltaT_ANT                                     = uniform_deltaT_ANT_config
    C%precip_CC_correction_NAM                               = precip_CC_correction_NAM_config
    C%precip_CC_correction_EAS                               = precip_CC_correction_EAS_config
    C%precip_CC_correction_GRL                               = precip_CC_correction_GRL_config
    C%precip_CC_correction_ANT                               = precip_CC_correction_ANT_config

    ! == Climate - snapshot plus a transient deltaT
    c%filename_climate_snapshot_trans_dT_NAM                 = filename_climate_snapshot_trans_dT_NAM_config
    c%filename_climate_snapshot_trans_dT_EAS                 = filename_climate_snapshot_trans_dT_EAS_config
    c%filename_climate_snapshot_trans_dT_GRL                 = filename_climate_snapshot_trans_dT_GRL_config
    c%filename_climate_snapshot_trans_dT_ANT                 = filename_climate_snapshot_trans_dT_ANT_config
    C%filename_atmosphere_dT_NAM                             = filename_atmosphere_dT_NAM_config
    C%filename_atmosphere_dT_EAS                             = filename_atmosphere_dT_EAS_config
    C%filename_atmosphere_dT_GRL                             = filename_atmosphere_dT_GRL_config
    C%filename_atmosphere_dT_ANT                             = filename_atmosphere_dT_ANT_config

    ! == Climate snapshot_plus_anomalies
    C%climate_snp_p_anml_filename_snapshot_NAM               = climate_snp_p_anml_filename_snapshot_NAM_config
    C%climate_snp_p_anml_filename_snapshot_EAS               = climate_snp_p_anml_filename_snapshot_EAS_config
    C%climate_snp_p_anml_filename_snapshot_GRL               = climate_snp_p_anml_filename_snapshot_GRL_config
    C%climate_snp_p_anml_filename_snapshot_ANT               = climate_snp_p_anml_filename_snapshot_ANT_config
    C%climate_snp_p_anml_filename_anomalies_NAM              = climate_snp_p_anml_filename_anomalies_NAM_config
    C%climate_snp_p_anml_filename_anomalies_EAS              = climate_snp_p_anml_filename_anomalies_EAS_config
    C%climate_snp_p_anml_filename_anomalies_GRL              = climate_snp_p_anml_filename_anomalies_GRL_config
    C%climate_snp_p_anml_filename_anomalies_ANT              = climate_snp_p_anml_filename_anomalies_ANT_config

    C%choice_insolation_forcing                              = choice_insolation_forcing_config
    C%filename_insolation                                    = filename_insolation_config
    C%static_insolation_time                                 = static_insolation_time_config

    ! Climate matrix
    C%choice_matrix_forcing                                  = choice_matrix_forcing_config
    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    C%climate_matrix_filename_PD_obs_climate                   = climate_matrix_filename_PD_obs_climate_config

    ! GCM snapshots in the matrix_warm_cold option
    C%climate_matrix_filename_climate_snapshot_PI              = climate_matrix_filename_climate_snapshot_PI_config
    C%climate_matrix_filename_climate_snapshot_warm            = climate_matrix_filename_climate_snapshot_warm_config
    C%climate_matrix_filename_climate_snapshot_cold            = climate_matrix_filename_climate_snapshot_cold_config

    C%climate_matrix_constant_lapserate                        = climate_matrix_constant_lapserate_config

    ! Scaling factor for CO2 vs ice weights
    C%climate_matrix_CO2vsice_NAM                              = climate_matrix_CO2vsice_NAM_config
    C%climate_matrix_CO2vsice_EAS                              = climate_matrix_CO2vsice_EAS_config
    C%climate_matrix_CO2vsice_GRL                              = climate_matrix_CO2vsice_GRL_config
    C%climate_matrix_CO2vsice_ANT                              = climate_matrix_CO2vsice_ANT_config

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    C%climate_matrix_high_CO2_level                            = climate_matrix_high_CO2_level_config
    C%climate_matrix_low_CO2_level                             = climate_matrix_low_CO2_level_config
    C%climate_matrix_warm_orbit_time                           = climate_matrix_warm_orbit_time_config
    C%climate_matrix_cold_orbit_time                           = climate_matrix_cold_orbit_time_config

    ! Whether or not to apply a bias correction to the GCM snapshots
    C%climate_matrix_biascorrect_warm                          = climate_matrix_biascorrect_warm_config
    C%climate_matrix_biascorrect_cold                          = climate_matrix_biascorrect_cold_config
    C%climate_matrix_switch_glacial_index_precip               = climate_matrix_switch_glacial_index_precip_config

    ! Settings for the ISMIP7 climate model
    C%climate_ISMIP7_choice_baseline                           = climate_ISMIP7_choice_baseline_config
    C%climate_ISMIP7_filename_baseline                         = climate_ISMIP7_filename_baseline_config
    C%climate_ISMIP7_choice_refgeo                             = climate_ISMIP7_choice_refgeo_config
    C%climate_ISMIP7_forcing_foldername                        = climate_ISMIP7_forcing_foldername_config
    C%climate_ISMIP7_forcing_version                           = climate_ISMIP7_forcing_version_config

    ! == Ocean
    ! ========

    ! Time step
    C%do_asynchronous_ocean                                  = do_asynchronous_ocean_config
    C%dt_ocean                                               = dt_ocean_config

    ! Vertical grid
    C%ocean_vertical_grid_max_depth                          = ocean_vertical_grid_max_depth_config
    C%ocean_vertical_grid_dz                                 = ocean_vertical_grid_dz_config

    ! Choice of ocean model
    C%choice_ocean_model_NAM                                 = choice_ocean_model_NAM_config
    C%choice_ocean_model_EAS                                 = choice_ocean_model_EAS_config
    C%choice_ocean_model_GRL                                 = choice_ocean_model_GRL_config
    C%choice_ocean_model_ANT                                 = choice_ocean_model_ANT_config

    ! Choice of idealised ocean model
    C%choice_ocean_model_idealised                           = choice_ocean_model_idealised_config
    C%choice_ocean_isomip_scenario                           = choice_ocean_isomip_scenario_config
    C%ocean_tanh_deep_temperature                            = ocean_tanh_deep_temperature_config
    C%ocean_tanh_thermocline_depth                           = ocean_tanh_thermocline_depth_config
    C%ocean_tanh_thermocline_scale_depth                     = ocean_tanh_thermocline_scale_depth_config
    C%ocean_linear_deep_temperature                          = ocean_linear_deep_temperature_config
    C%ocean_linear_deep_salinity                             = ocean_linear_deep_salinity_config
    C%ocean_linear_reference_depth                           = ocean_linear_reference_depth_config
    C%ocean_lin_therm_surf_salinity                          = ocean_lin_therm_surf_salinity_config
    C%ocean_lin_therm_deep_salinity                          = ocean_lin_therm_deep_salinity_config
    C%ocean_lin_therm_surf_temperature                       = ocean_lin_therm_surf_temperature_config
    C%ocean_lin_therm_deep_temperature                       = ocean_lin_therm_deep_temperature_config
    C%ocean_lin_therm_thermocline_top                        = ocean_lin_therm_thermocline_top_config
    C%ocean_lin_therm_thermocline_bottom                     = ocean_lin_therm_thermocline_bottom_config
    C%ocean_uniform_T                                        = ocean_uniform_T_config
    C%ocean_uniform_S                                        = ocean_uniform_S_config

    ! Choice of realistic ocean model
    C%choice_ocean_model_realistic                           = choice_ocean_model_realistic_config

    ! Paths to files containing fields for realistic ocean
    C%filename_ocean_snapshot_NAM                            = filename_ocean_snapshot_NAM_config
    C%filename_ocean_snapshot_EAS                            = filename_ocean_snapshot_EAS_config
    C%filename_ocean_snapshot_GRL                            = filename_ocean_snapshot_GRL_config
    C%filename_ocean_snapshot_ANT                            = filename_ocean_snapshot_ANT_config

    C%filename_ocean_warm_snapshot_NAM                       = filename_ocean_warm_snapshot_NAM_config
    C%filename_ocean_warm_snapshot_EAS                       = filename_ocean_warm_snapshot_EAS_config
    C%filename_ocean_warm_snapshot_GRL                       = filename_ocean_warm_snapshot_GRL_config
    C%filename_ocean_warm_snapshot_ANT                       = filename_ocean_warm_snapshot_ANT_config

    C%filename_ocean_cold_snapshot_NAM                       = filename_ocean_cold_snapshot_NAM_config
    C%filename_ocean_cold_snapshot_EAS                       = filename_ocean_cold_snapshot_EAS_config
    C%filename_ocean_cold_snapshot_GRL                       = filename_ocean_cold_snapshot_GRL_config
    C%filename_ocean_cold_snapshot_ANT                       = filename_ocean_cold_snapshot_ANT_config

    ! Choice of ocean deltaT
    C%ocean_uniform_deltaT_NAM                               = ocean_uniform_deltaT_NAM_config
    C%ocean_uniform_deltaT_EAS                               = ocean_uniform_deltaT_EAS_config
    C%ocean_uniform_deltaT_GRL                               = ocean_uniform_deltaT_GRL_config
    C%ocean_uniform_deltaT_ANT                               = ocean_uniform_deltaT_ANT_config

    ! Choice of extrapolation method
    C%choice_ocean_extrapolation_method                      = choice_ocean_extrapolation_method_config

    ! Choice of transient ocean model
    C%choice_ocean_model_transient                           = choice_ocean_model_transient_config

    ! Paths to files containing the deltaT record for the transient deltaT ocean model
    C%filename_ocean_dT_NAM                                  = filename_ocean_dT_NAM_config
    C%filename_ocean_dT_EAS                                  = filename_ocean_dT_EAS_config
    C%filename_ocean_dT_GRL                                  = filename_ocean_dT_GRL_config
    C%filename_ocean_dT_ANT                                  = filename_ocean_dT_ANT_config

    ! Paths to files containing the GI record for the transient GlacialIndex ocean model
    C%filename_ocean_GI_NAM                                  = filename_ocean_GI_NAM_config
    C%filename_ocean_GI_EAS                                  = filename_ocean_GI_EAS_config
    C%filename_ocean_GI_GRL                                  = filename_ocean_GI_GRL_config
    C%filename_ocean_GI_ANT                                  = filename_ocean_GI_ANT_config

    ! Settings for the snapshot_plus_anomalies ocean model
    C%ocean_snp_p_anml_filename_snapshot                     = ocean_snp_p_anml_filename_snapshot_config
    C%ocean_snp_p_anml_filename_anomalies                    = ocean_snp_p_anml_filename_anomalies_config

    ! Settings for the ISMIP7 ocean model
    C%ocean_ISMIP7_forcing_foldername                        = ocean_ISMIP7_forcing_foldername_config
    C%ocean_ISMIP7_forcing_version                           = ocean_ISMIP7_forcing_version_config

    ! == Surface mass balance
    ! =======================

    ! Time step
    C%do_asynchronous_SMB                                    = do_asynchronous_SMB_config
    C%dt_SMB                                                 = dt_SMB_config

    ! Choice of SMB model
    C%choice_SMB_model_NAM                                   = choice_SMB_model_NAM_config
    C%choice_SMB_model_EAS                                   = choice_SMB_model_EAS_config
    C%choice_SMB_model_GRL                                   = choice_SMB_model_GRL_config
    C%choice_SMB_model_ANT                                   = choice_SMB_model_ANT_config

    ! Value to be used for uniform SMB (no regional variants, only used for idealised-geometry experiments)
    C%uniform_SMB                                            = uniform_SMB_config

    ! Choice of idealised SMB model
    C%choice_SMB_model_idealised                             = choice_SMB_model_idealised_config

    ! Prescribed SMB forcing
    C%choice_SMB_prescribed_NAM                              = choice_SMB_prescribed_NAM_config
    C%choice_SMB_prescribed_EAS                              = choice_SMB_prescribed_EAS_config
    C%choice_SMB_prescribed_GRL                              = choice_SMB_prescribed_GRL_config
    C%choice_SMB_prescribed_ANT                              = choice_SMB_prescribed_ANT_config

    ! Files containing prescribed SMB forcing
    C%filename_SMB_prescribed_NAM                            = filename_SMB_prescribed_NAM_config
    C%filename_SMB_prescribed_EAS                            = filename_SMB_prescribed_EAS_config
    C%filename_SMB_prescribed_GRL                            = filename_SMB_prescribed_GRL_config
    C%filename_SMB_prescribed_ANT                            = filename_SMB_prescribed_ANT_config

    ! Timeframes for reading prescribed SMB forcing from file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_SMB_prescribed_NAM                           = timeframe_SMB_prescribed_NAM_config
    C%timeframe_SMB_prescribed_EAS                           = timeframe_SMB_prescribed_EAS_config
    C%timeframe_SMB_prescribed_GRL                           = timeframe_SMB_prescribed_GRL_config
    C%timeframe_SMB_prescribed_ANT                           = timeframe_SMB_prescribed_ANT_config

    ! IMAU-ITM SMB model
    C%choice_SMB_IMAUITM_init_firn_NAM                       = choice_SMB_IMAUITM_init_firn_NAM_config
    C%choice_SMB_IMAUITM_init_firn_EAS                       = choice_SMB_IMAUITM_init_firn_EAS_config
    C%choice_SMB_IMAUITM_init_firn_GRL                       = choice_SMB_IMAUITM_init_firn_GRL_config
    C%choice_SMB_IMAUITM_init_firn_ANT                       = choice_SMB_IMAUITM_init_firn_ANT_config

    ! Files containing the firn model (yearly firn depth and melt)
    C%filename_firn_IMAUITM_NAM                              = filename_firn_IMAUITM_NAM_config
    C%filename_firn_IMAUITM_EAS                              = filename_firn_IMAUITM_EAS_config
    C%filename_firn_IMAUITM_GRL                              = filename_firn_IMAUITM_GRL_config
    C%filename_firn_IMAUITM_ANT                              = filename_firn_IMAUITM_ANT_config

    ! timeframe for restarting from the firn model
    C%timeframe_restart_firn_IMAUITM_NAM                    = timeframe_restart_firn_IMAUITM_NAM_config
    C%timeframe_restart_firn_IMAUITM_EAS                    = timeframe_restart_firn_IMAUITM_EAS_config
    C%timeframe_restart_firn_IMAUITM_GRL                    = timeframe_restart_firn_IMAUITM_GRL_config
    C%timeframe_restart_firn_IMAUITM_ANT                    = timeframe_restart_firn_IMAUITM_ANT_config

    ! Tuning parameters for the IMAU-ITM SMB model
    C%SMB_IMAUITM_initial_firn_thickness                     = SMB_IMAUITM_initial_firn_thickness_config
    C%SMB_IMAUITM_C_abl_constant_NAM                         = SMB_IMAUITM_C_abl_constant_NAM_config
    C%SMB_IMAUITM_C_abl_constant_EAS                         = SMB_IMAUITM_C_abl_constant_EAS_config
    C%SMB_IMAUITM_C_abl_constant_GRL                         = SMB_IMAUITM_C_abl_constant_GRL_config
    C%SMB_IMAUITM_C_abl_constant_ANT                         = SMB_IMAUITM_C_abl_constant_ANT_config
    C%SMB_IMAUITM_C_abl_Ts_NAM                               = SMB_IMAUITM_C_abl_Ts_NAM_config
    C%SMB_IMAUITM_C_abl_Ts_EAS                               = SMB_IMAUITM_C_abl_Ts_EAS_config
    C%SMB_IMAUITM_C_abl_Ts_GRL                               = SMB_IMAUITM_C_abl_Ts_GRL_config
    C%SMB_IMAUITM_C_abl_Ts_ANT                               = SMB_IMAUITM_C_abl_Ts_ANT_config
    C%SMB_IMAUITM_C_abl_Q_NAM                                = SMB_IMAUITM_C_abl_Q_NAM_config
    C%SMB_IMAUITM_C_abl_Q_EAS                                = SMB_IMAUITM_C_abl_Q_EAS_config
    C%SMB_IMAUITM_C_abl_Q_GRL                                = SMB_IMAUITM_C_abl_Q_GRL_config
    C%SMB_IMAUITM_C_abl_Q_ANT                                = SMB_IMAUITM_C_abl_Q_ANT_config
    C%SMB_IMAUITM_C_refr_NAM                                 = SMB_IMAUITM_C_refr_NAM_config
    C%SMB_IMAUITM_C_refr_EAS                                 = SMB_IMAUITM_C_refr_EAS_config
    C%SMB_IMAUITM_C_refr_GRL                                 = SMB_IMAUITM_C_refr_GRL_config
    C%SMB_IMAUITM_C_refr_ANT                                 = SMB_IMAUITM_C_refr_ANT_config
    C%SMB_IMAUITM_albedo_water                               = SMB_IMAUITM_albedo_water_config
    C%SMB_IMAUITM_albedo_soil                                = SMB_IMAUITM_albedo_soil_config
    C%SMB_IMAUITM_albedo_ice                                 = SMB_IMAUITM_albedo_ice_config
    c%SMB_IMAUITM_albedo_snow                                = SMB_IMAUITM_albedo_snow_config

    ! Settings for the snapshot_plus_anomalies SMB model
    C%SMB_snp_p_anml_filename_snapshot_T2m                   = SMB_snp_p_anml_filename_snapshot_T2m_config
    C%SMB_snp_p_anml_filename_snapshot_SMB                   = SMB_snp_p_anml_filename_snapshot_SMB_config
    C%SMB_snp_p_anml_filename_anomalies                      = SMB_snp_p_anml_filename_anomalies_config

    ! Settings for the ISMIP7 SMB model
    C%SMB_ISMIP7_choice_SMB_baseline                         = SMB_ISMIP7_choice_SMB_baseline_config
    C%SMB_ISMIP7_filename_SMB_baseline_fixed                 = SMB_ISMIP7_filename_SMB_baseline_fixed_config
    C%SMB_ISMIP7_choice_SMB_offset                           = SMB_ISMIP7_choice_SMB_offset_config
    C%SMB_ISMIP7_filename_SMB_offset                         = SMB_ISMIP7_filename_SMB_offset_config
    C%SMB_ISMIP7_choice_refgeo                               = SMB_ISMIP7_choice_refgeo_config
    C%SMB_ISMIP7_forcing_foldername                          = SMB_ISMIP7_forcing_foldername_config
    C%SMB_ISMIP7_forcing_version                             = SMB_ISMIP7_forcing_version_config

    ! == Basal mass balance
    ! =====================

    ! Time step
    C%do_asynchronous_BMB                                    = do_asynchronous_BMB_config
    C%dt_BMB                                                 = dt_BMB_config
    C%dt_BMB_reinit                                          = dt_BMB_reinit_config

    ! Hard limits on melt/refreezing rates
    C%BMB_maximum_allowed_melt_rate                          = BMB_maximum_allowed_melt_rate_config
    C%BMB_maximum_allowed_refreezing_rate                    = BMB_maximum_allowed_refreezing_rate_config

    ! BMB transition phase
    C%do_BMB_transition_phase                                = do_BMB_transition_phase_config
    C%BMB_transition_phase_t_start                           = BMB_transition_phase_t_start_config
    C%BMB_transition_phase_t_end                             = BMB_transition_phase_t_end_config

    ! Grounding line treatment
    C%choice_BMB_subgrid                                     = choice_BMB_subgrid_config

    ! Choice of BMB model
    C%choice_BMB_model_NAM                                   = choice_BMB_model_NAM_config
    C%choice_BMB_model_EAS                                   = choice_BMB_model_EAS_config
    C%choice_BMB_model_GRL                                   = choice_BMB_model_GRL_config
    C%choice_BMB_model_ANT                                   = choice_BMB_model_ANT_config

    ! Choice of BMB model in ROI
    C%choice_BMB_model_NAM_ROI                               = choice_BMB_model_NAM_ROI_config
    C%choice_BMB_model_EAS_ROI                               = choice_BMB_model_EAS_ROI_config
    C%choice_BMB_model_GRL_ROI                               = choice_BMB_model_GRL_ROI_config
    C%choice_BMB_model_ANT_ROI                               = choice_BMB_model_ANT_ROI_config

    ! Prescribed BMB forcing
    C%choice_BMB_prescribed_NAM                              = choice_BMB_prescribed_NAM_config
    C%choice_BMB_prescribed_EAS                              = choice_BMB_prescribed_EAS_config
    C%choice_BMB_prescribed_GRL                              = choice_BMB_prescribed_GRL_config
    C%choice_BMB_prescribed_ANT                              = choice_BMB_prescribed_ANT_config

    ! Files containing prescribed BMB forcing
    C%filename_BMB_prescribed_NAM                            = filename_BMB_prescribed_NAM_config
    C%filename_BMB_prescribed_EAS                            = filename_BMB_prescribed_EAS_config
    C%filename_BMB_prescribed_GRL                            = filename_BMB_prescribed_GRL_config
    C%filename_BMB_prescribed_ANT                            = filename_BMB_prescribed_ANT_config

    ! Timeframes for reading prescribed BMB forcing from file (set to 1E9_dp if the file has no time dimension)
    C%timeframe_BMB_prescribed_NAM                           = timeframe_BMB_prescribed_NAM_config
    C%timeframe_BMB_prescribed_EAS                           = timeframe_BMB_prescribed_EAS_config
    C%timeframe_BMB_prescribed_GRL                           = timeframe_BMB_prescribed_GRL_config
    C%timeframe_BMB_prescribed_ANT                           = timeframe_BMB_prescribed_ANT_config

    ! Choice of idealised BMB model
    C%choice_BMB_model_idealised                             = choice_BMB_model_idealised_config

    ! Choice of parameterised BMB model
    C%choice_BMB_model_parameterised                         = choice_BMB_model_parameterised_config

    ! "uniform"
    C%uniform_BMB                                            = uniform_BMB_config
    C%uniform_BMB_ROI                                        = uniform_BMB_ROI_config
    C%uniform_BMB_t_start                                    = uniform_BMB_t_start_config

    ! "parameterised"
    C%BMB_Favier2019_gamma                                   = BMB_Favier2019_gamma_config
    C%BMB_Holland_Cmelt                                      = BMB_Holland_Cmelt_config

    ! "laddie"
    C%choice_BMB_laddie_system                               = choice_BMB_laddie_system_config
    C%filename_BMB_laddie_configname                         = filename_BMB_laddie_configname_config
    C%filename_BMB_laddie_initial_restart                    = filename_BMB_laddie_initial_restart_config
    C%filename_BMB_laddie_initial_output                     = filename_BMB_laddie_initial_output_config
    C%dir_BMB_laddie_model                                   = dir_BMB_laddie_model_config
    C%conda_activate_prompt                                  = conda_activate_prompt_config

    ! "inverted|
    C%BMB_inversion_t_start                                  = BMB_inversion_t_start_config
    C%BMB_inversion_t_end                                    = BMB_inversion_t_end_config

    ! == LADDIE model
    ! ===============

    ! Parallellisation
    C%do_repartition_laddie                                  = do_repartition_laddie_config

    ! Remapping
    C%choice_laddie_remapping_option                         = choice_laddie_remapping_option_config

    ! Output
    C%do_write_laddie_output_fields                          = do_write_laddie_output_fields_config
    C%do_write_laddie_output_scalar                          = do_write_laddie_output_scalar_config
    C%time_interval_scalar_output                            = time_interval_scalar_output_config
    ! Time step
    C%dt_laddie                                              = dt_laddie_config
    C%time_duration_laddie                                   = time_duration_laddie_config
    C%time_duration_laddie_init                              = time_duration_laddie_init_config

    ! Integration
    C%choice_laddie_integration_scheme                       = choice_laddie_integration_scheme_config
    C%laddie_fbrk3_beta1                                     = laddie_fbrk3_beta1_config
    C%laddie_fbrk3_beta2                                     = laddie_fbrk3_beta2_config
    C%laddie_fbrk3_beta3                                     = laddie_fbrk3_beta3_config
    C%laddie_lfra_nu                                         = laddie_lfra_nu_config

    ! Momentum advection
    C%choice_laddie_momentum_advection                       = choice_laddie_momentum_advection_config

    ! Initialisation
    C%choice_laddie_model_initialisation                     = choice_laddie_model_initialisation_config

    ! Initialisation 'read_from_file'
    C%filename_laddie_restart                                = filename_laddie_restart_config
    C%timeframe_laddie_restart                               = timeframe_laddie_restart_config

    ! Initialisation 'uniform'
    C%laddie_initial_thickness                               = laddie_initial_thickness_config
    C%laddie_initial_T_offset                                = laddie_initial_T_offset_config
    C%laddie_initial_S_offset                                = laddie_initial_S_offset_config

    ! Equation of state
    C%choice_laddie_equation_of_state                        = choice_laddie_equation_of_state_config
    C%uniform_laddie_eos_linear_alpha                        = uniform_laddie_eos_linear_alpha_config
    C%uniform_laddie_eos_linear_beta                         = uniform_laddie_eos_linear_beta_config

    ! Coriolis
    C%choice_laddie_coriolis                                 = choice_laddie_coriolis_config
    C%uniform_laddie_coriolis_parameter                      = uniform_laddie_coriolis_parameter_config

    ! Turbulent heat exchange
    C%choice_laddie_gamma                                    = choice_laddie_gamma_config
    C%uniform_laddie_gamma_T                                 = uniform_laddie_gamma_T_config

    ! Drag coefficients
    C%laddie_drag_coefficient_top                            = laddie_drag_coefficient_top_config
    C%laddie_drag_coefficient_mom                            = laddie_drag_coefficient_mom_config

    ! Viscosity and diffusivity
    C%laddie_viscosity                                       = laddie_viscosity_config
    C%laddie_diffusivity                                     = laddie_diffusivity_config

    ! Entrainment
    C%choice_laddie_entrainment                              = choice_laddie_entrainment_config
    C%laddie_Holland2006_cl                                  = laddie_Holland2006_cl_config
    C%laddie_Gaspar1988_mu                                   = laddie_Gaspar1988_mu_config

    ! Stability
    C%laddie_thickness_minimum                               = laddie_thickness_minimum_config
    C%laddie_thickness_maximum                               = laddie_thickness_maximum_config
    C%laddie_velocity_maximum                                = laddie_velocity_maximum_config
    C%laddie_buoyancy_minimum                                = laddie_buoyancy_minimum_config

    ! Tides
    C%choice_laddie_tides                                    = choice_laddie_tides_config
    C%uniform_laddie_tidal_velocity                          = uniform_laddie_tidal_velocity_config

    ! Subglacial discharge (SGD)
    C%choice_laddie_SGD                                      = choice_laddie_SGD_config
    C%choice_laddie_SGD_idealised                            = choice_laddie_SGD_idealised_config
    C%laddie_SGD_flux                                        = laddie_SGD_flux_config
    C%filename_laddie_mask_SGD                               = filename_laddie_mask_SGD_config
    C%start_time_of_applying_SGD                             = start_time_of_applying_SGD_config
    C%transects_SGD                                          = transects_SGD_config
    C%distribute_SGD                                         = distribute_SGD_config

    ! == Lateral mass balance
    ! =======================

    ! Time step
    C%dt_LMB                                                 = dt_LMB_config

    ! Choice of LMB model
    C%choice_LMB_model_NAM                                   = choice_LMB_model_NAM_config
    C%choice_LMB_model_EAS                                   = choice_LMB_model_EAS_config
    C%choice_LMB_model_GRL                                   = choice_LMB_model_GRL_config
    C%choice_LMB_model_ANT                                   = choice_LMB_model_ANT_config

    ! "uniform"
    C%uniform_LMB                                            = uniform_LMB_config

    ! "GlacialIndex"
    C%filename_LMB_GI_NAM                                    = filename_LMB_GI_NAM_config
    C%filename_LMB_GI_EAS                                    = filename_LMB_GI_EAS_config
    C%filename_LMB_GI_GRL                                    = filename_LMB_GI_GRL_config
    C%filename_LMB_GI_ANT                                    = filename_LMB_GI_ANT_config

    C%warm_LMB_NAM                                           = warm_LMB_NAM_config
    C%warm_LMB_EAS                                           = warm_LMB_EAS_config
    C%warm_LMB_GRL                                           = warm_LMB_GRL_config
    C%warm_LMB_ANT                                           = warm_LMB_ANT_config

    C%cold_LMB_NAM                                           = cold_LMB_NAM_config
    C%cold_LMB_EAS                                           = cold_LMB_EAS_config
    C%cold_LMB_GRL                                           = cold_LMB_GRL_config
    C%cold_LMB_ANT                                           = cold_LMB_ANT_config

    ! == Glacial isostatic adjustment
    ! ===============================

    ! General settings
    C%choice_GIA_model                                       = choice_GIA_model_config
    C%dt_GIA                                                 = dt_GIA_config
    C%dx_GIA                                                 = dx_GIA_config
    C%ELRA_lithosphere_flex_rigidity                         = ELRA_lithosphere_flex_rigidity_config
    C%choice_GIA_ELRA_relaxation_time                        = choice_GIA_ELRA_relaxation_time_config
    C%ELRA_bedrock_relaxation_time                           = ELRA_bedrock_relaxation_time_config
    C%filename_GIA_ELRA_relaxation_time                      = filename_GIA_ELRA_relaxation_time_config
    C%ELRA_mantle_density                                    = ELRA_mantle_density_config

    ! == Sea level
    ! ============

    C%choice_sealevel_model                                  = choice_sealevel_model_config
    C%fixed_sealevel                                         = fixed_sealevel_config
    C%filename_prescribed_sealevel                           = filename_prescribed_sealevel_config

    ! == SELEN
    ! ========

    C%SELEN_run_at_t_start                                   = SELEN_run_at_t_start_config
    C%SELEN_n_TDOF_iterations                                = SELEN_n_TDOF_iterations_config
    C%SELEN_n_recursion_iterations                           = SELEN_n_recursion_iterations_config
    C%SELEN_use_rotational_feedback                          = SELEN_use_rotational_feedback_config
    C%SELEN_n_harmonics                                      = SELEN_n_harmonics_config
    C%SELEN_display_progress                                 = SELEN_display_progress_config

    C%SELEN_dir                                              = SELEN_dir_config
    C%SELEN_global_topo_filename                             = SELEN_global_topo_filename_config
    C%SELEN_TABOO_init_filename                              = SELEN_TABOO_init_filename_config
    C%SELEN_LMJ_VALUES_filename                              = SELEN_LMJ_VALUES_filename_config

    C%SELEN_irreg_time_n                                     = SELEN_irreg_time_n_config
    C%SELEN_irreg_time_window                                = SELEN_irreg_time_window_config

    C%SELEN_lith_thickness                                   = SELEN_lith_thickness_config
    C%SELEN_visc_n                                           = SELEN_visc_n_config
    C%SELEN_visc_prof                                        = SELEN_visc_prof_config

    ! Settings for the TABOO Earth deformation model
    C%SELEN_TABOO_CDE                                        = SELEN_TABOO_CDE_config
    C%SELEN_TABOO_TLOVE                                      = SELEN_TABOO_TLOVE_config
    C%SELEN_TABOO_DEG1                                       = SELEN_TABOO_DEG1_config
    C%SELEN_TABOO_RCMB                                       = SELEN_TABOO_RCMB_config

    ! == Tracer tracking
    ! ==================

    C%choice_tracer_tracking_model                           = choice_tracer_tracking_model_config

    ! Settings for the particle-based tracer-tracking model
    C%tractrackpart_dt_coupling                              = tractrackpart_dt_coupling_config
    C%tractrackpart_dx_particle                              = tractrackpart_dx_particle_config
    C%tractrackpart_dt_particle_min                          = tractrackpart_dt_particle_min_config
    C%tractrackpart_dt_particle_max                          = tractrackpart_dt_particle_max_config
    C%tractrackpart_n_max_particles                          = tractrackpart_n_max_particles_config
    C%tractrackpart_dt_new_particles                         = tractrackpart_dt_new_particles_config
    C%tractrackpart_dx_new_particles                         = tractrackpart_dx_new_particles_config
    C%tractrackpart_remap_n_nearest                          = tractrackpart_remap_n_nearest_config
    C%tractrackpart_write_raw_output                         = tractrackpart_write_raw_output_config
    C%tractrackpart_dt_raw_output                            = tractrackpart_dt_raw_output_config

    ! == Output
    ! =========

    ! Basic settings
    C%do_create_netcdf_output                                = do_create_netcdf_output_config
    C%output_precision                                       = output_precision_config
    C%do_compress_output                                     = do_compress_output_config
    C%dt_output                                              = dt_output_config
    C%dt_output_restart                                      = dt_output_restart_config
    C%dt_output_grid                                         = dt_output_grid_config
    C%dx_output_grid_NAM                                     = dx_output_grid_NAM_config
    C%dx_output_grid_EAS                                     = dx_output_grid_EAS_config
    C%dx_output_grid_GRL                                     = dx_output_grid_GRL_config
    C%dx_output_grid_ANT                                     = dx_output_grid_ANT_config
    C%dx_output_grid_ROI_NAM                                 = dx_output_grid_ROI_NAM_config
    C%dx_output_grid_ROI_EAS                                 = dx_output_grid_ROI_EAS_config
    C%dx_output_grid_ROI_GRL                                 = dx_output_grid_ROI_GRL_config
    C%dx_output_grid_ROI_ANT                                 = dx_output_grid_ROI_ANT_config

    ! ISMIP
    C%do_create_ismip_output                                 = do_create_ismip_output_config
    C%ismip_scenario_name                                    = ismip_scenario_name_config
    C%ismip_group_name                                       = ismip_group_name_config
    C%ismip_model_name                                       = ismip_model_name_config
    C%ismip_member_id                                        = ismip_member_id_config
    C%ismip_forcing_member_id                                = ismip_forcing_member_id_config
    C%ismip_counter                                          = ismip_counter_config
    C%ismip_esm_name                                         = ismip_esm_name_config
    C%ismip_contact_name                                     = ismip_contact_name_config
    C%ismip_contact_email                                    = ismip_contact_email_config
    C%ismip_conventions                                      = ismip_conventions_config
    C%dt_output_ismip                                        = dt_output_ismip_config

    ! Transects
    C%transects_NAM                                          = transects_NAM_config
    C%transects_EAS                                          = transects_EAS_config
    C%transects_GRL                                          = transects_GRL_config
    C%transects_ANT                                          = transects_ANT_config

    ! Which data fields we want to write to the main NetCDF output files
    C%choice_output_field_01                                 = choice_output_field_01_config
    C%choice_output_field_02                                 = choice_output_field_02_config
    C%choice_output_field_03                                 = choice_output_field_03_config
    C%choice_output_field_04                                 = choice_output_field_04_config
    C%choice_output_field_05                                 = choice_output_field_05_config
    C%choice_output_field_06                                 = choice_output_field_06_config
    C%choice_output_field_07                                 = choice_output_field_07_config
    C%choice_output_field_08                                 = choice_output_field_08_config
    C%choice_output_field_09                                 = choice_output_field_09_config
    C%choice_output_field_10                                 = choice_output_field_10_config
    C%choice_output_field_11                                 = choice_output_field_11_config
    C%choice_output_field_12                                 = choice_output_field_12_config
    C%choice_output_field_13                                 = choice_output_field_13_config
    C%choice_output_field_14                                 = choice_output_field_14_config
    C%choice_output_field_15                                 = choice_output_field_15_config
    C%choice_output_field_16                                 = choice_output_field_16_config
    C%choice_output_field_17                                 = choice_output_field_17_config
    C%choice_output_field_18                                 = choice_output_field_18_config
    C%choice_output_field_19                                 = choice_output_field_19_config
    C%choice_output_field_20                                 = choice_output_field_20_config
    C%choice_output_field_21                                 = choice_output_field_21_config
    C%choice_output_field_22                                 = choice_output_field_22_config
    C%choice_output_field_23                                 = choice_output_field_23_config
    C%choice_output_field_24                                 = choice_output_field_24_config
    C%choice_output_field_25                                 = choice_output_field_25_config
    C%choice_output_field_26                                 = choice_output_field_26_config
    C%choice_output_field_27                                 = choice_output_field_27_config
    C%choice_output_field_28                                 = choice_output_field_28_config
    C%choice_output_field_29                                 = choice_output_field_29_config
    C%choice_output_field_30                                 = choice_output_field_30_config
    C%choice_output_field_31                                 = choice_output_field_31_config
    C%choice_output_field_32                                 = choice_output_field_32_config
    C%choice_output_field_33                                 = choice_output_field_33_config
    C%choice_output_field_34                                 = choice_output_field_34_config
    C%choice_output_field_35                                 = choice_output_field_35_config
    C%choice_output_field_36                                 = choice_output_field_36_config
    C%choice_output_field_37                                 = choice_output_field_37_config
    C%choice_output_field_38                                 = choice_output_field_38_config
    C%choice_output_field_39                                 = choice_output_field_39_config
    C%choice_output_field_40                                 = choice_output_field_40_config
    C%choice_output_field_41                                 = choice_output_field_41_config
    C%choice_output_field_42                                 = choice_output_field_42_config
    C%choice_output_field_43                                 = choice_output_field_43_config
    C%choice_output_field_44                                 = choice_output_field_44_config
    C%choice_output_field_45                                 = choice_output_field_45_config
    C%choice_output_field_46                                 = choice_output_field_46_config
    C%choice_output_field_47                                 = choice_output_field_47_config
    C%choice_output_field_48                                 = choice_output_field_48_config
    C%choice_output_field_49                                 = choice_output_field_49_config
    C%choice_output_field_50                                 = choice_output_field_50_config

    ! Finished copying the values of the _config variables to the C structure

    ! Total mask values (used only for diagnostic output)
    ! ===================================================

    C%type_icefree_land                        = 1
    C%type_icefree_ocean                       = 2
    C%type_grounded_ice                        = 3
    C%type_floating_ice                        = 4
    C%type_groundingline_gr                    = 5
    C%type_groundingline_fl                    = 6
    C%type_calvingfront_gr                     = 7
    C%type_calvingfront_fl                     = 8
    C%type_margin                              = 9
    C%type_coastline                           = 10

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine copy_config_variables_to_struct

end module model_configuration_type_and_namelist
