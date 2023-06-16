MODULE unit_tests_ice

  ! Unit tests for different ice velocity solvers.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE scalar_types                                           , ONLY: type_regional_scalars
  USE reference_geometries                                   , ONLY: type_reference_geometry, initialise_reference_geometries_raw, initialise_reference_geometries_on_model_mesh
  USE mesh_creation                                          , ONLY: create_mesh_from_gridded_geometry, write_mesh_success
  USE ice_model_main                                         , ONLY: initialise_ice_model, calc_zeta_gradients
  USE ice_velocity_main                                      , ONLY: solve_stress_balance
  USE mesh_operators                                         , ONLY: calc_3D_matrix_operators_mesh
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file
  USE netcdf_output                                          , ONLY: add_zeta_dimension_to_file, setup_mesh_in_netcdf_file, add_field_mesh_dp_3D_b_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_b_notime

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_ice_unit_tests
    ! Run all unit tests for the different ice velocity solvers

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_all_ice_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all unit tests for the different ice velocity solvers
    CALL test_ice_velocities_Halfar_dome

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_ice_unit_tests

  SUBROUTINE test_ice_velocities_Halfar_dome
    ! Generate a simple Halfar dome geometry and calculate ice velocities
    ! with the SIA, DIVA, and BPA to check if the results make sense.

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ice_velocities_Halfar_dome'
    TYPE(type_mesh)                                                    :: mesh
    TYPE(type_reference_geometry)                                      :: refgeo_init, refgeo_PD, refgeo_GIAeq
    CHARACTER(LEN=256)                                                 :: region_name, mesh_name
    REAL(dp)                                                           :: lambda_M, phi_M, beta_stereo, xmin, xmax, ymin, ymax
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_SIA , v_3D_b_SIA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_DIVA, v_3D_b_DIVA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_BPA , v_3D_b_BPA
    LOGICAL                                                            :: found_errors_SIA, found_errors_DIVA, found_errors_BPA
    REAL(dp), PARAMETER                                                :: R_check         = 440E3_dp  ! At this distance from the ice divide, uabs_surf should be equal to...
    REAL(dp), PARAMETER                                                :: uabs_surf_check = 105._dp   ! ...this value...
    REAL(dp), PARAMETER                                                :: uabs_surf_tol   = 20._dp    ! ...within this tolerance.
    INTEGER                                                            :: ti
    REAL(dp)                                                           :: R, uabs_surf_target, uabs_surf_SIA, uabs_surf_DIVA, uabs_surf_BPA
    CHARACTER(LEN=256)                                                 :: filename
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

  ! Set model configuration for this experiment
  ! ===========================================

    ! Model domain
    region_name = 'ANT'
    lambda_M                          = 0._dp                            ! Longitude of the pole of the stereographic projection for the Antarctica domain [degrees east]
    phi_M                             = -90._dp                          ! Latitude  of the pole of the stereographic projection for the Antarctica domain [degrees north]
    beta_stereo                       = 71._dp                           ! Standard parallel     of the stereographic projection for the Antarctica domain [degrees]
    xmin                              = -800E3_dp                        ! Western  boundary     of the Antarctica domain [m]
    xmax                              =  800E3_dp                        ! Eastern  boundary     of the Antarctica domain [m]
    ymin                              = -800E3_dp                        ! Southern boundary     of the Antarctica domain [m]
    ymax                              =  800E3_dp                        ! Northern boundary     of the Antarctica domain [m]

    ! == Reference geometries (initial, present-day, and GIA equilibrium)

    ! Some pre-processing stuff for reference ice geometry
    C%refgeo_Hi_min                         = 2._dp                            ! Remove ice thinner than this value in the reference ice geometry. Particularly useful for BedMachine Greenland, which somehow covers the entire tundra with half a meter of ice...
    C%remove_Lake_Vostok                    = .FALSE.

    ! == Initial geometry
    C%choice_refgeo_init_ANT                = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    C%choice_refgeo_init_idealised          = 'Halfar'                         ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_init_idealised              = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry

    ! == Present-day geometry
    C%choice_refgeo_PD_ANT                  = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    C%choice_refgeo_PD_idealised            = 'Halfar'                         ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_PD_idealised                = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry

    ! == GIA equilibrium geometry
    C%choice_refgeo_GIAeq_ANT               = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    C%choice_refgeo_GIAeq_idealised         = 'Halfar'                         ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_GIAeq_idealised             = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry

    ! == Parameters for idealised geometries
    C%refgeo_idealised_Halfar_H0            = 3000._dp                         ! Suggested value: 3000 m
    C%refgeo_idealised_Halfar_R0            = 500E3_dp                         ! Suggested value: 500E3 m

    ! == Mesh generation

    ! How to set up the initial mesh
    C%choice_initial_mesh_ANT               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'

    ! Resolutions for different parts of the ice sheet
    C%maximum_resolution_uniform            = 100e3_dp                         ! [m]          Maximum resolution for the entire domain
    C%maximum_resolution_grounded_ice       = 20e3_dp                          ! [m]          Maximum resolution for grounded ice
    C%maximum_resolution_floating_ice       = 100e3_dp                         ! [m]          Maximum resolution for floating ice
    C%maximum_resolution_grounding_line     = 100e3_dp                         ! [m]          Maximum resolution for the grounding line
    C%grounding_line_width                  = 100e3_dp                         ! [m]          Width of the band around the grounding line that should get this resolution
    C%maximum_resolution_calving_front      = 100e3_dp                         ! [m]          Maximum resolution for the calving front
    C%calving_front_width                   = 100e3_dp                         ! [m]          Width of the band around the calving front that should get this resolution
    C%maximum_resolution_ice_front          = 10e3_dp                          ! [m]          Maximum resolution for the ice front
    C%ice_front_width                       = 5e3_dp                           ! [m]          Width of the band around the ice front that should get this resolution
    C%maximum_resolution_coastline          = 100e3_dp                         ! [m]          Maximum resolution for the coastline
    C%coastline_width                       = 100e3_dp                         ! [m]          Width of the band around the coastline that should get this resolution

    ! Regions of interest
    C%choice_regions_of_interest            = ''                               ! Regions of interest where other (higher) resolutions apply. Separated by double vertical bars "||", e.g. "PineIsland||Thwaites"

    ! Advanced geometry parameters
    C%do_singlecore_mesh_creation           = .TRUE.                           !              Whether or not to use only a single core for mesh generation (for better reproducibility)
    C%alpha_min                             = 0.4363_dp                        ! [radians]    Smallest allowed internal triangle angle
    C%nit_Lloyds_algorithm                  = 3                                ! [-]          Number of iterations of Lloyds algorithm to be applied after refinement
    C%mesh_resolution_tolerance             = 1.25_dp                          ! [-]          Factors the target resolution for trangle-size requirement. 1=strict, use >1 to avoid unnecesarily high resolution

    ! Memory
    C%nC_mem                                = 32                               ! [-]          How many columns of memory should be allocated for connectivity lists

    ! == The scaled vertical coordinate zeta

    C%choice_zeta_grid                      = 'regular'                        ! The type of vertical grid to use; can be "regular", "irregular_log", "old_15_layer_zeta"
    C%nz                                    = 12                               ! The number of vertical layers to use

    ! == Ice dynamics - velocity

    ! General
    C%choice_stress_balance_approximation   = 'SIA'                            ! Choice of stress balance approximation: "none" (= no flow, though geometry can still change due to mass balance), "SIA", "SSA", "SIA/SSA", "DIVA", "BPA"
    C%n_flow                                = 3._dp                            ! Exponent in Glen's flow law
    C%m_enh_sheet                           = 1._dp                            ! Ice flow enhancement factor for grounded ice
    C%m_enh_shelf                           = 1._dp                            ! Ice flow enhancement factor for floating ice
    C%choice_hybrid_SIASSA_scheme           = 'add'                            ! Choice of scheme for combining SIA and SSA velocities in the hybrid approach
    C%do_GL_subgrid_friction                = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
    C%subgrid_friction_exponent             = 2._dp                            ! Exponent to which f_grnd should be raised before being used to scale beta
    C%do_include_SSADIVA_crossterms         = .TRUE.                           ! Whether or not to include the gradients of the effective viscosity (the "cross-terms") in the solution of the SSA/DIVA

    ! Initialisation
    C%choice_initial_velocity_ANT           = 'zero'

    ! Some parameters for numerically solving the stress balance
    C%SIA_maximum_diffusivity               = 1E5_dp                           ! Limit the diffusivity in the SIA to this value
    C%visc_it_norm_dUV_tol                  = 5E-6_dp                          ! Stop criterion for the viscosity iteration: the L2-norm of successive velocity solutions should be smaller than this number
    C%visc_it_nit                           = 500._dp                          ! Maximum number of effective viscosity iterations
    C%visc_it_relax                         = 0.4_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    C%epsilon_sq_0                          = 1E-10_dp                         ! Normalisation term so that zero velocity gives non-zero viscosity
    C%visc_eff_min                          = 1E0_dp                           ! Minimum value for effective viscosity
    C%beta_max                              = 1E20_dp                          ! Maximum value for basal friction coefficient
    C%vel_max                               = 5000._dp                         ! Velocities are limited to this value
    C%stress_balance_PETSc_rtol             = 1E-6_dp                          ! PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    C%stress_balance_PETSc_abstol           = 1E-4_dp                          ! PETSc solver - stop criterion, absolute difference

    ! == Ice dynamics - sliding law

    ! Sliding laws
    C%choice_sliding_law                    = 'no_sliding'                     ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Budd", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"

    ! == Ice dynamics - boundary conditions

    C%BC_u_west                             = 'zero'                           ! Boundary conditions for the ice velocity field at the domain border
    C%BC_u_east                             = 'zero'                           ! Allowed choices: "infinite", "zero", "zero"
    C%BC_u_south                            = 'zero'
    C%BC_u_north                            = 'zero'
    C%BC_v_west                             = 'zero'
    C%BC_v_east                             = 'zero'
    C%BC_v_south                            = 'zero'
    C%BC_v_north                            = 'zero'
    C%BC_H_west                             = 'zero'                           ! Boundary conditions for ice thickness at the domain boundary
    C%BC_H_east                             = 'zero'                           ! Allowed choices:  "infinite", "zero", "ISMIP_HOM_F"
    C%BC_H_south                            = 'zero'
    C%BC_H_north                            = 'zero'

    ! Rheological model (relating Glen's flow parameter to ice temperature)
    C%choice_ice_rheology                   = 'uniform'                        ! Choice of ice rheology model: "uniform", "Huybrechts1992", "MISMIP_mod"
    C%uniform_flow_factor                   = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")

  ! Initialise the model
  ! ====================

    ! Initialise all the reference geometries on their raw input grids
    CALL initialise_reference_geometries_raw( region_name, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Create mesh from gridded initial geometry data
    mesh_name = 'Halfar_mesh'
    CALL create_mesh_from_gridded_geometry( region_name, mesh_name, &
      refgeo_init%grid_raw, &
      refgeo_init%Hi_grid_raw, &
      refgeo_init%Hb_grid_raw, &
      refgeo_init%Hs_grid_raw, &
      refgeo_init%SL_grid_raw, &
      xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, &
      mesh)

    ! Write the mesh creation success message to the terminal
    CALL write_mesh_success( mesh)

    ! Remap reference geometries from their raw input grids to the model mesh
    CALL initialise_reference_geometries_on_model_mesh( region_name, mesh, refgeo_init, refgeo_PD, refgeo_GIAeq)

  ! Calculate ice velocities
  ! ========================

    ! Allocate memory
    ALLOCATE( u_3D_b_SIA(  mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_SIA(  mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_DIVA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_DIVA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_BPA(  mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_BPA(  mesh%ti1:mesh%ti2, mesh%nz))

    ! Calculate velocities
    CALL test_ice_velocities_Halfar_dome_SIA(  region_name, mesh, refgeo_init, refgeo_PD, u_3D_b_SIA , v_3D_b_SIA )
    CALL test_ice_velocities_Halfar_dome_DIVA( region_name, mesh, refgeo_init, refgeo_PD, u_3D_b_DIVA, v_3D_b_DIVA)
    CALL test_ice_velocities_Halfar_dome_BPA(  region_name, mesh, refgeo_init, refgeo_PD, u_3D_b_BPA , v_3D_b_BPA )

    ! Validate results - check if the surface velocity follows the expected
    ! linear increase away from the ice divide

    found_errors_SIA  = .FALSE.
    found_errors_DIVA = .FALSE.
    found_errors_BPA  = .FALSE.

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate distance from the ice divide
      R = NORM2( mesh%TriGC( ti,:))

      ! Only check up to a certain distance (near the margin, the SIA and BPA become a bit "noisy")
      IF (R > R_check) CYCLE

      ! Calculate target surface velocity
      uabs_surf_target = uabs_surf_check * R / R_check

      ! Calculate actual surface velocities
      uabs_surf_SIA  = SQRT( u_3D_b_SIA(  ti,1)**2 + v_3D_b_SIA(  ti,1)**2)
      uabs_surf_DIVA = SQRT( u_3D_b_DIVA( ti,1)**2 + v_3D_b_DIVA( ti,1)**2)
      uabs_surf_BPA  = SQRT( u_3D_b_BPA(  ti,1)**2 + v_3D_b_BPA(  ti,1)**2)

      ! Compare
      IF (ABS( uabs_surf_SIA  - uabs_surf_target) > uabs_surf_tol) found_errors_SIA  = .TRUE.
      IF (ABS( uabs_surf_DIVA - uabs_surf_target) > uabs_surf_tol) found_errors_DIVA = .TRUE.
      IF (ABS( uabs_surf_BPA  - uabs_surf_target) > uabs_surf_tol) found_errors_BPA  = .TRUE.

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_SIA , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_DIVA, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_BPA , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

    IF (.NOT. (found_errors_SIA .OR. found_errors_DIVA .OR. found_errors_BPA)) THEN
      IF (par%master) CALL happy('validated SIA, DIVA, and BPA solvers for diagnostic velocities in the Halfar dome geometry')
    END IF
    IF (found_errors_SIA) THEN
      IF (par%master) CALL warning('found errors in SIA solver for diagnostic velocities in the Halfar dome geometry')
    END IF
    IF (found_errors_DIVA) THEN
      IF (par%master) CALL warning('found errors in DIVA solver for diagnostic velocities in the Halfar dome geometry')
    END IF
    IF (found_errors_BPA) THEN
      IF (par%master) CALL warning('found errors in BPA solver for diagnostic velocities in the Halfar dome geometry')
    END IF

  ! Write results to NetCDF
  ! =======================

    ! Create a NetCDF output file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add a zeta dimension for the 3-D ice velocities
    CALL add_zeta_dimension_to_file( filename, ncid, mesh%zeta)

    ! Add all the ice velocities as fields
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_SIA' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_SIA' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_DIVA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_DIVA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_BPA' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_BPA' )

    ! Write the velocities to the file
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_SIA' , u_3D_b_SIA )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_SIA' , v_3D_b_SIA )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_DIVA', u_3D_b_DIVA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_DIVA', v_3D_b_DIVA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_BPA' , u_3D_b_BPA )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_BPA' , v_3D_b_BPA )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    DEALLOCATE( u_3D_b_SIA)
    DEALLOCATE( v_3D_b_SIA)
    DEALLOCATE( u_3D_b_DIVA)
    DEALLOCATE( v_3D_b_DIVA)
    DEALLOCATE( u_3D_b_BPA)
    DEALLOCATE( v_3D_b_BPA)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ice_velocities_Halfar_dome

  SUBROUTINE test_ice_velocities_Halfar_dome_SIA( region_name, mesh, refgeo_init, refgeo_PD, u_3D_b, v_3D_b)
    ! Generate a simple Halfar dome geometry and calculate ice velocities with the SIA

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),                                    INTENT(IN)    :: region_name
    TYPE(type_mesh),                                     INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),                       INTENT(IN)    :: refgeo_init, refgeo_PD
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,mesh%nz),      INTENT(OUT)   :: u_3D_b, v_3D_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ice_velocities_Halfar_dome_SIA'
    TYPE(type_ice_model)                                               :: ice
    TYPE(type_regional_scalars)                                        :: scalars

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set stres balance to SIA
    C%choice_stress_balance_approximation   = 'SIA'

    ! Initialise
    CALL initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ice%A_flow_3D = C%uniform_flow_factor
    ALLOCATE( ice%beta_b( mesh%vi1:mesh%vi2))

    ! Calculate velocities
    CALL solve_stress_balance( mesh, ice)

    ! Copy result
    u_3D_b = ice%u_3D_b
    v_3D_b = ice%v_3D_b

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ice_velocities_Halfar_dome_SIA

  SUBROUTINE test_ice_velocities_Halfar_dome_DIVA( region_name, mesh, refgeo_init, refgeo_PD, u_3D_b, v_3D_b)
    ! Generate a simple Halfar dome geometry and calculate ice velocities with the DIVA

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),                                    INTENT(IN)    :: region_name
    TYPE(type_mesh),                                     INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),                       INTENT(IN)    :: refgeo_init, refgeo_PD
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,mesh%nz),      INTENT(OUT)   :: u_3D_b, v_3D_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ice_velocities_Halfar_dome_DIVA'
    TYPE(type_ice_model)                                               :: ice
    TYPE(type_regional_scalars)                                        :: scalars

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set stres balance to SIA
    C%choice_stress_balance_approximation   = 'DIVA'

    ! Initialise
    CALL initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ice%A_flow_3D = C%uniform_flow_factor
    ALLOCATE( ice%beta_b( mesh%vi1:mesh%vi2))

    ! Calculate velocities
    CALL solve_stress_balance( mesh, ice)

    ! Copy result
    u_3D_b = ice%u_3D_b
    v_3D_b = ice%v_3D_b

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ice_velocities_Halfar_dome_DIVA

  SUBROUTINE test_ice_velocities_Halfar_dome_BPA( region_name, mesh, refgeo_init, refgeo_PD, u_3D_b, v_3D_b)
    ! Generate a simple Halfar dome geometry and calculate ice velocities with the BPA

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),                                    INTENT(IN)    :: region_name
    TYPE(type_mesh),                                     INTENT(INOUT) :: mesh
    TYPE(type_reference_geometry),                       INTENT(IN)    :: refgeo_init, refgeo_PD
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,mesh%nz),      INTENT(OUT)   :: u_3D_b, v_3D_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ice_velocities_Halfar_dome_BPA'
    TYPE(type_ice_model)                                               :: ice
    TYPE(type_regional_scalars)                                        :: scalars

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set stres balance to BPA
    C%choice_stress_balance_approximation   = 'BPA'

    ! Initialise
    CALL initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ice%A_flow_3D = C%uniform_flow_factor
    ALLOCATE( ice%beta_b( mesh%vi1:mesh%vi2))

    ! Calculate necessary 3-D matrix operators
    CALL calc_zeta_gradients( mesh, ice)
    CALL calc_3D_matrix_operators_mesh( mesh, ice)

    ! Calculate velocities
    CALL solve_stress_balance( mesh, ice)

    ! Copy result
    u_3D_b = ice%u_3D_b
    v_3D_b = ice%v_3D_b

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ice_velocities_Halfar_dome_BPA

END MODULE unit_tests_ice
