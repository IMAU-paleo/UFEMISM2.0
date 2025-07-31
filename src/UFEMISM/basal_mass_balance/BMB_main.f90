MODULE BMB_main

  ! The main BMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE laddie_forcing_types                                   , ONLY: type_laddie_forcing
  USE BMB_idealised                                          , ONLY: initialise_BMB_model_idealised, run_BMB_model_idealised
  USE BMB_prescribed                                         , ONLY: initialise_BMB_model_prescribed, run_BMB_model_prescribed
  USE BMB_parameterised                                      , ONLY: initialise_BMB_model_parameterised, run_BMB_model_parameterised
  USE BMB_laddie                                             , ONLY: initialise_BMB_model_laddie, run_BMB_model_laddie, remap_BMB_model_laddie
  use BMB_inverted, only: initialise_BMB_model_inverted, run_BMB_model_inverted
  use laddie_main, only: initialise_laddie_model, run_laddie_model, remap_laddie_model
  use laddie_utilities, only: allocate_laddie_forcing
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use ice_geometry_basics, only: is_floating
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
  use netcdf_io_main
  use checksum_mod, only: checksum
  use mpi_distributed_shared_memory, only: reallocate_dist_shared

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model( mesh, ice, ocean, refgeo, SMB, BMB, region_name, time, is_initial)
    ! Calculate the basal mass balance

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo
    TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time
    logical,                                intent(in)    :: is_initial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model
    CHARACTER(LEN=256)                                    :: choice_BMB_model_ROI
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which BMB model to run for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model      = C%choice_BMB_model_NAM
        choice_BMB_model_ROI  = C%choice_BMB_model_NAM_ROI
      CASE ('EAS')
        choice_BMB_model      = C%choice_BMB_model_EAS
        choice_BMB_model_ROI  = C%choice_BMB_model_EAS_ROI
      CASE ('GRL')
        choice_BMB_model      = C%choice_BMB_model_GRL
        choice_BMB_model_ROI  = C%choice_BMB_model_GRL_ROI
      CASE ('ANT')
        choice_BMB_model      = C%choice_BMB_model_ANT
        choice_BMB_model_ROI  = C%choice_BMB_model_ANT_ROI
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Check if we need to calculate a new BMB
    IF (C%do_asynchronous_BMB) THEN
      ! Asynchronous coupling: do not calculate a new BMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next BMB time step
      IF (time == BMB%t_next) THEN
        ! Go on to calculate a new BMB
        BMB%t_next = time + C%dt_BMB
      ELSEIF (time > BMB%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the BMB time step')
      ELSE
        ! It is not yet time to calculate a new BMB

        ! Apply subgrid scheme of old BMB to new mask
        SELECT CASE (choice_BMB_model)
          CASE ('inverted')
            ! No need to do anything
          CASE ('prescribed_fixed')
            ! No need to do anything
          CASE ('laddie')
            ! No need to do anything
          CASE DEFAULT
            CALL apply_BMB_subgrid_scheme( mesh, ice, BMB)
        END SELECT

        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_BMB) THEN
      ! Synchronous coupling: calculate a new BMB in every model loop
      BMB%t_next = time + C%dt_BMB
    END IF

    ! Run the chosen BMB model
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        BMB%BMB_shelf = 0._dp
        DO vi = mesh%vi1, mesh%vi2
          IF (ice%mask_floating_ice( vi) .OR. ice%mask_icefree_ocean( vi) .OR. ice%mask_gl_gr( vi)) THEN
            BMB%BMB_shelf( vi) = C%uniform_BMB
          END IF
        END DO
      CASE ('prescribed')
        CALL run_BMB_model_prescribed( mesh, ice, BMB, region_name, time)
      CASE ('prescribed_fixed')
        ! No need to do anything
      CASE ('idealised')
        CALL run_BMB_model_idealised( mesh, ice, BMB, time)
      CASE ('parameterised')
        CALL run_BMB_model_parameterised( mesh, ice, ocean, BMB)
      CASE ('inverted')
        CALL run_BMB_model_inverted( mesh, ice, BMB%inv, time)
        BMB%BMB = BMB%inv%BMB
      CASE ('laddie_py')
        CALL run_BMB_model_laddie( mesh, ice, BMB, time, .FALSE.)
      CASE ('laddie')
        call update_laddie_forcing( mesh, ice, ocean, BMB%forcing)
        CALL run_laddie_model( mesh, ice, ocean, BMB%laddie, BMB%forcing, time, is_initial, region_name)

        DO vi = mesh%vi1, mesh%vi2
          BMB%BMB( vi) = -BMB%laddie%melt( vi) * sec_per_year
        END DO
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Check hybrid_ROI_BMB
    SELECT CASE (choice_BMB_model_ROI)
      CASE ('identical_to_choice_BMB_model')
        ! No need to do anything
      CASE ('uniform')
        ! Update BMB only for cells in ROI
        DO vi = mesh%vi1, mesh%vi2
          IF (ice%mask_ROI(vi)) THEN
            IF (ice%mask_floating_ice( vi) .OR. ice%mask_icefree_ocean( vi) .OR. ice%mask_gl_gr( vi)) THEN
              BMB%BMB_shelf( vi) = C%uniform_BMB_ROI
            END IF
          END IF
        END DO
        CALL apply_BMB_subgrid_scheme_ROI( mesh, ice, BMB)
      CASE ('laddie_py')
        ! run_BMB_model_laddie and read BMB values only for region of interest
        CALL run_BMB_model_laddie( mesh, ice, BMB, time, .TRUE.)
        CALL apply_BMB_subgrid_scheme_ROI( mesh, ice, BMB)
      CASE ('prescribed', 'prescribed_fixed', 'idealised', 'parameterised', 'inverted', 'laddie')
        CALL crash('this BMB_model "' // TRIM( choice_BMB_model_ROI) // '" is not implemented for hybrid-BMB in ROI yet')
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_ROI "' // TRIM( choice_BMB_model_ROI) // '"')
    END SELECT

    ! Apply subgrid scheme of old BMB to new mask
    SELECT CASE (choice_BMB_model)
      CASE ('inverted')
        ! No need to do anything
      CASE ('prescribed_fixed')
        ! No need to do anything
      CASE ('laddie')
        ! No need to do anything
      CASE DEFAULT
        CALL apply_BMB_subgrid_scheme( mesh, ice, BMB)
    END SELECT

    ! save BMB in BMB_modelled if applying transition phase
    IF (C%do_BMB_transition_phase) THEN
      DO vi = mesh%vi1, mesh%vi2
        BMB%BMB_modelled( vi) = BMB%BMB( vi)
      END DO
    END IF

    ! Apply limits
    BMB%BMB = max( -C%BMB_maximum_allowed_melt_rate, min( C%BMB_maximum_allowed_refreezing_rate, BMB%BMB ))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model

  SUBROUTINE initialise_BMB_model( mesh, ice, ocean, BMB, refgeo_PD, refgeo_init, region_name)
    ! Initialise the BMB model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                   INTENT(OUT)   :: BMB
    type(type_reference_geometry),          intent(in   ) :: refgeo_PD, refgeo_init
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model
    CHARACTER(LEN=256)                                    :: choice_BMB_model_ROI

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising basal mass balance model...'

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model      = C%choice_BMB_model_NAM
        choice_BMB_model_ROI  = C%choice_BMB_model_NAM_ROI
      CASE ('EAS')
        choice_BMB_model      = C%choice_BMB_model_EAS
        choice_BMB_model_ROI  = C%choice_BMB_model_EAS_ROI
      CASE ('GRL')
        choice_BMB_model      = C%choice_BMB_model_GRL
        choice_BMB_model_ROI  = C%choice_BMB_model_GRL_ROI
      CASE ('ANT')
        choice_BMB_model      = C%choice_BMB_model_ANT
        choice_BMB_model_ROI  = C%choice_BMB_model_ANT_ROI
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Allocate memory for main variables
    ALLOCATE( BMB%BMB( mesh%vi1:mesh%vi2))
    BMB%BMB = 0._dp

    ! Allocate shelf BMB
    ALLOCATE( BMB%BMB_shelf( mesh%vi1:mesh%vi2))
    BMB%BMB_shelf = 0._dp

    ! Allocate inverted BMB
    ALLOCATE( BMB%BMB_inv( mesh%vi1:mesh%vi2))
    BMB%BMB_inv = 0._dp

    ! Allocate reference BMB
    ALLOCATE( BMB%BMB_ref( mesh%vi1:mesh%vi2))
    BMB%BMB_ref = 0._dp

    ! Allocate transition phase BMB
    ALLOCATE( BMB%BMB_transition_phase( mesh%vi1:mesh%vi2))
    BMB%BMB_transition_phase = 0._dp

    ! Allocate modelled BMB
    ALLOCATE( BMB%BMB_modelled( mesh%vi1:mesh%vi2))
    BMB%BMB_modelled = 0._dp

    ! Allocate mask for cavities
    ALLOCATE( BMB%mask_floating_ice( mesh%vi1:mesh%vi2))
    ALLOCATE( BMB%mask_gl_fl( mesh%vi1:mesh%vi2))
    ALLOCATE( BMB%mask_gl_gr( mesh%vi1:mesh%vi2))
    BMB%mask_floating_ice = ice%mask_floating_ice .AND. .NOT. ice%mask_gl_fl
    BMB%mask_gl_fl = ice%mask_gl_fl
    BMB%mask_gl_gr = ice%mask_gl_gr

    ! Set time of next calculation to start time
    BMB%t_next = C%start_time_of_run

    ! Determine which BMB model to initialise
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('prescribed')
        CALL initialise_BMB_model_prescribed( mesh, BMB, region_name)
      CASE ('prescribed_fixed')
        CALL initialise_BMB_model_prescribed( mesh, BMB, region_name)
        CALL apply_BMB_subgrid_scheme( mesh, ice, BMB)
      CASE ('idealised')
        CALL initialise_BMB_model_idealised( mesh, BMB)
      CASE ('parameterised')
        CALL initialise_BMB_model_parameterised( mesh, BMB)
      CASE ('inverted')
        call initialise_BMB_model_inverted( mesh, BMB%inv, refgeo_PD, refgeo_init)
      CASE ('laddie_py')
        CALL initialise_BMB_model_laddie( mesh, BMB)
      CASE ('laddie')
        call allocate_laddie_forcing( mesh, BMB%forcing)
        call update_laddie_forcing( mesh, ice, ocean, BMB%forcing)
        CALL initialise_laddie_model( mesh, BMB%laddie, BMB%forcing, ocean, ice, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Check hybrid_ROI_BMB
    SELECT CASE (choice_BMB_model_ROI)
      CASE ('identical_to_choice_BMB_model')
        ! No need to do anything
      CASE ('uniform')
        ! No need to do anything
      CASE ('laddie_py')
         CALL initialise_BMB_model_laddie( mesh, BMB)
      CASE ('prescribed', 'prescribed_fixed', 'idealised', 'parameterised', 'inverted', 'laddie')
        CALL crash('this BMB_model "' // TRIM( choice_BMB_model_ROI) // '" is not implemented for hybrid-BMB in ROI yet')
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_ROI "' // TRIM( choice_BMB_model_ROI) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model

  SUBROUTINE write_to_restart_file_BMB_model( mesh, BMB, region_name, time)
    ! Write to the restart file for the BMB model

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(IN)    :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Write to the restart file of the chosen BMB model
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE ('prescribed_fixed')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        ! No need to do anything
      CASE ('inverted')
        CALL write_to_restart_file_BMB_model_region( mesh, BMB, region_name, time)
      CASE ('laddie_py')
        ! No need to do anything
      CASE ('laddie')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_BMB_model

  SUBROUTINE write_to_restart_file_BMB_model_region( mesh, BMB, region_name, time)
    ! Write to the restart NetCDF file for the BMB model

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN) :: mesh
    TYPE(type_BMB_model),     INTENT(IN) :: BMB
    CHARACTER(LEN=3),         INTENT(IN) :: region_name
    REAL(dp),                 INTENT(IN) :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER        :: routine_name = 'write_to_restart_file_BMB_model_region'
    INTEGER                              :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Writing to BMB restart file "' // &
      colour_string( TRIM( BMB%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( BMB%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( BMB%restart_filename, ncid, time)

    ! ! Write the BMB fields to the file
    CALL write_to_field_multopt_mesh_dp_2D( mesh, BMB%restart_filename, ncid, 'BMB', BMB%BMB)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_BMB_model_region

  SUBROUTINE create_restart_file_BMB_model( mesh, BMB, region_name)
    ! Create the restart file for the BMB model

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Create the restart file of the chosen BMB model
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE ('prescribed_fixed')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        ! No need to do anything
      CASE ('inverted')
        CALL create_restart_file_BMB_model_region( mesh, BMB, region_name)
      CASE ('laddie_py')
        ! No need to do anything
      CASE ('laddie')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_BMB_model

  SUBROUTINE create_restart_file_BMB_model_region( mesh, BMB, region_name)
    ! Create a restart NetCDF file for the BMB submodel
    ! Includes generation of the procedural filename (e.g. "restart_BMB_00001.nc")

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)    :: mesh
    TYPE(type_BMB_model),     INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),         INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER           :: routine_name = 'create_restart_file_BMB_model_region'
    CHARACTER(LEN=256)                      :: filename_base
    INTEGER                                 :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_BMB_' // region_name
    CALL generate_filename_XXXXXdotnc( filename_base, BMB%restart_filename)

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Creating BMB model restart file "' // &
      colour_string( TRIM( BMB%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( BMB%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( BMB%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( BMB%restart_filename, ncid)

    ! Add the data fields to the file
    CALL add_field_mesh_dp_2D( BMB%restart_filename, ncid, 'BMB', long_name = 'Basal mass balance', units = 'm/yr')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_BMB_model_region

  SUBROUTINE remap_BMB_model( mesh_old, mesh_new, ice, ocean, BMB, region_name, time)
    ! Remap the BMB model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping basal mass balance model data to the new mesh...'

    ! Determine which BMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_BMB_model = C%choice_BMB_model_NAM
      CASE ('EAS')
        choice_BMB_model = C%choice_BMB_model_EAS
      CASE ('GRL')
        choice_BMB_model = C%choice_BMB_model_GRL
      CASE ('ANT')
        choice_BMB_model = C%choice_BMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Reallocate memory for main variables
    CALL reallocate_bounds( BMB%BMB, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( BMB%BMB_shelf, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( BMB%BMB_inv, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( BMB%BMB_ref, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( BMB%BMB_transition_phase, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( BMB%BMB_modelled, mesh_new%vi1, mesh_new%vi2)

    ! Reallocate memory for cavity mask
    CALL reallocate_bounds( BMB%mask_floating_ice, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( BMB%mask_gl_fl, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( BMB%mask_gl_gr, mesh_new%vi1, mesh_new%vi2)

    ! Re-initialise
    BMB%mask_floating_ice = ice%mask_floating_ice .AND. .NOT. ice%mask_gl_fl
    BMB%mask_gl_fl = ice%mask_gl_fl
    BMB%mask_gl_gr = ice%mask_gl_gr

    ! Determine which BMB model to initialise
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('prescribed')
        CALL initialise_BMB_model_prescribed( mesh_new, BMB, region_name)
      CASE ('prescribed_fixed')
        CALL initialise_BMB_model_prescribed( mesh_new, BMB, region_name)
        CALL apply_BMB_subgrid_scheme( mesh_new, ice, BMB)
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        ! we only need to run the BMB model again, considering the ocean model is remapped just before a call to this function
        CALL run_BMB_model_parameterised( mesh_new, ice, ocean, BMB)
      CASE ('inverted')
        ! No need to do anything
      CASE ('laddie_py')
        CALL remap_BMB_model_laddie( mesh_new, BMB)
      CASE ('laddie')
        call remap_laddie_forcing( mesh_old, mesh_new, BMB%forcing)
        call update_laddie_forcing( mesh_new, ice, ocean, BMB%forcing)
        CALL remap_laddie_model( mesh_old, mesh_new, ice, ocean, BMB%laddie, BMB%forcing, time, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BMB_model

! ===== Utilities =====
! =====================

  SUBROUTINE apply_BMB_subgrid_scheme( mesh, ice, BMB)
    ! Apply selected scheme for sub-grid shelf melt
    ! (see Leguy et al. 2021 for explanations of the three schemes)

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_BMB_subgrid_scheme'
    CHARACTER(LEN=256)                                    :: choice_BMB_subgrid
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Note: apply extrapolation_FCMP_to_PMP to non-laddie BMB models before applying sub-grid schemes

    BMB%BMB = 0._dp

    DO vi = mesh%vi1, mesh%vi2
      ! Different sub-grid schemes for sub-shelf melt
      CALL compute_subgrid_BMB(ice, BMB, vi)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_BMB_subgrid_scheme

  SUBROUTINE apply_BMB_subgrid_scheme_ROI( mesh, ice, BMB)
    ! Apply selected scheme for sub-grid shelf melt
    ! (see Leguy et al. 2021 for explanations of the three schemes)

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_BMB_subgrid_scheme_ROI'
    CHARACTER(LEN=256)                                    :: choice_BMB_subgrid
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Note: apply extrapolation_FCMP_to_PMP to non-laddie BMB models before applying sub-grid schemes

    DO vi = mesh%vi1, mesh%vi2
      ! Only for ROI cells
      IF (ice%mask_ROI(vi)) THEN
        CALL compute_subgrid_BMB(ice, BMB, vi)
      END IF
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_BMB_subgrid_scheme_ROI

  SUBROUTINE compute_subgrid_BMB(ice, BMB, vi)

    INTEGER                     , INTENT(IN)              :: vi
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Different sub-grid schemes for sub-shelf melt
        IF (C%do_subgrid_BMB_at_grounding_line) THEN
          IF     (C%choice_BMB_subgrid == 'FCMP') THEN
            ! Apply FCMP scheme
            IF (ice%mask_floating_ice( vi)) BMB%BMB( vi) = BMB%BMB_shelf( vi)

          ELSEIF (C%choice_BMB_subgrid == 'PMP') THEN
            ! Apply PMP scheme
            IF (ice%mask_floating_ice( vi) .OR. ice%mask_gl_gr( vi)) BMB%BMB( vi) = (1._dp - ice%fraction_gr( vi)) * BMB%BMB_shelf( vi)
          ELSE
            CALL crash('unknown choice_BMB_subgrid "' // TRIM(C%choice_BMB_subgrid) // '"!')
          END IF
        ELSE
          ! Apply NMP scheme
          IF (ice%fraction_gr( vi) == 0._dp) BMB%BMB( vi) = BMB%BMB_shelf( vi)
        END IF

  END SUBROUTINE compute_subgrid_BMB

  subroutine update_laddie_forcing( mesh, ice, ocean, forcing)

    ! In/output variables
    type(type_mesh),           intent(in   ) :: mesh 
    type(type_ice_model),      intent(in   ) :: ice
    type(type_ocean_model),    intent(in   ) :: ocean
    type(type_laddie_forcing), intent(inout) :: forcing

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_laddie_forcing'

    ! Add routine to path
    call init_routine( routine_name)

    forcing%Hi                ( mesh%vi1:mesh%vi2  ) = ice%Hi                ( mesh%vi1:mesh%vi2  )
    forcing%Hib               ( mesh%vi1:mesh%vi2  ) = ice%Hib               ( mesh%vi1:mesh%vi2  )
    forcing%dHib_dx_b         ( mesh%ti1:mesh%ti2  ) = ice%dHib_dx_b         ( mesh%ti1:mesh%ti2  )
    forcing%dHib_dy_b         ( mesh%ti1:mesh%ti2  ) = ice%dHib_dy_b         ( mesh%ti1:mesh%ti2  )
    forcing%mask_icefree_land ( mesh%vi1:mesh%vi2  ) = ice%mask_icefree_land ( mesh%vi1:mesh%vi2  )
    forcing%mask_icefree_ocean( mesh%vi1:mesh%vi2  ) = ice%mask_icefree_ocean( mesh%vi1:mesh%vi2  )
    forcing%mask_grounded_ice ( mesh%vi1:mesh%vi2  ) = ice%mask_grounded_ice ( mesh%vi1:mesh%vi2  )
    forcing%mask_floating_ice ( mesh%vi1:mesh%vi2  ) = ice%mask_floating_ice ( mesh%vi1:mesh%vi2  )

    forcing%mask_gl_fl        ( mesh%vi1:mesh%vi2  ) = ice%mask_gl_fl        ( mesh%vi1:mesh%vi2  )
    forcing%mask_SGD          ( mesh%vi1:mesh%vi2  ) = ice%mask_SGD          ( mesh%vi1:mesh%vi2  )
    
    forcing%Ti                ( mesh%vi1:mesh%vi2,:) = ice%Ti                ( mesh%vi1:mesh%vi2,:) - 273.15 ! [degC]
    forcing%T_ocean           ( mesh%vi1:mesh%vi2,:) = ocean%T               ( mesh%vi1:mesh%vi2,:)
    forcing%S_ocean           ( mesh%vi1:mesh%vi2,:) = ocean%S               ( mesh%vi1:mesh%vi2,:)

    call checksum( forcing%Hi                , 'forcing%Hi'                , mesh%pai_V)
    call checksum( forcing%Hib               , 'forcing%Hib'               , mesh%pai_V)
    call checksum( forcing%dHib_dx_b         , 'forcing%dHib_dx_b'         , mesh%pai_Tri)
    call checksum( forcing%dHib_dy_b         , 'forcing%dHib_dy_b'         , mesh%pai_Tri)
    call checksum( forcing%mask_icefree_land , 'forcing%mask_icefree_land' , mesh%pai_V)
    call checksum( forcing%mask_icefree_ocean, 'forcing%mask_icefree_ocean', mesh%pai_V)
    call checksum( forcing%mask_grounded_ice , 'forcing%mask_grounded_ice' , mesh%pai_V)
    call checksum( forcing%mask_floating_ice , 'forcing%mask_floating_ice' , mesh%pai_V)
    call checksum( forcing%mask_gl_fl        , 'forcing%mask_gl_fl'        , mesh%pai_V)
    call checksum( forcing%mask_SGD          , 'forcing%mask_SGD'          , mesh%pai_V)
    call checksum( forcing%Ti                , 'forcing%Ti'                , mesh%pai_V)
    call checksum( forcing%T_ocean           , 'forcing%T_ocean'           , mesh%pai_V)
    call checksum( forcing%S_ocean           , 'forcing%S_ocean'           , mesh%pai_V)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_laddie_forcing

  subroutine remap_laddie_forcing( mesh_old, mesh_new, forcing)
    ! Reallocate and remap laddie forcing

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh_old
    type(type_mesh),                        intent(in)    :: mesh_new
    type(type_laddie_forcing),              intent(inout) :: forcing

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'remap_laddie_forcing'

    ! Add routine to path
    call init_routine( routine_name)

    ! Forcing
    call reallocate_dist_shared( forcing%Hi                , forcing%wHi                , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%Hib               , forcing%wHib               , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%dHib_dx_b         , forcing%wdHib_dx_b         , mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( forcing%dHib_dy_b         , forcing%wdHib_dy_b         , mesh_new%pai_Tri%n_nih)
    call reallocate_dist_shared( forcing%mask_icefree_land , forcing%wmask_icefree_land , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%mask_icefree_ocean, forcing%wmask_icefree_ocean, mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%mask_grounded_ice , forcing%wmask_grounded_ice , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%mask_floating_ice , forcing%wmask_floating_ice , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%mask_gl_fl        , forcing%wmask_gl_fl        , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%mask_SGD          , forcing%wmask_SGD          , mesh_new%pai_V%n_nih)
    call reallocate_dist_shared( forcing%Ti                , forcing%wTi                , mesh_new%pai_V%n_nih, mesh_new%nz)
    call reallocate_dist_shared( forcing%T_ocean           , forcing%wT_ocean           , mesh_new%pai_V%n_nih, C%nz_ocean)
    call reallocate_dist_shared( forcing%S_ocean           , forcing%wS_ocean           , mesh_new%pai_V%n_nih, C%nz_ocean)
    forcing%Hi                ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%Hi
    forcing%Hib               ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%Hib
    forcing%dHib_dx_b         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih             ) => forcing%dHib_dx_b
    forcing%dHib_dy_b         ( mesh_new%pai_Tri%i1_nih:mesh_new%pai_Tri%i2_nih             ) => forcing%dHib_dy_b
    forcing%mask_icefree_land ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%mask_icefree_land
    forcing%mask_icefree_ocean( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%mask_icefree_ocean
    forcing%mask_grounded_ice ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%mask_grounded_ice
    forcing%mask_floating_ice ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%mask_floating_ice
    forcing%mask_gl_fl        ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%mask_gl_fl
    forcing%mask_SGD          ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih               ) => forcing%mask_SGD
    forcing%Ti                ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:mesh_new%nz) => forcing%Ti
    forcing%T_ocean           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:C%nz_ocean ) => forcing%T_ocean
    forcing%S_ocean           ( mesh_new%pai_V%i1_nih  :mesh_new%pai_V%i2_nih, 1:C%nz_ocean ) => forcing%S_ocean

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_laddie_forcing


  subroutine Frankas_BMB_transition

    ! ! == Total BMB
    ! ! ============

    ! ! Initialise
    ! region%BMB%BMB        = 0._dp
    ! region%BMB%BMB_transition_phase = 0._dp

    ! if (C%do_BMB_transition_phase) then

    !   ! Safety
    !   if (C%BMB_transition_phase_t_start < C%BMB_inversion_t_start .or. C%BMB_transition_phase_t_end > C%BMB_inversion_t_end ) then
    !     ! If the window of smoothing falls outside window of BMB inversion, crash.
    !     call crash(' The time window for BMB smoothing does not fall within the time window for BMB inversion. Make sure that "BMB_transition_phase_t_start" >= "BMB_inversion_t_start", and "BMB_transition_phase_t_end" <= "BMB_inversion_t_end".')

    !   elseif (C%BMB_transition_phase_t_start >= C%BMB_transition_phase_t_end) then
    !     ! If start and end time of smoothing window is equal or start > end, crash.
    !     call crash(' "BMB_transition_phase_t_start" is equivalent or larger than "BMB_transition_phase_t_end".')

    !   end if

    !   ! Compute smoothing weights for BMB inversion smoothing
    !   if (region%time < C%BMB_transition_phase_t_start) then
    !     w = 1.0_dp

    !   elseif (region%time >= C%BMB_transition_phase_t_start .and. &
    !     region%time <= C%BMB_transition_phase_t_end) then
    !     w = 1.0_dp - ((region%time - C%BMB_transition_phase_t_start)/(C%BMB_transition_phase_t_end - C%BMB_transition_phase_t_start))

    !   elseif (region%time > C%BMB_transition_phase_t_end) then
    !     w = 0.0_dp

    !   end if

    ! end if

    ! ! Compute total BMB
    ! do vi = region%mesh%vi1, region%mesh%vi2

    !   ! Skip vertices where BMB does not operate
    !   if (.not. region%ice%mask_gl_gr( vi) .and. &
    !       .not. region%ice%mask_floating_ice( vi) .and. &
    !       .not. region%ice%mask_cf_fl( vi)) cycle

    !   if (C%do_BMB_transition_phase) then
    !     ! If BMB_transition_phase is turned ON, use weight 'w' to compute BMB field
    !     region%BMB%BMB( vi) = w * region%BMB%BMB_inv( vi) + (1.0_dp - w) * region%BMB%BMB_modelled( vi)

    !     ! Save smoothed BMB field for diagnostic output
    !     region%BMB%BMB_transition_phase( vi) = region%BMB%BMB( vi)

    !   else
    !     ! If BMB_transition_phase is turned OFF, just apply inverted melt rates
    !     region%BMB%BMB( vi) = region%BMB%BMB_inv( vi)

    !   end if

    ! end do

  end subroutine Frankas_BMB_transition

END MODULE BMB_main
