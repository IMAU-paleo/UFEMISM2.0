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
  USE BMB_idealised                                          , ONLY: initialise_BMB_model_idealised, run_BMB_model_idealised
  USE BMB_prescribed                                         , ONLY: initialise_BMB_model_prescribed, run_BMB_model_prescribed
  USE BMB_parameterised                                      , ONLY: initialise_BMB_model_parameterised, run_BMB_model_parameterised
  USE BMB_laddie                                             , ONLY: initialise_BMB_model_laddie, run_BMB_model_laddie, remap_BMB_model_laddie
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE math_utilities                                         , ONLY: is_floating
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file, open_existing_netcdf_file_for_writing
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, add_time_dimension_to_file, &
                                                                     add_field_mesh_dp_2D, write_to_field_multopt_mesh_dp_2D, write_time_to_file, write_to_field_multopt_mesh_dp_3D_ocean

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model( mesh, ice, ocean, refgeo, SMB, BMB, region_name, time)
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

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which BMB model to run for this region
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
        CALL run_BMB_model_inverted( mesh, ice, BMB, time)
      CASE ('laddie')
        CALL run_BMB_model_laddie( mesh, BMB, time)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Apply subgrid scheme of old BMB to new mask
    SELECT CASE (choice_BMB_model)
      CASE ('inverted')
        ! No need to do anything
      CASE ('prescribed_fixed')
        ! No need to do anything
      CASE DEFAULT
        CALL apply_BMB_subgrid_scheme( mesh, ice, BMB)
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model

  SUBROUTINE initialise_BMB_model( mesh, ice, BMB, region_name)
    ! Initialise the BMB model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(OUT)   :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising basal mass balance model...'

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
        ! No need to do anything
      CASE ('laddie')
        CALL initialise_BMB_model_laddie( mesh, BMB)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
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
    IF (par%master) WRITE(0,'(A)') '   Writing to BMB restart file "' // &
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
    IF (par%master) WRITE(0,'(A)') '   Creating BMB model restart file "' // &
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

  SUBROUTINE remap_BMB_model( mesh_old, mesh_new, ice, BMB, region_name)
    ! Remap the BMB model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(OUT)   :: BMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_BMB_model'
    CHARACTER(LEN=256)                                    :: choice_BMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '    Remapping basal mass balance model data to the new mesh...'

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
        CALL crash('Remapping after mesh update not implemented yet for parameterised BMB')
      CASE ('inverted')
        ! No need to do anything
      CASE ('laddie')
        CALL remap_BMB_model_laddie( mesh_new, BMB)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BMB_model

! ===== Utilities =====
! =====================

  SUBROUTINE run_BMB_model_inverted( mesh, ice, BMB, time)
    ! Extrapolate inverted BMB values

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_inverted'
    INTEGER                                               :: vi, ci, vj
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: BMB_floating_ice,  BMB_gl_fl, BMB_gl_gr
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2)                :: mask_floating_ice, mask_gl_fl, mask_gl_gr
    REAL(dp)                                              :: max_cavity_size, min_neighbour

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Only run this routine if we are outside the core inversion
    ! which is performed within the ice dynamics model
    IF (time >= C%BMB_inversion_t_start .AND. &
        time <= C%BMB_inversion_t_end) THEN

      ! Save the current BMB field
      BMB%BMB_ref = BMB%BMB

      ! Save the current masks
      BMB%mask_floating_ice = ice%mask_floating_ice
      BMB%mask_gl_fl = ice%mask_gl_fl
      BMB%mask_gl_gr = ice%mask_gl_gr

      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! And exit
      RETURN

    END IF

    ! ! DENK DROM
    ! IF (time == C%start_time_of_run) THEN
    !   DO vi = mesh%vi1, mesh%vi2
    !     IF(ice%mask_floating_ice( vi)) BMB%BMB_ref( vi) = -.05_dp
    !     IF(ice%mask_gl_fl(        vi)) BMB%BMB_ref( vi) = -.2_dp
    !   END DO
    !   print*, ''
    !   print*, ':)'
    !   print*, ''
    ! END IF

    ! Set extrapolatable fields to the reference (inverted) BMB field
    BMB_floating_ice = BMB%BMB_ref
    BMB_gl_fl = BMB%BMB_ref
    BMB_gl_gr = BMB%BMB_ref

    ! Initialise extrapolation mask
    mask_floating_ice = 0
    mask_gl_fl = 0
    mask_gl_gr = 0

    ! Initialise maximum cavity size. The maximum size among
    ! all valid BMB cavities will be used as the search radius
    ! in the extrapolation later.
    max_cavity_size = MINVAL(mesh%R)

    ! Set extrapolation masks
    DO vi = mesh%vi1, mesh%vi2

      ! Skip vertices where BMB will not operate. These stay with extrapolation masks set to 0
      IF (.NOT. ice%mask_floating_ice( vi) .AND. .NOT. ice%mask_gl_gr( vi) .AND. .NOT. ice%mask_icefree_ocean( vi)) THEN
        CYCLE
      END IF

      ! Interior-shelf vertices. Values at floating grounding lines
      ! will be overwritten later after the extrapolations.
      IF (BMB%mask_floating_ice( vi)) THEN
        ! Inverted cavity: use as seed
        mask_floating_ice( vi) = 2
      ELSE
        ! New cavity: extrapolate here
        mask_floating_ice( vi) = 1
      END IF

      ! Floating-side grounding line vertices
      IF (BMB%mask_gl_fl( vi)) THEN
        ! Inverted cavity: use as seed
        mask_gl_fl( vi) = 2
      ELSE
        ! New cavity: extrapolate here
        mask_gl_fl( vi) = 1
      END IF

      ! Floating-side grounding line vertices
      IF (BMB%mask_gl_gr( vi)) THEN
        ! Inverted cavity: use as seed
        mask_gl_gr( vi) = 2
      ELSE
        ! New cavity: extrapolate here
        mask_gl_gr( vi) = 1
      END IF

      ! Check if this cavity has a lower resolution
      max_cavity_size = MAX( max_cavity_size, mesh%R( vi))

    END DO

    max_cavity_size = max_cavity_size / 3._dp

    ! == Extrapolate into new cavities
    ! ================================

    ! Perform the extrapolations - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    CALL extrapolate_Gaussian( mesh, mask_floating_ice, BMB_floating_ice, max_cavity_size)
    CALL extrapolate_Gaussian( mesh, mask_gl_fl, BMB_gl_fl, max_cavity_size)
    CALL extrapolate_Gaussian( mesh, mask_gl_gr, BMB_gl_gr, max_cavity_size)

    ! == Final BMB field
    ! ==================

    DO vi = mesh%vi1, mesh%vi2

      ! Floating-side grounding lines, or previously iced inverted ones
      IF (ice%mask_gl_fl( vi) .OR. (ice%mask_icefree_ocean( vi) .AND. BMB%mask_gl_fl( vi))) THEN
        BMB%BMB( vi) = BMB_gl_fl( vi)

      ! Interior shelves, or previously iced inverted ones
      ELSEIF (ice%mask_floating_ice( vi) .OR. (ice%mask_icefree_ocean( vi) .AND. BMB%mask_floating_ice( vi))) THEN
        BMB%BMB( vi) = BMB_floating_ice( vi)

      ! Grounded-side grounding lines, or previously iced inverted ones
      ELSEIF (ice%mask_gl_gr( vi) .OR. (ice%mask_icefree_ocean( vi) .AND. BMB%mask_gl_gr( vi))) THEN
        IF (BMB%mask_gl_gr( vi)) THEN
          ! Original grounded-side grounding line: apply full value
          BMB%BMB( vi) = BMB_gl_gr( vi)
        ELSE
          ! New grounded-side grounding line: scale value based on floating fraction
          BMB%BMB( vi) = BMB_gl_gr( vi) * (1._dp - ice%fraction_gr( vi))
        END IF

      ! Original ice-free ocean
      ELSEIF (ice%mask_icefree_ocean( vi)) THEN
        ! Use inverted value directly
        BMB%BMB( vi) = BMB%BMB_ref( vi)

      ! Not a place where we want BMB
      ELSE
        BMB%BMB( vi) = 0._dp

      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_inverted

  SUBROUTINE apply_BMB_subgrid_scheme( mesh, ice, BMB)
    ! Apply selected scheme for sub-grid shelf melt
    ! (see Leguy et al. 2021 for explanations of the three schemes)

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT)   :: BMB

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
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_BMB_subgrid_scheme

END MODULE BMB_main
