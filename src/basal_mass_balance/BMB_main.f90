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
  USE BMB_parameterised                                      , ONLY: initialise_BMB_model_parameterised, run_BMB_model_parameterised
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

    IMPLICIT NONE

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

    ! Add routine to path
    CALL init_routine( routine_name)

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
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_BMB) THEN
      ! Synchronous coupling: calculate a new BMB in every model loop
      BMB%t_next = time + C%dt_BMB
    END IF

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

    ! Run the chosen BMB model
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        BMB%BMB = C%uniform_BMB
      CASE ('idealised')
        CALL run_BMB_model_idealised( mesh, ice, BMB, time)
      CASE ('parameterised')
        CALL run_BMB_model_parameterised( mesh, ice, ocean, BMB)
      CASE ('inverted')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Extrapolate BMB
    CALL extrapolate_BMB_inland( mesh, ice, BMB)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model

  SUBROUTINE initialise_BMB_model( mesh, BMB, region_name)
    ! Initialise the BMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
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

    ! Set time of next calculation to start time
    BMB%t_next = C%start_time_of_run

    ! Determine which BMB model to initialise
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        CALL initialise_BMB_model_idealised( mesh, BMB)
      CASE ('parameterised')
        CALL initialise_BMB_model_parameterised( mesh, BMB)
      CASE ('inverted')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model

  SUBROUTINE write_to_restart_file_BMB_model( mesh, BMB, region_name, time)
    ! Write to the restart file for the BMB model

    IMPLICIT NONE

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
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        ! No need to do anything
      CASE ('inverted')
        CALL write_to_restart_file_BMB_model_region( mesh, BMB, region_name, time)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_BMB_model

  SUBROUTINE write_to_restart_file_BMB_model_region( mesh, BMB, region_name, time)
    ! Write to the restart NetCDF file for the BMB model

    IMPLICIT NONE

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

    ! ! Write the velocity fields to the file
    CALL write_to_field_multopt_mesh_dp_2D( mesh, BMB%restart_filename, ncid, 'BMB', BMB%BMB)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_BMB_model_region

  SUBROUTINE create_restart_file_BMB_model( mesh, BMB, region_name)
    ! Create the restart file for the BMB model

    IMPLICIT NONE

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
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        ! No need to do anything
      CASE ('inverted')
        CALL create_restart_file_BMB_model_region( mesh, BMB, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_BMB_model

  SUBROUTINE create_restart_file_BMB_model_region( mesh, BMB, region_name)
    ! Create a restart NetCDF file for the BMB submodel
    ! Includes generation of the procedural filename (e.g. "restart_BMB_00001.nc")

    IMPLICIT NONE

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

  SUBROUTINE remap_BMB_model( mesh_old, mesh_new, BMB, region_name)
    ! Remap the BMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
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

    ! Determine which BMB model to initialise
    SELECT CASE (choice_BMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('parameterised')
        CALL crash('Remapping after mesh update not implemented yet for parameterised BMB')
      CASE ('inverted')
        CALL crash('Remapping after mesh update not implemented yet for parameterised BMB')
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model "' // TRIM( choice_BMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BMB_model

! ===== Utilities =====
! =====================

  ! SUBROUTINE BMB_inversion( mesh, ice, SMB, BMB, dHi_dt_predicted, Hi_predicted, dt, time, region_name)
  !   ! Calculate the basal mass balance
  !   !
  !   ! Use an inversion based on the computed dHi_dt

  !   IMPLICIT NONE

  !   ! In/output variables:
  !   TYPE(type_mesh),                        INTENT(IN)    :: mesh
  !   TYPE(type_ice_model),                   INTENT(IN)    :: ice
  !   TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
  !   TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
  !   REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: dHi_dt_predicted
  !   REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi_predicted
  !   REAL(dp),                               INTENT(IN)    :: dt
  !   REAL(dp),                               INTENT(IN)    :: time
  !   CHARACTER(LEN=3)                                      :: region_name

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'BMB_inversion'
  !   INTEGER                                               :: vi
  !   CHARACTER(LEN=256)                                    :: choice_BMB_model
  !   REAL(dp)                                              :: BMB_change

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Determine filename for this model region
  !   SELECT CASE (region_name)
  !     CASE ('NAM')
  !       choice_BMB_model  = C%choice_BMB_model_NAM
  !     CASE ('EAS')
  !       choice_BMB_model  = C%choice_BMB_model_EAS
  !     CASE ('GRL')
  !       choice_BMB_model  = C%choice_BMB_model_GRL
  !     CASE ('ANT')
  !       choice_BMB_model  = C%choice_BMB_model_ANT
  !     CASE DEFAULT
  !       CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
  !   END SELECT

  !   ! Invert ocean BMB based on the full dHi_dt at each time step
  !   IF (.NOT. choice_BMB_model == 'inverted' .OR. &
  !       time < C%BMB_inversion_t_start .OR. &
  !       time > C%BMB_inversion_t_end) THEN
  !     ! Finalise routine path
  !     CALL finalise_routine( routine_name)
  !     RETURN
  !   END IF

  !   DO vi = mesh%vi1, mesh%vi2

  !     ! IF (ice%mask_gl_gr( vi) .OR. &
  !     !     ice%mask_cf_gr( vi)) THEN

  !     ! For these areas, use dHi_dt to get an "inversion" of equilibrium BMB.
  !     IF (ice%mask_cf_fl( vi)) THEN

  !       ! BMB will absorb all remaining change after calving did its thing
  !       BMB%BMB( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

  !       ! Adjust rate of ice thickness change dHi/dt to compensate the change
  !       dHi_dt_predicted( vi) = 0._dp

  !       ! Adjust corrected ice thickness to compensate the change
  !       Hi_predicted( vi) = ice%Hi_prev( vi)

  !     ELSEIF (ice%mask_floating_ice( vi)) THEN

  !       ! Basal melt will account for all change here
  !       BMB%BMB( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

  !       ! Adjust rate of ice thickness change dHi/dt to compensate the change
  !       dHi_dt_predicted( vi) = 0._dp

  !       ! Adjust corrected ice thickness to compensate the change
  !       Hi_predicted( vi) = ice%Hi_prev( vi)

  !     ELSEIF (ice%mask_icefree_ocean( vi)) THEN

  !       ! For open ocean, asume that all SMB melts
  !       BMB%BMB( vi) = -SMB%SMB( vi)

  !       ! Adjust rate of ice thickness change dHi/dt to compensate the change
  !       dHi_dt_predicted( vi) = 0._dp

  !       ! Adjust corrected ice thickness to compensate the change
  !       Hi_predicted( vi) = ice%Hi_prev( vi)

  !     ELSE
  !       ! Not a place where basal melt operates
  !       BMB%BMB( vi) = 0._dp
  !     END IF

  !   END DO ! vi = mesh%vi1, mesh%vi2

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE BMB_inversion

  SUBROUTINE extrapolate_BMB_inland( mesh, ice, BMB)
    ! Extrapolate basal rates into partially floating margin vertices

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),      INTENT(IN)    :: mesh
    TYPE(type_ice_model), INTENT(IN)    :: ice
    TYPE(type_BMB_model), INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER       :: routine_name = 'extrapolate_BMB_inland'
    INTEGER                             :: vi
    INTEGER, DIMENSION(:), ALLOCATABLE  :: mask

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Fill new floating vertices with values from surrounding ones
    ! ============================================================

    ! TBD

    ! Extrapolate inverted BMB into partially floating grounding line vertices
    ! ========================================================================

    ! If not desired, exit
    IF (.NOT. C%do_subgrid_BMB_at_grounding_line) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Allocate memory
    ALLOCATE( mask( mesh%vi1:mesh%vi2), source = 0)

    ! Extrapolate BMB from floating-side grounding line to
    ! grunded-side grounding line

    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_gl_fl( vi)) THEN
        mask( vi) = 2
      ELSEIF (ice%mask_gl_gr( vi) .AND. ice%Hib( vi) < ice%SL( vi)) THEN
        mask( vi) = 1
      ELSE
        mask( vi) = 0
      END IF
    END DO

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    CALL extrapolate_Gaussian( mesh, mask, BMB%BMB, 20000._dp)

    ! Now multiply the extrapolated values by each vertex's grounded fraction
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_gl_gr( vi) .AND. ice%Hib( vi) < ice%SL( vi)) THEN
        ! Subgrid basal melt rate
        BMB%BMB( vi) = (1._dp - ice%fraction_gr( vi)) * BMB%BMB( vi)
        ! Limit it to only melt (refreezing is tricky)
        BMB%BMB( vi) = MAX( BMB%BMB( vi), 0._dp)
      END IF
    END DO

    ! Clean up after yourself
    DEALLOCATE( mask)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extrapolate_BMB_inland

END MODULE BMB_main
