MODULE ocean_realistic

  ! Realistic ocean models

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
  USE netcdf_input                                           , ONLY: read_field_from_file_3D_ocean
  USE netcdf_basic                                           , ONLY: field_name_options_T_ocean, field_name_options_S_ocean

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ocean_model_realistic( mesh, ice, ocean)
    ! Calculate the ocean
    !
    ! Use an realistic ocean scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model_realistic'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen realistic ocean model
    IF     (C%choice_ocean_model_realistic == 'snapshot') THEN
      ! Do nothing
    ELSE
      CALL crash('unknown choice_ocean_model_realistic "' // TRIM( C%choice_ocean_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_realistic

  SUBROUTINE initialise_ocean_model_realistic( mesh, ocean, region_name)
    ! Initialise the ocean model
    !
    ! Use an realistic ocean scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model_realistic'
    CHARACTER(LEN=256)                                    :: filename_ocean_snapshot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising realistic ocean model "' // &
      colour_string( TRIM( C%choice_ocean_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic ocean model
    IF (C%choice_ocean_model_realistic == 'snapshot') THEN
      ! Read single-time data from external file

      ! Determine which ocean model to initialise for this region
      IF     (region_name == 'NAM') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_NAM
      ELSEIF (region_name == 'EAS') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_EAS
      ELSEIF (region_name == 'GRL') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_GRL
      ELSEIF (region_name == 'ANT') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_ANT
      ELSE
        CALL crash('unknown region_name "' // region_name // '"')
      END IF

      ! Fill in  main variables
      CALL read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_T_ocean, mesh, ocean%T)
      CALL read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_S_ocean, mesh, ocean%S)

    ELSE
      CALL crash('unknown choice_ocean_model_realistic "' // TRIM( C%choice_ocean_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model_realistic

END MODULE ocean_realistic
