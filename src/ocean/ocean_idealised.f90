MODULE ocean_idealised

  ! Idealised ocean models

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

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ocean_model_idealised( mesh, ice, ocean, time)
    ! Calculate the ocean
    !
    ! Use an idealised ocean scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised

  SUBROUTINE initialise_ocean_model_idealised( mesh, ocean)
    ! Initialise the ocean model
    !
    ! Use an idealised ocean scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising idealised ocean model "' // &
      colour_string( TRIM( C%choice_ocean_model_idealised),'light blue') // '"...'

    ! Run the chosen idealised climate model
    IF (C%choice_climate_model_idealised == 'ISOMIP_WARM' .OR. &
        C%choice_climate_model_idealised == 'ISOMIP_COLD') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_ocean_model_idealised "' // TRIM( C%choice_ocean_model_idealised) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model_idealised

  ! == ISOMIP ==
  ! ============

  SUBROUTINE run_ocean_model_idealised_ISOMIP( mesh, ice, ocean, time)
    ! 

    IMPLICIT NONE

    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_ocean_model),             INTENT(INOUT) :: ocean
    REAL(dp),                             INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_ocean_model_idealised_ISOMIP'
    INTEGER                                             :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised_ISOMIP

END MODULE ocean_idealised
