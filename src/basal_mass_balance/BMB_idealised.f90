MODULE BMB_idealised

  ! Idealised BMB models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE BMB_model_types                                        , ONLY: type_BMB_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model_idealised( mesh, ice, BMB, time)
    ! Calculate the basal mass balance
    !
    ! Use an idealised BMB scheme

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen idealised BMB model
    SELECT CASE (C%choice_BMB_model_idealised)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_idealised "' // TRIM( C%choice_BMB_model_idealised) // '"')
      CASE ('MISMIP+')
        CALL run_BMB_model_idealised_MISMIPplus( mesh, ice, BMB, time)
      CASE ('MISMIPplus')
        CALL run_BMB_model_idealised_MISMIPplus( mesh, ice, BMB, time)
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_idealised

  SUBROUTINE run_BMB_model_idealised_MISMIPplus( mesh, ice, BMB, time)
    ! The schematic basal melt used in the MISMIPplus experiments

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_idealised_MISMIPplus'
    INTEGER                                            :: vi
    REAL(dp)                                           :: zd, cavity_thickness

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    BMB%BMB_shelf = 0._dp

    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_floating_ice( vi)) THEN

        zd = ice%Hs( vi) - ice%Hi( vi)
        cavity_thickness = MAX( 0._dp, zd - ice%Hb( vi))

        ! Cornford et al. (2020), Eq. 7
        BMB%BMB_shelf( vi) = -0.2_dp * TANH( cavity_thickness / 75._dp) * MAX( -100._dp - zd, 0._dp)

      END IF ! IF (ice%mask_floating_ice( vi)) THEN
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_idealised_MISMIPplus

  SUBROUTINE initialise_BMB_model_idealised( mesh, BMB)
    ! Initialise the BMB model
    !
    ! Use an idealised BMB scheme

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising idealised BMB model "' // &
      colour_string( TRIM( C%choice_BMB_model_idealised),'light blue') // '"...'

    ! Initialise the chosen idealised BMB model
    SELECT CASE (C%choice_BMB_model_idealised)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_idealised "' // TRIM( C%choice_BMB_model_idealised) // '"')
      CASE ('MISMIP+')
        ! No need to do anything
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_idealised

END MODULE BMB_idealised
