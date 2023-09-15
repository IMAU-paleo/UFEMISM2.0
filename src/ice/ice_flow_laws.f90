MODULE ice_flow_laws

  ! Ice flow laws (currently only Glen's flow law)

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  ! ===== Glen's Flow Law =====
  ! ===========================

  ! The flow law itself, relating the effective viscosity to the effective strain rate,
  ! and a (temperature-dependent) "flow factor"

  PURE FUNCTION calc_effective_viscosity_Glen_2D( Glens_flow_law_epsilon_sq_0_applied, du_dx, du_dy, dv_dx, dv_dy, A) RESULT( eta)
    ! Calculate the effective viscosity eta as a function of the strain rates du/dx,
    ! du/dy, dv/dx, dv/dy, (so excluding the vertical shear strain rates du/dz, dv/dz,
    ! and the gradients of w), according to Glen's flow law.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                                      INTENT(IN)    :: Glens_flow_law_epsilon_sq_0_applied
    REAL(dp),                                      INTENT(IN)    :: du_dx, du_dy, dv_dx, dv_dy
    REAL(dp),                                      INTENT(IN)    :: A
    REAL(dp)                                                     :: eta

    ! Local variables:
    REAL(dp)                                                     :: epsilon_sq

    ! Calculate the square of the effective strain rate epsilon
    epsilon_sq = du_dx**2 + &
                 dv_dy**2 + &
                 du_dx * dv_dy + &
                 0.25_dp * (du_dy + dv_dx)**2 + &
                 Glens_flow_law_epsilon_sq_0_applied

    ! Calculate the effective viscosity eta
    eta = 0.5_dp * A**(-1._dp / C%Glens_flow_law_exponent) * (epsilon_sq)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

  END FUNCTION calc_effective_viscosity_Glen_2D

  PURE FUNCTION calc_effective_viscosity_Glen_3D_uv_only( Glens_flow_law_epsilon_sq_0_applied, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, A) RESULT( eta)
    ! Calculate the effective viscosity eta as a function of the strain rates du/dx,
    ! du/dy, du/dz, dv/dx, dv/dy, dv/dz (so excluding the gradients of w),
    ! according to Glen's flow law.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                                      INTENT(IN)    :: Glens_flow_law_epsilon_sq_0_applied
    REAL(dp),                                      INTENT(IN)    :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz
    REAL(dp),                                      INTENT(IN)    :: A
    REAL(dp)                                                     :: eta

    ! Local variables:
    REAL(dp)                                                     :: epsilon_sq

    ! Calculate the square of the effective strain rate epsilon
    epsilon_sq = du_dx**2 + &
                 dv_dy**2 + &
                 du_dx * dv_dy + &
                 0.25_dp * (du_dy + dv_dx)**2 + &
                 0.25_dp * (du_dz**2 + dv_dz**2) + &
                 Glens_flow_law_epsilon_sq_0_applied

    ! Calculate the effective viscosity eta
    eta = 0.5_dp * A**(-1._dp / C%Glens_flow_law_exponent) * (epsilon_sq)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

  END FUNCTION calc_effective_viscosity_Glen_3D_uv_only

  ! The calculation of the temperature-dependent flow factor

  SUBROUTINE calc_ice_rheology_Glen( mesh, ice)
    ! Calculate the flow factor A in Glen's flow law

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ice_rheology_Glen'
    INTEGER                                            :: vi,k
    REAL(dp), PARAMETER                                :: T_switch    = 263.15_dp     ! [K]           Temperature separating the "low" from the "high" temperature regime
    REAL(dp), PARAMETER                                :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1]    Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1]    Activation energy for creep in the Arrhenius relationship

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_rheology_Glen == 'uniform') THEN
      ! Apply a uniform value for the ice flow factor

      ice%A_flow = C%uniform_Glens_flow_factor

    ELSEIF (C%choice_ice_rheology_Glen == 'Huybrechts1992') THEN

      ! Calculate the ice flow factor as a function of the ice temperature according to the Arrhenius relationship (Huybrechts, 1992)
      DO vi = mesh%vi1, mesh%vi2
      DO k = 1, C%nz

        IF (ice%Ti( vi,k) < T_switch) THEN
          ice%A_flow( vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti( vi,k)))
        ELSE
          ice%A_flow( vi,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti( vi,k)))
        END IF

      END DO
      END DO

    ELSE
      CALL crash('unknown choice_ice_rheology_Glen "' // TRIM( C%choice_ice_rheology_Glen) // '"!')
    END IF

    ! Apply the flow enhancement factors
    DO vi = mesh%vi1, mesh%vi2

      IF (C%choice_enhancement_factor_transition == 'separate') THEN

        ! Totally separate values for grounded and floating areas
        IF     (ice%mask_grounded_ice( vi)) THEN
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_sheet
        ELSEIF (ice%mask_floating_ice( vi)) THEN
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_shelf
        END IF

      ELSEIF (C%choice_enhancement_factor_transition == 'interp') THEN

        IF(ice%Hi( vi) > 0._dp .AND. ice%Hib( vi) < ice%SL( vi)) THEN
          ! Interpolation between grounded and floating values depending on grounded fraction
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * (ice%fraction_gr( vi) * C%m_enh_sheet + (1._dp-ice%fraction_gr( vi)) * C%m_enh_shelf)
        ELSEIF( ice%mask_grounded_ice( vi)) THEN
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_sheet
        ELSEIF (ice%mask_floating_ice( vi)) THEN
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_shelf
        END IF

      ELSE
        CALL crash('unknown choice_enhancement_factor_transition "' // TRIM( C%choice_enhancement_factor_transition) // '"!')
      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_ice_rheology_Glen

END MODULE ice_flow_laws
