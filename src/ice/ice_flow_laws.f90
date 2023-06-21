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

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  PURE FUNCTION calc_effective_viscosity_Glen_2D( du_dx, du_dy, dv_dx, dv_dy, A) RESULT( eta)
    ! Calculate the effective viscosity eta as a function of the strain rates du/dx,
    ! du/dy, dv/dx, dv/dy, (so excluding the vertical shear strain rates du/dz, dv/dz,
    ! and the gradients of w), according to Glen's flow law.

    IMPLICIT NONE

    ! In/output variables:
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
                 C%Glens_flow_law_epsilon_sq_0

    ! Calculate the effective viscosity eta
    eta = 0.5_dp * A**(-1._dp / C%Glens_flow_law_exponent) * (epsilon_sq)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

  END FUNCTION calc_effective_viscosity_Glen_2D

  PURE FUNCTION calc_effective_viscosity_Glen_3D_uv_only( du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, A) RESULT( eta)
    ! Calculate the effective viscosity eta as a function of the strain rates du/dx,
    ! du/dy, du/dz, dv/dx, dv/dy, dv/dz (so excluding the gradients of w),
    ! according to Glen's flow law.

    IMPLICIT NONE

    ! In/output variables:
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
                 C%Glens_flow_law_epsilon_sq_0

    ! Calculate the effective viscosity eta
    eta = 0.5_dp * A**(-1._dp / C%Glens_flow_law_exponent) * (epsilon_sq)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

  END FUNCTION calc_effective_viscosity_Glen_3D_uv_only

END MODULE ice_flow_laws
