MODULE mesh_zeta

  ! Routines used in performing the vertical coordinate transformation z -> zeta

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE reallocate_mod                                         , ONLY: reallocate
  USE mesh_types                                             , ONLY: type_mesh

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Some operations on the scaled vertical coordinate

  PURE FUNCTION integrate_from_zeta_is_one_to_zeta_is_zetap( zeta, f) RESULT( integral_f)
    ! This subroutine integrates f from zeta( k=nz) = 1 (which corresponds to the ice base)
    ! to the level zetap = zeta( k) for all values of k
    !
    ! NOTE: if the integrand f is positive, the integral is negative because the integration is in
    ! the opposite zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the bottom up to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f( k+1) and f( k):
    !  integral_f( k) = integral_f( k+1) + 0.5 * (f( k+1) + f( k)) * (-dzeta)
    ! with dzeta = zeta( k+1) - zeta( k). So for f > 0,  integral_f < 0.
    !
    ! Heiko Goelzer (h.goelzer@uu.nl) Jan 2016

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp), DIMENSION( SIZE( zeta,1))                :: integral_f

    ! Local variables:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    integral_f( nz) = 0._dp

    DO k = nz-1, 1, -1
      integral_f( k) = integral_f( k+1) - 0.5_dp * (f( k+1) + f( k)) * (zeta( k+1) - zeta( k))
    END DO

  END FUNCTION integrate_from_zeta_is_one_to_zeta_is_zetap

  PURE FUNCTION integrate_from_zeta_is_zero_to_zeta_is_zetap( zeta, f) RESULT( integral_f)
    ! This subroutine integrates f from zeta( k=1) = 0 (which corresponds to the ice surface)
    ! to the level zetap = zeta( k) for all values of k
    !
    ! If the integrand f is positive, the integral is positive because the integration is in
    ! the zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the top down to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f( k) and f( k-1):
    ! integral_f( k) = integral_f( k-1) + 0.5 * (f( k) + f( k-1)) * (dzeta)
    ! with dzeta = zeta( k+1) - zeta( k).
    !
    ! Heiko Goelzer (h.goelzer@uu.nl) Jan 2016

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp), DIMENSION( SIZE( zeta,1))                :: integral_f

    ! Local variables:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    integral_f( 1) = 0._dp
    DO k = 2, nz, 1
      integral_f( k) = integral_f( k-1) + 0.5_dp * (f( k) + f( k-1)) * (zeta( k) - zeta( k-1))
    END DO

  END FUNCTION integrate_from_zeta_is_zero_to_zeta_is_zetap

  PURE FUNCTION integrate_over_zeta( zeta, f) RESULT( integral_f)
    ! Integrate f over zeta from zeta = 1 (the ice base) to zeta = 0 (the ice surface)
    !
    ! NOTE: defined so that if f is positive, then integral_f is positive too.

    IMPLICIT NONE

    ! In/output variable:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp)                                           :: integral_f

    ! Local variable:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    ! Initial value is zero
    integral_f = 0._dp

    ! Intermediate values include sum of all previous values
    ! Take current value as average between points

    DO k = 2, nz
       integral_f = integral_f + 0.5_dp * (f( k) + f( k-1)) * (zeta( k) - zeta( k-1))
    END DO

  END FUNCTION integrate_over_zeta

  PURE FUNCTION vertical_average( zeta, f) RESULT( average_f)
    ! Calculate the vertical average of any given function f defined at the vertical zeta grid.
    !
    ! The integration is in the direction of the positive zeta-axis from zeta( k=1) = 0 up to zeta( k=nz) = 1.
    ! Numerically: the average between layer k and k+1 is calculated and multiplied by the distance between those
    ! layers k and k+1, which is imediately the weight factor for this contribution because de total layer distance
    ! is scaled to 1. The sum of all the weighted contribution gives average_f the vertical average of f.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta
    REAL(dp), DIMENSION( SIZE( zeta,1)), INTENT(IN)    :: f
    REAL(dp)                                           :: average_f

    ! Local variables:
    INTEGER                                            :: nz, k

    nz = SIZE( zeta,1)

    average_f = 0._dp

    DO k = 1, nz-1
      average_f = average_f + 0.5_dp * (f( k+1) + f( k)) * (zeta( k+1) - zeta( k))
    END DO

  END FUNCTION vertical_average

END MODULE mesh_zeta