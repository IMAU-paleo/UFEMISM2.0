MODULE analytical_solutions

  ! Some known analytical solutions (Halfar dome, Bueler dome, Schoof ice stream)

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE Halfar_dome( A, n, H0, R0, x, y, t, H)
    ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity
    ! function with dome thickness H0 and margin radius R0 at t0.
    !
    ! Halfar, P.: On the Dynamics of Ice Sheets, Journal of Geophysical Research 86,
    !   11065-11072, 1981

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: A          ! [Pa^-3 yr^-1] Glen's flow law parameter
    REAL(dp),                            INTENT(IN)    :: n          ! [ ]           Glen's flow law exponent
    REAL(dp),                            INTENT(IN)    :: H0         ! [m]           Thickness at ice divide at t=0
    REAL(dp),                            INTENT(IN)    :: R0         ! [m]           Ice margin radius at t=0
    REAL(dp),                            INTENT(IN)    :: x          ! [m]           x-coordinate
    REAL(dp),                            INTENT(IN)    :: y          ! [m]           y-coordinate
    REAL(dp),                            INTENT(IN)    :: t          ! [yr]          time
    REAL(dp),                            INTENT(OUT)   :: H          ! [m]           Ice thickness at [x,y,t]

    ! Local variables
    REAL(dp)                                           :: Gamma, t0, r, f1, f2, f3, tp

    Gamma = (2._dp / 5._dp) * (A / sec_per_year) * (ice_density * grav)**n
    t0 = 1._dp / ((5._dp * n + 3._dp) * Gamma) * ((2._dp * n + 1._dp)/(n + 1._dp))**n * (R0**(n + 1._dp))/(H0**(2._dp * n  + 1))

    tp = (t * sec_per_year) + t0

    r = SQRT(x**2._dp + y**2._dp)

    f1 = (t0/tp)**(2._dp/(5._dp * n + 3._dp))
    f2 = (t0/tp)**(1._dp/(5._dp * n + 3._dp))
    f3 = (r/R0)

    H = H0 * f1 * MAX( 0._dp, (1._dp - (f2*f3)**((n + 1._dp)/n)))**(n/(2._dp * n + 1._dp))

  END SUBROUTINE Halfar_dome

  SUBROUTINE Bueler_dome( A, n, H0, R0, lambda, x, y, t, H, M)
    ! Extends the Halfar dome solution to include a positive accumulation rate
    !
    ! Bueler, E., Lingle, C. S., Kallen-Brown, J. A., Covey, D. N., and Bowman, L. N.:
    !   Exact solutions and verification of numerical models for isothermal ice sheets,
    !   Journal of Glaciology 51, 291-306, 2005



    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: A          ! [Pa^-3 yr^-1] Glen's flow law parameter
    REAL(dp),                            INTENT(IN)    :: n          ! [ ]           Glen's flow law exponent
    REAL(dp),                            INTENT(IN)    :: H0         ! [m]           Thickness at ice divide at t=0
    REAL(dp),                            INTENT(IN)    :: R0         ! [m]           Ice margin radius at t=0
    REAL(dp),                            INTENT(IN)    :: lambda     ! [ ]           Mass balance parameter (lambda = 5.0 gives a nice growing ice sheet)
    REAL(dp),                            INTENT(IN)    :: x          ! [m]           x-coordinate
    REAL(dp),                            INTENT(IN)    :: y          ! [m]           y-coordinate
    REAL(dp),                            INTENT(IN)    :: t          ! [yr]          time
    REAL(dp),                            INTENT(OUT)   :: H          ! [m]           Ice thickness at [x,y,t]
    REAL(dp),                            INTENT(OUT)   :: M          ! [m/yr]        Accumulation rate at [x,y,t]

    ! Local variables
    REAL(dp)                                           :: alpha, beta, Gamma, f1, f2, t0, tp, f3, f4

    alpha = (2._dp - (n+1._dp)*lambda) / ((5._dp*n)+3._dp)
    beta  = (1._dp + ((2._dp*n)+1._dp)*lambda) / ((5._dp*n)+3._dp)
    Gamma = 2._dp/5._dp * (A/sec_per_year) * (ice_density * grav)**n

    f1 = ((2._dp*n)+1)/(n+1._dp)
    f2 = (R0**(n+1._dp))/(H0**((2._dp*n)+1._dp))
    t0 = (beta / Gamma) * (f1**n) * f2

    tp = t * sec_per_year

    f1 = (tp / t0)**(-alpha)
    f2 = (tp / t0)**(-beta)
    f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
    f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
    H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))

    M = (lambda / tp) * H * sec_per_year

  END SUBROUTINE Bueler_dome

  SUBROUTINE Schoof2006_icestream( A, n, H, tantheta, L, m, y, u, tau_yield)
    ! Describes an ice stream flowing down an inclined plane. The plane slopes
    ! downward in the positive x-direction, and there is a band of increased
    ! bed slipperiness of width L running along y=0. The ice velocity u in the
    ! positive x-direction is then given by this solution.
    !
    ! Schoof, C.: A variational approach to ice stream flow, Journal of Fluid Mechanics
    !  556, 227-251, 2006

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: A          ! [Pa^-3 yr^-1] Glen's flow law parameter
    REAL(dp),                            INTENT(IN)    :: n          ! [ ]           Glen's flow law exponent
    REAL(dp),                            INTENT(IN)    :: H          ! [m]           Ice thickness
    REAL(dp),                            INTENT(IN)    :: tantheta   ! [m]           Surface slope in the along-stream direction
    REAL(dp),                            INTENT(IN)    :: L          ! [m]           Ice-stream width
    REAL(dp),                            INTENT(IN)    :: m          ! [ ]           Ice-stream margin half-width exponent
    REAL(dp),                            INTENT(IN)    :: y          ! [m]           Distance from central flowline in the cross-stream direction
    REAL(dp),                            INTENT(OUT)   :: u          ! [m/yr]        Ice velocity in the along-stream direction
    REAL(dp),                            INTENT(OUT)   :: tau_yield  !               Till yield stress for the Coulomb sliding law

    ! Local variables:
    REAL(dp)                                           :: B, f, W, ua, ub, uc, ud, ue

    ! Safety
    IF (n /= 3._dp) CALL crash('Schoof only derived a solution for the case of n=3!')

    ! Calculate the gravitational driving stress f
    f = -ice_density * grav * H * tantheta

    ! Calculate the ice hardness factor B
    B = A**(-1._dp/3._dp)

    ! Calculate the "ice stream half-width" W
    W = L * (m+1._dp)**(1._dp/m)

    ! Calculate the till yield stress across the stream
    tau_yield = f * ABS( y/L)**m

    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * H**3)
    ub = ( 1._dp / 4._dp                           ) * (   ( y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS( y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS( y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS( y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)

    ! Outside the ice-stream, velocity is zero
    IF (ABS( y) > w) u = 0._dp

  END SUBROUTINE Schoof2006_icestream

END MODULE analytical_solutions
