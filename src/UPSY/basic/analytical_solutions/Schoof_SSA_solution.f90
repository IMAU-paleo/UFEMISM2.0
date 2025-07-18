module Schoof_SSA_solution

  use precisions, only: dp
  use control_resources_and_error_messaging , only:  crash
  use parameters

  implicit none

  private

  public :: Schoof2006_icestream

contains

  subroutine Schoof2006_icestream( A, n, H, tantheta, L, m, y, u, tau_yield)
    ! Describes an ice stream flowing down an inclined plane. The plane slopes
    ! downward in the positive x-direction, and there is a band of increased
    ! bed slipperiness of width L running along y=0. The ice velocity u in the
    ! positive x-direction is then given by this solution.
    !
    ! Schoof, C.: A variational approach to ice stream flow, Journal of Fluid Mechanics
    !  556, 227-251, 2006

    ! In/output variables:
    real(dp), intent(in   ) :: A          ! [Pa^-3 yr^-1] Glen's flow law parameter
    real(dp), intent(in   ) :: n          ! [ ]           Glen's flow law exponent
    real(dp), intent(in   ) :: H          ! [m]           Ice thickness
    real(dp), intent(in   ) :: tantheta   ! [m]           Surface slope in the along-stream direction
    real(dp), intent(in   ) :: L          ! [m]           Ice-stream width
    real(dp), intent(in   ) :: m          ! [ ]           Ice-stream margin half-width exponent
    real(dp), intent(in   ) :: y          ! [m]           Distance from central flowline in the cross-stream direction
    real(dp), intent(  out) :: u          ! [m/yr]        Ice velocity in the along-stream direction
    real(dp), intent(  out) :: tau_yield  !               Till yield stress for the Coulomb sliding law

    ! Local variables:
    real(dp) :: B, f, W, ua, ub, uc, ud, ue

    ! Safety
    if (n /= 3._dp) call crash('Schoof only derived a solution for the case of n=3!')

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

  end subroutine Schoof2006_icestream

end module Schoof_SSA_solution