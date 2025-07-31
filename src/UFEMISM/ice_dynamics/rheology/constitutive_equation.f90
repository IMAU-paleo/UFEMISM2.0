module constitutive_equation

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model

  implicit none

  private

  public :: calc_effective_viscosity_Glen_2D, calc_effective_viscosity_Glen_3D_uv_only, &
    calc_ice_rheology_Glen

contains

  ! ===== Glen's Flow Law =====
  ! ===========================

  ! The flow law itself, relating the effective viscosity to the effective strain rate,
  ! and a (temperature-dependent) "flow factor"

  pure function calc_effective_viscosity_Glen_2D( Glens_flow_law_epsilon_sq_0_applied, &
    du_dx, du_dy, dv_dx, dv_dy, A) result( eta)
    !< Calculate the effective viscosity eta as a function of the strain rates du/dx,
    !< du/dy, dv/dx, dv/dy, (so excluding the vertical shear strain rates du/dz, dv/dz,
    !< and the gradients of w), according to Glen's flow law.

    ! In/output variables:
    real(dp), intent(in)    :: Glens_flow_law_epsilon_sq_0_applied
    real(dp), intent(in)    :: du_dx, du_dy, dv_dx, dv_dy
    real(dp), intent(in)    :: A
    real(dp)                :: eta

    ! Local variables:
    real(dp) :: epsilon_sq

    ! Calculate the square of the effective strain rate epsilon
    epsilon_sq = du_dx**2 + &
                 dv_dy**2 + &
                 du_dx * dv_dy + &
                 0.25_dp * (du_dy + dv_dx)**2 + &
                 Glens_flow_law_epsilon_sq_0_applied

    ! Calculate the effective viscosity eta
    eta = 0.5_dp * A**(-1._dp / C%Glens_flow_law_exponent) * &
      (epsilon_sq)**((1._dp - C%Glens_flow_law_exponent) / (2._dp*C%Glens_flow_law_exponent))

  end function calc_effective_viscosity_Glen_2D

  pure function calc_effective_viscosity_Glen_3D_uv_only( Glens_flow_law_epsilon_sq_0_applied, &
    du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, A) RESULT( eta)
    !< Calculate the effective viscosity eta as a function of the strain rates du/dx,
    !< du/dy, du/dz, dv/dx, dv/dy, dv/dz (so excluding the gradients of w),
    !< according to Glen's flow law.

    ! In/output variables:
    real(dp), intent(in)    :: Glens_flow_law_epsilon_sq_0_applied
    real(dp), intent(in)    :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz
    real(dp), intent(in)    :: A
    real(dp)                :: eta

    ! Local variables:
    real(dp) :: epsilon_sq

    ! Calculate the square of the effective strain rate epsilon
    epsilon_sq = du_dx**2 + &
                 dv_dy**2 + &
                 du_dx * dv_dy + &
                 0.25_dp * (du_dy + dv_dx)**2 + &
                 0.25_dp * (du_dz**2 + dv_dz**2) + &
                 Glens_flow_law_epsilon_sq_0_applied

    ! Calculate the effective viscosity eta
    eta = 0.5_dp * A**(-1._dp / C%Glens_flow_law_exponent) * &
      (epsilon_sq)**((1._dp - C%Glens_flow_law_exponent) / (2._dp*C%Glens_flow_law_exponent))

  end function calc_effective_viscosity_Glen_3D_uv_only

  ! The calculation of the temperature-dependent flow factor

  subroutine calc_ice_rheology_Glen( mesh, ice)
    !< Calculate the flow factor A in Glen's flow law

    ! In/output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_ice_rheology_Glen'
    integer                        :: vi,k
    real(dp), parameter            :: T_switch    = 263.15_dp     ! [K]           Temperature separating the "low" from the "high" temperature regime
    real(dp), parameter            :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    real(dp), parameter            :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    real(dp), parameter            :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1]    Activation energy for creep in the Arrhenius relationship
    real(dp), parameter            :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1]    Activation energy for creep in the Arrhenius relationship

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_ice_rheology_Glen)
    case default
      call crash('unknown choice_ice_rheology_Glen "' // trim( C%choice_ice_rheology_Glen) // '"!')
    case ('uniform')
      ! Apply a uniform value for the ice flow factor

      ice%A_flow = C%uniform_Glens_flow_factor

    case ('Huybrechts1992')

      ! Calculate the ice flow factor as a function of the ice temperature according to the Arrhenius relationship (Huybrechts, 1992)
      do vi = mesh%vi1, mesh%vi2
      do k = 1, C%nz

        if (ice%Ti( vi,k) < T_switch) then
          ice%A_flow( vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti( vi,k)))
        else
          ice%A_flow( vi,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti( vi,k)))
        end if

      end do
      end do

    end select

    ! Apply the flow enhancement factors
    do vi = mesh%vi1, mesh%vi2

      select case (C%choice_enhancement_factor_transition)
      case ('separate')

        ! Totally separate values for grounded and floating areas
        if (ice%mask_grounded_ice( vi)) then
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_sheet
        elseif (ice%mask_floating_ice( vi)) then
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_shelf
        end if

      case ('interp')

        if (ice%Hi( vi) > 0._dp .and. ice%Hib( vi) < ice%SL( vi)) then
          ! Interpolation between grounded and floating values depending on grounded fraction
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * &
            (       ice%fraction_gr( vi)  * C%m_enh_sheet + &
             (1._dp-ice%fraction_gr( vi)) * C%m_enh_shelf)
        elseif (ice%mask_grounded_ice( vi)) then
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_sheet
        elseif (ice%mask_floating_ice( vi)) then
          ice%A_flow( vi,:) = ice%A_flow( vi,:) * C%m_enh_shelf
        end if

      case default
        call crash('unknown choice_enhancement_factor_transition "' // trim( C%choice_enhancement_factor_transition) // '"!')
      end select

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_rheology_Glen

end module constitutive_equation
