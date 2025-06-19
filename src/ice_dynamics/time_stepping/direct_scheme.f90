module direct_scheme

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use region_types, only: type_model_region
  use conservation_of_momentum_main, only: solve_stress_balance
  use time_step_criteria, only: calc_critical_timestep_SIA, calc_critical_timestep_adv
  use conservation_of_mass_main, only: calc_dHi_dt
  use inversion_utilities, only: MB_inversion
  use ice_thickness_safeties, only: alter_ice_thickness

  implicit none

  private

  public :: run_ice_dynamics_model_direct

contains

  subroutine run_ice_dynamics_model_direct( region, dt_max)
    !< Calculate a new next modelled ice thickness

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in   ) :: dt_max

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_ice_dynamics_model_direct'
    real(dp)                       :: dt_crit_SIA, dt_crit_adv, dt
    integer                        :: vi
    integer                        :: n_visc_its
    integer                        :: n_Axb_its

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. (C%choice_stress_balance_approximation == 'SIA' .OR. &
              C%choice_stress_balance_approximation == 'SSA' .OR. &
              C%choice_stress_balance_approximation == 'SIA/SSA')) then
      call crash('direct timestepping only works for SIA, SSA, or SIA/SSA ice dynamics!')
    end if

    ! Store previous ice model state
    region%ice%t_Hi_prev  = region%ice%t_Hi_next
    region%ice%Hi_prev    = region%ice%Hi_next

    ! Calculate ice velocities
    call solve_stress_balance( region%mesh, region%ice, region%bed_roughness, &
      region%BMB%BMB, region%name, n_visc_its, n_Axb_its)

    ! Calculate time step

    ! Start with the maximum allowed time step
    dt = dt_max

    ! Limit to the SIA critical time step
    if (C%choice_stress_balance_approximation == 'SIA' .OR. &
        C%choice_stress_balance_approximation == 'SIA/SSA') then
      call calc_critical_timestep_SIA( region%mesh, region%ice, dt_crit_SIA)
      dt = MIN( dt, dt_crit_SIA)
    end if

    ! Limit to the advective critical time step
    if (C%choice_stress_balance_approximation == 'SSA' .OR. &
        C%choice_stress_balance_approximation == 'SIA/SSA') then
      call calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_adv)
      dt = MIN( dt, dt_crit_adv)
    end if

    ! Limit to the smallest allowed time step
    dt = MAX( C%dt_ice_min, dt)

    ! Calculate thinning rates and predicted geometry
    call calc_dHi_dt( region%mesh, region%ice%Hi, region%ice%Hb, region%ice%SL, region%ice%u_vav_b, region%ice%v_vav_b, region%SMB%SMB, region%BMB%BMB, region%LMB%LMB, region%AMB%AMB, region%ice%fraction_margin, &
                      region%ice%mask_noice, dt, region%ice%dHi_dt, region%ice%Hi_next, region%ice%divQ, region%ice%dHi_dt_target)

    ! if so desired, invert/adjust mass balance fluxes to get an equilibrium state
    call MB_inversion( region%mesh, region%ice, region%refgeo_PD, region%SMB, region%BMB, region%LMB, region%AMB, region%ice%dHi_dt, region%ice%Hi_next, dt, region%time, region%name)

    ! Save the "raw" dynamical dH/dt before any alterations
    region%ice%dHi_dt_raw = region%ice%dHi_dt

    ! Modify predicted ice thickness if desired
    call alter_ice_thickness( region%mesh, region%ice, region%ice%Hi_prev, region%ice%Hi_next, region%refgeo_PD, region%time)

    ! Compute residual between the "raw" and final thinning rates
    region%ice%dHi_dt_residual = region%ice%dHi_dt_raw - (region%ice%Hi_next - region%ice%Hi_prev) / dt

    ! Set next modelled ice thickness timestamp
    region%ice%t_Hi_next = region%ice%t_Hi_prev + dt

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_ice_dynamics_model_direct

end module direct_scheme
