module conservation_of_mass_main

  ! Contains all the routines needed to calculate ice thickness rates of change (dH/dt)

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use conservation_of_mass_utilities, only: calc_ice_flux_divergence_matrix_upwind, &
    apply_mask_noice_direct, calc_flux_limited_timestep
  use conservation_of_mass_explicit, only: calc_dHi_dt_explicit, apply_ice_thickness_BC_explicit
  use conservation_of_mass_semiimplicit, only: calc_dHi_dt_semiimplicit
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD

  implicit none

  private

  public :: calc_dHi_dt, apply_mask_noice_direct, apply_ice_thickness_BC_explicit

contains

  subroutine calc_dHi_dt( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
    fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    !< Calculate ice thickness at time t+dt

    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh                  ! [-]       The model mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hi                    ! [m]       Ice thickness at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hb                    ! [m]       Bedrock elevation at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SL                    ! [m]       Water surface elevation at time t
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SMB                   ! [m yr^-1] Surface mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: BMB                   ! [m yr^-1] Basal   mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: LMB                   ! [m yr^-1] Lateral mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: AMB                   ! [m yr^-1] Artificial mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   )           :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    real(dp),                               intent(inout)           :: dt                    ! [dt]      Time step
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_dHi_dt'
    integer                        :: vi
    logical                        :: found_negative_vals
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Reset artificial mass balance
    AMB = 0._dp

    ! Calculate Hi( t+dt) with the specified time discretisation scheme
    select case (C%choice_ice_integration_method)
    case default
      call crash('unknown choice_ice_integration_method "' // trim( C%choice_ice_integration_method) // '"!')
    case ('none')
      ! Unchanging ice geometry

      Hi_tplusdt = Hi
      dHi_dt     = 0._dp
      call finalise_routine( routine_name)
      return

    case ('explicit')
      call calc_dHi_dt_explicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
        fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    case ('semi-implicit')
      call calc_dHi_dt_semiimplicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
        fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    end select

    ! Limit Hi( t+dt) to zero; throw a warning if negative thickness are encountered
    found_negative_vals = .false.
    do vi = mesh%vi1, mesh%vi2
      if (Hi_tplusdt( vi) < 0._dp) then
        ! Implicit solvers sometimes give VERY small negative numbers (e.g. -2e-189),
        ! only throw a warning if things get properly negative. Also, ignore negative
        ! values over ice-free points and very thin ice, which can experience negative
        ! mass balance, but for which it is not feasible to find a dt that prevents a
        ! negative ice thickness.

        if (Hi_tplusdt( vi) < -0.1_dp .and. Hi( vi) > C%Hi_min) found_negative_vals = .TRUE.
        ! Limit to zero
        Hi_tplusdt( vi) = 0._dp
      end if
    end do ! do vi = mesh%vi1, mesh%vi2
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_negative_vals, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (found_negative_vals) then
      call warning('encountered negative values for Hi_tplusdt - time step too large?')
    end if

    ! Add difference between corrected and uncorrected dHi_dt to residual tracker
    AMB = AMB + (Hi_tplusdt - Hi) / dt - dHi_dt

    ! Recalculate dH/dt with adjusted values of H
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt

end module conservation_of_mass_main
