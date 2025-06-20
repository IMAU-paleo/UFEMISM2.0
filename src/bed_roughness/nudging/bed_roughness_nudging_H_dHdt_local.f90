module bed_roughness_nudging_H_dHdt_local

  use precisions, only: dp
  use control_resources_and_error_messaging, only: warning, crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use bed_roughness_model_types, only: type_bed_roughness_model, type_bed_roughness_nudging_model_H_dHdt_local
  use mesh_utilities, only: extrapolate_Gaussian
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D
  use nudging_utilities, only: calc_nudging_vs_extrapolation_masks

  implicit none

  private

  public :: initialise_bed_roughness_nudging_H_dHdt_local, run_bed_roughness_nudging_H_dHdt_local

contains

  subroutine run_bed_roughness_nudging_H_dHdt_local( mesh, grid_smooth, ice, &
    target_geometry, bed_roughness_prev, bed_roughness_next, nudge)
    ! Run the bed roughness nuding model based on local values of H and dH/dt (i.e. CISM method)

    ! In/output variables:
    type(type_mesh),                                     intent(in   ) :: mesh
    type(type_grid),                                     intent(in   ) :: grid_smooth
    type(type_ice_model),                                intent(in   ) :: ice
    type(type_reference_geometry),                       intent(in   ) :: target_geometry
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(in   ) :: bed_roughness_prev
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(  out) :: bed_roughness_next
    type(type_bed_roughness_nudging_model_H_dHdt_local), intent(inout) :: nudge

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'run_bed_roughness_nudging_H_dHdt_local'
    real(dp), dimension(mesh%ti1:mesh%ti2) :: dC_dx_b, dC_dy_b
    real(dp), dimension(mesh%vi1:mesh%vi2) :: d2C_dx2, d2C_dy2
    integer                                :: vi

    ! Add routine to path
    call init_routine( routine_name)

    call ddx_a_b_2D( mesh, bed_roughness_prev, dC_dx_b)
    call ddy_a_b_2D( mesh, bed_roughness_prev, dC_dy_b)

    call ddx_b_a_2D( mesh, dC_dx_b, d2C_dx2)
    call ddy_b_a_2D( mesh, dC_dy_b, d2C_dy2)

    nudge%Laplac_C = d2C_dx2 + d2C_dy2

    call calc_nudging_vs_extrapolation_masks( mesh, ice, &
      nudge%mask_calc_dCdt_from_nudging, &
      nudge%mask_calc_dCdt_from_extrapolation, &
      nudge%mask_extrapolation)

    nudge%dC_dt = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (nudge%mask_calc_dCdt_from_nudging( vi)) then

        ! Tim's big equation
        nudge%dC_dt( vi) = -bed_roughness_prev( vi) * (&
            (ice%Hs( vi) - target_geometry%Hs( vi)) / (C%bednudge_H_dHdt_local_H0 * C%bednudge_H_dHdt_local_tau) &
          + (2._dp / C%bednudge_H_dHdt_local_H0 * ice%dHs_dt( vi)) &
          ! - C%bednudge_H_dHdt_local_r / C%bednudge_H_dHdt_local_tau * log( bed_roughness_prev( vi) / bed_roughness_target( vi) &
          - C%bednudge_H_dHdt_local_L**2 / C%bednudge_H_dHdt_local_tau * nudge%Laplac_C( vi))

      end if

    end do

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, nudge%mask_extrapolation, nudge%dC_dt, 50e3_dp)

    ! Calculate predicted bed roughness at t+dt
    bed_roughness_next = max( C%generic_bed_roughness_min, min( C%generic_bed_roughness_max, &
      bed_roughness_prev + C%bed_roughness_nudging_dt * nudge%dC_dt ))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_bed_roughness_nudging_H_dHdt_local

  subroutine initialise_bed_roughness_nudging_H_dHdt_local( mesh, nudge)
    ! Initialise the bed roughness nudging model based on local values of H and dH/dt (i.e. CISM method)

    ! In/output variables:
    type(type_mesh),                                    intent(in   ) :: mesh
    type(type_bed_roughness_nudging_model_H_dHdt_local), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_nudging_H_dHdt_local'

    ! Add routine to path
    call init_routine( routine_name)

    ! Nudging masks
    allocate( nudge%mask_calc_dCdt_from_nudging      ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_calc_dCdt_from_extrapolation( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_extrapolation               ( mesh%vi1:mesh%vi2), source = 0)

    ! Intermediate terms
    allocate( nudge%C       ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%Laplac_C( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%dC_dt   ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_nudging_H_dHdt_local

end module bed_roughness_nudging_H_dHdt_local
