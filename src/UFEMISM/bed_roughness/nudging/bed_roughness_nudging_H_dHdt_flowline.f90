module bed_roughness_nudging_H_dHdt_flowline

  ! Contains all the routines for the basal inversion model
  ! based on flowline-averaged values of H and dH/dt

  use precisions, only: dp
  use control_resources_and_error_messaging, only: warning, crash, init_routine, finalise_routine
  use model_configuration, only: C
  use mpi_basic, only: par
  use parameters
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use bed_roughness_model_types, only: type_bed_roughness_model, type_bed_roughness_nudging_model_H_dHdt_flowline
  use mesh_utilities, only: extrapolate_Gaussian
  use mpi_distributed_memory, only: gather_to_all
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D
  use mesh_data_smoothing, only: smooth_Gaussian
  use nudging_utilities, only: calc_nudging_vs_extrapolation_masks, trace_flowline_upstream, &
    trace_flowline_downstream, map_from_vertices_to_half_flowline, calc_half_flowline_average

  implicit none

  private

  public :: initialise_bed_roughness_nudging_H_dHdt_flowline, run_bed_roughness_nudging_H_dHdt_flowline

contains

  subroutine run_bed_roughness_nudging_H_dHdt_flowline( mesh, grid_smooth, ice, target_geometry, bed_roughness)
    ! Run the bed roughness nuding model based on flowline-averaged values of H and dH/dt

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_grid),                     intent(in   ) :: grid_smooth
    type(type_ice_model),                intent(in   ) :: ice
    type(type_reference_geometry),       intent(in   ) :: target_geometry
    type(type_bed_roughness_model),      intent(inout) :: bed_roughness

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'run_bed_roughness_nudging_H_dHdt_flowline'

    ! Add routine to path
    call init_routine( routine_name)

    call calc_nudging_vs_extrapolation_masks( mesh, ice, &
      bed_roughness%nudging_H_dHdt_flowline%mask_calc_dCdt_from_nudging, &
      bed_roughness%nudging_H_dHdt_flowline%mask_calc_dCdt_from_extrapolation, &
      bed_roughness%nudging_H_dHdt_flowline%mask_extrapolation)

    call calc_flowline_averaged_deltaHs_dHsdt( mesh, ice, target_geometry, &
      bed_roughness%nudging_H_dHdt_flowline)

    call calc_dCdt( mesh, ice, grid_smooth, bed_roughness, &
      bed_roughness%nudging_H_dHdt_flowline)

    ! Calculate predicted bed roughness at t+dt
    bed_roughness%generic_bed_roughness_next = max( C%generic_bed_roughness_min, min( C%generic_bed_roughness_max, &
      bed_roughness%generic_bed_roughness_prev + C%bed_roughness_nudging_dt * &
      bed_roughness%nudging_H_dHdt_flowline%dC_dt ))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_bed_roughness_nudging_H_dHdt_flowline

  subroutine initialise_bed_roughness_nudging_H_dHdt_flowline( mesh, nudge)
    ! Initialise the bed roughness nudging model based on flowline-averaged values of H and dH/dt

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    type(type_bed_roughness_nudging_model_H_dHdt_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_nudging_H_dHdt_flowline'

    ! Add routine to path
    call init_routine( routine_name)

    ! Nudging masks
    allocate( nudge%mask_calc_dCdt_from_nudging      ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_calc_dCdt_from_extrapolation( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_extrapolation               ( mesh%vi1:mesh%vi2), source = 0)

    ! Half-flowline-averaged deltaHs and dHs/dt
    allocate( nudge%deltaHs_av_up  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%deltaHs_av_down( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%dHs_dt_av_up   ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%dHs_dt_av_down ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Intermediate terms
    allocate( nudge%R    ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%I_tot( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%dC_dt( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_nudging_H_dHdt_flowline

  subroutine calc_flowline_averaged_deltaHs_dHsdt( mesh, ice, target_geometry, nudge)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    type(type_ice_model),                                   intent(in   ) :: ice
    type(type_reference_geometry),                          intent(in   ) :: target_geometry
    type(type_bed_roughness_nudging_model_H_dHdt_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_flowline_averaged_deltaHs_dHsdt'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: deltaHs
    real(dp), dimension(mesh%nV)           :: Hi_tot
    real(dp), dimension(mesh%nV)           :: deltaHs_tot
    real(dp), dimension(mesh%nV)           :: dHs_dt_tot
    real(dp), dimension(mesh%nTri)         :: u_b_tot
    real(dp), dimension(mesh%nTri)         :: v_b_tot
    integer                                :: vi
    real(dp), dimension(2)                 :: p
    real(dp), dimension(mesh%nV, 2)        :: trace_up, trace_down
    integer                                :: n_up, n_down
    real(dp), dimension(mesh%nV   )        :: s_up, s_down
    real(dp), dimension(mesh%nV   )        :: deltaHs_up, deltaHs_down
    real(dp), dimension(mesh%nV   )        :: dHs_dt_up, dHs_dt_down

    ! Add routine to path
    call init_routine( routine_name)

    deltaHs = ice%Hs - target_geometry%Hs

    call gather_to_all( ice%Hi     , Hi_tot     )
    call gather_to_all( deltaHs    , deltaHs_tot)
    call gather_to_all( ice%dHs_dt , dHs_dt_tot )
    call gather_to_all( ice%u_vav_b, u_b_tot    )
    call gather_to_all( ice%v_vav_b, v_b_tot    )

    nudge%deltaHs_av_up   = 0._dp
    nudge%deltaHs_av_down = 0._dp
    nudge%dHs_dt_av_up    = 0._dp
    nudge%dHs_dt_av_down  = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (nudge%mask_calc_dCdt_from_nudging( vi)) then

        ! Trace both halves of the flowline
        p = [mesh%V( vi,1), mesh%V( vi,2)]
        call trace_flowline_upstream(   mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_up  , n_up  , s_up)
        call trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_down, n_down, s_down)

        ! Calculate thickness error and thinning rates on both halves of the flowline
        call map_from_vertices_to_half_flowline( mesh, deltaHs_tot, &
          vi, trace_up, n_up, deltaHs_up)
        call map_from_vertices_to_half_flowline( mesh, deltaHs_tot, &
          vi, trace_down, n_down, deltaHs_down)
        call map_from_vertices_to_half_flowline( mesh, dHs_dt_tot, &
          vi, trace_up, n_up, dHs_dt_up)
        call map_from_vertices_to_half_flowline( mesh, dHs_dt_tot, &
          vi, trace_down, n_down, dHs_dt_down)

        ! Calculate weighted average of thickness error and thinning rates on both halves of the flowline
        call calc_half_flowline_average( s_up  , n_up  , deltaHs_up  , nudge%deltaHs_av_up  ( vi))
        call calc_half_flowline_average( s_up  , n_up  , dHs_dt_up   , nudge%dHs_dt_av_up   ( vi))
        call calc_half_flowline_average( s_down, n_down, deltaHs_down, nudge%deltaHs_av_down( vi))
        call calc_half_flowline_average( s_down, n_down, dHs_dt_down , nudge%dHs_dt_av_down ( vi))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_flowline_averaged_deltaHs_dHsdt

  subroutine calc_dCdt( mesh, ice, grid_smooth, bed_roughness, nudge)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    type(type_ice_model),                                   intent(in   ) :: ice
    type(type_grid),                                        intent(in   ) :: grid_smooth
    type(type_bed_roughness_model),                         intent(in   ) :: bed_roughness
    type(type_bed_roughness_nudging_model_H_dHdt_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_dCdt'
    integer                                :: vi

    ! Add routine to path
    call init_routine( routine_name)

    nudge%R     = 0._dp
    nudge%I_tot = 0._dp
    nudge%dC_dt = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (nudge%mask_calc_dCdt_from_nudging( vi)) then

        nudge%R( vi) = max( 0._dp, min( 1._dp, &
          ((ice%uabs_vav( vi) * ice%Hi( vi)) / (C%bednudge_H_dHdt_flowline_u_scale * C%bednudge_H_dHdt_flowline_Hi_scale)) ))

        nudge%I_tot( vi) = (&
          (nudge%deltaHs_av_up( vi) - 0.25_dp * nudge%deltaHs_av_down( vi)) / C%bednudge_H_dHdt_flowline_dH0 + &
          (nudge%dHs_dt_av_up(  vi) - 0.25_dp * nudge%dHs_dt_av_down(  vi)) / C%bednudge_H_dHdt_flowline_dHdt0)

        nudge%dC_dt( vi) = -1._dp * (nudge%I_tot( vi) * bed_roughness%generic_bed_roughness( vi)) / C%bednudge_H_dHdt_flowline_t_scale

      end if

    end do

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, nudge%mask_extrapolation, nudge%dC_dt, 10e3_dp)

    call reduce_dCdt_on_steep_slopes( mesh, ice, nudge%dC_dt)

    call smooth_dCdt( mesh, grid_smooth, nudge%dC_dt)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dCdt

  subroutine reduce_dCdt_on_steep_slopes( mesh, ice, dC_dt)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: dC_dt

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'reduce_dCdt_on_steep_slopes'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dHs_dx, dHs_dy, abs_grad_Hs
    real(dp)                               :: fg_exp_mod
    integer                                :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate surface slopes
    call ddx_a_a_2D( mesh, ice%Hs, dHs_dx)
    call ddy_a_a_2D( mesh, ice%Hs, dHs_dy)

    ! Calculate absolute surface gradient
    abs_grad_Hs = sqrt( dHs_dx**2 + dHs_dy**2)

    ! Scale bed roughness rate of change for partially grounded, steep-sloped areas
    do vi = mesh%vi1, mesh%vi2

      ! Ice margin and grounding lines
      if (ice%mask_grounded_ice( vi)) then

        ! Strengthen the effect of grounded fractions for steep slopes
        fg_exp_mod = min( 1.0_dp, max( 0._dp, max( 0._dp, abs_grad_Hs( vi) - 0.02_dp) / (0.06_dp - 0.02_dp) ))

        ! Scale based on grounded fraction
        dC_dt( vi) = dC_dt( vi) * ice%fraction_gr( vi) ** (1._dp + fg_exp_mod)

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reduce_dCdt_on_steep_slopes

  subroutine smooth_dCdt( mesh, grid_smooth, dC_dt)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_grid),                        intent(in   ) :: grid_smooth
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: dC_dt

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'smooth_dCdt'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dC_dt_smoothed

    ! Add routine to path
    call init_routine( routine_name)

    dC_dt_smoothed = dC_dt
    call smooth_Gaussian( mesh, grid_smooth, C%output_dir, dC_dt_smoothed, C%bednudge_H_dHdt_flowline_r_smooth)

    dC_dt = (1._dp - C%bednudge_H_dHdt_flowline_w_smooth) * dC_dt + &
                     C%bednudge_H_dHdt_flowline_w_smooth  * dC_dt_smoothed

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_dCdt

end module bed_roughness_nudging_H_dHdt_flowline
