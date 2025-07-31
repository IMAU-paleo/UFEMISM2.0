module bed_roughness_nudging_H_u_flowline

  ! Contains all the routines for the basal inversion model
  ! based on flowline-averaged values of H and u

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use bed_roughness_model_types, only: type_bed_roughness_model, type_bed_roughness_nudging_model_H_u_flowline
  use netcdf_io_main, only: read_field_from_file_2D_b, find_last_output_file, find_last_timeframe
  use nudging_utilities, only: calc_nudging_vs_extrapolation_masks, trace_flowline_upstream, &
    trace_flowline_downstream, map_from_vertices_to_half_flowline, map_from_triangles_to_half_flowline, &
    calc_half_flowline_average
  use mpi_distributed_memory, only: gather_to_all
  use mesh_utilities, only: extrapolate_Gaussian
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D

  implicit none

  private

  public :: initialise_bed_roughness_nudging_H_u_flowline, run_bed_roughness_nudging_H_u_flowline

contains

  subroutine run_bed_roughness_nudging_H_u_flowline( mesh, ice, target_geometry, bed_roughness)
    ! Run the bed roughness nuding model based on flowline-averaged values of H and u

    ! In/output variables:
    type(type_mesh),                intent(in   ) :: mesh
    type(type_ice_model),           intent(in   ) :: ice
    type(type_reference_geometry),  intent(in   ) :: target_geometry
    type(type_bed_roughness_model), intent(inout) :: bed_roughness

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'run_bed_roughness_nudging_H_u_flowline'

    ! Add routine to path
    call init_routine( routine_name)

    call calc_nudging_vs_extrapolation_masks( mesh, ice, &
      bed_roughness%nudging_H_u_flowline%mask_calc_dCdt_from_nudging, &
      bed_roughness%nudging_H_u_flowline%mask_calc_dCdt_from_extrapolation, &
      bed_roughness%nudging_H_u_flowline%mask_extrapolation)

    call calc_flowline_averaged_deltaHs_deltau( mesh, ice, target_geometry, &
      bed_roughness%nudging_H_u_flowline)

    call calc_dCdt( mesh, ice, bed_roughness, &
      bed_roughness%nudging_H_u_flowline)

    ! Calculate predicted bed roughness at t+dt
    bed_roughness%generic_bed_roughness_next = max( C%generic_bed_roughness_min, min( C%generic_bed_roughness_max, &
      bed_roughness%generic_bed_roughness_prev + C%bed_roughness_nudging_dt * &
      bed_roughness%nudging_H_u_flowline%dC_dt ))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_bed_roughness_nudging_H_u_flowline

  subroutine initialise_bed_roughness_nudging_H_u_flowline( mesh, nudge)
    ! Initialise the bed roughness nudging model based on flowline-averaged values of H and u

    ! In/output variables:
    type(type_mesh),                                     intent(in   ) :: mesh
    type(type_bed_roughness_nudging_model_H_u_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_nudging_H_u_flowline'
    real(dp)                       :: bednudge_H_u_flowline_timeframe_u_target

    ! Add routine to path
    call init_routine( routine_name)

    ! Nudging masks
    allocate( nudge%mask_calc_dCdt_from_nudging      ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_calc_dCdt_from_extrapolation( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_extrapolation               ( mesh%vi1:mesh%vi2), source = 0)

    ! Target velocity
    allocate( nudge%uabs_surf_target_b               ( mesh%ti1:mesh%ti2), source = 0._dp)

    ! Half-flowline-averaged deltaHs and dHs/dt
    allocate( nudge%deltaHs_av_up  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%deltaHs_av_down( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%deltau_av_up   ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%deltau_av_down ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Intermediate terms
    allocate( nudge%R       ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%I_tot   ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%Laplac_C( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%dC_dt   ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Read target ice velocity
    if (index( C%bednudge_H_u_flowline_file_u_target,'_LAST.nc') > 1) then
      call find_last_output_file( C%bednudge_H_u_flowline_file_u_target)
      call find_last_timeframe(   C%bednudge_H_u_flowline_file_u_target, bednudge_H_u_flowline_timeframe_u_target)
      call read_field_from_file_2D_b( trim( C%bednudge_H_u_flowline_file_u_target), 'uabs_surf', &
        mesh, C%output_dir, nudge%uabs_surf_target_b, time_to_read = bednudge_H_u_flowline_timeframe_u_target)
    else
      ! Assume we're reading from a file with no time dimension
      call read_field_from_file_2D_b( trim( C%bednudge_H_u_flowline_file_u_target), 'uabs_surf', &
        mesh, C%output_dir, nudge%uabs_surf_target_b)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_nudging_H_u_flowline

  subroutine calc_flowline_averaged_deltaHs_deltau( mesh, ice, target_geometry, nudge)

    ! In/output variables:
    type(type_mesh),                                     intent(in   ) :: mesh
    type(type_ice_model),                                intent(in   ) :: ice
    type(type_reference_geometry),                       intent(in   ) :: target_geometry
    type(type_bed_roughness_nudging_model_H_u_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_flowline_averaged_deltaHs_deltau'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: deltaHs
    real(dp), dimension(mesh%ti1:mesh%ti2) :: deltau
    real(dp), dimension(mesh%nV)           :: Hi_tot
    real(dp), dimension(mesh%nV)           :: deltaHs_tot
    real(dp), dimension(mesh%nTri)         :: deltau_tot
    real(dp), dimension(mesh%nTri)         :: u_b_tot
    real(dp), dimension(mesh%nTri)         :: v_b_tot
    integer                                :: vi, ti
    real(dp), dimension(2)                 :: p
    real(dp), dimension(mesh%nV, 2)        :: trace_up, trace_down
    integer                                :: n_up, n_down
    real(dp), dimension(mesh%nV)           :: s_up, s_down
    real(dp), dimension(mesh%nV)           :: deltaHs_up, deltaHs_down
    real(dp), dimension(mesh%nTri)         :: deltau_up, deltau_down

    ! Add routine to path
    call init_routine( routine_name)

    deltaHs = ice%Hi - target_geometry%Hi

    deltau = 0._dp
    do ti = mesh%ti1, mesh%ti2
      if (.not. isnan( nudge%uabs_surf_target_b( ti))) then
        deltau( ti) = ice%uabs_surf_b( ti) - nudge%uabs_surf_target_b( ti)
      end if
    end do

    call gather_to_all( ice%Hi     , Hi_tot     )
    call gather_to_all( deltaHs    , deltaHs_tot)
    call gather_to_all( deltau     , deltau_tot )
    call gather_to_all( ice%u_vav_b, u_b_tot    )
    call gather_to_all( ice%v_vav_b, v_b_tot    )

    nudge%deltaHs_av_up   = 0._dp
    nudge%deltaHs_av_down = 0._dp
    nudge%deltau_av_up    = 0._dp
    nudge%deltau_av_down  = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (nudge%mask_calc_dCdt_from_nudging( vi)) then

        ! Trace both halves of the flowline
        p = [mesh%V( vi,1), mesh%V( vi,2)]
        call trace_flowline_upstream(   mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_up  , n_up  , s_up)
        call trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_down, n_down, s_down)

        ! Calculate thickness and velocity errors on both halves of the flowline
        call map_from_vertices_to_half_flowline( mesh, deltaHs_tot, &
          vi, trace_up, n_up, deltaHs_up)
        call map_from_vertices_to_half_flowline( mesh, deltaHs_tot, &
          vi, trace_down, n_down, deltaHs_down)
        call map_from_triangles_to_half_flowline( mesh, deltau_tot, &
          vi, trace_up, n_up, deltau_up)
        call map_from_triangles_to_half_flowline( mesh, deltau_tot, &
          vi, trace_down, n_down, deltau_down)

        ! Calculate weighted average of thickness error and thinning rates on both halves of the flowline
        call calc_half_flowline_average( s_up  , n_up  , deltaHs_up  , nudge%deltaHs_av_up  ( vi))
        call calc_half_flowline_average( s_down, n_down, deltaHs_down, nudge%deltaHs_av_down( vi))
        call calc_half_flowline_average( s_up  , n_up  , deltau_up   , nudge%deltau_av_up   ( vi))
        call calc_half_flowline_average( s_down, n_down, deltau_down , nudge%deltau_av_down ( vi))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_flowline_averaged_deltaHs_deltau

  subroutine calc_dCdt( mesh, ice, bed_roughness, nudge)

    ! In/output variables:
    type(type_mesh),                                     intent(in   ) :: mesh
    type(type_ice_model),                                intent(in   ) :: ice
    type(type_bed_roughness_model),                      intent(in   ) :: bed_roughness
    type(type_bed_roughness_nudging_model_H_u_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_dCdt'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: I1, I2, I3
    integer                                :: vi
    real(dp), dimension(mesh%ti1:mesh%ti2) :: dC_dx_b, dC_dy_b
    real(dp), dimension(mesh%vi1:mesh%vi2) :: d2C_dx2, d2C_dy2

    ! Add routine to path
    call init_routine( routine_name)

    nudge%R     = 0._dp
    I1          = 0._dp
    I2          = 0._dp
    I3          = 0._dp
    nudge%I_tot = 0._dp
    nudge%dC_dt = 0._dp

    ! Calculate Laplacian of the bed roughness field
    call ddx_a_b_2D( mesh, bed_roughness%generic_bed_roughness, dC_dx_b)
    call ddy_a_b_2D( mesh, bed_roughness%generic_bed_roughness, dC_dy_b)

    call ddx_b_a_2D( mesh, dC_dx_b, d2C_dx2)
    call ddy_b_a_2D( mesh, dC_dy_b, d2C_dy2)

    nudge%Laplac_C = d2C_dx2 + d2C_dy2
    do vi = mesh%vi1, mesh%vi2
      if (mesh%VBI( vi) > 0) nudge%Laplac_C( vi) = 0._dp
    end do

    ! Calculate dC/dt
    do vi = mesh%vi1, mesh%vi2

      if (nudge%mask_calc_dCdt_from_nudging( vi)) then

        nudge%R( vi) = max( 0._dp, min( 1._dp, &
          ((ice%uabs_vav( vi) * ice%Hi( vi)) / (C%bednudge_H_u_flowline_u_scale * C%bednudge_H_u_flowline_Hi_scale)) ))

        I1( vi) = -nudge%deltau_av_up  ( vi) / C%bednudge_H_u_flowline_u0
        I2( vi) = -nudge%deltau_av_down( vi) / C%bednudge_H_u_flowline_u0
        I3( vi) =  nudge%deltaHs_av_up ( vi) / C%bednudge_H_u_flowline_H0

        nudge%I_tot( vi) = (I1( vi) + I2( vi) + I3( vi)) * nudge%R( vi)

        nudge%dC_dt( vi) = -bed_roughness%generic_bed_roughness( vi) * &
          (nudge%I_tot( vi) / C%bednudge_H_u_flowline_t_scale &
          - C%bednudge_H_u_flowline_L**2 / C%bednudge_H_u_flowline_tau * nudge%Laplac_C( vi))

      end if

    end do

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, nudge%mask_extrapolation, nudge%dC_dt, 10e3_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dCdt

end module bed_roughness_nudging_H_u_flowline
