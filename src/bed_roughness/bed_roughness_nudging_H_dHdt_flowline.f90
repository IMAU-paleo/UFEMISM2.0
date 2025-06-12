module bed_roughness_nudging_H_dHdt_flowline

  ! Contains all the routines for the basal inversion model
  ! based on flowline-averaged values of H and dH/dt

  use precisions, only: dp
  use control_resources_and_error_messaging, only: warning, crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use bed_roughness_model_types, only: type_bed_roughness_model, type_bed_roughness_nudging_model_H_dHdt_flowline
  use mesh_utilities, only: find_containing_vertex, find_containing_triangle, extrapolate_Gaussian, &
    interpolate_to_point_dp_2D_singlecore
  use plane_geometry, only: triangle_area
  use mpi_distributed_memory, only: gather_to_all
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D
  use mesh_data_smoothing, only: smooth_Gaussian
  use netcdf_io_main

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

    call calc_nudging_vs_extrapolation_masks( mesh, ice, target_geometry, &
      bed_roughness%nudging_H_dHdt_flowline)

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
    allocate( nudge%mask_Hs_is_converging            ( mesh%vi1:mesh%vi2), source = .false.)
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

  subroutine calc_nudging_vs_extrapolation_masks( mesh, ice, target_geometry, nudge)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    type(type_ice_model),                                   intent(in   ) :: ice
    type(type_reference_geometry),                          intent(in   ) :: target_geometry
    type(type_bed_roughness_nudging_model_H_dHdt_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_nudging_vs_extrapolation_masks'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    nudge%mask_calc_dCdt_from_nudging       = .false.
    nudge%mask_calc_dCdt_from_extrapolation = .false.
    nudge%mask_Hs_is_converging             = .false.
    nudge%mask_extrapolation                = 0

    do vi = mesh%vi1, mesh%vi2

      ! Determine whether bed roughness should be
      ! updated by inversion or by extrapolation

      ! Only perform the inversion on fully grounded vertices
      if (ice%mask_grounded_ice( vi) .and. ice%Hi( vi) > 100._dp .and. &
        .not. (ice%mask_margin( vi) .or. ice%mask_gl_gr( vi) .or. ice%mask_cf_gr( vi))) then

        ! Perform the inversion here
        nudge%mask_calc_dCdt_from_nudging      ( vi) = .true.
        nudge%mask_calc_dCdt_from_extrapolation( vi) = .false.
        nudge%mask_extrapolation               ( vi) = 2

        ! If Hs is already converging to the target value, do not nudge bed roughness further
        if (ice%dHs_dt( vi) * (ice%Hs( vi) - target_geometry%Hs( vi)) < 0._dp) then
          ! nudge%mask_Hs_is_converging( vi) = .true.
        end if

      else

        ! Extrapolate here
        nudge%mask_calc_dCdt_from_nudging      ( vi) = .false.
        nudge%mask_calc_dCdt_from_extrapolation( vi) = .true.
        nudge%mask_extrapolation               ( vi) = 1

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_nudging_vs_extrapolation_masks

  subroutine calc_flowline_averaged_deltaHs_dHsdt( mesh, ice, target_geometry, nudge)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    type(type_ice_model),                                   intent(in   ) :: ice
    type(type_reference_geometry),                          intent(in   ) :: target_geometry
    type(type_bed_roughness_nudging_model_H_dHdt_flowline), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_flowline_averaged_deltaHs_dHsdt'
    real(dp), dimension(mesh%nV)    :: Hi_tot
    real(dp), dimension(mesh%nV)    :: Hs_tot
    real(dp), dimension(mesh%nV)    :: Hs_target_tot
    real(dp), dimension(mesh%nV)    :: dHs_dt_tot
    real(dp), dimension(mesh%nTri)  :: u_b_tot
    real(dp), dimension(mesh%nTri)  :: v_b_tot
    integer                         :: vi
    real(dp), dimension(2)          :: p
    real(dp), dimension(mesh%nV, 2) :: trace_up, trace_down
    integer                         :: n_up, n_down
    real(dp), dimension(mesh%nV   ) :: s_up, s_down
    real(dp), dimension(mesh%nV   ) :: deltaHs_up, deltaHs_down
    real(dp), dimension(mesh%nV   ) :: dHs_dt_up, dHs_dt_down

    ! Add routine to path
    call init_routine( routine_name)

    call gather_to_all( ice%Hi            , Hi_tot       )
    call gather_to_all( ice%Hs            , Hs_tot       )
    call gather_to_all( target_geometry%Hs, Hs_target_tot)
    call gather_to_all( ice%dHs_dt        , dHs_dt_tot   )
    call gather_to_all( ice%u_vav_b       , u_b_tot      )
    call gather_to_all( ice%v_vav_b       , v_b_tot      )

    nudge%deltaHs_av_up   = 0._dp
    nudge%deltaHs_av_down = 0._dp
    nudge%dHs_dt_av_up    = 0._dp
    nudge%dHs_dt_av_down  = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (nudge%mask_calc_dCdt_from_nudging( vi) .and. .not. nudge%mask_Hs_is_converging( vi)) then

        ! Trace both halves of the flowline
        p = [mesh%V( vi,1), mesh%V( vi,2)]
        call trace_flowline_upstream(   mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_up  , n_up  , s_up)
        call trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_down, n_down, s_down)

        ! If we couldn't trace the flowline here, extrapolate instead of inverting
        if (n_up < 3 .or. n_down < 3) then
          nudge%mask_calc_dCdt_from_nudging      ( vi) = .false.
          nudge%mask_calc_dCdt_from_extrapolation( vi) = .true.
          nudge%mask_extrapolation( vi) = 1
          cycle
        end if

        ! Calculate thickness error and thinning rates on both halves of the flowline
        call calc_deltaHs_dHdt_along_flowline( mesh, Hs_tot, Hs_target_tot, dHs_dt_tot, &
          vi, trace_up, n_up, deltaHs_up, dHs_dt_up)
        call calc_deltaHs_dHdt_along_flowline( mesh, Hs_tot, Hs_target_tot, dHs_dt_tot, &
          vi, trace_down, n_down, deltaHs_down, dHs_dt_down)

        ! Calculate weighted average of thickness error and thinning rates on both halves of the flowline
        call calc_flowline_average( s_up  , n_up  , deltaHs_up  , nudge%deltaHs_av_up  ( vi))
        call calc_flowline_average( s_up  , n_up  , dHs_dt_up   , nudge%dHs_dt_av_up   ( vi))
        call calc_flowline_average( s_down, n_down, deltaHs_down, nudge%deltaHs_av_down( vi))
        call calc_flowline_average( s_down, n_down, dHs_dt_down , nudge%dHs_dt_av_down ( vi))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_flowline_averaged_deltaHs_dHsdt

  subroutine calc_deltaHs_dHdt_along_flowline( mesh, Hs_tot, Hs_target_tot, dHs_dt_tot, &
    vi, T, n, deltaHs, dHs_dt)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    real(dp), dimension(mesh%nV), intent(in   ) :: Hs_tot
    real(dp), dimension(mesh%nV), intent(in   ) :: Hs_target_tot
    real(dp), dimension(mesh%nV), intent(in   ) :: dHs_dt_tot
    integer,                      intent(in   ) :: vi
    real(dp), dimension(:,:),     intent(in   ) :: T
    integer,                      intent(in   ) :: n
    real(dp), dimension(mesh%nV), intent(  out) :: deltaHs, dHs_dt

    ! Local variables:
    integer                :: ti, k
    real(dp), dimension(2) :: p
    real(dp)               :: Hs_mod, Hs_target, dHs_dt_mod

    deltaHs = 0._dp
    dHs_dt  = 0._dp
    ti      = mesh%iTri( vi,1)

    do k = 1, n

      p = T( k,:)

      call interpolate_to_point_dp_2D_singlecore( mesh, Hs_tot       , p, ti, Hs_mod)
      call interpolate_to_point_dp_2D_singlecore( mesh, Hs_target_tot, p, ti, Hs_target)
      call interpolate_to_point_dp_2D_singlecore( mesh, dHs_dt_tot   , p, ti, dHs_dt_mod)

      deltaHs( k) = Hs_mod - Hs_target
      dHs_dt(  k) = dHs_dt_mod

    end do

  end subroutine calc_deltaHs_dHdt_along_flowline

  subroutine calc_flowline_average( s, n, d, d_av)

    ! In/output variables:
    real(dp), dimension(:), intent(in   ) :: s
    integer,                intent(in   ) :: n
    real(dp), dimension(:), intent(in   ) :: d
    real(dp),               intent(  out) :: d_av

    ! Local variables:
    real(dp) :: int_w_d, int_w
    real(dp) :: s1, s2, ds, w1, w2, w_av, d1, d2, dd
    integer  :: k

    ! Trivial cases
    if (n == 0) then
      call crash('calc_flowline_average - flowline has length zero')
    elseif (n == 2) then
      d_av = d(1)
      return
    end if

    int_w_d = 0._dp
    int_w   = 0._dp

    do k = 2, n

      ! Distance of both points
      s1 = s( k-1)
      s2 = s( k  )
      ds = s2 - s1

      ! Weights for both points
      w1   = (2._dp / s( n)) * (1._dp - s1 / s( n))
      w2   = (2._dp / s( n)) * (1._dp - s2 / s( n))
      w_av = (w1 + w2) / 2._dp

      ! Thickness error and thinning rate for both points
      d1 = d( k-1)
      d2 = d( k)
      dd = (d1 + d2) / 2._dp

      ! Add to integrals
      int_w_d = int_w_d + (w_av * dd * ds)
      int_w   = int_w   + (w_av      * ds)

    end do

    d_av = int_w_d / int_w

  end subroutine calc_flowline_average

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

      if (nudge%mask_calc_dCdt_from_nudging( vi) .and. .not. nudge%mask_Hs_is_converging( vi)) then

        nudge%R( vi) = max( 0._dp, min( 1._dp, &
          ((ice%uabs_vav( vi) * ice%Hi( vi)) / (C%bednudge_H_dHdt_flowline_u_scale * C%bednudge_H_dHdt_flowline_Hi_scale)) ))

        ! nudge%I_tot( vi) = nudge%R( vi) * (&
        !   (nudge%deltaHs_av_up( vi)                             ) / C%bednudge_H_dHdt_flowline_dH0 + &
        !   (nudge%dHs_dt_av_up(  vi) + nudge%dHs_dt_av_down(  vi)) / C%bednudge_H_dHdt_flowline_dHdt0)
        ! nudge%I_tot( vi) = nudge%R( vi) * (&
        !   (nudge%deltaHs_av_up( vi)) / C%bednudge_H_dHdt_flowline_dH0 + &
        !   (nudge%dHs_dt_av_up(  vi)) / C%bednudge_H_dHdt_flowline_dHdt0)
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
    call smooth_Gaussian( mesh, grid_smooth, dC_dt_smoothed, C%bednudge_H_dHdt_flowline_r_smooth)

    dC_dt = (1._dp - C%bednudge_H_dHdt_flowline_w_smooth) * dC_dt + &
                     C%bednudge_H_dHdt_flowline_w_smooth * dC_dt_smoothed

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_dCdt

  subroutine trace_flowline_upstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n, s)
    ! Trace the flowline passing through point p upstream through
    ! the 2-D velocity field u_b,v_b.
    !
    ! Returns a list T of n points on the flowline
    !
    ! Stop the trace when it encounters the ice divide (defined as an ice velocity lower than 1 m/yr)

    ! Input variables:
    type(type_mesh),                     intent(in   ) :: mesh
    real(dp), dimension(mesh%nV),        intent(in   ) :: Hi_tot
    real(dp), dimension(mesh%nTri),      intent(in   ) :: u_b_tot, v_b_tot
    real(dp), dimension(2),              intent(in   ) :: p
    real(dp), dimension(:,:  ),          intent(  out) :: T
    integer,                             intent(  out) :: n
    real(dp), dimension(:),              intent(  out) :: s

    ! Local variables:
    real(dp), dimension(2) :: pt
    integer                :: vi, iti, ti, k
    real(dp)               :: dist, w, w_tot
    real(dp)               :: u_pt, v_pt, uabs_pt
    real(dp), dimension(2) :: u_hat_pt
    real(dp)               :: dist_prev

    ! Safety - if there's no ice, we can't do a trace
    vi = 1
    call find_containing_vertex( mesh, p, vi)
    if (Hi_tot( vi) < 1._dp) then
      T = 0._dp
      n = 1
      T( 1,:) = p
      return
    end if

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p

    ! Safety
    dist_prev = 0._dp

    do while (.true.)

      ! Find the vertex vi containing the tracer
      call find_containing_vertex( mesh, pt, vi)

      ! Interpolate between the surrounding triangles to find
      ! the velocities at the tracer's location

      w_tot = 0._dp
      u_pt  = 0._dp
      v_pt  = 0._dp

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        dist = norm2( mesh%TriGC( ti,:) - pt)
        w = 1._dp / dist**2
        w_tot = w_tot + w
        u_pt  = u_pt  + w * u_b_tot( ti)
        v_pt  = v_pt  + w * v_b_tot( ti)
      end do

      u_pt = u_pt / w_tot
      v_pt = v_pt / w_tot

      ! Calculate absolute velocity at the tracer's location
      uabs_pt = sqrt( u_pt**2 + v_pt**2)

      ! if we've reached the ice divide (defined as the place where
      ! we find velocities below 1 m/yr), end the trace
      if (uabs_pt < 1e-2_dp) exit

      ! Calculate the normalised velocity vector at the tracer's location
      u_hat_pt = [u_pt / uabs_pt, v_pt / uabs_pt]

      ! Add current position to the traces
      n = n + 1
      ! Safety
      if (n > size( T,1)) call crash('upstream flowline tracer got stuck!')
      T( n,:) = pt

      ! Save previous distance-to-origin
      dist_prev = norm2( pt - p)

      ! Move the tracer upstream by a distance of one local resolution
      pt = pt - u_hat_pt * mesh%R( vi)

      ! if the new distance-to-origin is shorter than the previous one, end the trace
      if (norm2( pt - p) < dist_prev) EXIT

      ! if the new tracer location is outside the domain, end the trace
      if (pt( 1) <= mesh%xmin .or. pt( 2) >= mesh%xmax .or. &
          pt( 2) <= mesh%ymin .or. pt( 2) >= mesh%ymax) exit

    end do

    ! Safety
    if (n == 0) then
      n = 1
      T( 1,:) = p
    end if

    ! Calculate distance along both halves of the flowline
    s = 0._dp
    do k = 2, n
      s( k) = s( k-1) + norm2( T( k,:) - T( k-1,:))
    end do

  end subroutine trace_flowline_upstream

  subroutine trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n, s)
    ! Trace the flowline passing through point p downstream through
    ! the 2-D velocity field u_b,v_b.
    !
    ! Returns a list T of n points on the flowline
    !
    ! Stop the trace when it encounters the ice margin

    ! Input variables:
    type(type_mesh),                intent(in   ) :: mesh
    real(dp), dimension(mesh%nV),   intent(in   ) :: Hi_tot
    real(dp), dimension(mesh%nTri), intent(in   ) :: u_b_tot, v_b_tot
    real(dp), dimension(2),         intent(in   ) :: p
    real(dp), dimension(:,:),       intent(  out) :: T
    integer,                        intent(  out) :: n
    real(dp), dimension(:),         intent(  out) :: s

    ! Local variables:
    real(dp), dimension(2) :: pt
    integer                :: vi, iti, ti, k
    real(dp)               :: dist, w, w_tot
    real(dp)               :: u_pt, v_pt, uabs_pt
    real(dp), dimension(2) :: u_hat_pt
    real(dp)               :: dist_prev

    ! Safety - if there's no ice, we can't do a trace
    vi = 1
    call find_containing_vertex( mesh, p, vi)
    if (Hi_tot( vi) < 1._dp) then
      T = 0._dp
      n = 1
      T( 1,:) = p
      return
    end if

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p

    ! Safety
    dist_prev = 0._dp

    do while (.true.)

      ! Find the vertex vi containing the tracer
      call find_containing_vertex( mesh, pt, vi)

      ! if ice thickness in this vertex is below 1 m, assume we've found the
      ! ice margin, and end the trace
      if (Hi_tot( vi) < 0.1_dp) exit

      ! Interpolate between the surrounding triangles to find
      ! the velocities at the tracer's location

      w_tot = 0._dp
      u_pt  = 0._dp
      v_pt  = 0._dp

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        dist = norm2( mesh%TriGC( ti,:) - pt)
        if (dist == 0._dp) call crash('whaa!')
        w = 1._dp / dist**2
        w_tot = w_tot + w
        u_pt  = u_pt  + w * u_b_tot( ti)
        v_pt  = v_pt  + w * v_b_tot( ti)
      end do

      u_pt = u_pt / w_tot
      v_pt = v_pt / w_tot

      ! Calculate absolute velocity at the tracer's location
      uabs_pt = sqrt( u_pt**2 + v_pt**2)

      ! if we're at the ice divide (defined as the place where
      ! we find velocities below 1 m/yr), we can't do the trace
      ! if (uabs_pt < 1._dp) exit

      ! Calculate the normalised velocity vector at the tracer's location
      u_hat_pt = [u_pt / uabs_pt, v_pt / uabs_pt]

      ! Add current position to the traces
      n = n + 1
      ! Safety
      if (n > size( T,1)) call crash('downstream flowline tracer got stuck!')
      T( n,:) = pt

      ! Save previous distance-to-origin
      dist_prev = norm2( pt - p)

      ! Move the tracer downstream by a distance of one local resolution
      pt = pt + u_hat_pt * mesh%R( vi)

      ! if the new distance-to-origin is shorter than the previous one, end the trace
      if (norm2( pt - p) < dist_prev) exit

      ! if the new tracer location is outside the domain, end the trace
      if (pt( 1) <= mesh%xmin .or. pt( 2) >= mesh%xmax .or. &
          pt( 2) <= mesh%ymin .or. pt( 2) >= mesh%ymax) exit

    end do

    ! Safety
    if (n == 0) then
      n = 1
      T( 1,:) = p
    end if

    ! Calculate distance along both halves of the flowline
    s = 0._dp
    do k = 2, n
      s( k) = s( k-1) + norm2( T( k,:) - T( k-1,:))
    end do

  end subroutine trace_flowline_downstream

end module bed_roughness_nudging_H_dHdt_flowline
