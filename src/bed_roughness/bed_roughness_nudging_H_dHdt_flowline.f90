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
  use bed_roughness_model_types, only: type_bed_roughness_model
  use mesh_utilities, only: find_containing_vertex, find_containing_triangle, extrapolate_Gaussian, &
    interpolate_to_point_dp_2D_singlecore
  use plane_geometry, only: triangle_area
  use mpi_distributed_memory, only: gather_to_all
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D
  use mesh_data_smoothing, only: smooth_Gaussian

  implicit none

  private

  public :: initialise_bed_roughness_nudging_H_dHdt_flowline, run_bed_roughness_nudging_H_dHdt_flowline

contains

  subroutine run_bed_roughness_nudging_H_dHdt_flowline( mesh, grid_smooth, ice, refgeo, bed_roughness)
    ! Run the bed roughness nuding model based on flowline-averaged values of H and dH/dt

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_grid),                     intent(in   ) :: grid_smooth
    type(type_ice_model),                intent(in   ) :: ice
    type(type_reference_geometry),       intent(in   ) :: refgeo
    type(type_bed_roughness_model),      intent(inout) :: bed_roughness

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'run_bed_roughness_nudging_H_dHdt_flowline'
    integer,  dimension(:    ), allocatable :: mask
    real(dp), dimension(mesh%nV)            :: Hi_tot
    real(dp), dimension(mesh%nV)            :: Hs_tot
    real(dp), dimension(mesh%nV)            :: Hs_target_tot
    real(dp), dimension(mesh%nV)            :: dHs_dt_tot
    real(dp), dimension(mesh%nTri)          :: u_b_tot
    real(dp), dimension(mesh%nTri)          :: v_b_tot
    logical,  dimension(mesh%nV)            :: mask_grounded_ice_tot
    logical,  dimension(mesh%nV)            :: mask_gl_gr_tot
    logical,  dimension(mesh%nV)            :: mask_margin_tot
    real(dp), dimension(mesh%nV)            :: fraction_gr_tot
    integer                                 :: vi
    real(dp), dimension(2)                  :: p
    real(dp), dimension(:,:  ), allocatable :: trace_up, trace_down
    real(dp), dimension(:    ), allocatable :: s_up, s_down
    real(dp), dimension(:    ), allocatable :: deltaHs_up, deltaHs_down
    real(dp), dimension(:    ), allocatable :: dHs_dt_up, dHs_dt_down
    integer                                 :: n_up,n_down
    integer                                 :: k
    real(dp), dimension(2)                  :: pt
    integer                                 :: ti
    real(dp)                                :: Hs_mod, Hs_target, dHs_dt_mod
    real(dp)                                :: s1, s2, w1, w2, deltaHs1, deltaHs2, dHs_dt1, dHs_dt2, w_av, deltaHs_av, dHs_dt_av, ds
    real(dp)                                :: int_w_deltaHs_up, int_w_dHs_dt_up, int_w_up
    real(dp)                                :: int_w_deltaHs_down, int_w_dHs_dt_down, int_w_down
    real(dp), dimension(:    ), allocatable :: deltaHs_av_up, deltaHs_av_down
    real(dp), dimension(:    ), allocatable :: dHs_dt_av_up, dHs_dt_av_down
    real(dp), dimension(:    ), allocatable :: I_tot, R
    real(dp), dimension(:    ), allocatable :: dC_dt
    real(dp), dimension(:    ), allocatable :: dHs_dx, dHs_dy, abs_grad_Hs
    real(dp)                                :: fg_exp_mod
    real(dp), dimension(:    ), allocatable :: dC_dt_smoothed
    real(dp)                                :: misfit

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    allocate( mask(            mesh%vi1:mesh%vi2), source = 0      )
    allocate( trace_up(        mesh%nV, 2       ), source = 0._dp  )
    allocate( trace_down(      mesh%nV, 2       ), source = 0._dp  )
    allocate( s_up(            mesh%nV          ), source = 0._dp  )
    allocate( s_down(          mesh%nV          ), source = 0._dp  )
    allocate( deltaHs_up(      mesh%nV          ), source = 0._dp  )
    allocate( deltaHs_down(    mesh%nV          ), source = 0._dp  )
    allocate( dHs_dt_up(       mesh%nV          ), source = 0._dp  )
    allocate( dHs_dt_down(     mesh%nV          ), source = 0._dp  )
    allocate( deltaHs_av_up(   mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( deltaHs_av_down( mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( dHs_dt_av_up(    mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( dHs_dt_av_down(  mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( R(               mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( I_tot(           mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( dC_dt(           mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( dHs_dx(          mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( dHs_dy(          mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( abs_grad_Hs(     mesh%vi1:mesh%vi2), source = 0._dp  )
    allocate( dC_dt_smoothed(  mesh%vi1:mesh%vi2), source = 0._dp  )

    ! Gather ice model data from all processes
    call gather_to_all( ice%Hi               , Hi_tot               )
    call gather_to_all( ice%Hs               , Hs_tot               )
    call gather_to_all( refgeo%Hs            , Hs_target_tot        )
    call gather_to_all( ice%dHs_dt           , dHs_dt_tot           )
    call gather_to_all( ice%u_vav_b          , u_b_tot              )
    call gather_to_all( ice%v_vav_b          , v_b_tot              )
    call gather_to_all( ice%mask_grounded_ice, mask_grounded_ice_tot)
    call gather_to_all( ice%mask_gl_gr       , mask_gl_gr_tot       )
    call gather_to_all( ice%mask_margin      , mask_margin_tot      )
    call gather_to_all( ice%fraction_gr      , fraction_gr_tot      )

    ! == Calculate bed roughness rates of changes
    ! ===========================================

    do vi = mesh%vi1, mesh%vi2

      ! Determine whether bed roughness should be
      ! updated by inversion or by extrapolation

      ! Only perform the inversion on fully grounded vertices
      if (ice%mask_grounded_ice( vi) .and. &
        .not. (ice%mask_margin( vi) .or. ice%mask_gl_gr( vi) .or. ice%mask_cf_gr( vi))) then

        ! Perform the inversion here
        mask( vi) = 2

        ! Surface elevation misfit
        misfit = ice%Hs( vi) - refgeo%Hs( vi)

        ! Is it improving already?
        if (ice%dHs_dt( vi)*misfit < 0._dp) then
          ! Yes, so leave this vertex alone
          cycle
        end if

      ELSE

        ! Extrapolate here
        mask( vi) = 1
        cycle

      end if

      ! Trace both halves of the flowline
      ! =================================

      ! The point p
      p = [mesh%V( vi,1), mesh%V( vi,2)]

      ! Trace both halves of the flowline
      call trace_flowline_upstream(   mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_up  , n_up  )
      call trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_down, n_down)

      ! if we couldn't trace the flowline here, extrapolate instead of inverting
      if (n_up < 3 .or. n_down < 3) then
        ! Mark for first extrapolation
        mask( vi) = 1
        ! Skip inversion and go to next vertex
        cycle
      end if

      ! Calculate distance along both halves of the flowline
      s_up = 0._dp
      do k = 2, n_up
        s_up( k) = s_up( k-1) + NORM2( trace_up( k,:) - trace_up( k-1,:))
      end do

      s_down = 0._dp
      do k = 2, n_down
        s_down( k) = s_down( k-1) + NORM2( trace_down( k,:) - trace_down( k-1,:))
      end do

      ! Calculate thickness error and thinning rates on both halves of the flowline
      ! ===========================================================================

      deltaHs_up = 0._dp
      dHs_dt_up  = 0._dp
      ti         = mesh%iTri( vi,1)

      do k = 1, n_up

        pt = trace_up( k,:)

        call interpolate_to_point_dp_2D_singlecore( mesh, Hs_tot       , pt, ti, Hs_mod)
        call interpolate_to_point_dp_2D_singlecore( mesh, Hs_target_tot, pt, ti, Hs_target)
        call interpolate_to_point_dp_2D_singlecore( mesh, dHs_dt_tot   , pt, ti, dHs_dt_mod)

        deltaHs_up( k) = Hs_mod - Hs_target
        dHs_dt_up(  k) = dHs_dt_mod

      end do

      deltaHs_down = 0._dp
      dHs_dt_down  = 0._dp
      ti           = mesh%iTri( vi,1)

      do k = 1, n_down

        pt = trace_down( k,:)

        call interpolate_to_point_dp_2D_singlecore( mesh, Hs_tot       , pt, ti, Hs_mod)
        call interpolate_to_point_dp_2D_singlecore( mesh, Hs_target_tot, pt, ti, Hs_target)
        call interpolate_to_point_dp_2D_singlecore( mesh, dHs_dt_tot   , pt, ti, dHs_dt_mod)

        deltaHs_down( k) = Hs_mod - Hs_target
        dHs_dt_down(  k) = dHs_dt_mod

      end do

      ! Calculate weighted average of thickness error and thinning rates on both halves of the flowline
      ! ===============================================================================================

      int_w_deltaHs_up = 0._dp
      int_w_dHs_dt_up  = 0._dp
      int_w_up         = 0._dp

      do k = 2, n_up

        ! Distance of both points
        s1 = s_up( k-1)
        s2 = s_up( k  )
        ds = s2 - s1

        ! Weights for both points
        w1 = (2._dp / s_up( n_up)) * (1._dp - s1 / s_up( n_up))
        w2 = (2._dp / s_up( n_up)) * (1._dp - s2 / s_up( n_up))
        w_av = (w1 + w2) / 2._dp

        ! Thickness error and thinning rate for both points
        deltaHs1   = deltaHs_up( k-1)
        deltaHs2   = deltaHs_up( k)
        deltaHs_av = (deltaHs1 + deltaHs2) / 2._dp
        dHs_dt1    = dHs_dt_up(  k-1)
        dHs_dt2    = dHs_dt_up(  k)
        dHs_dt_av  = (dHs_dt1 + dHs_dt2) / 2._dp

        ! Add to integrals
        int_w_deltaHs_up = int_w_deltaHs_up + (w_av * deltaHs_av * ds)
        int_w_dHs_dt_up  = int_w_dHs_dt_up  + (w_av * dHs_dt_av  * ds)
        int_w_up         = int_w_up         + (w_av              * ds)

      end do

      deltaHs_av_up( vi) = int_w_deltaHs_up / int_w_up
      dHs_dt_av_up(  vi) = int_w_dHs_dt_up  / int_w_up

      int_w_deltaHs_down = 0._dp
      int_w_dHs_dt_down  = 0._dp
      int_w_down         = 0._dp

      do k = 2, n_down

        ! Distance of both points
        s1 = s_down( k-1)
        s2 = s_down( k  )
        ds = s2 - s1

        ! Weights for both points
        w1 = (2._dp / s_down( n_down)) * (1._dp - s1 / s_down( n_down))
        w2 = (2._dp / s_down( n_down)) * (1._dp - s2 / s_down( n_down))
        w_av = (w1 + w2) / 2._dp

        ! Thickness error and thinning rate for both points
        deltaHs1   = deltaHs_down( k-1)
        deltaHs2   = deltaHs_down( k)
        deltaHs_av = (deltaHs1 + deltaHs2) / 2._dp
        dHs_dt1    = dHs_dt_down(  k-1)
        dHs_dt2    = dHs_dt_down(  k)
        dHs_dt_av  = (dHs_dt1 + dHs_dt2) / 2._dp

        ! Add to integrals
        int_w_deltaHs_down = int_w_deltaHs_down + (w_av * deltaHs_av * ds)
        int_w_dHs_dt_down  = int_w_dHs_dt_down  + (w_av * dHs_dt_av  * ds)
        int_w_down         = int_w_down         + (w_av              * ds)

      end do

      deltaHs_av_down( vi) = int_w_deltaHs_down / int_w_down
      dHs_dt_av_down(  vi) = int_w_dHs_dt_down  / int_w_down

      ! Calculate bed roughness rates of change
      ! =======================================

      R( vi) = max( 0._dp, min( 1._dp, &
        ((ice%uabs_vav( vi) * ice%Hi( vi)) / (C%bednudge_H_dHdt_flowline_u_scale * C%bednudge_H_dHdt_flowline_Hi_scale)) ))

      I_tot( vi) = R( vi) * (&
        (deltaHs_av_up( vi)                       ) / C%bednudge_H_dHdt_flowline_dH0 + &
        (dHs_dt_av_up(  vi) + dHs_dt_av_down(  vi)) / C%bednudge_H_dHdt_flowline_dHdt0)

      dC_dt( vi) = -1._dp * (I_tot( vi) * bed_roughness%generic_bed_roughness( vi)) / C%bednudge_H_dHdt_flowline_t_scale

    end do

    ! == Extrapolated inverted roughness rates of change to the whole domain
    ! ======================================================================

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, mask, dC_dt, C%bednudge_H_dHdt_flowline_r_smooth)

    ! Regularise tricky extrapolated areas

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

    ! Smoothing
    ! =========

    dC_dt_smoothed = dC_dt

    ! Smooth the local variable
    call smooth_Gaussian( mesh, grid_smooth, dC_dt_smoothed, C%bednudge_H_dHdt_flowline_r_smooth)

    do vi = mesh%vi1, mesh%vi2
      dC_dt( vi) = (1._dp - C%bednudge_H_dHdt_flowline_w_smooth) * dC_dt( vi) + C%bednudge_H_dHdt_flowline_w_smooth * dC_dt_smoothed( vi)
    end do ! do vi = mesh%vi1, mesh%vi2

    ! Final bed roughness field
    ! =========================

    ! Calculate predicted bed roughness at t+dt
    bed_roughness%generic_bed_roughness_next = MAX( C%generic_bed_roughness_min, MIN( C%generic_bed_roughness_max, &
      bed_roughness%generic_bed_roughness_prev + C%bed_roughness_nudging_dt * dC_dt ))

    ! Clean up after yourself
    deallocate( mask           )
    deallocate( trace_up       )
    deallocate( trace_down     )
    deallocate( s_up           )
    deallocate( s_down         )
    deallocate( deltaHs_up     )
    deallocate( deltaHs_down   )
    deallocate( dHs_dt_up      )
    deallocate( dHs_dt_down    )
    deallocate( deltaHs_av_up  )
    deallocate( deltaHs_av_down)
    deallocate( dHs_dt_av_up   )
    deallocate( dHs_dt_av_down )
    deallocate( R              )
    deallocate( I_tot          )
    deallocate( dC_dt          )
    deallocate( dHs_dx         )
    deallocate( dHs_dy         )
    deallocate( abs_grad_Hs    )
    deallocate( dC_dt_smoothed )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_bed_roughness_nudging_H_dHdt_flowline

  subroutine initialise_bed_roughness_nudging_H_dHdt_flowline( mesh, ice, bed_roughness, region_name)
    ! Initialise the bed roughness nudging model based on flowline-averaged values of H and dH/dt

    ! Input variables:
    type(type_mesh),                intent(in   ) :: mesh
    type(type_ice_model),           intent(in   ) :: ice
    type(type_bed_roughness_model), intent(inout) :: bed_roughness
    character(len=3),               intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_nudging_H_dHdt_flowline'
    real(dp)                       :: dummy_dp
    character                      :: dummy_char

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings
    dummy_dp = mesh%xmin
    dummy_dp = ice%Hi( mesh%vi1)
    dummy_dp = bed_roughness%generic_bed_roughness( mesh%vi1)
    dummy_char = region_name( 1:1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_nudging_H_dHdt_flowline

  ! == Flowline tracing
  ! ===================

  subroutine trace_flowline_upstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n)
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

    ! Local variables:
    real(dp), dimension(2) :: pt
    integer                :: vi, iti, ti
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
      if (uabs_pt < 1._dp) exit

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

  end subroutine trace_flowline_upstream

  subroutine trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n)
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
    real(dp), dimension(:,:  ),     intent(  out) :: T
    integer,                        intent(  out) :: n

    ! Local variables:
    real(dp), dimension(2) :: pt
    integer                :: vi, iti, ti
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
      if (Hi_tot( vi) < 1._dp) exit

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
      if (uabs_pt < 1._dp) exit

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

  end subroutine trace_flowline_downstream

end module bed_roughness_nudging_H_dHdt_flowline
