MODULE basal_inversion_H_dHdt_flowline

  ! Contains all the routines for the basal inversion model
  ! based on flowline-averaged values of H and dH/dt

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF, &
                                                                     save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D, &
                                                                     save_variable_as_netcdf_logical_1D
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE basal_inversion_types                                  , ONLY: type_basal_inversion
  USE mesh_utilities                                         , ONLY: find_containing_vertex, find_containing_triangle, extrapolate_Gaussian
  USE math_utilities                                         , ONLY: triangle_area
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D
  USE mesh_operators                                         , ONLY: ddx_a_a_2D, ddy_a_a_2D
  USE mesh_data_smoothing                                    , ONLY: smooth_Gaussian_2D


  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_basal_inversion_H_dHdt_flowline( mesh, grid_smooth, ice, refgeo, BIV)
    ! Run the basal inversion model based on flowline-averaged values of H and dH/dt

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_basal_inversion),          INTENT(INOUT) :: BIV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_basal_inversion_H_dHdt_flowline'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: mask
    REAL(dp), DIMENSION(mesh%nV)                       :: Hi_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: Hs_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: Hs_target_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: dHs_dt_tot
    REAL(dp), DIMENSION(mesh%nTri)                     :: u_b_tot
    REAL(dp), DIMENSION(mesh%nTri)                     :: v_b_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_grounded_ice_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_gl_gr_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_margin_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: fraction_gr_tot
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: p
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: trace_up, trace_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: s_up, s_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaHs_up, deltaHs_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHs_dt_up, dHs_dt_down
    INTEGER                                            :: n_up,n_down
    INTEGER                                            :: k
    REAL(dp), DIMENSION(2)                             :: pt
    INTEGER                                            :: ti,via,vib,vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    REAL(dp)                                           :: Atri_abp, Atri_bcp, Atri_cap, Atri_tot
    REAL(dp)                                           :: wa, wb, wc
    REAL(dp)                                           :: Hs_mod, Hs_target, dHs_dt_mod
    REAL(dp)                                           :: s1, s2, w1, w2, deltaHs1, deltaHs2, dHs_dt1, dHs_dt2, w_av, deltaHs_av, dHs_dt_av, ds
    REAL(dp)                                           :: int_w_deltaHs_up, int_w_dHs_dt_up, int_w_up
    REAL(dp)                                           :: int_w_deltaHs_down, int_w_dHs_dt_down, int_w_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaHs_av_up, deltaHs_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHs_dt_av_up, dHs_dt_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: I_tot, R
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt, dC2_dt
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHs_dx, dHs_dy, abs_grad_Hs
    REAL(dp)                                           :: fg_exp_mod
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt_smoothed, dC2_dt_smoothed
    REAL(dp)                                           :: misfit

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( mask(            mesh%vi1:mesh%vi2), source = 0      )
    ALLOCATE( trace_up(        mesh%nV, 2       ), source = 0._dp  )
    ALLOCATE( trace_down(      mesh%nV, 2       ), source = 0._dp  )
    ALLOCATE( s_up(            mesh%nV          ), source = 0._dp  )
    ALLOCATE( s_down(          mesh%nV          ), source = 0._dp  )
    ALLOCATE( deltaHs_up(      mesh%nV          ), source = 0._dp  )
    ALLOCATE( deltaHs_down(    mesh%nV          ), source = 0._dp  )
    ALLOCATE( dHs_dt_up(       mesh%nV          ), source = 0._dp  )
    ALLOCATE( dHs_dt_down(     mesh%nV          ), source = 0._dp  )
    ALLOCATE( deltaHs_av_up(   mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( deltaHs_av_down( mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHs_dt_av_up(    mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHs_dt_av_down(  mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( R(               mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( I_tot(           mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dC1_dt(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dC2_dt(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHs_dx(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHs_dy(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( abs_grad_Hs(     mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dC1_dt_smoothed( mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dC2_dt_smoothed( mesh%vi1:mesh%vi2), source = 0._dp  )

    ! Gather ice model data from all processes
    CALL gather_to_all_dp_1D(      ice%Hi               , Hi_tot               )
    CALL gather_to_all_dp_1D(      ice%Hs               , Hs_tot               )
    CALL gather_to_all_dp_1D(      refgeo%Hs            , Hs_target_tot        )
    CALL gather_to_all_dp_1D(      ice%dHs_dt           , dHs_dt_tot           )
    CALL gather_to_all_dp_1D(      ice%u_vav_b          , u_b_tot              )
    CALL gather_to_all_dp_1D(      ice%v_vav_b          , v_b_tot              )
    CALL gather_to_all_logical_1D( ice%mask_grounded_ice, mask_grounded_ice_tot)
    CALL gather_to_all_logical_1D( ice%mask_gl_gr       , mask_gl_gr_tot       )
    CALL gather_to_all_logical_1D( ice%mask_margin      , mask_margin_tot      )
    CALL gather_to_all_dp_1D(      ice%fraction_gr      , fraction_gr_tot      )

  ! == Calculate bed roughness rates of changes
  ! ===========================================

    DO vi = mesh%vi1, mesh%vi2

      ! Determine whether bed roughness should be
      ! updated by inversion or by extrapolation

      ! Only perform the inversion on fully grounded vertices
      IF (ice%mask_grounded_ice( vi) .AND. &
        .NOT. (ice%mask_margin( vi) .OR. ice%mask_gl_gr( vi) .OR. ice%mask_cf_gr( vi))) THEN

        ! Perform the inversion here
        mask( vi) = 2

        ! Surface elevation misfit
        misfit = ice%Hs( vi) - refgeo%Hs( vi)

        ! Is it improving already?
        IF (ice%dHs_dt( vi)*misfit < 0._dp) THEN
          ! Yes, so leave this vertex alone
          CYCLE
        END IF

      ELSE

        ! Extrapolate here
        mask( vi) = 1
        CYCLE

      END IF

      ! Trace both halves of the flowline
      ! =================================

      ! The point p
      p = [mesh%V( vi,1), mesh%V( vi,2)]

      ! Trace both halves of the flowline
      CALL trace_flowline_upstream(   mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_up  , n_up  )
      CALL trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_down, n_down)

      ! If we couldn't trace the flowline here, extrapolate instead of inverting
      IF (n_up < 3 .OR. n_down < 3) THEN
        ! Mark for first extrapolation
        mask( vi) = 1
        ! Skip inversion and go to next vertex
        CYCLE
      END IF

      ! Calculate distance along both halves of the flowline
      s_up = 0._dp
      DO k = 2, n_up
        s_up( k) = s_up( k-1) + NORM2( trace_up( k,:) - trace_up( k-1,:))
      END DO

      s_down = 0._dp
      DO k = 2, n_down
        s_down( k) = s_down( k-1) + NORM2( trace_down( k,:) - trace_down( k-1,:))
      END DO

      ! Calculate thickness error and thinning rates on both halves of the flowline
      ! ===========================================================================

      deltaHs_up = 0._dp
      dHs_dt_up  = 0._dp
      ti         = mesh%iTri( vi,1)

      DO k = 1, n_up

        ! The point along the flowline
        pt = trace_up( k,:)

        ! The mesh triangle containing the point
        CALL find_containing_triangle( mesh, pt, ti)

        ! The three vertices spanning ti
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        ! Trilinearly interpolate between a,b,c to find d_int
        pa = mesh%V( via,:)
        pb = mesh%V( vib,:)
        pc = mesh%V( vic,:)

        Atri_abp = triangle_area( pa, pb, p)
        Atri_bcp = triangle_area( pb, pc, p)
        Atri_cap = triangle_area( pc, pa, p)

        Atri_tot = Atri_abp + Atri_bcp + Atri_cap

        wc = Atri_abp / Atri_tot
        wa = Atri_bcp / Atri_tot
        wb = Atri_cap / Atri_tot

        Hs_mod     = Hs_tot(        via) * wa + Hs_tot(        vib) * wb + Hs_tot(        vic) * wc
        Hs_target  = Hs_target_tot( via) * wa + Hs_target_tot( vib) * wb + Hs_target_tot( vic) * wc
        dHs_dt_mod = dHs_dt_tot(    via) * wa + dHs_dt_tot(    vib) * wb + dHs_dt_tot(    vic) * wc

        deltaHs_up( k) = Hs_mod - Hs_target
        dHs_dt_up(  k) = dHs_dt_mod

      END DO !  DO k = 1, n_up

      deltaHs_down = 0._dp
      dHs_dt_down  = 0._dp
      ti           = mesh%iTri( vi,1)

      DO k = 1, n_down

        ! The point along the flowline
        pt = trace_down( k,:)

        ! The mesh triangle containing the point
        CALL find_containing_triangle( mesh, pt, ti)

        ! The three vertices spanning ti
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        ! Trilinearly interpolate between a,b,c to find d_int
        pa = mesh%V( via,:)
        pb = mesh%V( vib,:)
        pc = mesh%V( vic,:)

        Atri_abp = triangle_area( pa, pb, p)
        Atri_bcp = triangle_area( pb, pc, p)
        Atri_cap = triangle_area( pc, pa, p)

        Atri_tot = Atri_abp + Atri_bcp + Atri_cap

        wc = Atri_abp / Atri_tot
        wa = Atri_bcp / Atri_tot
        wb = Atri_cap / Atri_tot

        Hs_mod     = Hs_tot(        via) * wa + Hs_tot(        vib) * wb + Hs_tot(        vic) * wc
        Hs_target  = Hs_target_tot( via) * wa + Hs_target_tot( vib) * wb + Hs_target_tot( vic) * wc
        dHs_dt_mod = dHs_dt_tot(    via) * wa + dHs_dt_tot(    vib) * wb + dHs_dt_tot(    vic) * wc

        deltaHs_down( k) = Hs_mod - Hs_target
        dHs_dt_down(  k) = dHs_dt_mod

      END DO !  DO k = 1, n_down

      ! Calculate weighted average of thickness error and thinning rates on both halves of the flowline
      ! ===============================================================================================

      int_w_deltaHs_up = 0._dp
      int_w_dHs_dt_up  = 0._dp
      int_w_up         = 0._dp

      DO k = 2, n_up

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

      END DO ! DO k = 2, n_up

      deltaHs_av_up( vi) = int_w_deltaHs_up / int_w_up
      dHs_dt_av_up(  vi) = int_w_dHs_dt_up  / int_w_up

      int_w_deltaHs_down = 0._dp
      int_w_dHs_dt_down  = 0._dp
      int_w_down         = 0._dp

      DO k = 2, n_down

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

      END DO ! DO k = 2, n_down

      deltaHs_av_down( vi) = int_w_deltaHs_down / int_w_down
      dHs_dt_av_down(  vi) = int_w_dHs_dt_down  / int_w_down

      ! Calculate bed roughness rates of change
      ! =======================================

      R(     vi) = MAX( 0._dp, MIN( 1._dp, &
        ((ice%uabs_vav( vi) * ice%Hi( vi)) / (C%bednudge_H_dHdt_flowline_u_scale * C%bednudge_H_dHdt_flowline_Hi_scale)) ))

      I_tot( vi) = R( vi) * (&
        (deltaHs_av_up( vi)                       ) / C%bednudge_H_dHdt_flowline_dH0 + &
        (dHs_dt_av_up(  vi) + dHs_dt_av_down(  vi)) / C%bednudge_H_dHdt_flowline_dHdt0)

      dC1_dt( vi) = -1._dp * (I_tot( vi) * BIV%generic_bed_roughness_1( vi)) / C%bednudge_H_dHdt_flowline_t_scale
      dC2_dt( vi) = -1._dp * (I_tot( vi) * BIV%generic_bed_roughness_2( vi)) / C%bednudge_H_dHdt_flowline_t_scale

    END DO ! DO vi = mesh%vi1, mesh%vi2

  ! == Extrapolated inverted roughness rates of change to the whole domain
  ! ======================================================================

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    CALL extrapolate_Gaussian( mesh, mask, dC1_dt, C%bednudge_H_dHdt_flowline_r_smooth)
    CALL extrapolate_Gaussian( mesh, mask, dC2_dt, C%bednudge_H_dHdt_flowline_r_smooth)

    ! Regularise tricky extrapolated areas

    ! Calculate surface slopes
    CALL ddx_a_a_2D( mesh, ice%Hs, dHs_dx)
    CALL ddy_a_a_2D( mesh, ice%Hs, dHs_dy)

    ! Calculate absolute surface gradient
    abs_grad_Hs = SQRT( dHs_dx**2 + dHs_dy**2)

    ! Scale bed roughness rate of change for partially grounded, steep-sloped areas
    DO vi = mesh%vi1, mesh%vi2

      ! Ice margin and grounding lines
      IF (ice%mask_grounded_ice( vi)) THEN

        ! Strengthen the effect of grounded fractions for steep slopes
        fg_exp_mod = MIN( 1.0_dp, MAX( 0._dp, MAX( 0._dp, abs_grad_Hs( vi) - 0.02_dp) / (0.06_dp - 0.02_dp) ))

        ! Scale based on grounded fraction
        dC1_dt( vi) = dC1_dt( vi) * ice%fraction_gr( vi) ** (1._dp + fg_exp_mod)
        dC2_dt( vi) = dC2_dt( vi) * ice%fraction_gr( vi) ** (1._dp + fg_exp_mod)

      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Smoothing
    ! =========

    dC1_dt_smoothed = dC1_dt
    dC2_dt_smoothed = dC2_dt

    ! Smooth the local variable
    CALL smooth_Gaussian_2D( mesh, grid_smooth, dC1_dt_smoothed, C%bednudge_H_dHdt_flowline_r_smooth)
    CALL smooth_Gaussian_2D( mesh, grid_smooth, dC2_dt_smoothed, C%bednudge_H_dHdt_flowline_r_smooth)

    DO vi = mesh%vi1, mesh%vi2
      dC1_dt( vi) = (1._dp - C%bednudge_H_dHdt_flowline_w_smooth) * dC1_dt( vi) + C%bednudge_H_dHdt_flowline_w_smooth * dC1_dt_smoothed( vi)
      dC2_dt( vi) = (1._dp - C%bednudge_H_dHdt_flowline_w_smooth) * dC2_dt( vi) + C%bednudge_H_dHdt_flowline_w_smooth * dC2_dt_smoothed( vi)
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Final bed roughness field
    ! =========================

    ! Calculate predicted bed roughness at t+dt
    BIV%generic_bed_roughness_1_next = MAX( C%generic_bed_roughness_1_min, MIN( C%generic_bed_roughness_1_max, &
      BIV%generic_bed_roughness_1_prev + C%bed_roughness_nudging_dt * dC1_dt ))
    BIV%generic_bed_roughness_2_next = MAX( C%generic_bed_roughness_2_min, MIN( C%generic_bed_roughness_2_max, &
      BIV%generic_bed_roughness_2_prev + C%bed_roughness_nudging_dt * dC2_dt ))

    ! Clean up after yourself
    DEALLOCATE( mask           )
    DEALLOCATE( trace_up       )
    DEALLOCATE( trace_down     )
    DEALLOCATE( s_up           )
    DEALLOCATE( s_down         )
    DEALLOCATE( deltaHs_up     )
    DEALLOCATE( deltaHs_down   )
    DEALLOCATE( dHs_dt_up      )
    DEALLOCATE( dHs_dt_down    )
    DEALLOCATE( deltaHs_av_up  )
    DEALLOCATE( deltaHs_av_down)
    DEALLOCATE( dHs_dt_av_up   )
    DEALLOCATE( dHs_dt_av_down )
    DEALLOCATE( R              )
    DEALLOCATE( I_tot          )
    DEALLOCATE( dC1_dt         )
    DEALLOCATE( dC2_dt         )
    DEALLOCATE( dHs_dx         )
    DEALLOCATE( dHs_dy         )
    DEALLOCATE( abs_grad_Hs    )
    DEALLOCATE( dC1_dt_smoothed)
    DEALLOCATE( dC2_dt_smoothed)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_basal_inversion_H_dHdt_flowline

  SUBROUTINE initialise_basal_inversion_H_dHdt_flowline( mesh, ice, BIV, region_name)
    ! Initialise the basal inversion model based on flowline-averaged values of H and dH/dt

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_basal_inversion),          INTENT(INOUT) :: BIV
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion_H_dHdt_flowline'
    REAL(dp)                                           :: dummy_dp
    CHARACTER                                          :: dummy_char

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy_dp = mesh%xmin
    dummy_dp = ice%Hi( mesh%vi1)
    dummy_dp = BIV%generic_bed_roughness_1( mesh%vi1)
    dummy_char = region_name( 1:1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion_H_dHdt_flowline

  ! == Flowline tracing
  ! ===================

  SUBROUTINE trace_flowline_upstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n)
    ! Trace the flowline passing through point p upstream through
    ! the 2-D velocity field u_b,v_b.
    !
    ! Returns a list T of n points on the flowline
    !
    ! Stop the trace when it encounters the ice divide (defined as an ice velocity lower than 1 m/yr)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%nV),        INTENT(IN)    :: Hi_tot
    REAL(dp), DIMENSION(mesh%nTri),      INTENT(IN)    :: u_b_tot, v_b_tot
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: T
    INTEGER,                             INTENT(OUT)   :: n

    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: pt
    INTEGER                                            :: vi, iti, ti
    REAL(dp)                                           :: dist, w, w_tot
    REAL(dp)                                           :: u_pt, v_pt, uabs_pt
    REAL(dp), DIMENSION(2)                             :: u_hat_pt
    REAL(dp)                                           :: dist_prev

    ! Safety - if there's no ice, we can't do a trace
    vi = 1
    CALL find_containing_vertex( mesh, p, vi)
    IF (Hi_tot( vi) < 1._dp) THEN
      T = 0._dp
      n = 1
      T( 1,:) = p
      RETURN
    END IF

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p

    ! Safety
    dist_prev = 0._dp

    DO WHILE (.TRUE.)

      ! Find the vertex vi containing the tracer
      CALL find_containing_vertex( mesh, pt, vi)

      ! Interpolate between the surrounding triangles to find
      ! the velocities at the tracer's location

      w_tot = 0._dp
      u_pt  = 0._dp
      v_pt  = 0._dp

      DO iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        dist = NORM2( mesh%TriGC( ti,:) - pt)
        w = 1._dp / dist**2
        w_tot = w_tot + w
        u_pt  = u_pt  + w * u_b_tot( ti)
        v_pt  = v_pt  + w * v_b_tot( ti)
      END DO

      u_pt = u_pt / w_tot
      v_pt = v_pt / w_tot

      ! Calculate absolute velocity at the tracer's location
      uabs_pt = SQRT( u_pt**2 + v_pt**2)

      ! If we've reached the ice divide (defined as the place where
      ! we find velocities below 1 m/yr), end the trace
      IF (uabs_pt < 1._dp) EXIT

      ! Calculate the normalised velocity vector at the tracer's location
      u_hat_pt = [u_pt / uabs_pt, v_pt / uabs_pt]

      ! Add current position to the traces
      n = n + 1
      ! Safety
      IF (n > SIZE( T,1)) THEN
        ! DO iti = 1, MIN( SIZE(T,1), 1000)
        !   print*, T(iti,1), ',', T(iti,2)
        ! END DO
        CALL crash('upstream flowline tracer got stuck!')
      END IF
      T( n,:) = pt

      ! Save previous distance-to-origin
      dist_prev = NORM2( pt - p)

      ! Move the tracer upstream by a distance of one local resolution
      pt = pt - u_hat_pt * mesh%R( vi)

      ! If the new distance-to-origin is shorter than the previous one, end the trace
      IF (NORM2( pt - p) < dist_prev) EXIT

      ! If the new tracer location is outside the domain, end the trace
      IF (pt( 1) <= mesh%xmin .OR. pt( 2) >= mesh%xmax .OR. &
          pt( 2) <= mesh%ymin .OR. pt( 2) >= mesh%ymax) EXIT

    END DO ! DO WHILE (.TRUE.)

    ! Safety
    IF (n == 0) THEN
      n = 1
      T( 1,:) = p
    END IF

  END SUBROUTINE trace_flowline_upstream

  SUBROUTINE trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n)
    ! Trace the flowline passing through point p downstream through
    ! the 2-D velocity field u_b,v_b.
    !
    ! Returns a list T of n points on the flowline
    !
    ! Stop the trace when it encounters the ice margin

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%nV),        INTENT(IN)    :: Hi_tot
    REAL(dp), DIMENSION(mesh%nTri),      INTENT(IN)    :: u_b_tot, v_b_tot
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: T
    INTEGER,                             INTENT(OUT)   :: n

    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: pt
    INTEGER                                            :: vi, iti, ti
    REAL(dp)                                           :: dist, w, w_tot
    REAL(dp)                                           :: u_pt, v_pt, uabs_pt
    REAL(dp), DIMENSION(2)                             :: u_hat_pt
    REAL(dp)                                           :: dist_prev

    ! Safety - if there's no ice, we can't do a trace
    vi = 1
    CALL find_containing_vertex( mesh, p, vi)
    IF (Hi_tot( vi) < 1._dp) THEN
      T = 0._dp
      n = 1
      T( 1,:) = p
      RETURN
    END IF

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p

    ! Safety
    dist_prev = 0._dp

    DO WHILE (.TRUE.)

      ! Find the vertex vi containing the tracer
      CALL find_containing_vertex( mesh, pt, vi)

      ! If ice thickness in this vertex is below 1 m, assume we've found the
      ! ice margin, and end the trace
      IF (Hi_tot( vi) < 1._dp) EXIT

      ! Interpolate between the surrounding triangles to find
      ! the velocities at the tracer's location

      w_tot = 0._dp
      u_pt  = 0._dp
      v_pt  = 0._dp

      DO iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        dist = NORM2( mesh%TriGC( ti,:) - pt)
        IF (dist == 0._dp) CALL crash('whaa!')
        w = 1._dp / dist**2
        w_tot = w_tot + w
        u_pt  = u_pt  + w * u_b_tot( ti)
        v_pt  = v_pt  + w * v_b_tot( ti)
      END DO

      u_pt = u_pt / w_tot
      v_pt = v_pt / w_tot

      ! Calculate absolute velocity at the tracer's location
      uabs_pt = SQRT( u_pt**2 + v_pt**2)

      ! If we're at the ice divide (defined as the place where
      ! we find velocities below 1 m/yr), we can't do the trace
      IF (uabs_pt < 1._dp) EXIT

      ! Calculate the normalised velocity vector at the tracer's location
      u_hat_pt = [u_pt / uabs_pt, v_pt / uabs_pt]

      ! Add current position to the traces
      n = n + 1
      ! Safety
      IF (n > SIZE( T,1)) THEN
        ! DO iti = 1, MIN( SIZE(T,1), 1000)
        !   print*, T(iti,1), ',', T(iti,2)
        ! END DO
        CALL crash('downstream flowline tracer got stuck!')
      END IF
      T( n,:) = pt

      ! Save previous distance-to-origin
      dist_prev = NORM2( pt - p)

      ! Move the tracer downstream by a distance of one local resolution
      pt = pt + u_hat_pt * mesh%R( vi)

      ! If the new distance-to-origin is shorter than the previous one, end the trace
      IF (NORM2( pt - p) < dist_prev) EXIT

      ! If the new tracer location is outside the domain, end the trace
      IF (pt( 1) <= mesh%xmin .OR. pt( 2) >= mesh%xmax .OR. &
          pt( 2) <= mesh%ymin .OR. pt( 2) >= mesh%ymax) EXIT

    END DO ! DO WHILE (.TRUE.)

    ! Safety
    IF (n == 0) THEN
      n = 1
      T( 1,:) = p
    END IF

  END SUBROUTINE trace_flowline_downstream

END MODULE basal_inversion_H_dHdt_flowline
