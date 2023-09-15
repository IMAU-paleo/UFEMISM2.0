MODULE basal_hydrology

  ! Contains all the different basal hydrology models.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE basal_inversion_types                                  , ONLY: type_hydrology_inversion
  USE reallocate_mod                                         , ONLY: reallocate_clean_dp_1D
  USE mesh_utilities                                         , ONLY: find_containing_vertex, find_containing_triangle, extrapolate_Gaussian
  USE math_utilities                                         , ONLY: triangle_area, is_floating
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D
  USE mesh_remapping                                         , ONLY: smooth_Gaussian_2D
  USE mesh_operators                                         , ONLY: ddx_a_a_2D, ddy_a_a_2D

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_basal_hydrology_model( mesh, grid_smooth, ice, refgeo, HIV, time)
    ! Run the chosen basal hydrology model

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_hydrology_inversion),      INTENT(INOUT) :: HIV
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_basal_hydrology_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================

    IF     (C%choice_basal_hydrology_model == 'Martin2011') THEN
      ! The Martin et al. (2011) parameterisation of pore water pressure
      CALL calc_pore_water_pressure_Martin2011( mesh, ice)
    ELSEIF (C%choice_basal_hydrology_model == 'inversion') THEN
      ! Inversion of pore water pressure
      CALL run_pore_water_fraction_inversion( mesh, grid_smooth, ice, refgeo, HIV, time)
    ELSE
      CALL crash('unknown choice_basal_hydrology_model "' // TRIM( C%choice_basal_hydrology_model) // '"!')
    END IF

    ! Calculate overburden and effective pressure
    ! ===========================================

    DO vi = mesh%vi1, mesh%vi2
      ice%overburden_pressure( vi) = ice_density * grav * ice%Hi( vi)
      ice%effective_pressure(  vi) = MAX( 0._dp, ice%overburden_pressure( vi) - ice%pore_water_pressure( vi))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_basal_hydrology_model

  SUBROUTINE initialise_basal_hydrology_model( mesh, ice)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_hydrology_model'
    REAL(dp)                                           :: dummy1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy1 = mesh%xmin
    dummy1 = ice%Hi( mesh%vi1)

    ! Initialise the chosen basal hydrology model
    IF     (C%choice_basal_hydrology_model == 'Martin2011') THEN
      ! No need to do anything
    ELSEIF (C%choice_basal_hydrology_model == 'inversion') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_basal_hydrology_model "' // TRIM( C%choice_basal_hydrology_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_hydrology_model

! ===== Different basal hydrology models =====
! ===========================================

  ! == Martin et al. (2011)
  SUBROUTINE calc_pore_water_pressure_Martin2011( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Use the parameterisation from Martin et al. (2011)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_Martin2011'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      ice%pore_water_fraction( vi) = MIN( 1._dp, MAX( 0._dp, 1._dp - (ice%Hb( vi) - ice%SL( vi) - C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - C%Martin2011_hydro_Hb_min) ))

      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure( vi) = 0.96_dp * ice_density * grav * ice%Hi( vi) * ice%pore_water_fraction( vi)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_Martin2011

! ===== Inversion =====
! =====================

  ! == Main routines
  !=================

  SUBROUTINE run_pore_water_fraction_inversion( mesh, grid_smooth, ice, refgeo, HIV, time)
    ! Run the main hydrology inversion model

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_hydrology_inversion),      INTENT(INOUT) :: HIV
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_pore_water_fraction_inversion'
    INTEGER                                            :: vi
    REAL(dp)                                           :: wt_prev, wt_next

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Only do basal inversion within the specified time window
    IF (time < C%pore_water_nudging_t_start) THEN
      HIV%t_next = C%pore_water_nudging_t_start
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    IF (time == C%pore_water_nudging_t_end) THEN
      HIV%t_next = C%end_time_of_run
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! If the desired time is beyond the time of the next modelled pore water,
    ! run the basal inversion model to calculate a new next modelled pore water.
    ! ==========================================================================

    IF (time == HIV%t_next) THEN
      ! Need to calculate new predicted pore water

      ! Store previous modelled pore water
      HIV%pore_water_fraction_prev = HIV%pore_water_fraction_next
      HIV%t_prev = HIV%t_next
      HIV%t_next = HIV%t_prev + C%pore_water_nudging_dt

      call pore_water_fraction_inversion( mesh, grid_smooth, ice, refgeo, HIV)

    ELSEIF (time > HIV%t_next) THEN
      ! This should not be possible
      CALL crash('overshot the hydrology inversion time step')
    ELSE
      ! We're within the current HIV prediction window
    END IF ! IF (region%time == region%HIV%t_next) THEN

    ! Interpolate between previous and next modelled pore
    ! water fraction to find its value at the desired time
    ! ====================================================

    ! Calculate time interpolation weights
    wt_prev = (HIV%t_next - time) / (HIV%t_next - HIV%t_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled pore water fraction to desired time
    DO vi = mesh%vi1, mesh%vi2
      ice%pore_water_fraction( vi) = wt_prev * HIV%pore_water_fraction_prev( vi) + wt_next * HIV%pore_water_fraction_next( vi)
    END DO

    ! Compute effective pore water pressure
    ! =====================================

    DO vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure( vi) = ice%pore_water_fraction(vi) * ice_density * grav * ice%Hi( vi)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_pore_water_fraction_inversion

  SUBROUTINE initialise_pore_water_fraction_inversion( mesh, ice, HIV, region_name)
    ! Initialise pore water fraction inversion model

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_hydrology_inversion),      INTENT(OUT)   :: HIV
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_pore_water_fraction_inversion'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master) WRITE(0,*) '  Initialising basal hydrology inversion model...'

    ! Allocate memory for main variables
    ! ==================================

    ALLOCATE( HIV%pore_water_fraction_prev( mesh%vi1:mesh%vi2))
    ALLOCATE( HIV%pore_water_fraction_next( mesh%vi1:mesh%vi2))

    HIV%pore_water_fraction_prev = 0._dp
    HIV%pore_water_fraction_next = 0._dp

    ! Timeframes
    HIV%t_prev   = C%start_time_of_run
    HIV%t_next   = C%start_time_of_run

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_pore_water_fraction_inversion

  ! == Inversion methods
  ! ====================

  SUBROUTINE pore_water_fraction_inversion( mesh, grid_smooth, ice, refgeo, HIV)
    ! Invert the pore water fraction
    !
    ! Use the flowline approach to track potential for hydro stuff

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_hydrology_inversion),      INTENT(INOUT) :: HIV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'pore_water_fraction_inversion'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: pore_dryness
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: mask
    REAL(dp), DIMENSION(mesh%nV)                       :: Hi_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: Hi_target_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: dHi_dt_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: Ti_hom_tot
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
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaHi_up, deltaHi_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHi_dt_up, dHi_dt_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Ti_hom_up, Ti_hom_down
    INTEGER                                            :: n_up,n_down
    INTEGER                                            :: k
    REAL(dp), DIMENSION(2)                             :: pt
    INTEGER                                            :: ti,via,vib,vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    REAL(dp)                                           :: Atri_abp, Atri_bcp, Atri_cap, Atri_tot
    REAL(dp)                                           :: wa, wb, wc
    REAL(dp)                                           :: Hi_mod, Hi_target, dHi_dt_mod, Ti_hom_mod
    REAL(dp)                                           :: s1, s2, w1, w2, deltaHi1, deltaHi2, dHi_dt1, dHi_dt2, Ti_hom1, Ti_hom2, w_av, deltaHi_av, dHi_dt_av, Ti_hom_av, ds
    REAL(dp)                                           :: int_w_deltaHi_up, int_w_dHi_dt_up, int_w_Ti_hom_up, int_w_up
    REAL(dp)                                           :: int_w_deltaHi_down, int_w_dHi_dt_down, int_w_Ti_hom_down, int_w_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaHi_av_up, deltaHi_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHi_dt_av_up, dHi_dt_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Ti_hom_av_up, Ti_hom_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: I_tot, R
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt, dC2_dt
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHs_dx, dHs_dy, abs_grad_Hs
    REAL(dp)                                           :: fg_exp_mod, bf_exp_mod, basal_fric, hs_exp_mod, hi_exp_mod
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt_smoothed, dC2_dt_smoothed
    REAL(dp)                                           :: misfit
    REAL(dp)                                           :: is_nice, is_thick, is_hot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( pore_dryness(    mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( mask(            mesh%vi1:mesh%vi2), source = 0      )
    ALLOCATE( trace_up(        mesh%nV, 2       ), source = 0._dp  )
    ALLOCATE( trace_down(      mesh%nV, 2       ), source = 0._dp  )
    ALLOCATE( s_up(            mesh%nV          ), source = 0._dp  )
    ALLOCATE( s_down(          mesh%nV          ), source = 0._dp  )
    ALLOCATE( deltaHi_up(      mesh%nV          ), source = 0._dp  )
    ALLOCATE( deltaHi_down(    mesh%nV          ), source = 0._dp  )
    ALLOCATE( dHi_dt_up(       mesh%nV          ), source = 0._dp  )
    ALLOCATE( dHi_dt_down(     mesh%nV          ), source = 0._dp  )
    ALLOCATE( Ti_hom_up(       mesh%nV          ), source = 0._dp  )
    ALLOCATE( Ti_hom_down(     mesh%nV          ), source = 0._dp  )
    ALLOCATE( deltaHi_av_up(   mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( deltaHi_av_down( mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHi_dt_av_up(    mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHi_dt_av_down(  mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( Ti_hom_av_up(    mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( Ti_hom_av_down(  mesh%vi1:mesh%vi2), source = 0._dp  )
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
    CALL gather_to_all_dp_1D(      refgeo%Hi            , Hi_target_tot        )
    CALL gather_to_all_dp_1D(      ice%dHi_dt_predicted , dHi_dt_tot           )
    CALL gather_to_all_dp_1D(      ice%Ti_hom           , Ti_hom_tot           )
    CALL gather_to_all_dp_1D(      ice%u_vav_b          , u_b_tot              )
    CALL gather_to_all_dp_1D(      ice%v_vav_b          , v_b_tot              )
    CALL gather_to_all_logical_1D( ice%mask_grounded_ice, mask_grounded_ice_tot)
    CALL gather_to_all_logical_1D( ice%mask_gl_gr       , mask_gl_gr_tot       )
    CALL gather_to_all_logical_1D( ice%mask_margin      , mask_margin_tot      )
    CALL gather_to_all_dp_1D(      ice%fraction_gr      , fraction_gr_tot      )

    ! == Calculate pore water rates of changes
    ! ========================================

    ! Use the supplement of pore water fraction,
    ! as it scales better during the inversion
    pore_dryness = (1._dp - HIV%pore_water_fraction_prev)

    DO vi = mesh%vi1, mesh%vi2

      ! Determine whether pore water fraction should be
      ! updated by inversion or by extrapolation

      ! Only perform the inversion on fully grounded vertices
      IF (ice%mask_grounded_ice( vi) .AND. (.NOT. ice%mask_gl_gr( vi))) THEN
      ! IF (ice%mask_grounded_ice( vi)) THEN

        ! Perform the inversion here
        mask( vi) = 2

        ! Surface elevation misfit
        misfit = ice%Hi( vi) - refgeo%Hi( vi)

      ELSE

        ! Give non-fully grounded ice a default value (irrelevant)
        ice%pore_water_likelihood( vi) = 1._dp

        ! Extrapolate here
        mask( vi) = 1
        CYCLE

      END IF

      ! Trace both halves of the flowline
      ! =================================

      ! The point p
      p = [mesh%V( vi,1), mesh%V( vi,2)]

      ! Check whether we want a local- of flowline-type of inversion
      IF (    C%choice_pore_water_nudging_method == 'flowline') THEN

        ! Trace both halves of the flowline
        CALL trace_flowline_upstream(   mesh,      Hi_tot, u_b_tot, v_b_tot, p, trace_up  , n_up  )
        CALL trace_flowline_downstream( mesh, ice, Hi_tot, u_b_tot, v_b_tot, p, trace_down, n_down)

      ELSEIF (C%choice_pore_water_nudging_method == 'local') THEN

        ! Set these values so the check below always fails
        n_up   = 0
        n_down = 0

      ELSE
        CALL crash('unknown choice_pore_water_nudging_method "' // TRIM( C%choice_pore_water_nudging_method) // '"!')
      END IF

      ! If we couldn't trace the flowline here, perform a local adjustment
      IF (n_up < 3 .OR. n_down < 3) THEN

        ! Compute likelihood of local subglacial water
        ice%pore_water_likelihood( vi) = EXP(ice%Ti_hom( vi)/3._dp)

        ! Compute new adjustment for pore water fraction
        I_tot( vi) = ice%pore_water_likelihood( vi) * (misfit / C%porenudge_H_dHdt_flowline_dH0 + &
                                                       2._dp * (ice%dHi_dt_predicted( vi)) / C%porenudge_H_dHdt_flowline_dHdt0)

        ! Compute a rate of adjustment
        dC1_dt( vi) = -1._dp * (I_tot( vi) * pore_dryness( vi)) / C%porenudge_H_dHdt_flowline_t_scale

        ! We are done here, so go to next vertex
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

      deltaHi_up = 0._dp
      dHi_dt_up  = 0._dp
      Ti_hom_up  = 0._dp
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

        Hi_mod     = Hi_tot(        via) * wa + Hi_tot(        vib) * wb + Hi_tot(        vic) * wc
        Hi_target  = Hi_target_tot( via) * wa + Hi_target_tot( vib) * wb + Hi_target_tot( vic) * wc
        dHi_dt_mod = dHi_dt_tot(    via) * wa + dHi_dt_tot(    vib) * wb + dHi_dt_tot(    vic) * wc
        Ti_hom_mod = Ti_hom_tot(    via) * wa + Ti_hom_tot(    vib) * wb + Ti_hom_tot(    vic) * wc

        deltaHi_up( k) = Hi_mod - Hi_target
        dHi_dt_up(  k) = dHi_dt_mod
        Ti_hom_up(  k) = Ti_hom_mod


      END DO !  DO k = 1, n_up

      deltaHi_down = 0._dp
      dHi_dt_down  = 0._dp
      Ti_hom_down  = 0._dp
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

        Hi_mod     = Hi_tot(        via) * wa + Hi_tot(        vib) * wb + Hi_tot(        vic) * wc
        Hi_target  = Hi_target_tot( via) * wa + Hi_target_tot( vib) * wb + Hi_target_tot( vic) * wc
        dHi_dt_mod = dHi_dt_tot(    via) * wa + dHi_dt_tot(    vib) * wb + dHi_dt_tot(    vic) * wc
        Ti_hom_mod = Ti_hom_tot(    via) * wa + Ti_hom_tot(    vib) * wb + Ti_hom_tot(    vic) * wc

        deltaHi_down( k) = Hi_mod - Hi_target
        dHi_dt_down(  k) = dHi_dt_mod
        Ti_hom_down(  k) = Ti_hom_mod

      END DO !  DO k = 1, n_down

      ! Calculate weighted average of thickness error and thinning rates on both halves of the flowline
      ! ===============================================================================================

      int_w_deltaHi_up = 0._dp
      int_w_dHi_dt_up  = 0._dp
      int_w_Ti_hom_up  = 0._dp
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
        deltaHi1   = deltaHi_up( k-1)
        deltaHi2   = deltaHi_up( k)
        deltaHi_av = (deltaHi1 + deltaHi2) / 2._dp
        dHi_dt1    = dHi_dt_up(  k-1)
        dHi_dt2    = dHi_dt_up(  k)
        dHi_dt_av  = (dHi_dt1 + dHi_dt2) / 2._dp
        Ti_hom1    = Ti_hom_up(  k-1)
        Ti_hom2    = Ti_hom_up(  k)
        Ti_hom_av  = (Ti_hom1 + Ti_hom2) / 2._dp

        ! Add to integrals
        int_w_deltaHi_up = int_w_deltaHi_up + (w_av * deltaHi_av * ds)
        int_w_dHi_dt_up  = int_w_dHi_dt_up  + (w_av * dHi_dt_av  * ds)
        int_w_Ti_hom_up  = int_w_Ti_hom_up  + (w_av * Ti_hom_av  * ds)
        int_w_up         = int_w_up         + (w_av              * ds)

      END DO ! DO k = 2, n_up

      deltaHi_av_up( vi) = int_w_deltaHi_up / int_w_up
      dHi_dt_av_up(  vi) = int_w_dHi_dt_up  / int_w_up
      Ti_hom_av_up(  vi) = int_w_Ti_hom_up  / int_w_up

      int_w_deltaHi_down = 0._dp
      int_w_dHi_dt_down  = 0._dp
      int_w_Ti_hom_down  = 0._dp
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
        deltaHi1   = deltaHi_down( k-1)
        deltaHi2   = deltaHi_down( k)
        deltaHi_av = (deltaHi1 + deltaHi2) / 2._dp
        dHi_dt1    = dHi_dt_down(  k-1)
        dHi_dt2    = dHi_dt_down(  k)
        dHi_dt_av  = (dHi_dt1 + dHi_dt2) / 2._dp
        Ti_hom1    = Ti_hom_down(  k-1)
        Ti_hom2    = Ti_hom_down(  k)
        Ti_hom_av  = (Ti_hom1 + Ti_hom2) / 2._dp

        ! Add to integrals
        int_w_deltaHi_down = int_w_deltaHi_down + (w_av * deltaHi_av * ds)
        int_w_dHi_dt_down  = int_w_dHi_dt_down  + (w_av * dHi_dt_av  * ds)
        int_w_Ti_hom_down  = int_w_Ti_hom_down  + (w_av * Ti_hom_av  * ds)
        int_w_down         = int_w_down         + (w_av              * ds)

      END DO ! DO k = 2, n_down

      deltaHi_av_down( vi) = int_w_deltaHi_down / int_w_down
      dHi_dt_av_down(  vi) = int_w_dHi_dt_down  / int_w_down
      Ti_hom_av_down(  vi) = int_w_Ti_hom_down  / int_w_down

      ! Calculate pore water fraction rates of change
      ! =============================================

      ! Compute likelihood of subglacial water
      ice%pore_water_likelihood( vi) = EXP(Ti_hom_av_up( vi)/3._dp)

      ! Compute new adjustment for pore water fraction
      I_tot( vi) = ice%pore_water_likelihood( vi) * ((deltaHi_av_up( vi)                      ) / C%porenudge_H_dHdt_flowline_dH0 + &
                                                     (dHi_dt_av_up(  vi) + dHi_dt_av_down( vi)) / C%porenudge_H_dHdt_flowline_dHdt0)

      ! Compute a rate of adjustment
      dC1_dt( vi) = -1._dp * (I_tot( vi) * pore_dryness( vi)) / C%porenudge_H_dHdt_flowline_t_scale

    END DO ! vi = mesh%vi1, mesh%vi2

    ! Regularise tricky areas
    ! =======================

    ! Calculate surface slopes
    CALL ddx_a_a_2D( mesh, ice%Hs, dHs_dx)
    CALL ddy_a_a_2D( mesh, ice%Hs, dHs_dy)

    ! Calculate absolute surface gradient
    abs_grad_Hs = SQRT( dHs_dx**2 + dHs_dy**2)

    DO vi = mesh%vi1, mesh%vi2

        ! Grounded fractions
        fg_exp_mod = 1._dp - ice%fraction_gr( vi)**C%subgrid_friction_exponent

        ! Steep slopes
        hs_exp_mod = MIN( 1.0_dp, MAX( 0._dp, MAX( 0._dp, abs_grad_Hs( vi) - 0.003_dp) / (0.01_dp - 0.003_dp) ))

        ! Slippery regions
        bf_exp_mod = 1._dp - MAX( 0._dp, MIN( 1._dp, ice%basal_friction_coefficient( vi) / 10000._dp ))

        ! Ice thickness
        hi_exp_mod = MIN( 1.0_dp, MAX( 0._dp, ice%Hi( vi)/1000._dp))

        ! Prevent over-increase of sliding over partially grounded, steep-sloped areas
        IF (dC1_dt( vi) < 0._dp) THEN
          ! Scale based on friction-slope-slide-thickness modifiers
          dC1_dt( vi) = dC1_dt( vi) * (1._dp - .333_dp * (fg_exp_mod + bf_exp_mod + hs_exp_mod)) * hi_exp_mod
        END IF

    END DO

    ! Extrapolate over floating areas
    ! ===============================

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    CALL extrapolate_Gaussian( mesh, mask, dC1_dt, C%porenudge_H_dHdt_flowline_r_smooth)

    ! Smoothing
    ! =========

    IF (C%porenudge_H_dHdt_flowline_w_smooth > 0._dp) THEN

      dC1_dt_smoothed = dC1_dt

      ! Smooth the local variable
      CALL smooth_Gaussian_2D( mesh, grid_smooth, dC1_dt_smoothed, C%porenudge_H_dHdt_flowline_r_smooth)

      DO vi = mesh%vi1, mesh%vi2
        dC1_dt( vi) = (1._dp - C%porenudge_H_dHdt_flowline_w_smooth) * dC1_dt( vi) + C%porenudge_H_dHdt_flowline_w_smooth * dC1_dt_smoothed( vi)
      END DO

    END IF

    ! New pore dryness field
    ! ======================

    pore_dryness = pore_dryness + dC1_dt * C%pore_water_nudging_dt

    ! Final pore water fraction field
    ! ===============================

    DO vi = mesh%vi1, mesh%vi2
      HIV%pore_water_fraction_next( vi) = 1._dp - pore_dryness( vi)
      IF (.NOT. ice%mask_grounded_ice( vi)) THEN
        HIV%pore_water_fraction_next( vi) = MAX(HIV%pore_water_fraction_next( vi), 1._dp - ice%fraction_gr( vi)**3._dp)
      END IF
      HIV%pore_water_fraction_next( vi) = MIN( HIV%pore_water_fraction_next( vi), C%pore_water_fraction_max)
      HIV%pore_water_fraction_next( vi) = MAX( HIV%pore_water_fraction_next( vi), C%pore_water_fraction_min)
    END DO

    ! Clean up after yourself
    DEALLOCATE( pore_dryness   )
    DEALLOCATE( mask           )
    DEALLOCATE( trace_up       )
    DEALLOCATE( trace_down     )
    DEALLOCATE( s_up           )
    DEALLOCATE( s_down         )
    DEALLOCATE( deltaHi_up     )
    DEALLOCATE( deltaHi_down   )
    DEALLOCATE( dHi_dt_up      )
    DEALLOCATE( dHi_dt_down    )
    DEALLOCATE( Ti_hom_up      )
    DEALLOCATE( Ti_hom_down    )
    DEALLOCATE( deltaHi_av_up  )
    DEALLOCATE( deltaHi_av_down)
    DEALLOCATE( dHi_dt_av_up   )
    DEALLOCATE( dHi_dt_av_down )
    DEALLOCATE( Ti_hom_av_up   )
    DEALLOCATE( Ti_hom_av_down )
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

  END SUBROUTINE pore_water_fraction_inversion

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
    REAL(dp)                                           :: dist_prev, dist_tot

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
    dist_tot  = 0._dp

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

      ! Add moved distance to total
      dist_tot  = dist_tot + NORM2( pt - p)

      ! If the new distance-to-origin is shorter than the previous one, end the trace
      IF (NORM2( pt - p) < dist_prev) EXIT

      ! If the total distance-to-origin is larger than a chosen limit, end the trace
      IF (dist_tot > C%porenudge_H_dHdt_flowline_dist_max * 1e3_dp) EXIT

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

  SUBROUTINE trace_flowline_downstream( mesh, ice, Hi_tot, u_b_tot, v_b_tot, p, T, n)
    ! Trace the flowline passing through point p downstream through
    ! the 2-D velocity field u_b,v_b.
    !
    ! Returns a list T of n points on the flowline
    !
    ! Stop the trace when it encounters the ice margin

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
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
    REAL(dp)                                           :: dist_prev, dist_tot

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
    dist_tot  = 0._dp

    DO WHILE (.TRUE.)

      ! Find the vertex vi containing the tracer
      CALL find_containing_vertex( mesh, pt, vi)

      ! If ice thickness in this vertex is below 1 m, assume we've found the
      ! ice margin, and end the trace
      IF (Hi_tot( vi) < 1._dp) EXIT

      ! If this vertex is floating ice, assume that we are keeping
      ! those areas fixed. End trace to avoid no-error biases
      IF (ice%mask_floating_ice( vi)) EXIT

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

      ! Add moved distance to total
      dist_tot  = dist_tot + NORM2( pt - p)

      ! If the new distance-to-origin is shorter than the previous one, end the trace
      IF (NORM2( pt - p) < dist_prev) EXIT

      ! If the total distance-to-origin is larger than a chosen limit, end the trace
      IF (dist_tot > C%porenudge_H_dHdt_flowline_dist_max * 1e3_dp) EXIT

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

END MODULE basal_hydrology
