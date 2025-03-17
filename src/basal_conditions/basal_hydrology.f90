MODULE basal_hydrology

  ! Contains all the different basal hydrology models.

! ===== Preamble =====
! ====================

  use mpi_basic, only: par
  USE precisions                                             , ONLY: dp
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE basal_inversion_types                                  , ONLY: type_hydrology_inversion
  USE mesh_utilities                                         , ONLY: find_containing_vertex, find_containing_triangle, extrapolate_Gaussian
  use ice_geometry_basics, only: is_floating
  use plane_geometry, only: triangle_area
  use mpi_distributed_memory, only: gather_to_all
  use mesh_data_smoothing, only: smooth_Gaussian
  use netcdf_io_main

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
    TYPE(type_hydrology_inversion),      INTENT(IN)    :: HIV
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_basal_hydrology_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================

    IF     (C%choice_basal_hydrology_model == 'none') THEN
      ! No pore water pressure
      CALL calc_pore_water_pressure_none( mesh, ice)
    ELSEIF (C%choice_basal_hydrology_model == 'Martin2011') THEN
      ! The Martin et al. (2011) parameterisation of pore water pressure
      CALL calc_pore_water_pressure_Martin2011( mesh, ice)
    ELSEIF (C%choice_basal_hydrology_model == 'inversion') THEN
      ! Inversion of pore water pressure
      CALL calc_pore_water_pressure_inversion( mesh, ice, HIV)
    ELSEIF (C%choice_basal_hydrology_model == 'read_from_file') THEN
      ! Read from file
      CALL calc_pore_water_pressure_from_file( mesh, ice)
    ELSE
      CALL crash('unknown choice_basal_hydrology_model "' // TRIM( C%choice_basal_hydrology_model) // '"!')
    END IF

    ! Calculate overburden and effective pressure
    ! ===========================================

    DO vi = mesh%vi1, mesh%vi2
      ice%overburden_pressure( vi) = ice_density * grav * ice%Hi_eff( vi)
      ice%effective_pressure(  vi) = MAX( 0._dp, ice%overburden_pressure( vi) - ice%pore_water_pressure( vi))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_basal_hydrology_model

  SUBROUTINE initialise_basal_hydrology_model( mesh, ice, region_name)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_hydrology_model'
    REAL(dp)                                           :: dummy1
    CHARACTER(LEN=256)                                 :: pore_water_fraction_choice_initialise

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy1 = mesh%xmin
    dummy1 = ice%Hi( mesh%vi1)

    ! Determine choice of initial ice temperatures for this model region
    IF     (region_name == 'NAM') THEN
      pore_water_fraction_choice_initialise = C%pore_water_fraction_choice_initialise_NAM
    ELSEIF (region_name == 'EAS') THEN
      pore_water_fraction_choice_initialise = C%pore_water_fraction_choice_initialise_EAS
    ELSEIF (region_name == 'GRL') THEN
      pore_water_fraction_choice_initialise = C%pore_water_fraction_choice_initialise_GRL
    ELSEIF (region_name == 'ANT') THEN
      pore_water_fraction_choice_initialise = C%pore_water_fraction_choice_initialise_ANT
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    ! Initialise the chosen basal hydrology model
    IF     (pore_water_fraction_choice_initialise == 'zero') THEN
      ! Set values to zero
      ice%pore_water_fraction = 0._dp
    ELSEIF (pore_water_fraction_choice_initialise == 'read_from_file') THEN
      ! Read pore water fraction from file
      CALL initialise_pore_water_fraction_from_file( mesh, ice, region_name)
    ELSE
      CALL crash('unknown choice_basal_hydrology_model "' // TRIM( C%choice_basal_hydrology_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_hydrology_model

! ===== Different basal hydrology models ====
! ===========================================

  ! == No subglacial hydrology
  SUBROUTINE calc_pore_water_pressure_none( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Use the parameterisation from Martin et al. (2011)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_none'
    INTEGER                                            :: vi
    REAL(dp)                                           :: weight_gr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute effective pore water pressure
    ! =====================================

    ! Scale pore water fraction based on grounded area fractions
    CALL apply_grounded_fractions_to_pore_water_fraction( mesh, ice)

    ! Compute pore water pressure based on the pore water fraction as
    ! the fraction of the overburden pressure supported by basal water
    DO vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure( vi) = ice%pore_water_fraction(vi) * ice_density * grav * ice%Hi_eff( vi)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_none

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
      ice%pore_water_pressure( vi) = 0.96_dp * ice_density * grav * ice%Hi_eff( vi) * ice%pore_water_fraction( vi)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_Martin2011

  ! == Inverted values
  SUBROUTINE calc_pore_water_pressure_inversion( mesh, ice, HIV)
    ! Calculate the pore water pressure
    !
    ! Use the inverted values of pore water fraction

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_hydrology_inversion),      INTENT(IN)    :: HIV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_inversion'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Extrapolate inverted pore water fraction
    ! ========================================

    CALL apply_extrapolation_to_pore_water_fraction( mesh, ice, HIV)

    ! Apply ocean entrainment
    ! =======================

    CALL apply_grounded_fractions_to_pore_water_fraction( mesh, ice)

    ! Compute pore water pressure
    ! ===========================

    DO vi = mesh%vi1, mesh%vi2
      ! Compute pore water pressure based on the pore water fraction as
      ! the fraction of the overburden pressure supported by basal water
      ice%pore_water_pressure( vi) = ice%pore_water_fraction(vi) * ice_density * grav * ice%Hi_eff( vi)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_inversion

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

      call pore_water_fraction_inversion( mesh, grid_smooth, ice, refgeo, HIV, time)

    ELSEIF (time > HIV%t_next) THEN
      ! This should not be possible
      CALL crash('overshot the hydrology inversion time step')
    ELSE
      ! We're within the current HIV prediction window
    END IF

    ! Interpolate between previous and next modelled pore
    ! water fraction to find its value at the desired time
    ! ====================================================

    ! Calculate time interpolation weights
    wt_prev = (HIV%t_next - time) / (HIV%t_next - HIV%t_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled pore water fraction to desired time
    DO vi = mesh%vi1, mesh%vi2
      HIV%pore_water_fraction_app( vi) = wt_prev * HIV%pore_water_fraction_prev( vi) + wt_next * HIV%pore_water_fraction_next( vi)
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
    IF (par%primary) WRITE(0,*) '  Initialising basal hydrology inversion model...'

    ! Allocate memory for main variables
    ! ==================================

    ALLOCATE( HIV%pore_water_fraction_prev( mesh%vi1:mesh%vi2))
    ALLOCATE( HIV%pore_water_fraction_next( mesh%vi1:mesh%vi2))
    ALLOCATE( HIV%pore_water_fraction_app(  mesh%vi1:mesh%vi2))
    ALLOCATE( HIV%mask_inverted_point(      mesh%vi1:mesh%vi2))

    HIV%pore_water_fraction_prev = ice%pore_water_fraction
    HIV%pore_water_fraction_next = ice%pore_water_fraction
    HIV%pore_water_fraction_app  = ice%pore_water_fraction
    HIV%mask_inverted_point  = .FALSE.

    ! Timeframes
    HIV%t_prev   = C%start_time_of_run
    HIV%t_next   = C%start_time_of_run

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_pore_water_fraction_inversion

  ! == Inversion methods
  ! ====================

  SUBROUTINE pore_water_fraction_inversion( mesh, grid_smooth, ice, refgeo, HIV, time)
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
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'pore_water_fraction_inversion'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: pore_dryness
    REAL(dp), DIMENSION(mesh%nV)                       :: Hi_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: Hi_target_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: dHi_dt_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: U_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: U_target_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: Ti_hom_tot
    REAL(dp), DIMENSION(mesh%nTri)                     :: u_b_tot
    REAL(dp), DIMENSION(mesh%nTri)                     :: v_b_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_grounded_ice_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_gl_gr_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_cf_gr_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_margin_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: fraction_gr_tot
    REAL(dp), DIMENSION(mesh%nV)                       :: pore_water_fraction_next_tot
    INTEGER                                            :: vi,ci,vc
    REAL(dp), DIMENSION(2)                             :: p
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: trace_up, trace_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: s_up, s_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaHi_up, deltaHi_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHi_dt_up, dHi_dt_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaU_up, deltaU_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Ti_hom_up, Ti_hom_down
    INTEGER                                            :: n_up,n_down
    INTEGER                                            :: k
    REAL(dp), DIMENSION(2)                             :: pt
    INTEGER                                            :: ti,via,vib,vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    REAL(dp)                                           :: Atri_abp, Atri_bcp, Atri_cap, Atri_tot
    REAL(dp)                                           :: wa, wb, wc
    REAL(dp)                                           :: Hi_mod, Hi_target, dHi_dt_mod, U_mod, U_target, Ti_hom_mod
    REAL(dp)                                           :: s1, s2, w1, w2, deltaHi1, deltaHi2, dHi_dt1, dHi_dt2, deltaU1, deltaU2, Ti_hom1, Ti_hom2, w_av, deltaHi_av, dHi_dt_av, deltaU_av, Ti_hom_av, ds
    REAL(dp)                                           :: int_w_deltaHi_up,   int_w_dHi_dt_up,   int_w_deltaU_up,   int_w_Ti_hom_up,   int_w_up
    REAL(dp)                                           :: int_w_deltaHi_down, int_w_dHi_dt_down, int_w_deltaU_down, int_w_Ti_hom_down, int_w_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaHi_av_up, deltaHi_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHi_dt_av_up, dHi_dt_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: deltaU_av_up, deltaU_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Ti_hom_av_up, Ti_hom_av_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: I_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt
    REAL(dp)                                           :: Hi_misfit, uabs_surf_misfit
    REAL(dp)                                           :: fg_exp_mod, bf_exp_mod, hs_exp_mod, hi_exp_mod, max_neighbour, max_vertex_size, unstable_vertex, exponent_gr
    REAL(dp)                                           :: t_scale, porenudge_H_dHdt_flowline_t_scale, porenudge_H_dHdt_flowline_dHdt0, porenudge_H_dHdt_flowline_dH0, porenudge_H_dHdt_flowline_dU0
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt_smoothed, unstable_vertex_smoothed
    LOGICAL                                            :: found_grounded_neighbour
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( pore_dryness(             mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( trace_up(                 mesh%nV, 2       ), source = 0._dp )
    ALLOCATE( trace_down(               mesh%nV, 2       ), source = 0._dp )
    ALLOCATE( s_up(                     mesh%nV          ), source = 0._dp )
    ALLOCATE( s_down(                   mesh%nV          ), source = 0._dp )
    ALLOCATE( deltaHi_up(               mesh%nV          ), source = 0._dp )
    ALLOCATE( deltaHi_down(             mesh%nV          ), source = 0._dp )
    ALLOCATE( dHi_dt_up(                mesh%nV          ), source = 0._dp )
    ALLOCATE( dHi_dt_down(              mesh%nV          ), source = 0._dp )
    ALLOCATE( deltaU_up(                mesh%nV          ), source = 0._dp )
    ALLOCATE( deltaU_down(              mesh%nV          ), source = 0._dp )
    ALLOCATE( Ti_hom_up(                mesh%nV          ), source = 0._dp )
    ALLOCATE( Ti_hom_down(              mesh%nV          ), source = 0._dp )
    ALLOCATE( deltaHi_av_up(            mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( deltaHi_av_down(          mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( dHi_dt_av_up(             mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( dHi_dt_av_down(           mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( deltaU_av_up(             mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( deltaU_av_down(           mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( Ti_hom_av_up(             mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( Ti_hom_av_down(           mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( I_tot(                    mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( dC1_dt(                   mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( dC1_dt_smoothed(          mesh%vi1:mesh%vi2), source = 0._dp )
    ALLOCATE( unstable_vertex_smoothed( mesh%vi1:mesh%vi2), source = 0._dp )

    ! Gather ice model data from all processes
    CALL gather_to_all(      ice%Hi               , Hi_tot               )
    CALL gather_to_all(      refgeo%Hi            , Hi_target_tot        )
    CALL gather_to_all(      ice%dHi_dt           , dHi_dt_tot           )
    CALL gather_to_all(      ice%uabs_surf        , U_tot                )
    CALL gather_to_all(      ice%uabs_surf_target , U_target_tot         )
    CALL gather_to_all(      ice%Ti_hom           , Ti_hom_tot           )
    CALL gather_to_all(      ice%u_vav_b          , u_b_tot              )
    CALL gather_to_all(      ice%v_vav_b          , v_b_tot              )
    CALL gather_to_all( ice%mask_grounded_ice, mask_grounded_ice_tot)
    CALL gather_to_all( ice%mask_gl_gr       , mask_gl_gr_tot       )
    CALL gather_to_all( ice%mask_cf_gr       , mask_cf_gr_tot       )
    CALL gather_to_all( ice%mask_margin      , mask_margin_tot      )
    CALL gather_to_all(      ice%fraction_gr      , fraction_gr_tot      )

    ! == Reducing power
    ! =================

    porenudge_H_dHdt_flowline_t_scale = C%porenudge_H_dHdt_flowline_t_scale
    porenudge_H_dHdt_flowline_dHdt0   = C%porenudge_H_dHdt_flowline_dHdt0
    porenudge_H_dHdt_flowline_dH0     = C%porenudge_H_dHdt_flowline_dH0
    porenudge_H_dHdt_flowline_dU0     = C%porenudge_H_dHdt_flowline_dU0

    ! Increase time length scale if so desired
    ! Compute how much time has passed since start of equilibrium stage
    t_scale = (time - 1.0E4_dp) / (C%pore_water_nudging_t_end - 1.0E4_dp)
    ! Limit t_scale to [0 1]
    t_scale = max( 0._dp, min( t_scale, 1._dp))
    ! Curve t_scale a bit
    t_scale = t_scale ** 1._dp

    porenudge_H_dHdt_flowline_t_scale = t_scale * 2._dp * C%porenudge_H_dHdt_flowline_t_scale + (1._dp - t_scale) * C%porenudge_H_dHdt_flowline_t_scale

    ! == Calculate pore water rates of change
    ! =======================================

    ! Use the supplement of pore water fraction,
    ! as it scales better during the inversion
    pore_dryness = (1._dp - HIV%pore_water_fraction_prev)

    ! Initialise mask of inverted values
    HIV%mask_inverted_point = .FALSE.

    DO vi = mesh%vi1, mesh%vi2

      ! Ice thickness misfit
      Hi_misfit = ice%Hi( vi) - refgeo%Hi( vi)

      ! Surface velocity misfit
      IF (C%do_target_uabs_surf .AND. ice%uabs_surf_target( vi) > 0._dp) THEN
        uabs_surf_misfit = ice%uabs_surf_target( vi) - ice%uabs_surf( vi)
      ELSE
        uabs_surf_misfit = 0._dp
      END IF

      ! Only perform the inversion on fully grounded vertices
      IF ( .NOT. ice%mask_grounded_ice( vi) .OR. ice%mask_gl_gr( vi) .OR. ice%mask_cf_gr( vi) .OR. ice%mask_margin( vi)) THEN
        ! Give non-fully grounded ice a default value (irrelevant)
        ice%pore_water_likelihood( vi) = 0._dp
        ! Skip this vertex
        CYCLE
      END IF

      ! At this point, this is a valid inversion vertex
      HIV%mask_inverted_point( vi) = .TRUE.

      ! Skip non-evolving and already fit-enough vertices
      IF (ABS(ice%dHi_dt( vi)) <= 1._dp .AND. ABS(Hi_misfit) <= 10._dp) THEN
        ! Skip this vertex
        CYCLE
      END IF

      ! Check if velocity misfit is worth it
      IF (ABS(uabs_surf_misfit) < 10._dp) THEN
        ! Ignore this level of misfit
        uabs_surf_misfit = 0._dp
      END IF

      ! Trace both halves of the flowline
      ! =================================

      ! The point p
      p = [mesh%V( vi,1), mesh%V( vi,2)]

      ! Check whether we want a local-, flowline-, or mixed-type inversion
      IF (    C%choice_pore_water_nudging_method == 'flowline') THEN

        ! Trace both halves of the flowline
        CALL trace_flowline_upstream(   mesh,      Hi_tot, u_b_tot, v_b_tot, p, trace_up  , n_up  )
        CALL trace_flowline_downstream( mesh, ice, Hi_tot, u_b_tot, v_b_tot, p, trace_down, n_down)

      ELSEIF (C%choice_pore_water_nudging_method == 'mixed') THEN

        IF (Hi_misfit < 0._dp) THEN
          ! If misfit is negative, trace upstream half of the flowline
          CALL trace_flowline_upstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, trace_up, n_up)
          ! Downstream half is not used, so assign dummy values to save time
          trace_down = trace_up
          n_down = n_up
        ELSE
          ! Set these values so the check below always fails for positive misfits
          n_up   = 0
          n_down = 0
        END IF

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
        I_tot( vi) = ice%pore_water_likelihood( vi) * (Hi_misfit / porenudge_H_dHdt_flowline_dH0 + &
                                                       2._dp * (ice%dHi_dt( vi)) / porenudge_H_dHdt_flowline_dHdt0 + &
                                                       uabs_surf_misfit / porenudge_H_dHdt_flowline_dU0)

        ! Compute a rate of adjustment
        dC1_dt( vi) = -1._dp * (I_tot( vi) * pore_dryness( vi)) / porenudge_H_dHdt_flowline_t_scale

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
      deltaU_up  = 0._dp
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
        U_mod      = U_tot (        via) * wa + U_tot(         vib) * wb + U_tot(         vic) * wc
        U_target   = U_target_tot(  via) * wa + U_target_tot(  vib) * wb + U_target_tot(  vic) * wc
        Ti_hom_mod = Ti_hom_tot(    via) * wa + Ti_hom_tot(    vib) * wb + Ti_hom_tot(    vic) * wc

        deltaHi_up( k) = Hi_mod - Hi_target
        dHi_dt_up(  k) = dHi_dt_mod
        deltaU_up( k) = U_mod - U_target
        Ti_hom_up(  k) = Ti_hom_mod

        ! Ignore invalid target values
        IF (U_target <= 0._dp .OR. U_target_tot( via) <= 0._dp .OR. U_target_tot( vib) <= 0._dp .OR. U_target_tot( vic) <= 0._dp) THEN
          deltaU_up( k) = 0._dp
        END IF

        ! Ignore small misfits
        IF (ABS(deltaU_up( k)) < 10._dp) THEN
          deltaU_up( k) = 0._dp
        END IF

      END DO !  DO k = 1, n_up

      deltaHi_down = 0._dp
      dHi_dt_down  = 0._dp
      deltaU_down  = 0._dp
      Ti_hom_down  = 0._dp
      ti           = mesh%iTri( vi,1)

      DO k = 1, n_down

        ! Skip downstream flowline when using the "mixed" method
        IF (C%choice_pore_water_nudging_method == 'mixed') CYCLE

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
        U_mod      = U_tot(         via) * wa + U_tot(         vib) * wb + U_tot(         vic) * wc
        U_target   = U_target_tot(  via) * wa + U_target_tot(  vib) * wb + U_target_tot(  vic) * wc
        Ti_hom_mod = Ti_hom_tot(    via) * wa + Ti_hom_tot(    vib) * wb + Ti_hom_tot(    vic) * wc

        deltaHi_down( k) = Hi_mod - Hi_target
        dHi_dt_down(  k) = dHi_dt_mod
        deltaU_down(  k) = U_mod - U_target
        Ti_hom_down(  k) = Ti_hom_mod

        ! Ignore invalid target values
        IF (U_target <= 0._dp .OR. U_target_tot( via) <= 0._dp .OR. U_target_tot( vib) <= 0._dp .OR. U_target_tot( vic) <= 0._dp) THEN
          deltaU_down( k) = 0._dp
        END IF

        ! Ignore small misfits
        IF (ABS(deltaU_down( k)) < 10._dp) THEN
          deltaU_down( k) = 0._dp
        END IF

      END DO !  DO k = 1, n_down

      ! Calculate weighted average of thickness error and thinning rates on both halves of the flowline
      ! ===============================================================================================

      int_w_deltaHi_up = 0._dp
      int_w_dHi_dt_up  = 0._dp
      int_w_deltaU_up  = 0._dp
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
        deltaU1    = deltaU_up( k-1)
        deltaU2    = deltaU_up( k)
        deltaU_av  = (deltaU1 + deltaU2) / 2._dp
        Ti_hom1    = Ti_hom_up(  k-1)
        Ti_hom2    = Ti_hom_up(  k)
        Ti_hom_av  = (Ti_hom1 + Ti_hom2) / 2._dp

        ! Add to integrals
        int_w_deltaHi_up = int_w_deltaHi_up + (w_av * deltaHi_av * ds)
        int_w_dHi_dt_up  = int_w_dHi_dt_up  + (w_av * dHi_dt_av  * ds)
        int_w_deltaU_up  = int_w_deltaU_up  + (w_av * deltaU_av  * ds)
        int_w_Ti_hom_up  = int_w_Ti_hom_up  + (w_av * Ti_hom_av  * ds)
        int_w_up         = int_w_up         + (w_av              * ds)

      END DO ! DO k = 2, n_up

      deltaHi_av_up( vi) = int_w_deltaHi_up / int_w_up
      dHi_dt_av_up(  vi) = int_w_dHi_dt_up  / int_w_up
      deltaU_av_up(  vi) = int_w_deltaU_up  / int_w_up
      Ti_hom_av_up(  vi) = int_w_Ti_hom_up  / int_w_up

      int_w_deltaHi_down = 0._dp
      int_w_dHi_dt_down  = 0._dp
      int_w_deltaU_down  = 0._dp
      int_w_Ti_hom_down  = 0._dp
      int_w_down         = 0._dp

      DO k = 2, n_down

        ! Skip downstream flowline when using the "mixed" method
        IF (C%choice_pore_water_nudging_method == 'mixed') CYCLE

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
        deltaU1    = deltaU_down( k-1)
        deltaU2    = deltaU_down( k)
        deltaU_av  = (deltaU1 + deltaU2) / 2._dp
        Ti_hom1    = Ti_hom_down(  k-1)
        Ti_hom2    = Ti_hom_down(  k)
        Ti_hom_av  = (Ti_hom1 + Ti_hom2) / 2._dp

        ! Add to integrals
        int_w_deltaHi_down = int_w_deltaHi_down + (w_av * deltaHi_av * ds)
        int_w_dHi_dt_down  = int_w_dHi_dt_down  + (w_av * dHi_dt_av  * ds)
        int_w_deltaU_down  = int_w_deltaU_down  + (w_av * deltaU_av  * ds)
        int_w_Ti_hom_down  = int_w_Ti_hom_down  + (w_av * Ti_hom_av  * ds)
        int_w_down         = int_w_down         + (w_av              * ds)

      END DO ! DO k = 2, n_down

      IF (C%choice_pore_water_nudging_method == 'mixed') THEN
        deltaHi_av_down( vi) = 0._dp
        dHi_dt_av_down(  vi) = 0._dp
        deltaU_av_down(  vi) = 0._dp
        Ti_hom_av_down(  vi) = 0._dp
      ELSE
        deltaHi_av_down( vi) = int_w_deltaHi_down / int_w_down
        dHi_dt_av_down(  vi) = int_w_dHi_dt_down  / int_w_down
        deltaU_av_down(  vi) = int_w_deltaU_down  / int_w_down
        Ti_hom_av_down(  vi) = int_w_Ti_hom_down  / int_w_down
      END IF

      ! Calculate pore water fraction rates of change
      ! =============================================

      ! Compute likelihood of subglacial water
      ice%pore_water_likelihood( vi) = EXP(Ti_hom_av_up( vi)/3._dp)

      ! Check whether we want a local-, flowline-, or mixed-type inversion
      IF (    C%choice_pore_water_nudging_method == 'flowline') THEN

        ! Compute new adjustment for pore water fraction based on up and down flowlines
        I_tot( vi) = ice%pore_water_likelihood( vi) * ((deltaHi_av_up( vi)                      ) / porenudge_H_dHdt_flowline_dH0 + &
                                                       (dHi_dt_av_up(  vi) + dHi_dt_av_down( vi)) / porenudge_H_dHdt_flowline_dHdt0 + &
                                                       (deltaU_av_up(  vi)                      ) / porenudge_H_dHdt_flowline_dU0)

      ELSEIF (C%choice_pore_water_nudging_method == 'mixed') THEN

        ! Compute new adjustment for pore water fraction over too-thin regions based on upstream flowline
        I_tot( vi) = ice%pore_water_likelihood( vi) * ((deltaHi_av_up( vi)       ) / porenudge_H_dHdt_flowline_dH0 + &
                                                       (2._dp * dHi_dt_av_up( vi)) / porenudge_H_dHdt_flowline_dHdt0 + &
                                                       (deltaU_av_up(  vi)       ) / porenudge_H_dHdt_flowline_dU0)

      ELSEIF (C%choice_pore_water_nudging_method == 'local') THEN

        ! Nothing to do here, as this point will not be reached in this method

      ELSE
        CALL crash('unknown choice_pore_water_nudging_method "' // TRIM( C%choice_pore_water_nudging_method) // '"!')
      END IF

      ! Compute a rate of adjustment
      dC1_dt( vi) = -1._dp * (I_tot( vi) * pore_dryness( vi)) / porenudge_H_dHdt_flowline_t_scale

    END DO ! vi = mesh%vi1, mesh%vi2

    ! Regularise tricky areas
    ! =======================

    DO vi = mesh%vi1, mesh%vi2

      ! Grounded fractions
      fg_exp_mod = 1._dp - ice%fraction_gr( vi)
      ! Steep slopes
      hs_exp_mod = MIN( 1.0_dp, MAX( 0._dp, MAX( 0._dp, ice%Hs_slope( vi) - 0.003_dp) / (0.01_dp - 0.003_dp) ))

      hi_exp_mod = MIN( 1.0_dp, MAX( 0._dp, ice%Hi( vi) / 200._dp ))

      ! Prevent over-increase of sliding over partially grounded, steep-sloped areas
      IF (dC1_dt( vi) < 0._dp) THEN
        ! Scale based on friction-slope-slide-thickness modifiers
        dC1_dt( vi) = dC1_dt( vi) * (1._dp - .5_dp * (fg_exp_mod + hs_exp_mod)) * hi_exp_mod
      END IF

    END DO

    ! Smoothing
    ! =========

    IF (C%porenudge_H_dHdt_flowline_w_smooth > 0._dp) THEN

      dC1_dt_smoothed = dC1_dt

      ! Smooth the local variable
      CALL smooth_Gaussian( mesh, grid_smooth, dC1_dt_smoothed, C%porenudge_H_dHdt_flowline_r_smooth)

      DO vi = mesh%vi1, mesh%vi2
        dC1_dt( vi) = (1._dp - C%porenudge_H_dHdt_flowline_w_smooth) * dC1_dt( vi) + C%porenudge_H_dHdt_flowline_w_smooth * dC1_dt_smoothed( vi)
      END DO

    END IF

    ! Inverted pore dryness field
    ! ===========================

    pore_dryness = pore_dryness + dC1_dt * C%pore_water_nudging_dt

    ! New pore water fraction field
    ! =============================

    ! Get pore_water_fraction back from pore_dryness
    HIV%pore_water_fraction_next = 1._dp - pore_dryness

    ! Prevent fast change over unstable points
    ! ========================================

    ! This step takes the number of times each vertex has caused a
    ! reduction in the time step during the computation of the new
    ! ice thickness within the predictor-corrector method. The more
    ! times the pixel has triggered this reduction, the less its
    ! inverted pore water fraction field is allowed to evolve.
    ! However, this can prevent a recovery from an overshoot that
    ! caused too much basal sliding. Therefore, allow for increases
    ! of the basal friction over these problematic regions.

    IF (C%choice_timestepping == 'pc') THEN

      ! Smooth the instability field
      unstable_vertex_smoothed = 1._dp - EXP( REAL( MIN( 0, -ice%pc%tau_n_guilty + 0), dp) / 1._dp)
      CALL smooth_Gaussian( mesh, grid_smooth, unstable_vertex_smoothed, 40000._dp)

      ! Merge the smoothed and original instability fields: the guiltier
      ! the vertex, the stronger the original field dominates there.
      DO vi = mesh%vi1, mesh%vi2
        ! Original value
        unstable_vertex = 1._dp - EXP( REAL( MIN( 0, -ice%pc%tau_n_guilty( vi) + 0), dp) / 1._dp)
        ! Final smoothed field: a weighed average between the original and preliminary smoothed fields
        unstable_vertex_smoothed( vi) = (1._dp - unstable_vertex) * unstable_vertex_smoothed( vi) + unstable_vertex * unstable_vertex
      END DO

      ! Allow for some level of change even in the guiltiest of cases
      unstable_vertex_smoothed = MIN( 0.5_dp, MAX( 0._dp, unstable_vertex_smoothed))

      ! Only apply the reduction of adjustments in areas
      ! where sliding is expected to increase.
      DO vi = mesh%vi1, mesh%vi2
        IF (HIV%pore_water_fraction_next( vi) > HIV%pore_water_fraction_prev( vi)) THEN
          HIV%pore_water_fraction_next( vi) = (1._dp - unstable_vertex_smoothed( vi)) * HIV%pore_water_fraction_next( vi) + &
                                                       unstable_vertex_smoothed( vi)  * HIV%pore_water_fraction_prev( vi)
        END IF
      END DO

      ! Limit final pore_water_fraction_next values to hard limit
      HIV%pore_water_fraction_next = MIN( HIV%pore_water_fraction_next, 1._dp)
      HIV%pore_water_fraction_next = MAX( HIV%pore_water_fraction_next, 0._dp)

    END IF

    ! Extend onto neighbour vertices
    ! ==============================

    ! Gather ice model data from all processes
    CALL gather_to_all( HIV%pore_water_fraction_next, pore_water_fraction_next_tot)

    ! First, check at the grounded margins for the highest
    ! pore water fraction among inverted neighbours, to make
    ! sure that the ice flows smoothly

    ! Grounded margins
    DO vi = mesh%vi1, mesh%vi2

      ! Skip if not grounded margin
      IF (.NOT. (mask_cf_gr_tot( vi) .OR. mask_gl_gr_tot( vi) .OR. mask_margin_tot( vi))) CYCLE

      ! Initialise maximum neighbour
      found_grounded_neighbour = .FALSE.
      max_neighbour = 0._dp

      ! Check interior grounded vertices for greater values
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi, ci)
        IF (mask_grounded_ice_tot( vc) .AND. .NOT. mask_cf_gr_tot( vc) .AND. .NOT. mask_gl_gr_tot( vc) .AND. .NOT. mask_margin_tot( vc)) THEN
          max_neighbour = MAX( max_neighbour, pore_water_fraction_next_tot( vc))
          found_grounded_neighbour = .TRUE.
        END IF
      END DO

      IF (found_grounded_neighbour) THEN
        ! Use maximum value among neighbours, if any
        HIV%pore_water_fraction_next( vi) = max_neighbour
        HIV%mask_inverted_point( vi) = .TRUE.
      ELSE
        IF (ice%Hib( vi) < ice%SL( vi)) THEN
          ! Compute exponent for this vertex's weight based on ice thickness
          exponent_gr = MAX( .1_dp, LOG10( MAX( 1._dp, ice%Hi( vi))) - 2._dp)
          ! Get estimated pore water fraction based on grounded area fraction
          HIV%pore_water_fraction_next( vi) = 1._dp - ice%fraction_gr( vi)**exponent_gr
        ELSE
          ! Assume no hydrology
          HIV%pore_water_fraction_next( vi) = 0._dp
        END IF
      END IF

    END DO

    ! Then, check at the floating side of the grounding
    ! line for the highest pore water fraction among
    ! grounded grounding line neighbours

    ! Margins: floating side of grounding line
    DO vi = mesh%vi1, mesh%vi2

      ! Skip if not floating grounding line
      IF (.NOT. ice%mask_gl_fl( vi)) CYCLE

      ! Initialise maximum neighbour
      found_grounded_neighbour = .FALSE.
      max_neighbour = 0._dp

      ! Check grounded grounding lines for greater values
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi, ci)
        IF (mask_gl_gr_tot( vc)) THEN
          max_neighbour = MAX( max_neighbour, pore_water_fraction_next_tot( vc))
          found_grounded_neighbour = .TRUE.
        END IF
      END DO

      IF (found_grounded_neighbour) THEN
        ! Use maximum value among neighbours, if any
        HIV%pore_water_fraction_next( vi) = max_neighbour
        ! Consider this vertex a valid inversion point
        HIV%mask_inverted_point( vi) = .TRUE.
      ELSE
        ! Compute exponent for this vertex's weight based on ice thickness
        exponent_gr = MAX( .1_dp, LOG10( MAX( 1._dp, ice%Hi( vi))) - 2._dp)
        ! Get estimated pore water fraction based on grounded area fraction
        HIV%pore_water_fraction_next( vi) = 1._dp - ice%fraction_gr( vi)**exponent_gr
      END IF

    END DO

    ! Then, do the same for ice-free land next
    ! to grounded ice, to ensure a smooth transition

    ! Margins: land next to grounded front
    DO vi = mesh%vi1, mesh%vi2

      ! Skip if not ice-free land
      IF (.NOT. ice%mask_icefree_land( vi)) CYCLE

      ! Initialise maximum neighbour
      found_grounded_neighbour = .FALSE.
      max_neighbour = 0._dp

      ! Check grounded grounding lines for greater values
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi, ci)
        IF (mask_margin_tot( vc)) THEN
          max_neighbour = MAX( max_neighbour, pore_water_fraction_next_tot( vc))
          found_grounded_neighbour = .TRUE.
        END IF
      END DO

      IF (found_grounded_neighbour) THEN
        ! Use maximum value among neighbours, if any
        HIV%pore_water_fraction_next( vi) = max_neighbour
        ! Consider this vertex a valid inversion point
        HIV%mask_inverted_point( vi) = .TRUE.
      ELSE
        ! Assume no hydrology
        HIV%pore_water_fraction_next( vi) = 0._dp
      END IF

    END DO

    ! Finally, assign maximum values for interior ice
    ! shelf and ice-free ocean vertices. However, do
    ! not mark it as a valid inverted vertex, since
    ! this one should be overwritten later by the
    ! flowline extrapolation and ocean entrainment
    ! to account for the possibility of it not being
    ! floating in the future. If that extrapolation
    ! fails, then this value will be used though, so
    ! keep that in mind.

    ! Interior shelves and ocean
    DO vi = mesh%vi1, mesh%vi2

      IF (ice%mask_floating_ice( vi) .OR. ice%mask_icefree_ocean( vi)) THEN

        ! Skip floating side of grounding line
        IF (ice%mask_gl_fl( vi)) CYCLE

        ! Use maximum value
        HIV%pore_water_fraction_next( vi) = 1._dp

      END IF

    END DO

    ! Limit inverted values
    ! =====================

    ! Limit values to prescribed limits
    DO vi = mesh%vi1, mesh%vi2
      HIV%pore_water_fraction_next( vi) = MIN( HIV%pore_water_fraction_next( vi), C%pore_water_fraction_max)
      HIV%pore_water_fraction_next( vi) = MAX( HIV%pore_water_fraction_next( vi), C%pore_water_fraction_min)
    END DO

    ! ! DENK DROM : save smoothed field in a host variable for later inspection
    ! ice%pore_water_likelihood = unstable_vertex_smoothed

    ! ! DENK DROM
    ! ice%pore_water_likelihood = 0._dp
    ! DO vi = mesh%vi1, mesh%vi2
    !   IF (HIV%mask_inverted_point( vi)) ice%pore_water_likelihood( vi) = 1._dp
    ! END DO

    ! Finalise
    ! ========

    ! Clean up after yourself
    DEALLOCATE( pore_dryness   )
    DEALLOCATE( trace_up       )
    DEALLOCATE( trace_down     )
    DEALLOCATE( s_up           )
    DEALLOCATE( s_down         )
    DEALLOCATE( deltaHi_up     )
    DEALLOCATE( deltaHi_down   )
    DEALLOCATE( dHi_dt_up      )
    DEALLOCATE( dHi_dt_down    )
    DEALLOCATE( deltaU_up      )
    DEALLOCATE( deltaU_down    )
    DEALLOCATE( Ti_hom_up      )
    DEALLOCATE( Ti_hom_down    )
    DEALLOCATE( deltaHi_av_up  )
    DEALLOCATE( deltaHi_av_down)
    DEALLOCATE( dHi_dt_av_up   )
    DEALLOCATE( dHi_dt_av_down )
    DEALLOCATE( deltaU_av_up   )
    DEALLOCATE( deltaU_av_down )
    DEALLOCATE( Ti_hom_av_up   )
    DEALLOCATE( Ti_hom_av_down )
    DEALLOCATE( I_tot          )
    DEALLOCATE( dC1_dt         )
    DEALLOCATE( dC1_dt_smoothed)
    DEALLOCATE( unstable_vertex_smoothed)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE pore_water_fraction_inversion

! ===== Utilities =====
! =====================

  SUBROUTINE apply_extrapolation_to_pore_water_fraction( mesh, ice, HIV)
    ! Extrapolate inverted pore water fraction using a flowline approach

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_hydrology_inversion),      INTENT(IN)    :: HIV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_extrapolation_to_pore_water_fraction'
    INTEGER                                            :: vi
    REAL(dp)                                           :: weight_gr, exponent_gr, ocean_entrainment
    REAL(dp), DIMENSION(mesh%nV)                       :: pore_water_fraction_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_inverted_point_tot
    REAL(dp), DIMENSION(mesh%nTri)                     :: u_b_tot
    REAL(dp), DIMENSION(mesh%nTri)                     :: v_b_tot
    REAL(dp), DIMENSION(2)                             :: p
    LOGICAL                                            :: found_source

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Retrieve inverted field
    ice%pore_water_fraction = HIV%pore_water_fraction_app

    ! Gather inversion data, inversion mask, and horizontal velocities from all processes
    CALL gather_to_all(      ice%pore_water_fraction, pore_water_fraction_tot)
    CALL gather_to_all( HIV%mask_inverted_point, mask_inverted_point_tot)
    CALL gather_to_all(      ice%u_vav_b            , u_b_tot                )
    CALL gather_to_all(      ice%v_vav_b            , v_b_tot                )

    ! Extrapolate non-inverted areas using closest upstream inverted value
    DO vi = mesh%vi1, mesh%vi2

      ! Skip original inverted points
      IF (HIV%mask_inverted_point( vi)) CYCLE

      ! Skip ice-free ocean points
      IF (ice%mask_icefree_ocean( vi)) CYCLE

      ! Skip ice shelf points away from the grounding line
      IF (ice%mask_floating_ice( vi) .AND. .NOT. ice%mask_gl_fl( vi)) CYCLE

      ! Get x and y coordinates of this non-inverted vertex
      p = [mesh%V( vi,1), mesh%V( vi,2)]

      ! Extrpolate value from the closest upstream inverted vertex
      CALL extrapolate_from_upstream( mesh, ice, mask_inverted_point_tot, pore_water_fraction_tot, u_b_tot, v_b_tot, p, found_source)

      ! Check if an inverted value was found upstream
      IF (found_source) THEN
        ! If so, use that value for this vertex
        ice%pore_water_fraction( vi) = pore_water_fraction_tot( vi)
      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_extrapolation_to_pore_water_fraction

  SUBROUTINE apply_grounded_fractions_to_pore_water_fraction( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Use the inverted values of pore water fraction

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_grounded_fractions_to_pore_water_fraction'
    INTEGER                                            :: vi
    REAL(dp)                                           :: weight_gr, exponent_hi, exponent_hs, exponent_gr, ocean_entrainment

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this is actually wanted
    IF (.NOT. C%do_subgrid_friction_on_A_grid) THEN
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! And exit
      RETURN
    END IF

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise grounded area fraction weight
      weight_gr = 1._dp

      ! Compute exponent for this vertex's weight based on ice thickness
      exponent_hi = LOG10( MAX( 1._dp, ice%Hi( vi)))
      ! Compute exponent for this vertex's weight based on ice gradients
      exponent_hs = ice%Hs_slope( vi) / 0.005_dp
      ! Compute final exponent for this vertex's weight
      exponent_gr = MAX( 0._dp, exponent_hi - exponent_hs)

      ! Compute a weight based on the grounded area fractions
      IF (ice%mask_gl_gr( vi)) THEN
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      ELSEIF (ice%mask_cf_gr( vi)) THEN
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      ELSEIF (ice%mask_gl_fl( vi)) THEN
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      ELSEIF (ice%mask_grounded_ice( vi)) THEN
        weight_gr = 1._dp

      ELSEIF (ice%mask_floating_ice( vi)) THEN
        weight_gr = 0._dp

      ELSEIF (ice%mask_icefree_ocean( vi)) THEN
        weight_gr = 0._dp

      END IF

      ! Just in case
      weight_gr = MIN( 1._dp, MAX( 0._dp, weight_gr))

      ! DENK DROM : Test whether bathymetry is actually needed / helps with stability
      ! Compute entrainment factor based on bathymetry
      ocean_entrainment = 1._dp! - (ice%Hb( vi) - ice%SL( vi) + 500._dp) / 500._dp

      ! Limit ocean entrainment factor to a [0-1] range
      ocean_entrainment = MIN( 1._dp, MAX( 0._dp, ocean_entrainment))

      ! Add ocean entrainment to inverted pore water fraction through a weighed average
      ice%pore_water_fraction( vi) = weight_gr * ice%pore_water_fraction( vi) + (1._dp - weight_gr) * ocean_entrainment

      ! Limit final pore water fraction field to a valid range
      IF (ice%mask_grounded_ice( vi) .OR. ice%mask_icefree_land( vi)) THEN
        ! Grounded vertices
        ice%pore_water_fraction(vi) = MIN( ice%pore_water_fraction(vi), C%pore_water_fraction_max)
        ice%pore_water_fraction(vi) = MAX( ice%pore_water_fraction(vi), C%pore_water_fraction_min)
      ELSE
        ! Floating vertices
        ice%pore_water_fraction( vi) = MIN( 1._dp, MAX( 0._dp, ice%pore_water_fraction( vi)))
      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_grounded_fractions_to_pore_water_fraction

  SUBROUTINE extrapolate_from_upstream( mesh, ice, d_mask_tot, d_values_tot, u_b_tot, v_b_tot, p, found_source)
    ! Extrapolate an inverted field based on its values upstream of
    ! the point of interest based on its flowline, computed based on
    ! the 2-D velocity field u_b,v_b.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh                ! Mesh data
    TYPE(type_ice_model),                INTENT(IN)    :: ice                 ! Ice model structure
    LOGICAL,  DIMENSION(mesh%nV),        INTENT(IN)    :: d_mask_tot          ! Inversion mask: 1=inverted; 0=not-inverted
    REAL(dp), DIMENSION(mesh%nV),        INTENT(INOUT) :: d_values_tot        ! Inversion data: the field you want to extrapolate
    REAL(dp), DIMENSION(mesh%nTri),      INTENT(IN)    :: u_b_tot, v_b_tot    ! The u and v velocity fields
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p                   ! Your point of interest, i.e. where you want to extrapolate
    LOGICAL,                             INTENT(OUT)   :: found_source        ! Flag: 1=successful extrapolation; 0=failed extrapolation

    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: pt                  ! Current location of tracer
    INTEGER                                            :: vi, vp, iti, ti, n  ! Some indices
    REAL(dp)                                           :: dist, w, w_tot      ! For the interpolation of velocities at point pt
    REAL(dp)                                           :: u_pt, v_pt, uabs_pt ! Interpolated velocities at point pt
    REAL(dp), DIMENSION(2)                             :: u_hat_pt            ! Direction of the velocity vector at point pt
    REAL(dp)                                           :: dist_prev, dist_tot ! Safety: total travelled distance of tracer

    ! If the vertex vp containing this point is part of the
    ! original inverted field, no need to extrapolate there
    vp = 1
    CALL find_containing_vertex( mesh, p, vp)
    IF (d_mask_tot( vp)) THEN
      RETURN
    END IF

    ! Initialise
    n  = 0  ! Number of tracer movements
    pt = p  ! Current point of the tracer
    vi = vp ! Current vertex containing the tracer

    ! Initialise flag for successful extrapolation
    found_source = .FALSE.

    ! Safety
    dist_prev = 0._dp ! Current distance from origin
    dist_tot  = 0._dp ! Total distance travelled by tracer

    DO WHILE (.TRUE.)

      ! Check if tracer reached its destination
      ! =======================================

      ! Find the vertex vi currently containing the new tracer position
      CALL find_containing_vertex( mesh, pt, vi)

      ! If the current vertex is an inverted point,
      ! use it as the value for our origin vertex,
      ! mark it as a successfully extrapolated point,
      ! and exit
      IF (d_mask_tot( vi)) THEN
        d_values_tot( vp) = d_values_tot( vi)
        found_source = .TRUE.
        RETURN
      END IF

      ! If not, move the tracer upstream
      ! ================================

      ! Interpolate between the surrounding triangles
      ! to find the velocities at the tracer's location
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
      ! we find velocities below 0.01 m/yr), end the trace. This also
      ! prevents the tracer from getting stuck in place when the
      ! local velocity field is zero.
      IF (uabs_pt < .01_dp) EXIT

      ! Calculate the normalised velocity vector at the tracer's location
      u_hat_pt = [u_pt / uabs_pt, v_pt / uabs_pt]

      ! Add to counter of tracer movements
      n = n + 1
      ! Safety
      IF (n > SIZE( mesh%V,1)) THEN
        CALL crash('upstream flowline tracer got stuck!')
      END IF

      ! Save previous distance-to-origin
      dist_prev = NORM2( pt - p)

      ! Move the tracer upstream by a distance of 1/2 local resolution
      pt = pt - u_hat_pt * .5_dp * mesh%R( vi)

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

  END SUBROUTINE extrapolate_from_upstream

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

  SUBROUTINE initialise_pore_water_fraction_from_file( mesh, ice, region_name)
    ! Initialise the pore water fraction from file
    !
    ! pore water only, no time

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_pore_water_fraction_from_file'
    CHARACTER(LEN=256)                                    :: filename_pore_water_fraction
    REAL(dp)                                              :: timeframe_pore_water_fraction

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE ('NAM')
        filename_pore_water_fraction  = C%filename_pore_water_fraction_NAM
        timeframe_pore_water_fraction = C%timeframe_pore_water_fraction_NAM
      CASE ('EAS')
        filename_pore_water_fraction  = C%filename_pore_water_fraction_EAS
        timeframe_pore_water_fraction = C%timeframe_pore_water_fraction_EAS
      CASE ('GRL')
        filename_pore_water_fraction  = C%filename_pore_water_fraction_GRL
        timeframe_pore_water_fraction = C%timeframe_pore_water_fraction_GRL
      CASE ('ANT')
        filename_pore_water_fraction  = C%filename_pore_water_fraction_ANT
        timeframe_pore_water_fraction = C%timeframe_pore_water_fraction_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    if (index( filename_pore_water_fraction,'_LAST.nc') > 1) then
      call find_last_output_file( filename_pore_water_fraction)
      call find_last_timeframe(   filename_pore_water_fraction, timeframe_pore_water_fraction)
    end if

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising pore water fraction from file "' // colour_string( TRIM(filename_pore_water_fraction),'light blue') // '"...'

    ! Read pore water fraction from file
    IF (timeframe_pore_water_fraction == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_2D( filename_pore_water_fraction, 'pore_water_fraction', mesh, ice%pore_water_fraction)
    ELSE
      ! Assume the file has a time dimension, and read the specified timeframe
      CALL read_field_from_file_2D( filename_pore_water_fraction, 'pore_water_fraction', mesh, ice%pore_water_fraction, time_to_read = timeframe_pore_water_fraction)
    END IF

    CALL calc_pore_water_pressure_from_file( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_pore_water_fraction_from_file

  ! == Pore water fraction from file
  SUBROUTINE calc_pore_water_pressure_from_file( mesh, ice)
    ! Calculate the pore water pressure from a prescribed pore water fraction, read from file
    !
    ! Use the parameterisation from Martin et al. (2011)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_from_file'
    INTEGER                                            :: vi
    REAL(dp)                                           :: weight_gr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute effective pore water pressure
    ! =====================================

    ! Scale pore water fraction based on grounded area fractions
    CALL apply_grounded_fractions_to_pore_water_fraction( mesh, ice)

    ! Compute pore water pressure based on the pore water fraction as
    ! the fraction of the overburden pressure supported by basal water
    DO vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure( vi) = ice%pore_water_fraction(vi) * ice_density * grav * ice%Hi_eff( vi)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_from_file

END MODULE basal_hydrology
