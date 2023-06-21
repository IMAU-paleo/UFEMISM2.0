MODULE ice_model_main

  ! The main ice-dynamical model module.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE scalar_types                                           , ONLY: type_regional_scalars
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ice_model_memory                                       , ONLY: allocate_ice_model
  USE ice_model_utilities                                    , ONLY: determine_masks, calc_bedrock_CDFs, calc_grounded_fractions
  USE ice_model_scalars                                      , ONLY: calc_ice_model_scalars
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE math_utilities                                         , ONLY: ice_surface_elevation, thickness_above_floatation
  USE basal_conditions_main                                  , ONLY: initialise_basal_conditions
  USE ice_velocity_main                                      , ONLY: initialise_velocity_solver, solve_stress_balance, remap_velocity_solver, &
                                                                     create_restart_file_ice_velocity, write_to_restart_file_ice_velocity
  USE region_types                                           , ONLY: type_model_region

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ice_dynamics_model( region)
    ! Advance the ice dynamics model by one time step,
    ! which is determined by the ice dynamics itself.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ice_dynamics_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_model

! ===== Time stepping =====
! =========================

  SUBROUTINE calc_critical_timestep_SIA( mesh, ice, dt_crit_SIA)
    ! Calculate the critical time step for advective ice flow (CFL criterion)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_SIA'
    INTEGER                                            :: ti, via, vib, vic
    REAL(dp)                                           :: d_ab, d_bc, d_ca, d_min, Hi, D_SIA, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise time step with maximum allowed value
    dt_crit_SIA = C%dt_ice_max

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate shortest triangle side
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      d_ab = NORM2( mesh%V( vib,:) - mesh%V( via,:))
      d_bc = NORM2( mesh%V( vic,:) - mesh%V( vib,:))
      d_ca = NORM2( mesh%V( via,:) - mesh%V( vic,:))

      d_min = MINVAL([ d_ab, d_bc, d_ca])

      ! Find maximum diffusivity in the vertical column
      D_SIA = MAXVAL( ABS( ice%SIA%D_3D_b( ti,:)))

      ! Calculate critical timestep
      Hi = MAXVAL( [ice%Hi( via), ice%Hi( vib), ice%Hi( vic)])
      dt = d_min**2 / (6._dp * Hi * D_SIA)
      dt_crit_SIA = MIN( dt_crit_SIA, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_SIA = dt_crit_SIA * dt_correction_factor

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_SIA
!
!  SUBROUTINE calc_critical_timestep_adv( mesh, ice, dt_crit_adv)
!    ! Calculate the critical time step for advective ice flow (CFL criterion)
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(IN)    :: ice
!    REAL(dp),                            INTENT(OUT)   :: dt_crit_adv
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_adv'
!    INTEGER                                            :: aci, vi, vj
!    REAL(dp)                                           :: dist, dt
!    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    dt_crit_adv = 2._dp * C%dt_max
!
!    DO aci = mesh%ci1, mesh%ci2
!
!      ! Only check at ice-covered vertices
!      vi = mesh%Aci( aci,1)
!      vj = mesh%Aci( aci,2)
!      IF (ice%Hi_a( vi) == 0._dp .OR. ice%Hi_a( vj) == 0._dp) CYCLE
!
!      dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
!      dt = dist / (ABS( ice%u_vav_c( aci)) + ABS( ice%v_vav_c( aci)))
!      dt_crit_adv = MIN( dt_crit_adv, dt)
!
!    END DO
!
!    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
!    dt_crit_adv = MIN(C%dt_max, dt_crit_adv * dt_correction_factor)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_critical_timestep_adv
!
!  SUBROUTINE calc_pc_truncation_error( mesh, ice, dt, dt_prev)
!    ! Calculate the truncation error in the ice thickness rate of change (Robinson et al., 2020, Eq. 32)
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp),                            INTENT(IN)    :: dt, dt_prev
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pc_truncation_error'
!    INTEGER                                            :: vi, ci, vc
!    LOGICAL                                            :: has_GL_neighbour
!    REAL(dp)                                           :: zeta, eta_proc
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Ratio of time steps
!    zeta = dt / dt_prev
!
!    ! Find maximum truncation error
!    eta_proc = C%pc_eta_min
!
!    DO vi = mesh%vi1, mesh%vi2
!
!      ! Calculate truncation error (Robinson et al., 2020, Eq. 32)
!      ice%pc_tau( vi) = ABS( zeta * (ice%Hi_corr( vi) - ice%Hi_pred( vi)) / ((3._dp * zeta + 3._dp) * dt))
!
!      IF (ice%mask_sheet_a( vi) == 1) THEN
!
!        has_GL_neighbour = .FALSE.
!        DO ci = 1, mesh%nC( vi)
!          vc = mesh%C( vi,ci)
!          IF (ice%mask_gl_a( vc) == 1) THEN
!            has_GL_neighbour = .TRUE.
!            EXIT
!          END IF
!        END DO
!
!        IF (.NOT.has_GL_neighbour) eta_proc = MAX( eta_proc, ice%pc_tau( vi))
!
!      END IF
!
!    END DO
!    CALL MPI_REDUCE( eta_proc, ice%pc_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_pc_truncation_error
!
!  SUBROUTINE determine_timesteps_and_actions( region, t_end)
!    ! Determine how long we can run just ice dynamics before another "action" (thermodynamics,
!    ! GIA, output writing, inverse routine, etc.) has to be performed, and adjust the time step accordingly.
!
!    IMPLICIT NONE
!
!    ! Input variables:
!    TYPE(type_model_region),             INTENT(INOUT) :: region
!    REAL(dp),                            INTENT(IN)    :: t_end
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_timesteps_and_actions'
!    REAL(dp)                                           :: t_next
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (par%master) THEN
!
!      ! Determine when each model components should be updated
!
!      t_next = MIN(t_end, region%time + C%dt_max)
!
!      ! First the ice dynamics
!      ! ======================
!
!      IF     (C%choice_stress_balance_approximation == 'none') THEN
!        ! Just stick to the maximum time step
!      ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
!        t_next = MIN( t_next, region%t_next_SIA)
!      ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
!        t_next = MIN( t_next, region%t_next_SSA)
!      ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
!        t_next = MIN( t_next, region%t_next_SIA)
!        t_next = MIN( t_next, region%t_next_SSA)
!      ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
!        t_next = MIN( t_next, region%t_next_DIVA)
!      ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
!        t_next = MIN( t_next, region%t_next_BPA)
!      ELSE
!        CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
!      END IF ! IF (C%choice_stress_balance_approximation == 'SIA') THEN
!
!      ! Then the other model components
!      ! ===============================
!
!      region%do_thermo  = .FALSE.
!      IF (region%time == region%t_next_thermo) THEN
!        region%do_thermo      = .TRUE.
!        region%t_last_thermo  = region%time
!        region%t_next_thermo  = region%t_last_thermo + C%dt_thermo
!      END IF
!      t_next = MIN( t_next, region%t_next_thermo)
!
!      region%do_climate = .FALSE.
!      IF (region%time == region%t_next_climate) THEN
!        region%do_climate     = .TRUE.
!        region%t_last_climate = region%time
!        region%t_next_climate = region%t_last_climate + C%dt_climate
!      END IF
!      t_next = MIN( t_next, region%t_next_climate)
!
!      region%do_ocean   = .FALSE.
!      IF (region%time == region%t_next_ocean) THEN
!        region%do_ocean       = .TRUE.
!        region%t_last_ocean   = region%time
!        region%t_next_ocean   = region%t_last_ocean + C%dt_ocean
!      END IF
!      t_next = MIN( t_next, region%t_next_ocean)
!
!      region%do_SMB     = .FALSE.
!      IF (region%time == region%t_next_SMB) THEN
!        region%do_SMB         = .TRUE.
!        region%t_last_SMB     = region%time
!        region%t_next_SMB     = region%t_last_SMB + C%dt_SMB
!      END IF
!      t_next = MIN( t_next, region%t_next_SMB)
!
!      region%do_BMB     = .FALSE.
!      IF (region%time == region%t_next_BMB) THEN
!        region%do_BMB         = .TRUE.
!        region%t_last_BMB     = region%time
!        region%t_next_BMB     = region%t_last_BMB + C%dt_BMB
!      END IF
!      t_next = MIN( t_next, region%t_next_BMB)
!
!      region%do_ELRA    = .FALSE.
!      IF (C%choice_GIA_model == 'ELRA') THEN
!        IF (region%time == region%t_next_ELRA) THEN
!          region%do_ELRA      = .TRUE.
!          region%t_last_ELRA  = region%time
!          region%t_next_ELRA  = region%t_last_ELRA + C%dt_bedrock_ELRA
!        END IF
!        t_next = MIN( t_next, region%t_next_ELRA)
!      END IF
!
!      region%do_basal    = .FALSE.
!      IF (region%time == region%t_next_basal) THEN
!        region%do_basal       = .TRUE.
!        region%t_last_basal   = region%time
!        region%t_next_basal   = region%t_last_basal + C%BIVgeo_dt
!      END IF
!      t_next = MIN( t_next, region%t_next_basal)
!
!      region%do_SMB_inv = .FALSE.
!      IF (region%time == region%t_next_SMB_inv) THEN
!        region%do_SMB_inv       = .TRUE.
!        region%t_last_SMB_inv   = region%time
!        region%t_next_SMB_inv   = region%t_last_SMB_inv + C%dt_SMB_inv
!      END IF
!      t_next = MIN( t_next, region%t_next_SMB_inv)
!
!      region%do_output  = .FALSE.
!      IF (region%time == region%t_next_output) THEN
!        region%do_output      = .TRUE.
!        region%t_last_output  = region%time
!        region%t_next_output  = region%t_last_output + C%dt_output
!      END IF
!      t_next = MIN( t_next, region%t_next_output)
!
!      ! Set time step so that we move forward to the next action
!      region%dt = t_next - region%time
!
!    END IF ! IF (par%master) THEN
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE determine_timesteps_and_actions

! ===== Administration: allocation, initialisation, and remapping =====
! =====================================================================

  SUBROUTINE initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ! Initialise all data fields of the ice module

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    TYPE(type_ice_model),          INTENT(INOUT) :: ice
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo_init
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo_PD
    TYPE(type_regional_scalars),   INTENT(OUT)   :: scalars
    CHARACTER(LEN=3),              INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_ice_model'
    INTEGER                                      :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      WRITE(*,"(A)") '  Initialising ice dynamics model...'
    END IF
    CALL sync

    ! === Memory allocation ===
    ! =========================

    ! Allocate memory
    CALL allocate_ice_model( mesh, ice)

    ! === Value initialisation ===
    ! ============================

    ! Sea level
    ! =========

    SELECT CASE (C%choice_sealevel_model)

      CASE ('fixed')
        ! Fixed sea level
        ice%SL = C%fixed_sealevel

      CASE ('prescribed')
        ! Sea-level prescribed from external record file
        CALL crash('Sea level initialisation: prescribed method not implement yet!')
        ! ice%SL = forcing%sealevel_obs

      CASE ('eustatic')
        ! Eustatic sea level
        CALL crash('Sea level initialisation: eustatic method not implement yet!')
        ! ice%SL = C%initial_guess_sealevel

      CASE ('SELEN')
        ! Sea level from SELEN
        CALL crash('Sea level initialisation: SELEN method not implement yet!')
        ! ice%SL = C%initial_guess_sealevel

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_sealevel_model "' // &
                    TRIM( C%choice_sealevel_model) // '"!')

    END SELECT

    ! Initial topography
    ! ==================

    ! Basic topography
    DO vi = mesh%vi1, mesh%vi2

      ! Main quantities
      ice%Hi ( vi) = refgeo_init%Hi( vi)
      ice%Hb ( vi) = refgeo_init%Hb( vi)
      ice%Hs ( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))
      ice%Hib( vi) = ice%Hs( vi) - ice%Hi( vi)

      ice%TAF( vi)  = thickness_above_floatation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))

      ! Differences w.r.t. present-day
      ice%dHi ( vi)  = ice%Hi ( vi) - refgeo_PD%Hi ( vi)
      ice%dHb ( vi)  = ice%Hb ( vi) - refgeo_PD%Hb ( vi)
      ice%dHs ( vi)  = ice%Hs ( vi) - refgeo_PD%Hs ( vi)
      ice%dHib( vi)  = ice%Hib( vi) - (refgeo_PD%Hs ( vi) - refgeo_PD%Hi( vi))

    END DO

    ! Initialised predicted ice thickness
    ice%Hi_tplusdt = ice%Hi

    ! Initial masks
    ! =============

    ! Sector masks
    CALL determine_masks( mesh, ice)

    ! Initialise previous-time-step mask
    ice%mask_ice_prev = ice%mask_ice

    ! Initial rates of change
    ! =======================

    ice%dHi_dt  = 0._dp
    ice%dHb_dt  = 0._dp
    ice%dHs_dt  = 0._dp
    ice%dHib_dt = 0._dp

    ! Sub-grid fractions
    ! ==================

    ! Compute bedrock cumulative density function
    CALL calc_bedrock_CDFs( mesh, refgeo_PD, ice)
    ! Initialise sub-grid grounded-area fractions
    CALL calc_grounded_fractions( mesh, ice)

    ! Basal conditions
    ! ================

    ! Allocate and initialise basal conditions
    CALL initialise_basal_conditions( mesh, ice)

    ! Velocities
    ! ==========

    ! Initialise data and matrices for the velocity solver(s)
    CALL initialise_velocity_solver( mesh, ice, region_name)

    ! Ice-sheet-wide scalars
    ! ======================

    CALL calc_ice_model_scalars( mesh, ice, refgeo_PD, scalars)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_model

END MODULE ice_model_main
