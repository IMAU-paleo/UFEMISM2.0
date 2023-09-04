MODULE basal_inversion_main

  ! Contains all the routines for managing the basal inversion model

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
  USE basal_inversion_types                                  , ONLY: type_basal_inversion
  USE region_types                                           , ONLY: type_model_region
  USE basal_inversion_H_dHdt_flowline                        , ONLY: initialise_basal_inversion_H_dHdt_flowline, run_basal_inversion_H_dHdt_flowline
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
  USE mesh_remapping                                         , ONLY: smooth_Gaussian_2D
  USE mesh_operators                                         , ONLY: ddx_a_a_2D, ddy_a_a_2D

  IMPLICIT NONE

CONTAINS

  ! ===== Main routines =====
  ! =========================

  SUBROUTINE run_basal_inversion( region)
    ! Run the main basal inversion model

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_basal_inversion'
    INTEGER                                            :: vi
    REAL(dp)                                           :: wt_prev, wt_next

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Only do basal inversion within the specified time window
    IF (region%time < C%bed_roughness_nudging_t_start) THEN
      region%BIV%t_next = C%bed_roughness_nudging_t_start
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    IF (region%time == C%bed_roughness_nudging_t_end) THEN
      region%BIV%t_next = C%end_time_of_run
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! If the desired time is beyond the time of the next modelled bed roughness,
    ! run the basal inversion model to calculate a new next modelled bed roughness.
    ! =============================================================================

    IF (region%time == region%BIV%t_next) THEN
      ! Need to calculate new predicted bed roughness

      ! Store previous modelled bed roughness
      region%BIV%generic_bed_roughness_1_prev = region%BIV%generic_bed_roughness_1_next
      region%BIV%generic_bed_roughness_2_prev = region%BIV%generic_bed_roughness_2_next
      region%BIV%t_prev = region%BIV%t_next
      region%BIV%t_next = region%BIV%t_prev + C%bed_roughness_nudging_dt

      ! Run the basal inversion model to calculate a new next modelled bed roughness
      SELECT CASE (C%choice_bed_roughness_nudging_method)
        CASE ('H_dHdt_flowline')

          ! Run with the specified target geometry
          SELECT CASE (C%choice_inversion_target_geometry)
            CASE ('init')
              CALL run_basal_inversion_H_dHdt_flowline( region%mesh, region%grid_smooth, region%ice, region%refgeo_init, region%BIV)
            CASE ('PD')
              CALL run_basal_inversion_H_dHdt_flowline( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%BIV)
            CASE DEFAULT
              CALL crash('unknown choice_inversion_target_geometry "' // TRIM( C%choice_inversion_target_geometry) // '"!')
          END SELECT

        CASE ('Pien2023')

          ! Run with the specified target geometry
          SELECT CASE (C%choice_inversion_target_geometry)
            CASE ('init')
              CALL run_basal_inversion_Pien2023( region%mesh, region%grid_smooth, region%ice, region%refgeo_init, region%BIV)
            CASE ('PD')
              CALL run_basal_inversion_Pien2023( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%BIV)
            CASE DEFAULT
              CALL crash('unknown choice_inversion_target_geometry "' // TRIM( C%choice_inversion_target_geometry) // '"!')
          END SELECT

        CASE DEFAULT
          CALL crash('unknown choice_bed_roughness_nudging_method "' // TRIM( C%choice_bed_roughness_nudging_method) // '"!')
      END SELECT

    ELSEIF (region%time > region%BIV%t_next) THEN
      ! This should not be possible
      CALL crash('overshot the basal inversion time step')
    ELSE
      ! We're within the current BIV prediction window
    END IF ! IF (region%time == region%BIV%t_next) THEN

    ! Interpolate between previous and next modelled bed roughness
    ! to find the bed roughness at the desired time
    ! =================================================================

    ! Calculate time interpolation weights
    wt_prev = (region%BIV%t_next - region%time) / (region%BIV%t_next - region%BIV%t_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled bed roughness to desired time
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%BIV%generic_bed_roughness_1( vi) = wt_prev * region%BIV%generic_bed_roughness_1_prev( vi) + wt_next * region%BIV%generic_bed_roughness_1_next( vi)
      region%BIV%generic_bed_roughness_2( vi) = wt_prev * region%BIV%generic_bed_roughness_2_prev( vi) + wt_next * region%BIV%generic_bed_roughness_2_next( vi)
    END DO

    ! Update sliding law-specific bed roughness
    ! =========================================

    SELECT CASE (C%choice_sliding_law)
      CASE ('no_sliding')
        CALL crash('cannot run basal inversion for choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
      CASE ('idealised')
        CALL crash('cannot run basal inversion for choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
      CASE ('Weertman')
        ! Weertman sliding law; bed roughness is described by slid_beta_sq
        region%ice%slid_beta_sq = region%BIV%generic_bed_roughness_1
      CASE ('Coulomb')
        ! Coulomb sliding law; bed roughness is described by till_friction_angle
        region%ice%till_friction_angle = region%BIV%generic_bed_roughness_1
      CASE ('Budd')
        ! Budd-type sliding law; bed roughness is described by till_friction_angle
        region%ice%till_friction_angle = region%BIV%generic_bed_roughness_1
      CASE ('Tsai2015')
        ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
        region%ice%slid_alpha_sq = region%BIV%generic_bed_roughness_1
        region%ice%slid_beta_sq  = region%BIV%generic_bed_roughness_2
      CASE ('Schoof2005')
        ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
        region%ice%slid_alpha_sq = region%BIV%generic_bed_roughness_1
        region%ice%slid_beta_sq  = region%BIV%generic_bed_roughness_2
      CASE ('Zoet-Iverson')
        ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
        region%ice%till_friction_angle = region%BIV%generic_bed_roughness_1
      CASE DEFAULT
        CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_basal_inversion

  SUBROUTINE initialise_basal_inversion( mesh, ice, BIV, region_name)
    ! Initialise the main basal inversion model

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_basal_inversion),          INTENT(OUT)   :: BIV
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master) WRITE(0,*) ' Initialising basal inversion model "' // colour_string( TRIM( C%choice_bed_roughness_nudging_method),'light blue') // '"...'

    ! Allocate memory for main variables
    ! ==================================

    ALLOCATE( BIV%generic_bed_roughness_1( mesh%vi1:mesh%vi2))
    ALLOCATE( BIV%generic_bed_roughness_2( mesh%vi1:mesh%vi2))

    BIV%generic_bed_roughness_1 = 0._dp
    BIV%generic_bed_roughness_2 = 0._dp

    ALLOCATE( BIV%generic_bed_roughness_1_prev( mesh%vi1:mesh%vi2))
    ALLOCATE( BIV%generic_bed_roughness_2_prev( mesh%vi1:mesh%vi2))
    ALLOCATE( BIV%generic_bed_roughness_1_next( mesh%vi1:mesh%vi2))
    ALLOCATE( BIV%generic_bed_roughness_2_next( mesh%vi1:mesh%vi2))

    BIV%generic_bed_roughness_1_prev = 0._dp
    BIV%generic_bed_roughness_2_prev = 0._dp
    BIV%generic_bed_roughness_1_next = 0._dp
    BIV%generic_bed_roughness_2_next = 0._dp

    ! Timeframes
    BIV%t_prev   = C%start_time_of_run
    BIV%t_next   = C%start_time_of_run

    ! Get sliding law-specific bed roughness
    ! ======================================

    SELECT CASE (C%choice_sliding_law)
      CASE ('no_sliding')
        CALL crash('cannot run basal inversion for choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
      CASE ('idealised')
        CALL crash('cannot run basal inversion for choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
      CASE ('Weertman')
        ! Weertman sliding law; bed roughness is described by slid_beta_sq
        BIV%generic_bed_roughness_1      = ice%slid_beta_sq
        BIV%generic_bed_roughness_1_prev = ice%slid_beta_sq
        BIV%generic_bed_roughness_1_next = ice%slid_beta_sq
        BIV%generic_bed_roughness_2      = 0._dp
        BIV%generic_bed_roughness_2_prev = 0._dp
        BIV%generic_bed_roughness_2_next = 0._dp
      CASE ('Coulomb')
        ! Coulomb sliding law; bed roughness is described by till_friction_angle
        BIV%generic_bed_roughness_1      = ice%till_friction_angle
        BIV%generic_bed_roughness_1_prev = ice%till_friction_angle
        BIV%generic_bed_roughness_1_next = ice%till_friction_angle
        BIV%generic_bed_roughness_2      = 0._dp
        BIV%generic_bed_roughness_2_prev = 0._dp
        BIV%generic_bed_roughness_2_next = 0._dp
      CASE ('Budd')
        ! Budd-type sliding law; bed roughness is described by till_friction_angle
        BIV%generic_bed_roughness_1      = ice%till_friction_angle
        BIV%generic_bed_roughness_1_prev = ice%till_friction_angle
        BIV%generic_bed_roughness_1_next = ice%till_friction_angle
        BIV%generic_bed_roughness_2      = 0._dp
        BIV%generic_bed_roughness_2_prev = 0._dp
        BIV%generic_bed_roughness_2_next = 0._dp
      CASE ('Tsai2015')
        ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
        BIV%generic_bed_roughness_1      = ice%slid_alpha_sq
        BIV%generic_bed_roughness_1_prev = ice%slid_alpha_sq
        BIV%generic_bed_roughness_1_next = ice%slid_alpha_sq
        BIV%generic_bed_roughness_2      = ice%slid_beta_sq
        BIV%generic_bed_roughness_2_prev = ice%slid_beta_sq
        BIV%generic_bed_roughness_2_next = ice%slid_beta_sq
      CASE ('Schoof2005')
        ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
        BIV%generic_bed_roughness_1      = ice%slid_alpha_sq
        BIV%generic_bed_roughness_1_prev = ice%slid_alpha_sq
        BIV%generic_bed_roughness_1_next = ice%slid_alpha_sq
        BIV%generic_bed_roughness_2      = ice%slid_beta_sq
        BIV%generic_bed_roughness_2_prev = ice%slid_beta_sq
        BIV%generic_bed_roughness_2_next = ice%slid_beta_sq
      CASE ('Zoet-Iverson')
        ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
        BIV%generic_bed_roughness_1      = ice%till_friction_angle
        BIV%generic_bed_roughness_1_prev = ice%till_friction_angle
        BIV%generic_bed_roughness_1_next = ice%till_friction_angle
        BIV%generic_bed_roughness_2      = 0._dp
        BIV%generic_bed_roughness_2_prev = 0._dp
        BIV%generic_bed_roughness_2_next = 0._dp
      CASE DEFAULT
        CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END SELECT

    ! Initialise chosen basal inversion model
    ! =======================================

    SELECT CASE (C%choice_bed_roughness_nudging_method)
      CASE ('H_dHdt_flowline')
        CALL initialise_basal_inversion_H_dHdt_flowline( mesh, ice, BIV, region_name)
      CASE ('Pien2023')
        CALL initialise_basal_inversion_Pien2023( mesh, ice, BIV, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_bed_roughness_nudging_method "' // TRIM( C%choice_bed_roughness_nudging_method) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion

  !===

  SUBROUTINE run_basal_inversion_Pien2023( mesh, grid_smooth, ice, refgeo, BIV)
    ! Run the basal inversion model based on Pien2023

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_basal_inversion),          INTENT(INOUT) :: BIV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_basal_inversion_Pien2023'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: mask
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt, dC2_dt
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dC1_dt_smoothed, dC2_dt_smoothed
    REAL(dp)                                           :: misfit, reg_1st, reg_2nd, BIVgeo_Pien2023_H0, BIVgeo_Pien2023_tau, beta_min
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHs_dx, dHs_dy, abs_grad_Hs
    REAL(dp)                                           :: fg_exp_mod

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( mask(            mesh%vi1:mesh%vi2), source = 0      )
    ALLOCATE( dC1_dt(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dC2_dt(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dC1_dt_smoothed( mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dC2_dt_smoothed( mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHs_dx(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( dHs_dy(          mesh%vi1:mesh%vi2), source = 0._dp  )
    ALLOCATE( abs_grad_Hs(     mesh%vi1:mesh%vi2), source = 0._dp  )

  ! == Calculate bed roughness rates of changes
  ! ===========================================

    DO vi = mesh%vi1, mesh%vi2

      ! Determine whether bed roughness should be
      ! updated by inversion or by extrapolation

      ! Only perform the inversion on fully grounded vertices
      IF ( ice%fraction_gr( vi) == 1._dp) THEN

        ! Perform the inversion here
        mask( vi) = 2

        ! Surface elevation misfit
        misfit = ice%Hi( vi) - refgeo%Hi( vi)

        ! Is it improving already?
        IF (ice%dHi_dt( vi)*misfit < 0._dp) THEN
          ! Yes, so leave this vertex alone
          CYCLE
        END IF

      ELSE

        ! Extrapolate here
        mask( vi) = 1
        CYCLE

      END IF

      ! == Regularisation
      ! =================

      BIVgeo_Pien2023_H0  = 200._dp
      BIVgeo_Pien2023_tau = 500._dp

      ! First regularisation term
      reg_1st = misfit / BIVgeo_Pien2023_H0 / BIVgeo_Pien2023_tau

      ! Second regularisation term
      reg_2nd = ice%dHi_dt( vi) * 2._dp / BIVgeo_Pien2023_H0

      ! == Prevent over adjustment
      ! ==========================

      beta_min = 1000._dp * MIN( 1._dp, MAX( 0._dp, 1._dp - ice%Hi( vi) / 200._dp))

      IF ((reg_1st + reg_2nd) > 0._dp .AND. ice%basal_friction_coefficient( vi) <= beta_min) THEN
        ! It's okay, we'll get'em next time
        CYCLE
      END IF

      ! Calculate bed roughness rates of change
      ! =======================================

      dC1_dt( vi) = -1._dp * BIV%generic_bed_roughness_1( vi) * (reg_1st + reg_2nd)
      dC2_dt( vi) = -1._dp * BIV%generic_bed_roughness_2( vi) * (reg_1st + reg_2nd)

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

    ! Scale (reduce) bed roughness rate of change for partially grounded, steep-sloped areas
    DO vi = mesh%vi1, mesh%vi2

      ! Ice margin and grounding lines
      IF (ice%mask_grounded_ice( vi)) THEN

        ! Strengthen the effect of grounded fractions for steep slopes
        fg_exp_mod = MIN( 1.0_dp, MAX( 0._dp, MAX( 0._dp, abs_grad_Hs( vi) - 0.02_dp) / (0.06_dp - 0.02_dp) ))

        ! Scale based on grounded fraction
        dC1_dt( vi) = dC1_dt( vi) * ice%fraction_gr( vi) ** (1._dp + fg_exp_mod)
        dC2_dt( vi) = dC2_dt( vi) * ice%fraction_gr( vi) ** (1._dp + fg_exp_mod)

      END IF

    END DO

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
    DEALLOCATE( dC1_dt         )
    DEALLOCATE( dC2_dt         )
    DEALLOCATE( dC1_dt_smoothed)
    DEALLOCATE( dC2_dt_smoothed)
    DEALLOCATE( dHs_dx         )
    DEALLOCATE( dHs_dy         )
    DEALLOCATE( abs_grad_Hs    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_basal_inversion_Pien2023

  SUBROUTINE initialise_basal_inversion_Pien2023( mesh, ice, BIV, region_name)
    ! Initialise the basal inversion model based on Pien2023

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_basal_inversion),          INTENT(INOUT) :: BIV
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion_Pien2023'
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

  END SUBROUTINE initialise_basal_inversion_Pien2023

END MODULE basal_inversion_main
