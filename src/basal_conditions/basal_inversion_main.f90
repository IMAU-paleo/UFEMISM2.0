MODULE basal_inversion_main

  ! Contains all the routines for managing the basal inversion model

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE basal_inversion_types                                  , ONLY: type_basal_inversion
  USE region_types                                           , ONLY: type_model_region
  USE basal_inversion_H_dHdt_flowline                        , ONLY: initialise_basal_inversion_H_dHdt_flowline, run_basal_inversion_H_dHdt_flowline

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
      CASE DEFAULT
        CALL crash('unknown choice_bed_roughness_nudging_method "' // TRIM( C%choice_bed_roughness_nudging_method) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion

END MODULE basal_inversion_main
