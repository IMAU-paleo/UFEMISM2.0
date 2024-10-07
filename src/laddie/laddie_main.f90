MODULE laddie_main

  ! The main laddie model module.

! ===== Preamble =====
! ====================
    
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, allocate_laddie_model
  USE laddie_thickness                                       , ONLY: compute_H_np1 
  USE laddie_velocity                                        , ONLY: compute_UV_np1
  USE laddie_tracers                                         , ONLY: compute_TS_np1

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE run_laddie_model( mesh, ice, ocean, laddie, time)
    ! Run the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_laddie_model'
    INTEGER                                               :: vi
    REAL(dp)                                              :: tl               ! [s] Laddie time
    REAL(dp)                                              :: dt               ! [s] Laddie time step
 
    ! Add routine to path
    CALL init_routine( routine_name)
 
    ! == Preparation ==
    ! =================

    ! Get time step
    tl = 0.0_dp
    dt = C%dt_laddie

    ! Update masks
    ! TODO

    ! Extrapolate new cells
    ! TODO

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, laddie, ocean, ice)

    ! == Main time loop ==
    ! ====================

    DO WHILE (tl <= C%time_duration_laddie * sec_per_day)
      ! Integrate H 1 time step
      CALL compute_H_np1( mesh, ice, laddie, dt)

      ! Integrate U and V 1 time step
      CALL compute_UV_np1( mesh, ice, laddie, dt)

      ! Integrate T and S 1 time step
      CALL compute_TS_np1( mesh, ice, laddie, dt)

      ! Update secondary fields
      ! TODO

      ! == Move time ==
      ! Increase laddie time
      tl = tl + C%dt_laddie

      ! Move main variables by 1 time step
      DO vi = mesh%vi1, mesh%vi2
        laddie%H( vi) = laddie%H_next( vi)
      END DO

      ! Display or save fields
      ! TODO
      ! IF (par%master) THEN
      !   WRITE( *, "(F8.3)") MAXVAL(laddie%H)
      ! END IF     


    END DO !DO WHILE (tl <= C%time_duration_laddie)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_laddie_model

  SUBROUTINE initialise_laddie_model( mesh, laddie, ocean, ice)
    ! Initialise the laddie model

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_laddie_model'
    INTEGER                                               :: vi
    REAL(dp), PARAMETER                                   :: H_init = 10.0_dp ! [m]    Initial thickness
    REAL(dp), PARAMETER                                   :: T_off  = 0.0_dp  ! [degC] Initial temperature offset
    REAL(dp), PARAMETER                                   :: S_off  = -0.1_dp ! [PSU]  Initial salinity offset
 
    ! Add routine to path
    CALL init_routine( routine_name)
 
    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising LADDIE model...'

    ! Allocate variables
    CALL allocate_laddie_model( mesh, laddie)

    ! Layer thickness 
    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         laddie%H( vi)      = H_init
         laddie%H_prev( vi) = H_init
         laddie%H_next( vi) = H_init
       END IF
    END DO

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, laddie, ocean, ice)

    ! Initialise main T and S
    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         laddie%T( vi)      = laddie%T_amb( vi) + T_off
         laddie%T_prev( vi) = laddie%T_amb( vi) + T_off
         laddie%T_next( vi) = laddie%T_amb( vi) + T_off
         laddie%S( vi)      = laddie%S_amb( vi) + S_off
         laddie%S_prev( vi) = laddie%S_amb( vi) + S_off
         laddie%S_next( vi) = laddie%S_amb( vi) + S_off
       END IF
    END DO
 
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_laddie_model

END MODULE laddie_main

