MODULE laddie_tracers

  ! Tracer routines for the laddie model

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
  USE reallocate_mod                                         , ONLY: reallocate_bounds

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_TS_np1( mesh, ice, laddie, dt)
    ! Integrate T and S by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_TS_np1'
    INTEGER                                               :: vi
    REAL(dp)                                              :: dHTdt
    REAL(dp)                                              :: dHSdt
    REAL(dp)                                              :: HT_next
    REAL(dp)                                              :: HS_next
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Temperature integration ==

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_floating_ice( vi)) THEN

        ! Get dHT_dt
        dHTdt = -laddie%divQT( vi) &
              + laddie%melt( vi) * laddie%T_base( vi) &
              + laddie%entr( vi) * laddie%T_amb( vi) &
              + laddie%diffT( vi) 

        ! HT_n = HT_n + dHT_dt * dt
        HT_next = laddie%T( vi)*laddie%H( vi) + dHTdt * dt

        laddie%T_next( vi) = HT_next / laddie%H_next( vi)

      END IF !(ice%mask_floating_ice( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2

    ! == Salinity integration ==

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_floating_ice( vi)) THEN

        ! Get dHS_dt
        dHSdt = -laddie%divQS( vi) &
              + laddie%entr( vi) * laddie%S_amb( vi) &
              + laddie%diffS( vi)

        ! HS_n = HS_n + dHS_dt * dt
        HS_next = laddie%S( vi)*laddie%H( vi) + dHSdt * dt

        laddie%S_next( vi) = HS_next / laddie%H_next( vi)

      END IF !(ice%mask_floating_ice( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_TS_np1

  SUBROUTINE compute_diffTS( mesh, ice, laddie, Hstar)
    ! Compute horizontal diffusion of heat and salt

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_diffTS'
    INTEGER                                               :: vi
    REAL(dp), DIMENSION(mesh%nV)                          :: Hstar_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather
    CALL gather_to_all_dp_2D( Hstar, Hstar_tot)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_diffTS

END MODULE laddie_tracers

