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
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D

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
    INTEGER                                               :: vj
    INTEGER                                               :: ci
    REAL(dp), DIMENSION(mesh%nV)                          :: T_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: S_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_floating_ice_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather
    CALL gather_to_all_dp_1D( laddie%T, T_tot)
    CALL gather_to_all_dp_1D( laddie%S, S_tot)
    CALL gather_to_all_logical_1D( ice%mask_floating_ice, mask_floating_ice_tot)

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_floating_ice( vi)) THEN
        ! Get diffusivity parameter
        laddie%K_h( vi) = C%laddie_diffusivity
        ! TODO add scalable options

        ! Initialise at 0
        laddie%diffT( vi) = 0.0_dp
        laddie%diffS( vi) = 0.0_dp

        ! Loop over connected vertices
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi, ci)
          ! Can simply skip non-floating vertices to ensure d/dx = d/dy = 0 at boundaries
          IF (mask_floating_ice_tot( vj)) THEN
            laddie%diffT( vi) = laddie%diffT( vi) + (T_tot( vj)-laddie%T( vi)) * laddie%K_h( vi) * Hstar( vi) / mesh%A( vi)
            laddie%diffS( vi) = laddie%diffS( vi) + (S_tot( vj)-laddie%S( vi)) * laddie%K_h( vi) * Hstar( vi) / mesh%A( vi)
          END IF
        END DO

      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_diffTS

END MODULE laddie_tracers

