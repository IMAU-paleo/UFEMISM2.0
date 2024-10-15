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
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D
  USE math_utilities                                         , ONLY: check_for_NaN_dp_1D

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_TS_npx( mesh, ice, laddie, npx, dt, incldiff)
    ! Integrate T and S by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(IN)    :: laddie
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx
    REAL(dp),                               INTENT(IN)    :: dt
    LOGICAL,                                INTENT(IN)    :: incldiff

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_TS_npx'
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
      IF (laddie%mask_a( vi)) THEN

        ! Get dHT_dt
        dHTdt = -laddie%divQT( vi) &
              + laddie%melt( vi) * laddie%T_base( vi) &
              + MAX(0.0_dp,laddie%entr( vi)) * laddie%T_amb( vi) &
              + laddie%entr_dmin( vi) * laddie%T_amb( vi) &
              - laddie%detr( vi) * laddie%T_amb( vi)

        IF (incldiff) THEN
          dHTdt = dHTdt + laddie%diffT( vi)
        END IF

        ! HT_n = HT_n + dHT_dt * dt
        HT_next = laddie%now%T( vi)*laddie%now%H( vi) + dHTdt * dt

        npx%T( vi) = HT_next / npx%H( vi)

      END IF !(laddie%mask_a( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2

    CALL check_for_NaN_dp_1D( npx%T, 'T_lad')

    ! == Salinity integration ==

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN

        ! Get dHS_dt
        dHSdt = -laddie%divQS( vi) &
              + MAX(0.0_dp,laddie%entr( vi)) * laddie%S_amb( vi) &
              + laddie%entr_dmin( vi) * laddie%S_amb( vi) &
              - laddie%detr( vi) * laddie%S_amb( vi)

        IF (incldiff) THEN
          dHSdt = dHSdt + laddie%diffS( vi)
        END IF

        ! HS_n = HS_n + dHS_dt * dt
        HS_next = laddie%now%S( vi)*laddie%now%H( vi) + dHSdt * dt

        npx%S( vi) = HS_next / npx%H( vi)

      END IF !(laddie%mask_a( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2

    CALL check_for_NaN_dp_1D( npx%S, 'S_lad')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_TS_npx

  SUBROUTINE compute_diffTS( mesh, ice, laddie)
    ! Compute horizontal diffusion of heat and salt

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_diffTS'
    INTEGER                                               :: vi
    INTEGER                                               :: vj
    INTEGER                                               :: ci
    REAL(dp), DIMENSION(mesh%nV)                          :: T_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: S_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather
    CALL gather_to_all_dp_1D( laddie%now%T, T_tot)
    CALL gather_to_all_dp_1D( laddie%now%S, S_tot)
    CALL gather_to_all_logical_1D( laddie%mask_a, mask_a_tot)

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN
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
          IF (mask_a_tot( vj)) THEN
            laddie%diffT( vi) = laddie%diffT( vi) + (T_tot( vj)-laddie%now%T( vi)) * laddie%K_h( vi) * laddie%now%H( vi) / mesh%A( vi)
            laddie%diffS( vi) = laddie%diffS( vi) + (S_tot( vj)-laddie%now%S( vi)) * laddie%K_h( vi) * laddie%now%H( vi) / mesh%A( vi)
          END IF
        END DO

      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_diffTS

END MODULE laddie_tracers

