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
  USE mpi_distributed_memory                                 , ONLY: gather_to_all

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_TS_npx( mesh, ice, laddie, npxref, npx, Hstar, dt, include_diffusive_terms)
    ! Integrate T and S by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar
    REAL(dp),                               INTENT(IN)    :: dt
    LOGICAL,                                INTENT(IN)    :: include_diffusive_terms

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_TS_npx'
    INTEGER                                               :: vi
    REAL(dp)                                              :: dHTdt, dHSdt, HT_next, HS_next
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute divergence of heat and salt
    CALL compute_divQTS( mesh, laddie, npx, Hstar)

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN

        ! == Get time derivatives ==

        ! Time-derivative heat equation
        dHTdt = -laddie%divQT( vi) &
               + laddie%melt( vi) * laddie%T_base( vi) &
               + MAX(0.0_dp,laddie%entr( vi)) * laddie%T_amb( vi) &
               + laddie%entr_dmin( vi) * laddie%T_amb( vi) &
               - laddie%detr( vi) * npxref%T( vi)

        ! Time-derivative salt equation
        dHSdt = -laddie%divQS( vi) &
               + MAX(0.0_dp,laddie%entr( vi)) * laddie%S_amb( vi) &
               + laddie%entr_dmin( vi) * laddie%S_amb( vi) &
               - laddie%detr( vi) * npxref%S( vi)

        ! Add diffusive terms if requested
        IF (include_diffusive_terms) THEN
          dHTdt = dHTdt + laddie%diffT( vi)
          dHSdt = dHSdt + laddie%diffS( vi)
        END IF

        ! == Apply time-integration ==

        ! HT_n = HT + dHT_dt * dt
        HT_next = laddie%now%T( vi)*laddie%now%H( vi) + dHTdt * dt
        npx%T( vi) = HT_next / npx%H( vi)

        ! HS_n = HS + dHS_dt * dt
        HS_next = laddie%now%S( vi)*laddie%now%H( vi) + dHSdt * dt
        npx%S( vi) = HS_next / npx%H( vi)
      END IF !(laddie%mask_a( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_TS_npx

  SUBROUTINE compute_diffTS( mesh, ice, laddie, npxref)
    ! Compute horizontal diffusion of heat and salt

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_diffTS'
    INTEGER                                               :: vi, vj, ci, ei
    REAL(dp)                                              :: Kh
    REAL(dp), DIMENSION(mesh%nV)                          :: T_tot, S_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    REAL(dp), DIMENSION(mesh%nE)                          :: H_c_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather
    CALL gather_to_all( npxref%T, T_tot)
    CALL gather_to_all( npxref%S, S_tot)
    CALL gather_to_all( laddie%mask_a, mask_a_tot)
    CALL gather_to_all( npxref%H_c, H_c_tot)

    ! Initialise at 0
    laddie%diffT = 0.0_dp
    laddie%diffS = 0.0_dp

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2

      IF (laddie%mask_a( vi)) THEN

        ! Loop over connected vertices
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi, ci)
          ei = mesh%VE( vi, ci)
          ! Can simply skip non-floating vertices to ensure d/dx = d/dy = 0 at boundaries
          IF (mask_a_tot( vj)) THEN
            ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section

            Kh = C%laddie_diffusivity 

            laddie%diffT( vi) = laddie%diffT( vi) + (T_tot( vj) - T_tot( vi)) * Kh * H_c_tot( ei) / mesh%A( vi) * mesh%Cw( vi, ci)/mesh%D( vi, ci)
            laddie%diffS( vi) = laddie%diffS( vi) + (S_tot( vj) - S_tot( vi)) * Kh * H_c_tot( ei) / mesh%A( vi) * mesh%Cw( vi, ci)/mesh%D( vi, ci)
          END IF
        END DO

      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_diffTS

  SUBROUTINE compute_divQTS( mesh, laddie, npx, Hstar)
    ! Divergence of heat / salt

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQTS'
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: T_tot, S_tot, Hstar_tot
    INTEGER                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                               :: ti, ci, ei, tj, vi, vj, vi1, vi2, i, j, e, k
    REAL(dp)                                              :: u_perp
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_gr_a_tot, mask_oc_a_tot
    LOGICAL                                               :: isbound

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL gather_to_all( npx%U_c, U_c_tot)
    CALL gather_to_all( npx%V_c, V_c_tot)
    CALL gather_to_all( Hstar, Hstar_tot)
    CALL gather_to_all( npx%T, T_tot)
    CALL gather_to_all( npx%S, S_tot)
    CALL gather_to_all( laddie%mask_a, mask_a_tot)
    CALL gather_to_all( laddie%mask_gr_a, mask_gr_a_tot)
    CALL gather_to_all( laddie%mask_oc_a, mask_oc_a_tot)

    ! Initialise with zeros
    laddie%divQT = 0.0_dp
    laddie%divQS = 0.0_dp

    ! == Loop over vertices ==
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      IF (laddie%mask_a( vi)) THEN

        ! Loop over all connections of vertex vi
        DO ci = 1, mesh%nC( vi)

          ! Connection ci from vertex vi leads through edge ei to vertex vj
          vj = mesh%C(  vi,ci)

          ! Skip connection if neighbour is grounded. No flux across grounding line
          ! Can be made more flexible when accounting for partial cells (PMP instead of FCMP)
          IF (mask_gr_a_tot( vj)) CYCLE

          ei = mesh%VE( vi,ci)

          ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
          u_perp = U_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + V_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

          ! Calculate upwind momentum divergence
          ! =============================
          ! u_perp > 0: flow is exiting this vertex into vertex vj
          IF (u_perp > 0) THEN
            laddie%divQT( vi) = laddie%divQT( vi) + mesh%Cw( vi, ci) * u_perp * Hstar_tot( vi) * T_tot( vi) / mesh%A( vi)
            laddie%divQS( vi) = laddie%divQS( vi) + mesh%Cw( vi, ci) * u_perp * Hstar_tot( vi) * S_tot( vi) / mesh%A( vi)
          ! u_perp < 0: flow is entering this vertex from vertex vj
          ELSE
            IF (mask_oc_a_tot( vj)) THEN
              CYCLE ! no inflow
              ! TODO fix boundary condition inflow
              laddie%divQT( vi) = laddie%divQT( vi) + mesh%Cw( vi, ci) * u_perp * Hstar_tot( vi) * T_tot( vi) / mesh%A( vi)
              laddie%divQS( vi) = laddie%divQS( vi) + mesh%Cw( vi, ci) * u_perp * Hstar_tot( vi) * S_tot( vi) / mesh%A( vi)
            ELSE
              laddie%divQT( vi) = laddie%divQT( vi) + mesh%Cw( vi, ci) * u_perp * Hstar_tot( vj) * T_tot( vj) / mesh%A( vi)
              laddie%divQS( vi) = laddie%divQS( vi) + mesh%Cw( vi, ci) * u_perp * Hstar_tot( vj) * S_tot( vj) / mesh%A( vi)
            END IF
          END IF

        END DO ! DO ci = 1, mesh%nC( vi)
       
      END IF ! (laddie%mask_a( vi))

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQTS

END MODULE laddie_tracers

