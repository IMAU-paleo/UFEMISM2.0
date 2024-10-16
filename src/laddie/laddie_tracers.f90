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

  SUBROUTINE compute_TS_npx( mesh, ice, laddie, npxref, npx, dt, incldiff)
    ! Integrate T and S by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(IN)    :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref
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
              - laddie%detr( vi) * npxref%T( vi)

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
              - laddie%detr( vi) * npxref%S( vi)

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

  SUBROUTINE compute_divQTS( mesh, laddie, npx, U_c, V_c, Hstar, mask_a, mask_gr_a)
    ! Calculate the layer flux divergence matrix M_divQ using an upwind scheme
    !
    ! The vertically averaged ice flux divergence represents the net ice volume (which,
    ! assuming constant density, is proportional to the ice mass) entering each Voronoi
    ! cell per unit time. This is found by calculating the ice fluxes through each
    ! shared Voronoi cell boundary, using an upwind scheme: if ice flows from vertex vi
    ! to vertex vj, the flux is found by multiplying the velocity at their shared
    ! boundary u_c with the ice thickness at vi (and, of course, the length L_c of the
    ! shared boundary). If instead it flows from vj to vi, u_c is multiplied with the
    ! ice thickness at vj.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(IN)    :: U_c, V_c
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar
    LOGICAL, DIMENSION(mesh%vi1:mesh%vi2),  INTENT(IN)    :: mask_a, mask_gr_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQTS'
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: T_tot, S_tot, Hstar_tot
    INTEGER                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                               :: ti, ci, ei, tj, vi, vj, vi1, vi2, i, j, e, k
    REAL(dp)                                              :: A_i, L_c, D_x, D_y, D, u_perp
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_a_tot
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_gr_a_tot
    LOGICAL                                               :: isbound

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL gather_to_all_dp_1D( U_c, U_c_tot)
    CALL gather_to_all_dp_1D( V_c, V_c_tot)
    CALL gather_to_all_dp_1D( Hstar, Hstar_tot)
    CALL gather_to_all_dp_1D( npx%T, T_tot)
    CALL gather_to_all_dp_1D( npx%S, S_tot)
    CALL gather_to_all_logical_1D( mask_a, mask_a_tot)
    CALL gather_to_all_logical_1D( mask_gr_a, mask_gr_a_tot)

    ! Initialise with zeros
    laddie%divQT = 0.0_dp
    laddie%divQS = 0.0_dp

    ! == Loop over vertices ==
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      IF (mask_a( vi)) THEN
        ! Initialise

        ! Loop over all connections of vertex vi
        DO ci = 1, mesh%nC( vi)

          ! Connection ci from vertex vi leads through edge ei to vertex vj
          ei = mesh%VE( vi,ci)
          vj = mesh%C(  vi,ci)

          ! Skip connection if neighbour is grounded. No flux across grounding line
          ! Can be made more flexible when accounting for partial cells (PMP instead of FCMP)
          IF (mask_gr_a_tot( vj)) CYCLE

          ! The Voronoi cell of vertex vi has area A_i
          A_i = mesh%A( vi)

          ! The shared Voronoi cell boundary section between the Voronoi cells
          ! of vertices vi and vj has length L_c
          L_c = mesh%Cw( vi,ci)

          ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
          D_x = mesh%V( vj,1) - mesh%V( vi,1)
          D_y = mesh%V( vj,2) - mesh%V( vi,2)
          D   = SQRT( D_x**2 + D_y**2)
          u_perp = U_c_tot( ei) * D_x/D + V_c_tot( ei) * D_y/D

          ! Calculate momentum divergence
          ! =============================
          ! Upwind:
          ! u_perp > 0: flow is exiting this vertex into vertex vj
          laddie%divQT( vi) = laddie%divQT( vi) + L_c * MAX( 0._dp, u_perp) * Hstar_tot( vi) * T_tot( vi) / A_i
          laddie%divQS( vi) = laddie%divQS( vi) + L_c * MAX( 0._dp, u_perp) * Hstar_tot( vi) * S_tot( vi) / A_i
          ! u_perp < 0: flow is entering this vertex from vertex vj
          laddie%divQT( vi) = laddie%divQT( vi) + L_c * MIN( 0._dp, u_perp) * Hstar_tot( vj) * T_tot( vj) / A_i
          laddie%divQS( vi) = laddie%divQS( vi) + L_c * MIN( 0._dp, u_perp) * Hstar_tot( vj) * S_tot( vj) / A_i
          ! Centered:

        END DO ! DO ci = 1, mesh%nC( vi)
       
      END IF ! (mask_a( vi))

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQTS

END MODULE laddie_tracers

