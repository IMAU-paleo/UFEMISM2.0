MODULE laddie_thickness

  ! Thickness routines for the laddie model

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
  USE laddie_physics                                         , ONLY: compute_melt_rate, compute_entrainment, &
                                                                     compute_freezing_temperature, compute_buoyancy
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, map_H_a_b, map_H_a_c

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_H_npx( mesh, ice, ocean, laddie, npxref, npx, dt)
    ! Integrate H by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_H_npx'
    INTEGER                                               :: vi
    REAL(dp)                                              :: dHdt
 
    ! Add routine to path
    CALL init_routine( routine_name)

    CALL compute_divQH( mesh, laddie, npxref)

    ! Compute freezing temperature
    CALL compute_freezing_temperature( mesh, ice, laddie, npxref)

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, npxref%H)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, ice, laddie, npx, npxref%H)

    ! Compute melt rate
    CALL compute_melt_rate( mesh, ice, laddie, npxref, npxref%H)
     
    ! Compute entrainment                                    
    CALL compute_entrainment( mesh, ice, laddie, npxref, npxref%H)

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN

        ! Get first guess at dHdt
        dHdt = -laddie%divQH( vi) + laddie%melt( vi) + laddie%entr( vi)

        ! First guess at H_n
        npx%H( vi) = laddie%now%H( vi) + dHdt * dt

        ! If H_n < Hmin, enhance entrainment to ensure H_n >= Hmin
        laddie%entr_dmin( vi) = MAX( C%laddie_thickness_minimum - npx%H( vi), 0.0_dp) / dt

        ! If H_n > Hmax, suppress entrainment to ensure H_n <= Hmax
        laddie%entr( vi) = laddie%entr( vi) + MIN( C%laddie_thickness_maximum - npx%H( vi), 0.0_dp) / dt

        ! Update detrainment. Shouldn't matter but just in case
        laddie%detr( vi) = - MIN(laddie%entr( vi),0.0_dp)

        ! Get actual dHdt
        dHdt = -laddie%divQH( vi) + laddie%melt( vi) + laddie%entr( vi) + laddie%entr_dmin( vi)

        ! Get actual H_n
        npx%H( vi) = laddie%now%H( vi) + dHdt * dt

      END IF !(laddie%mask_a( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2

    ! Map H to b grid and c grid
    CALL map_H_a_b( mesh, laddie, npx%H, npx%H_b)
    CALL map_H_a_c( mesh, laddie, npx%H, npx%H_c)

    CALL check_for_NaN_dp_1D( npx%H, 'H_lad')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_H_npx

  SUBROUTINE compute_divQH( mesh, laddie, npx)
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

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQH'
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: H_tot
    INTEGER                                               :: vi, ci, vj, ei
    REAL(dp)                                              :: D_x, D_y, D, u_perp
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_gr_a_tot, mask_oc_a_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL gather_to_all_dp_1D( npx%U_c, U_c_tot)
    CALL gather_to_all_dp_1D( npx%V_c, V_c_tot)
    CALL gather_to_all_dp_1D( npx%H, H_tot)
    CALL gather_to_all_logical_1D( laddie%mask_gr_a, mask_gr_a_tot)
    CALL gather_to_all_logical_1D( laddie%mask_oc_a, mask_oc_a_tot)

    ! Initialise with zeros
    laddie%divQH = 0.0_dp

    ! == Loop over vertices ==
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      IF (laddie%mask_a( vi)) THEN
        ! Initialise

        ! Loop over all connections of vertex vi
        DO ci = 1, mesh%nC( vi)

          ! Connection ci from vertex vi leads through edge ei to vertex vj
          vj = mesh%C(  vi,ci)

          IF (vj == 0) CYCLE

          ! Skip connection if neighbour is grounded. No flux across grounding line
          ! Can be made more flexible when accounting for partial cells (PMP instead of FCMP)
          IF (mask_gr_a_tot( vj)) CYCLE

          ! Get edge
          ei = mesh%VE( vi,ci)

          ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
          D_x = mesh%V( vj,1) - mesh%V( vi,1)
          D_y = mesh%V( vj,2) - mesh%V( vi,2)
          D   = SQRT( D_x**2 + D_y**2)
          u_perp = U_c_tot( ei) * D_x/D + V_c_tot( ei) * D_y/D

          ! Calculate momentum divergence
          ! =============================
          ! Upwind:
          ! u_perp > 0: flow is exiting this vertex into vertex vj
          IF (u_perp > 0) THEN
            laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * H_tot( vi) / mesh%A( vi)
          ! u_perp < 0: flow is entering this vertex from vertex vj
          ELSE
            IF (mask_oc_a_tot( vj)) THEN
              ! Apply dH/dx = 0 in case of inflow from open ocean
              !laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * H_tot( vi) / mesh%A( vi)
              ! No inflow
              CYCLE
            ELSE
              laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * H_tot( vj) / mesh%A( vi)
            END IF
          END IF

        END DO ! DO ci = 1, mesh%nC( vi)
       
      END IF ! (laddie%mask_a( vi))

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQH

END MODULE laddie_thickness

