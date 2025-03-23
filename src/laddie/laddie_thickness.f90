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
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  USE laddie_physics                                         , ONLY: compute_melt_rate, compute_entrainment, &
                                                                     compute_freezing_temperature, compute_buoyancy
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, map_H_a_b, map_H_a_c
  use mesh_utilities, only: average_over_domain

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE compute_H_npx( mesh, ice, ocean, laddie, npxref, npx, time, dt)
    ! Integrate H by time step dt

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref    ! Reference time step as input
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx       ! New timestep as output
    REAL(dp),                               INTENT(IN)    :: time
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_H_npx'
    real(dp) :: d_av

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute thickness divergence
    CALL compute_divQH( mesh, laddie, npxref)

    ! Compute freezing temperature
    CALL compute_freezing_temperature( mesh, ice, laddie, npxref)

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, ice, ocean, laddie, npxref%H)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, ice, laddie, npx, npxref%H)

    ! Compute melt rate
    CALL compute_melt_rate( mesh, ice, laddie, npxref, npxref%H, time)

    ! Compute entrainment
    CALL compute_entrainment( mesh, ice, laddie, npxref, npxref%H)

    ! Do integration
    CALL integrate_H( mesh, ice, laddie, npx, dt)

    ! Map new values of H to b grid and c grid
    CALL map_H_a_b( mesh, laddie, npx%H, npx%H_b)
    CALL map_H_a_c( mesh, laddie, npx%H, npx%H_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_H_npx

  subroutine integrate_H( mesh, ice, laddie, npx, dt)
    ! Do the actual computation of npx%H

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_laddie_model),                intent(inout) :: laddie
    type(type_laddie_timestep),             intent(inout) :: npx       ! New timestep as output
    real(dp),                               intent(in)    :: dt

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'integrate_H'
    integer                                               :: vi
    real(dp)                                              :: dHdt
    real(dp) :: d_av

    ! Add routine to path
    call init_routine( routine_name)

    ! Loop over vertices
    do vi = mesh%vi1, mesh%vi2
      if (laddie%mask_a( vi)) then

        ! Get first guess at dHdt
        dHdt = -laddie%divQH( vi) + laddie%melt( vi) + laddie%entr( vi)

        ! First guess at H_n
        npx%H( vi) = laddie%now%H( vi) + dHdt * dt

        ! If H_n < Hmin, enhance entrainment to ensure H_n >= Hmin
        laddie%entr_dmin( vi) = MAX( C%laddie_thickness_minimum - npx%H( vi), 0.0_dp) / dt

        ! If H_n > Hmax, suppress entrainment to ensure H_n <= available water column thickness
        laddie%entr( vi) = laddie%entr( vi) + MIN( ice%Hib( vi)-ice%Hb( vi) - npx%H( vi), 0.0_dp) / dt

        ! Prevent strong entr_dmin and strong detrainment
        if (laddie%entr_dmin(vi) > 0) then
          laddie%entr( vi) = MAX(laddie%entr( vi), 0.0_dp)
        end if

        ! Update detrainment. Shouldn't matter but just in case
        laddie%detr( vi) = - MIN(laddie%entr( vi),0.0_dp)

        ! Get actual dHdt
        dHdt = -laddie%divQH( vi) + laddie%melt( vi) + laddie%entr( vi) + laddie%entr_dmin( vi)

        ! Get actual H_n
        npx%H( vi) = laddie%now%H( vi) + dHdt * dt

      end if !(laddie%mask_a( vi)) THEN
    end do !vi = mesh%vi, mesh%v2

    ! DENK DROM
    call average_over_domain( mesh, npx%H, d_av)
    if (par%primary) write(0,'(A,F12.8)') ' mean npx%H = ', d_av

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_H

  SUBROUTINE compute_divQH( mesh, laddie, npx)
    ! Divergence of layer thickness based on first order upwind scheme

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQH'
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: H_tot
    INTEGER                                               :: vi, ci, vj, ei
    REAL(dp)                                              :: u_perp
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_gr_a_tot, mask_oc_a_tot
    real(dp) :: d_av

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL gather_to_all( npx%U_c, U_c_tot)
    CALL gather_to_all( npx%V_c, V_c_tot)
    CALL gather_to_all( npx%H, H_tot)
    CALL gather_to_all( laddie%mask_gr_a, mask_gr_a_tot)
    CALL gather_to_all( laddie%mask_oc_a, mask_oc_a_tot)

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

          ! Skip connection if neighbour is grounded. No flux across grounding line
          ! Can be made more flexible when accounting for partial cells (PMP instead of FCMP)
          IF (mask_gr_a_tot( vj)) CYCLE

          ! Get edge
          ei = mesh%VE( vi,ci)

          ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
          u_perp = U_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + V_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

          ! Calculate upwind momentum divergence
          ! =============================
          ! u_perp > 0: flow is exiting this vertex into vertex vj
          IF (u_perp > 0) THEN
            laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * H_tot( vi) / mesh%A( vi)
          ! u_perp < 0: flow is entering this vertex from vertex vj
          ELSE
            IF (mask_oc_a_tot( vj)) THEN
              CYCLE ! No inflow
              ! TODO fix boundary condition for inflow
              ! laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * H_tot( vi) / mesh%A( vi)
            ELSE
              laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * H_tot( vj) / mesh%A( vi)
            END IF
          END IF

        END DO ! DO ci = 1, mesh%nC( vi)

      END IF ! (laddie%mask_a( vi))

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! DENK DROM
    call average_over_domain( mesh, laddie%divQH, d_av)
    if (par%primary) write(0,'(A,F12.8)') ' mean divQH = ', d_av

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQH

END MODULE laddie_thickness

