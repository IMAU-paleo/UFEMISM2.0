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
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  USE laddie_physics                                         , ONLY: compute_melt_rate, compute_entrainment, &
                                                                     compute_freezing_temperature, compute_buoyancy, compute_subglacial_discharge
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, map_H_a_b, map_H_a_c
  use mesh_halo_exchange, only: exchange_halos
  use checksum_mod, only: checksum

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE compute_H_npx( mesh, laddie, npx_old, npx_ref, npx_new, time, dt)
    ! Integrate H by time step dt

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx_old   ! Old time step
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx_ref   ! Reference time step for RHS terms
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx_new   ! New timestep as output
    REAL(dp),                               INTENT(IN)    :: time
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_H_npx'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute thickness divergence
    CALL compute_divQH( mesh, laddie, npx_ref)

    ! Compute freezing temperature
    CALL compute_freezing_temperature( mesh, laddie, npx_ref)

    ! Initialise ambient T and S
    CALL compute_ambient_TS( mesh, laddie, npx_ref%H)

    ! Compute buoyancy
    CALL compute_buoyancy( mesh, laddie, npx_ref, npx_ref%H)

    ! Compute melt rate
    CALL compute_melt_rate( mesh, laddie, npx_ref, npx_ref%H, time)

    ! Compute entrainment
    CALL compute_entrainment( mesh, laddie, npx_ref, npx_ref%H)

    ! Do integration
    CALL integrate_H( mesh, laddie, npx_old, npx_new, dt)

    ! Map new values of H to b grid and c grid
    CALL map_H_a_b( mesh, laddie, npx_new%H, npx_new%H_b)
    CALL map_H_a_c( mesh, laddie, npx_new%H, npx_new%H_c)
    call checksum( npx_new%H_b, 'npx_new%H_b', mesh%pai_Tri)
    call checksum( npx_new%H_c, 'npx_new%H_c', mesh%pai_E)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_H_npx

  subroutine integrate_H( mesh, laddie, npx_old, npx_new, dt)
    ! Do the actual computation of npx%H

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(inout) :: laddie
    type(type_laddie_timestep),             intent(in   ) :: npx_old   ! Old timestep as input
    type(type_laddie_timestep),             intent(inout) :: npx_new   ! New timestep as output
    real(dp),                               intent(in)    :: dt

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'integrate_H'
    integer                                               :: vi
    real(dp)                                              :: dHdt

    ! Add routine to path
    call init_routine( routine_name)

    ! Loop over vertices
    do vi = mesh%vi1, mesh%vi2
      if (laddie%mask_a( vi)) then

        ! Get first guess at dHdt
        dHdt = -laddie%divQH( vi) + laddie%melt( vi) + laddie%entr( vi) + laddie%SGD( vi)

        ! First guess at H_n
        npx_new%H( vi) = npx_old%H( vi) + dHdt * dt

        ! If H_n < Hmin, enhance entrainment to ensure H_n >= Hmin
        laddie%entr_dmin( vi) = MAX( C%laddie_thickness_minimum - npx_new%H( vi), 0.0_dp) / dt

        ! If H_n > Hmax, suppress entrainment to ensure H_n <= Hmax
        laddie%entr( vi) = laddie%entr( vi) + MIN( C%laddie_thickness_maximum - npx_new%H( vi), 0.0_dp) / dt

        ! Prevent strong entr_dmin and strong detrainment
        if (laddie%entr_dmin(vi) > 0) then
          laddie%entr( vi) = MAX(laddie%entr( vi), 0.0_dp)
        end if

        ! Update detrainment. Shouldn't matter but just in case
        laddie%detr( vi) = - MIN(laddie%entr( vi),0.0_dp)

        ! Get actual dHdt
        dHdt = -laddie%divQH( vi) + laddie%melt( vi) + laddie%entr( vi) + laddie%entr_dmin( vi) + laddie%SGD( vi)

        ! Get actual H_n
        npx_new%H( vi) = npx_old%H( vi) + dHdt * dt

      end if !(laddie%mask_a( vi)) THEN
    end do !vi = mesh%vi, mesh%v2
    call checksum( laddie%entr_dmin, 'laddie%entr_dmin', mesh%pai_V)
    call checksum( laddie%entr     , 'laddie%entr     ', mesh%pai_V)
    call checksum( laddie%detr     , 'laddie%detr     ', mesh%pai_V)
    call checksum( npx_new%H       , 'npx_new%H       ', mesh%pai_V)

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
    INTEGER                                               :: vi, ci, vj, ei
    REAL(dp)                                              :: u_perp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    call exchange_halos( mesh, npx%U_c)
    call exchange_halos( mesh, npx%V_c)
    call exchange_halos( mesh, npx%H  )
    call checksum( npx%U_c, 'npx%U_c', mesh%pai_E)
    call checksum( npx%V_c, 'npx%V_c', mesh%pai_E)
    call checksum( npx%H  , 'npx%H  ', mesh%pai_V)

    ! Initialise with zeros
    laddie%divQH( mesh%vi1:mesh%vi2) = 0.0_dp

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
          IF (laddie%mask_gr_a( vj)) CYCLE

          ! Get edge
          ei = mesh%VE( vi,ci)

          ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
          u_perp = npx%U_c( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + npx%V_c( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

          ! Calculate upwind momentum divergence
          ! =============================
          ! u_perp > 0: flow is exiting this vertex into vertex vj
          IF (u_perp > 0) THEN
            laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * npx%H( vi) / mesh%A( vi)
          ! u_perp < 0: flow is entering this vertex from vertex vj
          ELSE
            IF (laddie%mask_oc_a( vj)) THEN
              CYCLE ! No inflow
              ! TODO fix boundary condition for inflow
              ! laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * npx%H( vi) / mesh%A( vi)
            ELSE
              laddie%divQH( vi) = laddie%divQH( vi) + mesh%Cw( vi, ci) * u_perp * npx%H( vj) / mesh%A( vi)
            END IF
          END IF

        END DO ! DO ci = 1, mesh%nC( vi)

      END IF ! (laddie%mask_a( vi))

    END DO ! DO vi = mesh%vi1, mesh%vi2
    call checksum( laddie%divQH, 'laddie%divQH', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQH

END MODULE laddie_thickness

