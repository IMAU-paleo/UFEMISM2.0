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
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  use mesh_halo_exchange, only: exchange_halos
  use checksum_mod, only: checksum

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE compute_TS_npx( mesh, laddie, npx_old, npx_ref, npx_new, Hstar, dt, include_diffusive_terms)
    ! Integrate T and S by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx_old   ! Old time step
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx_ref   ! Reference time step for RHS terms
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx_new   ! New timestep as output
    REAL(dp), DIMENSION(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), INTENT(IN)    :: Hstar
    REAL(dp),                               INTENT(IN)    :: dt
    LOGICAL,                                INTENT(IN)    :: include_diffusive_terms

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_TS_npx'
    INTEGER                                               :: vi
    REAL(dp)                                              :: dHTdt, dHSdt, HT_next, HS_next

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute divergence of heat and salt
    CALL compute_divQTS( mesh, laddie, npx_ref, Hstar)

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN

        ! == Get time derivatives ==

        ! Time-derivative heat equation
        dHTdt = -laddie%divQT( vi) &
               + laddie%melt( vi) * laddie%T_base( vi) &
               + MAX(0.0_dp,laddie%entr( vi)) * laddie%T_amb( vi) &
               + laddie%entr_dmin( vi) * laddie%T_amb( vi) &
               - laddie%detr( vi) * npx_ref%T( vi) &
               + laddie%SGD( vi) * (freezing_lambda_2 + freezing_lambda_3*laddie%Hib( vi))

        ! Time-derivative salt equation
        dHSdt = -laddie%divQS( vi) &
               + MAX(0.0_dp,laddie%entr( vi)) * laddie%S_amb( vi) &
               + laddie%entr_dmin( vi) * laddie%S_amb( vi) &
               - laddie%detr( vi) * npx_ref%S( vi) &
               + laddie%SGD( vi) * 0._dp

        ! Add diffusive terms if requested
        IF (include_diffusive_terms) THEN
          dHTdt = dHTdt + laddie%diffT( vi)
          dHSdt = dHSdt + laddie%diffS( vi)
        END IF

        ! == Apply time-integration ==

        ! HT_n = HT + dHT_dt * dt
        HT_next = npx_old%T( vi)*npx_old%H( vi) + dHTdt * dt
        npx_new%T( vi) = HT_next / npx_new%H( vi)

        ! HS_n = HS + dHS_dt * dt
        HS_next = npx_old%S( vi)*npx_old%H( vi) + dHSdt * dt
        npx_new%S( vi) = HS_next / npx_new%H( vi)
      END IF !(laddie%mask_a( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2
    call checksum( npx_new%T, 'npx_new%T', mesh%pai_V)
    call checksum( npx_new%S, 'npx_new%S', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_TS_npx

  SUBROUTINE compute_diffTS( mesh, laddie, npxref)
    ! Compute horizontal diffusion of heat and salt

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npxref

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_diffTS'
    INTEGER                                               :: vi, vj, ci, ei
    REAL(dp)                                              :: Kh

    ! Add routine to path
    CALL init_routine( routine_name)

    call exchange_halos( mesh, npxref%T)
    call exchange_halos( mesh, npxref%S)
    call exchange_halos( mesh, npxref%H_c)

    ! Initialise at 0
    laddie%diffT( mesh%vi1:mesh%vi2) = 0.0_dp
    laddie%diffS( mesh%vi1:mesh%vi2) = 0.0_dp

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2

      IF (laddie%mask_a( vi)) THEN

        ! Loop over connected vertices
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi, ci)
          ei = mesh%VE( vi, ci)
          ! Can simply skip non-floating vertices to ensure d/dx = d/dy = 0 at boundaries
          IF (laddie%mask_a( vj)) THEN
            ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section

            Kh = C%laddie_diffusivity

            laddie%diffT( vi) = laddie%diffT( vi) + (npxref%T( vj) - npxref%T( vi)) * Kh * npxref%H_c( ei) / mesh%A( vi) * mesh%Cw( vi, ci)/mesh%D( vi, ci)
            laddie%diffS( vi) = laddie%diffS( vi) + (npxref%S( vj) - npxref%S( vi)) * Kh * npxref%H_c( ei) / mesh%A( vi) * mesh%Cw( vi, ci)/mesh%D( vi, ci)
          END IF
        END DO

      END IF
    END DO
    call checksum( laddie%diffT, 'laddie%diffT', mesh%pai_V)
    call checksum( laddie%diffS, 'laddie%diffS', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_diffTS

  SUBROUTINE compute_divQTS( mesh, laddie, npx, Hstar)
    ! Divergence of heat / salt

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx
    REAL(dp), DIMENSION(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_divQTS'
    INTEGER                                               :: ci, ei, vi, vj
    REAL(dp)                                              :: u_perp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! TODO: figure out which of these fields already had their halos exchanged
    call exchange_halos( mesh, npx%U_c)
    call exchange_halos( mesh, npx%V_c)
    ! call exchange_halos( mesh, Hstar) ! Definitely this one
    call exchange_halos( mesh, npx%T)
    call exchange_halos( mesh, npx%S)
    call exchange_halos( mesh, laddie%mask_a)
    call exchange_halos( mesh, laddie%mask_gr_a)
    call exchange_halos( mesh, laddie%mask_oc_a)

    call checksum( npx%U_c         , 'npx%U_c         ', mesh%pai_E)
    call checksum( npx%V_c         , 'npx%V_c         ', mesh%pai_E)
    call checksum( npx%T           , 'npx%T           ', mesh%pai_V)
    call checksum( npx%S           , 'npx%S           ', mesh%pai_V)
    call checksum( laddie%mask_a   , 'laddie%mask_a   ', mesh%pai_V)
    call checksum( laddie%mask_gr_a, 'laddie%mask_gr_a', mesh%pai_V)
    call checksum( laddie%mask_oc_a, 'laddie%mask_oc_a', mesh%pai_V)

    ! Initialise with zeros
    laddie%divQT( mesh%vi1:mesh%vi2) = 0.0_dp
    laddie%divQS( mesh%vi1:mesh%vi2) = 0.0_dp

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
          IF (laddie%mask_gr_a( vj)) CYCLE

          ei = mesh%VE( vi,ci)

          ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
          u_perp = npx%U_c( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + npx%V_c( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

          ! Calculate upwind momentum divergence
          ! =============================
          ! u_perp > 0: flow is exiting this vertex into vertex vj
          IF (u_perp > 0) THEN
            laddie%divQT( vi) = laddie%divQT( vi) + mesh%Cw( vi, ci) * u_perp * Hstar( vi) * npx%T( vi) / mesh%A( vi)
            laddie%divQS( vi) = laddie%divQS( vi) + mesh%Cw( vi, ci) * u_perp * Hstar( vi) * npx%S( vi) / mesh%A( vi)
          ! u_perp < 0: flow is entering this vertex from vertex vj
          ELSE
            IF (laddie%mask_oc_a( vj)) THEN
              CYCLE ! no inflow
              ! TODO fix boundary condition inflow
              laddie%divQT( vi) = laddie%divQT( vi) + mesh%Cw( vi, ci) * u_perp * Hstar( vi) * npx%T( vi) / mesh%A( vi)
              laddie%divQS( vi) = laddie%divQS( vi) + mesh%Cw( vi, ci) * u_perp * Hstar( vi) * npx%S( vi) / mesh%A( vi)
            ELSE
              laddie%divQT( vi) = laddie%divQT( vi) + mesh%Cw( vi, ci) * u_perp * Hstar( vj) * npx%T( vj) / mesh%A( vi)
              laddie%divQS( vi) = laddie%divQS( vi) + mesh%Cw( vi, ci) * u_perp * Hstar( vj) * npx%S( vj) / mesh%A( vi)
            END IF
          END IF

        END DO ! DO ci = 1, mesh%nC( vi)

      END IF ! (laddie%mask_a( vi))

    END DO ! DO vi = mesh%vi1, mesh%vi2
    call checksum( laddie%divQT, 'laddie%divQT', mesh%pai_V)
    call checksum( laddie%divQS, 'laddie%divQS', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_divQTS

END MODULE laddie_tracers

