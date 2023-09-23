MODULE ice_thickness

  ! Contains all the routines needed to calculate ice thickness rates of change (dH/dt)

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: save_variable_as_netcdf_dp_1D, write_CSR_matrix_to_NetCDF
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist, duplicate_matrix_CSR_dist
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ice_velocity_main                                      , ONLY: map_velocities_from_b_to_c_2D
  USE petsc_basic                                            , ONLY: multiply_CSR_matrix_with_vector_1D, solve_matrix_equation_CSR_PETSc
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D
  USE math_utilities                                         , ONLY: ice_surface_elevation, Hi_from_Hb_Hs_and_SL
  USE math_utilities                                         , ONLY: is_floating
  USE netcdf_input                                           , ONLY: read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

! == The main routines, to be called from the ice dynamics module

  SUBROUTINE calc_dHi_dt( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, dHi_dt_residual, BC_prescr_mask, BC_prescr_Hi)
    ! Calculate ice thickness at time t+dt

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh                  ! [-]       The model mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hi                    ! [m]       Ice thickness at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hb                    ! [m]       Bedrock elevation at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SL                    ! [m]       Water surface elevation at time t
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SMB                   ! [m yr^-1] Surface mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: BMB                   ! [m yr^-1] Basal   mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: LMB                   ! [m yr^-1] Lateral mass balance
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    REAL(dp),                               INTENT(INOUT)           :: dt                    ! [dt]      Time step
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT)           :: dHi_dt_residual       ! [m yr^-1] Residual ice thickness rate of change
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'calc_dHi_dt'
    INTEGER                                                         :: vi
    LOGICAL                                                         :: found_negative_vals

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise tracker for residual (imposed) mass changes
    dHi_dt_residual = 0._dp

    ! Calculate Hi( t+dt) with the specified time discretisation scheme
    IF     (C%choice_ice_integration_method == 'none') THEN
      ! Unchanging ice geometry

      Hi_tplusdt = Hi
      dHi_dt     = 0._dp
      CALL finalise_routine( routine_name)
      RETURN

    ELSEIF (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_dHi_dt_explicit(     mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, dHi_dt_residual, BC_prescr_mask, BC_prescr_Hi)
    ELSEIF (C%choice_ice_integration_method == 'implicit') THEN
      CALL calc_dHi_dt_implicit(     mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, dHi_dt_residual, BC_prescr_mask, BC_prescr_Hi)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      CALL calc_dHi_dt_semiimplicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, dHi_dt_residual, BC_prescr_mask, BC_prescr_Hi)
    ELSE
      CALL crash('unknown choice_ice_integration_method "' // TRIM( C%choice_ice_integration_method) // '"!')
    END IF

    ! Limit Hi( t+dt) to zero; throw a warning if negative thickness are encountered
    found_negative_vals = .FALSE.
    DO vi = mesh%vi1, mesh%vi2
      IF (Hi_tplusdt( vi) < 0._dp) THEN
        ! Implicit solvers sometimes give VERY small negative numbers (e.g. -2e-189),
        ! only throw a warning if things get properly negative
        IF (Hi_tplusdt( vi) < -0.1_dp) found_negative_vals = .TRUE.
        ! Limit to zero
        Hi_tplusdt( vi) = 0._dp
      END IF
    END DO ! DO vi = mesh%vi1, mesh%vi2
    IF (found_negative_vals) THEN
      CALL warning('encountered negative values for Hi_tplusdt - time step too large?')
    END IF

    ! Add difference between original and applied dHi_dt to residual tracker
    dHi_dt_residual = dHi_dt_residual + (dHi_dt - (Hi_tplusdt - Hi) / dt)

    ! Recalculate dH/dt with adjusted values of H
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt

  SUBROUTINE calc_dHi_dt_explicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, dHi_dt_residual, BC_prescr_mask, BC_prescr_Hi)
    ! Calculate ice thickness rates of change (dH/dt)
    !
    ! Use a time-explicit discretisation scheme for the ice fluxes
    !
    ! The ice continuity equation (alternatively known as the ice thickness equation,
    ! or just conservation of mass) reads:
    !
    !     [1] dH/dt = -div( Q) + m
    !
    ! Here, Q is the horizontal ice flux vector, and m is the net mass balance.
    !
    ! We define a matrix operator M_divQ that can be multiplied with the ice thickness
    ! to produce the flux divergence:
    !
    !     [2] div( Q) = M_divQ * H
    !
    ! Substituting [2] into [1] yields:
    !
    !     [3] dH/dt = -M_divQ H + m
    !
    ! Using a time-explicit discretisation scheme so that H on the right-hand side
    ! is defined at time t yields:
    !
    !     [4] (H( t+dt) - H( t)) / dt = -M_divQ H( t) + m

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh                  ! [-]       The model mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hi                    ! [m]       Ice thickness at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hb                    ! [m]       Bedrock elevation at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SL                    ! [m]       Water surface elevation at time t
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SMB                   ! [m yr^-1] Surface mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: BMB                   ! [m yr^-1] Basal   mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: LMB                   ! [m yr^-1] Lateral mass balance
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    REAL(dp),                               INTENT(INOUT)           :: dt                    ! [dt]      Time step
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT)           :: dHi_dt_residual       ! [m yr^-1] Residual ice thickness rate of change
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'calc_dHi_dt_explicit'
    TYPE(type_sparse_matrix_CSR_dp)                                 :: M_divQ
    REAL(dp)                                                        :: dt_max
    INTEGER                                                         :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    CALL calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, divQ)

    ! Calculate rate of ice thickness change dHi/dt
    dHi_dt = -divQ + SMB + BMB + LMB - dHi_dt_target

    ! Calculate largest time step possible based on dHi_dt
    CALL calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt, dt_max)

    ! Constrain dt based on new limit
    dt = MIN( dt, dt_max)

    ! Calculate ice thickness at t+dt
    Hi_tplusdt = MAX( 0._dp, Hi + dHi_dt * dt)

    ! Apply boundary conditions at the domain border
    CALL apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)

    ! Set predicted ice thickness to prescribed values where told to do so
    IF (PRESENT( BC_prescr_mask) .OR. PRESENT( BC_prescr_Hi)) THEN
      ! Safety
      IF (.NOT. (PRESENT( BC_prescr_mask) .AND. PRESENT( BC_prescr_Hi))) THEN
        CALL crash('need to provide prescribed both Hi and mask!')
      END IF
      DO vi = mesh%vi1, mesh%vi2
        IF (BC_prescr_mask( vi) == 1) THEN
          Hi_tplusdt( vi) = MAX( 0._dp, BC_prescr_Hi( vi))
        END IF
      END DO
    END IF

    ! Enforce Hi = 0 where told to do so
    CALL apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Add difference between original and applied dHi_dt to residual tracker
    dHi_dt_residual = dHi_dt_residual + (dHi_dt - (Hi_tplusdt - Hi) / dt)

    ! Recalculate dH/dt, accounting for limit of no negative ice thickness
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_explicit

  SUBROUTINE calc_dHi_dt_implicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, dHi_dt_residual, BC_prescr_mask, BC_prescr_Hi)
    ! Calculate ice thickness rates of change (dH/dt)
    !
    ! Use a time-implicit discretisation scheme for the ice fluxes
    !
    ! The ice continuity equation (alternatively known as the ice thickness equation,
    ! or just conservation of mass) reads:
    !
    !     [1] dH/dt = -div( Q) + m
    !
    ! Here, Q is the horizontal ice flux vector, and m is the net mass balance.
    !
    ! We define a matrix operator M_divQ that can be multiplied with the ice thickness
    ! to produce the flux divergence:
    !
    !     [2] div( Q) = M_divQ * H
    !
    ! Substituting [2] into [1] yields:
    !
    !     [3] dH/dt = -M_divQ H + m
    !
    ! Using a time-implicit discretisation scheme so that H on the right-hand side
    ! is defined at time t+dt yields:
    !
    !     [4] (H( t+dt) - H( t)) / dt = -M_divQ H( t+dt) + m
    !
    ! Rearranging to place all H( t+dt) terms on the left-hand side yields:
    !
    !     [5] (1/dt + M_divQ) H( t+dt) = H( t) / dt + m
    !
    ! Finally, multiplying both sides by dt yields:
    !
    !     [6] ( 1 + dt M_divQ) H( t+dt) = H( t) + m dt
    !
    ! This is a matrix equation, with the stiffness matrix A = 1 + dt M_divQ, the
    ! load vector b = H( t) + dt m, which can be solved for H( t+dt)
    !
    ! Multiplying by dt gives the advantage that the mass balance term on the right
    ! has units of meters, so we can more easily calculate the "applied" mass balance
    ! to limit melt to the available ice mass.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh                  ! [-]       The model mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hi                    ! [m]       Ice thickness at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hb                    ! [m]       Bedrock elevation at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SL                    ! [m]       Water surface elevation at time t
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SMB                   ! [m yr^-1] Surface mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: BMB                   ! [m yr^-1] Basal   mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: LMB                   ! [m yr^-1] Lateral mass balance
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    REAL(dp),                               INTENT(INOUT)           :: dt                    ! [dt]      Time step
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT)           :: dHi_dt_residual       ! [m yr^-1] Residual ice thickness rate of change
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'calc_dHi_dt_implicit'
    TYPE(type_sparse_matrix_CSR_dp)                                 :: M_divQ
    TYPE(type_sparse_matrix_CSR_dp)                                 :: AA
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                          :: bb
    INTEGER                                                         :: vi, k1, k2, k, vj
    REAL(dp)                                                        :: dt_max
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                          :: dHi_dt_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    CALL calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, divQ)

    ! Calculate an estimate of the rate of ice thickness change dHi/dt
    dHi_dt_dummy = -divQ + SMB + BMB + LMB - dHi_dt_target

    ! Calculate largest time step possible based on that estimate
    CALL calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt_dummy, dt_max)

    ! Constrain dt based on new limit
    dt = MIN( dt, dt_max)

    ! Calculate the stiffness matrix A and the load vector b

    ! Start by letting AA = M_divQ
    CALL duplicate_matrix_CSR_dist( M_divQ, AA)

    ! Multiply by dt so AA = dt M_divQ
    AA%val = AA%val * dt

    ! Add 1 to the diagonal
    DO vi = mesh%vi1, mesh%vi2
      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1)-1
      DO k = k1, k2
        vj = M_divQ%ind( k)
        IF (vj == vi) THEN
          AA%val( k) = AA%val( k) + 1._dp
        END IF
      END DO ! DO k = k1, k2
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Load vector
    DO vi = mesh%vi1, mesh%vi2
      bb( vi) = Hi( vi) + MAX( -1._dp * Hi( vi), dt * (SMB( vi) + BMB( vi) + LMB( vi) - dHi_dt_target( vi)))
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Take the current ice thickness plus the current thinning rate as the initial guess
    Hi_tplusdt = Hi + dt * dHi_dt

    ! Apply boundary conditions
    CALL apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt, BC_prescr_mask, BC_prescr_Hi)

    ! Solve for Hi_tplusdt
    CALL solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol)

    ! Save computed dHi_dt before enforced mass changes
    dHi_dt_dummy = (Hi_tplusdt - Hi) / dt

    ! Enforce Hi = 0 where told to do so
    CALL apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Apply boundary conditions at the domain border
    CALL apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Add difference between original and applied dHi_dt to residual tracker
    dHi_dt_residual = dHi_dt_residual + dHi_dt_dummy - dHi_dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)
    CALL deallocate_matrix_CSR_dist( AA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_implicit

  SUBROUTINE calc_dHi_dt_semiimplicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, dHi_dt_residual, BC_prescr_mask, BC_prescr_Hi)
    ! Calculate ice thickness rates of change (dH/dt)
    !
    ! Use a semi-implicit time discretisation scheme for the ice fluxes
    !
    ! The ice continuity equation (alternatively known as the ice thickness equation,
    ! or just conservation of mass) reads:
    !
    !     [1] dH/dt = -div( Q) + m
    !
    ! Here, Q is the horizontal ice flux vector, and m is the net mass balance.
    !
    ! We define a matrix operator M_divQ that can be multiplied with the ice thickness
    ! to produce the flux divergence:
    !
    !     [2] div( Q) = M_divQ * H
    !
    ! Substituting [2] into [1] yields:
    !
    !     [3] dH/dt = -M_divQ H + m
    !
    ! Using a semi-implicit discretisation scheme, so that H on the right-hand side
    ! is defined as a weighted average of H( t) and H( t+dt), yields:
    !
    !     [4] (H( t+dt) - H( t)) / dt = -M_divQ [f_s H( t+dt) + (1 - f_s) H( t)] + m
    !
    ! This implies that f_s = 0 is equivalent to the explicit scheme, f_s = 1 is the
    ! implicit scheme, 0 < f_s < 1 is called "semi-implicit", and f_s > 1 is called
    ! "over-implicit".
    !
    ! Rearranging to place all H( t+dt) terms on the left-hand side yields:
    !
    !         (1/dt + f_s M_divQ) H( t+dt) = (1 / dt - (1 - f_s) M_divQ) H( t) + m
    !     [5] (1/dt + f_s M_divQ) H( t+dt) = H( t) / dt - (1 - f_s) M_divQ H( t) + m
    !
    !
    ! Finally, multiplying both sides by dt yields:
    !
    !     [6] ( 1 + dt f_s M_divQ) H( t+dt) = H( t) - dt (1 - f_s) M_divQ H( t) + m dt
    !
    ! This is a matrix equation, with the stiffness matrix A = 1/dt + f_s M_divQ, the
    ! load vector b = H( t) / dt - (1 - f_s) M_divQ H( t) + m, which can be solved for H( t+dt).
    !
    ! Multiplying by dt gives the advantage that the mass balance term on the right
    ! has units of meters, so we can more easily calculate the "applied" mass balance
    ! to limit melt to the available ice mass.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh                  ! [-]       The model mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hi                    ! [m]       Ice thickness at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hb                    ! [m]       Bedrock elevation at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SL                    ! [m]       Water surface elevation at time t
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)              :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SMB                   ! [m yr^-1] Surface mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: BMB                   ! [m yr^-1] Basal   mass balance
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: LMB                   ! [m yr^-1] Lateral mass balance
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    REAL(dp),                               INTENT(INOUT)           :: dt                    ! [dt]      Time step
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)             :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT)           :: dHi_dt_residual       ! [m yr^-1] Residual ice thickness rate of change
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'calc_dHi_dt_semiimplicit'
    TYPE(type_sparse_matrix_CSR_dp)                                 :: M_divQ
    TYPE(type_sparse_matrix_CSR_dp)                                 :: AA
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                          :: M_divQ_H
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                          :: bb
    INTEGER                                                         :: vi, k1, k2, k, vj
    REAL(dp)                                                        :: dt_max
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                          :: dHi_dt_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    CALL calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, divQ)

    ! Calculate an estimate of the rate of ice thickness change dHi/dt
    dHi_dt_dummy = -divQ + SMB + BMB + LMB - dHi_dt_target

    ! Calculate largest time step possible based on that estimate
    CALL calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt, dt_max)

    ! Constrain dt based on new limit
    dt = MIN( dt, dt_max)

    ! Calculate the stiffness matrix A and the load vector b

    ! Start by letting AA = M_divQ
    CALL duplicate_matrix_CSR_dist( M_divQ, AA)

    ! Multiply by dt f_s
    AA%val = AA%val * dt * C%dHi_semiimplicit_fs

    ! Add 1 to the diagonal
    DO vi = mesh%vi1, mesh%vi2
      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1)-1
      DO k = k1, k2
        vj = M_divQ%ind( k)
        IF (vj == vi) THEN
          AA%val( k) = AA%val( k) + 1._dp
        END IF
      END DO ! DO k = k1, k2
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Load vector
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, M_divQ_H)
    DO vi = mesh%vi1, mesh%vi2
      bb( vi) = Hi( vi) - (dt * (1._dp - C%dHi_semiimplicit_fs) * M_divQ_H( vi)) + MAX( -1._dp * Hi( vi), dt * (SMB( vi) + BMB( vi) + LMB( vi) - dHi_dt_target( vi)))
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Take the current ice thickness plus the current thinning rate as the initial guess
    Hi_tplusdt = Hi + dt * dHi_dt

    ! Apply boundary conditions
    CALL apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt, BC_prescr_mask, BC_prescr_Hi)

    ! Solve for Hi_tplusdt
    CALL solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol)

    ! Save computed dHi_dt before enforced mass changes
    dHi_dt_dummy = (Hi_tplusdt - Hi) / dt

    ! Enforce Hi = 0 where told to do so
    CALL apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Apply boundary conditions at the domain border
    CALL apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Add difference between original and applied dHi_dt to residual tracker
    dHi_dt_residual = dHi_dt_residual + dHi_dt_dummy - dHi_dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)
    CALL deallocate_matrix_CSR_dist( AA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_semiimplicit

  SUBROUTINE calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)
    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
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
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti1), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti1), INTENT(IN)    :: v_vav_b
    TYPE(type_sparse_matrix_CSR_dp),        INTENT(OUT)   :: M_divQ

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_ice_flux_divergence_matrix_upwind'
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2)                :: u_vav_c, v_vav_c
    REAL(dp), DIMENSION(mesh%nE)                          :: u_vav_c_tot, v_vav_c_tot
    INTEGER                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                               :: vi, ci, ei, vj
    REAL(dp)                                              :: A_i, L_c
    REAL(dp)                                              :: D_x, D_y, D, u_perp
    REAL(dp), DIMENSION(0:mesh%nC_mem)                    :: cM_divQ

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL map_velocities_from_b_to_c_2D( mesh, u_vav_b, v_vav_b, u_vav_c, v_vav_c)
    CALL gather_to_all_dp_1D( u_vav_c, u_vav_c_tot)
    CALL gather_to_all_dp_1D( v_vav_c, v_vav_c_tot)

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV      ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_est_proc    = mesh%nV_loc + SUM( mesh%nC( mesh%vi1:mesh%vi2))

    CALL allocate_matrix_CSR_dist( M_divQ, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! == Calculate coefficients
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise
      cM_divQ = 0._dp

      ! Loop over all connections of vertex vi
      DO ci = 1, mesh%nC( vi)

        ! Connection ci from vertex vi leads through edge ei to vertex vj
        ei = mesh%VE( vi,ci)
        vj = mesh%C(  vi,ci)

        ! The Voronoi cell of vertex vi has area A_i
        A_i = mesh%A( vi)

        ! The shared Voronoi cell boundary section between the Voronoi cells
        ! of vertices vi and vj has length L_c
        L_c = mesh%Cw( vi,ci)

        ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
        D_x = mesh%V( vj,1) - mesh%V( vi,1)
        D_y = mesh%V( vj,2) - mesh%V( vi,2)
        D   = SQRT( D_x**2 + D_y**2)
        u_perp = u_vav_c_tot( ei) * D_x/D + v_vav_c_tot( ei) * D_y/D

        ! Calculate matrix coefficients
        cM_divQ( 0 ) = cM_divQ( 0) + L_c * MAX( 0._dp, u_perp) / A_i
        cM_divQ( ci) =               L_c * MIN( 0._dp, u_perp) / A_i

      END DO ! DO ci = 1, mesh%nC( vi)

      ! Add coefficients to matrix
      CALL add_entry_CSR_dist( M_divQ, vi, vi, cM_divQ( 0))
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C(  vi,ci)
        CALL add_entry_CSR_dist( M_divQ, vi, vj, cM_divQ( ci))
      END DO ! DO ci = 1, mesh%nC( vi)

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_ice_flux_divergence_matrix_upwind

  SUBROUTINE apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)
    ! Apply boundary conditions to the ice thickness on the domain border directly

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh                  ! [-]       The model mesh
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: Hb                    ! [m]       Bedrock elevation at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: SL                    ! [m]       Water surface elevation at time t
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT)           :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'apply_ice_thickness_BC_explicit'
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                          :: Hs_tplusdt
    INTEGER                                                         :: vi
    REAL(dp), DIMENSION(mesh%nV)                                    :: Hs_tplusdt_tot
    LOGICAL,  DIMENSION(mesh%nV)                                    :: mask_noice_tot
    CHARACTER(LEN=256)                                              :: BC_H
    INTEGER                                                         :: ci,vj,n_interior_neighbours
    REAL(dp)                                                        :: Hs_sum

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate Hs( t+dt)
    DO vi = mesh%vi1, mesh%vi2
      Hs_tplusdt( vi) = ice_surface_elevation( Hi_tplusdt( vi), Hb( vi), SL( vi))
    END DO

    ! Gather global data fields
    CALL gather_to_all_dp_1D(      Hs_tplusdt, Hs_tplusdt_tot)
    CALL gather_to_all_logical_1D( mask_noice, mask_noice_tot)

    ! == First pass: set values of border vertices to mean of interior neighbours
    !    ...for those border vertices that actually have interior neighbours.
    ! ===========================================================================

    DO vi = mesh%vi1, mesh%vi2

      IF     (mesh%VBI( vi) == 1 .OR. mesh%VBI( vi) == 2) THEN
        ! Northern domain border
        BC_H = C%BC_H_north
      ELSEIF (mesh%VBI( vi) == 3 .OR. mesh%VBI( vi) == 4) THEN
        ! Eastern domain border
        BC_H = C%BC_H_east
      ELSEIF (mesh%VBI( vi) == 5 .OR. mesh%VBI( vi) == 6) THEN
        ! Southern domain border
        BC_H = C%BC_H_south
      ELSEIF (mesh%VBI( vi) == 7 .OR. mesh%VBI( vi) == 8) THEN
        ! Western domain border
        BC_H = C%BC_H_west
      ELSE
        ! Free vertex
        CYCLE
      END IF

      SELECT CASE (BC_H)
        CASE ('zero')
          ! Set ice thickness to zero here

          Hi_tplusdt( vi) = 0._dp

        CASE ('infinite')
          ! Set H on this vertex equal to the average value on its neighbours

          n_interior_neighbours = 0
          Hs_sum = 0._dp

          DO ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            IF (mesh%VBI( vj) == 0 .AND. .NOT. mask_noice_tot( vj)) THEN
              n_interior_neighbours = n_interior_neighbours + 1
              Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
            END IF
          END DO ! DO ci = 1, mesh%nC( vi)

          IF (n_interior_neighbours > 0) THEN
            Hs_tplusdt( vi) = MAX( Hb( vi), Hs_sum / REAL( n_interior_neighbours,dp) )
            Hi_tplusdt( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))
          END IF

        CASE DEFAULT
          CALL crash('unknown BC_H "' // TRIM( BC_H) // '"')
      END SELECT

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! == Second pass: set values of border vertices to mean of all neighbours
    !    ...for those border vertices that have no interior neighbours.
    ! =======================================================================

    ! Gather global data fields again
    CALL gather_to_all_dp_1D( Hs_tplusdt, Hs_tplusdt_tot)

    DO vi = mesh%vi1, mesh%vi2

      IF     (mesh%VBI( vi) == 1 .OR. mesh%VBI( vi) == 2) THEN
        ! Northern domain border
        BC_H = C%BC_H_north
      ELSEIF (mesh%VBI( vi) == 3 .OR. mesh%VBI( vi) == 4) THEN
        ! Eastern domain border
        BC_H = C%BC_H_east
      ELSEIF (mesh%VBI( vi) == 5 .OR. mesh%VBI( vi) == 6) THEN
        ! Southern domain border
        BC_H = C%BC_H_south
      ELSEIF (mesh%VBI( vi) == 7 .OR. mesh%VBI( vi) == 8) THEN
        ! Western domain border
        BC_H = C%BC_H_west
      ELSE
        ! Free vertex
        CYCLE
      END IF

      SELECT CASE (BC_H)
        CASE ('zero')
          ! Set ice thickness to zero here

          Hi_tplusdt( vi) = 0._dp

        CASE ('infinite')
          ! Set H on this vertex equal to the average value on its neighbours

          n_interior_neighbours = 0
          Hs_sum = 0._dp

          DO ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
            IF (mesh%VBI( vj) == 0 .AND. .NOT. mask_noice_tot( vj)) THEN
              n_interior_neighbours = n_interior_neighbours + 1
            END IF
          END DO ! DO ci = 1, mesh%nC( vi)

          IF (n_interior_neighbours == 0) THEN
            Hs_tplusdt( vi) = MAX( Hb( vi), Hs_sum / REAL( mesh%nC( vi),dp) )
            Hi_tplusdt( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))
          END IF

        CASE DEFAULT
          CALL crash('unknown BC_H "' // TRIM( BC_H) // '"')
      END SELECT

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_ice_thickness_BC_explicit

  SUBROUTINE apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt, BC_prescr_mask, BC_prescr_Hi)
    ! Apply boundary conditions to the ice thickness matrix equation AA * Hi( t+dt) = bb

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)              :: mask_noice            ! Mask of vertices where no ice is allowed
    TYPE(type_sparse_matrix_CSR_dp),        INTENT(INOUT)           :: AA                    ! Stiffness matrix
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT)           :: bb                    ! Load vector
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT)           :: Hi_tplusdt            ! Initial guess
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_mask        ! Mask of vertices where thickness is prescribed
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)   , OPTIONAL :: BC_prescr_Hi          ! Prescribed thicknesses

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'apply_ice_thickness_BC_matrix'
    INTEGER                                                         :: vi,k1,k2,k,vj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Boundary conditions at the domain border
    ! ===========================================

    DO vi = mesh%vi1, mesh%vi2

      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1) - 1

      IF     (mesh%VBI( vi) == 1 .OR. mesh%VBI( vi) == 2) THEN
        ! Northern domain border

        SELECT CASE (C%BC_H_north)
          CASE ('zero')
            ! Set ice thickness to zero here

            ! Set diagonal element of A to 1, rest of row to 0
            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = 1._dp
              ELSE
                ! Off-diagonal element
                AA%val( k) = 0._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector and initial guess
            bb( vi) = 0._dp
            Hi_tplusdt( vi) = 0._dp

          CASE ('infinite')
            ! Set H on this vertex equal to the average value on its neighbours

            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = REAL( mesh%nC( vi), dp)
              ELSE
                ! Off-diagonal element
                AA%val( k) = -1._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector
            bb( vi) = 0._dp

          CASE DEFAULT
            CALL crash('unknown BC_H_north "' // TRIM( C%BC_H_north) // '"')
        END SELECT

      ELSEIF (mesh%VBI( vi) == 3 .OR. mesh%VBI( vi) == 4) THEN
        ! Eastern domain border

        SELECT CASE (C%BC_H_east)
          CASE ('zero')
            ! Set ice thickness to zero here

            ! Set diagonal element of A to 1, rest of row to 0
            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = 1._dp
              ELSE
                ! Off-diagonal element
                AA%val( k) = 0._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector and initial guess
            bb( vi) = 0._dp
            Hi_tplusdt( vi) = 0._dp

          CASE ('infinite')
            ! Set H on this vertex equal to the average value on its neighbours

            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = REAL( mesh%nC( vi), dp)
              ELSE
                ! Off-diagonal element
                AA%val( k) = -1._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector
            bb( vi) = 0._dp

          CASE DEFAULT
            CALL crash('unknown BC_H_east "' // TRIM( C%BC_H_east) // '"')
        END SELECT

      ELSEIF (mesh%VBI( vi) == 5 .OR. mesh%VBI( vi) == 6) THEN
        ! Southern domain border

        SELECT CASE (C%BC_H_south)
          CASE ('zero')
            ! Set ice thickness to zero here

            ! Set diagonal element of A to 1, rest of row to 0
            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = 1._dp
              ELSE
                ! Off-diagonal element
                AA%val( k) = 0._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector and initial guess
            bb( vi) = 0._dp
            Hi_tplusdt( vi) = 0._dp

          CASE ('infinite')
            ! Set H on this vertex equal to the average value on its neighbours

            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = REAL( mesh%nC( vi), dp)
              ELSE
                ! Off-diagonal element
                AA%val( k) = -1._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector
            bb( vi) = 0._dp

          CASE DEFAULT
            CALL crash('BC_H_south "' // TRIM( C%BC_H_south) // '"')
        END SELECT

      ELSEIF (mesh%VBI( vi) == 7 .OR. mesh%VBI( vi) == 8) THEN
        ! Western domain border

        SELECT CASE (C%BC_H_west)
          CASE ('zero')
            ! Set ice thickness to zero here

            ! Set diagonal element of A to 1, rest of row to 0
            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = 1._dp
              ELSE
                ! Off-diagonal element
                AA%val( k) = 0._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector and initial guess
            bb( vi) = 0._dp
            Hi_tplusdt( vi) = 0._dp

          CASE ('infinite')
            ! Set H on this vertex equal to the average value on its neighbours

            DO k = k1, k2
              vj = AA%ind( k)
              IF (vj == vi) THEN
                ! Diagonal element
                AA%val( k) = REAL( mesh%nC( vi), dp)
              ELSE
                ! Off-diagonal element
                AA%val( k) = -1._dp
              END IF
            END DO ! DO k = k1, k2

            ! Load vector
            bb( vi) = 0._dp

          CASE DEFAULT
            CALL crash('BC_H_west "' // TRIM( C%BC_H_west) // '"')
        END SELECT

      ELSE
        ! Free vertex
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Set predicted ice thickness to prescribed values where told to do so
    ! ====================================================================

    IF (PRESENT( BC_prescr_mask) .OR. PRESENT( BC_prescr_Hi)) THEN
      ! Safety
      IF (.NOT. (PRESENT( BC_prescr_mask) .AND. PRESENT( BC_prescr_Hi))) THEN
        CALL crash('need to provide prescribed both Hi and mask!')
      END IF

      DO vi = mesh%vi1, mesh%vi2
        IF (BC_prescr_mask( vi) == 1) THEN

          ! Set diagonal element of A to 1, rest of row to 0
          k1 = AA%ptr( vi)
          k2 = AA%ptr( vi+1) - 1
          DO k = k1, k2
            vj = AA%ind( k)
            IF (vj == vi) THEN
              ! Diagonal element
              AA%val( k) = 1._dp
            ELSE
              ! Off-diagonal element
              AA%val( k) = 0._dp
            END IF
          END DO ! DO k = k1, k2

          ! Load vector and initial guess
          bb        ( vi) = BC_prescr_Hi( vi)
          Hi_tplusdt( vi) = BC_prescr_Hi( vi)

        END IF ! IF (BC_prescr_mask( vi) == 1) THEN
      END DO ! DO vi = mesh%vi1, mesh%vi2

    END IF ! IF (PRESENT( BC_prescr_mask) .OR. PRESENT( BC_prescr_Hi)) THEN

    ! == No-ice mask
    ! ==============

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi)) THEN
        ! Set ice thickness to zero here

        ! Set diagonal element of A to 1, rest of row to 0
        k1 = AA%ptr( vi)
        k2 = AA%ptr( vi+1) - 1
        DO k = k1, k2
          vj = AA%ind( k)
          IF (vj == vi) THEN
            ! Diagonal element
            AA%val( k) = 1._dp
          ELSE
            ! Off-diagonal element
            AA%val( k) = 0._dp
          END IF
        END DO ! DO k = k1, k2

        ! Load vector and initial guess
        bb( vi) = 0._dp
        Hi_tplusdt( vi) = 0._dp

      END IF ! IF (mask_noice( vi)) THEN
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_ice_thickness_BC_matrix

  SUBROUTINE apply_mask_noice_direct( mesh, mask_noice, Hi)
    ! Enforce Hi = 0 where told to do so

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: mask_noice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_mask_noice_direct'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi)) Hi( vi) = 0._dp
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_mask_noice_direct

  SUBROUTINE calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt, dt_max)
    ! Calculate the largest time step that does not result in more
    ! ice flowing out of a cell than is contained within it.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)  :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)  :: Hi
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)  :: Hb
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)  :: SL
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)  :: dHi_dt
    REAL(dp),                               INTENT(OUT) :: dt_max

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'calc_flux_limited_timestep'
    INTEGER                                             :: vi
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)              :: dt_lim

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    dt_lim = C%dt_ice_max

    ! Loop over each mesh vertex within this process
    DO vi = mesh%vi1, mesh%vi2
      ! If there is [non-negligible grounded] ice, and there is mass loss
      IF (.NOT. is_floating( Hi( vi), Hb( vi), SL( vi)) .AND. Hi( vi) >= 10._dp .AND. dHi_dt( vi) < 0._dp) THEN

        ! Compute time step limit (in yr) based on
        ! available ice thickness and flux divergence
        dt_lim( vi) = Hi( vi) / MAX( dHi_dt( vi), 1E-9_dp)

      END IF
    END DO

    ! Get most strict time step limit for this process
    dt_max = MINVAL( dt_lim)

    ! Get most strict time step limit among all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_max, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Limit to minimum ice model time step
    dt_max = MAX( C%dt_ice_min, dt_max)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_flux_limited_timestep

! == Target dHi_dt initialisation

  SUBROUTINE initialise_dHi_dt_target( mesh, ice, region_name)
    ! Prescribe a target dHi_dt from a file without a time dimension

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_dHi_dt_target'
    CHARACTER(LEN=256)                                    :: filename_dHi_dt_prescribed
    REAL(dp)                                              :: timeframe_dHi_dt_prescribed

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE ('NAM')
        filename_dHi_dt_prescribed  = C%filename_dHi_dt_prescribed_NAM
        timeframe_dHi_dt_prescribed = C%timeframe_dHi_dt_prescribed_NAM
      CASE ('EAS')
        filename_dHi_dt_prescribed  = C%filename_dHi_dt_prescribed_EAS
        timeframe_dHi_dt_prescribed = C%timeframe_dHi_dt_prescribed_EAS
      CASE ('GRL')
        filename_dHi_dt_prescribed  = C%filename_dHi_dt_prescribed_GRL
        timeframe_dHi_dt_prescribed = C%timeframe_dHi_dt_prescribed_GRL
      CASE ('ANT')
        filename_dHi_dt_prescribed  = C%filename_dHi_dt_prescribed_ANT
        timeframe_dHi_dt_prescribed = C%timeframe_dHi_dt_prescribed_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising target dHi_dt from file "' // colour_string( TRIM( filename_dHi_dt_prescribed),'light blue') // '"...'

    ! Read dHi_dt from file
    IF (timeframe_dHi_dt_prescribed == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_2D( filename_dHi_dt_prescribed, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target)
    ELSE
      ! Assume the file has a time dimension, and read the specified timeframe
      CALL read_field_from_file_2D( filename_dHi_dt_prescribed, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target, time_to_read = timeframe_dHi_dt_prescribed)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_dHi_dt_target

END MODULE ice_thickness
