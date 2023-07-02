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
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D

  IMPLICIT NONE

CONTAINS

! == The main routines, to be called from the ice dynamics module

  SUBROUTINE calc_dHi_dt( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, mask_noice, dt, dHi_dt, Hi_tplusdt)
    ! Calculate ice thickness at time t+dt

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: mask_noice
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_dHi_dt'
    INTEGER                                               :: vi
    LOGICAL                                               :: found_negative_vals

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate Hi( t+dt) with the specified time discretisation scheme
    IF     (C%choice_ice_integration_method == 'none') THEN
      ! Unchanging ice geometry

      Hi_tplusdt = Hi
      dHi_dt     = 0._dp
      CALL finalise_routine( routine_name)
      RETURN

    ELSEIF (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_dHi_dt_explicit(     mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, mask_noice, dt, dHi_dt, Hi_tplusdt)
    ELSEIF (C%choice_ice_integration_method == 'implicit') THEN
      CALL calc_dHi_dt_implicit(     mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, mask_noice, dt, dHi_dt, Hi_tplusdt)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      CALL calc_dHi_dt_semiimplicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, mask_noice, dt, dHi_dt, Hi_tplusdt)
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

    ! Recalculate dH/dt with adjusted values of H
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt

  SUBROUTINE calc_dHi_dt_explicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, mask_noice, dt, dHi_dt, Hi_tplusdt)
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
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: mask_noice
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_dHi_dt_explicit'
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_divQ
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: divQ

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    CALL calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, divQ)

    ! Calculate rate of ice thickness change dHi/dt
    dHi_dt = -divQ + SMB + BMB

    ! Calculate ice thickness at t+dt
    Hi_tplusdt = MAX( 0._dp, Hi + dHi_dt * dt)

    ! Enfore Hi = 0 where told to do so
    CALL apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Recalculate dH/dt, accounting for limit of no negative ice thickness
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_explicit

  SUBROUTINE calc_dHi_dt_implicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, mask_noice, dt, dHi_dt, Hi_tplusdt)
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
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: mask_noice
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_dHi_dt_implicit'
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_divQ
    TYPE(type_sparse_matrix_CSR_dp)                       :: AA
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: bb
    INTEGER                                               :: vi, k1, k2, k, vj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    CALL calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)

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
      bb( vi) = Hi( vi) + MAX( -1._dp * Hi( vi), dt * (SMB( vi) + BMB( vi)))
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Take current ice thickness as the initial guess
    Hi_tplusdt = Hi

    ! Apply boundary conditions
    CALL apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt)

    ! Solve for Hi_tplusdt
    CALL solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol)

    ! Enfore Hi = 0 where told to do so
    CALL apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)
    CALL deallocate_matrix_CSR_dist( AA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_implicit

  SUBROUTINE calc_dHi_dt_semiimplicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, mask_noice, dt, dHi_dt, Hi_tplusdt)
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
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: mask_noice
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_dHi_dt_semiimplicit'
    TYPE(type_sparse_matrix_CSR_dp)                       :: M_divQ
    TYPE(type_sparse_matrix_CSR_dp)                       :: AA
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: M_divQ_H
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: bb
    INTEGER                                               :: vi, k1, k2, k, vj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    CALL calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)

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
      bb( vi) = Hi( vi) - (dt * (1._dp - C%dHi_semiimplicit_fs) * M_divQ_H( vi)) + MAX( -1._dp * Hi( vi), dt * (SMB( vi) + BMB( vi)))
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Take current ice thickness as the initial guess
    Hi_tplusdt = Hi

    ! Apply boundary conditions
    CALL apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt)

    ! Solve for Hi_tplusdt
    CALL solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol)

    ! Enfore Hi = 0 where told to do so
    CALL apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

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

      ! Loopover all connections of vertex vi
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

  SUBROUTINE apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt)
    ! Apply boundary conditions to the ice thickness matrix equation AA * Hi( t+dt) = bb

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: mask_noice
    TYPE(type_sparse_matrix_CSR_dp),        INTENT(INOUT) :: AA          ! Stiffness matrix
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: bb          ! Load vector
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi_tplusdt  ! Initial guess

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_ice_thickness_BC_matrix'
    INTEGER                                               :: vi,k1,k2,k,vj

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

  ! == No-ice mask
  ! ==============

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi)) THEN
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

      END IF ! IF (mask_noice( vi)) THEN
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_ice_thickness_BC_matrix

  SUBROUTINE apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)
    ! Enforce Hi = 0 where told to do so

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    LOGICAL,  DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: mask_noice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi_tplusdt  ! Initial guess

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_mask_noice_direct'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi)) Hi_tplusdt( vi) = 0._dp
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_mask_noice_direct

  SUBROUTINE calc_flux_limited_timestep( mesh, ice, dt_max)
    ! Calculate the largest time step that does not result in more
    ! ice flowing out of a cell than is contained within it.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_max

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_flux_limited_timestep'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL warning('Jorjo! Jorjo! Please fix me!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_flux_limited_timestep

END MODULE ice_thickness
