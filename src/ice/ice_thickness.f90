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

  SUBROUTINE calc_Hi_tplusdt( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ! Calculate ice thickness at time t+dt

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_Hi_tplusdt'
    INTEGER                                               :: vi
    LOGICAL                                               :: found_negative_vals

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate Hi( t+dt) with the specified time discretisation scheme
    IF     (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_Hi_tplusdt_explicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ELSEIF (C%choice_ice_integration_method == 'implicit') THEN
      CALL calc_Hi_tplusdt_implicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      CALL calc_Hi_tplusdt_semiimplicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ELSE
      CALL crash('unknown choice_ice_integration_method "' // TRIM( C%choice_ice_integration_method) // '"!')
    END IF

    ! Limit Hi( t+dt) to zero; throw a warning if negative thickness are encountered
    found_negative_vals = .FALSE.
    DO vi = mesh%vi1, mesh%vi2
      IF (Hi_tplusdt( vi) < 0._dp) THEN
        found_negative_vals = .TRUE.
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

  END SUBROUTINE calc_Hi_tplusdt

  SUBROUTINE calc_Hi_tplusdt_explicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ! Calculate ice thickness rates of change (dH/dt)
    !
    ! Use a time-explicit discretisation scheme for the ice fluxes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_Hi_tplusdt_explicit'
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
    Hi_tplusdt = Hi + dHi_dt * dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Hi_tplusdt_explicit

  SUBROUTINE calc_Hi_tplusdt_implicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ! Calculate ice thickness rates of change (dH/dt)
    !
    ! Use a time-implicit discretisation scheme for the ice fluxes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_Hi_tplusdt_implicit'
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

    ! Add 1/dt to the diagonal
    DO vi = mesh%vi1, mesh%vi2
      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1)-1
      DO k = k1, k2
        vj = M_divQ%ind( k)
        IF (vj == vi) THEN
          AA%val( k) = AA%val( k) + 1._dp / dt
        END IF
      END DO ! DO k = k1, k2
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Load vector
    DO vi = mesh%vi1, mesh%vi2
      bb( vi) = Hi( vi) / dt + SMB( vi) + BMB( vi)
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Take current ice thickness as the initial guess
    Hi_tplusdt = Hi

    ! Solve for Hi_tplusdt
    CALL solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)
    CALL deallocate_matrix_CSR_dist( AA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Hi_tplusdt_implicit

  SUBROUTINE calc_Hi_tplusdt_semiimplicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ! Calculate ice thickness rates of change (dH/dt)
    !
    ! Use a time-implicit discretisation scheme for the ice fluxes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_vav_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_vav_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: BMB
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: dHi_dt
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT)   :: Hi_tplusdt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_Hi_tplusdt_semiimplicit'
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

    ! Multiply by f_s
    AA%val = AA%val * C%dHi_semiimplicit_fs

    ! Add 1/dt to the diagonal
    DO vi = mesh%vi1, mesh%vi2
      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1)-1
      DO k = k1, k2
        vj = M_divQ%ind( k)
        IF (vj == vi) THEN
          AA%val( k) = AA%val( k) + 1._dp / dt
        END IF
      END DO ! DO k = k1, k2
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Load vector
    CALL multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, M_divQ_H)
    DO vi = mesh%vi1, mesh%vi2
      bb( vi) = Hi( vi) / dt - ((1._dp - C%dHi_semiimplicit_fs) * M_divQ_H( vi)) + SMB( vi) + BMB( vi)
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Take current ice thickness as the initial guess
    Hi_tplusdt = Hi

    ! Solve for Hi_tplusdt
    CALL solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)
    CALL deallocate_matrix_CSR_dist( AA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Hi_tplusdt_semiimplicit

  SUBROUTINE calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, M_divQ)
    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme

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
