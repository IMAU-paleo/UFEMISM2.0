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
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ice_velocity_main                                      , ONLY: map_velocities_from_b_to_c_2D
  USE petsc_basic                                            , ONLY: multiply_CSR_matrix_with_vector_1D
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

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_Hi_tplusdt_explicit( mesh, Hi, u_vav_b, v_vav_b, SMB, BMB, dt, dHi_dt, Hi_tplusdt)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      ! DENK DROM
      CALL crash('fixme!')
    ELSE
      CALL crash('unknown choice_ice_integration_method "' // TRIM( C%choice_ice_integration_method) // '"!')
    END IF

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
    INTEGER                                               :: vi
    LOGICAL                                               :: found_negative_vals

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

    ! Limit ice thickness to zero; throw a warning if negative thickness are encountered
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

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( M_divQ)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Hi_tplusdt_explicit

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
    REAL(dp)                                              :: cM_divQ

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
    nnz_est_proc    = SUM( mesh%nC( mesh%vi1:mesh%vi2))

    CALL allocate_matrix_CSR_dist( M_divQ, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

  ! == Calculate coefficients
  ! =========================

    DO vi = mesh%vi1, mesh%vi2
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

      ! Value of matrix coefficient
      cM_divQ = L_c * u_perp / A_i

      ! Upwind scheme
      IF (u_perp > 0._dp) THEN
        ! Ice flows from vi to vj
        CALL add_entry_CSR_dist( M_divQ, vi, vi, cM_divQ)
      ELSE
        ! Ice flows from vj to vi
        CALL add_entry_CSR_dist( M_divQ, vi, vj, cM_divQ)
      END IF

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
