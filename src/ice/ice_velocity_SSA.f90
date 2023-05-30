MODULE ice_velocity_SSA

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA)

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: save_variable_as_netcdf_dp_1D, save_variable_as_netcdf_dp_2D
  USE petsc_basic                                            , ONLY: solve_matrix_equation_CSR_PETSc
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_SSA
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate_clean
  USE mesh_operators                                         , ONLY: map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D
  USE mesh_zeta                                              , ONLY: vertical_average
  USE sliding_laws                                           , ONLY: calc_basal_friction_coefficient
  USE mesh_utilities                                         , ONLY: find_ti_copy_ISMIP_HOM_periodic
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist
  USE analytical_solutions                                   , ONLY: Schoof2006_icestream

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE initialise_SSA_solver( mesh, SSA)
    ! Initialise the SSA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(OUT)   :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SSA_solver'
    INTEGER :: tii, ti
    REAL(dp) :: y, tau_c

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    CALL allocate_SSA_solver( mesh, SSA)

    ! Set tolerances for PETSc matrix solver for the linearised SSA
    SSA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    SSA%PETSc_abstol = C%stress_balance_PETSc_abstol

!    ! DENK DROM
!    DO tii = 1, mesh%nTri_loc
!      ti = tii + mesh%ti1 - 1
!      y = mesh%TriGC( ti,2)
!      CALL Schoof2006_icestream( C%uniform_flow_factor, C%n_flow, C%refgeo_idealised_SSA_icestream_Hi, &
!        C%refgeo_idealised_SSA_icestream_dhdx, C%refgeo_idealised_SSA_icestream_L, C%refgeo_idealised_SSA_icestream_m, &
!        y, SSA%u_b( tii), tau_c)
!    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SSA_solver

  SUBROUTINE allocate_SSA_solver( mesh, SSA)
    ! Allocate memory the SSA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(OUT)   :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_SSA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Solution
    ALLOCATE( SSA%u_b(          mesh%nTri_loc), source = 0._dp)                   ! [m yr^-1] 2-D horizontal ice velocity
    ALLOCATE( SSA%v_b(          mesh%nTri_loc), source = 0._dp)

    ! Intermediate data fields
    ALLOCATE( SSA%A_flow_vav_a( mesh%nV_loc  ), source = 0._dp)                   ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter
    ALLOCATE( SSA%du_dx_a(      mesh%nV_loc  ), source = 0._dp)                   ! [yr^-1] 2-D horizontal strain rates
    ALLOCATE( SSA%du_dy_a(      mesh%nV_loc  ), source = 0._dp)
    ALLOCATE( SSA%dv_dx_a(      mesh%nV_loc  ), source = 0._dp)
    ALLOCATE( SSA%dv_dy_a(      mesh%nV_loc  ), source = 0._dp)
    ALLOCATE( SSA%eta_a(        mesh%nV_loc  ), source = 0._dp)                   ! Effective viscosity
    ALLOCATE( SSA%N_a(          mesh%nV_loc  ), source = 0._dp)                   ! Product term N = eta * H
    ALLOCATE( SSA%N_b(          mesh%nTri_loc), source = 0._dp)
    ALLOCATE( SSA%dN_dx_b(      mesh%nTri_loc), source = 0._dp)                   ! Gradients of N
    ALLOCATE( SSA%dN_dy_b(      mesh%nTri_loc), source = 0._dp)
    ALLOCATE( SSA%beta_b_b(     mesh%nTri_loc), source = 0._dp)                   ! Friction coefficient (tau_b = u * beta_b)
    ALLOCATE( SSA%tau_dx_b(     mesh%nTri_loc), source = 0._dp)                   ! Driving stress
    ALLOCATE( SSA%tau_dy_b(     mesh%nTri_loc), source = 0._dp)
    ALLOCATE( SSA%u_b_prev(     mesh%nTri_loc), source = 0._dp)                   ! Velocity solution from previous viscosity iteration
    ALLOCATE( SSA%v_b_prev(     mesh%nTri_loc), source = 0._dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_SSA_solver

  SUBROUTINE solve_SSA( mesh, ice, SSA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Calculate ice velocities by solving the Shallow Ice Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA
    INTEGER,  DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_SSA'
    TYPE(PetscErrorCode)                                         :: perr
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: BC_prescr_mask_b_applied
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: BC_prescr_u_b_applied
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: BC_prescr_v_b_applied
    INTEGER                                                      :: viscosity_iteration_i
    LOGICAL                                                      :: has_converged
    REAL(dp)                                                     :: resid_UV
    REAL(dp)                                                     :: uv_min, uv_max

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL warning('still need to calculate bed roughness and grounded fractions in main UFEMISM model!')

    ! If there is no grounded ice, or no sliding, no need to solve the SSA
    IF ((.NOT. ANY( ice%mask_sheet)) .OR. C%choice_sliding_law == 'no_sliding') THEN
      SSA%u_b = 0._dp
      SSA%v_b = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Handle the optional prescribed u,v boundary conditions
    ALLOCATE( BC_prescr_mask_b_applied( mesh%nTri_loc))
    ALLOCATE( BC_prescr_u_b_applied(    mesh%nTri_loc))
    ALLOCATE( BC_prescr_v_b_applied(    mesh%nTri_loc))
    IF (PRESENT( BC_prescr_mask_b) .OR. PRESENT( BC_prescr_u_b) .OR. PRESENT( BC_prescr_v_b)) THEN
      ! Safety
      IF (.NOT. (PRESENT( BC_prescr_mask_b) .AND. PRESENT( BC_prescr_u_b) .AND. PRESENT( BC_prescr_v_b))) THEN
        CALL crash('need to provide prescribed u,v fields and mask!')
      END IF
      BC_prescr_mask_b_applied = BC_prescr_mask_b
      BC_prescr_u_b_applied    = BC_prescr_u_b
      BC_prescr_v_b_applied    = BC_prescr_v_b
    ELSE
      BC_prescr_mask_b_applied = 0
      BC_prescr_u_b_applied    = 0._dp
      BC_prescr_v_b_applied    = 0._dp
    END IF

    ! Calculate the vertical average of Glen's flow parameter A
    CALL calc_vertically_averaged_flow_parameter( mesh, ice, SSA)

    ! Calculate the driving stress
    CALL calc_driving_stress( mesh, ice, SSA)

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the strain rates for the current velocity solution
      CALL calc_strain_rates( mesh, ice, SSA)

      ! Calculate the effective viscosity for the current velocity solution
      CALL calc_effective_viscosity( mesh, ice, SSA)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      CALL calc_applied_basal_friction_coefficient( mesh, ice, SSA)

      ! Solve the linearised SSA to calculate a new velocity solution
      CALL solve_SSA_linearised( mesh, SSA, BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Limit velocities for improved stability
      CALL apply_velocity_limits( mesh, SSA)

      ! Reduce the change between velocity solutions
      CALL relax_viscosity_iterations( mesh, SSA)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      CALL calc_visc_iter_UV_resid( mesh, SSA, resid_UV)

      ! DENK DROM
      uv_min = MINVAL( SSA%u_b)
      uv_max = MAXVAL( SSA%u_b)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
!      IF (par%master) WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

      ! If the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .FALSE.
      IF (resid_UV < C%visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      END IF

       ! If we've reached the maximum allowed number of iterations without converging, throw a warning
       IF (viscosity_iteration_i > C%visc_it_nit) THEN
         CALL warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
         EXIT viscosity_iteration
       END IF

    END DO viscosity_iteration

    ! Clean up after yourself
    DEALLOCATE( BC_prescr_mask_b_applied)
    DEALLOCATE( BC_prescr_u_b_applied   )
    DEALLOCATE( BC_prescr_v_b_applied   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SSA

  SUBROUTINE remap_SSA_solver( mesh_old, mesh_new, ice, SSA)
    ! Remap the SSA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT) :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_SSA_solver'
    REAL(dp)                                           :: dp_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SSA_solver

! == Assemble and solve the linearised SSA

  SUBROUTINE solve_SSA_linearised( mesh, SSA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Solve the linearised SSA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA
    INTEGER,  DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_SSA_linearised'
    INTEGER                                                      :: ncols, ncols_loc, nrows, nrows_loc, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                              :: A_CSR
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: bb
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: uv_buv
    INTEGER                                                      :: row_buv_glob,row_buv_loc,ti_glob,ti_loc,uv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Store the previous solution
    SSA%u_b_prev = SSA%u_b
    SSA%v_b_prev = SSA%v_b

  ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
  ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * 2      ! from
    ncols_loc       = mesh%nTri_loc * 2
    nrows           = mesh%nTri     * 2      ! to
    nrows_loc       = mesh%nTri_loc * 2
    nnz_per_row_est = mesh%nC_mem * 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    CALL allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector and the solution
    ALLOCATE( bb(     mesh%nTri_loc * 2))
    ALLOCATE( uv_buv( mesh%nTri_loc * 2))

    ! Fill in the current velocity solution
    DO ti_loc = 1, mesh%nTri_loc

      ti_glob = mesh%ti1 + ti_loc - 1

      ! u
      row_buv_glob = mesh%tiuv2n( ti_glob,1)
      row_buv_loc = row_buv_glob - A_CSR%i1 + 1
      uv_buv( row_buv_loc) = SSA%u_b( ti_loc)

      ! v
      row_buv_glob = mesh%tiuv2n( ti_glob,2)
      row_buv_loc = row_buv_glob - A_CSR%i1 + 1
      uv_buv( row_buv_loc) = SSA%v_b( ti_loc)

    END DO ! DO ti_loc = 1, mesh%nTri_loc

  ! == Construct the stiffness matrix for the linearised SSA
  ! ========================================================

    DO row_buv_glob = A_CSR%i1, A_CSR%i2

      row_buv_loc = row_buv_glob - A_CSR%i1 + 1
      ti_glob = mesh%n2tiuv( row_buv_glob,1)
      ti_loc = ti_glob - mesh%ti1 + 1
      uv = mesh%n2tiuv( row_buv_glob,2)

      IF (BC_prescr_mask_b( ti_loc) == 1) THEN
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector: prescribed velocity
        IF     (uv == 1) THEN
          bb( row_buv_loc) = BC_prescr_u_b( ti_loc)
        ELSEIF (uv == 2) THEN
          bb( row_buv_loc) = BC_prescr_v_b( ti_loc)
        ELSE
          CALL crash('uv can only be 1 or 2!')
        END IF

      ELSEIF (mesh%TriBI( ti_glob) == 1 .OR. mesh%TriBI( ti_glob) == 2) THEN
        ! Northern domain border

        CALL calc_SSA_stiffness_matrix_row_BC_north( mesh, SSA, A_CSR, bb, row_buv_glob)

      ELSEIF (mesh%TriBI( ti_glob) == 3 .OR. mesh%TriBI( ti_glob) == 4) THEN
        ! Eastern domain border

        CALL calc_SSA_stiffness_matrix_row_BC_east( mesh, SSA, A_CSR, bb, row_buv_glob)

      ELSEIF (mesh%TriBI( ti_glob) == 5 .OR. mesh%TriBI( ti_glob) == 6) THEN
        ! Northern domain border

        CALL calc_SSA_stiffness_matrix_row_BC_south( mesh, SSA, A_CSR, bb, row_buv_glob)

      ELSEIF (mesh%TriBI( ti_glob) == 7 .OR. mesh%TriBI( ti_glob) == 8) THEN
        ! Western domain border

        CALL calc_SSA_stiffness_matrix_row_BC_west( mesh, SSA, A_CSR, bb, row_buv_glob)

      ELSE
        ! No boundary conditions apply; solve the SSA

        CALL calc_SSA_stiffness_matrix_row_free( mesh, SSA, A_CSR, bb, row_buv_glob)

      END IF

    END DO

  ! == Solve the matrix equation
  ! ============================

    ! Use PETSc to solve the matrix equation
    CALL solve_matrix_equation_CSR_PETSc( A_CSR, bb, uv_buv, SSA%PETSc_rtol, SSA%PETSc_abstol)

    ! Disentangle the u and v components of the velocity solution
    DO ti_loc = 1, mesh%nTri_loc

      ti_glob = mesh%ti1 + ti_loc - 1

      ! u
      row_buv_glob = mesh%tiuv2n( ti_glob,1)
      row_buv_loc = row_buv_glob - A_CSR%i1 + 1
      SSA%u_b( ti_loc) = uv_buv( row_buv_loc)

      ! v
      row_buv_glob = mesh%tiuv2n( ti_glob,2)
      row_buv_loc = row_buv_glob - A_CSR%i1 + 1
      SSA%v_b( ti_loc) = uv_buv( row_buv_loc)

    END DO ! DO ti_loc = 1, mesh%nTri_loc

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( A_CSR)
    DEALLOCATE( bb)
    DEALLOCATE( uv_buv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SSA_linearised

  SUBROUTINE calc_SSA_stiffness_matrix_row_free( mesh, SSA, A_CSR, bb, row_buv_glob)
    ! Add coefficients to this matrix row to represent the linearised SSA
    !
    ! The SSA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_b u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_b v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_b u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_b v = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
    !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx
    !
    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
    !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy
    !
    ! We define the velocities u,v, the basal friction coefficient beta_b, and the driving
    ! stress tau_d on the b-grid (triangles), and the effective viscosity eta and the
    ! product term N = eta H on the a-grid (vertices).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)              :: SSA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT)           :: bb
    INTEGER,                             INTENT(IN)              :: row_buv_glob

    ! Local variables:
    INTEGER                                                      :: row_buv_loc, ti_glob, uv, row_b_glob, row_b_loc
    REAL(dp)                                                     :: N, dN_dx, dN_dy, beta_b, tau_dx, tau_dy
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: single_row_ind
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddx_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dx2_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dxdy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dy2_val
    INTEGER                                                      :: single_row_nnz
    REAL(dp)                                                     :: Au, Av
    INTEGER                                                      :: k, tj_glob, col_bu_glob, col_bv_glob

    ! Relevant indices for this triangle
    row_buv_loc = row_buv_glob - A_CSR%i1 + 1
    ti_glob     = mesh%n2tiuv( row_buv_glob,1)
    uv          = mesh%n2tiuv( row_buv_glob,2)
    row_b_glob  = mesh%ti2n( ti_glob)
    row_b_loc   = row_b_glob - mesh%ti1 + 1

    ! N, dN/dx, dN/dy, beta_b, tau_dx, and tau_dy on this triangle
    N      = SSA%N_b(      row_b_loc)
    dN_dx  = SSA%dN_dx_b(  row_b_loc)
    dN_dy  = SSA%dN_dy_b(  row_b_loc)
    beta_b = SSA%beta_b_b( row_b_loc)
    tau_dx = SSA%tau_dx_b( row_b_loc)
    tau_dy = SSA%tau_dy_b( row_b_loc)

    ! Allocate memory for single matrix rows
    ALLOCATE( single_row_ind(        mesh%nC_mem*2))
    ALLOCATE( single_row_ddx_val(    mesh%nC_mem*2))
    ALLOCATE( single_row_ddy_val(    mesh%nC_mem*2))
    ALLOCATE( single_row_d2dx2_val(  mesh%nC_mem*2))
    ALLOCATE( single_row_d2dxdy_val( mesh%nC_mem*2))
    ALLOCATE( single_row_d2dy2_val(  mesh%nC_mem*2))

    ! Read coefficients of the operator matrices
    CALL read_single_row_CSR_dist( mesh%M2_ddx_b_b   , row_b_glob, single_row_ind, single_row_ddx_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_ddy_b_b   , row_b_glob, single_row_ind, single_row_ddy_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , row_b_glob, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dxdy_b_b, row_b_glob, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , row_b_glob, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    IF (uv == 1) THEN
      ! x-component

      DO k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj_glob     = single_row_ind( k)
        col_bu_glob = mesh%tiuv2n( tj_glob,1)
        col_bv_glob = mesh%tiuv2n( tj_glob,2)

        !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
        !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * N     * single_row_d2dx2_val(  k) + &  ! 4  N    d2u/dx2
             4._dp * dN_dx * single_row_ddx_val(    k) + &  ! 4 dN/dx du/dx
                     N     * single_row_d2dy2_val(  k) + &  !    N    d2u/dy2
                     dN_dy * single_row_ddy_val(    k)      !   dN/dy du/dy
        IF (tj_glob == ti_glob) Au = Au - beta_b            ! - beta_b u

        Av = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2v/dxdy
             2._dp * dN_dx * single_row_ddy_val(    k) + &  ! 2 dN/dx dv/dy
                     dN_dy * single_row_ddx_val(    k)      !   dN/dy dv/dx

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_bu_glob, Au)
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_bv_glob, Av)

      END DO

      ! Load vector
      bb( row_buv_loc) = -SSA%tau_dx_b( row_b_loc)

    ELSEIF (uv == 2) THEN
      ! y-component

      DO k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj_glob     = single_row_ind( k)
        col_bu_glob = mesh%tiuv2n( tj_glob,1)
        col_bv_glob = mesh%tiuv2n( tj_glob,2)

        !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
        !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * N     * single_row_d2dy2_val(  k) + &  ! 4  N    d2v/dy2
             4._dp * dN_dy * single_row_ddy_val(    k) + &  ! 4 dN/dy dv/dy
                     N     * single_row_d2dx2_val(  k) + &  !    N    d2v/dx2
                     dN_dx * single_row_ddx_val(    k)      !   dN/dx dv/dx
        IF (tj_glob == ti_glob) Av = Av - beta_b            ! - beta_b v

        Au = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2u/dxdy
             2._dp * dN_dy * single_row_ddx_val(    k) + &  ! 2 dN/dy du/dx
                     dN_dx * single_row_ddy_val(    k)      !   dN/dx du/dy

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_bu_glob, Au)
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_bv_glob, Av)

      END DO

      ! Load vector
      bb( row_buv_loc) = -SSA%tau_dy_b( row_b_loc)

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

    ! Clean up after yourself
    DEALLOCATE( single_row_ind)
    DEALLOCATE( single_row_ddx_val)
    DEALLOCATE( single_row_ddy_val)
    DEALLOCATE( single_row_d2dx2_val)
    DEALLOCATE( single_row_d2dxdy_val)
    DEALLOCATE( single_row_d2dy2_val)

  END SUBROUTINE calc_SSA_stiffness_matrix_row_free

  SUBROUTINE calc_SSA_stiffness_matrix_row_BC_west( mesh, SSA, A_CSR, bb, row_buv_glob)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! western domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)              :: SSA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT)           :: bb
    INTEGER,                             INTENT(IN)              :: row_buv_glob

    ! Local variables:
    INTEGER                                                      :: row_buv_loc,ti,uv,row_b_glob
    INTEGER                                                      :: k, col_b_glob, tj, col_buv_glob
    INTEGER                                                      :: n, n_neighbours

    row_buv_loc = row_buv_glob - A_CSR%i1 + 1
    ti = mesh%n2tiuv( row_buv_glob,1)
    uv = mesh%n2tiuv( row_buv_glob,2)
    row_b_glob = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_west == 'infinite') THEN
        ! du/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_west == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_west == 'periodic_ISMIP_HOM') THEN
        ! u(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_u_west "' // TRIM( C%BC_u_west) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_west == 'infinite') THEN
        ! dv/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_west == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_west == 'periodic_ISMIP_HOM') THEN
        ! v(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_v_west "' // TRIM( C%BC_v_west) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_SSA_stiffness_matrix_row_BC_west

  SUBROUTINE calc_SSA_stiffness_matrix_row_BC_east( mesh, SSA, A_CSR, bb, row_buv_glob)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! eastern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)              :: SSA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT)           :: bb
    INTEGER,                             INTENT(IN)              :: row_buv_glob

    ! Local variables:
    INTEGER                                                      :: row_buv_loc,ti,uv,row_b_glob
    INTEGER                                                      :: k, col_b_glob, tj, col_buv_glob
    INTEGER                                                      :: n, n_neighbours

    row_buv_loc = row_buv_glob - A_CSR%i1 + 1
    ti = mesh%n2tiuv( row_buv_glob,1)
    uv = mesh%n2tiuv( row_buv_glob,2)
    row_b_glob = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_east == 'infinite') THEN
        ! du/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_east == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_east == 'periodic_ISMIP_HOM') THEN
        ! u(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_u_east "' // TRIM( C%BC_u_east) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_east == 'infinite') THEN
        ! dv/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_east == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_east == 'periodic_ISMIP_HOM') THEN
        ! v(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_v_east "' // TRIM( C%BC_v_east) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_SSA_stiffness_matrix_row_BC_east

  SUBROUTINE calc_SSA_stiffness_matrix_row_BC_south( mesh, SSA, A_CSR, bb, row_buv_glob)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! southern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)              :: SSA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT)           :: bb
    INTEGER,                             INTENT(IN)              :: row_buv_glob

    ! Local variables:
    INTEGER                                                      :: row_buv_loc,ti,uv,row_b_glob
    INTEGER                                                      :: k, col_b_glob, tj, col_buv_glob
    INTEGER                                                      :: n, n_neighbours

    row_buv_loc = row_buv_glob - A_CSR%i1 + 1
    ti = mesh%n2tiuv( row_buv_glob,1)
    uv = mesh%n2tiuv( row_buv_glob,2)
    row_b_glob = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_south == 'infinite') THEN
        ! du/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_south == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_south == 'periodic_ISMIP_HOM') THEN
        ! u(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_u_south "' // TRIM( C%BC_u_south) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_south == 'infinite') THEN
        ! dv/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_south == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_south == 'periodic_ISMIP_HOM') THEN
        ! v(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_v_south "' // TRIM( C%BC_v_south) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_SSA_stiffness_matrix_row_BC_south

  SUBROUTINE calc_SSA_stiffness_matrix_row_BC_north( mesh, SSA, A_CSR, bb, row_buv_glob)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! northern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)              :: SSA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT)           :: bb
    INTEGER,                             INTENT(IN)              :: row_buv_glob

    ! Local variables:
    INTEGER                                                      :: row_buv_loc,ti,uv,row_b_glob
    INTEGER                                                      :: k, col_b_glob, tj, col_buv_glob
    INTEGER                                                      :: n, n_neighbours

    row_buv_loc = row_buv_glob - A_CSR%i1 + 1
    ti = mesh%n2tiuv( row_buv_glob,1)
    uv = mesh%n2tiuv( row_buv_glob,2)
    row_b_glob = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_north == 'infinite') THEN
        ! du/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_north == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_u_north == 'periodic_ISMIP_HOM') THEN
        ! u(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_u_north "' // TRIM( C%BC_u_north) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_north == 'infinite') THEN
        ! dv/dy = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_buv_glob = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_buv_glob, col_buv_glob, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_north == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_buv_glob, row_buv_glob, 1._dp)

        ! Load vector
        bb( row_buv_loc) = 0._dp

      ELSEIF (C%BC_v_north == 'periodic_ISMIP_HOM') THEN
        ! v(x,y) = u(x+L/2,y+L/2)

        ! DENK DROM
        CALL crash('fixme!')

      ELSE
        CALL crash('unknown BC_v_north "' // TRIM( C%BC_v_north) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_SSA_stiffness_matrix_row_BC_north

! == Calculate several intermediate terms in the SSA

  SUBROUTINE calc_driving_stress( mesh, ice, SSA)
    ! Calculate the driving stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_driving_stress'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: Hi_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: dHs_dy_b
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    ALLOCATE( Hi_b(     mesh%nTri_loc))
    ALLOCATE( dHs_dx_b( mesh%nTri_loc))
    ALLOCATE( dHs_dy_b( mesh%nTri_loc))

    ! Calculate Hi, dHs/dx, and dHs/dy on the b-grid
    CALL map_a_b_2D( mesh, ice%Hi, Hi_b    )
    CALL ddx_a_b_2D( mesh, ice%Hs, dHs_dx_b)
    CALL ddy_a_b_2D( mesh, ice%Hs, dHs_dy_b)

    ! Calculate the driving stress
    DO ti = 1, mesh%nTri_loc
      SSA%tau_dx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      SSA%tau_dy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    END DO

    ! Clean up after yourself
    DEALLOCATE( Hi_b    )
    DEALLOCATE( dHs_dx_b)
    DEALLOCATE( dHs_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress

  SUBROUTINE calc_strain_rates( mesh, ice, SSA)
    ! Calculate the strain rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_strain_rates'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the strain rates
    CALL ddx_b_a_2D( mesh, SSA%u_b, SSA%du_dx_a)
    CALL ddy_b_a_2D( mesh, SSA%u_b, SSA%du_dy_a)
    CALL ddx_b_a_2D( mesh, SSA%v_b, SSA%dv_dx_a)
    CALL ddy_b_a_2D( mesh, SSA%v_b, SSA%dv_dy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_rates

  SUBROUTINE calc_vertically_averaged_flow_parameter( mesh, ice, SSA)
    ! Calculate the vertical average of Glen's flow parameter A

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_vertically_averaged_flow_parameter'
    INTEGER                                                      :: vi
    REAL(dp), DIMENSION( mesh%nz)                                :: A_prof

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the vertical average of Glen's flow parameter A
    DO vi = 1, mesh%nV_loc
      A_prof = ice%A_flow_3D( vi,:)
      SSA%A_flow_vav_a( vi) = vertical_average( mesh%zeta, A_prof)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertically_averaged_flow_parameter

  SUBROUTINE calc_effective_viscosity( mesh, ice, SSA)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_viscosity'
    INTEGER                                                      :: vi
    REAL(dp)                                                     :: epsilon_sq

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = 1, mesh%nV_loc

      ! Calculate the square of the effective strain rate epsilon
      epsilon_sq = SSA%du_dx_a( vi)**2 + &
                   SSA%dv_dy_a( vi)**2 + &
                   SSA%du_dx_a( vi) * SSA%dv_dy_a( vi) + &
                   0.25_dp * (SSA%du_dy_a( vi) + SSA%dv_dx_a( vi))**2 + &
                   C%epsilon_sq_0

      ! Calculate the effective viscosity eta
      SSA%eta_a( vi) = 0.5_dp * SSA%A_flow_vav_a( vi)**(-1._dp/  C%n_flow) * (epsilon_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

      ! Safety
      SSA%eta_a( vi) = MAX( SSA%eta_a( vi), C%visc_eff_min)

      ! Calculate the product term N = eta * H
      SSA%N_a( vi) = SSA%eta_a( vi) * MAX( 0.1_dp, ice%Hi( vi))

    END DO

    ! Calculate the product term N and its gradients on the b-grid
    CALL map_a_b_2D( mesh, SSA%N_a, SSA%N_b    )
    CALL ddx_a_b_2D( mesh, SSA%N_a, SSA%dN_dx_b)
    CALL ddy_a_b_2D( mesh, SSA%N_a, SSA%dN_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity

  SUBROUTINE calc_applied_basal_friction_coefficient( mesh, ice, SSA)
    ! Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    ! and scaled with the sub-grid grounded fraction
    !
    ! This is where the sliding law is called!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_applied_basal_friction_coefficient'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!
    CALL calc_basal_friction_coefficient( mesh, ice, SSA%u_b, SSA%v_b)

    ! Map basal friction coefficient beta_b to the b-grid
    CALL map_a_b_2D( mesh, ice%beta_b, SSA%beta_b_b)

    ! Apply the sub-grid grounded fraction
    IF (C%do_GL_subgrid_friction) THEN
      DO ti = 1, mesh%nTri_loc
        SSA%beta_b_b( ti) = SSA%beta_b_b( ti) * ice%fraction_gr_b( ti)**2
      END DO
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_applied_basal_friction_coefficient

! == Some useful tools for improving numerical stability of the viscosity iteration

  SUBROUTINE relax_viscosity_iterations( mesh, SSA)
    ! Reduce the change between velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'relax_viscosity_iterations'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = 1, mesh%nTri_loc
      SSA%u_b( ti) = (C%visc_it_relax * SSA%u_b( ti)) + ((1._dp - C%visc_it_relax) * SSA%u_b_prev( ti))
      SSA%v_b( ti) = (C%visc_it_relax * SSA%v_b( ti)) + ((1._dp - C%visc_it_relax) * SSA%v_b_prev( ti))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_viscosity_iterations

  SUBROUTINE calc_visc_iter_UV_resid( mesh, SSA, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)              :: SSA
    REAL(dp),                            INTENT(OUT)             :: resid_UV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_visc_iter_UV_resid'
    INTEGER                                                      :: ierr
    INTEGER                                                      :: ti
    REAL(dp)                                                     :: res1, res2

    ! Add routine to path
    CALL init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    DO ti = 1, mesh%nTri_loc

      res1 = res1 + (SSA%u_b( ti) - SSA%u_b_prev( ti))**2
      res1 = res1 + (SSA%v_b( ti) - SSA%v_b_prev( ti))**2

      res2 = res2 + (SSA%u_b( ti) + SSA%u_b_prev( ti))**2
      res2 = res2 + (SSA%v_b( ti) + SSA%v_b_prev( ti))**2

    END DO

    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_visc_iter_UV_resid

  SUBROUTINE apply_velocity_limits( mesh, SSA)
    ! Limit velocities for improved stability

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'apply_velocity_limits'
    INTEGER                                                      :: ti
    REAL(dp)                                                     :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = 1, mesh%nTri_loc

      ! Calculate absolute speed
      uabs = SQRT( SSA%u_b( ti)**2 + SSA%v_b( ti)**2)

      ! Reduce velocities if necessary
      IF (uabs > C%vel_max) THEN
        SSA%u_b( ti) = SSA%u_b( ti) * C%vel_max / uabs
        SSA%v_b( ti) = SSA%v_b( ti) * C%vel_max / uabs
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits

END MODULE ice_velocity_SSA
