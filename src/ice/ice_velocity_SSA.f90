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
  USE petsc_basic                                            , ONLY: mat_CSR2petsc, vec_double2petsc, multiply_PETSc_matrix_with_vector_1D, &
                                                                     multiply_PETSc_matrix_with_vector_2D, MatDestroy, VecDestroy
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_SSA
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate_clean
  USE mesh_operators                                         , ONLY: map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D
  USE mesh_zeta                                              , ONLY: vertical_average
  USE sliding_laws                                           , ONLY: calc_basal_friction_coefficient

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

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    CALL allocate_SSA_solver( mesh, SSA)

    ! Calculate combined mesh operators for efficient calculation of the stiffness matrix
    CALL calc_combined_mesh_operators_SSA( mesh, SSA)

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
    ALLOCATE( SSA%u_b(          mesh%nTri_loc))                   ! [m yr^-1] 2-D horizontal ice velocity
    ALLOCATE( SSA%v_b(          mesh%nTri_loc))

    ! Intermediate data fields
    ALLOCATE( SSA%A_flow_vav_a( mesh%nV_loc  ))                   ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter
    ALLOCATE( SSA%du_dx_a(      mesh%nV_loc  ))                   ! [yr^-1] 2-D horizontal strain rates
    ALLOCATE( SSA%du_dy_a(      mesh%nV_loc  ))
    ALLOCATE( SSA%dv_dx_a(      mesh%nV_loc  ))
    ALLOCATE( SSA%dv_dy_a(      mesh%nV_loc  ))
    ALLOCATE( SSA%eta_a(        mesh%nV_loc  ))                   ! Effective viscosity
    ALLOCATE( SSA%N_a(          mesh%nV_loc  ))                   ! Product term N = eta * H
    ALLOCATE( SSA%N_b(          mesh%nTri_loc))
    ALLOCATE( SSA%dN_dx_b(      mesh%nTri_loc))                   ! Gradients of N
    ALLOCATE( SSA%dN_dy_b(      mesh%nTri_loc))
    ALLOCATE( SSA%beta_b_b(     mesh%nTri_loc))                   ! Friction coefficient (tau_b = u * beta_b)
    ALLOCATE( SSA%tau_dx_b(     mesh%nTri_loc))                   ! Driving stress
    ALLOCATE( SSA%tau_dy_b(     mesh%nTri_loc))
    ALLOCATE( SSA%u_b_prev(     mesh%nTri_loc))                   ! Velocity solution from previous viscosity iteration
    ALLOCATE( SSA%v_b_prev(     mesh%nTri_loc))

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

    ! Add routine to path
    CALL init_routine( routine_name)

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
      IF (.NOT. (PRESENT(BC_prescr_mask_b) .AND. PRESENT(BC_prescr_u_b) .AND. PRESENT(BC_prescr_v_b))) THEN
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

!    ! Calculate boundary conditions mask matrices
!    CALL calc_SSA_BC_mask_matrices( mesh, SSA, BC_prescr_mask_b_applied)
!
!    ! Calculate the driving stress
!    CALL calc_driving_stress( mesh, ice, SSA)
!
!    ! Calculate the basal yield stress tau_c
!    CALL calc_basal_conditions( mesh, ice)
!
!    ! Determine sub-mesh grounded fractions for scaling the basal friction
!    CALL determine_grounded_fractions( mesh, ice)
!
!    ! The viscosity iteration
!    viscosity_iteration_i = 0
!    has_converged         = .FALSE.
!    viscosity_iteration: DO WHILE (.NOT. has_converged)
!      viscosity_iteration_i = viscosity_iteration_i + 1
!
!      ! Calculate the strain rates for the current velocity solution
!      CALL calc_strain_rates( mesh, ice, SSA)
!
!      ! Calculate the effective viscosity for the current velocity solution
!      CALL calc_effective_viscosity( mesh, ice, SSA)
!
!      ! Calculate the basal friction coefficient betab for the current velocity solution
!      CALL calc_applied_basal_friction_coefficient( mesh, ice, SSA)
!
!      ! Solve the linearised SSA to calculate a new velocity solution
!      CALL solve_SSA_linearised( mesh, SSA, BC_prescr_u_b_applied, BC_prescr_v_b_applied)
!
!      ! Limit velocities for improved stability
!      CALL apply_velocity_limits( mesh, SSA)
!
!      ! Reduce the change between velocity solutions
!      CALL relax_viscosity_iterations( mesh, SSA, C%DIVA_visc_it_relax)
!
!      ! Calculate the L2-norm of the two consecutive velocity solutions
!      CALL calc_visc_iter_UV_resid( mesh, SSA, resid_UV)
!
!!      ! DENK DROM
!!      IF (par%master) WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ', u = [', MINVAL( SSA%u_b), ' - ', MAXVAL( SSA%u_b), '], resid = ', resid_UV
!
!      ! If the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
!      has_converged = .FALSE.
!      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
!        has_converged = .TRUE.
!      ELSEIF (viscosity_iteration_i > C%DIVA_visc_it_nit) THEN
!        has_converged = .TRUE.
!      END IF
!
!    END DO viscosity_iteration

    ! Clean up after yourself
    DEALLOCATE( BC_prescr_mask_b_applied)
    DEALLOCATE( BC_prescr_u_b_applied   )
    DEALLOCATE( BC_prescr_v_b_applied   )
    CALL MatDestroy( SSA%m_free     , perr)
    CALL MatDestroy( SSA%m_BC_west  , perr)
    CALL MatDestroy( SSA%m_BC_east  , perr)
    CALL MatDestroy( SSA%m_BC_south , perr)
    CALL MatDestroy( SSA%m_BC_north , perr)
    CALL MatDestroy( SSA%m_BC_prescr, perr)

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

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the vertical average of Glen's flow parameter A
    DO vi = 1, mesh%nV_loc
      SSA%A_flow_vav_a( vi) = vertical_average( mesh%zeta, ice%A_flow_3D( vi,:))
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
    CALL crash('need to implement sliding laws!')
    CALL calc_basal_friction_coefficient( mesh, ice, SSA%u_b, SSA%v_b)

    ! Map basal friction coefficient beta_b to the b-grid
    CALL map_a_b_2D( mesh, ice%beta_b, SSA%beta_b_b)

    ! Apply the sub-grid grounded fraction
    IF (C%do_GL_subgrid_friction) THEN
      DO ti = 1, mesh%nTri_loc
        SSA%beta_b_b( ti) = SSA%beta_b_b( ti) * ice%f_grnd_b( ti)**2
      END DO
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_applied_basal_friction_coefficient

! == Calculate combined mesh operators for efficient calculation of the stiffness matrix

  SUBROUTINE calc_combined_mesh_operators_SSA( mesh, SSA)
    ! Calculate combined mesh operators for efficient calculation of the stiffness matrix

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_SSA),  INTENT(INOUT)           :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_combined_mesh_operators_SSA'
    TYPE(PetscErrorCode)                                         :: perr
    TYPE(tMat)                                                   :: M2_ddx_b_b
    TYPE(tMat)                                                   :: M2_ddy_b_b
    TYPE(tMat)                                                   :: M2_d2dx2_b_b
    TYPE(tMat)                                                   :: M2_d2dxdy_b_b
    TYPE(tMat)                                                   :: M2_d2dy2_b_b
    TYPE(tMat)                                                   :: M_dudx_buv_b
    TYPE(tMat)                                                   :: M_dudy_buv_b
    TYPE(tMat)                                                   :: M_d2udx2_buv_b
    TYPE(tMat)                                                   :: M_d2udxdy_buv_b
    TYPE(tMat)                                                   :: M_d2udy2_buv_b
    TYPE(tMat)                                                   :: M_dvdx_buv_b
    TYPE(tMat)                                                   :: M_dvdy_buv_b
    TYPE(tMat)                                                   :: M_d2vdx2_buv_b
    TYPE(tMat)                                                   :: M_d2vdxdy_buv_b
    TYPE(tMat)                                                   :: M_d2vdy2_buv_b
    TYPE(tMat)                                                   :: M1
    REAL(dp), DIMENSION(mesh%nTri_loc)                           :: d_4
    TYPE(tVec)                                                   :: v_4

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Convert gradient operators to PETSc format
    CALL mat_CSR2petsc( mesh%M2_ddx_b_b   , M2_ddx_b_b   )
    CALL mat_CSR2petsc( mesh%M2_ddy_b_b   , M2_ddy_b_b   )
    CALL mat_CSR2petsc( mesh%M2_d2dx2_b_b , M2_d2dx2_b_b )
    CALL mat_CSR2petsc( mesh%M2_d2dxdy_b_b, M2_d2dxdy_b_b)
    CALL mat_CSR2petsc( mesh%M2_d2dy2_b_b , M2_d2dy2_b_b )

    ! Vector containing only  and 4's
    d_4 = 4._dp
    CALL vec_double2petsc( d_4, v_4)

    ! Partial derivatives of u on the b-grid
    CALL MatMatMult( M2_ddx_b_b   , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx_buv_b   , perr)
    CALL MatMatMult( M2_ddy_b_b   , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudy_buv_b   , perr)
    CALL MatMatMult( M2_d2dx2_b_b , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udx2_buv_b , perr)
    CALL MatMatMult( M2_d2dxdy_b_b, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udxdy_buv_b, perr)
    CALL MatMatMult( M2_d2dy2_b_b , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udy2_buv_b , perr)

    ! Partial derivatives of v on the b-grid
    CALL MatMatMult( M2_ddx_b_b   , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx_buv_b   , perr)
    CALL MatMatMult( M2_ddy_b_b   , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdy_buv_b   , perr)
    CALL MatMatMult( M2_d2dx2_b_b , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdx2_buv_b , perr)
    CALL MatMatMult( M2_d2dxdy_b_b, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdxdy_buv_b, perr)
    CALL MatMatMult( M2_d2dy2_b_b , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdy2_buv_b , perr)

    ! M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b
    CALL MatDuplicate( M_d2udx2_buv_b, MAT_COPY_VALUES, SSA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, perr)
    CALL MatDiagonalScale( SSA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( SSA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, 3._dp, M_d2vdxdy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( SSA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, 1._dp, M_d2udy2_buv_b , DIFFERENT_NONZERO_PATTERN, perr)

    ! M_4_dudx_p_2_dvdy_buv_b
    CALL MatDuplicate( M_dudx_buv_b, MAT_COPY_VALUES, SSA%M_4_dudx_p_2_dvdy_buv_b, perr)
    CALL MatDiagonalScale( SSA%M_4_dudx_p_2_dvdy_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( SSA%M_4_dudx_p_2_dvdy_buv_b, 2._dp, M_dvdy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! M_dudy_p_dvdx_buv_b
    CALL MatDuplicate( M_dudy_buv_b, MAT_COPY_VALUES, SSA%M_dudy_p_dvdx_buv_b, perr)
    CALL MatAXPY( SSA%M_dudy_p_dvdx_buv_b, 1._dp, M_dvdx_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b
    CALL MatDuplicate( M_d2vdy2_buv_b, MAT_COPY_VALUES, SSA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, perr)
    CALL MatDiagonalScale( SSA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( SSA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, 3._dp, M_d2udxdy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( SSA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, 1._dp, M_d2vdx2_buv_b , DIFFERENT_NONZERO_PATTERN, perr)

    ! M_4_dvdy_p_2_dudx_buv_b
    CALL MatDuplicate( M_dvdy_buv_b, MAT_COPY_VALUES, SSA%M_4_dvdy_p_2_dudx_buv_b, perr)
    CALL MatDiagonalScale( SSA%M_4_dvdy_p_2_dudx_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( SSA%M_4_dvdy_p_2_dudx_buv_b, 2._dp, M_dudx_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! M_dvdx_p_dudy_buv_b
    CALL MatDuplicate( M_dvdx_buv_b, MAT_COPY_VALUES, SSA%M_dvdx_p_dudy_buv_b, perr)
    CALL MatAXPY( SSA%M_dvdx_p_dudy_buv_b, 1._dp, M_dudy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( M2_ddx_b_b     , perr)
    CALL MatDestroy( M2_ddy_b_b     , perr)
    CALL MatDestroy( M2_d2dx2_b_b   , perr)
    CALL MatDestroy( M2_d2dxdy_b_b  , perr)
    CALL MatDestroy( M2_d2dy2_b_b   , perr)
    CALL MatDestroy( M_dudx_buv_b   , perr)
    CALL MatDestroy( M_dudy_buv_b   , perr)
    CALL MatDestroy( M_d2udx2_buv_b , perr)
    CALL MatDestroy( M_d2udxdy_buv_b, perr)
    CALL MatDestroy( M_d2udy2_buv_b , perr)
    CALL MatDestroy( M_dvdx_buv_b   , perr)
    CALL MatDestroy( M_dvdy_buv_b   , perr)
    CALL MatDestroy( M_d2vdx2_buv_b , perr)
    CALL MatDestroy( M_d2vdxdy_buv_b, perr)
    CALL MatDestroy( M_d2vdy2_buv_b , perr)
    CALL VecDestroy( v_4            , perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_combined_mesh_operators_SSA

END MODULE ice_velocity_SSA
