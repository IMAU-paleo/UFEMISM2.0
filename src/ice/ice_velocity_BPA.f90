MODULE ice_velocity_BPA

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (BPA)

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: save_variable_as_netcdf_dp_1D, save_variable_as_netcdf_dp_2D
  USE petsc_basic                                            , ONLY: solve_matrix_equation_CSR_PETSc
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_BPA
  USE parameters
  USE mesh_operators                                         , ONLY: map_a_b_2D, map_a_b_3D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_3D, ddy_b_a_3D, &
                                                                     calc_3D_gradient_bk_ak, calc_3D_gradient_bk_bks, map_ak_bks, map_bks_ak, &
                                                                     calc_3D_gradient_ak_bk, calc_3D_gradient_bks_bk, calc_3D_matrix_operators_mesh, &
                                                                     map_b_a_3D
  USE mesh_zeta                                              , ONLY: vertical_average
  USE sliding_laws                                           , ONLY: calc_basal_friction_coefficient
  USE mesh_utilities                                         , ONLY: find_ti_copy_ISMIP_HOM_periodic
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file, open_existing_netcdf_file_for_writing
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, add_time_dimension_to_file, &
                                                                     add_zeta_dimension_to_file, add_field_mesh_dp_3D_b, write_time_to_file, write_to_field_multopt_mesh_dp_3D_b
  USE netcdf_input                                           , ONLY: read_field_from_mesh_file_3D_b
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_2D
  USE ice_flow_laws                                          , ONLY: calc_effective_viscosity_Glen_3D_uv_only, calc_ice_rheology_Glen
  USE ice_model_utilities                                    , ONLY: calc_zeta_gradients
  USE reallocate_mod                                         , ONLY: reallocate_bounds, reallocate_clean
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Main routines

  SUBROUTINE initialise_BPA_solver( mesh, BPA, region_name)
    ! Initialise the BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(OUT)   :: BPA
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BPA_solver'
    CHARACTER(LEN=256)                                 :: choice_initial_velocity

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    CALL allocate_BPA_solver( mesh, BPA)

    ! Determine the choice of initial velocities for this model region
    IF     (region_name == 'NAM') THEN
      choice_initial_velocity  = C%choice_initial_velocity_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_initial_velocity  = C%choice_initial_velocity_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_initial_velocity  = C%choice_initial_velocity_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_initial_velocity  = C%choice_initial_velocity_ANT
    ELSE
      CALL crash('unknown model region "' // region_name // '"!')
    END IF

    ! Initialise velocities according to the specified method
    IF     (choice_initial_velocity == 'zero') THEN
      BPA%u_bk = 0._dp
      BPA%v_bk = 0._dp
    ELSEIF (choice_initial_velocity == 'read_from_file') THEN
      CALL initialise_BPA_velocities_from_file( mesh, BPA, region_name)
    ELSE
      CALL crash('unknown choice_initial_velocity "' // TRIM( choice_initial_velocity) // '"!')
    END IF

    ! Set tolerances for PETSc matrix solver for the linearised BPA
    BPA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    BPA%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BPA_solver

  SUBROUTINE solve_BPA( mesh, ice, BPA, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    ! Calculate ice velocities by solving the Blatter-Pattyn Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_BPA'
    LOGICAL                                                      :: grounded_ice_exists
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE                      :: BC_prescr_mask_bk_applied
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      :: BC_prescr_u_bk_applied
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      :: BC_prescr_v_bk_applied
    INTEGER                                                      :: viscosity_iteration_i
    LOGICAL                                                      :: has_converged
    REAL(dp)                                                     :: resid_UV, resid_UV_prev
    REAL(dp)                                                     :: uv_min, uv_max
    REAL(dp)                                                     :: visc_it_relax_applied
    REAL(dp)                                                     :: Glens_flow_law_epsilon_sq_0_applied
    INTEGER                                                      :: nit_diverg_consec

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there is no grounded ice, or no sliding, no need to solve the BPA
    grounded_ice_exists = ANY( ice%mask_grounded_ice)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. grounded_ice_exists) THEN
      BPA%u_bk = 0._dp
      BPA%v_bk = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Handle the optional prescribed u,v boundary conditions
    ALLOCATE( BC_prescr_mask_bk_applied( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( BC_prescr_u_bk_applied(    mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( BC_prescr_v_bk_applied(    mesh%ti1:mesh%ti2,mesh%nz))
    IF (PRESENT( BC_prescr_mask_bk) .OR. PRESENT( BC_prescr_u_bk) .OR. PRESENT( BC_prescr_v_bk)) THEN
      ! Safety
      IF (.NOT. (PRESENT( BC_prescr_mask_bk) .AND. PRESENT( BC_prescr_u_bk) .AND. PRESENT( BC_prescr_v_bk))) THEN
        CALL crash('need to provide prescribed u,v fields and mask!')
      END IF
      BC_prescr_mask_bk_applied = BC_prescr_mask_bk
      BC_prescr_u_bk_applied    = BC_prescr_u_bk
      BC_prescr_v_bk_applied    = BC_prescr_v_bk
    ELSE
      BC_prescr_mask_bk_applied = 0
      BC_prescr_u_bk_applied    = 0._dp
      BC_prescr_v_bk_applied    = 0._dp
    END IF

    ! Calculate zeta gradients
    CALL calc_zeta_gradients( mesh, ice)

    ! Calculate 3-D matrix operators for the current ice geometry
    CALL calc_3D_matrix_operators_mesh( mesh, ice)

    ! Calculate the driving stress
    CALL calc_driving_stress( mesh, ice, BPA)

    ! Adaptive relaxation parameter for the viscosity iteration
    resid_UV                            = 1E9_dp
    nit_diverg_consec                   = 0
    visc_it_relax_applied               = C%visc_it_relax
    Glens_flow_law_epsilon_sq_0_applied = C%Glens_flow_law_epsilon_sq_0

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the strain rates for the current velocity solution
      CALL calc_strain_rates( mesh, BPA)

      ! Calculate the effective viscosity for the current velocity solution
      CALL calc_effective_viscosity( mesh, ice, BPA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      CALL calc_applied_basal_friction_coefficient( mesh, ice, BPA)

      ! Solve the linearised BPA to calculate a new velocity solution
      CALL solve_BPA_linearised( mesh, ice, BPA, BC_prescr_mask_bk_applied, BC_prescr_u_bk_applied, BC_prescr_v_bk_applied)

      ! Limit velocities for improved stability
      CALL apply_velocity_limits( mesh, BPA)

      ! Reduce the change between velocity solutions
      CALL relax_viscosity_iterations( mesh, BPA, visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      CALL calc_visc_iter_UV_resid( mesh, BPA, resid_UV)

      ! If the viscosity iteration diverges, lower the relaxation parameter
      IF (resid_UV > resid_UV_prev) THEN
        nit_diverg_consec = nit_diverg_consec + 1
      ELSE
        nit_diverg_consec = 0
      END IF
      IF (nit_diverg_consec > 2) THEN
        nit_diverg_consec = 0
        visc_it_relax_applied               = visc_it_relax_applied               * 0.9_dp
        Glens_flow_law_epsilon_sq_0_applied = Glens_flow_law_epsilon_sq_0_applied * 1.2_dp
      END IF
      IF (visc_it_relax_applied <= 0.05_dp .OR. Glens_flow_law_epsilon_sq_0_applied >= 1E-5_dp) THEN
        IF (visc_it_relax_applied < 0.05_dp) THEN
          CALL crash('viscosity iteration still diverges even with very low relaxation factor!')
        ELSEIF (Glens_flow_law_epsilon_sq_0_applied > 1E-5_dp) THEN
          CALL crash('viscosity iteration still diverges even with very high effective strain rate regularisation!')
        END IF
      END IF

      ! DENK DROM
      uv_min = MINVAL( BPA%u_bk)
      uv_max = MAXVAL( BPA%u_bk)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
!      IF (par%master) WRITE(0,*) '    BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

      ! If the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .FALSE.
      IF (resid_UV < C%visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      END IF

      ! If we've reached the maximum allowed number of iterations without converging, throw a warning
      IF (viscosity_iteration_i > C%visc_it_nit) THEN
        IF (par%master) CALL warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        EXIT viscosity_iteration
      END IF

    END DO viscosity_iteration

    ! Clean up after yourself
    DEALLOCATE( BC_prescr_mask_bk_applied)
    DEALLOCATE( BC_prescr_u_bk_applied   )
    DEALLOCATE( BC_prescr_v_bk_applied   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_BPA

  SUBROUTINE remap_BPA_solver( mesh_old, mesh_new, BPA)
    ! Remap the BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT) :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_BPA_solver'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: u_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: v_ak

    ! Add routine to path
    CALL init_routine( routine_name)

  ! Remap the fields that are re-used during the viscosity iteration
  ! ================================================================

    ! Allocate memory for velocities on the a-grid (vertices)
    ALLOCATE( u_ak( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    ALLOCATE( v_ak( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map velocities from the triangles of the old mesh to the vertices of the old mesh
    CALL map_b_a_3D( mesh_old, BPA%u_bk, u_ak)
    CALL map_b_a_3D( mesh_old, BPA%v_bk, v_ak)

    ! Remap velocities from the vertices of the old mesh to the vertices of the new mesh
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, u_ak, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, v_ak, '2nd_order_conservative')

    ! Reallocate memory for the velocities on the triangles
    CALL reallocate_bounds( BPA%u_bk                        , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( BPA%v_bk                        , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map velocities from the vertices of the new mesh to the triangles of the new mesh
    CALL map_a_b_3D( mesh_new, u_ak, BPA%u_bk)
    CALL map_a_b_3D( mesh_new, v_ak, BPA%v_bk)

    ! Clean up after yourself
    DEALLOCATE( u_ak)
    DEALLOCATE( v_ak)

  ! Reallocate everything else
  ! ==========================

    CALL reallocate_bounds( BPA%du_dx_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! [yr^-1] 2-D horizontal strain rates
    CALL reallocate_bounds( BPA%du_dy_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    CALL reallocate_bounds( BPA%du_dz_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    CALL reallocate_bounds( BPA%dv_dx_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! [yr^-1] 2-D horizontal strain rates
    CALL reallocate_bounds( BPA%dv_dy_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    CALL reallocate_bounds( BPA%dv_dz_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    CALL reallocate_bounds( BPA%du_dx_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)       ! [yr^-1] 2-D horizontal strain rates
    CALL reallocate_bounds( BPA%du_dy_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( BPA%du_dz_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( BPA%dv_dx_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)       ! [yr^-1] 2-D horizontal strain rates
    CALL reallocate_bounds( BPA%dv_dy_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( BPA%dv_dz_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( BPA%eta_ak                      , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! Effective viscosity
    CALL reallocate_bounds( BPA%eta_bks                     , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    CALL reallocate_bounds( BPA%eta_bk                      , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )
    CALL reallocate_bounds( BPA%deta_dx_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    CALL reallocate_bounds( BPA%deta_dy_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    CALL reallocate_bounds( BPA%deta_dz_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    CALL reallocate_bounds( BPA%basal_friction_coefficient_b, mesh_new%ti1 , mesh_new%ti2               )       ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    CALL reallocate_bounds( BPA%dh_dx_b                     , mesh_new%ti1 , mesh_new%ti2               )       ! Surface slope
    CALL reallocate_bounds( BPA%dh_dy_b                     , mesh_new%ti1 , mesh_new%ti2               )
    CALL reallocate_bounds( BPA%db_dx_b                     , mesh_new%ti1 , mesh_new%ti2               )       ! Basal slope
    CALL reallocate_bounds( BPA%db_dy_b                     , mesh_new%ti1 , mesh_new%ti2               )
    CALL reallocate_bounds( BPA%tau_dx_b                    , mesh_new%ti1 , mesh_new%ti2               )       ! Driving stress
    CALL reallocate_bounds( BPA%tau_dy_b                    , mesh_new%ti1 , mesh_new%ti2               )
    CALL reallocate_clean ( BPA%u_bk_prev                   , mesh_new%nTri              , mesh_new%nz  )       ! Velocity solution from previous viscosity iteration
    CALL reallocate_clean ( BPA%v_bk_prev                   , mesh_new%nTri              , mesh_new%nz  )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BPA_solver

! == Assemble and solve the linearised BPA

  SUBROUTINE solve_BPA_linearised( mesh, ice, BPA, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    ! Solve the linearised BPA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA
    INTEGER,  DIMENSION(mesh%ti1:mesh%ti2,mesh%nz),  INTENT(IN)  :: BC_prescr_mask_bk      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,mesh%nz),  INTENT(IN)  :: BC_prescr_u_bk         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,mesh%nz),  INTENT(IN)  :: BC_prescr_v_bk         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_BPA_linearised'
    INTEGER                                                      :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                              :: A_CSR
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: bb
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: uv_bkuv
    INTEGER                                                      :: row_tikuv,ti,k,uv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Store the previous solution
    CALL gather_to_all_dp_2D( BPA%u_bk, BPA%u_bk_prev)
    CALL gather_to_all_dp_2D( BPA%v_bk, BPA%v_bk_prev)

  ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
  ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * mesh%nz * 2      ! from
    ncols_loc       = mesh%nTri_loc * mesh%nz * 2
    nrows           = mesh%nTri     * mesh%nz * 2      ! to
    nrows_loc       = mesh%nTri_loc * mesh%nz * 2
    nnz_est_proc    = mesh%M2_ddx_bk_bk%nnz   * 4

    CALL allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector and the solution
    ALLOCATE( bb(      A_CSR%i1:A_CSR%i2))
    ALLOCATE( uv_bkuv( A_CSR%i1:A_CSR%i2))

    ! Fill in the current velocity solution
    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, mesh%nz

      ! u
      row_tikuv = mesh%tikuv2n( ti,k,1)
      uv_bkuv( row_tikuv) = BPA%u_bk( ti,k)

      ! v
      row_tikuv = mesh%tikuv2n( ti,k,2)
      uv_bkuv( row_tikuv) = BPA%v_bk( ti,k)

    END DO ! DO k  = 1, mesh%nz
    END DO ! DO ti = mesh%ti1, mesh%ti2

  ! == Construct the stiffness matrix for the linearised BPA
  ! ========================================================

    DO row_tikuv = A_CSR%i1, A_CSR%i2

      ti = mesh%n2tikuv( row_tikuv,1)
      k  = mesh%n2tikuv( row_tikuv,2)
      uv = mesh%n2tikuv( row_tikuv,3)

      IF (BC_prescr_mask_bk( ti,k) == 1) THEN
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector: prescribed velocity
        IF     (uv == 1) THEN
          bb( row_tikuv) = BC_prescr_u_bk( ti,k)
        ELSEIF (uv == 2) THEN
          bb( row_tikuv) = BC_prescr_v_bk( ti,k)
        ELSE
          CALL crash('uv can only be 1 or 2!')
        END IF

      ELSEIF (mesh%TriBI( ti) == 1 .OR. mesh%TriBI( ti) == 2) THEN
        ! Northern domain border

        CALL calc_BPA_stiffness_matrix_row_BC_north( mesh, BPA, A_CSR, bb, row_tikuv)

      ELSEIF (mesh%TriBI( ti) == 3 .OR. mesh%TriBI( ti) == 4) THEN
        ! Eastern domain border

        CALL calc_BPA_stiffness_matrix_row_BC_east( mesh, BPA, A_CSR, bb, row_tikuv)

      ELSEIF (mesh%TriBI( ti) == 5 .OR. mesh%TriBI( ti) == 6) THEN
        ! Southern domain border

        CALL calc_BPA_stiffness_matrix_row_BC_south( mesh, BPA, A_CSR, bb, row_tikuv)

      ELSEIF (mesh%TriBI( ti) == 7 .OR. mesh%TriBI( ti) == 8) THEN
        ! Western domain border

        CALL calc_BPA_stiffness_matrix_row_BC_west( mesh, BPA, A_CSR, bb, row_tikuv)

      ELSEIF (k == 1) THEN
        ! Ice surface

        CALL calc_BPA_stiffness_matrix_row_BC_surf( mesh, ice, BPA, A_CSR, bb, row_tikuv)

      ELSEIF (k == mesh%nz) THEN
        ! Ice base

        CALL calc_BPA_stiffness_matrix_row_BC_base( mesh, ice, BPA, A_CSR, bb, row_tikuv)

      ELSE
        ! No boundary conditions apply; solve the BPA

        CALL calc_BPA_stiffness_matrix_row_free( mesh, BPA, A_CSR, bb, row_tikuv)

      END IF

    END DO ! DO row_tikuv = A_CSR%i1, A_CSR%i2

  ! == Solve the matrix equation
  ! ============================

    ! Use PETSc to solve the matrix equation
    CALL solve_matrix_equation_CSR_PETSc( A_CSR, bb, uv_bkuv, BPA%PETSc_rtol, BPA%PETSc_abstol)

    ! Disentangle the u and v components of the velocity solution
    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, mesh%nz

      ! u
      row_tikuv = mesh%tikuv2n( ti,k,1)
      BPA%u_bk( ti,k) = uv_bkuv( row_tikuv)

      ! v
      row_tikuv = mesh%tikuv2n( ti,k,2)
      BPA%v_bk( ti,k) = uv_bkuv( row_tikuv)

    END DO ! DO k  = 1, mesh%nz
    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( A_CSR)
    DEALLOCATE( bb)
    DEALLOCATE( uv_bkuv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_BPA_linearised

  SUBROUTINE calc_BPA_stiffness_matrix_row_free( mesh, BPA, A_CSR, bb, row_tikuv)
    ! Add coefficients to this matrix row to represent the linearised BPA
    !
    ! The BPA reads;
    !
    !   d/dx [ 2 eta ( 2 du/dx + dv/dy )] + d/dy [ eta ( du/dy + dv/dx)] + d/dz [ eta du/dz] = -tau_dx
    !
    !   d/dy [ 2 eta ( 2 dv/dy + du/dx )] + d/dx [ eta ( dv/dx + du/dy)] + d/dz [ eta dv/dz] = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 eta d2u/dx2 + 4 deta/dx du/dx + 2 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !     eta d2u/dy2 +   deta/dy du/dy +   eta d2v/dxdy +   deta/dy dv/dx + ...
    !     eta d2u/dz2 +   deta/dz du/dz = -tau_dx
    !
    !   4 eta d2v/dy2 + 4 deta/dy dv/dy + 2 eta d2u/dxdy + 2 deta/dy du/dx + ...
    !     eta d2v/dx2 +   deta/dx dv/dx +   eta d2u/dxdy +   deta/dx du/dy + ...
    !     eta d2v/dz2 +   deta/dz dv/dz = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2 + deta/dy du/dy + eta d2u/dz2 + deta/dz du/dz + ...
    !   3 eta d2v/dxdy + 2 deta/dx dv/dy +               deta/dy dv/dx = -tau_dx
    !
    !   4 eta d2v/dy2  + 4 deta/dy dv/dy + eta d2v/dx2 + deta/dx dv/dx + eta d2v/dz2 + deta/dz dv/dz + ...
    !   3 eta d2u/dxdy + 2 deta/dy du/dx +               deta/dx du/dy = -tau_dy
    !
    ! We define the velocities u,v, and the driving stress tau_d on the bk-grid
    ! (triangles, regular vertical), and the effective viscosity eta on the ak-grid
    ! (vertices, regular vertical) and the bks-grid (triangles, staggered vertical).
    ! From this, we then calculate the gradients of eta on the bk-grid.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tikuv

    ! Local variables:
    INTEGER                                                      :: ti, k, uv
    REAL(dp)                                                     :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: single_row_ind
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddx_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddz_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dx2_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dxdy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dy2_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dz2_val
    INTEGER                                                      :: single_row_nnz
    INTEGER                                                      :: row_tik
    REAL(dp)                                                     :: Au, Av
    INTEGER                                                      :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta     = BPA%eta_bk(     ti,k)
    deta_dx = BPA%deta_dx_bk( ti,k)
    deta_dy = BPA%deta_dy_bk( ti,k)
    deta_dz = BPA%deta_dz_bk( ti,k)
    tau_dx  = BPA%tau_dx_b(   ti  )
    tau_dy  = BPA%tau_dy_b(   ti  )

    ! Allocate memory for single matrix rows
    ALLOCATE( single_row_ind(        mesh%nC_mem*3*2))
    ALLOCATE( single_row_ddx_val(    mesh%nC_mem*3*2))
    ALLOCATE( single_row_ddy_val(    mesh%nC_mem*3*2))
    ALLOCATE( single_row_ddz_val(    mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dx2_val(  mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dxdy_val( mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dy2_val(  mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dz2_val(  mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = mesh%tik2n( ti,k)
    CALL read_single_row_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_ddz_bk_bk   , row_tik, single_row_ind, single_row_ddz_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dz2_bk_bk , row_tik, single_row_ind, single_row_d2dz2_val , single_row_nnz)

    IF (uv == 1) THEN
      ! x-component

      DO n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        !   4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2 + deta/dy du/dy + eta d2u/dz2 + deta/dz du/dz + ...
        !   3 eta d2v/dxdy + 2 deta/dx dv/dy +               deta/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * eta     * single_row_d2dx2_val(  n) + &  ! 4  eta    d2u/dx2
             4._dp * deta_dx * single_row_ddx_val(    n) + &  ! 4 deta/dx  du/dx
                     eta     * single_row_d2dy2_val(  n) + &  !    eta    d2u/dy2
                     deta_dy * single_row_ddy_val(    n) + &  !   deta/dy  du/dy
                     eta     * single_row_d2dz2_val(  n) + &  !    eta    d2u/dz2
                     deta_dz * single_row_ddz_val(    n)      !   deta/dz  du/dz

        Av = 3._dp * eta     * single_row_d2dxdy_val( n) + &  ! 3  eta    d2v/dxdy
             2._dp * deta_dx * single_row_ddy_val(    n) + &  ! 2 deta/dx  dv/dy
                     deta_dy * single_row_ddx_val(    n)      !   deta/dy  dv/dx

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)

      END DO

      ! Load vector
      bb( row_tikuv) = -tau_dx

    ELSEIF (uv == 2) THEN
      ! y-component

      DO n = 1, single_row_nnz

        ! Releuant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        !   4 eta d2v/dy2  + 4 deta/dy dv/dy + eta d2v/dx2 + deta/dx dv/dx + eta d2v/dz2 + deta/dz dv/dz + ...
        !   3 eta d2u/dxdy + 2 deta/dy du/dx +               deta/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * eta     * single_row_d2dy2_val(  n) + &  ! 4  eta    d2v/dy2
             4._dp * deta_dy * single_row_ddy_val(    n) + &  ! 4 deta/dy  dv/dy
                     eta     * single_row_d2dx2_val(  n) + &  !    eta    d2v/dx2
                     deta_dx * single_row_ddx_val(    n) + &  !   deta/dx  dv/dx
                     eta     * single_row_d2dz2_val(  n) + &  !    eta    d2v/dz2
                     deta_dz * single_row_ddz_val(    n)      !   deta/dz  dv/dz

        Au = 3._dp * eta     * single_row_d2dxdy_val( n) + &  ! 3  eta    d2u/dxdy
             2._dp * deta_dy * single_row_ddx_val(    n) + &  ! 2 deta/dy  du/dx
                     deta_dx * single_row_ddy_val(    n)      !   deta/dx  du/dy

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)

      END DO

      ! Load uector
      bb( row_tikuv) = -tau_dy

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

  END SUBROUTINE calc_BPA_stiffness_matrix_row_free

  SUBROUTINE calc_BPA_stiffness_matrix_row_BC_surf( mesh, ice, BPA, A_CSR, bb, row_tikuv)
    ! Add coefficients to this matrix row to represent the boundary conditions to the BPA at the ice surface
    !
    ! At the ice surface (k=1), the zero-stress boundary condition implies that:
    !
    ! [1]     2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx) - du/dz = 0
    !
    ! The two-sided differencing schemes for the first and second derivatives du/dz, d2u/dz2 read:
    !
    ! [2]     du/dz  =  dzeta/dz    (u( k+1) - u( k-1)) / (2 dzeta)
    ! [3]    d2u/dz2 = (dzeta/dz)^2 (u( k+1) + u( k-1) - 2 u( k)) / dzeta^2
    !
    ! A the ice surface, u( k-1) doesn't actually exist, but we can treat it as a "ghost point".
    ! Since, at the ice surface, we know du/dz from the zero-stress boundary condition [1]:
    !
    ! [4]     du/dz = 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)
    !
    ! Substituting [4] into [2] yields:
    !
    !         dzeta/dz (u( k+1) - u( k-1)) / (2 dzeta) = 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)
    !         u( k+1) - u( k-1) = 2 dzeta / (dzeta/dz)  (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))
    ! [5]     u( k-1) = u( k+1) - 2 dzeta / (dzeta/dz)  (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))
    !
    ! Substituting [5] into [3] yields:
    !
    !         d2u/dz2 = (dzeta/dz)^2 1/dzeta^2 ( -2 u( k) + u( k+1) + u( k-1))
    !                 = (dzeta/dz)^2 1/dzeta^2 ( -2 u( k) + 2 u( k+1) - 2 dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)))
    ! [6]             = (dzeta/dz)^2 2/dzeta^2 ( u( k+1) - u( k) - dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)))
    !
    ! The product-rule-expanded form of the BPA reads:
    !
    ! [7]     4 eta d2u/dx2  + 4 deta/dx du/dx + ...
    !           eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !                            deta/dy dv/dx + ...
    !           eta d2u/dz2  +   deta/dz du/dz = - tau_d,x
    !
    ! Substituting [4] and [6] into [7] yields:
    !
    !         4 eta d2u/dx2  + 4 deta/dx du/dx + ...
    !           eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !                            deta/dy dv/dx + ...
    !           eta (dzeta/dz)^2 2/dzeta^2 ( u( k+1) - u( k) - dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))) + ...
    !          deta/dx 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx) = -tau_d,x
    !
    ! Rearranging to group the terms involving du/dx, du/dy, dv/dx, du/dy yields:
    !
    !         4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)) 4 dh/dx + 4 deta/dz dh/dx] + ...
    !         du/dy   [   deta/dy + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz))   dh/dy +   deta/dz dh/dy] + ...
    !         dv/dy   [ 2 deta/dx + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)) 2 dh/dx + 2 deta/dz dh/dx] + ...
    !         dv/dx   [   deta/dy + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)    dh/dy +   deta/dz dh/dy] + ...
    !         u( k+1) [             eta (dzeta/dz)^2 ( 2 / dzeta^2)] + ...
    !         u( k  ) [            -eta (dzeta/dz)^2 ( 2 / dzeta^2)] = -tau_d,x
    !
    ! Rearranging some of the cancelling dzeta/dz and dzeta terms yields:
    !
    ! [8]     4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx - 2 eta / dzeta    dzeta/dz   4 dh/dx + 4 deta/dz dh/dx] + ...
    !         du/dy   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
    !         dv/dy   [ 2 deta/dx - 2 eta / dzeta    dzeta/dz   2 dh/dx + 2 deta/dz dh/dx] + ...
    !         dv/dx   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
    !         u( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
    !         u( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,x

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tikuv

    ! Local variables:
    INTEGER                                                      :: ti, k, uv
    REAL(dp)                                                     :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy, dh_dx, dh_dy, dzeta_dz, dzeta
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: single_row_ind
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddx_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dx2_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dxdy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dy2_val
    INTEGER                                                      :: single_row_nnz
    INTEGER                                                      :: row_tik
    REAL(dp)                                                     :: cu_dudx, cu_dudy, cu_d2udx2, cu_d2udy2, cu_dvdx, cu_dvdy, cu_d2vdxdy, cu_uk, cu_ukp1
    REAL(dp)                                                     :: cv_dvdy, cv_dvdx, cv_d2vdy2, cv_d2vdx2, cv_dudy, cv_dudx, cv_d2udxdy, cv_vk, cv_vkp1
    REAL(dp)                                                     :: Au, Av
    INTEGER                                                      :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)

    ! Safety
    IF (k /= 1) CALL crash('Received k = {int_01}; only applicable at ice surface!', int_01 = k)

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta      = BPA%eta_bks(     ti,1)
    deta_dx  = BPA%deta_dx_bk(  ti,1)
    deta_dy  = BPA%deta_dy_bk(  ti,1)
    deta_dz  = BPA%deta_dz_bk(  ti,1)
    tau_dx   = BPA%tau_dx_b(    ti  )
    tau_dy   = BPA%tau_dy_b(    ti  )
    dh_dx    = BPA%dh_dx_b(     ti  )
    dh_dy    = BPA%dh_dy_b(     ti  )
    dzeta_dz = ice%dzeta_dz_bk( ti,1)
    dzeta    = mesh%zeta( 2) - mesh%zeta( 1)

    ! Calculate coefficients for the different gradients of u and v
    IF (uv == 1) THEN
      ! x-component

      ! 4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
      ! du/dx   [ 4 deta/dx - 2 eta / dzeta    dzeta/dz   4 dh/dx + 4 deta/dz dh/dx] + ...
      ! du/dy   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
      ! dv/dy   [ 2 deta/dx - 2 eta / dzeta    dzeta/dz   2 dh/dx + 2 deta/dz dh/dx] + ...
      ! dv/dx   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
      ! u( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
      ! u( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,x

      cu_dudx    =  4._dp * deta_dx - 2._dp * eta / dzeta * dzeta_dz * 4._dp * dh_dx + 4._dp * deta_dz * dh_dx
      cu_dudy    =          deta_dy - 2._dp * eta / dzeta * dzeta_dz *         dh_dy +         deta_dz * dh_dy
      cu_dvdy    =  2._dp * deta_dx - 2._dp * eta / dzeta * dzeta_dz * 2._dp * dh_dx + 2._dp * deta_dz * dh_dx
      cu_dvdx    =          deta_dy - 2._dp * eta / dzeta * dzeta_dz *         dh_dy +         deta_dz * dh_dy
      cu_d2udx2  =  4._dp * eta
      cu_d2udy2  =          eta
      cu_d2vdxdy =  3._dp * eta
      cu_uk      = -2._dp * eta / dzeta**2 * dzeta_dz**2
      cu_ukp1    =  2._dp * eta / dzeta**2 * dzeta_dz**2

    ELSEIF (uv == 2) THEN
      ! y-component

      ! 4 eta d2v/dy2 + eta d2v/dx2 + 3 eta d2u/dxdy + ...
      ! dv/dy   [ 4 deta/dy - 2 eta / dzeta    dzeta/dz   4 dh/dy + 4 deta/dz dh/dy] + ...
      ! dv/dx   [   deta/dx - 2 eta / dzeta    dzeta/dz     dh/dx +   deta/dz dh/dx] + ...
      ! du/dx   [ 2 deta/dy - 2 eta / dzeta    dzeta/dz   2 dh/dy + 2 deta/dz dh/dy] + ...
      ! du/dy   [   deta/dx - 2 eta / dzeta    dzeta/dz     dh/dx +   deta/dz dh/dx] + ...
      ! v( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
      ! v( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,y

      cv_dvdy    =  4._dp * deta_dy - 2._dp * eta / dzeta * dzeta_dz * 4._dp * dh_dy + 4._dp * deta_dz * dh_dy
      cv_dvdx    =          deta_dx - 2._dp * eta / dzeta * dzeta_dz *         dh_dx +         deta_dz * dh_dx
      cv_dudx    =  2._dp * deta_dy - 2._dp * eta / dzeta * dzeta_dz * 2._dp * dh_dy + 2._dp * deta_dz * dh_dy
      cv_dudy    =          deta_dx - 2._dp * eta / dzeta * dzeta_dz *         dh_dx +         deta_dz * dh_dx
      cv_d2vdy2  =  4._dp * eta
      cv_d2vdx2  =          eta
      cv_d2udxdy =  3._dp * eta
      cv_vk      = -2._dp * eta / dzeta**2 * dzeta_dz**2
      cv_vkp1    =  2._dp * eta / dzeta**2 * dzeta_dz**2

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

    ! Allocate memory for single matrix rows
    ALLOCATE( single_row_ind(        mesh%nC_mem*3*2))
    ALLOCATE( single_row_ddx_val(    mesh%nC_mem*3*2))
    ALLOCATE( single_row_ddy_val(    mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dx2_val(  mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dxdy_val( mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dy2_val(  mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = mesh%tik2n( ti,k)
    CALL read_single_row_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    IF (uv == 1) THEN
      ! x-component

      DO n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        ! Combine coefficients
        Au = cu_dudx    * single_row_ddx_val(    n) + &
             cu_dudy    * single_row_ddy_val(    n) + &
             cu_d2udx2  * single_row_d2dx2_val(  n) + &
             cu_d2udy2  * single_row_d2dy2_val(  n)
        IF (tj == ti .AND. kk == k  ) Au = Au + cu_uk
        IF (tj == ti .AND. kk == k+1) Au = Au + cu_ukp1

        Av = cu_dvdy    * single_row_ddy_val(    n) + &
             cu_dvdx    * single_row_ddx_val(    n) + &
             cu_d2vdxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)

      END DO

      ! Load vector
      bb( row_tikuv) = -tau_dx

    ELSEIF (uv == 2) THEN
      ! y-component

      DO n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)

        ! Combine coefficients
        Av = cv_dvdy    * single_row_ddy_val(    n) + &
             cv_dvdx    * single_row_ddx_val(    n) + &
             cv_d2vdy2  * single_row_d2dy2_val(  n) + &
             cv_d2vdx2  * single_row_d2dx2_val(  n)
        IF (tj == ti .AND. kk == k  ) Av = Av + cv_vk
        IF (tj == ti .AND. kk == k+1) Av = Av + cv_vkp1

        Au = cv_dudx    * single_row_ddx_val(    n) + &
             cv_dudy    * single_row_ddy_val(    n) + &
             cv_d2udxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)

      END DO

      ! Load uector
      bb( row_tikuv) = -tau_dy

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

  END SUBROUTINE calc_BPA_stiffness_matrix_row_BC_surf

  SUBROUTINE calc_BPA_stiffness_matrix_row_BC_base( mesh, ice, BPA, A_CSR, bb, row_tikuv)
    ! Add coefficients to this matrix row to represent the boundary conditions to the BPA at the ice base
    !
    ! At the ice surface (k=1), the zero-stress boundary condition implies that:
    !
    ! [1]     2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) - du/dz + beta_b/eta u = 0
    !
    ! As a shorthand, define:
    !
    ! [2]     P = du/dz = 2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) + beta_b/eta u
    !
    ! The two-sided differencing schemes for the first and second derivatives du/dz, d2u/dz2 read:
    !
    ! [3]     du/dz  =  dzeta/dz    (u( k+1) - u( k-1)          ) / (2 dzeta)
    ! [4]    d2u/dz2 = (dzeta/dz)^2 (u( k+1) + u( k-1) - 2 u( k)) /    dzeta^2
    !
    ! Substituting [3] into [2] and rearranging yields an expression for the ghost node u( k+1):
    !
    !         dzeta/dz (u( k+1) - u( k-1)) / (2 dzeta) = P
    !         u( k+1) - u( k-1) = 2 P dzeta / dzeta_dz
    ! [5]     u( k+1) = u( k-1) + 2 P dzeta / dzeta_dz
    !
    ! We then substitute [5] into [4] to find an expression for d2u/dz2 that no longer depends
    ! on the ghost node u( k+1):
    !
    !         d2u/dz2 = (dzeta/dz)^2 1 / dzeta^2 (u( k+1) + u( k-1) - 2 u( k))
    !                 = (dzeta/dz)^2 1 / dzeta^2 (u( k-1) + 2 P dzeta / dzeta_dz + u( k-1) - 2 u( k))
    !                 = (dzeta/dz)^2 1 / dzeta^2 (2 u( k-1) - 2 u( k) + 2 P dzeta / dzeta_dz)
    ! [6]             = (dzeta/dz)^2 2 / dzeta^2 (u( k-1) - u( k) + P dzeta / dzeta_dz)
    !
    ! The product-rule-expanded form of the BPA reads:
    !
    ! [7]     4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !           eta d2u/dz2  +   deta/dz du/dz = -tau_d,x
    !
    ! Substituting [6] and [2] into [7] yields:
    !
    ! [8]     4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !           eta [(dzeta/dz)^2 2 / dzeta^2 (u( k-1) - u( k) + P dzeta / dzeta_dz)] + ...
    !          deta/dz P = -tau_d,x
    !
    ! As more shorthand, define:
    !
    ! [9]     Q = eta (dzeta/dz)^2 2 / dzeta^2 = 2 eta / dzeta^2 (dzeta/dz)^2
    !
    ! Substituting [9] into [8] yields:
    !
    ! [10]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + P Q dzeta / dzeta_dz)] + deta/dz P = -tau_d,x
    !
    ! As even more shorthand, define:
    !
    ! [11]    R = Q dzeta / (dzeta/dz) + deta/dz = 2 eta / dzeta dzeta/dz + deta/dz
    !
    ! Substituting [11] into [10] yields:
    !
    ! [12]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + R P = -tau_d,x
    !
    ! Substituting [2] back into [12] yields:
    !
    ! [12]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + ...
    !         R [2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) + beta_b/eta u] = -tau_d,x
    !
    ! Rearranging to collect the terms involving du/dx, du/dy, dv/dx, dv/dy yields:
    !
    ! [13]    4 eta d2u/dx2  + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx + 4 R db/dx] + ...
    !         du/dy   [   deta/dy +   R db/dy] + ...
    !         dv/dy   [ 2 deta/dx + 2 R db/dx] + ...
    !         dv/dx   [   deta/dy +   R db/dy] + ...
    !         u( k  ) [R beta_b/eta - Q] + ...
    !         u( k-1) [Q] = -tau_d,x

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tikuv

    ! Local variables:
    INTEGER                                                      :: ti, k, uv
    REAL(dp)                                                     :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy, db_dx, db_dy, dzeta_dz, dzeta, basal_friction_coefficient, Q, R
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: single_row_ind
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddx_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dx2_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dxdy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dy2_val
    INTEGER                                                      :: single_row_nnz
    INTEGER                                                      :: row_tik
    REAL(dp)                                                     :: cu_dudx, cu_dudy, cu_d2udx2, cu_d2udy2, cu_dvdx, cu_dvdy, cu_d2vdxdy, cu_uk, cu_ukm1
    REAL(dp)                                                     :: cv_dvdy, cv_dvdx, cv_d2vdy2, cv_d2vdx2, cv_dudy, cv_dudx, cv_d2udxdy, cv_vk, cv_vkm1
    REAL(dp)                                                     :: Au, Av
    INTEGER                                                      :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)

    ! Safety
    IF (k /= mesh%nz) CALL crash('Received k = {int_01}; only applicable at ice base!', int_01 = k)

    ! Exception for the case of no sliding
    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! u = v = 0
      CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)
      bb( row_tikuv) = 0._dp
      RETURN
    END IF

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta                        = BPA%eta_bks(                      ti,mesh%nz-1)
    deta_dx                    = BPA%deta_dx_bk(                   ti,mesh%nz  )
    deta_dy                    = BPA%deta_dy_bk(                   ti,mesh%nz  )
    deta_dz                    = BPA%deta_dz_bk(                   ti,mesh%nz  )
    tau_dx                     = BPA%tau_dx_b(                     ti          )
    tau_dy                     = BPA%tau_dy_b(                     ti          )
    db_dx                      = BPA%db_dx_b(                      ti          )
    db_dy                      = BPA%db_dy_b(                      ti          )
    basal_friction_coefficient = BPA%basal_friction_coefficient_b( ti          )
    dzeta_dz                   = ice%dzeta_dz_bk(                  ti,mesh%nz  )
    dzeta                      = mesh%zeta( mesh%nz) - mesh%zeta( mesh%nz-1)

    ! [9]     Q = eta (dzeta/dz)^2 2 / dzeta^2 = 2 eta / dzeta^2 (dzeta/dz)^2
    Q = 2._dp * eta / dzeta**2 * dzeta_dz**2

    ! [11]    R = Q dzeta / (dzeta/dz) + deta/dz = 2 eta / dzeta dzeta/dz + deta/dz
    R = 2._dp * eta / dzeta    * dzeta_dz    + deta_dz

    ! Calculate coefficients for the different gradients of u and v
    IF (uv == 1) THEN
      ! x-component

      ! [13]    4 eta d2u/dx2  + eta d2u/dy2 + 3 eta d2v/dxdy + ...
      !         du/dx   [ 4 deta/dx + 4 R db/dx] + ...
      !         du/dy   [   deta/dy +   R db/dy] + ...
      !         dv/dy   [ 2 deta/dx + 2 R db/dx] + ...
      !         dv/dx   [   deta/dy +   R db/dy] + ...
      !         u( k  ) [R beta_b/eta - Q] + ...
      !         u( k-1) [Q] = -tau_d,x

      cu_dudx    =  4._dp * deta_dx + 4._dp * R * db_dx
      cu_dudy    =          deta_dy +         R * db_dy
      cu_dvdy    =  2._dp * deta_dx + 2._dp * R * db_dx
      cu_dvdx    =          deta_dy +         R * db_dy
      cu_d2udx2  =  4._dp * eta
      cu_d2udy2  =          eta
      cu_d2vdxdy =  3._dp * eta
      cu_uk      = ( R * basal_friction_coefficient / eta) - Q
      cu_ukm1    = Q

    ELSEIF (uv == 2) THEN
      ! y-component

      ! [13]    4 eta d2v/dy2  + eta d2v/dx2 + 3 eta d2u/dydx + ...
      !         dv/dy   [ 4 deta/dy + 4 R db/dy] + ...
      !         dv/dx   [   deta/dx +   R db/dx] + ...
      !         du/dx   [ 2 deta/dy + 2 R db/dy] + ...
      !         du/dy   [   deta/dx +   R db/dx] + ...
      !         v( k  ) [R beta_b/eta - Q] + ...
      !         v( k-1) [Q] = -tau_d,y

      cv_dvdy    =  4._dp * deta_dy + 4._dp * R * db_dy
      cv_dvdx    =          deta_dx +         R * db_dx
      cv_dudx    =  2._dp * deta_dy + 2._dp * R * db_dy
      cv_dudy    =          deta_dx +         R * db_dx
      cv_d2vdy2  =  4._dp * eta
      cv_d2vdx2  =          eta
      cv_d2udxdy =  3._dp * eta
      cv_vk      = ( R * basal_friction_coefficient / eta) - Q
      cv_vkm1    = Q

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

    ! Allocate memory for single matrix rows
    ALLOCATE( single_row_ind(        mesh%nC_mem*3*2))
    ALLOCATE( single_row_ddx_val(    mesh%nC_mem*3*2))
    ALLOCATE( single_row_ddy_val(    mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dx2_val(  mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dxdy_val( mesh%nC_mem*3*2))
    ALLOCATE( single_row_d2dy2_val(  mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = mesh%tik2n( ti,k)
    CALL read_single_row_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    IF (uv == 1) THEN
      ! x-component

      DO n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        ! Combine coefficients
        Au = cu_dudx    * single_row_ddx_val(    n) + &
             cu_dudy    * single_row_ddy_val(    n) + &
             cu_d2udx2  * single_row_d2dx2_val(  n) + &
             cu_d2udy2  * single_row_d2dy2_val(  n)
        IF (tj == ti .AND. kk == k  ) Au = Au + cu_uk
        IF (tj == ti .AND. kk == k-1) Au = Au + cu_ukm1

        Av = cu_dvdy    * single_row_ddy_val(    n) + &
             cu_dvdx    * single_row_ddx_val(    n) + &
             cu_d2vdxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)

      END DO

      ! Load vector
      bb( row_tikuv) = -tau_dx

    ELSEIF (uv == 2) THEN
      ! y-component

      DO n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)

        ! Combine coefficients
        Av = cv_dvdy    * single_row_ddy_val(    n) + &
             cv_dvdx    * single_row_ddx_val(    n) + &
             cv_d2vdy2  * single_row_d2dy2_val(  n) + &
             cv_d2vdx2  * single_row_d2dx2_val(  n)
        IF (tj == ti .AND. kk == k  ) Av = Av + cv_vk
        IF (tj == ti .AND. kk == k-1) Av = Av + cv_vkm1

        Au = cv_dudx    * single_row_ddx_val(    n) + &
             cv_dudy    * single_row_ddy_val(    n) + &
             cv_d2udxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)

      END DO

      ! Load uector
      bb( row_tikuv) = -tau_dy

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

  END SUBROUTINE calc_BPA_stiffness_matrix_row_BC_base

  SUBROUTINE calc_BPA_stiffness_matrix_row_BC_west( mesh, BPA, A_CSR, bb, row_tikuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! western domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tikuv

    ! Local variables:
    INTEGER                                                      :: ti,k,uv,row_ti
    INTEGER                                                      :: tj, col_tjkuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_west == 'infinite') THEN
        ! du/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_west == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_west == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      ELSE
        CALL crash('unknown BC_u_west "' // TRIM( C%BC_u_west) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_west == 'infinite') THEN
        ! dv/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_west == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_west == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_west "' // TRIM( C%BC_u_west) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_BPA_stiffness_matrix_row_BC_west

  SUBROUTINE calc_BPA_stiffness_matrix_row_BC_east( mesh, BPA, A_CSR, bb, row_tikuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! eastern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tikuv

    ! Local variables:
    INTEGER                                                      :: ti,k,uv,row_ti
    INTEGER                                                      :: tj, col_tjkuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_east == 'infinite') THEN
        ! du/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_east == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_east == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      ELSE
        CALL crash('unknown BC_u_east "' // TRIM( C%BC_u_east) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_east == 'infinite') THEN
        ! dv/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_east == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_east == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_east "' // TRIM( C%BC_u_east) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_BPA_stiffness_matrix_row_BC_east

  SUBROUTINE calc_BPA_stiffness_matrix_row_BC_south( mesh, BPA, A_CSR, bb, row_tikuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! southern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tikuv

    ! Local variables:
    INTEGER                                                      :: ti,k,uv,row_ti
    INTEGER                                                      :: tj, col_tjkuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_south == 'infinite') THEN
        ! du/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_south == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_south == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      ELSE
        CALL crash('unknown BC_u_south "' // TRIM( C%BC_u_south) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_south == 'infinite') THEN
        ! dv/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_south == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_south == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_south "' // TRIM( C%BC_u_south) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_BPA_stiffness_matrix_row_BC_south

  SUBROUTINE calc_BPA_stiffness_matrix_row_BC_north( mesh, BPA, A_CSR, bb, row_tikuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! northern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tikuv

    ! Local variables:
    INTEGER                                                      :: ti,k,uv,row_ti
    INTEGER                                                      :: tj, col_tjkuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    IF (uv == 1) THEN
      ! x-component

      IF     (C%BC_u_north == 'infinite') THEN
        ! du/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_north == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_u_north == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      ELSE
        CALL crash('unknown BC_u_north "' // TRIM( C%BC_u_north) // '"!')
      END IF

    ELSEIF (uv == 2) THEN
      ! y-component

      IF     (C%BC_v_north == 'infinite') THEN
        ! dv/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        DO n = 1, 3
          tj = mesh%TriC( ti,n)
          IF (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_north == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      ELSEIF (C%BC_v_north == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_north "' // TRIM( C%BC_u_north) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_BPA_stiffness_matrix_row_BC_north

! == Calculate several intermediate terms in the BPA

  SUBROUTINE calc_driving_stress( mesh, ice, BPA)
    ! Calculate the driving stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_driving_stress'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate dh/dx, dh/dy, db/dx, db/dy on the b-grid
    CALL ddx_a_b_2D( mesh, ice%Hs , BPA%dh_dx_b)
    CALL ddy_a_b_2D( mesh, ice%Hs , BPA%dh_dy_b)
    CALL ddx_a_b_2D( mesh, ice%Hib, BPA%db_dx_b)
    CALL ddy_a_b_2D( mesh, ice%Hib, BPA%db_dy_b)

    ! Calculate the driving stress
    DO ti = mesh%ti1, mesh%ti2
      BPA%tau_dx_b( ti) = -ice_density * grav * BPA%dh_dx_b( ti)
      BPA%tau_dy_b( ti) = -ice_density * grav * BPA%dh_dy_b( ti)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress

  SUBROUTINE calc_strain_rates( mesh, BPA)
    ! Calculate the strain rates
    !
    ! The velocities u and v are defined on the bk-grid (triangles, regular vertical)
    !
    ! The horizontal stretch/shear strain rates du/dx, du/dy, dv/dx, dv/dy are
    ! calculated on the ak-grid (vertices, regular vertical), and are then mapped
    ! to the bks-grid (triangles, staggered vertical)
    !
    ! The vertical shear strain rates du/dz, dv/dz are calculated on the bks-grid
    ! (triangles, staggered vertical), and are then mapped to the ak-grid (vertices,
    ! regular vertical).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_strain_rates'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate horizontal stretch/strain rates on the ak-grid
    CALL calc_3D_gradient_bk_ak(  mesh, mesh%M_ddx_bk_ak , BPA%u_bk, BPA%du_dx_ak )
    CALL calc_3D_gradient_bk_ak(  mesh, mesh%M_ddy_bk_ak , BPA%u_bk, BPA%du_dy_ak )
    CALL calc_3D_gradient_bk_ak(  mesh, mesh%M_ddx_bk_ak , BPA%v_bk, BPA%dv_dx_ak )
    CALL calc_3D_gradient_bk_ak(  mesh, mesh%M_ddy_bk_ak , BPA%v_bk, BPA%dv_dy_ak )

    ! Calculate vertical shear strain rates on the bks-grid
    CALL calc_3D_gradient_bk_bks( mesh, mesh%M_ddz_bk_bks, BPA%u_bk, BPA%du_dz_bks)
    CALL calc_3D_gradient_bk_bks( mesh, mesh%M_ddz_bk_bks, BPA%v_bk, BPA%dv_dz_bks)

    ! Map horizontal stretch/shear strain rates from the ak-grid to the bks-grid
    CALL map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%du_dx_ak, BPA%du_dx_bks)
    CALL map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%du_dy_ak, BPA%du_dy_bks)
    CALL map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%dv_dx_ak, BPA%dv_dx_bks)
    CALL map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%dv_dy_ak, BPA%dv_dy_bks)

    ! Map vertical shear strain rates from the bks-grid to the ak-grid
    CALL map_bks_ak( mesh, mesh%M_map_bks_ak, BPA%du_dz_bks, BPA%du_dz_ak)
    CALL map_bks_ak( mesh, mesh%M_map_bks_ak, BPA%dv_dz_bks, BPA%dv_dz_ak)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_rates

  SUBROUTINE calc_effective_viscosity( mesh, ice, BPA, Glens_flow_law_epsilon_sq_0_applied)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N
    !
    ! The effective viscosity eta is calculated separately on both the ak-grid (vertices, regular vertical)
    ! and on the bks-grid (triangles, staggered vertical), using the strain rates calculated in calc_strain_rates.
    !
    ! eta_bk, deta_dx_bk, and deta_dy_bk are calculated from eta_ak

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA
    REAL(dp),                            INTENT(IN)              :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_viscosity'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      ::  A_flow_bks
    INTEGER                                                      :: vi,ti,k,ks
    REAL(dp)                                                     :: A_min, eta_max
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      :: eta_bk_from_ak, eta_bk_from_bks
    REAL(dp)                                                     :: uabs_base, uabs_surf, R_shear

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! Allocate memory
    ALLOCATE( A_flow_bks( mesh%ti1:mesh%ti2, mesh%nz-1))

  ! == Calculate effective viscosity on the ak-grid

    ! Calculate the effective viscosity eta
    IF (C%choice_flow_law == 'Glen') THEN
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      CALL calc_ice_rheology_Glen( mesh, ice)

      ! Calculate effective viscosity
      DO vi = mesh%vi1, mesh%vi2
      DO k  = 1, mesh%nz
        BPA%eta_ak( vi,k) = calc_effective_viscosity_Glen_3D_uv_only( Glens_flow_law_epsilon_sq_0_applied, &
          BPA%du_dx_ak( vi,k), BPA%du_dy_ak( vi,k), BPA%du_dz_ak( vi,k), &
          BPA%dv_dx_ak( vi,k), BPA%dv_dy_ak( vi,k), BPA%dv_dz_ak( vi,k), ice%A_flow( vi,k))
      END DO
      END DO

    ELSE
      CALL crash('unknown choice_flow_law "' // TRIM( C%choice_flow_law) // '"!')
    END IF

    ! Safety
    BPA%eta_ak = MIN( MAX( BPA%eta_ak, C%visc_eff_min), eta_max)

  ! == Calculate effective viscosity on the bks-grid

    ! Calculate the effective viscosity eta
    IF (C%choice_flow_law == 'Glen') THEN
      ! Calculate the effective viscosity according to Glen's flow law

      ! Calculate flow factors: map ice flow factor from the ak-grid to the bks-grid
      CALL map_ak_bks( mesh, mesh%M_map_ak_bks, ice%A_flow, A_flow_bks)

      ! Calculate effective viscosity
      DO ti = mesh%ti1, mesh%ti2
      DO ks  = 1, mesh%nz-1
        BPA%eta_bks( ti,ks) = calc_effective_viscosity_Glen_3D_uv_only( C%Glens_flow_law_epsilon_sq_0, &
          BPA%du_dx_bks( ti,ks), BPA%du_dy_bks( ti,ks), BPA%du_dz_bks( ti,ks), &
          BPA%dv_dx_bks( ti,ks), BPA%dv_dy_bks( ti,ks), BPA%dv_dz_bks( ti,ks), A_flow_bks( ti,ks))
      END DO
      END DO

    ELSE
      CALL crash('unknown choice_flow_law "' // TRIM( C%choice_flow_law) // '"!')
    END IF

    ! Safety
    BPA%eta_bks = MIN( MAX( BPA%eta_bks, C%visc_eff_min), eta_max)

    ! Calculate the horizontal gradients of the effective viscosity from its value on the ak-grid
    CALL calc_3D_gradient_ak_bk(  mesh, mesh%M_ddx_ak_bk , BPA%eta_ak , BPA%deta_dx_bk)
    CALL calc_3D_gradient_ak_bk(  mesh, mesh%M_ddy_ak_bk , BPA%eta_ak , BPA%deta_dy_bk)

    ! Calculate the vertical gradients of the effective viscosity from its value on the bks-grid
    CALL calc_3D_gradient_bks_bk( mesh, mesh%M_ddz_bks_bk, BPA%eta_bks, BPA%deta_dz_bk)

    ! Map the effective viscosity from the ak- and bks-grids to the bk-grid

    ALLOCATE( eta_bk_from_ak(  mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( eta_bk_from_bks( mesh%ti1:mesh%ti2,mesh%nz))

    CALL map_a_b_3D( mesh, BPA%eta_ak, eta_bk_from_ak)
    CALL calc_3D_gradient_bks_bk( mesh, mesh%M_map_bks_bk, BPA%eta_bks, eta_bk_from_bks)

    ! Preliminary experiments suggest that in settings where ice flow is dominated by
    ! vertical shear (e.g. the Halfar dome with no sliding), the solver is only stable
    ! when using eta_bk_from_bks. But in settings with a lot of sliding and little
    ! vertical shear (e.g. ISMIP-HOM C), it needs eta_from_ak instead. The "shear factor"
    ! R_shear serves to provide a crude approximation to which flow mode dominates,
    ! which is then use to calculate a weighted average between the two versions of eta.

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate the shear factor R_shear
      uabs_surf = SQRT( 0.1_dp + BPA%u_bk( ti,1      )**2 + BPA%v_bk( ti,1      )**2)
      uabs_base = SQRT( 0.1_dp + BPA%u_bk( ti,mesh%nz)**2 + BPA%v_bk( ti,mesh%nz)**2)
      R_shear = uabs_base / uabs_surf

      ! By the nature of ice flow, uabs_base <= uabs_surf, so 0 <= R_shear <= 1,
      ! with 0 indicating no sliding and therefore full vertical shear, and
      ! 1 indicating full sliding.

      ! Weighted average
      BPA%eta_bk( ti,:) = R_shear * eta_bk_from_ak( ti,:) + (1._dp - R_shear) * eta_bk_from_bks( ti,:)

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    DEALLOCATE( A_flow_bks)
    DEALLOCATE( eta_bk_from_ak)
    DEALLOCATE( eta_bk_from_bks)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity

  SUBROUTINE calc_applied_basal_friction_coefficient( mesh, ice, BPA)
    ! Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    ! and scaled with the sub-grid grounded fraction
    !
    ! This is where the sliding law is called!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_applied_basal_friction_coefficient'
    INTEGER                                                      :: ti
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: u_base_b, v_base_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( u_base_b( mesh%ti1:mesh%ti2))
    ALLOCATE( v_base_b( mesh%ti1:mesh%ti2))

    ! Copy basal velocities from 3-D fields
    u_base_b = BPA%u_bk( :,mesh%nz)
    v_base_b = BPA%v_bk( :,mesh%nz)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!
    CALL calc_basal_friction_coefficient( mesh, ice, u_base_b, v_base_b)

    ! Map basal friction coefficient beta_b to the b-grid
    CALL map_a_b_2D( mesh, ice%basal_friction_coefficient, BPA%basal_friction_coefficient_b)

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    IF (C%do_GL_subgrid_friction) THEN
      DO ti = mesh%ti1, mesh%ti2
        BPA%basal_friction_coefficient_b( ti) = BPA%basal_friction_coefficient_b( ti) * ice%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      END DO
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_applied_basal_friction_coefficient

! == Some useful tools for improving numerical stability of the viscosity iteration

  SUBROUTINE relax_viscosity_iterations( mesh, BPA, visc_it_relax)
    ! Reduce the change between velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA
    REAL(dp),                            INTENT(IN)              :: visc_it_relax

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'relax_viscosity_iterations'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      BPA%u_bk( ti,:) = (visc_it_relax * BPA%u_bk( ti,:)) + ((1._dp - visc_it_relax) * BPA%u_bk_prev( ti,:))
      BPA%v_bk( ti,:) = (visc_it_relax * BPA%v_bk( ti,:)) + ((1._dp - visc_it_relax) * BPA%v_bk_prev( ti,:))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_viscosity_iterations

  SUBROUTINE calc_visc_iter_UV_resid( mesh, BPA, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    REAL(dp),                            INTENT(OUT)             :: resid_UV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_visc_iter_UV_resid'
    INTEGER                                                      :: ierr
    INTEGER                                                      :: ti,k
    REAL(dp)                                                     :: res1, res2

    ! Add routine to path
    CALL init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, mesh%nz

      res1 = res1 + (BPA%u_bk( ti,k) - BPA%u_bk_prev( ti,k))**2
      res1 = res1 + (BPA%v_bk( ti,k) - BPA%v_bk_prev( ti,k))**2

      res2 = res2 + (BPA%u_bk( ti,k) + BPA%u_bk_prev( ti,k))**2
      res2 = res2 + (BPA%v_bk( ti,k) + BPA%v_bk_prev( ti,k))**2

    END DO
    END DO

    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_visc_iter_UV_resid

  SUBROUTINE apply_velocity_limits( mesh, BPA)
    ! Limit velocities for improved stability

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'apply_velocity_limits'
    INTEGER                                                      :: ti,k
    REAL(dp)                                                     :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, mesh%nz

      ! Calculate absolute speed
      uabs = SQRT( BPA%u_bk( ti,k)**2 + BPA%v_bk( ti,k)**2)

      ! Reduce velocities if neceBPAry
      IF (uabs > C%vel_max) THEN
        BPA%u_bk( ti,k) = BPA%u_bk( ti,k) * C%vel_max / uabs
        BPA%v_bk( ti,k) = BPA%v_bk( ti,k) * C%vel_max / uabs
      END IF

    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits

! == Initialisation

  SUBROUTINE initialise_BPA_velocities_from_file( mesh, BPA, region_name)
    ! Initialise the velocities for the BPA solver from an external NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT) :: BPA
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BPA_velocities_from_file'
    REAL(dp)                                           :: dummy1
    CHARACTER(LEN=256)                                 :: filename
    REAL(dp)                                           :: timeframe

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy1 = mesh%xmin

    ! Determine the filename and timeframe to read for this model region
    IF     (region_name == 'NAM') THEN
      filename  = C%filename_initial_velocity_NAM
      timeframe = C%timeframe_initial_velocity_NAM
    ELSEIF (region_name == 'EAS') THEN
      filename  = C%filename_initial_velocity_EAS
      timeframe = C%timeframe_initial_velocity_EAS
    ELSEIF (region_name == 'GRL') THEN
      filename  = C%filename_initial_velocity_GRL
      timeframe = C%timeframe_initial_velocity_GRL
    ELSEIF (region_name == 'ANT') THEN
      filename  = C%filename_initial_velocity_ANT
      timeframe = C%timeframe_initial_velocity_ANT
    ELSE
      CALL crash('unknown model region "' // region_name // '"!')
    END IF

    ! Write to terminal
    IF (par%master) WRITE(0,*) '   Initialising BPA velocities from file "' // colour_string( TRIM( filename),'light blue') // '"...'

    ! Read velocities from the file
    IF (timeframe == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_mesh_file_3D_b( filename, 'u_bk', BPA%u_bk)
      CALL read_field_from_mesh_file_3D_b( filename, 'v_bk', BPA%v_bk)
    ELSE
      ! Read specified timeframe
      CALL read_field_from_mesh_file_3D_b( filename, 'u_bk', BPA%u_bk, time_to_read = timeframe)
      CALL read_field_from_mesh_file_3D_b( filename, 'v_bk', BPA%v_bk, time_to_read = timeframe)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BPA_velocities_from_file

  SUBROUTINE allocate_BPA_solver( mesh, BPA)
    ! Allocate memory the BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(OUT)   :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_BPA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Solution
    ALLOCATE( BPA%u_bk(          mesh%ti1:mesh%ti2,mesh%nz))                   ! [m yr^-1] 2-D horizontal ice velocity
    BPA%u_bk = 0._dp
    ALLOCATE( BPA%v_bk(          mesh%ti1:mesh%ti2,mesh%nz))
    BPA%v_bk = 0._dp

    ! Intermediate data fields
    ALLOCATE( BPA%du_dx_ak(      mesh%vi1:mesh%vi2,mesh%nz))                   ! [yr^-1] 2-D horizontal strain rates
    BPA%du_dx_ak = 0._dp
    ALLOCATE( BPA%du_dy_ak(      mesh%vi1:mesh%vi2,mesh%nz))
    BPA%du_dy_ak = 0._dp
    ALLOCATE( BPA%du_dz_ak(      mesh%vi1:mesh%vi2,mesh%nz))
    BPA%du_dz_ak = 0._dp
    ALLOCATE( BPA%dv_dx_ak(      mesh%vi1:mesh%vi2,mesh%nz))                   ! [yr^-1] 2-D horizontal strain rates
    BPA%dv_dx_ak = 0._dp
    ALLOCATE( BPA%dv_dy_ak(      mesh%vi1:mesh%vi2,mesh%nz))
    BPA%dv_dy_ak = 0._dp
    ALLOCATE( BPA%dv_dz_ak(      mesh%vi1:mesh%vi2,mesh%nz))
    BPA%dv_dz_ak = 0._dp
    ALLOCATE( BPA%du_dx_bks(     mesh%ti1:mesh%ti2,mesh%nz-1))                 ! [yr^-1] 2-D horizontal strain rates
    BPA%du_dx_bks = 0._dp
    ALLOCATE( BPA%du_dy_bks(     mesh%ti1:mesh%ti2,mesh%nz-1))
    BPA%du_dy_bks = 0._dp
    ALLOCATE( BPA%du_dz_bks(     mesh%ti1:mesh%ti2,mesh%nz-1))
    BPA%du_dz_bks = 0._dp
    ALLOCATE( BPA%dv_dx_bks(     mesh%ti1:mesh%ti2,mesh%nz-1))                 ! [yr^-1] 2-D horizontal strain rates
    BPA%dv_dx_bks = 0._dp
    ALLOCATE( BPA%dv_dy_bks(     mesh%ti1:mesh%ti2,mesh%nz-1))
    BPA%dv_dy_bks = 0._dp
    ALLOCATE( BPA%dv_dz_bks(     mesh%ti1:mesh%ti2,mesh%nz-1))
    BPA%dv_dz_bks = 0._dp
    ALLOCATE( BPA%eta_ak(        mesh%vi1:mesh%vi2,mesh%nz))                   ! Effective viscosity
    BPA%eta_ak = 0._dp
    ALLOCATE( BPA%eta_bks(       mesh%ti1:mesh%ti2,mesh%nz-1))
    BPA%eta_bks = 0._dp
    ALLOCATE( BPA%eta_bk(        mesh%ti1:mesh%ti2,mesh%nz))
    BPA%eta_bk = 0._dp
    ALLOCATE( BPA%deta_dx_bk(    mesh%ti1:mesh%ti2,mesh%nz))                   ! Gradients of eta
    BPA%deta_dx_bk = 0._dp
    ALLOCATE( BPA%deta_dy_bk(    mesh%ti1:mesh%ti2,mesh%nz))                   ! Gradients of eta
    BPA%deta_dy_bk = 0._dp
    ALLOCATE( BPA%deta_dz_bk(    mesh%ti1:mesh%ti2,mesh%nz))                   ! Gradients of eta
    BPA%deta_dz_bk = 0._dp
    ALLOCATE( BPA%basal_friction_coefficient_b(      mesh%ti1:mesh%ti2))       ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    BPA%basal_friction_coefficient_b = 0._dp
    ALLOCATE( BPA%dh_dx_b(       mesh%ti1:mesh%ti2))                           ! Surface slope
    BPA%dh_dx_b = 0._dp
    ALLOCATE( BPA%dh_dy_b(       mesh%ti1:mesh%ti2))
    BPA%dh_dy_b = 0._dp
    ALLOCATE( BPA%db_dx_b(       mesh%ti1:mesh%ti2))                           ! Basal slope
    BPA%db_dx_b = 0._dp
    ALLOCATE( BPA%db_dy_b(       mesh%ti1:mesh%ti2))
    BPA%db_dy_b = 0._dp
    ALLOCATE( BPA%tau_dx_b(      mesh%ti1:mesh%ti2))                           ! Driving stress
    BPA%tau_dx_b = 0._dp
    ALLOCATE( BPA%tau_dy_b(      mesh%ti1:mesh%ti2))
    BPA%tau_dy_b = 0._dp
    ALLOCATE( BPA%u_bk_prev(     mesh%nTri,mesh%nz))                   ! Velocity solution from previous viscosity iteration
    BPA%u_bk_prev = 0._dp
    ALLOCATE( BPA%v_bk_prev(     mesh%nTri,mesh%nz))
    BPA%v_bk_prev = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_BPA_solver

! == Restart NetCDF files

  SUBROUTINE write_to_restart_file_BPA( mesh, BPA, time)
    ! Write to the restart NetCDF file for the BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)              :: BPA
    REAL(dp),                            INTENT(IN)              :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'write_to_restart_file_BPA'
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Writing to BPA restart file "' // &
      colour_string( TRIM( BPA%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( BPA%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( BPA%restart_filename, ncid, time)

    ! Write the velocity fields to the file
    CALL write_to_field_multopt_mesh_dp_3D_b( mesh, BPA%restart_filename, ncid, 'u_bk', BPA%u_bk)
    CALL write_to_field_multopt_mesh_dp_3D_b( mesh, BPA%restart_filename, ncid, 'v_bk', BPA%v_bk)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_BPA

  SUBROUTINE create_restart_file_BPA( mesh, BPA)
    ! Create a restart NetCDF file for the BPA solver
    ! Includes generation of the procedural filename (e.g. "restart_BPA_00001.nc")

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'create_restart_file_BPA'
    CHARACTER(LEN=256)                                           :: filename_base
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_ice_velocity_BPA'
    CALL generate_filename_XXXXXdotnc( filename_base, BPA%restart_filename)

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Creating BPA restart file "' // &
      colour_string( TRIM( BPA%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( BPA%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( BPA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( BPA%restart_filename, ncid)

    ! Add a zeta dimension to the file
    CALL add_zeta_dimension_to_file( BPA%restart_filename, ncid, mesh%zeta)

    ! Add the velocity fields to the file
    CALL add_field_mesh_dp_3D_b( BPA%restart_filename, ncid, 'u_bk', long_name = '3-D horizontal ice velocity in the x-direction', units = 'm/yr')
    CALL add_field_mesh_dp_3D_b( BPA%restart_filename, ncid, 'v_bk', long_name = '3-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_BPA

END MODULE ice_velocity_BPA
