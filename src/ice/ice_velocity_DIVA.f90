MODULE ice_velocity_DIVA

  ! Routines for calculating ice velocities using the Depth-Integrated Viscosity Approximation (DIVA)

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF, &
                                                                     save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D
  USE parameters
  USE petsc_basic                                            , ONLY: solve_matrix_equation_CSR_PETSc
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_DIVA
  USE mesh_operators                                         , ONLY: map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D, map_b_a_2D, map_b_a_3D, map_a_b_3D
  USE mesh_zeta                                              , ONLY: vertical_average, integrate_from_zeta_is_one_to_zeta_is_zetap
  USE sliding_laws                                           , ONLY: calc_basal_friction_coefficient
  USE mesh_utilities                                         , ONLY: find_ti_copy_ISMIP_HOM_periodic
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file, open_existing_netcdf_file_for_writing
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, add_time_dimension_to_file, &
                                                                     add_field_mesh_dp_2D_b, add_field_mesh_dp_3D_b, write_time_to_file, write_to_field_multopt_mesh_dp_2D_b, &
                                                                     write_to_field_multopt_mesh_dp_3D_b, add_zeta_dimension_to_file
  USE netcdf_input                                           , ONLY: read_field_from_mesh_file_2D_b, read_field_from_mesh_file_3D_b
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D
  USE ice_flow_laws                                          , ONLY: calc_effective_viscosity_Glen_3D_uv_only, calc_ice_rheology_Glen
  USE reallocate_mod                                         , ONLY: reallocate_bounds, reallocate_clean
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Main routines

  SUBROUTINE initialise_DIVA_solver( mesh, DIVA, region_name)
    ! Initialise the DIVA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(OUT)   :: DIVA
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_DIVA_solver'
    CHARACTER(LEN=256)                                 :: choice_initial_velocity

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    CALL allocate_DIVA_solver( mesh, DIVA)

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
      DIVA%u_vav_b = 0._dp
      DIVA%v_vav_b = 0._dp
    ELSEIF (choice_initial_velocity == 'read_from_file') THEN
      CALL initialise_DIVA_velocities_from_file( mesh, DIVA, region_name)
    ELSE
      CALL crash('unknown choice_initial_velocity "' // TRIM( choice_initial_velocity) // '"!')
    END IF

    ! Set tolerances for PETSc matrix solver for the linearised DIVA
    DIVA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    DIVA%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_DIVA_solver

  SUBROUTINE solve_DIVA( mesh, ice, DIVA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Calculate ice velocities by solving the Depth-Integrated Viscosity Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA
    INTEGER,  DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_DIVA'
    LOGICAL                                                      :: grounded_ice_exists
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: BC_prescr_mask_b_applied
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: BC_prescr_u_b_applied
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: BC_prescr_v_b_applied
    INTEGER                                                      :: viscosity_iteration_i
    LOGICAL                                                      :: has_converged
    REAL(dp)                                                     :: resid_UV, resid_UV_prev
    REAL(dp)                                                     :: uv_min, uv_max
    REAL(dp)                                                     :: visc_it_relax_applied
    REAL(dp)                                                     :: Glens_flow_law_epsilon_sq_0_applied
    INTEGER                                                      :: nit_diverg_consec

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there is no grounded ice, no need (in fact, no way) to solve the DIVA
    grounded_ice_exists = ANY( ice%mask_grounded_ice)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. grounded_ice_exists) THEN
      DIVA%u_vav_b  = 0._dp
      DIVA%v_vav_b  = 0._dp
      DIVA%u_base_b = 0._dp
      DIVA%v_base_b = 0._dp
      DIVA%u_3D_b   = 0._dp
      DIVA%v_3D_b   = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Handle the optional prescribed u,v boundary conditions
    ALLOCATE( BC_prescr_mask_b_applied( mesh%ti1:mesh%ti2))
    ALLOCATE( BC_prescr_u_b_applied(    mesh%ti1:mesh%ti2))
    ALLOCATE( BC_prescr_v_b_applied(    mesh%ti1:mesh%ti2))
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

    ! Calculate the driving stress
    CALL calc_driving_stress( mesh, ice, DIVA)

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

      ! Calculate the horizontal strain rates for the current velocity solution
      CALL calc_horizontal_strain_rates( mesh, DIVA)

      ! Calculate the vertical shear strain rates
      CALL calc_vertical_shear_strain_rates( mesh, DIVA)

      ! Calculate the effective viscosity for the current velocity solution
      CALL calc_effective_viscosity( mesh, ice, DIVA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the F-integrals (Lipscomb et al. (2019), Eq. 30)
      CALL calc_F_integrals( mesh, ice, DIVA)

      ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)
      CALL calc_effective_basal_friction_coefficient( mesh, ice, DIVA)

      ! Solve the linearised DIVA to calculate a new velocity solution
      CALL solve_DIVA_linearised( mesh, DIVA, BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Limit velocities for improved stability
      CALL apply_velocity_limits( mesh, DIVA)

      ! Reduce the change between velocity solutions
      CALL relax_viscosity_iterations( mesh, DIVA, visc_it_relax_applied)

      ! Calculate basal velocities
      CALL calc_basal_velocities( mesh, DIVA)

      ! Calculate basal shear stress
      CALL calc_basal_shear_stress( mesh, DIVA)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      CALL calc_visc_iter_UV_resid( mesh, DIVA, resid_UV)

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
      uv_min = MINVAL( DIVA%u_vav_b)
      uv_max = MAXVAL( DIVA%u_vav_b)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      ! IF (par%master) WRITE(0,*) '    DIVA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

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

    ! Calculate 3-D ice velocities
    CALL calc_3D_velocities( mesh, DIVA)

    ! Clean up after yourself
    DEALLOCATE( BC_prescr_mask_b_applied)
    DEALLOCATE( BC_prescr_u_b_applied   )
    DEALLOCATE( BC_prescr_v_b_applied   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_DIVA

  SUBROUTINE remap_DIVA_solver( mesh_old, mesh_new, DIVA)
    ! Remap the DIVA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT) :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_DIVA_solver'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: u_vav_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: v_vav_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: tau_bx_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: tau_by_a
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: eta_3D_a
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: u_3D_a
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: v_3D_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! Allocate memory for velocities on the a-grid (vertices)
    ALLOCATE( u_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    ALLOCATE( v_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    ALLOCATE( tau_bx_a( mesh_old%vi1: mesh_old%vi2             ))
    ALLOCATE( tau_by_a( mesh_old%vi1: mesh_old%vi2             ))
    ALLOCATE( eta_3D_a( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    ALLOCATE( u_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    ALLOCATE( v_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    CALL map_b_a_2D( mesh_old, DIVA%u_vav_b , u_vav_a )
    CALL map_b_a_2D( mesh_old, DIVA%v_vav_b , v_vav_a )
    CALL map_b_a_2D( mesh_old, DIVA%tau_bx_b, tau_bx_a)
    CALL map_b_a_2D( mesh_old, DIVA%tau_by_b, tau_by_a)
    CALL map_b_a_3D( mesh_old, DIVA%eta_3D_b, eta_3D_a)
    CALL map_b_a_3D( mesh_old, DIVA%u_3D_b  , u_3D_a  )
    CALL map_b_a_3D( mesh_old, DIVA%v_3D_b  , v_3D_a  )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, u_vav_a , '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, v_vav_a , '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, tau_bx_a, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, tau_by_a, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, eta_3D_a, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, u_3D_a  , '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, v_3D_a  , '2nd_order_conservative')

    ! Reallocate memory for the data on the triangles
    CALL reallocate_bounds( DIVA%u_vav_b                     , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%v_vav_b                     , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%tau_bx_b                    , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%tau_by_b                    , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%eta_3D_b                    , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%u_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%v_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    CALL map_a_b_2D( mesh_new, u_vav_a , DIVA%u_vav_b )
    CALL map_a_b_2D( mesh_new, v_vav_a , DIVA%v_vav_b )
    CALL map_a_b_2D( mesh_new, tau_bx_a, DIVA%tau_bx_b)
    CALL map_a_b_2D( mesh_new, tau_by_a, DIVA%tau_by_b)
    CALL map_a_b_3D( mesh_new, eta_3D_a, DIVA%eta_3D_b)
    CALL map_a_b_3D( mesh_new, u_3D_a  , DIVA%u_3D_b  )
    CALL map_a_b_3D( mesh_new, v_3D_a  , DIVA%v_3D_b  )

    ! Clean up after yourself
    DEALLOCATE( u_vav_a )
    DEALLOCATE( v_vav_a )
    DEALLOCATE( tau_bx_a)
    DEALLOCATE( tau_by_a)
    DEALLOCATE( eta_3D_a)
    DEALLOCATE( u_3D_a  )
    DEALLOCATE( v_3D_a  )

    ! Reallocate everything else
    ! ==========================

    ! CALL reallocate_bounds( DIVA%u_vav_b                     , mesh_new%ti1, mesh_new%ti2             )           ! [m yr^-1] 2-D vertically averaged horizontal ice velocity
    ! CALL reallocate_bounds( DIVA%v_vav_b                     , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%u_base_b                    , mesh_new%ti1, mesh_new%ti2             )           ! [m yr^-1] 2-D horizontal ice velocity at the ice base
    CALL reallocate_bounds( DIVA%v_base_b                    , mesh_new%ti1, mesh_new%ti2             )
    ! CALL reallocate_bounds( DIVA%u_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)           ! [m yr^-1] 3-D horizontal ice velocity
    ! CALL reallocate_bounds( DIVA%v_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%du_dx_a                     , mesh_new%vi1, mesh_new%vi2             )           ! [yr^-1] 2-D horizontal strain rates
    CALL reallocate_bounds( DIVA%du_dy_a                     , mesh_new%vi1, mesh_new%vi2             )
    CALL reallocate_bounds( DIVA%dv_dx_a                     , mesh_new%vi1, mesh_new%vi2             )
    CALL reallocate_bounds( DIVA%dv_dy_a                     , mesh_new%vi1, mesh_new%vi2             )
    CALL reallocate_bounds( DIVA%du_dz_3D_a                  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! [yr^-1] 3-D vertical shear strain rates
    CALL reallocate_bounds( DIVA%dv_dz_3D_a                  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%eta_3D_a                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! Effective viscosity
    ! CALL reallocate_bounds( DIVA%eta_3D_b                    , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%eta_vav_a                   , mesh_new%vi1, mesh_new%vi2             )
    CALL reallocate_bounds( DIVA%N_a                         , mesh_new%vi1, mesh_new%vi2             )           ! Product term N = eta * H
    CALL reallocate_bounds( DIVA%N_b                         , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%dN_dx_b                     , mesh_new%ti1, mesh_new%ti2             )           ! Gradients of N
    CALL reallocate_bounds( DIVA%dN_dy_b                     , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%F1_3D_a                     , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! F-integrals
    CALL reallocate_bounds( DIVA%F2_3D_a                     , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%F1_3D_b                     , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%F2_3D_b                     , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( DIVA%basal_friction_coefficient_b, mesh_new%ti1, mesh_new%ti2             )           ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    CALL reallocate_bounds( DIVA%beta_eff_a                  , mesh_new%vi1, mesh_new%vi2             )           ! "Effective" friction coefficient (turning the SSA into the DIVA)
    CALL reallocate_bounds( DIVA%beta_eff_b                  , mesh_new%ti1, mesh_new%ti2             )
    ! CALL reallocate_bounds( DIVA%tau_bx_b                    , mesh_new%ti1, mesh_new%ti2             )           ! Basal shear stress
    ! CALL reallocate_bounds( DIVA%tau_by_b                    , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( DIVA%tau_dx_b                    , mesh_new%ti1, mesh_new%ti2             )           ! Driving stress
    CALL reallocate_bounds( DIVA%tau_dy_b                    , mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_clean ( DIVA%u_b_prev                    , mesh_new%nTri                          )           ! Velocity solution from previous viscosity iteration
    CALL reallocate_clean ( DIVA%v_b_prev                    , mesh_new%nTri                          )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_DIVA_solver

! == Assemble and solve the linearised DIVA

  SUBROUTINE solve_DIVA_linearised( mesh, DIVA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Solve the linearised DIVA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA
    INTEGER,  DIMENSION(mesh%ti1:mesh%ti2),          INTENT(IN)  :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2),          INTENT(IN)  :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2),          INTENT(IN)  :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_DIVA_linearised'
    INTEGER                                                      :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                              :: A_CSR
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: bb
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: uv_buv
    INTEGER                                                      :: row_tiuv,ti,uv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Store the previous solution
    CALL gather_to_all_dp_1D( DIVA%u_vav_b, DIVA%u_b_prev)
    CALL gather_to_all_dp_1D( DIVA%v_vav_b, DIVA%v_b_prev)

  ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
  ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * 2      ! from
    ncols_loc       = mesh%nTri_loc * 2
    nrows           = mesh%nTri     * 2      ! to
    nrows_loc       = mesh%nTri_loc * 2
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz * 4

    CALL allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector and the solution
    ALLOCATE( bb(     mesh%ti1*2-1: mesh%ti2*2))
    ALLOCATE( uv_buv( mesh%ti1*2-1: mesh%ti2*2))

    ! Fill in the current velocity solution
    DO ti = mesh%ti1, mesh%ti2

      ! u
      row_tiuv = mesh%tiuv2n( ti,1)
      uv_buv( row_tiuv) = DIVA%u_vav_b( ti)

      ! v
      row_tiuv = mesh%tiuv2n( ti,2)
      uv_buv( row_tiuv) = DIVA%v_vav_b( ti)

    END DO ! DO ti = mesh%ti1, mesh%ti2

  ! == Construct the stiffness matrix for the linearised DIVA
  ! ========================================================

    DO row_tiuv = A_CSR%i1, A_CSR%i2

      ti = mesh%n2tiuv( row_tiuv,1)
      uv = mesh%n2tiuv( row_tiuv,2)

      IF (BC_prescr_mask_b( ti) == 1) THEN
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector: prescribed velocity
        IF     (uv == 1) THEN
          bb( row_tiuv) = BC_prescr_u_b( ti)
        ELSEIF (uv == 2) THEN
          bb( row_tiuv) = BC_prescr_v_b( ti)
        ELSE
          CALL crash('uv can only be 1 or 2!')
        END IF

      ELSEIF (mesh%TriBI( ti) == 1 .OR. mesh%TriBI( ti) == 2) THEN
        ! Northern domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_north( mesh, DIVA, A_CSR, bb, row_tiuv)

      ELSEIF (mesh%TriBI( ti) == 3 .OR. mesh%TriBI( ti) == 4) THEN
        ! Eastern domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_east( mesh, DIVA, A_CSR, bb, row_tiuv)

      ELSEIF (mesh%TriBI( ti) == 5 .OR. mesh%TriBI( ti) == 6) THEN
        ! Southern domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_south( mesh, DIVA, A_CSR, bb, row_tiuv)

      ELSEIF (mesh%TriBI( ti) == 7 .OR. mesh%TriBI( ti) == 8) THEN
        ! Western domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_west( mesh, DIVA, A_CSR, bb, row_tiuv)

      ELSE
        ! No boundary conditions apply; solve the DIVA

        IF (C%do_include_SSADIVA_crossterms) THEN
          ! Calculate matrix coefficients for the full DIVA
          CALL calc_DIVA_stiffness_matrix_row_free( mesh, DIVA, A_CSR, bb, row_tiuv)
        ELSE
          ! Calculate matrix coefficients for the DIVA sans the gradients of the effective viscosity (the "cross-terms")
          CALL calc_DIVA_sans_stiffness_matrix_row_free( mesh, DIVA, A_CSR, bb, row_tiuv)
        END IF

      END IF

    END DO ! DO row_tiuv = A_CSR%i1, A_CSR%i2

  ! == Solve the matrix equation
  ! ============================

    ! Use PETSc to solve the matrix equation
    CALL solve_matrix_equation_CSR_PETSc( A_CSR, bb, uv_buv, DIVA%PETSc_rtol, DIVA%PETSc_abstol)

    ! Disentangle the u and v components of the velocity solution
    DO ti = mesh%ti1, mesh%ti2

      ! u
      row_tiuv = mesh%tiuv2n( ti,1)
      DIVA%u_vav_b( ti) = uv_buv( row_tiuv)

      ! v
      row_tiuv = mesh%tiuv2n( ti,2)
      DIVA%v_vav_b( ti) = uv_buv( row_tiuv)

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( A_CSR)
    DEALLOCATE( bb)
    DEALLOCATE( uv_buv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_DIVA_linearised

  SUBROUTINE calc_DIVA_stiffness_matrix_row_free( mesh, DIVA, A_CSR, bb, row_tiuv)
    ! Add coefficients to this matrix row to represent the linearised DIVA
    !
    ! The DIVA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_eff u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_eff v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_eff u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_eff v = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_eff u + ...
    !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx
    !
    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_eff v + ...
    !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy
    !
    ! We define the velocities u,v, the effective basal friction coefficient beta_eff, and the driving
    ! stress tau_d on the b-grid (triangles), and the effective viscosity eta and the
    ! product term N = eta H on the a-grid (vertices).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(IN)              :: DIVA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(mesh%ti1*2-1: mesh%ti2*2), INTENT(INOUT) :: bb
    INTEGER,                             INTENT(IN)              :: row_tiuv

    ! Local variables:
    INTEGER                                                      :: ti, uv
    REAL(dp)                                                     :: N, dN_dx, dN_dy, beta_eff, tau_dx, tau_dy
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: single_row_ind
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddx_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dx2_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dxdy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dy2_val
    INTEGER                                                      :: single_row_nnz
    REAL(dp)                                                     :: Au, Av
    INTEGER                                                      :: k, tj, col_tju, col_tjv

    ! Relevant indices for this triangle
    ti     = mesh%n2tiuv( row_tiuv,1)
    uv     = mesh%n2tiuv( row_tiuv,2)

    ! N, dN/dx, dN/dy, beta_eff, tau_dx, and tau_dy on this triangle
    N        = DIVA%N_b(        ti)
    dN_dx    = DIVA%dN_dx_b(    ti)
    dN_dy    = DIVA%dN_dy_b(    ti)
    beta_eff = DIVA%beta_eff_b( ti)
    tau_dx   = DIVA%tau_dx_b(   ti)
    tau_dy   = DIVA%tau_dy_b(   ti)

    ! Allocate memory for single matrix rows
    ALLOCATE( single_row_ind(        mesh%nC_mem*2))
    ALLOCATE( single_row_ddx_val(    mesh%nC_mem*2))
    ALLOCATE( single_row_ddy_val(    mesh%nC_mem*2))
    ALLOCATE( single_row_d2dx2_val(  mesh%nC_mem*2))
    ALLOCATE( single_row_d2dxdy_val( mesh%nC_mem*2))
    ALLOCATE( single_row_d2dy2_val(  mesh%nC_mem*2))

    ! Read coefficients of the operator matrices
    CALL read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ind, single_row_ddx_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ind, single_row_ddy_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dxdy_b_b, ti, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    IF (uv == 1) THEN
      ! x-component

      DO k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_eff u + ...
        !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * N     * single_row_d2dx2_val(  k) + &  ! 4  N    d2u/dx2
             4._dp * dN_dx * single_row_ddx_val(    k) + &  ! 4 dN/dx du/dx
                     N     * single_row_d2dy2_val(  k) + &  !    N    d2u/dy2
                     dN_dy * single_row_ddy_val(    k)      !   dN/dy du/dy
        IF (tj == ti) Au = Au - beta_eff                    ! - beta_eff u

        Av = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2v/dxdy
             2._dp * dN_dx * single_row_ddy_val(    k) + &  ! 2 dN/dx dv/dy
                     dN_dy * single_row_ddx_val(    k)      !   dN/dy dv/dx

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      END DO

      ! Load vector
      bb( row_tiuv) = -tau_dx

    ELSEIF (uv == 2) THEN
      ! y-component

      DO k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_eff v + ...
        !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * N     * single_row_d2dy2_val(  k) + &  ! 4  N    d2v/dy2
             4._dp * dN_dy * single_row_ddy_val(    k) + &  ! 4 dN/dy dv/dy
                     N     * single_row_d2dx2_val(  k) + &  !    N    d2v/dx2
                     dN_dx * single_row_ddx_val(    k)      !   dN/dx dv/dx
        IF (tj == ti) Av = Av - beta_eff                    ! - beta_eff v

        Au = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2u/dxdy
             2._dp * dN_dy * single_row_ddx_val(    k) + &  ! 2 dN/dy du/dx
                     dN_dx * single_row_ddy_val(    k)      !   dN/dx du/dy

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      END DO

      ! Load vector
      bb( row_tiuv) = -tau_dy

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

  END SUBROUTINE calc_DIVA_stiffness_matrix_row_free

  SUBROUTINE calc_DIVA_sans_stiffness_matrix_row_free( mesh, DIVA, A_CSR, bb, row_tiuv)
    ! Add coefficients to this matrix row to represent the linearised DIVA
    !
    ! The DIVA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_eff u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_eff v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_eff u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_eff v = -tau_dy
    !
    ! The "sans" approximation neglects the gradients dN/dx, dN/dy of N:
    !
    !   4 N d2u/dx2 + N d2u/dy2 + 3 N d2v/dxdy - beta_eff u = -tau_dx
    !   4 N d2v/dy2 + N d2v/dx2 + 3 N d2u/dxdy - beta_eff v = -tau_dy
    !
    ! Dividing both sides by N yields:
    !
    !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_eff u / N = -tau_dx / N
    !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_eff v / N = -tau_dy / N
    !
    ! Note that there is no clear mathematical or physical reason why this should be allowed.
    ! However, while I (Tijn Berends, 2023) have found a few cases where there are noticeable
    ! differences    ! in the solutions (e.g. ISMIP-HOM experiments with high strain rates),
    ! most of the time the difference with respect to the full SSA/DIVA is very small.
    ! The "sans" option makes the solver quite a lot more stable and therefore faster.
    ! Someone really ought to perform some proper experiments to determine whether or not
    ! this should be the default.
    !
    ! We define the velocities u,v, the effective basal friction coefficient beta_eff, and the driving
    ! stress tau_d on the b-grid (triangles), and the effective viscosity eta and the
    ! product term N = eta H on the a-grid (vertices).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(IN)              :: DIVA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(mesh%ti1*2-1: mesh%ti2*2), INTENT(INOUT) :: bb
    INTEGER,                             INTENT(IN)              :: row_tiuv

    ! Local variables:
    INTEGER                                                      :: ti, uv
    REAL(dp)                                                     :: N, beta_eff, tau_dx, tau_dy
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: single_row_ind
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddx_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_ddy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dx2_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dxdy_val
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: single_row_d2dy2_val
    INTEGER                                                      :: single_row_nnz
    REAL(dp)                                                     :: Au, Av
    INTEGER                                                      :: k, tj, col_tju, col_tjv

    ! Relevant indices for this triangle
    ti     = mesh%n2tiuv( row_tiuv,1)
    uv     = mesh%n2tiuv( row_tiuv,2)

    ! N, beta_eff, tau_dx, and tau_dy on this triangle
    N        = DIVA%N_b(        ti)
    beta_eff = DIVA%beta_eff_b( ti)
    tau_dx   = DIVA%tau_dx_b(   ti)
    tau_dy   = DIVA%tau_dy_b(   ti)

    ! Allocate memory for single matrix rows
    ALLOCATE( single_row_ind(        mesh%nC_mem*2))
    ALLOCATE( single_row_ddx_val(    mesh%nC_mem*2))
    ALLOCATE( single_row_ddy_val(    mesh%nC_mem*2))
    ALLOCATE( single_row_d2dx2_val(  mesh%nC_mem*2))
    ALLOCATE( single_row_d2dxdy_val( mesh%nC_mem*2))
    ALLOCATE( single_row_d2dy2_val(  mesh%nC_mem*2))

    ! Read coefficients of the operator matrices
    CALL read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ind, single_row_ddx_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ind, single_row_ddy_val   , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dxdy_b_b, ti, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    CALL read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    IF (uv == 1) THEN
      ! x-component

      DO k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_eff u / N = -tau_dx / N

        ! Combine the mesh operators
        Au = 4._dp * single_row_d2dx2_val(  k) + &  ! 4 d2u/dx2
                     single_row_d2dy2_val(  k)      !   d2u/dy2
        IF (tj == ti) Au = Au - beta_eff  / N       ! - beta_eff u / N

        Av = 3._dp * single_row_d2dxdy_val( k)      ! 3 d2v/dxdy

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      END DO

      ! Load vector
      bb( row_tiuv) = -tau_dx / N

    ELSEIF (uv == 2) THEN
      ! y-component

      DO k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_eff v / N = -tau_dy / N

        ! Combine the mesh operators
        Av = 4._dp * single_row_d2dy2_val(  k) + &  ! 4 d2v/dy2
                     single_row_d2dx2_val(  k)      !   d2v/dx2
        IF (tj == ti) Av = Av - beta_eff / N        ! - beta_eff v / N

        Au = 3._dp * single_row_d2dxdy_val( k)      ! 3 d2u/dxdy

        ! Add coefficients to the stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      END DO

      ! Load vector
      bb( row_tiuv) = -tau_dy / N

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

  END SUBROUTINE calc_DIVA_sans_stiffness_matrix_row_free

  SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_west( mesh, DIVA, A_CSR, bb, row_tiuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! western domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA),  INTENT(IN)              :: DIVA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tiuv

    ! Local variables:
    INTEGER                                                      :: ti,uv,row_ti
    INTEGER                                                      :: tj, col_tjuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_west == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_west == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * DIVA%u_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * DIVA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_west == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_west == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * DIVA%v_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * DIVA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_west "' // TRIM( C%BC_u_west) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_west

  SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_east( mesh, DIVA, A_CSR, bb, row_tiuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! eastern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA),  INTENT(IN)              :: DIVA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tiuv

    ! Local variables:
    INTEGER                                                      :: ti,uv,row_ti
    INTEGER                                                      :: tj, col_tjuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_east == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_east == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * DIVA%u_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * DIVA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_east == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_east == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * DIVA%v_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * DIVA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_east "' // TRIM( C%BC_u_east) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_east

  SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_south( mesh, DIVA, A_CSR, bb, row_tiuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! southern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA),  INTENT(IN)              :: DIVA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tiuv

    ! Local variables:
    INTEGER                                                      :: ti,uv,row_ti
    INTEGER                                                      :: tj, col_tjuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_south == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_south == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * DIVA%u_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * DIVA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_south == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_south == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * DIVA%v_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * DIVA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_south "' // TRIM( C%BC_u_south) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_south

  SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_north( mesh, DIVA, A_CSR, bb, row_tiuv)
    ! Add coefficients to this matrix row to represent boundary conditions at the
    ! northern domain border.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA),  INTENT(IN)              :: DIVA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: A_CSR
    REAL(dp), DIMENSION(A_CSR%i1:A_CSR%i2), INTENT(INOUT)        :: bb
    INTEGER,                             INTENT(IN)              :: row_tiuv

    ! Local variables:
    INTEGER                                                      :: ti,uv,row_ti
    INTEGER                                                      :: tj, col_tjuv
    INTEGER,  DIMENSION(mesh%nC_mem)                             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem)                             :: wti_copy
    REAL(dp)                                                     :: u_fixed, v_fixed
    INTEGER                                                      :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_north == 'zero') THEN
        ! u = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_u_north == 'periodic_ISMIP-HOM') THEN
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * DIVA%u_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * DIVA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

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
          col_tjuv = mesh%tiuv2n( tj,uv)
          CALL add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        END DO
        IF (n_neighbours == 0) CALL crash('whaa!')
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * REAL( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_north == 'zero') THEN
        ! v = 0

        ! Stiffness matrix
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      ELSEIF (C%BC_v_north == 'periodic_ISMIP-HOM') THEN
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        CALL add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        DO n = 1, mesh%nC_mem
          tj = ti_copy( n)
          IF (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * DIVA%v_b_prev( tj)
        END DO
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * DIVA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      ELSE
        CALL crash('unknown BC_u_north "' // TRIM( C%BC_u_north) // '"!')
      END IF

    ELSE
      CALL crash('uv can only be 1 or 2!')
    END IF

  END SUBROUTINE calc_DIVA_stiffness_matrix_row_BC_north

! == Calculate several intermediate terms in the DIVA

  SUBROUTINE calc_driving_stress( mesh, ice, DIVA)
    ! Calculate the driving stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_driving_stress'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: Hi_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: dHs_dy_b
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    ALLOCATE( Hi_b(     mesh%ti1:mesh%ti2))
    ALLOCATE( dHs_dx_b( mesh%ti1:mesh%ti2))
    ALLOCATE( dHs_dy_b( mesh%ti1:mesh%ti2))

    ! Calculate Hi, dHs/dx, and dHs/dy on the b-grid
    CALL map_a_b_2D( mesh, ice%Hi, Hi_b    )
    CALL ddx_a_b_2D( mesh, ice%Hs, dHs_dx_b)
    CALL ddy_a_b_2D( mesh, ice%Hs, dHs_dy_b)

    ! Calculate the driving stress
    DO ti = mesh%ti1, mesh%ti2
      DIVA%tau_dx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      DIVA%tau_dy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    END DO

    ! Clean up after yourself
    DEALLOCATE( Hi_b    )
    DEALLOCATE( dHs_dx_b)
    DEALLOCATE( dHs_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress

  SUBROUTINE calc_horizontal_strain_rates( mesh, DIVA)
    ! Calculate the horizontal strain rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_horizontal_strain_rates'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the strain rates
    CALL ddx_b_a_2D( mesh, DIVA%u_vav_b, DIVA%du_dx_a)
    CALL ddy_b_a_2D( mesh, DIVA%u_vav_b, DIVA%du_dy_a)
    CALL ddx_b_a_2D( mesh, DIVA%v_vav_b, DIVA%dv_dx_a)
    CALL ddy_b_a_2D( mesh, DIVA%v_vav_b, DIVA%dv_dy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_horizontal_strain_rates

  SUBROUTINE calc_vertical_shear_strain_rates( mesh, DIVA)
    ! Calculate the vertical shear strain rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_vertical_shear_strain_rates'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      :: du_dz_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      :: dv_dz_3D_b
    INTEGER                                                      :: ti,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( du_dz_3D_b( mesh%ti1:mesh%ti2,mesh%nz))
    ALLOCATE( dv_dz_3D_b( mesh%ti1:mesh%ti2,mesh%nz))

    ! Calculate (parameterised) vertical shear strain rates on the b-grid (Lipscomb et al., 2019, Eq. 36)
    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, mesh%nz
      du_dz_3D_b( ti,k) = DIVA%tau_bx_b( ti) * mesh%zeta( k) / MAX( C%visc_eff_min, DIVA%eta_3D_b( ti,k))
      dv_dz_3D_b( ti,k) = DIVA%tau_by_b( ti) * mesh%zeta( k) / MAX( C%visc_eff_min, DIVA%eta_3D_b( ti,k))
    END DO
    END DO

    ! Map vertical shear strain rates from the b-grid to the a-grid
    CALL map_b_a_3D( mesh, du_dz_3D_b, DIVA%du_dz_3D_a)
    CALL map_b_a_3D( mesh, dv_dz_3D_b, DIVA%dv_dz_3D_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_shear_strain_rates

  SUBROUTINE calc_effective_viscosity( mesh, ice, DIVA, Glens_flow_law_epsilon_sq_0_applied)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA
    REAL(dp),                            INTENT(IN)              :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_viscosity'
    INTEGER                                                      :: vi,k
    REAL(dp)                                                     :: A_min, eta_max
    REAL(dp), DIMENSION( mesh%nz)                                :: prof

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! Calculate the effective viscosity eta
    IF (C%choice_flow_law == 'Glen') THEN
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      CALL calc_ice_rheology_Glen( mesh, ice)

      ! Calculate effective viscosity
      DO vi = mesh%vi1, mesh%vi2
      DO k  = 1, mesh%nz
        DIVA%eta_3D_a( vi,k) = calc_effective_viscosity_Glen_3D_uv_only( &
          Glens_flow_law_epsilon_sq_0_applied, &
          DIVA%du_dx_a( vi), DIVA%du_dy_a( vi), DIVA%du_dz_3D_a( vi,k), &
          DIVA%dv_dx_a( vi), DIVA%dv_dy_a( vi), DIVA%dv_dz_3D_a( vi,k), ice%A_flow( vi,k))
      END DO
      END DO

    ELSE
      CALL crash('unknown choice_flow_law "' // TRIM( C%choice_flow_law) // '"!')
    END IF

    ! Safety
    DIVA%eta_3D_a = MIN( MAX( DIVA%eta_3D_a, C%visc_eff_min), eta_max)

    ! Map effective viscosity to the b-grid
    CALL map_a_b_3D( mesh, DIVA%eta_3D_a, DIVA%eta_3D_b)

    ! Calculate vertically averaged effective viscosity on the a-grid
    DO vi = mesh%vi1, mesh%vi2
      prof = DIVA%eta_3D_a( vi,:)
      DIVA%eta_vav_a( vi) = vertical_average( mesh%zeta, prof)
    END DO

    ! Calculate the product term N = eta * H on the a-grid
    DO vi = mesh%vi1, mesh%vi2
      DIVA%N_a( vi) = DIVA%eta_vav_a( vi) * MAX( 0.1, ice%Hi( vi))
    END DO

    ! Calculate the product term N and its gradients on the b-grid
    CALL map_a_b_2D( mesh, DIVA%N_a, DIVA%N_b    )
    CALL ddx_a_b_2D( mesh, DIVA%N_a, DIVA%dN_dx_b)
    CALL ddy_a_b_2D( mesh, DIVA%N_a, DIVA%dN_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity

  SUBROUTINE calc_F_integrals( mesh, ice, DIVA)
    ! Calculate the F-integrals on the a-grid (Lipscomb et al. (2019), Eq. 30)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_F_integrals'
    INTEGER                                                      :: vi,k
    REAL(dp), DIMENSION( mesh%nz)                                :: prof

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! F1
      DO k = 1, mesh%nz
        prof( k) = (mesh%zeta( k)    / DIVA%eta_3D_a( vi,k))
      END DO
      DIVA%F1_3D_a( vi,:) = -MAX( 0.1_dp, ice%Hi( vi)) * integrate_from_zeta_is_one_to_zeta_is_zetap( mesh%zeta, prof)

      ! F2
      DO k = 1, mesh%nz
        prof( k) = (mesh%zeta( k)**2 / DIVA%eta_3D_a( vi,k))
      END DO
      DIVA%F2_3D_a( vi,:) = -MAX( 0.1_dp, ice%Hi( vi)) * integrate_from_zeta_is_one_to_zeta_is_zetap( mesh%zeta, prof)

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Map F-integrals from the a-grid to the b-grid
    CALL map_a_b_3D( mesh, DIVA%F1_3D_a, DIVA%F1_3D_b)
    CALL map_a_b_3D( mesh, DIVA%F2_3D_a, DIVA%F2_3D_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_F_integrals

  SUBROUTINE calc_effective_basal_friction_coefficient( mesh, ice, DIVA)
    ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_basal_friction_coefficient'
    INTEGER                                                      :: vi,ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!
    CALL calc_basal_friction_coefficient( mesh, ice, DIVA%u_base_b, DIVA%v_base_b)

    ! Calculate beta_eff on the a-grid
    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Exception for the case of no sliding (Lipscomb et al., 2019, Eq. 35)

      DO vi = mesh%vi1, mesh%vi2
        DIVA%beta_eff_a( vi) = 1._dp / DIVA%F2_3D_a( vi,1)
      END DO

    ELSE ! IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Lipscomb et al., 2019, Eq. 33

      DO vi = mesh%vi1, mesh%vi2
        DIVA%beta_eff_a( vi) = ice%basal_friction_coefficient( vi) / (1._dp + ice%basal_friction_coefficient( vi) * DIVA%F2_3D_a( vi,1))
      END DO

    END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

    ! Map basal friction coefficient beta_b and effective basal friction coefficient beta_eff to the b-grid
    CALL map_a_b_2D( mesh, ice%basal_friction_coefficient, DIVA%basal_friction_coefficient_b)
    CALL map_a_b_2D( mesh, DIVA%beta_eff_a               , DIVA%beta_eff_b                  )

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    IF (C%do_GL_subgrid_friction) THEN
      ! On the b-grid
      DO ti = mesh%ti1, mesh%ti2
        DIVA%beta_eff_b( ti) = DIVA%beta_eff_b( ti) * ice%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      END DO
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_basal_friction_coefficient

  SUBROUTINE calc_basal_shear_stress( mesh, DIVA)
    ! Calculate the basal shear stress (Lipscomb et al., 2019, just above Eq. 33)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_basal_shear_stress'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      ! Lipscomb et al., 2019, just above Eq. 33
      DIVA%tau_bx_b( ti) = DIVA%u_vav_b( ti) * DIVA%beta_eff_b( ti)
      DIVA%tau_by_b( ti) = DIVA%v_vav_b( ti) * DIVA%beta_eff_b( ti)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_shear_stress

  SUBROUTINE calc_basal_velocities( mesh, DIVA)
    ! Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_basal_shear_stress'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Exception for the case of no sliding

      DIVA%u_base_b = 0._dp
      DIVA%v_base_b = 0._dp

    ELSE ! IF (C%choice_sliding_law == 'no_sliding') THEN

      ! Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)
      DO ti = mesh%ti1, mesh%ti2
        DIVA%u_base_b( ti) = DIVA%u_vav_b( ti) / (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F2_3D_b( ti,1))
        DIVA%v_base_b( ti) = DIVA%v_vav_b( ti) / (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F2_3D_b( ti,1))
      END DO

    END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_velocities

  SUBROUTINE calc_3D_velocities( mesh, DIVA)
    ! Calculate 3D velocities (Lipscomb et al., 2019, Eq. 29)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_3D_velocities'
    INTEGER                                                      :: ti,k

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Exception for the case of no sliding

      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, mesh%nz
        ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
        DIVA%u_3D_b( ti,k) = DIVA%tau_bx_b( ti) * DIVA%F1_3D_b( ti,k)
        DIVA%v_3D_b( ti,k) = DIVA%tau_by_b( ti) * DIVA%F1_3D_b( ti,k)
      END DO
      END DO

    ELSE ! IF (C%choice_sliding_law == 'no_sliding') THEN

      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, mesh%nz
        ! Lipscomb et al., 2019, Eq. 29
        DIVA%u_3D_b( ti,k) = DIVA%u_base_b( ti) * (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F1_3D_b( ti,k))
        DIVA%v_3D_b( ti,k) = DIVA%v_base_b( ti) * (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F1_3D_b( ti,k))
      END DO
      END DO

    END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_3D_velocities

! == Some useful tools for improving numerical stability of the viscosity iteration

  SUBROUTINE relax_viscosity_iterations( mesh, DIVA, visc_it_relax)
    ! Reduce the change between velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA
    REAL(dp),                            INTENT(IN)              :: visc_it_relax

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'relax_viscosity_iterations'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      DIVA%u_vav_b( ti) = (visc_it_relax * DIVA%u_vav_b( ti)) + ((1._dp - visc_it_relax) * DIVA%u_b_prev( ti))
      DIVA%v_vav_b( ti) = (visc_it_relax * DIVA%v_vav_b( ti)) + ((1._dp - visc_it_relax) * DIVA%v_b_prev( ti))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_viscosity_iterations

  SUBROUTINE calc_visc_iter_UV_resid( mesh, DIVA, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(IN)              :: DIVA
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

    DO ti = mesh%ti1, mesh%ti2

      res1 = res1 + (DIVA%u_vav_b( ti) - DIVA%u_b_prev( ti))**2
      res1 = res1 + (DIVA%v_vav_b( ti) - DIVA%v_b_prev( ti))**2

      res2 = res2 + (DIVA%u_vav_b( ti) + DIVA%u_b_prev( ti))**2
      res2 = res2 + (DIVA%v_vav_b( ti) + DIVA%v_b_prev( ti))**2

    END DO

    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_visc_iter_UV_resid

  SUBROUTINE apply_velocity_limits( mesh, DIVA)
    ! Limit velocities for improved stability

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'apply_velocity_limits'
    INTEGER                                                      :: ti
    REAL(dp)                                                     :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate absolute speed
      uabs = SQRT( DIVA%u_vav_b( ti)**2 + DIVA%v_vav_b( ti)**2)

      ! Reduce velocities if neceDIVAry
      IF (uabs > C%vel_max) THEN
        DIVA%u_vav_b( ti) = DIVA%u_vav_b( ti) * C%vel_max / uabs
        DIVA%v_vav_b( ti) = DIVA%v_vav_b( ti) * C%vel_max / uabs
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits

! == Initialisation

  SUBROUTINE initialise_DIVA_velocities_from_file( mesh, DIVA, region_name)
    ! Initialise the velocities for the DIVA solver from an external NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT) :: DIVA
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_DIVA_velocities_from_file'
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
    IF (par%master) WRITE(0,*) '   Initialising DIVA velocities from file "' // colour_string( TRIM( filename),'light blue') // '"...'

    ! Read velocities from the file
    IF (timeframe == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_mesh_file_2D_b( filename, 'u_vav_b' , DIVA%u_vav_b )
      CALL read_field_from_mesh_file_2D_b( filename, 'v_vav_b' , DIVA%v_vav_b )
      CALL read_field_from_mesh_file_2D_b( filename, 'tau_bx_b', DIVA%tau_bx_b)
      CALL read_field_from_mesh_file_2D_b( filename, 'tau_by_b', DIVA%tau_by_b)
      CALL read_field_from_mesh_file_3D_b( filename, 'eta_3D_b', DIVA%eta_3D_b)
    ELSE
      ! Read specified timeframe
      CALL read_field_from_mesh_file_2D_b( filename, 'u_vav_b' , DIVA%u_vav_b , time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D_b( filename, 'v_vav_b' , DIVA%v_vav_b , time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D_b( filename, 'tau_bx_b', DIVA%tau_bx_b, time_to_read = timeframe)
      CALL read_field_from_mesh_file_2D_b( filename, 'tau_by_b', DIVA%tau_by_b, time_to_read = timeframe)
      CALL read_field_from_mesh_file_3D_b( filename, 'eta_3D_b', DIVA%eta_3D_b, time_to_read = timeframe)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_DIVA_velocities_from_file

  SUBROUTINE allocate_DIVA_solver( mesh, DIVA)
    ! Allocate memory the DIVA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(OUT)   :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_DIVA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Solution
    ALLOCATE( DIVA%u_vav_b(      mesh%ti1:mesh%ti2))                   ! [m yr^-1] 2-D vertically averaged horizontal ice velocity
    DIVA%u_vav_b = 0._dp
    ALLOCATE( DIVA%v_vav_b(      mesh%ti1:mesh%ti2))
    DIVA%v_vav_b = 0._dp
    ALLOCATE( DIVA%u_base_b(     mesh%ti1:mesh%ti2))                   ! [m yr^-1] 2-D horizontal ice velocity at the ice base
    DIVA%u_base_b = 0._dp
    ALLOCATE( DIVA%v_base_b(     mesh%ti1:mesh%ti2))
    DIVA%v_base_b = 0._dp
    ALLOCATE( DIVA%u_3D_b(       mesh%ti1:mesh%ti2,mesh%nz))           ! [m yr^-1] 3-D horizontal ice velocity
    DIVA%u_3D_b = 0._dp
    ALLOCATE( DIVA%v_3D_b(       mesh%ti1:mesh%ti2,mesh%nz))
    DIVA%v_3D_b = 0._dp

    ! Intermediate data fields
    ALLOCATE( DIVA%du_dx_a(      mesh%vi1:mesh%vi2))                   ! [yr^-1] 2-D horizontal strain rates
    DIVA%du_dx_a = 0._dp
    ALLOCATE( DIVA%du_dy_a(      mesh%vi1:mesh%vi2))
    DIVA%du_dy_a = 0._dp
    ALLOCATE( DIVA%dv_dx_a(      mesh%vi1:mesh%vi2))
    DIVA%dv_dx_a = 0._dp
    ALLOCATE( DIVA%dv_dy_a(      mesh%vi1:mesh%vi2))
    DIVA%dv_dy_a = 0._dp
    ALLOCATE( DIVA%du_dz_3D_a(   mesh%vi1:mesh%vi2,mesh%nz))           ! [yr^-1] 3-D vertical shear strain rates
    DIVA%du_dz_3D_a = 0._dp
    ALLOCATE( DIVA%dv_dz_3D_a(   mesh%vi1:mesh%vi2,mesh%nz))
    DIVA%dv_dz_3D_a = 0._dp
    ALLOCATE( DIVA%eta_3D_a(     mesh%vi1:mesh%vi2,mesh%nz))           ! Effective viscosity
    DIVA%eta_3D_a = 0._dp
    ALLOCATE( DIVA%eta_3D_b(     mesh%ti1:mesh%ti2,mesh%nz))
    DIVA%eta_3D_b = 0._dp
    ALLOCATE( DIVA%eta_vav_a(    mesh%vi1:mesh%vi2))
    DIVA%eta_vav_a = 0._dp
    ALLOCATE( DIVA%N_a(          mesh%vi1:mesh%vi2))                   ! Product term N = eta * H
    DIVA%N_a = 0._dp
    ALLOCATE( DIVA%N_b(          mesh%ti1:mesh%ti2))
    DIVA%N_b = 0._dp
    ALLOCATE( DIVA%dN_dx_b(      mesh%ti1:mesh%ti2))                   ! Gradients of N
    DIVA%dN_dx_b = 0._dp
    ALLOCATE( DIVA%dN_dy_b(      mesh%ti1:mesh%ti2))
    DIVA%dN_dy_b = 0._dp
    ALLOCATE( DIVA%F1_3D_a(      mesh%vi1:mesh%vi2,mesh%nz))           ! F-integrals
    DIVA%F1_3D_a = 0._dp
    ALLOCATE( DIVA%F2_3D_a(      mesh%vi1:mesh%vi2,mesh%nz))
    DIVA%F2_3D_a = 0._dp
    ALLOCATE( DIVA%F1_3D_b(      mesh%ti1:mesh%ti2,mesh%nz))
    DIVA%F1_3D_b = 0._dp
    ALLOCATE( DIVA%F2_3D_b(      mesh%ti1:mesh%ti2,mesh%nz))
    DIVA%F2_3D_b = 0._dp
    ALLOCATE( DIVA%basal_friction_coefficient_b(     mesh%ti1:mesh%ti2)) ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    DIVA%basal_friction_coefficient_b = 0._dp
    ALLOCATE( DIVA%beta_eff_a(   mesh%vi1:mesh%vi2))                   ! "Effective" friction coefficient (turning the SSA into the DIVA)
    DIVA%beta_eff_a = 0._dp
    ALLOCATE( DIVA%beta_eff_b(   mesh%ti1:mesh%ti2))
    DIVA%beta_eff_b = 0._dp
    ALLOCATE( DIVA%tau_bx_b(     mesh%ti1:mesh%ti2))                   ! Basal shear stress
    DIVA%tau_bx_b = 0._dp
    ALLOCATE( DIVA%tau_by_b(     mesh%ti1:mesh%ti2))
    DIVA%tau_by_b = 0._dp
    ALLOCATE( DIVA%tau_dx_b(     mesh%ti1:mesh%ti2))                   ! Driving stress
    DIVA%tau_dx_b = 0._dp
    ALLOCATE( DIVA%tau_dy_b(     mesh%ti1:mesh%ti2))
    DIVA%tau_dy_b = 0._dp
    ALLOCATE( DIVA%u_b_prev(     mesh%nTri))                   ! Velocity solution from previous viscosity iteration
    DIVA%u_b_prev = 0._dp
    ALLOCATE( DIVA%v_b_prev(     mesh%nTri))
    DIVA%v_b_prev = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_DIVA_solver

! == Restart NetCDF files

  SUBROUTINE write_to_restart_file_DIVA( mesh, DIVA, time)
    ! Write to the restart NetCDF file for the DIVA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(IN)              :: DIVA
    REAL(dp),                            INTENT(IN)              :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'write_to_restart_file_DIVA'
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Writing to DIVA restart file "' // &
      colour_string( TRIM( DIVA%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( DIVA%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( DIVA%restart_filename, ncid, time)

    ! Write the velocity fields to the file
    CALL write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'u_vav_b' , DIVA%u_vav_b )
    CALL write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'v_vav_b' , DIVA%v_vav_b )
    CALL write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_bx_b', DIVA%tau_bx_b)
    CALL write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_by_b', DIVA%tau_by_b)
    CALL write_to_field_multopt_mesh_dp_3D_b( mesh, DIVA%restart_filename, ncid, 'eta_3D_b', DIVA%eta_3D_b)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_DIVA

  SUBROUTINE create_restart_file_DIVA( mesh, DIVA)
    ! Create a restart NetCDF file for the DIVA solver
    ! Includes generation of the procedural filename (e.g. "restart_DIVA_00001.nc")

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'create_restart_file_DIVA'
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
    filename_base = TRIM( C%output_dir) // 'restart_ice_velocity_DIVA'
    CALL generate_filename_XXXXXdotnc( filename_base, DIVA%restart_filename)

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Creating DIVA restart file "' // &
      colour_string( TRIM( DIVA%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( DIVA%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( DIVA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( DIVA%restart_filename, ncid)

    ! Add a zeta dimension to the file
    CALL add_zeta_dimension_to_file( DIVA%restart_filename, ncid, mesh%zeta)

    ! Add the velocity fields to the file
    CALL add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'u_vav_b' , long_name = 'Vertically averaged horizontal ice velocity in the x-direction', units = 'm/yr')
    CALL add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'v_vav_b' , long_name = 'Vertically averaged horizontal ice velocity in the y-direction', units = 'm/yr')
    CALL add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_bx_b', long_name = 'Basal shear stress in the x-direction', units = 'Pa')
    CALL add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_by_b', long_name = 'Basal shear stress in the y-direction', units = 'Pa')
    CALL add_field_mesh_dp_3D_b( DIVA%restart_filename, ncid, 'eta_3D_b', long_name = '3-D effective viscosity')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_DIVA

END MODULE ice_velocity_DIVA
