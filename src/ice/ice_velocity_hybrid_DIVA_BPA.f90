MODULE ice_velocity_hybrid_DIVA_BPA

  ! Routines for calculating ice velocities using the hybrid DIVA/BPA

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
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D, &
                                                                     save_variable_as_netcdf_logical_1D
  USE parameters
  USE petsc_basic                                            , ONLY: solve_matrix_equation_CSR_PETSc
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_DIVA, type_ice_velocity_solver_BPA, type_ice_velocity_solver_hybrid
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  USE ice_velocity_DIVA                                      , ONLY: allocate_DIVA_solver, remap_DIVA_solver, calc_DIVA_stiffness_matrix_row_free, &
                                                                     calc_DIVA_stiffness_matrix_row_BC_west, calc_DIVA_stiffness_matrix_row_BC_east, &
                                                                     calc_DIVA_stiffness_matrix_row_BC_south, calc_DIVA_stiffness_matrix_row_BC_north, &
                                                                     calc_DIVA_sans_stiffness_matrix_row_free
  USE ice_velocity_BPA                                       , ONLY: allocate_BPA_solver , remap_BPA_solver, calc_BPA_stiffness_matrix_row_free, &
                                                                     calc_BPA_stiffness_matrix_row_BC_west, calc_BPA_stiffness_matrix_row_BC_east, &
                                                                     calc_BPA_stiffness_matrix_row_BC_south, calc_BPA_stiffness_matrix_row_BC_north, &
                                                                     calc_BPA_stiffness_matrix_row_BC_base, calc_BPA_stiffness_matrix_row_BC_surf
  USE mesh_operators                                         , ONLY: calc_3D_matrix_operators_mesh, map_a_b_2D, map_b_a_2D, map_b_a_3D, map_a_b_3D
  use mesh_ROI_polygons
  USE math_utilities                                         , ONLY: is_in_polygon, is_in_polygons
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_logical_1D
  USE ice_model_utilities                                    , ONLY: calc_zeta_gradients
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist, add_empty_row_CSR_dist
  USE grid_basic                                             , ONLY: type_grid, gather_gridded_data_to_master_int_2D, calc_grid_mask_as_polygons
  USE netcdf_basic                                           , ONLY: open_existing_netcdf_file_for_reading, close_netcdf_file
  USE netcdf_input                                           , ONLY: setup_xy_grid_from_file, read_field_from_xy_file_2D_int

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Main routines

  SUBROUTINE initialise_hybrid_DIVA_BPA_solver( mesh, hybrid, region_name)
    ! Initialise the hybrid DIVA/BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_hybrid), INTENT(OUT)   :: hybrid
    CHARACTER(LEN=3),                      INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                        :: routine_name = 'initialise_hybrid_DIVA_BPA_solver'
    CHARACTER(LEN=256)                                   :: choice_initial_velocity

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    CALL allocate_hybrid_DIVA_BPA_solver( mesh, hybrid)

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
      hybrid%u_vav_b = 0._dp
      hybrid%v_vav_b = 0._dp
      hybrid%u_bk    = 0._dp
      hybrid%v_bk    = 0._dp
    ELSEIF (choice_initial_velocity == 'read_from_file') THEN
      CALL crash('restarting ice velocities not yet possible for the hybrid DIVA/BPA!')
    ELSE
      CALL crash('unknown choice_initial_velocity "' // TRIM( choice_initial_velocity) // '"!')
    END IF

    ! Set tolerances for PETSc matrix solver for the linearised hybrid DIVA/BPA
    hybrid%PETSc_rtol   = C%stress_balance_PETSc_rtol
    hybrid%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_hybrid_DIVA_BPA_solver

  SUBROUTINE solve_hybrid_DIVA_BPA( mesh, ice, hybrid, region_name, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Calculate ice velocities by solving the hybrid DIVA/BPA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                  INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_hybrid), INTENT(INOUT)           :: hybrid
    CHARACTER(LEN=3),                      INTENT(IN)              :: region_name
    INTEGER,  DIMENSION(:    ),            INTENT(IN)   , OPTIONAL :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:    ),            INTENT(IN)   , OPTIONAL :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:    ),            INTENT(IN)   , OPTIONAL :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'solve_hybrid_DIVA_BPA'
    LOGICAL                                                        :: grounded_ice_exists
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                        :: BC_prescr_mask_b_applied
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                        :: BC_prescr_u_b_applied
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                        :: BC_prescr_v_b_applied
    INTEGER                                                        :: viscosity_iteration_i
    LOGICAL                                                        :: has_converged
    REAL(dp)                                                       :: resid_UV, resid_UV_prev
    REAL(dp)                                                       :: uv_min, uv_max
    REAL(dp)                                                       :: visc_it_relax_applied
    REAL(dp)                                                       :: Glens_flow_law_epsilon_sq_0_applied
    INTEGER                                                        :: nit_diverg_consec

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there is no grounded ice, no need (in fact, no way) to solve the velocities
    grounded_ice_exists = ANY( ice%mask_grounded_ice)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. grounded_ice_exists) THEN
      hybrid%u_vav_b = 0._dp
      hybrid%v_vav_b = 0._dp
      hybrid%u_bk    = 0._dp
      hybrid%v_bk    = 0._dp
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

    ! Calculate zeta gradients
    CALL calc_zeta_gradients( mesh, ice)

    ! Calculate 3-D matrix operators for the current ice geometry
    CALL calc_3D_matrix_operators_mesh( mesh, ice)

    ! Calculate the driving stress
    CALL calc_driving_stress_DIVA( mesh, ice, hybrid%DIVA)
    CALL calc_driving_stress_BPA ( mesh, ice, hybrid%BPA )

    ! Calculate the solving masks for the hybrid solver
    CALL calc_hybrid_solver_masks_basic( mesh, ice, hybrid, region_name)

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

    ! == Calculate secondary terms in the DIVA
    ! ========================================

      ! Calculate the horizontal strain rates for the current velocity solution
      CALL calc_horizontal_strain_rates_DIVA( mesh, hybrid%DIVA)

      ! Calculate the vertical shear strain rates
      CALL calc_vertical_shear_strain_rates_DIVA( mesh, hybrid%DIVA)

      ! Calculate the effective viscosity for the current velocity solution
      CALL calc_effective_viscosity_DIVA( mesh, ice, hybrid%DIVA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the F-integrals
      CALL calc_F_integrals_DIVA( mesh, ice, hybrid%DIVA)

      ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)
      CALL calc_effective_basal_friction_coefficient_DIVA( mesh, ice, hybrid%DIVA)

    ! == Calculate secondary terms in the BPA
    ! =======================================

      ! Calculate the strain rates for the current velocity solution
      CALL calc_strain_rates_BPA( mesh, hybrid%BPA)

      ! Calculate the effective viscosity for the current velocity solution
      CALL calc_effective_viscosity_BPA( mesh, ice, hybrid%BPA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      CALL calc_applied_basal_friction_coefficient_BPA( mesh, ice, hybrid%BPA)

    ! == Solve the linearised hybrid DIVA/BPA
    ! =======================================

      ! Solve the linearised hybrid DIVA/BPA to calculate a new velocity solution
      CALL solve_hybrid_DIVA_BPA_linearised( mesh, ice, hybrid, BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Copy results to the DIVA and BPA structures
      hybrid%DIVA%u_vav_b = hybrid%u_vav_b
      hybrid%DIVA%v_vav_b = hybrid%v_vav_b
      hybrid%DIVA%u_3D_b  = hybrid%u_bk
      hybrid%DIVA%v_3D_b  = hybrid%v_bk
      hybrid%BPA%u_bk     = hybrid%u_bk
      hybrid%BPA%v_bk     = hybrid%v_bk

    ! == Calculate more secondary terms in the DIVA
    ! =============================================

      ! Calculate basal velocities
      CALL calc_basal_velocities_DIVA( mesh, hybrid%DIVA)

      ! Calculate basal shear stress
      CALL calc_basal_shear_stress_DIVA( mesh, hybrid%DIVA)

    ! == Improve stability and check for convergence
    ! ==============================================

      ! Limit velocities for improved stability
      CALL apply_velocity_limits( mesh, hybrid)

      ! Reduce the change between velocity solutions
      CALL relax_viscosity_iterations( mesh, hybrid, visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      CALL calc_visc_iter_UV_resid( mesh, hybrid, resid_UV)

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
!      uv_min = MINVAL( hybrid%u_bk)
!      uv_max = MAXVAL( hybrid%u_bk)
!      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
!      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
!      IF (par%master) WRITE(0,*) '    hybrid DIVA/BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

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
    DEALLOCATE( BC_prescr_mask_b_applied)
    DEALLOCATE( BC_prescr_u_b_applied   )
    DEALLOCATE( BC_prescr_v_b_applied   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_hybrid_DIVA_BPA

  SUBROUTINE remap_hybrid_DIVA_BPA_solver( mesh_old, mesh_new, hybrid)
    ! Remap the hybrid DIVA/BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                       INTENT(IN)    :: mesh_new
    TYPE(type_ice_velocity_solver_hybrid), INTENT(INOUT) :: hybrid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                        :: routine_name = 'remap_hybrid_DIVA_BPA_solver'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE              :: u_vav_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE              :: v_vav_a
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE              :: u_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE              :: v_ak

    ! Add routine to path
    CALL init_routine( routine_name)

  ! Remap the fields that are re-used during the viscosity iteration
  ! ================================================================

    ! Allocate memory for velocities on the a-grid (vertices)
    ALLOCATE( u_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    ALLOCATE( v_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    ALLOCATE( u_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    ALLOCATE( v_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    CALL map_b_a_2D( mesh_old, hybrid%u_vav_b, u_vav_a)
    CALL map_b_a_2D( mesh_old, hybrid%v_vav_b, v_vav_a)
    CALL map_b_a_3D( mesh_old, hybrid%u_bk   , u_ak   )
    CALL map_b_a_3D( mesh_old, hybrid%v_bk   , v_ak   )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, u_vav_a, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, v_vav_a, '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, u_ak   , '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, v_ak   , '2nd_order_conservative')

    ! Reallocate memory for the data on the triangles
    CALL reallocate_bounds( hybrid%u_vav_b, mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( hybrid%v_vav_b, mesh_new%ti1, mesh_new%ti2             )
    CALL reallocate_bounds( hybrid%u_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( hybrid%v_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    CALL map_a_b_2D( mesh_new, u_vav_a, hybrid%u_vav_b)
    CALL map_a_b_2D( mesh_new, v_vav_a, hybrid%v_vav_b)
    CALL map_a_b_3D( mesh_new, u_ak   , hybrid%u_bk   )
    CALL map_a_b_3D( mesh_new, v_ak   , hybrid%v_bk   )

    ! Clean up after yourself
    DEALLOCATE( u_vav_a)
    DEALLOCATE( v_vav_a)
    DEALLOCATE( u_ak   )
    DEALLOCATE( v_ak   )

  ! Remap data of the separate DIVA and BPA solvers
  ! ===============================================

    CALL remap_DIVA_solver( mesh_old, mesh_new, hybrid%DIVA)
    CALL remap_BPA_solver(  mesh_old, mesh_new, hybrid%BPA )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_hybrid_DIVA_BPA_solver

! == Basic masks and translation tables for the hybrid solver

  SUBROUTINE calc_hybrid_solver_masks_basic( mesh, ice, hybrid, region_name)
    ! Calculate the solving masks for the hybrid solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_model),                  INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_hybrid), INTENT(INOUT)           :: hybrid
    CHARACTER(LEN=3),                      INTENT(IN)              :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'calc_hybrid_solver_masks_basic'
    CHARACTER(LEN=256)                                             :: choice_hybrid_DIVA_BPA_mask

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    hybrid%mask_BPA_b  = .FALSE.
    hybrid%mask_DIVA_b = .FALSE.

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE DEFAULT
        CALL crash('unknown region name "' // TRIM( region_name) // '"!')
      CASE ('NAM')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_NAM
      CASE ('EAS')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_EAS
      CASE ('GRL')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_GRL
      CASE ('ANT')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_ANT
    END SELECT

    SELECT CASE (choice_hybrid_DIVA_BPA_mask)
      CASE DEFAULT
        CALL crash('unknown choice_hybrid_DIVA_BPA_mask "' // TRIM( choice_hybrid_DIVA_BPA_mask) // '"!')
      CASE ('ROI')
        CALL calc_hybrid_solver_masks_basic_ROI( mesh, ice, hybrid, region_name)
      CASE ('read_from_file')
        CALL calc_hybrid_solver_masks_basic_read_from_file( mesh, ice, hybrid, region_name)
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_hybrid_solver_masks_basic

  SUBROUTINE calc_hybrid_solver_masks_basic_ROI( mesh, ice, hybrid, region_name)
    ! Calculate the solving masks for the hybrid solver
    !
    ! Solve the BPA only in the specified regions of interest

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_model),                  INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_hybrid), INTENT(INOUT)           :: hybrid
    CHARACTER(LEN=3),                      INTENT(IN)              :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'calc_hybrid_solver_masks_basic_ROI'
    CHARACTER(LEN=256)                                             :: all_names_ROI, name_ROI
    INTEGER                                                        :: i
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                        :: poly_ROI
    INTEGER                                                        :: ti
    REAL(dp), DIMENSION(2)                                         :: p

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Go over all listed regions of interest
    all_names_ROI = C%choice_regions_of_interest

    DO WHILE (.TRUE.)

      ! Get the first region of interest from the list
      i = INDEX( all_names_ROI, '||')
      IF (i == 0) THEN
        ! There is only one left in the list
        name_ROI = TRIM( all_names_ROI)
        all_names_ROI = ''
      ELSE
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_TRIM( all_names_ROI))
      END IF


      SELECT CASE (region_name)
        CASE DEFAULT
          CALL crash('unknown region name "' // region_name // '"!')
        CASE ('NAM')
          ! North america

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              ! Don't need to do anything
              EXIT
            CASE ('Thwaites')
              ! Don't need to do anything
              EXIT
          END SELECT

        CASE ('EAS')
          ! Eurasia

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              ! Don't need to do anything
              EXIT
            CASE ('Thwaites')
              ! Don't need to do anything
              EXIT
          END SELECT

        CASE ('GRL')
          ! Greenland

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              ! Don't need to do anything
              EXIT
            CASE ('Thwaites')
              ! Don't need to do anything
              EXIT
          END SELECT

        CASE ('ANT')

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              CALL calc_polygon_Pine_Island_Glacier( poly_ROI)
            CASE ('Thwaites')
              CALL calc_polygon_Thwaites_Glacier( poly_ROI)
            CASE ('Tijn_test_ISMIP_HOM_A')
              CALL calc_polygon_Tijn_test_ISMIP_HOM_A( poly_ROI)
          END SELECT

      END SELECT

      ! Find all triangles that lie within this region of interest
      DO ti = mesh%ti1, mesh%ti2
        p = mesh%TriGC( ti,:)
        IF (is_in_polygon( poly_ROI, p)) THEN
          hybrid%mask_BPA_b(  ti) = .TRUE.
          hybrid%mask_DIVA_b( ti) = .FALSE.
        ELSE
          hybrid%mask_BPA_b(  ti) = .FALSE.
          hybrid%mask_DIVA_b( ti) = .TRUE.
        END IF
      END DO

      ! Clean up after yourself
      DEALLOCATE( poly_ROI)

      ! If no names are left, we are finished
      IF (all_names_ROI == '') EXIT

    END DO ! DO WHILE (.TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_hybrid_solver_masks_basic_ROI

  SUBROUTINE calc_hybrid_solver_masks_basic_read_from_file( mesh, ice, hybrid, region_name)
    ! Calculate the solving masks for the hybrid solver
    !
    ! Read the mask that determines where to solve the DIVA and where to solve the BPA from an external NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_model),                  INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_hybrid), INTENT(INOUT)           :: hybrid
    CHARACTER(LEN=3),                      INTENT(IN)              :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'calc_hybrid_solver_masks_basic_read_from_file'
    CHARACTER(LEN=256)                                             :: filename_hybrid_DIVA_BPA_mask
    INTEGER                                                        :: ncid
    TYPE(type_grid)                                                :: grid
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                        :: mask_int_grid_vec_partial
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE                        :: mask_int_grid
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE                        :: mask_grid
    INTEGER                                                        :: i,j
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                        :: poly_mult
    INTEGER                                                        :: ti
    REAL(dp), DIMENSION(2)                                         :: p

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE DEFAULT
        CALL crash('unknown region name "' // TRIM( region_name) // '"!')
      CASE ('NAM')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_NAM
      CASE ('EAS')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_EAS
      CASE ('GRL')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_GRL
      CASE ('ANT')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_ANT
    END SELECT

    ! Read grid from file
    CALL open_existing_netcdf_file_for_reading( filename_hybrid_DIVA_BPA_mask, ncid)
    CALL setup_xy_grid_from_file( filename_hybrid_DIVA_BPA_mask, ncid, grid)
    CALL close_netcdf_file( ncid)

    ! Read gridded mask from file
    ALLOCATE( mask_int_grid_vec_partial( grid%n1: grid%n2))
    CALL read_field_from_xy_file_2D_int( filename_hybrid_DIVA_BPA_mask, 'mask_BPA', mask_int_grid_vec_partial)

    ! Gather partial gridded data to the Master and broadcast the total field to all processes
    ALLOCATE( mask_int_grid( grid%nx, grid%ny))
    CALL gather_gridded_data_to_master_int_2D( grid, mask_int_grid_vec_partial, mask_int_grid)
    CALL MPI_BCAST( mask_int_grid, grid%nx * grid%ny, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Calculate logical mask (assumes data from file is integer 0 for FALSE and integer 1 for TRUE)
    ALLOCATE( mask_grid( grid%nx, grid%ny), source = .FALSE.)
    DO i = 1, grid%nx
    DO j = 1, grid%ny
      IF (mask_int_grid( i,j) == 1) mask_grid( i,j) = .TRUE.
    END DO
    END DO

    ! Calculate contour from gridded mask
    CALL calc_grid_mask_as_polygons( grid, mask_grid, poly_mult)

    ! Determine BPA solving masks on the mesh
    DO ti = mesh%ti1, mesh%ti2
      p = mesh%TriGC( ti,:)
      IF (is_in_polygons( poly_mult, p)) THEN
        hybrid%mask_BPA_b(  ti) = .TRUE.
        hybrid%mask_DIVA_b( ti) = .FALSE.
      ELSE
        hybrid%mask_BPA_b(  ti) = .FALSE.
        hybrid%mask_DIVA_b( ti) = .TRUE.
      END IF
    END DO

    ! Clean up after yourself
    DEALLOCATE( mask_int_grid_vec_partial)
    DEALLOCATE( mask_int_grid)
    DEALLOCATE( mask_grid)
    DEALLOCATE( poly_mult)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_hybrid_solver_masks_basic_read_from_file

! == Assemble and solve the linearised BPA

  SUBROUTINE solve_hybrid_DIVA_BPA_linearised( mesh, ice, hybrid, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Solve the linearised hybrid DIVA/BPA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_model),                  INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_hybrid), INTENT(INOUT)           :: hybrid
    INTEGER,  DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'solve_hybrid_DIVA_BPA_linearised'
    TYPE(type_sparse_matrix_CSR_dp)                                :: A_DIVA, A_BPA
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                        :: b_DIVA, b_BPA
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE                        :: tiuv2nh
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE                        :: tikuv2nh
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE                        :: nh2tiuv_tikuv
    INTEGER                                                        :: neq,i1,i2
    INTEGER                                                        :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                                :: A_combi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                        :: b_combi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                        :: uv_combi
    INTEGER                                                        :: neq_loc
    INTEGER                                                        :: row_nh,ti,k,uv,row_tiuv,row_tikuv,kk1,kk2,kk,col_tiuv,col_tikuv,tin,kn,uvn,col_nh
    REAL(dp)                                                       :: val, dzeta
    INTEGER                                                        :: nhu, nhv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Store the previous solution
    hybrid%u_bk_prev = hybrid%u_bk
    hybrid%v_bk_prev = hybrid%v_bk

    ! Calculate the stiffness matrix and load vector for the DIVA and the BPA
    CALL calc_masked_DIVA_stiffness_matrix_and_load_vector( mesh,      hybrid%DIVA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, hybrid%mask_DIVA_b, A_DIVA, b_DIVA)
    CALL calc_masked_BPA_stiffness_matrix_and_load_vector ( mesh, ice, hybrid%BPA , BC_prescr_mask_b                              , hybrid%mask_BPA_b , A_BPA , b_BPA )

    ! Calculate the "transition" solver masks
    CALL calc_hybrid_solver_masks_transition( mesh, hybrid, A_DIVA, A_BPA)

    ! Calculate combined DIVA/BPA translation tables
    CALL calc_hybrid_solver_translation_tables( mesh, hybrid, tiuv2nh, tikuv2nh, nh2tiuv_tikuv, neq, i1, i2)
    neq_loc = i2 + 1 - i1

  ! == Construct combined stiffness matrix and load vector
  ! ======================================================

    ! Initialise the stiffness matrix using the native UFEMISM CSR-matrix format

    ! Matrix size
    ncols           = neq      ! from
    ncols_loc       = neq_loc
    nrows           = neq      ! to
    nrows_loc       = neq_loc
    nnz_est_proc    = CEILING( 1.1_dp * REAL( A_DIVA%nnz + A_BPA%nnz, dp))

    CALL allocate_matrix_CSR_dist( A_combi, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector and the solution
    ALLOCATE( b_combi(  i1:i2))
    ALLOCATE( uv_combi( i1:i2))

    DO row_nh = i1, i2

      IF     (nh2tiuv_tikuv( row_nh,1) == 1) THEN
        ! This row represents a vertically averaged velocity

        ! This row in the combined matrix corresponds to this triangle and vertically averaged velocity component
        ti = nh2tiuv_tikuv( row_nh,2)
        uv = nh2tiuv_tikuv( row_nh,4)

        IF     (hybrid%mask_DIVA_b( ti)) THEN
          ! Copy the corresponding row from the DIVA stiffness matrix

          ! This triangle and vertically averaged velocity component correspond to this row in the DIVA matrix
          row_tiuv = mesh%tiuv2n( ti,uv)

          ! This row in the DIVA matrix contains these columns
          kk1 = A_DIVA%ptr( row_tiuv)
          kk2 = A_DIVA%ptr( row_tiuv+1)-1

          ! Loop over the columns of this row of the DIVA matrix
          DO kk = kk1, kk2
            ! This column index and coefficient of this entry in the DIVA matrix
            col_tiuv = A_DIVA%ind( kk)
            val      = A_DIVA%val( kk)
            ! This column in the DIVA matrix corresponds to this neighbouring triangle and vertically averaged velocity component
            tin = mesh%n2tiuv( col_tiuv,1)
            uvn = mesh%n2tiuv( col_tiuv,2)
            ! This neighbouring triangle and vertically averaged velocity component corresponds to this column in the combined matrix
            col_nh = tiuv2nh( tin,uvn)
            ! Add the coefficient from the DIVA matrix to the combined matrix
            CALL add_entry_CSR_dist( A_combi, row_nh, col_nh, val)
          END DO ! DO kk = kk1, kk2

          ! Copy the DIVA load vector
          b_combi( row_nh) = b_DIVA( row_tiuv)

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = hybrid%DIVA%u_vav_b( ti)

        ELSEIF (hybrid%mask_vav_from_BPA_b( ti)) THEN
          ! Define the vertically averaged velocities here from the 3-D BPA velocities
          !
          ! -u_vav + SUM_k [ u_3D( k) * dzeta( k)] = 0

          ! Add the coefficient of -1 for the vertically averaged velocity to the combined matrix
          CALL add_entry_CSR_dist( A_combi, row_nh, row_nh, -1._dp)

          ! Loop over the vertical column
          DO k = 1, mesh%nz

            ! Calculate the weight dzeta for the vertical average
            IF     (k == 1) THEN
              dzeta = mesh%zeta_stag( 1)
            ELSEIF (k == mesh%nz) THEN
              dzeta = 1._dp - mesh%zeta_stag( mesh%nz-1)
            ELSE
              dzeta = mesh%zeta_stag( k) - mesh%zeta_stag( k-1)
            END IF

            ! The 3-D velocity for this layer in this triangle corresponds to this column in the combined matrix
            col_nh = tikuv2nh( ti,k,uv)

            ! Add the coefficient to the combined matrix
            CALL add_entry_CSR_dist( A_combi, row_nh, col_nh, dzeta)

          END DO ! DO k = 1, mesh%nz

          ! The load vector is zero in this case
          b_combi( row_nh) = 0._dp

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = hybrid%DIVA%u_vav_b( ti)

        ELSE
          CALL crash('mask inconsistency; expected vertically averaged velocities, but both mask_DIVA_b and mask_vav_from_BPA_b are false!')
        END IF

      ELSEIF (nh2tiuv_tikuv( row_nh,1) == 2) THEN
        ! This row represents a 3-D velocity

        ! This row in the combined matrix corresponds to this triangle, layer, and 3-D velocity component
        ti = nh2tiuv_tikuv( row_nh,2)
        k  = nh2tiuv_tikuv( row_nh,3)
        uv = nh2tiuv_tikuv( row_nh,4)

        IF     (hybrid%mask_BPA_b( ti)) THEN
          ! Copy the corresponding row from the BPA stiffness matrix

          ! This triangle, layer, and 3-D velocity component correspond to this row in the BPA matrix
          row_tikuv = mesh%tikuv2n( ti,k,uv)

          ! This row in the BPA matrix contains these columns
          kk1 = A_BPA%ptr( row_tikuv)
          kk2 = A_BPA%ptr( row_tikuv+1)-1

          ! Loop over the columns of this row of the BPA matrix
          DO kk = kk1, kk2
            ! This column index and coefficient of this entry in the BPA matrix
            col_tikuv = A_BPA%ind( kk)
            val       = A_BPA%val( kk)
            ! This column in the BPA matrix corresponds to this neighbouring triangle, layer, and 3-D velocity component
            tin = mesh%n2tikuv( col_tikuv,1)
            kn  = mesh%n2tikuv( col_tikuv,2)
            uvn = mesh%n2tikuv( col_tikuv,3)
            ! This neighbouring triangle, layer, and 3-D velocity component corresponds to this column in the combined matrix
            col_nh = tikuv2nh( tin,kn,uvn)
            ! Add the coefficient from the DIVA matrix to the combined matrix
            CALL add_entry_CSR_dist( A_combi, row_nh, col_nh, val)
          END DO ! DO kk = kk1, kk2

          ! Copy the BPA load vector
          b_combi( row_nh) = b_BPA( row_tikuv)

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = hybrid%BPA%u_bk( ti,k)

        ELSEIF (hybrid%mask_3D_from_DIVA_b( ti)) THEN
          ! Define the 3-D velocities here from the vertically averaged DIVA velocities

          IF (C%choice_sliding_law == 'no_sliding') THEN
            ! Exception for the case of no sliding
            !
            ! According to Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34:
            !
            !   [1] u( z) = tau_bx * F1( z)
            !
            ! Also, according to the text just above Eq. 33:
            !
            !   [2] tau_bx = u_vav * beta_eff
            !
            ! Substituting [2] into [1] yields:
            !
            !   [3] u( z) = u_vav * beta_eff * F1( z)
            !
            ! This can be rearranged to read:
            !
            !   [4] u_vav * beta_eff * F1( z) - u( z) = 0

            ! u_vav term
            col_nh = tiuv2nh( ti,uv)
            val = hybrid%DIVA%beta_eff_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)
            CALL add_entry_CSR_dist( A_combi, row_nh, col_nh, val)

            ! u( z) term
            val = -1._dp
            CALL add_entry_CSR_dist( A_combi, row_nh, row_nh, val)

            ! The load vector is zero in this case
            b_combi( row_nh) = 0._dp

            ! Take the previous velocity solution as the initial guess
            uv_combi( row_nh) = hybrid%BPA%u_bk( ti,k)

          ELSE ! IF (C%choice_sliding_law == 'no_sliding') THEN
            ! The default case of finite sliding
            !
            ! According to Lipscomb et al., 2019, Eq. 29:
            !
            !   [1] u( z) = u_b * (1 + beta * F1( z))
            !
            ! Also, according to Eq. 32:
            !
            !   [2] u_b = u_vav / (1 + beta * F2( z=s))
            !
            ! Substituting [2] into [1] yields:
            !
            !   [3] u( z) = u_vav * (1 + beta * F1( z)) / (1 + beta * F2( z=s))
            !
            ! This can be rearranged to read:
            !
            !   [4] u_vav * (1 + beta * F1( z)) / (1 + beta * F2( z=s)) - u( z) = 0

            ! u_vav term
            col_nh = tiuv2nh( ti,uv)
            val = (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)) / &
                  (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F2_3D_b( ti,1))
            CALL add_entry_CSR_dist( A_combi, row_nh, col_nh, val)

            ! u( z) term
            val = -1._dp
            CALL add_entry_CSR_dist( A_combi, row_nh, row_nh, val)

            ! The load vector is zero in this case
            b_combi( row_nh) = 0._dp

            ! Take the previous velocity solution as the initial guess
            uv_combi( row_nh) = hybrid%BPA%u_bk( ti,k)

          END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

        ELSE
          CALL crash('mask inconsistency; expected 3-D velocities, but both mask_BPA_b and mask_3D_from_DIVA_b are false!')
        END IF

      ELSE
        CALL crash('nh2tiuv_tikuv( row_nh,1) = {int_01}, should be only 1 or 2!', int_01 = nh2tiuv_tikuv( row_nh,1))
      END IF ! IF     (nh2tiuv_tikuv( row_nh,1) == 1) THEN

    END DO ! DO row_nh = i1, i2

  ! == Solve the matrix equation
  ! ============================

    ! Use PETSc to solve the matrix equation
    CALL solve_matrix_equation_CSR_PETSc( A_combi, b_combi, uv_combi, hybrid%PETSc_rtol, hybrid%PETSc_abstol)

    ! Get velocities back from the combined vector
    DO ti = mesh%ti1, mesh%ti2

      IF (hybrid%mask_DIVA_b( ti)) THEN
        ! The DIVA was solved here

        ! Get vertically averaged DIVA velocities back from the combined vector
        nhu = tiuv2nh( ti,1)
        nhv = tiuv2nh( ti,2)
        hybrid%u_vav_b( ti) = uv_combi( nhu)
        hybrid%v_vav_b( ti) = uv_combi( nhv)

        ! Calculate 3-D velocities from the vertically averaged DIVA velocities
        IF (C%choice_sliding_law == 'no_sliding') THEN
          ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
          DO k = 1, mesh%nz
            hybrid%u_bk( ti,k) = hybrid%u_vav_b( ti) * hybrid%DIVA%beta_eff_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)
            hybrid%v_bk( ti,k) = hybrid%v_vav_b( ti) * hybrid%DIVA%beta_eff_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)
          END DO
        ELSE ! IF (C%choice_sliding_law == 'no_sliding') THEN
          ! Lipscomb et al., 2019, Eq. 29
          DO k = 1, mesh%nz
            hybrid%u_bk( ti,k) = hybrid%u_vav_b( ti) * &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)) / &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F2_3D_b( ti,1))
            hybrid%v_bk( ti,k) = hybrid%v_vav_b( ti) * &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)) / &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F2_3D_b( ti,1))
          END DO
        END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

      ELSEIF (hybrid%mask_BPA_b( ti)) THEN
        ! The BPA was solved here

        ! Get 3-D BPA velocities back from the combined vector
        DO k = 1, mesh%nz
          nhu = tikuv2nh( ti,k,1)
          nhv = tikuv2nh( ti,k,2)
          hybrid%u_bk( ti,k) = uv_combi( nhu)
          hybrid%v_bk( ti,k) = uv_combi( nhv)
        END DO

        ! Calculate vertically averaged velocities from the 3-D BPA velocities
        hybrid%u_vav_b( ti) = 0._dp
        hybrid%v_vav_b( ti) = 0._dp

        DO k = 1, mesh%nz
          IF     (k == 1) THEN
            dzeta = mesh%zeta_stag( 1)
          ELSEIF (k == mesh%nz) THEN
            dzeta = 1._dp - mesh%zeta_stag( mesh%nz-1)
          ELSE
            dzeta = mesh%zeta_stag( k) - mesh%zeta_stag( k-1)
          END IF
          hybrid%u_vav_b( ti) = hybrid%u_vav_b( ti) + dzeta * hybrid%u_bk( ti,k)
          hybrid%v_vav_b( ti) = hybrid%v_vav_b( ti) + dzeta * hybrid%v_bk( ti,k)
        END DO ! DO k = 1, mesh%nz

      ELSE
        ! Safety
        CALL crash('neither the DIVA nor the BPA was apparently solved here!')
      END IF ! IF (hybrid%mask_DIVA_b( ti)) THEN
    END DO ! DO row_nh = i1, i2

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( A_DIVA )
    CALL deallocate_matrix_CSR_dist( A_BPA  )
    CALL deallocate_matrix_CSR_dist( A_combi)
    DEALLOCATE( b_DIVA )
    DEALLOCATE( b_BPA  )
    DEALLOCATE( b_combi)
    DEALLOCATE( uv_combi)
    DEALLOCATE( tiuv2nh)
    DEALLOCATE( tikuv2nh)
    DEALLOCATE( nh2tiuv_tikuv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_hybrid_DIVA_BPA_linearised

  SUBROUTINE calc_masked_DIVA_stiffness_matrix_and_load_vector( mesh, DIVA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, mask_DIVA_b, A_DIVA, b_DIVA)
    ! Calculate the stiffness matrix for the masked DIVA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA),   INTENT(INOUT)           :: DIVA
    INTEGER,  DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    LOGICAL,  DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: mask_DIVA_b           ! T: solve the DIVA here, F: otherwise
    TYPE(type_sparse_matrix_CSR_dp),                   INTENT(OUT) :: A_DIVA
    REAL(dp), DIMENSION(:    ), ALLOCATABLE,           INTENT(OUT) :: b_DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'calc_masked_DIVA_stiffness_matrix_and_load_vector'
    INTEGER                                                        :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                                        :: row_tiuv,ti,uv

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
  ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * 2      ! from
    ncols_loc       = mesh%nTri_loc * 2
    nrows           = mesh%nTri     * 2      ! to
    nrows_loc       = mesh%nTri_loc * 2
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz * 4

    CALL allocate_matrix_CSR_dist( A_DIVA, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector
    ALLOCATE( b_DIVA( mesh%ti1*2-1: mesh%ti2*2))

  ! == Construct the stiffness matrix for the linearised DIVA
  ! ========================================================

    DO row_tiuv = A_DIVA%i1, A_DIVA%i2

      ti = mesh%n2tiuv( row_tiuv,1)
      uv = mesh%n2tiuv( row_tiuv,2)

      IF (BC_prescr_mask_b( ti) == 1) THEN
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        CALL add_entry_CSR_dist( A_DIVA, row_tiuv, row_tiuv, 1._dp)

        ! Load vector: prescribed velocity
        IF     (uv == 1) THEN
          b_DIVA( row_tiuv) = BC_prescr_u_b( ti)
        ELSEIF (uv == 2) THEN
          b_DIVA( row_tiuv) = BC_prescr_v_b( ti)
        ELSE
          CALL crash('uv can only be 1 or 2!')
        END IF

      ELSEIF (.NOT. mask_DIVA_b( ti)) THEN
        ! The BPA is solved here, not the DIVA

        CALL add_empty_row_CSR_dist( A_DIVA, row_tiuv)

      ELSEIF (mesh%TriBI( ti) == 1 .OR. mesh%TriBI( ti) == 2) THEN
        ! Northern domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_north( mesh, DIVA, A_DIVA, b_DIVA, row_tiuv)

      ELSEIF (mesh%TriBI( ti) == 3 .OR. mesh%TriBI( ti) == 4) THEN
        ! Eastern domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_east( mesh, DIVA, A_DIVA, b_DIVA, row_tiuv)

      ELSEIF (mesh%TriBI( ti) == 5 .OR. mesh%TriBI( ti) == 6) THEN
        ! Southern domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_south( mesh, DIVA, A_DIVA, b_DIVA, row_tiuv)

      ELSEIF (mesh%TriBI( ti) == 7 .OR. mesh%TriBI( ti) == 8) THEN
        ! Western domain border

        CALL calc_DIVA_stiffness_matrix_row_BC_west( mesh, DIVA, A_DIVA, b_DIVA, row_tiuv)

      ELSE
        ! No boundary conditions apply; solve the DIVA

        IF (C%do_include_SSADIVA_crossterms) THEN
          ! Calculate matrix coefficients for the full DIVA
          CALL calc_DIVA_stiffness_matrix_row_free( mesh, DIVA, A_DIVA, b_DIVA, row_tiuv)
        ELSE
          ! Calculate matrix coefficients for the DIVA sans the gradients of the effective viscosity (the "cross-terms")
          CALL calc_DIVA_sans_stiffness_matrix_row_free( mesh, DIVA, A_DIVA, b_DIVA, row_tiuv)
        END IF

      END IF

    END DO ! DO row_tiuv = A_DIVA%i1, A_DIVA%i2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_masked_DIVA_stiffness_matrix_and_load_vector

  SUBROUTINE calc_masked_BPA_stiffness_matrix_and_load_vector( mesh, ice, BPA, BC_prescr_mask_b, mask_BPA_b, A_BPA, b_BPA)
    ! Calculate the stiffness matrix for the masked BPA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_model),                  INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_BPA),    INTENT(INOUT)           :: BPA
    INTEGER,  DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    LOGICAL,  DIMENSION(mesh%ti1:mesh%ti2),            INTENT(IN)  :: mask_BPA_b            ! T: solve the BPA here, F: otherwise
    TYPE(type_sparse_matrix_CSR_dp),                   INTENT(OUT) :: A_BPA
    REAL(dp), DIMENSION(:    ), ALLOCATABLE,           INTENT(OUT) :: b_BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'calc_masked_BPA_stiffness_matrix_and_load_vector'
    INTEGER                                                        :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                                        :: row_tikuv,ti,k,uv

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
  ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * mesh%nz * 2      ! from
    ncols_loc       = mesh%nTri_loc * mesh%nz * 2
    nrows           = mesh%nTri     * mesh%nz * 2      ! to
    nrows_loc       = mesh%nTri_loc * mesh%nz * 2
    nnz_est_proc    = mesh%M2_ddx_bk_bk%nnz   * 4

    CALL allocate_matrix_CSR_dist( A_BPA, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector
    ALLOCATE( b_BPA( A_BPA%i1:A_BPA%i2))

  ! == Construct the stiffness matrix for the linearised BPA
  ! ========================================================

    DO row_tikuv = A_BPA%i1, A_BPA%i2

      ti = mesh%n2tikuv( row_tikuv,1)
      k  = mesh%n2tikuv( row_tikuv,2)
      uv = mesh%n2tikuv( row_tikuv,3)

      IF (BC_prescr_mask_b( ti) == 1 .OR. .NOT. mask_BPA_b( ti)) THEN
        ! The DIVA is solved here, not the BPA

        CALL add_empty_row_CSR_dist( A_BPA, row_tikuv)

      ELSEIF (mesh%TriBI( ti) == 1 .OR. mesh%TriBI( ti) == 2) THEN
        ! Northern domain border

        CALL calc_BPA_stiffness_matrix_row_BC_north( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      ELSEIF (mesh%TriBI( ti) == 3 .OR. mesh%TriBI( ti) == 4) THEN
        ! Eastern domain border

        CALL calc_BPA_stiffness_matrix_row_BC_east( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      ELSEIF (mesh%TriBI( ti) == 5 .OR. mesh%TriBI( ti) == 6) THEN
        ! Southern domain border

        CALL calc_BPA_stiffness_matrix_row_BC_south( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      ELSEIF (mesh%TriBI( ti) == 7 .OR. mesh%TriBI( ti) == 8) THEN
        ! Western domain border

        CALL calc_BPA_stiffness_matrix_row_BC_west( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      ELSEIF (k == 1) THEN
        ! Ice surface

        CALL calc_BPA_stiffness_matrix_row_BC_surf( mesh, ice, BPA, A_BPA, b_BPA, row_tikuv)

      ELSEIF (k == mesh%nz) THEN
        ! Ice base

        CALL calc_BPA_stiffness_matrix_row_BC_base( mesh, ice, BPA, A_BPA, b_BPA, row_tikuv)

      ELSE
        ! No boundary conditions apply; solve the BPA

        CALL calc_BPA_stiffness_matrix_row_free( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      END IF

    END DO ! DO row_tikuv = A_BPA%i1, A_BPA%i2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_masked_BPA_stiffness_matrix_and_load_vector

  SUBROUTINE calc_hybrid_solver_masks_transition( mesh, hybrid, A_DIVA, A_BPA)
    ! Calculate the "transition" solver masks

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_hybrid), INTENT(INOUT)           :: hybrid
    TYPE(type_sparse_matrix_CSR_dp),                   INTENT(IN)  :: A_DIVA
    TYPE(type_sparse_matrix_CSR_dp),                   INTENT(IN)  :: A_BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'calc_hybrid_solver_masks_transition'
    LOGICAL,  DIMENSION( mesh%nTri)                                :: mask_DIVA_halo_b, mask_BPA_halo_b
    INTEGER                                                        :: row_tiuv,ti,kk1,kk2,kk,col_tjuv,tj
    INTEGER                                                        :: row_tikuv,col_tjkuv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Mark all triangles where the DIVA solver needs vertically averaged velocities
    mask_DIVA_halo_b = .FALSE.
    DO row_tiuv = A_DIVA%i1, A_DIVA%i2

      ti = mesh%n2tiuv( row_tiuv,1)

      kk1 = A_DIVA%ptr( row_tiuv)
      kk2 = A_DIVA%ptr( row_tiuv+1) - 1

      IF (kk2 >= kk1) mask_DIVA_halo_b( ti) = .TRUE.

      DO kk = kk1, kk2

        col_tjuv = A_DIVA%ind( kk)
        tj = mesh%n2tiuv( col_tjuv,1)

        mask_DIVA_halo_b( tj) = .TRUE.

      END DO ! DO kk = kk1, kk2

    END DO ! DO row_tiuv = A_DIVA%i1, A_DIVA%i2
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, mask_DIVA_halo_b, mesh%nTri, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

    ! Mark all triangles where the BPA solver needs vertically averaged velocities
    mask_BPA_halo_b = .FALSE.
    DO row_tikuv = A_BPA%i1, A_BPA%i2

      ti = mesh%n2tikuv( row_tikuv,1)

      kk1 = A_BPA%ptr( row_tikuv)
      kk2 = A_BPA%ptr( row_tikuv+1) - 1

      IF (kk2 >= kk1) mask_BPA_halo_b( ti) = .TRUE.

      DO kk = kk1, kk2

        col_tjkuv = A_BPA%ind( kk)
        tj = mesh%n2tikuv( col_tjkuv,1)

        mask_BPA_halo_b( tj) = .TRUE.

      END DO ! DO kk = kk1, kk2

    END DO ! DO row_tiuv = A_BPA%i1, A_BPA%i2
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, mask_BPA_halo_b, mesh%nTri, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

    ! Mark all triangles where the DIVA is solved, but a nearby BPA triangle needs 3-D velocities
    hybrid%mask_3D_from_DIVA_b = .FALSE.
    DO ti = mesh%ti1, mesh%ti2
      IF (hybrid%mask_DIVA_b( ti) .AND. mask_BPA_halo_b( ti)) THEN
        hybrid%mask_3D_from_DIVA_b( ti) = .TRUE.
      END IF
    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Mark all triangles where the BPA is solved, but a nearby DIVA triangle needs vertically averaged velocities
    hybrid%mask_vav_from_BPA_b = .FALSE.
    DO ti = mesh%ti1, mesh%ti2
      IF (hybrid%mask_BPA_b( ti) .AND. mask_DIVA_halo_b( ti)) THEN
        hybrid%mask_vav_from_BPA_b( ti) = .TRUE.
      END IF
    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_hybrid_solver_masks_transition

  SUBROUTINE calc_hybrid_solver_translation_tables( mesh, hybrid, tiuv2nh, tikuv2nh, nh2tiuv_tikuv, neq, i1, i2)
    ! Calculate combined DIVA/BPA translation tables

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_hybrid), INTENT(IN)              :: hybrid
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)           :: tiuv2nh
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT)           :: tikuv2nh
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)           :: nh2tiuv_tikuv
    INTEGER,                                 INTENT(OUT)           :: neq,i1,i2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                  :: routine_name = 'calc_hybrid_solver_translation_tables'
    LOGICAL,  DIMENSION(mesh%nTri)                                 :: mask_DIVA_b_tot
    LOGICAL,  DIMENSION(mesh%nTri)                                 :: mask_BPA_b_tot
    LOGICAL,  DIMENSION(mesh%nTri)                                 :: mask_3D_from_DIVA_b_tot
    LOGICAL,  DIMENSION(mesh%nTri)                                 :: mask_vav_from_BPA_b_tot
    INTEGER                                                        :: ti,k,uv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather global masks
    CALL gather_to_all_logical_1D( hybrid%mask_DIVA_b        , mask_DIVA_b_tot        )
    CALL gather_to_all_logical_1D( hybrid%mask_BPA_b         , mask_BPA_b_tot         )
    CALL gather_to_all_logical_1D( hybrid%mask_3D_from_DIVA_b, mask_3D_from_DIVA_b_tot)
    CALL gather_to_all_logical_1D( hybrid%mask_vav_from_BPA_b, mask_vav_from_BPA_b_tot)

    ! Allocate memory
    ALLOCATE( tiuv2nh      ( mesh%nTri          ,  2   ))
    ALLOCATE( tikuv2nh     ( mesh%nTri,  mesh%nz,  2   ))
    ALLOCATE( nh2tiuv_tikuv( mesh%nTri * mesh%nz * 2, 4))

    neq = 0
    i1  = 0
    i2  = 0

    DO ti = 1, mesh%nTri

      IF (ti == mesh%ti1) THEN
        i1 = neq + 1
      END IF

      IF (mask_DIVA_b_tot( ti) .OR. mask_vav_from_BPA_b_tot( ti)) THEN
        ! Vertically averaged velocities must be defined here
        DO uv = 1, 2
          neq = neq + 1
          tiuv2nh( ti,uv) = neq
          nh2tiuv_tikuv( neq,:) = [1,ti,0,uv]
        END DO ! DO uv = 1, 2
      END IF

      IF (mask_BPA_b_tot( ti) .OR. mask_3D_from_DIVA_b_tot( ti)) THEN
        ! 3-D velocities must be defined here
        DO k = 1, mesh%nz
        DO uv = 1, 2
          neq = neq + 1
          tikuv2nh( ti,k,uv) = neq
          nh2tiuv_tikuv( neq,:) = [2,ti,k,uv]
        END DO
        END DO
      END IF

      IF (ti == mesh%ti2) THEN
        i2 = neq
      END IF

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_hybrid_solver_translation_tables

! == Interfaces with DIVA-specific routines

  SUBROUTINE calc_driving_stress_DIVA( mesh, ice, DIVA)
    ! Calculate the driving stress

    USE ice_velocity_DIVA, ONLY: calc_driving_stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_driving_stress_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_driving_stress( mesh, ice, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress_DIVA

  SUBROUTINE calc_horizontal_strain_rates_DIVA( mesh, DIVA)
    ! Calculate the horizontal strain rates

    USE ice_velocity_DIVA, ONLY: calc_horizontal_strain_rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_horizontal_strain_rates_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_horizontal_strain_rates( mesh, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_horizontal_strain_rates_DIVA

  SUBROUTINE calc_vertical_shear_strain_rates_DIVA( mesh, DIVA)
    ! Calculate the vertical shear strain rates

    USE ice_velocity_DIVA, ONLY: calc_vertical_shear_strain_rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_vertical_shear_strain_rates_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_vertical_shear_strain_rates( mesh, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_shear_strain_rates_DIVA

  SUBROUTINE calc_effective_viscosity_DIVA( mesh, ice, DIVA, Glens_flow_law_epsilon_sq_0_applied)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    USE ice_velocity_DIVA, ONLY: calc_effective_viscosity

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA
    REAL(dp),                            INTENT(IN)              :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_viscosity_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_effective_viscosity( mesh, ice, DIVA, Glens_flow_law_epsilon_sq_0_applied)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity_DIVA

  SUBROUTINE calc_F_integrals_DIVA( mesh, ice, DIVA)
    ! Calculate the F-integrals on the a-grid (Lipscomb et al. (2019), Eq. 30)

    USE ice_velocity_DIVA, ONLY: calc_F_integrals

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_F_integrals_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_F_integrals( mesh, ice, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_F_integrals_DIVA

  SUBROUTINE calc_effective_basal_friction_coefficient_DIVA( mesh, ice, DIVA)
    ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)

    USE ice_velocity_DIVA, ONLY: calc_effective_basal_friction_coefficient

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_basal_friction_coefficient_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_effective_basal_friction_coefficient( mesh, ice, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_basal_friction_coefficient_DIVA

  SUBROUTINE calc_basal_shear_stress_DIVA( mesh, DIVA)
    ! Calculate the basal shear stress (Lipscomb et al., 2019, just above Eq. 33)

    USE ice_velocity_DIVA, ONLY: calc_basal_shear_stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_basal_shear_stress_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_basal_shear_stress( mesh, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_shear_stress_DIVA

  SUBROUTINE calc_basal_velocities_DIVA( mesh, DIVA)
    ! Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)

    USE ice_velocity_DIVA, ONLY: calc_basal_velocities

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_DIVA), INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_basal_velocities_DIVA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_basal_velocities( mesh, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_velocities_DIVA

! == Interfaces with BPA-specific routines

  SUBROUTINE calc_driving_stress_BPA( mesh, ice, BPA)
    ! Calculate the driving stress

    USE ice_velocity_BPA, ONLY: calc_driving_stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_driving_stress_BPA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_driving_stress( mesh, ice, BPA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress_BPA

  SUBROUTINE calc_strain_rates_BPA( mesh, BPA)
    ! Calculate the strain rates

    USE ice_velocity_BPA, ONLY: calc_strain_rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_strain_rates_BPA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_strain_rates( mesh, BPA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_rates_BPA

  SUBROUTINE calc_effective_viscosity_BPA( mesh, ice, BPA, Glens_flow_law_epsilon_sq_0_applied)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    USE ice_velocity_BPA, ONLY: calc_effective_viscosity

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA
    REAL(dp),                            INTENT(IN)              :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_viscosity_BPA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_effective_viscosity( mesh, ice, BPA, Glens_flow_law_epsilon_sq_0_applied)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity_BPA

  SUBROUTINE calc_applied_basal_friction_coefficient_BPA( mesh, ice, BPA)
    ! Calculate the applied basal friction coefficient beta_b

    USE ice_velocity_BPA, ONLY: calc_applied_basal_friction_coefficient

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_applied_basal_friction_coefficient_BPA'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_applied_basal_friction_coefficient( mesh, ice, BPA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_applied_basal_friction_coefficient_BPA

! == Some useful tools for improving numerical stability of the viscosity iteration

  SUBROUTINE relax_viscosity_iterations( mesh, hybrid, visc_it_relax)
    ! Reduce the change between velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_hybrid),  INTENT(INOUT)           :: hybrid
    REAL(dp),                               INTENT(IN)              :: visc_it_relax

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'relax_viscosity_iterations'
    INTEGER                                                         :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      hybrid%u_bk( ti,:) = (visc_it_relax * hybrid%u_bk( ti,:)) + ((1._dp - visc_it_relax) * hybrid%u_bk_prev( ti,:))
      hybrid%v_bk( ti,:) = (visc_it_relax * hybrid%v_bk( ti,:)) + ((1._dp - visc_it_relax) * hybrid%v_bk_prev( ti,:))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_viscosity_iterations

  SUBROUTINE calc_visc_iter_UV_resid( mesh, hybrid, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    IMPLICIT NONE

    TYPE(type_mesh),                        INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_hybrid),  INTENT(INOUT)           :: hybrid
    REAL(dp),                               INTENT(OUT)             :: resid_UV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'calc_visc_iter_UV_resid'
    INTEGER                                                         :: ierr
    INTEGER                                                         :: ti,k
    REAL(dp)                                                        :: res1, res2

    ! Add routine to path
    CALL init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, mesh%nz

      res1 = res1 + (hybrid%u_bk( ti,k) - hybrid%u_bk_prev( ti,k))**2
      res1 = res1 + (hybrid%v_bk( ti,k) - hybrid%v_bk_prev( ti,k))**2

      res2 = res2 + (hybrid%u_bk( ti,k) + hybrid%u_bk_prev( ti,k))**2
      res2 = res2 + (hybrid%v_bk( ti,k) + hybrid%v_bk_prev( ti,k))**2

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

  SUBROUTINE apply_velocity_limits( mesh, hybrid)
    ! Limit velocities for improved stability

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_hybrid),  INTENT(INOUT)           :: hybrid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                   :: routine_name = 'apply_velocity_limits'
    INTEGER                                                         :: ti,k
    REAL(dp)                                                        :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, mesh%nz

      ! Calculate absolute speed
      uabs = SQRT( hybrid%u_bk( ti,k)**2 + hybrid%v_bk( ti,k)**2)

      ! Reduce velocities if neceBPAry
      IF (uabs > C%vel_max) THEN
        hybrid%u_bk( ti,k) = hybrid%u_bk( ti,k) * C%vel_max / uabs
        hybrid%v_bk( ti,k) = hybrid%v_bk( ti,k) * C%vel_max / uabs
      END IF

    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits

! == Initialisation

  SUBROUTINE allocate_hybrid_DIVA_BPA_solver( mesh, hybrid)
    ! Allocate memory the hybrid DIVA/BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                       INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_hybrid), INTENT(OUT)   :: hybrid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_hybrid_DIVA_BPA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Solution
    ALLOCATE( hybrid%u_vav_b            ( mesh%ti1:mesh%ti2        ))           ! [m yr^-1] 2-D horizontal ice velocity
    hybrid%u_vav_b = 0._dp
    ALLOCATE( hybrid%v_vav_b            ( mesh%ti1:mesh%ti2        ))
    hybrid%v_vav_b = 0._dp
    ALLOCATE( hybrid%u_bk               ( mesh%ti1:mesh%ti2,mesh%nz))           ! [m yr^-1] 2-D horizontal ice velocity
    hybrid%u_bk = 0._dp
    ALLOCATE( hybrid%v_bk               ( mesh%ti1:mesh%ti2,mesh%nz))
    hybrid%v_bk = 0._dp

    ! Separate DIVA/BPA solvers
    CALL allocate_DIVA_solver( mesh, hybrid%DIVA)
    CALL allocate_BPA_solver(  mesh, hybrid%BPA )

    ! Solver masks
    ALLOCATE( hybrid%mask_DIVA_b        ( mesh%ti1:mesh%ti2))                   ! T: solve the DIVA here, F: otherwise
    hybrid%mask_DIVA_b = .FALSE.
    ALLOCATE( hybrid%mask_BPA_b         ( mesh%ti1:mesh%ti2))                   ! T: solve the BPA  here, F: otherwise
    hybrid%mask_BPA_b = .FALSE.
    ALLOCATE( hybrid%mask_3D_from_DIVA_b( mesh%ti1:mesh%ti2))                   ! T: calculate 3-D velocities from the vertically averaged DIVA solution here, F: otherwise
    hybrid%mask_3D_from_DIVA_b = .FALSE.
    ALLOCATE( hybrid%mask_vav_from_BPA_b( mesh%ti1:mesh%ti2))                   ! T: calculate vertically averaged velocities from the 3-D BPA  solution here, F: otherwise
    hybrid%mask_vav_from_BPA_b = .FALSE.

    ! Intermediate data fields
    ALLOCATE( hybrid%u_bk_prev          ( mesh%nTri,mesh%nz))                   ! Velocity solution from previous viscosity iteration
    hybrid%u_bk_prev = 0._dp
    ALLOCATE( hybrid%v_bk_prev          ( mesh%nTri,mesh%nz))
    hybrid%v_bk_prev = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_hybrid_DIVA_BPA_solver

END MODULE ice_velocity_hybrid_DIVA_BPA
