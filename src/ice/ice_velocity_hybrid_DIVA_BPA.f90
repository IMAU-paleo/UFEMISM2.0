MODULE ice_velocity_hybrid_DIVA_BPA

  ! Routines for calculating ice velocities using the hybrid DIVA/BPA

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
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
  USE mesh_remapping                                         , ONLY: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  USE ice_velocity_DIVA                                      , ONLY: allocate_DIVA_solver, remap_DIVA_solver
  USE ice_velocity_BPA                                       , ONLY: allocate_BPA_solver , remap_BPA_solver
  USE mesh_operators                                         , ONLY: calc_3D_matrix_operators_mesh, map_a_b_2D, map_b_a_2D, map_b_a_3D, map_a_b_3D
  USE mesh_refinement                                        , ONLY: calc_polygon_Pine_Island_Glacier, calc_polygon_Thwaites_Glacier, &
                                                                     calc_polygon_Tijn_test_ISMIP_HOM_A
  USE math_utilities                                         , ONLY: is_in_polygon
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_logical_1D
  USE ice_model_utilities                                    , ONLY: calc_zeta_gradients

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
      hybrid%u_bk = 0._dp
      hybrid%v_bk = 0._dp
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
      hybrid%u_bk = 0._dp
      hybrid%v_bk = 0._dp
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
      uv_min = MINVAL( hybrid%u_bk)
      uv_max = MAXVAL( hybrid%u_bk)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      IF (par%master) WRITE(0,*) '    hybrid DIVA/BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

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

      ! DENK DROOM
      has_converged = .TRUE.
      CALL warning('DENK DROM - ending viscosity iteration after 1 loop!')

    END DO viscosity_iteration

    ! Clean up after yourself
    DEALLOCATE( BC_prescr_mask_b_applied)
    DEALLOCATE( BC_prescr_u_b_applied   )
    DEALLOCATE( BC_prescr_v_b_applied   )

    ! DENK DROM
    CALL crash('whoopsiedaisy!')

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
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE              :: u_ak
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE              :: v_ak

    ! Add routine to path
    CALL init_routine( routine_name)

  ! Remap the fields that are re-used during the viscosity iteration
  ! ================================================================

    ! Allocate memory for velocities on the a-grid (vertices)
    ALLOCATE( u_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    ALLOCATE( v_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    CALL map_b_a_3D( mesh_old, hybrid%u_bk   , u_ak   )
    CALL map_b_a_3D( mesh_old, hybrid%v_bk   , v_ak   )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, u_ak   , '2nd_order_conservative')
    CALL map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, v_ak   , '2nd_order_conservative')

    ! Reallocate memory for the data on the triangles
    CALL reallocate_bounds( hybrid%u_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    CALL reallocate_bounds( hybrid%v_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    CALL map_a_b_3D( mesh_new, u_ak   , hybrid%u_bk   )
    CALL map_a_b_3D( mesh_new, v_ak   , hybrid%v_bk   )

    ! Clean up after yourself
    DEALLOCATE( u_ak   )
    DEALLOCATE( v_ak   )

  ! Remap data of the separate DIVA and BPA solvers
  ! ===============================================

    CALL remap_DIVA_solver( mesh_old, mesh_new, hybrid%DIVA)
    CALL remap_BPA_solver(  mesh_old, mesh_new, hybrid%BPA )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_hybrid_DIVA_BPA_solver

! == Masks and translation tables for the hybrid solver

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

    ! Add routine to path
    CALL init_routine( routine_name)

    SELECT CASE (C%choice_hybrid_DIVA_BPA_mask)
      CASE DEFAULT
        CALL crash('unknown choice_hybrid_DIVA_BPA_mask "' // TRIM( C%choice_hybrid_DIVA_BPA_mask) // '"!')
      CASE ('ROI')
        CALL calc_hybrid_solver_masks_basic_ROI( mesh, ice, hybrid, region_name)
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

    ! DENK DROM
    CALL save_variable_as_netcdf_logical_1D( hybrid%mask_DIVA_b, 'mask_DIVA_b')
    CALL save_variable_as_netcdf_logical_1D( hybrid%mask_BPA_b , 'mask_BPA_b' )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_hybrid_solver_masks_basic_ROI

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

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_hybrid_DIVA_BPA_linearised

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
