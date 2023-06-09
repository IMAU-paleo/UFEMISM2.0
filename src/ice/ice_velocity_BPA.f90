MODULE ice_velocity_BPA

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (BPA)

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
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_BPA
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate_clean
  USE mesh_operators                                         , ONLY: map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_3D, ddy_b_a_3D
  USE mesh_zeta                                              , ONLY: vertical_average
  USE sliding_laws                                           , ONLY: calc_basal_friction_coefficient
  USE mesh_utilities                                         , ONLY: find_ti_copy_ISMIP_HOM_periodic
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, close_netcdf_file, open_existing_netcdf_file_for_writing
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, add_time_dimension_to_file, &
                                                                     add_zeta_dimension_to_file, add_field_mesh_dp_3D_b, write_time_to_file, write_to_field_multopt_mesh_dp_3D_b
  USE netcdf_input                                           , ONLY: read_field_from_mesh_file_3D_b

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
    ! Calculate ice velocities by solving the Shallow Ice Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_BPA'
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE                      :: BC_prescr_mask_bk_applied
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      :: BC_prescr_u_bk_applied
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                      :: BC_prescr_v_bk_applied
    INTEGER                                                      :: viscosity_iteration_i
    LOGICAL                                                      :: has_converged
    REAL(dp)                                                     :: resid_UV
    REAL(dp)                                                     :: uv_min, uv_max

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there is no grounded ice, or no sliding, no need to solve the BPA
    IF ((.NOT. ANY( ice%mask_sheet)) .OR. C%choice_sliding_law == 'no_sliding') THEN
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

    ! Calculate the driving stress
    CALL calc_driving_stress( mesh, ice, BPA)

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the strain rates for the current velocity solution
      CALL calc_strain_rates( mesh, BPA)

!      ! Calculate the effective viscosity for the current velocity solution
!      CALL calc_effective_viscosity( mesh, ice, BPA)
!
!      ! Calculate the basal friction coefficient betab for the current velocity solution
!      CALL calc_applied_basal_friction_coefficient( mesh, ice, BPA)
!
!      ! Solve the linearised BPA to calculate a new velocity solution
!      CALL solve_BPA_linearised( mesh, BPA, BC_prescr_mask_bk_applied, BC_prescr_u_bk_applied, BC_prescr_v_bk_applied)
!
!      ! Limit velocities for improved stability
!      CALL apply_velocity_limits( mesh, BPA)
!
!      ! Reduce the change between velocity solutions
!      CALL relax_viscosity_iterations( mesh, BPA)
!
!      ! Calculate the L2-norm of the two consecutive velocity solutions
!      CALL calc_visc_iter_UV_resid( mesh, BPA, resid_UV)

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
         CALL warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
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

  SUBROUTINE remap_BPA_solver( mesh_old, mesh_new, ice, BPA)
    ! Remap the BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT) :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_BPA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BPA_solver

! == Assemble and solve the linearised BPA

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
      BPA%tau_dx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      BPA%tau_dy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    END DO

    ! Clean up after yourself
    DEALLOCATE( Hi_b    )
    DEALLOCATE( dHs_dx_b)
    DEALLOCATE( dHs_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress

  SUBROUTINE calc_strain_rates( mesh, BPA)
    ! Calculate the strain rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_strain_rates'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the strain rates
!    CALL ddx_b_a_3D( mesh, BPA%u_bk, BPA%du_dx_ak)
!    CALL ddy_b_a_3D( mesh, BPA%u_bk, BPA%du_dy_ak)
!    CALL ddz_b_a_3D( mesh, BPA%u_bk, BPA%du_dz_ak)
!    CALL ddx_b_a_3D( mesh, BPA%v_bk, BPA%dv_dx_ak)
!    CALL ddy_b_a_3D( mesh, BPA%v_bk, BPA%dv_dy_ak)
!    CALL ddz_b_a_3D( mesh, BPA%v_bk, BPA%dv_dz_ak)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_rates

! == Some useful tools for improving numerical stability of the viscosity iteration

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
    CHARACTER(LEN=256)                                 :: filename
    REAL(dp)                                           :: timeframe

    ! Add routine to path
    CALL init_routine( routine_name)

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
    ALLOCATE( BPA%beta_b_b(      mesh%ti1:mesh%ti2))                           ! Basal friction coefficient (tau_b = u * beta_b)
    BPA%beta_b_b = 0._dp
    ALLOCATE( BPA%tau_dx_b(      mesh%ti1:mesh%ti2))                           ! Driving stress
    BPA%tau_dx_b = 0._dp
    ALLOCATE( BPA%tau_dy_b(      mesh%ti1:mesh%ti2))
    BPA%tau_dy_b = 0._dp
    ALLOCATE( BPA%u_bk_prev(     mesh%ti1:mesh%ti2,mesh%nz))                   ! Velocity solution from previous viscosity iteration
    BPA%u_bk_prev = 0._dp
    ALLOCATE( BPA%v_bk_prev(     mesh%ti1:mesh%ti2,mesh%nz))
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
    TYPE(type_ice_velocity_solver_BPA),  INTENT(INOUT)           :: BPA
    REAL(dp),                            INTENT(IN)              :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'write_to_restart_file_BPA'
    INTEGER                                                      :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_ice_velocity_BPA'
    CALL generate_filename_XXXXXdotnc( filename_base, BPA%restart_filename)

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( BPA%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( BPA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( BPA%restart_filename, ncid)

    ! Add a zeta dimension to the file
    CALL add_zeta_dimension_to_file( BPA%restart_filename, ncid, mesh%zeta)

    ! Add the velocity fields to the file
    CALL add_field_mesh_dp_3D_b( BPA%restart_filename, ncid, 'u_bk')
    CALL add_field_mesh_dp_3D_b( BPA%restart_filename, ncid, 'v_bk')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_BPA

END MODULE ice_velocity_BPA
