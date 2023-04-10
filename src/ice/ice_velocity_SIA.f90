MODULE ice_velocity_SIA

  ! Routines for calculating ice velocities using the Shallow Ice Approximation (SIA)

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_SIA
  USE mesh_operators                                         , ONLY: map_a_b_2D, map_a_b_3D, ddx_a_a_2D, ddy_a_a_2D, ddx_a_b_2D, ddy_a_b_2D
  USE mesh_zeta                                              , ONLY: integrate_from_zeta_is_one_to_zeta_is_zetap
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate_clean

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE initialise_SIA_solver( mesh, SIA)
    ! Initialise the SIA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_SIA),  INTENT(OUT)   :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SIA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    CALL allocate_SIA_solver( mesh, SIA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SIA_solver

  SUBROUTINE allocate_SIA_solver( mesh, SIA)
    ! Allocate memory the SIA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_velocity_solver_SIA),  INTENT(OUT)   :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_SIA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory

    ! Solution
    ALLOCATE( SIA%u_3D_b(     mesh%nTri_loc, mesh%nz), source = 0._dp)
    ALLOCATE( SIA%v_3D_b(     mesh%nTri_loc, mesh%nz), source = 0._dp)
    ALLOCATE( SIA%du_dz_3D_a( mesh%nV_loc  , mesh%nz), source = 0._dp)
    ALLOCATE( SIA%dv_dz_3D_a( mesh%nV_loc  , mesh%nz), source = 0._dp)

    ! Intermediate data fields
    ALLOCATE( SIA%D_3D_b(     mesh%nTri_loc, mesh%nz), source = 0._dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_SIA_solver

  SUBROUTINE solve_SIA( mesh, ice, SIA)
    ! Calculate ice velocities by solving the Shallow Ice Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ice_velocity_solver_SIA),  INTENT(INOUT) :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_SIA'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            ::  Hi_b,  Hs_b,  dHs_dx_a,  dHs_dy_a,  dHs_dx_b,  dHs_dy_b
    INTEGER                                            :: wHi_b, wHs_b, wdHs_dx_a, wdHs_dy_a, wdHs_dx_b, wdHs_dy_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            ::  A_flow_3D_b
    INTEGER                                            :: wA_flow_3D_b
    INTEGER                                            :: vi,ti,k
    REAL(dp)                                           :: abs_grad_Hs
    REAL(dp), DIMENSION(mesh%nz)                       :: z, int_A_hminzetan

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( Hi_b(        mesh%nTri_loc         ), source = 0._dp)
    ALLOCATE( Hs_b(        mesh%nTri_loc         ), source = 0._dp)
    ALLOCATE( dHs_dx_a(    mesh%nV_loc           ), source = 0._dp)
    ALLOCATE( dHs_dy_a(    mesh%nV_loc           ), source = 0._dp)
    ALLOCATE( dHs_dx_b(    mesh%nTri_loc         ), source = 0._dp)
    ALLOCATE( dHs_dy_b(    mesh%nTri_loc         ), source = 0._dp)
    ALLOCATE( A_flow_3D_b( mesh%nTri_loc, mesh%nz), source = 0._dp)

    ! Calculate ice thickness, surface elevation, surface slopes, and ice flow factor on the b-grid
    CALL map_a_b_2D( mesh, ice%Hi       , Hi_b       )
    CALL map_a_b_2D( mesh, ice%Hs       , Hs_b       )
    CALL ddx_a_a_2D( mesh, ice%Hs       , dHs_dx_a   )
    CALL ddy_a_a_2D( mesh, ice%Hs       , dHs_dy_a   )
    CALL ddx_a_b_2D( mesh, ice%Hs       , dHs_dx_b   )
    CALL ddy_a_b_2D( mesh, ice%Hs       , dHs_dy_b   )
    CALL map_a_b_3D( mesh, ice%A_flow_3D, A_flow_3D_b)

    ! Calculate velocities and strain rates according to the analytical solution of the SIA:
    ! (see also Bueler and Brown, 2009, Eqs. 12-13)
    !
    !   D( z) = -2 (rho g)^n (abs(grad H))^(n-1) int_b_z( A(T*) (h - zeta)^n ) dzeta
    !   u( z) = dh/dx D( z)
    !   v( z) = dh/dy D( z)
    !
    !   du/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dx
    !   dv/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dy

    ! Calculate velocities
    DO ti = 1, mesh%nTri_loc

      ! Calculate the integral from b to z of (A_flow * (h - zeta)^n) dzeta
      z = Hs_b( ti) - mesh%zeta * Hi_b( ti)
      int_A_hminzetan = integrate_from_zeta_is_one_to_zeta_is_zetap( z, A_flow_3D_b( ti,:) * (Hs_b( ti) - z)**C%n_flow)

      ! Calculate the diffusivity term
      abs_grad_Hs = SQRT( dHs_dx_b( ti)**2 + dHs_dy_b( ti)**2)
      SIA%D_3D_b( ti,:) = -2._dp * (ice_density * grav)**C%n_flow * abs_grad_Hs**(C%n_flow - 1._dp) * int_A_hminzetan

      ! Safety
      SIA%D_3D_b( ti,:) = MAX( -C%SIA_maximum_diffusivity, SIA%D_3D_b( ti,:))

      ! Calculate the velocities
      SIA%u_3D_b( ti,:) = SIA%D_3D_b( ti,:) * dHs_dx_b( ti)
      SIA%v_3D_b( ti,:) = SIA%D_3D_b( ti,:) * dHs_dy_b( ti)

    END DO ! DO ti = 1, mesh%nTri_loc

    ! Calculate vertical shear strain rates (needed later to calculate strain heating in thermodynamics)
    DO vi = 1, mesh%nV_loc

      abs_grad_Hs = SQRT( dHs_dx_a( vi)**2 + dHs_dy_a( vi)**2)
      z = ice%Hs( vi) - mesh%zeta * ice%Hi( vi)

      DO k = 1, mesh%nz
        SIA%du_dz_3D_a( vi,k) = -2._dp * (ice_density * grav)**C%n_flow * abs_grad_Hs**(C%n_flow - 1._dp) * &
          ice%A_flow_3D( vi,k) * (ice%Hs( vi) - z( k))**C%n_flow * dHs_dx_a( vi)
        SIA%dv_dz_3D_a( vi,k) = -2._dp * (ice_density * grav)**C%n_flow * abs_grad_Hs**(C%n_flow - 1._dp) * &
          ice%A_flow_3D( vi,k) * (ice%Hs( vi) - z( k))**C%n_flow * dHs_dy_a( vi)
      END DO

    END DO ! DO vi = 1, mesh%nV_loc

    ! Clean up after yourself
    DEALLOCATE( Hi_b       )
    DEALLOCATE( Hs_b       )
    DEALLOCATE( dHs_dx_a   )
    DEALLOCATE( dHs_dy_a   )
    DEALLOCATE( dHs_dx_b   )
    DEALLOCATE( dHs_dy_b   )
    DEALLOCATE( A_flow_3D_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SIA

  SUBROUTINE remap_SIA_solver( mesh_old, mesh_new, ice, SIA)
    ! Remap the SIA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ice_velocity_solver_SIA),  INTENT(INOUT) :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_SIA_solver'
    REAL(dp)                                           :: dp_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dp_dummy = mesh_old%V( 1,1)
    dp_dummy = mesh_new%V( 1,1)
    dp_dummy = ice%Hi( 1)

    ! Reallocate shared memory

    ! Solution
    CALL reallocate_clean( SIA%u_3D_b    , mesh_new%nTri_loc, mesh_new%nz)
    CALL reallocate_clean( SIA%v_3D_b    , mesh_new%nTri_loc, mesh_new%nz)
    CALL reallocate_clean( SIA%du_dz_3D_a, mesh_new%nV_loc  , mesh_new%nz)
    CALL reallocate_clean( SIA%dv_dz_3D_a, mesh_new%nV_loc  , mesh_new%nz)

    ! Intermediate data fields
    CALL reallocate_clean( SIA%D_3D_b    , mesh_new%nTri_loc, mesh_new%nz)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SIA_solver

END MODULE ice_velocity_SIA
