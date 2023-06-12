MODULE ice_velocity_main

  ! Contains all the routines needed to calculate instantaneous ice velocities for the current modelled ice-sheet geometry.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE petsc_basic                                            , ONLY: solve_matrix_equation_CSR_PETSc
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model, type_ice_velocity_solver_SIA, type_ice_velocity_solver_SSA, &
                                                                     type_ice_velocity_solver_DIVA, type_ice_velocity_solver_BPA
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate_clean
  USE mesh_operators                                         , ONLY: map_b_a_2D, map_b_a_3D, ddx_a_a_2D, ddy_a_a_2D
  USE ice_velocity_SIA                                       , ONLY: initialise_SIA_solver , solve_SIA , remap_SIA_solver
  USE ice_velocity_SSA                                       , ONLY: initialise_SSA_solver , solve_SSA , remap_SSA_solver
  USE ice_velocity_DIVA                                      , ONLY: initialise_DIVA_solver, solve_DIVA, remap_DIVA_solver
  USE ice_velocity_BPA                                       , ONLY: initialise_BPA_solver , solve_BPA , remap_BPA_solver
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_dp_2D
  USE mesh_zeta                                              , ONLY: vertical_average

  IMPLICIT NONE

CONTAINS

! == The main routines, to be called from the ice dynamics module

  SUBROUTINE initialise_velocity_solver( mesh, ice, region_name)
    ! Initialise the velocity solver for the chosen stress balance approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_velocity_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) then
      WRITE(*,"(A)") '   Initialising ice velocities using ' // &
                      colour_string( TRIM( C%choice_stress_balance_approximation),'light blue') // &
                     ' dynamics...'
    END IF

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      CALL initialise_SIA_solver(  mesh, ice%SIA)
    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      CALL initialise_SSA_solver(  mesh, ice%SSA, region_name)
    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      CALL initialise_SIA_solver(  mesh, ice%SIA)
      CALL initialise_SSA_solver(  mesh, ice%SSA, region_name)
    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      CALL initialise_DIVA_solver( mesh, ice%DIVA, region_name)
    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      CALL initialise_BPA_solver(  mesh, ice%BPA, region_name)
    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_velocity_solver

  SUBROUTINE solve_stress_balance( mesh, ice)
    ! Calculate all ice velocities based on the chosen stress balance approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_stress_balance'
    INTEGER                                            :: vi,ti
    REAL(dp), DIMENSION(mesh%nz)                       :: u_prof, v_prof

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      ! Calculate velocities according to the Shallow Ice Approximation

      CALL solve_SIA( mesh, ice, ice%SIA)
      CALL set_ice_velocities_to_SIA_results( mesh, ice, ice%SIA)

    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      ! Calculate velocities according to the Shallow Shelf Approximation

      CALL solve_SSA( mesh, ice, ice%SSA)
      CALL set_ice_velocities_to_SSA_results( mesh, ice, ice%SSA)

    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      ! Calculate velocities according to the hybrid SIA/SSA

      CALL solve_SIA( mesh, ice, ice%SIA)
      CALL solve_SSA( mesh, ice, ice%SSA)
      CALL set_ice_velocities_to_SIASSA_results( mesh, ice, ice%SIA, ice%SSA)

    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      ! Calculate velocities according to the Depth-Integrated Viscosity Approximation

      CALL solve_DIVA( mesh, ice, ice%DIVA)
      CALL set_ice_velocities_to_DIVA_results( mesh, ice, ice%DIVA)

    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      ! Calculate velocities according to the Depth-Integrated Viscosity Approximation

      CALL solve_BPA( mesh, ice, ice%BPA)
      CALL set_ice_velocities_to_BPA_results( mesh, ice, ice%BPA)

    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

  ! == Fill in derived velocity fields (surface, base, vertical average)

    DO ti = mesh%ti1, mesh%ti2

      ! Surface
      ice%u_surf_b(    ti) = ice%u_3D_b( ti,1)
      ice%v_surf_b(    ti) = ice%v_3D_b( ti,1)
      ice%uabs_surf_b( ti) = SQRT( ice%u_surf_b( ti)**2 + ice%v_surf_b( ti)**2)

      ! Base
      ice%u_base_b(    ti) = ice%u_3D_b( ti,C%nz)
      ice%v_base_b(    ti) = ice%v_3D_b( ti,C%nz)
      ice%uabs_base_b( ti) = SQRT( ice%u_base_b( ti)**2 + ice%v_base_b( ti)**2)

      ! Vertical average
      u_prof = ice%u_3D_b( ti,:)
      v_prof = ice%v_3D_b( ti,:)
      ice%u_vav_b( ti) = vertical_average( mesh%zeta, u_prof)
      ice%v_vav_b( ti) = vertical_average( mesh%zeta, v_prof)
      ice%uabs_vav_b( ti) = SQRT( ice%u_vav_b( ti)**2 + ice%v_vav_b( ti)**2)

    END DO

  ! == Calculate velocities on the a-grid (needed to calculate the vertical velocity w, and for writing to output)

    ! 3-D
    CALL map_b_a_3D( mesh, ice%u_3D_b  , ice%u_3D  )
    CALL map_b_a_3D( mesh, ice%v_3D_b  , ice%v_3D  )

    ! Surface
    CALL map_b_a_2D( mesh, ice%u_surf_b, ice%u_surf)
    CALL map_b_a_2D( mesh, ice%v_surf_b, ice%v_surf)

    ! Base
    CALL map_b_a_2D( mesh, ice%u_base_b, ice%u_base)
    CALL map_b_a_2D( mesh, ice%v_base_b, ice%v_base)

    ! Vertical average
    CALL map_b_a_2D( mesh, ice%u_vav_b , ice%u_vav )
    CALL map_b_a_2D( mesh, ice%v_vav_b , ice%v_vav )

    ! Absolute
    DO vi = mesh%vi1, mesh%vi2
      ice%uabs_surf( vi) = SQRT( ice%u_surf( vi)**2 + ice%v_surf( vi)**2)
      ice%uabs_base( vi) = SQRT( ice%u_base( vi)**2 + ice%v_base( vi)**2)
      ice%uabs_vav(  vi) = SQRT( ice%u_vav(  vi)**2 + ice%v_vav(  vi)**2)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_stress_balance

  SUBROUTINE calc_vertical_velocities( mesh, ice)!, BMB)
    ! Calculate vertical velocity w from conservation of mass
    !
    ! NOTE: since the vertical velocities for floating ice depend on
    !       the thinning rate dH/dt, this routine must be called
    !       after having calculated dHi_dt!
    !
    ! Derivation:
    !
    ! Conservation of mass, combined with the incompressibility
    ! condition (i.e. constant density) of ice, is described by:
    !
    !   du/dx + dv/dy + dw/dz = 0
    !
    ! Applying the zeta coordinate transformation yields:
    !
    !   du/dxp + dzeta/dx du/dzeta + dv/dxp + dzeta/dy dv/dzeta + dzeta/dz dw/dzeta = 0
    !
    ! The terms du/dxp + dv/dyp describe the two-dimensional divergence in scaled coordinates:
    !
    !   grad uv = du/dxp + dv/dyp
    !
    ! The average value over a single grid cell (Voronoi cell) of this divergence is:
    !
    !   grad uv = intint_Voronoi (grad uv) dA / intint dA = 1/A intint_Voronoi (grad uv) dA
    !
    ! By applying the divergence theorem, the surface integral over the Voronoi cell
    ! can be transformed into a loop integral over the boundary of that Voronoi cell:
    !
    !   grad uv = 1/A cint (uv * n_hat) dS
    !
    ! Here, n_hat is the outward unit normal to the Voronoi cell boundary. Substituting
    ! this into the equation for conservation of mass yields:
    !
    !   dw/dzeta = -1 / dzeta/dz [ 1/A cint (uv * n_hat) dS + dzeta/dx du/zeta + dzeta/dy dv/dzeta]
    !
    ! The vertical velocity w at the ice base is equal to the horizontal motion along
    ! the sloping ice base, plus the vertical motion of the ice base itself, plus the
    ! vertical motion of an ice particle with respect to the ice base (i.e. the basal melt rate):
    !
    !   w( z=b) = u( z=b) * dH_base/dx + v( z=b) * dH_base/dy + dH_base/dt + M_base
    !
    ! With this boundary condition, dw/dzeta can be integrated over zeta to yield w( z).

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    TYPE(type_BMB_model),                INTENT(IN)    :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_velocities'
    INTEGER                                            :: vi,ks,ci,vj,ei
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHib_dx
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHib_dy
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dHib_dt
    REAL(dp)                                           :: dzeta
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: u_3D_c
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: v_3D_c
    REAL(dp)                                           :: cint_un_dS, dS, u_ks, v_ks, un_dS, grad_uv_ks
    REAL(dp), DIMENSION(2)                             :: n_hat
    REAL(dp)                                           :: du_dzeta_ks, dv_dzeta_ks
    REAL(dp)                                           :: dzeta_dx_ks, dzeta_dy_ks, dzeta_dz_ks
    REAL(dp)                                           :: dw_dzeta_ks

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL warning('need to include basal melt rates here!')

    ! Allocate shared memory
    ALLOCATE( dHib_dx( mesh%nV_loc))
    ALLOCATE( dHib_dy( mesh%nV_loc))
    ALLOCATE( dHib_dt( mesh%nV_loc))
    ALLOCATE( u_3D_c(  mesh%nE_loc, mesh%nz))
    ALLOCATE( v_3D_c(  mesh%nE_loc, mesh%nz))

    DO vi = 1, mesh%nV_loc

      ! Calculate rate of change of ice base elevation
      IF     (ice%mask_sheet( vi)) THEN
        ! For grounded ice, the ice base simply moves with the bedrock
        dHib_dt( vi) =  ice%dHb_dt( vi)
      ELSEIF (ice%mask_shelf( vi)) THEN
        ! For floating ice, the ice base moves according to the thinning rate times the density fraction
        dHib_dt( vi) = -ice%dHi_dt( vi) * ice_density / seawater_density
      ELSE
        ! No ice, so no vertical velocity
        dHib_dt( vi) = 0._dp
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Calculate slopes of the ice base
    CALL ddx_a_a_2D( mesh, ice%Hib, dHib_dx)
    CALL ddy_a_a_2D( mesh, ice%Hib, dHib_dy)

    ! Calculate u,v on the c-grid (edges)
    CALL map_velocities_from_b_to_c_3D( mesh, ice%u_3D_b, ice%v_3D_b, u_3D_c, v_3D_c)

    ! Calculate vertical velocities by solving conservation of mass in each 3-D cell
    DO vi = 1, mesh%nV_loc

      ! No ice means no velocity
      IF (.NOT. ice%mask_ice( vi)) THEN
        ice%w_3D( vi,:) = 0._dp
        CYCLE
      END IF

      ! Calculate the vertical velocity at the ice base
      !
      ! NOTE: BMB is defined so that a positive number means accumulation of ice;
      !       at the ice base, that means that a positive BMB means a positive
      !       value of w

      ice%w_3D( vi,C%nz) = (ice%u_3D( vi,C%nz) * dHib_dx( vi)) + &
                           (ice%v_3D( vi,C%nz) * dHib_dy( vi)) + &
                            dHib_dt( vi)! + BMB%BMB( vi)

      ! Exception for very thin ice / ice margin: assume horizontal stretching
      ! is negligible, so that w( z) = w( z = b)
      IF (ice%Hi( vi) < 10._dp) THEN
        ice%w_3D( vi,:) = ice%w_3D( vi,C%nz)
        CYCLE
      END IF ! IF (ice%mask_margin_a( vi) == 1 .OR. ice%Hi_a( vi) < 10._dp) THEN

      ! Calculate vertical velocities by integrating dw/dz over the vertical column

      DO ks = mesh%nz-1, 1, -1

        dzeta = mesh%zeta( ks+1) - mesh%zeta( ks)

        ! Integrate u*n_hat around the Voronoi cell boundary
        cint_un_dS = 0._dp
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C(  vi,ci)
          ei = mesh%VE( vi,ci)
          ! Velocities at this section of the boundary
          u_ks = 0.5_dp * (u_3D_c( ei,ks) + u_3D_c( ei,ks+1))
          v_ks = 0.5_dp * (v_3D_c( ei,ks) + v_3D_c( ei,ks+1))
          ! Length of this section of the boundary
          dS = mesh%Cw( vi,ci)
          ! Outward normal vector to this section of the boundary
          n_hat = mesh%V( vj,:) - mesh%V( vi,:)
          n_hat = n_hat / NORM2( n_hat)
          ! Line integral over this section of the boundary
          un_dS = (u_ks * n_hat( 1) + v_ks * n_hat( 2)) * dS
          ! Add to loop integral
          cint_un_dS = cint_un_dS + un_dS
        END DO

        ! Calculate grad uv from the divergence theorem
        grad_uv_ks = cint_un_dS / mesh%A( vi)

        ! Calculate du/dzeta, dv/dzeta
        du_dzeta_ks = (ice%u_3D( vi,ks+1) - ice%u_3D( vi,ks)) / dzeta
        dv_dzeta_ks = (ice%v_3D( vi,ks+1) - ice%v_3D( vi,ks)) / dzeta

        ! Calculate dzeta/dx, dzeta/dy, dzeta/dz
        dzeta_dx_ks = 0.5_dp * (ice%dzeta_dx_ak( vi,ks) + ice%dzeta_dx_ak( vi,ks+1))
        dzeta_dy_ks = 0.5_dp * (ice%dzeta_dy_ak( vi,ks) + ice%dzeta_dy_ak( vi,ks+1))
        dzeta_dz_ks = 0.5_dp * (ice%dzeta_dz_ak( vi,ks) + ice%dzeta_dz_ak( vi,ks+1))

        ! Calculate dw/dzeta
        dw_dzeta_ks = -1._dp / dzeta_dz_ks * (grad_uv_ks + dzeta_dx_ks * du_dzeta_ks + dzeta_dy_ks * dv_dzeta_ks)

        ! Calculate w
        ice%w_3D( vi,ks) = ice%w_3D( vi,ks+1) - dzeta * dw_dzeta_ks

      END DO ! DO k = C%nz-1, 1, -1

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Clean up after yourself
    DEALLOCATE( dHib_dx)
    DEALLOCATE( dHib_dy)
    DEALLOCATE( dHib_dt)
    DEALLOCATE( u_3D_c)
    DEALLOCATE( v_3D_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_velocities

  SUBROUTINE remap_velocity_solver( mesh_old, mesh_new, ice)
    ! Remap the velocity solver for the chosen stress balance approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_velocity_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      CALL remap_SIA_solver(  mesh_old, mesh_new, ice, ice%SIA)
    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      CALL remap_SSA_solver(  mesh_old, mesh_new, ice, ice%SSA)
    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      CALL remap_SIA_solver(  mesh_old, mesh_new, ice, ice%SIA)
      CALL remap_SSA_solver(  mesh_old, mesh_new, ice, ice%SSA)
    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      CALL remap_DIVA_solver( mesh_old, mesh_new, ice, ice%DIVA)
    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      CALL crash('fixme!')
!      CALL remap_BPA_solver(  mesh_old, mesh_new, ice, ice%BPA)
    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_velocity_solver

! == Set applied ice model velocities to stress balance results

  SUBROUTINE set_ice_velocities_to_SIA_results( mesh, ice, SIA)
    ! Set applied ice model velocities and strain rates to SIA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_ice_velocity_solver_SIA),  INTENT(IN)    :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_SIA_results'
    INTEGER                                            :: vi,ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = SIA%u_3D_b( ti,:)
      ice%v_3D_b( ti,:) = SIA%v_3D_b( ti,:)
    END DO

    ! Strain rates
    DO vi = mesh%vi1, mesh%vi2
      ice%du_dz_3D( vi,:) = SIA%du_dz_3D( vi,:)
      ice%dv_dz_3D( vi,:) = SIA%dv_dz_3D( vi,:)
    END DO

    ! In the SIA, horizontal gradients of u,v, and all gradients of w, are neglected
    ice%du_dx_3D = 0._dp
    ice%du_dy_3D = 0._dp
    ice%dv_dx_3D = 0._dp
    ice%dv_dy_3D = 0._dp
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_SIA_results

  SUBROUTINE set_ice_velocities_to_SSA_results( mesh, ice, SSA)
    ! Set applied ice model velocities to SSA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)    :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_SSA_results'
    INTEGER                                            :: ti,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = SSA%u_b( ti)
      ice%v_3D_b( ti,:) = SSA%v_b( ti)
    END DO

    ! Strain rates
    DO vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D( vi,:) = SSA%du_dx_a( vi)
      ice%du_dy_3D( vi,:) = SSA%du_dy_a( vi)
      ice%dv_dx_3D( vi,:) = SSA%dv_dx_a( vi)
      ice%dv_dy_3D( vi,:) = SSA%dv_dy_a( vi)
    END DO

    ! In the SSA, vertical gradients of u,v, and all gradients of w, are neglected
    ice%du_dz_3D = 0._dp
    ice%dv_dz_3D = 0._dp
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_SSA_results

  SUBROUTINE set_ice_velocities_to_SIASSA_results( mesh, ice, SIA, SSA)
    ! Set applied ice model velocities to hybrid SIA/SSA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_ice_velocity_solver_SIA),  INTENT(IN)    :: SIA
    TYPE(type_ice_velocity_solver_SSA),  INTENT(IN)    :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_SIASSA_results'
    INTEGER                                            :: ti,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_hybrid_SIASSA_scheme == 'add') THEN
      ! u = u_SIA + u_SSA

      ! Velocities
      DO ti = mesh%ti1, mesh%ti2
        ice%u_3D_b( ti,:) = SIA%u_3D_b( ti,:) + SSA%u_b( ti)
        ice%v_3D_b( ti,:) = SIA%v_3D_b( ti,:) + SSA%v_b( ti)
      END DO

      ! Strain rates
      DO vi = mesh%vi1, mesh%vi2
        ice%du_dz_3D( vi,:) = SIA%du_dz_3D( vi,:)
        ice%dv_dz_3D( vi,:) = SIA%dv_dz_3D( vi,:)
        ice%du_dx_3D( vi,:) = SSA%du_dx_a(  vi  )
        ice%du_dy_3D( vi,:) = SSA%du_dy_a(  vi  )
        ice%dv_dx_3D( vi,:) = SSA%dv_dx_a(  vi  )
        ice%dv_dy_3D( vi,:) = SSA%dv_dy_a(  vi  )
      END DO

      ! In the hybrid SIA/SSA, gradients of w are neglected
      ice%dw_dx_3D = 0._dp
      ice%dw_dy_3D = 0._dp
      ice%dw_dz_3D = 0._dp

    ELSE
      CALL crash('unknown choice_hybrid_SIASSA_scheme_config "' // TRIM( C%choice_hybrid_SIASSA_scheme) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_SIASSA_results

  SUBROUTINE set_ice_velocities_to_DIVA_results( mesh, ice, DIVA)
    ! Set applied ice model velocities to DIVA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_ice_velocity_solver_DIVA), INTENT(IN)    :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_DIVA_results'
    INTEGER                                            :: ti,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = DIVA%u_3D_b( ti,:)
      ice%v_3D_b( ti,:) = DIVA%v_3D_b( ti,:)
    END DO

    ! Strain rates
    DO vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D( vi,:) = DIVA%du_dx_a(    vi  )
      ice%du_dy_3D( vi,:) = DIVA%du_dy_a(    vi  )
      ice%du_dz_3D( vi,:) = DIVA%du_dz_3D_a( vi,:)
      ice%dv_dx_3D( vi,:) = DIVA%dv_dx_a(    vi  )
      ice%dv_dy_3D( vi,:) = DIVA%dv_dy_a(    vi  )
      ice%dv_dz_3D( vi,:) = DIVA%dv_dz_3D_a( vi,:)
    END DO
    ! In the DIVA, gradients of w are neglected
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_DIVA_results

  SUBROUTINE set_ice_velocities_to_BPA_results( mesh, ice, BPA)
    ! Set applied ice model velocities and strain rates to BPA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_ice_velocity_solver_BPA),  INTENT(IN)    :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_BPA_results'
    INTEGER                                            :: ti,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = BPA%u_bk( ti,:)
      ice%v_3D_b( ti,:) = BPA%v_bk( ti,:)
    END DO

    ! Strain rates
    DO vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D( vi,:) = BPA%du_dx_ak( vi,:)
      ice%du_dy_3D( vi,:) = BPA%du_dy_ak( vi,:)
      ice%du_dz_3D( vi,:) = BPA%du_dz_ak( vi,:)
      ice%dv_dx_3D( vi,:) = BPA%dv_dx_ak( vi,:)
      ice%dv_dy_3D( vi,:) = BPA%dv_dy_ak( vi,:)
      ice%dv_dz_3D( vi,:) = BPA%dv_dz_ak( vi,:)
    END DO

    ! In the BPA, gradients of w are neglected
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_BPA_results

! == Calculate velocities on the c-grid for solving the ice thickness equation

  SUBROUTINE map_velocities_from_b_to_c_2D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    ! Calculate velocities on the c-grid for solving the ice thickness equation
    !
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b_partial
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b_partial
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: u_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_from_b_to_c_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: u_b_tot, v_b_tot
    INTEGER                                            :: eii, ei, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( u_b_tot( mesh%nTri))
    ALLOCATE( v_b_tot( mesh%nTri))

    ! Gather the full b-grid velocity fields to all processes
    CALL gather_to_all_dp_1D( u_b_partial, u_b_tot)
    CALL gather_to_all_dp_1D( v_b_partial, v_b_tot)

    ! Map velocities from the b-grid (triangles) to the c-grid (edges)
    DO eii = 1, mesh%nE_loc

      ei = mesh%ei1 + eii - 1
      til = mesh%ETri( ei,5)
      tir = mesh%ETri( ei,6)

      IF     (til == 0 .AND. tir > 0) THEN
        u_c( eii) = u_b_tot( tir)
        v_c( eii) = v_b_tot( tir)
      ELSEIF (tir == 0 .AND. til > 0) THEN
        u_c( eii) = u_b_tot( til)
        v_c( eii) = v_b_tot( til)
      ELSEIF (til >  0 .AND. tir > 0) THEN
        u_c( eii) = (u_b_tot( til) + u_b_tot( tir)) / 2._dp
        v_c( eii) = (v_b_tot( til) + v_b_tot( tir)) / 2._dp
      ELSE
        CALL crash('something is seriously wrong with the ETri array of this mesh!')
      END IF

    END DO

    ! Clean up after yourself
    DEALLOCATE( u_b_tot)
    DEALLOCATE( v_b_tot)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_from_b_to_c_2D

  SUBROUTINE map_velocities_from_b_to_c_3D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    ! Calculate velocities on the c-grid for solving the ice thickness equation
    !
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_b_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_b_partial
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: u_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_from_b_to_c_3D'
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: u_b_tot, v_b_tot
    INTEGER                                            :: eii, ei, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ALLOCATE( u_b_tot( mesh%nTri, mesh%nz))
    ALLOCATE( v_b_tot( mesh%nTri, mesh%nz))

    ! Gather the full b-grid velocity fields to all processes
    CALL gather_to_all_dp_2D( u_b_partial, u_b_tot)
    CALL gather_to_all_dp_2D( v_b_partial, v_b_tot)

    ! Map velocities from the b-grid (triangles) to the c-grid (edges)
    DO eii = 1, mesh%nE_loc

      ei = mesh%ei1 + eii - 1
      til = mesh%ETri( ei,5)
      tir = mesh%ETri( ei,6)

      IF     (til == 0 .AND. tir > 0) THEN
        u_c( eii,:) = u_b_tot( tir,:)
        v_c( eii,:) = v_b_tot( tir,:)
      ELSEIF (tir == 0 .AND. til > 0) THEN
        u_c( eii,:) = u_b_tot( til,:)
        v_c( eii,:) = v_b_tot( til,:)
      ELSEIF (til >  0 .AND. tir > 0) THEN
        u_c( eii,:) = (u_b_tot( til,:) + u_b_tot( tir,:)) / 2._dp
        v_c( eii,:) = (v_b_tot( til,:) + v_b_tot( tir,:)) / 2._dp
      ELSE
        CALL crash('something is seriously wrong with the ETri array of this mesh!')
      END IF

    END DO

    ! Clean up after yourself
    DEALLOCATE( u_b_tot)
    DEALLOCATE( v_b_tot)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_from_b_to_c_3D

END MODULE ice_velocity_main
