MODULE sliding_laws

  ! Contains all the different sliding laws.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE parameters
  use mesh_disc_apply_operators, only: map_b_a_2D
  USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
  use mpi_distributed_memory, only: gather_to_all

  IMPLICIT NONE

CONTAINS

! == Main routine
! ===============

  SUBROUTINE calc_basal_friction_coefficient( mesh, ice, u_b, v_b)
    ! Calculate the effective basal friction coefficient using the specified sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_friction_coefficient'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: u_a, v_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    ALLOCATE( u_a( mesh%vi1:mesh%vi2))
    ALLOCATE( v_a( mesh%vi1:mesh%vi2))

    ! Map velocities to the a-grid
    CALL map_b_a_2D( mesh, u_b, u_a)
    CALL map_b_a_2D( mesh, v_b, v_a)

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed (choice of beta is trivial)
      ice%basal_friction_coefficient = 0._dp
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
      CALL calc_sliding_law_idealised(   mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL calc_sliding_law_Weertman(    mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb(     mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Budd') THEN
      ! Regularised Coulomb-type sliding law
      CALL calc_sliding_law_Budd(        mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL calc_sliding_law_Tsai2015(    mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL calc_sliding_law_Schoof2005(  mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Limit basal friction coefficient to improve stability
    ice%basal_friction_coefficient = MIN( C%slid_beta_max, ice%basal_friction_coefficient)

    ! Clean up after yourself
    DEALLOCATE( u_a   )
    DEALLOCATE( v_a   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_friction_coefficient

! == Different sliding laws
! =========================

  SUBROUTINE calc_sliding_law_Weertman( mesh, ice, u_a, v_a)
    ! Weertman-type ("power law") sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_sliding_law_Weertman'
    INTEGER                                               :: vi
    REAL(dp)                                              :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Assume that bed roughness is represented by the slid_beta_sq field
    ice%bed_roughness = ice%slid_beta_sq

    ! Scale bed roughness based on grounded area fractions
    CALL apply_grounded_fractions_to_bed_roughness( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Davis et al. (2016), Eq. 6; replacing slid_beta_sq by the bed roughness field from above
      ice%basal_friction_coefficient( vi) = ice%bed_roughness( vi) * uabs ** (1._dp / C%slid_Weertman_m - 1._dp)

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Weertman

  SUBROUTINE calc_sliding_law_Coulomb( mesh, ice, u_a, v_a)
    ! Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_sliding_law_Coulomb'
    INTEGER                                               :: vi
    REAL(dp)                                              :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Compute bed roughness based on till friction angle
    ice%bed_roughness = TAN((pi / 180._dp) * ice%till_friction_angle)

    ! Scale bed roughness based on grounded area fractions
    CALL apply_grounded_fractions_to_bed_roughness( mesh, ice)

    ! == Till yield stress
    ! ====================

    ! Calculate the till yield stress from the effective pressure and bed roughness
    ice%till_yield_stress = ice%effective_pressure * ice%bed_roughness

    ! == Extend till yield stress over ice-free land neighbours
    ! =========================================================

    CALL extend_till_yield_stress_to_neighbours( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ice%basal_friction_coefficient( vi) = ice%till_yield_stress( vi) / uabs

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb

  SUBROUTINE calc_sliding_law_Budd( mesh, ice, u_a, v_a)
    ! Budd-type sliding law (formerly known as regularised Coulomb-type sliding law)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_sliding_law_Budd'
    INTEGER                                               :: vi
    REAL(dp)                                              :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Compute bed roughness based on till friction angle
    ice%bed_roughness = TAN((pi / 180._dp) * ice%till_friction_angle)

    ! Scale bed roughness based on grounded area fractions
    CALL apply_grounded_fractions_to_bed_roughness( mesh, ice)

    ! == Till yield stress
    ! ====================

    ! Calculate the till yield stress from the effective pressure and bed roughness
    ice%till_yield_stress = ice%effective_pressure * ice%bed_roughness

    ! == Extend till yield stress over ice-free land neighbours
    ! =========================================================

    CALL extend_till_yield_stress_to_neighbours( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ice%basal_friction_coefficient( vi) = ice%till_yield_stress( vi) * uabs ** (C%slid_Budd_q_plastic - 1._dp) / (C%slid_Budd_u_threshold ** C%slid_Budd_q_plastic)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Budd

  SUBROUTINE calc_sliding_law_Tsai2015(  mesh, ice, u_a, v_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    !
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Tsai et al.: Marine ice-sheet profiles and stability under Coulomb basal conditions,
    ! Journal of Glaciology 61, 205–215, doi:10.3189/2015JoG14J221, 2015.

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_sliding_law_Tsai2015'
    INTEGER                                               :: vi
    REAL(dp)                                              :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Assume that bed roughness is represented by the slid_beta_sq field
    ice%bed_roughness = ice%slid_beta_sq

    ! Scale bed roughness based on grounded area fractions
    CALL apply_grounded_fractions_to_bed_roughness( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Davis et al. (2016), Eq. 7; replacing slid_beta_sq by the bed roughness field from above
      ice%basal_friction_coefficient( vi) = MIN( ice%slid_alpha_sq( vi) * ice%effective_pressure( vi), ice%bed_roughness( vi) * uabs ** (1._dp / C%slid_Weertman_m)) * uabs**(-1._dp)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Tsai2015

  SUBROUTINE calc_sliding_law_Schoof2005(  mesh, ice, u_a, v_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    !
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Schoof: The effect of cvitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_sliding_law_Schoof2005'
    INTEGER                                               :: vi
    REAL(dp)                                              :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Assume that bed roughness is represented by the slid_beta_sq field
    ice%bed_roughness = ice%slid_beta_sq

    ! Scale bed roughness based on grounded area fractions
    CALL apply_grounded_fractions_to_bed_roughness( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Davis et al. (2016), Eq. 11; replacing slid_beta_sq by the bed roughness field from above
      ice%basal_friction_coefficient( vi) = ((ice%bed_roughness( vi) * uabs**(1._dp / C%slid_Weertman_m) * ice%slid_alpha_sq( vi) * ice%effective_pressure( vi)) / &
        ((ice%bed_roughness( vi)**C%slid_Weertman_m * uabs + (ice%slid_alpha_sq( vi) * ice%effective_pressure( vi))**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs**(-1._dp)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Schoof2005

  SUBROUTINE calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_sliding_law_ZoetIverson'
    INTEGER                                               :: vi, ci, vc
    REAL(dp)                                              :: uabs, min_neighbour
    LOGICAL,  DIMENSION(mesh%nV)                          :: mask_grounded_ice_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: till_yield_stress_tot
    LOGICAL                                               :: found_grounded_neighbour
    REAL(dp)                                              :: weight_gr, exponent_gr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Bed roughness
    ! ================

    ! Compute bed roughness based on till friction angle
    ice%bed_roughness = TAN((pi / 180._dp) * ice%till_friction_angle)

    ! Scale bed roughness based on grounded area fractions
    CALL apply_grounded_fractions_to_bed_roughness( mesh, ice)

    ! == Till yield stress
    ! ====================

    ! Calculate the till yield stress from the effective pressure and bed roughness
    ice%till_yield_stress = ice%effective_pressure * ice%bed_roughness

    ! == Extend till yield stress over ice-free land neighbours
    ! =========================================================

    CALL extend_till_yield_stress_to_neighbours( mesh, ice)

    ! == Basal friction field
    ! =======================

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
      ice%basal_friction_coefficient( vi) = ice%till_yield_stress( vi) * (uabs**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_ZoetIverson

  SUBROUTINE calc_sliding_law_idealised(  mesh, ice, u_a, v_a)
    ! Sliding laws for some idealised experiments

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_sliding_law_idealised'
    REAL(dp)                                              :: dummy1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dummy1 = u_a( mesh%vi1)
    dummy1 = v_a( mesh%vi1)

    ! Use the specified idealised sliding law
    IF     (C%choice_idealised_sliding_law == 'ISMIP-HOM_C') THEN
      ! ISMIP-HOM experiment C
      CALL calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP-HOM_D') THEN
      ! ISMIP-HOM experiment D
      CALL calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP-HOM_E') THEN
      ! ISMIP-HOM experiment E
      CALL crash('the Glacier Arolla experiment is not implemented in UFEMISM!')
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP-HOM_F') THEN
      ! ISMIP-HOM experiment F
      CALL calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)
    ELSE
      CALL crash('unknown choice_idealised_sliding_law "' // TRIM( C%choice_idealised_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised

  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment C

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_C'
    INTEGER                                            :: vi
    REAL(dp)                                           :: x,y

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      ice%basal_friction_coefficient( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%refgeo_idealised_ISMIP_HOM_L) * SIN( 2._dp * pi * y / C%refgeo_idealised_ISMIP_HOM_L)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C

  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment D

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_D'
    INTEGER                                            :: vi
    REAL(dp)                                           :: x

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      ice%basal_friction_coefficient( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%refgeo_idealised_ISMIP_HOM_L)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D

  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment F

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_F'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      ice%basal_friction_coefficient( vi) = (C%uniform_Glens_flow_factor * 1000._dp)**(-1._dp)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F

! == Utilities
! ============

  SUBROUTINE apply_grounded_fractions_to_bed_roughness( mesh, ice)
    ! Scale bed roughness based on grounded area fractions

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_grounded_fractions_to_bed_roughness'
    INTEGER                                               :: vi
    REAL(dp)                                              :: weight_gr, exponent_hi, exponent_hs, exponent_gr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this is actually wanted
    IF (.NOT. C%do_subgrid_friction_on_A_grid) THEN
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! And exit
      RETURN
    END IF

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise grounded area fraction weight
      weight_gr = 1._dp

      ! Compute exponent for this vertex's weight based on ice thickness
      exponent_hi = LOG10( MAX( 1._dp, ice%Hi( vi)))
      ! Compute exponent for this vertex's weight based on surface gradients
      exponent_hs = ice%Hs_slope( vi) / 0.005_dp
      ! Compute final exponent for this vertex's weight
      exponent_gr = MAX( 0._dp, exponent_hi - exponent_hs)

      ! Compute a weight based on the grounded area fractions
      IF (ice%mask_gl_gr( vi)) THEN
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      ELSEIF (ice%mask_cf_gr( vi)) THEN
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      ELSEIF (ice%mask_gl_fl( vi)) THEN
        weight_gr = ice%fraction_gr( vi)**exponent_gr

      ELSEIF (ice%mask_grounded_ice( vi)) THEN
        weight_gr = 1._dp

      ELSEIF (ice%mask_floating_ice( vi)) THEN
        weight_gr = 0._dp

      ELSEIF (ice%mask_icefree_ocean( vi)) THEN
        weight_gr = 0._dp

      END IF

      ! Just in case
      weight_gr = MIN( 1._dp, MAX( 0._dp, weight_gr))

      ! Compute till friction angle accounting for grounded area fractions
      ice%bed_roughness( vi) = ice%bed_roughness(vi) * weight_gr

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_grounded_fractions_to_bed_roughness

  SUBROUTINE extend_till_yield_stress_to_neighbours( mesh, ice)
    ! Extend till yield stress over ice-free land neighbours

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'extend_till_yield_stress_to_neighbours'
    INTEGER                                               :: vi, ci, vc
    REAL(dp)                                              :: min_neighbour
    LOGICAL,  DIMENSION(mesh%nV)                          :: mask_grounded_ice_tot
    REAL(dp), DIMENSION(mesh%nV)                          :: till_yield_stress_tot
    LOGICAL                                               :: found_grounded_neighbour

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather data from all processes
    CALL gather_to_all( ice%mask_grounded_ice, mask_grounded_ice_tot)
    CALL gather_to_all(      ice%till_yield_stress, till_yield_stress_tot)

    DO vi = mesh%vi1, mesh%vi2

      ! Skip if not ice-free land
      IF (.NOT. ice%mask_icefree_land( vi)) CYCLE

      ! Initialise
      found_grounded_neighbour = .FALSE.

      ! Initialise with default values and assume a till yield
      ! stress equivalent to a column of ice of 1000 metres
      ! with no hydrology and maximum bed roughness.
      min_neighbour = 1000._dp * ice_density * grav

      ! Check grounded neighbours
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi, ci)
        IF (mask_grounded_ice_tot( vc)) THEN
          min_neighbour = MIN( min_neighbour, till_yield_stress_tot( vc))
          found_grounded_neighbour = .TRUE.
        END IF
      END DO

      IF (found_grounded_neighbour) THEN
        ! Use the minimum value among neighbours
        ice%till_yield_stress( vi) = min_neighbour
      ELSE
        ! Use a default minimum value to avoid 0 friction
        ice%till_yield_stress( vi) = C%Hi_min * ice_density * grav
      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extend_till_yield_stress_to_neighbours

END MODULE sliding_laws
