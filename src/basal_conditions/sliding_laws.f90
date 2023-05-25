MODULE sliding_laws

  ! Contains all the different sliding laws.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE parameters
  USE mesh_operators                                         , ONLY: map_b_a_2D

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_basal_friction_coefficient( mesh, ice, u_b, v_b)
    ! Calculate the basal friction coefficient beta_b using the specified sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_friction_coefficient'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: u_a, v_a, uabs_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    ALLOCATE( u_a(    mesh%nV_loc))
    ALLOCATE( v_a(    mesh%nV_loc))
    ALLOCATE( uabs_a( mesh%nV_loc))

    ! Map velocities to the a-grid
    CALL map_b_a_2D( mesh, u_b, u_a)
    CALL map_b_a_2D( mesh, v_b, v_a)
    uabs_a = SQRT( u_a**2 + v_a**2)

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed (choice of beta is trivial)
      ice%beta_b = 0._dp
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
!      CALL calc_sliding_law_idealised(           mesh, ice, uabs_a)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
!      CALL calc_sliding_law_Weertman(            mesh, ice, uabs_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb-type sliding law
!      CALL calc_sliding_law_Coulomb(             mesh, ice, uabs_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
!      CALL calc_sliding_law_Coulomb_regularised( mesh, ice, uabs_a)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
!      CALL calc_sliding_law_Tsai2015(            mesh, ice, uabs_a)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
!      CALL calc_sliding_law_Schoof2005(          mesh, ice, uabs_a)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
!      CALL calc_sliding_law_ZoetIverson(         mesh, ice, uabs_a)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Clean up after yourself
    DEALLOCATE( u_a   )
    DEALLOCATE( v_a   )
    DEALLOCATE( uabs_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_friction_coefficient

!  SUBROUTINE calc_sliding_law_Weertman( mesh, ice, u_a, v_a)
!    ! Weertman-type ("power law") sliding law
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Weertman'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: uabs
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Calculate beta
!    DO vi = mesh%vi1, mesh%vi2
!
!      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
!      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)
!
!      ! Asay-Davis et al. (2016), Eq. 6
!      ice%beta_b_a( vi) = ice%beta_sq_a( vi) * uabs ** (1._dp / C%slid_Weertman_m - 1._dp)
!
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_Weertman
!
!  SUBROUTINE calc_sliding_law_Coulomb( mesh, ice, u_a, v_a)
!    ! Coulomb-type sliding law
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Coulomb'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: uabs
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Calculate the till yield stress from the till friction angle and the effective pressure
!    DO vi = mesh%vi1, mesh%vi2
!      ice%tauc_a( vi) = TAN((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
!    END DO
!    CALL sync
!
!    ! Calculate beta
!    DO vi = mesh%vi1, mesh%vi2
!
!      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
!      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)
!
!      ice%beta_b_a( vi) = ice%tauc_a( vi) / uabs
!
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_Coulomb
!
!  SUBROUTINE calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a)
!    ! Regularised Coulomb-type sliding law
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Coulomb_regularised'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: uabs
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Calculate the till yield stress from the till friction angle and the effective pressure
!    DO vi = mesh%vi1, mesh%vi2
!      ice%tauc_a( vi) = TAN((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
!    END DO
!    CALL sync
!
!    ! Calculate beta
!    DO vi = mesh%vi1, mesh%vi2
!
!      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
!      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)
!
!      ice%beta_b_a( vi) = ice%tauc_a( vi) * uabs ** (C%slid_Coulomb_reg_q_plastic - 1._dp) / (C%slid_Coulomb_reg_u_threshold ** C%slid_Coulomb_reg_q_plastic)
!
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_Coulomb_regularised
!
!  SUBROUTINE calc_sliding_law_Tsai2015(  mesh, ice, u_a, v_a)
!    ! Modified power-law relation according to Tsai et al. (2015)
!    ! (implementation based on equations provided by Asay-Davis et al., 2016)
!    !
!    ! Asay-Dvis et al.: Experimental design for three interrelated marine ice sheet and ocean model
!    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
!    ! Geoscientific Model Development 9, 2471-2497, 2016
!    !
!    ! Tsai et al.: Marine ice-sheet profiles and stability under Coulomb basal conditions,
!    ! Journal of Glaciology 61, 205–215, doi:10.3189/2015JoG14J221, 2015.
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Tsai2015'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: uabs
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Calculate beta
!    DO vi = mesh%vi1, mesh%vi2
!
!      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
!      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)
!
!      ! Asay-Dvis et al. (2016), Eq. 7
!      ice%beta_b_a( vi) = MIN( ice%alpha_sq_a( vi) * ice%Neff_a( vi), ice%beta_sq_a( vi) * uabs ** (1._dp / C%slid_Weertman_m)) * uabs**(-1._dp)
!
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_Tsai2015
!
!  SUBROUTINE calc_sliding_law_Schoof2005(  mesh, ice, u_a, v_a)
!    ! Modified power-law relation according to Tsai et al. (2015)
!    ! (implementation based on equations provided by Asay-Davis et al., 2016)
!    !
!    ! Asay-Dvis et al.: Experimental design for three interrelated marine ice sheet and ocean model
!    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
!    ! Geoscientific Model Development 9, 2471-2497, 2016
!    !
!    ! Schoof: The effect of cvitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Schoof2005'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: uabs
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Calculate beta
!    DO vi = mesh%vi1, mesh%vi2
!
!      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
!      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)
!
!      ! Asay-Dvis et al. (2016), Eq. 11
!      ice%beta_b_a( vi) = ((ice%beta_sq_a( vi) * uabs**(1._dp / C%slid_Weertman_m) * ice%alpha_sq_a( vi) * ice%Neff_a( vi)) / &
!        ((ice%beta_sq_a( vi)**C%slid_Weertman_m * uabs + (ice%alpha_sq_a( vi) * ice%Neff_a( vi))**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs**(-1._dp)
!
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_Schoof2005
!
!  SUBROUTINE calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a)
!    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_ZoetIverson'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: uabs
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Calculate the till yield stress from the till friction angle and the effective pressure
!    DO vi = mesh%vi1, mesh%vi2
!      ice%tauc_a( vi) = TAN((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
!    END DO
!    CALL sync
!
!    ! Calculate beta
!    DO vi = mesh%vi1, mesh%vi2
!
!      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
!      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)
!
!      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
!      ice%beta_b_a( vi) = ice%tauc_a( vi) * (uabs**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))
!
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_ZoetIverson
!
!  SUBROUTINE calc_sliding_law_idealised(  mesh, ice, u_a, v_a)
!    ! Sliding laws for some idealised experiments
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
!    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised'
!    REAL(dp) :: dummy_dp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! To prevent compiler warnings...
!    dummy_dp = u_a( 1)
!    dummy_dp = v_a( 1)
!
!    IF     (C%choice_idealised_sliding_law == 'ISMIP_HOM_C') THEN
!      ! ISMIP-HOM experiment C
!
!      CALL calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)
!
!    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_D') THEN
!      ! ISMIP-HOM experiment D
!
!      CALL calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)
!
!    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_E') THEN
!      ! ISMIP-HOM experiment E
!
!      CALL crash('the Glacier Arolla experiment is not implemented in UFEMISM!')
!
!    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_F') THEN
!      ! ISMIP-HOM experiment F
!
!      CALL calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)
!
!    ELSE
!      CALL crash('unknown choice_idealised_sliding_law "' // TRIM( C%choice_idealised_sliding_law) // '"!')
!    END IF
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_idealised
!
!  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)
!    ! Sliding laws for some idealised experiments
!    !
!    ! ISMIP-HOM experiment C
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_C'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: x,y
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    DO vi = mesh%vi1, mesh%vi2
!      x = mesh%V( vi,1)
!      y = mesh%V( vi,2)
!      ice%beta_b_a( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%ISMIP_HOM_L) * SIN( 2._dp * pi * y / C%ISMIP_HOM_L)
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C
!
!  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)
!    ! Sliding laws for some idealised experiments
!    !
!    ! ISMIP-HOM experiment D
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_D'
!    INTEGER                                            :: vi
!    REAL(dp)                                           :: x
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    DO vi = mesh%vi1, mesh%vi2
!      x = mesh%V( vi,1)
!      ice%beta_b_a( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%ISMIP_HOM_L)
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D
!
!  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)
!    ! Sliding laws for some idealised experiments
!    !
!    ! ISMIP-HOM experiment F
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    TYPE(type_mesh),                     INTENT(IN)    :: mesh
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_F'
!    INTEGER                                            :: vi
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    DO vi = mesh%vi1, mesh%vi2
!      ice%beta_b_a( vi) = (ice%A_flow_vav_a( vi) * 1000._dp)**(-1._dp)
!    END DO
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F


END MODULE sliding_laws
