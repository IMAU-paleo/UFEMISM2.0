MODULE laddie_physics

  ! Utilities for the laddie model

! ===== Preamble =====
! ====================
    
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE ocean_utilities                                        , ONLY: interpolate_ocean_depth

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_melt_rate( mesh, ice, ocean, laddie, Hstar)
    ! Compute melt rate using the three equations

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_melt_rate'
    INTEGER                                               :: vi
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get T and S at layer base
    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         laddie%melt( vi) = C%uniform_laddie_gamma_T * (laddie%T( vi) - ocean%T_freezing_point( vi))
       ELSE
         laddie%melt( vi) = 0.0_dp
       END IF
  
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_melt_rate

  SUBROUTINE compute_entrainment( mesh, ice, ocean, laddie, Hstar)
    ! Compute entrainment rate

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_entrainment'
    INTEGER                                               :: vi
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get T and S at layer base
    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         laddie%entr( vi) = 1.0E-6_dp
       ELSE
         laddie%entr( vi) = 0.0_dp
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_entrainment

  SUBROUTINE compute_freezing_temperature( mesh, ice, laddie)
    ! Compute freezing temperature at ice shelf base, based on Laddie salinity.
    ! TODO can maybe be merged with ice computation

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_freezing_temperature'
    INTEGER                                               :: vi
 
    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         laddie%T_freeze( vi) = freezing_lambda_1*laddie%S( vi) + freezing_lambda_2 - freezing_lambda_3*ice%Hib( vi)
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_freezing_temperature

  SUBROUTINE compute_buoyancy( mesh, ice, laddie, Hstar)
    ! Compute buoyancy = (rho_amb - rho)/rho_sw
    ! TODO update with Roquet EOS 

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_buoyancy'
    INTEGER                                               :: vi
 
    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         ! Get buoyancy
         laddie%drho_amb( vi) = C%uniform_laddie_eos_linear_beta  * (laddie%S_amb( vi)-laddie%S( vi)) &
                              - C%uniform_laddie_eos_linear_alpha * (laddie%T_amb( vi)-laddie%T( vi))
         ! Get depth-integrated buoyancy based on Hstar
         laddie%Hdrho_amb( vi) = Hstar( vi) * laddie%drho_amb( vi)
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_buoyancy

END MODULE laddie_physics

