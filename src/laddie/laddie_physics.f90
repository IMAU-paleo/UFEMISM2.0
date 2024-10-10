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
    REAL(dp)                                              :: That
    REAL(dp)                                              :: Chat
    REAL(dp)                                              :: Ctil
    REAL(dp)                                              :: Bval
    REAL(dp)                                              :: Cval
    REAL(dp)                                              :: Dval
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get friction velocity
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         laddie%u_star( vi) = C%laddie_drag_coefficient_mom * (laddie%U_a( vi)**2 + laddie%V_a( vi)**2 + C%uniform_laddie_tidal_velocity**2 )**.5
       END IF
    END DO

    ! Get gamma values. TODO add non-fixed, non-uniform option. If fixed,uniform, compute during initialisation and skip here
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         laddie%gamma_T( vi) = C%uniform_laddie_gamma_T
         laddie%gamma_S( vi) = C%uniform_laddie_gamma_T/35.0
       END IF
    END DO

    ! == Get melt rate ==
    ! ===================
    Chat = cp_ocean/L_fusion
    Ctil = cp_ice/cp_ocean

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         ! Solve three equations
         That = freezing_lambda_2 - freezing_lambda_3*ice%Hib( vi)
         ! Chat = cp_ocean / (L_fusion - cp_ice * ice%Ti( vi, 1)) TODO expand with proper Ti

         Bval = Chat*laddie%gamma_T( vi)*(That - laddie%T( vi)) + laddie%gamma_S( vi)*(1 + Chat*Ctil*(That + freezing_lambda_1*laddie%S( vi)))
         Cval = Chat*laddie%gamma_T( vi)*laddie%gamma_S( vi) * (That-laddie%T( vi) + freezing_lambda_1*laddie%S( vi))

         ! Get melt rate
         IF (4*Cval > Bval**2) THEN
           !Something wrong, set melt rate to zero. TODO check whether model needs to be crashed
           laddie%melt( vi) = 0.0
         ELSE
           laddie%melt( vi) = .5 * (-Bval + (Bval**2 - 4*Cval)**.5) 
         END IF

         ! Get temperature at ice base
         Dval = Chat*laddie%gamma_T( vi)+Chat*Ctil*laddie%melt( vi)
         IF (Dval == 0) THEN
           ! Seems like a very unlikely case, but better to be careful
           laddie%T_base( vi) = laddie%T_freeze( vi)
         ELSE
           laddie%T_base( vi) = Chat*laddie%gamma_T( vi)*laddie%T( vi)/Dval
         END IF

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
    REAL(dp), PARAMETER                                   :: mindrho = 0.000001_dp
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! TODO only Gaspar for now

    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         ! Get salinity at ice base
         laddie%S_base( vi) = (laddie%T_base( vi) - freezing_lambda_2 + freezing_lambda_3 * ice%Hib( vi)) / freezing_lambda_1
         
         ! Get buoyancy at ice base
         laddie%drho_base( vi) = C%uniform_laddie_eos_linear_beta  * (laddie%S( vi)-laddie%S_base( vi)) &
                               - C%uniform_laddie_eos_linear_alpha * (laddie%T( vi)-laddie%T_base( vi))

         ! Make sure buoyancy is positive
         laddie%drho_base( vi) = MAX(laddie%drho_base( vi),mindrho)

         ! Get entrainment
         laddie%entr( vi) = 2*C%laddie_Gaspar1988_mu/grav & 
                          * laddie%u_star( vi)**3 / (Hstar( vi) * laddie%drho_amb( vi)) &
                          - laddie%drho_base( vi) / laddie%drho_amb( vi) * laddie%melt( vi)

         ! Get detrainment = negative component of entrainment
         laddie%detr( vi) = - MIN(laddie%entr( vi),0.0_dp)
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
       IF (laddie%mask_a( vi)) THEN
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
    REAL(dp), PARAMETER                                   :: mindrho = 0.000001_dp
 
    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         ! Get buoyancy
         laddie%drho_amb( vi) = C%uniform_laddie_eos_linear_beta  * (laddie%S_amb( vi)-laddie%S( vi)) &
                              - C%uniform_laddie_eos_linear_alpha * (laddie%T_amb( vi)-laddie%T( vi))

         ! Make sure buoyancy is positive TODO expand with convection scheme
         laddie%drho_amb( vi) = MAX(laddie%drho_amb( vi),mindrho)

         ! Get depth-integrated buoyancy based on Hstar
         laddie%Hdrho_amb( vi) = Hstar( vi) * laddie%drho_amb( vi)
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_buoyancy

END MODULE laddie_physics

