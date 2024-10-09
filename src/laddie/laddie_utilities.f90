MODULE laddie_utilities

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

  SUBROUTINE compute_ambient_TS( mesh, laddie, ocean, ice)
    ! Compute T and S of ambient ocean water at the depth of LADDIE's layer bottom
    ! through vertical interpolation

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_ambient_TS'
    INTEGER                                               :: vi
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get T and S at layer base
    DO vi = mesh%vi1, mesh%vi2
       IF (ice%mask_floating_ice( vi)) THEN
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T( vi,:), laddie%H( vi) - ice%Hib( vi), laddie%T_amb( vi))
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S( vi,:), laddie%H( vi) - ice%Hib( vi), laddie%S_amb( vi))
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_ambient_TS

  SUBROUTINE allocate_laddie_model( mesh, laddie)
    ! Allocate variables of the laddie model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'allocate_laddie_model'
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Thickness
    ALLOCATE( laddie%H                  ( mesh%vi1:mesh%vi2              )) ! [m]             Layer thickness
    ALLOCATE( laddie%H_prev             ( mesh%vi1:mesh%vi2              )) ! [m]             previous timestep
    ALLOCATE( laddie%H_next             ( mesh%vi1:mesh%vi2              )) ! [m]             next timestep   
    ALLOCATE( laddie%dH_dt              ( mesh%vi1:mesh%vi2              )) ! [m]             change

    laddie%H              = 0._dp
    laddie%H_prev         = 0._dp
    laddie%H_next         = 0._dp
    laddie%dH_dt          = 0._dp

    ! Velocities
    ALLOCATE( laddie%U                  ( mesh%ti1:mesh%ti2              )) ! [m s^-1]        2D velocity
    ALLOCATE( laddie%V                  ( mesh%ti1:mesh%ti2              )) ! [m s^-1]  
    ALLOCATE( laddie%U_prev             ( mesh%ti1:mesh%ti2              )) ! [m s^-1]        2D velocity previous timestep
    ALLOCATE( laddie%V_prev             ( mesh%ti1:mesh%ti2              )) ! [m s^-1]  
    ALLOCATE( laddie%U_next             ( mesh%ti1:mesh%ti2              )) ! [m s^-1]        2D velocity next timestep
    ALLOCATE( laddie%V_next             ( mesh%ti1:mesh%ti2              )) ! [m s^-1]  

    laddie%U              = 0._dp
    laddie%V              = 0._dp
    laddie%U_prev         = 0._dp
    laddie%V_prev         = 0._dp
    laddie%U_next         = 0._dp
    laddie%V_next         = 0._dp

    ! Temperatures
    ALLOCATE( laddie%T                  ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature
    ALLOCATE( laddie%T_prev             ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature previous timestep
    ALLOCATE( laddie%T_next             ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature next timestep
    ALLOCATE( laddie%T_amb              ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature layer bottom
    ALLOCATE( laddie%T_base             ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature ice shelf base
    ALLOCATE( laddie%T_freeze           ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature freezing

    laddie%T              = 0._dp
    laddie%T_prev         = 0._dp
    laddie%T_next         = 0._dp
    laddie%T_amb          = 0._dp
    laddie%T_base         = 0._dp
    laddie%T_freeze       = 0._dp

    ! Salinities
    ALLOCATE( laddie%S                  ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity   
    ALLOCATE( laddie%S_prev             ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity previous timestep
    ALLOCATE( laddie%S_next             ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity next timestep
    ALLOCATE( laddie%S_amb              ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity layer bottom
    ALLOCATE( laddie%S_base             ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity ice shelf base

    laddie%S              = 0._dp
    laddie%S_prev         = 0._dp
    laddie%S_next         = 0._dp
    laddie%S_amb          = 0._dp
    laddie%S_base         = 0._dp

    ! Densities and buoyancies
    ALLOCATE( laddie%rho                ( mesh%vi1:mesh%vi2              )) ! [kg m^-3]       Layer density
    ALLOCATE( laddie%rho_amb            ( mesh%vi1:mesh%vi2              )) ! [kg m^-3]       Ambient water density
    ALLOCATE( laddie%drho_amb           ( mesh%vi1:mesh%vi2              )) ! []              Buoyancy at layer bottom
    ALLOCATE( laddie%drho_base          ( mesh%vi1:mesh%vi2              )) ! []              Buoyancy at ice base

    laddie%rho            = 0._dp
    laddie%rho_amb        = 0._dp
    laddie%drho_amb       = 0._dp
    laddie%drho_base      = 0._dp

    ! Friction velocity
    ALLOCATE( laddie%u_star             ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Friction velocity

    laddie%u_star         = 0._dp

    ! Physical parameter fields
    ALLOCATE( laddie%gamma_T            ( mesh%vi1:mesh%vi2              )) ! []              Turbulent heat exchange coefficient
    ALLOCATE( laddie%gamma_S            ( mesh%vi1:mesh%vi2              )) ! []              Turbulent salt exchange coefficient
    ALLOCATE( laddie%A_h                ( mesh%vi1:mesh%vi2              )) ! [m^2 s^-1]      Horizontal laplacian viscosity
    ALLOCATE( laddie%K_h                ( mesh%vi1:mesh%vi2              )) ! [m^2 s^-1]      Horizontal diffusivity

    laddie%gamma_T        = 0._dp
    laddie%gamma_S        = 0._dp
    laddie%A_h            = 0._dp
    laddie%K_h            = 0._dp

    ! Vertical rates
    ALLOCATE( laddie%melt               ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Melting / freezing rate
    ALLOCATE( laddie%entr               ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Entrainment
    ALLOCATE( laddie%entr_dmin          ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Entrainment for D_min
    ALLOCATE( laddie%detr               ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Detrainment
    ALLOCATE( laddie%entr_tot           ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Total (net) entrainment

    laddie%melt           = 0._dp
    laddie%entr           = 0._dp
    laddie%entr_dmin      = 0._dp
    laddie%detr           = 0._dp
    laddie%entr_tot       = 0._dp

    ! Horizontal fluxes
    ALLOCATE( laddie%divQ               ( mesh%vi1:mesh%vi2              )) ! [m^3 s^-1]      Divergence of layer thickness
    ALLOCATE( laddie%divQU              ( mesh%ti1:mesh%ti2              )) ! [m^4 s^-2]      Divergence of momentum
    ALLOCATE( laddie%divQV              ( mesh%ti1:mesh%ti2              )) ! [m^4 s^-2]   
    ALLOCATE( laddie%divQT              ( mesh%vi1:mesh%vi2              )) ! [degC m^3 s^-1] Divergence of heat
    ALLOCATE( laddie%divQS              ( mesh%vi1:mesh%vi2              )) ! [PSU m^3 s^-1]  Divergence of salt

    laddie%divQ           = 0._dp
    laddie%divQU          = 0._dp
    laddie%divQV          = 0._dp
    laddie%divQT          = 0._dp
    laddie%divQS          = 0._dp

    ! Viscosities
    ALLOCATE( laddie%viscU              ( mesh%ti1:mesh%ti2              )) ! [m^2 s^-2]      Horizontal viscosity term
    ALLOCATE( laddie%viscV              ( mesh%ti1:mesh%ti2              )) ! [m^2 s^-2]      

    laddie%viscU          = 0._dp
    laddie%viscV          = 0._dp

    ! Diffusivities
    ALLOCATE( laddie%diffT              ( mesh%ti1:mesh%ti2              )) ! [degC m s^-1]   Horizontal diffusivity of heat
    ALLOCATE( laddie%diffS              ( mesh%ti1:mesh%ti2              )) ! [PSU m s^-1]    Horizontal diffusivity of salt

    laddie%diffT          = 0._dp
    laddie%diffS          = 0._dp

    ! RHS terms
    ALLOCATE( laddie%ddrho_amb_dx_b     ( mesh%ti1:mesh%ti2              )) ! [m^-1]          Horizontal derivative of buoyancy
    ALLOCATE( laddie%ddrho_amb_dy_b     ( mesh%ti1:mesh%ti2              )) ! [m^-1]          
    ALLOCATE( laddie%dHib_dx_b          ( mesh%ti1:mesh%ti2              )) ! [m^-2]          Horizontal derivative of ice draft
    ALLOCATE( laddie%dHib_dy_b          ( mesh%ti1:mesh%ti2              )) ! [m^-2]          

    laddie%ddrho_amb_dx_b = 0._dp
    laddie%ddrho_amb_dy_b = 0._dp
    laddie%dHib_dx_b      = 0._dp
    laddie%dHib_dy_b      = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_laddie_model 

END MODULE laddie_utilities

