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
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE ocean_utilities                                        , ONLY: interpolate_ocean_depth
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_ambient_TS( mesh, ice, ocean, laddie, Hstar)
    ! Compute T and S of ambient ocean water at the depth of LADDIE's layer bottom
    ! through vertical interpolation

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_ambient_TS'
    INTEGER                                               :: vi
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get T and S at layer base
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T( vi,:), Hstar( vi) - ice%Hib( vi), laddie%T_amb( vi))
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S( vi,:), Hstar( vi) - ice%Hib( vi), laddie%S_amb( vi))
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_ambient_TS

  SUBROUTINE map_H_a_b( mesh, laddie, H_a, H_b)
    ! Map layer thickness from a to b grid, accounting for BCs

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(IN)    :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: H_a
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(INOUT) :: H_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'map_H_a_b'
    INTEGER                                               :: i, vi, ti, n
    REAL(dp), DIMENSION(mesh%nV)                          :: H_a_tot
 
    ! Add routine to path
    CALL init_routine( routine_name)


    CALL gather_to_all_dp_1D( H_a, H_a_tot)

    ! Get T and S at layer base
    DO ti = mesh%ti1, mesh%ti2
       H_b( ti) = 0.0_dp
       IF (laddie%mask_b( ti)) THEN
         
         ! Set zero
         n = 0

         ! Loop over vertices
         DO i = 1, 3
           vi = mesh%Tri( ti, i)
           IF (laddie%mask_a( vi)) THEN
             H_b( ti) = H_b( ti) + H_a_tot( vi)
             n = n + 1
           END IF
         END DO

         H_b( ti) = H_b( ti) / n

       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_H_a_b

  SUBROUTINE map_H_a_c( mesh, laddie, H_a, H_c)
    ! Map layer thickness from a to c grid, accounting for BCs

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(IN)    :: laddie
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: H_a
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(INOUT) :: H_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'map_H_a_c'
    INTEGER                                               :: i, vi, ei, n
    REAL(dp), DIMENSION(mesh%nV)                          :: H_a_tot
 
    ! Add routine to path
    CALL init_routine( routine_name)


    CALL gather_to_all_dp_1D( H_a, H_a_tot)

    DO ei = mesh%ei1, mesh%ei2
       H_c( ei) = 0.0_dp
       ! Set zero
       n = 0
       
       ! Loop over vertices
       DO i = 1, 2
         vi = mesh%EV( ei, i)
         IF (laddie%mask_a( vi)) THEN
           H_c( ei) = H_c( ei) + H_a_tot( vi)
           n = n+1
         END IF
       END DO

       IF (n==0) THEN
         H_c( ei) = 0
       ELSE
         H_c( ei) = H_c( ei) / n
       END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_H_a_c

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
    ALLOCATE( laddie%dH_dt              ( mesh%vi1:mesh%vi2              )) ! [m]             change

    laddie%dH_dt          = 0._dp

    ! Temperatures
    ALLOCATE( laddie%T_amb              ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature layer bottom
    ALLOCATE( laddie%T_base             ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature ice shelf base
    ALLOCATE( laddie%T_freeze           ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature freezing

    laddie%T_amb          = 0._dp
    laddie%T_base         = 0._dp
    laddie%T_freeze       = 0._dp

    ! Salinities
    ALLOCATE( laddie%S_amb              ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity layer bottom
    ALLOCATE( laddie%S_base             ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity ice shelf base

    laddie%S_amb          = 0._dp
    laddie%S_base         = 0._dp

    ! Densities and buoyancies
    ALLOCATE( laddie%rho                ( mesh%vi1:mesh%vi2              )) ! [kg m^-3]       Layer density
    ALLOCATE( laddie%rho_amb            ( mesh%vi1:mesh%vi2              )) ! [kg m^-3]       Ambient water density
    ALLOCATE( laddie%drho_amb           ( mesh%vi1:mesh%vi2              )) ! []              Buoyancy at layer bottom
    ALLOCATE( laddie%Hdrho_amb          ( mesh%vi1:mesh%vi2              )) ! []              Depth-integrated buoyancy at layer bottom
    ALLOCATE( laddie%Hdrho_amb_b        ( mesh%ti1:mesh%ti2              )) ! []              Depth-integrated buoyancy at layer bottom
    ALLOCATE( laddie%drho_base          ( mesh%vi1:mesh%vi2              )) ! []              Buoyancy at ice base

    laddie%rho            = 0._dp
    laddie%rho_amb        = 0._dp
    laddie%drho_amb       = 0._dp
    laddie%Hdrho_amb      = 0._dp
    laddie%Hdrho_amb_b    = 0._dp
    laddie%drho_base      = 0._dp

    ! Friction velocity
    ALLOCATE( laddie%u_star             ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Friction velocity

    laddie%u_star         = 0._dp

    ! Physical parameter fields
    ALLOCATE( laddie%gamma_T            ( mesh%vi1:mesh%vi2              )) ! []              Turbulent heat exchange coefficient
    ALLOCATE( laddie%gamma_S            ( mesh%vi1:mesh%vi2              )) ! []              Turbulent salt exchange coefficient
    ALLOCATE( laddie%A_h                ( mesh%ti1:mesh%ti2              )) ! [m^2 s^-1]      Horizontal laplacian viscosity
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
    ALLOCATE( laddie%divQH              ( mesh%vi1:mesh%vi2              )) ! [m^3 s^-1]      Divergence of layer thickness
    ALLOCATE( laddie%divQU              ( mesh%ti1:mesh%ti2              )) ! [m^4 s^-2]      Divergence of momentum
    ALLOCATE( laddie%divQV              ( mesh%ti1:mesh%ti2              )) ! [m^4 s^-2]   
    ALLOCATE( laddie%divQT              ( mesh%vi1:mesh%vi2              )) ! [degC m^3 s^-1] Divergence of heat
    ALLOCATE( laddie%divQS              ( mesh%vi1:mesh%vi2              )) ! [PSU m^3 s^-1]  Divergence of salt

    laddie%divQH          = 0._dp
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
    ALLOCATE( laddie%diffT              ( mesh%vi1:mesh%vi2              )) ! [degC m s^-1]   Horizontal diffusivity of heat
    ALLOCATE( laddie%diffS              ( mesh%vi1:mesh%vi2              )) ! [PSU m s^-1]    Horizontal diffusivity of salt

    laddie%diffT          = 0._dp
    laddie%diffS          = 0._dp

    ! RHS terms
    ALLOCATE( laddie%ddrho_amb_dx_b     ( mesh%ti1:mesh%ti2              )) ! [m^-1]          Horizontal derivative of buoyancy
    ALLOCATE( laddie%ddrho_amb_dy_b     ( mesh%ti1:mesh%ti2              )) ! [m^-1]          
    ALLOCATE( laddie%dHib_dx_b          ( mesh%ti1:mesh%ti2              )) ! [m^-2]          Horizontal derivative of ice draft
    ALLOCATE( laddie%dHib_dy_b          ( mesh%ti1:mesh%ti2              )) ! [m^-2]          
    ALLOCATE( laddie%dH_dx_b            ( mesh%ti1:mesh%ti2              )) ! [m^-2]          Horizontal derivative of thickness
    ALLOCATE( laddie%dH_dy_b            ( mesh%ti1:mesh%ti2              )) ! [m^-2]          
    ALLOCATE( laddie%detr_b             ( mesh%ti1:mesh%ti2              )) ! [m s^-1]        Detrainment on b grid

    laddie%ddrho_amb_dx_b = 0._dp
    laddie%ddrho_amb_dy_b = 0._dp
    laddie%dHib_dx_b      = 0._dp
    laddie%dHib_dy_b      = 0._dp
    laddie%dH_dx_b        = 0._dp
    laddie%dH_dy_b        = 0._dp
    laddie%detr_b         = 0._dp

    ! Mapped main variables
    ALLOCATE( laddie%H_c                ( mesh%ei1:mesh%ei2              )) ! [m]             Layer thickness on c grid 

    laddie%H_c            = 0._dp

    ! Masks

    ALLOCATE( laddie%mask_a             ( mesh%vi1:mesh%vi2              )) !                 Mask on a-grid
    ALLOCATE( laddie%mask_gr_a          ( mesh%vi1:mesh%vi2              )) !                 Grounded mask on a-grid
    ALLOCATE( laddie%mask_oc_a          ( mesh%vi1:mesh%vi2              )) !                 Icefree ocean mask on a-grid
    ALLOCATE( laddie%mask_b             ( mesh%ti1:mesh%ti2              )) !                 Mask on b-grid
    ALLOCATE( laddie%mask_gl_b          ( mesh%ti1:mesh%ti2              )) !                 Grounding line mask on b-grid
    ALLOCATE( laddie%mask_cf_b          ( mesh%ti1:mesh%ti2              )) !                 Calving front mask on b-grid
    ALLOCATE( laddie%mask_oc_b          ( mesh%ti1:mesh%ti2              )) !                 Icefree ocean mask on b-grid

    laddie%mask_a         = .false.
    laddie%mask_gr_a      = .false.
    laddie%mask_oc_a      = .false.
    laddie%mask_b         = .false.
    laddie%mask_gl_b      = .false.
    laddie%mask_cf_b      = .false.
    laddie%mask_oc_b      = .false.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_laddie_model 

  SUBROUTINE allocate_laddie_timestep( mesh, npx)
    ! Allocate variables of the laddie model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'allocate_laddie_timestep'
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( npx%H                  ( mesh%vi1:mesh%vi2              )) ! [m]             Layer thickness
    ALLOCATE( npx%H_b                ( mesh%ti1:mesh%ti2              )) ! [m]             Layer thickness on b grid
    ALLOCATE( npx%H_c                ( mesh%ei1:mesh%ei2              )) ! [m]             Layer thickness on c grid
    ALLOCATE( npx%U                  ( mesh%ti1:mesh%ti2              )) ! [m s^-1]        2D velocity
    ALLOCATE( npx%U_a                ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        2D velocity on a grid
    ALLOCATE( npx%U_c                ( mesh%ei1:mesh%ei2              )) ! [m s^-1]        2D velocity on b grid
    ALLOCATE( npx%V                  ( mesh%ti1:mesh%ti2              )) ! [m s^-1]  
    ALLOCATE( npx%V_a                ( mesh%vi1:mesh%vi2              )) ! [m s^-1]  
    ALLOCATE( npx%V_c                ( mesh%ei1:mesh%ei2              )) ! [m s^-1]  
    ALLOCATE( npx%T                  ( mesh%vi1:mesh%vi2              )) ! [degC]          Temperature
    ALLOCATE( npx%S                  ( mesh%vi1:mesh%vi2              )) ! [PSU]           Salinity   

    npx%H              = 0._dp
    npx%H_b            = 0._dp
    npx%H_c            = 0._dp
    npx%U              = 0._dp
    npx%U_a            = 0._dp
    npx%U_c            = 0._dp
    npx%V              = 0._dp
    npx%V_a            = 0._dp
    npx%V_c            = 0._dp
    npx%T              = 0._dp
    npx%S              = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_laddie_timestep

END MODULE laddie_utilities

