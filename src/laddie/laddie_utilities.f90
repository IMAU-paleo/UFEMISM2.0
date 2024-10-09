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
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, add_empty_row_CSR_dist   
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, gather_to_all_logical_1D

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
       IF (ice%mask_floating_ice( vi)) THEN
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T( vi,:), Hstar( vi) - ice%Hib( vi), laddie%T_amb( vi))
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S( vi,:), Hstar( vi) - ice%Hib( vi), laddie%S_amb( vi))
       END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_ambient_TS

  SUBROUTINE calc_laddie_flux_divergence_matrix_upwind( mesh, U, V, mask_floating, M_divQ)
    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    !
    ! The vertically averaged ice flux divergence represents the net ice volume (which,
    ! assuming constant density, is proportional to the ice mass) entering each Voronoi
    ! cell per unit time. This is found by calculating the ice fluxes through each
    ! shared Voronoi cell boundary, using an upwind scheme: if ice flows from vertex vi
    ! to vertex vj, the flux is found by multiplying the velocity at their shared
    ! boundary u_c with the ice thickness at vi (and, of course, the length L_c of the
    ! shared boundary). If instead it flows from vj to vi, u_c is multiplied with the
    ! ice thickness at vj.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: U
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: V
    LOGICAL, DIMENSION(mesh%vi1:mesh%vi2),  INTENT(IN)    :: mask_floating
    TYPE(type_sparse_matrix_CSR_dp),        INTENT(OUT)   :: M_divQ


    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_laddie_flux_divergence_matrix_upwind'
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2)                :: U_c, V_c
    REAL(dp), DIMENSION(mesh%nE)                          :: U_c_tot, V_c_tot
    INTEGER                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    INTEGER                                               :: vi, ci, ei, vj
    REAL(dp)                                              :: A_i, L_c
    REAL(dp)                                              :: D_x, D_y, D, u_perp
    REAL(dp), DIMENSION(0:mesh%nC_mem)                    :: cM_divQ
    LOGICAL, DIMENSION(mesh%nV)                           :: mask_floating_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    CALL map_laddie_velocities_from_b_to_c_2D( mesh, U, V, U_c, V_c)
    CALL gather_to_all_dp_1D( U_c, U_c_tot)
    CALL gather_to_all_dp_1D( V_c, V_c_tot)
    CALL gather_to_all_logical_1D( mask_floating, mask_floating_tot)

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV      ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_est_proc    = mesh%nV_loc + SUM( mesh%nC( mesh%vi1:mesh%vi2))

    CALL allocate_matrix_CSR_dist( M_divQ, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! == Calculate coefficients
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      IF (mask_floating( vi)) THEN
        ! Initialise
        cM_divQ = 0._dp

        ! Loop over all connections of vertex vi
        DO ci = 1, mesh%nC( vi)

          ! Connection ci from vertex vi leads through edge ei to vertex vj
          ei = mesh%VE( vi,ci)
          vj = mesh%C(  vi,ci)

          ! TODO skip flux into grounded ice. Need mask_grounded_tot( vj)

          ! The Voronoi cell of vertex vi has area A_i
          A_i = mesh%A( vi)

          ! The shared Voronoi cell boundary section between the Voronoi cells
          ! of vertices vi and vj has length L_c
          L_c = mesh%Cw( vi,ci)

          ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
          D_x = mesh%V( vj,1) - mesh%V( vi,1)
          D_y = mesh%V( vj,2) - mesh%V( vi,2)
          D   = SQRT( D_x**2 + D_y**2)
          u_perp = U_c_tot( ei) * D_x/D + V_c_tot( ei) * D_y/D

          ! Calculate matrix coefficients
          ! =============================

          ! u_perp > 0: flow is exiting this vertex into vertex vj
          cM_divQ( 0) = cM_divQ( 0) + L_c * MAX( 0._dp, u_perp) / A_i

          ! u_perp < 0: flow is entering this vertex from vertex vj
          cM_divQ( ci) = L_c * MIN( 0._dp, u_perp) / A_i

        END DO ! DO ci = 1, mesh%nC( vi)

        ! Add coefficients to matrix
        CALL add_entry_CSR_dist( M_divQ, vi, vi, cM_divQ( 0))
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C(  vi,ci)
          CALL add_entry_CSR_dist( M_divQ, vi, vj, cM_divQ( ci))
        END DO ! DO ci = 1, mesh%nC( vi)

      ELSE

        ! Add empty row
        CALL add_empty_row_CSR_dist( M_divQ, vi)

      END IF ! (ice%mask_floating_ice( vi))

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_laddie_flux_divergence_matrix_upwind

  SUBROUTINE map_laddie_velocities_from_b_to_c_2D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    ! Calculate velocities on the c-grid for solving the ice thickness equation
    ! 
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive
        
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_b_partial
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_b_partial
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(OUT)   :: u_c
    REAL(dp), DIMENSION(mesh%ei1:mesh%ei2), INTENT(OUT)   :: v_c
      
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'map_laddie_velocities_from_b_to_c_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE               :: u_b_tot, v_b_tot
    INTEGER                                               :: ei, til, tir
      
    ! Add routine to path
    CALL init_routine( routine_name)
        
    ! Allocate memory
    ALLOCATE( u_b_tot( mesh%nTri))
    ALLOCATE( v_b_tot( mesh%nTri))
        
    ! Gather the full b-grid velocity fields to all processes
    CALL gather_to_all_dp_1D( u_b_partial, u_b_tot)
    CALL gather_to_all_dp_1D( v_b_partial, v_b_tot)

    ! Map velocities from the b-grid (triangles) to the c-grid (edges)
    DO ei = mesh%ei1, mesh%ei2

      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      IF     (til == 0 .AND. tir > 0) THEN
        u_c( ei) = u_b_tot( tir)
        v_c( ei) = v_b_tot( tir)
      ELSEIF (tir == 0 .AND. til > 0) THEN
        u_c( ei) = u_b_tot( til)
        v_c( ei) = v_b_tot( til)
      ELSEIF (til >  0 .AND. tir > 0) THEN
        u_c( ei) = (u_b_tot( til) + u_b_tot( tir)) / 2._dp
        v_c( ei) = (v_b_tot( til) + v_b_tot( tir)) / 2._dp
      ELSE
        CALL crash('something is seriously wrong with the ETri array of this mesh!')
      END IF

    END DO

    ! Clean up after yourself
    DEALLOCATE( u_b_tot)
    DEALLOCATE( v_b_tot)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_laddie_velocities_from_b_to_c_2D


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
    ALLOCATE( laddie%H_b                ( mesh%ti1:mesh%ti2              )) ! [m]             Layer thickness on b grid
    ALLOCATE( laddie%H_b_next           ( mesh%ti1:mesh%ti2              )) ! [m]             Layer next thickness on b grid
    ALLOCATE( laddie%U_a                ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        Layer velocity on a grid  
    ALLOCATE( laddie%V_a                ( mesh%vi1:mesh%vi2              )) ! [m s^-1]        

    laddie%H_b            = 0._dp
    laddie%H_b_next       = 0._dp
    laddie%U_a            = 0._dp
    laddie%V_a            = 0._dp

    ! Masks

    ALLOCATE( laddie%mask_a             ( mesh%vi1:mesh%vi2              )) !                 Mask on a-grid
    ALLOCATE( laddie%mask_b             ( mesh%ti1:mesh%ti2              )) !                 Mask on b-grid

    laddie%mask_a         = .false.
    laddie%mask_b         = .false.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_laddie_model 

END MODULE laddie_utilities

