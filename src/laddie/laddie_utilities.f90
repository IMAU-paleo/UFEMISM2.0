MODULE laddie_utilities

  ! Utilities for the laddie model

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE ocean_utilities                                        , ONLY: interpolate_ocean_depth
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D
  use mesh_integrate_over_domain, only: average_over_domain, calc_and_print_min_mean_max
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD
  use mpi_distributed_shared_memory, only: allocate_dist_shared

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
    REAL(dp), DIMENSION(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), INTENT(IN)    :: Hstar

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
    ! call calc_and_print_min_mean_max( mesh, laddie%T_amb, 'laddie%T_amb')
    ! call calc_and_print_min_mean_max( mesh, laddie%S_amb, 'laddie%S_amb')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_ambient_TS

  subroutine map_H_a_b( mesh, laddie, H_a, H_b)
    ! Map layer thickness from a to b grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(:),                 intent(in)    :: H_a
    real(dp), dimension(:),                 intent(inout) :: H_b

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_b'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_b, &
      mesh%pai_V, H_a, mesh%pai_Tri, H_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_H_a_b

  subroutine map_H_a_c( mesh, laddie, H_a, H_c)
    ! Map layer thickness from a to c grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(:),                 intent(in)    :: H_a
    real(dp), dimension(:),                 intent(inout) :: H_c

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_c'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_c, &
      mesh%pai_V, H_a, mesh%pai_E, H_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_H_a_c

  SUBROUTINE allocate_laddie_model( mesh, laddie)
    ! Allocate variables of the laddie model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'allocate_laddie_model'
    INTEGER                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_per_row_est, nnz_est_proc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Thickness
    call allocate_dist_shared( laddie%dH_dt         , laddie%wdH_dt         , mesh%pai_V%n_nih  )    ! [m]             change
    laddie%dH_dt         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%dH_dt

    ! Temperatures
    call allocate_dist_shared( laddie%T_amb         , laddie%wT_amb         , mesh%pai_V%n_nih  )    ! [degC]          Temperature layer bottom
    call allocate_dist_shared( laddie%T_base        , laddie%wT_base        , mesh%pai_V%n_nih  )    ! [degC]          Temperature ice shelf base
    call allocate_dist_shared( laddie%T_freeze      , laddie%wT_freeze      , mesh%pai_V%n_nih  )    ! [degC]          Temperature freezing
    laddie%T_amb         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%T_amb
    laddie%T_base        ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%T_base
    laddie%T_freeze      ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%T_freeze

    ! Salinities
    call allocate_dist_shared( laddie%S_amb         , laddie%wS_amb         , mesh%pai_V%n_nih  )    ! [PSU]           Salinity layer bottom
    call allocate_dist_shared( laddie%S_base        , laddie%wS_base        , mesh%pai_V%n_nih  )    ! [PSU]           Salinity ice shelf base
    laddie%S_amb         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%S_amb
    laddie%S_base        ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%S_base

    ! Densities and buoyancies
    call allocate_dist_shared( laddie%rho           , laddie%wrho           , mesh%pai_V%n_nih  )    ! [kg m^-3]       Layer density
    call allocate_dist_shared( laddie%rho_amb       , laddie%wrho_amb       , mesh%pai_V%n_nih  )    ! [kg m^-3]       Ambient water density
    call allocate_dist_shared( laddie%drho_amb      , laddie%wdrho_amb      , mesh%pai_V%n_nih  )    ! []              Buoyancy at layer bottom
    call allocate_dist_shared( laddie%Hdrho_amb     , laddie%wHdrho_amb     , mesh%pai_V%n_nih  )    ! []              Depth-integrated buoyancy at layer bottom
    call allocate_dist_shared( laddie%Hdrho_amb_b   , laddie%wHdrho_amb_b   , mesh%pai_Tri%n_nih)    ! []              Depth-integrated buoyancy at layer bottom
    call allocate_dist_shared( laddie%drho_base     , laddie%wdrho_base     , mesh%pai_V%n_nih  )    ! []              Buoyancy at ice base
    laddie%rho           ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%rho
    laddie%rho_amb       ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%rho_amb
    laddie%drho_amb      ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%drho_amb
    laddie%Hdrho_amb     ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%Hdrho_amb
    laddie%Hdrho_amb_b   ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%Hdrho_amb_b
    laddie%drho_base     ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%drho_base

    ! Friction velocity
    call allocate_dist_shared( laddie%u_star        , laddie%wu_star        , mesh%pai_V%n_nih  )    ! [m s^-1]        Friction velocity
    laddie%u_star        ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%u_star

    ! Physical parameter fields
    call allocate_dist_shared( laddie%gamma_T       , laddie%wgamma_T       , mesh%pai_V%n_nih  )    ! []              Turbulent heat exchange coefficient
    call allocate_dist_shared( laddie%gamma_S       , laddie%wgamma_S       , mesh%pai_V%n_nih  )    ! []              Turbulent salt exchange coefficient
    call allocate_dist_shared( laddie%A_h           , laddie%wA_h           , mesh%pai_Tri%n_nih)    ! [m^2 s^-1]      Horizontal laplacian viscosity
    call allocate_dist_shared( laddie%K_h           , laddie%wK_h           , mesh%pai_V%n_nih  )    ! [m^2 s^-1]      Horizontal diffusivity
    laddie%gamma_T       ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%gamma_T
    laddie%gamma_S       ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%gamma_S
    laddie%A_h           ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%A_h
    laddie%K_h           ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%K_h

    ! Vertical rates
    call allocate_dist_shared( laddie%melt          , laddie%wmelt          , mesh%pai_V%n_nih  )    ! [m s^-1]        Melting / freezing rate
    call allocate_dist_shared( laddie%entr          , laddie%wentr          , mesh%pai_V%n_nih  )    ! [m s^-1]        Entrainment
    call allocate_dist_shared( laddie%entr_dmin     , laddie%wentr_dmin     , mesh%pai_V%n_nih  )    ! [m s^-1]        Entrainment for D_min
    call allocate_dist_shared( laddie%detr          , laddie%wdetr          , mesh%pai_V%n_nih  )    ! [m s^-1]        Detrainment
    call allocate_dist_shared( laddie%entr_tot      , laddie%wentr_tot      , mesh%pai_V%n_nih  )    ! [m s^-1]        Total (net) entrainment
    laddie%melt          ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%melt
    laddie%entr          ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%entr
    laddie%entr_dmin     ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%entr_dmin
    laddie%detr          ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%detr
    laddie%entr_tot      ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%entr_tot

    ! Horizontal fluxes
    call allocate_dist_shared( laddie%divQH         , laddie%wdivQH         , mesh%pai_V%n_nih  )    ! [m^3 s^-1]      Divergence of layer thickness
    call allocate_dist_shared( laddie%divQU         , laddie%wdivQU         , mesh%pai_Tri%n_nih)    ! [m^4 s^-2]      Divergence of momentum
    call allocate_dist_shared( laddie%divQV         , laddie%wdivQV         , mesh%pai_Tri%n_nih)    ! [m^4 s^-2]
    call allocate_dist_shared( laddie%divQT         , laddie%wdivQT         , mesh%pai_V%n_nih  )    ! [degC m^3 s^-1] Divergence of heat
    call allocate_dist_shared( laddie%divQS         , laddie%wdivQS         , mesh%pai_V%n_nih  )    ! [PSU m^3 s^-1]  Divergence of salt
    laddie%divQH         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%divQH
    laddie%divQU         ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%divQU
    laddie%divQV         ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%divQV
    laddie%divQT         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%divQT
    laddie%divQS         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%divQS

    ! Viscosities
    call allocate_dist_shared( laddie%viscU         , laddie%wviscU         , mesh%pai_Tri%n_nih)    ! [m^2 s^-2]      Horizontal viscosity term
    call allocate_dist_shared( laddie%viscV         , laddie%wviscV         , mesh%pai_Tri%n_nih)    ! [m^2 s^-2]
    laddie%viscU         ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%viscU
    laddie%viscV         ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%viscV

    ! Diffusivities
    call allocate_dist_shared( laddie%diffT         , laddie%wdiffT         , mesh%pai_V%n_nih  )    ! [degC m s^-1]   Horizontal diffusivity of heat
    call allocate_dist_shared( laddie%diffS         , laddie%wdiffS         , mesh%pai_V%n_nih  )    ! [PSU m s^-1]    Horizontal diffusivity of salt
    laddie%diffT         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%diffT
    laddie%diffS         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%diffS

    ! RHS terms
    call allocate_dist_shared( laddie%ddrho_amb_dx_b, laddie%wddrho_amb_dx_b, mesh%pai_Tri%n_nih)    ! [m^-1]          Horizontal derivative of buoyancy
    call allocate_dist_shared( laddie%ddrho_amb_dy_b, laddie%wddrho_amb_dy_b, mesh%pai_Tri%n_nih)    ! [m^-1]
    call allocate_dist_shared( laddie%dH_dx_b       , laddie%wdH_dx_b       , mesh%pai_Tri%n_nih)    ! [m^-2]          Horizontal derivative of thickness
    call allocate_dist_shared( laddie%dH_dy_b       , laddie%wdH_dy_b       , mesh%pai_Tri%n_nih)    ! [m^-2]
    call allocate_dist_shared( laddie%detr_b        , laddie%wdetr_b        , mesh%pai_Tri%n_nih)    ! [m s^-1]        Detrainment on b grid
    laddie%ddrho_amb_dx_b( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%ddrho_amb_dx_b
    laddie%ddrho_amb_dy_b( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%ddrho_amb_dy_b
    laddie%dH_dx_b       ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%dH_dx_b
    laddie%dH_dy_b       ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%dH_dy_b
    laddie%detr_b        ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%detr_b

    ! Forward-Backward Runge-Kutta 3 scheme
    call allocate_dist_shared( laddie%Hstar         , laddie%wHstar         , mesh%pai_V%n_nih  )    ! [m]               Intermediate layer thickness
    laddie%Hstar         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%Hstar

    ! Mapped variables
    call allocate_dist_shared( laddie%H_c           , laddie%wH_c           , mesh%pai_E%n_nih  )
    call allocate_dist_shared( laddie%Hstar_b       , laddie%wHstar_b       , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( laddie%Hstar_c       , laddie%wHstar_c       , mesh%pai_E%n_nih  )
    laddie%H_c           ( mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih  ) => laddie%H_c
    laddie%Hstar_b       ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%Hstar_b
    laddie%Hstar_c       ( mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih  ) => laddie%Hstar_c

    ! Forward-Backward Runge-Kutta 3 scheme
    ALLOCATE( laddie%Hstar              ( mesh%vi1:mesh%vi2), source=0._dp) ! [m]               Intermediate layer thickness

    ! Masks
    call allocate_dist_shared( laddie%mask_a        , laddie%wmask_a        , mesh%pai_V%n_nih  )    !                 Mask on a-grid
    call allocate_dist_shared( laddie%mask_gr_a     , laddie%wmask_gr_a     , mesh%pai_V%n_nih  )    !                 Grounded mask on a-grid
    call allocate_dist_shared( laddie%mask_oc_a     , laddie%wmask_oc_a     , mesh%pai_V%n_nih  )    !                 Icefree ocean mask on a-grid
    call allocate_dist_shared( laddie%mask_b        , laddie%wmask_b        , mesh%pai_Tri%n_nih)    !                 Mask on b-grid
    call allocate_dist_shared( laddie%mask_gl_b     , laddie%wmask_gl_b     , mesh%pai_Tri%n_nih)    !                 Grounding line mask on b-grid
    call allocate_dist_shared( laddie%mask_cf_b     , laddie%wmask_cf_b     , mesh%pai_Tri%n_nih)    !                 Calving front mask on b-grid
    call allocate_dist_shared( laddie%mask_oc_b     , laddie%wmask_oc_b     , mesh%pai_Tri%n_nih)    !                 Icefree ocean mask on b-grid
    laddie%mask_a        ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%mask_a
    laddie%mask_gr_a     ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%mask_gr_a
    laddie%mask_oc_a     ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%mask_oc_a
    laddie%mask_b        ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%mask_b
    laddie%mask_gl_b     ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%mask_gl_b
    laddie%mask_cf_b     ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%mask_cf_b
    laddie%mask_oc_b     ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => laddie%mask_oc_b

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nTri      ! to
    nrows_loc       = mesh%nTri_loc
    nnz_per_row_est = 3
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( laddie%M_map_H_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Matrix size
    ncols           = mesh%nV        ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nE      ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( laddie%M_map_H_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Matrix size
    ncols           = mesh%nTri        ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nE      ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( laddie%M_map_UV_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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

    call allocate_dist_shared( npx%H  , npx%wH  , mesh%pai_V%n_nih  )   ! [m]             Layer thickness
    call allocate_dist_shared( npx%H_b, npx%wH_b, mesh%pai_Tri%n_nih)   ! [m]             Layer thickness on b grid
    call allocate_dist_shared( npx%H_c, npx%wH_c, mesh%pai_E%n_nih  )   ! [m]             Layer thickness on c grid
    call allocate_dist_shared( npx%U  , npx%wU  , mesh%pai_Tri%n_nih)   ! [m s^-1]        2D velocity
    call allocate_dist_shared( npx%U_a, npx%wU_a, mesh%pai_V%n_nih  )   ! [m s^-1]        2D velocity on a grid
    call allocate_dist_shared( npx%U_c, npx%wU_c, mesh%pai_E%n_nih  )   ! [m s^-1]        2D velocity on b grid
    call allocate_dist_shared( npx%V  , npx%wV  , mesh%pai_Tri%n_nih)   ! [m s^-1]
    call allocate_dist_shared( npx%V_a, npx%wV_a, mesh%pai_V%n_nih  )   ! [m s^-1]
    call allocate_dist_shared( npx%V_c, npx%wV_c, mesh%pai_E%n_nih  )   ! [m s^-1]
    call allocate_dist_shared( npx%T  , npx%wT  , mesh%pai_V%n_nih  )   ! [degC]          Temperature
    call allocate_dist_shared( npx%S  , npx%wS  , mesh%pai_V%n_nih  )   ! [PSU]           Salinity

    npx%H  ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => npx%H
    npx%H_b( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => npx%H_b
    npx%H_c( mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih  ) => npx%H_c
    npx%U  ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => npx%U
    npx%U_a( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => npx%U_a
    npx%U_c( mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih  ) => npx%U_c
    npx%V  ( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => npx%V
    npx%V_a( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => npx%V_a
    npx%V_c( mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih  ) => npx%V_c
    npx%T  ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => npx%T
    npx%S  ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => npx%S

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_laddie_timestep

  subroutine print_diagnostics( mesh, laddie, tl)
    !< Print out diagnostics

    ! In- and output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(in   ) :: laddie
    real(dp),                intent(in   ) :: tl

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_diagnostics'
    real(dp)                       :: H_av, melt_max, U_max, T_max
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call average_over_domain( mesh, laddie%now%H, H_av)

    melt_max = maxval( laddie%melt) * sec_per_year
    call MPI_ALLREDUCE( MPI_IN_PLACE, melt_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    U_max = maxval( sqrt( laddie%now%U**2 + laddie%now%V**2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, U_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    T_max = maxval( laddie%now%T)
    call MPI_ALLREDUCE( MPI_IN_PLACE, T_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (par%primary) write( *, "(F8.3,A,F8.3,A,F8.2,A,F8.3,A,F8.3)") tl/sec_per_day, &
      '  Dmean ', H_av, '  Meltmax', melt_max, '   U', U_max, '   Tmax', T_max

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine print_diagnostics

END MODULE laddie_utilities

