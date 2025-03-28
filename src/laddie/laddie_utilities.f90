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
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D
  use mesh_utilities, only: average_over_domain
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
    REAL(dp), DIMENSION(mesh%vi1_node:mesh%vi2_node), INTENT(IN)    :: Hstar

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

  subroutine map_H_a_b( mesh, laddie, H_a, H_b)
    ! Map layer thickness from a to b grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(mesh%vi1_node:mesh%vi2_node), intent(in)    :: H_a
    real(dp), dimension(mesh%ti1_node:mesh%ti2_node), intent(inout) :: H_b

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_b'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_b, H_a, H_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      xx_tot_buf = mesh%d_a_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_H_a_b

  subroutine map_H_a_c( mesh, laddie, H_a, H_c)
    ! Map layer thickness from a to c grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(mesh%vi1_node:mesh%vi2_node), intent(in)    :: H_a
    real(dp), dimension(mesh%ei1_node:mesh%ei2_node), intent(inout) :: H_c

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_c'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_c, H_a, H_c, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      xx_tot_buf = mesh%d_a_tot)

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
    call allocate_dist_shared(  laddie%dH_dt         , laddie%wdH_dt         , mesh%nV_node  ) ! [m]             change

    laddie%dH_dt         (mesh%vi1_node:mesh%vi2_node) => laddie%dH_dt

    ! Temperatures
    call allocate_dist_shared(  laddie%T_amb         , laddie%wT_amb         , mesh%nV_node  ) ! [degC]          Temperature layer bottom
    call allocate_dist_shared(  laddie%T_base        , laddie%wT_base        , mesh%nV_node  ) ! [degC]          Temperature ice shelf base
    call allocate_dist_shared(  laddie%T_freeze      , laddie%wT_freeze      , mesh%nV_node  ) ! [degC]          Temperature freezing

    laddie%T_amb         (mesh%vi1_node:mesh%vi2_node) => laddie%T_amb
    laddie%T_base        (mesh%vi1_node:mesh%vi2_node) => laddie%T_base
    laddie%T_freeze      (mesh%vi1_node:mesh%vi2_node) => laddie%T_freeze

    ! Salinities
    call allocate_dist_shared(  laddie%S_amb         , laddie%wS_amb         , mesh%nV_node  ) ! [PSU]           Salinity layer bottom
    call allocate_dist_shared(  laddie%S_base        , laddie%wS_base        , mesh%nV_node  ) ! [PSU]           Salinity ice shelf base

    laddie%S_amb         (mesh%vi1_node:mesh%vi2_node) => laddie%S_amb
    laddie%S_base        (mesh%vi1_node:mesh%vi2_node) => laddie%S_base

    ! Densities and buoyancies
    call allocate_dist_shared(  laddie%rho           , laddie%wrho           , mesh%nV_node  ) ! [kg m^-3]       Layer density
    call allocate_dist_shared(  laddie%rho_amb       , laddie%wrho_amb       , mesh%nV_node  ) ! [kg m^-3]       Ambient water density
    call allocate_dist_shared(  laddie%drho_amb      , laddie%wdrho_amb      , mesh%nV_node  ) ! []              Buoyancy at layer bottom
    call allocate_dist_shared(  laddie%Hdrho_amb     , laddie%wHdrho_amb     , mesh%nV_node  ) ! []              Depth-integrated buoyancy at layer bottom
    call allocate_dist_shared(  laddie%Hdrho_amb_b   , laddie%wHdrho_amb_b   , mesh%nTri_node) ! []              Depth-integrated buoyancy at layer bottom
    call allocate_dist_shared(  laddie%drho_base     , laddie%wdrho_base     , mesh%nV_node  ) ! []              Buoyancy at ice base

    laddie%rho           (mesh%vi1_node:mesh%vi2_node) => laddie%rho
    laddie%rho_amb       (mesh%vi1_node:mesh%vi2_node) => laddie%rho_amb
    laddie%drho_amb      (mesh%vi1_node:mesh%vi2_node) => laddie%drho_amb
    laddie%Hdrho_amb     (mesh%vi1_node:mesh%vi2_node) => laddie%Hdrho_amb
    laddie%Hdrho_amb_b   (mesh%ti1_node:mesh%ti2_node) => laddie%Hdrho_amb_b
    laddie%drho_base     (mesh%vi1_node:mesh%vi2_node) => laddie%drho_base

    ! Friction velocity
    call allocate_dist_shared(  laddie%u_star        , laddie%wu_star        , mesh%nV_node  ) ! [m s^-1]        Friction velocity

    laddie%u_star        (mesh%vi1_node:mesh%vi2_node) => laddie%u_star

    ! Physical parameter fields
    call allocate_dist_shared(  laddie%gamma_T       , laddie%wgamma_T       , mesh%nV_node  ) ! []              Turbulent heat exchange coefficient
    call allocate_dist_shared(  laddie%gamma_S       , laddie%wgamma_S       , mesh%nV_node  ) ! []              Turbulent salt exchange coefficient
    call allocate_dist_shared(  laddie%A_h           , laddie%wA_h           , mesh%nTri_node) ! [m^2 s^-1]      Horizontal laplacian viscosity
    call allocate_dist_shared(  laddie%K_h           , laddie%wK_h           , mesh%nV_node  ) ! [m^2 s^-1]      Horizontal diffusivity

    laddie%gamma_T       (mesh%vi1_node:mesh%vi2_node) => laddie%gamma_T
    laddie%gamma_S       (mesh%vi1_node:mesh%vi2_node) => laddie%gamma_S
    laddie%A_h           (mesh%ti1_node:mesh%ti2_node) => laddie%A_h
    laddie%K_h           (mesh%vi1_node:mesh%vi2_node) => laddie%K_h

    ! Vertical rates
    call allocate_dist_shared(  laddie%melt          , laddie%wmelt          , mesh%nV_node  ) ! [m s^-1]        Melting / freezing rate
    call allocate_dist_shared(  laddie%entr          , laddie%wentr          , mesh%nV_node  ) ! [m s^-1]        Entrainment
    call allocate_dist_shared(  laddie%entr_dmin     , laddie%wentr_dmin     , mesh%nV_node  ) ! [m s^-1]        Entrainment for D_min
    call allocate_dist_shared(  laddie%detr          , laddie%wdetr          , mesh%nV_node  ) ! [m s^-1]        Detrainment
    call allocate_dist_shared(  laddie%entr_tot      , laddie%wentr_tot      , mesh%nV_node  ) ! [m s^-1]        Total (net) entrainment

    laddie%melt          (mesh%vi1_node:mesh%vi2_node) => laddie%melt
    laddie%entr          (mesh%vi1_node:mesh%vi2_node) => laddie%entr
    laddie%entr_dmin     (mesh%vi1_node:mesh%vi2_node) => laddie%entr_dmin
    laddie%detr          (mesh%vi1_node:mesh%vi2_node) => laddie%detr
    laddie%entr_tot      (mesh%vi1_node:mesh%vi2_node) => laddie%entr_tot

    ! Horizontal fluxes
    call allocate_dist_shared(  laddie%divQH         , laddie%wdivQH         , mesh%nV_node  ) ! [m^3 s^-1]      Divergence of layer thickness
    call allocate_dist_shared(  laddie%divQU         , laddie%wdivQU         , mesh%nTri_node) ! [m^4 s^-2]      Divergence of momentum
    call allocate_dist_shared(  laddie%divQV         , laddie%wdivQV         , mesh%nTri_node) ! [m^4 s^-2]
    call allocate_dist_shared(  laddie%divQT         , laddie%wdivQT         , mesh%nV_node  ) ! [degC m^3 s^-1] Divergence of heat
    call allocate_dist_shared(  laddie%divQS         , laddie%wdivQS         , mesh%nV_node  ) ! [PSU m^3 s^-1]  Divergence of salt

    laddie%divQH         (mesh%vi1_node:mesh%vi2_node) => laddie%divQH
    laddie%divQU         (mesh%ti1_node:mesh%ti2_node) => laddie%divQU
    laddie%divQV         (mesh%ti1_node:mesh%ti2_node) => laddie%divQV
    laddie%divQT         (mesh%vi1_node:mesh%vi2_node) => laddie%divQT
    laddie%divQS         (mesh%vi1_node:mesh%vi2_node) => laddie%divQS

    ! Viscosities
    call allocate_dist_shared(  laddie%viscU         , laddie%wviscU         , mesh%nTri_node) ! [m^2 s^-2]      Horizontal viscosity term
    call allocate_dist_shared(  laddie%viscV         , laddie%wviscV         , mesh%nTri_node) ! [m^2 s^-2]

    laddie%viscU         (mesh%ti1_node:mesh%ti2_node) => laddie%viscU
    laddie%viscV         (mesh%ti1_node:mesh%ti2_node) => laddie%viscV

    ! Diffusivities
    call allocate_dist_shared(  laddie%diffT         , laddie%wdiffT         , mesh%nV_node  ) ! [degC m s^-1]   Horizontal diffusivity of heat
    call allocate_dist_shared(  laddie%diffS         , laddie%wdiffS         , mesh%nV_node  ) ! [PSU m s^-1]    Horizontal diffusivity of salt

    laddie%diffT         (mesh%vi1_node:mesh%vi2_node) => laddie%diffT
    laddie%diffS         (mesh%vi1_node:mesh%vi2_node) => laddie%diffS

    ! RHS terms
    call allocate_dist_shared(  laddie%ddrho_amb_dx_b, laddie%wddrho_amb_dx_b, mesh%nTri_node) ! [m^-1]          Horizontal derivative of buoyancy
    call allocate_dist_shared(  laddie%ddrho_amb_dy_b, laddie%wddrho_amb_dy_b, mesh%nTri_node) ! [m^-1]
    call allocate_dist_shared(  laddie%dH_dx_b       , laddie%wdH_dx_b       , mesh%nTri_node) ! [m^-2]          Horizontal derivative of thickness
    call allocate_dist_shared(  laddie%dH_dy_b       , laddie%wdH_dy_b       , mesh%nTri_node) ! [m^-2]
    call allocate_dist_shared(  laddie%detr_b        , laddie%wdetr_b        , mesh%nTri_node) ! [m s^-1]        Detrainment on b grid

    laddie%ddrho_amb_dx_b(mesh%ti1_node:mesh%ti2_node) => laddie%ddrho_amb_dx_b
    laddie%ddrho_amb_dy_b(mesh%ti1_node:mesh%ti2_node) => laddie%ddrho_amb_dy_b
    laddie%dH_dx_b       (mesh%ti1_node:mesh%ti2_node) => laddie%dH_dx_b
    laddie%dH_dy_b       (mesh%ti1_node:mesh%ti2_node) => laddie%dH_dy_b
    laddie%detr_b        (mesh%ti1_node:mesh%ti2_node) => laddie%detr_b

    ! Masks
    call allocate_dist_shared(  laddie%mask_a        , laddie%wmask_a        , mesh%nV_node  ) !                 Mask on a-grid
    call allocate_dist_shared(  laddie%mask_gr_a     , laddie%wmask_gr_a     , mesh%nV_node  ) !                 Grounded mask on a-grid
    call allocate_dist_shared(  laddie%mask_oc_a     , laddie%wmask_oc_a     , mesh%nV_node  ) !                 Icefree ocean mask on a-grid
    call allocate_dist_shared(  laddie%mask_b        , laddie%wmask_b        , mesh%nTri_node) !                 Mask on b-grid
    call allocate_dist_shared(  laddie%mask_gl_b     , laddie%wmask_gl_b     , mesh%nTri_node) !                 Grounding line mask on b-grid
    call allocate_dist_shared(  laddie%mask_cf_b     , laddie%wmask_cf_b     , mesh%nTri_node) !                 Calving front mask on b-grid
    call allocate_dist_shared(  laddie%mask_oc_b     , laddie%wmask_oc_b     , mesh%nTri_node) !                 Icefree ocean mask on b-grid

    laddie%mask_a        (mesh%vi1_node:mesh%vi2_node) => laddie%mask_a
    laddie%mask_gr_a     (mesh%vi1_node:mesh%vi2_node) => laddie%mask_gr_a
    laddie%mask_oc_a     (mesh%vi1_node:mesh%vi2_node) => laddie%mask_oc_a
    laddie%mask_b        (mesh%ti1_node:mesh%ti2_node) => laddie%mask_b
    laddie%mask_gl_b     (mesh%ti1_node:mesh%ti2_node) => laddie%mask_gl_b
    laddie%mask_cf_b     (mesh%ti1_node:mesh%ti2_node) => laddie%mask_cf_b
    laddie%mask_oc_b     (mesh%ti1_node:mesh%ti2_node) => laddie%mask_oc_b

    ! Total data fields
    call allocate_dist_shared( laddie%U_tot        , laddie%wU_tot        , mesh%nTri)
    call allocate_dist_shared( laddie%V_tot        , laddie%wV_tot        , mesh%nTri)
    call allocate_dist_shared( laddie%H_b_tot      , laddie%wH_b_tot      , mesh%nTri)
    call allocate_dist_shared( laddie%U_c_tot      , laddie%wU_c_tot      , mesh%nE)
    call allocate_dist_shared( laddie%V_c_tot      , laddie%wV_c_tot      , mesh%nE)
    call allocate_dist_shared( laddie%mask_gl_b_tot, laddie%wmask_gl_b_tot, mesh%nTri)
    call allocate_dist_shared( laddie%mask_a_tot   , laddie%wmask_a_tot   , mesh%nV)
    call allocate_dist_shared( laddie%mask_b_tot   , laddie%wmask_b_tot   , mesh%nTri)
    call allocate_dist_shared( laddie%mask_oc_b_tot, laddie%wmask_oc_b_tot, mesh%nTri)
    call allocate_dist_shared( laddie%H_c_tot      , laddie%wH_c_tot      , mesh%nE)
    call allocate_dist_shared( laddie%H_tot        , laddie%wH_tot        , mesh%nV)
    call allocate_dist_shared( laddie%mask_gr_a_tot, laddie%wmask_gr_a_tot, mesh%nV)
    call allocate_dist_shared( laddie%mask_oc_a_tot, laddie%wmask_oc_a_tot, mesh%nV)
    call allocate_dist_shared( laddie%T_tot        , laddie%wT_tot        , mesh%nV)
    call allocate_dist_shared( laddie%S_tot        , laddie%wS_tot        , mesh%nV)
    call allocate_dist_shared( laddie%Hstar_tot    , laddie%wHstar_tot    , mesh%nV)

    call allocate_dist_shared( laddie%Hstar_b      , laddie%wHstar_b      , mesh%nTri_node)

    laddie%Hstar_b( mesh%ti1_node:mesh%ti2_node) => laddie%Hstar_b

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

    call allocate_dist_shared( npx%H  , npx%wH  , mesh%nV_node  ) ! [m]             Layer thickness
    call allocate_dist_shared( npx%H_b, npx%wH_b, mesh%nTri_node) ! [m]             Layer thickness on b grid
    call allocate_dist_shared( npx%H_c, npx%wH_c, mesh%nE_node  ) ! [m]             Layer thickness on c grid
    call allocate_dist_shared( npx%U  , npx%wU  , mesh%nTri_node) ! [m s^-1]        2D velocity
    call allocate_dist_shared( npx%U_a, npx%wU_a, mesh%nV_node  ) ! [m s^-1]        2D velocity on a grid
    call allocate_dist_shared( npx%U_c, npx%wU_c, mesh%nE_node  ) ! [m s^-1]        2D velocity on b grid
    call allocate_dist_shared( npx%V  , npx%wV  , mesh%nTri_node) ! [m s^-1]
    call allocate_dist_shared( npx%V_a, npx%wV_a, mesh%nV_node  ) ! [m s^-1]
    call allocate_dist_shared( npx%V_c, npx%wV_c, mesh%nE_node  ) ! [m s^-1]
    call allocate_dist_shared( npx%T  , npx%wT  , mesh%nV_node  ) ! [degC]          Temperature
    call allocate_dist_shared( npx%S  , npx%wS  , mesh%nV_node  ) ! [PSU]           Salinity

    npx%H  ( mesh%vi1_node:mesh%vi2_node) => npx%H
    npx%H_b( mesh%ti1_node:mesh%ti2_node) => npx%H_b
    npx%H_c( mesh%ei1_node:mesh%ei2_node) => npx%H_c
    npx%U  ( mesh%ti1_node:mesh%ti2_node) => npx%U
    npx%U_a( mesh%vi1_node:mesh%vi2_node) => npx%U_a
    npx%U_c( mesh%ei1_node:mesh%ei2_node) => npx%U_c
    npx%V  ( mesh%ti1_node:mesh%ti2_node) => npx%V
    npx%V_a( mesh%vi1_node:mesh%vi2_node) => npx%V_a
    npx%V_c( mesh%ei1_node:mesh%ei2_node) => npx%V_c
    npx%T  ( mesh%vi1_node:mesh%vi2_node) => npx%T
    npx%S  ( mesh%vi1_node:mesh%vi2_node) => npx%S

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
    real(dp)                       :: H_av, Meltmax, Umax, Tmax
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call average_over_domain( mesh, laddie%now%H, H_av, d_is_hybrid = .true.)

    Meltmax = maxval( laddie%melt) * sec_per_year
    call MPI_ALLREDUCE( MPI_IN_PLACE, Meltmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    Umax = maxval( sqrt( laddie%now%U**2 + laddie%now%V**2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, Umax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    Tmax = maxval( laddie%now%T)
    call MPI_ALLREDUCE( MPI_IN_PLACE, Tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (par%primary) then
      write( *, "(F8.3,A,F8.3,A,F8.2,A,F8.3,A,F8.3)") tl/sec_per_day, &
        '  Dmean ', H_av, '  Meltmax', Meltmax, '   U', Umax, '   Tmax', Tmax
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine print_diagnostics

END MODULE laddie_utilities

