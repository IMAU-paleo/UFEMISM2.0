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
    real(dp) :: d_av

    ! Add routine to path
    CALL init_routine( routine_name)

    call crash('fixme!')

    ! ! Get T and S at layer base
    ! DO vi = mesh%vi1, mesh%vi2
    !    IF (laddie%mask_a( vi)) THEN
    !      CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T( vi,:), Hstar( vi) - ice%Hib( vi), laddie%T_amb( vi))
    !      CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S( vi,:), Hstar( vi) - ice%Hib( vi), laddie%S_amb( vi))
    !    END IF
    ! END DO

    ! ! DENK DROM
    ! call average_over_domain( mesh, laddie%T_amb, d_av)
    ! if (par%primary) write(0,'(A,F12.8)') ' mean T_amb = ', d_av
    ! call average_over_domain( mesh, laddie%S_amb, d_av)
    ! if (par%primary) write(0,'(A,F12.8)') ' mean S_amb = ', d_av

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_ambient_TS

  subroutine map_H_a_b( mesh, laddie, H_a, H_b)
    ! Map layer thickness from a to b grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)    :: H_a
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: H_b

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_b'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme!')

    ! call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_b, H_a, H_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_H_a_b

  subroutine map_H_a_c( mesh, laddie, H_a, H_c)
    ! Map layer thickness from a to c grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)    :: H_a
    real(dp), dimension(mesh%ei1:mesh%ei2), intent(inout) :: H_c

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_c'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme!')

    ! call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_c, H_a, H_c)

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

    call crash('fixme!')

    ! ! Thickness
    ! ALLOCATE( laddie%dH_dt              ( mesh%vi1:mesh%vi2), source=0._dp) ! [m]             change

    ! ! Temperatures
    ! ALLOCATE( laddie%T_amb              ( mesh%vi1:mesh%vi2), source=0._dp) ! [degC]          Temperature layer bottom
    ! ALLOCATE( laddie%T_base             ( mesh%vi1:mesh%vi2), source=0._dp) ! [degC]          Temperature ice shelf base
    ! ALLOCATE( laddie%T_freeze           ( mesh%vi1:mesh%vi2), source=0._dp) ! [degC]          Temperature freezing

    ! ! Salinities
    ! ALLOCATE( laddie%S_amb              ( mesh%vi1:mesh%vi2), source=0._dp) ! [PSU]           Salinity layer bottom
    ! ALLOCATE( laddie%S_base             ( mesh%vi1:mesh%vi2), source=0._dp) ! [PSU]           Salinity ice shelf base

    ! ! Densities and buoyancies
    ! ALLOCATE( laddie%rho                ( mesh%vi1:mesh%vi2), source=0._dp) ! [kg m^-3]       Layer density
    ! ALLOCATE( laddie%rho_amb            ( mesh%vi1:mesh%vi2), source=0._dp) ! [kg m^-3]       Ambient water density
    ! ALLOCATE( laddie%drho_amb           ( mesh%vi1:mesh%vi2), source=0._dp) ! []              Buoyancy at layer bottom
    ! ALLOCATE( laddie%Hdrho_amb          ( mesh%vi1:mesh%vi2), source=0._dp) ! []              Depth-integrated buoyancy at layer bottom
    ! ALLOCATE( laddie%Hdrho_amb_b        ( mesh%ti1:mesh%ti2), source=0._dp) ! []              Depth-integrated buoyancy at layer bottom
    ! ALLOCATE( laddie%drho_base          ( mesh%vi1:mesh%vi2), source=0._dp) ! []              Buoyancy at ice base

    ! ! Friction velocity
    ! ALLOCATE( laddie%u_star             ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]        Friction velocity

    ! ! Physical parameter fields
    ! ALLOCATE( laddie%gamma_T            ( mesh%vi1:mesh%vi2), source=0._dp) ! []              Turbulent heat exchange coefficient
    ! ALLOCATE( laddie%gamma_S            ( mesh%vi1:mesh%vi2), source=0._dp) ! []              Turbulent salt exchange coefficient
    ! ALLOCATE( laddie%A_h                ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^2 s^-1]      Horizontal laplacian viscosity
    ! ALLOCATE( laddie%K_h                ( mesh%vi1:mesh%vi2), source=0._dp) ! [m^2 s^-1]      Horizontal diffusivity

    ! ! Vertical rates
    ! ALLOCATE( laddie%melt               ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]        Melting / freezing rate
    ! ALLOCATE( laddie%entr               ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]        Entrainment
    ! ALLOCATE( laddie%entr_dmin          ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]        Entrainment for D_min
    ! ALLOCATE( laddie%detr               ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]        Detrainment
    ! ALLOCATE( laddie%entr_tot           ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]        Total (net) entrainment

    ! ! Horizontal fluxes
    ! ALLOCATE( laddie%divQH              ( mesh%vi1:mesh%vi2), source=0._dp) ! [m^3 s^-1]      Divergence of layer thickness
    ! ALLOCATE( laddie%divQU              ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^4 s^-2]      Divergence of momentum
    ! ALLOCATE( laddie%divQV              ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^4 s^-2]
    ! ALLOCATE( laddie%divQT              ( mesh%vi1:mesh%vi2), source=0._dp) ! [degC m^3 s^-1] Divergence of heat
    ! ALLOCATE( laddie%divQS              ( mesh%vi1:mesh%vi2), source=0._dp) ! [PSU m^3 s^-1]  Divergence of salt

    ! ! Viscosities
    ! ALLOCATE( laddie%viscU              ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^2 s^-2]      Horizontal viscosity term
    ! ALLOCATE( laddie%viscV              ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^2 s^-2]

    ! ! Diffusivities
    ! ALLOCATE( laddie%diffT              ( mesh%vi1:mesh%vi2), source=0._dp) ! [degC m s^-1]   Horizontal diffusivity of heat
    ! ALLOCATE( laddie%diffS              ( mesh%vi1:mesh%vi2), source=0._dp) ! [PSU m s^-1]    Horizontal diffusivity of salt

    ! ! RHS terms
    ! ALLOCATE( laddie%ddrho_amb_dx_b     ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^-1]          Horizontal derivative of buoyancy
    ! ALLOCATE( laddie%ddrho_amb_dy_b     ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^-1]
    ! ALLOCATE( laddie%dH_dx_b            ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^-2]          Horizontal derivative of thickness
    ! ALLOCATE( laddie%dH_dy_b            ( mesh%ti1:mesh%ti2), source=0._dp) ! [m^-2]
    ! ALLOCATE( laddie%detr_b             ( mesh%ti1:mesh%ti2), source=0._dp) ! [m s^-1]        Detrainment on b grid

    ! ! Masks
    ! ALLOCATE( laddie%mask_a             ( mesh%vi1:mesh%vi2), source=.false.) !                 Mask on a-grid
    ! ALLOCATE( laddie%mask_gr_a          ( mesh%vi1:mesh%vi2), source=.false.) !                 Grounded mask on a-grid
    ! ALLOCATE( laddie%mask_oc_a          ( mesh%vi1:mesh%vi2), source=.false.) !                 Icefree ocean mask on a-grid
    ! ALLOCATE( laddie%mask_b             ( mesh%ti1:mesh%ti2), source=.false.) !                 Mask on b-grid
    ! ALLOCATE( laddie%mask_gl_b          ( mesh%ti1:mesh%ti2), source=.false.) !                 Grounding line mask on b-grid
    ! ALLOCATE( laddie%mask_cf_b          ( mesh%ti1:mesh%ti2), source=.false.) !                 Calving front mask on b-grid
    ! ALLOCATE( laddie%mask_oc_b          ( mesh%ti1:mesh%ti2), source=.false.) !                 Icefree ocean mask on b-grid

    ! ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ! ===================================================================

    ! ! Matrix size
    ! ncols           = mesh%nV        ! from
    ! ncols_loc       = mesh%nV_loc
    ! nrows           = mesh%nTri      ! to
    ! nrows_loc       = mesh%nTri_loc
    ! nnz_per_row_est = 3
    ! nnz_est_proc    = nrows_loc * nnz_per_row_est

    ! call allocate_matrix_CSR_dist( laddie%M_map_H_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! ! Matrix size
    ! ncols           = mesh%nV        ! from
    ! ncols_loc       = mesh%nV_loc
    ! nrows           = mesh%nE      ! to
    ! nrows_loc       = mesh%nE_loc
    ! nnz_per_row_est = 2
    ! nnz_est_proc    = nrows_loc * nnz_per_row_est

    ! call allocate_matrix_CSR_dist( laddie%M_map_H_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! ! Matrix size
    ! ncols           = mesh%nTri        ! from
    ! ncols_loc       = mesh%nTri_loc
    ! nrows           = mesh%nE      ! to
    ! nrows_loc       = mesh%nE_loc
    ! nnz_per_row_est = 2
    ! nnz_est_proc    = nrows_loc * nnz_per_row_est

    ! call allocate_matrix_CSR_dist( laddie%M_map_UV_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

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

    call crash('fixme!')

    ! ALLOCATE( npx%H                  ( mesh%vi1:mesh%vi2), source=0._dp) ! [m]             Layer thickness
    ! ALLOCATE( npx%H_b                ( mesh%ti1:mesh%ti2), source=0._dp) ! [m]             Layer thickness on b grid
    ! ALLOCATE( npx%H_c                ( mesh%ei1:mesh%ei2), source=0._dp) ! [m]             Layer thickness on c grid
    ! ALLOCATE( npx%U                  ( mesh%ti1:mesh%ti2), source=0._dp) ! [m s^-1]        2D velocity
    ! ALLOCATE( npx%U_a                ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]        2D velocity on a grid
    ! ALLOCATE( npx%U_c                ( mesh%ei1:mesh%ei2), source=0._dp) ! [m s^-1]        2D velocity on b grid
    ! ALLOCATE( npx%V                  ( mesh%ti1:mesh%ti2), source=0._dp) ! [m s^-1]
    ! ALLOCATE( npx%V_a                ( mesh%vi1:mesh%vi2), source=0._dp) ! [m s^-1]
    ! ALLOCATE( npx%V_c                ( mesh%ei1:mesh%ei2), source=0._dp) ! [m s^-1]
    ! ALLOCATE( npx%T                  ( mesh%vi1:mesh%vi2), source=0._dp) ! [degC]          Temperature
    ! ALLOCATE( npx%S                  ( mesh%vi1:mesh%vi2), source=0._dp) ! [PSU]           Salinity

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

    call crash('fixme!')

    ! call average_over_domain( mesh, laddie%now%H, H_av)

    ! Meltmax = maxval( laddie%melt) * sec_per_year
    ! call MPI_ALLREDUCE( MPI_IN_PLACE, Meltmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Umax = maxval( sqrt( laddie%now%U**2 + laddie%now%V**2))
    ! call MPI_ALLREDUCE( MPI_IN_PLACE, Umax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Tmax = maxval( laddie%now%T)
    ! call MPI_ALLREDUCE( MPI_IN_PLACE, Tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! if (par%primary) then
    !   write( *, "(F8.3,A,F8.3,A,F8.2,A,F8.3,A,F8.3)") tl/sec_per_day, &
    !     '  Dmean ', H_av, '  Meltmax', Meltmax, '   U', Umax, '   Tmax', Tmax
    ! end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine print_diagnostics

END MODULE laddie_utilities

