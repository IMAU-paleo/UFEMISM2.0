MODULE laddie_model_types

  ! The different data types used in the laddie modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE CSR_sparse_matrix_type                                 , ONLY: type_sparse_matrix_CSR_dp
  use mpi_f08, only: MPI_WIN

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_laddie_timestep
    ! Fields of (partial) timesteps

    ! Main data fields
    real(dp), dimension(:), contiguous, pointer :: H               => null()  ! [m]               Layer thickness
    real(dp), dimension(:), contiguous, pointer :: H_b             => null()  ! [m]               Layer thickness on b grid
    real(dp), dimension(:), contiguous, pointer :: H_c             => null()  ! [m]               Layer thickness on c grid
    real(dp), dimension(:), contiguous, pointer :: U               => null()  ! [m s^-1]          2D velocity
    real(dp), dimension(:), contiguous, pointer :: U_a             => null()  ! [m s^-1]          2D velocity on a grid
    real(dp), dimension(:), contiguous, pointer :: U_c             => null()  ! [m s^-1]          2D velocity on c grid
    real(dp), dimension(:), contiguous, pointer :: V               => null()  ! [m s^-1]
    real(dp), dimension(:), contiguous, pointer :: V_a             => null()  ! [m s^-1]
    real(dp), dimension(:), contiguous, pointer :: V_c             => null()  ! [m s^-1]
    real(dp), dimension(:), contiguous, pointer :: T               => null()  ! [degrees Celsius] Temperature
    real(dp), dimension(:), contiguous, pointer :: S               => null()  ! [PSU]             Salinity
    type(MPI_WIN) :: wH, wH_b, wH_c, wU, wU_a, wU_c, wV, wV_a, wV_c, wT, wS

  END TYPE type_laddie_timestep

  type type_scalar_output_buffer
    !< Memory for buffering scalar output (from every model time step) between output writing intervals

    integer :: n_mem         !< Number of timeframes for which memory has been allocated
    integer :: n             !< Number of timeframes that are currently buffered

    real(dp), dimension(:), allocatable :: time

    real(dp), dimension(:), allocatable :: layer_volume
    real(dp), dimension(:), allocatable :: area_a
    real(dp), dimension(:), allocatable :: area_b

    real(dp), dimension(:), allocatable :: thickness_mean
    real(dp), dimension(:), allocatable :: thickness_min
    real(dp), dimension(:), allocatable :: thickness_max

    real(dp), dimension(:), allocatable :: melt_mean
    real(dp), dimension(:), allocatable :: melt_max
    real(dp), dimension(:), allocatable :: melt_min
    real(dp), dimension(:), allocatable :: melt_tot

    real(dp), dimension(:), allocatable :: uabs_max

    real(dp), dimension(:), allocatable :: T_mean
    real(dp), dimension(:), allocatable :: T_max
    real(dp), dimension(:), allocatable :: T_min

    real(dp), dimension(:), allocatable :: S_mean
    real(dp), dimension(:), allocatable :: S_max
    real(dp), dimension(:), allocatable :: S_min

    real(dp), dimension(:), allocatable :: entr_tot
    real(dp), dimension(:), allocatable :: entr_dmin_tot
    real(dp), dimension(:), allocatable :: detr_tot
    real(dp), dimension(:), allocatable :: divQH_sum

  end type type_scalar_output_buffer

  TYPE type_laddie_model
    ! The laddie model structure

    ! Output
    character(len=1024)                         :: output_fields_filename
    character(len=1024)                         :: output_scalar_filename
    logical                                     :: output_fields_file_matches_current_mesh

    ! Time domain
    real(dp)                                    :: dt                          ! [s]               Time step
    real(dp)                                    :: tend                        ! [s]               Time end of Laddie cycle

    real(dp), dimension(:), contiguous, pointer :: dH_dt           => null()  ! [m s^-1]          Layer thickness change
    type(MPI_WIN) :: wdH_dt

    ! Forcing
    real(dp), dimension(:),   contiguous, pointer :: Hi                 => null()  ! [m]               Ice thickness
    real(dp), dimension(:),   contiguous, pointer :: Hib                => null()  ! [m]               Ice base elevation (w.r.t. PD sea level)
    real(dp), dimension(:),   contiguous, pointer :: dHib_dx_b          => null()  ! []                Horizontal derivative of ice draft on b-grid
    real(dp), dimension(:),   contiguous, pointer :: dHib_dy_b          => null()  ! []                Horizontal derivative of ice draft on b-grid
    logical,  dimension(:),   contiguous, pointer :: mask_icefree_land  => null()  ! []                T: ice-free land , F: otherwise
    logical,  dimension(:),   contiguous, pointer :: mask_icefree_ocean => null()  ! []                T: ice-free ocean, F: otherwise
    logical,  dimension(:),   contiguous, pointer :: mask_grounded_ice  => null()  ! []                T: grounded ice  , F: otherwise
    logical,  dimension(:),   contiguous, pointer :: mask_floating_ice  => null()  ! []                T: floating ice  , F: otherwise
    logical,  dimension(:),   contiguous, pointer :: mask_gl_fl         => null()  ! []                T: gl_fl ice     , F: otherwise
    logical,  dimension(:),   contiguous, pointer :: mask_SGD           => null()  ! []                T: potential subglacial discharge areas, F: otherwise
    real(dp), dimension(:,:), contiguous, pointer :: Ti                 => null()  ! [K]               Englacial temperature
    real(dp), dimension(:,:), contiguous, pointer :: T_ocean            => null()  ! [degrees Celsius] 3-D ocean temperature
    real(dp), dimension(:,:), contiguous, pointer :: S_ocean            => null()  ! [PSU]             3-D ocean salinity
    type(MPI_WIN) :: wHi, wHib, wdHib_dx_b, wdHib_dy_b
    type(MPI_WIN) :: wmask_icefree_land, wmask_icefree_ocean, wmask_grounded_ice, wmask_floating_ice, wmask_gl_fl, wmask_SGD
    type(MPI_WIN) :: wTi, wT_ocean, wS_ocean

    ! Ambient fields
    real(dp), dimension(:), contiguous, pointer :: T_amb           => null()  ! [degrees Celsius] Ambient temperature at layer base
    real(dp), dimension(:), contiguous, pointer :: S_amb           => null()  ! [PSU]             Ambient salinity at layer base
    real(dp), dimension(:), contiguous, pointer :: rho_amb         => null()  ! [kg m^-3]         Ambient density at layer base
    type(MPI_WIN) :: wT_amb, wS_amb, wrho_amb

    ! Physical variables
    real(dp), dimension(:), contiguous, pointer :: T_freeze        => null()  ! [degrees Celsius] Freezing temperature at ice shelf base
    real(dp), dimension(:), contiguous, pointer :: T_base          => null()  ! [degrees Celsius] Temperature at ice shelf base
    real(dp), dimension(:), contiguous, pointer :: S_base          => null()  ! [PSU]             Salinity at ice shelf base
    real(dp), dimension(:), contiguous, pointer :: rho             => null()  ! [kg m^-3]         Density of mixed layer water
    real(dp), dimension(:), contiguous, pointer :: drho_amb        => null()  ! []                Buoyancy at layer bottom (rho_amb-rho)/rho_sw
    real(dp), dimension(:), contiguous, pointer :: Hdrho_amb       => null()  ! [m]               Depth-integrated buoyancy
    real(dp), dimension(:), contiguous, pointer :: Hdrho_amb_b     => null()  ! [m]               Depth-integrated buoyancy
    real(dp), dimension(:), contiguous, pointer :: drho_base       => null()  ! []                Buoyancy at ice base (rho-rho_base)/rho_sw
    real(dp), dimension(:), contiguous, pointer :: u_star          => null()  ! [m s^-1]          Friction velocity
    type(MPI_WIN) :: wT_freeze, wT_base, wS_base, wrho
    type(MPI_WIN) :: wdrho_amb, wHdrho_amb, wHdrho_amb_b, wdrho_base, wu_star

    ! Physical parameter fields
    real(dp), dimension(:), contiguous, pointer :: gamma_T         => null()  ! []                Turbulent heat exchange coefficient
    real(dp), dimension(:), contiguous, pointer :: gamma_S         => null()  ! []                Turbulent salt exchange coefficient
    real(dp), dimension(:), contiguous, pointer :: A_h             => null()  ! [m^2 s^-1]        Horizontal laplacian viscosity
    real(dp), dimension(:), contiguous, pointer :: K_h             => null()  ! [m^2 s^-1]        Horizontal diffusivity
    type(MPI_WIN) :: wgamma_T, wgamma_S, wA_h, wK_h

    ! Vertical fluxes
    real(dp), dimension(:), contiguous, pointer :: melt            => null()  ! [m s^-1]          Basal melting or freezing rate
    real(dp), dimension(:), contiguous, pointer :: entr            => null()  ! [m s^-1]          Entrainment rate of ambient water
    real(dp), dimension(:), contiguous, pointer :: entr_dmin       => null()  ! [m s^-1]          Additional entrainment to retain D_min
    real(dp), dimension(:), contiguous, pointer :: detr            => null()  ! [m s^-1]          Detrainment rate into ambient water
    real(dp), dimension(:), contiguous, pointer :: entr_tot        => null()  ! [m s^-1]          Total (net) entrainment = entr+entr_dmin-detr
    real(dp), dimension(:), contiguous, pointer :: SGD             => null()  ! [m s^-1]          Subglacial discharge rate
    type(MPI_WIN) :: wmelt, wentr, wentr_dmin, wdetr, wentr_tot, wSGD

    ! Horizontal fluxes
    real(dp), dimension(:), contiguous, pointer :: divQH           => null()  ! [m^3 s^-1]        Divergence of layer thickness
    real(dp), dimension(:), contiguous, pointer :: divQU           => null()  ! [m^4 s^-2]        Divergence of momentum
    real(dp), dimension(:), contiguous, pointer :: divQV           => null()  ! [m^4 s^-2]
    real(dp), dimension(:), contiguous, pointer :: divQT           => null()  ! [degC m^3 s^-1]   Divergence of heat
    real(dp), dimension(:), contiguous, pointer :: divQS           => null()  ! [PSU m^3 s^-1]    Divergence of salt
    type(MPI_WIN) :: wdivQH, wdivQU, wdivQV, wdivQT, wdivQS

    ! Viscosities
    real(dp), dimension(:), contiguous, pointer :: viscU           => null()  ! [m^2 s^-2]        Horizontal viscosity term
    real(dp), dimension(:), contiguous, pointer :: viscV           => null()  ! [m^2 s^-2]
    type(MPI_WIN) :: wviscU, wviscV

    ! Diffusivities
    real(dp), dimension(:), contiguous, pointer :: diffT           => null()  ! [degC m s^-1]     Horizontal diffusivity of heat
    real(dp), dimension(:), contiguous, pointer :: diffS           => null()  ! [PSU m s^-1]      Horizontal diffusivity of salt
    type(MPI_WIN) :: wdiffT, wdiffS

    ! RHS terms
    real(dp), dimension(:), contiguous, pointer :: ddrho_amb_dx_b  => null()  ! [m^-1]            Horizontal derivative of buoyancy
    real(dp), dimension(:), contiguous, pointer :: ddrho_amb_dy_b  => null()  ! [m^-1]            Horizontal derivative of buoyancy
    real(dp), dimension(:), contiguous, pointer :: dH_dx_b         => null()  ! [m^-2]            Horizontal derivative of thickness
    real(dp), dimension(:), contiguous, pointer :: dH_dy_b         => null()  ! [m^-2]            Horizontal derivative of thickness
    real(dp), dimension(:), contiguous, pointer :: detr_b          => null()  ! [m s^-1]          Detrainment on b grid
    type(MPI_WIN) :: wddrho_amb_dx_b, wddrho_amb_dy_b, wdH_dx_b, wdH_dy_b, wdetr_b

    ! Forward-Backward Runge-Kutta 3 scheme
    real(dp), dimension(:), contiguous, pointer :: Hstar           => null()  ! [m]               Intermediate layer thickness
    type(MPI_WIN) :: wHstar

    ! Mapped variables
    real(dp), dimension(:), contiguous, pointer :: H_c             => null()  ! [m]               Layer thickness on c grid
    real(dp), dimension(:), contiguous, pointer :: Hstar_b         => null()
    real(dp), dimension(:), contiguous, pointer :: Hstar_c         => null()
    type(MPI_WIN) :: wH_c, wHstar_b, wHstar_c

    ! Masks
    logical,  dimension(:), contiguous, pointer :: mask_a          => null()  !                   Mask on a-grid on which to apply computation
    logical,  dimension(:), contiguous, pointer :: mask_gr_a       => null()  !                   Grounded mask on a-grid
    logical,  dimension(:), contiguous, pointer :: mask_oc_a       => null()  !                   Icefree ocean mask on a-grid
    logical,  dimension(:), contiguous, pointer :: mask_b          => null()  !                   Mask on b-grid on which to apply computation
    logical,  dimension(:), contiguous, pointer :: mask_gl_b       => null()  !                   Grounding line mask on b-grid
    logical,  dimension(:), contiguous, pointer :: mask_cf_b       => null()  !                   Calving front mask on b-grid
    logical,  dimension(:), contiguous, pointer :: mask_oc_b       => null()  !                   Icefree ocean mask on b-grid
    type(MPI_WIN) :: wmask_a, wmask_gr_a, wmask_oc_a, wmask_b, wmask_gl_b, wmask_cf_b, wmask_oc_b

    ! Domains and areas
    real(dp), dimension(:), contiguous, pointer :: domain_a        => null()  ! []                Floating domain on a grid
    real(dp), dimension(:), contiguous, pointer :: domain_b        => null()  ! []                Floating domain on b grid
    type(MPI_WIN) :: wdomain_a, wdomain_b
    real(dp)                                    :: area_a                     ! [m^2]             Integrated area on a grid
    real(dp)                                    :: area_b                     ! [m^2]             Integrated area on b grid

    ! Mapping operators
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_H_a_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_H_a_c
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_UV_b_c

    ! Timestepping types
    TYPE(type_laddie_timestep)              :: nm1                         !                   Timestep n minus 1
    TYPE(type_laddie_timestep)              :: now                         !                   Timestep now
    TYPE(type_laddie_timestep)              :: np1                         !                   Timestep n plus 1
    TYPE(type_laddie_timestep)              :: np12                        !                   Timestep n plus 1/2
    TYPE(type_laddie_timestep)              :: np13                        !                   Timestep n plus 1/3

    ! Scalar output buffer
    type(type_scalar_output_buffer)         :: buffer

  END TYPE type_laddie_model

CONTAINS

END MODULE laddie_model_types
