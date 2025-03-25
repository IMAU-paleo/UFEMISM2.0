module laddie_model_types

  ! The different data types used in the laddie modules

  use precisions, only: dp
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_laddie_timestep, type_laddie_model

  type type_laddie_timestep
    ! Fields of (partial) timesteps

    ! Main data fields
    real(dp), dimension(:), pointer :: H           ! [m]               Layer thickness
    real(dp), dimension(:), pointer :: H_b         ! [m]               Layer thickness on b grid
    real(dp), dimension(:), pointer :: H_c         ! [m]               Layer thickness on c grid
    real(dp), dimension(:), pointer :: U           ! [m s^-1]          2D velocity
    real(dp), dimension(:), pointer :: U_a         ! [m s^-1]          2D velocity on a grid
    real(dp), dimension(:), pointer :: U_c         ! [m s^-1]          2D velocity on c grid
    real(dp), dimension(:), pointer :: V           ! [m s^-1]
    real(dp), dimension(:), pointer :: V_a         ! [m s^-1]
    real(dp), dimension(:), pointer :: V_c         ! [m s^-1]
    real(dp), dimension(:), pointer :: T           ! [degrees Celsius] Temperature
    real(dp), dimension(:), pointer :: S           ! [PSU]             Salinity
    type(MPI_WIN) :: wH,wH_b,wH_c,wU,wU_a,wU_c,wV,wV_a,wV_c,wT,wS

  end type type_laddie_timestep

  type type_laddie_model
    ! The laddie model structure

    ! Time domain
    real(dp)                        :: dt                          ! [s]               Time step
    real(dp)                        :: tend                        ! [s]               Time end of Laddie cycle

    real(dp), dimension(:), pointer :: dH_dt                       ! [m s^-1]          Layer thickness change
    type(MPI_WIN) :: wdH_dt

    ! Ambient fields
    real(dp), dimension(:), pointer :: T_amb                       ! [degrees Celsius] Ambient temperature at layer base
    real(dp), dimension(:), pointer :: S_amb                       ! [PSU]             Ambient salinity at layer base
    real(dp), dimension(:), pointer :: rho_amb                     ! [kg m^-3]         Ambient density at layer base
    type(MPI_WIN) :: wT_amb, wS_amb, wrho_amb

    ! Physical variables
    real(dp), dimension(:), pointer :: T_freeze                    ! [degrees Celsius] Freezing temperature at ice shelf base
    real(dp), dimension(:), pointer :: T_base                      ! [degrees Celsius] Temperature at ice shelf base
    real(dp), dimension(:), pointer :: S_base                      ! [PSU]             Salinity at ice shelf base
    real(dp), dimension(:), pointer :: rho                         ! [kg m^-3]         Density of mixed layer water
    real(dp), dimension(:), pointer :: drho_amb                    ! []                Buoyancy at layer bottom (rho_amb-rho)/rho_sw
    real(dp), dimension(:), pointer :: Hdrho_amb                   ! [m]               Depth-integrated buoyancy
    real(dp), dimension(:), pointer :: Hdrho_amb_b                 ! [m]               Depth-integrated buoyancy
    real(dp), dimension(:), pointer :: drho_base                   ! []                Buoyancy at ice base (rho-rho_base)/rho_sw
    real(dp), dimension(:), pointer :: u_star                      ! [m s^-1]          Friction velocity
    type(MPI_WIN) :: wT_freeze, wT_base, wS_base, wrho, wdrho_amb, wHdrho_amb, wHdrho_amb_b, wdrho_base, wu_star

    ! Physical parameter fields
    real(dp), dimension(:), pointer :: gamma_T                     ! []                Turbulent heat exchange coefficient
    real(dp), dimension(:), pointer :: gamma_S                     ! []                Turbulent salt exchange coefficient
    real(dp), dimension(:), pointer :: A_h                         ! [m^2 s^-1]        Horizontal laplacian viscosity
    real(dp), dimension(:), pointer :: K_h                         ! [m^2 s^-1]        Horizontal diffusivity
    type(MPI_WIN) :: wgamma_T, wgamma_S, wA_h, wK_h

    ! Vertical fluxes
    real(dp), dimension(:), pointer :: melt                        ! [m s^-1]          Basal melting or freezing rate
    real(dp), dimension(:), pointer :: entr                        ! [m s^-1]          Entrainment rate of ambient water
    real(dp), dimension(:), pointer :: entr_dmin                   ! [m s^-1]          Additional entrainment to retain D_min
    real(dp), dimension(:), pointer :: detr                        ! [m s^-1]          Detrainment rate into ambient water
    real(dp), dimension(:), pointer :: entr_tot                    ! [m s^-1]          Total (net) entrainment = entr+entr_dmin-detr
    type(MPI_WIN) :: wmelt, wentr, wentr_dmin, wdetr, wentr_tot

    ! Horizontal fluxes
    real(dp), dimension(:), pointer :: divQH                       ! [m^3 s^-1]        Divergence of layer thickness
    real(dp), dimension(:), pointer :: divQU                       ! [m^4 s^-2]        Divergence of momentum
    real(dp), dimension(:), pointer :: divQV                       ! [m^4 s^-2]
    real(dp), dimension(:), pointer :: divQT                       ! [degC m^3 s^-1]   Divergence of heat
    real(dp), dimension(:), pointer :: divQS                       ! [PSU m^3 s^-1]    Divergence of salt
    type(MPI_WIN) :: wdivQH, wdivQU, wdivQV, wdivQT, wdivQS

    ! Viscosities
    real(dp), dimension(:), pointer :: viscU                       ! [m^2 s^-2]        Horizontal viscosity term
    real(dp), dimension(:), pointer :: viscV                       ! [m^2 s^-2]
    type(MPI_WIN) :: wviscU, wviscV

    ! Diffusivities
    real(dp), dimension(:), pointer :: diffT                       ! [degC m s^-1]     Horizontal diffusivity of heat
    real(dp), dimension(:), pointer :: diffS                       ! [PSU m s^-1]      Horizontal diffusivity of salt
    type(MPI_WIN) :: wdiffT, wdiffS

    ! RHS terms
    real(dp), dimension(:), pointer :: ddrho_amb_dx_b              ! [m^-1]            Horizontal derivative of buoyancy
    real(dp), dimension(:), pointer :: ddrho_amb_dy_b              ! [m^-1]            Horizontal derivative of buoyancy
    real(dp), dimension(:), pointer :: dH_dx_b                     ! [m^-2]            Horizontal derivative of thickness
    real(dp), dimension(:), pointer :: dH_dy_b                     ! [m^-2]            Horizontal derivative of thickness
    real(dp), dimension(:), pointer :: detr_b                      ! [m s^-1]          Detrainment on b grid
    type(MPI_WIN) :: wddrho_amb_dx_b, wddrho_amb_dy_b, wdH_dx_b, wdH_dy_b, wdetr_b

    ! Mapped variables
    real(dp), dimension(:), pointer :: H_c                         ! [m]               Layer thickness on c grid
    type(MPI_WIN) :: wH_c

    ! Masks
    logical,  dimension(:), pointer :: mask_a                      !                   Mask on a-grid on which to apply computation
    logical,  dimension(:), pointer :: mask_gr_a                   !                   Grounded mask on a-grid
    logical,  dimension(:), pointer :: mask_oc_a                   !                   Icefree ocean mask on a-grid
    logical,  dimension(:), pointer :: mask_b                      !                   Mask on b-grid on which to apply computation
    logical,  dimension(:), pointer :: mask_gl_b                   !                   Grounding line mask on b-grid
    logical,  dimension(:), pointer :: mask_cf_b                   !                   Calving front mask on b-grid
    logical,  dimension(:), pointer :: mask_oc_b                   !                   Icefree ocean mask on b-grid
    type(MPI_WIN) :: wmask_a, wmask_gr_a, wmask_oc_a, wmask_b, wmask_gl_b, wmask_cf_b, wmask_oc_b

    ! Mapping operators
    type(type_sparse_matrix_CSR_dp) :: M_map_H_a_b
    type(type_sparse_matrix_CSR_dp) :: M_map_H_a_c
    type(type_sparse_matrix_CSR_dp) :: M_map_UV_b_c

    ! Timestepping types
    type(type_laddie_timestep)      :: now                         !                   Timestep now
    type(type_laddie_timestep)      :: np1                         !                   Timestep n plus 1
    type(type_laddie_timestep)      :: np12                        !                   Timestep n plus 1/2
    type(type_laddie_timestep)      :: np13                        !                   Timestep n plus 1/3

  end type type_laddie_model

end module laddie_model_types
