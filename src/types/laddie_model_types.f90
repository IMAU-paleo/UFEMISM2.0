MODULE laddie_model_types

  ! The different data types used in the laddie modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_laddie_model
    ! The laddie model structure
  
    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: H                           ! [m]               Layer thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: U                           ! [m s^-1]          2D velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: V                           ! [m s^-1]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: S                           ! [PSU]             Salinity  
  
    ! Time domain
    REAL(dp)                                :: dt                          ! [s]               Time step
    REAL(dp)                                :: tend                        ! [s]               Time end of Laddie cycle
  
    ! Time stepping
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: H_prev                      ! [m]               Layer thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: U_prev                      ! [m s^-1]          2D velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: V_prev                      ! [m s^-1]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_prev                      ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: S_prev                      ! [PSU]             Salinity  
  
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: H_next                      ! [m]               Layer thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: U_next                      ! [m s^-1]          2D velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: V_next                      ! [m s^-1]
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_next                      ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: S_next                      ! [PSU]             Salinity   
  
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dH_dt                       ! [m s^-1]          Layer thickness change
  
    ! Ambient fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_amb                       ! [degrees Celsius] Ambient temperature at layer base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: S_amb                       ! [PSU]             Ambient salinity at layer base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: rho_amb                     ! [kg m^-3]         Ambient density at layer base
  
    ! Physical variables
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_freeze                    ! [degrees Celsius] Freezing temperature at ice shelf base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_base                      ! [degrees Celsius] Temperature at ice shelf base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: S_base                      ! [PSU]             Salinity at ice shelf base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: rho                         ! [kg m^-3]         Density of mixed layer water
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: drho_amb                    ! []                Buoyancy at layer bottom (rho_amb-rho)/rho_sw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: drho_base                   ! []                Buoyancy at ice base (rho-rho_base)/rho_sw
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_star                      ! [m s^-1]          Friction velocity
  
    ! Physical parameter fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: gamma_T                     ! []                Turbulent heat exchange coefficient
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: gamma_S                     ! []                Turbulent salt exchange coefficient
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: A_h                         ! [m^2 s^-1]        Horizontal laplacian viscosity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: K_h                         ! [m^2 s^-1]        Horizontal diffusivity
  
    ! Vertical fluxes
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: melt                        ! [m s^-1]          Basal melting or freezing rate
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: entr                        ! [m s^-1]          Entrainment rate of ambient water
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: entr_dmin                   ! [m s^-1]          Additional entrainment to retain D_min
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: detr                        ! [m s^-1]          Detrainment rate into ambient water
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: entr_tot                    ! [m s^-1]          Total (net) entrainment = entr+entr_dmin-detr
  
    ! Horizontal fluxes
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: divQ                        ! [m^3 s^-1]        Divergence of layer thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: divQU                       ! [m^4 s^-2]        Divergence of momentum
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: divQV                       ! [m^4 s^-2]        
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: divQT                       ! [degC m^3 s^-1]   Divergence of heat
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: divQS                       ! [PSU m^3 s^-1]    Divergence of salt
  
    ! Viscosities
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: viscU                       ! [m^2 s^-2]        Horizontal viscosity term
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: viscV                       ! [m^2 s^-2]              
  
    ! Diffusivities
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: diffT                       ! [degC m s^-1]     Horizontal diffusivity of heat 
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: diffS                       ! [PSU m s^-1]      Horizontal diffusivity of salt
  
    ! RHS terms
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: ddrho_amb_dx_b              ! [m^-1]            Horizontal derivative of buoyancy 
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: ddrho_amb_dy_b              ! [m^-1]            Horizontal derivative of buoyancy 
  
  END TYPE type_laddie_model


CONTAINS

END MODULE laddie_model_types
