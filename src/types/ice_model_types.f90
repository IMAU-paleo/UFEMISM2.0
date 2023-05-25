MODULE ice_model_types

  ! The different data types used in the ice modules

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE mpi_distributed_memory                                 , ONLY: partition_list
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  TYPE type_ice_velocity_solver_SIA
    ! Data fields needed to solve the Shallow Ice Approximation

    ! Solution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D_b                      ! [m yr^-1] 3-D ice velocity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D_a                  ! [yr^-1] Vertical shear strain rates
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D_a

    ! Intermediate data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: D_3D_b                      ! [m yr^-1] Diffusivity

  END TYPE type_ice_velocity_solver_SIA

  TYPE type_ice_velocity_solver_SSA
    ! Data fields needed to solve the Shallow Shelf Approximation

    ! Solution
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_b                         ! [m yr^-1] 2-D horizontal ice velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_b

    ! Intermediate data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: A_flow_vav_a                ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: du_dx_a                     ! [yr^-1] 2-D horizontal strain rates
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: du_dy_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dv_dx_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dv_dy_a
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: eta_a                       ! Effective viscosity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: N_a                         ! Product term N = eta * H
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: N_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dN_dx_b                     ! Gradients of N
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dN_dy_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: beta_b_b                    ! Friction coefficient (tau_b = u * beta_b)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dx_b                    ! Driving stress
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: tau_dy_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_b_prev                    ! Velocity solution from previous viscosity iteration
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_b_prev

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised SSA
    REAL(dp)                                :: PETSc_rtol
    REAL(dp)                                :: PETSc_abstol

    ! Stiffness matrices and load vectors
    TYPE(tMat)                              :: AA                          ! Total stiffness matrix
    TYPE(tMat)                              :: AA_free                     ! Stiffness matrix describing the free SSA
    TYPE(tMat)                              :: AA_BC_west                  ! Stiffness matrix describing boundary conditions on the western  domain boundary
    TYPE(tMat)                              :: AA_BC_east                  ! Stiffness matrix describing boundary conditions on the eastern  domain boundary
    TYPE(tMat)                              :: AA_BC_south                 ! Stiffness matrix describing boundary conditions on the southern domain boundary
    TYPE(tMat)                              :: AA_BC_north                 ! Stiffness matrix describing boundary conditions on the northern domain boundary
    TYPE(tMat)                              :: AA_BC_prescr                ! Stiffness matrix describing boundary conditions for prescribed velocities

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: bb                          ! Total load vector
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: bb_free                     ! Load vector for the free SSA
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: bb_BC_west                  ! Load vector for the boundary conditions on the western  domain boundary
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: bb_BC_east                  ! Load vector for the boundary conditions on the eastern  domain boundary
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: bb_BC_south                 ! Load vector for the boundary conditions on the southern domain boundary
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: bb_BC_north                 ! Load vector for the boundary conditions on the northern domain boundary
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: bb_BC_prescr                ! Load vector for the boundary conditions for prescribed velocities

    ! Matrices describing masks for the different components of the equation
    TYPE(tMat)                              :: m_free                      ! Matrix mask for the free SSA
    TYPE(tMat)                              :: m_BC_west                   ! Matrix mask for the boundary conditions on the western  domain boundary
    TYPE(tMat)                              :: m_BC_east                   ! Matrix mask for the boundary conditions on the eastern  domain boundary
    TYPE(tMat)                              :: m_BC_south                  ! Matrix mask for the boundary conditions on the southern domain boundary
    TYPE(tMat)                              :: m_BC_north                  ! Matrix mask for the boundary conditions on the northern domain boundary
    TYPE(tMat)                              :: m_BC_prescr                 ! Matrix mask for the boundary conditions for prescribed velocities

    ! Some useful combined mesh operators
    TYPE(tMat)                              :: M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b
    TYPE(tMat)                              :: M_4_dudx_p_2_dvdy_buv_b
    TYPE(tMat)                              :: M_dudy_p_dvdx_buv_b
    TYPE(tMat)                              :: M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b
    TYPE(tMat)                              :: M_4_dvdy_p_2_dudx_buv_b
    TYPE(tMat)                              :: M_dvdx_p_dudy_buv_b

  END TYPE type_ice_velocity_solver_SSA

  TYPE type_ice_model
    ! The ice dynamics model data structure.

    ! === Masks ===
    ! =============

    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_land                   ! T: land above water level (sea and/or lake), F: land below water level (idem)
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_ocean                  ! T: land below sea level, F: land above sea level
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_lake                   ! T: land below lake level, F: land above lake level
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_ice                    ! T: Hi > 0, F: Hi = 0
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_sheet                  ! T: grounded ice, F: floating ice or ice-free
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_shelf                  ! T: floating ice, F: grounded ice or ice-free
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_coast                  ! T: land next to sea/lake, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_margin                 ! T: ice next to ice-free, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_gr                  ! T: grounded ice next to floating ice, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_gl_fl                  ! T: floating ice next to grounded ice, F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_cf_gr                  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_cf_fl                  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE :: mask_ice_prev               ! T: ice-covered during previous time step, F: otherwise
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: mask                        ! Combined mask with integers representing different parts, only used for quick visual inspection of output

    ! === Basic geometry ===
    ! ======================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hi                          ! [m] Ice thickness
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hb                          ! [m] Bedrock elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hs                          ! [m] Surface elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: Hib                         ! [m] Ice base elevation (w.r.t. PD sea level)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: TAF                         ! [m] Thickness above flotation

    ! === Geometry changes ===
    ! ========================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHi                         ! [m] Ice thickness difference (w.r.t. to reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHb                         ! [m] Bedrock elevation difference (w.r.t. to reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHs                         ! [m] Surface elevation difference (w.r.t. to reference)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: dHib                        ! [m] Base elevation difference (w.r.t. to reference)

    ! === Ice change ===
    ! ==================

    ! TBA

    ! === Thermodynamics ===
    ! ======================

    ! Ice temperatures
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti                          ! [K] Englacial temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ti_pmp                      ! [K] Pressure melting point temperature

    ! Physical quantities
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Cpi                         ! [J kg^-1 K^-1] Specific heat capacity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Ki                          ! [J m^-1 K^-1 yr^-1] Thermal conductivity

    ! === Ice flow ===
    ! ================

    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: A_flow_3D                   ! [Pa^-3 y^-1] Glen's flow law parameter

    ! === Ice velocities ===
    ! ======================

    ! Velocity solvers
    TYPE(type_ice_velocity_solver_SIA)      :: SIA                         ! Shallow Ice Approximation
    TYPE(type_ice_velocity_solver_SSA)      :: SSA                         ! Shallow Shelf Approximation
    ! TYPE(type_ice_velocity_solver_DIVA)     :: DIVA                        ! Depth-Integrated Viscosity Approximation
    ! TYPE(type_ice_velocity_solver_BPA)      :: BPA                         ! Blatter-Pattyn Approximation

    ! 3-D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D                        ! [m yr^-1] 3-D ice velocity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: u_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: w_3D

    ! Vertically integrated
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav                       ! [m yr^-1] Vertically averaged ice velocity
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_vav_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_vav_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_vav
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_vav_b

    ! Surface
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_surf                      ! [m yr^-1] Ice velocity at the surface
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_surf_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_surf_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: w_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_surf
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_surf_b

    ! Basal
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_base                      ! [m yr^-1] Ice velocity at the base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: u_base_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: v_base_b
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: w_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: uabs_base_b

    ! === Strain rates ===
    ! ====================

    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dx_3D                    ! [yr^-1]
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: du_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dv_dz_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dx_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dy_3D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: dw_dz_3D

    ! == Basal conditions ==
    ! ======================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: beta_b                      ! Basal friction coefficient (tau_b = u * beta_b)

    ! === Sea level ===
    ! =================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SL                          ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)

    ! === Grounded area fractions ===
    ! ===============================

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: f_grnd                      ! [0-1] Grounded fractions of vertices
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: f_grnd_b                    ! [0-1] Grounded fractions of triangles

    ! === Extras ===
    ! ==============

    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: basin_ID                    ! The drainage basin to which each vertex belongs

  END TYPE type_ice_model

CONTAINS

! ===== Subroutines ======
! ========================

END MODULE ice_model_types