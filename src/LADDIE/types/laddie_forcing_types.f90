module laddie_forcing_types

  ! The different data types used in the laddie modules

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use CSR_sparse_matrix_type                                 , only: type_sparse_matrix_CSR_dp
  use mpi_f08, only: MPI_WIN

  implicit none

! ===== Types =====
! =================

  type type_laddie_forcing
    ! The laddie model structure

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

  end type type_laddie_forcing

contains

end module laddie_forcing_types
