MODULE ocean_utilities

  ! Realistic ocean models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model

  IMPLICIT NONE

CONTAINS

! ===== Regridding ======
! =======================

  subroutine initialise_ocean_vertical_grid
    ! Set up the vertical grid used for ocean data - regular grid

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_ocean_vertical_grid'
    INTEGER                       :: k

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the number of vertical layers to be used
    C%nz_ocean = 1 + floor( C%ocean_vertical_grid_max_depth / C%ocean_vertical_grid_dz)

    ! Allocate memory
    allocate( C%z_ocean( C%nz_ocean))

    ! Fill in the values
    do k = 1, C%nz_ocean
      C%z_ocean( k) = real(k-1,dp) * C%ocean_vertical_grid_dz
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_vertical_grid

END MODULE ocean_utilities
