MODULE UFEMISM_main_model

  ! The main regional ice-sheet model

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model

  IMPLICIT NONE

  TYPE type_model_region
    ! Gathers the different model components for a single model region

    ! Metadata
    CHARACTER(LEN=3)                        :: name                                           ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                      :: long_name                                      ! North America, Eurasia, Greenland, Antarctica

    ! The current time of this particular region.
    REAL(dp)                                :: time

    ! The mesh that all model components define their data on
    TYPE(type_mesh)                         :: mesh

    ! The ice dynamics model
    TYPE(type_ice_model)                    :: ice

  END TYPE type_model_region

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_model( region, t_end)
    ! Integrate this model region forward in time until t_end

    IMPLICIT NONE

    ! In/output variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    TYPE(type_model_region)                            , INTENT(INOUT) :: region
    REAL(dp)                                           , INTENT(IN)    :: t_end    ! [yr]

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! PIEP

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model


END MODULE UFEMISM_main_model