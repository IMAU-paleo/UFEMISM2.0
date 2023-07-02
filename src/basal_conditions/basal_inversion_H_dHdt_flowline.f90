MODULE basal_inversion_H_dHdt_flowline

  ! Contains all the routines for the basal inversion model
  ! based on flowline-averaged values of H and dH/dt

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE basal_inversion_types                                  , ONLY: type_basal_inversion_H_dHdt_flowline

  IMPLICIT NONE

CONTAINS

  ! ===== Main routines =====
  ! =========================

  SUBROUTINE run_basal_inversion_H_dHdt_flowline( mesh, ice, BIV)
    ! Initialise the basal inversion model based on flowline-averaged values of H and dH/dt

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)    :: ice
    TYPE(type_basal_inversion_H_dHdt_flowline), INTENT(INOUT)   :: BIV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_basal_inversion_H_dHdt_flowline'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_basal_inversion_H_dHdt_flowline

  SUBROUTINE initialise_basal_inversion_H_dHdt_flowline( mesh, ice, BIV, region_name)
    ! Initialise the basal inversion model based on flowline-averaged values of H and dH/dt

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_basal_inversion_H_dHdt_flowline), INTENT(OUT)   :: BIV
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion_H_dHdt_flowline'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion_H_dHdt_flowline

END MODULE basal_inversion_H_dHdt_flowline
