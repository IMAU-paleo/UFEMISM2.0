MODULE AMB_main

  ! The main AMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE AMB_model_types                                        , ONLY: type_AMB_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE initialise_AMB_model( mesh, AMB)
    ! Initialise the AMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_AMB_model),                   INTENT(OUT)   :: AMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_AMB_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising artificial mass balance model...'

    ! Allocate memory for main variables
    ALLOCATE( AMB%AMB( mesh%vi1:mesh%vi2))
    AMB%AMB = 0._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_AMB_model

  SUBROUTINE remap_AMB_model( mesh_old, mesh_new, AMB)
    ! Remap the AMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_AMB_model),                   INTENT(OUT)   :: AMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_AMB_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping artificial mass balance model data to the new mesh...'

    ! Reallocate memory for main variables
    CALL reallocate_bounds( AMB%AMB, mesh_new%vi1, mesh_new%vi2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_AMB_model

END MODULE AMB_main
