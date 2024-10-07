MODULE laddie_tracers

  ! Tracer routines for the laddie model

! ===== Preamble =====
! ====================
    
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE laddie_model_types                                     , ONLY: type_laddie_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds

  IMPLICIT NONE
    
CONTAINS
    
! ===== Main routines =====
! =========================

  SUBROUTINE compute_TS_np1( mesh, ice, laddie, dt)
    ! Integrate T and S by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_TS_np1'
    INTEGER                                               :: vi
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_TS_np1

END MODULE laddie_tracers

