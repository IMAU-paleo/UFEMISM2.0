MODULE laddie_thickness

  ! Thickness routines for the laddie model

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

  SUBROUTINE compute_H_np1( mesh, ice, laddie, dt)
    ! Integrate H by one time step

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_H_np1'
    INTEGER                                               :: vi
    REAL(dp)                                              :: dHdt
 
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
      IF (laddie%mask_a( vi)) THEN

        ! Get first guess at dHdt
        dHdt = -laddie%divQ( vi) + laddie%melt( vi) + laddie%entr( vi)

        ! First guess at H_n
        laddie%H_next( vi) = laddie%H( vi) + dHdt * dt

        ! If H_n < Hmin, enhance entrainment to ensure H_n >= Hmin
        laddie%entr_dmin( vi) = MAX( C%laddie_thickness_minimum - laddie%H_next( vi), 0.0_dp) / dt

        ! If H_n > Hmax, suppress entrainment to ensure H_n <= Hmax
        laddie%entr( vi) = laddie%entr( vi) + MIN( C%laddie_thickness_maximum - laddie%H_next( vi), 0.0_dp) / dt

        ! Update detrainment. Shouldn't matter but just in case
        laddie%detr( vi) = - MIN(laddie%entr( vi),0.0_dp)

        ! Get actual dHdt
        dHdt = -laddie%divQ( vi) + laddie%melt( vi) + laddie%entr( vi) + laddie%entr_dmin( vi)

        ! Get actual H_n
        laddie%H_next( vi) = laddie%H( vi) + dHdt * dt

      END IF !(laddie%mask_a( vi)) THEN
    END DO !vi = mesh%vi, mesh%v2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_H_np1

END MODULE laddie_thickness

