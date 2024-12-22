MODULE SMB_main

  ! The main SMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE SMB_idealised                                          , ONLY: initialise_SMB_model_idealised, run_SMB_model_idealised
  USE SMB_prescribed                                         , ONLY: initialise_SMB_model_prescribed, run_SMB_model_prescribed
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  use plane_geometry, only: is_in_polygon
  USE mesh_data_smoothing                                    , ONLY: smooth_Gaussian_2D

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_SMB_model( mesh, grid_smooth, ice, climate, SMB, region_name, time)
    ! Calculate the surface mass balance

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_grid),                        INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(IN)    :: climate
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new SMB
    IF (C%do_asynchronous_SMB) THEN
      ! Asynchronous coupling: do not calculate a new SMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next SMB time step
      IF (time == SMB%t_next) THEN
        ! Go on to calculate a new SMB
        SMB%t_next = time + C%dt_SMB
      ELSEIF (time > SMB%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the SMB time step')
      ELSE
        ! It is not yet time to calculate a new SMB
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_SMB) THEN
      ! Synchronous coupling: calculate a new SMB in every model loop
      SMB%t_next = time + C%dt_SMB
    END IF

    ! Determine which SMB model to run for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Run the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        SMB%SMB = C%uniform_SMB
      CASE ('idealised')
        CALL run_SMB_model_idealised( mesh, ice, SMB, time)
      CASE ('prescribed')
        CALL run_SMB_model_prescribed( mesh, ice, SMB, region_name, time)
      CASE ('reconstructed')
        CALL run_SMB_model_reconstructed( mesh, grid_smooth, ice, SMB, region_name, time)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model

  SUBROUTINE initialise_SMB_model( mesh, SMB, region_name)
    ! Initialise the SMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(OUT)   :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising surface mass balance model...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Allocate memory for main variables
    ALLOCATE( SMB%SMB( mesh%vi1:mesh%vi2))
    SMB%SMB = 0._dp

    ! Set time of next calculation to start time
    SMB%t_next = C%start_time_of_run

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        SMB%SMB = C%uniform_SMB
      CASE ('idealised')
        CALL initialise_SMB_model_idealised( mesh, SMB)
      CASE ('prescribed')
        CALL initialise_SMB_model_prescribed( mesh, SMB, region_name)
      CASE ('reconstructed')
        CALL initialise_SMB_model_reconstructed( mesh, SMB, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model

  SUBROUTINE write_to_restart_file_SMB_model( mesh, SMB, region_name, time)
    ! Write to the restart file for the SMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which SMB model to use for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Write to the restart file of the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE ('reconstructed')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_SMB_model

  SUBROUTINE create_restart_file_SMB_model( mesh, SMB, region_name)
    ! Create the restart file for the SMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which SMB model to use for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Create the restart file of the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE ('reconstructed')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_SMB_model

  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, SMB, region_name)
    ! Remap the SMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_SMB_model),                   INTENT(OUT)   :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '    Remapping surface mass balance model data to the new mesh...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Reallocate memory for main variables
    CALL reallocate_bounds( SMB%SMB, mesh_new%vi1, mesh_new%vi2)

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        CALL initialise_SMB_model_prescribed( mesh_new, SMB, region_name)
      CASE ('reconstructed')
        CALL crash('Remapping after mesh update not implemented yet for reconstructed SMB')
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SMB_model

! ===== Inversion =====
! =====================

  SUBROUTINE run_SMB_model_reconstructed( mesh, grid_smooth, ice, SMB, region_name, time)
    ! Calculate the surface mass balance
    !
    ! Use a reconstructed SMB approach

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_grid),                        INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model_reconstructed'
    INTEGER                                               :: vi
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE               :: poly_ROI             ! Polygon defining reconstructed area
    REAL(dp), DIMENSION(2)                                :: p                    ! Coordinates of a vertex
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: SMB_smoothed         ! Smoothed SMB field
    REAL(dp)                                              :: w_smooth             ! Weight of the smoothed SMB field
    REAL(dp), PARAMETER                                   :: r_smooth =  2.E4_dp  ! Radius used to smooth the SMB field
    REAL(dp), PARAMETER                                   :: Hs_ela   =  500._dp  ! Equilibrium line altitud: SMB becomes positive here
    REAL(dp), PARAMETER                                   :: Hs_tla   =  1500._dp ! Transitional line altitud: SMB reaches maximum here
    REAL(dp), PARAMETER                                   :: Hs_dla   =  2500._dp ! Desertification line altitude: SMB becomes zero here
    REAL(dp), PARAMETER                                   :: SMB_max  =  2._dp    ! Maximum SMB value allowed
    REAL(dp), PARAMETER                                   :: SMB_min  = -10._dp   ! Minimum SMB value allowed

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%choice_regions_of_interest == 'Patagonia') THEN
      CALL crash('reconstructed SMB method only implemented for C%choice_regions_of_interest == Patagonia')
    END IF

    ! Compute polygon for reconstruction
    CALL calc_polygon_Patagonia( poly_ROI)

    DO vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Check if point lies within our reconstruction polygon
      IF (is_in_polygon(poly_ROI, p)) THEN
        ! If yes, check whether point lies above or below estimated transitional line altitude
        IF (ice%Hs( vi) <= Hs_tla) THEN
          ! If below, SMB goes from 0 at the ELA to its estimated maximum at the TLA
          SMB%SMB( vi) = SMB_max * MAX( 0._dp, MIN( 1._dp, (ice%Hs( vi) - Hs_ela)/(Hs_tla - Hs_ela)))
        ELSE
          ! If above, SMB goes from estimated maximum at the TLA to 0 at the DLA
          SMB%SMB( vi) = SMB_max * (1._dp - MAX( 0._dp, MIN( 1._dp, (ice%Hs( vi) - Hs_tla)/(Hs_dla - Hs_tla))))
        END IF
      ELSE
        ! If vertex lies outside of the reconstructed polygon, assume a negative
        ! SMB that counters (and then some) the flux convergence there
        SMB%SMB( vi) = MIN( 0._dp, MAX( SMB_min, ice%divQ( vi) - .5_dp))
      END IF

    END DO

    ! Smooth the reconstructed field
    SMB_smoothed = SMB%SMB
    CALL smooth_Gaussian_2D( mesh, grid_smooth, SMB_smoothed, r_smooth)

    ! Only apply the smoothed field inside the reconstructed area
    ! to reduce the power of positive SMB there
    DO vi = mesh%vi1, mesh%vi2
      ! Our vextex coordinates
      p = mesh%V( vi,:)
      ! Check if point lies inside polygon
      IF (is_in_polygon(poly_ROI, p)) THEN
        ! Compute a weight based on Hs: the higher, the less smoothing
        w_smooth = MAX( 0._dp, MIN( 1._dp, ice%Hs( vi) / Hs_dla))
        ! Apply weighed smoothing
        SMB%SMB( vi) = w_smooth * SMB%SMB( vi) + (1._dp - w_smooth) * SMB_smoothed( vi)
      END IF
    END DO

    ! Smooth the field once more
    SMB_smoothed = SMB%SMB
    CALL smooth_Gaussian_2D( mesh, grid_smooth, SMB_smoothed, r_smooth)

    ! Apply this second smoothing only outside of the reconstructed area
    ! to conserve the power of negative SMB there
    DO vi = mesh%vi1, mesh%vi2
      p = mesh%V( vi,:)
      IF (.NOT. is_in_polygon(poly_ROI, p)) SMB%SMB( vi) = SMB_smoothed( vi)
    END DO

    ! Clean up after yourself
    DEALLOCATE( poly_ROI)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_reconstructed

  SUBROUTINE initialise_SMB_model_reconstructed( mesh, SMB, region_name)
    ! Initialise the SMB model
    !
    ! Use a reconstructed SMB approach

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model_reconstructed'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising reconstructed SMB model...'

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model_reconstructed

END MODULE SMB_main
