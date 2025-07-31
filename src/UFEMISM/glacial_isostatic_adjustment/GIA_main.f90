MODULE GIA_main

  ! The main GIA model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE GIA_model_types                                        , ONLY: type_GIA_model, type_ELRA_model
  USE region_types                                           , ONLY: type_model_region
  USE grid_basic                                             , ONLY: setup_square_grid
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use reference_geometry_types                               , only: type_reference_geometry
  USE GIA_ELRA                                               , only: run_ELRA_model, calculate_ELRA_bedrock_deformation_rate, initialise_ELRA_model, remap_ELRA_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_GIA_model( region)
    ! Calculate bedrock deformation at the desired time, and update
    ! predicted deformation if necessary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),                INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_GIA_model'
    REAL(dp)                                              :: wt_prev, wt_next
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the desired time is beyond the time of the next modelled bedrock deformation,
    ! run the GIA model to calculate a new next modelled bedrock deformation.
    ! ============================================================================

    IF (region%time == region%GIA%t_next) THEN
      ! Need to calculate new predicted bedrock deformation

      ! Store previous modelled bedrock deformation
      region%GIA%dHb_prev = region%GIA%dHb_next
      region%GIA%t_prev = region%GIA%t_next
      region%GIA%t_next = region%GIA%t_prev + C%dt_GIA

      ! Run the GIA model to calculate a new next modelled bedrock deformation
      IF     (C%choice_GIA_model == 'none') THEN
        ! No need to do anything
      ELSEIF (C%choice_GIA_model == 'ELRA') THEN
        CALL run_ELRA_model( region)
      ELSE
        CALL crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"!')
      END IF

    ELSEIF (region%time > region%GIA%t_next) THEN
      ! This should not be possible
      CALL crash('overshot the GIA time step')
    ELSE
      ! We're within the current GIA prediction window
    END IF ! IF (region%time == region%GIA%t_next) THEN

    ! Interpolate between previous and next modelled bedrock deformation
    ! to find the bedrock deformation and elevation at the desired time
    ! =================================================================

    ! Calculate time interpolation weights
    wt_prev = (region%GIA%t_next - region%time) / (region%GIA%t_next - region%GIA%t_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled bedrock deformation to desired time
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%dHb( vi) = wt_prev * region%GIA%dHb_prev( vi) + wt_next * region%GIA%dHb_next( vi)
    END DO

    ! Calculate all other GIA quantities
    ! ==================================

!    DO vi = region%mesh%vi1, region%mesh%vi2
!      region%ice%dHb_dt( vi) = (region%GIA%dHb_next( vi) - region%GIA%dHb_prev( vi)) / C%dt_GIA
!    END DO

    DO vi = region%mesh%vi1, region%mesh%vi2
  	  region%ice%Hb( vi) = region%refgeo_GIAeq%Hb( vi) + region%ice%dHb( vi)
	END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_GIA_model

  SUBROUTINE initialise_GIA_model( mesh, GIA, region_name, refgeo_GIAeq, ELRA)
    ! Initialise the GIA model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_GIA_model),                   INTENT(OUT)   :: GIA
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    TYPE(type_reference_geometry), optional, INTENT(IN)   :: refgeo_GIAeq
    TYPE(type_ELRA_model), optional,         INTENT(OUT)   :: ELRA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_GIA_model'
    CHARACTER(LEN=256)                                    :: grid_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising GIA model...'

    ! Create the square grid for the GIA model
    grid_name = 'square_grid_GIA_' // region_name
    CALL setup_square_grid( grid_name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, &
       C%dx_GIA, GIA%grid, &
       lambda_M = mesh%lambda_M, phi_M = mesh%phi_M, beta_stereo = mesh%beta_stereo)

    ! Allocate memory for main variables

    ! In: relative surface load
    ALLOCATE( GIA%relative_surface_load_mesh( mesh%vi1:mesh%vi2))
    ALLOCATE( GIA%relative_surface_load_grid( GIA%grid%n1:GIA%grid%n2))
    GIA%relative_surface_load_mesh = 0._dp
    GIA%relative_surface_load_grid = 0._dp

    ! Out: modelled bedrock deformation
    ALLOCATE( GIA%dHb_prev( mesh%vi1:mesh%vi2))
    ALLOCATE( GIA%dHb_next( mesh%vi1:mesh%vi2))
    GIA%dHb_prev = 0._dp
    GIA%dHb_next = 0._dp

    ! Model states for GIA model
    GIA%t_prev   = C%start_time_of_run
    GIA%t_next   = C%start_time_of_run
    GIA%dHb_prev = 0._dp
    GIA%dHb_next = 0._dp

    ! Determine which GIA model to initialise
    IF     (C%choice_GIA_model == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL initialise_ELRA_model( mesh, GIA%grid, ELRA, refgeo_GIAeq)
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_GIA_model

  SUBROUTINE write_to_restart_file_GIA_model( mesh, GIA, region_name, time)
    ! Write to the restart file for the GIA model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_GIA_model),                   INTENT(IN)    :: GIA
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_GIA_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write to the restart file of the chosen GIA model
    IF     (C%choice_GIA_model == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_GIA_model

  SUBROUTINE create_restart_file_GIA_model( mesh, GIA, region_name)
    ! Create the restart file for the GIA model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_GIA_model),                   INTENT(INOUT) :: GIA
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_GIA_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create the restart file of the chosen GIA model
    IF     (C%choice_GIA_model == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_GIA_model

  SUBROUTINE remap_GIA_model( mesh_old, mesh_new, GIA, refgeo_GIAeq, ELRA)
    ! Remap the GIA model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_GIA_model),                   INTENT(OUT)   :: GIA
    TYPE(type_reference_geometry), optional, INTENT(IN)    :: refgeo_GIAeq
    TYPE(type_ELRA_model), optional, INTENT(OUT) :: ELRA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_GIA_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping GIA model data to the new mesh...'

    ! Reallocate memory for main variables

    ! In: relative surface load
    CALL reallocate_bounds( GIA%relative_surface_load_mesh, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( GIA%dHb_prev                  , mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( GIA%dHb_next                  , mesh_new%vi1, mesh_new%vi2)

    ! Determine which GIA model to remap
    IF     (C%choice_GIA_model == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL remap_ELRA_model( mesh_old, mesh_new, ELRA, refgeo_GIAeq, GIA%grid)
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_GIA_model

END MODULE GIA_main
