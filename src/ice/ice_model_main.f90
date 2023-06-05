MODULE ice_model_main

  ! The main ice-dynamical model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE scalar_types                                           , ONLY: type_regional_scalars
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ice_model_memory                                       , ONLY: allocate_ice_model
  USE ice_model_utilities                                    , ONLY: determine_masks, calc_bedrock_CDFs, calc_grounded_fractions
  USE ice_model_scalars                                      , ONLY: calc_ice_model_scalars
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE math_utilities                                         , ONLY: ice_surface_elevation, thickness_above_floatation
  USE basal_conditions_main                                  , ONLY: initialise_basal_conditions
  USE ice_velocity_main                                      , ONLY: initialise_velocity_solver

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  subroutine initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ! Initialise all data fields of the ice module

    implicit none

    ! In- and output variables
    type(type_mesh),               intent(in)    :: mesh
    type(type_ice_model),          intent(inout) :: ice
    type(type_reference_geometry), intent(in)    :: refgeo_init
    type(type_reference_geometry), intent(in)    :: refgeo_PD
    type(type_regional_scalars),   intent(out)   :: scalars
    character(len=3),              intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'initialise_ice_model'
    integer                                      :: vi

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write(*,"(A)") '  Initialising ice dynamics model...'
    end if
    call sync

    ! === Memory allocation ===
    ! =========================

    ! Allocate memory
    call allocate_ice_model( mesh, ice)

    ! === Value initialisation ===
    ! ============================

    ! Sea level
    ! =========

    select case (C%choice_sealevel_model)

      case ('fixed')
        ! Fixed sea level
        ice%SL = C%fixed_sealevel

      case ('prescribed')
        ! Sea-level prescribed from external record file
        call crash('Sea level initialisation: prescribed method not implement yet!')
        ! ice%SL = forcing%sealevel_obs

      case ('eustatic')
        ! Eustatic sea level
        call crash('Sea level initialisation: eustatic method not implement yet!')
        ! ice%SL = C%initial_guess_sealevel

      case ('SELEN')
        ! Sea level from SELEN
        call crash('Sea level initialisation: SELEN method not implement yet!')
        ! ice%SL = C%initial_guess_sealevel

      case default
        ! Unknown case
        call crash('unknown choice_sealevel_model "' // &
                    TRIM( C%choice_sealevel_model) // '"!')

    end select

    ! Initial topography
    ! ==================

    ! Basic topography
    do vi = mesh%vi1, mesh%vi2

      ! Main quantities
      ice%Hi ( vi) = refgeo_init%Hi( vi)
      ice%Hb ( vi) = refgeo_init%Hb( vi)
      ice%Hs ( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))
      ice%Hib( vi) = ice%Hs( vi) - ice%Hi( vi)

      ice%TAF( vi)  = thickness_above_floatation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))

      ! Differences w.r.t. present-day
      ice%dHi ( vi)  = ice%Hi ( vi) - refgeo_PD%Hi ( vi)
      ice%dHb ( vi)  = ice%Hb ( vi) - refgeo_PD%Hb ( vi)
      ice%dHs ( vi)  = ice%Hs ( vi) - refgeo_PD%Hs ( vi)
      ice%dHib( vi)  = ice%Hib( vi) - (refgeo_PD%Hs ( vi) - refgeo_PD%Hi( vi))

    end do

    ! Initialised predicted ice thickness
    ice%Hi_tplusdt = ice%Hi

    ! Initial masks
    ! =============

    ! Sector masks
    call determine_masks( mesh, ice)

    ! Initialise previous-time-step mask
    ice%mask_ice_prev = ice%mask_ice

    ! Initial rates of change
    ! =======================

    ice%dHi_dt  = 0._dp
    ice%dHb_dt  = 0._dp
    ice%dHs_dt  = 0._dp
    ice%dHib_dt = 0._dp

    ! Sub-grid fractions
    ! ==================

    ! Compute bedrock cumulative density function
    CALL calc_bedrock_CDFs( mesh, refgeo_PD, ice)
    ! Initialise sub-grid grounded-area fractions
    CALL calc_grounded_fractions( mesh, ice)

    ! Basal conditions
    ! ================

    ! Allocate and initialise basal conditions
    call initialise_basal_conditions( mesh, ice)

    ! Velocities
    ! ==========

    ! Initialise data and matrices for the velocity solver(s)
    call initialise_velocity_solver( mesh, ice, region_name)

    ! Ice-sheet-wide scalars
    ! ======================

    CALL calc_ice_model_scalars( mesh, ice, refgeo_PD, scalars)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ice_model

END MODULE ice_model_main
