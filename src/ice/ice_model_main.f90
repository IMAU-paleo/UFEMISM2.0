MODULE ice_model_main

  ! The main ice-dynamical model module.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ice_model_memory                                       , ONLY: allocate_ice_model
  USE ice_model_utilities                                    , ONLY: determine_masks
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE math_utilities                                         , ONLY: ice_surface_elevation, thickness_above_floatation

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

subroutine initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD)
    ! Initialise all data fields of the ice module

    implicit none

    ! In- and output variables
    type(type_mesh),               intent(in)    :: mesh
    type(type_ice_model),          intent(inout) :: ice
    type(type_reference_geometry), intent(in)    :: refgeo_init
    type(type_reference_geometry), intent(in)    :: refgeo_PD

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

    do vi = 1, mesh%nV_loc
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

    ! Initial masks
    ! =============

    call determine_masks( mesh, ice)

    ! ! Initialise some numbers for the predictor/corrector ice thickness update method
    ! ice%pc_zeta        = 1._dp
    ! ice%pc_eta         = C%pc_epsilon
    ! ice%pc_eta_prev    = C%pc_epsilon

    ! ice%dHb_dt_a = 0._dp
    ! ice%dHi_dt_a = 0._dp
    ! ice%dHs_dt_a = 0._dp

    ! ! Initialise the "previous ice mask", so that the first call to thermodynamics works correctly
    ! ice%mask_ice_a_prev( mesh%vi1:mesh%vi2) = ice%mask_ice_a( mesh%vi1:mesh%vi2)
    ! call allgather_array(ice%mask_ice_a_prev)

    ! ! Allocate and initialise basal conditions
    ! call initialise_basal_conditions( mesh, ice)

    ! ! Geothermal heat flux
    ! select case (C%choice_geothermal_heat_flux)

    !   case ('constant')
    !     ! Uniform value over whole domain
    !     ice%GHF_a( mesh%vi1:mesh%vi2) = C%constant_geothermal_heat_flux

    !   case ('spatial')
    !     ! Spatially variable field
    !     call crash ('spatially variable GHF not yet implemented!')
    !     ! call map_geothermal_heat_flux_to_mesh( mesh, ice)

    !   case default
    !     ! Unknown case
    !     call crash('unknown choice_geothermal_heat_flux "' // &
    !                 trim( C%choice_geothermal_heat_flux) // '"!')

    ! end select

    ! ! Initialise data and matrices for the velocity solver(s)
    ! call initialise_velocity_solver( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ice_model

END MODULE ice_model_main
