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
  USE mesh_operators                                         , ONLY: ddx_a_a_2D, ddy_a_a_2D, map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, &
                                                                     ddx_b_a_2D, ddy_b_a_2D

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

  SUBROUTINE calc_zeta_gradients( mesh, ice)
    ! Calculate all the gradients of zeta, needed to perform the scaled vertical coordinate transformation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_zeta_gradients'
!    REAL(dp), DIMENSION(:    ), POINTER                :: Hi_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dx2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dxdy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dy2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: Hi_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dx2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dxdy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dy2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: Hs_b
!    REAL(dp), DIMENSION(:    ), POINTER                :: Hs_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dx2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dxdy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dy2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dx2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dxdy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dy2_b
    INTEGER                                            :: vi,ti,k,ks
    REAL(dp)                                           :: Hi, zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
   !ALLOCATE( Hi_a(        mesh%vi1:mesh%vi2))
    ALLOCATE( dHi_dx_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( dHi_dy_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hi_dx2_a(  mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hi_dxdy_a( mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hi_dy2_a(  mesh%vi1:mesh%vi2))

    ALLOCATE( Hi_b(        mesh%ti1:mesh%ti2))
    ALLOCATE( dHi_dx_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( dHi_dy_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hi_dx2_b(  mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hi_dxdy_b( mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hi_dy2_b(  mesh%ti1:mesh%ti2))

   !ALLOCATE( Hs_a(        mesh%vi1:mesh%vi2))
    ALLOCATE( dHs_dx_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( dHs_dy_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hs_dx2_a(  mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hs_dxdy_a( mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hs_dy2_a(  mesh%vi1:mesh%vi2))

    ALLOCATE( Hs_b(        mesh%ti1:mesh%ti2))
    ALLOCATE( dHs_dx_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( dHs_dy_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hs_dx2_b(  mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hs_dxdy_b( mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hs_dy2_b(  mesh%ti1:mesh%ti2))

  ! Calculate gradients of Hi and Hs on both grids

    !CALL map_a_a_2D( mesh, ice%Hi_a, Hi_a       )
    CALL ddx_a_a_2D( mesh, ice%Hi  , dHi_dx_a   )
    CALL ddy_a_a_2D( mesh, ice%Hi  , dHi_dy_a   )
    CALL map_a_b_2D( mesh, ice%Hi  , Hi_b       )
    CALL ddx_a_b_2D( mesh, ice%Hi  , dHi_dx_b   )
    CALL ddy_a_b_2D( mesh, ice%Hi  , dHi_dy_b   )
    CALL ddx_b_a_2D( mesh, dHi_dx_b, d2Hi_dx2_a )
    CALL ddy_b_a_2D( mesh, dHi_dx_b, d2Hi_dxdy_a)
    CALL ddy_b_a_2D( mesh, dHi_dy_b, d2Hi_dy2_a )
    CALL ddx_a_b_2D( mesh, dHi_dx_a, d2Hi_dx2_b )
    CALL ddy_a_b_2D( mesh, dHi_dx_a, d2Hi_dxdy_b)
    CALL ddy_a_b_2D( mesh, dHi_dy_a, d2Hi_dy2_b )

    !CALL map_a_a_2D( mesh, ice%Hs_a, Hs_a       )
    CALL ddx_a_a_2D( mesh, ice%Hs  , dHs_dx_a   )
    CALL ddy_a_a_2D( mesh, ice%Hs  , dHs_dy_a   )
    CALL map_a_b_2D( mesh, ice%Hs  , Hs_b       )
    CALL ddx_a_b_2D( mesh, ice%Hs  , dHs_dx_b   )
    CALL ddy_a_b_2D( mesh, ice%Hs  , dHs_dy_b   )
    CALL ddx_b_a_2D( mesh, dHs_dx_b, d2Hs_dx2_a )
    CALL ddy_b_a_2D( mesh, dHs_dx_b, d2Hs_dxdy_a)
    CALL ddy_b_a_2D( mesh, dHs_dy_b, d2Hs_dy2_a )
    CALL ddx_a_b_2D( mesh, dHs_dx_a, d2Hs_dx2_b )
    CALL ddy_a_b_2D( mesh, dHs_dx_a, d2Hs_dxdy_b)
    CALL ddy_a_b_2D( mesh, dHs_dy_a, d2Hs_dy2_b )

    ! Calculate zeta gradients on all grids

    ! ak
    DO vi = mesh%vi1, mesh%vi2

      Hi = MAX( 10._dp, ice%Hi( vi))

      DO k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dt_ak(    vi,k) = ( 1._dp / Hi) * (ice%dHs_dt( vi) - zeta * ice%dHi_dt( vi))

        ice%dzeta_dx_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dx_a( vi) - zeta * dHi_dx_a( vi))
        ice%dzeta_dy_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dy_a( vi) - zeta * dHi_dy_a( vi))
        ice%dzeta_dz_ak(    vi,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_ak(  vi,k) = (dHi_dx_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dx2_a(  vi) - zeta * d2Hi_dx2_a(  vi))
        ice%d2zeta_dxdy_ak( vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dxdy_a( vi) - zeta * d2Hi_dxdy_a( vi))
        ice%d2zeta_dy2_ak(  vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dy_ak( vi,k) + (1._dp / Hi) * (d2Hs_dy2_a(  vi) - zeta * d2Hi_dy2_a(  vi))

      END DO
    END DO

    ! bk
    DO ti = mesh%ti1, mesh%ti2

      Hi = MAX( 10._dp, Hi_b( ti))

      DO k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dx_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bk(    ti,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_bk(  ti,k) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bk( ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bk(  ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bk( ti,k) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      END DO
    END DO

    ! bks
    DO ti = mesh%ti1, mesh%ti2

      Hi = MAX( 10._dp, Hi_b( ti))

      DO ks = 1, mesh%nz-1

        zeta = mesh%zeta_stag( ks)

        ice%dzeta_dx_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bks(    ti,ks) = (-1._dp / Hi)

        ice%d2zeta_dx2_bks(  ti,ks) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bks( ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bks(  ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_zeta_gradients

END MODULE ice_model_main
