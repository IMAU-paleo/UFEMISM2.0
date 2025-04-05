MODULE BMB_parameterised

  ! parameterised BMB models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE BMB_model_types                                        , ONLY: type_BMB_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model_parameterised( mesh, ice, ocean, BMB)
    ! Calculate the basal mass balance
    !
    ! Use a parameterised BMB scheme

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_model),              INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_parameterised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen parameterised BMB model
    SELECT CASE (C%choice_BMB_model_parameterised)
      CASE ('Favier2019')
        CALL run_BMB_model_parameterised_Favier2019( mesh, ice, ocean, BMB)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_parameterised "' // TRIM( C%choice_BMB_model_parameterised) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_parameterised

  subroutine run_BMB_model_parameterised_Favier2019( mesh, ice, ocean, BMB)
    ! The basal melt parameterisation used in Favier et al. (2019)

    ! In/output variables
    type(type_mesh),                     intent(in)    :: mesh
    type(type_ice_model),                intent(in)    :: ice
    type(type_ocean_model),              intent(in)    :: ocean
    type(type_BMB_model),                intent(inout) :: BMB

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'run_BMB_model_parameterised_Favier2019'
    integer                                            :: vi
    real(dp)                                           :: dT

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    BMB%BMB_shelf = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if ((ocean%T_draft( vi) == ocean%T_draft( vi)) .and. (ocean%T_freezing_point( vi) == ocean%T_freezing_point( vi))) then
        ! Both parameters are non-NaN

        ! Temperature forcing
        dT = ocean%T_draft( vi) - ocean%T_freezing_point( vi)

        ! Favier et al. (2019), Eq. 4
        ! Altered to allow for negative basal melt (i.e. refreezing) when dT < 0
        BMB%BMB_shelf( vi) =  -1._dp * sec_per_year * C%BMB_Favier2019_gamma * sign(dT,1._dp) * (seawater_density * cp_ocean * dT / (ice_density * L_fusion))**2._dp

        ! Apply grounded fractions
        if (ice%mask_gl_gr( vi) .and. ice%Hib(vi) < ice%SL(vi)) then
          ! Subgrid basal melt rate
          ! BMB%BMB_shelf( vi) = (1._dp - ice%fraction_gr( vi)) * BMB%BMB_shelf( vi)
          ! Limit it to only melt (refreezing is tricky)
          BMB%BMB_shelf( vi) = max( BMB%BMB_shelf( vi), 0._dp)
        end if
      else
        ! Either parameter is NaN, output zero BMB
        BMB%BMB_shelf( vi) = 0._dp

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_BMB_model_parameterised_Favier2019

  SUBROUTINE initialise_BMB_model_parameterised( mesh, BMB)
    ! Initialise the BMB model
    !
    ! Use a parameterised BMB scheme

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_parameterised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising parameterised BMB model "' // &
      colour_string( TRIM( C%choice_BMB_model_parameterised),'light blue') // '"...'

    ! Initialise the chosen parameterised BMB model
    SELECT CASE (C%choice_BMB_model_parameterised)
      CASE ('Favier2019')
        ! No need to do anything
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_parameterised "' // TRIM( C%choice_BMB_model_parameterised) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_parameterised

END MODULE BMB_parameterised
