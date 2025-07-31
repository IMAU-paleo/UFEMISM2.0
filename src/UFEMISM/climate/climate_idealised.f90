MODULE climate_idealised

  ! Idealised climate models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model_idealised( mesh, ice, climate, time)
    ! Calculate the climate
    !
    ! Use an idealised climate scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen idealised climate model
    IF     (C%choice_climate_model_idealised == 'EISMINT1_A' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_C' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_D' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_E' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_F') THEN
      CALL run_climate_model_idealised_EISMINT1( mesh, ice, climate, time)
    ELSE
      CALL crash('unknown choice_climate_model_idealised "' // TRIM( C%choice_climate_model_idealised) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_idealised

  SUBROUTINE initialise_climate_model_idealised( mesh, climate)
    ! Initialise the climate model
    !
    ! Use an idealised climate scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(INOUT) :: climate

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '     Initialising idealised climate model "' // &
      colour_string( TRIM( C%choice_climate_model_idealised),'light blue') // '"...'

    ! Run the chosen idealised climate model
    IF (C%choice_climate_model_idealised == 'EISMINT1_A' .OR. &
        C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
        C%choice_climate_model_idealised == 'EISMINT1_C' .OR. &
        C%choice_climate_model_idealised == 'EISMINT1_D' .OR. &
        C%choice_climate_model_idealised == 'EISMINT1_E' .OR. &
        C%choice_climate_model_idealised == 'EISMINT1_F') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_climate_model_idealised "' // TRIM( C%choice_climate_model_idealised) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_idealised

  ! == EISMINT1
  ! ===========

  SUBROUTINE run_climate_model_idealised_EISMINT1( mesh, ice, climate, time)
    ! Temperature for the EISMINT1 experiments (Huybrechts et al., 1996)

    IMPLICIT NONE

    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_model),             INTENT(INOUT) :: climate
    REAL(dp),                             INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_climate_model_idealised_EISMINT1'
    REAL(dp), PARAMETER                                 :: x_summit = 0._dp      ! x-coordinate of ice divide [m]
    REAL(dp), PARAMETER                                 :: y_summit = 0._dp      ! y-coordinate of ice divide [m]
    INTEGER                                             :: vi
    REAL(dp)                                            :: x, y, d, T, dT

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set precipitation to zero - SMB is parameterised anyway...
    climate%Precip = 0._dp

    ! Baseline temperature
    IF     (C%choice_climate_model_idealised == 'EISMINT1_A' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_C') THEN
      ! Moving margin: Huybrechts et al., Eq. 11

      DO vi = mesh%vi1, mesh%vi2
        ! Calculate baseline temperature
        climate%T2m( vi,:) = 270._dp - 0.01_dp * ice%Hs( vi)
      END DO

    ELSEIF (C%choice_climate_model_idealised == 'EISMINT1_D' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_E' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_F') THEN
      ! Fixed margin: Huybrechts et al., Eq. 9

      DO vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for fixed margin experiments, use square distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = MAX( ABS( x - x_summit), ABS( y - y_summit)) / 1E3_dp  ! [km]

        ! Calculate baseline temperature
        climate%T2m( vi,:) = 239._dp + (8.0E-08_dp * d**3)

      END DO

    ELSE
      CALL crash('unknown choice_climate_model_idealised"' // TRIM( C%choice_climate_model_idealised) // '"!')
    END IF
    CALL sync

    ! Add temperature change for the glacial cycle experiments
    IF     (C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_E') THEN
      ! 20,000-yr cyclicity

      T = 20E3_dp

      IF (time > 0._dp) THEN
        dT = 10._dp * SIN( 2._dp * pi * time / T)
        climate%T2m = climate%T2m + dT
      END IF

    ELSEIF (C%choice_climate_model_idealised == 'EISMINT1_C' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_F') THEN
      ! 40,000-yr cyclicity

      T = 40E3_dp

      IF (time > 0._dp) THEN
        dT = 10._dp * SIN( 2._dp * pi * time / T)
        climate%T2m = climate%T2m + dT
      END IF

    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_idealised_EISMINT1

END MODULE climate_idealised
