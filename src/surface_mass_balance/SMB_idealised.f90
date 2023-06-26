MODULE SMB_idealised

  ! Idealised SMB models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE SMB_model_types                                        , ONLY: type_SMB_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_SMB_model_idealised( mesh, ice, SMB, time)
    ! Calculate the surface mass balance
    !
    ! Use an idealised SMB scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen idealised SMB model
    IF     (C%choice_SMB_model_idealised == 'EISMINT1_A' .OR. &
            C%choice_SMB_model_idealised == 'EISMINT1_B' .OR. &
            C%choice_SMB_model_idealised == 'EISMINT1_C' .OR. &
            C%choice_SMB_model_idealised == 'EISMINT1_D' .OR. &
            C%choice_SMB_model_idealised == 'EISMINT1_E' .OR. &
            C%choice_SMB_model_idealised == 'EISMINT1_F') THEN
      CALL run_SMB_model_idealised_EISMINT1( mesh, SMB, time)
    ELSE
      CALL crash('unknown choice_SMB_model_idealised "' // TRIM( C%choice_SMB_model_idealised) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_idealised

  SUBROUTINE run_SMB_model_idealised_EISMINT1( mesh, SMB, time)
    ! Calculate the surface mass balance
    !
    ! Use an idealised SMB scheme
    !
    ! SMB for the EISMINT1 experiments (Huybrechts et al., 1996)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_SMB_model_idealised_EISMINT1'
    INTEGER                                            :: vi
    REAL(dp), PARAMETER                                :: s = 1E-2_dp    ! Mass balance change with distance from divide [m yr^-1 km^-1]
    REAL(dp)                                           :: x, y, d, R_el, T

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_SMB_model_idealised == 'EISMINT1_A' ) THEN
      ! Moving margin, no cyclicity

      DO vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for moving margin experiments, use Euclidean distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = SQRT( x**2 + y**2) / 1E3_dp  ! [km]

        ! Calculate distance from equilibrium line to ice divide
        R_el = 450._dp

        ! Calculate SMB (Huybrechts et al., Eq. 10)
        SMB%SMB( vi) = MIN( 0.5_dp, s * (R_el - d))

      END DO

    ELSEIF (C%choice_SMB_model_idealised == 'EISMINT1_B' ) THEN
      ! Moving margin, 20,000-yr cyclicity

      T = 20E3_dp

      DO vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for moving margin experiments, use Euclidean distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = SQRT( x**2 + y**2) / 1E3_dp  ! [km]

        ! Calculate distance from equilibrium line to ice divide (Huybrechts et al., Eq. 14)
        R_el = 450._dp + 100._dp * SIN( 2 * pi * time / T)

        ! Calculate SMB (Huybrechts et al., Eq. 10)
        SMB%SMB( vi) = MIN( 0.5_dp, s * (R_el - d))

      END DO

    ELSEIF (C%choice_SMB_model_idealised == 'EISMINT1_C' ) THEN
      ! Moving margin, 40,000-yr cyclicity

      T = 40E3_dp

      DO vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for moving margin experiments, use Euclidean distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = SQRT( x**2 + y**2) / 1E3_dp  ! [km]

        ! Calculate distance from equilibrium line to ice divide (Huybrechts et al., Eq. 14)
        R_el = 450._dp + 100._dp * SIN( 2._dp * pi * time / T)

        ! Calculate SMB (Huybrechts et al., Eq. 10)
        SMB%SMB( vi) = MIN( 0.5_dp, s * (R_el - d))

      END DO

    ELSEIF (C%choice_SMB_model_idealised == 'EISMINT1_D' ) THEN
      ! Fixed margin, no cyclicity

      ! Calculate SMB (Huybrechts et al., Eq. 8)
      DO vi = mesh%vi1, mesh%vi2
        SMB%SMB( vi) = 0.3_dp
      END DO

    ELSEIF (C%choice_SMB_model_idealised == 'EISMINT1_E' ) THEN
      ! Fixed margin, 20,000-yr cyclicity

      T = 20E3_dp

      ! Calculate SMB (Huybrechts et al., Eq. 13)
      DO vi = mesh%vi1, mesh%vi2
        SMB%SMB( vi) = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / T)
      END DO

    ELSEIF (C%choice_SMB_model_idealised == 'EISMINT1_F' ) THEN
      ! Fixed margin, 40,000-yr cyclicity

      T = 40E3_dp

      ! Calculate SMB (Huybrechts et al., Eq. 13)
      DO vi = mesh%vi1, mesh%vi2
        SMB%SMB( vi) = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / T)
      END DO

    ELSE
      CALL crash('unknown choice_SMB_model_idealised "' // TRIM( C%choice_SMB_model_idealised) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

    END SUBROUTINE run_SMB_model_idealised_EISMINT1

  SUBROUTINE initialise_SMB_model_idealised( mesh, SMB)
    ! Initialise the SMB model
    !
    ! Use an idealised SMB scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(OUT)   :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising idealised surface mass balance model...'

    ! Run the chosen idealised SMB model
    IF (C%choice_SMB_model_idealised == 'EISMINT1_A' .OR. &
        C%choice_SMB_model_idealised == 'EISMINT1_B' .OR. &
        C%choice_SMB_model_idealised == 'EISMINT1_C' .OR. &
        C%choice_SMB_model_idealised == 'EISMINT1_D' .OR. &
        C%choice_SMB_model_idealised == 'EISMINT1_E' .OR. &
        C%choice_SMB_model_idealised == 'EISMINT1_F') THEN
      ! No need to do anything
    ELSE
      CALL crash('unknown choice_SMB_model_idealised "' // TRIM( C%choice_SMB_model_idealised) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model_idealised

END MODULE SMB_idealised
