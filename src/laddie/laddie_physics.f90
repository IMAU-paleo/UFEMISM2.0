MODULE laddie_physics

  ! Utilities for the laddie model

! ===== Preamble =====
! ====================
  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use checksum_mod, only: checksum

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE compute_melt_rate( mesh, laddie, npx, Hstar, time)
    ! Compute melt rate using the three equations

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx
    REAL(dp), DIMENSION(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), INTENT(IN)    :: Hstar
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_melt_rate'
    INTEGER                                               :: vi
    REAL(dp)                                              :: That, Chat, Ctil, Bval, Cval, Dval, AA
    REAL(dp), PARAMETER                                   :: nu0 = 1.95E-6
    REAL(dp), PARAMETER                                   :: eps = 1.0E-12 ! Some small parameter to prevent div. by zero
    REAL(dp), PARAMETER                                   :: tol = 1.0E-12 ! Some small parameter to prevent div. by zero

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get friction velocity
    do vi = mesh%vi1, mesh%vi2
      laddie%u_star( vi) = (C%laddie_drag_coefficient_top * (npx%U_a( vi)**2 + npx%V_a( vi)**2 + C%uniform_laddie_tidal_velocity**2 ))**.5
    end do
    call checksum( laddie%u_star, 'laddie%u_star', mesh%pai_V)

    ! Get gamma values
    SELECT CASE (C%choice_laddie_gamma)
      CASE DEFAULT
        CALL crash('unknown choice_laddie_gamma "' // TRIM( C%choice_laddie_gamma) // '"')
      CASE ('uniform')
        laddie%gamma_T( mesh%vi1:mesh%vi2) = C%uniform_laddie_gamma_T
        laddie%gamma_S( mesh%vi1:mesh%vi2) = C%uniform_laddie_gamma_T/35.0_dp
      CASE ('Jenkins1991')
        DO vi = mesh%vi1, mesh%vi2
           IF (laddie%mask_a( vi)) THEN
             AA = 2.12_dp*LOG(laddie%u_star( vi) * Hstar( vi)/nu0+eps)
             laddie%gamma_T( vi) = laddie%u_star( vi) / (AA + 12.5_dp * Prandtl_number**(2.0_dp/3) - 8.68_dp)
             laddie%gamma_S( vi) = laddie%u_star( vi) / (AA + 12.5_dp * Schmidt_number**(2.0_dp/3) - 8.68_dp)
           END IF
        END DO
    END SELECT
    call checksum( laddie%gamma_T, 'laddie%gamma_T', mesh%pai_V)
    call checksum( laddie%gamma_S, 'laddie%gamma_S', mesh%pai_V)

    ! == Get melt rate ==
    ! ===================
    Ctil = cp_ice/cp_ocean

    ! Loop over vertices
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         ! Solve three equations
         That = freezing_lambda_2 + freezing_lambda_3*laddie%Hib( vi)
         IF (time == C%start_time_of_run .OR. C%choice_thermo_model == 'none') THEN
           ! Ignore heat diffusion into ice
           Chat = cp_ocean / L_fusion
         ELSE
           Chat = cp_ocean / (L_fusion - cp_ice * laddie%Ti( vi, 1))
         END IF

         Bval = Chat*laddie%gamma_T( vi)*(That - npx%T( vi)) + laddie%gamma_S( vi)*(1 + Chat*Ctil*(That + freezing_lambda_1*npx%S( vi)))
         Cval = Chat*laddie%gamma_T( vi)*laddie%gamma_S( vi) * (That-npx%T( vi) + freezing_lambda_1*npx%S( vi))

         ! Get melt rate
         IF (4*Cval > Bval**2) THEN
           ! Probably not possible, but to prevent NaNs, set melt rate to zero
           laddie%melt( vi) = 0.0
         ELSE
           laddie%melt( vi) = 0.5_dp * (-Bval + SQRT(Bval**2 - 4.0_dp*Cval))
         END IF

         ! Get temperature at ice base
         Dval = laddie%melt( vi) * cp_ice - cp_ocean * laddie%gamma_T( vi)
         IF (ABS(Dval) < tol) THEN
           ! Seems like a very unlikely case, but better to be careful
           laddie%T_base( vi) = laddie%T_freeze( vi)
         ELSE
           IF (time == C%start_time_of_run .OR. C%choice_thermo_model == 'none') THEN
             ! Ignore heat diffusion into ice
             laddie%T_base( vi) = (laddie%melt( vi) * L_fusion - cp_ocean * laddie%gamma_T( vi) * npx%T( vi)) / Dval
           ELSE
             laddie%T_base( vi) = (laddie%melt( vi) * (L_fusion - cp_ice * laddie%Ti( vi, 1)) - cp_ocean * laddie%gamma_T( vi) * npx%T( vi)) / Dval
           END IF
         END IF

       END IF
    END DO
    call checksum( laddie%melt  , 'laddie%melt  ', mesh%pai_V)
    call checksum( laddie%T_base, 'laddie%T_base', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_melt_rate

  SUBROUTINE compute_entrainment( mesh, laddie, npx, Hstar)
    ! Compute entrainment rate

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx
    REAL(dp), DIMENSION(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_entrainment'
    INTEGER                                               :: vi
    REAL(dp), PARAMETER                                   :: maxdetr = 0.001_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Only Gaspar option for now

    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         ! Get salinity at ice base
         laddie%S_base( vi) = (laddie%T_base( vi) - freezing_lambda_2 - freezing_lambda_3 * laddie%Hib( vi)) / freezing_lambda_1

         ! Get buoyancy at ice base
         laddie%drho_base( vi) = C%uniform_laddie_eos_linear_beta  * (npx%S( vi)-laddie%S_base( vi)) &
                               - C%uniform_laddie_eos_linear_alpha * (npx%T( vi)-laddie%T_base( vi))

         ! Get entrainment
         laddie%entr( vi) = 2*C%laddie_Gaspar1988_mu/grav &
                          * laddie%u_star( vi)**3 / (Hstar( vi) * laddie%drho_amb( vi)) &
                          - laddie%drho_base( vi) / laddie%drho_amb( vi) * laddie%melt( vi)

         ! Prevent too strong detrainment
         laddie%entr( vi) = MAX(laddie%entr( vi), -maxdetr)

         ! Get detrainment = negative component of entrainment
         laddie%detr( vi) = - MIN(laddie%entr( vi),0.0_dp)
       END IF
    END DO
    call checksum( laddie%S_base   , 'laddie%S_base   ', mesh%pai_V)
    call checksum( laddie%drho_base, 'laddie%drho_base', mesh%pai_V)
    call checksum( laddie%entr     , 'laddie%entr     ', mesh%pai_V)
    call checksum( laddie%detr     , 'laddie%detr     ', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_entrainment

  SUBROUTINE compute_subglacial_discharge( mesh, laddie)
  ! Compute subglacial discharge (SGD)
  ! TODO clean up routine; avoid so many if statements
  ! TODO allow option to read in SGD mask from file

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_subglacial_discharge'
    INTEGER                                               :: vi, ierr
    REAL(dp)                                              :: total_area 

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise total_area (the area over which the SGD will be distributed) at zero
    total_area = 0._dp

    ! Initialise SGD at zero
    laddie%SGD = 0._dp

    ! Determine total_area by looping over the vertices
    DO vi = mesh%vi1, mesh%vi2 
      IF (laddie%mask_a( vi)) THEN
        IF (laddie%mask_gl_fl( vi) .and. laddie%mask_SGD( vi)) THEN
          total_area = total_area + mesh%A( vi) 
        END IF
      END IF 
    END DO

    ! Make sure processes have total_area 
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_area,      1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    IF (total_area == 0._dp) THEN
      ! No SGD is applied, issue a warning
      IF (par%primary) CALL warning('No subglacial discharge (SGD) is applied because the total_area where SGD is applied is zero!')
    ELSEIF (total_area > 0._dp) THEN
      ! Distribute SGD flux [m^3/s] over the total area to get the SGD in [m/s]
      DO vi = mesh%vi1, mesh%vi2 
        IF (laddie%mask_a( vi)) THEN
          IF (laddie%mask_gl_fl( vi) .and. laddie%mask_SGD( vi)) THEN
            laddie%SGD( vi) = C%laddie_SGD_flux / total_area
          END IF
        END IF 
      END DO  
    ELSE
      CALL crash('The total_area where SGD is applied is a negative value')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_subglacial_discharge

  SUBROUTINE compute_freezing_temperature( mesh, laddie, npx)
    ! Compute freezing temperature at ice shelf base, based on Laddie salinity.
    ! TODO can maybe be merged with ice computation

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_freezing_temperature'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         laddie%T_freeze( vi) = freezing_lambda_1*npx%S( vi) + freezing_lambda_2 + freezing_lambda_3*laddie%Hib( vi)
       END IF
    END DO
    call checksum( laddie%T_freeze, 'laddie%T_freeze', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_freezing_temperature

  SUBROUTINE compute_buoyancy( mesh, laddie, npx, Hstar)
    ! Compute buoyancy = (rho_amb - rho)/rho_sw
    ! TODO update with Roquet EOS

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    TYPE(type_laddie_timestep),             INTENT(IN)    :: npx
    REAL(dp), DIMENSION(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_buoyancy'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    laddie%drho_amb( mesh%vi1:mesh%vi2) = 0.0_dp

    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         ! Get buoyancy
         laddie%drho_amb( vi) = C%uniform_laddie_eos_linear_beta  * (laddie%S_amb( vi)-npx%S( vi)) &
                              - C%uniform_laddie_eos_linear_alpha * (laddie%T_amb( vi)-npx%T( vi))

         ! Make sure buoyancy is positive TODO expand with convection scheme
         laddie%drho_amb( vi) = MAX(laddie%drho_amb( vi),C%laddie_buoyancy_minimum/seawater_density)

         ! Get depth-integrated buoyancy based on Hstar
         laddie%Hdrho_amb( vi) = Hstar( vi) * laddie%drho_amb( vi)

       END IF
    END DO
    call checksum( laddie%drho_amb , 'laddie%drho_amb ', mesh%pai_V)
    call checksum( laddie%Hdrho_amb, 'laddie%Hdrho_amb', mesh%pai_V)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_buoyancy

END MODULE laddie_physics

