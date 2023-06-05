MODULE basal_hydrology

  ! Contains all the different basal hydrology models.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reallocate_mod                                         , ONLY: reallocate_clean_dp_1D

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_basal_hydrology( mesh, ice)
    ! Calculate the pore water pressure and effective basal pressure

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_hydrology'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================

    IF     (C%choice_basal_hydrology == 'saturated') THEN
      ! Assume all marine till is saturated (i.e. pore water pressure is equal to water pressure at depth everywhere)
      CALL calc_pore_water_pressure_saturated( mesh, ice)
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      ! The Martin et al. (2011) parameterisation of pore water pressure
      CALL calc_pore_water_pressure_Martin2011( mesh, ice)
    ELSE
      CALL crash('unknown choice_basal_hydrology "' // TRIM( C%choice_basal_hydrology) // '"!')
    END IF

    ! Calculate overburden and effective pressure
    ! ===========================================

    DO vi = mesh%vi1, mesh%vi2
      ice%overburden_pressure( vi) = ice_density * grav * ice%Hi( vi)
      ice%effective_pressure(  vi) = MAX( 0._dp, ice%overburden_pressure( vi) - ice%pore_water_pressure( vi))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_hydrology

  SUBROUTINE initialise_basal_hydrology( mesh, ice)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_hydrology'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    IF     (C%choice_basal_hydrology == 'saturated') THEN
      ALLOCATE( ice%pore_water_pressure( mesh%vi1:mesh%vi2))
      ALLOCATE( ice%overburden_pressure( mesh%vi1:mesh%vi2))
      ALLOCATE( ice%effective_pressure(  mesh%vi1:mesh%vi2))
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      ALLOCATE( ice%pore_water_pressure( mesh%vi1:mesh%vi2))
      ALLOCATE( ice%overburden_pressure( mesh%vi1:mesh%vi2))
      ALLOCATE( ice%effective_pressure(  mesh%vi1:mesh%vi2))
    ELSE
      CALL crash('unknown choice_basal_hydrology "' // TRIM( C%choice_basal_hydrology) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_hydrology

  SUBROUTINE calc_pore_water_pressure_saturated( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Assume all till is saturated, i.e. pore water pressure = -rho_w * g * Hb

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_saturated'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure( vi) = -seawater_density * grav * (ice%SL( vi) - ice%Hb( vi))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_saturated

  SUBROUTINE calc_pore_water_pressure_Martin2011( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Use the parameterisation from Martin et al. (2011)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_Martin2011'
    INTEGER                                            :: vi
    REAL(dp)                                           :: lambda_p

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      lambda_p = MIN( 1._dp, MAX( 0._dp, 1._dp - (ice%Hb( vi) - ice%SL( vi) - C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - C%Martin2011_hydro_Hb_min) ))

      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure( vi) = 0.96_dp * ice_density * grav * ice%Hi( vi) * lambda_p

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_Martin2011

  SUBROUTINE remap_basal_hydrology( mesh_old, mesh_new, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_basal_hydrology'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV

    ! Allocate shared memory
    IF     (C%choice_basal_hydrology == 'saturated') THEN
      CALL reallocate_clean_dp_1D( ice%pore_water_pressure, mesh_new%nV_loc)
      CALL reallocate_clean_dp_1D( ice%overburden_pressure, mesh_new%nV_loc)
      CALL reallocate_clean_dp_1D( ice%effective_pressure              , mesh_new%nV_loc)
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      CALL reallocate_clean_dp_1D( ice%pore_water_pressure, mesh_new%nV_loc)
      CALL reallocate_clean_dp_1D( ice%overburden_pressure, mesh_new%nV_loc)
      CALL reallocate_clean_dp_1D( ice%effective_pressure              , mesh_new%nV_loc)
    ELSE
      CALL crash('unknown choice_basal_hydrology "' // TRIM( C%choice_basal_hydrology) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_basal_hydrology

END MODULE basal_hydrology
