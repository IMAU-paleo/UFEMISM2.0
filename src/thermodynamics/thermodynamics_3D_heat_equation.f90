MODULE thermodynamics_3D_heat_equation

  ! All the routines to solve the 3-D heat equation in the ice

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF, &
                                                                     save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE SMB_types                                              , ONLY: type_SMB_model
  USE BMB_types                                              , ONLY: type_BMB_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE solve_3D_heat_equation( mesh, ice, SMB, BMB, dt)
    ! Solve the 3-D heat equation to calculate temperature at time t+dt

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                   INTENT(IN)    :: BMB
    REAL(dp),                               INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'solve_3D_heat_equation'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_3D_heat_equation

END MODULE thermodynamics_3D_heat_equation
