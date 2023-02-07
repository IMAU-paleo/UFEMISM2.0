MODULE netcdf_debug_module

! ===== Creating and writing to debug files =====
! ===============================================
!
! These routines create and write data to the NetCDF debug files.
! At pretty much any point in the model code, any data field
! on the mesh (both A-, B-, and C-grids) can be immediately
! written to the NetCDF debug file for inspection. This is by
! far the most useful developer's tool available in UFEMISM.
!
! Example: say we want to see what's going on with the field ice%Hi_a:
!
!   debug%dp_2D_a_01 = ice%Hi_a
!   CALL write_to_debug_file
!
! And voila, there you go.
! NOTE: write_to_debug_file, like all NetCDF routines, is shared; if called
!       by only a single process, the program will freeze.
!
! The debug structure is a global variable of this module; in order to
! be able to access it from other modules, import it from this module:
!
!   USE netcdf_debug_module, ONLY: debug, write_to_debug_file
!
! The fact that the different model regions have different meshes, and
! therefore need different output files, is already taken care of;
! the structure "debug" actually contains pointers that point to
! either debug_nam, debug_eas, debug_grl, or debug_ant; binding these
! pointers is done at the start of run_model.
!
! A new debug file is automatically created, and the old one discarded,
! when the mesh is updated.
!
! =====================================
!
! Additionally: single variables can now also be immediately saved as
!               NetCDF files. They will have no dimensions at all, so beware!
!
! Additionally: CSR and PETSc matrices can be written directly to NetCDF,
!               to allow for smooth comparison with Matlab prototype code.

! ===== Preamble =====
! ====================

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE petscksp
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
                                             allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                             allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                             allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                             allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                             allocate_shared_int_4D,   allocate_shared_dp_4D, &
                                             allocate_shared_bool_0D,  allocate_shared_bool_1D, &
                                             reallocate_shared_int_0D, reallocate_shared_dp_0D, &
                                             reallocate_shared_int_1D, reallocate_shared_dp_1D, &
                                             reallocate_shared_int_2D, reallocate_shared_dp_2D, &
                                             reallocate_shared_int_3D, reallocate_shared_dp_3D, &
                                             deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_debug_fields, type_model_region, type_mesh, type_sparse_matrix_CSR_dp, &
                                             type_grid
  USE netcdf,                          ONLY: NF90_INT, NF90_DOUBLE
  USE netcdf_basic_module,             ONLY: create_new_netcdf_file_for_writing, create_dimension, create_variable, &
                                             open_existing_netcdf_file_for_writing, close_netcdf_file, &
                                             write_var_int_1D, write_var_dp_1D, &
                                             write_var_int_2D, write_var_dp_2D, &
                                             write_var_int_3D, write_var_dp_3D
  USE netcdf_output_module,            ONLY: setup_mesh_in_netcdf_file, setup_xy_grid_in_netcdf_file, &
                                             add_zeta_dimension_to_file, add_month_dimension_to_file, &
                                             add_field_grid_int_2D_notime, add_field_grid_dp_2D_notime, &
                                             add_field_grid_dp_2D_monthly_notime, add_field_grid_dp_3D_notime, &
                                             add_field_mesh_int_2D_notime, add_field_mesh_int_2D_b_notime, add_field_mesh_int_2D_c_notime, &
                                             add_field_mesh_dp_2D_notime , add_field_mesh_dp_2D_b_notime , add_field_mesh_dp_2D_c_notime, &
                                             add_field_mesh_dp_2D_monthly_notime, add_field_mesh_dp_3D_notime, &
                                             write_to_field_multiple_options_mesh_int_2D_notime, &
                                             write_to_field_multiple_options_mesh_int_2D_b_notime, &
                                             write_to_field_multiple_options_mesh_int_2D_c_notime, &
                                             write_to_field_multiple_options_mesh_dp_2D_notime, &
                                             write_to_field_multiple_options_mesh_dp_2D_b_notime, &
                                             write_to_field_multiple_options_mesh_dp_2D_c_notime, &
                                             write_to_field_multiple_options_mesh_dp_2D_monthly_notime, &
                                             write_to_field_multiple_options_mesh_dp_3D_notime
  USE petsc_module,                    ONLY: mat_petsc2CSR
  USE sparse_matrix_module,            ONLY: deallocate_matrix_CSR

  IMPLICIT NONE

  TYPE(type_debug_fields) :: debug_NAM, debug_EAS, debug_GRL, debug_ANT, debug

CONTAINS

! ===== Debug NetCDF file =====
! =============================

  ! Main routines: create and write to NetCDF debug file
  SUBROUTINE write_to_debug_file
    ! Write the current set of debug data fields to the debug NetCDF file

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_debug_file'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( debug%filename, ncid)

    ! Write data
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_01', debug%int_2D_a_01)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_02', debug%int_2D_a_02)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_03', debug%int_2D_a_03)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_04', debug%int_2D_a_04)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_05', debug%int_2D_a_05)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_06', debug%int_2D_a_06)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_07', debug%int_2D_a_07)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_08', debug%int_2D_a_08)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_09', debug%int_2D_a_09)
    CALL write_to_field_multiple_options_mesh_int_2D_notime( debug%filename, ncid, 'int_2D_a_10', debug%int_2D_a_10)

    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_01', debug%int_2D_b_01)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_02', debug%int_2D_b_02)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_03', debug%int_2D_b_03)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_04', debug%int_2D_b_04)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_05', debug%int_2D_b_05)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_06', debug%int_2D_b_06)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_07', debug%int_2D_b_07)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_08', debug%int_2D_b_08)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_09', debug%int_2D_b_09)
    CALL write_to_field_multiple_options_mesh_int_2D_b_notime( debug%filename, ncid, 'int_2D_b_10', debug%int_2D_b_10)

    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_01', debug%int_2D_c_01)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_02', debug%int_2D_c_02)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_03', debug%int_2D_c_03)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_04', debug%int_2D_c_04)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_05', debug%int_2D_c_05)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_06', debug%int_2D_c_06)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_07', debug%int_2D_c_07)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_08', debug%int_2D_c_08)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_09', debug%int_2D_c_09)
    CALL write_to_field_multiple_options_mesh_int_2D_c_notime( debug%filename, ncid, 'int_2D_c_10', debug%int_2D_c_10)

    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_01', debug%dp_2D_a_01)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_02', debug%dp_2D_a_02)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_03', debug%dp_2D_a_03)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_04', debug%dp_2D_a_04)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_05', debug%dp_2D_a_05)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_06', debug%dp_2D_a_06)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_07', debug%dp_2D_a_07)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_08', debug%dp_2D_a_08)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_09', debug%dp_2D_a_09)
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( debug%filename, ncid, 'dp_2D_a_10', debug%dp_2D_a_10)

    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_01', debug%dp_2D_b_01)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_02', debug%dp_2D_b_02)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_03', debug%dp_2D_b_03)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_04', debug%dp_2D_b_04)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_05', debug%dp_2D_b_05)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_06', debug%dp_2D_b_06)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_07', debug%dp_2D_b_07)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_08', debug%dp_2D_b_08)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_09', debug%dp_2D_b_09)
    CALL write_to_field_multiple_options_mesh_dp_2D_b_notime( debug%filename, ncid, 'dp_2D_b_10', debug%dp_2D_b_10)

    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_01', debug%dp_2D_c_01)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_02', debug%dp_2D_c_02)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_03', debug%dp_2D_c_03)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_04', debug%dp_2D_c_04)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_05', debug%dp_2D_c_05)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_06', debug%dp_2D_c_06)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_07', debug%dp_2D_c_07)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_08', debug%dp_2D_c_08)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_09', debug%dp_2D_c_09)
    CALL write_to_field_multiple_options_mesh_dp_2D_c_notime( debug%filename, ncid, 'dp_2D_c_10', debug%dp_2D_c_10)

    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_01', debug%dp_2D_monthly_a_01)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_02', debug%dp_2D_monthly_a_02)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_03', debug%dp_2D_monthly_a_03)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_04', debug%dp_2D_monthly_a_04)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_05', debug%dp_2D_monthly_a_05)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_06', debug%dp_2D_monthly_a_06)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_07', debug%dp_2D_monthly_a_07)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_08', debug%dp_2D_monthly_a_08)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_09', debug%dp_2D_monthly_a_09)
    CALL write_to_field_multiple_options_mesh_dp_2D_monthly_notime( debug%filename, ncid, 'dp_2D_monthly_a_10', debug%dp_2D_monthly_a_10)

    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_01', debug%dp_3D_a_01)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_02', debug%dp_3D_a_02)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_03', debug%dp_3D_a_03)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_04', debug%dp_3D_a_04)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_05', debug%dp_3D_a_05)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_06', debug%dp_3D_a_06)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_07', debug%dp_3D_a_07)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_08', debug%dp_3D_a_08)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_09', debug%dp_3D_a_09)
    CALL write_to_field_multiple_options_mesh_dp_3D_notime( debug%filename, ncid, 'dp_3D_a_10', debug%dp_3D_a_10)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_debug_file

  SUBROUTINE create_debug_file( region)
    ! Create the NetCDF debug file for this model region.
    ! if one already exists, delete it (assume this is because we've just done a mesh update)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_debug_file'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (par%master) WRITE(0,*) '  Creating new debug file based on new mesh...'

    ! Determine NetCDF debug file name for this model region
    filename = 'debug_NAM.nc'
    filename(7:9) = region%name
    filename = TRIM( C%output_dir) // TRIM( filename)

    ! Delete existing debug file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the mesh in this file
    CALL setup_mesh_in_netcdf_file( filename, ncid, region%mesh)

    ! Add zeta and month dimensions
    CALL add_zeta_dimension_to_file(  filename, ncid)
    CALL add_month_dimension_to_file( filename, ncid)

    ! Add field variables
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_01')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_02')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_03')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_04')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_05')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_06')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_07')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_08')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_09')
    CALL add_field_mesh_int_2D_notime( filename, ncid, 'int_2D_a_10')

    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_01')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_02')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_03')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_04')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_05')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_06')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_07')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_08')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_09')
    CALL add_field_mesh_int_2D_b_notime( filename, ncid, 'int_2D_b_10')

    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_01')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_02')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_03')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_04')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_05')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_06')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_07')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_08')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_09')
    CALL add_field_mesh_int_2D_c_notime( filename, ncid, 'int_2D_c_10')

    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_01')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_02')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_03')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_04')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_05')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_06')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_07')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_08')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_09')
    CALL add_field_mesh_dp_2D_notime( filename, ncid, 'dp_2D_a_10')

    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_01')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_02')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_03')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_04')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_05')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_06')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_07')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_08')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_09')
    CALL add_field_mesh_dp_2D_b_notime( filename, ncid, 'dp_2D_b_10')

    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_01')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_02')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_03')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_04')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_05')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_06')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_07')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_08')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_09')
    CALL add_field_mesh_dp_2D_c_notime( filename, ncid, 'dp_2D_c_10')

    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_01')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_02')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_03')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_04')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_05')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_06')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_07')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_08')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_09')
    CALL add_field_mesh_dp_2D_monthly_notime( filename, ncid, 'dp_2D_monthly_a_10')

    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_01')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_02')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_03')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_04')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_05')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_06')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_07')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_08')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_09')
    CALL add_field_mesh_dp_3D_notime( filename, ncid, 'dp_3D_a_10')

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_debug_file

  ! Manage memory for the debug data fields
  SUBROUTINE associate_debug_fields( region)
    ! Since the dimensions vary, each region needs its own set of debug fields. However, if
    ! we make them part of the "region" TYPE, they need to be passed to every subroutine as an
    ! argument before they can be used, which is a lot of hassle. So instead they are saved as
    ! global variables of this module, where they can be accessed from anywhere. This is done
    ! via the "intermediary" set of pointers, which are bound to the region-specific debug structure
    ! with this here subroutine.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'associate_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set the filename
    IF     (region%name == 'NAM') THEN
      debug%filename = TRIM( C%output_dir) // 'debug_NAM.nc'
    ELSEIF (region%name == 'EAS') THEN
      debug%filename = TRIM( C%output_dir) // 'debug_EAS.nc'
    ELSEIF (region%name == 'GRL') THEN
      debug%filename = TRIM( C%output_dir) // 'debug_GRL.nc'
    ELSEIF (region%name == 'ANT') THEN
      debug%filename = TRIM( C%output_dir) // 'debug_ANT.nc'
    END IF

    ! If necessary (i.e. every time except the first ever time this subroutine is called), de-associate the intermediary pointers first.
    IF (ASSOCIATED(debug%dp_2D_a_01)) THEN

      NULLIFY( debug%int_2D_a_01)
      NULLIFY( debug%int_2D_a_02)
      NULLIFY( debug%int_2D_a_03)
      NULLIFY( debug%int_2D_a_04)
      NULLIFY( debug%int_2D_a_05)
      NULLIFY( debug%int_2D_a_06)
      NULLIFY( debug%int_2D_a_07)
      NULLIFY( debug%int_2D_a_08)
      NULLIFY( debug%int_2D_a_09)
      NULLIFY( debug%int_2D_a_10)

      NULLIFY( debug%int_2D_b_01)
      NULLIFY( debug%int_2D_b_02)
      NULLIFY( debug%int_2D_b_03)
      NULLIFY( debug%int_2D_b_04)
      NULLIFY( debug%int_2D_b_05)
      NULLIFY( debug%int_2D_b_06)
      NULLIFY( debug%int_2D_b_07)
      NULLIFY( debug%int_2D_b_08)
      NULLIFY( debug%int_2D_b_09)
      NULLIFY( debug%int_2D_b_10)

      NULLIFY( debug%int_2D_c_01)
      NULLIFY( debug%int_2D_c_02)
      NULLIFY( debug%int_2D_c_03)
      NULLIFY( debug%int_2D_c_04)
      NULLIFY( debug%int_2D_c_05)
      NULLIFY( debug%int_2D_c_06)
      NULLIFY( debug%int_2D_c_07)
      NULLIFY( debug%int_2D_c_08)
      NULLIFY( debug%int_2D_c_09)
      NULLIFY( debug%int_2D_c_10)

      NULLIFY( debug%dp_2D_a_01)
      NULLIFY( debug%dp_2D_a_02)
      NULLIFY( debug%dp_2D_a_03)
      NULLIFY( debug%dp_2D_a_04)
      NULLIFY( debug%dp_2D_a_05)
      NULLIFY( debug%dp_2D_a_06)
      NULLIFY( debug%dp_2D_a_07)
      NULLIFY( debug%dp_2D_a_08)
      NULLIFY( debug%dp_2D_a_09)
      NULLIFY( debug%dp_2D_a_10)

      NULLIFY( debug%dp_2D_b_01)
      NULLIFY( debug%dp_2D_b_02)
      NULLIFY( debug%dp_2D_b_03)
      NULLIFY( debug%dp_2D_b_04)
      NULLIFY( debug%dp_2D_b_05)
      NULLIFY( debug%dp_2D_b_06)
      NULLIFY( debug%dp_2D_b_07)
      NULLIFY( debug%dp_2D_b_08)
      NULLIFY( debug%dp_2D_b_09)
      NULLIFY( debug%dp_2D_b_10)

      NULLIFY( debug%dp_2D_c_01)
      NULLIFY( debug%dp_2D_c_02)
      NULLIFY( debug%dp_2D_c_03)
      NULLIFY( debug%dp_2D_c_04)
      NULLIFY( debug%dp_2D_c_05)
      NULLIFY( debug%dp_2D_c_06)
      NULLIFY( debug%dp_2D_c_07)
      NULLIFY( debug%dp_2D_c_08)
      NULLIFY( debug%dp_2D_c_09)
      NULLIFY( debug%dp_2D_c_10)

      NULLIFY( debug%dp_3D_a_01)
      NULLIFY( debug%dp_3D_a_02)
      NULLIFY( debug%dp_3D_a_03)
      NULLIFY( debug%dp_3D_a_04)
      NULLIFY( debug%dp_3D_a_05)
      NULLIFY( debug%dp_3D_a_06)
      NULLIFY( debug%dp_3D_a_07)
      NULLIFY( debug%dp_3D_a_08)
      NULLIFY( debug%dp_3D_a_09)
      NULLIFY( debug%dp_3D_a_10)

      NULLIFY( debug%dp_2D_monthly_a_01)
      NULLIFY( debug%dp_2D_monthly_a_02)
      NULLIFY( debug%dp_2D_monthly_a_03)
      NULLIFY( debug%dp_2D_monthly_a_04)
      NULLIFY( debug%dp_2D_monthly_a_05)
      NULLIFY( debug%dp_2D_monthly_a_06)
      NULLIFY( debug%dp_2D_monthly_a_07)
      NULLIFY( debug%dp_2D_monthly_a_08)
      NULLIFY( debug%dp_2D_monthly_a_09)
      NULLIFY( debug%dp_2D_monthly_a_10)

    END IF

    ! Bind to the actual memory for this region
    IF (region%name == 'NAM') THEN

      debug%int_2D_a_01 => debug_NAM%int_2D_a_01
      debug%int_2D_a_02 => debug_NAM%int_2D_a_02
      debug%int_2D_a_03 => debug_NAM%int_2D_a_03
      debug%int_2D_a_04 => debug_NAM%int_2D_a_04
      debug%int_2D_a_05 => debug_NAM%int_2D_a_05
      debug%int_2D_a_06 => debug_NAM%int_2D_a_06
      debug%int_2D_a_07 => debug_NAM%int_2D_a_07
      debug%int_2D_a_08 => debug_NAM%int_2D_a_08
      debug%int_2D_a_09 => debug_NAM%int_2D_a_09
      debug%int_2D_a_10 => debug_NAM%int_2D_a_10

      debug%int_2D_b_01 => debug_NAM%int_2D_b_01
      debug%int_2D_b_02 => debug_NAM%int_2D_b_02
      debug%int_2D_b_03 => debug_NAM%int_2D_b_03
      debug%int_2D_b_04 => debug_NAM%int_2D_b_04
      debug%int_2D_b_05 => debug_NAM%int_2D_b_05
      debug%int_2D_b_06 => debug_NAM%int_2D_b_06
      debug%int_2D_b_07 => debug_NAM%int_2D_b_07
      debug%int_2D_b_08 => debug_NAM%int_2D_b_08
      debug%int_2D_b_09 => debug_NAM%int_2D_b_09
      debug%int_2D_b_10 => debug_NAM%int_2D_b_10

      debug%int_2D_c_01 => debug_NAM%int_2D_c_01
      debug%int_2D_c_02 => debug_NAM%int_2D_c_02
      debug%int_2D_c_03 => debug_NAM%int_2D_c_03
      debug%int_2D_c_04 => debug_NAM%int_2D_c_04
      debug%int_2D_c_05 => debug_NAM%int_2D_c_05
      debug%int_2D_c_06 => debug_NAM%int_2D_c_06
      debug%int_2D_c_07 => debug_NAM%int_2D_c_07
      debug%int_2D_c_08 => debug_NAM%int_2D_c_08
      debug%int_2D_c_09 => debug_NAM%int_2D_c_09
      debug%int_2D_c_10 => debug_NAM%int_2D_c_10

      debug%dp_2D_a_01 => debug_NAM%dp_2D_a_01
      debug%dp_2D_a_02 => debug_NAM%dp_2D_a_02
      debug%dp_2D_a_03 => debug_NAM%dp_2D_a_03
      debug%dp_2D_a_04 => debug_NAM%dp_2D_a_04
      debug%dp_2D_a_05 => debug_NAM%dp_2D_a_05
      debug%dp_2D_a_06 => debug_NAM%dp_2D_a_06
      debug%dp_2D_a_07 => debug_NAM%dp_2D_a_07
      debug%dp_2D_a_08 => debug_NAM%dp_2D_a_08
      debug%dp_2D_a_09 => debug_NAM%dp_2D_a_09
      debug%dp_2D_a_10 => debug_NAM%dp_2D_a_10

      debug%dp_2D_b_01 => debug_NAM%dp_2D_b_01
      debug%dp_2D_b_02 => debug_NAM%dp_2D_b_02
      debug%dp_2D_b_03 => debug_NAM%dp_2D_b_03
      debug%dp_2D_b_04 => debug_NAM%dp_2D_b_04
      debug%dp_2D_b_05 => debug_NAM%dp_2D_b_05
      debug%dp_2D_b_06 => debug_NAM%dp_2D_b_06
      debug%dp_2D_b_07 => debug_NAM%dp_2D_b_07
      debug%dp_2D_b_08 => debug_NAM%dp_2D_b_08
      debug%dp_2D_b_09 => debug_NAM%dp_2D_b_09
      debug%dp_2D_b_10 => debug_NAM%dp_2D_b_10

      debug%dp_2D_c_01 => debug_NAM%dp_2D_c_01
      debug%dp_2D_c_02 => debug_NAM%dp_2D_c_02
      debug%dp_2D_c_03 => debug_NAM%dp_2D_c_03
      debug%dp_2D_c_04 => debug_NAM%dp_2D_c_04
      debug%dp_2D_c_05 => debug_NAM%dp_2D_c_05
      debug%dp_2D_c_06 => debug_NAM%dp_2D_c_06
      debug%dp_2D_c_07 => debug_NAM%dp_2D_c_07
      debug%dp_2D_c_08 => debug_NAM%dp_2D_c_08
      debug%dp_2D_c_09 => debug_NAM%dp_2D_c_09
      debug%dp_2D_c_10 => debug_NAM%dp_2D_c_10

      debug%dp_3D_a_01 => debug_NAM%dp_3D_a_01
      debug%dp_3D_a_02 => debug_NAM%dp_3D_a_02
      debug%dp_3D_a_03 => debug_NAM%dp_3D_a_03
      debug%dp_3D_a_04 => debug_NAM%dp_3D_a_04
      debug%dp_3D_a_05 => debug_NAM%dp_3D_a_05
      debug%dp_3D_a_06 => debug_NAM%dp_3D_a_06
      debug%dp_3D_a_07 => debug_NAM%dp_3D_a_07
      debug%dp_3D_a_08 => debug_NAM%dp_3D_a_08
      debug%dp_3D_a_09 => debug_NAM%dp_3D_a_09
      debug%dp_3D_a_10 => debug_NAM%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_NAM%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_NAM%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_NAM%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_NAM%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_NAM%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_NAM%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_NAM%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_NAM%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_NAM%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_NAM%dp_2D_monthly_a_10

    ELSEIF (region%name == 'EAS') THEN

      debug%int_2D_a_01 => debug_EAS%int_2D_a_01
      debug%int_2D_a_02 => debug_EAS%int_2D_a_02
      debug%int_2D_a_03 => debug_EAS%int_2D_a_03
      debug%int_2D_a_04 => debug_EAS%int_2D_a_04
      debug%int_2D_a_05 => debug_EAS%int_2D_a_05
      debug%int_2D_a_06 => debug_EAS%int_2D_a_06
      debug%int_2D_a_07 => debug_EAS%int_2D_a_07
      debug%int_2D_a_08 => debug_EAS%int_2D_a_08
      debug%int_2D_a_09 => debug_EAS%int_2D_a_09
      debug%int_2D_a_10 => debug_EAS%int_2D_a_10

      debug%int_2D_b_01 => debug_EAS%int_2D_b_01
      debug%int_2D_b_02 => debug_EAS%int_2D_b_02
      debug%int_2D_b_03 => debug_EAS%int_2D_b_03
      debug%int_2D_b_04 => debug_EAS%int_2D_b_04
      debug%int_2D_b_05 => debug_EAS%int_2D_b_05
      debug%int_2D_b_06 => debug_EAS%int_2D_b_06
      debug%int_2D_b_07 => debug_EAS%int_2D_b_07
      debug%int_2D_b_08 => debug_EAS%int_2D_b_08
      debug%int_2D_b_09 => debug_EAS%int_2D_b_09
      debug%int_2D_b_10 => debug_EAS%int_2D_b_10

      debug%int_2D_c_01 => debug_EAS%int_2D_c_01
      debug%int_2D_c_02 => debug_EAS%int_2D_c_02
      debug%int_2D_c_03 => debug_EAS%int_2D_c_03
      debug%int_2D_c_04 => debug_EAS%int_2D_c_04
      debug%int_2D_c_05 => debug_EAS%int_2D_c_05
      debug%int_2D_c_06 => debug_EAS%int_2D_c_06
      debug%int_2D_c_07 => debug_EAS%int_2D_c_07
      debug%int_2D_c_08 => debug_EAS%int_2D_c_08
      debug%int_2D_c_09 => debug_EAS%int_2D_c_09
      debug%int_2D_c_10 => debug_EAS%int_2D_c_10

      debug%dp_2D_a_01 => debug_EAS%dp_2D_a_01
      debug%dp_2D_a_02 => debug_EAS%dp_2D_a_02
      debug%dp_2D_a_03 => debug_EAS%dp_2D_a_03
      debug%dp_2D_a_04 => debug_EAS%dp_2D_a_04
      debug%dp_2D_a_05 => debug_EAS%dp_2D_a_05
      debug%dp_2D_a_06 => debug_EAS%dp_2D_a_06
      debug%dp_2D_a_07 => debug_EAS%dp_2D_a_07
      debug%dp_2D_a_08 => debug_EAS%dp_2D_a_08
      debug%dp_2D_a_09 => debug_EAS%dp_2D_a_09
      debug%dp_2D_a_10 => debug_EAS%dp_2D_a_10

      debug%dp_2D_b_01 => debug_EAS%dp_2D_b_01
      debug%dp_2D_b_02 => debug_EAS%dp_2D_b_02
      debug%dp_2D_b_03 => debug_EAS%dp_2D_b_03
      debug%dp_2D_b_04 => debug_EAS%dp_2D_b_04
      debug%dp_2D_b_05 => debug_EAS%dp_2D_b_05
      debug%dp_2D_b_06 => debug_EAS%dp_2D_b_06
      debug%dp_2D_b_07 => debug_EAS%dp_2D_b_07
      debug%dp_2D_b_08 => debug_EAS%dp_2D_b_08
      debug%dp_2D_b_09 => debug_EAS%dp_2D_b_09
      debug%dp_2D_b_10 => debug_EAS%dp_2D_b_10

      debug%dp_2D_c_01 => debug_EAS%dp_2D_c_01
      debug%dp_2D_c_02 => debug_EAS%dp_2D_c_02
      debug%dp_2D_c_03 => debug_EAS%dp_2D_c_03
      debug%dp_2D_c_04 => debug_EAS%dp_2D_c_04
      debug%dp_2D_c_05 => debug_EAS%dp_2D_c_05
      debug%dp_2D_c_06 => debug_EAS%dp_2D_c_06
      debug%dp_2D_c_07 => debug_EAS%dp_2D_c_07
      debug%dp_2D_c_08 => debug_EAS%dp_2D_c_08
      debug%dp_2D_c_09 => debug_EAS%dp_2D_c_09
      debug%dp_2D_c_10 => debug_EAS%dp_2D_c_10

      debug%dp_3D_a_01 => debug_EAS%dp_3D_a_01
      debug%dp_3D_a_02 => debug_EAS%dp_3D_a_02
      debug%dp_3D_a_03 => debug_EAS%dp_3D_a_03
      debug%dp_3D_a_04 => debug_EAS%dp_3D_a_04
      debug%dp_3D_a_05 => debug_EAS%dp_3D_a_05
      debug%dp_3D_a_06 => debug_EAS%dp_3D_a_06
      debug%dp_3D_a_07 => debug_EAS%dp_3D_a_07
      debug%dp_3D_a_08 => debug_EAS%dp_3D_a_08
      debug%dp_3D_a_09 => debug_EAS%dp_3D_a_09
      debug%dp_3D_a_10 => debug_EAS%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_EAS%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_EAS%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_EAS%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_EAS%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_EAS%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_EAS%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_EAS%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_EAS%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_EAS%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_EAS%dp_2D_monthly_a_10

    ELSEIF (region%name == 'GRL') THEN

      debug%int_2D_a_01 => debug_GRL%int_2D_a_01
      debug%int_2D_a_02 => debug_GRL%int_2D_a_02
      debug%int_2D_a_03 => debug_GRL%int_2D_a_03
      debug%int_2D_a_04 => debug_GRL%int_2D_a_04
      debug%int_2D_a_05 => debug_GRL%int_2D_a_05
      debug%int_2D_a_06 => debug_GRL%int_2D_a_06
      debug%int_2D_a_07 => debug_GRL%int_2D_a_07
      debug%int_2D_a_08 => debug_GRL%int_2D_a_08
      debug%int_2D_a_09 => debug_GRL%int_2D_a_09
      debug%int_2D_a_10 => debug_GRL%int_2D_a_10

      debug%int_2D_b_01 => debug_GRL%int_2D_b_01
      debug%int_2D_b_02 => debug_GRL%int_2D_b_02
      debug%int_2D_b_03 => debug_GRL%int_2D_b_03
      debug%int_2D_b_04 => debug_GRL%int_2D_b_04
      debug%int_2D_b_05 => debug_GRL%int_2D_b_05
      debug%int_2D_b_06 => debug_GRL%int_2D_b_06
      debug%int_2D_b_07 => debug_GRL%int_2D_b_07
      debug%int_2D_b_08 => debug_GRL%int_2D_b_08
      debug%int_2D_b_09 => debug_GRL%int_2D_b_09
      debug%int_2D_b_10 => debug_GRL%int_2D_b_10

      debug%int_2D_c_01 => debug_GRL%int_2D_c_01
      debug%int_2D_c_02 => debug_GRL%int_2D_c_02
      debug%int_2D_c_03 => debug_GRL%int_2D_c_03
      debug%int_2D_c_04 => debug_GRL%int_2D_c_04
      debug%int_2D_c_05 => debug_GRL%int_2D_c_05
      debug%int_2D_c_06 => debug_GRL%int_2D_c_06
      debug%int_2D_c_07 => debug_GRL%int_2D_c_07
      debug%int_2D_c_08 => debug_GRL%int_2D_c_08
      debug%int_2D_c_09 => debug_GRL%int_2D_c_09
      debug%int_2D_c_10 => debug_GRL%int_2D_c_10

      debug%dp_2D_a_01 => debug_GRL%dp_2D_a_01
      debug%dp_2D_a_02 => debug_GRL%dp_2D_a_02
      debug%dp_2D_a_03 => debug_GRL%dp_2D_a_03
      debug%dp_2D_a_04 => debug_GRL%dp_2D_a_04
      debug%dp_2D_a_05 => debug_GRL%dp_2D_a_05
      debug%dp_2D_a_06 => debug_GRL%dp_2D_a_06
      debug%dp_2D_a_07 => debug_GRL%dp_2D_a_07
      debug%dp_2D_a_08 => debug_GRL%dp_2D_a_08
      debug%dp_2D_a_09 => debug_GRL%dp_2D_a_09
      debug%dp_2D_a_10 => debug_GRL%dp_2D_a_10

      debug%dp_2D_b_01 => debug_GRL%dp_2D_b_01
      debug%dp_2D_b_02 => debug_GRL%dp_2D_b_02
      debug%dp_2D_b_03 => debug_GRL%dp_2D_b_03
      debug%dp_2D_b_04 => debug_GRL%dp_2D_b_04
      debug%dp_2D_b_05 => debug_GRL%dp_2D_b_05
      debug%dp_2D_b_06 => debug_GRL%dp_2D_b_06
      debug%dp_2D_b_07 => debug_GRL%dp_2D_b_07
      debug%dp_2D_b_08 => debug_GRL%dp_2D_b_08
      debug%dp_2D_b_09 => debug_GRL%dp_2D_b_09
      debug%dp_2D_b_10 => debug_GRL%dp_2D_b_10

      debug%dp_2D_c_01 => debug_GRL%dp_2D_c_01
      debug%dp_2D_c_02 => debug_GRL%dp_2D_c_02
      debug%dp_2D_c_03 => debug_GRL%dp_2D_c_03
      debug%dp_2D_c_04 => debug_GRL%dp_2D_c_04
      debug%dp_2D_c_05 => debug_GRL%dp_2D_c_05
      debug%dp_2D_c_06 => debug_GRL%dp_2D_c_06
      debug%dp_2D_c_07 => debug_GRL%dp_2D_c_07
      debug%dp_2D_c_08 => debug_GRL%dp_2D_c_08
      debug%dp_2D_c_09 => debug_GRL%dp_2D_c_09
      debug%dp_2D_c_10 => debug_GRL%dp_2D_c_10

      debug%dp_3D_a_01 => debug_GRL%dp_3D_a_01
      debug%dp_3D_a_02 => debug_GRL%dp_3D_a_02
      debug%dp_3D_a_03 => debug_GRL%dp_3D_a_03
      debug%dp_3D_a_04 => debug_GRL%dp_3D_a_04
      debug%dp_3D_a_05 => debug_GRL%dp_3D_a_05
      debug%dp_3D_a_06 => debug_GRL%dp_3D_a_06
      debug%dp_3D_a_07 => debug_GRL%dp_3D_a_07
      debug%dp_3D_a_08 => debug_GRL%dp_3D_a_08
      debug%dp_3D_a_09 => debug_GRL%dp_3D_a_09
      debug%dp_3D_a_10 => debug_GRL%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_GRL%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_GRL%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_GRL%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_GRL%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_GRL%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_GRL%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_GRL%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_GRL%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_GRL%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_GRL%dp_2D_monthly_a_10

    ELSEIF (region%name == 'ANT') THEN

      debug%int_2D_a_01 => debug_ANT%int_2D_a_01
      debug%int_2D_a_02 => debug_ANT%int_2D_a_02
      debug%int_2D_a_03 => debug_ANT%int_2D_a_03
      debug%int_2D_a_04 => debug_ANT%int_2D_a_04
      debug%int_2D_a_05 => debug_ANT%int_2D_a_05
      debug%int_2D_a_06 => debug_ANT%int_2D_a_06
      debug%int_2D_a_07 => debug_ANT%int_2D_a_07
      debug%int_2D_a_08 => debug_ANT%int_2D_a_08
      debug%int_2D_a_09 => debug_ANT%int_2D_a_09
      debug%int_2D_a_10 => debug_ANT%int_2D_a_10

      debug%int_2D_b_01 => debug_ANT%int_2D_b_01
      debug%int_2D_b_02 => debug_ANT%int_2D_b_02
      debug%int_2D_b_03 => debug_ANT%int_2D_b_03
      debug%int_2D_b_04 => debug_ANT%int_2D_b_04
      debug%int_2D_b_05 => debug_ANT%int_2D_b_05
      debug%int_2D_b_06 => debug_ANT%int_2D_b_06
      debug%int_2D_b_07 => debug_ANT%int_2D_b_07
      debug%int_2D_b_08 => debug_ANT%int_2D_b_08
      debug%int_2D_b_09 => debug_ANT%int_2D_b_09
      debug%int_2D_b_10 => debug_ANT%int_2D_b_10

      debug%int_2D_c_01 => debug_ANT%int_2D_c_01
      debug%int_2D_c_02 => debug_ANT%int_2D_c_02
      debug%int_2D_c_03 => debug_ANT%int_2D_c_03
      debug%int_2D_c_04 => debug_ANT%int_2D_c_04
      debug%int_2D_c_05 => debug_ANT%int_2D_c_05
      debug%int_2D_c_06 => debug_ANT%int_2D_c_06
      debug%int_2D_c_07 => debug_ANT%int_2D_c_07
      debug%int_2D_c_08 => debug_ANT%int_2D_c_08
      debug%int_2D_c_09 => debug_ANT%int_2D_c_09
      debug%int_2D_c_10 => debug_ANT%int_2D_c_10

      debug%dp_2D_a_01 => debug_ANT%dp_2D_a_01
      debug%dp_2D_a_02 => debug_ANT%dp_2D_a_02
      debug%dp_2D_a_03 => debug_ANT%dp_2D_a_03
      debug%dp_2D_a_04 => debug_ANT%dp_2D_a_04
      debug%dp_2D_a_05 => debug_ANT%dp_2D_a_05
      debug%dp_2D_a_06 => debug_ANT%dp_2D_a_06
      debug%dp_2D_a_07 => debug_ANT%dp_2D_a_07
      debug%dp_2D_a_08 => debug_ANT%dp_2D_a_08
      debug%dp_2D_a_09 => debug_ANT%dp_2D_a_09
      debug%dp_2D_a_10 => debug_ANT%dp_2D_a_10

      debug%dp_2D_b_01 => debug_ANT%dp_2D_b_01
      debug%dp_2D_b_02 => debug_ANT%dp_2D_b_02
      debug%dp_2D_b_03 => debug_ANT%dp_2D_b_03
      debug%dp_2D_b_04 => debug_ANT%dp_2D_b_04
      debug%dp_2D_b_05 => debug_ANT%dp_2D_b_05
      debug%dp_2D_b_06 => debug_ANT%dp_2D_b_06
      debug%dp_2D_b_07 => debug_ANT%dp_2D_b_07
      debug%dp_2D_b_08 => debug_ANT%dp_2D_b_08
      debug%dp_2D_b_09 => debug_ANT%dp_2D_b_09
      debug%dp_2D_b_10 => debug_ANT%dp_2D_b_10

      debug%dp_2D_c_01 => debug_ANT%dp_2D_c_01
      debug%dp_2D_c_02 => debug_ANT%dp_2D_c_02
      debug%dp_2D_c_03 => debug_ANT%dp_2D_c_03
      debug%dp_2D_c_04 => debug_ANT%dp_2D_c_04
      debug%dp_2D_c_05 => debug_ANT%dp_2D_c_05
      debug%dp_2D_c_06 => debug_ANT%dp_2D_c_06
      debug%dp_2D_c_07 => debug_ANT%dp_2D_c_07
      debug%dp_2D_c_08 => debug_ANT%dp_2D_c_08
      debug%dp_2D_c_09 => debug_ANT%dp_2D_c_09
      debug%dp_2D_c_10 => debug_ANT%dp_2D_c_10

      debug%dp_3D_a_01 => debug_ANT%dp_3D_a_01
      debug%dp_3D_a_02 => debug_ANT%dp_3D_a_02
      debug%dp_3D_a_03 => debug_ANT%dp_3D_a_03
      debug%dp_3D_a_04 => debug_ANT%dp_3D_a_04
      debug%dp_3D_a_05 => debug_ANT%dp_3D_a_05
      debug%dp_3D_a_06 => debug_ANT%dp_3D_a_06
      debug%dp_3D_a_07 => debug_ANT%dp_3D_a_07
      debug%dp_3D_a_08 => debug_ANT%dp_3D_a_08
      debug%dp_3D_a_09 => debug_ANT%dp_3D_a_09
      debug%dp_3D_a_10 => debug_ANT%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_ANT%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_ANT%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_ANT%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_ANT%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_ANT%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_ANT%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_ANT%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_ANT%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_ANT%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_ANT%dp_2D_monthly_a_10

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE associate_debug_fields

  SUBROUTINE initialise_debug_fields( region)
    ! Allocate all the fields of the debug structure

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (region%name == 'NAM') THEN
      CALL initialise_debug_fields_region( debug_NAM, region%mesh)
    ELSEIF (region%name == 'EAS') THEN
      CALL initialise_debug_fields_region( debug_EAS, region%mesh)
    ELSEIF (region%name == 'GRL') THEN
      CALL initialise_debug_fields_region( debug_GRL, region%mesh)
    ELSEIF (region%name == 'ANT') THEN
      CALL initialise_debug_fields_region( debug_ANT, region%mesh)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 80)

  END SUBROUTINE initialise_debug_fields

  SUBROUTINE reallocate_debug_fields( region)
    ! Deallocate and allocate all the fields of the debug structure

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'reallocate_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (region%name == 'NAM') THEN
      CALL deallocate_debug_fields_region( debug_NAM)
      CALL initialise_debug_fields_region( debug_NAM, region%mesh)
    ELSEIF (region%name == 'EAS') THEN
      CALL deallocate_debug_fields_region( debug_EAS)
      CALL initialise_debug_fields_region( debug_EAS, region%mesh)
    ELSEIF (region%name == 'GRL') THEN
      CALL deallocate_debug_fields_region( debug_GRL)
      CALL initialise_debug_fields_region( debug_GRL, region%mesh)
    ELSEIF (region%name == 'ANT') THEN
      CALL deallocate_debug_fields_region( debug_ANT)
      CALL initialise_debug_fields_region( debug_ANT, region%mesh)
    END IF
    CALL associate_debug_fields( region)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE reallocate_debug_fields

  SUBROUTINE initialise_debug_fields_region( debug, mesh)
    ! Allocate all the fields of the debug structure

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug
    TYPE(type_mesh),                 INTENT(IN)        :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_debug_fields_region'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_01, debug%wint_2D_a_01)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_02, debug%wint_2D_a_02)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_03, debug%wint_2D_a_03)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_04, debug%wint_2D_a_04)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_05, debug%wint_2D_a_05)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_06, debug%wint_2D_a_06)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_07, debug%wint_2D_a_07)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_08, debug%wint_2D_a_08)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_09, debug%wint_2D_a_09)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_10, debug%wint_2D_a_10)

    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_01, debug%wint_2D_b_01)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_02, debug%wint_2D_b_02)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_03, debug%wint_2D_b_03)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_04, debug%wint_2D_b_04)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_05, debug%wint_2D_b_05)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_06, debug%wint_2D_b_06)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_07, debug%wint_2D_b_07)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_08, debug%wint_2D_b_08)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_09, debug%wint_2D_b_09)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_10, debug%wint_2D_b_10)

    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_01, debug%wint_2D_c_01)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_02, debug%wint_2D_c_02)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_03, debug%wint_2D_c_03)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_04, debug%wint_2D_c_04)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_05, debug%wint_2D_c_05)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_06, debug%wint_2D_c_06)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_07, debug%wint_2D_c_07)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_08, debug%wint_2D_c_08)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_09, debug%wint_2D_c_09)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_10, debug%wint_2D_c_10)

    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_01, debug%wdp_2D_a_01)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_02, debug%wdp_2D_a_02)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_03, debug%wdp_2D_a_03)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_04, debug%wdp_2D_a_04)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_05, debug%wdp_2D_a_05)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_06, debug%wdp_2D_a_06)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_07, debug%wdp_2D_a_07)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_08, debug%wdp_2D_a_08)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_09, debug%wdp_2D_a_09)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_10, debug%wdp_2D_a_10)

    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_01, debug%wdp_2D_b_01)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_02, debug%wdp_2D_b_02)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_03, debug%wdp_2D_b_03)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_04, debug%wdp_2D_b_04)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_05, debug%wdp_2D_b_05)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_06, debug%wdp_2D_b_06)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_07, debug%wdp_2D_b_07)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_08, debug%wdp_2D_b_08)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_09, debug%wdp_2D_b_09)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_10, debug%wdp_2D_b_10)

    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_01, debug%wdp_2D_c_01)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_02, debug%wdp_2D_c_02)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_03, debug%wdp_2D_c_03)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_04, debug%wdp_2D_c_04)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_05, debug%wdp_2D_c_05)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_06, debug%wdp_2D_c_06)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_07, debug%wdp_2D_c_07)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_08, debug%wdp_2D_c_08)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_09, debug%wdp_2D_c_09)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_10, debug%wdp_2D_c_10)

    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_01, debug%wdp_3D_a_01)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_02, debug%wdp_3D_a_02)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_03, debug%wdp_3D_a_03)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_04, debug%wdp_3D_a_04)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_05, debug%wdp_3D_a_05)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_06, debug%wdp_3D_a_06)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_07, debug%wdp_3D_a_07)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_08, debug%wdp_3D_a_08)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_09, debug%wdp_3D_a_09)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_10, debug%wdp_3D_a_10)

    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_01, debug%wdp_2D_monthly_a_01)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_02, debug%wdp_2D_monthly_a_02)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_03, debug%wdp_2D_monthly_a_03)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_04, debug%wdp_2D_monthly_a_04)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_05, debug%wdp_2D_monthly_a_05)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_06, debug%wdp_2D_monthly_a_06)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_07, debug%wdp_2D_monthly_a_07)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_08, debug%wdp_2D_monthly_a_08)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_09, debug%wdp_2D_monthly_a_09)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_10, debug%wdp_2D_monthly_a_10)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 80)

  END SUBROUTINE initialise_debug_fields_region

  SUBROUTINE deallocate_debug_fields_region( debug)
    ! Deallocate all the fields of the debug structure

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_debug_fields_region'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL deallocate_shared( debug%wint_2D_a_01)
    CALL deallocate_shared( debug%wint_2D_a_02)
    CALL deallocate_shared( debug%wint_2D_a_03)
    CALL deallocate_shared( debug%wint_2D_a_04)
    CALL deallocate_shared( debug%wint_2D_a_05)
    CALL deallocate_shared( debug%wint_2D_a_06)
    CALL deallocate_shared( debug%wint_2D_a_07)
    CALL deallocate_shared( debug%wint_2D_a_08)
    CALL deallocate_shared( debug%wint_2D_a_09)
    CALL deallocate_shared( debug%wint_2D_a_10)

    CALL deallocate_shared( debug%wint_2D_b_01)
    CALL deallocate_shared( debug%wint_2D_b_02)
    CALL deallocate_shared( debug%wint_2D_b_03)
    CALL deallocate_shared( debug%wint_2D_b_04)
    CALL deallocate_shared( debug%wint_2D_b_05)
    CALL deallocate_shared( debug%wint_2D_b_06)
    CALL deallocate_shared( debug%wint_2D_b_07)
    CALL deallocate_shared( debug%wint_2D_b_08)
    CALL deallocate_shared( debug%wint_2D_b_09)
    CALL deallocate_shared( debug%wint_2D_b_10)

    CALL deallocate_shared( debug%wint_2D_c_01)
    CALL deallocate_shared( debug%wint_2D_c_02)
    CALL deallocate_shared( debug%wint_2D_c_03)
    CALL deallocate_shared( debug%wint_2D_c_04)
    CALL deallocate_shared( debug%wint_2D_c_05)
    CALL deallocate_shared( debug%wint_2D_c_06)
    CALL deallocate_shared( debug%wint_2D_c_07)
    CALL deallocate_shared( debug%wint_2D_c_08)
    CALL deallocate_shared( debug%wint_2D_c_09)
    CALL deallocate_shared( debug%wint_2D_c_10)

    CALL deallocate_shared( debug%wdp_2D_a_01)
    CALL deallocate_shared( debug%wdp_2D_a_02)
    CALL deallocate_shared( debug%wdp_2D_a_03)
    CALL deallocate_shared( debug%wdp_2D_a_04)
    CALL deallocate_shared( debug%wdp_2D_a_05)
    CALL deallocate_shared( debug%wdp_2D_a_06)
    CALL deallocate_shared( debug%wdp_2D_a_07)
    CALL deallocate_shared( debug%wdp_2D_a_08)
    CALL deallocate_shared( debug%wdp_2D_a_09)
    CALL deallocate_shared( debug%wdp_2D_a_10)

    CALL deallocate_shared( debug%wdp_2D_b_01)
    CALL deallocate_shared( debug%wdp_2D_b_02)
    CALL deallocate_shared( debug%wdp_2D_b_03)
    CALL deallocate_shared( debug%wdp_2D_b_04)
    CALL deallocate_shared( debug%wdp_2D_b_05)
    CALL deallocate_shared( debug%wdp_2D_b_06)
    CALL deallocate_shared( debug%wdp_2D_b_07)
    CALL deallocate_shared( debug%wdp_2D_b_08)
    CALL deallocate_shared( debug%wdp_2D_b_09)
    CALL deallocate_shared( debug%wdp_2D_b_10)

    CALL deallocate_shared( debug%wdp_2D_c_01)
    CALL deallocate_shared( debug%wdp_2D_c_02)
    CALL deallocate_shared( debug%wdp_2D_c_03)
    CALL deallocate_shared( debug%wdp_2D_c_04)
    CALL deallocate_shared( debug%wdp_2D_c_05)
    CALL deallocate_shared( debug%wdp_2D_c_06)
    CALL deallocate_shared( debug%wdp_2D_c_07)
    CALL deallocate_shared( debug%wdp_2D_c_08)
    CALL deallocate_shared( debug%wdp_2D_c_09)
    CALL deallocate_shared( debug%wdp_2D_c_10)

    CALL deallocate_shared( debug%wdp_3D_a_01)
    CALL deallocate_shared( debug%wdp_3D_a_02)
    CALL deallocate_shared( debug%wdp_3D_a_03)
    CALL deallocate_shared( debug%wdp_3D_a_04)
    CALL deallocate_shared( debug%wdp_3D_a_05)
    CALL deallocate_shared( debug%wdp_3D_a_06)
    CALL deallocate_shared( debug%wdp_3D_a_07)
    CALL deallocate_shared( debug%wdp_3D_a_08)
    CALL deallocate_shared( debug%wdp_3D_a_09)
    CALL deallocate_shared( debug%wdp_3D_a_10)

    CALL deallocate_shared( debug%wdp_2D_monthly_a_01)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_02)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_03)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_04)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_05)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_06)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_07)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_08)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_09)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_10)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_debug_fields_region

! ===== Single-variable NetCDF files =====
! ========================================

  SUBROUTINE save_variable_as_netcdf_int_1D( d, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:    ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_1D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_n1
    INTEGER                                            :: id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine file name
    filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', SIZE( d,1), id_dim_n1)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_INT, (/ id_dim_n1 /), id_var)

    ! Write data
    CALL write_var_int_1D(  filename, ncid, id_var, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_1D

  SUBROUTINE save_variable_as_netcdf_int_2D( d, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:,:  ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_2D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_n1, id_dim_n2
    INTEGER                                            :: id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine file name
    filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', SIZE( d,1), id_dim_n1)
    CALL create_dimension( filename, ncid, 'n2', SIZE( d,2), id_dim_n2)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_INT, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    CALL write_var_int_2D(  filename, ncid, id_var, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_2D

  SUBROUTINE save_variable_as_netcdf_int_3D( d, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:,:,:),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_3D'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1, id_dim_n2, id_dim_n3
    INTEGER                                            :: id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine file name
    filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', SIZE( d,1), id_dim_n1)
    CALL create_dimension( filename, ncid, 'n2', SIZE( d,2), id_dim_n2)
    CALL create_dimension( filename, ncid, 'n3', SIZE( d,3), id_dim_n3)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_INT, (/ id_dim_n1, id_dim_n2, id_dim_n3 /), id_var)

    ! Write data
    CALL write_var_int_3D(  filename, ncid, id_var, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_3D

  SUBROUTINE save_variable_as_netcdf_dp_1D( d, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:    ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_1D'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1
    INTEGER                                            :: id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine file name
    filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', SIZE( d,1), id_dim_n1)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_DOUBLE, (/ id_dim_n1 /), id_var)

    ! Write data
    CALL write_var_dp_1D(  filename, ncid, id_var, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_1D

  SUBROUTINE save_variable_as_netcdf_dp_2D( d, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:,:  ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_2D'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1, id_dim_n2
    INTEGER                                            :: id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine file name
    filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', SIZE( d,1), id_dim_n1)
    CALL create_dimension( filename, ncid, 'n2', SIZE( d,2), id_dim_n2)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_DOUBLE, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    CALL write_var_dp_2D(  filename, ncid, id_var, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_2D

  SUBROUTINE save_variable_as_netcdf_dp_3D( d, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:,:,:),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_3D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_n1, id_dim_n2, id_dim_n3
    INTEGER                                            :: id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine file name
    filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', SIZE( d,1), id_dim_n1)
    CALL create_dimension( filename, ncid, 'n2', SIZE( d,2), id_dim_n2)
    CALL create_dimension( filename, ncid, 'n3', SIZE( d,3), id_dim_n3)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_DOUBLE, (/ id_dim_n1, id_dim_n2, id_dim_n3 /), id_var)

    ! Write data
    CALL write_var_dp_3D(  filename, ncid, id_var, d)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_3D

! ===== Matrix NetCDF files =====
! ===============================

  SUBROUTINE write_PETSc_matrix_to_NetCDF( A, filename)
    ! Write a PETSc matrix to a NetCDF file

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(tMat),                          INTENT(IN)    :: A
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_PETSc_matrix_to_NetCDF'
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_CSR

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get matrix in CSR format using native Fortran arrays
    CALL mat_petsc2CSR( A, A_CSR)

    ! Write the CSR matrix to a file
    CALL write_CSR_matrix_to_NetCDF( A_CSR, filename)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR( A_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_PETSc_matrix_to_NetCDF

  SUBROUTINE write_CSR_matrix_to_NetCDF( A_CSR, filename)
    ! Write a CSR matrix to a NetCDF file

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: A_CSR
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_CSR_matrix_to_NetCDF'
    CHARACTER(LEN=256)                                 :: filename_applied
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz
    INTEGER                                            :: id_var_ptr, id_var_index, id_var_val

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Append output directory to filename
    filename_applied = TRIM( C%output_dir) // TRIM( filename) // '.nc'

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename_applied, ncid)

    ! Create dimensions
    CALL create_dimension( filename_applied, ncid, 'm'     , A_CSR%m  , id_dim_m  )
    CALL create_dimension( filename_applied, ncid, 'mplus1', A_CSR%m+1, id_dim_mp1)
    CALL create_dimension( filename_applied, ncid, 'n'     , A_CSR%n  , id_dim_n  )
    CALL create_dimension( filename_applied, ncid, 'nnz'   , A_CSR%nnz, id_dim_nnz)

    ! Create variables
    CALL create_variable( filename_applied, ncid, 'ptr'  , NF90_INT   , [id_dim_mp1], id_var_ptr  )
    CALL create_variable( filename_applied, ncid, 'index', NF90_INT   , [id_dim_nnz], id_var_index)
    CALL create_variable( filename_applied, ncid, 'val'  , NF90_DOUBLE, [id_dim_nnz], id_var_val  )

    ! Write data
    CALL write_var_int_1D( filename_applied, ncid, id_var_ptr  , A_CSR%ptr  )
    CALL write_var_int_1D( filename_applied, ncid, id_var_index, A_CSR%index( 1:A_CSR%nnz))
    CALL write_var_dp_1D(  filename_applied, ncid, id_var_val  , A_CSR%val(   1:A_CSR%nnz))

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_CSR_matrix_to_NetCDF

END MODULE netcdf_debug_module