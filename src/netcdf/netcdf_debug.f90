MODULE netcdf_debug

! Contains routines for quickly writing variables to NetCDF to check their content.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, deallocate_matrix_CSR_dist, gather_CSR_dist_to_master
  USE petsc_basic                                            , ONLY: mat_petsc2CSR
  use mpi_distributed_memory, only: gather_to_master

  USE netcdf,       ONLY: NF90_UNLIMITED, NF90_INT, NF90_FLOAT, NF90_DOUBLE
  USE netcdf_basic, ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, &
                          field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                          field_name_options_dim_nV, field_name_options_dim_nTri, field_name_options_dim_nC_mem, &
                          field_name_options_dim_nE, field_name_options_dim_two, field_name_options_dim_three, &
                          field_name_options_dim_four, field_name_options_V, field_name_options_Tri, field_name_options_nC, &
                          field_name_options_C, field_name_options_niTri, field_name_options_iTri, &
                          field_name_options_VBI, field_name_options_Tricc, field_name_options_TriC, &
                          field_name_options_TriBI, field_name_options_E, field_name_options_VE, field_name_options_EV, &
                          field_name_options_ETri, field_name_options_EBI, field_name_options_A, field_name_options_R, &
                          field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                          field_name_options_SL, field_name_options_Ti, get_first_option_from_list, &
                          open_existing_netcdf_file_for_reading, close_netcdf_file, create_new_netcdf_file_for_writing, &
                          inquire_dim_multopt, inquire_var_multopt, &
                          check_x, check_y, check_lon, check_lat, check_mesh_dimensions, check_zeta, find_timeframe, &
                          check_xy_grid_field_int_2D, check_xy_grid_field_dp_2D, check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, &
                          check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, check_lonlat_grid_field_dp_2D_monthly, check_lonlat_grid_field_dp_3D, &
                          check_mesh_field_int_2D, check_mesh_field_int_2D_b, check_mesh_field_int_2D_c, &
                          check_mesh_field_dp_2D, check_mesh_field_dp_2D_b, check_mesh_field_dp_2D_c, &
                          check_mesh_field_dp_2D_monthly, check_mesh_field_dp_3D, &
                          inquire_xy_grid, inquire_lonlat_grid, inquire_mesh, &
                          write_var_master_int_0D, write_var_master_int_1D, write_var_master_int_2D, write_var_master_int_3D, write_var_master_int_4D, &
                          write_var_master_dp_0D, write_var_master_dp_1D, write_var_master_dp_2D, write_var_master_dp_3D, write_var_master_dp_4D, &
                          add_attribute_char, check_month, check_time, create_dimension, create_variable, inquire_var, &
                          open_existing_netcdf_file_for_writing

  IMPLICIT NONE

CONTAINS

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
    TYPE(type_sparse_matrix_CSR_dp)                    :: AA

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get matrix in CSR format using native Fortran arrays
    CALL mat_petsc2CSR( A, AA)

    ! Write the CSR matrix to a file
    CALL write_CSR_matrix_to_NetCDF( AA, filename)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( AA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_PETSc_matrix_to_NetCDF

  SUBROUTINE write_CSR_matrix_to_NetCDF( AA, filename)
    ! Write a CSR matrix to a NetCDF file

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_CSR_matrix_to_NetCDF'
    CHARACTER(LEN=256)                                 :: filename_applied
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz
    INTEGER                                            :: id_var_ptr, id_var_ind, id_var_val
    TYPE(type_sparse_matrix_CSR_dp)                    :: AA_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather distributed matrix to the master
    CALL gather_CSR_dist_to_master( AA, AA_tot)

    ! Append output directory to filename
    filename_applied = TRIM( C%output_dir) // TRIM( filename) // '.nc'

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename_applied))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename_applied)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename_applied, ncid)

    ! Create dimensions
    CALL create_dimension( filename_applied, ncid, 'm'     , AA_tot%m  , id_dim_m  )
    CALL create_dimension( filename_applied, ncid, 'mplus1', AA_tot%m+1, id_dim_mp1)
    CALL create_dimension( filename_applied, ncid, 'n'     , AA_tot%n  , id_dim_n  )
    CALL create_dimension( filename_applied, ncid, 'nnz'   , AA_tot%nnz, id_dim_nnz)

    ! Create variables
    CALL create_variable( filename_applied, ncid, 'ptr', NF90_INT   , [id_dim_mp1], id_var_ptr)
    CALL create_variable( filename_applied, ncid, 'ind', NF90_INT   , [id_dim_nnz], id_var_ind)
    CALL create_variable( filename_applied, ncid, 'val', NF90_DOUBLE, [id_dim_nnz], id_var_val)

    ! Write to NetCDF
    CALL write_var_master_int_1D( filename_applied, ncid, id_var_ptr, AA_tot%ptr               )
    CALL write_var_master_int_1D( filename_applied, ncid, id_var_ind, AA_tot%ind( 1:AA_tot%nnz))
    CALL write_var_master_dp_1D(  filename_applied, ncid, id_var_val, AA_tot%val( 1:AA_tot%nnz))

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    IF (par%master) THEN
      CALL deallocate_matrix_CSR_dist( AA_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_CSR_matrix_to_NetCDF

! ===== Single-variable NetCDF files =====
! ========================================

  SUBROUTINE save_variable_as_netcdf_logical_1D( d_partial, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    LOGICAL,  DIMENSION(:    ),      INTENT(IN)        :: d_partial
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_logical_1D'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: d_partial_int

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for integer version
    ALLOCATE( d_partial_int( SIZE( d_partial,1)))

    ! Fill integer version
    WHERE (d_partial)
      d_partial_int = 1
    ELSEWHERE
      d_partial_int = 0
    END WHERE

    ! Write integer version to NetCDF
    CALL save_variable_as_netcdf_int_1D( d_partial_int, field_name)

    ! Clean up after yourself
    DEALLOCATE( d_partial_int)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_logical_1D

  SUBROUTINE save_variable_as_netcdf_int_1D( d_partial, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:    ),      INTENT(IN)        :: d_partial
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_1D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: n_partial, n_tot
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: d_tot
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

    ! Gather data to the master
    n_partial = SIZE( d_partial,1)
    CALL MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot))
      CALL gather_to_master( d_partial, d_tot)
    ELSE
      CALL gather_to_master( d_partial)
    END IF

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_INT, (/ id_dim_n1 /), id_var)

    ! Write data
    CALL write_var_master_int_1D( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_1D

  SUBROUTINE save_variable_as_netcdf_int_2D( d_partial, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:,:  ),      INTENT(IN)        :: d_partial
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_2D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: n_partial, n_tot, n2
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: d_tot
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

    ! Gather data to the master
    n_partial = SIZE( d_partial,1)
    n2        = SIZE( d_partial,2)
    CALL MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2))
      CALL gather_to_master( d_partial, d_tot)
    ELSE
      CALL gather_to_master( d_partial)
    END IF

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)
    CALL create_dimension( filename, ncid, 'n2', n_tot, id_dim_n2)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_INT, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    CALL write_var_master_int_2D( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_2D

  SUBROUTINE save_variable_as_netcdf_dp_1D( d_partial, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:    ),      INTENT(IN)        :: d_partial
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_1D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: n_partial, n_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_tot
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

    ! Gather data to the master
    n_partial = SIZE( d_partial,1)
    CALL MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot))
      CALL gather_to_master( d_partial, d_tot)
    ELSE
      CALL gather_to_master( d_partial)
    END IF

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_DOUBLE, (/ id_dim_n1 /), id_var)

    ! Write data
    CALL write_var_master_dp_1D( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_1D

  SUBROUTINE save_variable_as_netcdf_dp_2D( d_partial, field_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:,:  ),      INTENT(IN)        :: d_partial
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_2D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: n_partial, n_tot, n2
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot
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

    ! Gather data to the master
    n_partial = SIZE( d_partial,1)
    n2        = SIZE( d_partial,2)
    CALL MPI_ALLREDUCE( n_partial, n_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (par%master) THEN
      ALLOCATE( d_tot( n_tot, n2))
      CALL gather_to_master( d_partial, d_tot)
    ELSE
      CALL gather_to_master( d_partial)
    END IF

    ! Create dimensions
    CALL create_dimension( filename, ncid, 'n1', n_tot, id_dim_n1)
    CALL create_dimension( filename, ncid, 'n2', n2   , id_dim_n2)

    ! Create variable
    CALL create_variable( filename, ncid, field_name, NF90_DOUBLE, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    CALL write_var_master_dp_2D( filename, ncid, id_var, d_tot)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_2D

END MODULE netcdf_debug