MODULE netcdf_output

! ===== Creating and writing to output files =====
! ================================================
!
! These routines create and write data to output files, both on meshes and on grids. For
! grid files, remapping of data from the provided model mesh to the output grid is done
! automatically.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  use grid_basic, only: type_grid
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_master
  USE grid_lonlat_basic                                      , ONLY: type_grid_lonlat, gather_lonlat_gridded_data_to_master_dp_2D, &
                                                                     gather_lonlat_gridded_data_to_master_dp_3D
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  use mpi_distributed_memory, only: gather_to_master

  use netcdf, only: NF90_UNLIMITED, NF90_INT, NF90_FLOAT, NF90_DOUBLE, NF90_MAX_VAR_DIMS, NF90_DEF_GRP
  use netcdf_field_name_options
  use netcdf_inquire_grid_mesh
  use netcdf_write_var_master
  use netcdf_basic_wrappers
  use netcdf_check_dimensions
  use netcdf_check_fields
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_sparse_matrix_utilities, only: gather_CSR_dist_to_master

  IMPLICIT NONE

CONTAINS

  ! ==== Write data to flexibly-defined fields =====
  ! ================================================

  ! Write data to a gridded output file
  SUBROUTINE write_to_field_multopt_grid_dp_2D(                       grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D data field defined on a grid to a NetCDF file variable on the same grid
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_grid_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nx, grid%ny))
    CALL gather_gridded_data_to_master( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nx, grid%ny,1))
      d_grid_with_time( :,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
      DEALLOCATE( d_grid_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_grid_dp_2D

  SUBROUTINE write_to_field_multopt_grid_dp_2D_monthly(               grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D monthly data field defined on a grid to a NetCDF file variable on the same grid
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_grid_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE          :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE          :: d_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nx, grid%ny, 12))
    CALL gather_gridded_data_to_master( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nx, grid%ny,12,1))
      d_grid_with_time( :,:,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, 12, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
      DEALLOCATE( d_grid_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_grid_dp_2D_monthly

  SUBROUTINE write_to_field_multopt_grid_dp_3D        (               grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 3-D data field defined on a grid to a NetCDF file variable on the same grid
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_grid_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti, nz
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE          :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE          :: d_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( d_grid_vec_partial,2)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nx, grid%ny, nz))
    CALL gather_gridded_data_to_master( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nx, grid%ny,nz,1))
      d_grid_with_time( :,:,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nz, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
      DEALLOCATE( d_grid_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_grid_dp_3D

  SUBROUTINE write_to_field_multopt_grid_dp_2D_notime(                grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D data field defined on a grid to a NetCDF file variable on the same grid
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_grid_dp_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nx, grid%ny))
    CALL gather_gridded_data_to_master( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_grid_dp_2D_notime

  SUBROUTINE write_to_field_multopt_grid_dp_2D_monthly_notime(        grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D monthly data field defined on a grid to a NetCDF file variable on the same grid
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nx, grid%ny, 12))
    CALL gather_gridded_data_to_master( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_grid_dp_2D_monthly_notime

  SUBROUTINE write_to_field_multopt_grid_dp_3D_notime(                grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 3-D data field defined on a grid to a NetCDF file variable on the same grid
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_grid_dp_3D_notime'
    INTEGER                                            :: id_var, nz
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( d_grid_vec_partial,2)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nx, grid%ny, nz))
    CALL gather_gridded_data_to_master( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_grid_dp_3D_notime

  ! Write data to a lon/lat-gridded output file
  SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D(                grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D data field defined on a lon/lat-grid to a NetCDF file variable on the same grid
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nlon, grid%nlat))
    CALL gather_lonlat_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nlon, grid%nlat,1))
      d_grid_with_time( :,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
      DEALLOCATE( d_grid_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D

  SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D_monthly(        grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D monthly data field defined on a lon/lat-grid to a NetCDF file variable on the same grid
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE          :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE          :: d_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nlon, grid%nlat, 12))
    CALL gather_lonlat_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nlon, grid%nlat,12,1))
      d_grid_with_time( :,:,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 12, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
      DEALLOCATE( d_grid_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D_monthly

  SUBROUTINE write_to_field_multopt_lonlat_grid_dp_3D(                grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 3-D data field defined on a lon/lat-grid to a NetCDF file variable on the same grid
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti, nz
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE          :: d_grid
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE          :: d_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( d_grid_vec_partial,2)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nlon, grid%nlat, nz))
    CALL gather_lonlat_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nlon, grid%nlat,nz,1))
      d_grid_with_time( :,:,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, nz, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
      DEALLOCATE( d_grid_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_lonlat_grid_dp_3D

  SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D_notime(         grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D data field defined on a lon/lat-grid to a NetCDF file variable on the same grid
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nlon, grid%nlat))
    CALL gather_lonlat_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D_notime

  SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D_monthly_notime( grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 2-D monthly data field defined on a lon/lat-grid to a NetCDF file variable on the same grid
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nlon, grid%nlat, 12))
    CALL gather_lonlat_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_lonlat_grid_dp_2D_monthly_notime

  SUBROUTINE write_to_field_multopt_lonlat_grid_dp_3D_notime(         grid, filename, ncid, field_name_options, d_grid_vec_partial)
    ! Write a 3-D data field defined on a lon/lat-grid to a NetCDF file variable on the same grid
    !
    ! d is stored in vector form, distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_3D_notime'
    INTEGER                                            :: id_var, nz
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( d_grid_vec_partial,2)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_grid( grid%nlon, grid%nlat, nz))
    CALL gather_lonlat_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_grid)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_grid)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_lonlat_grid_dp_3D_notime

  ! Write data to a mesh output file
  SUBROUTINE write_to_field_multopt_mesh_int_2D(                      mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_int_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: d_tot
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: d_tot_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV))
    CALL gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,1))
      d_tot_with_time( :,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_int_2D

  SUBROUTINE write_to_field_multopt_mesh_dp_2D(                       mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_tot
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV))
    CALL gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,1))
      d_tot_with_time( :,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D

  SUBROUTINE write_to_field_multopt_mesh_dp_2D_b(                     mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_2D_b'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_tot
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nTri))
    CALL gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nTri,1))
      d_tot_with_time( :,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nTri, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D_b

  SUBROUTINE write_to_field_multopt_mesh_dp_2D_monthly(               mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D monthly in the physical sense, so a 2-D array!)
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_tot_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV, 12))
    CALL gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,12,1))
      d_tot_with_time( :,:,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, 12, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D_monthly

  SUBROUTINE write_to_field_multopt_mesh_dp_3D(                       mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 3-D in the physical sense, so a 2-D array!)
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_tot_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV, mesh%nz))
    CALL gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,mesh%nz,1))
      d_tot_with_time( :,:,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, mesh%nz, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_3D

  SUBROUTINE write_to_field_multopt_mesh_dp_3D_b(                     mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 3-D in the physical sense, so a 2-D array!)
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_3D_b'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_tot_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nTri, mesh%nz))
    CALL gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nTri,mesh%nz,1))
      d_tot_with_time( :,:,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nTri, mesh%nz, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_3D_b

  SUBROUTINE write_to_field_multopt_mesh_dp_3D_ocean(                 mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 3-D in the physical sense, so a 2-D array!)
    !
    ! Write to the last time frame of the variable
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_3D_ocean'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE            :: d_tot_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV, C%nz_ocean))
    CALL gather_to_master( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,C%nz_ocean,1))
      d_tot_with_time( :,:,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, C%nz_ocean, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_3D_ocean

  SUBROUTINE write_to_field_multopt_mesh_int_2D_notime(               mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_int_2D_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_int_2D_notime

  SUBROUTINE write_to_field_multopt_mesh_int_2D_b_notime(             mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_int_2D_b_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_int_2D_b( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nTri))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_int_2D_b_notime

  SUBROUTINE write_to_field_multopt_mesh_int_2D_c_notime(             mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_int_2D_c_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_int_2D_c( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nE))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_int_2D_c_notime

  SUBROUTINE write_to_field_multopt_mesh_dp_2D_notime(                mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_2D_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D_notime

  SUBROUTINE write_to_field_multopt_mesh_dp_2D_b_notime(              mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_2D_b_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nTri))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D_b_notime

  SUBROUTINE write_to_field_multopt_mesh_dp_2D_c_notime(              mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_2D_c_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D_c( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nE))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D_c_notime

  SUBROUTINE write_to_field_multopt_mesh_dp_2D_monthly_notime(        mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D monthly in the physical sense, so a 2-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_2D_monthly_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV, 12))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D_monthly_notime

  SUBROUTINE write_to_field_multopt_mesh_dp_3D_notime(                mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 3-D in the physical sense, so a 2-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_3D_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nV, mesh%nz))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_3D_notime

  SUBROUTINE write_to_field_multopt_mesh_dp_3D_b_notime(              mesh, filename, ncid, field_name_options, d_partial)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 3-D in the physical sense, so a 2-D array!)
    !
    ! d is stored distributed over the processes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_partial

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_mesh_dp_3D_b_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Gather data to the master
    IF (par%master) ALLOCATE( d_tot( mesh%nTri, mesh%nz))
    CALL gather_to_master( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_3D_b_notime

  ! Write a scalar variable
  SUBROUTINE write_to_field_multopt_int_0D( filename, ncid, field_name_options, d)
    ! Write a 0-D data field to a NetCDF file variable
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    integer,                             INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_int_0D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if the file has a time dimension and variable
    CALL check_time( filename, ncid)

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Inquire file time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Check if the variable has time as a dimension
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, (/ d /), start = (/ ti /), count = (/ 1 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_int_0D

  SUBROUTINE write_to_field_multopt_dp_0D( filename, ncid, field_name_options, d)
    ! Write a 0-D data field to a NetCDF file variable
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp),                            INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_to_field_multopt_dp_0D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=1024)                                :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if the file has a time dimension and variable
    CALL check_time( filename, ncid)

    ! Inquire variable info
    CALL inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Inquire file time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Check if the variable has time as a dimension
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master( filename, ncid, id_var, (/ d /), start = (/ ti /), count = (/ 1 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_dp_0D

  ! Write new time value to file
  SUBROUTINE write_time_to_file( filename, ncid, time)
    ! Write new time value to file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'write_time_to_file'
    INTEGER                                            :: id_dim_time
    INTEGER                                            :: id_var_time
    INTEGER                                            :: nt

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check time dimension and variable
    CALL check_time( filename, ncid)

    ! Determine current length of time dimension in file
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! Inquire variable id
    CALL inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! Write time
    nt = nt + 1
    CALL write_var_master( filename, ncid, id_var_time, (/ time /), start = (/ nt /), count = (/ 1 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_time_to_file

  ! ===== Set up mesh/grid, time, and field variables in a NetCDF file =====
  ! ========================================================================

  ! Set up x/y-grid and gridded variables
  SUBROUTINE setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    ! Set up a regular x/y-grid in an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'setup_xy_grid_in_netcdf_file'
    INTEGER                                            :: id_dim_x
    INTEGER                                            :: id_dim_y
    INTEGER                                            :: id_var_x
    INTEGER                                            :: id_var_y
    INTEGER                                            :: id_var_lon
    INTEGER                                            :: id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create x/y dimensions
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_x), grid%nx, id_dim_x)
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_y), grid%ny, id_dim_y)

    ! Create and write x/y variables

    ! x
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_x), NF90_DOUBLE, (/ id_dim_x /), id_var_x)
    CALL add_attribute_char( filename, ncid, id_var_x, 'long_name', 'x-coordinate')
    CALL add_attribute_char( filename, ncid, id_var_x, 'units'    , 'm'           )
    CALL write_var_master( filename, ncid, id_var_x, grid%x)

    ! y
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_y), NF90_DOUBLE, (/ id_dim_y /), id_var_y)
    CALL add_attribute_char( filename, ncid, id_var_y, 'long_name', 'y-coordinate')
    CALL add_attribute_char( filename, ncid, id_var_y, 'units'    , 'm'           )
    CALL write_var_master( filename, ncid, id_var_y, grid%y)

    ! lon/lat-coordinates
    IF (ALLOCATED( grid%lon) .OR. ALLOCATED( grid%lat)) THEN

      ! Safety
      IF (.NOT. ALLOCATED( grid%lon)) CALL crash('grid has lat but no lon coordinates!')
      IF (.NOT. ALLOCATED( grid%lat)) CALL crash('grid has lon but no lat coordinates!')

      ! lon
      CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_lon), NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var_lon)
      CALL add_attribute_char( filename, ncid, id_var_lon, 'long_name', 'Longitude')
      CALL add_attribute_char( filename, ncid, id_var_lon, 'units'    , 'degrees east')
      CALL write_var_master( filename, ncid, id_var_lon, grid%lon)

      ! lat
      CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_lat), NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var_lat)
      CALL add_attribute_char( filename, ncid, id_var_lat, 'long_name', 'Latitude')
      CALL add_attribute_char( filename, ncid, id_var_lat, 'units'    , 'degrees north')
      CALL write_var_master( filename, ncid, id_var_lat, grid%lat)

    END IF ! IF (ALLOCATED( grid%lon)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_xy_grid_in_netcdf_file

  SUBROUTINE add_field_grid_int_2D(                      filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_int_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename, ncid)
    CALL check_y(    filename, ncid)
    CALL check_time( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_x, id_dim_y, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_int_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_int_2D

  SUBROUTINE add_field_grid_dp_2D(                       filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_dp_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename, ncid)
    CALL check_y(    filename, ncid)
    CALL check_time( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D

  SUBROUTINE add_field_grid_dp_2D_monthly(               filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_dp_2D_monthly'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(     filename, ncid)
    CALL check_y(     filename, ncid)
    CALL check_month( filename, ncid)
    CALL check_time(  filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time , id_dim_time )

    ! Safety
    IF (id_dim_x     == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y     == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time  == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_month, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_monthly

  SUBROUTINE add_field_grid_dp_3D(                       filename, ncid, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_dp_3D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename, ncid)
    CALL check_y(    filename, ncid)
    CALL check_zeta( filename, ncid)
    CALL check_time( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_3D

  SUBROUTINE add_field_grid_int_2D_notime(               filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_int_2D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y dimensions and variables are there
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! Safety
    IF (id_dim_x == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_x, id_dim_y /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_int_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_int_2D_notime

  SUBROUTINE add_field_grid_dp_2D_notime(                filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_dp_2D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x( filename, ncid)
    CALL check_y( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)

    ! Safety
    IF (id_dim_x == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_notime

  SUBROUTINE add_field_grid_dp_2D_monthly_notime(        filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(     filename, ncid)
    CALL check_y(     filename, ncid)
    CALL check_month( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! Safety
    IF (id_dim_x     == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y     == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_month /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_monthly_notime

  SUBROUTINE add_field_grid_dp_3D_notime(                filename, ncid, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_grid_dp_3D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename, ncid)
    CALL check_y(    filename, ncid)
    CALL check_zeta( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_3D_notime

  ! Set up lon/lat-grid and gridded variables
  SUBROUTINE setup_lonlat_grid_in_netcdf_file( filename, ncid, grid)
    ! Set up a regular lon/lat-grid in an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'setup_lonlat_grid_in_netcdf_file'
    INTEGER                                            :: id_dim_lon
    INTEGER                                            :: id_dim_lat
    INTEGER                                            :: id_var_lon
    INTEGER                                            :: id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create lon/lat dimensions
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_lon), grid%nlon, id_dim_lon)
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_lat), grid%nlat, id_dim_lat)

    ! Create and write lon/lat variables

    ! lon
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_lon), NF90_DOUBLE, (/ id_dim_lon /), id_var_lon)
    CALL add_attribute_char( filename, ncid, id_var_lon, 'long_name', 'Longitude')
    CALL add_attribute_char( filename, ncid, id_var_lon, 'units'    , 'degrees east')
    CALL write_var_master( filename, ncid, id_var_lon, grid%lon)

    ! lat
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_lat), NF90_DOUBLE, (/ id_dim_lat /), id_var_lat)
    CALL add_attribute_char( filename, ncid, id_var_lat, 'long_name', 'Latitude')
    CALL add_attribute_char( filename, ncid, id_var_lat, 'units'    , 'degrees north')
    CALL write_var_master( filename, ncid, id_var_lat, grid%lat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_lonlat_grid_in_netcdf_file

  SUBROUTINE add_field_lonlat_grid_int_2D(               filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_int_2D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(  filename, ncid)
    CALL check_lat(  filename, ncid)
    CALL check_time( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_lon  == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat  == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_lon, id_dim_lat, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_int_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_int_2D

  SUBROUTINE add_field_lonlat_grid_dp_2D(                filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_dp_2D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(  filename, ncid)
    CALL check_lat(  filename, ncid)
    CALL check_time( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_lon  == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat  == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_FLOAT, (/ id_dim_lon, id_dim_lat, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_dp_2D

  SUBROUTINE add_field_lonlat_grid_dp_2D_monthly(        filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_dp_2D_monthly'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_month, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(   filename, ncid)
    CALL check_lat(   filename, ncid)
    CALL check_month( filename, ncid)
    CALL check_time(  filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time , id_dim_time )

    ! Safety
    IF (id_dim_lon   == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat   == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time  == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_FLOAT, (/ id_dim_lon, id_dim_lat, id_dim_month, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_dp_2D_monthly

  SUBROUTINE add_field_lonlat_grid_dp_3D(                filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_dp_3D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(  filename, ncid)
    CALL check_lat(  filename, ncid)
    CALL check_zeta( filename, ncid)
    CALL check_time( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_lon  == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat  == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_FLOAT, (/ id_dim_lon, id_dim_lat, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_dp_3D

  SUBROUTINE add_field_lonlat_grid_int_2D_notime(        filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_int_2D_notime'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(  filename, ncid)
    CALL check_lat(  filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )

    ! Safety
    IF (id_dim_lon  == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat  == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_lon, id_dim_lat /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_int_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_int_2D_notime

  SUBROUTINE add_field_lonlat_grid_dp_2D_notime(         filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_dp_2D_notime'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(  filename, ncid)
    CALL check_lat(  filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )

    ! Safety
    IF (id_dim_lon  == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat  == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_FLOAT, (/ id_dim_lon, id_dim_lat /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_dp_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_dp_2D_notime

  SUBROUTINE add_field_lonlat_grid_dp_2D_monthly_notime( filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_month, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(   filename, ncid)
    CALL check_lat(   filename, ncid)
    CALL check_month( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon  , id_dim_lon  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat  , id_dim_lat  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month, id_dim_month)

    ! Safety
    IF (id_dim_lon   == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat   == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_FLOAT, (/ id_dim_lon, id_dim_lat, id_dim_month /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_dp_2D_monthly_notime

  SUBROUTINE add_field_lonlat_grid_dp_3D_notime(         filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a lon/lat-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_lonlat_grid_dp_3D_notime'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if lon, lat, and time dimensions and variables are there
    CALL check_lon(  filename, ncid)
    CALL check_lat(  filename, ncid)
    CALL check_zeta( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_lat , id_dim_lat )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta, id_dim_zeta)

    ! Safety
    IF (id_dim_lon  == -1) CALL crash('no lon dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_lat  == -1) CALL crash('no lat dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_FLOAT, (/ id_dim_lon, id_dim_lat, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_lonlat_grid_field_dp_3D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_lonlat_grid_dp_3D_notime

  ! Set up mesh and meshed variables
  SUBROUTINE setup_mesh_in_netcdf_file( filename, ncid, mesh)
    ! Set up a mesh in an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(INOUT) :: ncid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'setup_mesh_in_netcdf_file'

    INTEGER                                            :: id_dim_vi
    INTEGER                                            :: id_dim_ti
    INTEGER                                            :: id_dim_ci
    INTEGER                                            :: id_dim_ei
    INTEGER                                            :: id_dim_vori
    INTEGER                                            :: id_dim_two
    INTEGER                                            :: id_dim_three
    INTEGER                                            :: id_dim_four

    INTEGER                                            :: id_var_xmin
    INTEGER                                            :: id_var_xmax
    INTEGER                                            :: id_var_ymin
    INTEGER                                            :: id_var_ymax
    INTEGER                                            :: id_var_tol_dist
    INTEGER                                            :: id_var_lambda_M
    INTEGER                                            :: id_var_phi_M
    INTEGER                                            :: id_var_beta_stereo

    INTEGER                                            :: id_var_V
    INTEGER                                            :: id_var_nC
    INTEGER                                            :: id_var_C
    INTEGER                                            :: id_var_niTri
    INTEGER                                            :: id_var_iTri
    INTEGER                                            :: id_var_VBI

    INTEGER                                            :: id_var_Tri
    INTEGER                                            :: id_var_Tricc
    INTEGER                                            :: id_var_TriC
    INTEGER                                            :: id_var_TriBI

    INTEGER                                            :: id_var_E
    INTEGER                                            :: id_var_VE
    INTEGER                                            :: id_var_EV
    INTEGER                                            :: id_var_ETri
    INTEGER                                            :: id_var_EBI

    INTEGER                                            :: id_var_vi2vori
    INTEGER                                            :: id_var_ti2vori
    INTEGER                                            :: id_var_ei2vori
    INTEGER                                            :: id_var_vori2vi
    INTEGER                                            :: id_var_vori2ti
    INTEGER                                            :: id_var_vori2ei
    INTEGER                                            :: id_var_Vor
    INTEGER                                            :: id_var_VornC
    INTEGER                                            :: id_var_VorC
    INTEGER                                            :: id_var_nVVor
    INTEGER                                            :: id_var_VVor

    INTEGER                                            :: id_var_TriGC
    INTEGER                                            :: id_var_R
    INTEGER                                            :: id_var_A
    INTEGER                                            :: id_var_lon
    INTEGER                                            :: id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create mesh dimensions
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nV    ), mesh%nV    , id_dim_vi   )
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nTri  ), mesh%nTri  , id_dim_ti   )
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nC_mem), mesh%nC_mem, id_dim_ci   )
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nE    ), mesh%nE    , id_dim_ei   )
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nVor  ), mesh%nVor  , id_dim_vori )
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_two   ), 2          , id_dim_two  )
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_three ), 3          , id_dim_three)
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_four  ), 4          , id_dim_four )

  ! == Create mesh variables - metadata
  ! ===================================

    ! xmin
    CALL create_scalar_variable( filename, ncid, 'xmin', NF90_DOUBLE, id_var_xmin)
    CALL add_attribute_char( filename, ncid, id_var_xmin, 'long_name'  , 'Location of western domain border')
    CALL add_attribute_char( filename, ncid, id_var_xmin, 'units', 'm')

    ! xmax
    CALL create_scalar_variable( filename, ncid, 'xmax', NF90_DOUBLE, id_var_xmax)
    CALL add_attribute_char( filename, ncid, id_var_xmax, 'long_name'  , 'Location of eastern domain border')
    CALL add_attribute_char( filename, ncid, id_var_xmax, 'units', 'm')

    ! ymin
    CALL create_scalar_variable( filename, ncid, 'ymin', NF90_DOUBLE, id_var_ymin)
    CALL add_attribute_char( filename, ncid, id_var_ymin, 'long_name'  , 'Location of southern domain border')
    CALL add_attribute_char( filename, ncid, id_var_ymin, 'units', 'm')

    ! ymax
    CALL create_scalar_variable( filename, ncid, 'ymax', NF90_DOUBLE, id_var_ymax)
    CALL add_attribute_char( filename, ncid, id_var_ymax, 'long_name'  , 'Location of northern domain border')
    CALL add_attribute_char( filename, ncid, id_var_ymax, 'units', 'm')

    ! tol_dist
    CALL create_scalar_variable( filename, ncid, 'tol_dist', NF90_DOUBLE, id_var_tol_dist)
    CALL add_attribute_char( filename, ncid, id_var_tol_dist, 'long_name'  , 'Spatial tolerance (points within this distance are assumed identical)')
    CALL add_attribute_char( filename, ncid, id_var_tol_dist, 'units', 'm')

    ! lambda_M
    CALL create_scalar_variable( filename, ncid, 'lambda_M', NF90_DOUBLE, id_var_lambda_M)
    CALL add_attribute_char( filename, ncid, id_var_lambda_M, 'long_name'  , 'Longitude of the pole of the oblique stereographic projection')
    CALL add_attribute_char( filename, ncid, id_var_lambda_M, 'units', 'degrees east')

    ! phi_M
    CALL create_scalar_variable( filename, ncid, 'phi_M', NF90_DOUBLE, id_var_phi_M)
    CALL add_attribute_char( filename, ncid, id_var_phi_M, 'long_name'  , 'Latitude of the pole of the oblique stereographic projection')
    CALL add_attribute_char( filename, ncid, id_var_phi_M, 'units', 'degrees north')

    ! beta_stereo
    CALL create_scalar_variable( filename, ncid, 'beta_stereo', NF90_DOUBLE, id_var_beta_stereo)
    CALL add_attribute_char( filename, ncid, id_var_beta_stereo, 'long_name'  , 'Standard parallel of the oblique stereographic projection')
    CALL add_attribute_char( filename, ncid, id_var_beta_stereo, 'units', 'degrees')

  ! == Create mesh variables - vertex data
  ! ======================================

    ! V
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_V             ), NF90_DOUBLE, (/ id_dim_vi, id_dim_two   /), id_var_V             )
    CALL add_attribute_char( filename, ncid, id_var_V             , 'long_name'  , 'Vertex coordinates'         )
    CALL add_attribute_char( filename, ncid, id_var_V             , 'units'      , 'm'                          )
    ! nC
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_nC            ), NF90_INT   , (/ id_dim_vi               /), id_var_nC            )
    CALL add_attribute_char( filename, ncid, id_var_nC            , 'long_name'  , 'Number of vertex-vertex connections')
    ! C
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_C             ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_C             )
    CALL add_attribute_char( filename, ncid, id_var_C             , 'long_name'  , 'Vertex-vertex connections')
    CALL add_attribute_char( filename, ncid, id_var_C             , 'orientation', 'counter-clockwise'          )
    ! niTri
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_niTri         ), NF90_INT   , (/ id_dim_vi               /), id_var_niTri         )
    CALL add_attribute_char( filename, ncid, id_var_niTri         , 'long_name'  , 'Number of vertex-triangle connections')
    ! iTri
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_iTri          ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_iTri          )
    CALL add_attribute_char( filename, ncid, id_var_iTri          , 'long_name'  , 'Vertex-triangle connections')
    CALL add_attribute_char( filename, ncid, id_var_iTri          , 'orientation', 'counter-clockwise'          )
    ! VBI
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_VBI           ), NF90_INT   , (/ id_dim_vi               /), id_var_VBI           )
    CALL add_attribute_char( filename, ncid, id_var_VBI           , 'long_name'  , 'Vertex boundary index')
    CALL add_attribute_char( filename, ncid, id_var_VBI           , 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')

  ! == Create mesh variables - triangle data
  ! ========================================

    ! Tri
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_Tri           ), NF90_INT   , (/ id_dim_ti, id_dim_three /), id_var_Tri           )
    CALL add_attribute_char( filename, ncid, id_var_Tri           , 'long_name'  , 'Vertex indices per triangle')
    CALL add_attribute_char( filename, ncid, id_var_Tri           , 'orientation', 'counter-clockwise'          )
    ! Tricc
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_Tricc         ), NF90_DOUBLE, (/ id_dim_ti, id_dim_two   /), id_var_Tricc         )
    CALL add_attribute_char( filename, ncid, id_var_Tricc         , 'long_name'  , 'Triangle circumcentre coordinates')
    CALL add_attribute_char( filename, ncid, id_var_Tricc         , 'units'      , 'm'                          )
    ! TriC
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_TriC          ), NF90_INT   , (/ id_dim_ti, id_dim_three /), id_var_TriC          )
    CALL add_attribute_char( filename, ncid, id_var_TriC          , 'long_name'  , 'Triangle-triangle connections')
    CALL add_attribute_char( filename, ncid, id_var_TriC          , 'orientation', 'counter-clockwise, opposite from constituent vertices (i.e. first entry is opposite from first vertex)')
    ! TriBI
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_TriBI         ), NF90_INT   , (/ id_dim_ti               /), id_var_TriBI)
    CALL add_attribute_char( filename, ncid, id_var_TriBI         , 'long_name'  , 'Triangle boundary index')
    CALL add_attribute_char( filename, ncid, id_var_TriBI         , 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')

  ! == Create mesh variables - edge data
  ! ====================================

    ! E
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_E             ), NF90_DOUBLE, (/ id_dim_ei, id_dim_two   /), id_var_E             )
    CALL add_attribute_char( filename, ncid, id_var_E             , 'long_name'  , 'Edge midpoint coordinates')
    CALL add_attribute_char( filename, ncid, id_var_E             , 'units'      , 'm'                           )
    ! VE
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_VE            ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_VE            )
    CALL add_attribute_char( filename, ncid, id_var_VE            , 'long_name'  , 'Vertex-to-edge connectivity list')
    CALL add_attribute_char( filename, ncid, id_var_VE            , 'orientation', 'same as vertex-vertex connectivity list')
    ! EV
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_EV            ), NF90_INT   , (/ id_dim_ei,  id_dim_four /), id_var_EV            )
    CALL add_attribute_char( filename, ncid, id_var_EV            , 'long_name'  , 'Edge-to-vertex connectivity list')
    CALL add_attribute_char( filename, ncid, id_var_EV            , 'orientation', 'vi,vj,vl,vr (start,end,left,right)')
    ! ETri
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_ETri          ), NF90_INT   , (/ id_dim_ei,  id_dim_two  /), id_var_ETri          )
    CALL add_attribute_char( filename, ncid, id_var_ETri          , 'long_name'  , 'Edge-to-triangle connectivity list')
    CALL add_attribute_char( filename, ncid, id_var_ETri          , 'orientation', 'tl,tr (left,right)')
    ! EBI
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_EBI           ), NF90_INT   , (/ id_dim_ei               /), id_var_EBI           )
    CALL add_attribute_char( filename, ncid, id_var_EBI           , 'long_name'  , 'Edge boundary index')
    CALL add_attribute_char( filename, ncid, id_var_EBI           , 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')

  ! == Create mesh variables - Voronoi mesh data
  ! ============================================

    ! vi2vori
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vi2vori), NF90_INT, (/ id_dim_vi /), id_var_vi2vori)
    call add_attribute_char( filename, ncid, id_var_vi2vori, 'long_name' , 'Translation table from regular vertices to Voronoi vertices')
    ! ti2vori
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_ti2vori), NF90_INT, (/ id_dim_ti /), id_var_ti2vori)
    call add_attribute_char( filename, ncid, id_var_ti2vori, 'long_name' , 'Translation table from triangles to Voronoi vertices')
    ! ei2vori
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_ei2vori), NF90_INT, (/ id_dim_ei /), id_var_ei2vori)
    call add_attribute_char( filename, ncid, id_var_ei2vori, 'long_name' , 'Translation table from edges to Voronoi vertices')
    ! vori2vi
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vori2vi), NF90_INT, (/ id_dim_vori /), id_var_vori2vi)
    call add_attribute_char( filename, ncid, id_var_vori2vi, 'long_name' , 'Translation table from Voronoi vertices to regular vertices')
    ! vori2ti
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vori2ti), NF90_INT, (/ id_dim_vori /), id_var_vori2ti)
    call add_attribute_char( filename, ncid, id_var_vori2ti, 'long_name' , 'Translation table from Voronoi vertices to triangles')
    ! vori2ei
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_vori2ei), NF90_INT, (/ id_dim_vori /), id_var_vori2ei)
    call add_attribute_char( filename, ncid, id_var_vori2ei, 'long_name' , 'Translation table from Voronoi vertices to edges')
    ! Vor
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_Vor), NF90_DOUBLE, (/ id_dim_vori, id_dim_two /), id_var_Vor)
    call add_attribute_char( filename, ncid, id_var_Vor, 'long_name' , 'Voronoi vertex coordinates')
    CALL add_attribute_char( filename, ncid, id_var_Vor, 'units', 'm')
    ! VornC
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VornC), NF90_INT, (/ id_dim_vori /), id_var_VornC)
    call add_attribute_char( filename, ncid, id_var_VornC, 'long_name' , 'Number of Voronoi vertex connections')
    ! VorC
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VorC), NF90_INT, (/ id_dim_vori, id_dim_three /), id_var_VorC)
    call add_attribute_char( filename, ncid, id_var_VorC, 'long_name' , 'Indices of Voronoi vertex connections')
    ! nVVor
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_nVVor), NF90_INT, (/ id_dim_vi /), id_var_nVVor)
    call add_attribute_char( filename, ncid, id_var_nVVor, 'long_name' , 'Number of Voronoi vertices spanning each Voronoi cell')
    ! VVor
    call create_variable( filename, ncid, get_first_option_from_list( field_name_options_VVor), NF90_INT, (/ id_dim_vi, id_dim_ci /), id_var_VVor)
    call add_attribute_char( filename, ncid, id_var_VVor, 'long_name' , 'Indices of Voronoi vertices spanning each Voronoi cell')

  ! == Create mesh variables - secondary geometry data
  ! ==================================================

    ! TriGC
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_TriGC         ), NF90_DOUBLE, (/ id_dim_ti, id_dim_two   /), id_var_TriGC         )
    CALL add_attribute_char( filename, ncid, id_var_TriGC         , 'long_name'  , 'Triangle geometric centre coordinates')
    CALL add_attribute_char( filename, ncid, id_var_TriGC         , 'units'      , 'm')
    ! R
    CALL add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_R             ), long_name = 'Resolution', units = 'm')
    CALL inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_R             ), id_var_R)
    ! A
    CALL add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_A             ), long_name = 'Voronoi cell area', units = 'm^2')
    CALL inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_A             ), id_var_A)
    ! lon
    CALL add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_lon           ), long_name = 'Longitude', units = 'degrees east')
    CALL inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_lon           ), id_var_lon)
    ! lat
    CALL add_field_mesh_dp_2D_notime( filename, ncid, get_first_option_from_list( field_name_options_lat           ), long_name = 'Latitude' , units = 'degrees north')
    CALL inquire_var(                 filename, ncid, get_first_option_from_list( field_name_options_lat           ), id_var_lat)

  ! == Write mesh data to file
  ! ==========================

    ! Metadata
    CALL write_var_master(    filename, ncid, id_var_xmin       , mesh%xmin       )
    CALL write_var_master(    filename, ncid, id_var_xmax       , mesh%xmax       )
    CALL write_var_master(    filename, ncid, id_var_ymin       , mesh%ymin       )
    CALL write_var_master(    filename, ncid, id_var_ymax       , mesh%ymax       )
    CALL write_var_master(    filename, ncid, id_var_tol_dist   , mesh%tol_dist   )
    CALL write_var_master(    filename, ncid, id_var_lambda_M   , mesh%lambda_M   )
    CALL write_var_master(    filename, ncid, id_var_phi_M      , mesh%phi_M      )
    CALL write_var_master(    filename, ncid, id_var_beta_stereo, mesh%beta_stereo)

    ! Vertex data
    CALL write_var_master(    filename, ncid, id_var_V          , mesh%V          )
    CALL write_var_master(   filename, ncid, id_var_nC         , mesh%nC         )
    CALL write_var_master(   filename, ncid, id_var_C          , mesh%C          )
    CALL write_var_master(   filename, ncid, id_var_niTri      , mesh%niTri      )
    CALL write_var_master(   filename, ncid, id_var_iTri       , mesh%iTri       )
    CALL write_var_master(   filename, ncid, id_var_VBI        , mesh%VBI        )

    ! Triangle data
    CALL write_var_master(   filename, ncid, id_var_Tri        , mesh%Tri        )
    CALL write_var_master(    filename, ncid, id_var_Tricc      , mesh%Tricc      )
    CALL write_var_master(   filename, ncid, id_var_TriC       , mesh%TriC       )
    CALL write_var_master(   filename, ncid, id_var_TriBI      , mesh%TriBI      )

    ! Edge data
    CALL write_var_master(    filename, ncid, id_var_E          , mesh%E          )
    CALL write_var_master(   filename, ncid, id_var_VE         , mesh%VE         )
    CALL write_var_master(   filename, ncid, id_var_EV         , mesh%EV         )
    CALL write_var_master(   filename, ncid, id_var_ETri       , mesh%ETri       )
    CALL write_var_master(   filename, ncid, id_var_EBI        , mesh%EBI        )

    ! Voronoi mesh data
    call write_var_master(   filename, ncid, id_var_vi2vori    , mesh%vi2vori    )
    call write_var_master(   filename, ncid, id_var_ti2vori    , mesh%ti2vori    )
    call write_var_master(   filename, ncid, id_var_ei2vori    , mesh%ei2vori    )
    call write_var_master(   filename, ncid, id_var_vori2vi    , mesh%vori2vi    )
    call write_var_master(   filename, ncid, id_var_vori2ti    , mesh%vori2ti    )
    call write_var_master(   filename, ncid, id_var_vori2ei    , mesh%vori2ei    )
    call write_var_master(    filename, ncid, id_var_Vor        , mesh%Vor        )
    call write_var_master(   filename, ncid, id_var_VornC      , mesh%VornC      )
    call write_var_master(   filename, ncid, id_var_VorC       , mesh%VorC       )
    call write_var_master(   filename, ncid, id_var_nVVor      , mesh%nVVor      )
    call write_var_master(   filename, ncid, id_var_VVor       , mesh%VVor       )

    ! Secondary geometry data
    CALL write_var_master(    filename, ncid, id_var_TriGC      , mesh%TriGC      )
    CALL write_var_master(    filename, ncid, id_var_R          , mesh%R          )
    CALL write_var_master(    filename, ncid, id_var_A          , mesh%A          )
    CALL write_var_master(    filename, ncid, id_var_lon        , mesh%lon        )
    CALL write_var_master(    filename, ncid, id_var_lat        , mesh%lat        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_mesh_in_netcdf_file

  subroutine write_matrix_operators_to_netcdf_file( filename, ncid, mesh)
    !< Write all the matrix operators to the netcdf output file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(inout) :: ncid
    type(type_mesh),  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_matrix_operators_to_netcdf_file'

    call init_routine( routine_name)

    call write_mesh_translation_tables_to_netcdf_file( filename, ncid, mesh)

    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_a_a, 'M_ddx_a_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_a_b, 'M_ddx_a_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_a_c, 'M_ddx_a_c')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_b_a, 'M_ddx_b_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_b_b, 'M_ddx_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_b_c, 'M_ddx_b_c')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_c_a, 'M_ddx_c_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_c_b, 'M_ddx_c_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddx_c_c, 'M_ddx_c_c')

    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_a_a, 'M_ddy_a_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_a_b, 'M_ddy_a_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_a_c, 'M_ddy_a_c')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_b_a, 'M_ddy_b_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_b_b, 'M_ddy_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_b_c, 'M_ddy_b_c')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_c_a, 'M_ddy_c_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_c_b, 'M_ddy_c_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_ddy_c_c, 'M_ddy_c_c')

    ! call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_a_a, 'M_map_a_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_a_b, 'M_map_a_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_a_c, 'M_map_a_c')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_b_a, 'M_map_b_a')
    ! call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_b_b, 'M_map_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_b_c, 'M_map_b_c')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_c_a, 'M_map_c_a')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_c_b, 'M_map_c_b')
    ! call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M_map_c_c, 'M_map_c_c')

    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_ddx_b_b   , 'M2_ddx_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_ddy_b_b   , 'M2_ddy_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_d2dx2_b_b , 'M2_d2dx2_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_d2dxdy_b_b, 'M2_d2dxdy_b_b')
    call write_matrix_operator_to_netcdf_file( filename, ncid, mesh%M2_d2dy2_b_b , 'M2_d2dy2_b_b')

    call finalise_routine( routine_name)

  end subroutine write_matrix_operators_to_netcdf_file

  subroutine write_mesh_translation_tables_to_netcdf_file( filename, ncid, mesh)
    !< Write the mesh translation tables to the netcdf output file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(inout) :: ncid
    type(type_mesh),  intent(in   ) :: mesh

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER  :: routine_name = 'write_mesh_translation_tables_to_netcdf_file'
    integer :: id_dim_nz, id_dim_nzp1, id_dim_vi, id_dim_ti, id_dim_ei, id_dim_two, id_dim_three
    integer :: ierr
    integer :: grp_ncid
    integer :: id_dim_nna, id_dim_nnauv, id_dim_nnak, id_dim_nnaks, id_dim_nnakuv, id_dim_nnaksuv
    integer :: id_var_n2vi, id_var_n2viuv, id_var_n2vik, id_var_n2vikuv, id_var_n2viks, id_var_n2viksuv
    integer :: id_var_vi2n, id_var_viuv2n, id_var_vik2n, id_var_vikuv2n, id_var_viks2n, id_var_viksuv2n
    integer :: id_dim_nnb, id_dim_nnbuv, id_dim_nnbk, id_dim_nnbks, id_dim_nnbkuv, id_dim_nnbksuv
    integer :: id_var_n2ti, id_var_n2tiuv, id_var_n2tik, id_var_n2tikuv, id_var_n2tiks, id_var_n2tiksuv
    integer :: id_var_ti2n, id_var_tiuv2n, id_var_tik2n, id_var_tikuv2n, id_var_tiks2n, id_var_tiksuv2n
    integer :: id_dim_nnc, id_dim_nncuv, id_dim_nnck, id_dim_nncks, id_dim_nnckuv, id_dim_nncksuv
    integer :: id_var_n2ei, id_var_n2eiuv, id_var_n2eik, id_var_n2eikuv, id_var_n2eiks, id_var_n2eiksuv
    integer :: id_var_ei2n, id_var_eiuv2n, id_var_eik2n, id_var_eikuv2n, id_var_eiks2n, id_var_eiksuv2n

    call init_routine( routine_name)

    ! Create a group for them
    ierr = NF90_DEF_GRP( ncid, 'mesh_translation_tables', grp_ncid)

    call create_dimension( filename, grp_ncid, 'nz'  , mesh%nz  , id_dim_nz)
    call create_dimension( filename, grp_ncid, 'nzp1', mesh%nz+1, id_dim_nzp1)

    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_nV    ), id_dim_vi)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_nTri  ), id_dim_ti)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_nE    ), id_dim_ei)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_two   ), id_dim_two)
    call inquire_dim_multopt( filename, ncid, get_first_option_from_list( field_name_options_dim_three ), id_dim_three)

    ! a-grid (vertices)
    ! =================

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'nna'    , mesh%nna    , id_dim_nna)
    call create_dimension( filename, grp_ncid, 'nna'    , mesh%nna    , id_dim_nna)
    call create_dimension( filename, grp_ncid, 'nnauv'  , mesh%nnauv  , id_dim_nnauv)
    call create_dimension( filename, grp_ncid, 'nnak'   , mesh%nnak   , id_dim_nnak)
    call create_dimension( filename, grp_ncid, 'nnaks'  , mesh%nnaks  , id_dim_nnaks)
    call create_dimension( filename, grp_ncid, 'nnakuv' , mesh%nnakuv , id_dim_nnakuv)
    call create_dimension( filename, grp_ncid, 'nnaksuv', mesh%nnaksuv, id_dim_nnaksuv)

    ! Create variables
    call create_variable( filename, grp_ncid, 'n2vi'    , NF90_INT, [id_dim_nna                             ], id_var_n2vi)
    call create_variable( filename, grp_ncid, 'n2viuv'  , NF90_INT, [id_dim_nnauv  , id_dim_two             ], id_var_n2viuv)
    call create_variable( filename, grp_ncid, 'n2vik'   , NF90_INT, [id_dim_nnak   , id_dim_two             ], id_var_n2vik)
    call create_variable( filename, grp_ncid, 'n2vikuv' , NF90_INT, [id_dim_nnakuv , id_dim_three           ], id_var_n2vikuv)
    call create_variable( filename, grp_ncid, 'n2viks'  , NF90_INT, [id_dim_nnaks  , id_dim_two             ], id_var_n2viks)
    call create_variable( filename, grp_ncid, 'n2viksuv', NF90_INT, [id_dim_nnaksuv, id_dim_three           ], id_var_n2viksuv)
    call create_variable( filename, grp_ncid, 'vi2n'    , NF90_INT, [id_dim_vi                              ], id_var_vi2n)
    call create_variable( filename, grp_ncid, 'viuv2n'  , NF90_INT, [id_dim_vi                  , id_dim_two], id_var_viuv2n)
    call create_variable( filename, grp_ncid, 'vik2n'   , NF90_INT, [id_dim_vi     , id_dim_nz              ], id_var_vik2n)
    call create_variable( filename, grp_ncid, 'vikuv2n' , NF90_INT, [id_dim_vi     , id_dim_nz,   id_dim_two], id_var_vikuv2n)
    call create_variable( filename, grp_ncid, 'viks2n'  , NF90_INT, [id_dim_vi     , id_dim_nzp1            ], id_var_viks2n)
    call create_variable( filename, grp_ncid, 'viksuv2n', NF90_INT, [id_dim_vi     , id_dim_nzp1, id_dim_two], id_var_viksuv2n)

    ! Write variables
    call write_var_master( filename, grp_ncid, id_var_n2vi    , mesh%n2vi)
    call write_var_master( filename, grp_ncid, id_var_n2viuv  , mesh%n2viuv)
    call write_var_master( filename, grp_ncid, id_var_n2vik   , mesh%n2vik)
    call write_var_master( filename, grp_ncid, id_var_n2vikuv , mesh%n2vikuv)
    call write_var_master( filename, grp_ncid, id_var_n2viks  , mesh%n2viks)
    call write_var_master( filename, grp_ncid, id_var_n2viksuv, mesh%n2viksuv)
    call write_var_master( filename, grp_ncid, id_var_vi2n    , mesh%vi2n)
    call write_var_master( filename, grp_ncid, id_var_viuv2n  , mesh%viuv2n)
    call write_var_master( filename, grp_ncid, id_var_vik2n   , mesh%vik2n)
    call write_var_master( filename, grp_ncid, id_var_vikuv2n , mesh%vikuv2n)
    call write_var_master( filename, grp_ncid, id_var_viks2n  , mesh%viks2n)
    call write_var_master( filename, grp_ncid, id_var_viksuv2n, mesh%viksuv2n)

    ! b-grid (triangles)
    ! ==================

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'nnb'    , mesh%nnb    , id_dim_nnb)
    call create_dimension( filename, grp_ncid, 'nnb'    , mesh%nnb    , id_dim_nnb)
    call create_dimension( filename, grp_ncid, 'nnbuv'  , mesh%nnbuv  , id_dim_nnbuv)
    call create_dimension( filename, grp_ncid, 'nnbk'   , mesh%nnbk   , id_dim_nnbk)
    call create_dimension( filename, grp_ncid, 'nnbks'  , mesh%nnbks  , id_dim_nnbks)
    call create_dimension( filename, grp_ncid, 'nnbkuv' , mesh%nnbkuv , id_dim_nnbkuv)
    call create_dimension( filename, grp_ncid, 'nnbksuv', mesh%nnbksuv, id_dim_nnbksuv)

    ! Create variables
    call create_variable( filename, grp_ncid, 'n2ti'    , NF90_INT, [id_dim_nnb                             ], id_var_n2ti)
    call create_variable( filename, grp_ncid, 'n2tiuv'  , NF90_INT, [id_dim_nnbuv  , id_dim_two             ], id_var_n2tiuv)
    call create_variable( filename, grp_ncid, 'n2tik'   , NF90_INT, [id_dim_nnbk   , id_dim_two             ], id_var_n2tik)
    call create_variable( filename, grp_ncid, 'n2tikuv' , NF90_INT, [id_dim_nnbkuv , id_dim_three           ], id_var_n2tikuv)
    call create_variable( filename, grp_ncid, 'n2tiks'  , NF90_INT, [id_dim_nnbks  , id_dim_two             ], id_var_n2tiks)
    call create_variable( filename, grp_ncid, 'n2tiksuv', NF90_INT, [id_dim_nnbksuv, id_dim_three           ], id_var_n2tiksuv)
    call create_variable( filename, grp_ncid, 'ti2n'    , NF90_INT, [id_dim_ti                              ], id_var_ti2n)
    call create_variable( filename, grp_ncid, 'tiuv2n'  , NF90_INT, [id_dim_ti                  , id_dim_two], id_var_tiuv2n)
    call create_variable( filename, grp_ncid, 'tik2n'   , NF90_INT, [id_dim_ti     , id_dim_nz              ], id_var_tik2n)
    call create_variable( filename, grp_ncid, 'tikuv2n' , NF90_INT, [id_dim_ti     , id_dim_nz,   id_dim_two], id_var_tikuv2n)
    call create_variable( filename, grp_ncid, 'tiks2n'  , NF90_INT, [id_dim_ti     , id_dim_nzp1            ], id_var_tiks2n)
    call create_variable( filename, grp_ncid, 'tiksuv2n', NF90_INT, [id_dim_ti     , id_dim_nzp1, id_dim_two], id_var_tiksuv2n)

    ! Write variables
    call write_var_master( filename, grp_ncid, id_var_n2ti    , mesh%n2ti)
    call write_var_master( filename, grp_ncid, id_var_n2tiuv  , mesh%n2tiuv)
    call write_var_master( filename, grp_ncid, id_var_n2tik   , mesh%n2tik)
    call write_var_master( filename, grp_ncid, id_var_n2tikuv , mesh%n2tikuv)
    call write_var_master( filename, grp_ncid, id_var_n2tiks  , mesh%n2tiks)
    call write_var_master( filename, grp_ncid, id_var_n2tiksuv, mesh%n2tiksuv)
    call write_var_master( filename, grp_ncid, id_var_ti2n    , mesh%ti2n)
    call write_var_master( filename, grp_ncid, id_var_tiuv2n  , mesh%tiuv2n)
    call write_var_master( filename, grp_ncid, id_var_tik2n   , mesh%tik2n)
    call write_var_master( filename, grp_ncid, id_var_tikuv2n , mesh%tikuv2n)
    call write_var_master( filename, grp_ncid, id_var_tiks2n  , mesh%tiks2n)
    call write_var_master( filename, grp_ncid, id_var_tiksuv2n, mesh%tiksuv2n)

    ! c-grid (edges)
    ! ==============

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'nnc'    , mesh%nnc    , id_dim_nnc)
    call create_dimension( filename, grp_ncid, 'nnc'    , mesh%nnc    , id_dim_nnc)
    call create_dimension( filename, grp_ncid, 'nncuv'  , mesh%nncuv  , id_dim_nncuv)
    call create_dimension( filename, grp_ncid, 'nnck'   , mesh%nnck   , id_dim_nnck)
    call create_dimension( filename, grp_ncid, 'nncks'  , mesh%nncks  , id_dim_nncks)
    call create_dimension( filename, grp_ncid, 'nnckuv' , mesh%nnckuv , id_dim_nnckuv)
    call create_dimension( filename, grp_ncid, 'nncksuv', mesh%nncksuv, id_dim_nncksuv)

    ! Create variables
    call create_variable( filename, grp_ncid, 'n2ei'    , NF90_INT, [id_dim_nnc                             ], id_var_n2ei)
    call create_variable( filename, grp_ncid, 'n2eiuv'  , NF90_INT, [id_dim_nncuv  , id_dim_two             ], id_var_n2eiuv)
    call create_variable( filename, grp_ncid, 'n2eik'   , NF90_INT, [id_dim_nnck   , id_dim_two             ], id_var_n2eik)
    call create_variable( filename, grp_ncid, 'n2eikuv' , NF90_INT, [id_dim_nnckuv , id_dim_three           ], id_var_n2eikuv)
    call create_variable( filename, grp_ncid, 'n2eiks'  , NF90_INT, [id_dim_nncks  , id_dim_two             ], id_var_n2eiks)
    call create_variable( filename, grp_ncid, 'n2eiksuv', NF90_INT, [id_dim_nncksuv, id_dim_three           ], id_var_n2eiksuv)
    call create_variable( filename, grp_ncid, 'ei2n'    , NF90_INT, [id_dim_ei                              ], id_var_ei2n)
    call create_variable( filename, grp_ncid, 'eiuv2n'  , NF90_INT, [id_dim_ei                  , id_dim_two], id_var_eiuv2n)
    call create_variable( filename, grp_ncid, 'eik2n'   , NF90_INT, [id_dim_ei     , id_dim_nz              ], id_var_eik2n)
    call create_variable( filename, grp_ncid, 'eikuv2n' , NF90_INT, [id_dim_ei     , id_dim_nz,   id_dim_two], id_var_eikuv2n)
    call create_variable( filename, grp_ncid, 'eiks2n'  , NF90_INT, [id_dim_ei     , id_dim_nzp1            ], id_var_eiks2n)
    call create_variable( filename, grp_ncid, 'eiksuv2n', NF90_INT, [id_dim_ei     , id_dim_nzp1, id_dim_two], id_var_eiksuv2n)

    ! Write variables
    call write_var_master( filename, grp_ncid, id_var_n2ei    , mesh%n2ei)
    call write_var_master( filename, grp_ncid, id_var_n2eiuv  , mesh%n2eiuv)
    call write_var_master( filename, grp_ncid, id_var_n2eik   , mesh%n2eik)
    call write_var_master( filename, grp_ncid, id_var_n2eikuv , mesh%n2eikuv)
    call write_var_master( filename, grp_ncid, id_var_n2eiks  , mesh%n2eiks)
    call write_var_master( filename, grp_ncid, id_var_n2eiksuv, mesh%n2eiksuv)
    call write_var_master( filename, grp_ncid, id_var_ei2n    , mesh%ei2n)
    call write_var_master( filename, grp_ncid, id_var_eiuv2n  , mesh%eiuv2n)
    call write_var_master( filename, grp_ncid, id_var_eik2n   , mesh%eik2n)
    call write_var_master( filename, grp_ncid, id_var_eikuv2n , mesh%eikuv2n)
    call write_var_master( filename, grp_ncid, id_var_eiks2n  , mesh%eiks2n)
    call write_var_master( filename, grp_ncid, id_var_eiksuv2n, mesh%eiksuv2n)

    call finalise_routine( routine_name)

  end subroutine write_mesh_translation_tables_to_netcdf_file

  subroutine write_matrix_operator_to_netcdf_file( filename, ncid, A, name)
    !< Write a single matrix operator to the netcdf output file

    ! In/output variables:
    character(len=*),                 intent(in   ) :: filename
    integer,                          intent(inout) :: ncid
    type(type_sparse_matrix_CSR_dp),  intent(in   ) :: A
    character(len=*),                 intent(in   ) :: name

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER  :: routine_name = 'write_matrix_operator_to_netcdf_file'
    type(type_sparse_matrix_CSR_dp) :: A_tot
    integer                         :: ierr
    integer                         :: grp_ncid, id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz
    integer                         :: id_var_ptr, id_var_ind, id_var_val

    call init_routine( routine_name)

    ! Gather distributed matrix to the master
    call gather_CSR_dist_to_master( A, A_tot)

    ! Create a new NetCDF group for this matrix operator
    ierr = NF90_DEF_GRP( ncid, name, grp_ncid)

    ! Create dimensions
    call create_dimension( filename, grp_ncid, 'm'     , A_tot%m  , id_dim_m  )
    call create_dimension( filename, grp_ncid, 'mplus1', A_tot%m+1, id_dim_mp1)
    call create_dimension( filename, grp_ncid, 'n'     , A_tot%n  , id_dim_n  )
    call create_dimension( filename, grp_ncid, 'nnz'   , A_tot%nnz, id_dim_nnz)

    ! Create variables
    call create_variable( filename, grp_ncid, 'ptr', NF90_INT   , [id_dim_mp1], id_var_ptr)
    call create_variable( filename, grp_ncid, 'ind', NF90_INT   , [id_dim_nnz], id_var_ind)
    call create_variable( filename, grp_ncid, 'val', NF90_DOUBLE, [id_dim_nnz], id_var_val)

    ! Write to NetCDF
    call write_var_master( filename, grp_ncid, id_var_ptr, A_tot%ptr              )
    call write_var_master( filename, grp_ncid, id_var_ind, A_tot%ind( 1:A_tot%nnz))
    call write_var_master(  filename, grp_ncid, id_var_val, A_tot%val( 1:A_tot%nnz))

    call finalise_routine( routine_name)

  end subroutine write_matrix_operator_to_netcdf_file

  ! Set up mesh and meshed variables
  SUBROUTINE setup_CDF_in_netcdf_file( filename, ncid, ice)
    ! Set up a bedrock CDF in an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(INOUT) :: ncid
    TYPE(type_ice_model),                INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'setup_CDF_in_netcdf_file'

    INTEGER                                            :: id_dim_vi, id_dim_ti, id_dim_bin
    INTEGER                                            :: id_var_cdf, id_var_cdf_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create CDF bin dimension
    CALL create_dimension( filename, ncid, 'bin', C%subgrid_bedrock_cdf_nbins, id_dim_bin)

    ! Inquire mesh dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV,   id_dim_vi)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! Vertex data
    CALL create_variable( filename, ncid, 'bedrock_cdf', NF90_DOUBLE, (/ id_dim_vi, id_dim_bin /), id_var_cdf)
    CALL add_attribute_char( filename, ncid, id_var_cdf, 'long_name', 'Bedrock CDF of vertices')
    CALL add_attribute_char( filename, ncid, id_var_cdf, 'units'    , '%'                 )
    CALL write_var_master( filename, ncid, id_var_cdf, ice%bedrock_cdf)

    ! Triangle data
    CALL create_variable( filename, ncid, 'bedrock_cdf_b', NF90_DOUBLE, (/ id_dim_ti, id_dim_bin /), id_var_cdf_b)
    CALL add_attribute_char( filename, ncid, id_var_cdf_b, 'long_name', 'Bedrock CDF of triangles')
    CALL add_attribute_char( filename, ncid, id_var_cdf_b, 'units', '%')
    CALL write_var_master( filename, ncid, id_var_cdf_b, ice%bedrock_cdf_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_CDF_in_netcdf_file

  SUBROUTINE add_field_mesh_int_2D(                      filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_int_2D'
    INTEGER                                            :: id_dim_vi, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_vi, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_int_2D

  SUBROUTINE add_field_mesh_dp_2D(                       filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_2D'
    INTEGER                                            :: id_dim_vi, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D

  SUBROUTINE add_field_mesh_dp_2D_b(                     filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_2D_b'
    INTEGER                                            :: id_dim_ti, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time    , id_dim_time)

    ! Safety
    IF (id_dim_ti   == -1) CALL crash('no ti dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_b

  SUBROUTINE add_field_mesh_dp_2D_monthly(               filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_2D_monthly'
    INTEGER                                            :: id_dim_vi, id_dim_month, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_month(           filename, ncid)
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month , id_dim_month)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time )

    ! Safety
    IF (id_dim_vi    == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time  == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_month, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_monthly

  SUBROUTINE add_field_mesh_dp_3D(                       filename, ncid, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_3D'
    INTEGER                                            :: id_dim_vi, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_zeta(            filename, ncid)
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta  , id_dim_zeta)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_3D

  SUBROUTINE add_field_mesh_dp_3D_b(                     filename, ncid, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_3D_b'
    INTEGER                                            :: id_dim_ti, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_zeta(            filename, ncid)
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta    , id_dim_zeta)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time    , id_dim_time)

    ! Safety
    IF (id_dim_ti   == -1) CALL crash('no ti dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_3D_b

  SUBROUTINE add_field_mesh_dp_3D_ocean(                 filename, ncid, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_3D_ocean'
    INTEGER                                            :: id_dim_vi, id_dim_depth, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_depth(           filename, ncid)
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_depth , id_dim_depth)
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time )

    ! Safety
    IF (id_dim_vi    == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_depth == -1) CALL crash('no depth dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time  == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_depth, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_3D_ocean( filename, ncid, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_3D_ocean

  SUBROUTINE add_field_mesh_int_2D_notime(               filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_int_2D_notime'
    INTEGER                                            :: id_dim_vi, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_vi /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_int_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_int_2D_notime

  SUBROUTINE add_field_mesh_int_2D_b_notime(             filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_int_2D_b_notime'
    INTEGER                                            :: id_dim_ti, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )

    ! Safety
    IF (id_dim_ti   == -1) CALL crash('no ti dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_ti /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_int_2D_b( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_int_2D_b_notime

  SUBROUTINE add_field_mesh_int_2D_c_notime(             filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_int_2D_c_notime'
    INTEGER                                            :: id_dim_ei, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei  )

    ! Safety
    IF (id_dim_ei   == -1) CALL crash('no ei dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_ei /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_int_2D_c( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_int_2D_c_notime

  SUBROUTINE add_field_mesh_dp_2D_notime(                filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_2D_notime'
    INTEGER                                            :: id_dim_vi, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_notime

  SUBROUTINE add_field_mesh_dp_2D_b_notime(              filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_2D_b_notime'
    INTEGER                                            :: id_dim_ti, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )

    ! Safety
    IF (id_dim_ti   == -1) CALL crash('no ti dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D_b( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_b_notime

  SUBROUTINE add_field_mesh_dp_2D_c_notime(              filename, ncid, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_2D_c_notime'
    INTEGER                                            :: id_dim_ei, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nE, id_dim_ei  )

    ! Safety
    IF (id_dim_ei   == -1) CALL crash('no ei dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ei /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D_c( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_c_notime

  SUBROUTINE add_field_mesh_dp_2D_monthly_notime(        filename, ncid, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_2D_monthly_notime'
    INTEGER                                            :: id_dim_vi, id_dim_month, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_month(           filename, ncid)
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi   )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_month , id_dim_month)

    ! Safety
    IF (id_dim_vi    == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_month /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D_monthly( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_monthly_notime

  SUBROUTINE add_field_mesh_dp_3D_notime(                filename, ncid, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_3D_notime'
    INTEGER                                            :: id_dim_vi, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_zeta(            filename, ncid)
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta  , id_dim_zeta)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_3D( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_3D_notime

  SUBROUTINE add_field_mesh_dp_3D_b_notime(              filename, ncid, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_mesh_dp_3D_b_notime'
    INTEGER                                            :: id_dim_ti, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_zeta(            filename, ncid)
    CALL check_mesh_dimensions( filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti  )
    CALL inquire_dim_multopt( filename, ncid, field_name_options_zeta    , id_dim_zeta)

    ! Safety
    IF (id_dim_ti   == -1) CALL crash('no ti dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_ti, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_3D_b( filename, ncid, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_3D_b_notime

  ! Add extra dimensions
  SUBROUTINE add_time_dimension_to_file( filename, ncid)
    ! Add a time dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_time_dimension_to_file'
    INTEGER                                            :: id_dim_time
    INTEGER                                            :: id_var_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create time dimension
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_time), NF90_UNLIMITED, id_dim_time)

    ! Create time variable
    CALL create_variable(  filename, ncid, get_first_option_from_list( field_name_options_time), NF90_DOUBLE, (/ id_dim_time /), id_var_time)
    CALL add_attribute_char( filename, ncid, id_var_time, 'long_name', 'Time')
    CALL add_attribute_char( filename, ncid, id_var_time, 'units', 'years')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_time_dimension_to_file

  SUBROUTINE add_month_dimension_to_file( filename, ncid)
    ! Add a month dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_month_dimension_to_file'
    INTEGER                                            :: id_dim_month
    INTEGER                                            :: id_var_month

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create month dimension
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_month), 12, id_dim_month)

    ! Create month variable
    CALL create_variable(  filename, ncid, get_first_option_from_list( field_name_options_month), NF90_INT, (/ id_dim_month /), id_var_month)
    CALL add_attribute_char( filename, ncid, id_var_month, 'long_name', 'Month')
    CALL add_attribute_char( filename, ncid, id_var_month, 'units', '1-12')
    CALL add_attribute_char( filename, ncid, id_var_month, 'description', '1 = Jan, 2 = Feb, ..., 12 = Dec')

    ! Write month variable
    CALL write_var_master( filename, ncid, id_var_month, (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_month_dimension_to_file

  SUBROUTINE add_zeta_dimension_to_file( filename, ncid, zeta)
    ! Add a zeta dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: zeta

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_zeta_dimension_to_file'
    INTEGER                                            :: nz
    INTEGER                                            :: id_dim_zeta
    INTEGER                                            :: id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( zeta,1)

    ! Create month dimension
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_zeta), nz, id_dim_zeta)

    ! Create month variable
    CALL create_variable(  filename, ncid, get_first_option_from_list( field_name_options_zeta), NF90_DOUBLE, (/ id_dim_zeta /), id_var_zeta)
    CALL add_attribute_char( filename, ncid, id_var_zeta, 'long_name', 'Scaled vertical coordinate')
    CALL add_attribute_char( filename, ncid, id_var_zeta, 'units', '0-1')
    CALL add_attribute_char( filename, ncid, id_var_zeta, 'transformation', 'zeta = (h - z) / H; zeta = 0 at the ice surface; zeta = 1 at the ice base')

    ! Write month variable
    CALL write_var_master( filename, ncid, id_var_zeta, zeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_zeta_dimension_to_file

  SUBROUTINE add_depth_dimension_to_file( filename, ncid, depth)
    ! Add a depth dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: depth

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_depth_dimension_to_file'
    INTEGER                                            :: nz
    INTEGER                                            :: id_dim_depth
    INTEGER                                            :: id_var_depth

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( depth,1)

    ! Create month dimension
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_depth), nz, id_dim_depth)

    ! Create month variable
    CALL create_variable(  filename, ncid, get_first_option_from_list( field_name_options_depth), NF90_DOUBLE, (/ id_dim_depth /), id_var_depth)
    CALL add_attribute_char( filename, ncid, id_var_depth, 'long_name', 'Depth')
    CALL add_attribute_char( filename, ncid, id_var_depth, 'units', 'meters')

    ! Write month variable
    CALL write_var_master( filename, ncid, id_var_depth, depth)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_depth_dimension_to_file

  SUBROUTINE add_cdf_dimension_to_file( filename, ncid)
    ! Add a bin dimension and a CDF to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_cdf_dimension_to_file'
    INTEGER                                            :: id_dim_bins
    INTEGER                                            :: id_var_cdf, k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create month dimension
    CALL create_dimension( filename, ncid, 'bin', C%subgrid_bedrock_cdf_nbins, id_dim_bins)

    ! Create month variable
    CALL create_variable(  filename, ncid, 'bin', NF90_INT, (/ id_dim_bins /), id_var_cdf)
    CALL add_attribute_char( filename, ncid, id_var_cdf, 'long_name', 'CDF bin')
    CALL add_attribute_char( filename, ncid, id_var_cdf, 'units', 'Bin Nr.')
    CALL add_attribute_char( filename, ncid, id_var_cdf, 'description', 'Each of the bins that form the bedrock CDF')

    ! Write month variable
    CALL write_var_master( filename, ncid, id_var_cdf, (/ (k, k = 1, C%subgrid_bedrock_cdf_nbins) /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_cdf_dimension_to_file

  ! Add scalar variables
  SUBROUTINE add_field_dp_0D( filename, ncid, var_name, long_name, units)
    ! Add a 0-D variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_dp_0D'
    INTEGER                                            :: id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_dp_0D

  SUBROUTINE add_field_int_0D( filename, ncid, var_name, long_name, units)
    ! Add a 0-D variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'add_field_int_0D'
    INTEGER                                            :: id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire dimensions
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_int_0D

  ! ===== Generate procedural file names =====
  ! ==========================================

  SUBROUTINE generate_filename_XXXXXdotnc( filename_base, filename_base_XXXXXdotnc)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename_base
    CHARACTER(LEN=*),                    INTENT(OUT)   :: filename_base_XXXXXdotnc

    ! Local variables:
    CHARACTER(LEN=1024), PARAMETER                     :: routine_name = 'generate_filename_XXXXXdotnc'
    INTEGER                                            :: i, ierr
    CHARACTER(LEN=5)                                   :: i_str
    LOGICAL                                            :: ex

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Let the Master do this, and broadcast the result to the other processes,
    ! to prevent racing conditions (i.e. one processes creating file _00001 before
    ! the others get a chance to look, so they see _00001 already exists and try
    ! to create _00002 instead.)

    IF (par%master) THEN

      i = 1
      filename_base_XXXXXdotnc = TRIM( filename_base) // '_00001.nc'

      INQUIRE( FILE = filename_base_XXXXXdotnc, EXIST = ex)

      DO WHILE (ex)

        i = i+1

        IF     (i < 10) THEN
          WRITE( i_str,'(A,I1)') '0000',i
        ELSEIF (i < 100) THEN
          WRITE( i_str,'(A,I2)') '000',i
        ELSEIF (i < 1000) THEN
          WRITE( i_str,'(A,I3)') '00',i
        ELSEIF (i < 10000) THEN
          WRITE( i_str,'(A,I4)') '0',i
        ELSEIF (i < 100000) THEN
          WRITE( i_str,'(A,I5)') i
        ELSE
          CALL crash('10000 files of base name "' // TRIM( filename_base) // '" already exist!')
        END IF

        filename_base_XXXXXdotnc = TRIM( filename_base) // '_' // i_str // '.nc'

        INQUIRE( FILE = filename_base_XXXXXdotnc, EXIST = ex)

      END DO ! DO WHILE (ex)

    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( filename_base_XXXXXdotnc, LEN( filename_base_XXXXXdotnc), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE generate_filename_XXXXXdotnc

END MODULE netcdf_output
