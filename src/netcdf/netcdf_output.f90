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
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE main_configuration                                     , ONLY: C
  USE grid_basic                                             , ONLY: type_grid, gather_gridded_data_to_master_dp_2D, gather_gridded_data_to_master_dp_3D
  USE grid_lonlat_basic                                      , ONLY: type_grid_lonlat, gather_lonlat_gridded_data_to_master_dp_2D, &
                                                                     gather_lonlat_gridded_data_to_master_dp_3D
  USE math_utilities                                         , ONLY: permute_2D_int, permute_2D_dp, permute_3D_int, permute_3D_dp, &
                                                                     flip_1D_dp, flip_2D_x1_dp, flip_2D_x2_dp, flip_3D_x1_dp, flip_3D_x2_dp, flip_3D_x3_dp, &
                                                                     inverse_oblique_sg_projection
  USE mesh_types                                             , ONLY: type_mesh
  USE mpi_distributed_memory                                 , ONLY: gather_to_master_int_1D, gather_to_master_int_2D, gather_to_master_dp_1D, &
                                                                     gather_to_master_dp_2D

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
                          open_existing_netcdf_file_for_reading, close_netcdf_file, &
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
                          add_attribute_char, check_month, check_time, create_dimension, create_variable, create_scalar_variable, inquire_var, &
                          open_existing_netcdf_file_for_writing, switch_to_data_mode

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_grid_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nx, grid%ny,1))
      d_grid_with_time( :,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_grid_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nx, grid%ny,12,1))
      d_grid_with_time( :,:,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, 12, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_grid_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti, nz
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_grid_with_time( grid%nx, grid%ny,nz,1))
      d_grid_with_time( :,:,:,1) = d_grid
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nz, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_grid_dp_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_gridded_data_to_master_dp_2D( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master_dp_2D( filename, ncid, id_var, d_grid)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_grid)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_grid_dp_3D_notime'
    INTEGER                                            :: id_var, nz
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_gridded_data_to_master_dp_3D( grid, d_grid_vec_partial, d_grid)

    ! Write data to the variable
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_grid)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL write_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 12, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti, nz
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL write_var_master_dp_4D( filename, ncid, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, nz, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL write_var_master_dp_2D( filename, ncid, id_var, d_grid)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_grid)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_lonlat_grid_dp_3D_notime'
    INTEGER                                            :: id_var, nz
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_grid)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_int_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_int_1D( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,1))
      d_tot_with_time( :,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_int_2D( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_1D( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,1))
      d_tot_with_time( :,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_2D( filename, ncid, id_var, d_tot_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_2D

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_2D( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,12,1))
      d_tot_with_time( :,:,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, 12, 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_2D( d_partial, d_tot)

    ! Add "pretend" time dimension
    IF (par%master) THEN
      ALLOCATE( d_tot_with_time( mesh%nV,mesh%nz,1))
      d_tot_with_time( :,:,1) = d_tot
    END IF

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_3D( filename, ncid, id_var, d_tot_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, mesh%nz, 1 /) )

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
      DEALLOCATE( d_tot_with_time)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_3D

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_int_2D_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_int_1D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_int_1D( filename, ncid, id_var, d_tot)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_int_2D_b_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_int_1D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_int_1D( filename, ncid, id_var, d_tot)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_int_2D_c_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_int_1D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_int_1D( filename, ncid, id_var, d_tot)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_1D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_1D( filename, ncid, id_var, d_tot)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_b_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_1D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_1D( filename, ncid, id_var, d_tot)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_c_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_1D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_1D( filename, ncid, id_var, d_tot)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_2D_monthly_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_2D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_2D( filename, ncid, id_var, d_tot)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multopt_mesh_dp_3D_notime'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
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
    CALL gather_to_master_dp_2D( d_partial, d_tot)

    ! Inquire length of time dimension
    CALL inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    CALL write_var_master_dp_2D( filename, ncid, id_var, d_tot)

    ! Clean up after yourself
    IF (par%master) THEN
      DEALLOCATE( d_tot)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multopt_mesh_dp_3D_notime

  ! Write new time value to file
  SUBROUTINE write_time_to_file( filename, ncid, time)
    ! Write new time value to file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_time_to_file'
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
    CALL write_var_master_dp_1D( filename, ncid, id_var_time, (/ time /), start = (/ nt /), count = (/ 1 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_xy_grid_in_netcdf_file'
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
    CALL write_var_master_dp_1D( filename, ncid, id_var_x, grid%x)

    ! y
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_y), NF90_DOUBLE, (/ id_dim_y /), id_var_y)
    CALL add_attribute_char( filename, ncid, id_var_y, 'long_name', 'y-coordinate')
    CALL add_attribute_char( filename, ncid, id_var_y, 'units'    , 'm'           )
    CALL write_var_master_dp_1D( filename, ncid, id_var_y, grid%y)

    ! lon/lat-coordinates
    IF (ALLOCATED( grid%lon) .OR. ALLOCATED( grid%lat)) THEN

      ! Safety
      IF (.NOT. ALLOCATED( grid%lon)) CALL crash('grid has lat but no lon coordinates!')
      IF (.NOT. ALLOCATED( grid%lat)) CALL crash('grid has lon but no lat coordinates!')

      ! lon
      CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_lon), NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var_lon)
      CALL add_attribute_char( filename, ncid, id_var_lon, 'long_name', 'Longitude')
      CALL add_attribute_char( filename, ncid, id_var_lon, 'units'    , 'degrees east')
      CALL write_var_master_dp_2D( filename, ncid, id_var_lon, grid%lon)

      ! lat
      CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_lat), NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var_lat)
      CALL add_attribute_char( filename, ncid, id_var_lat, 'long_name', 'Latitude')
      CALL add_attribute_char( filename, ncid, id_var_lat, 'units'    , 'degrees north')
      CALL write_var_master_dp_2D( filename, ncid, id_var_lat, grid%lat)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_int_2D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_monthly'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_3D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_int_2D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_monthly_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_3D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_lonlat_grid_in_netcdf_file'
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
    CALL write_var_master_dp_1D( filename, ncid, id_var_lon, grid%lon)

    ! lat
    CALL create_variable( filename, ncid, get_first_option_from_list( field_name_options_lat), NF90_DOUBLE, (/ id_dim_lat /), id_var_lat)
    CALL add_attribute_char( filename, ncid, id_var_lat, 'long_name', 'Latitude')
    CALL add_attribute_char( filename, ncid, id_var_lat, 'units'    , 'degrees north')
    CALL write_var_master_dp_1D( filename, ncid, id_var_lat, grid%lat)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_int_2D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_dp_2D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_dp_2D_monthly'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_dp_3D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_int_2D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_dp_2D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_dp_2D_monthly_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_lonlat_grid_dp_3D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_mesh_in_netcdf_file'

    INTEGER                                            :: id_dim_vi
    INTEGER                                            :: id_dim_ti
    INTEGER                                            :: id_dim_ci
    INTEGER                                            :: id_dim_ei
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
    CALL create_dimension( filename, ncid, get_first_option_from_list( field_name_options_dim_nE    ), mesh%nE    , id_dim_ei  )
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

  ! == Create mesh variables - secondary geometry data
  ! ==================================================

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
    CALL write_var_master_dp_0D(    filename, ncid, id_var_xmin       , mesh%xmin       )
    CALL write_var_master_dp_0D(    filename, ncid, id_var_xmax       , mesh%xmax       )
    CALL write_var_master_dp_0D(    filename, ncid, id_var_ymin       , mesh%ymin       )
    CALL write_var_master_dp_0D(    filename, ncid, id_var_ymax       , mesh%ymax       )
    CALL write_var_master_dp_0D(    filename, ncid, id_var_tol_dist   , mesh%tol_dist   )
    CALL write_var_master_dp_0D(    filename, ncid, id_var_lambda_M   , mesh%lambda_M   )
    CALL write_var_master_dp_0D(    filename, ncid, id_var_phi_M      , mesh%phi_M      )
    CALL write_var_master_dp_0D(    filename, ncid, id_var_beta_stereo, mesh%beta_stereo)

    ! Vertex data
    CALL write_var_master_dp_2D(    filename, ncid, id_var_V          , mesh%V          )
    CALL write_var_master_int_1D(   filename, ncid, id_var_nC         , mesh%nC         )
    CALL write_var_master_int_2D(   filename, ncid, id_var_C          , mesh%C          )
    CALL write_var_master_int_1D(   filename, ncid, id_var_niTri      , mesh%niTri      )
    CALL write_var_master_int_2D(   filename, ncid, id_var_iTri       , mesh%iTri       )
    CALL write_var_master_int_1D(   filename, ncid, id_var_VBI        , mesh%VBI        )

    ! Triangle data
    CALL write_var_master_int_2D(   filename, ncid, id_var_Tri        , mesh%Tri        )
    CALL write_var_master_dp_2D(    filename, ncid, id_var_Tricc      , mesh%Tricc      )
    CALL write_var_master_int_2D(   filename, ncid, id_var_TriC       , mesh%TriC       )
    CALL write_var_master_int_1D(   filename, ncid, id_var_TriBI      , mesh%TriBI      )

    ! Edge data
    CALL write_var_master_dp_2D(    filename, ncid, id_var_E          , mesh%E          )
    CALL write_var_master_int_2D(   filename, ncid, id_var_VE         , mesh%VE         )
    CALL write_var_master_int_2D(   filename, ncid, id_var_EV         , mesh%EV         )
    CALL write_var_master_int_2D(   filename, ncid, id_var_ETri       , mesh%ETri       )
    CALL write_var_master_int_1D(   filename, ncid, id_var_EBI        , mesh%EBI        )

    ! Secondary geometry data
    CALL write_var_master_dp_1D(    filename, ncid, id_var_R          , mesh%R          )
    CALL write_var_master_dp_1D(    filename, ncid, id_var_A          , mesh%A          )
    CALL write_var_master_dp_1D(    filename, ncid, id_var_lon        , mesh%lon        )
    CALL write_var_master_dp_1D(    filename, ncid, id_var_lat        , mesh%lat        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_mesh_in_netcdf_file

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_int_2D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_monthly'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_3D'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_int_2D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_int_2D_b_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_int_2D_c_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_b_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_c_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_monthly_notime'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_3D_notime'
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

  ! Add extra dimensions
  SUBROUTINE add_time_dimension_to_file( filename, ncid)
    ! Add a time dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_time_dimension_to_file'
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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_month_dimension_to_file'
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
    CALL write_var_master_int_1D( filename, ncid, id_var_month, (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /) )

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_zeta_dimension_to_file'
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
    CALL write_var_master_dp_1D( filename, ncid, id_var_zeta, zeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_zeta_dimension_to_file

END MODULE netcdf_output