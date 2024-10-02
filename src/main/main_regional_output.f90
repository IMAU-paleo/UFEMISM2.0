MODULE main_regional_output

  ! Creating and writing to the main regional output files

! ===== Preamble =====
! ====================

  ! USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: happy, warning, crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE grid_basic                                             , ONLY: type_grid
  USE region_types                                           , ONLY: type_model_region
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, open_existing_netcdf_file_for_writing, close_netcdf_file
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, setup_xy_grid_in_netcdf_file, setup_CDF_in_netcdf_file, &
                                                                     add_time_dimension_to_file, add_zeta_dimension_to_file, add_month_dimension_to_file, &
                                                                     write_time_to_file, &
                                                                     add_field_mesh_int_2D, add_field_mesh_int_2D_notime, &
                                                                     add_field_mesh_dp_2D, add_field_mesh_dp_2D_notime, &
                                                                     add_field_mesh_dp_2D_monthly, add_field_mesh_dp_2D_monthly_notime, &
                                                                     add_field_mesh_dp_3D, add_field_mesh_dp_3D_notime, &
                                                                     add_field_mesh_dp_2D_b, add_field_mesh_dp_2D_b_notime, &
                                                                     add_field_mesh_dp_3D_b, add_field_mesh_dp_3D_b_notime, &
                                                                     add_field_grid_dp_2D, add_field_grid_dp_2D_notime, &
                                                                     add_field_grid_dp_3D, add_field_grid_dp_3D_notime, &
                                                                     add_field_grid_dp_2D_monthly, add_field_grid_dp_2D_monthly_notime, &
                                                                     add_field_dp_0D, &
                                                                     write_to_field_multopt_mesh_int_2D, write_to_field_multopt_mesh_int_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D, write_to_field_multopt_mesh_dp_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_monthly, write_to_field_multopt_mesh_dp_2D_monthly_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D, write_to_field_multopt_mesh_dp_3D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_b, write_to_field_multopt_mesh_dp_2D_b_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_b, write_to_field_multopt_mesh_dp_3D_b_notime, &
                                                                     write_to_field_multopt_grid_dp_2D, write_to_field_multopt_grid_dp_2D_notime, &
                                                                     write_to_field_multopt_grid_dp_2D_monthly, write_to_field_multopt_grid_dp_2D_monthly_notime, &
                                                                     write_to_field_multopt_grid_dp_3D, write_to_field_multopt_grid_dp_3D_notime, &
                                                                     write_to_field_multopt_dp_0D
  use remapping_main, only: map_from_mesh_to_xy_grid_2D, map_from_mesh_to_xy_grid_3D, map_from_mesh_to_xy_grid_2D_minval

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  ! == Write to main regional output files
  ! ======================================

  SUBROUTINE write_to_main_regional_output_file_mesh( region)
    ! Write to the main regional output NetCDF file - mesh version

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_main_regional_output_file_mesh'
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Writing to mesh output file "' // colour_string( TRIM( region%output_filename_mesh), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( region%output_filename_mesh, ncid)

    ! Write the time to the file
    CALL write_time_to_file( region%output_filename_mesh, ncid, region%time)

    ! Write the default data fields to the file
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'Hi')
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'Hb')
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'Hs')
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'SL')
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'u_surf')
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'v_surf')
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'uabs_surf')

    ! Write all user-defined data fields to the file
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_01)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_02)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_03)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_04)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_05)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_06)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_07)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_08)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_09)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_10)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_11)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_12)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_13)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_14)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_15)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_16)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_17)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_18)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_19)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_20)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_21)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_22)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_23)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_24)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_25)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_26)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_27)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_28)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_29)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_30)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_31)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_32)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_33)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_34)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_35)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_36)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_37)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_38)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_39)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_40)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_41)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_42)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_43)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_44)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_45)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_46)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_47)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_48)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_49)
    CALL write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_50)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_main_regional_output_file_mesh

  SUBROUTINE write_to_main_regional_output_file_grid( region)
    ! Write to the main regional output NetCDF file - grid version

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_main_regional_output_file_mesh'
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Writing to grid output file "' // colour_string( TRIM( region%output_filename_grid), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( region%output_filename_grid, ncid)

    ! Write the time to the file
    CALL write_time_to_file( region%output_filename_grid, ncid, region%time)

    ! Write the default data fields to the file
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'Hi')
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'Hb')
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'Hs')
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'SL')
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'u_surf')
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'v_surf')
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'uabs_surf')

    ! Write all user-defined data fields to the file
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_01)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_02)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_03)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_04)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_05)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_06)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_07)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_08)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_09)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_10)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_11)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_12)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_13)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_14)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_15)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_16)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_17)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_18)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_19)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_20)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_21)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_22)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_23)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_24)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_25)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_26)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_27)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_28)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_29)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_30)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_31)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_32)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_33)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_34)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_35)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_36)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_37)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_38)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_39)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_40)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_41)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_42)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_43)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_44)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_45)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_46)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_47)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_48)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_49)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_50)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_main_regional_output_file_grid

  SUBROUTINE write_to_main_regional_output_file_grid_ROI( region, grid, filename)
    ! Write to the main regional output NetCDF file - grid version

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region
    TYPE(type_grid)                                    , INTENT(IN)    :: grid
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_main_regional_output_file_grid_ROI'
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Writing to grid output file "' // colour_string( TRIM( filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( filename, ncid, region%time)

    ! Write the default data fields to the file
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'Hi')
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'Hb')
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'Hs')
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'SL')
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'u_surf')
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'v_surf')
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'uabs_surf')

    ! Write all user-defined data fields to the file
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_01)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_02)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_03)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_04)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_05)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_06)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_07)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_08)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_09)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_10)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_11)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_12)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_13)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_14)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_15)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_16)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_17)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_18)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_19)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_20)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_21)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_22)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_23)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_24)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_25)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_26)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_27)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_28)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_29)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_30)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_31)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_32)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_33)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_34)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_35)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_36)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_37)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_38)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_39)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_40)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_41)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_42)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_43)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_44)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_45)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_46)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_47)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_48)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_49)
    CALL write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_50)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_main_regional_output_file_grid_ROI

  SUBROUTINE write_to_main_regional_output_file_mesh_field( region, filename, ncid, choice_output_field)
    ! Write to the main regional output NetCDF file - mesh version
    !
    ! Write a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: choice_output_field

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_main_regional_output_file_mesh_field'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                            :: mask_int

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Allocate memory
    ALLOCATE( mask_int( region%mesh%vi1:region%mesh%vi2), source = 0)

    ! Add the specified data field to the file
    SELECT CASE (choice_output_field)
      CASE ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      CASE ('resolution')
        ! Do nothing - this is already part of the regular mesh data; only write this to the square grid output

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      CASE ('Hi_init')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hi_init', region%refgeo_init%Hi)
      CASE ('Hb_init')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hb_init', region%refgeo_init%Hb)
      CASE ('Hs_init')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hs_init', region%refgeo_init%Hs)
      CASE ('SL_init')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'SL_init', region%refgeo_init%SL)

      ! Present-day ice-sheet geometry
      CASE ('Hi_PD')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hi_PD', region%refgeo_PD%Hi)
      CASE ('Hb_PD')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hb_PD', region%refgeo_PD%Hb)
      CASE ('Hs_PD')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hs_PD', region%refgeo_PD%Hs)
      CASE ('SL_PD')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'SL_PD', region%refgeo_PD%SL)

      ! GIA equilibrium ice-sheet geometry
      CASE ('Hi_GIAeq')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hi_GIAeq', region%refgeo_GIAeq%Hi)
      CASE ('Hb_GIAeq')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hb_GIAeq', region%refgeo_GIAeq%Hb)
      CASE ('Hs_GIAeq')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hs_GIAeq', region%refgeo_GIAeq%Hs)
      CASE ('SL_GIAeq')
        CALL write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'SL_GIAeq', region%refgeo_GIAeq%SL)

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      CASE ('Hi')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hi', region%ice%Hi)
      CASE ('Hb')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hb', region%ice%Hb)
      CASE ('Hs')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hs', region%ice%Hs)
      CASE ('Hib')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hib', region%ice%Hib)
      CASE ('SL')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'SL', region%ice%SL)
      CASE ('TAF')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'TAF', region%ice%TAF)
      CASE ('Hi_eff')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hi_eff', region%ice%Hi_eff)
      CASE ('Hs_slope')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hs_slope', region%ice%Hs_slope)

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      CASE ('dHi')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi', region%ice%dHi)
      CASE ('dHb')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHb', region%ice%dHb)
      CASE ('dHs')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHs', region%ice%dHs)
      CASE ('dHib')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHib', region%ice%dHib)

    ! ===== Geometry rates of changes =====
    ! =====================================

      CASE ('dHi_dt')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt', region%ice%dHi_dt)
      CASE ('dHb_dt')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHb_dt', region%ice%dHb_dt)
      CASE ('dHs_dt')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHs_dt', region%ice%dHs_dt)
      CASE ('dHib_dt')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHib_dt', region%ice%dHib_dt)
      CASE ('dHi_dt_raw')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt_raw', region%ice%dHi_dt_raw)
      CASE ('dHi_dt_residual')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt_residual', region%ice%dHi_dt_residual)

    ! ===== Target quantities =====
    ! =============================

      CASE ('dHi_dt_target')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt_target', region%ice%dHi_dt_target)
      CASE ('uabs_surf_target')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'uabs_surf_target', region%ice%uabs_surf_target)

    ! ===== Masks =====
    ! =================

      CASE ('mask_icefree_land')
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_icefree_land', mask_int)
      CASE ('mask_icefree_ocean')
        WHERE (region%ice%mask_icefree_ocean .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_icefree_ocean', mask_int)
      CASE ('mask_grounded_ice')
        WHERE (region%ice%mask_grounded_ice .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_grounded_ice', mask_int)
      CASE ('mask_floating_ice')
        WHERE (region%ice%mask_floating_ice .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_floating_ice', mask_int)
      CASE ('mask_margin')
        WHERE (region%ice%mask_margin .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_margin', mask_int)
      CASE ('mask_gl_gr')
        WHERE (region%ice%mask_gl_gr .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_gl_gr', mask_int)
      CASE ('mask_gl_fl')
        WHERE (region%ice%mask_gl_fl .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_gl_fl', mask_int)
      CASE ('mask_cf_gr')
        WHERE (region%ice%mask_cf_gr .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_cf_gr', mask_int)
      CASE ('mask_cf_fl')
        WHERE (region%ice%mask_cf_fl .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_cf_fl', mask_int)
      CASE ('mask_coastline')
        WHERE (region%ice%mask_coastline .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_coastline', mask_int)
      CASE ('mask')
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask', region%ice%mask)
      CASE ('basin_ID')
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'basin_ID', region%ice%basin_ID)

    ! ===== Area fractions =====
    ! ==========================

      CASE ('fraction_gr')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'fraction_gr', region%ice%fraction_gr)
      CASE ('fraction_gr_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'fraction_gr_b', region%ice%fraction_gr_b)
      CASE ('fraction_margin')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'fraction_margin', region%ice%fraction_margin)

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ti', region%ice%Ti)
      CASE ('Ti_pmp')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ti_pmp', region%ice%Ti_pmp)
      CASE ('Ti_hom')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Ti_hom', region%ice%Ti_hom)
      CASE ('Cpi')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Cpi', region%ice%Cpi)
      CASE ('Ki')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ki', region%ice%Ki)
      CASE ('internal_heating')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'internal_heating', region%ice%internal_heating)
      CASE ('frictional_heating')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'frictional_heating', region%ice%frictional_heating)
      CASE ('A_flow')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'A_flow', region%ice%A_flow)

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      CASE ('u_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'u_3D', region%ice%u_3D)
      CASE ('v_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'v_3D', region%ice%v_3D)
      CASE ('u_3D_b')
        CALL write_to_field_multopt_mesh_dp_3D_b( region%mesh, filename, ncid, 'u_3D_b', region%ice%u_3D_b)
      CASE ('v_3D_b')
        CALL write_to_field_multopt_mesh_dp_3D_b( region%mesh, filename, ncid, 'v_3D_b', region%ice%v_3D_b)
      CASE ('w_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'w_3D', region%ice%w_3D)

      ! Vertically integrated
      CASE ('u_vav')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'u_vav', region%ice%u_vav)
      CASE ('v_vav')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'v_vav', region%ice%v_vav)
      CASE ('u_vav_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'u_vav_b', region%ice%u_vav_b)
      CASE ('v_vav_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'v_vav_b', region%ice%v_vav_b)
      CASE ('uabs_vav')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'uabs_vav', region%ice%uabs_vav)
      CASE ('uabs_vav_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'uabs_vav_b', region%ice%uabs_vav_b)

      ! Surface
      CASE ('u_surf')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'u_surf', region%ice%u_surf)
      CASE ('v_surf')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'v_surf', region%ice%v_surf)
      CASE ('u_surf_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'u_surf_b', region%ice%u_surf_b)
      CASE ('v_surf_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'v_surf_b', region%ice%v_surf_b)
      CASE ('w_surf')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'w_surf', region%ice%w_surf)
      CASE ('uabs_surf')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'uabs_surf', region%ice%uabs_surf)
      CASE ('uabs_surf_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'uabs_surf_b', region%ice%uabs_surf_b)

      ! Base
      CASE ('u_base')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'u_base', region%ice%u_base)
      CASE ('v_base')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'v_base', region%ice%v_base)
      CASE ('u_base_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'u_base_b', region%ice%u_base_b)
      CASE ('v_base_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'v_base_b', region%ice%v_base_b)
      CASE ('w_base')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'w_base', region%ice%w_base)
      CASE ('uabs_base')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'uabs_base', region%ice%uabs_base)
      CASE ('uabs_base_b')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'uabs_base_b', region%ice%uabs_base_b)

    ! === Strain rates ===
    ! ====================

      CASE ('du_dx_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'du_dx_3D', region%ice%du_dx_3D)
      CASE ('du_dy_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'du_dy_3D', region%ice%du_dy_3D)
      CASE ('du_dz_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'du_dz_3D', region%ice%du_dz_3D)
      CASE ('dv_dx_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dv_dx_3D', region%ice%dv_dx_3D)
      CASE ('dv_dy_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dv_dy_3D', region%ice%dv_dy_3D)
      CASE ('dv_dz_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dv_dz_3D', region%ice%dv_dz_3D)
      CASE ('dw_dx_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dw_dx_3D', region%ice%dw_dx_3D)
      CASE ('dw_dy_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dw_dy_3D', region%ice%dw_dy_3D)
      CASE ('dw_dz_3D')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dw_dz_3D', region%ice%dw_dz_3D)

    ! == Ice flow regime ==
    ! =====================

      CASE ('divQ')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'divQ', region%ice%divQ)
      CASE ('R_shear')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'R_shear', region%ice%R_shear)

    ! == Ice P/C time stepping ==
    ! ===========================

      CASE ('pc_truncation_error')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pc_truncation_error', region%ice%pc%tau_np1)
      CASE ('pc_untolerated_events')
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'pc_untolerated_events', region%ice%pc%tau_n_guilty)

    ! == Basal hydrology ==
    ! =====================

      CASE ('pore_water_pressure')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pore_water_pressure', region%ice%pore_water_pressure)
      CASE ('overburden_pressure')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'overburden_pressure', region%ice%overburden_pressure)
      CASE ('effective_pressure')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'effective_pressure', region%ice%effective_pressure)
      CASE ('pore_water_likelihood')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pore_water_likelihood', region%ice%pore_water_likelihood)
      CASE ('pore_water_fraction')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pore_water_fraction', region%ice%pore_water_fraction)

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      CASE ('till_friction_angle')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'till_friction_angle', region%ice%till_friction_angle)
      CASE ('bed_roughness')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'bed_roughness', region%ice%bed_roughness)
      CASE ('till_yield_stress')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'till_yield_stress', region%ice%till_yield_stress)
      CASE ('slid_alpha_sq')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'slid_alpha_sq', region%ice%slid_alpha_sq)
      CASE ('slid_beta_sq')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'slid_beta_sq', region%ice%slid_beta_sq)

      ! Basal friction and shear stress
      CASE ('basal_friction_coefficient')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'basal_friction_coefficient', region%ice%basal_friction_coefficient)
      CASE ('basal_shear_stress')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'basal_shear_stress', region%ice%basal_shear_stress)

    ! == Geothermal heat ==
    ! =====================

      CASE ('geothermal_heat_flux')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'geothermal_heat_flux', region%ice%geothermal_heat_flux)

    ! == Climate ==
    ! =============

      ! Main climate variables
      CASE ('T2m')
        CALL write_to_field_multopt_mesh_dp_2D_monthly( region%mesh, filename, ncid, 'T2m', region%climate%T2m)
      CASE ('Precip')
        CALL write_to_field_multopt_mesh_dp_2D_monthly( region%mesh, filename, ncid, 'Precip', region%climate%Precip)

    ! == Ocean ==
    ! ===========

      ! Main ocean variables
      CASE ('T_ocean')
        CALL warning('ocean temperature not implemented yet!')
      CASE ('S_ocean')
        CALL warning('ocean salinity not implemented yet!')

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      CASE ('SMB')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'SMB', region%SMB%SMB)

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      CASE ('BMB')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'BMB', region%BMB%BMB)

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      CASE ('H_lad')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'H_lad', region%BMB%laddie%H)
      CASE ('U_lad')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'U_lad', region%BMB%laddie%U)
      CASE ('V_lad')
        CALL write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'V_lad', region%BMB%laddie%V)
      CASE ('T_lad')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'T_lad', region%BMB%laddie%T)
      CASE ('S_lad')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'S_lad', region%BMB%laddie%S)

    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      CASE ('LMB')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'LMB', region%LMB%LMB)

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      CASE ('AMB')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'AMB', region%AMB%AMB)

    ! ===== End of user-defined output fields =====
    ! =============================================

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_output_field "' // TRIM( choice_output_field) // '"!')
    END SELECT

    ! Clean up after yourself
    DEALLOCATE( mask_int)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_main_regional_output_file_mesh_field

  SUBROUTINE write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, choice_output_field)
    ! Write to the main regional output NetCDF file - grid version
    !
    ! Write a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region
    TYPE(type_grid)                                    , INTENT(IN)    :: grid
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: choice_output_field

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_main_regional_output_file_grid_field'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                            :: d_mesh_vec_partial_2D
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                            :: d_grid_vec_partial_2D
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: d_grid_vec_partial_2D_monthly
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: d_grid_vec_partial_3D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Allocate memory
    ALLOCATE( d_mesh_vec_partial_2D(         region%mesh%vi1:region%mesh%vi2         ))
    ALLOCATE( d_grid_vec_partial_2D(         grid%n_loc                ))
    ALLOCATE( d_grid_vec_partial_2D_monthly( grid%n_loc, 12            ))
    ALLOCATE( d_grid_vec_partial_3D(         grid%n_loc, region%mesh%nz))

    ! Add the specified data field to the file
    SELECT CASE (choice_output_field)
      CASE ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      CASE ('resolution')
        d_mesh_vec_partial_2D = region%mesh%R( region%mesh%vi1:region%mesh%vi2)
        CALL map_from_mesh_to_xy_grid_2D_minval( region%mesh, grid, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'resolution', d_grid_vec_partial_2D)

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      CASE ('Hi_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_init%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hi_init', d_grid_vec_partial_2D)
      CASE ('Hb_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_init%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hb_init', d_grid_vec_partial_2D)
      CASE ('Hs_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_init%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hs_init', d_grid_vec_partial_2D)
      CASE ('SL_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_init%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'SL_init', d_grid_vec_partial_2D)

      ! Present-day ice-sheet geometry
      CASE ('Hi_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_PD%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hi_PD', d_grid_vec_partial_2D)
      CASE ('Hb_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_PD%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hb_PD', d_grid_vec_partial_2D)
      CASE ('Hs_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_PD%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hs_PD', d_grid_vec_partial_2D)
      CASE ('SL_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_PD%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'SL_PD', d_grid_vec_partial_2D)

      ! GIA equilibrium ice-sheet geometry
      CASE ('Hi_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_GIAeq%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hi_GIAeq', d_grid_vec_partial_2D)
      CASE ('Hb_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_GIAeq%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hb_GIAeq', d_grid_vec_partial_2D)
      CASE ('Hs_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_GIAeq%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hs_GIAeq', d_grid_vec_partial_2D)
      CASE ('SL_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%refgeo_GIAeq%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'SL_GIAeq', d_grid_vec_partial_2D)

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      CASE ('Hi')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hi', d_grid_vec_partial_2D)
      CASE ('Hb')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hb', d_grid_vec_partial_2D)
      CASE ('Hs')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hs', d_grid_vec_partial_2D)
      CASE ('Hib')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%Hib, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hib', d_grid_vec_partial_2D)
      CASE ('SL')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'SL', d_grid_vec_partial_2D)
      CASE ('TAF')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%TAF, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'TAF', d_grid_vec_partial_2D)
      CASE ('Hi_eff')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%Hi_eff, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hi_eff', d_grid_vec_partial_2D)
      CASE ('Hs_slope')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%Hs_slope, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hs_slope', d_grid_vec_partial_2D)

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      CASE ('dHi')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi', d_grid_vec_partial_2D)
      CASE ('dHb')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHb', d_grid_vec_partial_2D)
      CASE ('dHs')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHs', d_grid_vec_partial_2D)
      CASE ('dHib')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHib, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHib', d_grid_vec_partial_2D)

    ! ===== Geometry rates of change =====
    ! ====================================

      CASE ('dHi_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHi_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt', d_grid_vec_partial_2D)
      CASE ('dHb_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHb_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHb_dt', d_grid_vec_partial_2D)
      CASE ('dHs_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHs_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHs_dt', d_grid_vec_partial_2D)
      CASE ('dHib_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHib_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHib_dt', d_grid_vec_partial_2D)
      CASE ('dHi_dt_raw')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHi_dt_raw, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt_raw', d_grid_vec_partial_2D)
      CASE ('dHi_dt_residual')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHi_dt_residual, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt_residual', d_grid_vec_partial_2D)

    ! ===== Target quantities =====
    ! =============================

      CASE ('dHi_dt_target')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%dHi_dt_target, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt_target', d_grid_vec_partial_2D)

      CASE ('uabs_surf_target')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%uabs_surf_target, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'uabs_surf_target', d_grid_vec_partial_2D)

    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('mask_icefree_land')
      CASE ('mask_icefree_ocean')
      CASE ('mask_grounded_ice')
      CASE ('mask_floating_ice')
      CASE ('mask_margin')
      CASE ('mask_gl_gr')
      CASE ('mask_gl_fl')
      CASE ('mask_cf_gr')
      CASE ('mask_cf_fl')
      CASE ('mask_coastline')
      CASE ('mask')
      CASE ('basin_ID')

    ! ===== Area fractions =====
    ! ==========================

      ! NOTE: sub-grid area fractions cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('fraction_gr')
      CASE ('fraction_gr_b')
      CASE ('fraction_margin')

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%Ti, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Ti', d_grid_vec_partial_3D)
      CASE ('Ti_pmp')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%Ti_pmp, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Ti_pmp', d_grid_vec_partial_3D)
      CASE ('Ti_hom')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%Ti_hom, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Ti_hom', d_grid_vec_partial_2D)
      CASE ('Cpi')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%Cpi, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Cpi', d_grid_vec_partial_3D)
      CASE ('Ki')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%Ki, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Ki', d_grid_vec_partial_3D)
      CASE ('internal_heating')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%internal_heating, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'internal_heating', d_grid_vec_partial_3D)
      CASE ('frictional_heating')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%frictional_heating, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'frictional_heating', d_grid_vec_partial_2D)
      CASE ('A_flow')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%A_flow, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'A_flow', d_grid_vec_partial_3D)

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      CASE ('u_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%u_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'u_3D', d_grid_vec_partial_3D)
      CASE ('v_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%v_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'v_3D', d_grid_vec_partial_3D)
      CASE ('u_3D_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_3D_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%w_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'w_3D', d_grid_vec_partial_3D)

      ! Vertically integrated
      CASE ('u_vav')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%u_vav, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'u_vav', d_grid_vec_partial_2D)
      CASE ('v_vav')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%v_vav, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'v_vav', d_grid_vec_partial_2D)
      CASE ('u_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('uabs_vav')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%uabs_vav, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'uabs_vav', d_grid_vec_partial_2D)
      CASE ('uabs_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

      ! Surface
      CASE ('u_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%u_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'u_surf', d_grid_vec_partial_2D)
      CASE ('v_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%v_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'v_surf', d_grid_vec_partial_2D)
      CASE ('u_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%w_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'w_surf', d_grid_vec_partial_2D)
      CASE ('uabs_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%uabs_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'uabs_surf', d_grid_vec_partial_2D)
      CASE ('uabs_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

      ! Base
      CASE ('u_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%u_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'u_base', d_grid_vec_partial_2D)
      CASE ('v_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%v_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'v_base', d_grid_vec_partial_2D)
      CASE ('u_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%w_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'w_base', d_grid_vec_partial_2D)
      CASE ('uabs_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%uabs_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'uabs_base', d_grid_vec_partial_2D)
      CASE ('uabs_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

    ! === Strain rates ===
    ! ====================

      CASE ('du_dx_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%du_dx_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'du_dx_3D', d_grid_vec_partial_3D)
      CASE ('du_dy_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%du_dy_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'du_dy_3D', d_grid_vec_partial_3D)
      CASE ('du_dz_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%du_dz_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'du_dz_3D', d_grid_vec_partial_3D)
      CASE ('dv_dx_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%dv_dx_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dv_dx_3D', d_grid_vec_partial_3D)
      CASE ('dv_dy_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%dv_dy_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dv_dy_3D', d_grid_vec_partial_3D)
      CASE ('dv_dz_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%dv_dz_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dv_dz_3D', d_grid_vec_partial_3D)
      CASE ('dw_dx_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%dw_dx_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dw_dx_3D', d_grid_vec_partial_3D)
      CASE ('dw_dy_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%dw_dy_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dw_dy_3D', d_grid_vec_partial_3D)
      CASE ('dw_dz_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%ice%dw_dz_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dw_dz_3D', d_grid_vec_partial_3D)

    ! == Ice flow regime ==
    ! =====================

      CASE ('divQ')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%divQ, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'divQ', d_grid_vec_partial_2D)
      CASE ('R_shear')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%R_shear, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'R_shear', d_grid_vec_partial_2D)

    ! == Ice P/C time stepping ==
    ! ===========================

      CASE ('pc_truncation_error')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%pc%tau_np1, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pc_truncation_error', d_grid_vec_partial_2D)
      CASE ('pc_untolerated_events')
        ! DENK DROM : Not gridable

    ! == Basal hydrology ==
    ! =====================

      CASE ('pore_water_pressure')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%pore_water_pressure, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pore_water_pressure', d_grid_vec_partial_2D)
      CASE ('overburden_pressure')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%overburden_pressure, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'overburden_pressure', d_grid_vec_partial_2D)
      CASE ('effective_pressure')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%effective_pressure, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'effective_pressure', d_grid_vec_partial_2D)
      CASE ('pore_water_likelihood')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%pore_water_likelihood, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pore_water_likelihood', d_grid_vec_partial_2D)
      CASE ('pore_water_fraction')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%pore_water_fraction, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pore_water_fraction', d_grid_vec_partial_2D)

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      CASE ('till_friction_angle')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%till_friction_angle, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'till_friction_angle', d_grid_vec_partial_2D)
      CASE ('bed_roughness')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%bed_roughness, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness', d_grid_vec_partial_2D)
      CASE ('till_yield_stress')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%till_yield_stress, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'till_yield_stress', d_grid_vec_partial_2D)
      CASE ('slid_alpha_sq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%slid_alpha_sq, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'slid_alpha_sq', d_grid_vec_partial_2D)
      CASE ('slid_beta_sq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%slid_beta_sq, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'slid_beta_sq', d_grid_vec_partial_2D)

      ! Basal friction and shear stress
      CASE ('basal_friction_coefficient')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%basal_friction_coefficient, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'basal_friction_coefficient', d_grid_vec_partial_2D)
      CASE ('basal_shear_stress')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%basal_shear_stress, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'basal_shear_stress', d_grid_vec_partial_2D)

    ! == Geothermal heat ==
    ! =====================

      CASE ('geothermal_heat_flux')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%ice%geothermal_heat_flux, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'geothermal_heat_flux', d_grid_vec_partial_2D)

    ! == Climate ==
    ! =============

      ! Main climate variables
      CASE ('T2m')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%climate%T2m, d_grid_vec_partial_2D_monthly)
        CALL write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, 'T2m', d_grid_vec_partial_2D_monthly)
      CASE ('Precip')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, grid, region%climate%Precip, d_grid_vec_partial_2D_monthly)
        CALL write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, 'Precip', d_grid_vec_partial_2D_monthly)

    ! == Ocean ==
    ! ==========================

      ! Main ocean variables
      CASE ('T_ocean')
        CALL warning('ocean temperature not implemented yet!')
      CASE ('S_ocean')
        CALL warning('ocean salinity not implemented yet!')

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      CASE ('SMB')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%SMB%SMB, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'SMB', d_grid_vec_partial_2D)

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      CASE ('BMB')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%BMB%BMB, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'BMB', d_grid_vec_partial_2D)

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      CASE ('H_lad')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%BMB%laddie%H, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'H_lad', d_grid_vec_partial_2D)
      CASE ('U_lad')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('V_lad')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('T_lad')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%BMB%laddie%T, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'T_lad', d_grid_vec_partial_2D)
      CASE ('S_lad')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%BMB%laddie%S, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'S_lad', d_grid_vec_partial_2D)

    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      CASE ('LMB')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%LMB%LMB, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'LMB', d_grid_vec_partial_2D)

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      CASE ('AMB')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, grid, region%AMB%AMB, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'AMB', d_grid_vec_partial_2D)

    ! ===== End of user-defined output fields =====
    ! =============================================

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_output_field "' // TRIM( choice_output_field) // '"!')
    END SELECT

    ! Clean up after yourself
    DEALLOCATE( d_mesh_vec_partial_2D)
    DEALLOCATE( d_grid_vec_partial_2D)
    DEALLOCATE( d_grid_vec_partial_2D_monthly)
    DEALLOCATE( d_grid_vec_partial_3D)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_main_regional_output_file_grid_field

  SUBROUTINE write_to_scalar_regional_output_file( region)
    ! Write to the scalar regional output NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_scalar_regional_output_file'
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Writing to scalar output file "' // colour_string( TRIM( region%output_filename_scalar), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( region%output_filename_scalar, ncid)

    ! Write the time to the file
    CALL write_time_to_file( region%output_filename_scalar, ncid, region%time)

    ! Write the default data fields to the file
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_area',          region%scalars%ice_area)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume',        region%scalars%ice_volume)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume_af',     region%scalars%ice_volume_af)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_area_PD',       region%scalars%ice_area_PD)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume_PD',     region%scalars%ice_volume_PD)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'ice_volume_af_PD',  region%scalars%ice_volume_af_PD)

    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_total',         region%scalars%SMB_total)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_gr',            region%scalars%SMB_gr)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_fl',            region%scalars%SMB_fl)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_land',          region%scalars%SMB_land)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'SMB_ocean',         region%scalars%SMB_ocean)

    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_total',         region%scalars%BMB_total)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_gr',            region%scalars%BMB_gr)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_fl',            region%scalars%BMB_fl)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_land',          region%scalars%BMB_land)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'BMB_ocean',         region%scalars%BMB_ocean)

    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'LMB_total',         region%scalars%LMB_total)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'LMB_gr',            region%scalars%LMB_gr)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'LMB_fl',            region%scalars%LMB_fl)

    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_total',         region%scalars%AMB_total)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_gr',            region%scalars%AMB_gr)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_fl',            region%scalars%AMB_fl)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_land',          region%scalars%AMB_land)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'AMB_ocean',         region%scalars%AMB_ocean)

    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'gl_flux',           region%scalars%gl_flux)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'cf_gr_flux',        region%scalars%cf_gr_flux)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'cf_fl_flux',        region%scalars%cf_fl_flux)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'margin_land_flux',  region%scalars%margin_land_flux)
    CALL write_to_field_multopt_dp_0D( region%output_filename_scalar, ncid, 'margin_ocean_flux', region%scalars%margin_ocean_flux)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_scalar_regional_output_file

  ! == Create main regional output files
  ! ====================================

  SUBROUTINE create_main_regional_output_file_mesh( region)
    ! Create the main regional output NetCDF file - mesh version

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'create_main_regional_output_file_mesh'
    CHARACTER(LEN=256)                                                 :: filename_base
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'main_output_' // region%name
    CALL generate_filename_XXXXXdotnc( filename_base, region%output_filename_mesh)

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Creating mesh output file "' // colour_string( TRIM( region%output_filename_mesh), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( region%output_filename_mesh, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( region%output_filename_mesh, ncid, region%mesh)

    IF (C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .OR. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') THEN
      ! Set up bedrock CDF in the file
      CALL setup_CDF_in_netcdf_file( region%output_filename_mesh, ncid, region%ice)
    END IF

    ! Add time, zeta, and month dimensions+variables to the file
    CALL add_time_dimension_to_file(  region%output_filename_mesh, ncid)
    CALL add_zeta_dimension_to_file(  region%output_filename_mesh, ncid, region%mesh%zeta)
    CALL add_month_dimension_to_file( region%output_filename_mesh, ncid)

    ! Add the default data fields to the file
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'Hi')
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'Hb')
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'Hs')
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'SL')
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'u_surf')
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'v_surf')
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'uabs_surf')

    ! Add all user-defined data fields to the file
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_01)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_02)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_03)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_04)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_05)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_06)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_07)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_08)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_09)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_10)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_11)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_12)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_13)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_14)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_15)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_16)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_17)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_18)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_19)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_20)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_21)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_22)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_23)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_24)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_25)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_26)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_27)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_28)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_29)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_30)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_31)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_32)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_33)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_34)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_35)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_36)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_37)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_38)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_39)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_40)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_41)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_42)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_43)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_44)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_45)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_46)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_47)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_48)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_49)
    CALL create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_50)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_main_regional_output_file_mesh

  SUBROUTINE create_main_regional_output_file_grid( region)
    ! Create the main regional output NetCDF file - grid version

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'create_main_regional_output_file_grid'
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    region%output_filename_grid = TRIM( C%output_dir) // 'main_output_' // region%name // '_grid.nc'

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Creating grid output file "' // colour_string( TRIM( region%output_filename_grid), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( region%output_filename_grid, ncid)

    ! Set up the grid in the file
    CALL setup_xy_grid_in_netcdf_file( region%output_filename_grid, ncid, region%output_grid)

    ! Add time, zeta, and month dimensions+variables to the file
    CALL add_time_dimension_to_file(  region%output_filename_grid, ncid)
    CALL add_zeta_dimension_to_file(  region%output_filename_grid, ncid, region%mesh%zeta)
    CALL add_month_dimension_to_file( region%output_filename_grid, ncid)

    ! Add the default data fields to the file
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'Hi')
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'Hb')
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'Hs')
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'SL')
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'u_surf')
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'v_surf')
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'uabs_surf')

    ! Add all user-defined data fields to the file
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_01)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_02)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_03)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_04)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_05)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_06)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_07)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_08)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_09)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_10)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_11)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_12)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_13)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_14)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_15)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_16)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_17)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_18)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_19)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_20)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_21)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_22)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_23)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_24)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_25)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_26)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_27)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_28)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_29)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_30)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_31)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_32)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_33)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_34)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_35)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_36)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_37)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_38)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_39)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_40)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_41)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_42)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_43)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_44)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_45)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_46)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_47)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_48)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_49)
    CALL create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_50)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_main_regional_output_file_grid

  SUBROUTINE create_main_regional_output_file_grid_ROI( region, grid, filename)
    ! Create the main regional output NetCDF file - grid version

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region
    TYPE(type_grid)                                    , INTENT(IN)    :: grid
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'create_main_regional_output_file_grid_ROI'
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Creating ROI output file "' // colour_string( TRIM( filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the grid in the file
    CALL setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    ! Add time, zeta, and month dimensions+variables to the file
    CALL add_time_dimension_to_file(  filename, ncid)
    CALL add_zeta_dimension_to_file(  filename, ncid, region%mesh%zeta)
    CALL add_month_dimension_to_file( filename, ncid)

    ! Add the default data fields to the file
    CALL create_main_regional_output_file_grid_field( filename, ncid, 'Hi')
    CALL create_main_regional_output_file_grid_field( filename, ncid, 'Hb')
    CALL create_main_regional_output_file_grid_field( filename, ncid, 'Hs')
    CALL create_main_regional_output_file_grid_field( filename, ncid, 'SL')
    CALL create_main_regional_output_file_grid_field( filename, ncid, 'u_surf')
    CALL create_main_regional_output_file_grid_field( filename, ncid, 'v_surf')
    CALL create_main_regional_output_file_grid_field( filename, ncid, 'uabs_surf')

    ! Add all user-defined data fields to the file
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_01)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_02)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_03)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_04)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_05)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_06)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_07)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_08)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_09)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_10)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_11)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_12)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_13)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_14)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_15)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_16)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_17)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_18)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_19)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_20)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_21)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_22)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_23)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_24)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_25)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_26)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_27)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_28)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_29)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_30)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_31)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_32)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_33)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_34)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_35)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_36)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_37)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_38)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_39)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_40)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_41)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_42)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_43)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_44)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_45)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_46)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_47)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_48)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_49)
    CALL create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_50)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_main_regional_output_file_grid_ROI

  SUBROUTINE create_main_regional_output_file_mesh_field( filename, ncid, choice_output_field)
    ! Create the main regional output NetCDF file - mesh version
    !
    ! Add a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: choice_output_field

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'create_main_regional_output_file_mesh_field'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Add the specified data field to the file
    SELECT CASE (choice_output_field)
      CASE ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      CASE ('resolution')
        ! Do nothing - this is already part of the regular mesh data; only write this to the square grid output

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      CASE ('Hi_init')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hi_init', long_name = 'Initial ice thickness', units = 'm')
      CASE ('Hb_init')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hb_init', long_name = 'Initial bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs_init')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hs_init', long_name = 'Initial surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL_init')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'SL_init', long_name = 'Initial geoid elevation', units = 'm w.r.t. PD sea level')

      ! Present-day ice-sheet geometry
      CASE ('Hi_PD')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hi_PD', long_name = 'Present-day ice thickness', units = 'm')
      CASE ('Hb_PD')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hb_PD', long_name = 'Present-day bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs_PD')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hs_PD', long_name = 'Present-day surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL_PD')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'SL_PD', long_name = 'Present-day geoid elevation', units = 'm w.r.t. PD sea level')

      ! GIA equilibrium ice-sheet geometry
      CASE ('Hi_GIAeq')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hi_GIAeq', long_name = 'GIA equilibrium ice thickness', units = 'm')
      CASE ('Hb_GIAeq')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hb_GIAeq', long_name = 'GIA equilibrium bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs_GIAeq')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'Hs_GIAeq', long_name = 'GIA equilibrium surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL_GIAeq')
        CALL add_field_mesh_dp_2D_notime( filename, ncid, 'SL_GIAeq', long_name = 'GIA equilibrium geoid elevation', units = 'm w.r.t. PD sea level')

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      CASE ('Hi')
        CALL add_field_mesh_dp_2D( filename, ncid, 'Hi', long_name = 'Ice thickness', units = 'm')
      CASE ('Hb')
        CALL add_field_mesh_dp_2D( filename, ncid, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs')
        CALL add_field_mesh_dp_2D( filename, ncid, 'Hs', long_name = 'Surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hib')
        CALL add_field_mesh_dp_2D( filename, ncid, 'Hib', long_name = 'Ice base elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL')
        CALL add_field_mesh_dp_2D( filename, ncid, 'SL', long_name = 'Geoid elevation', units = 'm w.r.t. PD sea level')
      CASE ('TAF')
        CALL add_field_mesh_dp_2D( filename, ncid, 'TAF', long_name = 'Thickness above floatation', units = 'm')
      CASE ('Hi_eff')
        CALL add_field_mesh_dp_2D( filename, ncid, 'Hi_eff', long_name = 'Effective ice thickness', units = 'm')
      CASE ('Hs_slope')
        CALL add_field_mesh_dp_2D( filename, ncid, 'Hs_slope', long_name = 'Absolute surface gradient', units = '-')

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      CASE ('dHi')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHi', long_name = 'Ice thickness difference w.r.t. reference', units = 'm')
      CASE ('dHb')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHb', long_name = 'Bedrock elevation difference w.r.t. reference', units = 'm')
      CASE ('dHs')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHs', long_name = 'Surface elevation difference w.r.t. reference', units = 'm')
      CASE ('dHib')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHib', long_name = 'Ice base elevation difference w.r.t. reference', units = 'm')

    ! ===== Geometry rates of change =====
    ! ====================================

      CASE ('dHi_dt')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHi_dt', long_name = 'Ice thickness rate of change', units = 'm yr^-1')
      CASE ('dHb_dt')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHb_dt', long_name = 'Bedrock elevation rate of change', units = 'm yr^-1')
      CASE ('dHs_dt')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHs_dt', long_name = 'Surface elevation rate of change', units = 'm yr^-1')
      CASE ('dHib_dt')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHib_dt', long_name = 'Ice base elevation rate of change', units = 'm yr^-1')
      CASE ('dHi_dt_raw')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHi_dt_raw', long_name = 'Ice thickness rate of change before any modifications', units = 'm yr^-1')
      CASE ('dHi_dt_residual')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHi_dt_residual', long_name = 'Residual ice thickness rate of change during model calibration', units = 'm yr^-1')

    ! ===== Target quantities =====
    ! =============================

      CASE ('dHi_dt_target')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHi_dt_target', long_name = 'Target ice thickness rate of change during model calibration', units = 'm yr^-1')
      CASE ('uabs_surf_target')
        CALL add_field_mesh_dp_2D( filename, ncid, 'uabs_surf_target', long_name = 'Target ice surface speed during model calibration', units = 'm yr^-1')

    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('mask_icefree_land')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_icefree_land', long_name = 'Mask indicating ice-free land')
      CASE ('mask_icefree_ocean')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_icefree_ocean', long_name = 'Mask indicating ice-free ocean')
      CASE ('mask_grounded_ice')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_grounded_ice', long_name = 'Mask indicating grounded ice')
      CASE ('mask_floating_ice')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_floating_ice', long_name = 'Mask indicating floating ice')
      CASE ('mask_margin')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_margin', long_name = 'Mask indicating ice next to ice-free')
      CASE ('mask_gl_gr')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_gl_gr', long_name = 'Mask indicating grounded side of grounding line')
      CASE ('mask_gl_fl')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_gl_fl', long_name = 'Mask indicating floating side of grounding line')
      CASE ('mask_cf_gr')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_cf_gr', long_name = 'Mask indicating grounded calving front')
      CASE ('mask_cf_fl')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_cf_fl', long_name = 'Mask indicating floating calving front')
      CASE ('mask_coastline')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_coastline', long_name = 'Mask indicating ice-free land next to ice-free ocean')
      CASE ('mask')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask', long_name = 'General mask')
      CASE ('basin_ID')
        CALL add_field_mesh_int_2D( filename, ncid, 'basin_ID', long_name = 'Drainage basin ID', units = 'ID code')

    ! ===== Area fractions =====
    ! ==========================

      CASE ('fraction_gr')
        CALL add_field_mesh_dp_2D( filename, ncid, 'fraction_gr', long_name = 'Grounded area fractions of vertices', units = '0-1')
      CASE ('fraction_gr_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'fraction_gr_b', long_name = 'Grounded area fractions of triangles', units = '0-1')
      CASE ('fraction_margin')
        CALL add_field_mesh_dp_2D( filename, ncid, 'fraction_margin', long_name = 'Ice-covered area fractions of ice margins', units = '0-1')

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL add_field_mesh_dp_3D( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')
      CASE ('Ti_pmp')
        CALL add_field_mesh_dp_3D( filename, ncid, 'Ti_pmp', long_name = 'Pressure melting point temperature', units = 'K')
      CASE ('Ti_hom')
        CALL add_field_mesh_dp_2D( filename, ncid, 'Ti_hom', long_name = 'Temperature at base w.r.t. pressure melting point', units = 'K')
      CASE ('Cpi')
        CALL add_field_mesh_dp_3D( filename, ncid, 'Cpi', long_name = 'Specific heat capacity', units = 'J kg^-1 K^-1')
      CASE ('Ki')
        CALL add_field_mesh_dp_3D( filename, ncid, 'Ki', long_name = 'Thermal conductivity', units = 'J m^-1 K^-1 yr^-1')
      CASE ('internal_heating')
        CALL add_field_mesh_dp_3D( filename, ncid, 'internal_heating', long_name = 'Internal heating', units = '?')
      CASE ('frictional_heating')
        CALL add_field_mesh_dp_2D( filename, ncid, 'frictional_heating', long_name = 'Frictional heating', units = '?')
      CASE ('A_flow')
        CALL add_field_mesh_dp_3D( filename, ncid, 'A_flow', long_name = 'Glens flow law factor', units = 'Pa^-3 y^-1')

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      CASE ('u_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'u_3D', long_name = '3-D ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'v_3D', long_name = '3-D ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_3D_b')
        CALL add_field_mesh_dp_3D_b( filename, ncid, 'u_3D_b', long_name = '3-D ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_3D_b')
        CALL add_field_mesh_dp_3D_b( filename, ncid, 'v_3D_b', long_name = '3-D ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('w_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'w_3D', long_name = '3-D ice velocity in the z-direction', units = 'm yr^-1')

      ! Vertically integrated
      CASE ('u_vav')
        CALL add_field_mesh_dp_2D( filename, ncid, 'u_vav', long_name = 'Vertically averaged ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_vav')
        CALL add_field_mesh_dp_2D( filename, ncid, 'v_vav', long_name = 'Vertically averaged ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_vav_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'u_vav_b', long_name = 'Vertically averaged ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_vav_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'v_vav_b', long_name = 'Vertically averaged ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('uabs_vav')
        CALL add_field_mesh_dp_2D( filename, ncid, 'uabs_vav', long_name = 'Vertically averaged absolute ice velocity', units = 'm yr^-1')
      CASE ('uabs_vav_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'uabs_vav_b', long_name = 'Vertically averaged absolute ice velocity', units = 'm yr^-1')

      ! Surface
      CASE ('u_surf')
        CALL add_field_mesh_dp_2D( filename, ncid, 'u_surf', long_name = 'Surface ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_surf')
        CALL add_field_mesh_dp_2D( filename, ncid, 'v_surf', long_name = 'Surface ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_surf_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'u_surf_b', long_name = 'Surface ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_surf_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'v_surf_b', long_name = 'Surface ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('w_surf')
        CALL add_field_mesh_dp_2D( filename, ncid, 'w_surf', long_name = 'Surface ice velocity in the z-direction', units = 'm yr^-1')
      CASE ('uabs_surf')
        CALL add_field_mesh_dp_2D( filename, ncid, 'uabs_surf', long_name = 'Absolute surface ice velocity', units = 'm yr^-1')
      CASE ('uabs_surf_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'uabs_surf_b', long_name = 'Absolute surface ice velocity', units = 'm yr^-1')

      ! Base
      CASE ('u_base')
        CALL add_field_mesh_dp_2D( filename, ncid, 'u_base', long_name = 'Basal ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_base')
        CALL add_field_mesh_dp_2D( filename, ncid, 'v_base', long_name = 'Basal ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_base_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'u_base_b', long_name = 'Basal ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_base_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'v_base_b', long_name = 'Basal ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('w_base')
        CALL add_field_mesh_dp_2D( filename, ncid, 'w_base', long_name = 'Basal ice velocity in the z-direction', units = 'm yr^-1')
      CASE ('uabs_base')
        CALL add_field_mesh_dp_2D( filename, ncid, 'uabs_base', long_name = 'Absolute basal ice velocity', units = 'm yr^-1')
      CASE ('uabs_base_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'uabs_base_b', long_name = 'Absolute basal ice velocity', units = 'm yr^-1')

    ! === Strain rates ===
    ! ====================

      CASE ('du_dx_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'du_dx_3D', long_name = '3-D xx strain rate', units = 'yr^-1')
      CASE ('du_dy_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'du_dy_3D', long_name = '3-D xy strain rate', units = 'yr^-1')
      CASE ('du_dz_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'du_dz_3D', long_name = '3-D xz strain rate', units = 'yr^-1')
      CASE ('dv_dx_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'dv_dx_3D', long_name = '3-D yx strain rate', units = 'yr^-1')
      CASE ('dv_dy_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'dv_dy_3D', long_name = '3-D yy strain rate', units = 'yr^-1')
      CASE ('dv_dz_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'dv_dz_3D', long_name = '3-D yz strain rate', units = 'yr^-1')
      CASE ('dw_dx_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'dw_dx_3D', long_name = '3-D zx strain rate', units = 'yr^-1')
      CASE ('dw_dy_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'dw_dy_3D', long_name = '3-D zy strain rate', units = 'yr^-1')
      CASE ('dw_dz_3D')
        CALL add_field_mesh_dp_3D( filename, ncid, 'dw_dz_3D', long_name = '3-D zz strain rate', units = 'yr^-1')

    ! == Ice flow regime ==
    ! =====================

      CASE ('divQ')
        CALL add_field_mesh_dp_2D( filename, ncid, 'divQ', long_name = 'Horizontal ice flux divergence', units = 'm yr^-1')
      CASE ('R_shear')
        CALL add_field_mesh_dp_2D( filename, ncid, 'R_shear', long_name = 'Slide/shear ratio', units = '0-1')

    ! == Ice P/C time stepping ==
    ! ===========================

      CASE ('pc_truncation_error')
        CALL add_field_mesh_dp_2D( filename, ncid, 'pc_truncation_error', long_name = 'Ice P/C truncation error tau', units = 'm')
      CASE ('pc_untolerated_events')
        CALL add_field_mesh_int_2D( filename, ncid, 'pc_untolerated_events', long_name = 'Ice P/C number of events above error tolerance', units = '-')

    ! == Basal hydrology ==
    ! =====================

      CASE ('pore_water_pressure')
        CALL add_field_mesh_dp_2D( filename, ncid, 'pore_water_pressure', long_name = 'Till pore water pressure', units = 'Pa')
      CASE ('overburden_pressure')
        CALL add_field_mesh_dp_2D( filename, ncid, 'overburden_pressure', long_name = 'Ice overburden pressure', units = 'Pa')
      CASE ('effective_pressure')
        CALL add_field_mesh_dp_2D( filename, ncid, 'effective_pressure', long_name = 'Effective basal pressure', units = 'Pa')
      CASE ('pore_water_likelihood')
        CALL add_field_mesh_dp_2D( filename, ncid, 'pore_water_likelihood', long_name = 'Till pore water likelihood', units = '0-1')
      CASE ('pore_water_fraction')
        CALL add_field_mesh_dp_2D( filename, ncid, 'pore_water_fraction', long_name = 'Fraction of overburden pressure reduced by pore water ', units = '0-1')

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      CASE ('till_friction_angle')
        CALL add_field_mesh_dp_2D( filename, ncid, 'till_friction_angle', long_name = 'Till friction angle', units = 'degrees')
      CASE ('bed_roughness')
        CALL add_field_mesh_dp_2D( filename, ncid, 'bed_roughness', long_name = 'Bed roughness', units = '0-1')
      CASE ('till_yield_stress')
        CALL add_field_mesh_dp_2D( filename, ncid, 'till_yield_stress', long_name = 'Till yield stress', units = 'Pa')
      CASE ('slid_alpha_sq')
        CALL add_field_mesh_dp_2D( filename, ncid, 'slid_alpha_sq', long_name = 'Coulomb-law friction coefficientn', units = 'dimensionless')
      CASE ('slid_beta_sq')
        CALL add_field_mesh_dp_2D( filename, ncid, 'slid_beta_sq', long_name = 'Power-law friction coefficient', units = 'Pa m^−1/m yr^1/m')

      ! Basal friction and shear stress
      CASE ('basal_friction_coefficient')
        CALL add_field_mesh_dp_2D( filename, ncid, 'basal_friction_coefficient', long_name = 'Basal friction coefficient', units = 'Pa yr m^-1')
      CASE ('basal_shear_stress')
        CALL add_field_mesh_dp_2D( filename, ncid, 'basal_shear_stress', long_name = 'Basal shear stress', units = 'Pa')

    ! == Geothermal heat ==
    ! =====================

      CASE ('geothermal_heat_flux')
        CALL add_field_mesh_dp_2D( filename, ncid, 'geothermal_heat_flux', long_name = 'Geothermal heat flux', units = 'J m^-2 yr^-1')

    ! == Climate ==
    ! =============

      ! Main climate variables
      CASE ('T2m')
        CALL add_field_mesh_dp_2D_monthly( filename, ncid, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
      CASE ('Precip')
        CALL add_field_mesh_dp_2D_monthly( filename, ncid, 'Precip', long_name = 'Monthly total precipitation', units = 'm.w.e.')

    ! == Ocean ==
    ! ===========

      ! Main ocean variables
      CASE ('T_ocean')
        CALL warning('ocean temperature not implemented yet!')
      CASE ('S_ocean')
        CALL warning('ocean salinity not implemented yet!')

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      CASE ('SMB')
        CALL add_field_mesh_dp_2D( filename, ncid, 'SMB', long_name = 'Surface mass balance', units = 'm yr^-1')

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      CASE ('BMB')
        CALL add_field_mesh_dp_2D( filename, ncid, 'BMB', long_name = 'Basal mass balance', units = 'm yr^-1')

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      CASE ('H_lad')
        CALL add_field_mesh_dp_2D( filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
      CASE ('U_lad')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
      CASE ('V_lad')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
      CASE ('T_lad')
        CALL add_field_mesh_dp_2D( filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
      CASE ('S_lad')
        CALL add_field_mesh_dp_2D( filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      CASE ('LMB')
        CALL add_field_mesh_dp_2D( filename, ncid, 'LMB', long_name = 'Lateral mass balance', units = 'm yr^-1')

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      CASE ('AMB')
        CALL add_field_mesh_dp_2D( filename, ncid, 'AMB', long_name = 'Artificial mass balance', units = 'm yr^-1')

    ! ===== End of user-defined output fields =====
    ! =============================================

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_output_field "' // TRIM( choice_output_field) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_main_regional_output_file_mesh_field

  SUBROUTINE create_main_regional_output_file_grid_field( filename, ncid, choice_output_field)
    ! Create the main regional output NetCDF file - grid version
    !
    ! Add a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: choice_output_field

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'create_main_regional_output_file_grid_field'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Add the specified data field to the file
    SELECT CASE (choice_output_field)
      CASE ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      CASE ('resolution')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'resolution', long_name = 'Mesh resolution (distance to nearest neighbour)', units = 'm')

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      CASE ('Hi_init')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hi_init', long_name = 'Initial ice thickness', units = 'm')
      CASE ('Hb_init')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hb_init', long_name = 'Initial bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs_init')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hs_init', long_name = 'Initial surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL_init')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'SL_init', long_name = 'Initial geoid elevation', units = 'm w.r.t. PD sea level')

      ! Present-day ice-sheet geometry
      CASE ('Hi_PD')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hi_PD', long_name = 'Present-day ice thickness', units = 'm')
      CASE ('Hb_PD')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hb_PD', long_name = 'Present-day bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs_PD')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hs_PD', long_name = 'Present-day surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL_PD')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'SL_PD', long_name = 'Present-day geoid elevation', units = 'm w.r.t. PD sea level')

      ! GIA equilibrium ice-sheet geometry
      CASE ('Hi_GIAeq')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hi_GIAeq', long_name = 'GIA equilibrium ice thickness', units = 'm')
      CASE ('Hb_GIAeq')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hb_GIAeq', long_name = 'GIA equilibrium bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs_GIAeq')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'Hs_GIAeq', long_name = 'GIA equilibrium surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL_GIAeq')
        CALL add_field_grid_dp_2D_notime( filename, ncid, 'SL_GIAeq', long_name = 'GIA equilibrium geoid elevation', units = 'm w.r.t. PD sea level')

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      CASE ('Hi')
        CALL add_field_grid_dp_2D( filename, ncid, 'Hi', long_name = 'Ice thickness', units = 'm')
      CASE ('Hb')
        CALL add_field_grid_dp_2D( filename, ncid, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hs')
        CALL add_field_grid_dp_2D( filename, ncid, 'Hs', long_name = 'Surface elevation', units = 'm w.r.t. PD sea level')
      CASE ('Hib')
        CALL add_field_grid_dp_2D( filename, ncid, 'Hib', long_name = 'Ice base elevation', units = 'm w.r.t. PD sea level')
      CASE ('SL')
        CALL add_field_grid_dp_2D( filename, ncid, 'SL', long_name = 'Geoid elevation', units = 'm w.r.t. PD sea level')
      CASE ('TAF')
        CALL add_field_grid_dp_2D( filename, ncid, 'TAF', long_name = 'Thickness above floatation', units = 'm')
      CASE ('Hi_eff')
        CALL add_field_grid_dp_2D( filename, ncid, 'Hi_eff', long_name = 'Effective ice thickness', units = 'm')
      CASE ('Hs_slope')
        CALL add_field_grid_dp_2D( filename, ncid, 'Hs_slope', long_name = 'Absolute surface gradient', units = '-')

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      CASE ('dHi')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHi', long_name = 'Ice thickness difference w.r.t. reference', units = 'm')
      CASE ('dHb')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHb', long_name = 'Bedrock elevation difference w.r.t. reference', units = 'm')
      CASE ('dHs')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHs', long_name = 'Surface elevation difference w.r.t. reference', units = 'm')
      CASE ('dHib')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHib', long_name = 'Ice base elevation difference w.r.t. reference', units = 'm')

    ! ===== Geometry rates of change =====
    ! ====================================

      CASE ('dHi_dt')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHi_dt', long_name = 'Ice thickness rate of change', units = 'm yr^-1')
      CASE ('dHb_dt')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHb_dt', long_name = 'Bedrock elevation rate of change', units = 'm yr^-1')
      CASE ('dHs_dt')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHs_dt', long_name = 'Surface elevation rate of change', units = 'm yr^-1')
      CASE ('dHib_dt')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHib_dt', long_name = 'Ice base elevation rate of change', units = 'm yr^-1')
      CASE ('dHi_dt_raw')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHi_dt_raw', long_name = 'Ice thickness rate of change before any modifications', units = 'm yr^-1')
      CASE ('dHi_dt_residual')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHi_dt_residual', long_name = 'Residual ice thickness rate of change during model calibration', units = 'm yr^-1')

    ! ===== Target quantities =====
    ! =============================

      CASE ('dHi_dt_target')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHi_dt_target', long_name = 'Target ice thickness rate of change during model calibration', units = 'm yr^-1')
      CASE ('uabs_surf_target')
        CALL add_field_grid_dp_2D( filename, ncid, 'uabs_surf_target', long_name = 'Target ice surface speed during model calibration', units = 'm yr^-1')

    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('mask_icefree_land')
      CASE ('mask_icefree_ocean')
      CASE ('mask_grounded_ice')
      CASE ('mask_floating_ice')
      CASE ('mask_margin')
      CASE ('mask_gl_gr')
      CASE ('mask_gl_fl')
      CASE ('mask_cf_gr')
      CASE ('mask_cf_fl')
      CASE ('mask_coastline')
      CASE ('mask')
      CASE ('basin_ID')

    ! ===== Area fractions =====
    ! ==========================

      ! NOTE: sub-grid area fractions cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('fraction_gr')
      CASE ('fraction_gr_b')
      CASE ('fraction_margin')

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL add_field_grid_dp_3D( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')
      CASE ('Ti_pmp')
        CALL add_field_grid_dp_3D( filename, ncid, 'Ti_pmp', long_name = 'Pressure melting point temperature', units = 'K')
      CASE ('Ti_hom')
        CALL add_field_grid_dp_2D( filename, ncid, 'Ti_hom', long_name = 'Temperature at base w.r.t. pressure melting point', units = 'K')
      CASE ('Cpi')
        CALL add_field_grid_dp_3D( filename, ncid, 'Cpi', long_name = 'Specific heat capacity', units = 'J kg^-1 K^-1')
      CASE ('Ki')
        CALL add_field_grid_dp_3D( filename, ncid, 'Ki', long_name = 'Thermal conductivity', units = 'J m^-1 K^-1 yr^-1')
      CASE ('internal_heating')
        CALL add_field_grid_dp_3D( filename, ncid, 'internal_heating', long_name = 'Internal heating', units = '?')
      CASE ('frictional_heating')
        CALL add_field_grid_dp_2D( filename, ncid, 'frictional_heating', long_name = 'Frictional heating', units = '?')
      CASE ('A_flow')
        CALL add_field_grid_dp_3D( filename, ncid, 'A_flow', long_name = 'Glens flow law factor', units = 'Pa^-3 y^-1')

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      CASE ('u_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'u_3D', long_name = '3-D ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'v_3D', long_name = '3-D ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_3D_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_3D_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'w_3D', long_name = '3-D ice velocity in the z-direction', units = 'm yr^-1')

      ! Vertically integrated
      CASE ('u_vav')
        CALL add_field_grid_dp_2D( filename, ncid, 'u_vav', long_name = 'Vertically averaged ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_vav')
        CALL add_field_grid_dp_2D( filename, ncid, 'v_vav', long_name = 'Vertically averaged ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('uabs_vav')
        CALL add_field_grid_dp_2D( filename, ncid, 'uabs_vav', long_name = 'Vertically averaged absolute ice velocity', units = 'm yr^-1')
      CASE ('uabs_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

      ! Surface
      CASE ('u_surf')
        CALL add_field_grid_dp_2D( filename, ncid, 'u_surf', long_name = 'Surface ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_surf')
        CALL add_field_grid_dp_2D( filename, ncid, 'v_surf', long_name = 'Surface ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_surf')
        CALL add_field_grid_dp_2D( filename, ncid, 'w_surf', long_name = 'Surface ice velocity in the z-direction', units = 'm yr^-1')
      CASE ('uabs_surf')
        CALL add_field_grid_dp_2D( filename, ncid, 'uabs_surf', long_name = 'Absolute surface ice velocity', units = 'm yr^-1')
      CASE ('uabs_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

      ! Base
      CASE ('u_base')
        CALL add_field_grid_dp_2D( filename, ncid, 'u_base', long_name = 'Basal ice velocity in the x-direction', units = 'm yr^-1')
      CASE ('v_base')
        CALL add_field_grid_dp_2D( filename, ncid, 'v_base', long_name = 'Basal ice velocity in the y-direction', units = 'm yr^-1')
      CASE ('u_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_base')
        CALL add_field_grid_dp_2D( filename, ncid, 'w_base', long_name = 'Basal ice velocity in the z-direction', units = 'm yr^-1')
      CASE ('uabs_base')
        CALL add_field_grid_dp_2D( filename, ncid, 'uabs_base', long_name = 'Absolute basal ice velocity', units = 'm yr^-1')
      CASE ('uabs_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

    ! === Strain rates ===
    ! ====================

      CASE ('du_dx_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'du_dx_3D', long_name = '3-D xx strain rate', units = 'yr^-1')
      CASE ('du_dy_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'du_dy_3D', long_name = '3-D xy strain rate', units = 'yr^-1')
      CASE ('du_dz_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'du_dz_3D', long_name = '3-D xz strain rate', units = 'yr^-1')
      CASE ('dv_dx_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'dv_dx_3D', long_name = '3-D yx strain rate', units = 'yr^-1')
      CASE ('dv_dy_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'dv_dy_3D', long_name = '3-D yy strain rate', units = 'yr^-1')
      CASE ('dv_dz_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'dv_dz_3D', long_name = '3-D yz strain rate', units = 'yr^-1')
      CASE ('dw_dx_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'dw_dx_3D', long_name = '3-D zx strain rate', units = 'yr^-1')
      CASE ('dw_dy_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'dw_dy_3D', long_name = '3-D zy strain rate', units = 'yr^-1')
      CASE ('dw_dz_3D')
        CALL add_field_grid_dp_3D( filename, ncid, 'dw_dz_3D', long_name = '3-D zz strain rate', units = 'yr^-1')

    ! == Ice flow regime ==
    ! =====================

      CASE ('divQ')
        CALL add_field_grid_dp_2D( filename, ncid, 'divQ', long_name = 'Horizontal ice flux divergence', units = 'm yr^-1')
      CASE ('R_shear')
        CALL add_field_grid_dp_2D( filename, ncid, 'R_shear', long_name = 'Slide/shear ratio', units = '0-1')

    ! == Ice P/C time stepping ==
    ! ===========================

      CASE ('pc_truncation_error')
        CALL add_field_grid_dp_2D( filename, ncid, 'pc_truncation_error', long_name = 'Ice P/C truncation error tau', units = 'm')
      CASE ('pc_untolerated_events')
        ! DENK DROM : not gridable

    ! == Basal hydrology ==
    ! =====================

      CASE ('pore_water_pressure')
        CALL add_field_grid_dp_2D( filename, ncid, 'pore_water_pressure', long_name = 'Till pore water pressure', units = 'Pa')
      CASE ('overburden_pressure')
        CALL add_field_grid_dp_2D( filename, ncid, 'overburden_pressure', long_name = 'Ice overburden pressure', units = 'Pa')
      CASE ('effective_pressure')
        CALL add_field_grid_dp_2D( filename, ncid, 'effective_pressure', long_name = 'Effective basal pressure', units = 'Pa')
      CASE ('pore_water_likelihood')
        CALL add_field_grid_dp_2D( filename, ncid, 'pore_water_likelihood', long_name = 'Till pore water likelihood', units = '0-1')
      CASE ('pore_water_fraction')
        CALL add_field_grid_dp_2D( filename, ncid, 'pore_water_fraction', long_name = 'Fraction of overburden pressure reduced by pore water', units = '0-1')

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      CASE ('till_friction_angle')
        CALL add_field_grid_dp_2D( filename, ncid, 'till_friction_angle', long_name = 'Till friction angle', units = 'degrees')
      CASE ('bed_roughness')
        CALL add_field_grid_dp_2D( filename, ncid, 'bed_roughness', long_name = 'Bed roughness', units = '0-1')
      CASE ('till_yield_stress')
        CALL add_field_grid_dp_2D( filename, ncid, 'till_yield_stress', long_name = 'Till yield stress', units = 'Pa')
      CASE ('slid_alpha_sq')
        CALL add_field_grid_dp_2D( filename, ncid, 'slid_alpha_sq', long_name = 'Coulomb-law friction coefficientn', units = 'dimensionless')
      CASE ('slid_beta_sq')
        CALL add_field_grid_dp_2D( filename, ncid, 'slid_beta_sq', long_name = 'Power-law friction coefficient', units = 'Pa m^−1/m yr^1/m')

      ! Basal friction and shear stress
      CASE ('basal_friction_coefficient')
        CALL add_field_grid_dp_2D( filename, ncid, 'basal_friction_coefficient', long_name = 'Basal friction coefficient', units = 'Pa yr m^-1')
      CASE ('basal_shear_stress')
        CALL add_field_grid_dp_2D( filename, ncid, 'basal_shear_stress', long_name = 'Basal shear stress', units = 'Pa')

    ! == Geothermal heat ==
    ! =====================

      CASE ('geothermal_heat_flux')
        CALL add_field_grid_dp_2D( filename, ncid, 'geothermal_heat_flux', long_name = 'Geothermal heat flux', units = 'J m^-2 yr^-1')

    ! == Climate ==
    ! =============

      ! Main climate variables
      CASE ('T2m')
        CALL add_field_grid_dp_2D_monthly( filename, ncid, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
      CASE ('Precip')
        CALL add_field_grid_dp_2D_monthly( filename, ncid, 'Precip', long_name = 'Monthly total precipitation', units = 'm.w.e.')

    ! == Ocean ==
    ! ==========================

      ! Main ocean variables
      CASE ('T_ocean')
        CALL warning('ocean temperature not implemented yet!')
      CASE ('S_ocean')
        CALL warning('ocean salinity not implemented yet!')

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      CASE ('SMB')
        CALL add_field_grid_dp_2D( filename, ncid, 'SMB', long_name = 'Surface mass balance', units = 'm yr^-1')

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      CASE ('BMB')
        CALL add_field_grid_dp_2D( filename, ncid, 'BMB', long_name = 'Basal mass balance', units = 'm yr^-1')

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      CASE ('H_lad')
        CALL add_field_grid_dp_2D( filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
      CASE ('U_lad')
        CALL add_field_grid_dp_2D( filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
      CASE ('V_lad')
        CALL add_field_grid_dp_2D( filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
      CASE ('T_lad')
        CALL add_field_grid_dp_2D( filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
      CASE ('S_lad')
        CALL add_field_grid_dp_2D( filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      CASE ('LMB')
        CALL add_field_grid_dp_2D( filename, ncid, 'LMB', long_name = 'Lateral mass balance', units = 'm yr^-1')

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      CASE ('AMB')
        CALL add_field_grid_dp_2D( filename, ncid, 'AMB', long_name = 'Artificial mass balance', units = 'm yr^-1')

    ! ===== End of user-defined output fields =====
    ! =============================================

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_output_field "' // TRIM( choice_output_field) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_main_regional_output_file_grid_field

  SUBROUTINE create_scalar_regional_output_file( region)
    ! Create the scalar regional output NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'create_scalar_regional_output_file'
    CHARACTER(LEN=256)                                                 :: filename_base
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'scalar_output_' // region%name
    CALL generate_filename_XXXXXdotnc( filename_base, region%output_filename_scalar)

    ! Print to terminal
    IF (par%master) WRITE(0,'(A)') '   Creating scalar output file "' // colour_string( TRIM( region%output_filename_scalar), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( region%output_filename_scalar, ncid)

    ! Add time, zeta, and month dimensions+variables to the file
    CALL add_time_dimension_to_file(  region%output_filename_scalar, ncid)

    ! Add the default data fields to the file
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_area')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume_af')

    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_area_PD')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume_PD')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'ice_volume_af_PD')

    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_total')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_gr')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_fl')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_land')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'SMB_ocean')

    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_total')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_gr')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_fl')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_land')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'BMB_ocean')

    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'LMB_total')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'LMB_gr')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'LMB_fl')

    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_total')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_gr')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_fl')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_land')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'AMB_ocean')

    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'gl_flux')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'cf_gr_flux')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'cf_fl_flux')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'margin_land_flux')
    CALL create_scalar_regional_output_file_field( region%output_filename_scalar, ncid, 'margin_ocean_flux')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_scalar_regional_output_file

  SUBROUTINE create_scalar_regional_output_file_field( filename, ncid, choice_output_field)
    ! Create the main regional output NetCDF file - mesh version
    !
    ! Add a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: choice_output_field

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'create_scalar_regional_output_file_field'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Add the specified data field to the file
    SELECT CASE (choice_output_field)
      CASE ('none')
        ! Do nothing

      ! Total ice sheet area
      CASE ('ice_area')
        CALL add_field_dp_0D( filename, ncid, 'ice_area', long_name = 'Total ice area', units = 'm^2')

        ! Total ice sheet volume in metres of sea level equivalent
      CASE ('ice_volume')
        CALL add_field_dp_0D( filename, ncid, 'ice_volume', long_name = 'Total ice volume', units = 'm s.l.e.')

      ! Total ice sheet volume above floatation in metres of sea level equivalent
      CASE ('ice_volume_af')
        CALL add_field_dp_0D( filename, ncid, 'ice_volume_af', long_name = 'Total ice volume above floatation', units = 'm s.l.e.')

      ! Total ice sheet area for present-day
      CASE ('ice_area_PD')
        CALL add_field_dp_0D( filename, ncid, 'ice_area_PD', long_name = 'Total ice area for present-day', units = 'm^2')

        ! Total ice sheet volume in metres of sea level equivalent for present-day
      CASE ('ice_volume_PD')
        CALL add_field_dp_0D( filename, ncid, 'ice_volume_PD', long_name = 'Total ice volume for present-day', units = 'm s.l.e.')

      ! Total ice sheet volume above floatation in metres of sea level equivalent for present-day
      CASE ('ice_volume_af_PD')
        CALL add_field_dp_0D( filename, ncid, 'ice_volume_af_PD', long_name = 'Total ice volume above floatation for present-day', units = 'm s.l.e.')

      ! Total SMB integrated over the entire domain
      CASE ('SMB_total')
        CALL add_field_dp_0D( filename, ncid, 'SMB_total', long_name = 'Area-integrated total SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice sheet area
      CASE ('SMB_gr')
        CALL add_field_dp_0D( filename, ncid, 'SMB_gr', long_name = 'Area-integrated ice sheet SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice shelf area
      CASE ('SMB_fl')
        CALL add_field_dp_0D( filename, ncid, 'SMB_fl', long_name = 'Area-integrated ice shelf SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice-free land area
      CASE ('SMB_land')
        CALL add_field_dp_0D( filename, ncid, 'SMB_land', long_name = 'Area-integrated ice-free land SMB', units = 'Gt yr^-1')

      ! Total SMB integrated over the entire ice-free ocean area
      CASE ('SMB_ocean')
        CALL add_field_dp_0D( filename, ncid, 'SMB_ocean', long_name = 'Area-integrated ice-free ocean SMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire domain
      CASE ('BMB_total')
        CALL add_field_dp_0D( filename, ncid, 'BMB_total', long_name = 'Area-integrated total BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice sheet area
      CASE ('BMB_gr')
        CALL add_field_dp_0D( filename, ncid, 'BMB_gr', long_name = 'Area-integrated ice sheet BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice shelf area
      CASE ('BMB_fl')
        CALL add_field_dp_0D( filename, ncid, 'BMB_fl', long_name = 'Area-integrated ice shelf BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice-free land area
      CASE ('BMB_land')
        CALL add_field_dp_0D( filename, ncid, 'BMB_land', long_name = 'Area-integrated ice-free land BMB', units = 'Gt yr^-1')

      ! Total BMB integrated over the entire ice-free ocean area
      CASE ('BMB_ocean')
        CALL add_field_dp_0D( filename, ncid, 'BMB_ocean', long_name = 'Area-integrated ice-free ocean BMB', units = 'Gt yr^-1')

      ! Total LMB integrated over the entire domain
      CASE ('LMB_total')
        CALL add_field_dp_0D( filename, ncid, 'LMB_total', long_name = 'Area-integrated total LMB', units = 'Gt yr^-1')

      ! Total LMB integrated over the entire ice sheet area
      CASE ('LMB_gr')
        CALL add_field_dp_0D( filename, ncid, 'LMB_gr', long_name = 'Area-integrated ice sheet LMB', units = 'Gt yr^-1')

      ! Total LMB integrated over the entire ice shelf area
      CASE ('LMB_fl')
        CALL add_field_dp_0D( filename, ncid, 'LMB_fl', long_name = 'Area-integrated ice shelf LMB', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire domain
      CASE ('AMB_total')
        CALL add_field_dp_0D( filename, ncid, 'AMB_total', long_name = 'Area-integrated total additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice sheet area
      CASE ('AMB_gr')
        CALL add_field_dp_0D( filename, ncid, 'AMB_gr', long_name = 'Area-integrated ice sheet additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice shelf area
      CASE ('AMB_fl')
        CALL add_field_dp_0D( filename, ncid, 'AMB_fl', long_name = 'Area-integrated ice shelf additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice-free land area
      CASE ('AMB_land')
        CALL add_field_dp_0D( filename, ncid, 'AMB_land', long_name = 'Area-integrated ice-free land additional MB from other sources', units = 'Gt yr^-1')

      ! Total additional MB from other sources integrated over the entire ice-free ocean area
      CASE ('AMB_ocean')
        CALL add_field_dp_0D( filename, ncid, 'AMB_ocean', long_name = 'Area-integrated ice-free ocean additional MB from other sources', units = 'Gt yr^-1')

      ! Total flux through the grounding line
      CASE ('gl_flux')
        CALL add_field_dp_0D( filename, ncid, 'gl_flux', long_name = 'Total lateral grounding line flux', units = 'Gt yr^-1')

      ! Total flux through grounded calving fronts
      CASE ('cf_gr_flux')
        CALL add_field_dp_0D( filename, ncid, 'cf_gr_flux', long_name = 'Total lateral grounded calving front flux', units = 'Gt yr^-1')

      ! Total flux through floating calving fronts
      CASE ('cf_fl_flux')
        CALL add_field_dp_0D( filename, ncid, 'cf_fl_flux', long_name = 'Total lateral floating calving front flux', units = 'Gt yr^-1')

      ! Total flux exiting ice margins into grounded areas
      CASE ('margin_land_flux')
        CALL add_field_dp_0D( filename, ncid, 'margin_land_flux', long_name = 'Total lateral flux exiting the ice margin into ground', units = 'Gt yr^-1')

      ! Total flux exiting ice margins into marine areas
      CASE ('margin_ocean_flux')
        CALL add_field_dp_0D( filename, ncid, 'margin_ocean_flux', long_name = 'Total lateral flux exiting the ice margin into water', units = 'Gt yr^-1')

    ! ===== End of user-defined output fields =====
    ! =============================================

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_output_field "' // TRIM( choice_output_field) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_scalar_regional_output_file_field

END MODULE main_regional_output
