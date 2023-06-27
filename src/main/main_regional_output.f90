MODULE main_regional_output

  ! Creating and writing to the main regional output files

! ===== Preamble =====
! ====================

  ! USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: happy, warning, crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE region_types                                           , ONLY: type_model_region
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, open_existing_netcdf_file_for_writing, close_netcdf_file
  USE netcdf_output                                          , ONLY: generate_filename_XXXXXdotnc, setup_mesh_in_netcdf_file, setup_xy_grid_in_netcdf_file, &
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
                                                                     write_to_field_multopt_mesh_int_2D, write_to_field_multopt_mesh_int_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D, write_to_field_multopt_mesh_dp_2D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_monthly, write_to_field_multopt_mesh_dp_2D_monthly_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D, write_to_field_multopt_mesh_dp_3D_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_b, write_to_field_multopt_mesh_dp_2D_b_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_b, write_to_field_multopt_mesh_dp_3D_b_notime, &
                                                                     write_to_field_multopt_grid_dp_2D, write_to_field_multopt_grid_dp_2D_notime, &
                                                                     write_to_field_multopt_grid_dp_2D_monthly, write_to_field_multopt_grid_dp_2D_monthly_notime, &
                                                                     write_to_field_multopt_grid_dp_3D, write_to_field_multopt_grid_dp_3D_notime
  USE mesh_remapping                                         , ONLY: map_from_mesh_to_xy_grid_2D, map_from_mesh_to_xy_grid_3D

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

    ! Add all data fields to the file
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

    ! Add all data fields to the file
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_01)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_02)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_03)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_04)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_05)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_06)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_07)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_08)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_09)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_10)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_11)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_12)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_13)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_14)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_15)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_16)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_17)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_18)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_19)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_20)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_21)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_22)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_23)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_24)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_25)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_26)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_27)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_28)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_29)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_30)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_31)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_32)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_33)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_34)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_35)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_36)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_37)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_38)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_39)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_40)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_41)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_42)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_43)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_44)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_45)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_46)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_47)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_48)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_49)
    CALL write_to_main_regional_output_file_grid_field( region, region%output_filename_grid, ncid, C%choice_output_field_50)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_main_regional_output_file_grid

  SUBROUTINE write_to_main_regional_output_file_mesh_field( region, filename, ncid, choice_output_field)
    ! Write to the main regional output NetCDF file - mesh version
    !
    ! Write a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_output_field

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
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_icefree_ocean', mask_int)
      CASE ('mask_grounded_ice')
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_grounded_ice', mask_int)
      CASE ('mask_floating_ice')
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_floating_ice', mask_int)
      CASE ('mask_gl_gr')
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_gl_gr', mask_int)
      CASE ('mask_gl_fl')
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_gl_fl', mask_int)
      CASE ('mask_cf_gr')
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_cf_gr', mask_int)
      CASE ('mask_cf_fl')
        WHERE (region%ice%mask_icefree_land .EQV. .TRUE.)
          mask_int = 1
        ELSEWHERE
          mask_int = 0
        END WHERE
        CALL write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_cf_fl', mask_int)
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
      CASE ('fraction_cf')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'fraction_cf', region%ice%fraction_cf)

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ti', region%ice%Ti)
      CASE ('Ti_pmp')
        CALL write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ti_pmp', region%ice%Ti_pmp)
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

    ! == Basal conditions ==
    ! ======================

      ! Basal hydrology
      CASE ('pore_water_pressure')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pore_water_pressure', region%ice%pore_water_pressure)
      CASE ('overburden_pressure')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'overburden_pressure', region%ice%overburden_pressure)
      CASE ('effective_pressure')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'effective_pressure', region%ice%effective_pressure)

      ! Basal roughness / friction
      CASE ('phi_fric')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'phi_fric', region%ice%phi_fric)
      CASE ('tau_c')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'tau_c', region%ice%tau_c)
      CASE ('alpha_sq')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'alpha_sq', region%ice%alpha_sq)
      CASE ('beta_sq')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'beta_sq', region%ice%beta_sq)
      CASE ('beta_b')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'beta_b', region%ice%beta_b)

      ! Basal sliding
      CASE ('friction_coef_1')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'friction_coef_1', region%ice%friction_coef_1)
      CASE ('friction_coef_2')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'friction_coef_2', region%ice%friction_coef_2)
      CASE ('basal_shear_stress')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'basal_shear_stress', region%ice%basal_shear_stress)

      ! Geothermal heat
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
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'SMB', region%SMB%SMB)

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      CASE ('BMB')
        CALL write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'BMB', region%BMB%BMB)

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

  SUBROUTINE write_to_main_regional_output_file_grid_field( region, filename, ncid, choice_output_field)
    ! Write to the main regional output NetCDF file - grid version
    !
    ! Write a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region)                            , INTENT(IN)    :: region
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_output_field

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'write_to_main_regional_output_file_grid_field'
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
    ALLOCATE( d_grid_vec_partial_2D(         region%output_grid%n_loc                ))
    ALLOCATE( d_grid_vec_partial_2D_monthly( region%output_grid%n_loc, 12            ))
    ALLOCATE( d_grid_vec_partial_3D(         region%output_grid%n_loc, region%mesh%nz))

    ! Add the specified data field to the file
    SELECT CASE (choice_output_field)
      CASE ('none')
        ! Do nothing

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      CASE ('Hi_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_init%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hi_init', d_grid_vec_partial_2D)
      CASE ('Hb_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_init%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hb_init', d_grid_vec_partial_2D)
      CASE ('Hs_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_init%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hs_init', d_grid_vec_partial_2D)
      CASE ('SL_init')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_init%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'SL_init', d_grid_vec_partial_2D)

      ! Present-day ice-sheet geometry
      CASE ('Hi_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_PD%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hi_PD', d_grid_vec_partial_2D)
      CASE ('Hb_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_PD%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hb_PD', d_grid_vec_partial_2D)
      CASE ('Hs_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_PD%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hs_PD', d_grid_vec_partial_2D)
      CASE ('SL_PD')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_PD%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'SL_PD', d_grid_vec_partial_2D)

      ! GIA equilibrium ice-sheet geometry
      CASE ('Hi_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_GIAeq%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hi_GIAeq', d_grid_vec_partial_2D)
      CASE ('Hb_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_GIAeq%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hb_GIAeq', d_grid_vec_partial_2D)
      CASE ('Hs_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_GIAeq%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hs_GIAeq', d_grid_vec_partial_2D)
      CASE ('SL_GIAeq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%refgeo_GIAeq%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'SL_GIAeq', d_grid_vec_partial_2D)

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      CASE ('Hi')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%Hi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hi', d_grid_vec_partial_2D)
      CASE ('Hb')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%Hb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hb', d_grid_vec_partial_2D)
      CASE ('Hs')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%Hs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hs', d_grid_vec_partial_2D)
      CASE ('Hib')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%Hib, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'Hib', d_grid_vec_partial_2D)
      CASE ('SL')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%SL, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'SL', d_grid_vec_partial_2D)
      CASE ('TAF')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%TAF, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'TAF', d_grid_vec_partial_2D)

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      CASE ('dHi')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHi, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHi', d_grid_vec_partial_2D)
      CASE ('dHb')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHb, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHb', d_grid_vec_partial_2D)
      CASE ('dHs')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHs, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHs', d_grid_vec_partial_2D)
      CASE ('dHib')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHib, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHib', d_grid_vec_partial_2D)

    ! ===== Geometry rates of change =====
    ! ====================================

      CASE ('dHi_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHi_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHi_dt', d_grid_vec_partial_2D)
      CASE ('dHb_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHb_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHb_dt', d_grid_vec_partial_2D)
      CASE ('dHs_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHs_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHs_dt', d_grid_vec_partial_2D)
      CASE ('dHib_dt')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%dHib_dt, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'dHib_dt', d_grid_vec_partial_2D)

    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('mask_icefree_land')
      CASE ('mask_icefree_ocean')
      CASE ('mask_grounded_ice')
      CASE ('mask_floating_ice')
      CASE ('mask_gl_gr')
      CASE ('mask_gl_fl')
      CASE ('mask_cf_gr')
      CASE ('mask_cf_fl')
      CASE ('mask')
      CASE ('basin_ID')

    ! ===== Area fractions =====
    ! ==========================

      ! NOTE: sub-grid area fractions cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('fraction_gr')
      CASE ('fraction_gr_b')
      CASE ('fraction_cf')

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%Ti, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'Ti', d_grid_vec_partial_3D)
      CASE ('Ti_pmp')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%Ti_pmp, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'Ti_pmp', d_grid_vec_partial_3D)
      CASE ('Cpi')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%Cpi, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'Cpi', d_grid_vec_partial_3D)
      CASE ('Ki')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%Ki, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'Ki', d_grid_vec_partial_3D)
      CASE ('internal_heating')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%internal_heating, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'internal_heating', d_grid_vec_partial_3D)
      CASE ('frictional_heating')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%frictional_heating, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'frictional_heating', d_grid_vec_partial_2D)
      CASE ('A_flow')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%A_flow, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'A_flow', d_grid_vec_partial_3D)

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      CASE ('u_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%u_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'u_3D', d_grid_vec_partial_3D)
      CASE ('v_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%v_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'v_3D', d_grid_vec_partial_3D)
      CASE ('u_3D_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_3D_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%w_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'w_3D', d_grid_vec_partial_3D)

      ! Vertically integrated
      CASE ('u_vav')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%u_vav, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'u_vav', d_grid_vec_partial_2D)
      CASE ('v_vav')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%v_vav, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'v_vav', d_grid_vec_partial_2D)
      CASE ('u_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('uabs_vav')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%uabs_vav, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'uabs_vav', d_grid_vec_partial_2D)
      CASE ('uabs_vav_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

      ! Surface
      CASE ('u_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%u_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'u_surf', d_grid_vec_partial_2D)
      CASE ('v_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%v_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'v_surf', d_grid_vec_partial_2D)
      CASE ('u_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%w_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'w_surf', d_grid_vec_partial_2D)
      CASE ('uabs_surf')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%uabs_surf, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'uabs_surf', d_grid_vec_partial_2D)
      CASE ('uabs_surf_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

      ! Base
      CASE ('u_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%u_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'u_base', d_grid_vec_partial_2D)
      CASE ('v_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%v_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'v_base', d_grid_vec_partial_2D)
      CASE ('u_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('v_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!
      CASE ('w_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%w_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'w_base', d_grid_vec_partial_2D)
      CASE ('uabs_base')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%uabs_base, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'uabs_base', d_grid_vec_partial_2D)
      CASE ('uabs_base_b')
        ! NOTE: mapping from mesh triangles to square grid is not (yet) available!

    ! === Strain rates ===
    ! ====================

      CASE ('du_dx_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%du_dx_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'du_dx_3D', d_grid_vec_partial_3D)
      CASE ('du_dy_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%du_dy_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'du_dy_3D', d_grid_vec_partial_3D)
      CASE ('du_dz_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%du_dz_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'du_dz_3D', d_grid_vec_partial_3D)
      CASE ('dv_dx_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%dv_dx_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'dv_dx_3D', d_grid_vec_partial_3D)
      CASE ('dv_dy_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%dv_dy_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'dv_dy_3D', d_grid_vec_partial_3D)
      CASE ('dv_dz_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%dv_dz_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'dv_dz_3D', d_grid_vec_partial_3D)
      CASE ('dw_dx_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%dw_dx_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'dw_dx_3D', d_grid_vec_partial_3D)
      CASE ('dw_dy_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%dw_dy_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'dw_dy_3D', d_grid_vec_partial_3D)
      CASE ('dw_dz_3D')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%ice%dw_dz_3D, d_grid_vec_partial_3D)
        CALL write_to_field_multopt_grid_dp_3D( region%output_grid, filename, ncid, 'dw_dz_3D', d_grid_vec_partial_3D)

    ! == Basal conditions ==
    ! ======================

      ! Basal hydrology
      CASE ('pore_water_pressure')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%pore_water_pressure, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'pore_water_pressure', d_grid_vec_partial_2D)
      CASE ('overburden_pressure')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%overburden_pressure, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'overburden_pressure', d_grid_vec_partial_2D)
      CASE ('effective_pressure')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%effective_pressure, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'effective_pressure', d_grid_vec_partial_2D)

      ! Basal roughness / friction
      CASE ('phi_fric')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%phi_fric, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'phi_fric', d_grid_vec_partial_2D)
      CASE ('tau_c')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%tau_c, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'tau_c', d_grid_vec_partial_2D)
      CASE ('alpha_sq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%alpha_sq, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'alpha_sq', d_grid_vec_partial_2D)
      CASE ('beta_sq')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%beta_sq, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'beta_sq', d_grid_vec_partial_2D)
      CASE ('beta_b')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%beta_b, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'beta_b', d_grid_vec_partial_2D)

       ! Basal sliding
      CASE ('friction_coef_1')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%friction_coef_1, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'friction_coef_1', d_grid_vec_partial_2D)
      CASE ('friction_coef_2')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%friction_coef_2, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'friction_coef_2', d_grid_vec_partial_2D)
      CASE ('basal_shear_stress')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%basal_shear_stress, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'basal_shear_stress', d_grid_vec_partial_2D)

      ! Geothermal heat
      CASE ('geothermal_heat_flux')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%ice%geothermal_heat_flux, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'geothermal_heat_flux', d_grid_vec_partial_2D)

    ! == Climate ==
    ! =============

      ! Main climate variables
      CASE ('T2m')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%climate%T2m, d_grid_vec_partial_2D_monthly)
        CALL write_to_field_multopt_grid_dp_2D_monthly( region%output_grid, filename, ncid, 'T2m', d_grid_vec_partial_2D_monthly)
      CASE ('Precip')
        CALL map_from_mesh_to_xy_grid_3D( region%mesh, region%output_grid, region%climate%Precip, d_grid_vec_partial_2D_monthly)
        CALL write_to_field_multopt_grid_dp_2D_monthly( region%output_grid, filename, ncid, 'Precip', d_grid_vec_partial_2D_monthly)

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
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%SMB%SMB, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'SMB', d_grid_vec_partial_2D)

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      CASE ('BMB')
        CALL map_from_mesh_to_xy_grid_2D( region%mesh, region%output_grid, region%BMB%BMB, d_grid_vec_partial_2D)
        CALL write_to_field_multopt_grid_dp_2D( region%output_grid, filename, ncid, 'BMB', d_grid_vec_partial_2D)

    ! ===== End of user-defined output fields =====
    ! =============================================

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_output_field "' // TRIM( choice_output_field) // '"!')
    END SELECT

    ! Clean up after yourself
    DEALLOCATE( d_grid_vec_partial_2D)
    DEALLOCATE( d_grid_vec_partial_2D_monthly)
    DEALLOCATE( d_grid_vec_partial_3D)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_main_regional_output_file_grid_field

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
    IF (par%master) WRITE(0,'(A)') '  Creating mesh output file "' // colour_string( TRIM( region%output_filename_mesh), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( region%output_filename_mesh, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( region%output_filename_mesh, ncid, region%mesh)

    ! Add time, zeta, and month dimensions+variables to the file
    CALL add_time_dimension_to_file(  region%output_filename_mesh, ncid)
    CALL add_zeta_dimension_to_file(  region%output_filename_mesh, ncid, region%mesh%zeta)
    CALL add_month_dimension_to_file( region%output_filename_mesh, ncid)

    ! Add all data fields to the file
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
    IF (par%master) WRITE(0,'(A)') '  Creating grid output file "' // colour_string( TRIM( region%output_filename_grid), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( region%output_filename_grid, ncid)

    ! Set up the grid in the file
    CALL setup_xy_grid_in_netcdf_file( region%output_filename_grid, ncid, region%output_grid)

    ! Add time, zeta, and month dimensions+variables to the file
    CALL add_time_dimension_to_file(  region%output_filename_grid, ncid)
    CALL add_zeta_dimension_to_file(  region%output_filename_grid, ncid, region%mesh%zeta)
    CALL add_month_dimension_to_file( region%output_filename_grid, ncid)

    ! Add all data fields to the file
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

  SUBROUTINE create_main_regional_output_file_mesh_field( filename, ncid, choice_output_field)
    ! Create the main regional output NetCDF file - mesh version
    !
    ! Add a single field to the file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_output_field

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

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      CASE ('dHi')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHi', long_name = 'Ice thickness difference (w.r.t. to reference)', units = 'm')
      CASE ('dHb')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHb', long_name = 'Bedrock elevation difference (w.r.t. to reference)', units = 'm')
      CASE ('dHs')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHs', long_name = 'Surface elevation difference (w.r.t. to reference)', units = 'm')
      CASE ('dHib')
        CALL add_field_mesh_dp_2D( filename, ncid, 'dHib', long_name = 'Ice base elevation difference (w.r.t. to reference)', units = 'm')

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
      CASE ('mask_gl_gr')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_gl_gr', long_name = 'Mask indicating grounded side of grounding line')
      CASE ('mask_gl_fl')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_gl_fl', long_name = 'Mask indicating floating side of grounding line')
      CASE ('mask_cf_gr')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_cf_gr', long_name = 'Mask indicating grounded calving front')
      CASE ('mask_cf_fl')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask_cf_fl', long_name = 'Mask indicating floating calving front')
      CASE ('mask')
        CALL add_field_mesh_int_2D( filename, ncid, 'mask', long_name = 'General mask')
      CASE ('basin_ID')
        CALL add_field_mesh_int_2D( filename, ncid, 'basin_ID', long_name = 'Drainage basin ID')

    ! ===== Area fractions =====
    ! ==========================

      CASE ('fraction_gr')
        CALL add_field_mesh_dp_2D( filename, ncid, 'fraction_gr', long_name = 'Grounded area fractions of vertices')
      CASE ('fraction_gr_b')
        CALL add_field_mesh_dp_2D_b( filename, ncid, 'fraction_gr_b', long_name = 'Grounded area fractions of triangles')
      CASE ('fraction_cf')
        CALL add_field_mesh_dp_2D( filename, ncid, 'fraction_gr', long_name = 'Ice-covered area fractions of calving fronts')

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL add_field_mesh_dp_3D( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')
      CASE ('Ti_pmp')
        CALL add_field_mesh_dp_3D( filename, ncid, 'Ti_pmp', long_name = 'Pressure melting point temperature', units = 'K')
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

    ! == Basal conditions ==
    ! ======================

      ! Basal hydrology
      CASE ('pore_water_pressure')
        CALL add_field_mesh_dp_2D( filename, ncid, 'pore_water_pressure', long_name = 'Basal pore water pressure', units = 'Pa')
      CASE ('overburden_pressure')
        CALL add_field_mesh_dp_2D( filename, ncid, 'overburden_pressure', long_name = 'Basal overburden pressure', units = 'Pa')
      CASE ('effective_pressure')
        CALL add_field_mesh_dp_2D( filename, ncid, 'effective_pressure', long_name = 'Basal effective pressure', units = 'Pa')

      ! Basal roughness / friction
      CASE ('phi_fric')
        CALL add_field_mesh_dp_2D( filename, ncid, 'phi_fric', long_name = 'Till friction angle', units = 'degrees')
      CASE ('tau_c')
        CALL add_field_mesh_dp_2D( filename, ncid, 'tau_c', long_name = 'Till yield stress', units = 'Pa')
      CASE ('alpha_sq')
        CALL add_field_mesh_dp_2D( filename, ncid, 'alpha_sq', long_name = 'Coulomb-law friction coefficient')
      CASE ('beta_sq')
        CALL add_field_mesh_dp_2D( filename, ncid, 'beta_sq', long_name = 'Power-law friction coefficient', units = 'Pa m^−1/3 yr^1/3')
      CASE ('beta_b')
        CALL add_field_mesh_dp_2D( filename, ncid, 'beta_b', long_name = 'Basal friction coefficient', units = 'Pa m^-1 yr')

      ! Basal sliding
      CASE ('friction_coef_1')
        CALL add_field_mesh_dp_2D( filename, ncid, 'friction_coef_1', long_name = 'Generic basal friction coefficient 1')
      CASE ('friction_coef_2')
        CALL add_field_mesh_dp_2D( filename, ncid, 'friction_coef_2', long_name = 'Generic basal friction coefficient 2')
      CASE ('basal_shear_stress')
        CALL add_field_mesh_dp_2D( filename, ncid, 'basal_shear_stress', long_name = 'Basal shear stress', units = 'Pa')

      ! Geothermal heat
      CASE ('geothermal_heat_flux')
        CALL add_field_mesh_dp_2D( filename, ncid, 'geothermal_heat_flux', long_name = 'Geothermal heat flux', units = 'J yr^-1 m^-2')

    ! == Climate ==
    ! =============

      ! Main climate variables
      CASE ('T2m')
        CALL add_field_mesh_dp_2D_monthly( filename, ncid, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
      CASE ('Precip')
        CALL add_field_mesh_dp_2D_monthly( filename, ncid, 'Precip', long_name = 'Monthly total precipitation', units = 'm.w.e.')

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
        CALL add_field_mesh_dp_2D( filename, ncid, 'SMB', long_name = 'Surface mass balance', units = 'm yr^-1')

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      CASE ('BMB')
        CALL add_field_mesh_dp_2D( filename, ncid, 'BMB', long_name = 'Basal mass balance', units = 'm yr^-1')

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
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: filename
    INTEGER                                            , INTENT(IN)    :: ncid
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: choice_output_field

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

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      CASE ('dHi')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHi', long_name = 'Ice thickness difference (w.r.t. to reference)', units = 'm')
      CASE ('dHb')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHb', long_name = 'Bedrock elevation difference (w.r.t. to reference)', units = 'm')
      CASE ('dHs')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHs', long_name = 'Surface elevation difference (w.r.t. to reference)', units = 'm')
      CASE ('dHib')
        CALL add_field_grid_dp_2D( filename, ncid, 'dHib', long_name = 'Ice base elevation difference (w.r.t. to reference)', units = 'm')

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

    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('mask_icefree_land')
      CASE ('mask_icefree_ocean')
      CASE ('mask_grounded_ice')
      CASE ('mask_floating_ice')
      CASE ('mask_gl_gr')
      CASE ('mask_gl_fl')
      CASE ('mask_cf_gr')
      CASE ('mask_cf_fl')
      CASE ('mask')
      CASE ('basin_ID')

    ! ===== Area fractions =====
    ! ==========================

      ! NOTE: sub-grid area fractions cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      CASE ('fraction_gr')
      CASE ('fraction_gr_b')
      CASE ('fraction_cf')

    ! === Thermodynamics and rheology ===
    ! ===================================

      CASE ('Ti')
        CALL add_field_grid_dp_3D( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')
      CASE ('Ti_pmp')
        CALL add_field_grid_dp_3D( filename, ncid, 'Ti_pmp', long_name = 'Pressure melting point temperature', units = 'K')
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

    ! == Basal conditions ==
    ! ======================

      ! Basal hydrology
      CASE ('pore_water_pressure')
        CALL add_field_grid_dp_2D( filename, ncid, 'pore_water_pressure', long_name = 'Basal pore water pressure', units = 'Pa')
      CASE ('overburden_pressure')
        CALL add_field_grid_dp_2D( filename, ncid, 'overburden_pressure', long_name = 'Basal overburden pressure', units = 'Pa')
      CASE ('effective_pressure')
        CALL add_field_grid_dp_2D( filename, ncid, 'effective_pressure', long_name = 'Basal effective pressure', units = 'Pa')

      ! Basal roughness / friction
      CASE ('phi_fric')
        CALL add_field_grid_dp_2D( filename, ncid, 'phi_fric', long_name = 'Till friction angle', units = 'degrees')
      CASE ('tau_c')
        CALL add_field_grid_dp_2D( filename, ncid, 'tau_c', long_name = 'Till yield stress', units = 'Pa')
      CASE ('alpha_sq')
        CALL add_field_grid_dp_2D( filename, ncid, 'alpha_sq', long_name = 'Coulomb-law friction coefficient')
      CASE ('beta_sq')
        CALL add_field_grid_dp_2D( filename, ncid, 'beta_sq', long_name = 'Power-law friction coefficient', units = 'Pa m^−1/3 yr^1/3')
      CASE ('beta_b')
        CALL add_field_grid_dp_2D( filename, ncid, 'beta_b', long_name = 'Basal friction coefficient', units = 'Pa m^-1 yr')

      ! Basal sliding
      CASE ('friction_coef_1')
        CALL add_field_grid_dp_2D( filename, ncid, 'friction_coef_1', long_name = 'Generic basal friction coefficient 1')
      CASE ('friction_coef_2')
        CALL add_field_grid_dp_2D( filename, ncid, 'friction_coef_2', long_name = 'Generic basal friction coefficient 2')
      CASE ('basal_shear_stress')
        CALL add_field_grid_dp_2D( filename, ncid, 'basal_shear_stress', long_name = 'Basal shear stress', units = 'Pa')

      ! Geothermal heat
      CASE ('geothermal_heat_flux')
        CALL add_field_grid_dp_2D( filename, ncid, 'geothermal_heat_flux', long_name = 'Geothermal heat flux', units = 'J yr^-1 m^-2')

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

    ! ===== End of user-defined output fields =====
    ! =============================================

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_output_field "' // TRIM( choice_output_field) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_main_regional_output_file_grid_field

END MODULE main_regional_output