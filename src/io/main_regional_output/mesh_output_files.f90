module mesh_output_files

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use grid_basic, only: type_grid
  use region_types, only: type_model_region
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use netcdf_io_main
  use netcdf_bedrock_CDF
  use netcdf, only: NF90_DOUBLE
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan
  use mesh_contour, only: calc_mesh_contour

  implicit none

  private

  public :: create_main_regional_output_file_mesh, write_to_main_regional_output_file_mesh

contains

  subroutine write_to_main_regional_output_file_mesh( region)
    !< Write to the main regional output NetCDF file - mesh version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_main_regional_output_file_mesh'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to mesh output file "' // colour_string( trim( region%output_filename_mesh), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( region%output_filename_mesh, ncid)

    ! write the time to the file
    call write_time_to_file( region%output_filename_mesh, ncid, region%time)

    ! write the default data fields to the file
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'Hi')
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'Hb')
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'Hs')
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'SL')
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'u_surf')
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'v_surf')
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, 'uabs_surf')

    ! write all user-defined data fields to the file
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_01)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_02)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_03)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_04)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_05)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_06)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_07)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_08)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_09)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_10)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_11)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_12)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_13)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_14)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_15)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_16)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_17)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_18)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_19)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_20)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_21)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_22)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_23)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_24)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_25)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_26)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_27)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_28)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_29)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_30)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_31)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_32)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_33)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_34)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_35)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_36)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_37)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_38)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_39)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_40)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_41)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_42)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_43)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_44)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_45)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_46)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_47)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_48)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_49)
    call write_to_main_regional_output_file_mesh_field( region, region%output_filename_mesh, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_main_regional_output_file_mesh

  subroutine write_to_main_regional_output_file_mesh_field( region, filename, ncid, choice_output_field)
    !< Write a specific field to the main regional output NetCDF file - mesh version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid
    character(len=*),        intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_main_regional_output_file_mesh_field'
    integer, dimension(region%mesh%vi1:region%mesh%vi2) :: mask_int

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Add the specified data field to the file
    select case (choice_output_field)
      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')
      case ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      case ('resolution')
        ! Do nothing - this is already part of the regular mesh data; only write this to the square grid output

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      case ('Hi_init')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hi_init', region%refgeo_init%Hi)
      case ('Hb_init')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hb_init', region%refgeo_init%Hb)
      case ('Hs_init')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hs_init', region%refgeo_init%Hs)
      case ('SL_init')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'SL_init', region%refgeo_init%SL)

      ! Present-day ice-sheet geometry
      case ('Hi_PD')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hi_PD', region%refgeo_PD%Hi)
      case ('Hb_PD')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hb_PD', region%refgeo_PD%Hb)
      case ('Hs_PD')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hs_PD', region%refgeo_PD%Hs)
      case ('SL_PD')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'SL_PD', region%refgeo_PD%SL)

      ! GIA equilibrium ice-sheet geometry
      case ('Hi_GIAeq')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hi_GIAeq', region%refgeo_GIAeq%Hi)
      case ('Hb_GIAeq')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hb_GIAeq', region%refgeo_GIAeq%Hb)
      case ('Hs_GIAeq')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'Hs_GIAeq', region%refgeo_GIAeq%Hs)
      case ('SL_GIAeq')
        call write_to_field_multopt_mesh_dp_2D_notime( region%mesh, filename, ncid, 'SL_GIAeq', region%refgeo_GIAeq%SL)

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      case ('Hi')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hi', region%ice%Hi)
      case ('Hb')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hb', region%ice%Hb)
      case ('Hs')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hs', region%ice%Hs)
      case ('Hib')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hib', region%ice%Hib)
      case ('SL')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'SL', region%ice%SL)
      case ('TAF')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'TAF', region%ice%TAF)
      case ('Hi_eff')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hi_eff', region%ice%Hi_eff)
      case ('Hs_slope')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Hs_slope', region%ice%Hs_slope)
      case ('grounding_line')
        call write_grounding_line_to_file( filename, ncid, region%mesh, region%ice)
      case ('ice_margin')
        call write_ice_margin_to_file( filename, ncid, region%mesh, region%ice)
      case ('calving_front')
        call write_calving_front_to_file( filename, ncid, region%mesh, region%ice)
      case ('coastline')
        call write_coastline_to_file( filename, ncid, region%mesh, region%ice)
      case ('grounded_ice_contour')
        call write_grounded_ice_contour_to_file( filename, ncid, region%mesh, region%ice)

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      case ('dHi')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi', region%ice%dHi)
      case ('dHb')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHb', region%ice%dHb)
      case ('dHs')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHs', region%ice%dHs)
      case ('dHib')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHib', region%ice%dHib)

    ! ===== Geometry rates of changes =====
    ! =====================================

      case ('dHi_dt')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt', region%ice%dHi_dt)
      case ('dHb_dt')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHb_dt', region%ice%dHb_dt)
      case ('dHs_dt')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHs_dt', region%ice%dHs_dt)
      case ('dHib_dt')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHib_dt', region%ice%dHib_dt)
      case ('dHi_dt_raw')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt_raw', region%ice%dHi_dt_raw)
      case ('dHi_dt_residual')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt_residual', region%ice%dHi_dt_residual)

    ! ===== Target quantities =====
    ! =============================

      case ('dHi_dt_target')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHi_dt_target', region%ice%dHi_dt_target)

    ! ===== Masks =====
    ! =================

      case ('mask_icefree_land')
        where (region%ice%mask_icefree_land)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_icefree_land', mask_int)
      case ('mask_icefree_ocean')
        where (region%ice%mask_icefree_ocean)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_icefree_ocean', mask_int)
      case ('mask_grounded_ice')
        where (region%ice%mask_grounded_ice)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_grounded_ice', mask_int)
      case ('mask_floating_ice')
        where (region%ice%mask_floating_ice)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_floating_ice', mask_int)
      case ('mask_margin')
        where (region%ice%mask_margin)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_margin', mask_int)
      case ('mask_gl_gr')
        where (region%ice%mask_gl_gr)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_gl_gr', mask_int)
      case ('mask_gl_fl')
        where (region%ice%mask_gl_fl)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_gl_fl', mask_int)
      case ('mask_cf_gr')
        where (region%ice%mask_cf_gr)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_cf_gr', mask_int)
      case ('mask_cf_fl')
        where (region%ice%mask_cf_fl)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_cf_fl', mask_int)
      case ('mask_coastline')
        where (region%ice%mask_coastline)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_coastline', mask_int)
      case ('mask_ROI')
        where (region%ice%mask_ROI)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_ROI', mask_int)
      case ('mask_SGD')
        where (region%ice%mask_SGD)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask_SGD', mask_int)
      case ('mask')
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'mask', region%ice%mask)
      case ('basin_ID')
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'basin_ID', region%ice%basin_ID)

    ! ===== Area fractions =====
    ! ==========================

      case ('fraction_gr')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'fraction_gr', region%ice%fraction_gr)
      case ('fraction_gr_b')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'fraction_gr_b', region%ice%fraction_gr_b)
      case ('fraction_margin')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'fraction_margin', region%ice%fraction_margin)

    ! === Thermodynamics and rheology ===
    ! ===================================

      case ('Ti')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ti', region%ice%Ti)
      case ('Ti_pmp')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ti_pmp', region%ice%Ti_pmp)
      case ('Ti_hom')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'Ti_hom', region%ice%Ti_hom)
      case ('Cpi')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Cpi', region%ice%Cpi)
      case ('Ki')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'Ki', region%ice%Ki)
      case ('internal_heating')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'internal_heating', region%ice%internal_heating)
      case ('frictional_heating')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'frictional_heating', region%ice%frictional_heating)
      case ('A_flow')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'A_flow', region%ice%A_flow)

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      case ('u_3D')
        call write_to_field_multopt_mesh_dp_3D_b( region%mesh, filename, ncid, 'u_3D', region%ice%u_3D_b)
      case ('v_3D')
        call write_to_field_multopt_mesh_dp_3D_b( region%mesh, filename, ncid, 'v_3D', region%ice%v_3D_b)
      case ('u_3D_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_3D_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('w_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'w_3D', region%ice%w_3D)

      ! Vertically integrated
      case ('u_vav')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'u_vav', region%ice%u_vav_b)
      case ('v_vav')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'v_vav', region%ice%v_vav_b)
      case ('u_vav_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_vav_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('uabs_vav')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'uabs_vav', region%ice%uabs_vav_b)
      case ('uabs_vav_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')

      ! Surface
      case ('u_surf')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'u_surf', region%ice%u_surf_b)
      case ('v_surf')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'v_surf', region%ice%v_surf_b)
      case ('u_surf_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_surf_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('w_surf')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'w_surf', region%ice%w_surf)
      case ('uabs_surf')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'uabs_surf', region%ice%uabs_surf_b)
      case ('uabs_surf_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')

      ! Base
      case ('u_base')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'u_base', region%ice%u_base_b)
      case ('v_base')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'v_base', region%ice%v_base_b)
      case ('u_base_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_base_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('w_base')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'w_base', region%ice%w_base)
      case ('uabs_base')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'uabs_base', region%ice%uabs_base_b)
      case ('uabs_base_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')

    ! === Strain rates ===
    ! ====================

      case ('du_dx_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'du_dx_3D', region%ice%du_dx_3D)
      case ('du_dy_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'du_dy_3D', region%ice%du_dy_3D)
      case ('du_dz_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'du_dz_3D', region%ice%du_dz_3D)
      case ('dv_dx_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dv_dx_3D', region%ice%dv_dx_3D)
      case ('dv_dy_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dv_dy_3D', region%ice%dv_dy_3D)
      case ('dv_dz_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dv_dz_3D', region%ice%dv_dz_3D)
      case ('dw_dx_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dw_dx_3D', region%ice%dw_dx_3D)
      case ('dw_dy_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dw_dy_3D', region%ice%dw_dy_3D)
      case ('dw_dz_3D')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'dw_dz_3D', region%ice%dw_dz_3D)

    ! == Ice flow regime ==
    ! =====================

      case ('divQ')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'divQ', region%ice%divQ)
      case ('R_shear')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'R_shear', region%ice%R_shear)

    ! == Ice P/C time stepping ==
    ! ===========================

      case ('pc_truncation_error')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pc_truncation_error', region%ice%pc%tau_np1)
      case ('pc_untolerated_events')
        call write_to_field_multopt_mesh_int_2D( region%mesh, filename, ncid, 'pc_untolerated_events', region%ice%pc%tau_n_guilty)

    ! == Basal hydrology ==
    ! =====================

      case ('pore_water_pressure')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pore_water_pressure', region%ice%pore_water_pressure)
      case ('overburden_pressure')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'overburden_pressure', region%ice%overburden_pressure)
      case ('effective_pressure')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'effective_pressure', region%ice%effective_pressure)
      case ('pore_water_likelihood')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pore_water_likelihood', region%ice%pore_water_likelihood)
      case ('pore_water_fraction')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'pore_water_fraction', region%ice%pore_water_fraction)

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      case ('till_friction_angle')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'till_friction_angle', region%bed_roughness%till_friction_angle)
      case ('alpha_sq')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'alpha_sq', region%bed_roughness%alpha_sq)
      case ('beta_sq')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'beta_sq', region%bed_roughness%beta_sq)

      ! Basal friction and shear stress
      case ('till_yield_stress')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'till_yield_stress', region%ice%till_yield_stress)
      case ('basal_friction_coefficient')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'basal_friction_coefficient', region%ice%basal_friction_coefficient)
      case ('basal_shear_stress')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'basal_shear_stress', region%ice%basal_shear_stress)

      ! Bed roughness nudging - H, dH/dt, flowline
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up', &
          region%bed_roughness%nudging_H_dHdt_flowline%deltaHs_av_up)
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down', &
          region%bed_roughness%nudging_H_dHdt_flowline%deltaHs_av_down)
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up', &
          region%bed_roughness%nudging_H_dHdt_flowline%dHs_dt_av_up)
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down', &
          region%bed_roughness%nudging_H_dHdt_flowline%dHs_dt_av_down)
      case ('bed_roughness_nudge_H_dHdt_flowline_R')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_R', &
          region%bed_roughness%nudging_H_dHdt_flowline%R)
      case ('bed_roughness_nudge_H_dHdt_flowline_I_tot')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_I_tot', &
          region%bed_roughness%nudging_H_dHdt_flowline%I_tot)
      case ('bed_roughness_nudge_H_dHdt_flowline_dC_dt')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dC_dt', &
          region%bed_roughness%nudging_H_dHdt_flowline%dC_dt)

      ! Bed roughness nudging - H, u, flowline
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_up')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltaHs_av_up', &
          region%bed_roughness%nudging_H_u_flowline%deltaHs_av_up)
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_down')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltaHs_av_down', &
          region%bed_roughness%nudging_H_u_flowline%deltaHs_av_down)
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_up')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltau_av_up', &
          region%bed_roughness%nudging_H_u_flowline%deltau_av_up)
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_down')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltau_av_down', &
          region%bed_roughness%nudging_H_u_flowline%deltau_av_down)
      case ('bed_roughness_nudge_H_u_flowline_R')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_R', &
          region%bed_roughness%nudging_H_u_flowline%R)
      case ('bed_roughness_nudge_H_u_flowline_I_tot')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_I_tot', &
          region%bed_roughness%nudging_H_u_flowline%I_tot)
      case ('bed_roughness_nudge_H_u_flowline_dC_dt')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_dC_dt', &
          region%bed_roughness%nudging_H_u_flowline%dC_dt)
      case ('bed_roughness_nudge_H_u_target_velocity')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, &
          'bed_roughness_nudge_H_u_target_velocity', &
          region%bed_roughness%nudging_H_u_flowline%uabs_surf_target_b)

    ! == Geothermal heat ==
    ! =====================

      case ('geothermal_heat_flux')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'geothermal_heat_flux', region%ice%geothermal_heat_flux)

    ! == Climate ==
    ! =============

      ! Main climate variables
      case ('T2m')
        call write_to_field_multopt_mesh_dp_2D_monthly( region%mesh, filename, ncid, 'T2m', region%climate%T2m)
      case ('Precip')
        call write_to_field_multopt_mesh_dp_2D_monthly( region%mesh, filename, ncid, 'Precip', region%climate%Precip)
      case ('Q_TOA')
        call write_to_field_multopt_mesh_dp_2D_monthly( region%mesh, filename, ncid, 'Q_TOA', region%climate%snapshot%Q_TOA)

    ! == Ocean ==
    ! ===========

      ! Main ocean variables
      case ('T_ocean')
        call write_to_field_multopt_mesh_dp_3D_ocean( region%mesh, filename, ncid, 'T_ocean', region%ocean%T)
      case ('S_ocean')
        call write_to_field_multopt_mesh_dp_3D_ocean( region%mesh, filename, ncid, 'S_ocean', region%ocean%S)
      case ('T_draft')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'T_draft', region%ocean%T_draft)
      case ('T_freezing_point')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'T_freezing_point', region%ocean%T_freezing_point)

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      case ('SMB')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'SMB', region%SMB%SMB)
      case ('Albedo')
        call write_to_field_multopt_mesh_dp_2D_monthly( region%mesh, filename, ncid, 'Albedo', region%SMB%IMAUITM%Albedo)
      CASE ('FirnDepth')
        call write_to_field_multopt_mesh_dp_2D_monthly( region%mesh, filename, ncid, 'FirnDepth', region%SMB%IMAUITM%FirnDepth)
      CASE ('MeltPreviousYear')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'MeltPreviousYear', region%SMB%IMAUITM%MeltPreviousYear)

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      case ('BMB')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'BMB', region%BMB%BMB)
      case ('BMB_inv')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'BMB_inv', region%BMB%BMB_inv)
      case ('BMB_transition_phase')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'BMB_transition_phase', region%BMB%BMB_transition_phase)
      case ('BMB_modelled')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'BMB_modelled', region%BMB%BMB_modelled)

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'H_lad', region%BMB%laddie%now%H, d_is_hybrid = .true.)
      case ('U_lad')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'U_lad', region%BMB%laddie%now%U, d_is_hybrid = .true.)
      case ('V_lad')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'V_lad', region%BMB%laddie%now%V, d_is_hybrid = .true.)
      case ('T_lad')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'T_lad', region%BMB%laddie%now%T, d_is_hybrid = .true.)
      case ('S_lad')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'S_lad', region%BMB%laddie%now%S, d_is_hybrid = .true.)

      ! Useful laddie fields
      case ('drho_amb')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'drho_amb', region%BMB%laddie%drho_amb, d_is_hybrid = .true.)
      case ('drho_base')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'drho_base', region%BMB%laddie%drho_base, d_is_hybrid = .true.)
      case ('entr')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'entr', region%BMB%laddie%entr, d_is_hybrid = .true.)
      case ('entr_dmin')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'entr_dmin', region%BMB%laddie%entr_dmin, d_is_hybrid = .true.)
      case ('SGD')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'SGD', region%BMB%laddie%SGD, d_is_hybrid = .true.)
      case ('melt')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'melt', region%BMB%laddie%melt, d_is_hybrid = .true.)
      case ('divQH')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'divQH', region%BMB%laddie%divQH, d_is_hybrid = .true.)
      case ('divQT')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'divQT', region%BMB%laddie%divQT, d_is_hybrid = .true.)
      case ('divQS')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'divQS', region%BMB%laddie%divQS, d_is_hybrid = .true.)
      case ('diffT')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'diffT', region%BMB%laddie%diffT, d_is_hybrid = .true.)
      case ('diffS')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'diffS', region%BMB%laddie%diffS, d_is_hybrid = .true.)
      case ('viscU')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'viscU', region%BMB%laddie%viscU, d_is_hybrid = .true.)
      case ('viscV')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'viscV', region%BMB%laddie%viscV, d_is_hybrid = .true.)
      case ('T_base')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'T_base', region%BMB%laddie%T_base, d_is_hybrid = .true.)
      case ('T_amb')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'T_amb', region%BMB%laddie%T_amb, d_is_hybrid = .true.)
      case ('u_star')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'u_star', region%BMB%laddie%u_star, d_is_hybrid = .true.)
      case ('gamma_T')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'gamma_T', region%BMB%laddie%gamma_T, d_is_hybrid = .true.)
      case ('divQU')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'divQU', region%BMB%laddie%divQU, d_is_hybrid = .true.)
      case ('divQV')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'divQV', region%BMB%laddie%divQV, d_is_hybrid = .true.)
      case ('HU_lad')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'HU_lad', region%BMB%laddie%now%H_b*region%BMB%laddie%now%U, d_is_hybrid = .true.)
      case ('HV_lad')
        call write_to_field_multopt_mesh_dp_2D_b( region%mesh, filename, ncid, 'HV_lad', region%BMB%laddie%now%H_b*region%BMB%laddie%now%V, d_is_hybrid = .true.)

    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      case ('LMB')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'LMB', region%LMB%LMB)

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      case ('AMB')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'AMB', region%AMB%AMB)

    ! == Glacial isostatic adjustment ==
    ! ==================================

      ! Main GIA variables
      case ('dHb_next')
        call write_to_field_multopt_mesh_dp_2D( region%mesh, filename, ncid, 'dHb_next', region%GIA%dHb_next)

    ! == Tracer tracking ==
    ! =====================

      case ('age')
        call write_to_field_multopt_mesh_dp_3D( region%mesh, filename, ncid, 'age', region%tracer_tracking%age)

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_main_regional_output_file_mesh_field

  subroutine create_main_regional_output_file_mesh( region)
    !< Create the main regional output NetCDF file - mesh version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_main_regional_output_file_mesh'
    character(len=1024)            :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = trim( C%output_dir) // 'main_output_' // region%name
    call generate_filename_XXXXXdotnc( filename_base, region%output_filename_mesh)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating mesh output file "' // colour_string( trim( region%output_filename_mesh), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( region%output_filename_mesh, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( region%output_filename_mesh, ncid, region%mesh)

    if (C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .or. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') then
      ! Set up bedrock CDF in the file
      call setup_bedrock_CDF_in_netcdf_file( region%output_filename_mesh, ncid, region%ice)
    end if

    ! Add time, zeta, and month dimensions+variables to the file
    call add_time_dimension_to_file(  region%output_filename_mesh, ncid)
    call add_zeta_dimension_to_file(  region%output_filename_mesh, ncid, region%mesh%zeta)
    call add_month_dimension_to_file( region%output_filename_mesh, ncid)
    call add_depth_dimension_to_file( region%output_filename_mesh, ncid, C%z_ocean)

    ! Operator matrices
    if (C%do_write_matrix_operators) then
      call write_matrix_operators_to_netcdf_file( region%output_filename_mesh, ncid, region%mesh)
    end if

    ! Add the default data fields to the file
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'Hi')
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'Hb')
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'Hs')
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'SL')
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'u_surf')
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'v_surf')
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, 'uabs_surf')

    ! Add all user-defined data fields to the file
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_01)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_02)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_03)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_04)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_05)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_06)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_07)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_08)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_09)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_10)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_11)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_12)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_13)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_14)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_15)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_16)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_17)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_18)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_19)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_20)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_21)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_22)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_23)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_24)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_25)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_26)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_27)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_28)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_29)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_30)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_31)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_32)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_33)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_34)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_35)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_36)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_37)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_38)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_39)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_40)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_41)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_42)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_43)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_44)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_45)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_46)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_47)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_48)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_49)
    call create_main_regional_output_file_mesh_field( region%output_filename_mesh, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_main_regional_output_file_mesh

  subroutine create_main_regional_output_file_mesh_field( filename, ncid, choice_output_field)
    !< Create a single field in the main regional output NetCDF file - mesh version

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_main_regional_output_file_mesh_field'
    integer                        :: int_dummy, id_dim_ei, id_dim_two, id_dim_time
    integer                        :: id_var_grounding_line
    integer                        :: id_var_calving_front
    integer                        :: id_var_ice_margin
    integer                        :: id_var_coastline
    integer                        :: id_var_grounded_ice_contour

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Add the specified data field to the file
    select case (choice_output_field)

      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')

      case ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      case ('resolution')
        ! Do nothing - this is already part of the regular mesh data; only write this to the square grid output

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      case ('Hi_init')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hi_init', long_name = 'Initial ice thickness', units = 'm')
      case ('Hb_init')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hb_init', long_name = 'Initial bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs_init')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hs_init', long_name = 'Initial surface elevation', units = 'm w.r.t. PD sea level')
      case ('SL_init')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'SL_init', long_name = 'Initial geoid elevation', units = 'm w.r.t. PD sea level')

      ! Present-day ice-sheet geometry
      case ('Hi_PD')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hi_PD', long_name = 'Present-day ice thickness', units = 'm')
      case ('Hb_PD')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hb_PD', long_name = 'Present-day bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs_PD')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hs_PD', long_name = 'Present-day surface elevation', units = 'm w.r.t. PD sea level')
      case ('SL_PD')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'SL_PD', long_name = 'Present-day geoid elevation', units = 'm w.r.t. PD sea level')

      ! GIA equilibrium ice-sheet geometry
      case ('Hi_GIAeq')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hi_GIAeq', long_name = 'GIA equilibrium ice thickness', units = 'm')
      case ('Hb_GIAeq')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hb_GIAeq', long_name = 'GIA equilibrium bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs_GIAeq')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hs_GIAeq', long_name = 'GIA equilibrium surface elevation', units = 'm w.r.t. PD sea level')
      case ('SL_GIAeq')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'SL_GIAeq', long_name = 'GIA equilibrium geoid elevation', units = 'm w.r.t. PD sea level')

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      case ('Hi')
        call add_field_mesh_dp_2D( filename, ncid, 'Hi', long_name = 'Ice thickness', units = 'm')
      case ('Hb')
        call add_field_mesh_dp_2D( filename, ncid, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs')
        call add_field_mesh_dp_2D( filename, ncid, 'Hs', long_name = 'Surface elevation', units = 'm w.r.t. PD sea level')
      case ('Hib')
        call add_field_mesh_dp_2D( filename, ncid, 'Hib', long_name = 'Ice base elevation', units = 'm w.r.t. PD sea level')
      case ('SL')
        call add_field_mesh_dp_2D( filename, ncid, 'SL', long_name = 'Geoid elevation', units = 'm w.r.t. PD sea level')
      case ('TAF')
        call add_field_mesh_dp_2D( filename, ncid, 'TAF', long_name = 'Thickness above floatation', units = 'm')
      case ('Hi_eff')
        call add_field_mesh_dp_2D( filename, ncid, 'Hi_eff', long_name = 'Effective ice thickness', units = 'm')
      case ('Hs_slope')
        call add_field_mesh_dp_2D( filename, ncid, 'Hs_slope', long_name = 'Absolute surface gradient', units = '-')
      case ('grounding_line')
        call inquire_dim( filename, ncid, 'ei', int_dummy, id_dim_ei)
        call inquire_dim( filename, ncid, 'two', int_dummy, id_dim_two)
        call inquire_dim( filename, ncid, 'time', int_dummy, id_dim_time)
        call create_variable( filename, ncid, 'grounding_line', NF90_DOUBLE, (/ id_dim_ei, id_dim_two, id_dim_time /), id_var_grounding_line)
        call add_attribute_char( filename, ncid, id_var_grounding_line, 'long_name', 'Grounding-line coordinates')
        call add_attribute_char( filename, ncid, id_var_grounding_line, 'units', 'm')
        call add_attribute_char( filename, ncid, id_var_grounding_line, 'format', 'Matlab contour format')
      case ('ice_margin')
        call inquire_dim( filename, ncid, 'ei', int_dummy, id_dim_ei)
        call inquire_dim( filename, ncid, 'two', int_dummy, id_dim_two)
        call inquire_dim( filename, ncid, 'time', int_dummy, id_dim_time)
        call create_variable( filename, ncid, 'ice_margin', NF90_DOUBLE, (/ id_dim_ei, id_dim_two, id_dim_time /), id_var_ice_margin)
        call add_attribute_char( filename, ncid, id_var_ice_margin, 'long_name', 'Ice margin coordinates')
        call add_attribute_char( filename, ncid, id_var_ice_margin, 'units', 'm')
        call add_attribute_char( filename, ncid, id_var_ice_margin, 'format', 'Matlab contour format')
      case ('calving_front')
        call inquire_dim( filename, ncid, 'ei', int_dummy, id_dim_ei)
        call inquire_dim( filename, ncid, 'two', int_dummy, id_dim_two)
        call inquire_dim( filename, ncid, 'time', int_dummy, id_dim_time)
        call create_variable( filename, ncid, 'calving_front', NF90_DOUBLE, (/ id_dim_ei, id_dim_two, id_dim_time /), id_var_calving_front)
        call add_attribute_char( filename, ncid, id_var_calving_front, 'long_name', 'Calving-front coordinates')
        call add_attribute_char( filename, ncid, id_var_calving_front, 'units', 'm')
        call add_attribute_char( filename, ncid, id_var_calving_front, 'format', 'Matlab contour format')
      case ('coastline')
        call inquire_dim( filename, ncid, 'ei', int_dummy, id_dim_ei)
        call inquire_dim( filename, ncid, 'two', int_dummy, id_dim_two)
        call inquire_dim( filename, ncid, 'time', int_dummy, id_dim_time)
        call create_variable( filename, ncid, 'coastline', NF90_DOUBLE, (/ id_dim_ei, id_dim_two, id_dim_time /), id_var_coastline)
        call add_attribute_char( filename, ncid, id_var_coastline, 'long_name', 'Coastline coordinates')
        call add_attribute_char( filename, ncid, id_var_coastline, 'units', 'm')
        call add_attribute_char( filename, ncid, id_var_coastline, 'format', 'Matlab contour format')
      case ('grounded_ice_contour')
        call inquire_dim( filename, ncid, 'ei', int_dummy, id_dim_ei)
        call inquire_dim( filename, ncid, 'two', int_dummy, id_dim_two)
        call inquire_dim( filename, ncid, 'time', int_dummy, id_dim_time)
        call create_variable( filename, ncid, 'grounded_ice_contour', NF90_DOUBLE, (/ id_dim_ei, id_dim_two, id_dim_time /), id_var_grounded_ice_contour)
        call add_attribute_char( filename, ncid, id_var_grounded_ice_contour, 'long_name', 'Grounded ice contour coordinates')
        call add_attribute_char( filename, ncid, id_var_grounded_ice_contour, 'units', 'm')
        call add_attribute_char( filename, ncid, id_var_grounded_ice_contour, 'format', 'Matlab contour format')

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      case ('dHi')
        call add_field_mesh_dp_2D( filename, ncid, 'dHi', long_name = 'Ice thickness difference w.r.t. reference', units = 'm')
      case ('dHb')
        call add_field_mesh_dp_2D( filename, ncid, 'dHb', long_name = 'Bedrock elevation difference w.r.t. reference', units = 'm')
      case ('dHs')
        call add_field_mesh_dp_2D( filename, ncid, 'dHs', long_name = 'Surface elevation difference w.r.t. reference', units = 'm')
      case ('dHib')
        call add_field_mesh_dp_2D( filename, ncid, 'dHib', long_name = 'Ice base elevation difference w.r.t. reference', units = 'm')

    ! ===== Geometry rates of change =====
    ! ====================================

      case ('dHi_dt')
        call add_field_mesh_dp_2D( filename, ncid, 'dHi_dt', long_name = 'Ice thickness rate of change', units = 'm yr^-1')
      case ('dHb_dt')
        call add_field_mesh_dp_2D( filename, ncid, 'dHb_dt', long_name = 'Bedrock elevation rate of change', units = 'm yr^-1')
      case ('dHs_dt')
        call add_field_mesh_dp_2D( filename, ncid, 'dHs_dt', long_name = 'Surface elevation rate of change', units = 'm yr^-1')
      case ('dHib_dt')
        call add_field_mesh_dp_2D( filename, ncid, 'dHib_dt', long_name = 'Ice base elevation rate of change', units = 'm yr^-1')
      case ('dHi_dt_raw')
        call add_field_mesh_dp_2D( filename, ncid, 'dHi_dt_raw', long_name = 'Ice thickness rate of change before any modifications', units = 'm yr^-1')
      case ('dHi_dt_residual')
        call add_field_mesh_dp_2D( filename, ncid, 'dHi_dt_residual', long_name = 'Residual ice thickness rate of change during model calibration', units = 'm yr^-1')

    ! ===== Target quantities =====
    ! =============================

      case ('dHi_dt_target')
        call add_field_mesh_dp_2D( filename, ncid, 'dHi_dt_target', long_name = 'Target ice thickness rate of change during model calibration', units = 'm yr^-1')
      case ('uabs_surf_target')
        call add_field_mesh_dp_2D( filename, ncid, 'uabs_surf_target', long_name = 'Target ice surface speed during model calibration', units = 'm yr^-1')

    ! ===== Masks =====
    ! =================

      ! notE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      case ('mask_icefree_land')
        call add_field_mesh_int_2D( filename, ncid, 'mask_icefree_land', long_name = 'Mask indicating ice-free land')
      case ('mask_icefree_ocean')
        call add_field_mesh_int_2D( filename, ncid, 'mask_icefree_ocean', long_name = 'Mask indicating ice-free ocean')
      case ('mask_grounded_ice')
        call add_field_mesh_int_2D( filename, ncid, 'mask_grounded_ice', long_name = 'Mask indicating grounded ice')
      case ('mask_floating_ice')
        call add_field_mesh_int_2D( filename, ncid, 'mask_floating_ice', long_name = 'Mask indicating floating ice')
      case ('mask_margin')
        call add_field_mesh_int_2D( filename, ncid, 'mask_margin', long_name = 'Mask indicating ice next to ice-free')
      case ('mask_gl_gr')
        call add_field_mesh_int_2D( filename, ncid, 'mask_gl_gr', long_name = 'Mask indicating grounded side of grounding line')
      case ('mask_gl_fl')
        call add_field_mesh_int_2D( filename, ncid, 'mask_gl_fl', long_name = 'Mask indicating floating side of grounding line')
      case ('mask_cf_gr')
        call add_field_mesh_int_2D( filename, ncid, 'mask_cf_gr', long_name = 'Mask indicating grounded calving front')
      case ('mask_cf_fl')
        call add_field_mesh_int_2D( filename, ncid, 'mask_cf_fl', long_name = 'Mask indicating floating calving front')
      case ('mask_coastline')
        call add_field_mesh_int_2D( filename, ncid, 'mask_coastline', long_name = 'Mask indicating ice-free land next to ice-free ocean')
      case ('mask_ROI')
        call add_field_mesh_int_2D( filename, ncid, 'mask_ROI', long_name = 'Mask indicating ROI')
      case ('mask_SGD')
        call add_field_mesh_int_2D( filename, ncid, 'mask_SGD', long_name = 'Mask indicating potential subglacial discharge cells')
      case ('mask')
        call add_field_mesh_int_2D( filename, ncid, 'mask', long_name = 'General mask')
      case ('basin_ID')
        call add_field_mesh_int_2D( filename, ncid, 'basin_ID', long_name = 'Drainage basin ID', units = 'ID code')

    ! ===== Area fractions =====
    ! ==========================

      case ('fraction_gr')
        call add_field_mesh_dp_2D( filename, ncid, 'fraction_gr', long_name = 'Grounded area fractions of vertices', units = '0-1')
      case ('fraction_gr_b')
        call add_field_mesh_dp_2D_b( filename, ncid, 'fraction_gr_b', long_name = 'Grounded area fractions of triangles', units = '0-1')
      case ('fraction_margin')
        call add_field_mesh_dp_2D( filename, ncid, 'fraction_margin', long_name = 'Ice-covered area fractions of ice margins', units = '0-1')

    ! === Thermodynamics and rheology ===
    ! ===================================

      case ('Ti')
        call add_field_mesh_dp_3D( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')
      case ('Ti_pmp')
        call add_field_mesh_dp_3D( filename, ncid, 'Ti_pmp', long_name = 'Pressure melting point temperature', units = 'K')
      case ('Ti_hom')
        call add_field_mesh_dp_2D( filename, ncid, 'Ti_hom', long_name = 'Temperature at base w.r.t. pressure melting point', units = 'K')
      case ('Cpi')
        call add_field_mesh_dp_3D( filename, ncid, 'Cpi', long_name = 'Specific heat capacity', units = 'J kg^-1 K^-1')
      case ('Ki')
        call add_field_mesh_dp_3D( filename, ncid, 'Ki', long_name = 'Thermal conductivity', units = 'J m^-1 K^-1 yr^-1')
      case ('internal_heating')
        call add_field_mesh_dp_3D( filename, ncid, 'internal_heating', long_name = 'Internal heating', units = '?')
      case ('frictional_heating')
        call add_field_mesh_dp_2D( filename, ncid, 'frictional_heating', long_name = 'Frictional heating', units = '?')
      case ('A_flow')
        call add_field_mesh_dp_3D( filename, ncid, 'A_flow', long_name = 'Glens flow law factor', units = 'Pa^-3 y^-1')

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      case ('u_3D')
        call add_field_mesh_dp_3D_b( filename, ncid, 'u_3D', long_name = '3-D ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_3D')
        call add_field_mesh_dp_3D_b( filename, ncid, 'v_3D', long_name = '3-D ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_3D_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_3D_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('w_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'w_3D', long_name = '3-D ice velocity in the z-direction', units = 'm yr^-1')

      ! Vertically integrated
      case ('u_vav')
        call add_field_mesh_dp_2D_b( filename, ncid, 'u_vav', long_name = 'Vertically averaged ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_vav')
        call add_field_mesh_dp_2D_b( filename, ncid, 'v_vav', long_name = 'Vertically averaged ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_vav_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_vav_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('uabs_vav')
        call add_field_mesh_dp_2D_b( filename, ncid, 'uabs_vav', long_name = 'Vertically averaged absolute ice velocity', units = 'm yr^-1')
      case ('uabs_vav_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')

      ! Surface
      case ('u_surf')
        call add_field_mesh_dp_2D_b( filename, ncid, 'u_surf', long_name = 'Surface ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_surf')
        call add_field_mesh_dp_2D_b( filename, ncid, 'v_surf', long_name = 'Surface ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_surf_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_surf_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('w_surf')
        call add_field_mesh_dp_2D( filename, ncid, 'w_surf', long_name = 'Surface ice velocity in the z-direction', units = 'm yr^-1')
      case ('uabs_surf')
        call add_field_mesh_dp_2D_b( filename, ncid, 'uabs_surf', long_name = 'Absolute surface ice velocity', units = 'm yr^-1')
      case ('uabs_surf_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')

      ! Base
      case ('u_base')
        call add_field_mesh_dp_2D_b( filename, ncid, 'u_base', long_name = 'Basal ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_base')
        call add_field_mesh_dp_2D_b( filename, ncid, 'v_base', long_name = 'Basal ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_base_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('v_base_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')
      case ('w_base')
        call add_field_mesh_dp_2D( filename, ncid, 'w_base', long_name = 'Basal ice velocity in the z-direction', units = 'm yr^-1')
      case ('uabs_base')
        call add_field_mesh_dp_2D_b( filename, ncid, 'uabs_base', long_name = 'Absolute basal ice velocity', units = 'm yr^-1')
      case ('uabs_base_b')
        call crash( trim(choice_output_field)//' no longer an option; horizontal velocities are always returned on the b-grid')

    ! === Strain rates ===
    ! ====================

      case ('du_dx_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'du_dx_3D', long_name = '3-D xx strain rate', units = 'yr^-1')
      case ('du_dy_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'du_dy_3D', long_name = '3-D xy strain rate', units = 'yr^-1')
      case ('du_dz_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'du_dz_3D', long_name = '3-D xz strain rate', units = 'yr^-1')
      case ('dv_dx_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'dv_dx_3D', long_name = '3-D yx strain rate', units = 'yr^-1')
      case ('dv_dy_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'dv_dy_3D', long_name = '3-D yy strain rate', units = 'yr^-1')
      case ('dv_dz_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'dv_dz_3D', long_name = '3-D yz strain rate', units = 'yr^-1')
      case ('dw_dx_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'dw_dx_3D', long_name = '3-D zx strain rate', units = 'yr^-1')
      case ('dw_dy_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'dw_dy_3D', long_name = '3-D zy strain rate', units = 'yr^-1')
      case ('dw_dz_3D')
        call add_field_mesh_dp_3D( filename, ncid, 'dw_dz_3D', long_name = '3-D zz strain rate', units = 'yr^-1')

    ! == Ice flow regime ==
    ! =====================

      case ('divQ')
        call add_field_mesh_dp_2D( filename, ncid, 'divQ', long_name = 'Horizontal ice flux divergence', units = 'm yr^-1')
      case ('R_shear')
        call add_field_mesh_dp_2D( filename, ncid, 'R_shear', long_name = 'Slide/shear ratio', units = '0-1')

    ! == Ice P/C time stepping ==
    ! ===========================

      case ('pc_truncation_error')
        call add_field_mesh_dp_2D( filename, ncid, 'pc_truncation_error', long_name = 'Ice P/C truncation error tau', units = 'm')
      case ('pc_untolerated_events')
        call add_field_mesh_int_2D( filename, ncid, 'pc_untolerated_events', long_name = 'Ice P/C number of events above error tolerance', units = '-')

    ! == Basal hydrology ==
    ! =====================

      case ('pore_water_pressure')
        call add_field_mesh_dp_2D( filename, ncid, 'pore_water_pressure', long_name = 'Till pore water pressure', units = 'Pa')
      case ('overburden_pressure')
        call add_field_mesh_dp_2D( filename, ncid, 'overburden_pressure', long_name = 'Ice overburden pressure', units = 'Pa')
      case ('effective_pressure')
        call add_field_mesh_dp_2D( filename, ncid, 'effective_pressure', long_name = 'Effective basal pressure', units = 'Pa')
      case ('pore_water_likelihood')
        call add_field_mesh_dp_2D( filename, ncid, 'pore_water_likelihood', long_name = 'Till pore water likelihood', units = '0-1')
      case ('pore_water_fraction')
        call add_field_mesh_dp_2D( filename, ncid, 'pore_water_fraction', long_name = 'Fraction of overburden pressure reduced by pore water ', units = '0-1')

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      case ('till_friction_angle')
        call add_field_mesh_dp_2D( filename, ncid, 'till_friction_angle', long_name = 'Till friction angle', units = 'degrees')
      case ('alpha_sq')
        call add_field_mesh_dp_2D( filename, ncid, 'alpha_sq', long_name = 'Coulomb-law friction coefficientn', units = 'dimensionless')
      case ('beta_sq')
        call add_field_mesh_dp_2D( filename, ncid, 'beta_sq', long_name = 'Power-law friction coefficient', units = 'Pa m^−1/m yr^1/m')

        ! Basal friction and shear stress
      case ('till_yield_stress')
        call add_field_mesh_dp_2D( filename, ncid, 'till_yield_stress', long_name = 'Till yield stress', units = 'Pa')
      case ('basal_friction_coefficient')
        call add_field_mesh_dp_2D( filename, ncid, 'basal_friction_coefficient', long_name = 'Basal friction coefficient', units = 'Pa yr m^-1')
      case ('basal_shear_stress')
        call add_field_mesh_dp_2D( filename, ncid, 'basal_shear_stress', long_name = 'Basal shear stress', units = 'Pa')

      ! Bed roughness nudging - H, dH/dt, flowline
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up', &
          long_name = 'Upstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down', &
          long_name = 'Downstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up', &
          long_name = 'Upstream flowline-averaged thinning rate', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down', &
          long_name = 'Downstream flowline-averaged thinning rate', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_dHdt_flowline_R')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_R', &
          long_name = 'Ice flux-based scaling factor')
      case ('bed_roughness_nudge_H_dHdt_flowline_I_tot')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_I_tot', &
          long_name = 'Weighted average of flowline-averaged terms')
      case ('bed_roughness_nudge_H_dHdt_flowline_dC_dt')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dC_dt', &
          long_name = 'Bed roughness rate of change')

      ! Bed roughness nudging - H, u, flowline
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_up')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltaHs_av_up', &
          long_name = 'Upstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_down')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltaHs_av_down', &
          long_name = 'Downstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_up')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltau_av_up', &
          long_name = 'Upstream flowline-averaged velocity error', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_down')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltau_av_down', &
          long_name = 'Downstream flowline-averaged velocity error', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_u_flowline_R')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_R', &
          long_name = 'Ice flux-based scaling factor')
      case ('bed_roughness_nudge_H_u_flowline_I_tot')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_I_tot', &
          long_name = 'Weighted average of flowline-averaged terms')
      case ('bed_roughness_nudge_H_u_flowline_dC_dt')
        call add_field_mesh_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_dC_dt', &
          long_name = 'Bed roughness rate of change')
      case ('bed_roughness_nudge_H_u_target_velocity')
        call add_field_mesh_dp_2D_b( filename, ncid, &
          'bed_roughness_nudge_H_u_target_velocity', &
          long_name = 'Target velocity', units = 'm yr^-1')

    ! == Geothermal heat ==
    ! =====================

      case ('geothermal_heat_flux')
        call add_field_mesh_dp_2D( filename, ncid, 'geothermal_heat_flux', long_name = 'Geothermal heat flux', units = 'J m^-2 yr^-1')

    ! == Climate ==
    ! =============

      ! Main climate variables
      case ('T2m')
        call add_field_mesh_dp_2D_monthly( filename, ncid, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
      case ('Precip')
        call add_field_mesh_dp_2D_monthly( filename, ncid, 'Precip', long_name = 'Monthly total precipitation', units = 'm.w.e.')
      case ('Q_TOA')
        CALL add_field_mesh_dp_2D_monthly( filename, ncid, 'Q_TOA', long_name = 'Monthly insolation at the top of the atmosphere', units = 'W m^-2')

    ! == Ocean ==
    ! ===========

      ! Main ocean variables
      case ('T_ocean')
        call add_field_mesh_dp_3D_ocean( filename, ncid, 'T_ocean', long_name = 'Ocean temperature', units = 'deg C')
      case ('S_ocean')
        call add_field_mesh_dp_3D_ocean( filename, ncid, 'S_ocean', long_name = 'Ocean salinity', units = 'PSU')
      case ('T_draft')
        call add_field_mesh_dp_2D( filename, ncid, 'T_draft', long_name = 'Ocean temperature at ice draft', units = 'deg C')
      case ('T_freezing_point')
        call add_field_mesh_dp_2D( filename, ncid, 'T_freezing_point', long_name = 'Ocean freezing temperature at ice draft', units = 'deg C')

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      case ('SMB')
        call add_field_mesh_dp_2D( filename, ncid, 'SMB', long_name = 'Surface mass balance', units = 'm yr^-1')
      CASE ('Albedo')
        CALL add_field_mesh_dp_2D_monthly( filename, ncid, 'Albedo', long_name = 'Surface albedo', units = '0-1')
      CASE ('FirnDepth')
        CALL add_field_mesh_dp_2D_monthly( filename, ncid, 'FirnDepth', long_name = 'Monthly firn layer depth', units = 'm')
      CASE ('MeltPreviousYear')
        CALL add_field_mesh_dp_2D( filename, ncid, 'MeltPreviousYear', long_name = 'Total ice melt from previous year', units = 'm')

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      case ('BMB')
        call add_field_mesh_dp_2D( filename, ncid, 'BMB', long_name = 'Basal mass balance', units = 'm yr^-1')
      case ('BMB_inv')
        call add_field_mesh_dp_2D( filename, ncid, 'BMB_inv', long_name = 'Basal mass balance - inverted', units = 'm yr^-1')
      case ('BMB_transition_phase')
        call add_field_mesh_dp_2D( filename, ncid, 'BMB_transition_phase', long_name = 'Basal mass balance - transition phase', units = 'm yr^-1')
      case ('BMB_modelled')
        call add_field_mesh_dp_2D( filename, ncid, 'BMB_modelled', long_name = 'Basal mass balance - modelled', units = 'm yr^-1')

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call add_field_mesh_dp_2D( filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
      case ('U_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
      case ('V_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
      case ('T_lad')
        call add_field_mesh_dp_2D( filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
      case ('S_lad')
        call add_field_mesh_dp_2D( filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

      ! Useful laddie fields
      case ('drho_amb')
        call add_field_mesh_dp_2D( filename, ncid, 'drho_amb', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('drho_base')
        call add_field_mesh_dp_2D( filename, ncid, 'drho_base', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('entr')
        call add_field_mesh_dp_2D( filename, ncid, 'entr', long_name = 'Entrainment rate', units = 'm s^-1')
      case ('entr_dmin')
        call add_field_mesh_dp_2D( filename, ncid, 'entr_dmin', long_name = 'Entrainment rate for Dmin', units = 'm s^-1')
      case ('SGD')
        call add_field_mesh_dp_2D( filename, ncid, 'SGD', long_name = 'Subglacial discharge rate', units = 'm s^-1')
      case ('melt')
        call add_field_mesh_dp_2D( filename, ncid, 'melt', long_name = 'Melt rate', units = 'm s^-1')
      case ('divQH')
        call add_field_mesh_dp_2D( filename, ncid, 'divQH', long_name = 'Thickness divergence', units = 'm s^-1')
      case ('divQT')
        call add_field_mesh_dp_2D( filename, ncid, 'divQT', long_name = 'Heat divergence', units = 'degC m s^-1')
      case ('divQS')
        call add_field_mesh_dp_2D( filename, ncid, 'divQS', long_name = 'Salt divergence', units = 'PSU m s^-1')
      case ('diffT')
        call add_field_mesh_dp_2D( filename, ncid, 'diffT', long_name = 'Heat diffusion', units = 'degC m s^-1')
      case ('diffS')
        call add_field_mesh_dp_2D( filename, ncid, 'diffS', long_name = 'Salt diffusion', units = 'PSU m s^-1')
      case ('viscU')
        call add_field_mesh_dp_2D_b( filename, ncid, 'viscU', long_name = 'Laddie U viscosity', units = 'm^2 s^-2')
      case ('viscV')
        call add_field_mesh_dp_2D_b( filename, ncid, 'viscV', long_name = 'Laddie V viscosity', units = 'm^2 s^-2')
      case ('T_base')
        call add_field_mesh_dp_2D( filename, ncid, 'T_base', long_name = 'Temperature at ice/ocean interface', units = 'deg C')
      case ('T_amb')
        call add_field_mesh_dp_2D( filename, ncid, 'T_amb', long_name = 'Temperature at interface with ambient ocean', units = 'deg C')
      case ('u_star')
        call add_field_mesh_dp_2D( filename, ncid, 'u_star', long_name = 'Friction velocity', units = 'm s^-1')
      case ('gamma_T')
        call add_field_mesh_dp_2D( filename, ncid, 'gamma_T', long_name = 'Heat exchange coefficient', units = 'm s^-1')
      case ('divQU')
        call add_field_mesh_dp_2D_b( filename, ncid, 'divQU', long_name = 'Laddie U divergence', units = 'm^2 s^-2')
      case ('divQV')
        call add_field_mesh_dp_2D_b( filename, ncid, 'divQV', long_name = 'Laddie V divergence', units = 'm^2 s^-2')
      case ('HU_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'HU_lad', long_name = 'Laddie HU ', units = 'm^2 s^-1')
      case ('HV_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'HV_lad', long_name = 'Laddie HV ', units = 'm^2 s^-1')

    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      case ('LMB')
        call add_field_mesh_dp_2D( filename, ncid, 'LMB', long_name = 'Lateral mass balance', units = 'm yr^-1')

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      case ('AMB')
        call add_field_mesh_dp_2D( filename, ncid, 'AMB', long_name = 'Artificial mass balance', units = 'm yr^-1')

    ! == Glacial isostatic adjustment ==
    ! ==================================

      ! Main GIA variables
      case ('dHb_next')
        call add_field_mesh_dp_2D( filename, ncid, 'dHb_next', long_name = 'Bedrock elevation difference from ELRA', units = 'm')

    ! == Tracer tracking ==
    ! =====================

      case ('age')
        call add_field_mesh_dp_3D( filename, ncid, 'age', long_name = 'Age of ice', units = 'yr')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_main_regional_output_file_mesh_field

  subroutine write_grounding_line_to_file( filename, ncid, mesh, ice)

    ! In/output variables:
    character(len=*),     intent(in   ) :: filename
    integer,              intent(in   ) :: ncid
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_grounding_line_to_file'
    real(dp)                                :: NaN
    real(dp), dimension(mesh%vi1:mesh%vi2)  :: TAF_for_GL
    integer                                 :: vi
    real(dp), dimension(:,:  ), allocatable :: CC

    ! Add routine to path
    call init_routine( routine_name)

    NaN = ieee_value( NaN, ieee_signaling_nan)

    ! Replace thickness above floatation with NaN in ice-free vertices so GL wont be found there
    do vi = mesh%vi1, mesh%vi2
      if (ice%Hi( vi) > 0.1_dp) then
        TAF_for_GL( vi) = ice%TAF( vi)
      else
        TAF_for_GL( vi) = NaN
      end if
    end do

    ! Calculate grounding line contour
    if (par%primary) allocate( CC( mesh%nE,2))
    call calc_mesh_contour( mesh, TAF_for_GL, 0._dp, CC)

    ! Write to NetCDF
    call write_contour_to_file( filename, ncid, mesh, CC, 'grounding_line')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_grounding_line_to_file

  subroutine write_calving_front_to_file( filename, ncid, mesh, ice)

    ! In/output variables:
    character(len=*),     intent(in   ) :: filename
    integer,              intent(in   ) :: ncid
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_calving_front_to_file'
    real(dp)                                :: NaN
    real(dp), dimension(mesh%vi1:mesh%vi2)  :: Hi_for_GL
    integer                                 :: vi
    real(dp), dimension(:,:  ), allocatable :: CC

    ! Add routine to path
    call init_routine( routine_name)

    NaN = ieee_value( NaN, ieee_signaling_nan)

    ! Replace ice thickness with NaN in grounded vertices so CF wont be found there
    do vi = mesh%vi1, mesh%vi2
      if (ice%TAF( vi) < 0._dp) then
        Hi_for_GL( vi) = ice%Hi( vi)
      else
        Hi_for_GL( vi) = NaN
      end if
    end do

    ! Calculate calving front contour
    if (par%primary) allocate( CC( mesh%nE,2))
    call calc_mesh_contour( mesh, Hi_for_GL, 0.05_dp, CC)

    ! Write to NetCDF
    call write_contour_to_file( filename, ncid, mesh, CC, 'calving_front')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_calving_front_to_file

  subroutine write_ice_margin_to_file( filename, ncid, mesh, ice)

    ! In/output variables:
    character(len=*),     intent(in   ) :: filename
    integer,              intent(in   ) :: ncid
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_ice_margin_to_file'
    real(dp), dimension(:,:  ), allocatable :: CC

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate ice margin contour
    if (par%primary) allocate( CC( mesh%nE,2))
    call calc_mesh_contour( mesh, ice%Hi, 0.05_dp, CC)

    ! Write to NetCDF
    call write_contour_to_file( filename, ncid, mesh, CC, 'ice_margin')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_ice_margin_to_file

  subroutine write_coastline_to_file( filename, ncid, mesh, ice)

    ! In/output variables:
    character(len=*),     intent(in   ) :: filename
    integer,              intent(in   ) :: ncid
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_coastline_to_file'
    real(dp)                                :: NaN
    real(dp), dimension(mesh%vi1:mesh%vi2)  :: water_depth_for_coastline
    integer                                 :: vi
    real(dp), dimension(:,:  ), allocatable :: CC

    ! Add routine to path
    call init_routine( routine_name)

    NaN = ieee_value( NaN, ieee_signaling_nan)

    ! Replace water depth with NaN in ice-covered vertices so coastline wont be found there
    do vi = mesh%vi1, mesh%vi2
      if (ice%Hi( vi) > 0.05_dp) then
        water_depth_for_coastline( vi) = NaN
      else
        water_depth_for_coastline( vi) = ice%SL( vi) - ice%Hb( vi)
      end if
    end do

    ! Calculate coastline contour
    if (par%primary) allocate( CC( mesh%nE,2))
    call calc_mesh_contour( mesh, water_depth_for_coastline, 0._dp, CC)

    ! Write to NetCDF
    call write_contour_to_file( filename, ncid, mesh, CC, 'coastline')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_coastline_to_file

  subroutine write_grounded_ice_contour_to_file( filename, ncid, mesh, ice)

    ! In/output variables:
    character(len=*),     intent(in   ) :: filename
    integer,              intent(in   ) :: ncid
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_grounded_ice_contour_to_file'
    integer                                 :: vi
    real(dp), dimension(mesh%vi1:mesh%vi2)  :: Hi_grounded_only
    real(dp), dimension(:,:  ), allocatable :: CC

    ! Add routine to path
    call init_routine( routine_name)

    ! Remove floating ice
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_grounded_ice( vi)) then
        Hi_grounded_only( vi) = ice%Hi( vi)
      else
        Hi_grounded_only( vi) = 0._dp
      end if
    end do

    ! Calculate grounding ice contour
    if (par%primary) allocate( CC( mesh%nE,2))
    call calc_mesh_contour( mesh, Hi_grounded_only, 0.05_dp, CC)

    ! Write to NetCDF
    call write_contour_to_file( filename, ncid, mesh, CC, 'grounded_ice_contour')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_grounded_ice_contour_to_file

  subroutine write_contour_to_file( filename, ncid, mesh, CC, var_name)

    ! In/output variables:
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(in   ) :: CC
    character(len=*),         intent(in   ) :: var_name

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_contour_to_file'
    real(dp), dimension(:,:,:), allocatable :: CC_with_time
    integer                                 :: id_dim_time, ti, id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Add "pretend" time dimension
    if (par%primary) then
      allocate( CC_with_time( mesh%nE,2,1))
      CC_with_time( :,:,1) = CC
    end if

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write to NetCDF
    call inquire_var( filename, ncid, var_name, id_var)
    call write_var_primary( filename, ncid, id_var, CC_with_time, &
      start = (/ 1, 1, ti /), count = (/ mesh%nE, 2, 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_contour_to_file

end module mesh_output_files
