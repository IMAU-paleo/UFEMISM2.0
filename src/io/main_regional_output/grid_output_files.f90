module grid_output_files

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, warning, crash
  use model_configuration, only: C
  use region_types, only: type_model_region
  use grid_types, only: type_grid
  use netcdf_io_main
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, &
    map_from_mesh_vertices_to_xy_grid_3D, map_from_mesh_vertices_to_xy_grid_2D_minval, &
    map_from_mesh_triangles_to_xy_grid_2D, map_from_mesh_triangles_to_xy_grid_3D

  implicit none

  private

  public :: create_main_regional_output_file_grid, write_to_main_regional_output_file_grid, &
            create_main_regional_output_file_grid_ROI, write_to_main_regional_output_file_grid_ROI

contains

  subroutine write_to_main_regional_output_file_grid( region)
    !< Write to the main regional output NetCDF file - grid version

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
    if (par%primary) write(0,'(A)') '   Writing to grid output file "' // colour_string( trim( region%output_filename_grid), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( region%output_filename_grid, ncid)

    ! write the time to the file
    call write_time_to_file( region%output_filename_grid, ncid, region%time)

    ! write the default data fields to the file
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'Hi')
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'Hb')
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'Hs')
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'SL')
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'u_surf')
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'v_surf')
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, 'uabs_surf')

    ! write all user-defined data fields to the file
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_01)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_02)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_03)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_04)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_05)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_06)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_07)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_08)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_09)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_10)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_11)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_12)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_13)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_14)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_15)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_16)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_17)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_18)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_19)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_20)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_21)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_22)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_23)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_24)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_25)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_26)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_27)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_28)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_29)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_30)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_31)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_32)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_33)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_34)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_35)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_36)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_37)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_38)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_39)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_40)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_41)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_42)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_43)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_44)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_45)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_46)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_47)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_48)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_49)
    call write_to_main_regional_output_file_grid_field( region, region%output_grid, region%output_filename_grid, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_main_regional_output_file_grid

  subroutine write_to_main_regional_output_file_grid_ROI( region, grid, filename)
    !< Write to the gridded output NetCDF file for a region-of-interest

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region
    type(type_grid),         intent(in   ) :: grid
    character(len=*),        intent(in   ) :: filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_main_regional_output_file_grid_ROI'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to grid output file "' // colour_string( trim( filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( filename, ncid)

    ! write the time to the file
    call write_time_to_file( filename, ncid, region%time)

    ! write the default data fields to the file
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'Hi')
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'Hb')
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'Hs')
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'SL')
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'u_surf')
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'v_surf')
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, 'uabs_surf')

    ! write all user-defined data fields to the file
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_01)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_02)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_03)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_04)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_05)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_06)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_07)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_08)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_09)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_10)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_11)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_12)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_13)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_14)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_15)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_16)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_17)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_18)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_19)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_20)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_21)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_22)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_23)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_24)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_25)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_26)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_27)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_28)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_29)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_30)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_31)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_32)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_33)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_34)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_35)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_36)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_37)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_38)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_39)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_40)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_41)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_42)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_43)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_44)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_45)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_46)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_47)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_48)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_49)
    call write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_main_regional_output_file_grid_ROI

  subroutine write_to_main_regional_output_file_grid_field( region, grid, filename, ncid, choice_output_field)
    !< Write a single field to the main regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region
    type(type_grid),         intent(in   ) :: grid
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid
    character(len=*),        intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_main_regional_output_file_grid_field'
    real(dp), dimension(:),   allocatable :: d_mesh_vec_partial_2D
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_2D_monthly
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_3D
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_3D_ocean
    real(dp), dimension(:),   allocatable :: mask_int

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! allocate memory
    allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
    allocate( d_grid_vec_partial_2D(         grid%n_loc                ))
    allocate( d_grid_vec_partial_2D_monthly( grid%n_loc, 12            ))
    allocate( d_grid_vec_partial_3D(         grid%n_loc, region%mesh%nz))
    allocate( d_grid_vec_partial_3D_ocean(   grid%n_loc, C%nz_ocean    ))
    allocate( mask_int( region%mesh%vi1:region%mesh%vi2), source = 0._dp)

    ! Add the specified data field to the file
    select case (choice_output_field)
      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')
      case ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      case ('resolution')
        d_mesh_vec_partial_2D = region%mesh%R( region%mesh%vi1:region%mesh%vi2)
        call map_from_mesh_vertices_to_xy_grid_2D_minval( region%mesh, grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'resolution', d_grid_vec_partial_2D)

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      case ('Hi_init')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_init%Hi, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hi_init', d_grid_vec_partial_2D)
      case ('Hb_init')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_init%Hb, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hb_init', d_grid_vec_partial_2D)
      case ('Hs_init')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_init%Hs, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hs_init', d_grid_vec_partial_2D)
      case ('SL_init')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_init%SL, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'SL_init', d_grid_vec_partial_2D)

      ! Present-day ice-sheet geometry
      case ('Hi_PD')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_PD%Hi, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hi_PD', d_grid_vec_partial_2D)
      case ('Hb_PD')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_PD%Hb, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hb_PD', d_grid_vec_partial_2D)
      case ('Hs_PD')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_PD%Hs, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hs_PD', d_grid_vec_partial_2D)
      case ('SL_PD')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_PD%SL, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'SL_PD', d_grid_vec_partial_2D)

      ! GIA equilibrium ice-sheet geometry
      case ('Hi_GIAeq')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_GIAeq%Hi, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hi_GIAeq', d_grid_vec_partial_2D)
      case ('Hb_GIAeq')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_GIAeq%Hb, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hb_GIAeq', d_grid_vec_partial_2D)
      case ('Hs_GIAeq')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_GIAeq%Hs, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'Hs_GIAeq', d_grid_vec_partial_2D)
      case ('SL_GIAeq')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%refgeo_GIAeq%SL, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'SL_GIAeq', d_grid_vec_partial_2D)

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      case ('Hi')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Hi, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hi', d_grid_vec_partial_2D)
      case ('Hb')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Hb, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hb', d_grid_vec_partial_2D)
      case ('Hs')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Hs, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hs', d_grid_vec_partial_2D)
      case ('Hib')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Hib, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hib', d_grid_vec_partial_2D)
      case ('SL')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%SL, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'SL', d_grid_vec_partial_2D)
      case ('TAF')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%TAF, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'TAF', d_grid_vec_partial_2D)
      case ('Hi_eff')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Hi_eff, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hi_eff', d_grid_vec_partial_2D)
      case ('Hs_slope')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Hs_slope, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Hs_slope', d_grid_vec_partial_2D)
      case ('grounding_line')
        ! Do nothing; only written to mesh files
      case ('ice_margin')
        ! Do nothing; only written to mesh files
      case ('calving_front')
        ! Do nothing; only written to mesh files
      case ('coastline')
        ! Do nothing; only written to mesh files
      case ('grounded_ice_contour')
        ! Do nothing; only written to mesh files

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      case ('dHi')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHi, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi', d_grid_vec_partial_2D)
      case ('dHb')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHb, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHb', d_grid_vec_partial_2D)
      case ('dHs')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHs, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHs', d_grid_vec_partial_2D)
      case ('dHib')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHib, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHib', d_grid_vec_partial_2D)

    ! ===== Geometry rates of change =====
    ! ====================================

      case ('dHi_dt')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHi_dt, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt', d_grid_vec_partial_2D)
      case ('dHb_dt')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHb_dt, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHb_dt', d_grid_vec_partial_2D)
      case ('dHs_dt')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHs_dt, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHs_dt', d_grid_vec_partial_2D)
      case ('dHib_dt')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHib_dt, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHib_dt', d_grid_vec_partial_2D)
      case ('dHi_dt_raw')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHi_dt_raw, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt_raw', d_grid_vec_partial_2D)
      case ('dHi_dt_residual')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHi_dt_residual, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt_residual', d_grid_vec_partial_2D)

    ! ===== Target quantities =====
    ! =============================

      case ('dHi_dt_target')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%dHi_dt_target, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHi_dt_target', d_grid_vec_partial_2D)

    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      case ('mask_icefree_land')
      case ('mask_icefree_ocean')
      case ('mask_grounded_ice')
      case ('mask_floating_ice')
      case ('mask_margin')
      case ('mask_gl_gr')
      case ('mask_gl_fl')
      case ('mask_cf_gr')
      case ('mask_cf_fl')
      case ('mask_coastline')
      case ('mask')
      case ('basin_ID')
      case ('mask_SGD')
      case ('mask_ROI')
        ! Exception for mask_ROI, needed for offline laddie coupling
        where (region%ice%mask_ROI .eqv. .TRUE.)
          mask_int = 1.0_dp
        elsewhere
          mask_int = 0.0_dp
        end where
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, mask_int, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'mask_ROI', d_grid_vec_partial_2D)

    ! ===== Area fractions =====
    ! ==========================

      ! notE: sub-grid area fractions cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      case ('fraction_gr')
      case ('fraction_gr_b')
      case ('fraction_margin')

    ! === Thermodynamics and rheology ===
    ! ===================================

      case ('Ti')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%Ti, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Ti', d_grid_vec_partial_3D)
      case ('Ti_pmp')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%Ti_pmp, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Ti_pmp', d_grid_vec_partial_3D)
      case ('Ti_hom')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Ti_hom, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'Ti_hom', d_grid_vec_partial_2D)
      case ('Cpi')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%Cpi, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Cpi', d_grid_vec_partial_3D)
      case ('Ki')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%Ki, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'Ki', d_grid_vec_partial_3D)
      case ('internal_heating')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%internal_heating, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'internal_heating', d_grid_vec_partial_3D)
      case ('frictional_heating')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%frictional_heating, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'frictional_heating', d_grid_vec_partial_2D)
      case ('A_flow')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%A_flow, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'A_flow', d_grid_vec_partial_3D)

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      case ('u_3D')
        call map_from_mesh_triangles_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%u_3D_b, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'u_3D', d_grid_vec_partial_3D)
      case ('v_3D')
        call map_from_mesh_triangles_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%v_3D_b, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'v_3D', d_grid_vec_partial_3D)
      case ('w_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%w_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'w_3D', d_grid_vec_partial_3D)

      ! Vertically integrated
      case ('u_vav')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%u_vav_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'u_vav', d_grid_vec_partial_2D)
      case ('v_vav')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%v_vav_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'v_vav', d_grid_vec_partial_2D)
      case ('uabs_vav')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%uabs_vav_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'uabs_vav', d_grid_vec_partial_2D)

      ! Surface
      case ('u_surf')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%u_surf_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'u_surf', d_grid_vec_partial_2D)
      case ('v_surf')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%v_surf_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'v_surf', d_grid_vec_partial_2D)
      case ('w_surf')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%w_surf, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'w_surf', d_grid_vec_partial_2D)
      case ('uabs_surf')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%uabs_surf_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'uabs_surf', d_grid_vec_partial_2D)

      ! Base
      case ('u_base')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%u_base_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'u_base', d_grid_vec_partial_2D)
      case ('v_base')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%v_base_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'v_base', d_grid_vec_partial_2D)
      case ('w_base')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%w_base, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'w_base', d_grid_vec_partial_2D)
      case ('uabs_base')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%uabs_base_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'uabs_base', d_grid_vec_partial_2D)

    ! === Strain rates ===
    ! ====================

      case ('du_dx_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%du_dx_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'du_dx_3D', d_grid_vec_partial_3D)
      case ('du_dy_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%du_dy_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'du_dy_3D', d_grid_vec_partial_3D)
      case ('du_dz_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%du_dz_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'du_dz_3D', d_grid_vec_partial_3D)
      case ('dv_dx_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%dv_dx_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dv_dx_3D', d_grid_vec_partial_3D)
      case ('dv_dy_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%dv_dy_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dv_dy_3D', d_grid_vec_partial_3D)
      case ('dv_dz_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%dv_dz_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dv_dz_3D', d_grid_vec_partial_3D)
      case ('dw_dx_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%dw_dx_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dw_dx_3D', d_grid_vec_partial_3D)
      case ('dw_dy_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%dw_dy_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dw_dy_3D', d_grid_vec_partial_3D)
      case ('dw_dz_3D')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ice%dw_dz_3D, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'dw_dz_3D', d_grid_vec_partial_3D)

    ! == Ice flow regime ==
    ! =====================

      case ('divQ')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%divQ, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'divQ', d_grid_vec_partial_2D)
      case ('R_shear')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%R_shear, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'R_shear', d_grid_vec_partial_2D)

    ! == Ice P/C time stepping ==
    ! ===========================

      case ('pc_truncation_error')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%pc%tau_np1, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pc_truncation_error', d_grid_vec_partial_2D)
      case ('pc_untolerated_events')
        ! DENK DROM : not gridable

    ! == Basal hydrology ==
    ! =====================

      case ('pore_water_pressure')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%pore_water_pressure, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pore_water_pressure', d_grid_vec_partial_2D)
      case ('overburden_pressure')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%overburden_pressure, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'overburden_pressure', d_grid_vec_partial_2D)
      case ('effective_pressure')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%effective_pressure, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'effective_pressure', d_grid_vec_partial_2D)
      case ('pore_water_likelihood')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%pore_water_likelihood, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pore_water_likelihood', d_grid_vec_partial_2D)
      case ('pore_water_fraction')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%pore_water_fraction, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'pore_water_fraction', d_grid_vec_partial_2D)

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      case ('till_friction_angle')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%till_friction_angle, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'till_friction_angle', d_grid_vec_partial_2D)
      case ('alpha_sq')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%alpha_sq, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'alpha_sq', d_grid_vec_partial_2D)
      case ('beta_sq')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%beta_sq, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'beta_sq', d_grid_vec_partial_2D)

      ! Basal friction and shear stress
      case ('till_yield_stress')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%till_yield_stress, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'till_yield_stress', d_grid_vec_partial_2D)
      case ('basal_friction_coefficient')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%basal_friction_coefficient, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'basal_friction_coefficient', d_grid_vec_partial_2D)
      case ('basal_shear_stress')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%basal_shear_stress, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'basal_shear_stress', d_grid_vec_partial_2D)

      ! Bed roughness nudging - H, dH/dt, flowline
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_dHdt_flowline%deltaHs_av_up, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_dHdt_flowline%deltaHs_av_down, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_dHdt_flowline%dHs_dt_av_up, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_dHdt_flowline%dHs_dt_av_down, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_dHdt_flowline_R')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_dHdt_flowline%R, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_dHdt_flowline_R', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_dHdt_flowline_I_tot')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_dHdt_flowline%I_tot, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_dHdt_flowline_I_tot', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_dHdt_flowline_dC_dt')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_dHdt_flowline%dC_dt, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_dHdt_flowline_dC_dt', d_grid_vec_partial_2D)

      ! Bed roughness nudging - H, u, flowline
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_up')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%deltaHs_av_up, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_flowline_deltaHs_av_up', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_down')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%deltaHs_av_down, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_flowline_deltaHs_av_down', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_up')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%deltau_av_up, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_flowline_deltau_av_up', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_down')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%deltau_av_down, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_flowline_deltau_av_down', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_u_flowline_R')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%R, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_flowline_R', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_u_flowline_I_tot')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%I_tot, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_flowline_I_tot', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_u_flowline_dC_dt')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%dC_dt, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_flowline_dC_dt', d_grid_vec_partial_2D)
      case ('bed_roughness_nudge_H_u_target_velocity')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%bed_roughness%nudging_H_u_flowline%uabs_surf_target_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'bed_roughness_nudge_H_u_target_velocity', d_grid_vec_partial_2D)

    ! == Geothermal heat ==
    ! =====================

      case ('geothermal_heat_flux')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%geothermal_heat_flux, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'geothermal_heat_flux', d_grid_vec_partial_2D)

    ! == Climate ==
    ! =============

      ! Main climate variables
      case ('T2m')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%climate%T2m, d_grid_vec_partial_2D_monthly)
        call write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, 'T2m', d_grid_vec_partial_2D_monthly)
      case ('Precip')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%climate%Precip, d_grid_vec_partial_2D_monthly)
        call write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, 'Precip', d_grid_vec_partial_2D_monthly)
      CASE ('Q_TOA')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%climate%snapshot%Q_TOA, d_grid_vec_partial_2D_monthly)
        call write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, 'Q_TOA', d_grid_vec_partial_2D_monthly)

    ! == Ocean ==
    ! ==========================

      ! Main ocean variables
      case ('T_ocean')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ocean%T, d_grid_vec_partial_3D_ocean)
        call write_to_field_multopt_grid_dp_3D_ocean( grid, filename, ncid, 'T_ocean', d_grid_vec_partial_3D_ocean)
      case ('S_ocean')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%ocean%S, d_grid_vec_partial_3D_ocean)
        call write_to_field_multopt_grid_dp_3D_ocean( grid, filename, ncid, 'S_ocean', d_grid_vec_partial_3D_ocean)
      case ('T_draft')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ocean%T_draft, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'T_draft', d_grid_vec_partial_2D)
      case ('T_freezing_point')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ocean%T_freezing_point, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'T_freezing_point', d_grid_vec_partial_2D)

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      case ('SMB')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%SMB%SMB, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'SMB', d_grid_vec_partial_2D)
      case ('Albedo')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%SMB%IMAUITM%Albedo, d_grid_vec_partial_2D_monthly)
        call write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, 'Albedo', d_grid_vec_partial_2D_monthly)
      case ('FirnDepth')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%SMB%IMAUITM%FirnDepth, d_grid_vec_partial_2D_monthly)
        call write_to_field_multopt_grid_dp_2D_monthly( grid, filename, ncid, 'FirnDepth', d_grid_vec_partial_2D_monthly)
      case ('MeltPreviousYear')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%SMB%IMAUITM%MeltPreviousYear, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'MeltPreviousYear', d_grid_vec_partial_2D)

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      case ('BMB')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%BMB, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'BMB', d_grid_vec_partial_2D)
      case ('BMB_inv')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%BMB_inv, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'BMB_inv', d_grid_vec_partial_2D)
      case ('BMB_transition_phase')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%BMB_transition_phase, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'BMB_transition_phase', d_grid_vec_partial_2D)
      case ('BMB_modelled')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%BMB_modelled, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'BMB_modelled', d_grid_vec_partial_2D)

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%now%H, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'H_lad', d_grid_vec_partial_2D)
      case ('U_lad')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%now%U, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'U_lad', d_grid_vec_partial_2D)
      case ('V_lad')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%now%V, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'V_lad', d_grid_vec_partial_2D)
      case ('T_lad')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%now%T, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'T_lad', d_grid_vec_partial_2D)
      case ('S_lad')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%now%S, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'S_lad', d_grid_vec_partial_2D)

      ! Useful laddie fields
      case ('drho_amb')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%drho_amb, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'drho_amb', d_grid_vec_partial_2D)
      case ('drho_base')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%drho_base, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'drho_base', d_grid_vec_partial_2D)
      case ('entr')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%entr, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'entr', d_grid_vec_partial_2D)
      case ('entr_dmin')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%entr_dmin, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'entr_dmin', d_grid_vec_partial_2D)
      case ('SGD')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%SGD, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'SGD', d_grid_vec_partial_2D)
      case ('melt')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%melt, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'melt', d_grid_vec_partial_2D)
      case ('divQH')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%divQH, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'divQH', d_grid_vec_partial_2D)
      case ('divQT')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%divQT, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'divQT', d_grid_vec_partial_2D)
      case ('divQS')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%divQS, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'divQS', d_grid_vec_partial_2D)
      case ('diffT')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%diffT, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'diffT', d_grid_vec_partial_2D)
      case ('diffS')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%diffS, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'diffS', d_grid_vec_partial_2D)
      case ('viscU')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%viscU, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'viscU', d_grid_vec_partial_2D)
      case ('viscV')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%viscV, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'viscV', d_grid_vec_partial_2D)
      case ('T_base')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%T_base, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'T_base', d_grid_vec_partial_2D)
      case ('T_amb')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%T_amb, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'T_amb', d_grid_vec_partial_2D)
      case ('u_star')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%u_star, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'u_star', d_grid_vec_partial_2D)
      case ('gamma_T')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%gamma_T, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'gamma_T', d_grid_vec_partial_2D)
      case ('divQU')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%divQU, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'divQU', d_grid_vec_partial_2D)
      case ('divQV')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%BMB%laddie%divQV, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'divQV', d_grid_vec_partial_2D)
      case ('HU_lad')
        ! Not implemented
      case ('HV_lad')
        ! Not implemented

    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      case ('LMB')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%LMB%LMB, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'LMB', d_grid_vec_partial_2D)

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      case ('AMB')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%AMB%AMB, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'AMB', d_grid_vec_partial_2D)

    ! == Glacial isostatic adjustment ==
    ! ==================================

      ! Main GIA variables
      case ('dHb_next')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%GIA%dHb_next, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'dHb_next', d_grid_vec_partial_2D)

    ! == Tracer tracking ==
    ! =====================

      case ('age')
        call map_from_mesh_vertices_to_xy_grid_3D( region%mesh, grid, C%output_dir, region%tracer_tracking%age, d_grid_vec_partial_3D)
        call write_to_field_multopt_grid_dp_3D( grid, filename, ncid, 'age', d_grid_vec_partial_3D)

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_main_regional_output_file_grid_field

  subroutine create_main_regional_output_file_grid( region)
    !< Create the main regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_main_regional_output_file_grid'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    region%output_filename_grid = trim( C%output_dir) // 'main_output_' // region%name // '_grid.nc'

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating grid output file "' // colour_string( trim( region%output_filename_grid), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( region%output_filename_grid, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( region%output_filename_grid, ncid, region%output_grid)

    ! Add time, zeta, and month dimensions+variables to the file
    call add_time_dimension_to_file(  region%output_filename_grid, ncid)
    call add_zeta_dimension_to_file(  region%output_filename_grid, ncid, region%mesh%zeta)
    call add_month_dimension_to_file( region%output_filename_grid, ncid)
    call add_depth_dimension_to_file( region%output_filename_grid, ncid, C%z_ocean)

    ! Add the default data fields to the file
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'Hi')
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'Hb')
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'Hs')
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'SL')
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'u_surf')
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'v_surf')
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, 'uabs_surf')

    ! Add all user-defined data fields to the file
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_01)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_02)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_03)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_04)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_05)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_06)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_07)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_08)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_09)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_10)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_11)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_12)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_13)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_14)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_15)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_16)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_17)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_18)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_19)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_20)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_21)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_22)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_23)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_24)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_25)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_26)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_27)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_28)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_29)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_30)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_31)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_32)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_33)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_34)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_35)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_36)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_37)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_38)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_39)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_40)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_41)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_42)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_43)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_44)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_45)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_46)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_47)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_48)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_49)
    call create_main_regional_output_file_grid_field( region%output_filename_grid, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_main_regional_output_file_grid

  subroutine create_main_regional_output_file_grid_ROI( region, grid, filename)
    !< Create the main regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region
    type(type_grid),         intent(in   ) :: grid
    character(len=*),        intent(in   ) :: filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_main_regional_output_file_grid_ROI'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating ROI output file "' // colour_string( trim( filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    ! Add time, zeta, and month dimensions+variables to the file
    call add_time_dimension_to_file(  filename, ncid)
    call add_zeta_dimension_to_file(  filename, ncid, region%mesh%zeta)
    call add_month_dimension_to_file( filename, ncid)
    call add_depth_dimension_to_file( filename, ncid, C%z_ocean)

    ! Add the default data fields to the file
    call create_main_regional_output_file_grid_field( filename, ncid, 'Hi')
    call create_main_regional_output_file_grid_field( filename, ncid, 'Hb')
    call create_main_regional_output_file_grid_field( filename, ncid, 'Hs')
    call create_main_regional_output_file_grid_field( filename, ncid, 'SL')
    call create_main_regional_output_file_grid_field( filename, ncid, 'u_surf')
    call create_main_regional_output_file_grid_field( filename, ncid, 'v_surf')
    call create_main_regional_output_file_grid_field( filename, ncid, 'uabs_surf')

    ! Add all user-defined data fields to the file
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_01)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_02)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_03)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_04)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_05)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_06)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_07)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_08)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_09)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_10)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_11)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_12)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_13)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_14)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_15)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_16)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_17)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_18)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_19)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_20)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_21)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_22)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_23)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_24)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_25)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_26)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_27)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_28)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_29)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_30)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_31)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_32)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_33)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_34)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_35)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_36)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_37)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_38)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_39)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_40)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_41)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_42)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_43)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_44)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_45)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_46)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_47)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_48)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_49)
    call create_main_regional_output_file_grid_field( filename, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_main_regional_output_file_grid_ROI

  subroutine create_main_regional_output_file_grid_field( filename, ncid, choice_output_field)
    !< Create a single field in the main regional output NetCDF file - grid version

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_main_regional_output_file_grid_field'

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
        call add_field_grid_dp_2D_notime( filename, ncid, 'resolution', long_name = 'Mesh resolution (distance to nearest neighbour)', units = 'm')

    ! ===== Reference geometries =====
    ! ================================

      ! Initial ice-sheet geometry
      case ('Hi_init')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hi_init', long_name = 'Initial ice thickness', units = 'm')
      case ('Hb_init')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hb_init', long_name = 'Initial bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs_init')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hs_init', long_name = 'Initial surface elevation', units = 'm w.r.t. PD sea level')
      case ('SL_init')
        call add_field_grid_dp_2D_notime( filename, ncid, 'SL_init', long_name = 'Initial geoid elevation', units = 'm w.r.t. PD sea level')

      ! Present-day ice-sheet geometry
      case ('Hi_PD')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hi_PD', long_name = 'Present-day ice thickness', units = 'm')
      case ('Hb_PD')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hb_PD', long_name = 'Present-day bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs_PD')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hs_PD', long_name = 'Present-day surface elevation', units = 'm w.r.t. PD sea level')
      case ('SL_PD')
        call add_field_grid_dp_2D_notime( filename, ncid, 'SL_PD', long_name = 'Present-day geoid elevation', units = 'm w.r.t. PD sea level')

      ! GIA equilibrium ice-sheet geometry
      case ('Hi_GIAeq')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hi_GIAeq', long_name = 'GIA equilibrium ice thickness', units = 'm')
      case ('Hb_GIAeq')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hb_GIAeq', long_name = 'GIA equilibrium bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs_GIAeq')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hs_GIAeq', long_name = 'GIA equilibrium surface elevation', units = 'm w.r.t. PD sea level')
      case ('SL_GIAeq')
        call add_field_grid_dp_2D_notime( filename, ncid, 'SL_GIAeq', long_name = 'GIA equilibrium geoid elevation', units = 'm w.r.t. PD sea level')

    ! ===== Basic ice-sheet geometry =====
    ! ====================================

      case ('Hi')
        call add_field_grid_dp_2D( filename, ncid, 'Hi', long_name = 'Ice thickness', units = 'm')
      case ('Hb')
        call add_field_grid_dp_2D( filename, ncid, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. PD sea level')
      case ('Hs')
        call add_field_grid_dp_2D( filename, ncid, 'Hs', long_name = 'Surface elevation', units = 'm w.r.t. PD sea level')
      case ('Hib')
        call add_field_grid_dp_2D( filename, ncid, 'Hib', long_name = 'Ice base elevation', units = 'm w.r.t. PD sea level')
      case ('SL')
        call add_field_grid_dp_2D( filename, ncid, 'SL', long_name = 'Geoid elevation', units = 'm w.r.t. PD sea level')
      case ('TAF')
        call add_field_grid_dp_2D( filename, ncid, 'TAF', long_name = 'Thickness above floatation', units = 'm')
      case ('Hi_eff')
        call add_field_grid_dp_2D( filename, ncid, 'Hi_eff', long_name = 'Effective ice thickness', units = 'm')
      case ('Hs_slope')
        call add_field_grid_dp_2D( filename, ncid, 'Hs_slope', long_name = 'Absolute surface gradient', units = '-')
      case ('grounding_line')
        ! Do nothing; only written to mesh files
      case ('ice_margin')
        ! Do nothing; only written to mesh files
      case ('calving_front')
        ! Do nothing; only written to mesh files
      case ('coastline')
        ! Do nothing; only written to mesh files
      case ('grounded_ice_contour')
        ! Do nothing; only written to mesh files

    ! ===== Geometry changes w.r.t. reference =====
    ! =============================================

      case ('dHi')
        call add_field_grid_dp_2D( filename, ncid, 'dHi', long_name = 'Ice thickness difference w.r.t. reference', units = 'm')
      case ('dHb')
        call add_field_grid_dp_2D( filename, ncid, 'dHb', long_name = 'Bedrock elevation difference w.r.t. reference', units = 'm')
      case ('dHs')
        call add_field_grid_dp_2D( filename, ncid, 'dHs', long_name = 'Surface elevation difference w.r.t. reference', units = 'm')
      case ('dHib')
        call add_field_grid_dp_2D( filename, ncid, 'dHib', long_name = 'Ice base elevation difference w.r.t. reference', units = 'm')

    ! ===== Geometry rates of change =====
    ! ====================================

      case ('dHi_dt')
        call add_field_grid_dp_2D( filename, ncid, 'dHi_dt', long_name = 'Ice thickness rate of change', units = 'm yr^-1')
      case ('dHb_dt')
        call add_field_grid_dp_2D( filename, ncid, 'dHb_dt', long_name = 'Bedrock elevation rate of change', units = 'm yr^-1')
      case ('dHs_dt')
        call add_field_grid_dp_2D( filename, ncid, 'dHs_dt', long_name = 'Surface elevation rate of change', units = 'm yr^-1')
      case ('dHib_dt')
        call add_field_grid_dp_2D( filename, ncid, 'dHib_dt', long_name = 'Ice base elevation rate of change', units = 'm yr^-1')
      case ('dHi_dt_raw')
        call add_field_grid_dp_2D( filename, ncid, 'dHi_dt_raw', long_name = 'Ice thickness rate of change before any modifications', units = 'm yr^-1')
      case ('dHi_dt_residual')
        call add_field_grid_dp_2D( filename, ncid, 'dHi_dt_residual', long_name = 'Residual ice thickness rate of change during model calibration', units = 'm yr^-1')

    ! ===== Target quantities =====
    ! =============================

      case ('dHi_dt_target')
        call add_field_grid_dp_2D( filename, ncid, 'dHi_dt_target', long_name = 'Target ice thickness rate of change during model calibration', units = 'm yr^-1')
      case ('uabs_surf_target')
        call add_field_grid_dp_2D( filename, ncid, 'uabs_surf_target', long_name = 'Target ice surface speed during model calibration', units = 'm yr^-1')

    ! ===== Masks =====
    ! =================

      ! notE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      case ('mask_icefree_land')
      case ('mask_icefree_ocean')
      case ('mask_grounded_ice')
      case ('mask_floating_ice')
      case ('mask_margin')
      case ('mask_gl_gr')
      case ('mask_gl_fl')
      case ('mask_cf_gr')
      case ('mask_cf_fl')
      case ('mask_coastline')
      case ('mask')
      case ('basin_ID')
      case ('mask_SGD')
      case ('mask_ROI')
        ! Exception for mask_ROI, needed for offline laddie computation
        call add_field_grid_dp_2D( filename, ncid, 'mask_ROI', long_name = 'ROI mask', units = '')


    ! ===== Area fractions =====
    ! ==========================

      ! notE: sub-grid area fractions cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      case ('fraction_gr')
      case ('fraction_gr_b')
      case ('fraction_margin')

    ! === Thermodynamics and rheology ===
    ! ===================================

      case ('Ti')
        call add_field_grid_dp_3D( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')
      case ('Ti_pmp')
        call add_field_grid_dp_3D( filename, ncid, 'Ti_pmp', long_name = 'Pressure melting point temperature', units = 'K')
      case ('Ti_hom')
        call add_field_grid_dp_2D( filename, ncid, 'Ti_hom', long_name = 'Temperature at base w.r.t. pressure melting point', units = 'K')
      case ('Cpi')
        call add_field_grid_dp_3D( filename, ncid, 'Cpi', long_name = 'Specific heat capacity', units = 'J kg^-1 K^-1')
      case ('Ki')
        call add_field_grid_dp_3D( filename, ncid, 'Ki', long_name = 'Thermal conductivity', units = 'J m^-1 K^-1 yr^-1')
      case ('internal_heating')
        call add_field_grid_dp_3D( filename, ncid, 'internal_heating', long_name = 'Internal heating', units = '?')
      case ('frictional_heating')
        call add_field_grid_dp_2D( filename, ncid, 'frictional_heating', long_name = 'Frictional heating', units = '?')
      case ('A_flow')
        call add_field_grid_dp_3D( filename, ncid, 'A_flow', long_name = 'Glens flow law factor', units = 'Pa^-3 y^-1')

    ! === Ice velocities ===
    ! ======================

      ! 3-D
      case ('u_3D')
        call add_field_grid_dp_3D( filename, ncid, 'u_3D', long_name = '3-D ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_3D')
        call add_field_grid_dp_3D( filename, ncid, 'v_3D', long_name = '3-D ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_3D_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('v_3D_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('w_3D')
        call add_field_grid_dp_3D( filename, ncid, 'w_3D', long_name = '3-D ice velocity in the z-direction', units = 'm yr^-1')

      ! Vertically integrated
      case ('u_vav')
        call add_field_grid_dp_2D( filename, ncid, 'u_vav', long_name = 'Vertically averaged ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_vav')
        call add_field_grid_dp_2D( filename, ncid, 'v_vav', long_name = 'Vertically averaged ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_vav_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('v_vav_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('uabs_vav')
        call add_field_grid_dp_2D( filename, ncid, 'uabs_vav', long_name = 'Vertically averaged absolute ice velocity', units = 'm yr^-1')
      case ('uabs_vav_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!

      ! Surface
      case ('u_surf')
        call add_field_grid_dp_2D( filename, ncid, 'u_surf', long_name = 'Surface ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_surf')
        call add_field_grid_dp_2D( filename, ncid, 'v_surf', long_name = 'Surface ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_surf_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('v_surf_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('w_surf')
        call add_field_grid_dp_2D( filename, ncid, 'w_surf', long_name = 'Surface ice velocity in the z-direction', units = 'm yr^-1')
      case ('uabs_surf')
        call add_field_grid_dp_2D( filename, ncid, 'uabs_surf', long_name = 'Absolute surface ice velocity', units = 'm yr^-1')
      case ('uabs_surf_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!

      ! Base
      case ('u_base')
        call add_field_grid_dp_2D( filename, ncid, 'u_base', long_name = 'Basal ice velocity in the x-direction', units = 'm yr^-1')
      case ('v_base')
        call add_field_grid_dp_2D( filename, ncid, 'v_base', long_name = 'Basal ice velocity in the y-direction', units = 'm yr^-1')
      case ('u_base_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('v_base_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!
      case ('w_base')
        call add_field_grid_dp_2D( filename, ncid, 'w_base', long_name = 'Basal ice velocity in the z-direction', units = 'm yr^-1')
      case ('uabs_base')
        call add_field_grid_dp_2D( filename, ncid, 'uabs_base', long_name = 'Absolute basal ice velocity', units = 'm yr^-1')
      case ('uabs_base_b')
        ! notE: mapping from mesh triangles to square grid is not (yet) available!

    ! === Strain rates ===
    ! ====================

      case ('du_dx_3D')
        call add_field_grid_dp_3D( filename, ncid, 'du_dx_3D', long_name = '3-D xx strain rate', units = 'yr^-1')
      case ('du_dy_3D')
        call add_field_grid_dp_3D( filename, ncid, 'du_dy_3D', long_name = '3-D xy strain rate', units = 'yr^-1')
      case ('du_dz_3D')
        call add_field_grid_dp_3D( filename, ncid, 'du_dz_3D', long_name = '3-D xz strain rate', units = 'yr^-1')
      case ('dv_dx_3D')
        call add_field_grid_dp_3D( filename, ncid, 'dv_dx_3D', long_name = '3-D yx strain rate', units = 'yr^-1')
      case ('dv_dy_3D')
        call add_field_grid_dp_3D( filename, ncid, 'dv_dy_3D', long_name = '3-D yy strain rate', units = 'yr^-1')
      case ('dv_dz_3D')
        call add_field_grid_dp_3D( filename, ncid, 'dv_dz_3D', long_name = '3-D yz strain rate', units = 'yr^-1')
      case ('dw_dx_3D')
        call add_field_grid_dp_3D( filename, ncid, 'dw_dx_3D', long_name = '3-D zx strain rate', units = 'yr^-1')
      case ('dw_dy_3D')
        call add_field_grid_dp_3D( filename, ncid, 'dw_dy_3D', long_name = '3-D zy strain rate', units = 'yr^-1')
      case ('dw_dz_3D')
        call add_field_grid_dp_3D( filename, ncid, 'dw_dz_3D', long_name = '3-D zz strain rate', units = 'yr^-1')

    ! == Ice flow regime ==
    ! =====================

      case ('divQ')
        call add_field_grid_dp_2D( filename, ncid, 'divQ', long_name = 'Horizontal ice flux divergence', units = 'm yr^-1')
      case ('R_shear')
        call add_field_grid_dp_2D( filename, ncid, 'R_shear', long_name = 'Slide/shear ratio', units = '0-1')

    ! == Ice P/C time stepping ==
    ! ===========================

      case ('pc_truncation_error')
        call add_field_grid_dp_2D( filename, ncid, 'pc_truncation_error', long_name = 'Ice P/C truncation error tau', units = 'm')
      case ('pc_untolerated_events')
        ! DENK DROM : not gridable

    ! == Basal hydrology ==
    ! =====================

      case ('pore_water_pressure')
        call add_field_grid_dp_2D( filename, ncid, 'pore_water_pressure', long_name = 'Till pore water pressure', units = 'Pa')
      case ('overburden_pressure')
        call add_field_grid_dp_2D( filename, ncid, 'overburden_pressure', long_name = 'Ice overburden pressure', units = 'Pa')
      case ('effective_pressure')
        call add_field_grid_dp_2D( filename, ncid, 'effective_pressure', long_name = 'Effective basal pressure', units = 'Pa')
      case ('pore_water_likelihood')
        call add_field_grid_dp_2D( filename, ncid, 'pore_water_likelihood', long_name = 'Till pore water likelihood', units = '0-1')
      case ('pore_water_fraction')
        call add_field_grid_dp_2D( filename, ncid, 'pore_water_fraction', long_name = 'Fraction of overburden pressure reduced by pore water', units = '0-1')

    ! == Basal sliding ==
    ! ===================

      ! Sliding law coefficients
      case ('till_friction_angle')
        call add_field_grid_dp_2D( filename, ncid, 'till_friction_angle', long_name = 'Till friction angle', units = 'degrees')
      case ('alpha_sq')
        call add_field_grid_dp_2D( filename, ncid, 'alpha_sq', long_name = 'Coulomb-law friction coefficientn', units = 'dimensionless')
      case ('beta_sq')
        call add_field_grid_dp_2D( filename, ncid, 'beta_sq', long_name = 'Power-law friction coefficient', units = 'Pa m^−1/m yr^1/m')

        ! Basal friction and shear stress
      case ('till_yield_stress')
        call add_field_grid_dp_2D( filename, ncid, 'till_yield_stress', long_name = 'Till yield stress', units = 'Pa')
      case ('basal_friction_coefficient')
        call add_field_grid_dp_2D( filename, ncid, 'basal_friction_coefficient', long_name = 'Basal friction coefficient', units = 'Pa yr m^-1')
      case ('basal_shear_stress')
        call add_field_grid_dp_2D( filename, ncid, 'basal_shear_stress', long_name = 'Basal shear stress', units = 'Pa')

      ! Bed roughness nudging - H, dH/dt, flowline
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_up', &
          long_name = 'Upstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_deltaHs_av_down', &
          long_name = 'Downstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_up', &
          long_name = 'Upstream flowline-averaged thinning rate', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dHs_dt_av_down', &
          long_name = 'Downstream flowline-averaged thinning rate', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_dHdt_flowline_R')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_R', &
          long_name = 'Ice flux-based scaling factor')
      case ('bed_roughness_nudge_H_dHdt_flowline_I_tot')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_I_tot', &
          long_name = 'Weighted average of flowline-averaged terms')
      case ('bed_roughness_nudge_H_dHdt_flowline_dC_dt')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_dHdt_flowline_dC_dt', &
          long_name = 'Bed roughness rate of change')

      ! Bed roughness nudging - H, u, flowline
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_up')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltaHs_av_up', &
          long_name = 'Upstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_u_flowline_deltaHs_av_down')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltaHs_av_down', &
          long_name = 'Downstream flowline-averaged thickness error', units = 'm')
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_up')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltau_av_up', &
          long_name = 'Upstream flowline-averaged velocity error', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_u_flowline_deltau_av_down')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_deltau_av_down', &
          long_name = 'Downstream flowline-averaged velocity error', units = 'm yr^-1')
      case ('bed_roughness_nudge_H_u_flowline_R')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_R', &
          long_name = 'Ice flux-based scaling factor')
      case ('bed_roughness_nudge_H_u_flowline_I_tot')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_I_tot', &
          long_name = 'Weighted average of flowline-averaged terms')
      case ('bed_roughness_nudge_H_u_flowline_dC_dt')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_flowline_dC_dt', &
          long_name = 'Bed roughness rate of change')
      case ('bed_roughness_nudge_H_u_target_velocity')
        call add_field_grid_dp_2D( filename, ncid, &
          'bed_roughness_nudge_H_u_target_velocity', &
          long_name = 'Target velocity', units = 'm yr^-1')

    ! == Geothermal heat ==
    ! =====================

      case ('geothermal_heat_flux')
        call add_field_grid_dp_2D( filename, ncid, 'geothermal_heat_flux', long_name = 'Geothermal heat flux', units = 'J m^-2 yr^-1')

    ! == Climate ==
    ! =============

      ! Main climate variables
      case ('T2m')
        call add_field_grid_dp_2D_monthly( filename, ncid, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
      case ('Precip')
        call add_field_grid_dp_2D_monthly( filename, ncid, 'Precip', long_name = 'Monthly total precipitation', units = 'm.w.e.')
      CASE ('Q_TOA')
        CALL add_field_grid_dp_2D_monthly( filename, ncid, 'Q_TOA', long_name = 'Monthly insolation at the top of the atmosphere', units = 'W m^-2')

    ! == Ocean ==
    ! ==========================

      ! Main ocean variables
      case ('T_ocean')
        call add_field_grid_dp_3D_ocean( filename, ncid, 'T_ocean', long_name = 'Ocean temperature', units = 'deg C')
      case ('S_ocean')
        call add_field_grid_dp_3D_ocean( filename, ncid, 'S_ocean', long_name = 'Ocean salinity', units = 'psu')
      case ('T_draft')
        call add_field_grid_dp_2D( filename, ncid, 'T_draft', long_name = 'Ocean temperature at ice draft', units = 'deg C')
      case ('T_freezing_point')
        call add_field_grid_dp_2D( filename, ncid, 'T_freezing_point', long_name = 'Ocean freezing temperature at ice draft', units = 'deg C')

    ! == Surface mass balance ==
    ! ==========================

      ! Main SMB variables
      case ('SMB')
        call add_field_grid_dp_2D( filename, ncid, 'SMB', long_name = 'Surface mass balance', units = 'm yr^-1')
      case ('Albedo')
        call add_field_grid_dp_2D_monthly( filename, ncid, 'Albedo', long_name = 'Surface albedo', units = '0-1')
      case ('FirnDepth')
        call add_field_grid_dp_2D_monthly( filename, ncid, 'FirnDepth', long_name = 'Monthly firn layer depth', units = 'm')
      case ('MeltPreviousYear')
        call add_field_grid_dp_2D_monthly( filename, ncid, 'MeltPreviousYear', long_name = 'Total ice melt from previous year', units = 'm')

    ! == Basal mass balance ==
    ! ========================

      ! Main BMB variables
      case ('BMB')
        call add_field_grid_dp_2D( filename, ncid, 'BMB', long_name = 'Basal mass balance', units = 'm yr^-1')
      case ('BMB_inv')
        call add_field_grid_dp_2D( filename, ncid, 'BMB_inv', long_name = 'Basal mass balance - inverted', units = 'm yr^-1')
      case ('BMB_transition_phase')
        call add_field_grid_dp_2D( filename, ncid, 'BMB_transition_phase', long_name = 'Basal mass balance - transition phase', units = 'm yr^-1')
      case ('BMB_modelled')
        call add_field_grid_dp_2D( filename, ncid, 'BMB_modelled', long_name = 'Basal mass balance - modelled', units = 'm yr^-1')

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call add_field_grid_dp_2D( filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
      case ('U_lad')
        call add_field_grid_dp_2D( filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
      case ('V_lad')
        call add_field_grid_dp_2D( filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
      case ('T_lad')
        call add_field_grid_dp_2D( filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
      case ('S_lad')
        call add_field_grid_dp_2D( filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

      ! Useful laddie fields
      case ('drho_amb')
        call add_field_grid_dp_2D( filename, ncid, 'drho_amb', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('drho_base')
        call add_field_grid_dp_2D( filename, ncid, 'drho_base', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('entr')
        call add_field_grid_dp_2D( filename, ncid, 'entr', long_name = 'Entrainment rate', units = 'm s^-1')
      case ('entr_dmin')
        call add_field_grid_dp_2D( filename, ncid, 'entr_dmin', long_name = 'Entrainment rate for Dmin', units = 'm s^-1')
      case ('SGD')
        call add_field_grid_dp_2D( filename, ncid, 'SGD', long_name = 'Subglacial discharge rate', units = 'm s^-1')
      case ('melt')
        call add_field_grid_dp_2D( filename, ncid, 'melt', long_name = 'melt rate', units = 'm s^-1')
      case ('divQH')
        call add_field_grid_dp_2D( filename, ncid, 'divQH', long_name = 'Thickness divergence', units = 'm s^-1')
      case ('divQT')
        call add_field_grid_dp_2D( filename, ncid, 'divQT', long_name = 'Heat divergence', units = 'degC m s^-1')
      case ('divQS')
        call add_field_grid_dp_2D( filename, ncid, 'divQS', long_name = 'Salt divergence', units = 'PSU m s^-1')
      case ('diffT')
        call add_field_grid_dp_2D( filename, ncid, 'diffT', long_name = 'Heat diffusion', units = 'degC m s^-1')
      case ('diffS')
        call add_field_grid_dp_2D( filename, ncid, 'diffS', long_name = 'Salt diffusion', units = 'PSU m s^-1')
      case ('viscU')
        ! not implemented
      case ('viscV')
        ! not implemented
      case ('T_base')
        call add_field_grid_dp_2D( filename, ncid, 'T_base', long_name = 'Temperature at ice/ocean interface', units = 'deg C')
      case ('T_amb')
        call add_field_grid_dp_2D( filename, ncid, 'T_amb', long_name = 'Temperature at interface with ambient ocean', units = 'deg C')
      case ('u_star')
        call add_field_grid_dp_2D( filename, ncid, 'u_star', long_name = 'Friction velocity', units = 'm s^-1')
      case ('gamma_T')
        call add_field_grid_dp_2D( filename, ncid, 'gamma_T', long_name = 'Heat exchange coefficient', units = 'm s^-1')
      case ('divQU')
        ! not implemented
      case ('divQV')
        ! not implemented
      case ('HU_lad')
        ! not implemented
      case ('HV_lad')
        ! not implemented


    ! == Lateral mass balance ==
    ! ==========================

      ! Main LMB variables
      case ('LMB')
        call add_field_grid_dp_2D( filename, ncid, 'LMB', long_name = 'Lateral mass balance', units = 'm yr^-1')

    ! == Artificial mass balance ==
    ! =============================

      ! Main AMB variables
      case ('AMB')
        call add_field_grid_dp_2D( filename, ncid, 'AMB', long_name = 'Artificial mass balance', units = 'm yr^-1')

    ! == Glacial isostatic adjustment ==
    ! ==================================

      ! Main GIA variables
      case ('dHb_next')
        call add_field_grid_dp_2D( filename, ncid, 'dHb_next', long_name = 'Bedrock elevation difference from ELRA', units = 'm')

    ! == Tracer tracking ==
    ! =====================

      case ('age')
        call add_field_grid_dp_3D( filename, ncid, 'age', long_name = 'Age of ice', units = 'yr')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_main_regional_output_file_grid_field

end module grid_output_files
