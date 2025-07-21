module bed_roughness_main

  ! Contains all the routines for calculating the parameters
  ! that determine bed roughness.
  !
  ! NOTE: the variable ice%bed_roughness itself is computed
  ! from these parameters within each sliding law in their
  ! respective module

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use bed_roughness_model_types, only: type_bed_roughness_model
  use reference_geometry_types, only: type_reference_geometry
  use netcdf_io_main
  use reallocate_mod, only: reallocate_bounds
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D

  implicit none

  private

  public :: initialise_bed_roughness_model, remap_bed_roughness_model

contains

  subroutine initialise_bed_roughness_model( mesh, ice, bed_roughness, region_name)
    ! Initialise the bed roughness

    ! Input variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(in   ) :: ice
    type(type_bed_roughness_model),      intent(  out) :: bed_roughness
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_model'

    ! Add routine to path
    call init_routine( routine_name)

    allocate( bed_roughness%till_friction_angle( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( bed_roughness%alpha_sq           ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( bed_roughness%beta_sq            ( mesh%vi1:mesh%vi2), source = 0._dp)

    if (C%choice_bed_roughness == 'uniform') then
      ! Apply a uniform bed roughness

      call calc_bed_roughness_uniform( bed_roughness)

    elseif (C%choice_bed_roughness == 'parameterised') then
      ! Apply the chosen parameterisation of bed roughness

      call calc_bed_roughness_parameterised( mesh, ice, bed_roughness)

    elseif (C%choice_bed_roughness == 'read_from_file') then
      ! Initialise bed roughness from a NetCDF file

      call calc_bed_roughness_from_file( mesh, bed_roughness, region_name)

    else
      call crash('unknown choice_bed_roughness "' // trim( C%choice_bed_roughness) // '"')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_model

  subroutine remap_bed_roughness_model( mesh_old, mesh_new, ice, bed_roughness, region_name)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh_old
    type(type_mesh),                     intent(in   ) :: mesh_new
    type(type_ice_model),                intent(in   ) :: ice
    type(type_bed_roughness_model),      intent(inout) :: bed_roughness
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_bed_roughness_model'

    ! Add routine to path
    call init_routine( routine_name)

    if (C%do_bed_roughness_nudging) then
      ! Remap the existing bed roughness fields so they can continue to be nudged

      ! Bed roughness as described in different sliding laws
      call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, bed_roughness%till_friction_angle)
      call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, bed_roughness%alpha_sq)
      call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, bed_roughness%beta_sq)

      ! Main data fields
      call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, bed_roughness%generic_bed_roughness)

      ! Timestepping
      call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, bed_roughness%generic_bed_roughness_prev)
      call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, bed_roughness%generic_bed_roughness_next)

    else
      ! Simply re-initialise from whatever the user said

      deallocate( bed_roughness%till_friction_angle)
      deallocate( bed_roughness%alpha_sq           )
      deallocate( bed_roughness%beta_sq            )

      call initialise_bed_roughness_model( mesh_new, ice, bed_roughness, region_name)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_bed_roughness_model

! ===== Different bed roughness models ====
! =========================================

  subroutine calc_bed_roughness_uniform( bed_roughness)

    ! Input variables:
    type(type_bed_roughness_model), intent(inout) :: bed_roughness

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bed_roughness_uniform'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise field based on chosen sliding law
    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('no_sliding')
      ! No need to do anything
    case ('idealised')
      ! No need to do anything

    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by beta_sq
      bed_roughness%beta_sq = C%slid_Weertman_beta_sq_uniform

    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle
      bed_roughness%till_friction_angle = C%slid_Coulomb_phi_fric_uniform

    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle
      bed_roughness%till_friction_angle = C%slid_Coulomb_phi_fric_uniform

    case ('Tsai2015')
      ! Tsai2015 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part
      bed_roughness%alpha_sq = C%slid_Tsai2015_alpha_sq_uniform
      bed_roughness%beta_sq  = C%slid_Tsai2015_beta_sq_uniform

    case ('Schoof2005')
      ! Schoof2005 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part
      bed_roughness%alpha_sq = C%slid_Schoof2005_alpha_sq_uniform
      bed_roughness%beta_sq  = C%slid_Schoof2005_beta_sq_uniform

    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
      bed_roughness%till_friction_angle = C%slid_ZI_phi_fric_uniform

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_uniform

  subroutine calc_bed_roughness_parameterised( mesh, ice, bed_roughness)

    ! Input variables:
    type(type_mesh),                intent(in   ) :: mesh
    type(type_ice_model),           intent(in   ) :: ice
    type(type_bed_roughness_model), intent(inout) :: bed_roughness

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bed_roughness_parameterised'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_bed_roughness_parameterised)
    case default
      call crash('unknown choice_bed_roughness_parameterised "' // trim( C%choice_bed_roughness_parameterised) // '"')
    case ('Martin2011')
      call calc_bed_roughness_Martin2011( mesh, ice, bed_roughness)
    case ('MISMIPplus','MISMIP+')
      call calc_bed_roughness_MISMIPplus( bed_roughness)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_parameterised

  subroutine calc_bed_roughness_Martin2011( mesh, ice, bed_roughness)
    ! Calculate the till friction angle using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"

    ! Input variables:
    type(type_mesh),                intent(in   ) :: mesh
    type(type_ice_model),           intent(in   ) :: ice
    type(type_bed_roughness_model), intent(inout) :: bed_roughness

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bed_roughness_Martin2011'
    integer                        :: vi
    real(dp)                       :: weight_Hb

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. (&
      C%choice_sliding_law == 'Coulomb' .or. &
      C%choice_sliding_law == 'Budd' .or. &
      C%choice_sliding_law == 'Zoet-Iverson')) then
      call crash('only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"')
    end if

    do vi = mesh%vi1, mesh%vi2

      ! Compute till friction angle based on Martin et al. (2011) Eq. 10
      weight_Hb = min( 1._dp, max( 0._dp, &
        (ice%Hb( vi) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))

      bed_roughness%till_friction_angle( vi) = (1._dp - weight_Hb) * C%Martin2011till_phi_min + weight_Hb * C%Martin2011till_phi_max

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_Martin2011

  subroutine calc_bed_roughness_MISMIPplus( bed_roughness)

    ! Input variables:
    type(type_bed_roughness_model), intent(inout) :: bed_roughness

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bed_roughness_MISMIPplus'
    real(dp), parameter            :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    real(dp), parameter            :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_sliding_law)
    case default
      call crash('only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"')

    case ('Weertman')
      ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option

      bed_roughness%beta_sq  = MISMIPplus_beta_sq

    case ('Tsai2015')
      ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option

      bed_roughness%alpha_sq = MISMIPplus_alpha_sq
      bed_roughness%beta_sq  = MISMIPplus_beta_sq

    case ('Schoof2005')
      ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option

      bed_roughness%alpha_sq = MISMIPplus_alpha_sq
      bed_roughness%beta_sq  = MISMIPplus_beta_sq

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_MISMIPplus

  subroutine calc_bed_roughness_from_file( mesh, bed_roughness, region_name)
    ! Fill bed roughness with data from an external NetCDF file

    ! Input variables:
    type(type_mesh),                intent(in   ) :: mesh
    type(type_bed_roughness_model), intent(inout) :: bed_roughness
    character(len=3),               intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bed_roughness_from_file'
    character(len=1024)            :: filename_bed_roughness
    real(dp)                       :: timeframe_bed_roughness

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename and timeframe for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"')
    case ('NAM')
      filename_bed_roughness  = C%filename_bed_roughness_NAM
      timeframe_bed_roughness = C%timeframe_bed_roughness_NAM
    case ('EAS')
      filename_bed_roughness  = C%filename_bed_roughness_EAS
      timeframe_bed_roughness = C%timeframe_bed_roughness_EAS
    case ('GRL')
      filename_bed_roughness  = C%filename_bed_roughness_GRL
      timeframe_bed_roughness = C%timeframe_bed_roughness_GRL
    case ('ANT')
      filename_bed_roughness  = C%filename_bed_roughness_ANT
      timeframe_bed_roughness = C%timeframe_bed_roughness_ANT
    end select

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    if (index( filename_bed_roughness,'_LAST.nc') > 1) then
      call find_last_output_file( filename_bed_roughness)
      call find_last_timeframe(   filename_bed_roughness, timeframe_bed_roughness)
    end if

    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('no_sliding')
      ! No need to do anything
    case ('idealised')
      ! No need to do anything
    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by beta_sq

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'beta_sq', mesh, C%output_dir, bed_roughness%beta_sq)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'beta_sq', mesh, C%output_dir, bed_roughness%beta_sq, timeframe_bed_roughness)
      end if

    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, C%output_dir, bed_roughness%till_friction_angle)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, C%output_dir, bed_roughness%till_friction_angle, timeframe_bed_roughness)
      end if

    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, C%output_dir, bed_roughness%till_friction_angle)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, C%output_dir, bed_roughness%till_friction_angle, timeframe_bed_roughness)
      end if

    case ('Tsai2015')
      ! Tsai2015 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, C%output_dir, bed_roughness%alpha_sq)
        call read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, C%output_dir, bed_roughness%beta_sq)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, C%output_dir, bed_roughness%alpha_sq, timeframe_bed_roughness)
        call read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, C%output_dir, bed_roughness%beta_sq, timeframe_bed_roughness)
      end if

    case ('Schoof2005')
      ! Schoof2005 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, C%output_dir, bed_roughness%alpha_sq)
        call read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, C%output_dir, bed_roughness%beta_sq)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'alpha_sq', mesh, C%output_dir, bed_roughness%alpha_sq, timeframe_bed_roughness)
        call read_field_from_file_2D( filename_bed_roughness, 'beta_sq' , mesh, C%output_dir, bed_roughness%beta_sq, timeframe_bed_roughness)
      end if

    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, C%output_dir, bed_roughness%till_friction_angle)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, C%output_dir, bed_roughness%till_friction_angle, timeframe_bed_roughness)
      end if

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_from_file

end module bed_roughness_main
