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
  use reference_geometry_types, only: type_reference_geometry
  use analytical_solutions, only: Schoof2006_icestream
  use netcdf_io_main

  implicit none

  private

  public :: run_bed_roughness_model, initialise_bed_roughness

contains

  subroutine run_bed_roughness_model( mesh, ice, refgeo, time)
    ! Run the chosen bed roughness model

    ! Input variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_reference_geometry),       intent(in   ) :: refgeo
    real(dp),                            intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_bed_roughness_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate bed roughness using the chosen model
    ! ==============================================

    select case (C%choice_bed_roughness)
    case default
      call crash('unknown choice_bed_roughness "' // trim( C%choice_bed_roughness) // '"')

    case ('uniform')
      ! Apply a uniform bed roughness

      ! No need to do anything

    case ('parameterised')
      ! Apply the chosen parameterisation of bed roughness

      call calc_bed_roughness_parameterised( mesh, ice)

    case ('read_from_file')
      ! Initialise bed roughness from a NetCDF file

      ! No need to do anything

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_bed_roughness_model

  subroutine initialise_bed_roughness( mesh, ice, region_name)
    ! Initialise the bed roughness

    ! Input variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialisation
    ! ==============

    if (C%choice_bed_roughness == 'uniform') then
      ! Apply a uniform bed roughness

      call initialise_bed_roughness_uniform( mesh, ice)

    elseif (C%choice_bed_roughness == 'parameterised') then
      ! Apply the chosen parameterisation of bed roughness

      call initialise_bed_roughness_parameterised( mesh, ice)

    elseif (C%choice_bed_roughness == 'read_from_file') then
      ! Initialise bed roughness from a NetCDF file

      call initialise_bed_roughness_from_file( mesh, ice, region_name)

    else
      call crash('unknown choice_bed_roughness "' // trim( C%choice_bed_roughness) // '"')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness

  subroutine remap_bed_roughness( mesh_old, mesh_new, ice)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh_old
    type(type_mesh),                     intent(in   ) :: mesh_new
    type(type_ice_model),                intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_bed_roughness'

    ! Add routine to path
    call init_routine( routine_name)

    ! DENK DROM
    call crash('fixme')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_bed_roughness

! ===== Different bed roughness models ====
! =========================================

  ! == Uniform bed roughness
  subroutine initialise_bed_roughness_uniform( mesh, ice)
    ! Initialise bed roughness
    !
    ! Use a uniform value over the whole domain

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_uniform'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise bed roughness
    ! ========================

    ! Initialise field based on chosen sliding law
    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('no_sliding')
      ! No need to do anything
    case ('idealised')
      ! No need to do anything

    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by slid_beta_sq
      ice%slid_beta_sq = C%slid_Weertman_beta_sq_uniform

    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle
      ice%till_friction_angle = C%slid_Coulomb_phi_fric_uniform

    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle
      ice%till_friction_angle = C%slid_Coulomb_phi_fric_uniform

    case ('Tsai2015')
      ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
      ice%slid_alpha_sq = C%slid_Tsai2015_alpha_sq_uniform
      ice%slid_beta_sq  = C%slid_Tsai2015_beta_sq_uniform

    case ('Schoof2005')
      ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part
      ice%slid_alpha_sq = C%slid_Schoof2005_alpha_sq_uniform
      ice%slid_beta_sq  = C%slid_Schoof2005_beta_sq_uniform

    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle
      ice%till_friction_angle = C%slid_ZI_phi_fric_uniform

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_uniform

  ! == parameterised bed roughness
  subroutine calc_bed_roughness_parameterised( mesh, ice)
    ! Compute bed roughness
    ! Use a simple parameterisation to calculate bed roughness

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bed_roughness_parameterised'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_bed_roughness_parameterised)
    case default
      call crash('unknown choice_bed_roughness_parameterised "' // trim( C%choice_bed_roughness_parameterised) // '"')
    case ('Martin2011')
      call calc_bed_roughness_Martin2011( mesh, ice)
    case ('SSA_icestream')
      ! No need to do anything
    case ('MISMIPplus')
      ! No need to do anything
    case ('MISMIP+')
      ! No need to do anything
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_parameterised

  subroutine initialise_bed_roughness_parameterised( mesh, ice)
    ! Initialise the bed roughness
    ! Use a simple parameterisation to calculate bed roughness

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_parameterised'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_bed_roughness_parameterised)
    case default
      call crash('unknown choice_bed_roughness_parameterised "' // trim( C%choice_bed_roughness_parameterised) // '"')
    case ('Martin2011')
      call initialise_bed_roughness_Martin2011( mesh, ice)
    case ('SSA_icestream')
      call initialise_bed_roughness_SSA_icestream( mesh, ice)
    case ('MISMIPplus')
      call initialise_bed_roughness_MISMIPplus( mesh, ice)
    case ('MISMIP+')
      call initialise_bed_roughness_MISMIPplus( mesh, ice)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_parameterised

  ! The Martin et al. (2011) till parameterisation
  subroutine calc_bed_roughness_Martin2011( mesh, ice)
    ! Calculate the till friction angle using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

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

      ice%till_friction_angle( vi) = (1._dp - weight_Hb) * C%Martin2011till_phi_min + weight_Hb * C%Martin2011till_phi_max

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_Martin2011

  subroutine initialise_bed_roughness_Martin2011( mesh, ice)
    ! Calculate the till friction angle using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_Martin2011'
    integer                        :: vi
    real(dp)                       :: weight_Hb

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Budd' .OR. C%choice_sliding_law == 'Zoet-Iverson')) then
      call crash('only applicable when choice_sliding_law = "Coulomb", "Budd", or "Zoet-Iverson"')
    end if

    do vi = mesh%vi1, mesh%vi2

      ! Compute till friction angle based on Martin et al. (2011) Eq. 10
      weight_Hb = min( 1._dp, max( 0._dp, &
        (ice%Hb( vi) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))

      ice%till_friction_angle( vi) = (1._dp - weight_Hb) * C%Martin2011till_phi_min + weight_Hb * C%Martin2011till_phi_max

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_Martin2011

  ! Idealised cases
  subroutine initialise_bed_roughness_SSA_icestream( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)

    ! Input variables:
    type(type_mesh),     intent(in   ) :: mesh
    type(type_ice_model),intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_SSA_icestream'
    integer                        :: vi
    real(dp)                       :: y, u

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      y = mesh%V( vi,2)
      call Schoof2006_icestream( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, C%refgeo_idealised_SSA_icestream_Hi, &
        C%refgeo_idealised_SSA_icestream_dhdx, C%refgeo_idealised_SSA_icestream_L, C%refgeo_idealised_SSA_icestream_m, &
        y, u, ice%till_yield_stress( vi))
    end DO

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_SSA_icestream

  subroutine initialise_bed_roughness_MISMIPplus( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: MISMIP+ (i.e. just a uniform value)

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_MISMIPplus'
    real(dp), parameter            :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    real(dp), parameter            :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_sliding_law)
    case default
      call crash('only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"')

    case ('Weertman')
      ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option

      ice%slid_beta_sq  = MISMIPplus_beta_sq

    case ('Tsai2015')
      ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option

      ice%slid_alpha_sq = MISMIPplus_alpha_sq
      ice%slid_beta_sq  = MISMIPplus_beta_sq

    case ('Schoof2005')
      ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option

      ice%slid_alpha_sq = MISMIPplus_alpha_sq
      ice%slid_beta_sq  = MISMIPplus_beta_sq

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_MISMIPplus

  ! == Bed roughness from an external file
  subroutine initialise_bed_roughness_from_file( mesh, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice
    character(len=3),     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_from_file'
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

    select case (C%choice_sliding_law)
    case default
      call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"')
    case ('no_sliding')
      ! No need to do anything
    case ('idealised')
      ! No need to do anything
    case ('Weertman')
      ! Weertman sliding law; bed roughness is described by slid_beta_sq

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq', mesh, ice%slid_beta_sq)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq', mesh, ice%slid_beta_sq, timeframe_bed_roughness)
      end if

    case ('Coulomb')
      ! Coulomb sliding law; bed roughness is described by till_friction_angle

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle, timeframe_bed_roughness)
      end if

    case ('Budd')
      ! Budd-type sliding law; bed roughness is described by till_friction_angle

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle, timeframe_bed_roughness)
      end if

    case ('Tsai2015')
      ! Tsai2015 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq)
        call read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq, timeframe_bed_roughness)
        call read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq, timeframe_bed_roughness)
      end if

    case ('Schoof2005')
      ! Schoof2005 sliding law; bed roughness is described by slid_alpha_sq for the Coulomb part, and slid_beta_sq for the Weertman part

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq)
        call read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'slid_alpha_sq', mesh, ice%slid_alpha_sq, timeframe_bed_roughness)
        call read_field_from_file_2D( filename_bed_roughness, 'slid_beta_sq' , mesh, ice%slid_beta_sq, timeframe_bed_roughness)
      end if

    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law; bed roughness is described by till_friction_angle

      if (timeframe_bed_roughness == 1E9_dp) then
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle)
      else
        call read_field_from_file_2D( filename_bed_roughness, 'till_friction_angle||phi_fric', mesh, ice%till_friction_angle, timeframe_bed_roughness)
      end if

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_from_file

end module bed_roughness_main
