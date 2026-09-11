module climate_snapshot_lapse

  ! A snapshot option for climate forcing, typically used for spinups.
  ! The climate forcing can optionally be lapse-rate corrected for elevation changes.
  !
  ! An input file from, e.g., RACMO, is needed containing:
  !
  ! - Pr         : [m] Total precipitation accumulated of a month
  ! - T2m        : [K] 2m air temperature averaged over a month
  ! - Hs         : [m] Surface elevation at which Pr and T2m were determined
  !
  ! Both Pr and T2m should have a dimension [month]
  !
  ! An example file is provided in the repository:
  ! external/data/atmosphere/Antarctica/RACMO_clim_1979_2021.nc

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning, crash
  use mpi_f08, only: MPI_WIN
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mesh_types, only: type_mesh
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use reference_geometry_types, only: type_reference_geometry
  use netcdf_basic
  use netcdf_io_main, only: read_field_from_file_2D_monthly, read_field_from_file_2D, read_field_from_file_0D
  use climate_model_basic, only: atype_climate_model

  implicit none

  private

  public :: type_climate_model_snapshot_lapse

  type, extends(atype_climate_model) :: type_climate_model_snapshot_lapse
    !< Variables and functions that are specific to the snapshot_lapse climate model

      ! Baseline climate and surface elevation
      real(dp), dimension(:,:), contiguous, pointer :: T2m_baseline    => null()   !< [K]                      Baseline monthly mean 2-m air temperature
      real(dp), dimension(:,:), contiguous, pointer :: Precip_baseline => null()   !< [m.w.e. month^-1]        Baseline monthly total precipitation
      real(dp), dimension(:  ), contiguous, pointer :: Hs_baseline     => null()   !< [m w.r.t. PD sea level]  Baseline surface elevation
      type(MPI_WIN) :: wT2m_baseline, wPrecip_baseline, wHs_baseline

      ! Insolation
      real(dp), dimension(:,:), contiguous, pointer :: Q_TOA           => null()
      type(MPI_WIN) :: wQ_TOA

      ! Wind
      real(dp), dimension(:,:), contiguous, pointer :: Wind_LR         => null()   ! [m s^-1]     10m wind velocity in x-direction
      real(dp), dimension(:,:), contiguous, pointer :: Wind_DU         => null()   ! [m s^-1]     10m wind velocity in y-direction
      type(MPI_WIN) :: wWind_LR, wWind_DU

      ! Region-specific info
      character(len=1024)                           :: filename_climate_snapshot
      logical                                       :: do_lapse_rate_corrections
      real(dp)                                      :: lapse_rate_temp
      logical                                       :: has_insolation

    contains

      procedure, public :: allocate_climate_model   => climate_model_snapshot_lapse_allocate
      procedure, public :: deallocate_climate_model => climate_model_snapshot_lapse_deallocate
      procedure, public :: initialise_climate_model => climate_model_snapshot_lapse_initialise
      procedure, public :: run_climate_model        => climate_model_snapshot_lapse_run
      procedure, public :: remap_climate_model      => climate_model_snapshot_lapse_remap

      procedure, public :: get_climate_model_name

      procedure, private :: initialise_insolation_forcing
      procedure, private :: apply_geometry_downscaling_corrections

  end type type_climate_model_snapshot_lapse

contains

  subroutine climate_model_snapshot_lapse_allocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot_lapse), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_lapse_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the snapshot_lapse climate model

    ! Allocate baseline climate
    call self%create_field( self%T2m_baseline, self%wT2m_baseline, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m_baseline', &
      long_name = 'Baseline monthly mean 2-m air temperature', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip_baseline, self%wPrecip_baseline, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip_baseline', &
      long_name = 'Baseline monthly total precipitation', &
      units     = 'm.w.e. month^-1', &
      remap_method = 'reallocate')

    ! Elevation-based temperature correction
    call self%create_field( self%Hs_baseline, self%wHs_baseline, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hs_baseline', &
      long_name = 'Baseline surface elevation', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%Q_TOA, self%wQ_TOA, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Q_TOA', &
      long_name = 'Monthly insolation', &
      units     = 'W m^-2', &
      remap_method = 'reallocate')

    call self%create_field( self%Wind_LR, self%wWind_LR, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Wind_LR', &
      long_name = 'Monthly 10m wind velocity in x-direction', &
      units     = 'm s^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%Wind_DU, self%wWind_DU, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Wind_DU', &
      long_name = 'Monthly 10m wind velocity in y-direction', &
      units     = 'm s^-1', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_lapse_allocate

  subroutine climate_model_snapshot_lapse_deallocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot_lapse), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_lapse_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to climate model snapshot_lapse

    ! Baseline climate and surface elevation
    nullify( self%T2m_baseline)
    nullify( self%Precip_baseline)
    nullify( self%Hs_baseline)

    ! Insolation
    nullify( self%Q_TOA)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_lapse_deallocate

  subroutine climate_model_snapshot_lapse_initialise( self, geom, refgeo_PD, refgeo_init, region_name)

    ! In/output variables:
    class(type_climate_model_snapshot_lapse), intent(inout) :: self
    class(atype_ice_geometry_model_data),     intent(in   ) :: geom
    type(type_reference_geometry),            intent(in   ) :: refgeo_PD
    type(type_reference_geometry),            intent(in   ) :: refgeo_init
    character(len=3),                         intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'climate_model_snapshot_lapse_initialise'
    character(len=1024)            :: filename_climate_snapshot

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to climate model snapshot_lapse

    self%has_insolation = .false. ! Initialise

    ! Determine which climate model to initialise for this region
    if     (region_name == 'NAM') then
      filename_climate_snapshot      = C%filename_climate_snapshot_NAM
      self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_NAM
      self%lapse_rate_temp           = C%lapse_rate_temp_NAM
      if (C%choice_SMB_model_NAM == 'IMAU-ITM' .or. C%choice_SMB_model_NAM == 'ITM_v2') then
         self%has_insolation = .true.
      end if
    elseif (region_name == 'EAS') then
      filename_climate_snapshot      = C%filename_climate_snapshot_EAS
      self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_EAS
      self%lapse_rate_temp           = C%lapse_rate_temp_EAS
      if (C%choice_SMB_model_EAS == 'IMAU-ITM' .or. C%choice_SMB_model_EAS == 'ITM_v2') then
         self%has_insolation = .true.
      end if
    elseif (region_name == 'GRL') then
      filename_climate_snapshot      = C%filename_climate_snapshot_GRL
      self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_GRL
      self%lapse_rate_temp           = C%lapse_rate_temp_GRL
      if (C%choice_SMB_model_GRL == 'IMAU-ITM' .or. C%choice_SMB_model_GRL == 'ITM_v2') then
         self%has_insolation = .true.
      end if
    elseif (region_name == 'ANT') then
      filename_climate_snapshot      = C%filename_climate_snapshot_ANT
      self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_ANT
      self%lapse_rate_temp           = C%lapse_rate_temp_ANT
      if (C%choice_SMB_model_ANT == 'IMAU-ITM' .or. C%choice_SMB_model_ANT == 'ITM_v2') then
         self%has_insolation = .true.
      end if
    else
      call crash('unknown region_name "' // region_name // '"')
    end if

    if (par%primary) then
      write(0,*) '   Reading climate baseline from file: ', &
        UPSY%stru%colour_string( trim( filename_climate_snapshot), 'light blue')
    end if

    ! Read the fixed baseline climate
    call read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , self%mesh, C%output_dir, self%Hs_baseline)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , self%mesh, C%output_dir, self%T2m_baseline)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', self%mesh, C%output_dir, self%Precip_baseline)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'uas'   , self%mesh, C%output_dir, self%Wind_LR)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'vas'   , self%mesh, C%output_dir, self%Wind_DU)

    call self%apply_geometry_downscaling_corrections( geom)

    ! Initialise insolation
    if (self%has_insolation) then
      call self%initialise_insolation_forcing()
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_lapse_initialise

  subroutine climate_model_snapshot_lapse_run( self, geom, time)

    ! In/output variables:
    class(type_climate_model_snapshot_lapse),  intent(inout) :: self
    class(atype_ice_geometry_model_data),      intent(in   ) :: geom
    real(dp),                                  intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_lapse_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to climate model snapshot_lapse

    ! TODO update insolation here in case time-variant

    ! Determine and apply downscaling corrections
    call self%apply_geometry_downscaling_corrections( geom)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_lapse_run

  subroutine climate_model_snapshot_lapse_remap( self, mesh_new)

    ! In/output variables:
    class(type_climate_model_snapshot_lapse), intent(inout) :: self
    type(type_mesh), target,          intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_lapse_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to climate model snapshot_lapse
    call crash('remapping not yet supported for snapshot_lapse climate forcing')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_lapse_remap

  subroutine initialise_insolation_forcing( self)

    ! In/output variables:
    class(type_climate_model_snapshot_lapse), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_insolation_forcing'
    real(dp)                    :: t0 = 0._dp

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary)   write(0,*) ' Initialising insolation data from ', trim(C%filename_insolation), '...'

    select case( C%choice_insolation_forcing)
    case ('none')
      call crash('Chosen climate or SMB model cannot be used with choice_insolation_forcing = "none"!')
    case ('static')
      call read_field_from_file_0D( C%filename_insolation, field_name_options_time, t0, time_to_read = C%static_insolation_time)
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, self%mesh, C%output_dir, self%Q_TOA, time_to_read = t0)
    case ('realistic')
      ! TODO to be implemented
    case default
      call crash('unknown choice insolation forcing "' // trim(C%choice_insolation_forcing) // '"')
    end select 

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_insolation_forcing

  subroutine apply_geometry_downscaling_corrections( self, geom)

    ! In/output variables:
    class(type_climate_model_snapshot_lapse), intent(inout) :: self
    class(atype_ice_geometry_model_data),     intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter :: routine_name = 'apply_geometry_downscaling_corrections'
    integer                     :: vi, mi
    real(dp)                    :: delta_T2m

    ! Add routine to call stack
    call init_routine( routine_name)

      do vi = self%mesh%vi1, self%mesh%vi2

        if (geom%mask_icefree_ocean( vi)) then
          ! Don't apply corrections over open ocean, just inherit baseline
          delta_T2m = 0._dp
        else
          delta_T2m = (self%Hs_baseline( vi) - geom%Hs( vi)) * self%lapse_rate_temp
        end if

        do mi = 1, 12
          self%T2m( vi, mi)    = self%T2m_baseline(    vi, mi) + delta_T2m
          self%Precip( vi, mi) = self%Precip_baseline( vi, mi) * C%precip_CC_correction_ANT**delta_T2m
        end do
      end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine apply_geometry_downscaling_corrections

  function get_climate_model_name( self) result( climate_model_name)
    class(type_climate_model_snapshot_lapse), intent(in) :: self
    character(len=:), allocatable :: climate_model_name
    climate_model_name = 'snapshot_lapse'
  end function get_climate_model_name

end module climate_snapshot_lapse
