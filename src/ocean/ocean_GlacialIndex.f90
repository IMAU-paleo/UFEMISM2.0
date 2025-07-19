module ocean_GlacialIndex

  ! Glacial Index ocean model (i.e., T,S = T_warm,S_warm + GI * (T_cold,S_cold - T_warm,S_warm)

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use ocean_model_types                                      , only: type_ocean_model
  use netcdf_io_main
  use ocean_extrapolation                                    , only: extrapolate_ocean_forcing
  use series_utilities

  implicit none

contains

! ===== Main routines =====
! =========================

subroutine initialise_ocean_model_GlacialIndex( mesh, ice, ocean, region_name, start_time_of_run)
    ! Initialise the ocean model
    !
    ! Use a realistic ocean scheme

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean
    character(len=3),                       intent(in)    :: region_name
    real(dp),                               intent(in)    :: start_time_of_run

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_ocean_model_GlacialIndex'
    character(len=256)                                    :: filename_ocean_snapshot_warm, filename_ocean_snapshot_cold, filename_ocean_GI

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary)  write(*,"(A)") '     Initialising transient ocean model "' // &
            colour_string( trim( C%choice_ocean_model_transient),'light blue') // '"...'
    ! We need the snapshot and the dT to apply to it
    select case (region_name)
        case ('NAM')
        filename_ocean_snapshot_warm = C%filename_ocean_warm_snapshot_NAM
        filename_ocean_snapshot_cold = C%filename_ocean_cold_snapshot_NAM
        filename_ocean_GI            = C%filename_ocean_GI_NAM
        case ('EAS')
        filename_ocean_snapshot_warm = C%filename_ocean_warm_snapshot_EAS
        filename_ocean_snapshot_cold = C%filename_ocean_cold_snapshot_EAS
        filename_ocean_GI            = C%filename_ocean_GI_EAS
        case ('GRL')
        filename_ocean_snapshot_warm = C%filename_ocean_warm_snapshot_GRL
        filename_ocean_snapshot_cold = C%filename_ocean_cold_snapshot_GRL
        filename_ocean_GI            = C%filename_ocean_GI_GRL
        case ('ANT')
        filename_ocean_snapshot_warm = C%filename_ocean_warm_snapshot_ANT
        filename_ocean_snapshot_cold = C%filename_ocean_cold_snapshot_ANT
        filename_ocean_GI            = C%filename_ocean_GI_ANT
        case default
        call crash('unknown region_name "' // region_name // '"')
    end select

    ! Allocating timeframe variables; the series itself is allocated in the read function below
    allocate(ocean%GI%GI_t0)
    allocate(ocean%GI%GI_t1)
    allocate(ocean%GI%GI_at_t0)
    allocate(ocean%GI%GI_at_t1)
    allocate( ocean%GI%T0_warm( mesh%vi1:mesh%vi2,C%nz_ocean))
    allocate( ocean%GI%S0_warm( mesh%vi1:mesh%vi2,C%nz_ocean))
    allocate( ocean%GI%T0_cold( mesh%vi1:mesh%vi2,C%nz_ocean))
    allocate( ocean%GI%S0_cold( mesh%vi1:mesh%vi2,C%nz_ocean))
    ocean%GI%T0_warm = 0._dp
    ocean%GI%S0_warm = 0._dp
    ocean%GI%T0_cold = 0._dp
    ocean%GI%S0_cold = 0._dp

    ! Fill in  main variables
    call read_field_from_file_3D_ocean( filename_ocean_snapshot_warm, field_name_options_T_ocean,  mesh, C%output_dir, C%z_ocean, ocean%GI%T0_warm)
    call read_field_from_file_3D_ocean( filename_ocean_snapshot_warm, field_name_options_S_ocean,  mesh, C%output_dir, C%z_ocean, ocean%GI%S0_warm)
    call read_field_from_file_3D_ocean( filename_ocean_snapshot_cold, field_name_options_T_ocean,  mesh, C%output_dir, C%z_ocean, ocean%GI%T0_cold)
    call read_field_from_file_3D_ocean( filename_ocean_snapshot_cold, field_name_options_S_ocean,  mesh, C%output_dir, C%z_ocean, ocean%GI%S0_cold)

    call read_field_from_series_file(   filename_ocean_GI,       field_name_options_GI, ocean%GI%GI_series, ocean%GI%GI_series_time)
    call update_timeframes_from_record(ocean%GI%GI_series_time, ocean%GI%GI_series, ocean%GI%GI_t0, ocean%GI%GI_t1, ocean%GI%GI_at_t0, ocean%GI%GI_at_t1, start_time_of_run)

    ! Apply extrapolation method if required
    select case (C%choice_ocean_extrapolation_method)
        case('initialisation')
        call extrapolate_ocean_forcing( mesh, ice, ocean%GI%T0_warm)
        call extrapolate_ocean_forcing( mesh, ice, ocean%GI%S0_warm)
        call extrapolate_ocean_forcing( mesh, ice, ocean%GI%T0_cold)
        call extrapolate_ocean_forcing( mesh, ice, ocean%GI%S0_cold)
        case default
        call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
    end select

    call finalise_routine(routine_name)

end subroutine initialise_ocean_model_GlacialIndex

subroutine run_ocean_model_GlacialIndex(mesh, ocean, time)
  ! Runs the Glacial Index transient ocean model
    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ocean_model),                 intent(inout) :: ocean
    real(dp),                               intent(in)    :: time

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_GlacialIndex'
    REAL(dp)                                           :: GI_at_time
    INTEGER                                            :: vi, z

    call init_routine( routine_name)

    IF (time < ocean%GI%GI_t0 .OR. time > ocean%GI%GI_t1) THEN
          !IF (par%primary)  WRITE(0,*) '   Model time is out of the current GI timeframes. Updating timeframes...'
          call update_timeframes_from_record(ocean%GI%GI_series_time, ocean%GI%GI_series, ocean%GI%GI_t0, ocean%GI%GI_t1, ocean%GI%GI_at_t0, ocean%GI%GI_at_t1, time)
        END IF

        call interpolate_value_from_forcing_record(ocean%GI%GI_t0, ocean%GI%GI_t1, ocean%GI%GI_at_t0, ocean%GI%GI_at_t1, time, GI_at_time)

        do vi = mesh%vi1, mesh%vi2
          do z = 1, C%nz_ocean
            ocean%T(vi, z) = ocean%GI%T0_warm(vi, z) + GI_at_time * (ocean%GI%T0_cold(vi, z) - ocean%GI%T0_warm(vi, z))
            ocean%S(vi, z) = ocean%GI%S0_warm(vi, z) + GI_at_time * (ocean%GI%S0_cold(vi, z) - ocean%GI%S0_warm(vi, z))
          end do
        end do


    call finalise_routine(routine_name)

end subroutine run_ocean_model_GlacialIndex


end module ocean_GlacialIndex