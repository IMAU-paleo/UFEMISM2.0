MODULE LMB_GlacialIndex

  ! The Glacial Index-dependent LMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE LMB_model_types                                        , ONLY: type_LMB_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use netcdf_io_main
  use series_utilities

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

subroutine run_LMB_model_GlacialIndex(mesh, ice, LMB, time)
! Calculate the lateral mass balance based on a Glacial Index

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_LMB_model),                   INTENT(INOUT) :: LMB
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_LMB_model_GlacialIndex'
    REAL(dp)                                              :: GI_at_time, LMB_at_time
    INTEGER                                               :: vi

  ! Add routine to path
  CALL init_routine( routine_name)

  IF (time < LMB%GI%GI_t0 .OR. time > LMB%GI%GI_t1) THEN
    !IF (par%primary)  WRITE(0,*) '   Model time is out of the current GI timeframes. Updating timeframes...'
    call update_timeframes_from_record(LMB%GI%GI_series_time, LMB%GI%GI_series, LMB%GI%GI_t0, LMB%GI%GI_t1, LMB%GI%GI_at_t0, LMB%GI%GI_at_t1, time)
  END IF

  call interpolate_value_from_forcing_record(LMB%GI%GI_t0, LMB%GI%GI_t1, LMB%GI%GI_at_t0, LMB%GI%GI_at_t1, time, GI_at_time)
  LMB_at_time = LMB%GI%LMB_warm + GI_at_time * (LMB%GI%LMB_cold - LMB%GI%LMB_warm)

  do vi = mesh%vi1, mesh%vi2
    if (ice%mask_cf_fl( vi) .OR. ice%mask_cf_gr( vi)) then
      LMB%LMB(vi) = LMB_at_time
    else
      LMB%LMB(vi) = 0._dp
    end if
  end do

  ! Finalise routine path
  CALL finalise_routine( routine_name)


end subroutine run_LMB_model_GlacialIndex

subroutine initialise_LMB_model_GlacialIndex(mesh, LMB, region_name, start_time_of_run)
! Initialise the lateral mass balance based on a Glacial Index

  IMPLICIT NONE

  ! In- and output variables
  type(type_mesh),                        intent(in)    :: mesh
  type(type_LMB_model),                   intent(inout) :: LMB
  character(len=3),                       intent(in)    :: region_name
  real(dp),                               intent(in)    :: start_time_of_run

  ! Local variables:
  character(len=256), parameter                         :: routine_name = 'initialise_LMB_model_GlacialIndex'
  character(len=256)                                    :: filename_LMB_GI
  real(dp)                                              :: LMB_cold, LMB_warm


  ! Add routine to path
  CALL init_routine( routine_name)
    select case (region_name)
        case ('NAM')
          filename_LMB_GI              = C%filename_LMB_GI_NAM
          LMB_warm                     = C%warm_LMB_NAM
          LMB_cold                     = C%cold_LMB_NAM
        case ('EAS')
          filename_LMB_GI              = C%filename_LMB_GI_EAS
          LMB_warm                     = C%warm_LMB_EAS
          LMB_cold                     = C%cold_LMB_EAS
        case ('GRL')
          filename_LMB_GI              = C%filename_LMB_GI_GRL
          LMB_warm                     = C%warm_LMB_GRL
          LMB_cold                     = C%cold_LMB_GRL
        case ('ANT')
          filename_LMB_GI              = C%filename_LMB_GI_ANT
          LMB_warm                     = C%warm_LMB_ANT
          LMB_cold                     = C%cold_LMB_ANT
        case default
          call crash('unknown region_name "' // region_name // '"')
    end select

    ! Allocating timeframe variables; the series itself is allocated in the read function below
    allocate(LMB%GI%GI_t0)
    allocate(LMB%GI%GI_t1)
    allocate(LMB%GI%GI_at_t0)
    allocate(LMB%GI%GI_at_t1)
    allocate(LMB%GI%LMB_warm)
    allocate(LMB%GI%LMB_cold)
    LMB%GI%LMB_warm = LMB_warm
    LMB%GI%LMB_cold = LMB_cold

    ! Fill in  main variables
    call read_field_from_series_file(   filename_LMB_GI,       field_name_options_GI, LMB%GI%GI_series, LMB%GI%GI_series_time)
    call update_timeframes_from_record(LMB%GI%GI_series_time, LMB%GI%GI_series, LMB%GI%GI_t0, LMB%GI%GI_t1, LMB%GI%GI_at_t0, LMB%GI%GI_at_t1, start_time_of_run)

  ! Finalise routine path
  CALL finalise_routine( routine_name)

end subroutine initialise_LMB_model_GlacialIndex




END MODULE LMB_GlacialIndex