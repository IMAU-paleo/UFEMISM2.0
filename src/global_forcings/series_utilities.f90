MODULE series_utilities

  ! Utility functions used to update values from time series

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning, insert_val_into_string_int,insert_val_into_string_dp
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE global_forcing_types                                   , ONLY: type_global_forcing
  
  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

SUBROUTINE update_timeframes_from_record(time_axis, record, time_at_t0, time_at_t1, val_at_t0, val_at_t1, time)
    ! Update the timeframes so we can interpolate between two points in a time series record

    IMPLICIT NONE

    REAL(dp), dimension(:),         INTENT(IN)    :: time_axis
    REAL(dp), dimension(:),         INTENT(IN)    :: record
    REAL(dp),                       INTENT(INOUT) :: time_at_t0
    REAL(dp),                       INTENT(INOUT) :: time_at_t1
    REAL(dp),                       INTENT(INOUT) :: val_at_t0
    REAL(dp),                       INTENT(INOUT) :: val_at_t1
    REAL(dp),                       INTENT(IN   ) :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_timeframes_from_record'
    INTEGER                                            :: ti0, ti1, tii, ncid, nt
    CHARACTER(LEN=256)                                 :: str
    REAL(dp)                                           :: dt_min 

    ! Add routine to path
    CALL init_routine( routine_name)

    
    ! Update sea level
    ! Find timeframe closest to desired time
    nt = size(time_axis)
    if (time_axis( 1) > time) then
      ! Desired time beyond lower limit
      !call warning('desired timeframe at t = {dp_01} before start of sea level record time; reading data from t = {dp_02} instead!', &
      !  dp_01 = time, dp_02 = time_axis( 1))
      ti0 = 1
    elseif (time_axis( nt) < time) then
      ! Desired time beyond upper limit
      !call warning('desired timeframe at t = {dp_01} after end of sea level record time; reading data from t = {dp_02} instead!', &
      !  dp_01 = time, dp_02 = time_axis( nt))
      ti0 = nt
    else
      ! Desired time is within the file time
      dt_min = huge( 1._dp)
      do tii = 1, nt
        if (abs( time_axis( tii) - time) < dt_min) then
          ti0 = tii
          dt_min = abs( time_axis( tii) - time)
        end if
      end do
      if (dt_min > 0._dp) then
        !call warning('desired timeframe at t = {dp_01} not present in sea level record; reading data from closest match at t = {dp_02} instead!', &
        !  dp_01 = time, dp_02 = time_axis( ti0))
      end if
    end if
      
    
    time_at_t0 = time_axis(ti0)
    val_at_t0  = record(ti0)
      
    ! if the desired time is after t0, we take one record after for t1
    if (time >= time_at_t0) then
      if (ti0 == size(time_axis)) then
        IF (par%primary) WRITE(0,*) 'desired timeframe is at or beyond the last record. Using last available value for both timeframes...'
        time_at_t1 = time_axis(ti0)
        val_at_t1  = record(ti0)
      else
        time_at_t1 = time_axis(ti0+1)
        val_at_t1  = record(ti0+1)
      end if
    else
      ! otherwise we read one record before for t0, and that record is t1
      if (ti0 == 1) then
        IF (par%primary) WRITE(0,*) 'desired timeframe is at or before the first record. Using first available value for both timeframes...'
        time_at_t1 = time_axis(ti0)
        val_at_t1  = record(ti0)
      else
        time_at_t1  = time_axis(ti0)
        val_at_t1   = record(ti0)
        time_at_t0  = time_axis(ti0-1)
        val_at_t0   = record(ti0-1)
      end if
    end if

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_timeframes_from_record

  SUBROUTINE interpolate_value_from_forcing_record(time_at_t0, time_at_t1, val_at_t0, val_at_t1, time, interpolated_value)
  ! Interpolate between time frames

    IMPLICIT NONE

    REAL(dp),                          INTENT(IN   )   :: val_at_t0
    REAL(dp),                          INTENT(IN   )   :: val_at_t1
    REAL(dp),                          INTENT(IN   )   :: time_at_t0
    REAL(dp),                          INTENT(IN   )   :: time_at_t1
    REAL(dp),                          INTENT(IN   )   :: time
    REAL(dp),                          INTENT(  OUT)   :: interpolated_value
    

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'interpolate_value_from_forcing_record'
    REAL(dp)                                           :: wt0,wt1

    ! Add routine to path
    CALL init_routine( routine_name)


    if (time_at_t1 == time_at_t0) then
      wt0 = 0._dp
      wt1 = 1._dp
    else
      if (time > time_at_t1) then
        wt0 = 0._dp
      elseif (time < time_at_t0) then
        wt0 = 1._dp
      else
        wt0 = (time_at_t1 - time) / (time_at_t1 - time_at_t0)
      end if
      wt1 = 1._dp - wt0
    end if

    ! Interpolate the two timeframes - constant sea level over the entire region
    interpolated_value = wt0 * val_at_t0 + wt1 * val_at_t1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolate_value_from_forcing_record



END MODULE series_utilities