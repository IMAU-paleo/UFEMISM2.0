MODULE control_resources_and_error_messaging

  ! Keep track of which subroutine has control, so that error messages can actually tell
  ! you where they originated.
  !
  ! Also keep track of how much computation time each subroutine uses.

! ===== Preamble =====
! ====================

  use mpi_f08, only: MPI_WTIME
  USE precisions                                             , ONLY: dp
  use mpi_basic, only: par, sync

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  LOGICAL             :: do_colour_strings = .TRUE.
  LOGICAL             :: do_display_path   = .FALSE.

  CHARACTER(LEN=1024) :: routine_path
  integer             :: n_MPI_windows_used

  TYPE subroutine_resource_tracker
    CHARACTER(LEN = 2048) :: routine_path
    REAL(dp)              :: tstart, tcomp
    integer               :: n_MPI_windows_used_initial
    integer               :: n_MPI_windows_used_final
  END TYPE subroutine_resource_tracker

  TYPE( subroutine_resource_tracker), DIMENSION(:), ALLOCATABLE :: resource_tracker

CONTAINS

! ===== Control & resource tracking ======
! ========================================

  SUBROUTINE initialise_control_and_resource_tracker
    ! Initialise the control and resource tracker

    IMPLICIT NONE

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    INTEGER                                                            :: i,n

    ! Allocate space to track up to 2,000 subroutines, that should be enough for a while...
    n = 5000
    ALLOCATE( resource_tracker( n))

    ! Initialise values
    DO i = 1, n
      resource_tracker( i)%routine_path = 'subroutine_placeholder'
      resource_tracker( i)%tstart       = 0._dp
      resource_tracker( i)%tcomp        = 0._dp
    END DO

#endif

    ! Initialise the routine path
    routine_path = 'UFEMISM_program'

    ! Initialise the shared memory leak tracker
    n_MPI_windows_used = 0

  END SUBROUTINE initialise_control_and_resource_tracker

  SUBROUTINE init_routine( routine_name, do_track_resource_use)
    ! Initialise a subroutine; update the routine path

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: routine_name
    LOGICAL                                  , OPTIONAL, INTENT(IN)    :: do_track_resource_use

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    INTEGER                                                            :: len_path_tot, len_path_used, len_name
    INTEGER                                                            :: i
    LOGICAL                                                            :: do_track_resource_use_loc

    ! Check if routine_name has enough memory
    len_path_tot  = LEN(      routine_path)
    len_path_used = LEN_TRIM( routine_path)
    len_name      = LEN_TRIM( routine_name)

    IF (len_path_used + 1 + len_name > len_path_tot) THEN
      CALL crash('init_routine - ERROR: routine_path = "' // TRIM( routine_path) // '", no more space to append routine_name = "' // TRIM( routine_name) // '"!')
    END IF

    ! Append this routine to the routine path
    routine_path = TRIM( routine_path) // '/' // TRIM( routine_name)

    ! If so specified, print the current routine path to the terminal (useful for debugging)
    IF (do_display_path) THEN
      IF (par%primary) WRITE(0,'(A)') '   Initialising ' // TRIM( routine_path)
    END IF

    ! Check if resource use for this subroutine should be tracked
    ! (especially for the NetCDF routines we don't want to do this, as there are
    ! a great many of them and the resource tracker output file will become annoyingly big)

    IF (PRESENT( do_track_resource_use)) THEN
      do_track_resource_use_loc = do_track_resource_use
    ELSE
      do_track_resource_use_loc = .TRUE.
    END IF

    IF (do_track_resource_use_loc) THEN

      ! Initialise the computation time tracker
      CALL find_subroutine_in_resource_tracker( i)
      resource_tracker( i)%tstart = MPI_WTIME()

      ! Check how many MPI windows are in use when this subroutine is initialised
      resource_tracker( i)%n_MPI_windows_used_initial = n_MPI_windows_used

    ELSE

      routine_path = TRIM( routine_path) // '_NOTRACK'

    END IF

#endif

  END SUBROUTINE init_routine

  SUBROUTINE finalise_routine( routine_name, n_extra_MPI_windows_expected)
    ! Finalise; remove the current routine name from the routine path

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=256)                                 , INTENT(IN)    :: routine_name
    integer, optional, intent(in) :: n_extra_MPI_windows_expected

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    LOGICAL                                                            :: do_track_resource_use
    INTEGER                                                            :: len_path_tot, i
    REAL(dp)                                                           :: dt
    integer :: n_extra_MPI_windows_expected_

    ! If so specified, print the current routine path to the terminal (useful for debugging)
    IF (do_display_path) THEN
      IF (par%primary) WRITE(0,'(A)') '   Finalising   ' // TRIM( routine_path)
    END IF

    ! Check if resource use should be tracked for this subroutine
    i = INDEX( routine_path, '_NOTRACK')
    IF ( i == 0) THEN
      do_track_resource_use = .TRUE.
    ELSE
      do_track_resource_use = .FALSE.
    END IF

    IF (do_track_resource_use) THEN
      ! Resource use for this subroutine should be tracked

      ! Add computation time to the resource tracker
      CALL find_subroutine_in_resource_tracker( i)
      dt = MPI_WTIME() - resource_tracker( i)%tstart
      resource_tracker( i)%tcomp = resource_tracker( i)%tcomp + dt

      ! Check how many MPI windows are in use when this subroutine is finalised
      resource_tracker( i)%n_MPI_windows_used_final = n_MPI_windows_used

      if (present( n_extra_MPI_windows_expected)) then
        n_extra_MPI_windows_expected_ = n_extra_MPI_windows_expected
      else
        n_extra_MPI_windows_expected_ = 0
      end if

      ! If it is larger than at the start, mention this
      if (resource_tracker( i)%n_MPI_windows_used_final > &
        resource_tracker( i)%n_MPI_windows_used_initial + n_extra_MPI_windows_expected_) then

        ! Some exceptions where extra used MPI windows should not be interpreted as a memory leak
        if (index( routine_path, 'UFEMISM_program/initialise') == 0 .and. &
            index( routine_path, 'allocate_dist_shared') == 0 .and. &
            index( routine_path, 'UPSY_unit_tests_program') == 0 .and. &
            index( routine_path, 'run_all_unit_tests') == 0 .and. &
            index( routine_path, 'run_all_multinode_unit_tests') == 0 .and. &
            index( routine_path, 'run_all_component_tests') == 0) then

          call warning('shared memory was allocated but not freed, possibly memory leak ' // &
            '(n_init = {int_01}, n_final = {int_02})', &
            int_01 = resource_tracker( i)%n_MPI_windows_used_initial, &
            int_02 = resource_tracker( i)%n_MPI_windows_used_final)
        end if
      end if

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, routine_name)

      IF (i == 0) THEN
        CALL crash('finalise_routine - ERROR: routine_name = "' // TRIM( routine_name) // '" not found in routine_path = "' // TRIM( routine_path) // '"!')
      END IF

      ! Remove the current routine name from the routine path
      len_path_tot = LEN( routine_path)
      routine_path( i-1:len_path_tot) = ' '

    ELSE ! IF (do_track_resource_use) THEN
      ! Resource use for this subroutine should not be tracked

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, TRIM( routine_name) // '_NOTRACK')

      IF (i == 0) THEN
        CALL crash('finalise_routine - ERROR: routine_name = "' // TRIM( routine_name) // '" not found in routine_path = "' // TRIM( routine_path) // '"!')
      END IF

      ! Remove the current routine name from the routine path
      len_path_tot = LEN( routine_path)
      routine_path( i-1:len_path_tot) = ' '

    END IF ! IF (do_track_resource_use) THEN

#endif

  END SUBROUTINE finalise_routine

  SUBROUTINE find_subroutine_in_resource_tracker( i)
    ! Find the current subroutine in the resource tracker. If it's not there yet, add it.

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER                                            , INTENT(OUT)   :: i

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    INTEGER                                                            :: n

    n = SIZE( resource_tracker)

    DO i = 1, n
      IF     (resource_tracker( i)%routine_path == routine_path) THEN
        ! The current subroutine is listed at this position in the resource tracker
        RETURN
      ELSEIF (resource_tracker( i)%routine_path == 'subroutine_placeholder') THEN
        ! We've checked all listed subroutines and haven't found the current one; add it
        resource_tracker( i)%routine_path = routine_path
        RETURN
      END IF
    END DO

    ! If we've reached this point, then the resource tracker is overflowing
    CALL crash('Resource tracker overflows! Allocate more memory for it in control_resources_and_error_messaging/initialise_control_and_resource_tracker!')

#endif

  END SUBROUTINE find_subroutine_in_resource_tracker

  SUBROUTINE reset_resource_tracker
    ! Reset the computation times and maximum memory use for all subroutines in the resource tracker

    IMPLICIT NONE

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER                                                            :: i,n


    n = SIZE( resource_tracker)

    DO i = 1, n
      resource_tracker( i)%tstart      = 0._dp
      resource_tracker( i)%tcomp       = 0._dp
    END DO

#endif

  END SUBROUTINE reset_resource_tracker

! ===== Error messaging =====
! ===========================

  subroutine print_UFEMISM_start

    ! Local variables:
    character(len=1024) :: str1, str2
    integer             :: i

    str1 = ' '
    if (par%n_nodes == 1) then
      str1 = '===== Running UFEMISM on {int_01} cores ====='
      call insert_val_into_string_int( str1, '{int_01}', par%n)
    else
      str1 = '===== Running UFEMISM on {int_01} cores ({int_02} nodes) ====='
      call insert_val_into_string_int( str1, '{int_01}', par%n)
      call insert_val_into_string_int( str1, '{int_02}', par%n_nodes)
    end if

    str2 = ' '
    do i = 1, len_trim( str1)
      str2( i:i) = '='
    end do

    if (par%primary) then
      write(0,'(A)') ''
      write(0,'(A)') trim( colour_string( trim( str2),'green'))
      write(0,'(A)') trim( colour_string( trim( str1),'green'))
      write(0,'(A)') trim( colour_string( trim( str2),'green'))
    end if
    call sync

  end subroutine print_UFEMISM_start

  SUBROUTINE print_UFEMISM_end( tcomp)
    ! Print the UFEMISM end message to the screen

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: tcomp

    ! Local variables:
    CHARACTER(LEN=128)                                                 :: str1, str2
    INTEGER                                                            :: n,i
    INTEGER                                                            :: nr, ns, nm, nh, nd

    ! Calculate number of elapsed days, hours, minutes, and seconds since this run started
    ns = CEILING(tcomp)

    nr = MOD(ns, 60*60*24)
    nd = (ns - nr) / (60*60*24)
    ns = ns - (nd*60*60*24)

    nr = MOD(ns, 60*60)
    nh = (ns - nr) / (60*60)
    ns = ns - (nh*60*60)

    nr = MOD(ns, 60)
    nm = (ns - nr) / (60)
    ns = ns - (nm*60)

    ! Print to screen
    str1 = ' '
    str1 = '===== Finished running UFEMISM in {int_01} days, {int_02} hours, {int_03} minutes, and {int_04} seconds ====='
    CALL insert_val_into_string_int( str1, '{int_01}', nd)
    CALL insert_val_into_string_int( str1, '{int_02}', nh)
    CALL insert_val_into_string_int( str1, '{int_03}', nm)
    CALL insert_val_into_string_int( str1, '{int_04}', ns)

    n = LEN_TRIM( str1)
    str2 = ' '
    DO i = 1, n
      str2( i:i) = '='
    END DO

    IF (par%primary) WRITE(0,'(A)') ''
    IF (par%primary) WRITE(0,'(A)') TRIM( colour_string( str2,'green'))
    IF (par%primary) WRITE(0,'(A)') TRIM( colour_string( str1,'green'))
    IF (par%primary) WRITE(0,'(A)') TRIM( colour_string( str2,'green'))
    IF (par%primary) WRITE(0,'(A)') ''
    CALL sync

  END SUBROUTINE print_UFEMISM_end

  SUBROUTINE crash( err_msg, int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
                              dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    ! Crash the model, write the error message to the screen

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: err_msg
    INTEGER                                  , OPTIONAL, INTENT(IN)    :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    REAL(dp)                                 , OPTIONAL, INTENT(IN)    ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    CHARACTER(LEN=1024)                                                :: err_msg_loc
    INTEGER                                                            :: nc
    CHARACTER(LEN=9)                                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                                      :: process_str

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Set the process string (e.g. "05/16")
    IF     (par%n < 10) THEN
      nc = 1
    ELSEIF (par%n < 100) THEN
      nc = 2
    ELSEIF (par%n < 1000) THEN
      nc = 3
    ELSEIF (par%n < 10000) THEN
      nc = 4
    ELSE
      nc = 5
    END IF

    WRITE( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    ALLOCATE(CHARACTER(2*nc+1) :: process_str)
    WRITE( process_str,fmt) par%i, '/', par%n

    ! Insert numbers into string if needed
    IF (PRESENT( int_01)) CALL insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    IF (PRESENT( int_02)) CALL insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    IF (PRESENT( int_03)) CALL insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    IF (PRESENT( int_04)) CALL insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    IF (PRESENT( int_05)) CALL insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    IF (PRESENT( int_06)) CALL insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    IF (PRESENT( int_07)) CALL insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    IF (PRESENT( int_08)) CALL insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    IF (PRESENT( int_09)) CALL insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    IF (PRESENT( int_10)) CALL insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    IF (PRESENT( dp_01 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    IF (PRESENT( dp_02 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    IF (PRESENT( dp_03 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    IF (PRESENT( dp_04 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    IF (PRESENT( dp_05 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    IF (PRESENT( dp_06 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    IF (PRESENT( dp_07 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    IF (PRESENT( dp_08 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    IF (PRESENT( dp_09 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    IF (PRESENT( dp_10 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    WRITE(0,'(A,A,A,A,A,A)') colour_string(' ERROR: ' // TRIM( err_msg_loc),'red') // ' in ' // colour_string( TRIM(routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = primary)'

    ! Stop the program
    error stop
  END SUBROUTINE crash

  SUBROUTINE warning( err_msg, int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
                                dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    ! Write the warning message to the screen, but don't crash the model

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: err_msg
    INTEGER                                  , OPTIONAL, INTENT(IN)    :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    REAL(dp)                                 , OPTIONAL, INTENT(IN)    ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    CHARACTER(LEN=1024)                                                :: err_msg_loc
    INTEGER                                                            :: nc
    CHARACTER(LEN=9)                                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                                      :: process_str

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Set the process string (e.g. "05/16")
    IF     (par%n < 10) THEN
      nc = 1
    ELSEIF (par%n < 100) THEN
      nc = 2
    ELSEIF (par%n < 1000) THEN
      nc = 3
    ELSEIF (par%n < 10000) THEN
      nc = 4
    ELSE
      nc = 5
    END IF

    WRITE( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    ALLOCATE(CHARACTER(2*nc+1) :: process_str)
    WRITE( process_str,fmt) par%i, '/', par%n

    ! Insert numbers into string if needed
    IF (PRESENT( int_01)) CALL insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    IF (PRESENT( int_02)) CALL insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    IF (PRESENT( int_03)) CALL insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    IF (PRESENT( int_04)) CALL insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    IF (PRESENT( int_05)) CALL insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    IF (PRESENT( int_06)) CALL insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    IF (PRESENT( int_07)) CALL insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    IF (PRESENT( int_08)) CALL insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    IF (PRESENT( int_09)) CALL insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    IF (PRESENT( int_10)) CALL insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    IF (PRESENT( dp_01 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    IF (PRESENT( dp_02 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    IF (PRESENT( dp_03 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    IF (PRESENT( dp_04 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    IF (PRESENT( dp_05 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    IF (PRESENT( dp_06 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    IF (PRESENT( dp_07 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    IF (PRESENT( dp_08 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    IF (PRESENT( dp_09 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    IF (PRESENT( dp_10 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    WRITE(0,'(A,A,A,A,A,A)') colour_string(' WARNING: ' // TRIM( err_msg_loc),'yellow') // ' in ' // colour_string( TRIM(routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = primary)'

    ! Clean up after yourself
    DEALLOCATE( process_str)

  END SUBROUTINE warning

  SUBROUTINE happy( err_msg, int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
                                dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    ! Write a happy message to the screen

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: err_msg
    INTEGER                                  , OPTIONAL, INTENT(IN)    :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    REAL(dp)                                 , OPTIONAL, INTENT(IN)    ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    CHARACTER(LEN=1024)                                                :: err_msg_loc
    INTEGER                                                            :: nc
    CHARACTER(LEN=9)                                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                                      :: process_str

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Set the process string (e.g. "05/16")
    IF     (par%n < 10) THEN
      nc = 1
    ELSEIF (par%n < 100) THEN
      nc = 2
    ELSEIF (par%n < 1000) THEN
      nc = 3
    ELSEIF (par%n < 10000) THEN
      nc = 4
    ELSE
      nc = 5
    END IF

    WRITE( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    ALLOCATE(CHARACTER(2*nc+1) :: process_str)
    WRITE( process_str,fmt) par%i, '/', par%n

    ! Insert numbers into string if needed
    IF (PRESENT( int_01)) CALL insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    IF (PRESENT( int_02)) CALL insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    IF (PRESENT( int_03)) CALL insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    IF (PRESENT( int_04)) CALL insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    IF (PRESENT( int_05)) CALL insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    IF (PRESENT( int_06)) CALL insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    IF (PRESENT( int_07)) CALL insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    IF (PRESENT( int_08)) CALL insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    IF (PRESENT( int_09)) CALL insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    IF (PRESENT( int_10)) CALL insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    IF (PRESENT( dp_01 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    IF (PRESENT( dp_02 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    IF (PRESENT( dp_03 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    IF (PRESENT( dp_04 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    IF (PRESENT( dp_05 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    IF (PRESENT( dp_06 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    IF (PRESENT( dp_07 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    IF (PRESENT( dp_08 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    IF (PRESENT( dp_09 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    IF (PRESENT( dp_10 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    WRITE(0,'(A,A,A,A,A,A)') colour_string(' SUCCESS: ' // TRIM( err_msg_loc),'green') // ' in ' // colour_string( TRIM(routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = primary)'

    ! Clean up after yourself
    DEALLOCATE( process_str)

  END SUBROUTINE happy

  FUNCTION colour_string( str, col) RESULT( str_col)
    ! Add colour to a string for writing to the terminal

    IMPLICIT NONE

    ! Input variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*),                                    INTENT(IN)    :: str, col

    ! Result variables:
    CHARACTER(LEN=:), ALLOCATABLE                                      :: str_col

    ALLOCATE(CHARACTER(LEN(str)+9) :: str_col)       ! The +9 is just enough to store the color characters

    ! Safety: optionally don't do colours
    IF (.NOT. do_colour_strings) THEN
      str_col = str
      RETURN
    END IF

    ! The 91m gives red, 0m sets the default back
    ! Available colors: 90:gray, 91:red, 92:green, 93:yellow, 94:blue, 95:pink, 96:light blue
    IF     (col == 'gray') THEN
      str_col = achar(27)//'[90m'//str//achar(27)//'[0m'
    ELSEIF (col == 'red') THEN
      str_col = achar(27)//'[91m'//str//achar(27)//'[0m'
    ELSEIF (col == 'green') THEN
      str_col = achar(27)//'[92m'//str//achar(27)//'[0m'
    ELSEIF (col == 'yellow') THEN
      str_col = achar(27)//'[93m'//str//achar(27)//'[0m'
    ELSEIF (col == 'blue') THEN
      str_col = achar(27)//'[94m'//str//achar(27)//'[0m'
    ELSEIF (col == 'pink') THEN
      str_col = achar(27)//'[95m'//str//achar(27)//'[0m'
    ELSEIF (col == 'light blue') THEN
      str_col = achar(27)//'[96m'//str//achar(27)//'[0m'
    ELSE
      CALL crash('unknown colour "' // TRIM( col) // '"')
    END IF

  END FUNCTION colour_string

  SUBROUTINE insert_val_into_string_int( str, marker, val)
    ! Replace marker in str with val (where val is an integer)
    !
    ! Example: str    = 'Johnny has {int_01} apples.'
    !          marker = '{int_01}'
    !          val    = 5
    !
    ! This returns: str = 'Johnny has 5 apples'

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                   , INTENT(INOUT) :: str
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: marker
    INTEGER                                            , INTENT(IN)    :: val

    ! Local variables:
    INTEGER                                                            :: ci
    INTEGER                                                            :: nc
    CHARACTER(LEN=4)                                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                                      :: val_str
    INTEGER                                                            :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = INDEX( str, marker)

    ! Safety
    IF (ci == 0) CALL crash('insert_val_into_string_int: couldnt find marker "' // TRIM( marker) // '" in string "' // TRIM( str) // '"!')

    ! Write val to a string
    IF     (ABS( val) < 10) THEN
      nc = 1
    ELSEIF (ABS( val) < 100) THEN
      nc = 2
    ELSEIF (ABS( val) < 1000) THEN
      nc = 3
    ELSEIF (ABS( val) < 10000) THEN
      nc = 4
    ELSEIF (ABS( val) < 100000) THEN
      nc = 5
    ELSEIF (ABS( val) < 1000000) THEN
      nc = 6
    ELSEIF (ABS( val) < 10000000) THEN
      nc = 7
    ELSEIF (ABS( val) < 100000000) THEN
      nc = 8
    ELSE
      nc = 9
    END IF
    ! Add room for a minus sign if needed
    IF (val < 0) nc = nc + 1

    WRITE( fmt,'(A,I1,A)') '(I', nc, ')'
    ALLOCATE(CHARACTER(nc) :: val_str)
    WRITE( val_str,fmt) val

    ! Find total string length right now
    len_str    = LEN( str)
    len_marker = LEN( marker)

    ! Insert the integer string into the string
    str = str(1:ci-1) // val_str // str(ci+len_marker:len_str)

    ! Clean up after yourself
    DEALLOCATE( val_str)

  END SUBROUTINE insert_val_into_string_int

  SUBROUTINE insert_val_into_string_dp( str, marker, val)
    ! Replace marker in str with val (where val is a double-precision number)
    !
    ! Example: str    = 'Johnny weighs {dp_01} kg.'
    !          marker = '{dp_01}'
    !          val    = 57.098
    !
    ! This returns: str = 'Johnny weighs 57.098 kg'

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                   , INTENT(INOUT) :: str
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: marker
    REAL(dp)                                           , INTENT(IN)    :: val

    ! Local variables:
    INTEGER                                                            :: ci
    CHARACTER(LEN=20)                                                  :: val_str
    INTEGER                                                            :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = INDEX( str, marker)

    ! Safety
    IF (ci == 0) CALL crash('insert_val_into_string_dp: couldnt find marker "' // TRIM( marker) // '" in string "' // TRIM( str) // '"!')

    ! Write val to a string
    WRITE( val_str,'(E20.14)') val

    ! Find total string length right now
    len_str    = LEN( str)
    len_marker = LEN( marker)

    ! Insert the integer string into the string
    str = str(1:ci-1) // val_str // str(ci+len_marker:len_str)

  END SUBROUTINE insert_val_into_string_dp

  SUBROUTINE capitalise_string( str)
    ! Capitalise all letters in a string

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                   , INTENT(INOUT) :: str

    ! Local variables:
    INTEGER                                                            :: i, index_cap
    CHARACTER(26), PARAMETER                                           :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER                                           :: low = 'abcdefghijklmnopqrstuvwxyz'

    DO i = 1, LEN_TRIM( str)
      index_cap = INDEX( low, str( i:i))
      IF (index_cap > 0) str( i:i) = cap( index_cap:index_cap)
    END DO

  END SUBROUTINE capitalise_string

  SUBROUTINE remove_leading_spaces( str)
    ! Remove leading spaces from a character string

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                   , INTENT(INOUT) :: str

    ! Local variables:
    INTEGER                                                            :: lstr

    DO WHILE (str( 1:1) == ' ' .AND. LEN_TRIM( str) > 0)
      lstr = LEN_TRIM( str)
      str( 1:lstr-1) = str( 2:lstr)
      str( lstr:lstr) = ' '
    END DO

  END SUBROUTINE remove_leading_spaces

  SUBROUTINE str2int( str,int,stat)
    ! Convert a string containing an integer number to an actual integer

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*)                                   , INTENT(IN)    :: str
    INTEGER                                            , INTENT(OUT)   :: int
    INTEGER                                            , INTENT(OUT)   :: stat

    READ( str, *, IOSTAT = stat) int

  END SUBROUTINE str2int

  SUBROUTINE int2str( int,str)
    ! Convert integer to a character string (with leading zeros)

    IMPLICIT NONE

    ! In/output variables:
    INTEGER                                            , INTENT(IN)    :: int
    CHARACTER(LEN=*)                                   , INTENT(OUT)   :: str

    ! Local variables
    INTEGER                                                            :: n
    CHARACTER(LEN=4)                                                   :: fmt
    INTEGER                                                            :: i

    n = LEN( str)

    WRITE( fmt,'(A,I1,A)') '(I', n, ')'

    WRITE( str,fmt) int

    DO i = 1, n
      IF (str( i:i) == ' ') str( i:i) = '0'
    END DO

  END SUBROUTINE int2str

  function strrep( str, old, new) result( str_new)
    !< Replace all occurences in [str] of [old] with [new]

    ! In/output variables
    character(len=*), intent(inout) :: str
    character(len=*), intent(in   ) :: old, new
    character(len=len(str))         :: str_new

    ! Local variables:
    integer                         :: i, j, nit

    str_new = str
    nit = 0
    do while (index( str_new, old) > 0 .and. nit < len( str))
      nit = nit + 1
      i = index( str_new,old)
      j = i + len( old)
      str_new = str_new( 1:i-1) // trim(new) // str( j: len_trim( str_new))
    end do

  end function strrep

END MODULE control_resources_and_error_messaging
