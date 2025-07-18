module netcdf_find_timeframe

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_DOUBLE_PRECISION
  use mpi_basic, only: par, sync
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use netcdf_field_name_options
  use netcdf_check_dimensions
  use netcdf_read_var_primary
  use netcdf_basic_wrappers

  implicit none

  private

  public :: find_timeframe, find_last_output_file, find_last_timeframe

contains

  subroutine find_timeframe( filename, ncid, time, ti)
    ! Find the timeframe in the file that is closest to the desired time.
    ! if the file has no time dimension or variable, throw an error.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    real(dp),         intent(in   ) :: time
    integer,          intent(  out) :: ti

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'find_timeframe'
    integer                             :: nt, id_dim_time, id_var_time
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ierr
    integer                             :: tii
    real(dp)                            :: dt_min

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the file contains a valid time dimension and variable
    call check_time( filename, ncid)

    ! inquire size of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! inquire time variable ID
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! allocate memory
    allocate( time_from_file( nt))

    ! Read time from file
    call read_var_primary( filename, ncid, id_var_time, time_from_file)
    call MPI_BCAST( time_from_file(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Find timeframe closest to desired time
    if (time_from_file( 1) > time) then
      ! Desired time beyond lower limit
      call warning('desired timeframe at t = {dp_01} before start of file time for file "' &
        // trim( filename) // '"; reading data from t = {dp_02} instead!', &
        dp_01 = time, dp_02 = time_from_file( 1))
      ti = 1
    elseif (time_from_file( nt) < time) then
      ! Desired time beyond upper limit
      call warning('desired timeframe at t = {dp_01} after end of file time for file "' &
        // trim( filename) // '"; reading data from t = {dp_02} instead!', &
        dp_01 = time, dp_02 = time_from_file( nt))
      ti = nt
    else
      ! Desired time is within the file time
      dt_min = huge( 1._dp)
      do tii = 1, nt
        if (abs( time_from_file( tii) - time) < dt_min) then
          ti = tii
          dt_min = abs( time_from_file( tii) - time)
        end if
      end do
      if (dt_min > 0._dp) then
        call warning('desired timeframe at t = {dp_01} not present in file "' &
          // trim( filename) // '"; reading data from closest match at t = {dp_02} instead!', &
          dp_01 = time, dp_02 = time_from_file( ti))
      end if
    end if

    ! Clean up after yourself
    deallocate( time_from_file)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_timeframe

  subroutine find_last_output_file( filename)
    !< For a filename of the form 'dir_a/dir_b/dir_c/main_output_ANT_LAST.nc',
    !< find the last UFEMISM mesh output file in dir_a/dir_b/dir_c.

    ! In/output variables:
    character(len=*), intent(inout) :: filename

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'find_last_output_file'
    integer                             :: i, i_file
    character(len=5)                    :: i_str
    logical                             :: file_exists
    character(len=1024)                 :: filename_next

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    i = index( filename,'_LAST.nc')
    if (i <= 1) return

    filename = filename (1:i-1) // '_00001.nc'
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) call crash('couldnt find file ' // trim( filename))

    i_file = 1
    do while (.true.)

      i_file = i_file + 1

      if (i_file < 10) then
        write( i_str,'(A,I1)') '0000',i_file
      elseif (i_file < 100) then
        write( i_str,'(A,I2)') '000',i_file
      elseif (i_file < 1000) then
        write( i_str,'(A,I3)') '00',i_file
      elseif (i_file < 10000) then
        write( i_str,'(A,I4)') '0',i_file
      elseif (i_file < 100000) then
        write( i_str,'(A,I5)') i_file
      else
        call crash('whaa!')
      end if

      filename_next = filename(1:i-1) // '_' // i_str // '.nc'

      inquire( exist = file_exists, file = trim( filename_next))

      if (file_exists) then
        filename = filename_next
      else
        exit
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_last_output_file

  subroutine find_last_timeframe( filename, time_to_read)
    !< Find the last timeframe in the file

    ! In/output variables:
    character(len=*), intent(inout) :: filename
    real(dp),         intent(inout) :: time_to_read

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'find_last_timeframe'
    integer                             :: ncid, id_dim_time, id_var_time
    integer                             :: nt
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ierr

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    call open_existing_netcdf_file_for_reading( filename, ncid)

    ! Check if the file contains a valid time dimension and variable
    call check_time( filename, ncid)

    ! inquire size of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! inquire time variable ID
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! allocate memory
    allocate( time_from_file( nt))

    ! Read time from file
    call read_var_primary( filename, ncid, id_var_time, time_from_file)
    call MPI_BCAST( time_from_file(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Set the timeframe to the last one in this file
    time_to_read = time_from_file( nt)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_last_timeframe

end module netcdf_find_timeframe
