module checksum_mod

  use precisions, only: dp
  use control_resources_and_error_messaging, only: routine_path, crash, &
    insert_val_into_string_int, insert_val_into_string_dp
  use mpi_basic, only: par
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUM, &
    MPI_MIN, MPI_MAX, MPI_COMM_WORLD

  implicit none

  private

  public :: create_checksum_logfile, checksum

  interface checksum
    procedure :: checksum_logical_1D
    procedure :: checksum_logical_2D
    procedure :: checksum_int_0D
    procedure :: checksum_int_1D
    procedure :: checksum_int_2D
    procedure :: checksum_dp_0D
    procedure :: checksum_dp_1D
    procedure :: checksum_dp_2D
  end interface checksum

  logical             :: do_write_checksum_log
  logical             :: checksum_logfile_exists
  character(len=1024) :: filename_checksum_logfile
  integer             :: unit_checksum_logfile

contains

  subroutine checksum_logical_1D( d, var_name, pai)
    logical,  dimension(:),            intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d, ierr

    if (.not. do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = count( d)
    else
      sum_d = count( d( pai%i1:pai%i2))
    end if

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    call log_checksum_int( sum_d, 0, 0, 0, var_name)

  end subroutine checksum_logical_1D

  subroutine checksum_logical_2D( d, var_name, pai)
    logical,  dimension(:,:),          intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d, ierr

    if (.not. do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = count( d)
    else
      sum_d = count( d( pai%i1:pai%i2,:))
    end if

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    call log_checksum_int( sum_d, 0, 0, 0, var_name)

  end subroutine checksum_logical_2D

  subroutine checksum_int_0D( d, var_name)
    ! For 0-D values, verify that all processes have the same value,
    ! and if so, only write the primary's value to the log

    ! In/output variables:
    integer,          intent(in) :: d
    character(len=*), intent(in) :: var_name

    ! Local variables:
    integer :: min_d, max_d
    integer :: ierr

    if (.not. do_write_checksum_log) return

    call MPI_ALLREDUCE( d, min_d, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( d, max_d, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (min_d /= max_d) call crash('checksum - scalar variable has different values on different processes')

    call log_checksum_int( d, 0, 0, 0, var_name)

  end subroutine checksum_int_0D

  subroutine checksum_int_1D( d, var_name, pai)
    integer,  dimension(:),            intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d, sum_abs_d, min_d, max_d, ierr

    if (.not. do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d     = sum   ( d)
      sum_abs_d = sum( abs( d))
      min_d     = minval( d)
      max_d     = maxval( d)
    else
      sum_d     = sum     ( d( pai%i1:pai%i2))
      sum_abs_d = sum( abs( d( pai%i1:pai%i2)))
      min_d     = minval  ( d( pai%i1:pai%i2))
      max_d     = maxval  ( d( pai%i1:pai%i2))
    end if

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d    , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_abs_d, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, min_d    , 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, max_d    , 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    call log_checksum_int( sum_d, sum_abs_d, min_d, max_d, var_name)

  end subroutine checksum_int_1D

  subroutine checksum_int_2D( d, var_name, pai)
    integer,  dimension(:,:),          intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d, sum_abs_d, min_d, max_d, ierr

    if (.not. do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d     = sum   ( d)
      sum_abs_d = sum( abs( d))
      min_d     = minval( d)
      max_d     = maxval( d)
    else
      sum_d     = sum     ( d( pai%i1:pai%i2,:))
      sum_abs_d = sum( abs( d( pai%i1:pai%i2,:)))
      min_d     = minval  ( d( pai%i1:pai%i2,:))
      max_d     = maxval  ( d( pai%i1:pai%i2,:))
    end if

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d    , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_abs_d, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, min_d    , 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, max_d    , 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    call log_checksum_int( sum_d, sum_abs_d, min_d, max_d, var_name)

  end subroutine checksum_int_2D

  subroutine checksum_dp_0D( d, var_name)
    ! For 0-D values, verify that all processes have the same value,
    ! and if so, only write the primary's value to the log

    ! In/output variables:
    real(dp),         intent(in   ) :: d
    character(len=*), intent(in   ) :: var_name

    ! Local variables:
    real(dp) :: min_d, max_d
    integer  :: ierr

    if (.not. do_write_checksum_log) return

    call MPI_ALLREDUCE( d, min_d, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( d, max_d, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (min_d /= max_d) call crash('checksum - scalar variable has different values on different processes')

    call log_checksum_dp( d, 0._dp, 0._dp, 0._dp, var_name)

  end subroutine checksum_dp_0D

  subroutine checksum_dp_1D( d, var_name, pai)
    real(dp), dimension(:),            intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    real(dp) :: sum_d, sum_abs_d, min_d, max_d
    integer  :: ierr

    if (.not. do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d     = sum   ( d)
      sum_abs_d = sum( abs( d))
      min_d     = minval( d)
      max_d     = maxval( d)
    else
      sum_d     = sum     ( d( pai%i1:pai%i2))
      sum_abs_d = sum( abs( d( pai%i1:pai%i2)))
      min_d     = minval  ( d( pai%i1:pai%i2))
      max_d     = maxval  ( d( pai%i1:pai%i2))
    end if

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_abs_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, min_d    , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, max_d    , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call log_checksum_dp( sum_d, sum_abs_d, min_d, max_d, var_name)

  end subroutine checksum_dp_1D

  subroutine checksum_dp_2D( d, var_name, pai)
    real(dp), dimension(:,:),          intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    real(dp) :: sum_d, sum_abs_d, min_d, max_d
    integer  :: ierr

    if (.not. do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d     = sum   ( d)
      sum_abs_d = sum( abs( d))
      min_d     = minval( d)
      max_d     = maxval( d)
    else
      sum_d     = sum     ( d( pai%i1:pai%i2,:))
      sum_abs_d = sum( abs( d( pai%i1:pai%i2,:)))
      min_d     = minval  ( d( pai%i1:pai%i2,:))
      max_d     = maxval  ( d( pai%i1:pai%i2,:))
    end if

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_abs_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, min_d    , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, max_d    , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call log_checksum_dp( sum_d, sum_abs_d, min_d, max_d, var_name)

  end subroutine checksum_dp_2D

  subroutine log_checksum_int( sum_d, sum_abs_d, min_d, max_d, var_name)

    ! In/output variables:
    integer,          intent(in) :: sum_d, sum_abs_d, min_d, max_d
    character(len=*), intent(in) :: var_name

    ! Local variables:
    character(len=1024) :: str

    str = trim( var_name) // ': sum = {int_01}, sum(abs) = {int_02}, min = {int_03}, max = {int_04} [' &
      // trim( routine_path) // ']'
    call insert_val_into_string_int( str, '{int_01}', sum_d)
    call insert_val_into_string_int( str, '{int_02}', sum_abs_d)
    call insert_val_into_string_int( str, '{int_03}', min_d)
    call insert_val_into_string_int( str, '{int_04}', max_d)

    call log_checksum( str)

  end subroutine log_checksum_int

  subroutine log_checksum_dp( sum_d, sum_abs_d, min_d, max_d, var_name)

    ! In/output variables:
    real(dp),         intent(in) :: sum_d, sum_abs_d, min_d, max_d
    character(len=*), intent(in) :: var_name

    ! Local variables:
    character(len=1024) :: str

    str = trim( var_name) // ': sum = {dp_01}, sum(abs) = {dp_02}, min = {dp_03}, max = {dp_04} [' &
      // trim( routine_path) // ']'
    call insert_val_into_string_dp( str, '{dp_01}', sum_d)
    call insert_val_into_string_dp( str, '{dp_02}', sum_abs_d)
    call insert_val_into_string_dp( str, '{dp_03}', min_d)
    call insert_val_into_string_dp( str, '{dp_04}', max_d)

    call log_checksum( str)

  end subroutine log_checksum_dp

  subroutine log_checksum( str)
    character(len=*), intent(in) :: str
    if (par%primary) then
      ! Write to terminal
      write(0,*) trim(str)
      ! Write to logfile
      write(unit_checksum_logfile,*) trim(str)
    end if
  end subroutine log_checksum

  subroutine create_checksum_logfile( output_dir)
    character(len=*), intent(in) :: output_dir
    filename_checksum_logfile = trim( output_dir) // '/checksum_logfile.txt'
    if (par%primary) then
      open( file = trim( filename_checksum_logfile), newunit = unit_checksum_logfile)
      write( unit_checksum_logfile,*) '-= Checksum logfile =-'
    end if
    checksum_logfile_exists = .true.
  end subroutine create_checksum_logfile

end module checksum_mod
