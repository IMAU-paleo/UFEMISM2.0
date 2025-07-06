module checksum_mod

  use precisions, only: dp
  use model_configuration, only: C
  use control_resources_and_error_messaging, only: routine_path, crash, &
    insert_val_into_string_int, insert_val_into_string_dp
  use mpi_basic, only: par
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUM, &
    MPI_MIN, MPI_MAX, MPI_COMM_WORLD

  implicit none

  private

  public :: checksum

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

  logical             :: checksum_logfile_exists
  character(len=1024) :: filename_checksum_logfile
  integer             :: unit_checksum_logfile

contains

  subroutine checksum_logical_1D( d, var_name, pai)
    logical,  dimension(:),            intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d

    if (.not. C%do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = count( d)
    else
      sum_d = count( d( pai%i1_nih:pai%i2_nih))
    end if

    call log_checksum_int( sum_d, var_name)

  end subroutine checksum_logical_1D

  subroutine checksum_logical_2D( d, var_name, pai)
    logical,  dimension(:,:),          intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d

    if (.not. C%do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = count( d)
    else
      sum_d = count( d( pai%i1_nih:pai%i2_nih,:))
    end if

    call log_checksum_int( sum_d, var_name)

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

    if (.not. C%do_write_checksum_log) return

    call MPI_ALLREDUCE( d, min_d, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( d, max_d, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (min_d /= max_d) call crash('checksum - scalar variable has different values on different processes')

    call log_checksum_int( d, var_name)

  end subroutine checksum_int_0D

  subroutine checksum_int_1D( d, var_name, pai)
    integer,  dimension(:),            intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d

    if (.not. C%do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = sum( d)
    else
      sum_d = sum( d( pai%i1_nih:pai%i2_nih))
    end if

    call log_checksum_int( sum_d, var_name)

  end subroutine checksum_int_1D

  subroutine checksum_int_2D( d, var_name, pai)
    integer,  dimension(:,:),          intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    integer :: sum_d

    if (.not. C%do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = sum( d)
    else
      sum_d = sum( d( pai%i1_nih:pai%i2_nih,:))
    end if

    call log_checksum_int( sum_d, var_name)

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

    if (.not. C%do_write_checksum_log) return

    call MPI_ALLREDUCE( d, min_d, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( d, max_d, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (min_d /= max_d) call crash('checksum - scalar variable has different values on different processes')

    call log_checksum_dp( d, var_name)

  end subroutine checksum_dp_0D

  subroutine checksum_dp_1D( d, var_name, pai)
    real(dp), dimension(:),            intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    real(dp) :: sum_d

    if (.not. C%do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = sum( d)
    else
      sum_d = sum( d( pai%i1_nih:pai%i2_nih))
    end if

    call log_checksum_dp( sum_d, var_name)

  end subroutine checksum_dp_1D

  subroutine checksum_dp_2D( d, var_name, pai)
    real(dp), dimension(:,:),          intent(in) :: d
    character(len=*),                  intent(in) :: var_name
    type(type_par_arr_info), optional, intent(in) :: pai

    real(dp) :: sum_d

    if (.not. C%do_write_checksum_log) return

    if (.not. present( pai)) then
      sum_d = sum( d)
    else
      sum_d = sum( d( pai%i1_nih:pai%i2_nih,:))
    end if

    call log_checksum_dp( sum_d, var_name)

  end subroutine checksum_dp_2D

  subroutine log_checksum_int( sum_d, var_name)

    ! In/output variables:
    integer,          intent(in) :: sum_d
    character(len=*), intent(in) :: var_name

    ! Local variables:
    character(len=1024) :: str
    integer             :: ierr

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    str = 'sum(' // trim( var_name) // ') = {int_01} [' // trim( routine_path) // ']'
    call insert_val_into_string_int( str, '{int_01}', sum_d)

    call log_checksum( str)

  end subroutine log_checksum_int

  subroutine log_checksum_dp( sum_d, var_name)

    ! In/output variables:
    real(dp),         intent(in) :: sum_d
    character(len=*), intent(in) :: var_name

    ! Local variables:
    character(len=1024) :: str
    integer             :: ierr

    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    str = 'sum(' // trim( var_name) // ') = {int_01} [' // trim( routine_path) // ']'
    call insert_val_into_string_dp( str, '{int_01}', sum_d)

    call log_checksum( str)

  end subroutine log_checksum_dp

  subroutine log_checksum( str)
    character(len=*), intent(in) :: str
    if (par%primary) then
      ! Write to terminal
      write(0,*) trim(str)
      ! Write to logfile
      if (.not. checksum_logfile_exists) call create_checksum_logfile
      write(unit_checksum_logfile,*) trim(str)
    end if
  end subroutine log_checksum

  subroutine create_checksum_logfile
    filename_checksum_logfile = trim( C%output_dir) // '/checksum_logfile.txt'
    open( file = trim( filename_checksum_logfile), newunit = unit_checksum_logfile)
    write( unit_checksum_logfile,*) '-= Checksum logfile =-'
  end subroutine create_checksum_logfile

end module checksum_mod
