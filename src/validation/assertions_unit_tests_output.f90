module assertions_unit_tests_output

  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string
  use model_configuration, only: C

  implicit none

  private

  public :: create_unit_tests_output_folder, create_unit_tests_output_file, write_unit_test_result

  ! Module variables
  character(len=1024) :: filename_unit_tests_output

contains

  !> Write unit test output to the terminal and to an output file.
  subroutine write_unit_test_result( test_result, test_name)

    ! In/output variables:
    logical,          intent(in   ) :: test_result
    character(len=*), intent(in   ) :: test_name
    ! Local variables:
    character(len=1024) :: str_terminal, str_file
    integer             :: io_unit_test_file
    integer             :: stat
    character(len=512)  :: msg

    if (test_result .eqv. .true.) then
      str_file     = 'Unit test passed:'
      str_terminal = colour_string( 'Unit test passed:', 'green')
    else
      str_file     = 'Unit test failed:'
      str_terminal = colour_string( 'Unit test failed:', 'red')
    end if

    str_file     = trim( str_file)     // ' ' // trim( test_name)
    str_terminal = trim( str_terminal) // ' ' // trim( test_name)

    if (par%master) then

      ! Write to terminal
      write(0,*) trim( str_terminal)

      ! Write to file
      open(newunit = io_unit_test_file, file = filename_unit_tests_output, status = "old", action = "write", position = "append", &
        iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not open unit test output file, error message "' // trim(msg) // '"')
      end if
      write(io_unit_test_file,*) trim( str_file)
      close(io_unit_test_file)

    end if

  end subroutine write_unit_test_result

  !> Create the unit test output file
  subroutine create_unit_tests_output_file

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_unit_tests_output_folder'
    integer                        :: io_unit_test_file
    integer                        :: stat
    character(len=512)             :: msg

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      ! Create filename
      filename_unit_tests_output = trim(C%output_dir) // '/unit_tests_output.txt'
      ! Create file
      open(newunit = io_unit_test_file, file = filename_unit_tests_output, status = "new", action = "write", &
        iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not create unit test output file, error message "' // trim(msg) // '"')
      end if
      close(io_unit_test_file)
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_unit_tests_output_file

  !> Create the unit test output folder
  subroutine create_unit_tests_output_folder

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_unit_tests_output_folder'
    logical                        :: ex
    character(len=1024)            :: cwd

    ! Add routine to path
    call init_routine( routine_name)

    C%output_dir = 'results_unit_tests'

    ! Create the directory
    if (par%master) then

      ! Remove existing folder if necessary
      inquire( file = trim( C%output_dir) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( C%output_dir))
      end if

      ! Create output directory
      CALL system('mkdir ' // trim( C%output_dir))

      ! Tell the user where it is
      call getcwd( cwd)
      write(0,'(A)') ''
      write(0,'(A)') ' Output directory: ' // colour_string( trim(cwd)//'/'//trim( C%output_dir), 'light blue')
      write(0,'(A)') ''

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_unit_tests_output_folder

end module assertions_unit_tests_output