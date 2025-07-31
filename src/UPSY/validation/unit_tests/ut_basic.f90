module ut_basic

  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string

  implicit none

  private

  public :: foldername_unit_tests_output, filename_unit_tests_output, unit_test
  public :: create_unit_tests_output_folder, create_unit_tests_output_file

  ! Module variables
  character(len=1024) :: foldername_unit_tests_output
  character(len=1024) :: filename_unit_tests_output
  logical             :: all_unit_tests_passed = .true.

contains

  !> The basic unit test subroutine: test if the provided condition is true, and write the outcome
  !> to both the terminal and the unit test output file, without crashing the model
  subroutine unit_test( test_result, test_name)

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
      if (all_unit_tests_passed) then
        str_terminal = colour_string( 'Unit test passed:', 'green')
      else
        str_terminal = colour_string( 'Unit test passed:', 'yellow')
      end if
    else
      all_unit_tests_passed = .false.
      str_file     = 'Unit test failed:'
      str_terminal = colour_string( 'Unit test failed:', 'red')
    end if

    str_file     = trim( str_file)     // ' ' // trim( test_name)
    str_terminal = trim( str_terminal) // ' ' // trim( test_name)

    if (par%primary) then

      ! Write to terminal
      write(0,*) trim( str_terminal)

      ! Write to file
      open(newunit = io_unit_test_file, file = trim( filename_unit_tests_output), status = "old", action = "write", position = "append", &
        iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not open unit test output file, error message "' // trim(msg) // '"')
      end if
      write(io_unit_test_file,*) trim( str_file)
      close(io_unit_test_file)

    end if

  end subroutine unit_test

  !> Create the unit test output file
  subroutine create_unit_tests_output_file

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_unit_tests_output_folder'
    integer                        :: io_unit_test_file
    integer                        :: stat
    character(len=512)             :: msg

    ! Add routine to path
    call init_routine( routine_name)

    filename_unit_tests_output = trim( foldername_unit_tests_output) // '/unit_tests_output.txt'

    if (par%primary) then
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
  subroutine create_unit_tests_output_folder( foldername)

    ! In/output variables:
    character(len=*), intent(in) :: foldername

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_unit_tests_output_folder'
    logical                        :: ex
    character(len=1024)            :: cwd

    ! Add routine to path
    call init_routine( routine_name)

    foldername_unit_tests_output = foldername

    ! Create the directory
    if (par%primary) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_unit_tests_output) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_unit_tests_output))
      end if

      ! Create output directory
      CALL system('mkdir ' // trim( foldername_unit_tests_output))

      ! Tell the user where it is
      call getcwd( cwd)
      write(0,'(A)') ''
      write(0,'(A)') ' Output directory: ' // colour_string( trim(cwd)//'/'//trim( foldername_unit_tests_output), 'light blue')
      write(0,'(A)') ''

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_unit_tests_output_folder

end module ut_basic