module ut_basic

  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string

  implicit none

  private

  public :: foldername_unit_tests_output, filename_unit_tests_output, unit_test

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

end module ut_basic