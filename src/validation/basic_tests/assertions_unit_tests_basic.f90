module assertions_unit_tests_basic

  use control_resources_and_error_messaging, only: crash
  use assertions_unit_tests_output, only: write_unit_test_result

  implicit none

  private

  public :: ASSERTION, UNIT_TEST, process_test_result

  ! Module variables
  integer, parameter :: ASSERTION = 1
  integer, parameter :: UNIT_TEST = 2

contains

  ! ===== Process test results =====
  ! ================================

    !> Process the result of an assertion/unit test
    subroutine process_test_result( test_mode, test_result, message)
      ! In/output variables:
      integer,          intent(in   ) :: test_mode
      logical,          intent(in   ) :: test_result
      character(len=*), intent(in   ) :: message

      if (test_mode == ASSERTION) then
        call process_test_result_assertion( test_result, message)
      elseif (test_mode == UNIT_TEST) then
        call process_test_result_unit_test( test_result, message)
      else
        call crash('validation/assertions_unit_tests_basic/process_test_result: test_mode should be either ASSERTION or UNIT_TEST!')
      end if

    end subroutine process_test_result

    !> Process the result of an assertion
    subroutine process_test_result_assertion( test_result, message)
      ! In/output variables:
      logical,          intent(in   ) :: test_result
      character(len=*), intent(in   ) :: message

      if (test_result .eqv. .false.) then
        call crash('Failed assertion: "' // trim(message) // '"')
      end if

    end subroutine process_test_result_assertion

    !> Process the result of a unit test
    subroutine process_test_result_unit_test( test_result, test_name)

      ! In/output variables:
      logical,          intent(in   ) :: test_result
      character(len=*), intent(in   ) :: test_name
      ! Local variables:

      call write_unit_test_result( test_result, test_name)

    end subroutine process_test_result_unit_test

end module assertions_unit_tests_basic
