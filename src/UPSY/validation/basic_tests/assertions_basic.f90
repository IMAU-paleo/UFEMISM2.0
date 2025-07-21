module assertions_basic

  use control_resources_and_error_messaging, only: crash

  implicit none

  private

  public :: assert

contains

  !> The basic assertion subroutine: assert that the provided condition is true
  subroutine assert( condition, fail_message)

    !> In/output variables:
    logical,          intent(in) :: condition
    character(len=*), intent(in) :: fail_message

    if (.not. condition) then
      call crash('Assertion failed: '//trim(fail_message))
    end if

  end subroutine assert

end module assertions_basic