module tests_logical

  ! The assertions/unit tests for logical values.

  use assertions_unit_tests_basic, only: ASSERTION, UNIT_TEST, process_test_result

  implicit none

contains

  ! ===== 0-D =====
  ! ===============

  !> Test if a == b
  subroutine test_eqv_logical_0D( a, b, test_mode, message)
    ! In/output variables:
    logical,          intent(in   ) :: a
    logical,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a .eqv. b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_0D

  !> Test if a /= b
  subroutine test_neqv_logical_0D( a, b, test_mode, message)
    ! In/output variables:
    logical,          intent(in   ) :: a
    logical,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a .neqv. b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_0D

  ! ===== 1-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:) == b
  subroutine test_eqv_logical_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:), intent(in   ) :: a
    logical,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_1D_scalar

  !> Test if a(:) /= b
  subroutine test_neqv_logical_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:), intent(in   ) :: a
    logical,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_1D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:) == b(:)
  subroutine test_eqv_logical_1D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:), intent(in   ) :: a, b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_1D_array

  !> Test if a(:) /= b
  subroutine test_neqv_logical_1D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:), intent(in   ) :: a, b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_1D_array

  ! ===== 2-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:,:) == b
  subroutine test_eqv_logical_2D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:), intent(in   ) :: a
    logical,                          intent(in   ) :: b
    integer,                          intent(in   ) :: test_mode
    character(len=*),                 intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_2D_scalar

  !> Test if a(:,:) /= b
  subroutine test_neqv_logical_2D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:), intent(in   ) :: a
    logical,                          intent(in   ) :: b
    integer,                          intent(in   ) :: test_mode
    character(len=*),                 intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_2D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:,:) == b(:,:)
  subroutine test_eqv_logical_2D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:), intent(in   ) :: a, b
    integer,                          intent(in   ) :: test_mode
    character(len=*),                 intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_2D_array

  !> Test if a(:,:) /= b
  subroutine test_neqv_logical_2D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:), intent(in   ) :: a, b
    integer,                          intent(in   ) :: test_mode
    character(len=*),                 intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_2D_array

  ! ===== 3-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:,:,:) == b
  subroutine test_eqv_logical_3D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:), intent(in   ) :: a
    logical,                            intent(in   ) :: b
    integer,                            intent(in   ) :: test_mode
    character(len=*),                   intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_3D_scalar

  !> Test if a(:,:,:) /= b
  subroutine test_neqv_logical_3D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:), intent(in   ) :: a
    logical,                            intent(in   ) :: b
    integer,                            intent(in   ) :: test_mode
    character(len=*),                   intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_3D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:,:,:) == b(:,:,:)
  subroutine test_eqv_logical_3D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:), intent(in   ) :: a, b
    integer,                            intent(in   ) :: test_mode
    character(len=*),                   intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_3D_array

  !> Test if a(:,:,:) /= b
  subroutine test_neqv_logical_3D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:), intent(in   ) :: a, b
    integer,                            intent(in   ) :: test_mode
    character(len=*),                   intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_3D_array

  ! ===== 3-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:,:,:,:) == b
  subroutine test_eqv_logical_4D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:,:), intent(in   ) :: a
    logical,                              intent(in   ) :: b
    integer,                              intent(in   ) :: test_mode
    character(len=*),                     intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_4D_scalar

  !> Test if a(:,:,:,:) /= b
  subroutine test_neqv_logical_4D_scalar( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:,:), intent(in   ) :: a
    logical,                              intent(in   ) :: b
    integer,                              intent(in   ) :: test_mode
    character(len=*),                     intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_4D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:,:,:,:) == b(:,:,:,:)
  subroutine test_eqv_logical_4D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:,:), intent(in   ) :: a, b
    integer,                              intent(in   ) :: test_mode
    character(len=*),                     intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .eqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eqv_logical_4D_array

  !> Test if a(:,:,:,:) /= b
  subroutine test_neqv_logical_4D_array( a, b, test_mode, message)
    ! In/output variables:
    logical,          dimension(:,:,:,:), intent(in   ) :: a, b
    integer,                              intent(in   ) :: test_mode
    character(len=*),                     intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all( a .neqv. b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neqv_logical_4D_array

end module tests_logical
