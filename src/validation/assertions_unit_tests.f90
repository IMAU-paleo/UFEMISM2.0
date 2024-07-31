module assertions_unit_tests

  ! The main assertions/unit tests module.
  !
  ! The "test_XX" generic interfaces choose the appropriate subroutines automatically.
  !
  ! The user must provide the "test_mode" inoput variable, which should equal either
  ! ASSERTION or UNIT_TEST (which can be imported from thsi module). ASSERTION will
  ! crash the model if the test fails; UNIT_TEST will not crash the model, but will
  ! write the test result (both pass and fail) to output.

  use precisions                                             , only: dp
  use control_resources_and_error_messaging                  , only: crash

  implicit none

  private

  public :: ASSERTION, UNIT_TEST
  public :: test_eqv, test_neqv, test_eq, test_neq, test_gt, test_lt, test_ge, test_le, test_ge_le, test_tol

  ! Module variables
  integer, parameter :: ASSERTION = 1
  integer, parameter :: UNIT_TEST = 2

  ! Generic interfaces to specific test routines
  ! ============================================

  !> Test if a == b
  interface test_eqv
    procedure test_eqv_logical_0D
    procedure test_eqv_logical_1D_scalar, test_eqv_logical_1D_array
    procedure test_eqv_logical_2D_scalar, test_eqv_logical_2D_array
    procedure test_eqv_logical_3D_scalar, test_eqv_logical_3D_array
    procedure test_eqv_logical_4D_scalar, test_eqv_logical_4D_array
  end interface test_eqv

  !> Test if a /= b
  interface test_neqv
    procedure test_neqv_logical_0D
    procedure test_neqv_logical_1D_scalar, test_neqv_logical_1D_array
    procedure test_neqv_logical_2D_scalar, test_neqv_logical_2D_array
    procedure test_neqv_logical_3D_scalar, test_neqv_logical_3D_array
    procedure test_neqv_logical_4D_scalar, test_neqv_logical_4D_array
  end interface test_neqv

  !> Test if a == b
  interface test_eq
    procedure test_eq_int_0D
    procedure test_eq_int_1D_scalar, test_eq_int_1D_array
    procedure test_eq_int_2D_scalar, test_eq_int_2D_array
    procedure test_eq_int_3D_scalar, test_eq_int_3D_array
    procedure test_eq_int_4D_scalar, test_eq_int_4D_array
    procedure test_eq_dp_0D
    procedure test_eq_dp_1D_scalar, test_eq_dp_1D_array
    procedure test_eq_dp_2D_scalar, test_eq_dp_2D_array
    procedure test_eq_dp_3D_scalar, test_eq_dp_3D_array
    procedure test_eq_dp_4D_scalar, test_eq_dp_4D_array
  end interface test_eq

  !> Test if a /= b
  interface test_neq
    procedure test_neq_int_0D
    procedure test_neq_int_1D_scalar, test_neq_int_1D_array
    procedure test_neq_int_2D_scalar, test_neq_int_2D_array
    procedure test_neq_int_3D_scalar, test_neq_int_3D_array
    procedure test_neq_int_4D_scalar, test_neq_int_4D_array
    procedure test_neq_dp_0D
    procedure test_neq_dp_1D_scalar, test_neq_dp_1D_array
    procedure test_neq_dp_2D_scalar, test_neq_dp_2D_array
    procedure test_neq_dp_3D_scalar, test_neq_dp_3D_array
    procedure test_neq_dp_4D_scalar, test_neq_dp_4D_array
  end interface test_neq

  !> Test if a > b
  interface test_gt
    procedure test_gt_int_0D
    procedure test_gt_int_1D_scalar, test_gt_int_1D_array
    procedure test_gt_int_2D_scalar, test_gt_int_2D_array
    procedure test_gt_int_3D_scalar, test_gt_int_3D_array
    procedure test_gt_int_4D_scalar, test_gt_int_4D_array
    procedure test_gt_dp_0D
    procedure test_gt_dp_1D_scalar, test_gt_dp_1D_array
    procedure test_gt_dp_2D_scalar, test_gt_dp_2D_array
    procedure test_gt_dp_3D_scalar, test_gt_dp_3D_array
    procedure test_gt_dp_4D_scalar, test_gt_dp_4D_array
  end interface test_gt

  !> Test if a < b
  interface test_lt
    procedure test_lt_int_0D
    procedure test_lt_int_1D_scalar, test_lt_int_1D_array
    procedure test_lt_int_2D_scalar, test_lt_int_2D_array
    procedure test_lt_int_3D_scalar, test_lt_int_3D_array
    procedure test_lt_int_4D_scalar, test_lt_int_4D_array
    procedure test_lt_dp_0D
    procedure test_lt_dp_1D_scalar, test_lt_dp_1D_array
    procedure test_lt_dp_2D_scalar, test_lt_dp_2D_array
    procedure test_lt_dp_3D_scalar, test_lt_dp_3D_array
    procedure test_lt_dp_4D_scalar, test_lt_dp_4D_array
  end interface test_lt

  !> Test if a >= b
  interface test_ge
    procedure test_ge_int_0D
    procedure test_ge_int_1D_scalar, test_ge_int_1D_array
    procedure test_ge_int_2D_scalar, test_ge_int_2D_array
    procedure test_ge_int_3D_scalar, test_ge_int_3D_array
    procedure test_ge_int_4D_scalar, test_ge_int_4D_array
    procedure test_ge_dp_0D
    procedure test_ge_dp_1D_scalar, test_ge_dp_1D_array
    procedure test_ge_dp_2D_scalar, test_ge_dp_2D_array
    procedure test_ge_dp_3D_scalar, test_ge_dp_3D_array
    procedure test_ge_dp_4D_scalar, test_ge_dp_4D_array
  end interface test_ge

  !> Test if a <= b
  interface test_le
    procedure test_le_int_0D
    procedure test_le_int_1D_scalar, test_le_int_1D_array
    procedure test_le_int_2D_scalar, test_le_int_2D_array
    procedure test_le_int_3D_scalar, test_le_int_3D_array
    procedure test_le_int_4D_scalar, test_le_int_4D_array
    procedure test_le_dp_0D
    procedure test_le_dp_1D_scalar, test_le_dp_1D_array
    procedure test_le_dp_2D_scalar, test_le_dp_2D_array
    procedure test_le_dp_3D_scalar, test_le_dp_3D_array
    procedure test_le_dp_4D_scalar, test_le_dp_4D_array
  end interface test_le

  !> Test if a >= b1 && a <= b2
  interface test_ge_le
    procedure test_ge_le_int_0D
    procedure test_ge_le_int_1D_scalar, test_ge_le_int_1D_array
    procedure test_ge_le_int_2D_scalar, test_ge_le_int_2D_array
    procedure test_ge_le_int_3D_scalar, test_ge_le_int_3D_array
    procedure test_ge_le_int_4D_scalar, test_ge_le_int_4D_array
    procedure test_ge_le_dp_0D
    procedure test_ge_le_dp_1D_scalar, test_ge_le_dp_1D_array
    procedure test_ge_le_dp_2D_scalar, test_ge_le_dp_2D_array
    procedure test_ge_le_dp_3D_scalar, test_ge_le_dp_3D_array
    procedure test_ge_le_dp_4D_scalar, test_ge_le_dp_4D_array
  end interface test_ge_le

  !> Test if a >= (b - tol) && a <= (b + tol)
  interface test_tol
    procedure test_tol_int_0D
    procedure test_tol_int_1D_scalar, test_tol_int_1D_array
    procedure test_tol_int_2D_scalar, test_tol_int_2D_array
    procedure test_tol_int_3D_scalar, test_tol_int_3D_array
    procedure test_tol_int_4D_scalar, test_tol_int_4D_array
    procedure test_tol_dp_0D
    procedure test_tol_dp_1D_scalar, test_tol_dp_1D_array
    procedure test_tol_dp_2D_scalar, test_tol_dp_2D_array
    procedure test_tol_dp_3D_scalar, test_tol_dp_3D_array
    procedure test_tol_dp_4D_scalar, test_tol_dp_4D_array
  end interface test_tol

contains

! ===== Logical =====
! ===================

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

! ===== Integer =====
! ===================

  ! ===== 0-D =====
  ! ===============

  !> Test if a == b
  subroutine test_eq_int_0D( a, b, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a == b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eq_int_0D

  !> Test if a /= b
  subroutine test_neq_int_0D( a, b, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a /= b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neq_int_0D

  !> Test if a > b
  subroutine test_gt_int_0D( a, b, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a > b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_gt_int_0D

  !> Test if a < b
  subroutine test_lt_int_0D( a, b, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a < b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_lt_int_0D

  !> Test if a >= b
  subroutine test_ge_int_0D( a, b, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a >= b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_int_0D

  !> Test if a <= b
  subroutine test_le_int_0D( a, b, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a <= b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_le_int_0D

  !> Test if a >= b1 && a <= b2
  subroutine test_ge_le_int_0D( a, b1, b2, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b1, b2
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that b2 >= b1
    call test_ge( b2, b1, test_mode, message)

    test_result = a >= b1 .and. a <= b2

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_le_int_0D

  !> Test if a >= (b - tol) && a <= (b + tol)
  subroutine test_tol_int_0D( a, b, tol, test_mode, message)
    ! In/output variables:
    integer,          intent(in   ) :: a
    integer,          intent(in   ) :: b, tol
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that tol >= 0
    call test_ge( tol, 0, test_mode, message)

    test_result = a >= (b - tol) .and. a <= (b + tol)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_tol_int_0D

  ! ===== 1-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:) == b
  subroutine test_eq_int_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a == b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eq_int_1D_scalar

  !> Test if a(:) /= b
  subroutine test_neq_int_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a /= b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neq_int_1D_scalar

  !> Test if a(:) > b
  subroutine test_gt_int_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a > b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_gt_int_1D_scalar

  !> Test if a(:) < b
  subroutine test_lt_int_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a < b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_lt_int_1D_scalar

  !> Test if a(:) >= b
  subroutine test_ge_int_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a >= b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_int_1D_scalar

  !> Test if a(:) <= b
  subroutine test_le_int_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a <= b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_le_int_1D_scalar

  !> Test if a(:) >= b1 && a(:) <= b2
  subroutine test_ge_le_int_1D_scalar( a, b1, b2, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b1, b2
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that b2 >= b1
    call test_ge( b2, b1, test_mode, message)

    test_result = all(a >= b1 .and. a <= b2)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_le_int_1D_scalar

  !> Test if a(:) >= (b - tol) && a(:) <= (b + tol)
  subroutine test_tol_int_1D_scalar( a, b, tol, test_mode, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   ) :: a
    integer,                        intent(in   ) :: b, tol
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that tol >= 0
    call test_ge( tol, 0, test_mode, message)

    test_result = all(a >= (b - tol) .and. a <= (b + tol))

    call process_test_result( test_mode, test_result, message)

  end subroutine test_tol_int_1D_scalar

  ! ===== Array =====
  ! ==================

!> Test if a(:) == b(:)
subroutine test_eq_int_1D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_int_1D_array

!> Test if a(:) /= b(:)
subroutine test_neq_int_1D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_int_1D_array

!> Test if a(:) > b(:)
subroutine test_gt_int_1D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_int_1D_array

!> Test if a(:) < b(:)
subroutine test_lt_int_1D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_int_1D_array

!> Test if a(:) >= b(:)
subroutine test_ge_int_1D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_int_1D_array

!> Test if a(:) <= b(:)
subroutine test_le_int_1D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_int_1D_array

!> Test if a(:) >= b1(:) && a(:) <= b2(:)
subroutine test_ge_le_int_1D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b1, b2
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_int_1D_array

!> Test if a(:) >= (b(:) - tol(:)) && a(:) <= (b(:) + tol(:))
subroutine test_tol_int_1D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  integer,          dimension(:), intent(in   ) :: a, b, tol
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_int_1D_array

! ===== 2-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

!> Test if a(:,:) == b
subroutine test_eq_int_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_int_2D_scalar

!> Test if a(:,:) /= b
subroutine test_neq_int_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_int_2D_scalar

!> Test if a(:,:) > b
subroutine test_gt_int_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_int_2D_scalar

!> Test if a(:,:) < b
subroutine test_lt_int_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_int_2D_scalar

!> Test if a(:,:) >= b
subroutine test_ge_int_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_int_2D_scalar

!> Test if a(:,:) <= b
subroutine test_le_int_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_int_2D_scalar

!> Test if a(:,:) >= b1 && a(:,:) <= b2
subroutine test_ge_le_int_2D_scalar( a, b1, b2, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b1, b2
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_int_2D_scalar

!> Test if a(:,:) >= (b - tol) && a(:,:) <= (b + tol)
subroutine test_tol_int_2D_scalar( a, b, tol, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a
  integer,                          intent(in   ) :: b, tol
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_int_2D_scalar

    ! ===== Array =====
    ! ==================

!> Test if a(:,:) == b(:,:)
subroutine test_eq_int_2D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_int_2D_array

!> Test if a(:,:) /= b(:,:)
subroutine test_neq_int_2D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_int_2D_array

!> Test if a(:,:) > b(:,:)
subroutine test_gt_int_2D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_int_2D_array

!> Test if a(:,:) < b(:,:)
subroutine test_lt_int_2D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_int_2D_array

!> Test if a(:,:) >= b(:,:)
subroutine test_ge_int_2D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_int_2D_array

!> Test if a(:,:) <= b(:,:)
subroutine test_le_int_2D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_int_2D_array

!> Test if a(:,:) >= b1(:,:) && a(:,:) <= b2(:,:)
subroutine test_ge_le_int_2D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b1, b2
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_int_2D_array

!> Test if a(:,:) >= (b(:,:) - tol(:,:)) && a(:,:) <= (b(:,:) + tol(:,:))
subroutine test_tol_int_2D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:), intent(in   ) :: a, b, tol
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_int_2D_array

! ===== 3-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

!> Test if a(:,:,:) == b
subroutine test_eq_int_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_int_3D_scalar

!> Test if a(:,:,:) /= b
subroutine test_neq_int_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_int_3D_scalar

!> Test if a(:,:,:) > b
subroutine test_gt_int_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_int_3D_scalar

!> Test if a(:,:,:) < b
subroutine test_lt_int_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_int_3D_scalar

!> Test if a(:,:,:) >= b
subroutine test_ge_int_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_int_3D_scalar

!> Test if a(:,:,:) <= b
subroutine test_le_int_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_int_3D_scalar

!> Test if a(:,:,:) >= b1 && a(:,:,:) <= b2
subroutine test_ge_le_int_3D_scalar( a, b1, b2, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b1, b2
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_int_3D_scalar

!> Test if a(:,:,:) >= (b - tol) && a(:,:,:) <= (b + tol)
subroutine test_tol_int_3D_scalar( a, b, tol, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a
  integer,                            intent(in   ) :: b, tol
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_int_3D_scalar

    ! ===== Array =====
    ! ==================

!> Test if a(:,:,:) == b(:,:,:)
subroutine test_eq_int_3D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_int_3D_array

!> Test if a(:,:,:) /= b(:,:,:)
subroutine test_neq_int_3D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_int_3D_array

!> Test if a(:,:,:) > b(:,:,:)
subroutine test_gt_int_3D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_int_3D_array

!> Test if a(:,:,:) < b(:,:,:)
subroutine test_lt_int_3D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_int_3D_array

!> Test if a(:,:,:) >= b(:,:,:)
subroutine test_ge_int_3D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_int_3D_array

!> Test if a(:,:,:) <= b(:,:,:)
subroutine test_le_int_3D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_int_3D_array

!> Test if a(:,:,:) >= b1(:,:,:) && a(:,:,:) <= b2(:,:,:)
subroutine test_ge_le_int_3D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b1, b2
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_int_3D_array

!> Test if a(:,:,:) >= (b(:,:,:) - tol(:,:,:)) && a(:,:,:) <= (b(:,:,:) + tol(:,:,:))
subroutine test_tol_int_3D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:), intent(in   ) :: a, b, tol
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_int_3D_array

! ===== 4-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

!> Test if a(:,:,:,:) == b
subroutine test_eq_int_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_int_4D_scalar

!> Test if a(:,:,:,:) /= b
subroutine test_neq_int_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_int_4D_scalar

!> Test if a(:,:,:,:) > b
subroutine test_gt_int_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_int_4D_scalar

!> Test if a(:,:,:,:) < b
subroutine test_lt_int_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_int_4D_scalar

!> Test if a(:,:,:,:) >= b
subroutine test_ge_int_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_int_4D_scalar

!> Test if a(:,:,:,:) <= b
subroutine test_le_int_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_int_4D_scalar

!> Test if a(:,:,:,:) >= b1 && a(:,:,:,:) <= b2
subroutine test_ge_le_int_4D_scalar( a, b1, b2, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b1, b2
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_int_4D_scalar

!> Test if a(:,:,:,:) >= (b - tol) && a(:,:,:,:) <= (b + tol)
subroutine test_tol_int_4D_scalar( a, b, tol, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a
  integer,                              intent(in   ) :: b, tol
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_int_4D_scalar

    ! ===== Array =====
    ! ==================

!> Test if a(:,:,:,:) == b(:,:,:,:)
subroutine test_eq_int_4D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_int_4D_array

!> Test if a(:,:,:,:) /= b(:,:,:,:)
subroutine test_neq_int_4D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_int_4D_array

!> Test if a(:,:,:,:) > b(:,:,:,:)
subroutine test_gt_int_4D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_int_4D_array

!> Test if a(:,:,:,:) < b(:,:,:,:)
subroutine test_lt_int_4D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_int_4D_array

!> Test if a(:,:,:,:) >= b(:,:,:,:)
subroutine test_ge_int_4D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_int_4D_array

!> Test if a(:,:,:,:) <= b(:,:,:,:)
subroutine test_le_int_4D_array( a, b, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_int_4D_array

!> Test if a(:,:,:,:) >= b1(:,:,:,:) && a(:,:,:,:) <= b2(:,:,:,:)
subroutine test_ge_le_int_4D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b1, b2
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_int_4D_array

!> Test if a(:,:,:,:) >= (b(:,:,:,:) - tol(:,:,:,:)) && a(:,:,:,:) <= (b(:,:,:,:) + tol(:,:,:,:))
subroutine test_tol_int_4D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  integer,          dimension(:,:,:,:), intent(in   ) :: a, b, tol
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_int_4D_array

! ===== Real =====
! ================

  ! ===== 0-D =====
  ! ===============

  !> Test if a == b
  subroutine test_eq_dp_0D( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a == b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eq_dp_0D

  !> Test if a /= b
  subroutine test_neq_dp_0D( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a /= b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neq_dp_0D

  !> Test if a > b
  subroutine test_gt_dp_0D( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a > b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_gt_dp_0D

  !> Test if a < b
  subroutine test_lt_dp_0D( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a < b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_lt_dp_0D

  !> Test if a >= b
  subroutine test_ge_dp_0D( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a >= b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_dp_0D

  !> Test if a <= b
  subroutine test_le_dp_0D( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = a <= b

    call process_test_result( test_mode, test_result, message)

  end subroutine test_le_dp_0D

  !> Test if a >= b1 && a <= b2
  subroutine test_ge_le_dp_0D( a, b1, b2, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b1, b2
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that b2 >= b1
    call test_ge( b2, b1, test_mode, message)

    test_result = a >= b1 .and. a <= b2

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_le_dp_0D

  !> Test if a >= (b - tol) && a <= (b + tol)
  subroutine test_tol_dp_0D( a, b, tol, test_mode, message)
    ! In/output variables:
    real(dp),         intent(in   ) :: a
    real(dp),         intent(in   ) :: b, tol
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that tol >= 0
    call test_ge( tol, 0._dp, test_mode, message)

    test_result = a >= (b - tol) .and. a <= (b + tol)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_tol_dp_0D

  ! ===== 1-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:) == b
  subroutine test_eq_dp_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a == b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_eq_dp_1D_scalar

  !> Test if a(:) /= b
  subroutine test_neq_dp_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a /= b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_neq_dp_1D_scalar

  !> Test if a(:) > b
  subroutine test_gt_dp_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a > b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_gt_dp_1D_scalar

  !> Test if a(:) < b
  subroutine test_lt_dp_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a < b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_lt_dp_1D_scalar

  !> Test if a(:) >= b
  subroutine test_ge_dp_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a >= b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_dp_1D_scalar

  !> Test if a(:) <= b
  subroutine test_le_dp_1D_scalar( a, b, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = all(a <= b)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_le_dp_1D_scalar

  !> Test if a(:) >= b1 && a(:) <= b2
  subroutine test_ge_le_dp_1D_scalar( a, b1, b2, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b1, b2
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that b2 >= b1
    call test_ge( b2, b1, test_mode, message)

    test_result = all(a >= b1 .and. a <= b2)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_ge_le_dp_1D_scalar

  !> Test if a(:) >= (b - tol) && a(:) <= (b + tol)
  subroutine test_tol_dp_1D_scalar( a, b, tol, test_mode, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   ) :: a
    real(dp),                       intent(in   ) :: b, tol
    integer,                        intent(in   ) :: test_mode
    character(len=*),               intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    ! First test that tol >= 0
    call test_ge( tol, 0._dp, test_mode, message)

    test_result = all(a >= (b - tol) .and. a <= (b + tol))

    call process_test_result( test_mode, test_result, message)

  end subroutine test_tol_dp_1D_scalar

  ! ===== Array =====
  ! ==================

!> Test if a(:) == b(:)
subroutine test_eq_dp_1D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_dp_1D_array

!> Test if a(:) /= b(:)
subroutine test_neq_dp_1D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_dp_1D_array

!> Test if a(:) > b(:)
subroutine test_gt_dp_1D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_dp_1D_array

!> Test if a(:) < b(:)
subroutine test_lt_dp_1D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_dp_1D_array

!> Test if a(:) >= b(:)
subroutine test_ge_dp_1D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_dp_1D_array

!> Test if a(:) <= b(:)
subroutine test_le_dp_1D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_dp_1D_array

!> Test if a(:) >= b1(:) && a(:) <= b2(:)
subroutine test_ge_le_dp_1D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b1, b2
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_1D_array

!> Test if a(:) >= (b(:) - tol(:)) && a(:) <= (b(:) + tol(:))
subroutine test_tol_dp_1D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b, tol
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0._dp, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_1D_array

! ===== 2-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

!> Test if a(:,:) == b
subroutine test_eq_dp_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_dp_2D_scalar

!> Test if a(:,:) /= b
subroutine test_neq_dp_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_dp_2D_scalar

!> Test if a(:,:) > b
subroutine test_gt_dp_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_dp_2D_scalar

!> Test if a(:,:) < b
subroutine test_lt_dp_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_dp_2D_scalar

!> Test if a(:,:) >= b
subroutine test_ge_dp_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_dp_2D_scalar

!> Test if a(:,:) <= b
subroutine test_le_dp_2D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_dp_2D_scalar

!> Test if a(:,:) >= b1 && a(:,:) <= b2
subroutine test_ge_le_dp_2D_scalar( a, b1, b2, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b1, b2
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_2D_scalar

!> Test if a(:,:) >= (b - tol) && a(:,:) <= (b + tol)
subroutine test_tol_dp_2D_scalar( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a
  real(dp),                         intent(in   ) :: b, tol
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0._dp, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_2D_scalar

    ! ===== Array =====
    ! ==================

!> Test if a(:,:) == b(:,:)
subroutine test_eq_dp_2D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_dp_2D_array

!> Test if a(:,:) /= b(:,:)
subroutine test_neq_dp_2D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_dp_2D_array

!> Test if a(:,:) > b(:,:)
subroutine test_gt_dp_2D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_dp_2D_array

!> Test if a(:,:) < b(:,:)
subroutine test_lt_dp_2D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_dp_2D_array

!> Test if a(:,:) >= b(:,:)
subroutine test_ge_dp_2D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_dp_2D_array

!> Test if a(:,:) <= b(:,:)
subroutine test_le_dp_2D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_dp_2D_array

!> Test if a(:,:) >= b1(:,:) && a(:,:) <= b2(:,:)
subroutine test_ge_le_dp_2D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b1, b2
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_2D_array

!> Test if a(:,:) >= (b(:,:) - tol(:,:)) && a(:,:) <= (b(:,:) + tol(:,:))
subroutine test_tol_dp_2D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b, tol
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0._dp, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_2D_array

! ===== 3-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

!> Test if a(:,:,:) == b
subroutine test_eq_dp_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_dp_3D_scalar

!> Test if a(:,:,:) /= b
subroutine test_neq_dp_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_dp_3D_scalar

!> Test if a(:,:,:) > b
subroutine test_gt_dp_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_dp_3D_scalar

!> Test if a(:,:,:) < b
subroutine test_lt_dp_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_dp_3D_scalar

!> Test if a(:,:,:) >= b
subroutine test_ge_dp_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_dp_3D_scalar

!> Test if a(:,:,:) <= b
subroutine test_le_dp_3D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_dp_3D_scalar

!> Test if a(:,:,:) >= b1 && a(:,:,:) <= b2
subroutine test_ge_le_dp_3D_scalar( a, b1, b2, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b1, b2
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_3D_scalar

!> Test if a(:,:,:) >= (b - tol) && a(:,:,:) <= (b + tol)
subroutine test_tol_dp_3D_scalar( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a
  real(dp),                           intent(in   ) :: b, tol
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0._dp, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_3D_scalar

    ! ===== Array =====
    ! ==================

!> Test if a(:,:,:) == b(:,:,:)
subroutine test_eq_dp_3D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_dp_3D_array

!> Test if a(:,:,:) /= b(:,:,:)
subroutine test_neq_dp_3D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_dp_3D_array

!> Test if a(:,:,:) > b(:,:,:)
  subroutine test_gt_dp_3D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_dp_3D_array

!> Test if a(:,:,:) < b(:,:,:)
subroutine test_lt_dp_3D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_dp_3D_array

!> Test if a(:,:,:) >= b(:,:,:)
subroutine test_ge_dp_3D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_dp_3D_array

!> Test if a(:,:,:) <= b(:,:,:)
subroutine test_le_dp_3D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_dp_3D_array

!> Test if a(:,:,:) >= b1(:,:,:) && a(:,:,:) <= b2(:,:,:)
subroutine test_ge_le_dp_3D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b1, b2
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_3D_array

!> Test if a(:,:,:) >= (b(:,:,:) - tol(:,:,:)) && a(:,:,:) <= (b(:,:,:) + tol(:,:,:))
subroutine test_tol_dp_3D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b, tol
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0._dp, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_3D_array

! ===== 4-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

!> Test if a(:,:,:,:) == b
subroutine test_eq_dp_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_dp_4D_scalar

!> Test if a(:,:,:,:) /= b
subroutine test_neq_dp_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_dp_4D_scalar

!> Test if a(:,:,:,:) > b
subroutine test_gt_dp_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_dp_4D_scalar

!> Test if a(:,:,:,:) < b
subroutine test_lt_dp_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_dp_4D_scalar

!> Test if a(:,:,:,:) >= b
subroutine test_ge_dp_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_dp_4D_scalar

!> Test if a(:,:,:,:) <= b
subroutine test_le_dp_4D_scalar( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_dp_4D_scalar

!> Test if a(:,:,:,:) >= b1 && a(:,:,:,:) <= b2
subroutine test_ge_le_dp_4D_scalar( a, b1, b2, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b1, b2
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_4D_scalar

!> Test if a(:,:,:,:) >= (b - tol) && a(:,:,:,:) <= (b + tol)
subroutine test_tol_dp_4D_scalar( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a
  real(dp),                             intent(in   ) :: b, tol
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0._dp, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_4D_scalar

    ! ===== Array =====
    ! ==================

!> Test if a(:,:,:,:) == b(:,:,:,:)
subroutine test_eq_dp_4D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a == b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_eq_dp_4D_array

!> Test if a(:,:,:,:) /= b(:,:,:,:)
subroutine test_neq_dp_4D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a /= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_neq_dp_4D_array

!> Test if a(:,:,:,:) > b(:,:,:,:)
subroutine test_gt_dp_4D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a > b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_gt_dp_4D_array

!> Test if a(:,:,:,:) < b(:,:,:,:)
subroutine test_lt_dp_4D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a < b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_lt_dp_4D_array

!> Test if a(:,:,:,:) >= b(:,:,:,:)
subroutine test_ge_dp_4D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_dp_4D_array

!> Test if a(:,:,:,:) <= b(:,:,:,:)
subroutine test_le_dp_4D_array( a, b, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a <= b)

  call process_test_result( test_mode, test_result, message)

end subroutine test_le_dp_4D_array

!> Test if a(:,:,:,:) >= b1(:,:,:,:) && a(:,:,:,:) <= b2(:,:,:,:)
subroutine test_ge_le_dp_4D_array( a, b1, b2, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b1, b2
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that b2 >= b1
  call test_ge( b2, b1, test_mode, message)

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_4D_array

!> Test if a(:,:,:,:) >= (b(:,:,:,:) - tol(:,:,:,:)) && a(:,:,:,:) <= (b(:,:,:,:) + tol(:,:,:,:))
subroutine test_tol_dp_4D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b, tol
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  ! First test that tol >= 0
  call test_ge( tol, 0._dp, test_mode, message)

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_4D_array

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
      call crash('test_mode should be either ASSERTION or UNIT_TEST!')
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
    use unit_tests_output, only: write_unit_test_result

    ! In/output variables:
    logical,          intent(in   ) :: test_result
    character(len=*), intent(in   ) :: test_name
    ! Local variables:

    call write_unit_test_result( test_result, test_name)

  end subroutine process_test_result_unit_test

end module assertions_unit_tests
