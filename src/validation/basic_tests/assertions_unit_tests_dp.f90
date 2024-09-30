module assertions_unit_tests_dp

  ! The assertions/unit tests for double precision values.

  use assertions_unit_tests_basic, only: ASSERTION, UNIT_TEST, process_test_result
  use precisions, only: dp

  implicit none

contains

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

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_1D_array

!> Test if a(:) >= (b(:) - tol) && a(:) <= (b(:) + tol)
subroutine test_tol_dp_1D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:), intent(in   ) :: a, b
  real(dp),                       intent(in   ) :: tol
  integer,                        intent(in   ) :: test_mode
  character(len=*),               intent(in   ) :: message
  ! Local variables:
  logical :: test_result

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

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_2D_array

!> Test if a(:,:) >= (b(:,:) - tol) && a(:,:) <= (b(:,:) + tol)
subroutine test_tol_dp_2D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:), intent(in   ) :: a, b
  real(dp),                         intent(in   ) :: tol
  integer,                          intent(in   ) :: test_mode
  character(len=*),                 intent(in   ) :: message
  ! Local variables:
  logical :: test_result

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

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_3D_array

!> Test if a(:,:,:) >= (b(:,:,:) - tol) && a(:,:,:) <= (b(:,:,:) + tol)
subroutine test_tol_dp_3D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:), intent(in   ) :: a, b
  real(dp),                           intent(in   ) :: tol
  integer,                            intent(in   ) :: test_mode
  character(len=*),                   intent(in   ) :: message
  ! Local variables:
  logical :: test_result

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

  test_result = all(a >= b1 .and. a <= b2)

  call process_test_result( test_mode, test_result, message)

end subroutine test_ge_le_dp_4D_array

!> Test if a(:,:,:,:) >= (b(:,:,:,:) - tol) && a(:,:,:,:) <= (b(:,:,:,:) + tol)
subroutine test_tol_dp_4D_array( a, b, tol, test_mode, message)
  ! In/output variables:
  real(dp),         dimension(:,:,:,:), intent(in   ) :: a, b
  real(dp),                             intent(in   ) :: tol
  integer,                              intent(in   ) :: test_mode
  character(len=*),                     intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = all(a >= (b - tol) .and. a <= (b + tol))

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_dp_4D_array

end module assertions_unit_tests_dp
