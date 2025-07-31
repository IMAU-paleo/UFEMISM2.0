module tests_dp

  ! Basic tests for double precision values.

  use precisions, only: dp

  implicit none

contains

  ! ===== 0-D =====
  ! ===============

  !> Test if a == b
  pure function test_eq_dp_0D( a, b) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: res

    res = a == b

  end function test_eq_dp_0D

  !> Test if a /= b
  pure function test_neq_dp_0D( a, b) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: res

    res = a /= b

  end function test_neq_dp_0D

  !> Test if a > b
  pure function test_gt_dp_0D( a, b) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: res

    res = a > b

  end function test_gt_dp_0D

  !> Test if a < b
  pure function test_lt_dp_0D( a, b) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: res

    res = a < b

  end function test_lt_dp_0D

  !> Test if a >= b
  pure function test_ge_dp_0D( a, b) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: res

    res = a >= b

  end function test_ge_dp_0D

  !> Test if a <= b
  pure function test_le_dp_0D( a, b) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: res

    res = a <= b

  end function test_le_dp_0D

  !> Test if a >= b1 && a <= b2
  pure function test_ge_le_dp_0D( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b1, b2
    logical :: res

    res = a >= b1 .and. a <= b2

  end function test_ge_le_dp_0D

  !> Test if a >= (b - tol) && a <= (b + tol)
  pure function test_tol_dp_0D( a, b, tol) result( res)
    ! In/output variables:
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b, tol
    logical :: res

    res = a >= (b - tol) .and. a <= (b + tol)

  end function test_tol_dp_0D

  ! ===== 1-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:) == b
  pure function test_eq_dp_1D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_1D_scalar

  !> Test if a(:) /= b
  pure function test_neq_dp_1D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_1D_scalar

  !> Test if a(:) > b
  pure function test_gt_dp_1D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_1D_scalar

  !> Test if a(:) < b
  pure function test_lt_dp_1D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_1D_scalar

  !> Test if a(:) >= b
  pure function test_ge_dp_1D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_1D_scalar

  !> Test if a(:) <= b
  pure function test_le_dp_1D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_1D_scalar

  !> Test if a(:) >= b1 && a(:) <= b2
  pure function test_ge_le_dp_1D_scalar( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_1D_scalar

  !> Test if a(:) >= (b - tol) && a(:) <= (b + tol)
  pure function test_tol_dp_1D_scalar( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a
    real(dp),               intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_1D_scalar

  ! ===== Array =====
  ! ==================

  !> Test if a(:) == b(:)
  pure function test_eq_dp_1D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_1D_array

  !> Test if a(:) /= b(:)
  pure function test_neq_dp_1D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_1D_array

  !> Test if a(:) > b(:)
  pure function test_gt_dp_1D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_1D_array

  !> Test if a(:) < b(:)
  pure function test_lt_dp_1D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_1D_array

  !> Test if a(:) >= b(:)
  pure function test_ge_dp_1D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_1D_array

  !> Test if a(:) <= b(:)
  pure function test_le_dp_1D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_1D_array

  !> Test if a(:) >= b1(:) && a(:) <= b2(:)
  pure function test_ge_le_dp_1D_array( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_1D_array

  !> Test if a(:) >= (b(:) - tol) && a(:) <= (b(:) + tol)
  pure function test_tol_dp_1D_array( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:), intent(in) :: a, b
    real(dp),               intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_1D_array

! ===== 2-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

  !> Test if a(:,:) == b
  pure function test_eq_dp_2D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_2D_scalar

  !> Test if a(:,:) /= b
  pure function test_neq_dp_2D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_2D_scalar

  !> Test if a(:,:) > b
  pure function test_gt_dp_2D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_2D_scalar

  !> Test if a(:,:) < b
  pure function test_lt_dp_2D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_2D_scalar

  !> Test if a(:,:) >= b
  pure function test_ge_dp_2D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_2D_scalar

  !> Test if a(:,:) <= b
  pure function test_le_dp_2D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_2D_scalar

  !> Test if a(:,:) >= b1 && a(:,:) <= b2
  pure function test_ge_le_dp_2D_scalar( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_2D_scalar

  !> Test if a(:,:) >= (b - tol) && a(:,:) <= (b + tol)
  pure function test_tol_dp_2D_scalar( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a
    real(dp),                 intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_2D_scalar

    ! ===== Array =====
    ! ==================

  !> Test if a(:,:) == b(:,:)
  pure function test_eq_dp_2D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_2D_array

  !> Test if a(:,:) /= b(:,:)
  pure function test_neq_dp_2D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_2D_array

  !> Test if a(:,:) > b(:,:)
  pure function test_gt_dp_2D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_2D_array

  !> Test if a(:,:) < b(:,:)
  pure function test_lt_dp_2D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_2D_array

  !> Test if a(:,:) >= b(:,:)
  pure function test_ge_dp_2D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_2D_array

  !> Test if a(:,:) <= b(:,:)
  pure function test_le_dp_2D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_2D_array

  !> Test if a(:,:) >= b1(:,:) && a(:,:) <= b2(:,:)
  pure function test_ge_le_dp_2D_array( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_2D_array

  !> Test if a(:,:) >= (b(:,:) - tol) && a(:,:) <= (b(:,:) + tol)
  pure function test_tol_dp_2D_array( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:,:), intent(in) :: a, b
    real(dp),                 intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_2D_array

! ===== 3-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

  !> Test if a(:,:,:) == b
  pure function test_eq_dp_3D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_3D_scalar

  !> Test if a(:,:,:) /= b
  pure function test_neq_dp_3D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_3D_scalar

  !> Test if a(:,:,:) > b
  pure function test_gt_dp_3D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_3D_scalar

  !> Test if a(:,:,:) < b
  pure function test_lt_dp_3D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_3D_scalar

  !> Test if a(:,:,:) >= b
  pure function test_ge_dp_3D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_3D_scalar

  !> Test if a(:,:,:) <= b
  pure function test_le_dp_3D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_3D_scalar

  !> Test if a(:,:,:) >= b1 && a(:,:,:) <= b2
  pure function test_ge_le_dp_3D_scalar( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_3D_scalar

  !> Test if a(:,:,:) >= (b - tol) && a(:,:,:) <= (b + tol)
  pure function test_tol_dp_3D_scalar( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a
    real(dp),                   intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_3D_scalar

      ! ===== Array =====
      ! ==================

  !> Test if a(:,:,:) == b(:,:,:)
  pure function test_eq_dp_3D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_3D_array

  !> Test if a(:,:,:) /= b(:,:,:)
  pure function test_neq_dp_3D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_3D_array

  !> Test if a(:,:,:) > b(:,:,:)
    pure function test_gt_dp_3D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_3D_array

  !> Test if a(:,:,:) < b(:,:,:)
  pure function test_lt_dp_3D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_3D_array

  !> Test if a(:,:,:) >= b(:,:,:)
  pure function test_ge_dp_3D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_3D_array

  !> Test if a(:,:,:) <= b(:,:,:)
  pure function test_le_dp_3D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_3D_array

  !> Test if a(:,:,:) >= b1(:,:,:) && a(:,:,:) <= b2(:,:,:)
  pure function test_ge_le_dp_3D_array( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_3D_array

  !> Test if a(:,:,:) >= (b(:,:,:) - tol) && a(:,:,:) <= (b(:,:,:) + tol)
  pure function test_tol_dp_3D_array( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:), intent(in) :: a, b
    real(dp),                   intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_3D_array

! ===== 4-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

  !> Test if a(:,:,:,:) == b
  pure function test_eq_dp_4D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_4D_scalar

  !> Test if a(:,:,:,:) /= b
  pure function test_neq_dp_4D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_4D_scalar

  !> Test if a(:,:,:,:) > b
  pure function test_gt_dp_4D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_4D_scalar

  !> Test if a(:,:,:,:) < b
  pure function test_lt_dp_4D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_4D_scalar

  !> Test if a(:,:,:,:) >= b
  pure function test_ge_dp_4D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_4D_scalar

  !> Test if a(:,:,:,:) <= b
  pure function test_le_dp_4D_scalar( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_4D_scalar

  !> Test if a(:,:,:,:) >= b1 && a(:,:,:,:) <= b2
  pure function test_ge_le_dp_4D_scalar( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_4D_scalar

  !> Test if a(:,:,:,:) >= (b - tol) && a(:,:,:,:) <= (b + tol)
  pure function test_tol_dp_4D_scalar( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a
    real(dp),                     intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_4D_scalar

    ! ===== Array =====
    ! ==================

  !> Test if a(:,:,:,:) == b(:,:,:,:)
  pure function test_eq_dp_4D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_dp_4D_array

  !> Test if a(:,:,:,:) /= b(:,:,:,:)
  pure function test_neq_dp_4D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_dp_4D_array

  !> Test if a(:,:,:,:) > b(:,:,:,:)
  pure function test_gt_dp_4D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_dp_4D_array

  !> Test if a(:,:,:,:) < b(:,:,:,:)
  pure function test_lt_dp_4D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_dp_4D_array

  !> Test if a(:,:,:,:) >= b(:,:,:,:)
  pure function test_ge_dp_4D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_dp_4D_array

  !> Test if a(:,:,:,:) <= b(:,:,:,:)
  pure function test_le_dp_4D_array( a, b) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_dp_4D_array

  !> Test if a(:,:,:,:) >= b1(:,:,:,:) && a(:,:,:,:) <= b2(:,:,:,:)
  pure function test_ge_le_dp_4D_array( a, b1, b2) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_dp_4D_array

  !> Test if a(:,:,:,:) >= (b(:,:,:,:) - tol) && a(:,:,:,:) <= (b(:,:,:,:) + tol)
  pure function test_tol_dp_4D_array( a, b, tol) result( res)
    ! In/output variables:
    real(dp), dimension(:,:,:,:), intent(in) :: a, b
    real(dp),                     intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_dp_4D_array

end module tests_dp
