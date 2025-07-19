module tests_int

  ! Basic tests for integer values.

  implicit none

contains

  ! ===== 0-D =====
  ! ===============

  !> Test if a == b
  pure function test_eq_int_0D( a, b) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b
    logical :: res

    res = a == b

  end function test_eq_int_0D

  !> Test if a /= b
  pure function test_neq_int_0D( a, b) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b
    logical :: res

    res = a /= b

  end function test_neq_int_0D

  !> Test if a > b
  pure function test_gt_int_0D( a, b) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b
    logical :: res

    res = a > b

  end function test_gt_int_0D

  !> Test if a < b
  pure function test_lt_int_0D( a, b) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b
    logical :: res

    res = a < b

  end function test_lt_int_0D

  !> Test if a >= b
  pure function test_ge_int_0D( a, b) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b
    logical :: res

    res = a >= b

  end function test_ge_int_0D

  !> Test if a <= b
  pure function test_le_int_0D( a, b) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b
    logical :: res

    res = a <= b

  end function test_le_int_0D

  !> Test if a >= b1 && a <= b2
  pure function test_ge_le_int_0D( a, b1, b2) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b1, b2
    logical :: res

    res = a >= b1 .and. a <= b2

  end function test_ge_le_int_0D

  !> Test if a >= (b - tol) && a <= (b + tol)
  pure function test_tol_int_0D( a, b, tol) result( res)
    integer, intent(in) :: a
    integer, intent(in) :: b, tol
    logical :: res

    res = a >= (b - tol) .and. a <= (b + tol)

  end function test_tol_int_0D

  ! ===== 1-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:) == b
  pure function test_eq_int_1D_scalar( a, b) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_int_1D_scalar

  !> Test if a(:) /= b
  pure function test_neq_int_1D_scalar( a, b) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_1D_scalar

  !> Test if a(:) > b
  pure function test_gt_int_1D_scalar( a, b) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_int_1D_scalar

  !> Test if a(:) < b
  pure function test_lt_int_1D_scalar( a, b) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_int_1D_scalar

  !> Test if a(:) >= b
  pure function test_ge_int_1D_scalar( a, b) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_1D_scalar

  !> Test if a(:) <= b
  pure function test_le_int_1D_scalar( a, b) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_int_1D_scalar

  !> Test if a(:) >= b1 && a(:) <= b2
  pure function test_ge_le_int_1D_scalar( a, b1, b2) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_1D_scalar

  !> Test if a(:) >= (b - tol) && a(:) <= (b + tol)
  pure function test_tol_int_1D_scalar( a, b, tol) result( res)
    integer, dimension(:), intent(in) :: a
    integer,               intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_1D_scalar

  ! ===== Array =====
  ! ==================

  !> Test if a(:) == b(:)
  pure function test_eq_int_1D_array( a, b) result( res)
    integer, dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_int_1D_array

  !> Test if a(:) /= b(:)
  pure function test_neq_int_1D_array( a, b) result( res)
    integer, dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_1D_array

  !> Test if a(:) > b(:)
  pure function test_gt_int_1D_array( a, b) result( res)
    integer, dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_int_1D_array

  !> Test if a(:) < b(:)
  pure function test_lt_int_1D_array( a, b) result( res)
    integer, dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_int_1D_array

  !> Test if a(:) >= b(:)
  pure function test_ge_int_1D_array( a, b) result( res)
    integer, dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_1D_array

  !> Test if a(:) <= b(:)
  pure function test_le_int_1D_array( a, b) result( res)
    integer, dimension(:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_int_1D_array

  !> Test if a(:) >= b1(:) && a(:) <= b2(:)
  pure function test_ge_le_int_1D_array( a, b1, b2) result( res)
    integer, dimension(:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_1D_array

  !> Test if a(:) >= (b(:) - tol) && a(:) <= (b(:) + tol)
  pure function test_tol_int_1D_array( a, b, tol) result( res)
    integer, dimension(:), intent(in) :: a, b
    integer,                        intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_1D_array

  ! ===== 2-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:,:) == b
  pure function test_eq_int_2D_scalar( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_int_2D_scalar

  !> Test if a(:,:) /= b
  pure function test_neq_int_2D_scalar( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_2D_scalar

  !> Test if a(:,:) > b
  pure function test_gt_int_2D_scalar( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_int_2D_scalar

  !> Test if a(:,:) < b
  pure function test_lt_int_2D_scalar( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_int_2D_scalar

  !> Test if a(:,:) >= b
  pure function test_ge_int_2D_scalar( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_2D_scalar

  !> Test if a(:,:) <= b
  pure function test_le_int_2D_scalar( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_int_2D_scalar

  !> Test if a(:,:) >= b1 && a(:,:) <= b2
  pure function test_ge_le_int_2D_scalar( a, b1, b2) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_2D_scalar

  !> Test if a(:,:) >= (b - tol) && a(:,:) <= (b + tol)
  pure function test_tol_int_2D_scalar( a, b, tol) result( res)
    integer, dimension(:,:), intent(in) :: a
    integer,                 intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_2D_scalar

    ! ===== Array =====
    ! ==================

  !> Test if a(:,:) == b(:,:)
  pure function test_eq_int_2D_array( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_int_2D_array

  !> Test if a(:,:) /= b(:,:)
  pure function test_neq_int_2D_array( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_2D_array

  !> Test if a(:,:) > b(:,:)
  pure function test_gt_int_2D_array( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_int_2D_array

  !> Test if a(:,:) < b(:,:)
  pure function test_lt_int_2D_array( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_int_2D_array

  !> Test if a(:,:) >= b(:,:)
  pure function test_ge_int_2D_array( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_2D_array

  !> Test if a(:,:) <= b(:,:)
  pure function test_le_int_2D_array( a, b) result( res)
    integer, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_int_2D_array

  !> Test if a(:,:) >= b1(:,:) && a(:,:) <= b2(:,:)
  pure function test_ge_le_int_2D_array( a, b1, b2) result( res)
    integer, dimension(:,:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_2D_array

  !> Test if a(:,:) >= (b(:,:) - tol) && a(:,:) <= (b(:,:) + tol)
  pure function test_tol_int_2D_array( a, b, tol) result( res)
    integer, dimension(:,:), intent(in) :: a, b
    integer,                 intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_2D_array

! ===== 3-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

  !> Test if a(:,:,:) == b
  pure function test_eq_int_3D_scalar( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_int_3D_scalar

  !> Test if a(:,:,:) /= b
  pure function test_neq_int_3D_scalar( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_3D_scalar

  !> Test if a(:,:,:) > b
  pure function test_gt_int_3D_scalar( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_int_3D_scalar

  !> Test if a(:,:,:) < b
  pure function test_lt_int_3D_scalar( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_int_3D_scalar

  !> Test if a(:,:,:) >= b
  pure function test_ge_int_3D_scalar( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_3D_scalar

  !> Test if a(:,:,:) <= b
  pure function test_le_int_3D_scalar( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_int_3D_scalar

  !> Test if a(:,:,:) >= b1 && a(:,:,:) <= b2
  pure function test_ge_le_int_3D_scalar( a, b1, b2) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_3D_scalar

  !> Test if a(:,:,:) >= (b - tol) && a(:,:,:) <= (b + tol)
  pure function test_tol_int_3D_scalar( a, b, tol) result( res)
    integer, dimension(:,:,:), intent(in) :: a
    integer,                   intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_3D_scalar

    ! ===== Array =====
    ! ==================

  !> Test if a(:,:,:) == b(:,:,:)
  pure function test_eq_int_3D_array( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_int_3D_array

  !> Test if a(:,:,:) /= b(:,:,:)
  pure function test_neq_int_3D_array( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_3D_array

  !> Test if a(:,:,:) > b(:,:,:)
  pure function test_gt_int_3D_array( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_int_3D_array

  !> Test if a(:,:,:) < b(:,:,:)
  pure function test_lt_int_3D_array( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_int_3D_array

  !> Test if a(:,:,:) >= b(:,:,:)
  pure function test_ge_int_3D_array( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_3D_array

  !> Test if a(:,:,:) <= b(:,:,:)
  pure function test_le_int_3D_array( a, b) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_int_3D_array

  !> Test if a(:,:,:) >= b1(:,:,:) && a(:,:,:) <= b2(:,:,:)
  pure function test_ge_le_int_3D_array( a, b1, b2) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_3D_array

  !> Test if a(:,:,:) >= (b(:,:,:) - tol) && a(:,:,:) <= (b(:,:,:) + tol)
  pure function test_tol_int_3D_array( a, b, tol) result( res)
    integer, dimension(:,:,:), intent(in) :: a, b
    integer,                   intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_3D_array

! ===== 4-D =====
! ===============

  ! ===== Scalar =====
  ! ==================

  !> Test if a(:,:,:,:) == b
  pure function test_eq_int_4D_scalar( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b
    logical :: res

    res = all(a == b)

  end function test_eq_int_4D_scalar

  !> Test if a(:,:,:,:) /= b
  pure function test_neq_int_4D_scalar( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_4D_scalar

  !> Test if a(:,:,:,:) > b
  pure function test_gt_int_4D_scalar( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b
    logical :: res

    res = all(a > b)

  end function test_gt_int_4D_scalar

  !> Test if a(:,:,:,:) < b
  pure function test_lt_int_4D_scalar( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b
    logical :: res

    res = all(a < b)

  end function test_lt_int_4D_scalar

  !> Test if a(:,:,:,:) >= b
  pure function test_ge_int_4D_scalar( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_4D_scalar

  !> Test if a(:,:,:,:) <= b
  pure function test_le_int_4D_scalar( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b
    logical :: res

    res = all(a <= b)

  end function test_le_int_4D_scalar

  !> Test if a(:,:,:,:) >= b1 && a(:,:,:,:) <= b2
  pure function test_ge_le_int_4D_scalar( a, b1, b2) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_4D_scalar

  !> Test if a(:,:,:,:) >= (b - tol) && a(:,:,:,:) <= (b + tol)
  pure function test_tol_int_4D_scalar( a, b, tol) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a
    integer,                     intent(in) :: b, tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_4D_scalar

    ! ===== Array =====
    ! ==================

  !> Test if a(:,:,:,:) == b(:,:,:,:)
  pure function test_eq_int_4D_array( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a == b)

  end function test_eq_int_4D_array

  !> Test if a(:,:,:,:) /= b(:,:,:,:)
  pure function test_neq_int_4D_array( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a /= b)

  end function test_neq_int_4D_array

  !> Test if a(:,:,:,:) > b(:,:,:,:)
  pure function test_gt_int_4D_array( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a > b)

  end function test_gt_int_4D_array

  !> Test if a(:,:,:,:) < b(:,:,:,:)
  pure function test_lt_int_4D_array( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a < b)

  end function test_lt_int_4D_array

  !> Test if a(:,:,:,:) >= b(:,:,:,:)
  pure function test_ge_int_4D_array( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a >= b)

  end function test_ge_int_4D_array

  !> Test if a(:,:,:,:) <= b(:,:,:,:)
  pure function test_le_int_4D_array( a, b) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all(a <= b)

  end function test_le_int_4D_array

  !> Test if a(:,:,:,:) >= b1(:,:,:,:) && a(:,:,:,:) <= b2(:,:,:,:)
  pure function test_ge_le_int_4D_array( a, b1, b2) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b1, b2
    logical :: res

    res = all(a >= b1 .and. a <= b2)

  end function test_ge_le_int_4D_array

  !> Test if a(:,:,:,:) >= (b(:,:,:,:) - tol) && a(:,:,:,:) <= (b(:,:,:,:) + tol)
  pure function test_tol_int_4D_array( a, b, tol) result( res)
    integer, dimension(:,:,:,:), intent(in) :: a, b
    integer,                     intent(in) :: tol
    logical :: res

    res = all(a >= (b - tol) .and. a <= (b + tol))

  end function test_tol_int_4D_array

! ===== Permutations =====
! ========================

  !> Test if a(:) == b(:) for any cyclical permutation of a(:)
  !> (used e.g. for comparing triangle-vertex lists between meshes)
  pure function test_eq_permute_int_1D( a, b) result( res)
    ! In/output variables:
    integer, dimension(:), intent(in) :: a, b
    logical :: res
    ! Local variables:
    integer :: i
    integer, dimension(size(a)) :: a_

    res = .false.

    a_ = a
    do i = 1, size( a)
      a_ = [a_(size(a)), a_(1:size(a)-1)]
      res = res .or. all( a_ == b)
    end do

  end function test_eq_permute_int_1D

end module tests_int
