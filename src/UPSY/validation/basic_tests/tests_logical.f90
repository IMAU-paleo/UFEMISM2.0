module tests_logical

  ! Basic tests for logical values.

  implicit none

contains

  ! ===== 0-D =====
  ! ===============

  !> Test if a == b
  pure function test_eqv_logical_0D( a, b) result( res)
    logical, intent(in) :: a
    logical, intent(in) :: b
    logical :: res

    res = a .eqv. b

  end function test_eqv_logical_0D

  !> Test if a /= b
  pure function test_neqv_logical_0D( a, b) result( res)
    logical, intent(in) :: a
    logical, intent(in) :: b
    logical :: res

    res = a .neqv. b

  end function test_neqv_logical_0D

  ! ===== 1-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:) == b
  pure function test_eqv_logical_1D_scalar( a, b) result( res)
    logical, dimension(:), intent(in) :: a
    logical,               intent(in) :: b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_1D_scalar

  !> Test if a(:) /= b
  pure function test_neqv_logical_1D_scalar( a, b) result( res)
    logical, dimension(:), intent(in) :: a
    logical,               intent(in) :: b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_1D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:) == b(:)
  pure function test_eqv_logical_1D_array( a, b) result( res)
    logical, dimension(:), intent(in) :: a, b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_1D_array

  !> Test if a(:) /= b
  pure function test_neqv_logical_1D_array( a, b) result( res)
    logical, dimension(:), intent(in) :: a, b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_1D_array

  ! ===== 2-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:,:) == b
  pure function test_eqv_logical_2D_scalar( a, b) result( res)
    logical, dimension(:,:), intent(in) :: a
    logical,                 intent(in) :: b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_2D_scalar

  !> Test if a(:,:) /= b
  pure function test_neqv_logical_2D_scalar( a, b) result( res)
    logical, dimension(:,:), intent(in) :: a
    logical,                 intent(in) :: b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_2D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:,:) == b(:,:)
  pure function test_eqv_logical_2D_array( a, b) result( res)
    logical, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_2D_array

  !> Test if a(:,:) /= b
  pure function test_neqv_logical_2D_array( a, b) result( res)
    logical, dimension(:,:), intent(in) :: a, b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_2D_array

  ! ===== 3-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:,:,:) == b
  pure function test_eqv_logical_3D_scalar( a, b) result( res)
    logical, dimension(:,:,:), intent(in) :: a
    logical,                   intent(in) :: b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_3D_scalar

  !> Test if a(:,:,:) /= b
  pure function test_neqv_logical_3D_scalar( a, b) result( res)
    logical, dimension(:,:,:), intent(in) :: a
    logical,                   intent(in) :: b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_3D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:,:,:) == b(:,:,:)
  pure function test_eqv_logical_3D_array( a, b) result( res)
    logical, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_3D_array

  !> Test if a(:,:,:) /= b
  pure function test_neqv_logical_3D_array( a, b) result( res)
    logical, dimension(:,:,:), intent(in) :: a, b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_3D_array

  ! ===== 3-D =====
  ! ===============

    ! ===== Scalar =====
    ! ==================

  !> Test if a(:,:,:,:) == b
  pure function test_eqv_logical_4D_scalar( a, b) result( res)
    logical, dimension(:,:,:,:), intent(in) :: a
    logical,                     intent(in) :: b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_4D_scalar

  !> Test if a(:,:,:,:) /= b
  pure function test_neqv_logical_4D_scalar( a, b) result( res)
    logical, dimension(:,:,:,:), intent(in) :: a
    logical,                     intent(in) :: b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_4D_scalar

    ! ===== Array =====
    ! =================

  !> Test if a(:,:,:,:) == b(:,:,:,:)
  pure function test_eqv_logical_4D_array( a, b) result( res)
    logical, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all( a .eqv. b)

  end function test_eqv_logical_4D_array

  !> Test if a(:,:,:,:) /= b
  pure function test_neqv_logical_4D_array( a, b) result( res)
    logical, dimension(:,:,:,:), intent(in) :: a, b
    logical :: res

    res = all( a .neqv. b)

  end function test_neqv_logical_4D_array

end module tests_logical
