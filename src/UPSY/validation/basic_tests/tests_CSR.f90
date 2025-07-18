module tests_CSR

  ! Basic tests for CSR-type matrices.

  use precisions, only: dp
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp

  implicit none

contains

  pure function test_tol_CSR( a, b, tol) result( res)
    type(type_sparse_matrix_CSR_dp), intent(in) :: a, b
    real(dp),                        intent(in) :: tol
    logical :: res

    res = .true.

    res = res .and. a%m     == b%m
    res = res .and. a%m_loc == b%m_loc
    res = res .and. a%n     == b%n
    res = res .and. a%n_loc == b%n_loc
    res = res .and. a%nnz   == b%nnz

    if (size(a%ptr) == size(b%ptr)) then
      res = res .and. all(a%ptr == b%ptr)
    else
      res = .false.
    end if

    if (size(a%ptr) == size(b%ptr)) then
      res = res .and. all(a%ind == b%ind)
    else
      res = .false.
    end if

    if (size(a%val) == size(b%val)) then
      res = res .and. all(a%val >= b%val - tol) .and. all(a%val <= b%val + tol)
    else
      res = .false.
    end if

  end function test_tol_CSR

end module tests_CSR
