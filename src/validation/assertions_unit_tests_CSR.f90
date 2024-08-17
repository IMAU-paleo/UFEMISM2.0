module assertions_unit_tests_CSR

  ! The assertions/unit tests for CSR-type matrices.

  use assertions_unit_tests_basic, only: ASSERTION, UNIT_TEST, process_test_result
  use precisions, only: dp

  implicit none

contains

subroutine test_tol_CSR( a, b, tol, test_mode, message)
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  ! In/output variables:
  type(type_sparse_matrix_CSR_dp), intent(in   ) :: a, b
  real(dp),                        intent(in   ) :: tol
  integer,                         intent(in   ) :: test_mode
  character(len=*),                intent(in   ) :: message
  ! Local variables:
  logical :: test_result

  test_result = .true.

  test_result = test_result .and. a%m     == b%m
  test_result = test_result .and. a%m_loc == b%m_loc
  test_result = test_result .and. a%n     == b%n
  test_result = test_result .and. a%n_loc == b%n_loc
  test_result = test_result .and. a%nnz   == b%nnz

  if (size(a%ptr) == size(b%ptr)) then
    test_result = test_result .and. all(a%ptr == b%ptr)
  else
    test_result = .false.
  end if

  if (size(a%ptr) == size(b%ptr)) then
    test_result = test_result .and. all(a%ind == b%ind)
  else
    test_result = .false.
  end if

  if (size(a%val) == size(b%val)) then
    test_result = test_result .and. all(a%val >= b%val - tol) .and. all(a%val <= b%val + tol)
  else
    test_result = .false.
  end if

  call process_test_result( test_mode, test_result, message)

end subroutine test_tol_CSR

end module assertions_unit_tests_CSR
