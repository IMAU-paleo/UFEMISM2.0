module CSR_sparse_matrix_type

  ! The Compressed Sparse Row matrix type

  use precisions, only: dp

  implicit none

  private

  public :: type_sparse_matrix_CSR_dp

  ! The basic CSR matrix type
  type type_sparse_matrix_CSR_dp
    ! Compressed Sparse Row (CSR) format matrix

    integer                             :: m,n         ! A = [m-by-n]
    integer                             :: m_loc,n_loc ! number of rows and columns owned by each process
    integer                             :: i1,i2,j1,j2 ! rows and columns owned by each process
    integer                             :: nnz_max     ! Maximum number of non-zero entries in A
    integer                             :: nnz         ! Actual  number of non-zero entries in A
    integer,  dimension(:), allocatable :: ptr         ! Row start indices
    integer,  dimension(:), allocatable :: ind         ! Column indices
    real(dp), dimension(:), allocatable :: val         ! Values

  end type type_sparse_matrix_CSR_dp

contains

end module CSR_sparse_matrix_type
