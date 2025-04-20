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
    integer                             :: nnz_max     ! Maximum number of non-zero entries in A
    integer                             :: nnz         ! Actual  number of non-zero entries in A
    integer,  dimension(:), allocatable :: ptr         ! Row start indices
    integer,  dimension(:), allocatable :: ind         ! Column indices
    real(dp), dimension(:), allocatable :: val         ! Values

    ! Parallelisation
    integer :: m_loc,  i1,      i2      ! Rows    owned by each process
    integer :: m_node, i1_node, i2_node ! Rows    owned by each node
    integer :: n_loc,  j1,      j2      ! Columns owned by each process
    integer :: n_node, j1_node, j2_node ! Columns owned by each node

    logical :: is_finalised = .false.
    integer :: j_min_node, j_max_node   ! Range of rows of x needed by this node to compute y = A*x

  end type type_sparse_matrix_CSR_dp

contains

end module CSR_sparse_matrix_type
