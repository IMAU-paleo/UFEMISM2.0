MODULE CSR_sparse_matrix_utilities

  ! Subroutines to work with Compressed Sparse Row formatted matrices

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate
  USE mpi_distributed_memory                                 , ONLY: partition_list

  IMPLICIT NONE

  ! The basic CSR matrix type
  TYPE type_sparse_matrix_CSR_dp
    ! Compressed Sparse Row (CSR) format matrix

    INTEGER                                 :: m,n                           ! A = [m-by-n]
    INTEGER                                 :: m_loc,n_loc                   ! number of rows and columns owned by each process
    INTEGER                                 :: i1,i2,j1,j2                   ! rows and columns owned by each process
    INTEGER                                 :: nnz_max                       ! Maximum number of non-zero entries in A
    INTEGER                                 :: nnz                           ! Actual  number of non-zero entries in A
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: ptr                           ! Row start indices
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: ind                           ! Column indices
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: val                           ! Values

  END TYPE type_sparse_matrix_CSR_dp

CONTAINS

! ===== Subroutinea ======
! ========================

  SUBROUTINE allocate_matrix_CSR_dist( A, m_glob, n_glob, m_loc, n_loc, nnz_max_proc)
    ! Allocate memory for a CSR-format sparse m-by-n matrix A

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: m_glob, n_glob, m_loc, n_loc, nnz_max_proc

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_matrix_CSR_dist'
    INTEGER,  DIMENSION(par%n)                         :: m_loc_all, n_loc_all

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Matrix dimensions
    A%m       = m_glob
    A%n       = n_glob
    A%m_loc   = m_loc
    A%n_loc   = n_loc
    A%nnz_max = nnz_max_proc
    A%nnz     = 0

    ! Partition rows and columns over the processes
    CALL MPI_ALLGATHER( m_loc, 1, MPI_INTEGER, m_loc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER( n_loc, 1, MPI_INTEGER, n_loc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Safety
    IF (SUM( m_loc_all) /= m_glob) CALL crash('sum of numbers of local rows doesnt match number of global rows!')
    IF (SUM( n_loc_all) /= n_glob) CALL crash('sum of numbers of local columns doesnt match number of global columns!')

    A%i1 = 1 + SUM( m_loc_all( 1:par%i  ))
    A%i2 = 1 + SUM( m_loc_all( 1:par%i+1))-1
    A%j1 = 1 + SUM( n_loc_all( 1:par%i  ))
    A%j2 = 1 + SUM( n_loc_all( 1:par%i+1))-1

    ! Allocate memory
    ALLOCATE( A%ptr( A%m+1    ), source = 1    )
    ALLOCATE( A%ind( A%nnz_max), source = 0    )
    ALLOCATE( A%val( A%nnz_max), source = 0._dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_matrix_CSR_dist

  SUBROUTINE deallocate_matrix_CSR_dist( A)
    ! Deallocate memory for a CSR-format sparse matrix A

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: A

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_matrix_CSR_dist'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Matrix dimensions
    A%m       = 0
    A%n       = 0
    A%nnz_max = 0
    A%nnz     = 0

    IF (ALLOCATED( A%ptr)) DEALLOCATE( A%ptr)
    IF (ALLOCATED( A%ind)) DEALLOCATE( A%ind)
    IF (ALLOCATED( A%val)) DEALLOCATE( A%val)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_matrix_CSR_dist

  SUBROUTINE add_entry_CSR_dist( A, i, j, v)
    ! Add value v to row i, column j of CSR-formatted matrix A
    !
    ! NOTE: assumes all rows before i are finished and nothing exists yet for rows after i!

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(IN)    :: v

    ! Safety
    IF (i < A%i1 .OR. i > A%i2) CALL crash('out of ownership range!')

    ! Increase number of non-zeros
    A%nnz = A%nnz + 1

    ! List entry
    A%ind( A%nnz) = j
    A%val( A%nnz) = v

    ! Update pointer list
    A%ptr( i+1 : A%m+1) = A%nnz+1

    ! Extend memory if necessary
    IF (A%nnz > A%nnz_max - 10) CALL extend_matrix_CSR_dist( A, 1000)

  END SUBROUTINE add_entry_CSR_dist

  SUBROUTINE extend_matrix_CSR_dist( A, nnz_extra)
    ! Extend memory for a CSR-format sparse m-by-n matrix A

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: nnz_extra

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_matrix_CSR_dist'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: index_temp
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: val_temp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Extend memory
    CALL reallocate( A%ind, A%nnz + nnz_extra)
    CALL reallocate( A%val, A%nnz + nnz_extra)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extend_matrix_CSR_dist

END MODULE CSR_sparse_matrix_utilities
