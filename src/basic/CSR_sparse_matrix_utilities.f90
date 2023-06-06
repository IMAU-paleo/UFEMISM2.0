MODULE CSR_sparse_matrix_utilities

  ! Subroutines to work with Compressed Sparse Row formatted matrices

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
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
    ALLOCATE( A%ptr( A%i1: A%i2+1    ), source = 1    )
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
    A%m_loc   = 0
    A%i1      = 0
    A%i2      = 0
    A%n       = 0
    A%n_loc   = 0
    A%j1      = 0
    A%j2      = 0
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
    A%ptr( i+1) = A%nnz+1

    ! Extend memory if necessary
    IF (A%nnz > A%nnz_max - 10) CALL extend_matrix_CSR_dist( A, 1000)

  END SUBROUTINE add_entry_CSR_dist

  SUBROUTINE add_empty_row_CSR_dist( A, i)
    ! Add an empty row i to CSR-formatted matrix A
    !
    ! NOTE: assumes all rows before i are finished and nothing exists yet for rows after i!

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: i

    ! Safety
    IF (i < A%i1 .OR. i > A%i2) CALL crash('out of ownership range!')

    ! Update pointer list
    A%ptr( i+1) = A%nnz+1

  END SUBROUTINE add_empty_row_CSR_dist

  SUBROUTINE extend_matrix_CSR_dist( A, nnz_extra)
    ! Extend memory for a CSR-format sparse m-by-n matrix A

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: nnz_extra

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_matrix_CSR_dist'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Extend memory
    A%nnz_max = A%nnz + nnz_extra
    CALL reallocate( A%ind, A%nnz_max)
    CALL reallocate( A%val, A%nnz_max)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extend_matrix_CSR_dist

  SUBROUTINE gather_CSR_dist_to_master( A, A_tot)
    ! Gather a CSR-format sparse m-by-n matrix A that is distributed over the processes, to the master

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: A
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(OUT)   :: A_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'gather_CSR_dist_to_master'
    INTEGER,  DIMENSION(par%n)                         :: m_glob_all, n_glob_all, m_loc_all, n_loc_all
    INTEGER                                            :: nnz_tot
    INTEGER                                            :: p
    INTEGER                                            :: row, k1, k2, k, col
    REAL(dp)                                           :: val
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_proc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather dimensions
    CALL MPI_ALLGATHER( A%m    , 1, MPI_INTEGER, m_glob_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER( A%n    , 1, MPI_INTEGER, n_glob_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER( A%m_loc, 1, MPI_INTEGER, m_loc_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER( A%n_loc, 1, MPI_INTEGER, n_loc_all , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    CALL MPI_ALLREDUCE( A%nnz, nnz_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Safety - check if dimensions match
    IF (ANY( m_glob_all /= A%m)) CALL crash('global numbers of rows do not match across the processes!')
    IF (ANY( n_glob_all /= A%n)) CALL crash('global numbers of columns do not match across the processes!')
    IF (SUM( m_loc_all) /= A%m ) CALL crash('local numbers of rows do not add up across the processes!')
    IF (SUM( n_loc_all) /= A%n ) CALL crash('local numbers of columns do not add up across the processes!')

    ! Allocate memory
    IF (par%master) THEN
      A_tot%m       = A%m
      A_tot%m_loc   = A%m
      A_tot%i1      = 1
      A_tot%i2      = A%m
      A_tot%n       = A%n
      A_tot%n_loc   = A%n
      A_tot%j1      = 1
      A_tot%j2      = A%n
      A_tot%nnz     = 0
      A_tot%nnz_max = nnz_tot
      ALLOCATE( A_tot%ptr( A%m+1)  , source = 1    )
      ALLOCATE( A_tot%ind( nnz_tot), source = 0    )
      ALLOCATE( A_tot%val( nnz_tot), source = 0._dp)
    ELSE
      A_tot%m       = A%m
      A_tot%m_loc   = 0
      A_tot%i1      = 1
      A_tot%i2      = 0
      A_tot%n       = A%n
      A_tot%n_loc   = 0
      A_tot%j1      = 1
      A_tot%j2      = 0
      A_tot%nnz     = 0
      A_tot%nnz_max = 0
    END IF

    ! Start with the master's own data
    IF (par%master) THEN
      DO row = A%i1, A%i2
        k1 = A%ptr( row)
        k2 = A%ptr( row+1) - 1
        DO k = k1, k2
          col = A%ind( k)
          val = A%val( k)
          CALL add_entry_CSR_dist( A_tot, row, col, val)
        END DO
      END DO
    END IF

    ! Collect data from the other processes
    DO p = 1, par%n-1

      IF     (par%i == p) THEN

        ! Send matrix metadata to master
        CALL MPI_SEND( A%m      , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%m_loc  , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%i1     , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%i2     , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%n      , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%n_loc  , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%j1     , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%j2     , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%nnz    , 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%nnz_max, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)

        ! Send matrix data to master
        CALL MPI_SEND( A%ptr, A%m_loc+1, MPI_INTEGER         , 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%ind, A%nnz_max, MPI_INTEGER         , 0, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND( A%val, A%nnz_max, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)

      ELSEIF (par%master) THEN

        ! Receive matrix metadata from process
        CALL MPI_RECV( A_proc%m      , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%m_loc  , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%i1     , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%i2     , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%n      , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%n_loc  , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%j1     , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%j2     , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%nnz    , 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%nnz_max, 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)

        ! Allocate memory
        ALLOCATE( A_proc%ptr( A_proc%i1: A_proc%i2+1), source = 0    )
        ALLOCATE( A_proc%ind( A_proc%nnz_max), source = 0    )
        ALLOCATE( A_proc%val( A_proc%nnz_max), source = 0._dp)

        ! Receive matrix data from process
        CALL MPI_RECV( A_proc%ptr, A_proc%m_loc+1, MPI_INTEGER         , p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%ind, A_proc%nnz_max, MPI_INTEGER         , p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)
        CALL MPI_RECV( A_proc%val, A_proc%nnz_max, MPI_DOUBLE_PRECISION, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_status, ierr)

        ! Write to total matrix
        DO row = A_proc%i1, A_proc%i2
          k1 = A_proc%ptr( row)
          k2 = A_proc%ptr( row+1) - 1
          DO k = k1, k2
            col = A_proc%ind( k)
            val = A_proc%val( k)
            CALL add_entry_CSR_dist( A_tot, row, col, val)
          END DO
        END DO

        ! Clean up after yourself
        DEALLOCATE( A_proc%ptr)
        DEALLOCATE( A_proc%ind)
        DEALLOCATE( A_proc%val)

      END IF ! IF     (par%i == p) THEN
      CALL sync

    END DO ! DO p = 1, par%n-1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_CSR_dist_to_master

  SUBROUTINE read_single_row_CSR_dist( A, i, ind, val, nnz)
    ! Read the coefficients of a single row of A

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: A
    INTEGER,                             INTENT(IN)    :: i
    INTEGER,  DIMENSION(:    ),          INTENT(INOUT) :: ind
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: val
    INTEGER,                             INTENT(OUT)   :: nnz

    ! Local variables:
    INTEGER                                            :: k1,k2

    ! Safety
    IF (i < A%i1 .OR. i > A%i2) CALL crash('row {int_01} is not owned by process {int_02}!', int_01 = i, int_02 = par%i)

    k1 = A%ptr( i)
    k2 = A%ptr( i+1) - 1

    nnz = k2 + 1 - k1

    ind( 1:nnz) = A%ind( k1:k2)
    val( 1:nnz) = A%val( k1:k2)

  END SUBROUTINE read_single_row_CSR_dist

END MODULE CSR_sparse_matrix_utilities
