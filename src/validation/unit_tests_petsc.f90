MODULE unit_tests_petsc

  ! Unit tests for different PETSc routines
  !
  ! Convention: xx = Fortran, x = PETSc

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, deallocate_matrix_CSR_dist
  USE petsc_basic                                            , ONLY: perr, mat_CSR2petsc, multiply_CSR_matrix_with_vector_1D, multiply_petsc_matrix_with_vector_1D, MatDestroy, &
                                                                     mat_petsc2CSR
  USE netcdf_debug                                           , ONLY: write_CSR_matrix_to_NetCDF

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_petsc_unit_tests
    ! Run all unit tests for the PETSc subroutines

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_all_petsc_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all unit tests for the PETSc subroutines
    CALL test_multiply_matrix_with_vector
    CALL test_matrix_PETSc_CSR_conversion

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_petsc_unit_tests

  SUBROUTINE test_multiply_matrix_with_vector
    ! Test all multiply_matrix_with_vector routine

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_multiply_matrix_with_vector'
    TYPE(type_sparse_matrix_CSR_dp)                                    :: AA
    TYPE(tMat)                                                         :: A
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                            :: xx, yy
    LOGICAL                                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - should be run on at least two cores
    IF (par%n < 2) CALL crash('should be run on at least two cores')
    IF (par%i > 1) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    found_errors = .FALSE.

  ! == multiply_CSR_matrix_with_vector_1D

    ! Calculate the following matrix-vector multiplication:
    !
    ! [  1,   ,   ,   ,   ,   ] [ 1] = [  1]
    ! [  2,  3,   ,   ,   ,   ] [ 2] = [  8]
    ! [   ,  4,   ,  5,   ,   ] [ 3] = [ 28]
    ! [   ,   ,   ,  6,  7,   ] [ 4] = [ 59]
    ! [   ,   ,   ,   ,  8,   ] [ 5] = [ 40]
    ! [   ,   ,   ,   ,  9, 10] [ 6] = [105]

    ! Let process 0 own rows 1-2, and let process 1 own rows 3-6
    IF     (par%i == 0) THEN

      ! A
      AA%m       = 6
      AA%m_loc   = 2
      AA%i1      = 1
      AA%i2      = 2
      AA%n       = 6
      AA%n_loc   = 2
      AA%j1      = 1
      AA%j2      = 2
      AA%nnz     = 3
      AA%nnz_max = AA%nnz
      ALLOCATE( AA%ptr( AA%m+1))
      ALLOCATE( AA%ind( AA%nnz_max))
      ALLOCATE( AA%val( AA%nnz_max))

      AA%ptr = [1, 2, 4, 4, 4, 4, 4]
      AA%ind = [1, 1, 2]
      AA%val = [1._dp, 2._dp, 3._dp]

      ! x,y
      ALLOCATE( xx( 2), source = 0._dp)
      ALLOCATE( yy( 2), source = 0._dp)

      xx = [1._dp, 2._dp]

    ELSEIF (par%i == 1) THEN

      ! A
      AA%m       = 6
      AA%m_loc   = 4
      AA%i1      = 3
      AA%i2      = 6
      AA%n       = 6
      AA%n_loc   = 4
      AA%j1      = 3
      AA%j2      = 6
      AA%nnz     = 7
      AA%nnz_max = AA%nnz
      ALLOCATE( AA%ptr( AA%m+1))
      ALLOCATE( AA%ind( AA%nnz_max))
      ALLOCATE( AA%val( AA%nnz_max))

      AA%ptr = [1, 1, 1, 3, 5, 6, 8]
      AA%ind = [2, 4, 4, 5, 5, 5, 6]
      AA%val = [4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp, 10._dp]

      ! x,y
      ALLOCATE( xx( 4), source = 0._dp)
      ALLOCATE( yy( 4), source = 0._dp)

      xx = [3._dp, 4._dp, 5._dp, 6._dp]

    END IF

    ! Perform the multiplication
    CALL multiply_CSR_matrix_with_vector_1D( AA, xx, yy)

    ! Check results
    IF     (par%i == 0) THEN
      IF (ANY( yy /= [1._dp, 8._dp])) found_errors = .TRUE.
    ELSEIF (par%i == 1) THEN
      IF (ANY( yy /= [28._dp, 59._dp, 40._dp, 105._dp])) found_errors = .TRUE.
    END IF

  ! == multiply_PETSc_matrix_with_vector_1D

    ! Turn AA into a PETSc matrix
    CALL mat_CSR2petsc( AA, A)

    ! Perform the multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( A, xx, yy)

    ! Check results
    IF     (par%i == 0) THEN
      IF (ANY( yy /= [1._dp, 8._dp])) found_errors = .TRUE.
    ELSEIF (par%i == 1) THEN
      IF (ANY( yy /= [28._dp, 59._dp, 40._dp, 105._dp])) found_errors = .TRUE.
    END IF

    ! Clean up after yourself
    CALL deallocate_matrix_CSR_dist( AA)
    CALL MatDestroy( A, perr)
    DEALLOCATE( xx)
    DEALLOCATE( yy)

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all multiply_matrix_with_vector routines')
    ELSE
      IF (par%master) CALL warning('found errors in multiply_matrix_with_vector routines')
    END IF

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_multiply_matrix_with_vector

  SUBROUTINE test_matrix_PETSc_CSR_conversion
    ! Test matrix conversion between PETSc and CSR formats

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_matrix_PETSc_CSR_conversion'
    TYPE(type_sparse_matrix_CSR_dp)                                    :: AA, AA2
    TYPE(tMat)                                                         :: A
    LOGICAL                                                            :: found_errors
    INTEGER                                                            :: i,k1,k2,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - should be run on at least two cores
    IF (par%n < 2) CALL crash('should be run on at least two cores')
    IF (par%i > 1) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    found_errors = .FALSE.

  ! == Initialise CSR matrix

    ! Set up the following matrix
    !
    ! [  1,   ,   ,   ,   ,   ]
    ! [  2,  3,   ,   ,   ,   ]
    ! [   ,  4,   ,  5,   ,   ]
    ! [   ,   ,   ,  6,  7,   ]
    ! [   ,   ,   ,   ,  8,   ]
    ! [   ,   ,   ,   ,  9, 10]

    ! Let process 0 own rows 1-2, and let process 1 own rows 3-6
    IF     (par%i == 0) THEN

      ! A
      AA%m       = 6
      AA%m_loc   = 2
      AA%i1      = 1
      AA%i2      = 2
      AA%n       = 6
      AA%n_loc   = 2
      AA%j1      = 1
      AA%j2      = 2
      AA%nnz     = 3
      AA%nnz_max = AA%nnz
      ALLOCATE( AA%ptr( AA%m+1))
      ALLOCATE( AA%ind( AA%nnz_max))
      ALLOCATE( AA%val( AA%nnz_max))

      AA%ptr = [1, 2, 4, 4, 4, 4, 4]
      AA%ind = [1, 1, 2]
      AA%val = [1._dp, 2._dp, 3._dp]

    ELSEIF (par%i == 1) THEN

      ! A
      AA%m       = 6
      AA%m_loc   = 4
      AA%i1      = 3
      AA%i2      = 6
      AA%n       = 6
      AA%n_loc   = 4
      AA%j1      = 3
      AA%j2      = 6
      AA%nnz     = 7
      AA%nnz_max = AA%nnz
      ALLOCATE( AA%ptr( AA%m+1))
      ALLOCATE( AA%ind( AA%nnz_max))
      ALLOCATE( AA%val( AA%nnz_max))

      AA%ptr = [1, 1, 1, 3, 5, 6, 8]
      AA%ind = [2, 4, 4, 5, 5, 5, 6]
      AA%val = [4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp, 10._dp]

    END IF

  ! == Convert to PETSc format and back to CSR

    CALL mat_CSR2petsc( AA, A )
    CALL mat_petsc2CSR( A , AA2)

    ! Check if everything worked
    IF (AA2%m     /= AA%m     .OR. AA2%n     /= AA%n     .OR. &
        AA2%m_loc /= AA%m_loc .OR. AA2%n_loc /= AA%n_loc .OR. &
        AA2%nnz /= AA%nnz) THEN
      found_errors = .TRUE.
    ELSE
      ! At least the sizes match, now check the entries

      DO i = AA%i1, AA%i2

        IF (AA%ptr( i) /= AA2%ptr( i)) THEN
          found_errors = .TRUE.
          EXIT
        END IF

        k1 = AA%ptr( i)
        k2 = AA%ptr( i+1) - 1

        DO k = k1, k2
          IF (AA%ind( k) /= AA2%ind( k)) found_errors = .TRUE.
          IF (AA%val( k) /= AA2%val( k)) found_errors = .TRUE.
        END DO

      END DO ! DO ii = 1, AA%m_loc

    END IF

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated matrix conversion between PETSc and CSR formats')
    ELSE
      IF (par%master) CALL warning('found errors in matrix conversion between PETSc and CSR formats')
    END IF

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_matrix_PETSc_CSR_conversion

END MODULE unit_tests_petsc
