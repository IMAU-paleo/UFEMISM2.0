module unit_tests_petsc

  ! Unit tests for different PETSc routines
  !
  ! Convention: xx = Fortran, x = PETSc

#include <petsc/finclude/petscksp.h>
  use petscksp
  use mpi
  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging                  , only: warning, crash, happy, init_routine, finalise_routine
  use model_configuration                                    , only: C
  use CSR_sparse_matrix_utilities                            , only: type_sparse_matrix_CSR_dp, deallocate_matrix_CSR_dist
  use petsc_basic                                            , only: perr, mat_CSR2petsc, multiply_CSR_matrix_with_vector_1D, multiply_petsc_matrix_with_vector_1D, MatDestroy, &
                                                                     mat_petsc2CSR
  use netcdf_debug                                           , only: write_CSR_matrix_to_NetCDF

  implicit none

contains

  subroutine run_all_petsc_unit_tests
    ! Run all unit tests for the PETSc subroutines

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_petsc_unit_tests'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all unit tests for the PETSc subroutines
    call test_multiply_matrix_with_vector
    call test_matrix_PETSc_CSR_conversion

    ! Add routine to call stack
    call finalise_routine( routine_name)

  end subroutine run_all_petsc_unit_tests

  subroutine test_multiply_matrix_with_vector
    ! Test all multiply_matrix_with_vector routine

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'test_multiply_matrix_with_vector'
    type(type_sparse_matrix_CSR_dp)         :: AA
    type(tMat)                              :: A
    real(dp), dimension(:    ), allocatable :: xx, yy
    logical                                 :: found_errors

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on at least two cores
    if (par%n < 2) call crash('should be run on at least two cores')
    if (par%i > 1) then
      call finalise_routine( routine_name)
      return
    end if

    found_errors = .false.

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
    if (par%i == 0) then

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
      allocate( AA%ptr( AA%m+1))
      allocate( AA%ind( AA%nnz_max))
      allocate( AA%val( AA%nnz_max))

      AA%ptr = [1, 2, 4, 4, 4, 4, 4]
      AA%ind = [1, 1, 2]
      AA%val = [1._dp, 2._dp, 3._dp]

      ! x,y
      allocate( xx( 2), source = 0._dp)
      allocate( yy( 2), source = 0._dp)

      xx = [1._dp, 2._dp]

    elseif (par%i == 1) then

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
      allocate( AA%ptr( AA%m+1))
      allocate( AA%ind( AA%nnz_max))
      allocate( AA%val( AA%nnz_max))

      AA%ptr = [1, 1, 1, 3, 5, 6, 8]
      AA%ind = [2, 4, 4, 5, 5, 5, 6]
      AA%val = [4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp, 10._dp]

      ! x,y
      allocate( xx( 4), source = 0._dp)
      allocate( yy( 4), source = 0._dp)

      xx = [3._dp, 4._dp, 5._dp, 6._dp]

    end if

    ! Perform the multiplication
    call multiply_CSR_matrix_with_vector_1D( AA, xx, yy)

    ! Check results
    if     (par%i == 0) then
      if (ANY( yy /= [1._dp, 8._dp])) found_errors = .true.
    elseif (par%i == 1) then
      if (ANY( yy /= [28._dp, 59._dp, 40._dp, 105._dp])) found_errors = .true.
    end if

  ! == multiply_PETSc_matrix_with_vector_1D

    ! Turn AA into a PETSc matrix
    call mat_CSR2petsc( AA, A)

    ! Perform the multiplication
    call multiply_PETSc_matrix_with_vector_1D( A, xx, yy)

    ! Check results
    if     (par%i == 0) then
      if (ANY( yy /= [1._dp, 8._dp])) found_errors = .true.
    elseif (par%i == 1) then
      if (ANY( yy /= [28._dp, 59._dp, 40._dp, 105._dp])) found_errors = .true.
    end if

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( AA)
    call MatDestroy( A, perr)
    deallocate( xx)
    deallocate( yy)

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.NOT. found_errors) then
      if (par%master) call happy('validated all multiply_matrix_with_vector routines')
    else
      if (par%master) call warning('found errors in multiply_matrix_with_vector routines')
    end if

    ! Add routine to call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_matrix_with_vector

  subroutine test_matrix_PETSc_CSR_conversion
    ! Test matrix conversion between PETSc and CSR formats

    ! Local variables:
    character(len=1024), parameter                                      :: routine_name = 'test_matrix_PETSc_CSR_conversion'
    type(type_sparse_matrix_CSR_dp)                                    :: AA, AA2
    type(tMat)                                                         :: A
    logical                                                            :: found_errors
    integer                                                            :: i,k1,k2,k

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on at least two cores
    if (par%n < 2) call crash('should be run on at least two cores')
    if (par%i > 1) then
      call finalise_routine( routine_name)
      return
    end if

    found_errors = .false.

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
    if     (par%i == 0) then

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
      allocate( AA%ptr( AA%m+1))
      allocate( AA%ind( AA%nnz_max))
      allocate( AA%val( AA%nnz_max))

      AA%ptr = [1, 2, 4, 4, 4, 4, 4]
      AA%ind = [1, 1, 2]
      AA%val = [1._dp, 2._dp, 3._dp]

    elseif (par%i == 1) then

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
      allocate( AA%ptr( AA%m+1))
      allocate( AA%ind( AA%nnz_max))
      allocate( AA%val( AA%nnz_max))

      AA%ptr = [1, 1, 1, 3, 5, 6, 8]
      AA%ind = [2, 4, 4, 5, 5, 5, 6]
      AA%val = [4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp, 10._dp]

    end if

  ! == Convert to PETSc format and back to CSR

    call mat_CSR2petsc( AA, A )
    call mat_petsc2CSR( A , AA2)

    ! Check if everything worked
    if (AA2%m     /= AA%m     .OR. AA2%n     /= AA%n     .OR. &
        AA2%m_loc /= AA%m_loc .OR. AA2%n_loc /= AA%n_loc .OR. &
        AA2%nnz /= AA%nnz) then
      found_errors = .true.
    else
      ! At least the sizes match, now check the entries

      DO i = AA%i1, AA%i2

        if (AA%ptr( i) /= AA2%ptr( i)) then
          found_errors = .true.
          EXIT
        end if

        k1 = AA%ptr( i)
        k2 = AA%ptr( i+1) - 1

        DO k = k1, k2
          if (AA%ind( k) /= AA2%ind( k)) found_errors = .true.
          if (AA%val( k) /= AA2%val( k)) found_errors = .true.
        end DO

      end DO ! DO ii = 1, AA%m_loc

    end if

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    call MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.NOT. found_errors) then
      if (par%master) call happy('validated matrix conversion between PETSc and CSR formats')
    else
      if (par%master) call warning('found errors in matrix conversion between PETSc and CSR formats')
    end if

    ! Add routine to call stack
    call finalise_routine( routine_name)

  end subroutine test_matrix_PETSc_CSR_conversion

end module unit_tests_petsc
