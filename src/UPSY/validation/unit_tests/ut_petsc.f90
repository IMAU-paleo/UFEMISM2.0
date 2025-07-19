module ut_petsc

  ! Unit tests for different PETSc routines
  !
  ! Convention: xx = Fortran, x = PETSc

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: deallocate_matrix_CSR_dist, finalise_matrix_CSR_dist
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D
  use petsc_basic, only: mat_CSR2petsc, multiply_petsc_matrix_with_vector_1D, mat_petsc2CSR
  use tests_main
  use assertions_basic
  use ut_basic

  implicit none

  private

  public :: unit_tests_petsc_main

contains

  subroutine unit_tests_petsc_main( test_name_parent)
    ! Run all unit tests for the PETSc subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_petsc_main'
    character(len=1024), parameter :: test_name_local = 'petsc'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Run all unit tests for the PETSc subroutines
    call test_multiply_PETSc_matrix_with_vector_1D( test_name)
    call test_matrix_PETSc_CSR_conversion( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_petsc_main

  subroutine test_multiply_PETSc_matrix_with_vector_1D( test_name_parent)
    ! Test mmultiply_PETSc_matrix_with_vector_1D

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'test_multiply_PETSc_matrix_with_vector_1D'
    character(len=1024), parameter          :: test_name_local = 'multiply_PETSc_matrix_with_vector_1D'
    character(len=1024)                     :: test_name
    type(type_sparse_matrix_CSR_dp)         :: AA
    type(tMat)                              :: A
    real(dp), dimension(:    ), allocatable :: xx, yy, yy_correct
    integer                                 :: perr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! == multiply_CSR_matrix_with_vector_1D

    ! Perform the following matrix-vector multiplication:
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
      allocate( yy_correct( 2), source = 0._dp)

      xx = [1._dp, 2._dp]
      yy_correct = [1._dp, 8._dp]

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
      allocate( yy_correct( 4), source = 0._dp)

      xx = [3._dp, 4._dp, 5._dp, 6._dp]
      yy_correct = [28._dp, 59._dp, 40._dp, 105._dp]

    end if

    call finalise_matrix_CSR_dist( AA)

    ! Turn AA into a PETSc matrix
    call mat_CSR2petsc( AA, A)

    ! Perform the multiplication
    call multiply_PETSc_matrix_with_vector_1D( A, xx, yy)

    ! Check results
    call unit_test( test_eq( yy, yy_correct), test_name)

    ! Clean up after yourself
    call deallocate_matrix_CSR_dist( AA)
    call MatDestroy( A, perr)
    deallocate( xx)
    deallocate( yy)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_multiply_PETSc_matrix_with_vector_1D

  subroutine test_matrix_PETSc_CSR_conversion( test_name_parent)
    ! Test matrix conversion between PETSc and CSR formats

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_matrix_PETSc_CSR_conversion'
    character(len=1024), parameter  :: test_name_local = 'matrix_PETSc_CSR_conversion'
    character(len=1024)             :: test_name
    type(type_sparse_matrix_CSR_dp) :: AA, AA2
    type(tMat)                      :: A
    logical                         :: found_errors

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

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
      allocate( AA%ptr( AA%i1: AA%i2+1), source = 1)
      allocate( AA%ind( AA%nnz_max))
      allocate( AA%val( AA%nnz_max))

      AA%ptr = [1, 2, 4]
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
      allocate( AA%ptr( AA%i1: AA%i2+1), source = 1)
      allocate( AA%ind( AA%nnz_max))
      allocate( AA%val( AA%nnz_max))

      AA%ptr = [1, 3, 5, 6, 8]
      AA%ind = [2, 4, 4, 5, 5, 5, 6]
      AA%val = [4._dp, 5._dp, 6._dp, 7._dp, 8._dp, 9._dp, 10._dp]

    end if

    call finalise_matrix_CSR_dist( AA)

    ! == Convert to PETSc format and back to CSR

    call mat_CSR2petsc( AA, A )
    call mat_petsc2CSR( A , AA2)

    ! Check results
    call unit_test( test_tol ( AA, AA2, 1e-12_dp), test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_matrix_PETSc_CSR_conversion

end module ut_petsc
