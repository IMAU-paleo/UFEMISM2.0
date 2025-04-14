module CSR_matrix_vector_multiplication

  ! Subroutines to work with Compressed Sparse Row formatted matrices

  use assertions_basic, only: assert
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mpi_distributed_memory, only: gather_to_all

  implicit none

  private

  public :: multiply_CSR_matrix_with_vector_1D, multiply_CSR_matrix_with_vector_2D

contains

  subroutine multiply_CSR_matrix_with_vector_1D( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),  intent(in   ) :: AA
    real(dp), dimension(AA%j1:AA%j2), intent(in   ) :: xx
    real(dp), dimension(AA%i1:AA%i2), intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'multiply_CSR_matrix_with_vector_1D'
    real(dp), dimension(:), allocatable :: xx_tot
    integer                             :: i,k1,k2,k,j

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( size(xx,1) == AA%n_loc, 'matrix and vector sizes dont match!')
    call assert( size(yy,1) == AA%m_loc, 'matrix and vector sizes dont match!')
#endif

    allocate( xx_tot( AA%n))
    call gather_to_all( xx, xx_tot)

    ! Perform CSR matrix multiplication
    do i = AA%i1, AA%i2

      yy( i) = 0._dp

      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1)-1

      do k = k1, k2
        j = AA%ind( k)
        yy( i) = yy( i) + AA%val( k) * xx_tot( j)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D

  subroutine multiply_CSR_matrix_with_vector_2D( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA, xx, and yy are stored as distributed memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: AA
    real(dp), dimension(:,:),        intent(in   ) :: xx
    real(dp), dimension(:,:),        intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'multiply_CSR_matrix_with_vector_2D'
    integer                             :: j
    real(dp), dimension(:), allocatable :: xx_1D, yy_1D

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( size(xx,1) == AA%n_loc, 'matrix and vector sizes dont match!')
    call assert( size(yy,1) == AA%m_loc, 'matrix and vector sizes dont match!')
    call assert( size(xx,2) == size(yy,2), 'vector sizes dont match!')
#endif

    ! Allocate memory
    allocate( xx_1D( AA%n_loc), source = 0._dp)
    allocate( yy_1D( AA%m_loc), source = 0._dp)

    ! Calculate each column separately
    do j = 1, size(xx,2)
      xx_1D = xx( :,j)
      call multiply_CSR_matrix_with_vector_1D( AA, xx_1D, yy_1D)
      yy( :,j) = yy_1D
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_2D

end module CSR_matrix_vector_multiplication
