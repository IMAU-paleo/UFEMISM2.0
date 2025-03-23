module CSR_matrix_vector_multiplication

  ! Subroutines to work with Compressed Sparse Row formatted matrices

  use assertions_basic, only: assert
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_SEND, MPI_RECV, &
    MPI_STATUS, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_ALLREDUCE, MPI_BCAST, MPI_WIN
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use parameters
  use reallocate_mod, only: reallocate
  use mpi_distributed_memory, only: partition_list, gather_to_all
  use mpi_distributed_shared_memory, only: allocate_dist_shared, gather_dist_shared_to_all, &
    deallocate_dist_shared

  implicit none

  private

  public :: multiply_CSR_matrix_with_vector_1D, multiply_CSR_matrix_with_vector_2D

contains

  subroutine multiply_CSR_matrix_with_vector_1D( AA, xx, yy, xx_is_hybrid, yy_is_hybrid)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: AA
    real(dp), dimension(:),          intent(in   ) :: xx
    real(dp), dimension(:),          intent(  out) :: yy
    logical, optional,               intent(in   ) :: xx_is_hybrid, yy_is_hybrid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'multiply_CSR_matrix_with_vector_1D'
    logical                        :: xx_is_hybrid_, yy_is_hybrid_

    ! Add routine to path
    call init_routine( routine_name)

    if (present( xx_is_hybrid)) then
      xx_is_hybrid_ = xx_is_hybrid
    else
      xx_is_hybrid_ = .false.
    end if
    if (present( yy_is_hybrid)) then
      yy_is_hybrid_ = yy_is_hybrid
    else
      yy_is_hybrid_ = .false.
    end if

    if ((.not. xx_is_hybrid_) .and. (.not. yy_is_hybrid_)) then
      call multiply_CSR_matrix_with_vector_1D_dist_dist( AA, xx, yy)
    elseif ((.not. xx_is_hybrid_) .and. yy_is_hybrid_) then
      call multiply_CSR_matrix_with_vector_1D_dist_hybrid( AA, xx, yy)
    elseif (xx_is_hybrid_ .and. (.not. yy_is_hybrid_)) then
      call multiply_CSR_matrix_with_vector_1D_hybrid_dist( AA, xx, yy)
    elseif (xx_is_hybrid_ .and. yy_is_hybrid_) then
      call multiply_CSR_matrix_with_vector_1D_hybrid_hybrid( AA, xx, yy)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D

  subroutine multiply_CSR_matrix_with_vector_1D_dist_dist( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA, xx, and yy are stored as distributed memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),  intent(in   ) :: AA
    real(dp), dimension(AA%j1:AA%j2), intent(in   ) :: xx
    real(dp), dimension(AA%i1:AA%i2), intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'multiply_CSR_matrix_with_vector_1D_dist_dist'
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

  end subroutine multiply_CSR_matrix_with_vector_1D_dist_dist

  subroutine multiply_CSR_matrix_with_vector_1D_dist_hybrid( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA and xx are stored as distributed memory, yy as hybrid distributed/shared memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),            intent(in   ) :: AA
    real(dp), dimension(AA%j1     :AA%j2     ), intent(in   ) :: xx
    real(dp), dimension(AA%i1_node:AA%i2_node), intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'multiply_CSR_matrix_with_vector_1D_dist_hybrid'
    real(dp), dimension(:), allocatable :: xx_tot
    integer                             :: i,k1,k2,k,j

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( size(xx,1) == AA%n_loc, 'matrix and vector sizes dont match!')
    call assert( size(yy,1) == AA%m_node, 'matrix and vector sizes dont match!')
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

  end subroutine multiply_CSR_matrix_with_vector_1D_dist_hybrid

  subroutine multiply_CSR_matrix_with_vector_1D_hybrid_dist( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA and yy are stored as distributed memory, xx as hybrid distributed/shared memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),            intent(in   ) :: AA
    real(dp), dimension(AA%j1_node:AA%j2_node), intent(in   ) :: xx
    real(dp), dimension(AA%i1     :AA%i2     ), intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'multiply_CSR_matrix_with_vector_1D_hybrid_dist'
    real(dp), dimension(:), pointer :: xx_tot => null()
    type(MPI_WIN)                   :: wxx_tot
    integer                         :: i,k1,k2,k,j

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( size(xx,1) == AA%n_node, 'matrix and vector sizes dont match!')
    call assert( size(yy,1) == AA%m_loc, 'matrix and vector sizes dont match!')
#endif

    ! Gather xx to all nodes
    call allocate_dist_shared( xx_tot, wxx_tot, AA%n)
    call gather_dist_shared_to_all( xx, xx_tot)

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

    ! Clean up after yourself
    call deallocate_dist_shared( xx_tot, wxx_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D_hybrid_dist

  subroutine multiply_CSR_matrix_with_vector_1D_hybrid_hybrid( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA is as distributed memory, xx and yy as hybrid distributed/shared memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),            intent(in   ) :: AA
    real(dp), dimension(AA%j1_node:AA%j2_node), intent(in   ) :: xx
    real(dp), dimension(AA%i1_node:AA%i2_node), intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'multiply_CSR_matrix_with_vector_1D_hybrid_hybrid'
    real(dp), dimension(:), pointer :: xx_tot => null()
    type(MPI_WIN)                   :: wxx_tot
    integer                         :: i,k1,k2,k,j

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( size(xx,1) == AA%n_node, 'matrix and vector sizes dont match!')
    call assert( size(yy,1) == AA%m_node, 'matrix and vector sizes dont match!')
#endif

    ! Gather xx to all nodes
    call allocate_dist_shared( xx_tot, wxx_tot, AA%n)
    call gather_dist_shared_to_all( xx, xx_tot)

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

    ! Clean up after yourself
    call deallocate_dist_shared( xx_tot, wxx_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D_hybrid_hybrid

  subroutine multiply_CSR_matrix_with_vector_2D( AA, xx, yy, xx_is_hybrid, yy_is_hybrid)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA, xx, and yy are stored as distributed memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: AA
    real(dp), dimension(:,:),        intent(in   ) :: xx
    real(dp), dimension(:,:),        intent(  out) :: yy
    logical, optional,               intent(in   ) :: xx_is_hybrid, yy_is_hybrid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'multiply_CSR_matrix_with_vector_2D'
    logical                        :: xx_is_hybrid_, yy_is_hybrid_

    ! Add routine to path
    call init_routine( routine_name)

    if (present( xx_is_hybrid)) then
      xx_is_hybrid_ = xx_is_hybrid
    else
      xx_is_hybrid_ = .false.
    end if
    if (present( yy_is_hybrid)) then
      yy_is_hybrid_ = yy_is_hybrid
    else
      yy_is_hybrid_ = .false.
    end if

    if ((.not. xx_is_hybrid_) .and. (.not. yy_is_hybrid_)) then
      call multiply_CSR_matrix_with_vector_2D_dist_dist( AA, xx, yy)
    elseif ((.not. xx_is_hybrid_) .and. yy_is_hybrid_) then
      call crash('xx dist, yy hybrid not implemented yet!')
      ! call multiply_CSR_matrix_with_vector_2D_dist_hybrid( AA, xx, yy)
    elseif (xx_is_hybrid_ .and. (.not. yy_is_hybrid_)) then
      call multiply_CSR_matrix_with_vector_2D_hybrid_dist( AA, xx, yy)
    elseif (xx_is_hybrid_ .and. yy_is_hybrid_) then
      call crash('xx hybrid, yy hybrid not implemented yet!')
      ! call multiply_CSR_matrix_with_vector_2D_hybrid_hybrid( AA, xx, yy)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_2D

  subroutine multiply_CSR_matrix_with_vector_2D_dist_dist( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA, xx, and yy are stored as distributed memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: AA
    real(dp), dimension(:,:),        intent(in   ) :: xx
    real(dp), dimension(:,:),        intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'multiply_CSR_matrix_with_vector_2D_dist_dist'
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

  end subroutine multiply_CSR_matrix_with_vector_2D_dist_dist

  subroutine multiply_CSR_matrix_with_vector_2D_hybrid_dist( AA, xx, yy)
    !< Multiply a CSR matrix with a FORTRAN vector: yy = AA*xx

    ! NOTE: AA and yy are stored as distributed memory, xx as hybrid distributed/shared memory

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: AA
    real(dp), dimension(:,:),        intent(in   ) :: xx
    real(dp), dimension(:,:),        intent(  out) :: yy

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'multiply_CSR_matrix_with_vector_2D_hybrid_dist'
    integer                             :: j
    real(dp), dimension(:), pointer     :: xx_1D => null()
    type(MPI_WIN)                       :: wxx_1D
    real(dp), dimension(:), allocatable :: yy_1D

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( size(xx,1) == AA%n_node, 'matrix and vector sizes dont match!')
    call assert( size(yy,1) == AA%m_loc, 'matrix and vector sizes dont match!')
    call assert( size(xx,2) == size(yy,2), 'vector sizes dont match!')
#endif

    ! Allocate memory
    call allocate_dist_shared( xx_1D, wxx_1D, AA%n_node)
    allocate( yy_1D( AA%m_node), source = 0._dp)

    ! Calculate each column separately
    do j = 1, size(xx,2)
      xx_1D( AA%j1_node:AA%j2_node) = xx( AA%j1_node:AA%j2_node,j)
      call sync
      call multiply_CSR_matrix_with_vector_1D( AA, xx_1D, yy_1D)
      yy( :,j) = yy_1D
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( xx_1D, wxx_1D)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_2D_hybrid_dist

end module CSR_matrix_vector_multiplication
