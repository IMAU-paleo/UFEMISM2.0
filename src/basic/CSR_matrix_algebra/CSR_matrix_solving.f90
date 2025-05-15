module CSR_matrix_solving

  use assertions_basic, only: assert
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD
  use mpi_distributed_shared_memory, only: allocate_dist_shared, gather_dist_shared_to_all, &
    deallocate_dist_shared, dist_to_hybrid, hybrid_to_dist
  use halo_exchange_mod, only: basic_halo_exchange

  implicit none

  private

  public :: solve_matrix_equation_CSR_Jacobi_wrapper, solve_matrix_equation_CSR_Jacobi

contains

  subroutine solve_matrix_equation_CSR_Jacobi_wrapper( AA, pai_x, xx, pai_b, bb, nit, tol, &
    xx_is_hybrid, bb_is_hybrid, buffer_xx_nih, buffer_bb_nih)
    !< Interface between the old, purely distributed memory architecture,
    !< and the new, hybrid distributed/shared memory architecture.

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),          intent(in   ) :: AA
    type(type_par_arr_info),                  intent(in   ) :: pai_x
    real(dp), dimension(:), target,           intent(inout) :: xx
    type(type_par_arr_info),                  intent(in   ) :: pai_b
    real(dp), dimension(:), target,           intent(in   ) :: bb
    integer,                                  intent(in   ) :: nit
    real(dp),                                 intent(in   ) :: tol
    logical,                        optional, intent(in   ) :: xx_is_hybrid, bb_is_hybrid
    real(dp), dimension(:), target, optional, intent(in   ) :: buffer_xx_nih, buffer_bb_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_Jacobi_wrapper'
    logical                         :: xx_is_hybrid_, bb_is_hybrid_
    real(dp), dimension(:), pointer :: xx_nih => null()
    real(dp), dimension(:), pointer :: bb_nih => null()
    type(MPI_WIN)                   :: wxx_nih, wbb_nih

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
#if (DO_ASSERTIONS)
    call assert( AA%m == AA%n, 'matrix A is not square')
    call assert( pai_x%n == AA%m .and. pai_b%n == AA%n, 'vector sizes dont match')
#endif

    if (present( xx_is_hybrid)) then
      xx_is_hybrid_ = xx_is_hybrid
    else
      xx_is_hybrid_ = .false.
    end if

    if (present( bb_is_hybrid)) then
      bb_is_hybrid_ = bb_is_hybrid
    else
      bb_is_hybrid_ = .false.
    end if

    if (xx_is_hybrid_) then
      xx_nih => xx
    else
      if (present( buffer_xx_nih)) then
        xx_nih => buffer_xx_nih
      else
        call allocate_dist_shared( xx_nih, wxx_nih, pai_x%n_nih)
      end if
      call dist_to_hybrid( pai_x, xx, xx_nih)
      call basic_halo_exchange( pai_x, xx_nih)
    end if

    if (bb_is_hybrid_) then
      bb_nih => bb
    else
      if (present( buffer_bb_nih)) then
        bb_nih => buffer_bb_nih
      else
        call allocate_dist_shared( bb_nih, wbb_nih, pai_b%n_nih)
      end if
      call dist_to_hybrid( pai_b, bb, bb_nih)
    end if

    call solve_matrix_equation_CSR_Jacobi( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol)

    if (xx_is_hybrid_) then
      nullify( xx_nih)
    else
      call hybrid_to_dist( pai_x, xx_nih, xx)
      if (.not. present( buffer_xx_nih)) then
        call deallocate_dist_shared( xx_nih, wxx_nih)
      else
        nullify( xx_nih)
      end if
    end if

    if (bb_is_hybrid_) then
      nullify( bb_nih)
    else
      if (.not. present( buffer_bb_nih)) then
        call deallocate_dist_shared( bb_nih, wbb_nih)
      else
        nullify( bb_nih)
      end if
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_Jacobi_wrapper

  subroutine solve_matrix_equation_CSR_Jacobi( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol)
    !< Solve the matrix equation Ax = b using the Jacobi algorithm

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(inout) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_b
    real(dp), dimension(pai_b%i1_nih:pai_b%i2_nih), intent(in   ) :: bb_nih
    integer,                                        intent(in   ) :: nit
    real(dp),                                       intent(in   ) :: tol

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_Jacobi'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
#if (DO_ASSERTIONS)
    call assert( AA%m == AA%n, 'matrix A is not square')
    call assert( pai_x%n == AA%m .and. pai_b%n == AA%n, 'vector sizes dont match')
#endif

    if (.not. AA%is_finalised) call crash('A is not finalised')

    if (AA%needs_x_tot == 1) then
      call solve_matrix_equation_CSR_Jacobi_x_tot( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol)
    elseif (AA%needs_x_tot == 0) then
      call solve_matrix_equation_CSR_Jacobi_x_nih( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol)
    else
      call crash('needs_x_tot not initialised')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_Jacobi

  subroutine solve_matrix_equation_CSR_Jacobi_x_tot( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(inout) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_b
    real(dp), dimension(pai_b%i1_nih:pai_b%i2_nih), intent(in   ) :: bb_nih
    integer,                                        intent(in   ) :: nit
    real(dp),                                       intent(in   ) :: tol

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_Jacobi_x_tot'
    real(dp), dimension(:), pointer :: AA_diag => null()
    type(MPI_WIN)                   :: wAA_diag
    integer                         :: i,k,j
    real(dp), dimension(:), pointer :: xx_old_tot => null()
    type(MPI_WIN)                   :: wxx_old_tot
    integer                         :: it
    real(dp)                        :: lhs, residual, max_abs_residual
    integer                         :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Save diagonal elements of A
    call allocate_dist_shared( AA_diag, wAA_diag, AA%n_node)
    AA_diag( AA%i1_node:AA%i2_node) => AA_diag

    do i = AA%i1, AA%i2
      do k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%ind( k)
        if (j == i) AA_diag( i) = AA%val( k)
      end do
    end do

    ! Allocate memory for x_old
    call allocate_dist_shared( xx_old_tot, wxx_old_tot, pai_x%n)

    ! Run the Jacobi iteration until it converges
    Jacobi_iterate: do it = 1, nit

      ! Update x_old
      call gather_dist_shared_to_all( pai_x, xx_nih, xx_old_tot)

      max_abs_residual = 0._dp

      do i = AA%i1, AA%i2

        lhs = 0._dp
        do k = AA%ptr( i), AA%ptr( i+1)-1
          j = AA%ind( k)
          lhs = lhs + AA%val( k) * xx_old_tot( j)
        end do

        residual = (lhs - bb_nih( i)) / AA_diag( i)
        xx_nih( i) = xx_old_tot( i) - residual

        max_abs_residual = max( max_abs_residual, abs( residual))

      end do
      call MPI_ALLREDUCE( MPI_IN_PLACE, max_abs_residual, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

      ! ! DENK DROM
      ! if (par%primary) write(0,'(A,I5,A,E12.5)') '    Jacobi iteration ', &
      !   it, ': max_abs_residual = ', max_abs_residual

      ! If the iteration has converged, stop
      if (max_abs_residual < tol) exit Jacobi_iterate

    end do Jacobi_iterate

    ! Clean up after yourself
    call deallocate_dist_shared( AA_diag   , wAA_diag   )
    call deallocate_dist_shared( xx_old_tot, wxx_old_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_Jacobi_x_tot

  subroutine solve_matrix_equation_CSR_Jacobi_x_nih( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(inout) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_b
    real(dp), dimension(pai_b%i1_nih:pai_b%i2_nih), intent(in   ) :: bb_nih
    integer,                                        intent(in   ) :: nit
    real(dp),                                       intent(in   ) :: tol

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_Jacobi_x_nih'
    real(dp), dimension(:), pointer :: AA_diag => null()
    type(MPI_WIN)                   :: wAA_diag
    integer                         :: i,k,j
    real(dp), dimension(:), pointer :: xx_old_nih => null()
    type(MPI_WIN)                   :: wxx_old_nih
    integer                         :: it
    real(dp)                        :: lhs, residual, max_abs_residual
    integer                         :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Save diagonal elements of A
    call allocate_dist_shared( AA_diag, wAA_diag, AA%n_node)
    AA_diag( AA%i1_node:AA%i2_node) => AA_diag

    do i = AA%i1, AA%i2
      do k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%ind( k)
        if (j == i) AA_diag( i) = AA%val( k)
      end do
    end do

    ! Allocate memory for x_old
    call allocate_dist_shared( xx_old_nih, wxx_old_nih, pai_x%n)
    xx_old_nih( pai_x%i1_nih:pai_x%i2_nih) => xx_old_nih

    ! Run the Jacobi iteration until it converges
    Jacobi_iterate: do it = 1, nit

      ! Update x_old
      xx_old_nih( pai_x%i1:pai_x%i2) = xx_nih( pai_x%i1:pai_x%i2)
      call basic_halo_exchange( pai_x, xx_old_nih)

      max_abs_residual = 0._dp

      do i = AA%i1, AA%i2

        lhs = 0._dp
        do k = AA%ptr( i), AA%ptr( i+1)-1
          j = AA%ind( k)
          lhs = lhs + AA%val( k) * xx_old_nih( j)
        end do

        residual = (lhs - bb_nih( i)) / AA_diag( i)
        xx_nih( i) = xx_old_nih( i) - residual

        max_abs_residual = max( max_abs_residual, abs( residual))

      end do
      call MPI_ALLREDUCE( MPI_IN_PLACE, max_abs_residual, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

      ! ! DENK DROM
      ! if (par%primary) write(0,'(A,I5,A,E12.5)') '    Jacobi iteration ', &
      !   it, ': max_abs_residual = ', max_abs_residual

      ! If the iteration has converged, stop
      if (max_abs_residual < tol) exit Jacobi_iterate

    end do Jacobi_iterate

    ! Clean up after yourself
    call deallocate_dist_shared( AA_diag   , wAA_diag   )
    call deallocate_dist_shared( xx_old_nih, wxx_old_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_Jacobi_x_nih

end module CSR_matrix_solving
