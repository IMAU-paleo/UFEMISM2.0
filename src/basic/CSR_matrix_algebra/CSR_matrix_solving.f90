module CSR_matrix_solving

  use assertions_basic, only: assert
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_WIN
  use mpi_distributed_shared_memory, only: allocate_dist_shared, gather_dist_shared_to_all, &
    deallocate_dist_shared, dist_to_hybrid, hybrid_to_dist
  use halo_exchange_mod, only: basic_halo_exchange

  implicit none

  private

  public :: solve_matrix_equation_CSR_SOR_wrapper

contains

  subroutine solve_matrix_equation_CSR_SOR_wrapper( AA, pai_x, xx, pai_b, bb, nit, tol, omega, &
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
    real(dp),                                 intent(in   ) :: omega
    logical,                        optional, intent(in   ) :: xx_is_hybrid, bb_is_hybrid
    real(dp), dimension(:), target, optional, intent(in   ) :: buffer_xx_nih, buffer_bb_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_SOR_wrapper'
    logical                         :: xx_is_hybrid_, bb_is_hybrid_
    real(dp), dimension(:), pointer :: xx_nih => null()
    real(dp), dimension(:), pointer :: bb_nih => null()
    type(MPI_WIN)                   :: wxx_nih, wbb_nih

    ! Add routine to path
    call init_routine( routine_name)

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

    call solve_matrix_equation_CSR_SOR( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol, omega)

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

  end subroutine solve_matrix_equation_CSR_SOR_wrapper

  subroutine solve_matrix_equation_CSR_SOR( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol, omega)
    !< Solve the matrix equation Ax = b using successive over-relaxation (SOR)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(inout) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_b
    real(dp), dimension(pai_b%i1_nih:pai_b%i2_nih), intent(in   ) :: bb_nih
    integer,                                        intent(in   ) :: nit
    real(dp),                                       intent(in   ) :: tol
    real(dp),                                       intent(in   ) :: omega

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_SOR'

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. AA%is_finalised) call crash('A is not finalised')

    if (AA%needs_x_tot == 1) then
      call solve_matrix_equation_CSR_SOR_x_tot( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol, omega)
    elseif (AA%needs_x_tot == 0) then
      call solve_matrix_equation_CSR_SOR_x_nih( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol, omega)
    else
      call crash('needs_x_tot not initialised')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_SOR

  subroutine solve_matrix_equation_CSR_SOR_x_tot( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol, omega)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(inout) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_b
    real(dp), dimension(pai_b%i1_nih:pai_b%i2_nih), intent(in   ) :: bb_nih
    integer,                                        intent(in   ) :: nit
    real(dp),                                       intent(in   ) :: tol
    real(dp),                                       intent(in   ) :: omega

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_SOR_x_tot'
    real(dp), dimension(:), pointer :: xx_tot => null()
    type(MPI_WIN)                   :: wxx_tot

    ! Add routine to path
    call init_routine( routine_name)

    call warning('boop')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_SOR_x_tot

  subroutine solve_matrix_equation_CSR_SOR_x_nih( AA, pai_x, xx_nih, pai_b, bb_nih, nit, tol, omega)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(inout) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_b
    real(dp), dimension(pai_b%i1_nih:pai_b%i2_nih), intent(in   ) :: bb_nih
    integer,                                        intent(in   ) :: nit
    real(dp),                                       intent(in   ) :: tol
    real(dp),                                       intent(in   ) :: omega

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'solve_matrix_equation_CSR_SOR_x_nih'

    ! Add routine to path
    call init_routine( routine_name)

    call warning('baap')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_SOR_x_nih

  ! omega_dyn = omega

  ! res_max = tol * 2._dp
  ! it = 0
  ! SOR_iterate: DO WHILE (res_max > tol .AND. it < nit)
  !   it = it+1

  !   res_max = 0._dp

  !   DO i = i1, i2

  !     lhs = 0._dp
  !     cij = 0._dp
  !     DO k = CSR%A_ptr( i), CSR%A_ptr( i+1)-1
  !       j = CSR%A_index( k)
  !       lhs = lhs + CSR%A_val( k) * CSR%x( j)
  !       IF (j == i) cij = CSR%A_val( k)
  !     END DO

  !     res = (lhs - CSR%b( i)) / cij
  !     res_max = MAX( res_max, ABS(res))

  !     CSR%x( i) = CSR%x( i) - omega_dyn * res

  !   END DO ! DO i = i1, i2
  !   CALL sync

  !   ! Check if we've reached a stable solution
  !   CALL MPI_ALLREDUCE( MPI_IN_PLACE, res_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  !   !IF (par%master) WRITE(0,*) '      SOR iteration ', it, ': res_max = ', res_max

  !   IF (it > 100 .AND. res_max > 1E3_dp ) THEN

  !     ! Divergence detected - decrease omega, reset solution to zero, restart SOR.
  !     IF (par%master) WRITE(0,*) '  solve_matrix_equation_CSR_SOR - divergence detected; decrease omega, reset solution to zero, restart SOR'
  !     omega_dyn = omega_dyn - 0.1_dp
  !     it = 0
  !     CSR%x( i1:i2) = 0._dp
  !     CALL sync

  !     IF (omega_dyn <= 0.1_dp) THEN
  !       CALL crash('divergence detected even with extremely low relaxation parameter!')
  !     END IF
  !   END IF

  ! END DO SOR_iterate

end module CSR_matrix_solving
