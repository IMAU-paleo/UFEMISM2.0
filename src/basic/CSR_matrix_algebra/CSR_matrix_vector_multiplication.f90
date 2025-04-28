module CSR_matrix_vector_multiplication

  use assertions_basic, only: assert
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, MPI_WIN
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_distributed_shared_memory, only: allocate_dist_shared, gather_dist_shared_to_all, &
    deallocate_dist_shared, dist_to_hybrid, hybrid_to_dist
  use halo_exchange_mod, only: basic_halo_exchange

  implicit none

  private

  public :: multiply_CSR_matrix_with_vector_1D, multiply_CSR_matrix_with_vector_2D, &
    multiply_CSR_matrix_with_vector_1D_wrapper, multiply_CSR_matrix_with_vector_2D_wrapper

contains

  subroutine multiply_CSR_matrix_with_vector_1D_wrapper( AA, pai_x, xx, pai_y, yy, &
    xx_is_hybrid, yy_is_hybrid, buffer_xx_nih, buffer_yy_nih)
    !< Interface between the old, purely distributed memory architecture,
    !< and the new, hybrid distributed/shared memory architecture.

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(in   ) :: AA
    type(type_par_arr_info),         intent(in   ) :: pai_x
    real(dp), dimension(:), target,  intent(in   ) :: xx
    type(type_par_arr_info),         intent(in   ) :: pai_y
    real(dp), dimension(:), target,  intent(  out) :: yy
    logical, optional,               intent(in   ) :: xx_is_hybrid, yy_is_hybrid
    real(dp), dimension(:), target, optional, intent(in   ) :: buffer_xx_nih, buffer_yy_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'multiply_CSR_matrix_with_vector_1D_wrapper'
    logical                         :: xx_is_hybrid_, yy_is_hybrid_
    real(dp), dimension(:), pointer :: xx_nih => null()
    real(dp), dimension(:), pointer :: yy_nih => null()
    type(MPI_WIN)                   :: wxx_nih, wyy_nih

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

    if (yy_is_hybrid_) then
      yy_nih => yy
    else
      if (present( buffer_yy_nih)) then
        yy_nih => buffer_yy_nih
      else
        call allocate_dist_shared( yy_nih, wyy_nih, pai_y%n_nih)
      end if
      call dist_to_hybrid( pai_y, yy, yy_nih)
    end if

    call multiply_CSR_matrix_with_vector_1D( AA, pai_x, xx_nih, pai_y, yy_nih)

    if (xx_is_hybrid_) then
      nullify( xx_nih)
    else
      if (.not. present( buffer_xx_nih)) then
        call deallocate_dist_shared( xx_nih, wxx_nih)
      else
        nullify( xx_nih)
      end if
    end if

    if (yy_is_hybrid_) then
      nullify( yy_nih)
    else
      call hybrid_to_dist( pai_y, yy_nih, yy)
      if (.not. present( buffer_yy_nih)) then
        call deallocate_dist_shared( yy_nih, wyy_nih)
      else
        nullify( yy_nih)
      end if
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D_wrapper

  subroutine multiply_CSR_matrix_with_vector_2D_wrapper( AA, pai_x, xx, pai_y, yy, &
    xx_is_hybrid, yy_is_hybrid, buffer_xx_nih, buffer_yy_nih)
    !< Interface between the old, purely distributed memory architecture,
    !< and the new, hybrid distributed/shared memory architecture.

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),  intent(in   ) :: AA
    type(type_par_arr_info),          intent(in   ) :: pai_x
    real(dp), dimension(:,:), target, intent(in   ) :: xx
    type(type_par_arr_info),          intent(in   ) :: pai_y
    real(dp), dimension(:,:), target, intent(  out) :: yy
    logical, optional,                intent(in   ) :: xx_is_hybrid, yy_is_hybrid
    real(dp), dimension(:,:), target, optional, intent(in   ) :: buffer_xx_nih, buffer_yy_nih

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'multiply_CSR_matrix_with_vector_2D_wrapper'
    logical                           :: xx_is_hybrid_, yy_is_hybrid_
    real(dp), dimension(:,:), pointer :: xx_nih => null()
    real(dp), dimension(:,:), pointer :: yy_nih => null()
    type(MPI_WIN)                     :: wxx_nih, wyy_nih

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

    if (xx_is_hybrid_) then
      xx_nih => xx
    else
      if (present( buffer_xx_nih)) then
        xx_nih => buffer_xx_nih
      else
        call allocate_dist_shared( xx_nih, wxx_nih, pai_x%n_nih, size( xx,2))
      end if
      call dist_to_hybrid( pai_x, size( xx,2), xx, xx_nih)
      call basic_halo_exchange( pai_x, size( xx,2), xx_nih)
    end if

    if (yy_is_hybrid_) then
      yy_nih => yy
    else
      if (present( buffer_yy_nih)) then
        yy_nih => buffer_yy_nih
      else
        call allocate_dist_shared( yy_nih, wyy_nih, pai_y%n_nih, size( xx,2))
      end if
      call dist_to_hybrid( pai_y, size( xx,2), yy, yy_nih)
    end if

    call multiply_CSR_matrix_with_vector_2D( AA, pai_x, xx_nih, pai_y, yy_nih, size( xx,2))

    if (xx_is_hybrid_) then
      nullify( xx_nih)
    else
      if (.not. present( buffer_xx_nih)) then
        call deallocate_dist_shared( xx_nih, wxx_nih)
      else
        nullify( xx_nih)
      end if
    end if

    if (yy_is_hybrid_) then
      nullify( yy_nih)
    else
      call hybrid_to_dist( pai_y, size( xx,2), yy_nih, yy)
      if (.not. present( buffer_yy_nih)) then
        call deallocate_dist_shared( yy_nih, wyy_nih)
      else
        nullify( yy_nih)
      end if
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_2D_wrapper

  subroutine multiply_CSR_matrix_with_vector_1D( AA, pai_x, xx_nih, pai_y, yy_nih)
    !< Evaluate the matrix-vector product yy = AA*xx

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(in   ) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_y
    real(dp), dimension(pai_y%i1_nih:pai_y%i2_nih), intent(  out) :: yy_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'multiply_CSR_matrix_with_vector_1D'
    real(dp), dimension(:), pointer :: xx_tot => null()
    type(MPI_WIN)                   :: wxx_tot

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. AA%is_finalised) call crash('A is not finalised')

    if (AA%needs_x_tot == 1) then
      call allocate_dist_shared( xx_tot, wxx_tot, AA%n)
      call gather_dist_shared_to_all( pai_x, xx_nih, xx_tot)
      call multiply_CSR_matrix_with_vector_1D_x_tot( AA, pai_x, xx_tot, pai_y, yy_nih)
      call deallocate_dist_shared( xx_tot, wxx_tot)
    elseif (AA%needs_x_tot == 0) then
      call multiply_CSR_matrix_with_vector_1D_x_nih( AA, pai_x, xx_nih, pai_y, yy_nih)
    else
      call crash('needs_x_tot not initialised')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D

  subroutine multiply_CSR_matrix_with_vector_1D_x_tot( AA, pai_x, xx_tot, pai_y, yy_nih)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(1:pai_x%n),                 intent(in   ) :: xx_tot
    type(type_par_arr_info),                        intent(in   ) :: pai_y
    real(dp), dimension(pai_y%i1_nih:pai_y%i2_nih), intent(  out) :: yy_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'multiply_CSR_matrix_with_vector_1D_x_tot'
    integer                        :: i,k1,k2,k,j

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( &
      (AA%m == pai_y%n) .and. &
      (AA%i1 == pai_y%i1) .and. &
      (AA%i2 == pai_y%i2) .and. &
      (AA%i1_node == pai_y%i1_node) .and. &
      (AA%i2_node == pai_y%i2_node) .and. &
      (AA%n == pai_x%n) .and.  &
      (AA%j1 == pai_x%i1) .and. &
      (AA%j2 == pai_x%i2) .and. &
      (AA%j1_node == pai_x%i1_node) .and. &
      (AA%j2_node == pai_x%i2_node), &
      'size of x doesnt meet expectations of A')
#endif

    ! Perform CSR matrix multiplication
    do i = AA%i1, AA%i2

      yy_nih( i) = 0._dp

      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1)-1

      do k = k1, k2
        j = AA%ind( k)
        yy_nih( i) = yy_nih( i) + AA%val( k) * xx_tot( j)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D_x_tot

  subroutine multiply_CSR_matrix_with_vector_1D_x_nih( AA, pai_x, xx_nih, pai_y, yy_nih)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                intent(in   ) :: AA
    type(type_par_arr_info),                        intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih), intent(in   ) :: xx_nih
    type(type_par_arr_info),                        intent(in   ) :: pai_y
    real(dp), dimension(pai_y%i1_nih:pai_y%i2_nih), intent(  out) :: yy_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'multiply_CSR_matrix_with_vector_1D_x_nih'
    integer                        :: i,k1,k2,k,j

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( &
      (AA%m == pai_y%n) .and. &
      (AA%i1 == pai_y%i1) .and. &
      (AA%i2 == pai_y%i2) .and. &
      (AA%i1_node == pai_y%i1_node) .and. &
      (AA%i2_node == pai_y%i2_node) .and. &
      (AA%n == pai_x%n) .and.  &
      (AA%j1 == pai_x%i1) .and. &
      (AA%j2 == pai_x%i2) .and. &
      (AA%j1_node == pai_x%i1_node) .and. &
      (AA%j2_node == pai_x%i2_node), &
      'size of x doesnt meet expectations of A')
#endif

    ! Perform CSR matrix multiplication
    do i = AA%i1, AA%i2

      yy_nih( i) = 0._dp

      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1)-1

      do k = k1, k2
        j = AA%ind( k)
        yy_nih( i) = yy_nih( i) + AA%val( k) * xx_nih( j)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_1D_x_nih

  subroutine multiply_CSR_matrix_with_vector_2D( AA, pai_x, xx_nih, pai_y, yy_nih, nz)
    !< Evaluate the matrix-vector product yy = AA*xx

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),                             intent(in   ) :: AA
    type(type_par_arr_info),                                     intent(in   ) :: pai_x
    real(dp), dimension(pai_x%i1_nih:pai_x%i2_nih,1:nz), target, intent(in   ) :: xx_nih
    type(type_par_arr_info),                                     intent(in   ) :: pai_y
    real(dp), dimension(pai_y%i1_nih:pai_y%i2_nih,1:nz), target, intent(  out) :: yy_nih
    integer,                                                     intent(in   ) :: nz

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'multiply_CSR_matrix_with_vector_2D'
    real(dp), dimension(:), pointer :: xx_1D_nih, yy_1D_nih
    integer                         :: k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      xx_1D_nih => xx_nih( :,k)
      yy_1D_nih => yy_nih( :,k)
      call multiply_CSR_matrix_with_vector_1D( AA, pai_x, xx_1D_nih, pai_y, yy_1D_nih)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_CSR_matrix_with_vector_2D

end module CSR_matrix_vector_multiplication
