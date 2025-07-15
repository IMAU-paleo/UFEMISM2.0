module CSR_matrix_basics

  ! Subroutines to work with Compressed Sparse Row formatted matrices

  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_SEND, MPI_RECV, &
    MPI_STATUS, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_ALLREDUCE, MPI_MIN, MPI_MAX, MPI_IN_PLACE, &
    MPI_LOGICAL, MPI_LOR
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use parameters
  use reallocate_mod, only: reallocate
  use mpi_distributed_memory, only: partition_list, gather_to_all
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  private

  public :: allocate_matrix_CSR_dist, deallocate_matrix_CSR_dist, duplicate_matrix_CSR_dist, &
    add_entry_CSR_dist, add_empty_row_CSR_dist, extend_matrix_CSR_dist, finalise_matrix_CSR_dist, &
    gather_CSR_dist_to_primary, read_single_row_CSR_dist, allocate_matrix_CSR_loc, &
    set_diagonal_to_one_and_rest_of_row_to_zero, set_row_to_value, set_row_diag_to_val

contains

  ! ===== CSR matrices in distributed memory =====
  ! ==============================================

  subroutine allocate_matrix_CSR_dist( A, m_glob, n_glob, m_loc, n_loc, nnz_max_proc, pai_x, pai_y)
    ! Allocate memory for a CSR-format sparse m-by-n matrix A

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),   intent(inout) :: A
    integer,                           intent(in   ) :: m_glob, n_glob, m_loc, n_loc, nnz_max_proc
    type(type_par_arr_info), optional, intent(in   ) :: pai_x, pai_y

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_matrix_CSR_dist'
    integer                        :: ierr
    integer, dimension(par%n)      :: m_loc_all, n_loc_all

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Matrix dimensions
    A%m       = m_glob
    A%n       = n_glob
    A%m_loc   = m_loc
    A%n_loc   = n_loc
    A%nnz_max = nnz_max_proc
    A%nnz     = 0

    ! Partition rows and columns over the processes
    call MPI_ALLGATHER( m_loc, 1, MPI_integer, m_loc_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_loc, 1, MPI_integer, n_loc_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)

    ! Safety
    if (sum( m_loc_all) /= m_glob) call crash('sum of numbers of local rows doesnt match number of global rows!')
    if (sum( n_loc_all) /= n_glob) call crash('sum of numbers of local columns doesnt match number of global columns!')

    A%i1 = 1 + sum( m_loc_all( 1:par%i  ))
    A%i2 = 1 + sum( m_loc_all( 1:par%i+1))-1
    A%j1 = 1 + sum( n_loc_all( 1:par%i  ))
    A%j2 = 1 + sum( n_loc_all( 1:par%i+1))-1

    ! Range owned by this node
    call MPI_ALLREDUCE( A%i1, A%i1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( A%i2, A%i2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    A%m_node = A%i2_node + 1 - A%i1_node

    call MPI_ALLREDUCE( A%j1, A%j1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( A%j2, A%j2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    A%n_node = A%j2_node + 1 - A%j1_node

    ! Allocate memory
    allocate( A%ptr( A%i1: A%i2+1), source = 1)
    allocate( A%ind( A%nnz_max), source = 0    )
    allocate( A%val( A%nnz_max), source = 0._dp)

    ! Parallel array info
    if (present( pai_x)) then
      A%pai_x = pai_x
    else
      A%pai_x%n       = A%n

      A%pai_x%n_loc   = A%n_loc
      A%pai_x%i1      = A%j1
      A%pai_x%i2      = A%j2

      A%pai_x%n_node  = A%n_node
      A%pai_x%i1_node = A%j1_node
      A%pai_x%i2_node = A%j2_node

      A%pai_x%n_nih   = A%pai_x%n_node
      A%pai_x%i1_nih  = A%pai_x%i1_node
      A%pai_x%i2_nih  = A%pai_x%i2_node

      A%pai_x%n_hle   = 0
      A%pai_x%i1_hle  = 0
      A%pai_x%i2_hle  = -1

      A%pai_x%n_hli   = 0
      A%pai_x%i1_hli  = 0
      A%pai_x%i2_hli  = -1

      A%pai_x%n_hre   = 0
      A%pai_x%i1_hre  = 0
      A%pai_x%i2_hre  = -1

      A%pai_x%n_hri   = 0
      A%pai_x%i1_hri  = 0
      A%pai_x%i2_hri  = -1
    end if

    if (present( pai_y)) then
      A%pai_y = pai_y
    else
      A%pai_y%n       = A%m

      A%pai_y%n_loc   = A%m_loc
      A%pai_y%i1      = A%i1
      A%pai_y%i2      = A%i2

      A%pai_y%n_node  = A%m_node
      A%pai_y%i1_node = A%i1_node
      A%pai_y%i2_node = A%i2_node

      A%pai_y%n_nih   = A%pai_y%n_node
      A%pai_y%i1_nih  = A%pai_y%i1_node
      A%pai_y%i2_nih  = A%pai_y%i2_node

      A%pai_y%n_hle   = 0
      A%pai_y%i1_hle  = 0
      A%pai_y%i2_hle  = -1

      A%pai_y%n_hli   = 0
      A%pai_y%i1_hli  = 0
      A%pai_y%i2_hli  = -1

      A%pai_y%n_hre   = 0
      A%pai_y%i1_hre  = 0
      A%pai_y%i2_hre  = -1

      A%pai_y%n_hri   = 0
      A%pai_y%i1_hri  = 0
      A%pai_y%i2_hri  = -1
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_matrix_CSR_dist

  subroutine deallocate_matrix_CSR_dist( A)
    !< Deallocate memory for a CSR-format sparse matrix A

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(inout) :: A

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_matrix_CSR_dist'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Matrix dimensions
    A%m       = 0
    A%n       = 0
    A%nnz_max = 0
    A%nnz     = 0

    ! Parallelisation ranges
    A%m_loc   = 0
    A%i1      = 0
    A%i2      = 0

    A%n_loc   = 0
    A%j1      = 0
    A%j2      = 0

    A%m_node  = 0
    A%i1_node = 0
    A%i2_node = 0

    A%n_node  = 0
    A%j1_node = 0
    A%j2_node = 0

    if (allocateD( A%ptr)) deallocate( A%ptr)
    if (allocateD( A%ind)) deallocate( A%ind)
    if (allocateD( A%val)) deallocate( A%val)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_matrix_CSR_dist

  subroutine duplicate_matrix_CSR_dist( A, B)
    ! Duplicate the CSR-format sparse matrix A

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(in)    :: A
    type(type_sparse_matrix_CSR_dp),     intent(OUT)   :: B

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'duplicate_matrix_CSR_dist'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Matrix dimensions
    B%m       = A%m
    B%n       = A%n
    B%nnz_max = A%nnz_max
    B%nnz     = A%nnz

    ! Parallelisation ranges
    B%m_loc   = A%m_loc
    B%i1      = A%i1
    B%i2      = A%i2

    B%n_loc   = A%n_loc
    B%j1      = A%j1
    B%j2      = A%j2

    B%m_node  = A%m_node
    B%i1_node = A%i1_node
    B%i2_node = A%i2_node

    B%n_node  = A%n_node
    B%j1_node = A%j1_node
    B%j2_node = A%j2_node

    ! Allocate memory
    allocate( B%ptr( B%i1: B%i2+1    ), source = 1    )
    allocate( B%ind( B%nnz_max), source = 0    )
    allocate( B%val( B%nnz_max), source = 0._dp)

    ! Copy data
    B%ptr = A%ptr
    B%ind = A%ind
    B%val = A%val

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine duplicate_matrix_CSR_dist

  subroutine add_entry_CSR_dist( A, i, j, v)
    ! Add value v to row i, column j of CSR-formatted matrix A
    !
    ! NOTE: assumes all rows before i are finished and nothing exists yet for rows after i!

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(inout) :: A
    integer,                             intent(in)    :: i,j
    real(dp),                            intent(in)    :: v

    ! Safety
    if (i < A%i1 .OR. i > A%i2) call crash('out of ownership range!')

    ! Increase number of non-zeros
    A%nnz = A%nnz + 1

    ! List entry
    A%ind( A%nnz) = j
    A%val( A%nnz) = v

    ! Update pointer list
    A%ptr( i+1) = A%nnz+1

    ! Extend memory if necessary
    if (A%nnz > A%nnz_max - 10) call extend_matrix_CSR_dist( A, 1000)

  end subroutine add_entry_CSR_dist

  subroutine add_empty_row_CSR_dist( A, i)
    ! Add an empty row i to CSR-formatted matrix A
    !
    ! NOTE: assumes all rows before i are finished and nothing exists yet for rows after i!

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(inout) :: A
    integer,                             intent(in)    :: i

    ! Safety
    if (i < A%i1 .OR. i > A%i2) call crash('out of ownership range!')

    ! Update pointer list
    A%ptr( i+1) = A%nnz+1

  end subroutine add_empty_row_CSR_dist

  subroutine extend_matrix_CSR_dist( A, nnz_extra)
    ! Extend memory for a CSR-format sparse m-by-n matrix A

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(inout) :: A
    integer,                             intent(in)    :: nnz_extra

    ! Local variables:

    ! Extend memory
    A%nnz_max = A%nnz + nnz_extra
    call reallocate( A%ind, A%nnz_max)
    call reallocate( A%val, A%nnz_max)

  end subroutine extend_matrix_CSR_dist

  subroutine crop_matrix_CSR_dist( A)
    ! Crop memory for a CSR-format sparse m-by-n matrix A

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(inout) :: A

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'crop_matrix_CSR_dist'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Crop memory
    A%nnz_max = A%nnz
    call reallocate( A%ind, A%nnz_max)
    call reallocate( A%val, A%nnz_max)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine crop_matrix_CSR_dist

  subroutine gather_CSR_dist_to_primary( A, A_tot)
    ! Gather a CSR-format sparse m-by-n matrix A that is distributed over the processes, to the primary

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(in)    :: A
    type(type_sparse_matrix_CSR_dp),     intent(OUT)   :: A_tot

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'gather_CSR_dist_to_primary'
    integer                                            :: ierr
    type(MPI_STATUS)                                   :: recv_status
    integer,  dimension(par%n)                         :: m_glob_all, n_glob_all, m_loc_all, n_loc_all
    integer                                            :: nnz_tot
    integer                                            :: p
    integer                                            :: row, k1, k2, k, col
    real(dp)                                           :: val
    type(type_sparse_matrix_CSR_dp)                    :: A_proc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Gather dimensions
    call MPI_ALLGATHER( A%m    , 1, MPI_integer, m_glob_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( A%n    , 1, MPI_integer, n_glob_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( A%m_loc, 1, MPI_integer, m_loc_all , 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( A%n_loc, 1, MPI_integer, n_loc_all , 1, MPI_integer, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( A%nnz, nnz_tot, 1, MPI_integer, MPI_sum, MPI_COMM_WORLD, ierr)

    ! Safety - check if dimensions match
    if (any( m_glob_all /= A%m)) call crash('global numbers of rows do not match across the processes!')
    if (any( n_glob_all /= A%n)) call crash('global numbers of columns do not match across the processes!')
    if (sum( m_loc_all) /= A%m ) call crash('local numbers of rows do not add up across the processes!')
    if (sum( n_loc_all) /= A%n ) call crash('local numbers of columns do not add up across the processes!')

    ! Allocate memory
    if (par%primary) then
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
      allocate( A_tot%ptr( A%m+1)  , source = 1    )
      allocate( A_tot%ind( nnz_tot), source = 0    )
      allocate( A_tot%val( nnz_tot), source = 0._dp)
    else
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
    end if

    ! Start with the primary's own data
    if (par%primary) then
      do row = A%i1, A%i2
        k1 = A%ptr( row)
        k2 = A%ptr( row+1) - 1
        do k = k1, k2
          col = A%ind( k)
          val = A%val( k)
          call add_entry_CSR_dist( A_tot, row, col, val)
        end do
      end do
    end if

    ! Collect data from the other processes
    do p = 1, par%n-1

      if     (par%i == p) then

        ! Send matrix metadata to primary
        call MPI_Send( A%m      , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%m_loc  , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%i1     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%i2     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%n      , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%n_loc  , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%j1     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%j2     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%nnz    , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%nnz_max, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)

        ! Send matrix data to primary
        call MPI_Send( A%ptr, A%m_loc+1, MPI_integer         , 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%ind, A%nnz_max, MPI_integer         , 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%val, A%nnz_max, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)

      elseif (par%primary) then

        ! Receive matrix metadata from process
        call MPI_RECV( A_proc%m      , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%m_loc  , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%i1     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%i2     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%n      , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%n_loc  , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%j1     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%j2     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%nnz    , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%nnz_max, 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)

        ! Allocate memory
        allocate( A_proc%ptr( A_proc%i1: A_proc%i2+1), source = 0    )
        allocate( A_proc%ind( A_proc%nnz_max), source = 0    )
        allocate( A_proc%val( A_proc%nnz_max), source = 0._dp)

        ! Receive matrix data from process
        call MPI_RECV( A_proc%ptr, A_proc%m_loc+1, MPI_integer         , p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%ind, A_proc%nnz_max, MPI_integer         , p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%val, A_proc%nnz_max, MPI_DOUBLE_PRECISION, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)

        ! Write to total matrix
        do row = A_proc%i1, A_proc%i2
          k1 = A_proc%ptr( row)
          k2 = A_proc%ptr( row+1) - 1
          do k = k1, k2
            col = A_proc%ind( k)
            val = A_proc%val( k)
            call add_entry_CSR_dist( A_tot, row, col, val)
          end do
        end do

        ! Clean up after yourself
        deallocate( A_proc%ptr)
        deallocate( A_proc%ind)
        deallocate( A_proc%val)

      end if ! if     (par%i == p) then
      call sync

    end do ! do p = 1, par%n-1

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_CSR_dist_to_primary

  subroutine read_single_row_CSR_dist( A, i, ind, val, nnz)
    ! Read the coefficients of a single row of A

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(in)    :: A
    integer,                             intent(in)    :: i
    integer,  dimension(:    ),          intent(inout) :: ind
    real(dp), dimension(:    ),          intent(inout) :: val
    integer,                             intent(OUT)   :: nnz

    ! Local variables:
    integer                                            :: k1,k2

    ! Safety
    if (i < A%i1 .OR. i > A%i2) call crash('row {int_01} is not owned by process {int_02}!', int_01 = i, int_02 = par%i)

    k1 = A%ptr( i)
    k2 = A%ptr( i+1) - 1

    nnz = k2 + 1 - k1

    ind( 1:nnz) = A%ind( k1:k2)
    val( 1:nnz) = A%val( k1:k2)

  end subroutine read_single_row_CSR_dist

  subroutine finalise_matrix_CSR_dist( A)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(inout) :: A

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'finalise_matrix_CSR_dist'

    ! Add routine to call stack
    call init_routine( routine_name)

    call crop_matrix_CSR_dist( A)
    call calc_j_node_range( A)

    A%is_finalised = .true.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine finalise_matrix_CSR_dist

  subroutine calc_j_node_range( A)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(inout) :: A

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_j_node_range'
    integer                        :: i, k1, k2, k, j, ierr
    logical                        :: needs_x_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    A%j_min_node =  huge( A%j_min_node)
    A%j_max_node = -huge( A%j_max_node)

    do i = A%i1, A%i2

      k1 = A%ptr( i)
      k2 = A%ptr( i+1)-1

      do k = k1, k2
        j = A%ind( k)
        A%j_min_node = min( A%j_min_node, j)
        A%j_max_node = max( A%j_max_node, j)
      end do

    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, A%j_min_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, A%j_max_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)

    needs_x_tot = A%j_min_node < A%pai_x%i1_nih .or. A%j_max_node > A%pai_x%i2_nih
    call MPI_ALLREDUCE( MPI_IN_PLACE, needs_x_tot, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (needs_x_tot) then
      A%needs_x_tot = 1
    else
      A%needs_x_tot = 0
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_j_node_range

  subroutine set_diagonal_to_one_and_rest_of_row_to_zero( A, i)

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(inout) :: A
    integer,                         intent(in   ) :: i

    call set_row_to_value(    A, i, 0._dp)
    call set_row_diag_to_val( A, i, 1._dp)

  end subroutine set_diagonal_to_one_and_rest_of_row_to_zero

  subroutine set_row_to_value( A, i, val)
    !< Set A(i,:) to val

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(inout) :: A
    integer,                         intent(in   ) :: i
    real(dp),                        intent(in   ) :: val

    ! Local variables:
    integer :: k

    do k = A%ptr( i), A%ptr( i+1) - 1
      A%val( k) = val
    end do

  end subroutine set_row_to_value

  subroutine set_row_diag_to_val( A, i, val)
    !< Set A(i,i) to val

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp), intent(inout) :: A
    integer,                         intent(in   ) :: i
    real(dp),                        intent(in   ) :: val

    ! Local variables:
    integer :: k
    logical :: found_diag

    found_diag = .false.
    do k = A%ptr( i), A%ptr( i+1) - 1
      if (A%ind( k) == i) then
        A%val( k) = val
        found_diag = .true.
      end if
    end do

    ! Safety
    if (.not. found_diag) call crash('set_row_diag_to_val - row has no element on the diagonal')

  end subroutine set_row_diag_to_val

  ! ===== CSR matrices in local memory =====
  ! ========================================

  subroutine allocate_matrix_CSR_loc( A, m, n, nnz_max)
    ! Allocate memory for a CSR-format sparse m-by-n matrix A

    ! In- and output variables:
    type(type_sparse_matrix_CSR_dp),     intent(inout) :: A
    integer,                             intent(in)    :: m, n, nnz_max

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'allocate_matrix_CSR_loc'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Matrix dimensions
    A%m       = m
    A%n       = n
    A%m_loc   = m
    A%n_loc   = n
    A%nnz_max = nnz_max
    A%nnz     = 0

    A%i1 = 1
    A%i2 = m
    A%j1 = 1
    A%j2 = n

    ! Allocate memory
    allocate( A%ptr( A%m+1    ), source = 1    )
    allocate( A%ind( A%nnz_max), source = 0    )
    allocate( A%val( A%nnz_max), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_matrix_CSR_loc

end module CSR_matrix_basics
