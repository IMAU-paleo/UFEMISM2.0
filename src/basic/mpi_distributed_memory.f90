module mpi_distributed_memory

  ! Some routine to work with distributed memory in the MPI parallelised architecture.

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLGATHER, MPI_INTEGER, MPI_GATHERV, MPI_BCAST, &
    MPI_RECV, MPI_ALLGATHERV, MPI_ANY_TAG, MPI_DOUBLE_PRECISION, MPI_SEND, MPI_SCATTERV, &
    MPI_STATUS, MPI_LOGICAL
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine

  implicit none

  private

  public :: partition_list, gather_to_primary, gather_to_all, distribute_from_primary

  interface gather_to_primary
    procedure gather_to_primary_int_1D
    procedure gather_to_primary_int_2D
    procedure gather_to_primary_dp_1D
    procedure gather_to_primary_dp_2D
  end interface gather_to_primary

  interface gather_to_all
    procedure gather_to_all_logical_1D
    procedure gather_to_all_int_1D
    procedure gather_to_all_int_2D
    procedure gather_to_all_dp_1D
    procedure gather_to_all_dp_2D
  end interface gather_to_all

  interface distribute_from_primary
    procedure distribute_from_primary_int_1D
    procedure distribute_from_primary_int_2D
    procedure distribute_from_primary_dp_1D
    procedure distribute_from_primary_dp_2D
  end interface distribute_from_primary

contains

  subroutine partition_list( ntot, i, n, i1, i2)
    !< Partition a list into parallel ranges (e.g. vertex domains)

    ! In/output variables:
    integer, intent(in   ) :: ntot, i, n
    integer, intent(  out) :: i1, i2
    integer                :: remainder, slice_size

    !it must be calculated exactly as petsc does
    !TOdo it would be better to get this ranges from PetsC then ...
    if (ntot > n*2) then
      remainder  = mod( ntot,n)
      slice_size = ntot/n

      i1 = slice_size*i     + min( i  , remainder) + 1
      i2 = slice_size*(i+1) + min( i+1, remainder)
    else
      if (i==0) then
        i1 = 1
        i2 = ntot
      else
        i1 = 1
        i2 = 0
      end if
    end if

  end subroutine partition_list

! ===== Gather distributed variables to the primary =====
! ======================================================

  subroutine gather_to_primary_int_1D( d_partial, d_tot)
    !< Gather a distributed 1-D integer variable to the primary

    ! In/output variables:
    integer , dimension(:),           intent(in   ) :: d_partial
    integer , dimension(:), optional, intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter :: routine_name = 'gather_to_primary_int_1D'
    integer                       :: ierr,n1,i
    integer                       :: n_tot
    integer,  dimension(1:par%n)  :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)

    ! Determine total size of distributed array
    call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    n_tot = sum( counts)
    call MPI_BCAST( n_tot, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    if (par%primary) then
      if (present(d_tot)) then
        if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      else
        call crash('d_tot must be present on primary process')
      endif
    endif

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to the primary
    call MPI_GATHERV( d_partial, n1, MPI_integer, d_tot, counts, displs, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_primary_int_1D

  subroutine gather_to_primary_int_2D( d_partial, d_tot)
    !< Gather a distributed 2-D integer variable to the primary

    ! Input variables:
    integer , dimension(:,:),           intent(in   ) :: d_partial
    integer , dimension(:,:), optional, intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter :: routine_name = 'gather_to_primary_int_2D'
    integer                       :: ierr,n1,n2,i,n2_proc
    integer                       :: j
    integer                       :: dummy(1)
    type(MPI_STATUS)              :: recv_status

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)
    n2 = size( d_partial,2)

    ! Check sizes
    do i = 1, par%n-1
      if (par%i == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do

    if (par%primary) then
      if (.not. present(d_tot)) call crash('d_tot must be present on primary process')
    endif

    do j = 1, n2
      if (par%primary) then
        call gather_to_primary_int_1D( d_partial(:,j), d_tot( :, j))
      else
        call gather_to_primary_int_1D( d_partial(:,j), dummy)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_primary_int_2D

  subroutine gather_to_primary_dp_1D( d_partial, d_tot)
    !< Gather a distributed 1-D dp variable to the primary

    ! Input variables:
    real(dp), dimension(:),           intent(in   ) :: d_partial
    real(dp), dimension(:), optional, intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter :: routine_name = 'gather_to_primary_dp_1D'
    integer                       :: ierr,n1,i
    integer                       :: n_tot
    integer,  dimension(1:par%n)  :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)

    ! Determine total size of distributed array
    call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    n_tot = sum( counts)
    call MPI_BCAST( n_tot, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    if (par%primary) then
      if (present(d_tot)) then
        if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      else
        call crash('d_tot must be present on primary process')
      endif
    endif

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to the primary
    call MPI_GATHERV( d_partial, n1, MPI_DOUBLE_PRECISION, d_tot, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_primary_dp_1D

  subroutine gather_to_primary_dp_2D( d_partial, d_tot)
    !< Gather a distributed 2-D dp variable to the primary

    ! Input variables:
    real(dp), dimension(:,:),           intent(in   ) :: d_partial
    real(dp), dimension(:,:), optional, intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter :: routine_name = 'gather_to_primary_dp_2D'
    integer                       :: ierr,n1,n2,i,n2_proc
    integer                       :: j
    real(dp)                      :: dummy(1)
    type(MPI_STATUS)              :: recv_status

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)
    n2 = size( d_partial,2)

    ! Check sizes
    do i = 1, par%n-1
      if (par%i == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do

    if (par%primary) then
      if (.not. present(d_tot)) call crash('d_tot must be present on primary process')
    endif

    do j = 1, n2
      if (par%primary) then
        call gather_to_primary_dp_1D( d_partial(:,j), d_tot( :, j))
      else
        call gather_to_primary_dp_1D( d_partial(:,j), dummy)
      end if
    end do


    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_primary_dp_2D

! ===== Gather distributed variables to all processes =====
! =========================================================

  subroutine gather_to_all_logical_1D( d_partial, d_tot)
    !< Gather a distributed 1-D logical variable to all processes

    ! Input variables:
    logical, dimension(:), intent(in   ) :: d_partial
    logical, dimension(:), intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter :: routine_name = 'gather_to_all_logical_1D'
    integer                       :: ierr,n1,i
    integer                       :: n_tot
    integer,  dimension(1:par%n)  :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)

    ! Determine total size of distributed array
    call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    n_tot = sum( counts)
    call MPI_BCAST( n_tot, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    if (n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to all processes
    call MPI_ALLGATHERV( d_partial, n1, MPI_logical, d_tot, counts, displs, MPI_logical, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_all_logical_1D

  subroutine gather_to_all_int_1D( d_partial, d_tot)
    !< Gather a distributed 1-D integer variable to all processes

    ! Input variables:
    integer, dimension(:), intent(in   ) :: d_partial
    integer, dimension(:), intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter :: routine_name = 'gather_to_all_int_1D'
    integer                       :: ierr,n1,i
    integer                       :: n_tot
    integer,  dimension(1:par%n)  :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)

    ! Determine total size of distributed array
    call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    n_tot = sum( counts)
    call MPI_BCAST( n_tot, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    if (n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to all processes
    call MPI_ALLGATHERV( d_partial, n1, MPI_integer, d_tot, counts, displs, MPI_integer, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_all_int_1D

  subroutine gather_to_all_int_2D( d_partial, d_tot)
    !< Gather a distributed 2-D integer variable to all processes

    ! Input variables:
    integer, dimension(:,:), intent(in   ) :: d_partial
    integer, dimension(:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter             :: routine_name = 'gather_to_all_int_2D'
    integer                                   :: ierr,n1,n2,i,n2_proc
    integer , dimension(:      ), allocatable :: d_partial_1D
    integer , dimension(:      ), allocatable :: d_tot_1D
    integer                                   :: n1_tot,j
    type(MPI_STATUS)              :: recv_status

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)
    n2 = size( d_partial,2)

    ! Check sizes
    do i = 1, par%n-1
      if (par%i == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do

    ! Gather 1 column at a time
    ALLOCATE( d_partial_1D( n1))
    n1_tot = size( d_tot,1)
    ALLOCATE( d_tot_1D( n1_tot))

    do j = 1, n2
      d_partial_1D = d_partial( :,j)
      call gather_to_all_int_1D( d_partial_1D, d_tot_1D)
      d_tot( :,j) = d_tot_1D
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_all_int_2D

  subroutine gather_to_all_dp_1D( d_partial, d_tot)
    !< Gather a distributed 1-D dp variable to all processes

    ! Input variables:
    real(dp), dimension(:), intent(in   ) :: d_partial
    real(dp), dimension(:), intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter :: routine_name = 'gather_to_all_dp_1D'
    integer                       :: ierr,n1,i
    integer                       :: n_tot
    integer,  dimension(1:par%n)  :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)

    ! Determine total size of distributed array
    call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    n_tot = sum( counts)
    call MPI_BCAST( n_tot, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    if (n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to all processes
    call MPI_ALLGATHERV( d_partial, n1, MPI_DOUBLE_PRECISION, d_tot, counts, displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_all_dp_1D

  subroutine gather_to_all_dp_2D( d_partial, d_tot)
    !< Gather a distributed 2-D dp variable to all processes

    ! Input variables:
    real(dp), dimension(:,:), intent(in   ) :: d_partial
    real(dp), dimension(:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=256), parameter             :: routine_name = 'gather_to_all_dp_2D'
    integer                                   :: ierr,n1,n2,i,n2_proc
    real(dp), dimension(:      ), allocatable :: d_partial_1D
    real(dp), dimension(:      ), allocatable :: d_tot_1D
    integer                                   :: n1_tot,j
    type(MPI_STATUS)              :: recv_status

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)
    n2 = size( d_partial,2)

    ! Check sizes
    do i = 1, par%n-1
      if (par%i == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do

    ! Gather 1 column at a time
    ALLOCATE( d_partial_1D( n1))
    n1_tot = size( d_tot,1)
    ALLOCATE( d_tot_1D( n1_tot))

    do j = 1, n2
      d_partial_1D = d_partial( :,j)
      call gather_to_all_dp_1D( d_partial_1D, d_tot_1D)
      d_tot( :,j) = d_tot_1D
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_to_all_dp_2D

! ===== Distribute variables from the primary =====
! ================================================

  subroutine distribute_from_primary_int_1D( d_tot, d_partial)
    !< Distribute a 1-D integer variable from the primary
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    integer, dimension(:), optional, intent(in   ) :: d_tot
    integer, dimension(:),           intent(  out) :: d_partial

    ! Local variables:
    character(len=256), parameter :: routine_name = 'distribute_from_primary_int_1D'
    integer                       :: ierr,n1,n_tot,i
    integer,  dimension(1:par%n)  :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the partial array owned by this process
    n1 = size( d_partial,1)

    ! Determine total size of distributed array
    call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    n_tot = sum( counts)
    call MPI_BCAST( n_tot, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    if (par%primary) then
      if( .not. present( d_tot)) call crash('d_tot must be present on primary')
      if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
    end if

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Scatter data to all the processes
    call MPI_SCATTERV( d_tot, counts, displs, MPI_integer, d_partial, n1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_from_primary_int_1D

  subroutine distribute_from_primary_int_2D( d_tot, d_partial)
    !< Distribute a 2-D integer variable from the primary
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    integer, dimension(:,:), optional, intent(in   ) :: d_tot
    integer, dimension(:,:),           intent(  out) :: d_partial

    ! Local variables:
    character(len=256), parameter :: routine_name = 'distribute_from_primary_int_2D'
    integer                       :: ierr,n1,n2,i,n2_proc,j
    type(MPI_STATUS)              :: recv_status

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)
    n2 = size( d_partial,2)

    ! Check sizes
    do i = 1, par%n-1
      if (par%i == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do

    do j = 1, n2
      if (par%primary) then
        call distribute_from_primary_int_1D( d_tot( :, j), d_partial( : ,j))
      else
        call distribute_from_primary_int_1D( d_partial=d_partial( : ,j))
      endif
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_from_primary_int_2D

  subroutine distribute_from_primary_dp_1D( d_tot, d_partial)
    !< Distribute a 1-D dp variable from the primary
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    real(dp), dimension(:), optional, intent(in   ) :: d_tot
    real(dp), dimension(:),           intent(  out) :: d_partial

    ! Local variables:
    character(len=256), parameter :: routine_name = 'distribute_from_primary_dp_1D'
    integer                       :: ierr,n1,n_tot,i
    integer,  dimension(1:par%n)  :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the partial array owned by this process
    n1 = size( d_partial,1)

    ! Determine total size of distributed array
    call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    n_tot = sum( counts)
    call MPI_BCAST( n_tot, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Safety
    if (par%primary) then
      if( .not. present( d_tot)) call crash('d_tot must be present on primary')
      if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
    end if

    ! Calculate displacements for MPI_SCATTERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Scatter data to all the processes
    call MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_PRECISION, d_partial, n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_from_primary_dp_1D

  subroutine distribute_from_primary_dp_2D( d_tot, d_partial)
    !< Distribute a 2-D dp variable from the primary
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    real(dp), dimension(:,:), optional, intent(in   ) :: d_tot
    real(dp), dimension(:,:),           intent(  out) :: d_partial

    ! Local variables:
    character(len=256), parameter :: routine_name = 'distribute_from_primary_dp_2D'
    integer                       :: ierr,n1,n2,i,n2_proc,j
    type(MPI_STATUS)              :: recv_status

    ! Add routine to path
    call init_routine( routine_name)

    ! Size of the array owned by this process
    n1 = size( d_partial,1)
    n2 = size( d_partial,2)

    ! Check sizes
    do i = 1, par%n-1
      if (par%i == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do

    ! Distribute 1 column at a time
    do j = 1, n2
      if (par%primary) then
        call distribute_from_primary_dp_1D( d_tot( :, j), d_partial( : ,j))
      else
        call distribute_from_primary_dp_1D( d_partial=d_partial( : ,j))
      endif
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_from_primary_dp_2D

end module mpi_distributed_memory
