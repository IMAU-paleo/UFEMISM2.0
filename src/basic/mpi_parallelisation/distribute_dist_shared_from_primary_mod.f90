module distribute_dist_shared_from_primary_mod

  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mpi_f08, only: MPI_INTEGER, MPI_LOGICAL, MPI_ALLGATHER, MPI_SCATTERV, MPI_SEND, MPI_RECV, &
    MPI_ANY_TAG, MPI_STATUS, MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX

  implicit none

  private

  public :: distribute_dist_shared_from_primary

  interface distribute_dist_shared_from_primary
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    procedure :: distribute_dist_shared_from_primary_logical_1D
    procedure :: distribute_dist_shared_from_primary_logical_2D
    procedure :: distribute_dist_shared_from_primary_logical_3D
    procedure :: distribute_dist_shared_from_primary_int_1D
    procedure :: distribute_dist_shared_from_primary_int_2D
    procedure :: distribute_dist_shared_from_primary_int_3D
    procedure :: distribute_dist_shared_from_primary_dp_1D
    procedure :: distribute_dist_shared_from_primary_dp_2D
    procedure :: distribute_dist_shared_from_primary_dp_3D
    procedure :: distribute_dist_shared_from_primary_complex_1D
    procedure :: distribute_dist_shared_from_primary_complex_2D
    procedure :: distribute_dist_shared_from_primary_complex_3D
  end interface distribute_dist_shared_from_primary

contains

  subroutine distribute_dist_shared_from_primary_logical_1D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    logical, dimension(:), intent(in   ) :: d_tot
    logical, dimension(:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'distribute_dist_shared_from_primary_logical_1D'
    integer                            :: ierr,n1,n_tot,i
    integer,  dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the partial array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)

      ! Safety
      if (par%primary) then
        if (n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      end if

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data to all the processes
      call MPI_SCATTERV( d_tot, counts, displs, MPI_LOGICAL, d_partial, n1, MPI_LOGICAL, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_1D

  subroutine distribute_dist_shared_from_primary_logical_2D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    logical, dimension(:,:), intent(in   ) :: d_tot
    logical, dimension(:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'distribute_dist_shared_from_primary_logical_2D'
    integer                               :: ierr,n2,i,n2_proc
    integer                               :: j
    type(MPI_STATUS)                      :: recv_status
    logical, dimension(size(d_tot    ,1)) :: d_tot_1D
    logical, dimension(size(d_partial,1)) :: d_partial_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n2, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
        end if
      end do
#endif

      do j = 1, n2
        if (par%primary) d_tot_1D = d_tot(:,j)
        call distribute_dist_shared_from_primary_logical_1D( d_tot_1D, d_partial_1D)
        d_partial(:,j) = d_partial_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_2D

  subroutine distribute_dist_shared_from_primary_logical_3D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    logical, dimension(:,:,:), intent(in   ) :: d_tot
    logical, dimension(:,:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter                           :: routine_name = 'distribute_dist_shared_from_primary_logical_3D'
    integer                                                  :: ierr,n3,i,n3_proc
    integer                                                  :: k
    type(MPI_STATUS)                                         :: recv_status
    logical, dimension(size(d_tot    ,1), size(d_tot    ,2)) :: d_tot_2D
    logical, dimension(size(d_partial,1), size(d_partial,2)) :: d_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n3, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n3_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
        end if
      end do
#endif

      do k = 1, n3
        if (par%primary) d_tot_2D = d_tot(:,:,k)
        call distribute_dist_shared_from_primary_logical_2D( d_tot_2D, d_partial_2D)
        d_partial(:,:,k) = d_partial_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_3D

  subroutine distribute_dist_shared_from_primary_int_1D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    integer, dimension(:), intent(in   ) :: d_tot
    integer, dimension(:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_int_1D'
    integer                           :: ierr,n1,n_tot,i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the partial array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)

      ! Safety
      if (par%primary) then
        if (n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      end if

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data to all the processes
      call MPI_SCATTERV( d_tot, counts, displs, MPI_INTEGER, d_partial, n1, MPI_INTEGER, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_1D

  subroutine distribute_dist_shared_from_primary_int_2D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    integer, dimension(:,:), intent(in   ) :: d_tot
    integer, dimension(:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'distribute_dist_shared_from_primary_int_2D'
    integer                               :: ierr,n2,i,n2_proc
    integer                               :: j
    type(MPI_STATUS)                      :: recv_status
    integer, dimension(size(d_tot    ,1)) :: d_tot_1D
    integer, dimension(size(d_partial,1)) :: d_partial_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
        end if
      end do
#endif

      do j = 1, n2
        if (par%primary) d_tot_1D = d_tot(:,j)
        call distribute_dist_shared_from_primary_int_1D( d_tot_1D, d_partial_1D)
        d_partial(:,j) = d_partial_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_2D

  subroutine distribute_dist_shared_from_primary_int_3D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    integer, dimension(:,:,:), intent(in   ) :: d_tot
    integer, dimension(:,:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter                           :: routine_name = 'distribute_dist_shared_from_primary_int_3D'
    integer                                                  :: ierr,n3,i,n3_proc
    integer                                                  :: k
    type(MPI_STATUS)                                         :: recv_status
    integer, dimension(size(d_tot    ,1), size(d_tot    ,2)) :: d_tot_2D
    integer, dimension(size(d_partial,1), size(d_partial,2)) :: d_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n3, 1, MPI_INTEGER, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
        end if
      end do
#endif

      do k = 1, n3
        if (par%primary) d_tot_2D = d_tot(:,:,k)
        call distribute_dist_shared_from_primary_int_2D( d_tot_2D, d_partial_2D)
        d_partial(:,:,k) = d_partial_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_3D

  subroutine distribute_dist_shared_from_primary_dp_1D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    real(dp), dimension(:), intent(in   ) :: d_tot
    real(dp), dimension(:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_dp_1D'
    integer                           :: ierr,n1,n_tot,i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the partial array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)

      ! Safety
      if (par%primary) then
        if (n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      end if

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data to all the processes
      call MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_PRECISION, d_partial, n1, MPI_DOUBLE_PRECISION, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_1D

  subroutine distribute_dist_shared_from_primary_dp_2D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    real(dp), dimension(:,:), intent(in   ) :: d_tot
    real(dp), dimension(:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'distribute_dist_shared_from_primary_dp_2D'
    integer                               :: ierr,n2,i,n2_proc
    integer                               :: j
    type(MPI_STATUS)                      :: recv_status
    real(dp), dimension(size(d_tot    ,1)) :: d_tot_1D
    real(dp), dimension(size(d_partial,1)) :: d_partial_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
        end if
      end do
#endif

      do j = 1, n2
        if (par%primary) d_tot_1D = d_tot(:,j)
        call distribute_dist_shared_from_primary_dp_1D( d_tot_1D, d_partial_1D)
        d_partial(:,j) = d_partial_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_2D

  subroutine distribute_dist_shared_from_primary_dp_3D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    real(dp), dimension(:,:,:), intent(in   ) :: d_tot
    real(dp), dimension(:,:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter                           :: routine_name = 'distribute_dist_shared_from_primary_dp_3D'
    integer                                                  :: ierr,n3,i,n3_proc
    integer                                                  :: k
    type(MPI_STATUS)                                         :: recv_status
    real(dp), dimension(size(d_tot    ,1), size(d_tot    ,2)) :: d_tot_2D
    real(dp), dimension(size(d_partial,1), size(d_partial,2)) :: d_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n3, 1, MPI_INTEGER, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
        end if
      end do
#endif

      do k = 1, n3
        if (par%primary) d_tot_2D = d_tot(:,:,k)
        call distribute_dist_shared_from_primary_dp_2D( d_tot_2D, d_partial_2D)
        d_partial(:,:,k) = d_partial_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_3D

  subroutine distribute_dist_shared_from_primary_complex_1D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    complex*16, dimension(:), intent(in   ) :: d_tot
    complex*16, dimension(:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_complex_1D'
    integer                           :: ierr,n1,n_tot,i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the partial array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)

      ! Safety
      if (par%primary) then
        if (n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      end if

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data to all the processes
      call MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_COMPLEX, d_partial, n1, MPI_DOUBLE_COMPLEX, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_1D

  subroutine distribute_dist_shared_from_primary_complex_2D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    complex*16, dimension(:,:), intent(in   ) :: d_tot
    complex*16, dimension(:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter           :: routine_name = 'distribute_dist_shared_from_primary_complex_2D'
    integer                                  :: ierr,n2,i,n2_proc
    integer                                  :: j
    type(MPI_STATUS)                         :: recv_status
    complex*16, dimension(size(d_tot    ,1)) :: d_tot_1D
    complex*16, dimension(size(d_partial,1)) :: d_partial_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n2, 1, MPI_INTEGER, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n2_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
        end if
      end do
#endif

      do j = 1, n2
        if (par%primary) d_tot_1D = d_tot(:,j)
        call distribute_dist_shared_from_primary_complex_1D( d_tot_1D, d_partial_1D)
        d_partial(:,j) = d_partial_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_2D

  subroutine distribute_dist_shared_from_primary_complex_3D( d_tot, d_partial)
    !< Distribute a data field from the primary to hybrid distributed/shared memory
    !< (e.g. after reading from to NetCDF)

    ! Input variables:
    complex*16, dimension(:,:,:), intent(in   ) :: d_tot
    complex*16, dimension(:,:,:), intent(  out) :: d_partial

    ! Local variables:
    character(len=1024), parameter                              :: routine_name = 'distribute_dist_shared_from_primary_complex_3D'
    integer                                                     :: ierr,n3,i,n3_proc
    integer                                                     :: k
    type(MPI_STATUS)                                            :: recv_status
    complex*16, dimension(size(d_tot    ,1), size(d_tot    ,2)) :: d_tot_2D
    complex*16, dimension(size(d_partial,1), size(d_partial,2)) :: d_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_partial = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
      ! Check sizes
      do i = 1, par%n_nodes-1
        if (par%node_ID == i) then
          call MPI_SEND( n3, 1, MPI_INTEGER, 0, 0, par%mpi_comm_node_primaries, ierr)
        elseif (par%primary) then
          call MPI_RECV( n3_proc, 1, MPI_INTEGER, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
          if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
        end if
      end do
#endif

      do k = 1, n3
        if (par%primary) d_tot_2D = d_tot(:,:,k)
        call distribute_dist_shared_from_primary_complex_2D( d_tot_2D, d_partial_2D)
        d_partial(:,:,k) = d_partial_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_3D

end module distribute_dist_shared_from_primary_mod
