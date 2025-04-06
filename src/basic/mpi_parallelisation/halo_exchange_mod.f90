module halo_exchange_mod

  use assertions_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mpi_basic, only: par, sync
  use mpi_f08, only: MPI_ISEND, MPI_IRECV, MPI_INTEGER, MPI_REQUEST, MPI_WAITALL, MPI_STATUSES_IGNORE, &
    MPI_WAIT, MPI_STATUS_IGNORE, MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX

  implicit none

  private

  public :: exchange_halos

  interface exchange_halos
    procedure :: exchange_halos_int_1D
    procedure :: exchange_halos_int_2D
    procedure :: exchange_halos_int_3D
    procedure :: exchange_halos_dp_1D
    procedure :: exchange_halos_dp_2D
    procedure :: exchange_halos_dp_3D
    procedure :: exchange_halos_complex_1D
    procedure :: exchange_halos_complex_2D
    procedure :: exchange_halos_complex_3D
  end interface exchange_halos

contains

  subroutine exchange_halos_int_1D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)

    ! In/output variables:
    integer, dimension(i1_nih:i2_nih), target, intent(inout) :: d_nih
    integer,                                   intent(in   ) :: i1_nih, i2_nih
    integer,                                   intent(in   ) :: i1_hle, i2_hle
    integer,                                   intent(in   ) :: i1_hli, i2_hli
    integer,                                   intent(in   ) :: i1_hre, i2_hre
    integer,                                   intent(in   ) :: i1_hri, i2_hri

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'exchange_halos_int_1D'
    integer                         :: node_ID_left, node_ID_right
    integer                         :: n_hle, n_hli, n_hre, n_hri
    integer, dimension(:), pointer  :: d_hle, d_hli, d_hre, d_hri
    integer, dimension(2)           :: range_send, range_recv
    type(MPI_REQUEST)               :: req
    type(MPI_REQUEST), dimension(2) :: reqs
    integer                         :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Left halos
    if (par%node_ID == 0) then
      ! There is no node to the left
    else
      node_ID_left = par%node_ID - 1
      d_hle( i1_hle:i2_hle) => d_nih( i1_hle:i2_hle)
      d_hli( i1_hli:i2_hli) => d_nih( i1_hli:i2_hli)
      n_hle = i2_hle + 1 - i1_hle
      n_hli = i2_hli + 1 - i1_hli
    end if

    ! Right halos
    if (par%node_ID == par%n_nodes-1) then
      ! There is no node to the right
    else
      node_ID_right = par%node_ID + 1
      d_hre( i1_hre:i2_hre) => d_nih( i1_hre:i2_hre)
      d_hri( i1_hri:i2_hri) => d_nih( i1_hri:i2_hri)
      n_hre = i2_hre + 1 - i1_hre
      n_hri = i2_hri + 1 - i1_hri
    end if

#if (DO_ASSERTIONS)
    ! Safety: check that sending/receiving nodes agree on where the halos are

    ! Send halo ranges from left to right
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      range_send = [i1_hri, i2_hri]
      call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_right, &
        14, par%mpi_comm_node_primaries, req, ierr)
    end if
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_left, &
        14, par%mpi_comm_node_primaries, req, ierr)
      call MPI_WAIT( req, MPI_STATUS_IGNORE)
      call assert( (range_recv(1) == i1_hle .and. range_recv(2) == i2_hle), &
        'unmatched sending/receiving halo sizes')
    end if

    ! Send halo ranges from right to left
    if (par%node_ID > 0 .and. par%node_primary) then
      range_send = [i1_hli, i2_hli]
      call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_left, &
        15, par%mpi_comm_node_primaries, req, ierr)
    end if
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_right, &
        15, par%mpi_comm_node_primaries, req, ierr)
      call MPI_WAIT( req, MPI_STATUS_IGNORE)
      call assert( (range_recv(1) == i1_hre .and. range_recv(2) == i2_hre), &
        'unmatched sending/receiving halo sizes')
    end if

#endif

    ! Send halos from left to right
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_ISEND( d_hri, n_hri, MPI_INTEGER, node_ID_right, &
        13, par%mpi_comm_node_primaries, reqs(1), ierr)
    end if
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_IRECV( d_hle, n_hle, MPI_INTEGER, node_ID_left, &
        13, par%mpi_comm_node_primaries, reqs(1), ierr)
    end if

    ! Send halos from right to left
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_ISEND( d_hli, n_hli, MPI_INTEGER, node_ID_left, &
        37, par%mpi_comm_node_primaries, reqs(2), ierr)
    end if
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_IRECV( d_hre, n_hre, MPI_INTEGER, node_ID_right, &
        37, par%mpi_comm_node_primaries, reqs(2), ierr)
    end if

    if (par%node_primary) then
      call MPI_WAITALL( 2, reqs, MPI_STATUSES_IGNORE)
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_int_1D

  subroutine exchange_halos_int_2D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, n2)

    ! In/output variables:
    integer, dimension(i1_nih:i2_nih,n2), target, intent(inout) :: d_nih
    integer,                                      intent(in   ) :: i1_nih, i2_nih
    integer,                                      intent(in   ) :: i1_hle, i2_hle
    integer,                                      intent(in   ) :: i1_hli, i2_hli
    integer,                                      intent(in   ) :: i1_hre, i2_hre
    integer,                                      intent(in   ) :: i1_hri, i2_hri
    integer,                                      intent(in   ) :: n2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_int_2D'
    integer                        :: k
    integer, dimension(:), pointer :: d_nih_1D

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, n2
      d_nih_1D( i1_nih:i2_nih) => d_nih( i1_nih:i2_nih,k)
      call exchange_halos_int_1D( d_nih_1D, i1_nih, i2_nih, &
        i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)
      nullify( d_nih_1D)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_int_2D

  subroutine exchange_halos_int_3D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, n2, n3)

    ! In/output variables:
    integer, dimension(i1_nih:i2_nih,n2,n3), target, intent(inout) :: d_nih
    integer,                                         intent(in   ) :: i1_nih, i2_nih
    integer,                                         intent(in   ) :: i1_hle, i2_hle
    integer,                                         intent(in   ) :: i1_hli, i2_hli
    integer,                                         intent(in   ) :: i1_hre, i2_hre
    integer,                                         intent(in   ) :: i1_hri, i2_hri
    integer,                                         intent(in   ) :: n2, n3

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_int_3D'
    integer                        :: k,l
    integer, dimension(:), pointer :: d_nih_1D

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, n2
      do l = 1, n3
        d_nih_1D( i1_nih:i2_nih) => d_nih( i1_nih:i2_nih,k,l)
        call exchange_halos_int_1D( d_nih_1D, i1_nih, i2_nih, &
          i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)
        nullify( d_nih_1D)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_int_3D

  subroutine exchange_halos_dp_1D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)

    ! In/output variables:
    real(dp), dimension(i1_nih:i2_nih), target, intent(inout) :: d_nih
    integer,                                    intent(in   ) :: i1_nih, i2_nih
    integer,                                    intent(in   ) :: i1_hle, i2_hle
    integer,                                    intent(in   ) :: i1_hli, i2_hli
    integer,                                    intent(in   ) :: i1_hre, i2_hre
    integer,                                    intent(in   ) :: i1_hri, i2_hri

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'exchange_halos_dp_1D'
    integer                         :: node_ID_left, node_ID_right
    integer                         :: n_hle, n_hli, n_hre, n_hri
    real(dp), dimension(:), pointer :: d_hle, d_hli, d_hre, d_hri
    integer, dimension(2)           :: range_send, range_recv
    type(MPI_REQUEST)               :: req
    type(MPI_REQUEST), dimension(2) :: reqs
    integer                         :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Left halos
    if (par%node_ID == 0) then
      ! There is no node to the left
    else
      node_ID_left = par%node_ID - 1
      d_hle( i1_hle:i2_hle) => d_nih( i1_hle:i2_hle)
      d_hli( i1_hli:i2_hli) => d_nih( i1_hli:i2_hli)
      n_hle = i2_hle + 1 - i1_hle
      n_hli = i2_hli + 1 - i1_hli
    end if

    ! Right halos
    if (par%node_ID == par%n_nodes-1) then
      ! There is no node to the right
    else
      node_ID_right = par%node_ID + 1
      d_hre( i1_hre:i2_hre) => d_nih( i1_hre:i2_hre)
      d_hri( i1_hri:i2_hri) => d_nih( i1_hri:i2_hri)
      n_hre = i2_hre + 1 - i1_hre
      n_hri = i2_hri + 1 - i1_hri
    end if

#if (DO_ASSERTIONS)
    ! Safety: check that sending/receiving nodes agree on where the halos are

    ! Send halo ranges from left to right
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      range_send = [i1_hri, i2_hri]
      call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_right, &
        14, par%mpi_comm_node_primaries, req, ierr)
    end if
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_left, &
        14, par%mpi_comm_node_primaries, req, ierr)
      call MPI_WAIT( req, MPI_STATUS_IGNORE)
      call assert( (range_recv(1) == i1_hle .and. range_recv(2) == i2_hle), &
        'unmatched sending/receiving halo sizes')
    end if

    ! Send halo ranges from right to left
    if (par%node_ID > 0 .and. par%node_primary) then
      range_send = [i1_hli, i2_hli]
      call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_left, &
        15, par%mpi_comm_node_primaries, req, ierr)
    end if
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_right, &
        15, par%mpi_comm_node_primaries, req, ierr)
      call MPI_WAIT( req, MPI_STATUS_IGNORE)
      call assert( (range_recv(1) == i1_hre .and. range_recv(2) == i2_hre), &
        'unmatched sending/receiving halo sizes')
    end if

#endif

    ! Send halos from left to right
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_ISEND( d_hri, n_hri, MPI_DOUBLE_PRECISION, node_ID_right, &
        13, par%mpi_comm_node_primaries, reqs(1), ierr)
    end if
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_IRECV( d_hle, n_hle, MPI_DOUBLE_PRECISION, node_ID_left, &
        13, par%mpi_comm_node_primaries, reqs(1), ierr)
    end if

    ! Send halos from right to left
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_ISEND( d_hli, n_hli, MPI_DOUBLE_PRECISION, node_ID_left, &
        37, par%mpi_comm_node_primaries, reqs(2), ierr)
    end if
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_IRECV( d_hre, n_hre, MPI_DOUBLE_PRECISION, node_ID_right, &
        37, par%mpi_comm_node_primaries, reqs(2), ierr)
    end if

    if (par%node_primary) then
      call MPI_WAITALL( 2, reqs, MPI_STATUSES_IGNORE)
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp_1D

  subroutine exchange_halos_dp_2D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, n2)

    ! In/output variables:
    real(dp), dimension(i1_nih:i2_nih,n2), target, intent(inout) :: d_nih
    integer,                                       intent(in   ) :: i1_nih, i2_nih
    integer,                                       intent(in   ) :: i1_hle, i2_hle
    integer,                                       intent(in   ) :: i1_hli, i2_hli
    integer,                                       intent(in   ) :: i1_hre, i2_hre
    integer,                                       intent(in   ) :: i1_hri, i2_hri
    integer,                                       intent(in   ) :: n2

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'exchange_halos_dp_2D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_nih_1D

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, n2
      d_nih_1D( i1_nih:i2_nih) => d_nih( i1_nih:i2_nih,k)
      call exchange_halos_dp_1D( d_nih_1D, i1_nih, i2_nih, &
        i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)
      nullify( d_nih_1D)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp_2D

  subroutine exchange_halos_dp_3D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, n2, n3)

    ! In/output variables:
    real(dp), dimension(i1_nih:i2_nih,n2,n3), target, intent(inout) :: d_nih
    integer,                                          intent(in   ) :: i1_nih, i2_nih
    integer,                                          intent(in   ) :: i1_hle, i2_hle
    integer,                                          intent(in   ) :: i1_hli, i2_hli
    integer,                                          intent(in   ) :: i1_hre, i2_hre
    integer,                                          intent(in   ) :: i1_hri, i2_hri
    integer,                                          intent(in   ) :: n2, n3

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'exchange_halos_dp_3D'
    integer                         :: k,l
    real(dp), dimension(:), pointer :: d_nih_1D

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, n2
      do l = 1, n3
        d_nih_1D( i1_nih:i2_nih) => d_nih( i1_nih:i2_nih,k,l)
        call exchange_halos_dp_1D( d_nih_1D, i1_nih, i2_nih, &
          i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)
        nullify( d_nih_1D)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp_3D

  subroutine exchange_halos_complex_1D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)

    ! In/output variables:
    complex*16, dimension(i1_nih:i2_nih), target, intent(inout) :: d_nih
    integer,                                      intent(in   ) :: i1_nih, i2_nih
    integer,                                      intent(in   ) :: i1_hle, i2_hle
    integer,                                      intent(in   ) :: i1_hli, i2_hli
    integer,                                      intent(in   ) :: i1_hre, i2_hre
    integer,                                      intent(in   ) :: i1_hri, i2_hri

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'exchange_halos_complex_1D'
    integer                           :: node_ID_left, node_ID_right
    integer                           :: n_hle, n_hli, n_hre, n_hri
    complex*16, dimension(:), pointer :: d_hle, d_hli, d_hre, d_hri
    integer, dimension(2)             :: range_send, range_recv
    type(MPI_REQUEST)                 :: req
    type(MPI_REQUEST), dimension(2)   :: reqs
    integer                           :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Left halos
    if (par%node_ID == 0) then
      ! There is no node to the left
    else
      node_ID_left = par%node_ID - 1
      d_hle( i1_hle:i2_hle) => d_nih( i1_hle:i2_hle)
      d_hli( i1_hli:i2_hli) => d_nih( i1_hli:i2_hli)
      n_hle = i2_hle + 1 - i1_hle
      n_hli = i2_hli + 1 - i1_hli
    end if

    ! Right halos
    if (par%node_ID == par%n_nodes-1) then
      ! There is no node to the right
    else
      node_ID_right = par%node_ID + 1
      d_hre( i1_hre:i2_hre) => d_nih( i1_hre:i2_hre)
      d_hri( i1_hri:i2_hri) => d_nih( i1_hri:i2_hri)
      n_hre = i2_hre + 1 - i1_hre
      n_hri = i2_hri + 1 - i1_hri
    end if

#if (DO_ASSERTIONS)
    ! Safety: check that sending/receiving nodes agree on where the halos are

    ! Send halo ranges from left to right
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      range_send = [i1_hri, i2_hri]
      call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_right, &
        14, par%mpi_comm_node_primaries, req, ierr)
    end if
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_left, &
        14, par%mpi_comm_node_primaries, req, ierr)
      call MPI_WAIT( req, MPI_STATUS_IGNORE)
      call assert( (range_recv(1) == i1_hle .and. range_recv(2) == i2_hle), &
        'unmatched sending/receiving halo sizes')
    end if

    ! Send halo ranges from right to left
    if (par%node_ID > 0 .and. par%node_primary) then
      range_send = [i1_hli, i2_hli]
      call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_left, &
        15, par%mpi_comm_node_primaries, req, ierr)
    end if
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_right, &
        15, par%mpi_comm_node_primaries, req, ierr)
      call MPI_WAIT( req, MPI_STATUS_IGNORE)
      call assert( (range_recv(1) == i1_hre .and. range_recv(2) == i2_hre), &
        'unmatched sending/receiving halo sizes')
    end if

#endif

    ! Send halos from left to right
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_ISEND( d_hri, n_hri, MPI_DOUBLE_COMPLEX, node_ID_right, &
        13, par%mpi_comm_node_primaries, reqs(1), ierr)
    end if
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_IRECV( d_hle, n_hle, MPI_DOUBLE_COMPLEX, node_ID_left, &
        13, par%mpi_comm_node_primaries, reqs(1), ierr)
    end if

    ! Send halos from right to left
    if (par%node_ID > 0 .and. par%node_primary) then
      call MPI_ISEND( d_hli, n_hli, MPI_DOUBLE_COMPLEX, node_ID_left, &
        37, par%mpi_comm_node_primaries, reqs(2), ierr)
    end if
    if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
      call MPI_IRECV( d_hre, n_hre, MPI_DOUBLE_COMPLEX, node_ID_right, &
        37, par%mpi_comm_node_primaries, reqs(2), ierr)
    end if

    if (par%node_primary) then
      call MPI_WAITALL( 2, reqs, MPI_STATUSES_IGNORE)
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_complex_1D

  subroutine exchange_halos_complex_2D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, n2)

    ! In/output variables:
    complex*16, dimension(i1_nih:i2_nih,n2), target, intent(inout) :: d_nih
    integer,                                         intent(in   ) :: i1_nih, i2_nih
    integer,                                         intent(in   ) :: i1_hle, i2_hle
    integer,                                         intent(in   ) :: i1_hli, i2_hli
    integer,                                         intent(in   ) :: i1_hre, i2_hre
    integer,                                         intent(in   ) :: i1_hri, i2_hri
    integer,                                         intent(in   ) :: n2

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'exchange_halos_complex_2D'
    integer                         :: k
    complex*16, dimension(:), pointer :: d_nih_1D

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, n2
      d_nih_1D( i1_nih:i2_nih) => d_nih( i1_nih:i2_nih,k)
      call exchange_halos_complex_1D( d_nih_1D, i1_nih, i2_nih, &
        i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)
      nullify( d_nih_1D)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_complex_2D

  subroutine exchange_halos_complex_3D( d_nih, i1_nih, i2_nih, &
    i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri, n2, n3)

    ! In/output variables:
    complex*16, dimension(i1_nih:i2_nih,n2,n3), target, intent(inout) :: d_nih
    integer,                                            intent(in   ) :: i1_nih, i2_nih
    integer,                                            intent(in   ) :: i1_hle, i2_hle
    integer,                                            intent(in   ) :: i1_hli, i2_hli
    integer,                                            intent(in   ) :: i1_hre, i2_hre
    integer,                                            intent(in   ) :: i1_hri, i2_hri
    integer,                                            intent(in   ) :: n2, n3

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'exchange_halos_complex_3D'
    integer                           :: k,l
    complex*16, dimension(:), pointer :: d_nih_1D

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, n2
      do l = 1, n3
        d_nih_1D( i1_nih:i2_nih) => d_nih( i1_nih:i2_nih,k,l)
        call exchange_halos_complex_1D( d_nih_1D, i1_nih, i2_nih, &
          i1_hle, i2_hle, i1_hli, i2_hli, i1_hre, i2_hre, i1_hri, i2_hri)
        nullify( d_nih_1D)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_complex_3D

end module halo_exchange_mod
