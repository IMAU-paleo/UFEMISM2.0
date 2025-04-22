module halo_exchange_mod

  use assertions_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par, sync
  use mpi_f08, only: MPI_ISEND, MPI_IRECV, MPI_INTEGER, MPI_REQUEST, MPI_WAITALL, MPI_STATUSES_IGNORE, &
    MPI_WAIT, MPI_STATUS_IGNORE, MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_LOGICAL
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  private

  public :: basic_halo_exchange

  interface basic_halo_exchange
    procedure :: exchange_halos_logical_1D
    procedure :: exchange_halos_logical_2D
    procedure :: exchange_halos_logical_3D
    procedure :: exchange_halos_int_1D
    procedure :: exchange_halos_int_2D
    procedure :: exchange_halos_int_3D
    procedure :: exchange_halos_dp_1D
    procedure :: exchange_halos_dp_2D
    procedure :: exchange_halos_dp_3D
    procedure :: exchange_halos_complex_1D
    procedure :: exchange_halos_complex_2D
    procedure :: exchange_halos_complex_3D
  end interface basic_halo_exchange

contains

subroutine exchange_halos_logical_1D( pai, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                           intent(in   ) :: pai
  logical, dimension(pai%i1_nih:pai%i2_nih), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter  :: routine_name = 'exchange_halos_logical_1D'
  integer                         :: node_ID_left, node_ID_right
  logical, dimension(:), pointer  :: d_hle, d_hli, d_hre, d_hri
  integer, dimension(2)           :: range_send, range_recv
  type(MPI_REQUEST)               :: req
  type(MPI_REQUEST), dimension(2) :: reqs
  integer                         :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  ! Left halos
  if (par%node_ID == 0) then
    ! There is no node to the left
  else
    node_ID_left = par%node_ID - 1
    d_hle( pai%i1_hle:pai%i2_hle) => d_nih( pai%i1_hle:pai%i2_hle)
    d_hli( pai%i1_hli:pai%i2_hli) => d_nih( pai%i1_hli:pai%i2_hli)
  end if

  ! Right halos
  if (par%node_ID == par%n_nodes-1) then
    ! There is no node to the right
  else
    node_ID_right = par%node_ID + 1
    d_hre( pai%i1_hre:pai%i2_hre) => d_nih( pai%i1_hre:pai%i2_hre)
    d_hri( pai%i1_hri:pai%i2_hri) => d_nih( pai%i1_hri:pai%i2_hri)
  end if

#if (DO_ASSERTIONS)
  ! Safety: check that sending/receiving nodes agree on where the halos are

  ! Send halo ranges from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    range_send = [pai%i1_hri, pai%i2_hri]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_right, &
      14, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_left, &
      14, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hle .and. range_recv(2) == pai%i2_hle), &
      'unmatched sending/receiving halo sizes')
  end if

  ! Send halo ranges from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    range_send = [pai%i1_hli, pai%i2_hli]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_left, &
      15, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_right, &
      15, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hre .and. range_recv(2) == pai%i2_hre), &
      'unmatched sending/receiving halo sizes')
  end if

#endif

  ! Send halos from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_ISEND( d_hri, pai%n_hri, MPI_LOGICAL, node_ID_right, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( d_hle, pai%n_hle, MPI_LOGICAL, node_ID_left, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if

  ! Send halos from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_ISEND( d_hli, pai%n_hli, MPI_LOGICAL, node_ID_left, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( d_hre, pai%n_hre, MPI_LOGICAL, node_ID_right, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if

  if (par%node_primary) then
    call MPI_WAITALL( 2, reqs, MPI_STATUSES_IGNORE)
  end if
  call sync

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_logical_1D

subroutine exchange_halos_logical_2D( pai, nz, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                intent(in   ) :: pai
  integer,                                                intent(in   ) :: nz
  logical, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'exchange_halos_logical_2D'
  integer                        :: k
  logical, dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k)
    call exchange_halos_logical_1D( pai, d_nih_1D)
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_logical_2D

subroutine exchange_halos_logical_3D( pai, nz, nl, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                     intent(in   ) :: pai
  integer,                                                     intent(in   ) :: nz, nl
  logical, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'exchange_halos_logical_3D'
  integer                        :: k,l
  logical, dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    do l = 1, nl
      d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k,l)
      call exchange_halos_logical_1D( pai, d_nih_1D)
    end do
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_logical_3D

subroutine exchange_halos_int_1D( pai, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                           intent(in   ) :: pai
  integer, dimension(pai%i1_nih:pai%i2_nih), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter  :: routine_name = 'exchange_halos_int_1D'
  integer                         :: node_ID_left, node_ID_right
  integer, dimension(:), pointer  :: d_hle, d_hli, d_hre, d_hri
  integer, dimension(2)           :: range_send, range_recv
  type(MPI_REQUEST)               :: req
  type(MPI_REQUEST), dimension(2) :: reqs
  integer                         :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  ! Left halos
  if (par%node_ID == 0) then
    ! There is no node to the left
  else
    node_ID_left = par%node_ID - 1
    d_hle( pai%i1_hle:pai%i2_hle) => d_nih( pai%i1_hle:pai%i2_hle)
    d_hli( pai%i1_hli:pai%i2_hli) => d_nih( pai%i1_hli:pai%i2_hli)
  end if

  ! Right halos
  if (par%node_ID == par%n_nodes-1) then
    ! There is no node to the right
  else
    node_ID_right = par%node_ID + 1
    d_hre( pai%i1_hre:pai%i2_hre) => d_nih( pai%i1_hre:pai%i2_hre)
    d_hri( pai%i1_hri:pai%i2_hri) => d_nih( pai%i1_hri:pai%i2_hri)
  end if

#if (DO_ASSERTIONS)
  ! Safety: check that sending/receiving nodes agree on where the halos are

  ! Send halo ranges from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    range_send = [pai%i1_hri, pai%i2_hri]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_right, &
      14, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_left, &
      14, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hle .and. range_recv(2) == pai%i2_hle), &
      'unmatched sending/receiving halo sizes')
  end if

  ! Send halo ranges from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    range_send = [pai%i1_hli, pai%i2_hli]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_left, &
      15, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_right, &
      15, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hre .and. range_recv(2) == pai%i2_hre), &
      'unmatched sending/receiving halo sizes')
  end if

#endif

  ! Send halos from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_ISEND( d_hri, pai%n_hri, MPI_INTEGER, node_ID_right, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( d_hle, pai%n_hle, MPI_INTEGER, node_ID_left, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if

  ! Send halos from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_ISEND( d_hli, pai%n_hli, MPI_INTEGER, node_ID_left, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( d_hre, pai%n_hre, MPI_INTEGER, node_ID_right, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if

  if (par%node_primary) then
    call MPI_WAITALL( 2, reqs, MPI_STATUSES_IGNORE)
  end if
  call sync

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_int_1D

subroutine exchange_halos_int_2D( pai, nz, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                intent(in   ) :: pai
  integer,                                                intent(in   ) :: nz
  integer, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'exchange_halos_int_2D'
  integer                        :: k
  integer, dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k)
    call exchange_halos_int_1D( pai, d_nih_1D)
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_int_2D

subroutine exchange_halos_int_3D( pai, nz, nl, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                     intent(in   ) :: pai
  integer,                                                     intent(in   ) :: nz, nl
  integer, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'exchange_halos_int_3D'
  integer                        :: k,l
  integer, dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    do l = 1, nl
      d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k,l)
      call exchange_halos_int_1D( pai, d_nih_1D)
    end do
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_int_3D

subroutine exchange_halos_dp_1D( pai, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                            intent(in   ) :: pai
  real(dp), dimension(pai%i1_nih:pai%i2_nih), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter  :: routine_name = 'exchange_halos_dp_1D'
  integer                         :: node_ID_left, node_ID_right
  real(dp), dimension(:), pointer :: d_hle, d_hli, d_hre, d_hri
  integer, dimension(2)           :: range_send, range_recv
  type(MPI_REQUEST)               :: req
  type(MPI_REQUEST), dimension(2) :: reqs
  integer                         :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  ! Left halos
  if (par%node_ID == 0) then
    ! There is no node to the left
  else
    node_ID_left = par%node_ID - 1
    d_hle( pai%i1_hle:pai%i2_hle) => d_nih( pai%i1_hle:pai%i2_hle)
    d_hli( pai%i1_hli:pai%i2_hli) => d_nih( pai%i1_hli:pai%i2_hli)
  end if

  ! Right halos
  if (par%node_ID == par%n_nodes-1) then
    ! There is no node to the right
  else
    node_ID_right = par%node_ID + 1
    d_hre( pai%i1_hre:pai%i2_hre) => d_nih( pai%i1_hre:pai%i2_hre)
    d_hri( pai%i1_hri:pai%i2_hri) => d_nih( pai%i1_hri:pai%i2_hri)
  end if

#if (DO_ASSERTIONS)
  ! Safety: check that sending/receiving nodes agree on where the halos are

  ! Send halo ranges from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    range_send = [pai%i1_hri, pai%i2_hri]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_right, &
      14, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_left, &
      14, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hle .and. range_recv(2) == pai%i2_hle), &
      'unmatched sending/receiving halo sizes')
  end if

  ! Send halo ranges from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    range_send = [pai%i1_hli, pai%i2_hli]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_left, &
      15, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_right, &
      15, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hre .and. range_recv(2) == pai%i2_hre), &
      'unmatched sending/receiving halo sizes')
  end if

#endif

  ! Send halos from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_ISEND( d_hri, pai%n_hri, MPI_DOUBLE_PRECISION, node_ID_right, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( d_hle, pai%n_hle, MPI_DOUBLE_PRECISION, node_ID_left, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if

  ! Send halos from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_ISEND( d_hli, pai%n_hli, MPI_DOUBLE_PRECISION, node_ID_left, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( d_hre, pai%n_hre, MPI_DOUBLE_PRECISION, node_ID_right, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if

  if (par%node_primary) then
    call MPI_WAITALL( 2, reqs, MPI_STATUSES_IGNORE)
  end if
  call sync

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_dp_1D

subroutine exchange_halos_dp_2D( pai, nz, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                 intent(in   ) :: pai
  integer,                                                 intent(in   ) :: nz
  real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter  :: routine_name = 'exchange_halos_dp_2D'
  integer                         :: k
  real(dp), dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k)
    call exchange_halos_dp_1D( pai, d_nih_1D)
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_dp_2D

subroutine exchange_halos_dp_3D( pai, nz, nl, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                      intent(in   ) :: pai
  integer,                                                      intent(in   ) :: nz, nl
  real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter  :: routine_name = 'exchange_halos_dp_3D'
  integer                         :: k,l
  real(dp), dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    do l = 1, nl
      d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k,l)
      call exchange_halos_dp_1D( pai, d_nih_1D)
    end do
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_dp_3D

subroutine exchange_halos_complex_1D( pai, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                              intent(in   ) :: pai
  complex*16, dimension(pai%i1_nih:pai%i2_nih), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter    :: routine_name = 'exchange_halos_complex_1D'
  integer                           :: node_ID_left, node_ID_right
  complex*16, dimension(:), pointer :: d_hle, d_hli, d_hre, d_hri
  integer, dimension(2)             :: range_send, range_recv
  type(MPI_REQUEST)                 :: req
  type(MPI_REQUEST), dimension(2)   :: reqs
  integer                           :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  ! Left halos
  if (par%node_ID == 0) then
    ! There is no node to the left
  else
    node_ID_left = par%node_ID - 1
    d_hle( pai%i1_hle:pai%i2_hle) => d_nih( pai%i1_hle:pai%i2_hle)
    d_hli( pai%i1_hli:pai%i2_hli) => d_nih( pai%i1_hli:pai%i2_hli)
  end if

  ! Right halos
  if (par%node_ID == par%n_nodes-1) then
    ! There is no node to the right
  else
    node_ID_right = par%node_ID + 1
    d_hre( pai%i1_hre:pai%i2_hre) => d_nih( pai%i1_hre:pai%i2_hre)
    d_hri( pai%i1_hri:pai%i2_hri) => d_nih( pai%i1_hri:pai%i2_hri)
  end if

#if (DO_ASSERTIONS)
  ! Safety: check that sending/receiving nodes agree on where the halos are

  ! Send halo ranges from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    range_send = [pai%i1_hri, pai%i2_hri]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_right, &
      14, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_left, &
      14, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hle .and. range_recv(2) == pai%i2_hle), &
      'unmatched sending/receiving halo sizes')
  end if

  ! Send halo ranges from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    range_send = [pai%i1_hli, pai%i2_hli]
    call MPI_ISEND( range_send, 2, MPI_INTEGER, node_ID_left, &
      15, par%mpi_comm_node_primaries, req, ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( range_recv, 2, MPI_INTEGER, node_ID_right, &
      15, par%mpi_comm_node_primaries, req, ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE)
    call assert( (range_recv(1) == pai%i1_hre .and. range_recv(2) == pai%i2_hre), &
      'unmatched sending/receiving halo sizes')
  end if

#endif

  ! Send halos from left to right
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_ISEND( d_hri, pai%n_hri, MPI_DOUBLE_COMPLEX, node_ID_right, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_IRECV( d_hle, pai%n_hle, MPI_DOUBLE_COMPLEX, node_ID_left, &
      13, par%mpi_comm_node_primaries, reqs(1), ierr)
  end if

  ! Send halos from right to left
  if (par%node_ID > 0 .and. par%node_primary) then
    call MPI_ISEND( d_hli, pai%n_hli, MPI_DOUBLE_COMPLEX, node_ID_left, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if
  if (par%node_ID < par%n_nodes-1 .and. par%node_primary) then
    call MPI_IRECV( d_hre, pai%n_hre, MPI_DOUBLE_COMPLEX, node_ID_right, &
      37, par%mpi_comm_node_primaries, reqs(2), ierr)
  end if

  if (par%node_primary) then
    call MPI_WAITALL( 2, reqs, MPI_STATUSES_IGNORE)
  end if
  call sync

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_complex_1D

subroutine exchange_halos_complex_2D( pai, nz, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                   intent(in   ) :: pai
  integer,                                                   intent(in   ) :: nz
  complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter    :: routine_name = 'exchange_halos_complex_2D'
  integer                           :: k
  complex*16, dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k)
    call exchange_halos_complex_1D( pai, d_nih_1D)
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_complex_2D

subroutine exchange_halos_complex_3D( pai, nz, nl, d_nih)

  ! In/output variables:
  type(type_par_arr_info),                                        intent(in   ) :: pai
  integer,                                                        intent(in   ) :: nz, nl
  complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(inout) :: d_nih

  ! Local variables:
  character(len=1024), parameter    :: routine_name = 'exchange_halos_complex_3D'
  integer                           :: k,l
  complex*16, dimension(:), pointer :: d_nih_1D

  ! Add routine to path
  call init_routine( routine_name)

  ! If running on one node, do nothing
  if (par%n_nodes == 1) then
    call finalise_routine( routine_name)
    return
  end if

  do k = 1, nz
    do l = 1, nl
      d_nih_1D( pai%i1_nih:pai%i2_nih) => d_nih( pai%i1_nih:pai%i2_nih,k,l)
      call exchange_halos_complex_1D( pai, d_nih_1D)
    end do
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine exchange_halos_complex_3D

end module halo_exchange_mod
