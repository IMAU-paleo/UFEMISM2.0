module gather_dist_shared_to_all_mod

  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mpi_f08, only: MPI_INTEGER, MPI_ALLGATHER, MPI_ALLGATHERV, &
    MPI_STATUS, MPI_ANY_TAG, MPI_SEND, MPI_RECV, MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_DOUBLE_COMPLEX
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  private

  public :: gather_dist_shared_to_all

  interface gather_dist_shared_to_all
    procedure :: gather_dist_shared_to_all_logical_1D
    procedure :: gather_dist_shared_to_all_logical_2D
    procedure :: gather_dist_shared_to_all_logical_3D
    procedure :: gather_dist_shared_to_all_int_1D
    procedure :: gather_dist_shared_to_all_int_2D
    procedure :: gather_dist_shared_to_all_int_3D
    procedure :: gather_dist_shared_to_all_dp_1D
    procedure :: gather_dist_shared_to_all_dp_2D
    procedure :: gather_dist_shared_to_all_dp_3D
    procedure :: gather_dist_shared_to_all_complex_1D
    procedure :: gather_dist_shared_to_all_complex_2D
    procedure :: gather_dist_shared_to_all_complex_3D
  end interface gather_dist_shared_to_all

contains

  subroutine gather_dist_shared_to_all_logical_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    logical, dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    logical, dimension(1:pai%n), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_all_logical_1D'
    logical, dimension(:), pointer    :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! We only need to gather the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_interior
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( pai%n_node, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= pai%n) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_ALLGATHERV( d_interior, pai%n_node, MPI_LOGICAL, &
        d_tot, counts, displs, MPI_LOGICAL, par%mpi_comm_node_primaries, ierr)

    end if
    call sync_node

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_logical_1D

  subroutine gather_dist_shared_to_all_logical_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    logical, dimension(1:pai%n,1:nz), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_dist_shared_to_all_logical_2D'
    logical, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      d_nih_1D => d_nih(:,k)
      d_tot_1D => d_tot(:,k)
      call gather_dist_shared_to_all_logical_1D( pai, d_nih_1D, d_tot_1D)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_logical_2D

  subroutine gather_dist_shared_to_all_logical_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                     intent(in   ) :: pai
    integer,                                                     intent(in   ) :: nz, nl
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    logical, dimension(1:pai%n,1:nz,1:nl), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_dist_shared_to_all_logical_3D'
    logical, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      do l = 1, nl
        d_nih_1D => d_nih(:,k,l)
        d_tot_1D => d_tot(:,k,l)
        call gather_dist_shared_to_all_logical_1D( pai, d_nih_1D, d_tot_1D)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_logical_3D

  subroutine gather_dist_shared_to_all_int_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    integer, dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    integer, dimension(1:pai%n), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_all_int_1D'
    integer, dimension(:), pointer    :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! We only need to gather the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_interior
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( pai%n_node, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= pai%n) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_ALLGATHERV( d_interior, pai%n_node, MPI_INTEGER, &
        d_tot, counts, displs, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)

    end if
    call sync_node

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_int_1D

  subroutine gather_dist_shared_to_all_int_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    integer, dimension(1:pai%n,1:nz), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_dist_shared_to_all_int_2D'
    integer, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      d_nih_1D => d_nih(:,k)
      d_tot_1D => d_tot(:,k)
      call gather_dist_shared_to_all_int_1D( pai, d_nih_1D, d_tot_1D)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_int_2D

  subroutine gather_dist_shared_to_all_int_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                     intent(in   ) :: pai
    integer,                                                     intent(in   ) :: nz, nl
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    integer, dimension(1:pai%n,1:nz,1:nl), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_dist_shared_to_all_int_3D'
    integer, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      do l = 1, nl
        d_nih_1D => d_nih(:,k,l)
        d_tot_1D => d_tot(:,k,l)
        call gather_dist_shared_to_all_int_1D( pai, d_nih_1D, d_tot_1D)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_int_3D

  subroutine gather_dist_shared_to_all_dp_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                            intent(in   ) :: pai
    real(dp), dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    real(dp), dimension(1:pai%n), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_all_dp_1D'
    real(dp), dimension(:), pointer   :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! We only need to gather the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_interior
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( pai%n_node, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= pai%n) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_ALLGATHERV( d_interior, pai%n_node, MPI_DOUBLE_PRECISION, &
        d_tot, counts, displs, MPI_DOUBLE_PRECISION, par%mpi_comm_node_primaries, ierr)

    end if
    call sync_node

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_dp_1D

  subroutine gather_dist_shared_to_all_dp_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                 intent(in   ) :: pai
    integer,                                                 intent(in   ) :: nz
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    real(dp), dimension(1:pai%n,1:nz), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'gather_dist_shared_to_all_dp_2D'
    real(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                         :: k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      d_nih_1D => d_nih(:,k)
      d_tot_1D => d_tot(:,k)
      call gather_dist_shared_to_all_dp_1D( pai, d_nih_1D, d_tot_1D)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_dp_2D

  subroutine gather_dist_shared_to_all_dp_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                      intent(in   ) :: pai
    integer,                                                      intent(in   ) :: nz, nl
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    real(dp), dimension(1:pai%n,1:nz,1:nl), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'gather_dist_shared_to_all_dp_3D'
    real(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                         :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      do l = 1, nl
        d_nih_1D => d_nih(:,k,l)
        d_tot_1D => d_tot(:,k,l)
        call gather_dist_shared_to_all_dp_1D( pai, d_nih_1D, d_tot_1D)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_dp_3D

  subroutine gather_dist_shared_to_all_complex_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                              intent(in   ) :: pai
    complex*16, dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    complex*16, dimension(1:pai%n), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_all_complex_1D'
    complex*16, dimension(:), pointer :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! We only need to gather the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_interior
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( pai%n_node, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= pai%n) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_ALLGATHERV( d_interior, pai%n_node, MPI_DOUBLE_COMPLEX, &
        d_tot, counts, displs, MPI_DOUBLE_COMPLEX, par%mpi_comm_node_primaries, ierr)

    end if
    call sync_node

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_complex_1D

  subroutine gather_dist_shared_to_all_complex_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                   intent(in   ) :: pai
    integer,                                                   intent(in   ) :: nz
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    complex*16, dimension(1:pai%n,1:nz), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_all_complex_2D'
    complex*16, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                           :: k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      d_nih_1D => d_nih(:,k)
      d_tot_1D => d_tot(:,k)
      call gather_dist_shared_to_all_complex_1D( pai, d_nih_1D, d_tot_1D)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_complex_2D

  subroutine gather_dist_shared_to_all_complex_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                        intent(in   ) :: pai
    integer,                                                        intent(in   ) :: nz, nl
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    complex*16, dimension(1:pai%n,1:nz,1:nl), target,               intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_all_complex_3D'
    complex*16, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                           :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, nz
      do l = 1, nl
        d_nih_1D => d_nih(:,k,l)
        d_tot_1D => d_tot(:,k,l)
        call gather_dist_shared_to_all_complex_1D( pai, d_nih_1D, d_tot_1D)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_all_complex_3D

end module gather_dist_shared_to_all_mod
