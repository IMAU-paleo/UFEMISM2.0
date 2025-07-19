module distribute_dist_shared_from_primary_mod

  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mpi_f08, only: MPI_INTEGER, MPI_LOGICAL, MPI_ALLGATHER, MPI_SCATTERV, MPI_SEND, MPI_RECV, &
    MPI_ANY_TAG, MPI_STATUS, MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX
  use parallel_array_info_type, only: type_par_arr_info

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

  subroutine distribute_dist_shared_from_primary_logical_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    logical, dimension(pai%i1_nih:pai%i2_nih), target, intent(  out) :: d_nih
    logical, dimension(1:pai%n), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_logical_1D'
    logical, dimension(:), pointer    :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to distribute the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
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

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_LOGICAL, &
        d_interior, pai%n_node, MPI_LOGICAL, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( pai%i1_nih:pai%i1_node-1) = .false.
      d_nih( pai%i2_node+1:pai%i2_nih) = .false.

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_1D

  subroutine distribute_dist_shared_from_primary_logical_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(  out) :: d_nih
    logical, dimension(1:pai%n,1:nz), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'distribute_dist_shared_from_primary_logical_2D'
    logical, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call distribute_dist_shared_from_primary_logical_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_logical_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_2D

  subroutine distribute_dist_shared_from_primary_logical_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                     intent(in   ) :: pai
    integer,                                                     intent(in   ) :: nz, nl
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(  out) :: d_nih
    logical, dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'distribute_dist_shared_from_primary_logical_3D'
    logical, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call distribute_dist_shared_from_primary_logical_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_logical_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_3D

  subroutine distribute_dist_shared_from_primary_int_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    integer, dimension(pai%i1_nih:pai%i2_nih), target, intent(  out) :: d_nih
    integer, dimension(1:pai%n), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_int_1D'
    integer, dimension(:), pointer    :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to distribute the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
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

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_integer, &
        d_interior, pai%n_node, MPI_integer, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( pai%i1_nih:pai%i1_node-1) = 0
      d_nih( pai%i2_node+1:pai%i2_nih) = 0

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_1D

  subroutine distribute_dist_shared_from_primary_int_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(  out) :: d_nih
    integer, dimension(1:pai%n,1:nz), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'distribute_dist_shared_from_primary_int_2D'
    integer, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call distribute_dist_shared_from_primary_int_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_int_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_2D

  subroutine distribute_dist_shared_from_primary_int_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                     intent(in   ) :: pai
    integer,                                                     intent(in   ) :: nz, nl
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(  out) :: d_nih
    integer, dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'distribute_dist_shared_from_primary_int_3D'
    integer, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call distribute_dist_shared_from_primary_int_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_int_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_3D

  subroutine distribute_dist_shared_from_primary_dp_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                            intent(in   ) :: pai
    real(dp), dimension(pai%i1_nih:pai%i2_nih), target, intent(  out) :: d_nih
    real(dp), dimension(1:pai%n), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_dp_1D'
    real(dp), dimension(:), pointer   :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to distribute the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
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

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_PRECISION, &
        d_interior, pai%n_node, MPI_DOUBLE_PRECISION, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( pai%i1_nih:pai%i1_node-1) = 0._dp
      d_nih( pai%i2_node+1:pai%i2_nih) = 0._dp

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_1D

  subroutine distribute_dist_shared_from_primary_dp_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                 intent(in   ) :: pai
    integer,                                                 intent(in   ) :: nz
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(  out) :: d_nih
    real(dp), dimension(1:pai%n,1:nz), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'distribute_dist_shared_from_primary_dp_2D'
    real(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                         :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call distribute_dist_shared_from_primary_dp_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_dp_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_2D

  subroutine distribute_dist_shared_from_primary_dp_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                      intent(in   ) :: pai
    integer,                                                      intent(in   ) :: nz, nl
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(  out) :: d_nih
    real(dp), dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'distribute_dist_shared_from_primary_dp_3D'
    real(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                         :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call distribute_dist_shared_from_primary_dp_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_dp_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_3D

  subroutine distribute_dist_shared_from_primary_complex_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                              intent(in   ) :: pai
    complex*16, dimension(pai%i1_nih:pai%i2_nih), target, intent(  out) :: d_nih
    complex*16, dimension(1:pai%n), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_complex_1D'
    complex*16, dimension(:), pointer :: d_interior
    integer                           :: ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to distribute the interior of each node
    d_interior( pai%i1_node:pai%i2_node) => d_nih( pai%i1_node:pai%i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
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

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_COMPLEX, &
        d_interior, pai%n_node, MPI_DOUBLE_COMPLEX, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( pai%i1_nih:pai%i1_node-1) = 0._dp
      d_nih( pai%i2_node+1:pai%i2_nih) = 0._dp

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_1D

  subroutine distribute_dist_shared_from_primary_complex_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                   intent(in   ) :: pai
    integer,                                                   intent(in   ) :: nz
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(  out) :: d_nih
    complex*16, dimension(1:pai%n,1:nz), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_complex_2D'
    complex*16, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                           :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call distribute_dist_shared_from_primary_complex_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_complex_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_2D

  subroutine distribute_dist_shared_from_primary_complex_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                        intent(in   ) :: pai
    integer,                                                        intent(in   ) :: nz, nl
    complex*16, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(  out) :: d_nih
    complex*16, dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_complex_3D'
    complex*16, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                           :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call distribute_dist_shared_from_primary_complex_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_complex_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_3D

end module distribute_dist_shared_from_primary_mod
