module distribute_dist_shared_from_primary_mod

  use assertions_basic
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

  subroutine distribute_dist_shared_from_primary_logical_1D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, d_tot)

    ! In/output variables:
    logical, dimension(i1_nih:i2_nih), target, intent(  out) :: d_nih
    integer,                                   intent(in   ) :: i1_node, i2_node
    integer,                                   intent(in   ) :: i1_nih, i2_nih
    integer,                                   intent(in   ) :: n_tot
    logical, dimension(1:n_tot), optional,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_logical_1D'
    logical, dimension(:), pointer    :: d_interior
    integer                           :: n_interior, ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to gather the interior of each node
    n_interior = i2_node + 1 - i1_node
    d_interior( i1_node:i2_node) => d_nih( i1_node:i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( n_interior, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= n_tot) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_LOGICAL, &
        d_interior, n_interior, MPI_LOGICAL, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( i1_nih:i1_node-1) = .false.
      d_nih( i2_node+1:i2_nih) = .false.

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_1D

  subroutine distribute_dist_shared_from_primary_logical_2D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, d_tot)

    ! In/output variables:
    logical, dimension(i1_nih:i2_nih,1:nz), target,     intent(  out) :: d_nih
    integer,                                            intent(in   ) :: i1_node, i2_node
    integer,                                            intent(in   ) :: i1_nih, i2_nih
    integer,                                            intent(in   ) :: n_tot, nz
    logical, dimension(1:n_tot,1:nz), optional, target, intent(in   ) :: d_tot

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
        call distribute_dist_shared_from_primary_logical_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_logical_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_2D

  subroutine distribute_dist_shared_from_primary_logical_3D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, nl, d_tot)

    ! In/output variables:
    logical, dimension(i1_nih:i2_nih,1:nz,1:nl), target,     intent(  out) :: d_nih
    integer,                                                 intent(in   ) :: i1_node, i2_node
    integer,                                                 intent(in   ) :: i1_nih, i2_nih
    integer,                                                 intent(in   ) :: n_tot, nz, nl
    logical, dimension(1:n_tot,1:nz,1:nl), optional, target, intent(in   ) :: d_tot

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
          call distribute_dist_shared_from_primary_logical_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_logical_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_logical_3D

  subroutine distribute_dist_shared_from_primary_int_1D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, d_tot)

    ! In/output variables:
    integer, dimension(i1_nih:i2_nih), target, intent(  out) :: d_nih
    integer,                                   intent(in   ) :: i1_node, i2_node
    integer,                                   intent(in   ) :: i1_nih, i2_nih
    integer,                                   intent(in   ) :: n_tot
    integer, dimension(1:n_tot), optional,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_int_1D'
    integer, dimension(:), pointer    :: d_interior
    integer                           :: n_interior, ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to gather the interior of each node
    n_interior = i2_node + 1 - i1_node
    d_interior( i1_node:i2_node) => d_nih( i1_node:i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( n_interior, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= n_tot) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_INTEGER, &
        d_interior, n_interior, MPI_INTEGER, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( i1_nih:i1_node-1) = 0
      d_nih( i2_node+1:i2_nih) = 0

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_1D

  subroutine distribute_dist_shared_from_primary_int_2D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, d_tot)

    ! In/output variables:
    integer, dimension(i1_nih:i2_nih,1:nz), target,     intent(  out) :: d_nih
    integer,                                            intent(in   ) :: i1_node, i2_node
    integer,                                            intent(in   ) :: i1_nih, i2_nih
    integer,                                            intent(in   ) :: n_tot, nz
    integer, dimension(1:n_tot,1:nz), optional, target, intent(in   ) :: d_tot

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
        call distribute_dist_shared_from_primary_int_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_int_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_2D

  subroutine distribute_dist_shared_from_primary_int_3D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, nl, d_tot)

    ! In/output variables:
    integer, dimension(i1_nih:i2_nih,1:nz,1:nl), target,     intent(  out) :: d_nih
    integer,                                                 intent(in   ) :: i1_node, i2_node
    integer,                                                 intent(in   ) :: i1_nih, i2_nih
    integer,                                                 intent(in   ) :: n_tot, nz, nl
    integer, dimension(1:n_tot,1:nz,1:nl), optional, target, intent(in   ) :: d_tot

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
          call distribute_dist_shared_from_primary_int_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_int_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_int_3D

  subroutine distribute_dist_shared_from_primary_dp_1D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, d_tot)

    ! In/output variables:
    real(dp), dimension(i1_nih:i2_nih), target, intent(  out) :: d_nih
    integer,                                    intent(in   ) :: i1_node, i2_node
    integer,                                    intent(in   ) :: i1_nih, i2_nih
    integer,                                    intent(in   ) :: n_tot
    real(dp), dimension(1:n_tot), optional,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_dp_1D'
    real(dp), dimension(:), pointer   :: d_interior
    integer                           :: n_interior, ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to gather the interior of each node
    n_interior = i2_node + 1 - i1_node
    d_interior( i1_node:i2_node) => d_nih( i1_node:i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( n_interior, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= n_tot) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_PRECISION, &
        d_interior, n_interior, MPI_DOUBLE_PRECISION, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( i1_nih:i1_node-1) = 0._dp
      d_nih( i2_node+1:i2_nih) = 0._dp

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_1D

  subroutine distribute_dist_shared_from_primary_dp_2D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, d_tot)

    ! In/output variables:
    real(dp), dimension(i1_nih:i2_nih,1:nz), target,     intent(  out) :: d_nih
    integer,                                             intent(in   ) :: i1_node, i2_node
    integer,                                             intent(in   ) :: i1_nih, i2_nih
    integer,                                             intent(in   ) :: n_tot, nz
    real(dp), dimension(1:n_tot,1:nz), optional, target, intent(in   ) :: d_tot

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
        call distribute_dist_shared_from_primary_dp_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_dp_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_2D

  subroutine distribute_dist_shared_from_primary_dp_3D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, nl, d_tot)

    ! In/output variables:
    real(dp), dimension(i1_nih:i2_nih,1:nz,1:nl), target,     intent(  out) :: d_nih
    integer,                                                  intent(in   ) :: i1_node, i2_node
    integer,                                                  intent(in   ) :: i1_nih, i2_nih
    integer,                                                  intent(in   ) :: n_tot, nz, nl
    real(dp), dimension(1:n_tot,1:nz,1:nl), optional, target, intent(in   ) :: d_tot

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
          call distribute_dist_shared_from_primary_dp_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_dp_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_dp_3D

  subroutine distribute_dist_shared_from_primary_complex_1D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, d_tot)

    ! In/output variables:
    complex*16, dimension(i1_nih:i2_nih), target, intent(  out) :: d_nih
    integer,                                      intent(in   ) :: i1_node, i2_node
    integer,                                      intent(in   ) :: i1_nih, i2_nih
    integer,                                      intent(in   ) :: n_tot
    complex*16, dimension(1:n_tot), optional,     intent(in   ) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'distribute_dist_shared_from_primary_complex_1D'
    complex*16, dimension(:), pointer :: d_interior
    integer                           :: n_interior, ierr, i
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')
#endif

    ! We only need to gather the interior of each node
    n_interior = i2_node + 1 - i1_node
    d_interior( i1_node:i2_node) => d_nih( i1_node:i2_node)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_interior = d_tot
      call sync
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Determine ranges owned by each process
      call MPI_ALLGATHER( n_interior, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)

#if (DO_ASSERTIONS)
      if( sum( counts) /= n_tot) call crash('combined sizes of d_partial dont match size of d_tot')
#endif

      ! Calculate displacements for MPI_SCATTERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Scatter data from the primary
      call MPI_SCATTERV( d_tot, counts, displs, MPI_DOUBLE_COMPLEX, &
        d_interior, n_interior, MPI_DOUBLE_COMPLEX, 0, par%mpi_comm_node_primaries, ierr)

      ! Leave exterior halos empty
      d_nih( i1_nih:i1_node-1) = complex( 0._dp, 0._dp)
      d_nih( i2_node+1:i2_nih) = complex( 0._dp, 0._dp)

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_1D

  subroutine distribute_dist_shared_from_primary_complex_2D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, d_tot)

    ! In/output variables:
    complex*16, dimension(i1_nih:i2_nih,1:nz), target,     intent(  out) :: d_nih
    integer,                                               intent(in   ) :: i1_node, i2_node
    integer,                                               intent(in   ) :: i1_nih, i2_nih
    integer,                                               intent(in   ) :: n_tot, nz
    complex*16, dimension(1:n_tot,1:nz), optional, target, intent(in   ) :: d_tot

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
        call distribute_dist_shared_from_primary_complex_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call distribute_dist_shared_from_primary_complex_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
          n_tot)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_2D

  subroutine distribute_dist_shared_from_primary_complex_3D( d_nih, i1_node, i2_node, i1_nih, i2_nih, &
    n_tot, nz, nl, d_tot)

    ! In/output variables:
    complex*16, dimension(i1_nih:i2_nih,1:nz,1:nl), target,     intent(  out) :: d_nih
    integer,                                                    intent(in   ) :: i1_node, i2_node
    integer,                                                    intent(in   ) :: i1_nih, i2_nih
    integer,                                                    intent(in   ) :: n_tot, nz, nl
    complex*16, dimension(1:n_tot,1:nz,1:nl), optional, target, intent(in   ) :: d_tot

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
          call distribute_dist_shared_from_primary_complex_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call distribute_dist_shared_from_primary_complex_1D( d_nih_1D, i1_node, i2_node, i1_nih, i2_nih, &
            n_tot)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine distribute_dist_shared_from_primary_complex_3D

end module distribute_dist_shared_from_primary_mod
