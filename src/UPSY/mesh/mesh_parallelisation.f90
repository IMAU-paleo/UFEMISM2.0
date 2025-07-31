module mesh_parallelisation

  ! Imagine a mesh of 100 vertices, with the model running on 2 shared-memory nodes of 2 processes
  ! each, so 4 processes in total. Each process will own 25 vertices, each node owns 50 vertices. So:
  !
  !   Node         0    0    1    1
  !   Process      0    1    2    3
  !
  !   vi1          1   26   51   76        => Each process "owns" vertices vi1-vi2
  !   vi2         25   50   75  100
  !   vi1_node     1    1   51   51        => Each node "owns" vertices vi1_node-vi2_node
  !   vi2_node    50   50  100  100
  !
  ! Because the vertices are ordered in the x-direction, the first 50 vertices, owned by
  ! node 0, cover pretty much the left half of the domain, while the second 50 vertices, owned
  ! by node 1, cover the right half. Let's say vertices 41-50 lie along the "border", i.e. they
  ! are owned by node 0, but all have at least one neighbour that is owned by node 1. On the other
  ! side of the border, vertices 51-60 are owned by node 1, but all have at least one neighbour
  ! that is owned by node 0. These groups of vertices constitute the halos. From the point of
  ! view of node 0, vertices 41-50 form the right interior halo, while vertices 51-60 form
  ! the right exterior halo. Conversely, from the point of view of node 1, vertices 41-50 form
  ! the left exterior halo, while vertices 51-60 form the left interior halo. So:
  !
  !   Node         0    0    1    1
  !   Process      0    1    2    3
  !
  !   vi1_nih      1    1   41   41        => nih = "node including halo", so the 50 vertices owned
  !   vi2_nih     60   60  100  100           by each node, plus the exterior halos
  !   vi1_hle      0    0   41   41        => hle = "halo left exterior"
  !   vi2_hle     -1   -1   50   50
  !   vi1_hli      0    0   51   51        => hli = "halo left interior"
  !   vi2_hli     -1   -1   60   60
  !   vi1_hre     51   51    0    0        => hre = "halo right exterior"
  !   vi2_hre     60   60   -1   -1
  !   vi1_hri     41   41    0    0        => hri = "halo right interior"
  !   vi2_hri     50   50   -1   -1
  !
  ! Each node allocates [vi1_nih:vi2_nih] of node-shared memory to store a data field,
  ! but "operates" only on vi1_node-vi2_node. Because of the left-right ordering of vertices,
  ! the data looks like:
  !
  !  [  vi1_nih                                :                                 vi2_nih  ]
  !  [ [vi1_hle:vi2_hle] [  vi1_node           :            vi2_node  ] [vi1_hre:vi2_hre] ]
  !  [ [vi1_hle:vi2_hle] [ [vi1_hli:vi2_hri] [...] [vi1_hri:vi2_hri ] ] [vi1_hre:vi2_hre] ]
  !
  ! When a gradient of the field needs to be
  ! calculated, the nodes perform a "halo exchange", in order to obtain the values of the
  ! field in their own exterior halos. In this example, node 0 will send the values on
  ! vertices 41-50 to process 1, while process 1 will send the values on vertices 51-60
  ! to node 0.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mpi_basic, only: par, sync
  use mpi_distributed_memory, only: partition_list
  use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_MIN, MPI_MAX, MPI_IN_PLACE, MPI_COMM_WORLD, &
    MPI_ALLGATHER, MPI_ALLGATHERV, MPI_WIN
  use mpi_distributed_shared_memory, only: allocate_dist_shared, gather_dist_shared_to_all, deallocate_dist_shared

  implicit none

  private

  public :: setup_mesh_parallelisation

contains

  subroutine setup_mesh_parallelisation( mesh, mask_active_a_tot, mask_active_b_tot)
    !< Setup the parallelisation of memory and halos on them mesh

    ! In/output variables:
    type(type_mesh),                         intent(inout) :: mesh
    logical, dimension(mesh%nV),   optional, intent(in   ) :: mask_active_a_tot
    logical, dimension(mesh%nTri), optional, intent(in   ) :: mask_active_b_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_mesh_parallelisation'

    ! Add routine to path
    call init_routine( routine_name)

    allocate( mesh%V_owning_process  ( mesh%nV))
    allocate( mesh%V_owning_node     ( mesh%nV))
    allocate( mesh%Tri_owning_process( mesh%nTri))
    allocate( mesh%Tri_owning_node   ( mesh%nTri))
    allocate( mesh%E_owning_process  ( mesh%nE))
    allocate( mesh%E_owning_node     ( mesh%nE))

    if (.not. present( mask_active_a_tot)) then
      ! Divide all vertices equally over the processes
      call determine_ownership_ranges_equal( mesh%nV, mesh%vi1, mesh%vi2, mesh%nV_loc, &
        mesh%pai_V%i1_node, mesh%pai_V%i2_node, mesh%pai_V%n_node, &
        mesh%V_owning_process, mesh%V_owning_node)
    else
      ! Divide vertices over the processes so each process has the same number of active vertices
      call determine_ownership_ranges_balanced( mesh%nV, mesh%vi1, mesh%vi2, mesh%nV_loc, &
        mesh%pai_V%i1_node, mesh%pai_V%i2_node, mesh%pai_V%n_node, &
        mesh%V_owning_process, mesh%V_owning_node, mask_active_a_tot)
    end if

    if (.not. present( mask_active_b_tot)) then
      ! Divide all triangles equally over the processes
      call determine_ownership_ranges_equal( mesh%nTri, mesh%ti1, mesh%ti2, mesh%nTri_loc, &
        mesh%pai_Tri%i1_node, mesh%pai_Tri%i2_node, mesh%pai_Tri%n_node, &
        mesh%Tri_owning_process, mesh%Tri_owning_node)
    else
      call determine_ownership_ranges_balanced( mesh%nTri, mesh%ti1, mesh%ti2, mesh%nTri_loc, &
        mesh%pai_Tri%i1_node, mesh%pai_Tri%i2_node, mesh%pai_Tri%n_node, &
        mesh%Tri_owning_process, mesh%Tri_owning_node, mask_active_b_tot)
    end if

    call determine_ownership_ranges_equal( mesh%nE, mesh%ei1, mesh%ei2, mesh%nE_loc, &
      mesh%pai_E%i1_node, mesh%pai_E%i2_node, mesh%pai_E%n_node, &
      mesh%E_owning_process, mesh%E_owning_node)

    mesh%pai_V%i1 = mesh%vi1
    mesh%pai_V%i2 = mesh%vi2
    mesh%pai_V%n_loc = mesh%nV_loc

    mesh%pai_Tri%i1 = mesh%ti1
    mesh%pai_Tri%i2 = mesh%ti2
    mesh%pai_Tri%n_loc = mesh%nTri_loc

    mesh%pai_E%i1 = mesh%ei1
    mesh%pai_E%i2 = mesh%ei2
    mesh%pai_E%n_loc = mesh%nE_loc

    ! Determine all halos
    call determine_halos( mesh)

    ! Allocate buffer shared memory for e.g. matrix multiplications
    call allocate_dist_shared( mesh%buffer1_d_a_nih , mesh%wbuffer1_d_a_nih , mesh%pai_V%n_nih)
    call allocate_dist_shared( mesh%buffer2_d_a_nih , mesh%wbuffer2_d_a_nih , mesh%pai_V%n_nih)
    call allocate_dist_shared( mesh%buffer1_d_ak_nih, mesh%wbuffer1_d_ak_nih, mesh%pai_V%n_nih,   mesh%nz)
    call allocate_dist_shared( mesh%buffer2_d_ak_nih, mesh%wbuffer2_d_ak_nih, mesh%pai_V%n_nih,   mesh%nz)
    mesh%buffer1_d_a_nih(  mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih             ) => mesh%buffer1_d_a_nih
    mesh%buffer2_d_a_nih(  mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih             ) => mesh%buffer2_d_a_nih
    mesh%buffer1_d_ak_nih( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  , 1:mesh%nz) => mesh%buffer1_d_ak_nih
    mesh%buffer2_d_ak_nih( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  , 1:mesh%nz) => mesh%buffer2_d_ak_nih

    call allocate_dist_shared( mesh%buffer1_d_b_nih , mesh%wbuffer1_d_b_nih , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( mesh%buffer2_d_b_nih , mesh%wbuffer2_d_b_nih , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( mesh%buffer1_d_bk_nih, mesh%wbuffer1_d_bk_nih, mesh%pai_Tri%n_nih, mesh%nz)
    call allocate_dist_shared( mesh%buffer2_d_bk_nih, mesh%wbuffer2_d_bk_nih, mesh%pai_Tri%n_nih, mesh%nz)
    mesh%buffer1_d_b_nih(  mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih           ) => mesh%buffer1_d_b_nih
    mesh%buffer2_d_b_nih(  mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih           ) => mesh%buffer2_d_b_nih
    mesh%buffer1_d_bk_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih, 1:mesh%nz) => mesh%buffer1_d_bk_nih
    mesh%buffer2_d_bk_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih, 1:mesh%nz) => mesh%buffer2_d_bk_nih

    call allocate_dist_shared( mesh%buffer1_d_c_nih , mesh%wbuffer1_d_c_nih , mesh%pai_E%n_nih)
    call allocate_dist_shared( mesh%buffer2_d_c_nih , mesh%wbuffer2_d_c_nih , mesh%pai_E%n_nih)
    call allocate_dist_shared( mesh%buffer1_d_ck_nih, mesh%wbuffer1_d_ck_nih, mesh%pai_E%n_nih,   mesh%nz)
    call allocate_dist_shared( mesh%buffer2_d_ck_nih, mesh%wbuffer2_d_ck_nih, mesh%pai_E%n_nih,   mesh%nz)
    mesh%buffer1_d_c_nih(  mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih             ) => mesh%buffer1_d_c_nih
    mesh%buffer2_d_c_nih(  mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih             ) => mesh%buffer2_d_c_nih
    mesh%buffer1_d_ck_nih( mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih  , 1:mesh%nz) => mesh%buffer1_d_ck_nih
    mesh%buffer2_d_ck_nih( mesh%pai_E%i1_nih  :mesh%pai_E%i2_nih  , 1:mesh%nz) => mesh%buffer2_d_ck_nih

    ! call print_parallelisation_info( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 12)

  end subroutine setup_mesh_parallelisation

  subroutine determine_ownership_ranges_equal( n, n1, n2, n_loc, n1_node, n2_node, n_node, &
    owning_process, owning_node)

    ! In/output variables:
    integer,               intent(in   ) :: n
    integer,               intent(  out) :: n1, n2, n_loc, n1_node, n2_node, n_node
    integer, dimension(n), intent(  out) :: owning_process, owning_node

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_ownership_ranges_equal'
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Range owned by this process
    call partition_list( n, par%i, par%n, n1, n2)
    n_loc = n2 + 1 - n1

    ! Range owned by this node
    call MPI_ALLREDUCE( n1, n1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( n2, n2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    n_node = n2_node + 1 - n1_node

    ! Find out which process/node owns each element
    owning_process = -1
    owning_process( n1:n2) = par%i
    call MPI_ALLREDUCE( MPI_IN_PLACE, owning_process, size( owning_process), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    owning_node = -1
    owning_node( n1_node:n2_node) = par%node_ID
    call MPI_ALLREDUCE( MPI_IN_PLACE, owning_node, size( owning_node), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_ownership_ranges_equal

  subroutine determine_ownership_ranges_balanced( n, n1, n2, n_loc, n1_node, n2_node, n_node, &
    owning_process, owning_node, mask_active_tot)

    ! In/output variables:
    integer,                         intent(in   ) :: n
    integer,                         intent(  out) :: n1, n2, n_loc, n1_node, n2_node, n_node
    integer, dimension(n),           intent(  out) :: owning_process, owning_node
    logical, dimension(n), optional, intent(in   ) :: mask_active_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_ownership_ranges_balanced'
    integer                        :: n_active_tot, i1, i2, n_active_proc
    integer, dimension(0:par%n-1)  :: n_active_proc_all
    integer                        :: ierr
    integer                        :: ip, n_active_proc_

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine range of vertices/triangles owned by each process, so that each of them
    ! owns the same number of masked ones

    n_active_tot = count( mask_active_tot)
    call partition_list( n_active_tot, par%i, par%n, i1, i2)
    n_active_proc = i2 + 1 - i1
    call MPI_ALLGATHER( n_active_proc, 1, MPI_INTEGER, n_active_proc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    i1 = -1
    i2 = 0

    do ip = 0, par%n-1

      i1 = i2+1
      i2 = i1-1
      n_active_proc_ = 0

      do while (n_active_proc_ < n_active_proc_all( ip))
        i2 = i2 + 1
        if (mask_active_tot( i2)) n_active_proc_ = n_active_proc_ + 1
      end do

      if (ip == par%n-1) i2 = n

      if (ip == par%i) then
        n1 = i1
        n2 = i2
        n_loc = n2 + 1 - n1
      end if

    end do

    ! Range owned by this node
    call MPI_ALLREDUCE( n1, n1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( n2, n2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    n_node = n2_node + 1 - n1_node

    ! Find out which process/node owns each element
    owning_process = -1
    owning_process( n1:n2) = par%i
    call MPI_ALLREDUCE( MPI_IN_PLACE, owning_process, size( owning_process), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    owning_node = -1
    owning_node( n1_node:n2_node) = par%node_ID
    call MPI_ALLREDUCE( MPI_IN_PLACE, owning_node, size( owning_node), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_ownership_ranges_balanced

  subroutine determine_halos( mesh)
    !< Determine all the halos

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_halos'
    integer                        :: node_ID_left, node_ID_right

    ! Add routine to path
    call init_routine( routine_name)

    ! Left ("west")
    if (par%node_ID == 0) then
      ! There is no node to the left

      mesh%pai_V%i1_nih  = mesh%pai_V%i1_node
      mesh%pai_V%n_hle   = 0
      mesh%pai_V%i1_hle  =  0
      mesh%pai_V%i2_hle  = -1
      mesh%pai_V%n_hli   = 0
      mesh%pai_V%i1_hli  =  0
      mesh%pai_V%i2_hli  = -1

      mesh%pai_Tri%i1_nih  = mesh%pai_Tri%i1_node
      mesh%pai_Tri%n_hle = 0
      mesh%pai_Tri%i1_hle  =  0
      mesh%pai_Tri%i2_hle  = -1
      mesh%pai_Tri%n_hli = 0
      mesh%pai_Tri%i1_hli  =  0
      mesh%pai_Tri%i2_hli  = -1

      mesh%pai_E%i1_nih  = mesh%pai_E%i1_node
      mesh%pai_E%n_hle   = 0
      mesh%pai_E%i1_hle  =  0
      mesh%pai_E%i2_hle  = -1
      mesh%pai_E%n_hli   = 0
      mesh%pai_E%i1_hli  =  0
      mesh%pai_E%i2_hli  = -1

    else
      node_ID_left = par%node_ID - 1

      call determine_halo_range_a( mesh, node_ID_left, par%node_ID , mesh%pai_V%i1_hle, mesh%pai_V%i2_hle)
      call determine_halo_range_a( mesh, par%node_ID , node_ID_left, mesh%pai_V%i1_hli, mesh%pai_V%i2_hli)
      mesh%pai_V%n_hle = mesh%pai_V%i2_hle + 1 - mesh%pai_V%i1_hle
      mesh%pai_V%n_hli = mesh%pai_V%i2_hli + 1 - mesh%pai_V%i1_hli
      mesh%pai_V%i1_nih = mesh%pai_V%i1_hle

      call determine_halo_range_b( mesh, node_ID_left, par%node_ID , mesh%pai_Tri%i1_hle, mesh%pai_Tri%i2_hle)
      call determine_halo_range_b( mesh, par%node_ID , node_ID_left, mesh%pai_Tri%i1_hli, mesh%pai_Tri%i2_hli)
      mesh%pai_Tri%n_hle = mesh%pai_Tri%i2_hle + 1 - mesh%pai_Tri%i1_hle
      mesh%pai_Tri%n_hli = mesh%pai_Tri%i2_hli + 1 - mesh%pai_Tri%i1_hli
      mesh%pai_Tri%i1_nih = mesh%pai_Tri%i1_hle

      call determine_halo_range_c( mesh, node_ID_left, par%node_ID , mesh%pai_E%i1_hle, mesh%pai_E%i2_hle)
      call determine_halo_range_c( mesh, par%node_ID , node_ID_left, mesh%pai_E%i1_hli, mesh%pai_E%i2_hli)
      mesh%pai_E%n_hle = mesh%pai_E%i2_hle + 1 - mesh%pai_E%i1_hle
      mesh%pai_E%n_hli = mesh%pai_E%i2_hli + 1 - mesh%pai_E%i1_hli
      mesh%pai_E%i1_nih = mesh%pai_E%i1_hle

    end if

    ! Right ("east")
    if (par%node_ID == par%n_nodes-1) then
      ! There is no node to the right

      mesh%pai_V%i2_nih  = mesh%pai_V%i2_node
      mesh%pai_V%n_hre   = 0
      mesh%pai_V%i1_hre  =  0
      mesh%pai_V%i2_hre  = -1
      mesh%pai_V%n_hri   = 0
      mesh%pai_V%i1_hri  =  0
      mesh%pai_V%i2_hri  = -1

      mesh%pai_Tri%i2_nih  = mesh%pai_Tri%i2_node
      mesh%pai_Tri%n_hre = 0
      mesh%pai_Tri%i1_hre  =  0
      mesh%pai_Tri%i2_hre  = -1
      mesh%pai_Tri%n_hri = 0
      mesh%pai_Tri%i1_hri  =  0
      mesh%pai_Tri%i2_hri  = -1

      mesh%pai_E%i2_nih  = mesh%pai_E%i2_node
      mesh%pai_E%n_hre   = 0
      mesh%pai_E%i1_hre  =  0
      mesh%pai_E%i2_hre  = -1
      mesh%pai_E%n_hri   = 0
      mesh%pai_E%i1_hri  =  0
      mesh%pai_E%i2_hri  = -1

    else
      node_ID_right = par%node_ID + 1

      call determine_halo_range_a( mesh, node_ID_right, par%node_ID  , mesh%pai_V%i1_hre, mesh%pai_V%i2_hre)
      call determine_halo_range_a( mesh, par%node_ID  , node_ID_right, mesh%pai_V%i1_hri, mesh%pai_V%i2_hri)
      mesh%pai_V%n_hre = mesh%pai_V%i2_hre + 1 - mesh%pai_V%i1_hre
      mesh%pai_V%n_hri = mesh%pai_V%i2_hri + 1 - mesh%pai_V%i1_hri
      mesh%pai_V%i2_nih = mesh%pai_V%i2_hre

      call determine_halo_range_b( mesh, node_ID_right, par%node_ID  , mesh%pai_Tri%i1_hre, mesh%pai_Tri%i2_hre)
      call determine_halo_range_b( mesh, par%node_ID  , node_ID_right, mesh%pai_Tri%i1_hri, mesh%pai_Tri%i2_hri)
      mesh%pai_Tri%n_hre = mesh%pai_Tri%i2_hre + 1 - mesh%pai_Tri%i1_hre
      mesh%pai_Tri%n_hri = mesh%pai_Tri%i2_hri + 1 - mesh%pai_Tri%i1_hri
      mesh%pai_Tri%i2_nih = mesh%pai_Tri%i2_hre

      call determine_halo_range_c( mesh, node_ID_right, par%node_ID  , mesh%pai_E%i1_hre, mesh%pai_E%i2_hre)
      call determine_halo_range_c( mesh, par%node_ID  , node_ID_right, mesh%pai_E%i1_hri, mesh%pai_E%i2_hri)
      mesh%pai_E%n_hre = mesh%pai_E%i2_hre + 1 - mesh%pai_E%i1_hre
      mesh%pai_E%n_hri = mesh%pai_E%i2_hri + 1 - mesh%pai_E%i1_hri
      mesh%pai_E%i2_nih = mesh%pai_E%i2_hre

    end if

    mesh%pai_V%n       = mesh%nV
    mesh%pai_Tri%n     = mesh%nTri
    mesh%pai_E%n       = mesh%nE

    mesh%pai_V%n_nih   = mesh%pai_V%i2_nih   + 1 - mesh%pai_V%i1_nih
    mesh%pai_Tri%n_nih = mesh%pai_Tri%i2_nih + 1 - mesh%pai_Tri%i1_nih
    mesh%pai_E%n_nih   = mesh%pai_E%i2_nih   + 1 - mesh%pai_E%i1_nih

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_halos

  subroutine determine_halo_range_a( mesh, node_in, node_next_to, vi1, vi2)
    !< Find all vertices that are owned by node_in but adjacent to node_next_to,
    !< and return their range vi1-vi2

    ! In/output variables:
    type(type_mesh), intent(in   ) :: mesh
    integer,         intent(in   ) :: node_in, node_next_to
    integer,         intent(  out) :: vi1, vi2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_halo_range_a'
    integer                        :: vi, ci, vj, iti, ti, ei
    logical                        :: has_vertices_in_node_next_to
    logical                        :: has_triangles_in_node_next_to
    logical                        :: has_edges_in_node_next_to

    ! Add routine to path
    call init_routine( routine_name)

    vi1 = mesh%nV
    vi2 = 1
    do vi = 1, mesh%nV
      if (mesh%V_owning_node( vi) == node_in) then
        ! This vertex is owned by node_in

        has_vertices_in_node_next_to  = .false.
        has_triangles_in_node_next_to = .false.
        has_edges_in_node_next_to     = .false.

        ! Check if this vertex has any neighbouring vertices owned by node_next_to
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mesh%V_owning_node( vj) == node_next_to) then
            has_vertices_in_node_next_to = .true.
            exit
          end if
        end do

        ! Check if this vertex has any neighbouring triangles owned by node_next_to
        do iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          if (mesh%Tri_owning_node( ti) == node_next_to) then
            has_triangles_in_node_next_to = .true.
            exit
          end if
        end do

        ! Check if this vertex has any neighbouring edges owned by node_next_to
        do ci = 1, mesh%nC( vi)
          ei = mesh%VE( vi,ci)
          if (mesh%E_owning_node( ei) == node_next_to) then
            has_edges_in_node_next_to = .true.
            exit
          end if
        end do

        if (has_vertices_in_node_next_to .or. &
            has_triangles_in_node_next_to .or. &
            has_edges_in_node_next_to) then
          vi1 = min( vi1, vi)
          vi2 = max( vi2, vi)
        end if

      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_halo_range_a

  subroutine determine_halo_range_b( mesh, node_in, node_next_to, ti1, ti2)
    !< Find all triangles that are owned by node_in but adjacent to node_next_to,
    !< and return their range ti1-ti2

    ! In/output variables:
    type(type_mesh), intent(in   ) :: mesh
    integer,         intent(in   ) :: node_in, node_next_to
    integer,         intent(  out) :: ti1, ti2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_halo_range_b'
    integer                        :: ti, n, vi, tj, ei
    logical                        :: has_vertices_in_node_next_to
    logical                        :: has_triangles_in_node_next_to
    logical                        :: has_edges_in_node_next_to

    ! Add routine to path
    call init_routine( routine_name)

    ti1 = mesh%nTri
    ti2 = 1
    do ti = 1, mesh%nTri
      if (mesh%Tri_owning_node( ti) == node_in) then
        ! This triangle is owned by node_in

        has_vertices_in_node_next_to  = .false.
        has_triangles_in_node_next_to = .false.
        has_edges_in_node_next_to     = .false.

        ! Check if this triangle has any neighbouring vertices owned by node_next_to
        do n = 1, 3
          vi = mesh%Tri( ti,n)
          if (mesh%V_owning_node( vi) == node_next_to) then
            has_vertices_in_node_next_to = .true.
            exit
          end if
        end do

        ! Check if this triangle has any neighbouring triangles owned by node_next_to
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj > 0) then
            if (mesh%Tri_owning_node( tj) == node_next_to) then
              has_triangles_in_node_next_to = .true.
              exit
            end if
          end if
        end do

        ! Check if this triangle has any neighbouring edges owned by node_next_to
        do n = 1, 3
          ei = mesh%TriE( ti,n)
          if (mesh%E_owning_node( ei) == node_next_to) then
            has_edges_in_node_next_to = .true.
            exit
          end if
        end do

        if (has_vertices_in_node_next_to .or. &
            has_triangles_in_node_next_to .or. &
            has_edges_in_node_next_to) then
          ti1 = min( ti1, ti)
          ti2 = max( ti2, ti)
        end if

      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_halo_range_b

  subroutine determine_halo_range_c( mesh, node_in, node_next_to, ei1, ei2)
    !< Find all edges that are owned by node_in but adjacent to node_next_to,
    !< and return their range ei1-ei2

    ! In/output variables:
    type(type_mesh), intent(in   ) :: mesh
    integer,         intent(in   ) :: node_in, node_next_to
    integer,         intent(  out) :: ei1, ei2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_halo_range_c'
    integer                        :: ei, vi, vj, vil, vir, til, tir
    logical                        :: has_vertices_in_node_next_to
    logical                        :: has_triangles_in_node_next_to
    logical                        :: has_edges_in_node_next_to

    ! Add routine to path
    call init_routine( routine_name)

    ei1 = mesh%nE
    ei2 = 1
    do ei = 1, mesh%nE
      if (mesh%E_owning_node( ei) == node_in) then
        ! This edge is owned by node_in

        vi  = mesh%EV( ei,1)
        vj  = mesh%EV( ei,2)
        vil = mesh%EV( ei,3)
        vir = mesh%EV( ei,4)

        til = mesh%ETri( ei,1)
        tir = mesh%ETri( ei,2)

        has_vertices_in_node_next_to  = .false.
        has_triangles_in_node_next_to = .false.
        has_edges_in_node_next_to     = .false.

        ! Check if this edge has any neighbouring vertices owned by node_next_to
        if (mesh%V_owning_node( vi) == node_next_to) has_vertices_in_node_next_to = .true.
        if (mesh%V_owning_node( vj) == node_next_to) has_vertices_in_node_next_to = .true.
        if (vil > 0) then
          if (mesh%V_owning_node( vil) == node_next_to) has_vertices_in_node_next_to = .true.
        end if
        if (vir > 0) then
          if (mesh%V_owning_node( vir) == node_next_to) has_vertices_in_node_next_to = .true.
        end if

        ! ! Check if this edge has any neighbouring triangles owned by node_next_to
        if (til > 0) then
          if (mesh%Tri_owning_node( til) == node_next_to) has_triangles_in_node_next_to = .true.
        end if
        if (tir > 0) then
          if (mesh%Tri_owning_node( tir) == node_next_to) has_triangles_in_node_next_to = .true.
        end if

        ! ! Check if this edge has any neighbouring edges owned by node_next_to
        ! FIXME!

        if (has_vertices_in_node_next_to .or. &
            has_triangles_in_node_next_to .or. &
            has_edges_in_node_next_to) then
          ei1 = min( ei1, ei)
          ei2 = max( ei2, ei)
        end if

      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_halo_range_c

  subroutine print_parallelisation_info( mesh)

    ! In/output variables:
    type(type_mesh), intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_parallelisation_info'
    integer, dimension(0:par%n-1) :: node, process, &
      vi1, vi2, vi1_node, vi2_node, ti1, ti2, ti1_node, ti2_node, ei1, ei2, ei1_node, ei2_node, &
      vi1_nih, vi2_nih, vi1_hle, vi2_hle, nV_hle, vi1_hli, vi2_hli, nV_hli, vi1_hre, vi2_hre, nV_hre, vi1_hri, vi2_hri, nV_hri, &
      ti1_nih, ti2_nih, ti1_hle, ti2_hle, nTri_hle, ti1_hli, ti2_hli, nTri_hli, ti1_hre, ti2_hre, nTri_hre, ti1_hri, ti2_hri, nTri_hri, &
      ei1_nih, ei2_nih, ei1_hle, ei2_hle, nE_hle, ei1_hli, ei2_hli, nE_hli, ei1_hre, ei2_hre, nE_hre, ei1_hri, ei2_hri, nE_hri
    integer :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call MPI_ALLGATHER( par%node_ID  , 1, MPI_INTEGER, node    , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( par%i        , 1, MPI_INTEGER, process , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( mesh%vi1     , 1, MPI_INTEGER, vi1     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%vi2     , 1, MPI_INTEGER, vi2     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( mesh%pai_V%i1_node, 1, MPI_INTEGER, vi1_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i2_node, 1, MPI_INTEGER, vi2_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i1_nih , 1, MPI_INTEGER, vi1_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i2_nih , 1, MPI_INTEGER, vi2_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%n_hle  , 1, MPI_INTEGER, nV_hle  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i1_hle , 1, MPI_INTEGER, vi1_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i2_hle , 1, MPI_INTEGER, vi2_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%n_hli  , 1, MPI_INTEGER, nV_hli  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i1_hli , 1, MPI_INTEGER, vi1_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i2_hli , 1, MPI_INTEGER, vi2_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%n_hre  , 1, MPI_INTEGER, nV_hre  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i1_hre , 1, MPI_INTEGER, vi1_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i2_hre , 1, MPI_INTEGER, vi2_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%n_hri  , 1, MPI_INTEGER, nV_hri  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i1_hri , 1, MPI_INTEGER, vi1_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_V%i2_hri , 1, MPI_INTEGER, vi2_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( mesh%pai_Tri%i1     , 1, MPI_INTEGER, ti1     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i2     , 1, MPI_INTEGER, ti2     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i1_node, 1, MPI_INTEGER, ti1_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i2_node, 1, MPI_INTEGER, ti2_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i1_nih , 1, MPI_INTEGER, ti1_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i2_nih , 1, MPI_INTEGER, ti2_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%n_hle  , 1, MPI_INTEGER, nTri_hle, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i1_hle , 1, MPI_INTEGER, ti1_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i2_hle , 1, MPI_INTEGER, ti2_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%n_hli  , 1, MPI_INTEGER, nTri_hli, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i1_hli , 1, MPI_INTEGER, ti1_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i2_hli , 1, MPI_INTEGER, ti2_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%n_hre  , 1, MPI_INTEGER, nTri_hre, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i1_hre , 1, MPI_INTEGER, ti1_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i2_hre , 1, MPI_INTEGER, ti2_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%n_hri  , 1, MPI_INTEGER, nTri_hri, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i1_hri , 1, MPI_INTEGER, ti1_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_Tri%i2_hri , 1, MPI_INTEGER, ti2_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( mesh%pai_E%i1     , 1, MPI_INTEGER, ei1     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i2     , 1, MPI_INTEGER, ei2     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i1_node, 1, MPI_INTEGER, ei1_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i2_node, 1, MPI_INTEGER, ei2_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i1_nih , 1, MPI_INTEGER, ei1_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i2_nih , 1, MPI_INTEGER, ei2_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%n_hle  , 1, MPI_INTEGER, nE_hle  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i1_hle , 1, MPI_INTEGER, ei1_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i2_hle , 1, MPI_INTEGER, ei2_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%n_hli  , 1, MPI_INTEGER, nE_hli  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i1_hli , 1, MPI_INTEGER, ei1_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i2_hli , 1, MPI_INTEGER, ei2_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%n_hre  , 1, MPI_INTEGER, nE_hre  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i1_hre , 1, MPI_INTEGER, ei1_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i2_hre , 1, MPI_INTEGER, ei2_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%n_hri  , 1, MPI_INTEGER, nE_hri  , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i1_hri , 1, MPI_INTEGER, ei1_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( mesh%pai_E%i2_hri , 1, MPI_INTEGER, ei2_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if (par%primary) then
      write(0,'(A,7I8)') 'Node     :', node
      write(0,'(A,7I8)') 'Process  :', process
      write(0,*) ''
      write(0,*) '-= VERTICES =-'
      write(0,*) ''
      write(0,'(A,I8)')  'nV = ', mesh%nV
      write(0,*) ''
      write(0,'(A,7I8)') 'vi1      :', vi1
      write(0,'(A,7I8)') 'vi2      :', vi2
      write(0,'(A,7I8)') 'vi1_node :', vi1_node
      write(0,'(A,7I8)') 'vi2_node :', vi2_node
      write(0,'(A,7I8)') 'vi1_nih  :', vi1_nih
      write(0,'(A,7I8)') 'vi2_nih  :', vi2_nih
      write(0,*) ''
      write(0,'(A,7I8)') 'nV_hle   :', nV_hle
      write(0,'(A,7I8)') 'vi1_hle  :', vi1_hle
      write(0,'(A,7I8)') 'vi2_hle  :', vi2_hle
      write(0,*) ''
      write(0,'(A,7I8)') 'nV_hli   :', nV_hli
      write(0,'(A,7I8)') 'vi1_hli  :', vi1_hli
      write(0,'(A,7I8)') 'vi2_hli  :', vi2_hli
      write(0,*) ''
      write(0,'(A,7I8)') 'nV_hre   :', nV_hre
      write(0,'(A,7I8)') 'vi1_hre  :', vi1_hre
      write(0,'(A,7I8)') 'vi2_hre  :', vi2_hre
      write(0,*) ''
      write(0,'(A,7I8)') 'nV_hri   :', nV_hri
      write(0,'(A,7I8)') 'vi1_hri  :', vi1_hri
      write(0,'(A,7I8)') 'vi2_hri  :', vi2_hri
      write(0,*) ''
      write(0,*) '-= TRIANGLES =-'
      write(0,*) ''
      write(0,'(A,I8)')  'nTri = ', mesh%nTri
      write(0,*) ''
      write(0,'(A,7I8)') 'ti1      :', ti1
      write(0,'(A,7I8)') 'ti2      :', ti2
      write(0,'(A,7I8)') 'ti1_node :', ti1_node
      write(0,'(A,7I8)') 'ti2_node :', ti2_node
      write(0,'(A,7I8)') 'ti1_nih  :', ti1_nih
      write(0,'(A,7I8)') 'ti2_nih  :', ti2_nih
      write(0,*) ''
      write(0,'(A,7I8)') 'nTri_hle :', nTri_hle
      write(0,'(A,7I8)') 'ti1_hle  :', ti1_hle
      write(0,'(A,7I8)') 'ti2_hle  :', ti2_hle
      write(0,*) ''
      write(0,'(A,7I8)') 'nTri_hli :', nTri_hli
      write(0,'(A,7I8)') 'ti1_hli  :', ti1_hli
      write(0,'(A,7I8)') 'ti2_hli  :', ti2_hli
      write(0,*) ''
      write(0,'(A,7I8)') 'nTri_hre :', nTri_hre
      write(0,'(A,7I8)') 'ti1_hre  :', ti1_hre
      write(0,'(A,7I8)') 'ti2_hre  :', ti2_hre
      write(0,*) ''
      write(0,'(A,7I8)') 'nTri_hri :', nTri_hri
      write(0,'(A,7I8)') 'ti1_hri  :', ti1_hri
      write(0,'(A,7I8)') 'ti2_hri  :', ti2_hri
      write(0,*) ''
      write(0,*) '-= EDGES =-'
      write(0,*) ''
      write(0,'(A,I8)')  'nE = ', mesh%nE
      write(0,*) ''
      write(0,'(A,7I8)') 'ei1      :', ei1
      write(0,'(A,7I8)') 'ei2      :', ei2
      write(0,'(A,7I8)') 'ei1_node :', ei1_node
      write(0,'(A,7I8)') 'ei2_node :', ei2_node
      write(0,'(A,7I8)') 'ei1_nih  :', ei1_nih
      write(0,'(A,7I8)') 'ei2_nih  :', ei2_nih
      write(0,*) ''
      write(0,'(A,7I8)') 'nE_hle   :', nE_hle
      write(0,'(A,7I8)') 'ei1_hle  :', ei1_hle
      write(0,'(A,7I8)') 'ei2_hle  :', ei2_hle
      write(0,*) ''
      write(0,'(A,7I8)') 'nE_hli   :', nE_hli
      write(0,'(A,7I8)') 'ei1_hli  :', ei1_hli
      write(0,'(A,7I8)') 'ei2_hli  :', ei2_hli
      write(0,*) ''
      write(0,'(A,7I8)') 'nE_hre   :', nE_hre
      write(0,'(A,7I8)') 'ei1_hre  :', ei1_hre
      write(0,'(A,7I8)') 'ei2_hre  :', ei2_hre
      write(0,*) ''
      write(0,'(A,7I8)') 'nE_hri   :', nE_hri
      write(0,'(A,7I8)') 'ei1_hri  :', ei1_hri
      write(0,'(A,7I8)') 'ei2_hri  :', ei2_hri
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine print_parallelisation_info

end module mesh_parallelisation
