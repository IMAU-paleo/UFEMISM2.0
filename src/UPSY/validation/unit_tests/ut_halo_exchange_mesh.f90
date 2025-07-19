module ut_halo_exchange_mesh

  ! Unit tests for mesh halo exchange code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mesh_halo_exchange, only: exchange_halos
  use mesh_types, only: type_mesh
  use parameters, only: pi
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_contiguous_domains, only: enforce_contiguous_process_domains

  implicit none

  private

  public :: test_mesh_halo_exchange_main

contains

  subroutine test_mesh_halo_exchange_main( test_name_parent)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_mesh_halo_exchange_main'
    character(len=1024), parameter :: test_name_local = 'mesh_halo_exchange'
    character(len=1024)            :: test_name
    real(dp)                       :: xmin, xmax, ymin, ymax, alpha_min, res_max
    character(len=1024)            :: name
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    name = 'test_mesh'
    xmin = 0._dp
    xmax = pi
    ymin = 0._dp
    ymax = pi
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 20._dp

    call allocate_mesh_primary( mesh, name, 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call enforce_contiguous_process_domains( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    ! Test halo exchange on this mesh
    call test_halo_exchange_mesh_a_logical( test_name, mesh)
    call test_halo_exchange_mesh_a_int    ( test_name, mesh)
    call test_halo_exchange_mesh_a_int_3D ( test_name, mesh)
    call test_halo_exchange_mesh_a_dp     ( test_name, mesh)
    call test_halo_exchange_mesh_a_dp_3D  ( test_name, mesh)

    call test_halo_exchange_mesh_b_logical( test_name, mesh)
    call test_halo_exchange_mesh_b_int    ( test_name, mesh)
    call test_halo_exchange_mesh_b_int_3D ( test_name, mesh)
    call test_halo_exchange_mesh_b_dp     ( test_name, mesh)
    call test_halo_exchange_mesh_b_dp_3D  ( test_name, mesh)

    call test_halo_exchange_mesh_c_logical( test_name, mesh)
    call test_halo_exchange_mesh_c_int    ( test_name, mesh)
    call test_halo_exchange_mesh_c_int_3D ( test_name, mesh)
    call test_halo_exchange_mesh_c_dp     ( test_name, mesh)
    call test_halo_exchange_mesh_c_dp_3D  ( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_halo_exchange_main

  ! == a-grid

  subroutine test_halo_exchange_mesh_a_logical( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_halo_exchange_mesh_a_logical'
    character(len=1024), parameter :: test_name_local = 'a/logical'
    character(len=1024)            :: test_name
    logical, dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                  :: wd_nih
    integer                        :: vi
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_V%n_nih)
    d_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = .false.
      do vi = mesh%pai_V%i1_node, mesh%pai_V%i2_node
        d_nih( vi) = .true.
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do vi = mesh%pai_V%i1_nih, mesh%pai_V%i2_nih
      test_result = test_result .and. d_nih( vi)
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_a_logical

  subroutine test_halo_exchange_mesh_a_int( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_halo_exchange_mesh_a_int'
    character(len=1024), parameter :: test_name_local = 'a/int'
    character(len=1024)            :: test_name
    integer, dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                  :: wd_nih
    integer                        :: vi
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_V%n_nih)
    d_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do vi = mesh%pai_V%i1_node, mesh%pai_V%i2_node
        d_nih( vi) = 1
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do vi = mesh%pai_V%i1_nih, mesh%pai_V%i2_nih
      test_result = test_result .and. d_nih( vi) == 1
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_a_int

  subroutine test_halo_exchange_mesh_a_int_3D( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_halo_exchange_mesh_a_int_3D'
    character(len=1024), parameter               :: test_name_local = 'a/int_3D'
    character(len=1024)                          :: test_name
    integer, dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                :: wd_nih
    integer                                      :: nz,vi,k
    logical                                      :: test_result
    integer                                      :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    nz = 13
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_V%n_nih, nz)
    d_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih,1:nz) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do vi = mesh%pai_V%i1_node, mesh%pai_V%i2_node
        do k = 1, nz
          d_nih( vi,k) = k
        end do
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do vi = mesh%pai_V%i1_nih, mesh%pai_V%i2_nih
      do k = 1, nz
        test_result = test_result .and. d_nih( vi,k) == k
      end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_a_int_3D

  subroutine test_halo_exchange_mesh_a_dp( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_halo_exchange_mesh_a_dp'
    character(len=1024), parameter  :: test_name_local = 'a/dp'
    character(len=1024)             :: test_name
    real(dp), dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                   :: wd_nih
    integer                         :: vi
    logical                         :: test_result
    integer                         :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_V%n_nih)
    d_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do vi = mesh%pai_V%i1_node, mesh%pai_V%i2_node
        d_nih( vi) = 1._dp
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do vi = mesh%pai_V%i1_nih, mesh%pai_V%i2_nih
      test_result = test_result .and. d_nih( vi) == 1._dp
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_a_dp

  subroutine test_halo_exchange_mesh_a_dp_3D( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_halo_exchange_mesh_a_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'a/dp_3D'
    character(len=1024)                           :: test_name
    real(dp), dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                 :: wd_nih
    integer                                       :: nz,vi,k
    logical                                       :: test_result
    integer                                       :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    nz = 13
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_V%n_nih, nz)
    d_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih,1:nz) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do vi = mesh%pai_V%i1_node, mesh%pai_V%i2_node
        do k = 1, nz
          d_nih( vi,k) = real( k,dp)
        end do
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do vi = mesh%pai_V%i1_nih, mesh%pai_V%i2_nih
      do k = 1, nz
        test_result = test_result .and. d_nih( vi,k) == real( k,dp)
      end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_a_dp_3D

  ! == b-grid

  subroutine test_halo_exchange_mesh_b_logical( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_halo_exchange_mesh_b_logical'
    character(len=1024), parameter :: test_name_local = 'b/logical'
    character(len=1024)            :: test_name
    logical, dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                  :: wd_nih
    integer                        :: ti
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_Tri%n_nih)
    d_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = .false.
      do ti = mesh%pai_Tri%i1_node, mesh%pai_Tri%i2_node
        d_nih( ti) = .true.
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ti = mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih
      test_result = test_result .and. d_nih( ti)
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_b_logical

  subroutine test_halo_exchange_mesh_b_int( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_halo_exchange_mesh_b_int'
    character(len=1024), parameter :: test_name_local = 'b/int'
    character(len=1024)            :: test_name
    integer, dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                  :: wd_nih
    integer                        :: ti
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_Tri%n_nih)
    d_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ti = mesh%pai_Tri%i1_node, mesh%pai_Tri%i2_node
        d_nih( ti) = 1
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ti = mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih
      test_result = test_result .and. d_nih( ti) == 1
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_b_int

  subroutine test_halo_exchange_mesh_b_int_3D( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_halo_exchange_mesh_b_int_3D'
    character(len=1024), parameter               :: test_name_local = 'b/int_3D'
    character(len=1024)                          :: test_name
    integer, dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                :: wd_nih
    integer                                      :: nz,ti,k
    logical                                      :: test_result
    integer                                      :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    nz = 13
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_Tri%n_nih, nz)
    d_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih,1:nz) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ti = mesh%pai_Tri%i1_node, mesh%pai_Tri%i2_node
        do k = 1, nz
          d_nih( ti,k) = k
        end do
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ti = mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih
      do k = 1, nz
        test_result = test_result .and. d_nih( ti,k) == k
      end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_b_int_3D

  subroutine test_halo_exchange_mesh_b_dp( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_halo_exchange_mesh_b_dp'
    character(len=1024), parameter  :: test_name_local = 'b/dp'
    character(len=1024)             :: test_name
    real(dp), dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                   :: wd_nih
    integer                         :: ti
    logical                         :: test_result
    integer                         :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_Tri%n_nih)
    d_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ti = mesh%pai_Tri%i1_node, mesh%pai_Tri%i2_node
        d_nih( ti) = 1._dp
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ti = mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih
      test_result = test_result .and. d_nih( ti) == 1._dp
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_b_dp

  subroutine test_halo_exchange_mesh_b_dp_3D( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_halo_exchange_mesh_b_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'b/dp_3D'
    character(len=1024)                           :: test_name
    real(dp), dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                 :: wd_nih
    integer                                       :: nz,ti,k
    logical                                       :: test_result
    integer                                       :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    nz = 13
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_Tri%n_nih, nz)
    d_nih( mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih,1:nz) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ti = mesh%pai_Tri%i1_node, mesh%pai_Tri%i2_node
        do k = 1, nz
          d_nih( ti,k) = real( k,dp)
        end do
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ti = mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih
      do k = 1, nz
        test_result = test_result .and. d_nih( ti,k) == real( k,dp)
      end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_b_dp_3D

  ! == c-grid

  subroutine test_halo_exchange_mesh_c_logical( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_halo_exchange_mesh_c_logical'
    character(len=1024), parameter :: test_name_local = 'c/logical'
    character(len=1024)            :: test_name
    logical, dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                  :: wd_nih
    integer                        :: ei
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_E%n_nih)
    d_nih( mesh%pai_E%i1_nih:mesh%pai_E%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = .false.
      do ei = mesh%pai_E%i1_node, mesh%pai_E%i2_node
        d_nih( ei) = .true.
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ei = mesh%pai_E%i1_nih, mesh%pai_E%i2_nih
      test_result = test_result .and. d_nih( ei)
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_c_logical

  subroutine test_halo_exchange_mesh_c_int( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_halo_exchange_mesh_c_int'
    character(len=1024), parameter :: test_name_local = 'c/int'
    character(len=1024)            :: test_name
    integer, dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                  :: wd_nih
    integer                        :: ei
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_E%n_nih)
    d_nih( mesh%pai_E%i1_nih:mesh%pai_E%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ei = mesh%pai_E%i1_node, mesh%pai_E%i2_node
        d_nih( ei) = 1
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ei = mesh%pai_E%i1_nih, mesh%pai_E%i2_nih
      test_result = test_result .and. d_nih( ei) == 1
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_c_int

  subroutine test_halo_exchange_mesh_c_int_3D( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_halo_exchange_mesh_c_int_3D'
    character(len=1024), parameter               :: test_name_local = 'c/int_3D'
    character(len=1024)                          :: test_name
    integer, dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                :: wd_nih
    integer                                      :: nz,ei,k
    logical                                      :: test_result
    integer                                      :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    nz = 13
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_E%n_nih, nz)
    d_nih( mesh%pai_E%i1_nih:mesh%pai_E%i2_nih,1:nz) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ei = mesh%pai_E%i1_node, mesh%pai_E%i2_node
        do k = 1, nz
          d_nih( ei,k) = k
        end do
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ei = mesh%pai_E%i1_nih, mesh%pai_E%i2_nih
      do k = 1, nz
        test_result = test_result .and. d_nih( ei,k) == k
      end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_c_int_3D

  subroutine test_halo_exchange_mesh_c_dp( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_halo_exchange_mesh_c_dp'
    character(len=1024), parameter  :: test_name_local = 'c/dp'
    character(len=1024)             :: test_name
    real(dp), dimension(:), pointer :: d_nih => null()
    type(MPI_WIN)                   :: wd_nih
    integer                         :: ei
    logical                         :: test_result
    integer                         :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_E%n_nih)
    d_nih( mesh%pai_E%i1_nih:mesh%pai_E%i2_nih) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ei = mesh%pai_E%i1_node, mesh%pai_E%i2_node
        d_nih( ei) = 1._dp
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ei = mesh%pai_E%i1_nih, mesh%pai_E%i2_nih
      test_result = test_result .and. d_nih( ei) == 1._dp
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_c_dp

  subroutine test_halo_exchange_mesh_c_dp_3D( test_name_parent, mesh)
    ! Test the mesh halo exchange subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_halo_exchange_mesh_c_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'c/dp_3D'
    character(len=1024)                           :: test_name
    real(dp), dimension(:,:), contiguous, pointer :: d_nih => null()
    type(MPI_WIN)                                 :: wd_nih
    integer                                       :: nz,ei,k
    logical                                       :: test_result
    integer                                       :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on 7 cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate some node-shared memory
    nz = 13
    call allocate_dist_shared( d_nih, wd_nih, mesh%pai_E%n_nih, nz)
    d_nih( mesh%pai_E%i1_nih:mesh%pai_E%i2_nih,1:nz) => d_nih

    ! Fill in data
    if (par%node_primary) then
      d_nih = 0
      do ei = mesh%pai_E%i1_node, mesh%pai_E%i2_node
        do k = 1, nz
          d_nih( ei,k) = real( k,dp)
        end do
      end do
    end if
    call sync_node

    ! Exchange halos
    call exchange_halos( mesh, d_nih)

    ! Verify results
    test_result = .true.
    do ei = mesh%pai_E%i1_nih, mesh%pai_E%i2_nih
      do k = 1, nz
        test_result = test_result .and. d_nih( ei,k) == real( k,dp)
      end do
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call deallocate_dist_shared( d_nih, wd_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_halo_exchange_mesh_c_dp_3D

end module ut_halo_exchange_mesh
