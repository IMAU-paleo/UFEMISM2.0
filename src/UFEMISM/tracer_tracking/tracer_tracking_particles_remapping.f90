module tracer_tracking_model_particles_remapping

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use mpi_f08, only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INT, MPI_IRECV, MPI_REQUEST, &
    MPI_RSEND, MPI_STATUS_IGNORE, MPI_WAIT
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles, type_map_particles_to_mesh
  use model_configuration, only: C
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex, &
    interpolate_to_point_dp_2D, interpolate_to_point_dp_3D
  use reallocate_mod, only: reallocate
  use netcdf, only: NF90_UNLIMITED, NF90_INT64, NF90_DOUBLE
  use netcdf_io_main
  use mpi_distributed_memory, only: partition_list, gather_to_all

  implicit none

  private

  public :: calc_particles_to_mesh_map, map_tracer_to_mesh

contains

  subroutine map_tracer_to_mesh( mesh, particles, f_particles, f_mesh)
    !< Map a data field defined on the tracer particles to the model mesh

    ! In/output variables
    type(type_mesh),                             intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles),  intent(in   ) :: particles
    real(dp), dimension(particles%n_max),        intent(in   ) :: f_particles
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz), intent(  out) :: f_mesh

    ! Local variables
    character(len=1024), parameter               :: routine_name = 'map_tracer_to_mesh'
    real(dp), dimension(particles%n_max * par%n) :: f_particles_tot
    integer                                      :: vi, k, ii, ip
    real(dp)                                     :: w_sum, f_sum, w, dist
    logical                                      :: coincides_with_particle

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather data from all particles across the processes
    call gather_to_all( f_particles, f_particles_tot)

    f_mesh = 0._dp

    do vi = mesh%vi1, mesh%vi2

      do k = 1, C%nz

        w_sum = 0._dp
        f_sum = 0._dp
        coincides_with_particle = .false.

        do ii = 1, particles%map%n
          ip   = particles%map%ip( vi,k,ii)
          dist = particles%map%d(  vi,k,ii)
          if (dist < mesh%tol_dist) then
            coincides_with_particle = .true.
            f_mesh( vi,k) = f_particles_tot( ip)
            exit
          end if
          w = 1._dp / dist**2
          w_sum = w_sum + w
          f_sum = f_sum + w * f_particles_tot( ip)
        end do
        if (.not. coincides_with_particle) then
          f_mesh( vi,k) = f_sum / w_sum
        end if

      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_tracer_to_mesh

  subroutine calc_particles_to_mesh_map( mesh, particles)
    !< Calculate the interpolation weights describing the
    !< mapping operation from the tracer particles to the model mesh.

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables
    character(len=1024), parameter   :: routine_name = 'calc_particles_to_mesh_map'
    integer                          :: p
    type(type_map_particles_to_mesh) :: map_other

    ! Add routine to path
    call init_routine( routine_name)

    call calc_particles_to_mesh_map_for_this_process( mesh, particles)

    do p = 1, par%n-1
      call exchange_map_data_between_processes( mesh, particles, p, map_other)
      call combine_received_map_data_with_own( mesh, particles, map_other)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_particles_to_mesh_map

  subroutine calc_particles_to_mesh_map_for_this_process( mesh, particles)
    !< Calculate the interpolation weights describing the
    !< mapping operation from this process' particles to the model mesh.

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables
    character(len=1024), parameter          :: routine_name = 'calc_particles_to_mesh_map_for_this_process'
    integer                                 :: ip_offset
    real(dp), dimension(mesh%nV,C%nz,3)     :: rs_mesh
    real(dp), dimension(particles%n_max ,3) :: rs_particles
    integer                                 :: vi, k, ip
    real(dp)                                :: dist_max
    integer                                 :: vi_in
    integer,  dimension(mesh%nV)            :: map, stack
    integer                                 :: stackN
    integer                                 :: ii, jj, ci, vj
    integer                                 :: n_cycles
    logical                                 :: is_one_of_n_nearest_particles
    real(dp)                                :: dist

    ! Add routine to path
    call init_routine( routine_name)

    ! Offset particle index (to assign unique indices to particles across the processes)
    ip_offset = par%i * particles%n_max

    ! Calculate scaled coordinates for mesh vertex-layers and for particles
    call calc_scaled_coordinates( mesh, particles, rs_mesh, rs_particles)

    ! Initialise map data
    if (allocated( particles%map%ip)) deallocate( particles%map%ip)
    if (allocated( particles%map%d )) deallocate( particles%map%d )
    allocate( particles%map%ip( mesh%nV, C%nz, particles%map%n))
    allocate( particles%map%d ( mesh%nV, C%nz, particles%map%n))

    dist_max = norm2( [mesh%xmin,mesh%ymin] - [mesh%xmax,mesh%ymax])
    particles%map%ip = 0
    particles%map%d  = dist_max

    do ip = 1, particles%n_max

      if (.not. particles%is_in_use( ip)) cycle

      vi_in = particles%vi_in( ip)

      ! Clean up map and stack
      map    = 0
      stack  = 0
      stackN = 0

      ! Initialise map and stack with the vertex whose Voronoi cell contains particle ip
      map( vi_in) = 1
      stackN = 1
      stack( 1) = vi_in

      ! Expand outward until we've covered all vertices that are nearest to particle ip
      n_cycles = 0
      do while (stackN > 0)

        ! Safety
        n_cycles = n_cycles + 1
        if (n_cycles > mesh%nV) call crash('Flood-fill expansion got stuck!')

        ! Take the last vertex from the stack
        vi = stack( stackN)
        stackN = stackN - 1
        ! Mark it as checked
        map( vi) = 2

        is_one_of_n_nearest_particles = .false.

        do k = 1, C%nz
          ! Calculate the scaled distance between particle ip and vertex-layer [vi,k]
          dist = norm2( rs_particles( ip,:) - rs_mesh( vi,k,:))
          ! If this particles is closer to [vi,k] than listed nearest particle ii,
          ! move ii:end down a slot and insert this particle in between.
          do ii = 1, particles%map%n
            if (dist < particles%map%d( vi,k,ii)) then
              is_one_of_n_nearest_particles = .true.
              do jj = particles%map%n, ii+1, -1
                particles%map%d(  vi,k,jj) = particles%map%d(  vi,k,jj-1)
                particles%map%ip( vi,k,jj) = particles%map%ip( vi,k,jj-1)
              end do
              particles%map%d(  vi,k,ii) = dist
              particles%map%ip( vi,k,ii) = ip + ip_offset
              exit
            end if
          end do
        end do

        ! If this particle was nearest to any of the layers of this vertex, add
        ! this vertex' unchecked neighbours to the stack
        if (is_one_of_n_nearest_particles) then
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (map( vj) == 0) then
              ! Add this unchecked neighbour to the stack
              stackN = stackN + 1
              stack( stackN) = vj
              map( vj) = 1
            end if
          end do
        end if

      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_particles_to_mesh_map_for_this_process

  subroutine calc_scaled_coordinates( mesh, particles, rs_mesh, rs_particles)

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(in   ) :: particles
    real(dp), dimension(mesh%nV,C%nz,3),        intent(  out) :: rs_mesh
    real(dp), dimension(particles%n_max ,3),    intent(  out) :: rs_particles

    ! Local variables
    character(len=1024), parameter :: routine_name = 'calc_scaled_coordinates'
    integer                        :: vi, k, ip

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate scaled coordinates for mesh vertex-layers and for particles
    do vi = 1, mesh%nV
      do k = 1, C%nz
        rs_mesh( vi,k,1) = (mesh%V( vi,1) - mesh%xmin) / (mesh%xmax - mesh%xmin)
        rs_mesh( vi,k,2) = (mesh%V( vi,2) - mesh%ymin) / (mesh%ymax - mesh%ymin)
        rs_mesh( vi,k,3) = mesh%zeta( k)
      end do
    end do

    do ip = 1, particles%n_max
      if (.not. particles%is_in_use( ip)) cycle
      rs_particles( ip,1) = (particles%r( ip,1) - mesh%xmin) / (mesh%xmax - mesh%xmin)
      rs_particles( ip,2) = (particles%r( ip,2) - mesh%ymin) / (mesh%ymax - mesh%ymin)
      rs_particles( ip,3) = particles%zeta( ip)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_scaled_coordinates

  subroutine exchange_map_data_between_processes( mesh, particles, p, map_other)
    !< Send map data for their vertices to process par%i+p,
    !< receive data for your vertices from process par%i-p)

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(in   ) :: particles
    integer,                                    intent(in   ) :: p
    type(type_map_particles_to_mesh),           intent(  out) :: map_other

    ! Local variables
    character(len=1024), parameter          :: routine_name = 'exchange_map_data_between_processes'
    integer                                 :: par_i_send, par_i_recv
    integer                                 :: vi1_send, vi2_send, nV_send, nel_send, nel_recv
    integer,  dimension(:,:,:), allocatable :: ip_send, ip_recv
    real(dp), dimension(:,:,:), allocatable :: d_send, d_recv
    integer                                 :: ierr
    type(MPI_REQUEST)                       :: req

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine which processes we want to exchange data with
    par_i_send = par%i + p
    if (par_i_send > par%n-1) par_i_send = par_i_send - par%n

    par_i_recv = par%i - p
    if (par_i_recv < 0) par_i_recv = par_i_recv + par%n

    ! Determine the vertices owned by the process we want to send to
    call partition_list( mesh%nV, par_i_send, par%n, vi1_send, vi2_send)
    nV_send = vi2_send + 1 - vi1_send

    nel_send = nV_send     * C%nz * particles%map%n
    nel_recv = mesh%nV_loc * C%nz * particles%map%n

    allocate( ip_send( vi1_send:vi2_send, C%nz, particles%map%n))
    allocate( d_send ( vi1_send:vi2_send, C%nz, particles%map%n))
    allocate( ip_recv( mesh%vi1:mesh%vi2, C%nz, particles%map%n))
    allocate( d_recv ( mesh%vi1:mesh%vi2, C%nz, particles%map%n))

    ! Collect data to send
    ip_send = particles%map%ip( vi1_send:vi2_send,:,:)
    d_send  = particles%map%d ( vi1_send:vi2_send,:,:)

    ! Use ready-mode send and non-blocking receive to make sure that
    ! all processes can send/receive at the same time

    call MPI_IRECV( ip_recv, nel_recv, MPI_INT, par_i_recv, 13, MPI_COMM_WORLD, req, ierr)
    call MPI_RSEND( ip_send, nel_send, MPI_INT, par_i_send, 13, MPI_COMM_WORLD,      ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE, ierr)

    call MPI_IRECV( d_recv, nel_recv, MPI_DOUBLE_PRECISION, par_i_recv, 13, MPI_COMM_WORLD, req, ierr)
    call MPI_RSEND( d_send, nel_send, MPI_DOUBLE_PRECISION, par_i_send, 13, MPI_COMM_WORLD,      ierr)
    call MPI_WAIT( req, MPI_STATUS_IGNORE, ierr)

    ! Collect received data into map_other
    map_other%n = particles%map%n
    allocate( map_other%ip( mesh%vi1:mesh%vi2, C%nz, map_other%n))
    allocate( map_other%d ( mesh%vi1:mesh%vi2, C%nz, map_other%n))
    map_other%ip = ip_recv
    map_other%d  = d_recv

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_map_data_between_processes

  subroutine combine_received_map_data_with_own( mesh, particles, map_other)

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    type(type_map_particles_to_mesh),           intent(in   ) :: map_other

    ! Local variables
    character(len=1024), parameter :: routine_name = 'combine_received_map_data_with_own'
    integer                        :: vi,k,ii_other,ip_other,ii,jj
    real(dp)                       :: dist_other

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
    do k = 1, C%nz

      ! For all n particles in map_other, check if they are closer than
      ! any of the n particles in this process' own map. If so, insert them

      do ii_other = 1, map_other%n

        ip_other   = map_other%ip( vi,k,ii_other)
        dist_other = map_other%d ( vi,k,ii_other)

        ! If this particle from the other map is closer to [vi,k] than listed nearest particle ii,
        ! move ii:end down a slot and insert this particle in between.
        do ii = 1, particles%map%n
          if (dist_other < particles%map%d( vi,k,ii)) then
            do jj = particles%map%n, ii+1, -1
              particles%map%d(  vi,k,jj) = particles%map%d(  vi,k,jj-1)
              particles%map%ip( vi,k,jj) = particles%map%ip( vi,k,jj-1)
            end do
            particles%map%d(  vi,k,ii) = dist_other
            particles%map%ip( vi,k,ii) = ip_other
            exit
          end if
        end do

      end do

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine combine_received_map_data_with_own

end module tracer_tracking_model_particles_remapping