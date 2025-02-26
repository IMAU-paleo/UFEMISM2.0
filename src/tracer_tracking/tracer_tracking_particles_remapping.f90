module tracer_tracking_model_particles_remapping

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par
  use mpi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles
  use model_configuration, only: C
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex, &
    interpolate_to_point_dp_2D, interpolate_to_point_dp_3D
  use reallocate_mod, only: reallocate
  use netcdf, only: NF90_UNLIMITED, NF90_INT64, NF90_DOUBLE
  use netcdf_io_main

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
    character(len=1024), parameter :: routine_name = 'map_tracer_to_mesh'
    integer                        :: vi, k, ii, ip
    real(dp)                       :: w_sum, f_sum, w, dist
    logical                        :: coincides_with_particle

    ! Add routine to path
    call init_routine( routine_name)

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
            f_mesh( vi,k) = f_particles( ip)
            exit
          end if
          w = 1._dp / dist**2
          w_sum = w_sum + w
          f_sum = f_sum + w * f_particles( ip)
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
    character(len=1024), parameter      :: routine_name = 'calc_particles_to_mesh_map'
    real(dp), dimension(mesh%nV,C%nz,3) :: rs_mesh
    real(dp), dimension(particles%n_max ,3) :: rs_particles
    integer                             :: vi, k, ip
    real(dp)                            :: dist_max
    integer                             :: vi_in
    integer,  dimension(mesh%nV)        :: map, stack
    integer                             :: stackN
    integer                             :: ii, jj, ci, vj
    integer                             :: n_cycles
    logical                             :: is_one_of_n_nearest_particles
    real(dp)                            :: dist

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    dist_max = norm2( [mesh%xmin,mesh%ymin] - [mesh%xmax,mesh%ymax])
    particles%map%ip = 0
    particles%map%d  = dist_max

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

    ! Initialise map and stack
    map    = 0
    stack  = 0
    stackN = 0

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
              particles%map%ip( vi,k,ii) = ip
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

  end subroutine calc_particles_to_mesh_map

end module tracer_tracking_model_particles_remapping