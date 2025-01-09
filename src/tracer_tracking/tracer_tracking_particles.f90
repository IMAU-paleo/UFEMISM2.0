module tracer_tracking_model_particles

  ! The main tracer tracking model module.

  use tests_main
  use assertions_basic
  use precisions, only: dp, int8
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
  use netcdf_output, only: write_time_to_file
  use netcdf_field_name_options
  use netcdf_basic_wrappers
  use netcdf_write_var_master

  implicit none

  private

  public :: initialise_tracer_tracking_model_particles, calc_particle_zeta, &
    interpolate_3d_velocities_to_3D_point, calc_particles_to_mesh_map, add_particle, &
    move_and_remove_particle, create_particles_netcdf_file, write_to_netcdf_file

  integer, parameter :: n_tracers         = 1
  integer, parameter :: n_nearest_to_find = 4

contains

  subroutine initialise_tracer_tracking_model_particles( mesh, ice, particles, n_max)
    !< Initialise the particle-based tracer-tracking model.
    !< Allocates memory for a number of particles.

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(in   ) :: ice
    type(type_tracer_tracking_model_particles), intent(  out) :: particles
    integer                                   , intent(in   ) :: n_max

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_tracer_tracking_model_particles'

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%master)  write(*,'(a)') '     Initialising particle-based tracer tracking model...'

    particles%n_max  = n_max
    particles%id_max = 0_int8
    allocate( particles%is_in_use( particles%n_max           ), source = .false.)
    allocate( particles%id       ( particles%n_max           ), source = 0_int8 )
    allocate( particles%r        ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%zeta     ( particles%n_max           ), source = 0._dp  )
    allocate( particles%vi_in    ( particles%n_max           ), source = 1      )
    allocate( particles%ti_in    ( particles%n_max           ), source = 1      )
    allocate( particles%u        ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%r_origin ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%t_origin ( particles%n_max           ), source = 0._dp  )
    allocate( particles%tracers  ( particles%n_max, n_tracers), source = 0._dp  )

    particles%map%n = n_nearest_to_find
    allocate( particles%map%ip( mesh%nV, C%nz, particles%map%n), source = 0)
    allocate( particles%map%d ( mesh%nV, C%nz, particles%map%n), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_tracer_tracking_model_particles

  subroutine move_and_remove_particle( mesh, ice, particles, ip, dt)
    !< Move particle ip for time dt to a new position.
    !< Assumes velocity has alread been calculated for the current position.
    !< Calculates velocities for the new position.
    !< Removes the particle if its new position lies outside the ice sheet.

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(in   ) :: ice
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    integer,                                    intent(in   ) :: ip
    real(dp),                                   intent(in   ) :: dt

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'move_and_remove_particle'
    real(dp), dimension(2)         :: p
    real(dp)                       :: Hi_interp

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( test_ge_le( ip, 1, particles%n_max), 'ip out of bounds')
    call assert( particles%is_in_use( ip), 'particle is not in use')
#endif

    ! Using the interpolated velocities from the last time
    ! this model was run, move the particle to its new position
    particles%r( ip,:) = particles%r( ip,:) + particles%u( ip,:) * dt

    ! If the new position is outside the mesh domain, remove the particle
    if ((.not. test_ge_le( particles%r( ip,1), mesh%xmin, mesh%xmax)) .or. &
        (.not. test_ge_le( particles%r( ip,2), mesh%ymin, mesh%ymax))) then
      call remove_particle( particles, ip)
      call finalise_routine( routine_name)
      return
    end if

    ! Find the vertex and triangle containing the new position
    p = particles%r( ip,1:2)
    call find_containing_triangle( mesh, p, particles%ti_in( ip))
    call find_containing_vertex  ( mesh, p, particles%vi_in( ip))

    ! Calculate the zeta coordinate of the new position
    call calc_particle_zeta( mesh, ice%Hi, ice%Hs, &
      particles%r( ip,1), particles%r( ip,2), particles%r( ip,3), particles%ti_in( ip), &
      particles%zeta( ip), Hi_interp_ = Hi_interp)

    ! If the new position is outside the ice sheet, remove the particle
    if (particles%zeta( ip) < 0._dp .or. particles%zeta( ip) > 1._dp .or. Hi_interp < 0.1_dp) then
      call remove_particle( particles, ip)
      call finalise_routine( routine_name)
      return
    end if

    ! Interpolate the current ice velocity solution to the new position
    call interpolate_3d_velocities_to_3D_point( mesh, ice%u_3D_b, ice%v_3D_b, ice%w_3D, &
      particles%r( ip,1), particles%r( ip,2), particles%zeta( ip), &
      particles%vi_in( ip), particles%ti_in( ip), &
      particles%u( ip,1), particles%u( ip,2), particles%u( ip,3))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine move_and_remove_particle

  subroutine add_particle( mesh, ice, particles, x, y, z, time)
    !< Add a new particle to the particle-based tracer tracking model.
    !< Inserts it into the first available memory slot.
    !< If no slots are available, allocates additional memory.

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(in   ) :: ice
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    real(dp),                                   intent(in   ) :: x, y, z, time

    ! Local variables
    character(len=1024), parameter :: routine_name = 'add_particle'
    real(dp), dimension(2)         :: p
    integer                        :: ti_in, vi_in
    real(dp)                       :: zeta, Hi_interp
    real(dp)                       :: u,v,w
    integer                        :: ipp, ip
    logical                        :: out_of_memory

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( test_ge_le( x, mesh%xmin, mesh%xmax) .and. test_ge_le( y, mesh%ymin, mesh%ymax), &
      '[x,y] lies outside the mesh domain')
#endif

    ! Find where on the mesh the new particle is located
    p = [x,y]
    vi_in = 1
    call find_containing_vertex  ( mesh, p, vi_in)
    ti_in = mesh%iTri( vi_in,1)
    call find_containing_triangle( mesh, p, ti_in)

    call calc_particle_zeta( mesh, ice%Hi, ice%Hs, x, y, z, &
      ti_in, zeta, Hi_interp_ = Hi_interp)

#if (DO_ASSERTIONS)
    call assert( test_ge_le( zeta, 0._dp, 1._dp) .and. Hi_interp > 0.1_dp, &
      '[x,y,z] does not lie inside the ice sheet')
#endif

    call interpolate_3d_velocities_to_3D_point( mesh, ice%u_3D_b, ice%v_3D_b, ice%w_3D, &
      x, y, zeta, vi_in, ti_in, u, v, w)

    ! Find the first empty memory slot. If none can be found, throw an error.
    out_of_memory = .true.
    do ipp = 1, particles%n_max
      if (.not. particles%is_in_use( ipp)) then
        ! Place the new particle here
        out_of_memory = .false.
        ip = ipp
        exit
      end if
    end do
    if (out_of_memory) then
      call crash('Out of memory - exceeded maximum number of particles')
    end if

    ! Add the new particle data to the list
    particles%id_max = particles%id_max + 1_int8

    particles%is_in_use( ip  ) = .true.
    particles%id       ( ip  ) = particles%id_max
    particles%r        ( ip,:) = [x,y,z]
    particles%zeta     ( ip  ) = zeta
    particles%vi_in    ( ip  ) = vi_in
    particles%ti_in    ( ip  ) = ti_in
    particles%u        ( ip,:) = [u,v,w]
    particles%r_origin ( ip,:) = [x,y,z]
    particles%t_origin ( ip  ) = time

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_particle

  subroutine remove_particle( particles, ip)
    !< Removes particle ip from the particle-based tracer-tracking model,
    !< and marks its memory slot as available.

    ! In/output variables
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    integer,                                    intent(in   ) :: ip

    ! Local variables
    character(len=1024), parameter :: routine_name = 'remove_particle'

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( particles%is_in_use( ip), 'particle is not in use')
#endif

    particles%is_in_use( ip   ) = .false.
    particles%id       ( ip   ) = 0_int8
    particles%r        ( ip, :) = 0._dp
    particles%zeta     ( ip   ) = 0._dp
    particles%vi_in    ( ip   ) = 1
    particles%ti_in    ( ip   ) = 1
    particles%u        ( ip, :) = 0._dp
    particles%r_origin ( ip, :) = 0._dp
    particles%t_origin ( ip   ) = 0._dp
    particles%tracers  ( ip, :) = 0._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_particle

  subroutine calc_particle_zeta( mesh, Hi, Hs, x, y, z, ti_in, zeta, Hi_interp_, Hs_interp_)
    !< Calculate the zeta coordinate of a particle located at position [x,y,z]
    !< NOTE: allows the particle to be located outside the ice sheet
    !< (in which case, zeta will be < 0 or > 1).

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi, Hs
    real(dp),                               intent(in   ) :: x,y,z
    integer,                                intent(inout) :: ti_in
    real(dp),                               intent(  out) :: zeta
    real(dp), optional,                     intent(  out) :: Hi_interp_, Hs_interp_

    ! Local variables
    real(dp), dimension(2) :: p
    real(dp)               :: Hi_interp, Hs_interp

    ! Interpolate Hi, Hs horizontally to [x,y]
    p = [x,y]
    call interpolate_to_point_dp_2D( mesh, Hi, p, ti_in, Hi_interp)
    call interpolate_to_point_dp_2D( mesh, Hs, p, ti_in, Hs_interp)

    Hi_interp = max( 0.1_dp, Hi_interp)

    ! Calculate zeta
    zeta = (Hs_interp - z) / Hi_interp

    if (present( Hi_interp_)) Hi_interp_ = Hi_interp
    if (present( Hs_interp_)) Hs_interp_ = Hs_interp

  end subroutine calc_particle_zeta

  subroutine interpolate_3d_velocities_to_3D_point( mesh, u_3D_b, v_3D_b, w_3D, &
    x, y, zeta, vi_in, ti_in, u, v, w)
    !< Interpolate the ice velocity fields u,v,w to a particle located at [x,y,zeta]

    ! In- and output variables
    type(type_mesh),                             intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2,C%nz), intent(in   ) :: u_3D_b, v_3D_b
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz), intent(in   ) :: w_3D
    real(dp),                                    intent(in   ) :: x,y,zeta
    integer,                                     intent(in   ) :: vi_in, ti_in
    real(dp),                                    intent(  out) :: u,v,w

    ! Local variables
    real(dp), dimension(C%nz) :: u_col, v_col, w_col
    real(dp)                  :: zeta_limited, wwk1, wwk2
    integer                   :: k1, k2

    call interpolate_3d_velocities_to_3D_point_uv( mesh, u_3D_b, v_3D_b, &
      x, y, vi_in, u_col, v_col)

    call interpolate_3d_velocities_to_3D_point_w( mesh, w_3D, &
      x, y, ti_in, w_col)

    ! Interpolate u,v,w vertically
    zeta_limited = min( 1._dp, max( 0._dp, zeta ))
    k1 = 1
    k2 = 2
    do while (mesh%zeta( k2) < zeta_limited)
      k1 = k2
      k2 = k2 + 1
    end do

    wwk1 = (mesh%zeta( k2) - zeta) / (mesh%zeta( k2) - mesh%zeta( k1))
    wwk2 = 1._dp - wwk1

    u = wwk1 * u_col( k1) + wwk2 * u_col( k2)
    v = wwk1 * v_col( k1) + wwk2 * v_col( k2)
    w = wwk1 * w_col( k1) + wwk2 * w_col( k2)

  end subroutine interpolate_3d_velocities_to_3D_point

  subroutine interpolate_3d_velocities_to_3D_point_uv( mesh, u_3D_b, v_3D_b, &
    x, y, vi_in, u_col, v_col)
    !< Interpolate the horizontal ice velocity fields u,v to a particle located at [x,y]
    !< NOTE: does not include vertical interpolation; instead, returns
    !< velocity columns u_col( nz), v_col( nz)

    ! In- and output variables
    type(type_mesh),                             intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2,C%nz), intent(in   ) :: u_3D_b, v_3D_b
    real(dp),                                    intent(in   ) :: x,y
    integer,                                     intent(in   ) :: vi_in
    real(dp), dimension(C%nz),                   intent(  out) :: u_col, v_col

    ! Local variables
    integer                   :: iti, ti, ti_nearest, ierr
    real(dp)                  :: ww_sum, dist_min, dist, ww
    real(dp), dimension(C%nz) :: u, v, u_nearest, v_nearest, u_col_sum, v_col_sum

    ! u,v are defined on the triangles, so average over itriangles around vi_in

    ti_nearest = 0
    ww_sum     = 0._dp
    u_col_sum  = 0._dp
    v_col_sum  = 0._dp
    dist_min   = 1._dp / mesh%tol_dist

    do iti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,iti)

      if (test_ge_le( ti, mesh%ti1, mesh%ti2)) then
        u = u_3D_b( ti,:)
        v = v_3D_b( ti,:)
      else
        u = -huge( u)
        v = -huge( v)
      end if
      call MPI_ALLREDUCE( MPI_IN_PLACE, u, C%nz, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, v, C%nz, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

      dist = norm2( mesh%Tricc( ti,:) - [x,y])
      ww = 1._dp / dist**2
      ww_sum    = ww_sum    + ww
      u_col_sum = u_col_sum + ww * u
      v_col_sum = v_col_sum + ww * v
      if (dist < dist_min) then
        dist_min   = dist
        ti_nearest = ti
        u_nearest = u
        v_nearest = v
      end if
    end do
    if (dist_min < mesh%tol_dist) then
      ! p lies on a triangle circumcentre; just use u,v from that triangle
      u_col = u_nearest
      v_col = v_nearest
    else
      u_col = u_col_sum / ww_sum
      v_col = v_col_sum / ww_sum
    end if

  end subroutine interpolate_3d_velocities_to_3D_point_uv

  subroutine interpolate_3d_velocities_to_3D_point_w( mesh, w_3D, &
    x, y, ti_in, w_col)
    !< Interpolate the vertical ice velocity field w to a particle located at [x,y]
    !< NOTE: does not include vertical interpolation; instead, returns
    !< velocity column w_col( nz)

    ! In- and output variables
    type(type_mesh),                             intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz), intent(in   ) :: w_3D
    real(dp),                                    intent(in   ) :: x,y
    integer,                                     intent(in   ) :: ti_in
    real(dp), dimension(C%nz),                   intent(  out) :: w_col

    ! Local variables
    real(dp), dimension(2) :: p
    integer                :: ti_in_

    ! w is defined on the vertices
    p = [x,y]
    ti_in_ = ti_in
    call interpolate_to_point_dp_3D( mesh, w_3D, p, ti_in_, w_col)

  end subroutine interpolate_3d_velocities_to_3D_point_w

  subroutine map_tracer_to_mesh( mesh, particles, f_particles, f_mesh)
    !< Map a data field defined on the tracer particles to the model mesh

    ! In/output variables
    type(type_mesh),                             intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles),  intent(in   ) :: particles
    real(dp), dimension(particles%n_max),            intent(in   ) :: f_particles
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

  subroutine write_to_netcdf_file( particles, time)
    !< Write particle data to NetCDF

    ! In/output variables:
    type(type_tracer_tracking_model_particles), intent(in) :: particles
    real(dp),                                   intent(in) :: time

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'write_to_netcdf_file'
    character(len=1024)                          :: filename
    integer                                      :: ncid, id_dim_time, ti
    integer(int8), dimension(:,:  ), allocatable :: id_with_time
    real(dp),      dimension(:,:,:), allocatable :: r_with_time
    real(dp),      dimension(:,:  ), allocatable :: t_origin_with_time

    ! Add routine to path
    call init_routine( routine_name)

    filename = particles%nc%filename

    call open_existing_netcdf_file_for_writing( filename, ncid)

    ! Write time and find timeframe
    call write_time_to_file( filename, ncid, time)
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Add "pretend" time dimension
    if (par%master) then

      allocate( id_with_time      ( particles%n_max,    1))
      allocate( r_with_time       ( particles%n_max, 3, 1))
      allocate( t_origin_with_time( particles%n_max,    1))

      id_with_time      ( :  ,1) = particles%id
      r_with_time       ( :,:,1) = particles%r
      t_origin_with_time( :  ,1) = particles%t_origin

    end if

    ! Write data
    call write_var_master( filename, ncid, particles%nc%id_var_id       , id_with_time      , &
      start = [1,ti], count = [particles%n_max,1])
    call write_var_master( filename, ncid, particles%nc%id_var_r        , r_with_time       , &
      start = [1,1,ti], count = [particles%n_max,3,1])
    call write_var_master( filename, ncid, particles%nc%id_var_t_origin , t_origin_with_time, &
      start = [1,ti], count = [particles%n_max,1])

    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_file

  subroutine create_particles_netcdf_file( filename, particles)
    !< Create a NetCDF output file for the raw particle data

    ! In/output variables:
    character(len=*),                           intent(in   ) :: filename
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables
    character(len=1024), parameter :: routine_name = 'create_particles_netcdf_file'
    integer                        :: ncid, n, three, time

    ! Add routine to path
    call init_routine( routine_name)

    ! Create and open new file
    call create_new_netcdf_file_for_writing( filename, ncid)

    particles%nc%filename = trim(filename)

    ! Define dimensions
    call create_dimension( filename, ncid, 'n',     particles%n_max   , particles%nc%id_dim_n    )
    call create_dimension( filename, ncid, 'three', 3             , particles%nc%id_dim_three)
    call create_dimension( filename, ncid, 'time',  NF90_UNLIMITED, particles%nc%id_dim_time )

    ! Abbreviations for shorter code
    n     = particles%nc%id_dim_n
    three = particles%nc%id_dim_three
    time  = particles%nc%id_dim_time

    ! Define variables
    call create_variable( filename, ncid, 'time', &
      NF90_DOUBLE, [time], particles%nc%id_var_time)
    call create_variable( filename, ncid, 'id', &
      NF90_INT64, [n, time], particles%nc%id_var_id)
    call create_variable( filename, ncid, 'r', &
      NF90_DOUBLE, [n, three, time], particles%nc%id_var_r)
    call create_variable( filename, ncid, 't_origin', &
      NF90_DOUBLE, [n, time], particles%nc%id_var_t_origin)

    ! Close file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_particles_netcdf_file

end module tracer_tracking_model_particles