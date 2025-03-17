module tracer_tracking_model_particles_basic

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex, &
    interpolate_to_point_dp_2D_singlecore, interpolate_to_point_dp_3D_singlecore
  use reallocate_mod, only: reallocate
  use netcdf, only: NF90_UNLIMITED, NF90_INT64, NF90_DOUBLE
  use netcdf_io_main

  implicit none

  private

  public :: create_particle_at_ice_surface, destroy_particle, update_particle_velocity

contains

  subroutine update_particle_velocity( mesh, Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
    particles, ip)
    !< Update the velocity and position timeframes of particle ip,
    !< and remove the particle if it exits the ice sheet or the region domain

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    real(dp), dimension(mesh%nV),               intent(in   ) :: Hi_tot, Hs_tot
    real(dp), dimension(mesh%nTri,C%nz),        intent(in   ) :: u_3D_b_tot, v_3D_b_tot
    real(dp), dimension(mesh%nV  ,C%nz),        intent(in   ) :: w_3D_tot
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    integer,                                    intent(in   ) :: ip

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_particle_velocity'
    real(dp), dimension(2)         :: p
    real(dp)                       :: Hi_interp
    real(dp)                       :: dt

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( test_ge_le( ip, 1, particles%n_max), 'ip out of bounds')
    call assert( particles%is_in_use( ip), 'particle is not in use')
#endif

    ! Locate the particle inside the ice sheet
    p = particles%r( ip,1:2)
    call find_containing_triangle( mesh, p, particles%ti_in( ip))
    call find_containing_vertex  ( mesh, p, particles%vi_in( ip))
    call calc_particle_zeta( mesh, Hi_tot, Hs_tot, &
      particles%r( ip,1), particles%r( ip,2), particles%r( ip,3), particles%ti_in( ip), &
      particles%zeta( ip), Hi_interp_ = Hi_interp)

    ! If the particle now lies outside the ice sheet, remove it
    if (particles%zeta( ip) < 0._dp .or. particles%zeta( ip) > 1._dp .or. Hi_interp < 0.1_dp) then
      call destroy_particle( particles, ip)
      call finalise_routine( routine_name)
      return
    end if

    ! Interpolate the current ice velocity solution to the new position
    call interpolate_3d_velocities_to_3D_point( mesh, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
      particles%r( ip,1), particles%r( ip,2), particles%zeta( ip), &
      particles%vi_in( ip), particles%ti_in( ip), &
      particles%u( ip,1), particles%u( ip,2), particles%u( ip,3))

    ! Calculate particle time step
    dt = calc_particle_dt( particles%u( ip,1), particles%u( ip,2), particles%u( ip,3))

    ! Update timeframes
    particles%t0  ( ip  ) = particles%t1( ip)
    particles%t1  ( ip  ) = particles%t1( ip) + dt
    particles%r_t0( ip,:) = particles%r_t1( ip,:)
    particles%r_t1( ip,:) = particles%r_t1( ip,:) + dt * particles%u( ip,:)

    ! If the particle's new position takes it ouside the mesh domain, remove it
    if (particles%r_t1( ip,1) < mesh%xmin + C%tractrackpart_dx_particle .or. &
        particles%r_t1( ip,1) > mesh%xmax - C%tractrackpart_dx_particle .or. &
        particles%r_t1( ip,2) < mesh%ymin + C%tractrackpart_dx_particle .or. &
        particles%r_t1( ip,2) > mesh%ymax - C%tractrackpart_dx_particle) then
      call destroy_particle( particles, ip)
      call finalise_routine( routine_name)
      return
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_particle_velocity

  subroutine create_particle_at_ice_surface( mesh, Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
    particles, x, y, time)
    !< Create a new particle at position [x,y,zeta=0]

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    real(dp), dimension(mesh%nV),               intent(in   ) :: Hi_tot, Hs_tot
    real(dp), dimension(mesh%nTri,C%nz),        intent(in   ) :: u_3D_b_tot, v_3D_b_tot
    real(dp), dimension(mesh%nV  ,C%nz),        intent(in   ) :: w_3D_tot
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    real(dp),                                   intent(in   ) :: x, y, time

    ! Local variables
    character(len=1024), parameter :: routine_name = 'create_particle_at_ice_surface'
    integer                        :: ipp, ip
    logical                        :: out_of_memory
    real(dp), dimension(2)         :: p
    integer                        :: ti_in, vi_in
    real(dp)                       :: zeta, z, Hi_interp
    real(dp)                       :: u,v,w,dt

    ! Add routine to path
    call init_routine( routine_name)

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

    ! Locate the particle inside the ice sheet
#if (DO_ASSERTIONS)
    call assert( test_ge_le( x, mesh%xmin, mesh%xmax) .and. test_ge_le( y, mesh%ymin, mesh%ymax), &
      '[x,y] lies outside the mesh domain')
#endif

    p = [x,y]
    vi_in = 1
    call find_containing_vertex  ( mesh, p, vi_in)
    ti_in = mesh%iTri( vi_in,1)
    call find_containing_triangle( mesh, p, ti_in)

    zeta = 0._dp
    call calc_particle_z( mesh, Hi_tot, Hs_tot, x, y, zeta, ti_in, z, Hi_interp_ = Hi_interp)

    ! Interpolate ice velocity to particle position
    call interpolate_3d_velocities_to_3D_point( mesh, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
      x, y, zeta, vi_in, ti_in, u, v, w)

    ! Calculate particle time step
    dt = calc_particle_dt( u, v, w)

    ! Add the new particle data to the list
    particles%id_max = particles%id_max + par%n

    particles%is_in_use( ip  ) = .true.
    particles%id       ( ip  ) = particles%id_max
    particles%r        ( ip,:) = [x,y,z]
    particles%zeta     ( ip  ) = zeta
    particles%vi_in    ( ip  ) = vi_in
    particles%ti_in    ( ip  ) = ti_in
    particles%u        ( ip,:) = [u,v,w]
    particles%r_origin ( ip,:) = [x,y,z]
    particles%t_origin ( ip  ) = time

    particles%t0       ( ip  ) = time
    particles%t1       ( ip  ) = time + dt
    particles%r_t0     ( ip,:) = [x,y,z]
    particles%r_t1     ( ip,:) = [x,y,z] + dt * [u,v,w]

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_particle_at_ice_surface

  subroutine destroy_particle( particles, ip)
    !< Destroy particle ip

    ! In/output variables
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    integer,                                    intent(in   ) :: ip

    ! Local variables
    character(len=1024), parameter :: routine_name = 'destroy_particle'

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    call assert( particles%is_in_use( ip), 'particle is not in use')
#endif

    particles%is_in_use( ip   ) = .false.
    particles%id       ( ip   ) = 0
    particles%r        ( ip, :) = 0._dp
    particles%zeta     ( ip   ) = 0._dp
    particles%vi_in    ( ip   ) = 1
    particles%ti_in    ( ip   ) = 1
    particles%u        ( ip, :) = 0._dp
    particles%r_origin ( ip, :) = 0._dp
    particles%t_origin ( ip   ) = 0._dp
    ! particles%tracers  ( ip, :) = 0._dp

    particles%t0       ( ip  ) = 0._dp
    particles%t1       ( ip  ) = 0._dp
    particles%r_t0     ( ip,:) = 0._dp
    particles%r_t1     ( ip,:) = 0._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine destroy_particle

  function calc_particle_dt( u, v, w) result( dt)
    !< Calculate particle time step

    ! In- and output variables
    real(dp), intent(in   ) :: u,v,w
    real(dp)                :: dt

    ! Local variables
    real(dp)                :: uabs

    uabs = sqrt( u**2 + v**2 + w**2)
    dt = C%tractrackpart_dx_particle / uabs
    dt = min( C%tractrackpart_dt_particle_max, max( C%tractrackpart_dt_particle_min, dt))

  end function calc_particle_dt

  subroutine calc_particle_zeta( mesh, Hi_tot, Hs_tot, x, y, z, ti_in, zeta, Hi_interp_, Hs_interp_)
    !< Calculate the zeta coordinate of a particle located at position [x,y,z]
    !< NOTE: allows the particle to be located outside the ice sheet
    !< (in which case, zeta will be < 0 or > 1).

    ! In- and output variables
    type(type_mesh),              intent(in   ) :: mesh
    real(dp), dimension(mesh%nV), intent(in   ) :: Hi_tot, Hs_tot
    real(dp),                     intent(in   ) :: x,y,z
    integer,                      intent(inout) :: ti_in
    real(dp),                     intent(  out) :: zeta
    real(dp), optional,           intent(  out) :: Hi_interp_, Hs_interp_

    ! Local variables
    real(dp), dimension(2) :: p
    real(dp)               :: Hi_interp, Hs_interp

    ! Interpolate Hi, Hs horizontally to [x,y]
    p = [x,y]
    call interpolate_to_point_dp_2D_singlecore( mesh, Hi_tot, p, ti_in, Hi_interp)
    call interpolate_to_point_dp_2D_singlecore( mesh, Hs_tot, p, ti_in, Hs_interp)

    Hi_interp = max( 0.1_dp, Hi_interp)

    ! Calculate zeta
    zeta = (Hs_interp - z) / Hi_interp

    if (present( Hi_interp_)) Hi_interp_ = Hi_interp
    if (present( Hs_interp_)) Hs_interp_ = Hs_interp

  end subroutine calc_particle_zeta

  subroutine calc_particle_z( mesh, Hi_tot, Hs_tot, x, y, zeta, ti_in, z, Hi_interp_, Hs_interp_)
    !< Calculate the z coordinate of a particle located at position [x,y,zeta]
    !< NOTE: allows the particle to be located outside the ice sheet
    !< (in which case, zeta will be < 0 or > 1).

    ! In- and output variables
    type(type_mesh),              intent(in   ) :: mesh
    real(dp), dimension(mesh%nV), intent(in   ) :: Hi_tot, Hs_tot
    real(dp),                     intent(in   ) :: x,y,zeta
    integer,                      intent(inout) :: ti_in
    real(dp),                     intent(  out) :: z
    real(dp), optional,           intent(  out) :: Hi_interp_, Hs_interp_

    ! Local variables
    real(dp), dimension(2) :: p
    real(dp)               :: Hi_interp, Hs_interp

    ! Interpolate Hi, Hs horizontally to [x,y]
    p = [x,y]
    call interpolate_to_point_dp_2D_singlecore( mesh, Hi_tot, p, ti_in, Hi_interp)
    call interpolate_to_point_dp_2D_singlecore( mesh, Hs_tot, p, ti_in, Hs_interp)

    Hi_interp = max( 0.1_dp, Hi_interp)

    ! Calculate z
    z = Hs_interp - zeta * Hi_interp

    if (present( Hi_interp_)) Hi_interp_ = Hi_interp
    if (present( Hs_interp_)) Hs_interp_ = Hs_interp

  end subroutine calc_particle_z

  subroutine interpolate_3d_velocities_to_3D_point( mesh, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
    x, y, zeta, vi_in, ti_in, u, v, w)
    !< Interpolate the ice velocity fields u,v,w to a particle located at [x,y,zeta]

    ! In- and output variables
    type(type_mesh),                     intent(in   ) :: mesh
    real(dp), dimension(mesh%nTri,C%nz), intent(in   ) :: u_3D_b_tot, v_3D_b_tot
    real(dp), dimension(mesh%nV  ,C%nz), intent(in   ) :: w_3D_tot
    real(dp),                            intent(in   ) :: x,y,zeta
    integer,                             intent(in   ) :: vi_in, ti_in
    real(dp),                            intent(  out) :: u,v,w

    ! Local variables
    real(dp), dimension(C%nz) :: u_col, v_col, w_col
    real(dp)                  :: zeta_limited, wwk1, wwk2
    integer                   :: k1, k2

    call interpolate_3d_velocities_to_3D_point_uv( mesh, u_3D_b_tot, v_3D_b_tot, &
      x, y, vi_in, u_col, v_col)

    call interpolate_3d_velocities_to_3D_point_w( mesh, w_3D_tot, &
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

  subroutine interpolate_3d_velocities_to_3D_point_uv( mesh, u_3D_b_tot, v_3D_b_tot, &
    x, y, vi_in, u_col, v_col)
    !< Interpolate the horizontal ice velocity fields u,v to a particle located at [x,y]
    !< NOTE: does not include vertical interpolation; instead, returns
    !< velocity columns u_col( nz), v_col( nz)

    ! In- and output variables
    type(type_mesh),                     intent(in   ) :: mesh
    real(dp), dimension(mesh%ntri,C%nz), intent(in   ) :: u_3D_b_tot, v_3D_b_tot
    real(dp),                            intent(in   ) :: x,y
    integer,                             intent(in   ) :: vi_in
    real(dp), dimension(C%nz),           intent(  out) :: u_col, v_col

    ! Local variables
    integer                   :: iti, ti
    real(dp)                  :: ww_sum, dist_min, dist, ww
    real(dp), dimension(C%nz) :: u, v, u_nearest, v_nearest, u_col_sum, v_col_sum

    ! u,v are defined on the triangles, so average over itriangles around vi_in

    ww_sum     = 0._dp
    u_col_sum  = 0._dp
    v_col_sum  = 0._dp
    dist_min   = 1._dp / mesh%tol_dist

    do iti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,iti)

      u = u_3D_b_tot( ti,:)
      v = v_3D_b_tot( ti,:)

      dist = norm2( mesh%Tricc( ti,:) - [x,y])
      ww = 1._dp / dist**2
      ww_sum    = ww_sum    + ww
      u_col_sum = u_col_sum + ww * u
      v_col_sum = v_col_sum + ww * v

      if (dist < dist_min) then
        dist_min  = dist
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

  subroutine interpolate_3d_velocities_to_3D_point_w( mesh, w_3D_tot, &
    x, y, ti_in, w_col)
    !< Interpolate the vertical ice velocity field w to a particle located at [x,y]
    !< NOTE: does not include vertical interpolation; instead, returns
    !< velocity column w_col( nz)

    ! In- and output variables
    type(type_mesh),                   intent(in   ) :: mesh
    real(dp), dimension(mesh%nV,C%nz), intent(in   ) :: w_3D_tot
    real(dp),                          intent(in   ) :: x,y
    integer,                           intent(in   ) :: ti_in
    real(dp), dimension(C%nz),         intent(  out) :: w_col

    ! Local variables
    real(dp), dimension(2) :: p
    integer                :: ti_in_

    ! w is defined on the vertices
    p = [x,y]
    ti_in_ = ti_in
    call interpolate_to_point_dp_3D_singlecore( mesh, w_3D_tot, p, ti_in_, w_col)

  end subroutine interpolate_3d_velocities_to_3D_point_w

end module tracer_tracking_model_particles_basic