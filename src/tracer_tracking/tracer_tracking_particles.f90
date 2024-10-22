module tracer_tracking_model_particles

  ! The main tracer tracking model module.

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles
  use model_configuration, only: C
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex

  implicit none

  private

  public :: initialise_tracer_tracking_model_particles, interpolate_particles_velocities, &
    calc_particle_zeta

  integer, parameter :: n_particles_init = 100
  integer, parameter :: n_tracers        = 16

contains

  subroutine initialise_tracer_tracking_model_particles( mesh, particles)

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_tracer_tracking_model_particles'

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%master)  write(*,'(a)') '     Initialising particle-based tracer tracking model...'

    allocate( particles%is_in_use( n_particles_init           ), source = .false.)
    allocate( particles%r        ( n_particles_init, 3        ), source = 0._dp  )
    allocate( particles%vi_in    ( n_particles_init           ), source = 1      )
    allocate( particles%ti_in    ( n_particles_init           ), source = 1      )
    allocate( particles%u        ( n_particles_init, 3        ), source = 0._dp  )
    allocate( particles%r_origin ( n_particles_init, 3        ), source = 0._dp  )
    allocate( particles%t_origin ( n_particles_init           ), source = 0._dp  )
    allocate( particles%tracers  ( n_particles_init, n_tracers), source = 0._dp  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_tracer_tracking_model_particles

  subroutine update_particles_vi_ti_in( mesh, particles)

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables
    character(len=1024), parameter :: routine_name = 'update_particles_vi_ti_in'
    integer                        :: ip
    real(dp), dimension(2)         :: p

    ! Add routine to path
    call init_routine( routine_name)

    do ip = 1, size( particles%r, 1)
      if (particles%is_in_use( ip)) then

        p = [particles%r( ip,1), particles%r( ip,2)]
        call find_containing_triangle( mesh, p, particles%ti_in( ip))
        call find_containing_vertex  ( mesh, p, particles%vi_in( ip))

      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_particles_vi_ti_in

  subroutine interpolate_particles_velocities( mesh, ice, particles)

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(in   ) :: ice
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'interpolate_particles_velocities'
    integer                        :: ip
    real(dp)                       :: zeta

    ! Add routine to path
    call init_routine( routine_name)

    do ip = 1, size( particles%r, 1)
      if (particles%is_in_use( ip)) then

        call calc_particle_zeta( mesh, ice%Hi, ice%Hs, &
          particles%r( ip,1), particles%r( ip,2), particles%r( ip,3), &
          particles%ti_in( ip), zeta)

        call interpolate_3d_velocities_to_3D_point( mesh, ice%u_3D_b, ice%v_3D_b, ice%w_3D, &
          particles%r( ip,1), particles%r( ip,2), zeta, &
          particles%vi_in( ip), particles%ti_in( ip), &
          particles%u( ip,1), particles%u( ip,2), particles%u( ip,3))

      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine interpolate_particles_velocities

  subroutine calc_particle_zeta( mesh, Hi, Hs, x, y, z, ti_in, zeta)

    ! In- and output variables
    type(type_mesh),        intent(in   ) :: mesh
    real(dp), dimension(:), intent(in   ) :: Hi, Hs
    real(dp),               intent(in   ) :: x,y,z
    integer,                intent(in   ) :: ti_in
    real(dp),               intent(  out) :: zeta

    ! Local variables
    integer  :: n, vi, vi_nearest
    real(dp) :: ww_sum, dist_min, dist, ww, Hi_sum, Hs_sum, Hi_interp, Hs_interp

    ! Interpolate Hib,Hs horizontally to [x,y]
    vi_nearest = 0
    ww_sum     = 0._dp
    Hi_sum     = 0._dp
    Hs_sum     = 0._dp
    dist_min   = 1._dp / mesh%tol_dist

    do n = 1, 3
      vi = mesh%Tri( ti_in,n)
      dist = norm2( mesh%V( vi,:) - [x,y])
      ww = 1._dp / dist**2
      ww_sum = ww_sum + ww
      Hi_sum = Hi_sum + ww * Hi( vi)
      Hs_sum = Hs_sum + ww * Hs( vi)
      if (dist < dist_min) then
        dist_min   = dist
        vi_nearest = vi
      end if
    end do
    if (dist_min < mesh%tol_dist) then
      ! p lies on a vertex; just use Hib,Hs from that triangle
      Hi_interp = Hi( vi_nearest)
      Hs_interp = Hs( vi_nearest)
    else
      Hi_interp = Hi_sum / ww_sum
      Hs_interp = Hs_sum / ww_sum
    end if

    Hi_interp = max( 0.1_dp, Hi_interp)

#if (DO_ASSERTIONS)
    ! Safety
    call assert( test_ge_le( z, Hs_interp - Hi_interp, Hs_interp), 'particle lies outside of the ice')
#endif

    ! Calculate zeta
    zeta = (Hs_interp - z) / Hi_interp

  end subroutine calc_particle_zeta

  subroutine interpolate_3d_velocities_to_3D_point( mesh, u_3D_b, v_3D_b, w_3D, &
    x, y, zeta, vi_in, ti_in, u, v, w)

    ! In- and output variables
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(in   ) :: u_3D_b, v_3D_b, w_3D
    real(dp),                 intent(in   ) :: x,y,zeta
    integer,                  intent(in   ) :: vi_in, ti_in
    real(dp),                 intent(  out) :: u,v,w

    ! Local variables
    real(dp), dimension(size(u_3D_b,2)) :: u_col, v_col, w_col
    real(dp)                            :: zeta_limited, wwk1, wwk2
    integer                             :: k1, k2

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

    ! In- and output variables
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(in   ) :: u_3D_b, v_3D_b
    real(dp),                 intent(in   ) :: x,y
    integer,                  intent(in   ) :: vi_in
    real(dp), dimension(:),   intent(  out) :: u_col, v_col

    ! Local variables
    integer                             :: iti, ti, ti_nearest
    real(dp)                            :: ww_sum, dist_min, dist, ww
    real(dp), dimension(size(u_3D_b,2)) :: u_col_sum, v_col_sum

    ! u,v are defined on the triangles, so average over itriangles around vi_in

    ti_nearest = 0
    ww_sum     = 0._dp
    u_col_sum  = 0._dp
    v_col_sum  = 0._dp
    dist_min   = 1._dp / mesh%tol_dist

    do iti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,iti)
      dist = norm2( mesh%Tricc( ti,:) - [x,y])
      ww = 1._dp / dist**2
      ww_sum    = ww_sum    + ww
      u_col_sum = u_col_sum + ww * u_3D_b( ti,:)
      v_col_sum = v_col_sum + ww * v_3D_b( ti,:)
      if (dist < dist_min) then
        dist_min   = dist
        ti_nearest = ti
      end if
    end do
    if (dist_min < mesh%tol_dist) then
      ! p lies on a triangle circumcentre; just use u,v from that triangle
      u_col = u_3D_b( ti_nearest,:)
      v_col = v_3D_b( ti_nearest,:)
    else
      u_col = u_col_sum / ww_sum
      v_col = v_col_sum / ww_sum
    end if

  end subroutine interpolate_3d_velocities_to_3D_point_uv

  subroutine interpolate_3d_velocities_to_3D_point_w( mesh, w_3D, &
    x, y, ti_in, w_col)

    ! In- and output variables
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(in   ) :: w_3D
    real(dp),                 intent(in   ) :: x,y
    integer,                  intent(in   ) :: ti_in
    real(dp), dimension(:),   intent(  out) :: w_col

    ! Local variables
    integer                           :: n, vi, vi_nearest
    real(dp)                          :: ww_sum, dist_min, dist, ww
    real(dp), dimension(size(w_3D,2)) :: w_col_sum

    ! w is defined on the vertices, so average over the vertices of ti_in

    vi_nearest = 0
    ww_sum     = 0._dp
    w_col_sum  = 0._dp
    dist_min   = 1._dp / mesh%tol_dist

    do n = 1, 3
      vi = mesh%Tri( ti_in,n)
      dist = norm2( mesh%V( vi,:) - [x,y])
      ww = 1._dp / dist**2
      ww_sum    = ww_sum    + ww
      w_col_sum = w_col_sum + ww * w_3D( vi,:)
      if (dist < dist_min) then
        dist_min   = dist
        vi_nearest = vi
      end if
    end do
    if (dist_min < mesh%tol_dist) then
      ! p lies on a vertex; just use Hib,Hs from that triangle
      w_col = w_3D( vi_nearest,:)
    else
      w_col = w_col_sum / ww_sum
    end if

  end subroutine interpolate_3d_velocities_to_3D_point_w

  subroutine move_particle( x_t, y_t, z_t, u_t, v_t, w_t, dt, x_tplusdt, y_tplusdt, z_tplusdt)

    ! In/output variables:
    real(dp), intent(in ) :: x_t, y_t, z_t                   !< [m]    Particle position at time t
    real(dp), intent(in ) :: u_t, v_t, w_t                   !< [m/yr] Particle velocity at time t
    real(dp), intent(in ) :: dt                              !< [yr]   Time step
    real(dp), intent(out) :: x_tplusdt, y_tplusdt, z_tplusdt !< [m] Particle position at time t + dt

    x_tplusdt = x_t + dt * u_t
    y_tplusdt = y_t + dt * v_t
    z_tplusdt = z_t + dt * w_t

  end subroutine move_particle

end module tracer_tracking_model_particles