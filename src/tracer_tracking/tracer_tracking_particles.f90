module tracer_tracking_model_particles

  ! The main tracer tracking model module.

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles
  use model_configuration, only: C
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex

  implicit none

  private

  public :: initialise_tracer_tracking_model_particles, interpolate_particles_velocities, &
    calc_particle_zeta, interpolate_3d_velocities_to_3D_point, update_particles_vi_ti_in, &
    create_particles_to_mesh_map

  integer, parameter :: n_particles_init  = 100
  integer, parameter :: n_tracers         = 16

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

    particles%n = n_particles_init
    allocate( particles%is_in_use( particles%n           ), source = .false.)
    allocate( particles%r        ( particles%n, 3        ), source = 0._dp  )
    allocate( particles%zeta     ( particles%n           ), source = 0._dp  )
    allocate( particles%vi_in    ( particles%n           ), source = 1      )
    allocate( particles%ti_in    ( particles%n           ), source = 1      )
    allocate( particles%u        ( particles%n, 3        ), source = 0._dp  )
    allocate( particles%r_origin ( particles%n, 3        ), source = 0._dp  )
    allocate( particles%t_origin ( particles%n           ), source = 0._dp  )
    allocate( particles%tracers  ( particles%n, n_tracers), source = 0._dp  )

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

    do ip = 1, particles%n
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

    do ip = 1, particles%n
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

  subroutine map_particle_tracer_to_mesh( mesh, particles, f_particles, f_interp)

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(in   ) :: particles
    real(dp), dimension(particles%n),           intent(in   ) :: f_particles
    real(dp), dimension(mesh%nV,C%nz),          intent(  out) :: f_interp

    ! Local variables
    character(len=1024), parameter :: routine_name = 'map_particle_tracer_to_mesh'
    integer                        :: vi, k, ii, ip
    real(dp)                       :: w_sum, f_sum, w, dist
    logical                        :: coincides_with_particle

    ! Add routine to path
    call init_routine( routine_name)

    f_interp = 0._dp

    do vi = 1, mesh%nV

      do k = 1, C%nz

        w_sum = 0._dp
        f_sum = 0._dp
        coincides_with_particle = .false.

        do ii = 1, particles%map%n
          ip   = particles%map%ip( vi,k,ii)
          dist = particles%map%d(  vi,k,ii)
          if (dist < mesh%tol_dist) then
            coincides_with_particle = .true.
            f_interp( vi,k) = f_particles( ip)
            exit
          end if
          w = 1._dp / dist**2
          w_sum = w_sum + w
          f_sum = f_sum + w * f_particles( ip)
        end do
        if (.not. coincides_with_particle) then
          f_interp( vi,k) = f_sum / w_sum
        end if

      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_particle_tracer_to_mesh

  subroutine create_particles_to_mesh_map( mesh, n_nearest_to_find, particles)

    ! In/output variables
    type(type_mesh),                            intent(in   ) :: mesh
    integer,                                    intent(in   ) :: n_nearest_to_find
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables
    character(len=1024), parameter      :: routine_name = 'create_particles_to_mesh_map'
    real(dp), dimension(mesh%nV,C%nz,3) :: rs_mesh
    real(dp), dimension(particles%n ,3) :: rs_particles
    integer                             :: vi,k,ip
    real(dp)                            :: dist_max
    integer                             :: vi_in
    integer,  dimension(mesh%nV)        :: map, stack
    integer                             :: stackN
    integer                             :: ii,jj,ci,vj
    integer                             :: n_cycles
    logical                             :: is_one_of_n_nearest_particles
    real(dp)                            :: dist

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

    do ip = 1, particles%n
      if (.not. particles%is_in_use( ip)) cycle
      rs_particles( ip,1) = (particles%r( ip,1) - mesh%xmin) / (mesh%xmax - mesh%xmin)
      rs_particles( ip,2) = (particles%r( ip,2) - mesh%ymin) / (mesh%ymax - mesh%ymin)
      rs_particles( ip,3) = particles%zeta( ip)
    end do

    ! Allocate memory
    dist_max = norm2( [mesh%xmin,mesh%ymin] - [mesh%xmax,mesh%ymax])
    particles%map%n = n_nearest_to_find
    allocate( particles%map%ip( mesh%nV, C%nz, particles%map%n), source = 0)
    allocate( particles%map%d ( mesh%nV, C%nz, particles%map%n), source = dist_max)

    ! Do a flood-fill style expansion around each particle

    ! Initialise map and stack
    map    = 0
    stack  = 0
    stackN = 0

    do ip = 1, particles%n

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

  end subroutine create_particles_to_mesh_map

end module tracer_tracking_model_particles