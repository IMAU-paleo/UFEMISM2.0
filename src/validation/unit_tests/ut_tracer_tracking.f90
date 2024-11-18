module ut_tracer_tracking

  ! Unit tests for the tracer tracking module

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par
  use model_configuration, only: C
  use parameters, only: pi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles, type_map_particles_to_mesh
  use tracer_tracking_model_particles, only: calc_particle_zeta, interpolate_3d_velocities_to_3D_point, &
    calc_particles_to_mesh_map
  use mpi
  use mpi_distributed_memory, only: gather_to_all_dp_1D

  implicit none

  private

  public :: unit_tests_tracer_tracking_main

contains

  subroutine unit_tests_tracer_tracking_main( test_name_parent)
    ! Test the tracer tracking subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_tracer_tracking_main'
    character(len=1024), parameter :: test_name_local = 'tracer_tracking'
    character(len=1024)            :: test_name
    real(dp)                       :: xmin, xmax, ymin, ymax, alpha_min, res_max
    character(len=1024)            :: name
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    name = 'test_mesh'
    xmin = -500e3_dp
    xmax =  500e3_dp
    ymin = -500e3_dp
    ymax =  500e3_dp
    alpha_min = 25._dp * pi / 180._dp
    res_max = 50e3_dp

    call allocate_mesh_primary( mesh, name, 100, 200, C%nC_mem)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call calc_all_secondary_mesh_data( mesh, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT)

    ! Run unit tests on this test mesh
    call test_calc_particle_zeta                   ( test_name, mesh)
    call test_interpolate_3d_velocities_to_3D_point( test_name, mesh)
    call test_calc_particles_to_mesh_map           ( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_tracer_tracking_main

  subroutine test_calc_particle_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_calc_particle_zeta'
    character(len=1024), parameter         :: test_name_local = 'calc_particle_zeta'
    character(len=1024)                    :: test_name
    real(dp), dimension(mesh%vi1:mesh%vi2) :: Hi,Hb,Hs
    integer                                :: vi
    integer                                :: nx,ny,nz,i,j,k
    real(dp)                               :: zmin,zmax
    integer                                :: ierr
    real(dp)                               :: x,y,z
    real(dp), dimension(2)                 :: p
    integer                                :: ti_in
    real(dp)                               :: Hi_ex, Hb_ex, Hs_ex, zeta_ex
    real(dp)                               :: Hi_interp, Hs_interp, zeta_interp, max_err_zeta
    logical                                :: verified

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Set up a simple test geometry
    do vi = mesh%vi1, mesh%vi2
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%V( vi,1), mesh%V( vi,2), &
        Hi = Hi( vi), Hb = Hb( vi), Hs = Hs( vi))
    end do

    zmin = minval( Hb)
    zmax = maxval( Hs)
    call MPI_ALLREDUCE( MPI_IN_PLACE, zmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, zmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Determine zeta for a grid of test particles
    nx = 100
    ny = 100
    nz = 100
    ti_in = 1
    max_err_zeta = 0._dp
    do i = 2, nx-1
    do j = 2, ny-1
    do k = 1, nz

      x = mesh%xmin + (mesh%xmax - mesh%xmin) * real( i-1,dp) / real( nx-1,dp)
      y = mesh%ymin + (mesh%ymax - mesh%ymin) * real( j-1,dp) / real( ny-1,dp)
      z =      zmin + (     zmax -      zmin) * real( k-1,dp) / real( nz-1,dp)

      p = [x,y]
      call find_containing_triangle( mesh, p, ti_in)

      ! Exact solution
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, x, y, &
        Hi = Hi_ex, Hb = Hb_ex, Hs = Hs_ex, &
        z = z, zeta = zeta_ex)

      ! Skip particles that lie outside of the ice
      if (zeta_ex < 0.05_dp .or. zeta_ex > 0.95_dp) cycle

      ! Interpolation
      call calc_particle_zeta( mesh, Hi, Hs, x, y, z, ti_in, zeta_interp, &
        Hi_interp_ = Hi_interp, Hs_interp_ = Hs_interp)
      max_err_zeta = max( max_err_zeta, abs( zeta_interp - zeta_ex))

    end do
    end do
    end do

    verified = max_err_zeta <= 0.005_dp

    call unit_test( verified, trim( test_name) // '/zeta')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_calc_particle_zeta

  subroutine test_interpolate_3d_velocities_to_3D_point( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_interpolate_3d_velocities_to_3D_point'
    character(len=1024), parameter              :: test_name_local = 'interpolate_3d_velocities_to_3D_point'
    character(len=1024)                         :: test_name
    real(dp), dimension(mesh%ti1:mesh%ti2,C%nz) :: u_3D_b, v_3D_b
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz) :: w_3D
    integer                                     :: vi,ti,k
    integer                                     :: n_particles,i,j
    real(dp)                                    :: x,y,zeta
    real(dp), dimension(2)                      :: p
    integer                                     :: vi_in, ti_in
    real(dp)                                    :: u_ex, v_ex, w_ex
    real(dp)                                    :: u, v, w, max_diff_u, max_diff_v, max_diff_w
    logical                                     :: verified_u, verified_v, verified_w

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_u = .true.
    verified_v = .true.
    verified_w = .true.

    ! Set up a simple test velocity field
    do ti = mesh%ti1, mesh%ti2
    do k = 1, C%nz
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%Tricc( ti,1), mesh%Tricc( ti,2), &
        zeta_q = mesh%zeta( k), u = u_3D_b( ti,k), v = v_3D_b( ti,k))
    end do
    end do

    do vi = mesh%vi1, mesh%vi2
    do k = 1, C%nz
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%V( vi,1), mesh%V( vi,2), &
        zeta_q = mesh%zeta( k), w = w_3D( vi,k))
    end do
    end do

    ! Determine velocities for a grid of test particles
    n_particles = 100
    vi_in = 1
    ti_in = 1

    max_diff_u = 0._dp
    max_diff_v = 0._dp
    max_diff_w = 0._dp

    do i = 2, n_particles-1
    do j = 2, n_particles-1
    do k = 1, C%nz

      x    = mesh%xmin + (mesh%xmax - mesh%xmin) * real( i-1,dp) / real( n_particles-1,dp)
      y    = mesh%ymin + (mesh%ymax - mesh%ymin) * real( j-1,dp) / real( n_particles-1,dp)
      zeta = mesh%zeta( k)

      p = [x,y]
      call find_containing_vertex  ( mesh, p, vi_in)
      call find_containing_triangle( mesh, p, ti_in)

      ! Exact solution
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, x, y, &
        zeta_q = zeta, u = u_ex, v = v_ex, w = w_ex)

      ! Interpolation
      call interpolate_3d_velocities_to_3D_point( mesh, u_3D_b, v_3D_b, w_3D, &
        x, y, zeta, vi_in, ti_in, u, v, w)

      max_diff_u = max( max_diff_u, abs( u - u_ex))
      max_diff_v = max( max_diff_v, abs( v - v_ex))
      max_diff_w = max( max_diff_w, abs( w - w_ex))

    end do
    end do
    end do

    verified_u = max_diff_u <= 5._dp
    verified_v = max_diff_v <= 5._dp
    verified_w = max_diff_w <= 0.001_dp

    call unit_test( verified_u, trim( test_name) // '/u')
    call unit_test( verified_v, trim( test_name) // '/v')
    call unit_test( verified_w, trim( test_name) // '/w')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_interpolate_3d_velocities_to_3D_point

  subroutine test_calc_particles_to_mesh_map( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_calc_particles_to_mesh_map'
    character(len=1024), parameter             :: test_name_local = 'calc_particles_to_mesh_map'
    character(len=1024)                        :: test_name
    real(dp), dimension(mesh%nV)               :: Hi, Hb, Hs
    integer                                    :: nx,ny,nz,vi
    type(type_tracer_tracking_model_particles) :: particles
    integer                                    :: i,j,k,ip
    real(dp), dimension(2)                     :: p
    logical                                    :: verified
    real(dp), dimension(mesh%nV,C%nz,3)        :: rs_mesh
    real(dp), dimension(:,:), allocatable      :: rs_particles
    type(type_map_particles_to_mesh)           :: map_bruteforce
    real(dp)                                   :: dist_max, dist
    integer                                    :: ii, jj

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Initialise a simple ice geometry
    do vi = 1, mesh%nV
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%V( vi,1), mesh%V( vi,2), &
        Hi = Hi( vi), Hb = Hb( vi), Hs = Hs( vi))
    end do

    ! Initialise a set of particles
    nx = 25
    ny = 25
    nz = 16
    particles%n = nx * ny * nz
    allocate( particles%is_in_use( particles%n   ), source = .true.)
    allocate( particles%r        ( particles%n, 3), source = 0._dp)
    allocate( particles%zeta     ( particles%n   ), source = 0._dp)
    allocate( particles%vi_in    ( particles%n   ), source = 1    )
    allocate( particles%ti_in    ( particles%n   ), source = 1    )

    ip = 0
    do i = 1, nx
    do j = 1, ny
    do k = 1, nz
      ip = ip + 1
      particles%r( ip,1) = mesh%xmin + (mesh%xmax - mesh%xmin) * real( i-1,dp) / real( nx-1,dp)
      particles%r( ip,2) = mesh%ymin + (mesh%ymax - mesh%ymin) * real( j-1,dp) / real( nx-1,dp)
      particles%zeta( ip) = real( k-1,dp) / real(nz-1,dp)
      p = [particles%r( ip,1), particles%r( ip,2)]
      call find_containing_triangle( mesh, p, particles%ti_in( ip))
      call find_containing_vertex  ( mesh, p, particles%vi_in( ip))
    end do
    end do
    end do

    particles%map%n = 4
    allocate( particles%map%ip( mesh%nV, C%nz, particles%map%n), source = 0)
    allocate( particles%map%d ( mesh%nV, C%nz, particles%map%n), source = 0._dp)

    call calc_particles_to_mesh_map( mesh, particles)

    ! Compare to a brute-force search
    ! ===============================

    ! Calculate scaled coordinates for mesh vertex-layers and for particles
    do vi = 1, mesh%nV
      do k = 1, C%nz
        rs_mesh( vi,k,1) = (mesh%V( vi,1) - mesh%xmin) / (mesh%xmax - mesh%xmin)
        rs_mesh( vi,k,2) = (mesh%V( vi,2) - mesh%ymin) / (mesh%ymax - mesh%ymin)
        rs_mesh( vi,k,3) = mesh%zeta( k)
      end do
    end do

    allocate(rs_particles( particles%n, 3))
    do ip = 1, particles%n
      rs_particles( ip,1) = (particles%r( ip,1) - mesh%xmin) / (mesh%xmax - mesh%xmin)
      rs_particles( ip,2) = (particles%r( ip,2) - mesh%ymin) / (mesh%ymax - mesh%ymin)
      rs_particles( ip,3) = particles%zeta( ip)
    end do

    ! Allocate memory
    dist_max = norm2( [mesh%xmin,mesh%ymin] - [mesh%xmax,mesh%ymax])
    map_bruteforce%n = particles%map%n
    allocate( map_bruteforce%ip( mesh%nV, C%nz, map_bruteforce%n), source = 0)
    allocate( map_bruteforce%d ( mesh%nV, C%nz, map_bruteforce%n), source = dist_max)

    do vi = 1, mesh%nV
      do k = 1, C%nz
        do ip = 1, particles%n
          ! Calculate the scaled distance between particle ip and vertex-layer [vi,k]
          dist = norm2( rs_particles( ip,:) - rs_mesh( vi,k,:))
          do ii = 1, map_bruteforce%n
            if (dist < map_bruteforce%d( vi,k,ii)) then
              do jj = map_bruteforce%n, ii+1, -1
                map_bruteforce%d(  vi,k,jj) = map_bruteforce%d(  vi,k,jj-1)
                map_bruteforce%ip( vi,k,jj) = map_bruteforce%ip( vi,k,jj-1)
              end do
              map_bruteforce%d(  vi,k,ii) = dist
              map_bruteforce%ip( vi,k,ii) = ip
              exit
            end if
          end do
        end do
      end do
    end do

    ! Compare results
    verified = .true.
    do vi = 1, mesh%nV
      do k = 1, C%nz
        do ii = 1, particles%map%n
          verified = verified .and. any( particles%map%ip( vi,k,ii) == map_bruteforce%ip( vi,k,:))
        end do
      end do
    end do

    call unit_test( verified, trim( test_name) // '/map')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_calc_particles_to_mesh_map

  subroutine calc_simple_test_geometry( xmin, xmax, ymin, ymax, x, y, &
    Hi, Hb, Hs, z, zeta, zeta_q, u, v, w)

    use parameters, only: ice_density, grav
    use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap

    ! In/output variables
    real(dp), intent(in   )           :: xmin, xmax, ymin, ymax, x, y
    real(dp), intent(  out), optional :: Hi, Hb, Hs
    real(dp), intent(in   ), optional :: z
    real(dp), intent(  out), optional :: zeta
    real(dp), intent(in   ), optional :: zeta_q
    real(dp), intent(  out), optional :: u, v, w

    ! Local variables
    real(dp) :: cx, cy

    cx = 1._dp / (xmax - xmin)
    cy = 1._dp / (ymax - ymin)

    if (present( Hi)) then
      call assert( present( Hs) .and. present( Hs), 'need both Hi, Hb and Hs')
      Hi = 2000._dp + 500._dp * (cos( 1._dp * pi * cx * (x - xmin)) + cos( 2._dp * pi * cy * (y - ymin)))
      Hb =            200._dp * (sin( 2._dp * pi * cx * (x - xmin)) + sin( 1._dp * pi * cy * (y - ymin)))
      Hs = Hb + Hi
    end if

    if (present( z)) then
      call assert( present(zeta), 'need both z and zeta')
      zeta = (Hs - z) / Hi
    end if

    if (present( u)) then
      call assert( present( zeta_q), 'need both zeta_q and u')
      u = (1._dp - zeta_q) * (100._dp + 50._dp * (sin( 1.5_dp * pi * cx * (x - xmin)) + cos( 0.5_dp * pi * cy * (y - ymin))))
    end if

    if (present( v)) then
      call assert( present( zeta_q), 'need both zeta_q and v')
      v = (1._dp - zeta_q) * (100._dp + 50._dp * (cos( 0.5_dp * pi * cx * (x - xmin)) + sin( 1.5_dp * pi * cy * (y - ymin))))
    end if

    if (present( w)) then
      call assert( present( zeta_q), 'need both zeta_q and w')
      w = (1._dp - zeta_q) * (          -1._dp * (sin( 0.5_dp * pi * cx * (x - xmin)) + sin( 0.2_dp * pi * cy * (y - ymin))))
    end if

  end subroutine calc_simple_test_geometry

end module ut_tracer_tracking