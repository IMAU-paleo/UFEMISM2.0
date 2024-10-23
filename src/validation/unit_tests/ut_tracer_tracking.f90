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
  use tracer_tracking_model_particles, only: calc_particle_zeta, interpolate_3d_velocities_to_3D_point

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

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_tracer_tracking_main

  subroutine test_calc_particle_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_calc_particle_zeta'
    character(len=1024), parameter :: test_name_local = 'calc_particle_zeta'
    character(len=1024)            :: test_name
    real(dp), dimension(mesh%nV)   :: Hi,Hb,Hs
    integer                        :: vi
    integer                        :: n_particles,i,j,k
    real(dp)                       :: zmin,zmax
    real(dp)                       :: x,y,z
    real(dp), dimension(2)         :: p
    integer                        :: ti_in
    real(dp)                       :: Hi_ex, Hb_ex, Hs_ex, zeta_ex
    real(dp)                       :: zeta
    logical                        :: verified

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified = .true.

    ! Set up a simple test geometry
    do vi = 1, mesh%nV
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%V( vi,1), mesh%V( vi,2), &
        Hi = Hi( vi), Hb = Hb( vi), Hs = Hs( vi))
    end do

    zmin = minval( Hb)
    zmax = maxval( Hs)

    ! Determine zeta for a grid of test particles
    n_particles = 100
    ti_in = 1
    do i = 2, n_particles-1
    do j = 2, n_particles-1
    do k = 1, n_particles

      x = mesh%xmin + (mesh%xmax - mesh%xmin) * real( i-1,dp) / real( n_particles-1,dp)
      y = mesh%ymin + (mesh%ymax - mesh%ymin) * real( j-1,dp) / real( n_particles-1,dp)
      z =      zmin + (     zmax -      zmin) * real( k-1,dp) / real( n_particles-1,dp)

      p = [x,y]
      call find_containing_triangle( mesh, p, ti_in)

      ! Exact solution
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, x, y, &
        Hi = Hi_ex, Hb = Hb_ex, Hs = Hs_ex, &
        z = z, zeta = zeta_ex)

      ! Skip particles that lie outside of the ice
      if (zeta_ex < 0.05_dp .or. zeta_ex > 0.95_dp) cycle

      ! Interpolation
      call calc_particle_zeta( mesh, Hi, Hs, x, y, z, ti_in, zeta)

      ! Check result
      verified = verified .and. test_tol( zeta, zeta_ex, 0.05_dp)

    end do
    end do
    end do

    call unit_test( verified, trim( test_name) // '/zeta')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_calc_particle_zeta

  subroutine test_interpolate_3d_velocities_to_3D_point( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_interpolate_3d_velocities_to_3D_point'
    character(len=1024), parameter      :: test_name_local = 'interpolate_3d_velocities_to_3D_point'
    character(len=1024)                 :: test_name
    real(dp), dimension(mesh%nTri,C%nz) :: u_3D_b, v_3D_b
    real(dp), dimension(mesh%nV  ,C%nz) :: w_3D
    integer                             :: vi,ti,k
    integer                             :: n_particles,i,j
    real(dp)                            :: x,y,zeta
    real(dp), dimension(2)              :: p
    integer                             :: vi_in, ti_in
    real(dp)                            :: u_ex, v_ex, w_ex
    real(dp)                            :: u,v,w
    logical                             :: verified_u, verified_v, verified_w

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_u = .true.
    verified_v = .true.
    verified_w = .true.

    ! Set up a simple test velocity field
    do ti = 1, mesh%nTri
    do k = 1, C%nz
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%Tricc( ti,1), mesh%Tricc( ti,2), &
        zeta_q = mesh%zeta( k), u = u_3D_b( ti,k), v = v_3D_b( ti,k))
    end do
    end do

    do vi = 1, mesh%nV
    do k = 1, C%nz
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%V( vi,1), mesh%V( vi,2), &
        zeta_q = mesh%zeta( k), w = w_3D( vi,k))
    end do
    end do

    ! Determine velocities for a grid of test particles
    n_particles = 100
    vi_in = 1
    ti_in = 1
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

      ! Check result
      verified_u = verified_u .and. test_tol( u, u_ex, 10._dp)
      verified_v = verified_v .and. test_tol( v, v_ex, 10._dp)
      verified_w = verified_w .and. test_tol( w, w_ex, 0.1_dp)

    end do
    end do
    end do

    call unit_test( verified_u, trim( test_name) // '/u')
    call unit_test( verified_v, trim( test_name) // '/v')
    call unit_test( verified_w, trim( test_name) // '/w')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_interpolate_3d_velocities_to_3D_point

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
      Hi = 2000._dp + 500._dp * (cos( 3._dp * pi * cx * (x - xmin)) + cos( 2._dp * pi * cy * (y - ymin)))
      Hb =            200._dp * (sin( 2._dp * pi * cx * (x - xmin)) + sin( 3._dp * pi * cy * (y - ymin)))
      Hs = Hb + Hi
    end if

    if (present( z)) then
      call assert( present(zeta), 'need both z and zeta')
      zeta = (Hs - z) / Hi
    end if

    if (present( u)) then
      call assert( present( zeta_q), 'need both zeta_q and u')
      u = (1._dp - zeta_q) * (100._dp + 50._dp * (sin( 2.5_dp * pi * cx * (x - xmin)) + cos( 3.5_dp * pi * cy * (y - ymin))))
    end if

    if (present( v)) then
      call assert( present( zeta_q), 'need both zeta_q and v')
      v = (1._dp - zeta_q) * (100._dp + 50._dp * (cos( 3.5_dp * pi * cx * (x - xmin)) + sin( 2.5_dp * pi * cy * (y - ymin))))
    end if

    if (present( w)) then
      call assert( present( zeta_q), 'need both zeta_q and w')
      w = (1._dp - zeta_q) * (          -1._dp * (sin( 1.5_dp * pi * cx * (x - xmin)) + sin( 2.2_dp * pi * cy * (y - ymin))))
    end if

  end subroutine calc_simple_test_geometry

end module ut_tracer_tracking