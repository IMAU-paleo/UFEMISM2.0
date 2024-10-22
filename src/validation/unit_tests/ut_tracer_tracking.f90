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
  use mesh_utilities, only: find_containing_triangle
  use tracer_tracking_model_particles, only: calc_particle_zeta

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
    call test_calc_particle_zeta( test_name, mesh)

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
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, &
        mesh%V( vi,1), mesh%V( vi,2), Hi( vi), Hb( vi), Hs( vi))
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
      call calc_simple_test_geometry( mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, &
        x, y, Hi_ex, Hb_ex, Hs_ex, z, zeta_ex)

      ! Skip particles that lie outside of the ice
      if (zeta_ex < 0.05_dp .or. zeta_ex > 0.95_dp) cycle

      ! Interpolation
      call calc_particle_zeta( mesh, Hi, Hs, x, y, z, ti_in, zeta)

      ! Check result
      verified = verified .and. test_tol( zeta, zeta_ex, 0.05_dp)

    end do
    end do
    end do

    call unit_test( verified, test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_calc_particle_zeta

  subroutine calc_simple_test_geometry( xmin, xmax, ymin, ymax, x, y, Hi, Hb, Hs, z, zeta)
    ! In/output variables
    real(dp), intent(in   )           :: xmin, xmax, ymin, ymax, x, y
    real(dp), intent(  out)           :: Hi, Hb, Hs
    real(dp), intent(in   ), optional :: z
    real(dp), intent(  out), optional :: zeta
    ! Local variables
    real(dp) :: xp, yp

    xp = (x - xmin) / (xmax - xmin)
    yp = (y - ymin) / (ymax - ymin)

    Hb =            200._dp * (sin( 2._dp * pi * xp) + sin( 3._dp * pi * yp))
    Hi = 2000._dp + 500._dp * (cos( 3._dp * pi * xp) + cos( 2._dp * pi * yp))
    Hs = Hb + Hi

    if (present( z)) then
      zeta = (Hs - z) / Hi
    end if

  end subroutine calc_simple_test_geometry

end module ut_tracer_tracking