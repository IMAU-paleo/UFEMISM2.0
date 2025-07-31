module tracer_tracking_model_particles_main

  use mpi_basic, only: par, sync
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles
  use model_configuration, only: C
  use tracer_tracking_model_particles_basic, only: update_particle_velocity, create_particle_at_ice_surface
  use SMB_model_types, only: type_SMB_model
  use grid_basic, only: setup_square_grid
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_all
  use tracer_tracking_model_particles_io, only: create_particles_netcdf_file, write_to_particles_netcdf_file
  use tracer_tracking_model_particles_remapping, only: calc_particles_to_mesh_map, map_tracer_to_mesh
  use mesh_utilities, only: find_containing_vertex, find_containing_triangle

  implicit none

  private

  public :: initialise_tracer_tracking_model_particles, run_tracer_tracking_model_particles, &
    remap_tracer_tracking_model_particles

  ! integer,  parameter :: n_tracers         = 1

contains

  subroutine initialise_tracer_tracking_model_particles( mesh, ice, particles)
    !< Initialise the particle-based tracer-tracking model

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(in   ) :: ice
    type(type_tracer_tracking_model_particles), intent(  out) :: particles

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_tracer_tracking_model_particles'
    character(len=256)             :: grid_name
    character(len=1024)            :: filename

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  write(*,'(a)') '     Initialising particle-based tracer-tracking model...'

    ! Basic data
    particles%n_max  = C%tractrackpart_n_max_particles
    particles%id_max = par%i
    allocate( particles%is_in_use    ( particles%n_max           ), source = .false.)
    allocate( particles%id           ( particles%n_max           ), source = 0      )
    allocate( particles%r            ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%zeta         ( particles%n_max           ), source = 0._dp  )
    allocate( particles%vi_in        ( particles%n_max           ), source = 1      )
    allocate( particles%ti_in        ( particles%n_max           ), source = 1      )
    allocate( particles%u            ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%r_origin     ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%t_origin     ( particles%n_max           ), source = 0._dp  )
    ! allocate( particles%tracers      ( particles%n_max, n_tracers), source = 0._dp  )

    ! Position timeframes (to enable asynchronous time-stepping)
    allocate( particles%t0           ( particles%n_max           ), source = 0._dp  )
    allocate( particles%t1           ( particles%n_max           ), source = 0._dp  )
    allocate( particles%r_t0         ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%r_t1         ( particles%n_max, 3        ), source = 0._dp  )

    ! Particles-to-mesh mapping
    particles%map%n = C%tractrackpart_remap_n_nearest
    allocate( particles%map%ip( mesh%nV, C%nz, particles%map%n), source = 0)
    allocate( particles%map%d ( mesh%nV, C%nz, particles%map%n), source = 0._dp)

    ! Grid for creating new particles
    particles%t_add_new_particles = C%start_time_of_run
    grid_name = 'particle_creation_grid'
    call setup_square_grid( grid_name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, &
      C%tractrackpart_dx_new_particles, particles%grid_new_particles, &
      mesh%lambda_M, mesh%phi_M, mesh%beta_stereo)

    ! Raw particle data output file
    if (C%tractrackpart_write_raw_output) then
      particles%t_write_raw_output = C%start_time_of_run
      filename = trim( C%output_dir) // '/tracer_tracking_particles.nc'
      call create_particles_netcdf_file( filename, particles)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_tracer_tracking_model_particles

  subroutine run_tracer_tracking_model_particles( mesh, ice, SMB, particles, time, age)
    !< Run the particle-based tracer-tracking model

    ! In- and output variables
    type(type_mesh),                             intent(in   ) :: mesh
    type(type_ice_model),                        intent(in   ) :: ice
    type(type_SMB_model),                        intent(in   ) :: SMB
    type(type_tracer_tracking_model_particles),  intent(inout) :: particles
    real(dp),                                    intent(in   ) :: time
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz), intent(  out) :: age

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'run_tracer_tracking_model_particles'
    real(dp), dimension(mesh%nV)          :: Hi_tot, Hs_tot
    real(dp), dimension(mesh%nTri,C%nz)   :: u_3D_b_tot, v_3D_b_tot
    real(dp), dimension(mesh%nV  ,C%nz)   :: w_3D_tot
    real(dp), dimension(:,:), allocatable :: Hi_grid_tot
    real(dp), dimension(:,:), allocatable :: SMB_grid_tot
    real(dp), dimension(particles%n_max)  :: age_p

    ! Add routine to path
    call init_routine( routine_name)

    allocate( Hi_grid_tot         ( particles%grid_new_particles%nx, particles%grid_new_particles%ny))
    allocate( SMB_grid_tot        ( particles%grid_new_particles%nx, particles%grid_new_particles%ny))

    call gather_ice_model_data( mesh, ice, SMB, particles, Hi_tot, Hs_tot, &
      u_3D_b_tot, v_3D_b_tot, w_3D_tot, Hi_grid_tot, SMB_grid_tot)

    ! If the time is right, add a new batch of particles
    if (time >= particles%t_add_new_particles) then
      particles%t_add_new_particles = particles%t_add_new_particles + C%tractrackpart_dt_new_particles
      call add_new_particles_from_SMB( mesh, &
        Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
        Hi_grid_tot, SMB_grid_tot, particles, time)
    end if

    ! Move and remove particles
    call move_and_remove_particles( mesh, particles, time, &
      Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot)

    ! Map tracers to the model mesh
    call calc_particles_to_mesh_map( mesh, particles)
    age_p = time - particles%t_origin
    call map_tracer_to_mesh( mesh, particles, age_p, age)

    ! If the time is right, write raw particle data to output file
    if (C%tractrackpart_write_raw_output) then
      if (time >= particles%t_write_raw_output) then
        particles%t_write_raw_output = particles%t_write_raw_output + C%tractrackpart_dt_raw_output
        call write_to_particles_netcdf_file( particles, time)
      end if
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_tracer_tracking_model_particles

  subroutine gather_ice_model_data( mesh, ice, SMB, particles, Hi_tot, Hs_tot, &
      u_3D_b_tot, v_3D_b_tot, w_3D_tot, Hi_grid_tot, SMB_grid_tot)
    !< Gather distributed ice model data for interpolating to the particles

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(in   ) :: ice
    type(type_SMB_model),                       intent(in   ) :: SMB
    type(type_tracer_tracking_model_particles), intent(in   ) :: particles
    real(dp), dimension(mesh%nV),               intent(  out) :: Hi_tot, Hs_tot
    real(dp), dimension(mesh%nTri,C%nz),        intent(  out) :: u_3D_b_tot, v_3D_b_tot
    real(dp), dimension(mesh%nV  ,C%nz),        intent(  out) :: w_3D_tot
    real(dp), dimension(:,:),                   intent(  out) :: Hi_grid_tot
    real(dp), dimension(:,:),                   intent(  out) :: SMB_grid_tot

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'gather_ice_model_data'
    real(dp), dimension(:  ), allocatable :: Hi_grid_vec_partial
    real(dp), dimension(:  ), allocatable :: SMB_grid_vec_partial

    ! Add routine to path
    call init_routine( routine_name)

    allocate( Hi_grid_vec_partial ( particles%grid_new_particles%n1: particles%grid_new_particles%n2))
    allocate( SMB_grid_vec_partial( particles%grid_new_particles%n1: particles%grid_new_particles%n2))

    ! Map ice thickness and SMB to the particle creation grid
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, particles%grid_new_particles, C%output_dir, &
      ice%Hi , Hi_grid_vec_partial)
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, particles%grid_new_particles, C%output_dir, &
      SMB%SMB, SMB_grid_vec_partial)

    ! Gather data to all processes, so they can be interpolated to the particle positions
    ! (necessary, as a particle owned by process n will generally not be located in the
    ! domain of that process)
    call gather_to_all( ice%Hi    , Hi_tot)
    call gather_to_all( ice%Hs    , Hs_tot)
    call gather_to_all( ice%u_3D_b, u_3D_b_tot)
    call gather_to_all( ice%v_3D_b, v_3D_b_tot)
    call gather_to_all( ice%w_3D  , w_3D_tot)
    call gather_gridded_data_to_all( particles%grid_new_particles, Hi_grid_vec_partial,  Hi_grid_tot)
    call gather_gridded_data_to_all( particles%grid_new_particles, SMB_grid_vec_partial, SMB_grid_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_ice_model_data

  subroutine add_new_particles_from_SMB( mesh, &
    Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
    Hi_grid_tot, SMB_grid_tot, particles, time)
    !< Add a batch of regularly-spaced particles in the accumulation zone

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    real(dp), dimension(mesh%nV),               intent(in   ) :: Hi_tot, Hs_tot
    real(dp), dimension(mesh%nTri,C%nz),        intent(in   ) :: u_3D_b_tot, v_3D_b_tot
    real(dp), dimension(mesh%nV  ,C%nz),        intent(in   ) :: w_3D_tot
    real(dp), dimension(:,:),                   intent(in   ) :: Hi_grid_tot
    real(dp), dimension(:,:),                   intent(in   ) :: SMB_grid_tot
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    real(dp),                                   intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_new_particles_from_SMB'
    integer                        :: n,i,j,n_added
    real(dp)                       :: x,y

    ! Add routine to path
    call init_routine( routine_name)

    ! For each grid point with a positive ice thickness and SMB, add a new particle
    n_added = 0
    do n = par%i+1, particles%grid_new_particles%n, par%n

      i = particles%grid_new_particles%n2ij( n,1)
      j = particles%grid_new_particles%n2ij( n,2)

      x = particles%grid_new_particles%x( i)
      y = particles%grid_new_particles%y( j)

      if (Hi_grid_tot( i,j) > 10._dp .and. SMB_grid_tot( i,j) > 0._dp .and. &
          x > mesh%xmin + C%tractrackpart_dx_particle .and. &
          x < mesh%xmax - C%tractrackpart_dx_particle .and. &
          y > mesh%ymin + C%tractrackpart_dx_particle .and. &
          y < mesh%ymax - C%tractrackpart_dx_particle) then
        n_added = n_added + 1
        call create_particle_at_ice_surface( mesh, Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
          particles, x, y, time)
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_new_particles_from_SMB

  subroutine move_and_remove_particles( mesh, particles, time, &
    Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot)

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    real(dp),                                   intent(in   ) :: time
    real(dp), dimension(mesh%nV),               intent(in   ) :: Hi_tot, Hs_tot
    real(dp), dimension(mesh%nTri,C%nz),        intent(in   ) :: u_3D_b_tot, v_3D_b_tot
    real(dp), dimension(mesh%nV  ,C%nz),        intent(in   ) :: w_3D_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'move_and_remove_particles'
    integer                        :: ip
    real(dp)                       :: w0, w1

    ! Add routine to path
    call init_routine( routine_name)

    do ip = 1, particles%n_max

      if (.not. particles%is_in_use( ip)) cycle

      ! If needed, update particle velocity and timeframes
      do while (particles%t1( ip) < time)

        particles%r( ip,:) = particles%r_t1( ip,:)
        call update_particle_velocity( mesh, Hi_tot, Hs_tot, u_3D_b_tot, v_3D_b_tot, w_3D_tot, &
          particles, ip)

        ! If the particle has been removed, stop tracing it
        if (.not. particles%is_in_use( ip)) exit

      end do

      ! If the particle has been removed, stop tracing it
      if (.not. particles%is_in_use( ip)) cycle

      ! Interpolate particle position between timeframes
      w0 = (particles%t1( ip) - time) / (particles%t1( ip) - particles%t0( ip))
      w1 = 1._dp - w0
      particles%r( ip,:) = (w0 * particles%r_t0( ip,:)) + (w1 * particles%r_t1( ip,:))

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine move_and_remove_particles

  subroutine remap_tracer_tracking_model_particles( mesh_old, mesh_new, particles, time, age)
    !< Run the particle-based tracer-tracking model

    ! In- and output variables
    type(type_mesh),                                     intent(in   ) :: mesh_old, mesh_new
    type(type_tracer_tracking_model_particles),          intent(inout) :: particles
    real(dp),                                            intent(in   ) :: time
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2,C%nz), intent(  out) :: age

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'remap_tracer_tracking_model_particles'
    integer                               :: ip, vi, ti
    real(dp), dimension(2)                :: p
    real(dp), dimension(particles%n_max)  :: age_p

    ! Add routine to path
    call init_routine( routine_name)

    particles%vi_in = 1
    particles%ti_in = 1
    do ip = 1, particles%n_max
      if (.not. particles%is_in_use( ip)) cycle
      p = particles%r( ip,1:2)
      call find_containing_vertex  ( mesh_new, p, particles%vi_in( ip))
      particles%ti_in( ip) = mesh_new%iTri( particles%vi_in( ip),1)
      call find_containing_triangle( mesh_new, p, particles%ti_in( ip))
    end do

    ! Map tracers from the particles to the new model mesh
    call calc_particles_to_mesh_map( mesh_new, particles)
    age_p = time - particles%t_origin
    call map_tracer_to_mesh( mesh_new, particles, age_p, age)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_tracer_tracking_model_particles

end module tracer_tracking_model_particles_main