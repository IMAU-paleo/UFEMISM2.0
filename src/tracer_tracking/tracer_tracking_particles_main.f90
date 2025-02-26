module tracer_tracking_model_particles_main

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model_particles
  use model_configuration, only: C
  use tracer_tracking_model_particles_basic, only: update_particle_velocity, create_particle
  use SMB_model_types, only: type_SMB_model

  implicit none

  private

  public :: initialise_tracer_tracking_model_particles, run_tracer_tracking_model_particles

  integer,  parameter :: n_tracers         = 1
  integer,  parameter :: n_nearest_to_find = 4
  real(dp), parameter :: dt_tracer_tracking_add_new_particles = 100._dp

contains

  subroutine initialise_tracer_tracking_model_particles( mesh, ice, particles, n_max)
    !< Initialise the particle-based tracer-tracking model

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

    ! Basic data
    particles%n_max  = n_max
    particles%id_max = 0
    allocate( particles%is_in_use    ( particles%n_max           ), source = .false.)
    allocate( particles%id           ( particles%n_max           ), source = 0      )
    allocate( particles%r            ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%zeta         ( particles%n_max           ), source = 0._dp  )
    allocate( particles%vi_in        ( particles%n_max           ), source = 1      )
    allocate( particles%ti_in        ( particles%n_max           ), source = 1      )
    allocate( particles%u            ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%r_origin     ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%t_origin     ( particles%n_max           ), source = 0._dp  )
    allocate( particles%tracers      ( particles%n_max, n_tracers), source = 0._dp  )

    ! Position timeframes (to enable asynchronous time-stepping)
    allocate( particles%t0           ( particles%n_max           ), source = 0._dp  )
    allocate( particles%t1           ( particles%n_max           ), source = 0._dp  )
    allocate( particles%r_t0         ( particles%n_max, 3        ), source = 0._dp  )
    allocate( particles%r_t1         ( particles%n_max, 3        ), source = 0._dp  )

    ! Particles-to-mesh mapping
    particles%map%n = n_nearest_to_find
    allocate( particles%map%ip( mesh%nV, C%nz, particles%map%n), source = 0)
    allocate( particles%map%d ( mesh%nV, C%nz, particles%map%n), source = 0._dp)


    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_tracer_tracking_model_particles

  subroutine run_tracer_tracking_model_particles( mesh, ice, SMB, particles, time)
    !< Run the particle-based tracer-tracking model

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(in   ) :: ice
    type(type_SMB_model),                       intent(in   ) :: SMB
    type(type_tracer_tracking_model_particles), intent(inout) :: particles
    real(dp),                                   intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_tracer_tracking_model_particles'
    integer                        :: ip
    real(dp)                       :: w0, w1

    ! Add routine to path
    call init_routine( routine_name)


    ! Move and remove particles
    do ip = 1, particles%n_max

      if (.not. particles%is_in_use( ip)) cycle

      ! If needed, update particle velocity and timeframes
      do while (particles%t1( ip) < time)

        particles%r( ip,:) = particles%r_t1( ip,:)
        call update_particle_velocity( mesh, ice, particles, ip)

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

  end subroutine run_tracer_tracking_model_particles

end module tracer_tracking_model_particles_main