module tracer_tracking_model_main

  ! The main tracer tracking model module.

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use tracer_tracking_model_types, only: type_tracer_tracking_model, &
    type_tracer_tracking_model_particles
  use model_configuration, only: C

  implicit none

  private

  public :: initialise_tracer_tracking_model

  integer, parameter :: n_particles_init = 100
  integer, parameter :: n_tracers        = 16

contains

  subroutine initialise_tracer_tracking_model( mesh, tracer_tracking)

    ! In- and output variables
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_tracer_tracking_model), intent(inout) :: tracer_tracking

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_tracer_tracking_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%master)  write(*,'(a)') '   Initialising tracer tracking model...'

    allocate( tracer_tracking%age    ( mesh%nV, C%nz           ))
    allocate( tracer_tracking%tracers( mesh%nV, C%nz, n_tracers))

    call initialise_tracer_tracking_model_particles( mesh, tracer_tracking%particles)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_tracer_tracking_model

  subroutine initialise_tracer_tracking_model_particles( mesh, particles)

    ! In- and output variables
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_tracer_tracking_model_particles), intent(inout) :: particles

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_tracer_tracking_model_particles'

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%master)  write(*,'(a)') '     Initialising particle-based tracer tracking model...'

    allocate( particles%is_in_use( n_particles_init           ))
    allocate( particles%r        ( n_particles_init, 3        ))
    allocate( particles%r_origin ( n_particles_init, 3        ))
    allocate( particles%t_origin ( n_particles_init           ))
    allocate( particles%tracers  ( n_particles_init, n_tracers))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_tracer_tracking_model_particles

end module tracer_tracking_model_main