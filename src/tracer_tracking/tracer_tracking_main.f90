module tracer_tracking_model_main

  ! The main tracer tracking model module.

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, warning, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use SMB_model_types, only: type_SMB_model
  use tracer_tracking_model_types, only: type_tracer_tracking_model
  use tracer_tracking_model_particles_main, only: initialise_tracer_tracking_model_particles, &
    run_tracer_tracking_model_particles, remap_tracer_tracking_model_particles
  use reallocate_mod

  implicit none

  private

  public :: initialise_tracer_tracking_model, run_tracer_tracking_model, remap_tracer_tracking_model

  ! integer, parameter :: n_tracers        = 16

contains

  subroutine run_tracer_tracking_model( mesh, ice, SMB, tracer_tracking, time)

    ! In- and output variables
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_ice_model),             intent(in   ) :: ice
    type(type_SMB_model),             intent(in   ) :: SMB
    type(type_tracer_tracking_model), intent(inout) :: tracer_tracking
    real(dp),                         intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_tracer_tracking_model'

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_tracer_tracking_model == 'none') then
      call finalise_routine( routine_name)
      return
    end if

    if (time == tracer_tracking%t_next) then

      tracer_tracking%t_prev = tracer_tracking%t_next
      tracer_tracking%t_next = tracer_tracking%t_prev + C%tractrackpart_dt_coupling

      select case (C%choice_tracer_tracking_model)
      case default
        call crash('unknown choice_tracer_tracking_model "' // trim(C%choice_tracer_tracking_model) // '"')
      case ('particles')
        call run_tracer_tracking_model_particles( mesh, ice, SMB, &
          tracer_tracking%particles, time, tracer_tracking%age)
      end select

    elseif (time > tracer_tracking%t_next) THEN
      ! This should not be possible
      call crash('overshot the tracer tracking model coupling time step')
    else
      ! No need to do anything
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_tracer_tracking_model

  subroutine initialise_tracer_tracking_model( mesh, ice, tracer_tracking)

    ! In- and output variables
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_ice_model),             intent(in   ) :: ice
    type(type_tracer_tracking_model), intent(inout) :: tracer_tracking

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_tracer_tracking_model'

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_tracer_tracking_model == 'none') then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary)  write(*,'(a)') '   Initialising tracer-tracking model ' // &
      colour_string( trim(C%choice_tracer_tracking_model), 'light blue') // '...'

    ! Allocate memory for the model-independent tracer-tracking data
    allocate( tracer_tracking%age( mesh%vi1:mesh%vi2, C%nz))
    ! allocate( tracer_tracking%tracers( mesh%vi1:mesh%vi2, C%nz, n_tracers))

    ! Initialise coupling times
    tracer_tracking%t_prev = -huge( tracer_tracking%t_prev)
    tracer_tracking%t_next = C%start_time_of_run

    ! Initialise the chosen tracer-tracking model
    select case (C%choice_tracer_tracking_model)
    case default
      call crash('unknown choice_tracer_tracking_model "' // trim(C%choice_tracer_tracking_model) // '"')
    case ('particles')
      call initialise_tracer_tracking_model_particles( mesh, ice, tracer_tracking%particles)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_tracer_tracking_model

  subroutine remap_tracer_tracking_model( mesh_old, mesh_new, tracer_tracking, time)

    ! In- and output variables
    type(type_mesh),                  intent(in   ) :: mesh_old, mesh_new
    type(type_tracer_tracking_model), intent(inout) :: tracer_tracking
    real(dp),                         intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_tracer_tracking_model'

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_tracer_tracking_model == 'none') then
      call finalise_routine( routine_name)
      return
    end if

    call reallocate_bounds( tracer_tracking%age, mesh_new%vi1, mesh_new%vi2, mesh_new%nz)

    select case (C%choice_tracer_tracking_model)
    case default
      call crash('unknown choice_tracer_tracking_model "' // trim(C%choice_tracer_tracking_model) // '"')
    case ('particles')
      call remap_tracer_tracking_model_particles( mesh_old, mesh_new, &
        tracer_tracking%particles, time, tracer_tracking%age)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_tracer_tracking_model

end module tracer_tracking_model_main