module BMB_inverted

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use BMB_model_types, only: type_BMB_model_inverted
  use reference_geometry_types, only: type_reference_geometry

  implicit none

  private

  public :: initialise_BMB_model_inverted, run_BMB_model_inverted

contains

  subroutine run_BMB_model_inverted( mesh, ice, BMB_inv, time)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_ice_model),          intent(in   ) :: ice
    type(type_BMB_model_inverted), intent(inout) :: BMB_inv
    real(dp),                      intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_BMB_model_inverted'

    ! Add routine to path
    call init_routine( routine_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_BMB_model_inverted

  subroutine initialise_BMB_model_inverted( mesh, BMB_inv, refgeo_PD, refgeo_init)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_BMB_model_inverted), intent(inout) :: BMB_inv
    type(type_reference_geometry), intent(in   ) :: refgeo_PD, refgeo_init

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_BMB_model_inverted'

    ! Add routine to path
    call init_routine( routine_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BMB_model_inverted

end module BMB_inverted
