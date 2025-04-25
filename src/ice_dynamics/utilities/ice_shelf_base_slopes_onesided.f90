module ice_shelf_base_slopes_onesided

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D
  use mpi_distributed_memory, only: gather_to_all
  use ice_geometry_basics, only: ice_surface_elevation

  implicit none

  private

  public :: calc_ice_shelf_base_slopes_onesided

contains

  subroutine calc_ice_shelf_base_slopes_onesided( mesh, ice)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_ice_shelf_base_slopes_onesided'
    logical,  dimension(mesh%nV)           :: mask_floating_ice_tot, mask_cf_fl_tot
    integer                                :: ti, via, vib, vic
    logical                                :: all_are_shelf, none_are_cf

    ! Add routine to path
    call init_routine( routine_name)

    ! Straightforward calculation everywhere
    call ddx_a_b_2D( mesh, ice%Hib, ice%dHib_dx_b)
    call ddy_a_b_2D( mesh, ice%Hib, ice%dHib_dy_b)

    ! Gather global mask
    call gather_to_all( ice%mask_floating_ice, mask_floating_ice_tot)
    call gather_to_all( ice%mask_cf_fl       , mask_cf_fl_tot       )

    ! Set slopes to zero on triangles that are only partly shelf
    do ti = mesh%ti1, mesh%ti2

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      all_are_shelf = mask_floating_ice_tot( via) .and. mask_floating_ice_tot( vib) .and. mask_floating_ice_tot( vic)
      none_are_cf = (.not. mask_cf_fl_tot( via)) .and. (.not. mask_cf_fl_tot( vib)) .and. (.not. mask_cf_fl_tot( vic))

      if (.not. (all_are_shelf .and. none_are_cf)) then
        ice%dHib_dx_b( ti) = 0._dp
        ice%dHib_dy_b( ti) = 0._dp
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_shelf_base_slopes_onesided

end module ice_shelf_base_slopes_onesided
