module BMB_inverted

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use mpi_basic, only: par
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use BMB_model_types, only: type_BMB_model_inverted
  use reference_geometry_types, only: type_reference_geometry
  use reference_geometries_main, only: reallocate_reference_geometry_on_mesh
  use ice_geometry_basics, only: is_floating
  use mpi_distributed_memory, only: gather_to_all

  implicit none

  private

  public :: initialise_BMB_model_inverted, run_BMB_model_inverted

contains

  subroutine run_BMB_model_inverted( mesh, ice, BMB_inv, time)
    !< Nudge basal melt rate based on mismatch between modelled and target ice geometry

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_ice_model),          intent(in   ) :: ice
    type(type_BMB_model_inverted), intent(inout) :: BMB_inv
    real(dp),                      intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_BMB_model_inverted'
    integer                        :: vi, ci, vj
    logical,  dimension(mesh%nV)   :: mask_floating_ice_tot, mask_cf_fl_tot
    real(dp), dimension(mesh%nV)   :: Hi_target_tot
    real(dp)                       :: w_sum, wH_sum, deltaH, dHdt, dBMBdt
    real(dp), parameter            :: c_H     = -0.003_dp
    real(dp), parameter            :: c_dHdt  = -0.03_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception: target ice thickness at the floating calving front
    ! is often wrong (because of the difficulty of remapping a discontinuous
    ! field), so instead use the mean of the neighbouring non-front shelf
    ! vertices.
    call gather_to_all( ice%mask_floating_ice     , mask_floating_ice_tot)
    call gather_to_all( ice%mask_cf_fl            , mask_cf_fl_tot)
    call gather_to_all( BMB_inv%target_geometry%Hi, Hi_target_tot)

    do vi = mesh%vi1, mesh% vi2
      if (mask_cf_fl_tot( vi)) then
        w_sum  = 0._dp
        wH_sum = 0._dp
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mask_floating_ice_tot( vj) .and. .not. mask_cf_fl_tot( vj)) then
            w_sum = w_sum + 1._dp
            wH_sum = wH_sum + Hi_target_tot( vj)
          end if
        end do
        if (w_sum > 0._dp) then
          Hi_target_tot( vi) = wH_sum / w_sum
        end if
      end if
    end do

    ! Only nudge during the user-defined time window
    if (time >= C%BMB_inversion_t_start .and.  time <= C%BMB_inversion_t_end) then

      do vi = mesh%vi1, mesh%vi2

        ! Only nudge vertices that are shelf in the modelled or reference geometry

        ! (This means that, during the inversion (which happens in pseudo-time), it can happen
        ! that we apply basal melt to grounded ice, if the grounding line temporarily advances,
        ! or refreezing to shelf vertices when the grounding line temporarily retreats. However,
        ! this happens in pseudo-time, and will stop once the geometry converges to the target)

        if (BMB_inv%target_mask_shelf( vi) .or. ice%mask_floating_ice( vi)) then

          deltaH = ice%Hi( vi) - Hi_target_tot( vi)
          dHdt   = ice%dHi_dt( vi)

          dBMBdt = c_H * deltaH + c_dHdt * dHdt

          BMB_inv%BMB( vi) = BMB_inv%BMB( vi) + C%dt_BMB * dBMBdt

        end if
      end do

      ! Apply limits
      BMB_inv%BMB =  max( -C%BMB_maximum_allowed_melt_rate, min( C%BMB_maximum_allowed_refreezing_rate, BMB_inv%BMB ))

    end if

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

    if (par%primary) write(0,*) '   Initialising basal mass balance model ' // colour_string('inverted','light blue') // &
      ' with ' // colour_string(trim(C%choice_inversion_target_geometry),'light blue') // ' target geometry...'

    allocate( BMB_inv%BMB( mesh%vi1:mesh%vi2), source = 0._dp)

    call initialise_BMB_model_inverted_set_target_geometry( mesh, BMB_inv, refgeo_PD, refgeo_init)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BMB_model_inverted

  subroutine initialise_BMB_model_inverted_set_target_geometry( mesh, BMB_inv, refgeo_PD, refgeo_init)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_BMB_model_inverted), intent(inout) :: BMB_inv
    type(type_reference_geometry), intent(in   ) :: refgeo_PD, refgeo_init

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_BMB_model_inverted_set_target_geometry'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_reference_geometry_on_mesh( mesh, BMB_inv%target_geometry)

    select case (C%choice_inversion_target_geometry)
    case default
      call crash('unknown choice_inversion_target_geometry "' // trim(C%choice_inversion_target_geometry) // '"')
    case ('init')
      BMB_inv%target_geometry%Hi = refgeo_init%Hi
      BMB_inv%target_geometry%Hb = refgeo_init%Hb
      BMB_inv%target_geometry%Hs = refgeo_init%Hs
      BMB_inv%target_geometry%SL = refgeo_init%SL
    case ('PD')
      BMB_inv%target_geometry%Hi = refgeo_PD%Hi
      BMB_inv%target_geometry%Hb = refgeo_PD%Hb
      BMB_inv%target_geometry%Hs = refgeo_PD%Hs
      BMB_inv%target_geometry%SL = refgeo_PD%SL
    end select

    ! Determine the shelf mask of the target geometry
    allocate( BMB_inv%target_mask_shelf( mesh%vi1:mesh%vi2), source = .false.)
    do vi = mesh%vi1, mesh%vi2
      if (BMB_inv%target_geometry%Hi( vi) > 0.1_dp) then
        BMB_inv%target_mask_shelf( vi) = is_floating( &
          BMB_inv%target_geometry%Hi( vi), &
          BMB_inv%target_geometry%Hb( vi), &
          BMB_inv%target_geometry%SL( vi))
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BMB_model_inverted_set_target_geometry

end module BMB_inverted
