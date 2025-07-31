module ice_thickness_safeties
  !< Different kinds of "safeties" to keep the ice sheet stable during nudging-based initialisation runs

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use parameters, only: ice_density, seawater_density
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use subgrid_ice_margin, only: calc_effective_thickness
  use ice_geometry_basics, only: is_floating

  implicit none

  private

  public :: alter_ice_thickness

contains

  subroutine alter_ice_thickness( mesh, ice, Hi_old, Hi_new, refgeo, time)
    !< Modify the predicted ice thickness in some sneaky way

    ! In- and output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi_old
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_new
    type(type_reference_geometry),          intent(in   ) :: refgeo
    real(dp),                               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'alter_ice_thickness'
    integer                                :: vi
    real(dp)                               :: decay_start, decay_end
    real(dp)                               :: fixiness, limitness, fix_H_applied, limit_H_applied
    real(dp), dimension(mesh%vi1:mesh%vi2) :: modiness_up, modiness_down
    real(dp), dimension(mesh%vi1:mesh%vi2) :: Hi_save, Hi_eff_new, fraction_margin_new
    real(dp)                               :: floating_area, calving_area, mass_lost

    ! Add routine to path
    call init_routine( routine_name)

    ! Save predicted ice thickness for future reference
    Hi_save = Hi_new

    ! Calculate would-be effective thickness
    call calc_effective_thickness( mesh, ice, Hi_new, Hi_eff_new, fraction_margin_new)

    ! == Mask conservation
    ! ====================

    ! if desired, don't let grounded ice cross the floatation threshold
    if (C%do_protect_grounded_mask .and. time <= C%protect_grounded_mask_t_end) then
      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_grounded_ice( vi)) then
          Hi_new( vi) = max( Hi_new( vi), (ice%SL( vi) - ice%Hb( vi)) * seawater_density/ice_density + .1_dp)
        end if
      end do
    end if

    ! General cases
    ! =============

    ! if so specified, remove very thin ice
    do vi = mesh%vi1, mesh%vi2
      if (Hi_eff_new( vi) < C%Hi_min) then
        Hi_new( vi) = 0._dp
      end if
    end do

    ! if so specified, remove thin floating ice
    if (C%choice_calving_law == 'threshold_thickness') then
      do vi = mesh%vi1, mesh%vi2
        if (is_floating( Hi_eff_new( vi), ice%Hb( vi), ice%SL( vi)) .and. Hi_eff_new( vi) < C%calving_threshold_thickness_shelf) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! DENK DROM
    if (C%remove_ice_absent_at_PD) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo%Hi( vi) == 0._dp) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! if so specified, remove all floating ice
    if (C%do_remove_shelves) then
      do vi = mesh%vi1, mesh%vi2
        if (is_floating( Hi_eff_new( vi), ice%Hb( vi), ice%SL( vi))) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! if so specified, remove all floating ice beyond the present-day calving front
    if (C%remove_shelves_larger_than_PD) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo%Hi( vi) == 0._dp .and. refgeo%Hb( vi) < 0._dp) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! if so specified, remove all floating ice crossing the continental shelf edge
    if (C%continental_shelf_calving) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo%Hi( vi) == 0._dp .and. refgeo%Hb( vi) < C%continental_shelf_min_height) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! === Fixiness ===
    ! ================

    ! Intial value
    fixiness = 1._dp

    ! Make sure that the start and end times make sense
    decay_start = C%fixiness_t_start
    decay_end   = C%fixiness_t_end

    ! Compute decaying fixiness
    if (decay_start >= decay_end) then
      ! Fix interval makes no sense: ignore fixiness
      fixiness = 0._dp
    elseif (time <= decay_start) then
      ! Before fix interval: check chosen option
      if (C%do_fixiness_before_start) then
        fixiness = 1._dp
      else
        fixiness = 0._dp
      end if
    elseif (time >= decay_end) then
      ! After fix interval: remove any fix/delay
      fixiness = 0._dp
    else
      ! We're within the fix interval: fixiness decreases with time
      fixiness = 1._dp - (time - decay_start) / (decay_end - decay_start)
    end if

    ! Just in case
    fixiness = min( 1._dp, max( 0._dp, fixiness))

    ! === Limitness ===
    ! =================

    ! Intial value
    limitness = 1._dp

    ! Make sure that the start and end times make sense
    decay_start = C%limitness_t_start
    decay_end   = C%limitness_t_end

    ! Compute decaying limitness
    if (decay_start >= decay_end) then
      ! Limit interval makes no sense: ignore limitness
      limitness = 0._dp
    elseif (time <= decay_start) then
      ! Before limit interval: check chosen option
      if (C%do_limitness_before_start) then
        limitness = 1._dp
      else
        limitness = 0._dp
      end if
    elseif (time >= decay_end) then
      ! After limit interval: remove any limits
      limitness = 0._dp
    else
      ! Limitness decreases with time
      limitness = 1._dp - (time - decay_start) / (decay_end - decay_start)
    end if

    ! Just in case
    limitness = min( 1._dp, max( 0._dp, limitness))

    ! === Modifier ===
    ! ================

    ! Intial value
    modiness_up   = 0._dp
    modiness_down = 0._dp

    select case (C%modiness_H_style)
    case default
      call crash('unknown modiness_H_choice "' // trim( C%modiness_H_style) // '"')
    case ('none')
      modiness_up   = 0._dp
      modiness_down = 0._dp

    case ('Ti_hom')
      modiness_up   = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)
      modiness_down = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)

    case ('Ti_hom_up')
      modiness_up   = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)
      modiness_down = 0._dp

    case ('Ti_hom_down')
      modiness_up   = 0._dp
      modiness_down = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)

    case ('no_thick_inland')
      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_grounded_ice( vi) .and. .not. ice%mask_gl_gr( vi)) then
          modiness_up( vi) = 1._dp
        end if
      end do
      modiness_down = 0._dp

    case ('no_thin_inland')
      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_grounded_ice( vi) .and. .not. ice%mask_gl_gr( vi)) then
          modiness_down( vi) = 1._dp
        end if
      end do
      modiness_up = 0._dp

    end select

    ! Just in case
    modiness_up   = min( 1._dp, max( 0._dp, modiness_up  ))
    modiness_down = min( 1._dp, max( 0._dp, modiness_down))

    ! === Fix, delay, limit ===
    ! =========================

    do vi = mesh%vi1, mesh%vi2

      ! Initialise
      fix_H_applied   = 0._dp
      limit_H_applied = 0._dp

      if (    ice%mask_gl_gr( vi)) then
        fix_H_applied   = C%fixiness_H_gl_gr * fixiness
        limit_H_applied = C%limitness_H_gl_gr * limitness

      elseif (ice%mask_gl_fl( vi)) then
        fix_H_applied   = C%fixiness_H_gl_fl * fixiness
        limit_H_applied = C%limitness_H_gl_fl * limitness

      elseif (ice%mask_grounded_ice( vi)) then
        fix_H_applied   = C%fixiness_H_grounded * fixiness
        limit_H_applied = C%limitness_H_grounded * limitness

      elseif (ice%mask_floating_ice( vi)) then
        fix_H_applied   = C%fixiness_H_floating * fixiness
        limit_H_applied = C%limitness_H_floating * limitness

      elseif (ice%mask_icefree_land( vi)) then
        if (C%fixiness_H_freeland .and. fixiness > 0._dp) fix_H_applied = 1._dp
        limit_H_applied = C%limitness_H_grounded * limitness

      elseif (ice%mask_icefree_ocean( vi)) then
        if (C%fixiness_H_freeocean .and. fixiness > 0._dp) fix_H_applied = 1._dp
        limit_H_applied = C%limitness_H_floating * limitness
      else
        ! if we reached this point, vertex is neither grounded, floating, nor ice free. That's a problem
        call crash('vertex neither grounded, floating, nor ice-free?')
      end if

      ! Apply fixiness
      Hi_new( vi) = Hi_old( vi) * fix_H_applied + Hi_new( vi) * (1._dp - fix_H_applied)

      ! Apply limitness
      Hi_new( vi) = min( Hi_new( vi), refgeo%Hi( vi) + (1._dp - modiness_up(   vi)) * limit_H_applied &
                                                    + (1._dp - limitness         ) * (Hi_new( vi) - refgeo%Hi( vi)) )

      Hi_new( vi) = max( Hi_new( vi), refgeo%Hi( vi) - (1._dp - modiness_down( vi)) * limit_H_applied &
                                                    - (1._dp - limitness         ) * (refgeo%Hi( vi) - Hi_new( vi)) )

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine alter_ice_thickness

end module ice_thickness_safeties
