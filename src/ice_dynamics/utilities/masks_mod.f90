module masks_mod
  !< Calculating masks (e.g. mask_ice, mask_shelf, etc.)

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mpi_distributed_memory, only: gather_to_all
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use ice_geometry_basics, only: is_floating
  use projections, only: oblique_sg_projection

  implicit none

  private

  public :: determine_masks, calc_mask_noice, calc_mask_noice_remove_Ellesmere

contains

  subroutine determine_masks( mesh, ice)
    !< Determine the different masks

    ! mask_icefree_land       ! T: ice-free land , F: otherwise
    ! mask_icefree_ocean      ! T: ice-free ocean, F: otherwise
    ! mask_grounded_ice       ! T: grounded ice  , F: otherwise
    ! mask_floating_ice       ! T: floating ice  , F: otherwise
    ! mask_icefree_land_prev  ! T: ice-free land , F: otherwise (during previous time step)
    ! mask_icefree_ocean_prev ! T: ice-free ocean, F: otherwise (during previous time step)
    ! mask_grounded_ice_prev  ! T: grounded ice  , F: otherwise (during previous time step)
    ! mask_floating_ice_prev  ! T: floating ice  , F: otherwise (during previous time step)
    ! mask_margin             ! T: ice next to ice-free, F: otherwise
    ! mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
    ! mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
    ! mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    ! mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    ! mask_coastline          ! T: ice-free land next to ice-free ocean, F: otherwise

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_masks'
    integer                        :: vi, ci, vj
    logical, dimension(mesh%nV)    :: mask_icefree_land_tot
    logical, dimension(mesh%nV)    :: mask_icefree_ocean_tot
    logical, dimension(mesh%nV)    :: mask_grounded_ice_tot
    logical, dimension(mesh%nV)    :: mask_floating_ice_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! === Basic masks ===
    ! ===================

    ! Store previous basic masks
    ice%mask_icefree_land_prev  = ice%mask_icefree_land
    ice%mask_icefree_ocean_prev = ice%mask_icefree_ocean
    ice%mask_grounded_ice_prev  = ice%mask_grounded_ice
    ice%mask_floating_ice_prev  = ice%mask_floating_ice

    ! Initialise basic masks
    ice%mask_icefree_land  = .false.
    ice%mask_icefree_ocean = .false.
    ice%mask_grounded_ice  = .false.
    ice%mask_floating_ice  = .false.
    ice%mask               = 0

    ! Calculate
    do vi = mesh%vi1, mesh%vi2

      if (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) then
        ! Ice thickness is below the floatation thickness; either floating ice, or ice-free ocean

        if (ice%Hi( vi) > 0._dp) then
          ! Floating ice

          ice%mask_floating_ice( vi) = .true.
          ice%mask( vi) = C%type_floating_ice

        else
          ! Ice-free ocean

          ice%mask_icefree_ocean( vi) = .true.
          ice%mask( vi) = C%type_icefree_ocean

        end if

      else
        ! Ice thickness is above the floatation thickness; either grounded ice, or ice-free land

        if (ice%Hi( vi) > 0._dp) then
          ! Grounded ice

          ice%mask_grounded_ice( vi) = .true.
          ice%mask( vi) = C%type_grounded_ice

        else
          ! Ice-free land

          ice%mask_icefree_land( vi) = .true.
          ice%mask( vi) = C%type_icefree_land

        end if

      end if

    end do

    ! === Transitional masks ===
    ! ==========================

    ! Gather basic masks to all processes
    call gather_to_all( ice%mask_icefree_land , mask_icefree_land_tot )
    call gather_to_all( ice%mask_icefree_ocean, mask_icefree_ocean_tot)
    call gather_to_all( ice%mask_grounded_ice , mask_grounded_ice_tot )
    call gather_to_all( ice%mask_floating_ice , mask_floating_ice_tot )

    ! Initialise transitional masks
    ice%mask_margin    = .false.
    ice%mask_gl_gr     = .false.
    ice%mask_gl_fl     = .false.
    ice%mask_cf_gr     = .false.
    ice%mask_cf_fl     = .false.
    ice%mask_coastline = .false.

    do vi = mesh%vi1, mesh%vi2

      ! Ice margin
      if (mask_grounded_ice_tot( vi) .OR. mask_floating_ice_tot( vi)) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (.not. (mask_grounded_ice_tot( vj) .OR. mask_floating_ice_tot( vj))) then
            ice%mask_margin( vi) = .true.
            ice%mask( vi) = C%type_margin
          end if
        end do
      end if

      ! Grounding line (grounded side)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mask_floating_ice_tot( vj)) then
            ice%mask_gl_gr( vi) = .true.
            ice%mask( vi) = C%type_groundingline_gr
          end if
        end do
      end if

      ! Grounding line (floating side)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mask_grounded_ice_tot( vj)) then
            ice%mask_gl_fl( vi) = .true.
            ice%mask( vi) = C%type_groundingline_fl
          end if
        end do
      end if

      ! Calving front (grounded)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            ice%mask_cf_gr( vi) = .true.
            ice%mask( vi) = C%type_calvingfront_gr
          end if
        end do
      end if

      ! Calving front (floating)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            ice%mask_cf_fl( vi) = .true.
            ice%mask( vi) = C%type_calvingfront_fl
          end if
        end do
      end if

      ! Coastline
      if (mask_icefree_land_tot( vi)) then
        do ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            ice%mask_coastline( vi) = .true.
            ice%mask( vi) = C%type_coastline
          end if
        end do
      end if

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_masks

  subroutine calc_mask_noice( mesh, ice)
    !< Calculate the no-ice mask

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_noice'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    ! ==========

    ice%mask_noice = .false.

    ! domain-specific cases (mutually exclusive)
    ! ==========================================

    select case (C%choice_mask_noice)
    case default
      call crash('unknown choice_mask_noice "' // trim( C%choice_mask_noice) // '"')
      case ('none')
        ! Ice is (in principle) allowed everywhere

        ice%mask_noice = .false.

      case ('MISMIP_mod')
        ! Kill all ice when r > 900 km

        do vi = mesh%vi1, mesh%vi2
          if (NORM2( mesh%V( vi,:)) > 900E3_dp) then
            ice%mask_noice( vi) = .true.
          else
            ice%mask_noice( vi) = .false.
          end if
        end do

      case ('MISMIP+')
        ! Kill all ice when x > 640 km

        do vi = mesh%vi1, mesh%vi2
          if (mesh%V( vi,1) > 640E3_dp) then
            ice%mask_noice( vi) = .true.
          else
            ice%mask_noice( vi) = .false.
          end if
        end do

      case ('remove_Ellesmere')
        ! Prevent ice growth in the Ellesmere Island part of the Greenland domain

        call calc_mask_noice_remove_Ellesmere( mesh, ice%mask_noice)

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_noice

  subroutine calc_mask_noice_remove_Ellesmere( mesh, mask_noice)
    !< Prevent ice growth in the Ellesmere Island part of the Greenland domain

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(inout) :: mask_noice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_noice_remove_Ellesmere'
    integer                        :: vi
    real(dp), dimension(2)         :: pa_latlon, pb_latlon, pa, pb
    real(dp)                       :: xa, ya, xb, yb, yl_ab

    ! Add routine to path
    call init_routine( routine_name)

    ! The two endpoints in lat,lon
    pa_latlon = [76.74_dp, -74.79_dp]
    pb_latlon = [82.19_dp, -60.00_dp]

    ! The two endpoints in x,y
    call oblique_sg_projection( pa_latlon(2), pa_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%beta_stereo, xa, ya)
    call oblique_sg_projection( pb_latlon(2), pb_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%beta_stereo, xb, yb)

    pa = [xa,ya]
    pb = [xb,yb]

    do vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1)) * (pb(2)-pa(2)) / (pb(1)-pa(1))
      if (mesh%V( vi,2) > pa(2) .and. mesh%V( vi,2) > yl_ab .and. mesh%V( vi,1) < pb(1)) then
        mask_noice( vi) = .true.
      else
        mask_noice( vi) = .false.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_noice_remove_Ellesmere

end module masks_mod
