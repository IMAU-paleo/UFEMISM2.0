module subgrid_ice_margin

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mpi_distributed_memory, only: gather_to_all
  use ice_geometry_basics, only: is_floating

  implicit none

  private

  public :: calc_effective_thickness

contains

  subroutine calc_effective_thickness( mesh, ice, Hi, Hi_eff, fraction_margin)
    !< Determine the ice-filled fraction and effective ice thickness of floating margin pixels

    ! In- and output variables
    type(type_mesh),      intent(in   )                 :: mesh
    type(type_ice_model), intent(in   )                 :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)  :: Hi
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out) :: Hi_eff
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out) :: fraction_margin

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_effective_thickness'
    integer                                :: vi, ci, vc
    real(dp)                               :: Hi_neighbour_max
    real(dp), dimension(mesh%nV)           :: Hi_tot, Hb_tot, SL_tot
    logical,  dimension(mesh%vi1:mesh%vi2) :: mask_margin, mask_floating
    logical,  dimension(mesh%nV)           :: mask_margin_tot, mask_floating_tot

    ! == Initialisation
    ! =================

    ! Add routine to path
    call init_routine( routine_name)

    ! Collect Hi from all processes
    call gather_to_all( Hi,     Hi_tot)
    call gather_to_all( ice%Hb, Hb_tot)
    call gather_to_all( ice%SL, SL_tot)

    ! == Margin mask
    ! ==============

    ! Initialise
    mask_margin = .false.

    do vi = mesh%vi1, mesh%vi2
      do ci = 1, mesh%nC( vi)
        vc = mesh%C( vi,ci)
        if (Hi_tot( vi) > 0._dp .and. Hi_tot( vc) == 0._dp) then
          mask_margin( vi) = .true.
        end if
      end do
    end do

    ! Gather mask values from all processes
    call gather_to_all( mask_margin, mask_margin_tot)

    ! == Floating mask
    ! ================

    ! Initialise
    mask_floating = .false.

    do vi = mesh%vi1, mesh%vi2
      if (is_floating( Hi_tot( vi), ice%Hb( vi), ice%SL( vi))) then
        mask_floating( vi) = .true.
      end if
    end do

    ! Gather mask values from all processes
    call gather_to_all( mask_floating, mask_floating_tot)

    ! == default values
    ! =================

    ! Initialise values
    do vi = mesh%vi1, mesh%vi2
      ! Grounded regions: ice-free land and grounded ice
      ! NOTE: important to let ice-free land have a non-zero
      ! margin fraction to let it get SMB in the ice thickness equation
      if (.not. mask_floating( vi)) then
        fraction_margin( vi) = 1._dp
        Hi_eff( vi) = Hi_tot( vi)
      ! Old and new floating regions
      elseif (Hi_tot( vi) > 0._dp) then
        fraction_margin( vi) = 1._dp
        Hi_eff( vi) = Hi_tot( vi)
      ! New ice-free ocean regions
      else
        fraction_margin( vi) = 0._dp
        Hi_eff( vi) = 0._dp
      end if
    end do

    ! === Compute ===
    ! ===============

    do vi = mesh%vi1, mesh%vi2

      ! Only check margin vertices
      if (.not. mask_margin_tot( vi)) then
        ! Simply use initialised values
        cycle
      end if

      ! === Max neighbour thickness ===
      ! ===============================

      ! Find the max ice thickness among non-margin neighbours
      Hi_neighbour_max = 0._dp

      do ci = 1, mesh%nC( vi)
        vc = mesh%C( vi,ci)

        ! Ignore margin neighbours
        if (mask_margin_tot( vc)) then
          cycle
        end if

        ! Floating margins check for neighbours
        if (mask_floating( vi)) then
          Hi_neighbour_max = max( Hi_neighbour_max, Hi_tot( vc))
        end if

      end do

      ! === Effective ice thickness ===
      ! ===============================

      ! Only apply if the thickest non-margin neighbour is thicker than
      ! this vertex. Otherwise, simply use initialised values.
      if (Hi_neighbour_max > Hi_tot( vi)) then
        ! Calculate sub-grid ice-filled fraction
        Hi_eff( vi) = Hi_neighbour_max
        fraction_margin( vi) = Hi_tot( vi) / Hi_neighbour_max
      end if

    end do

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_thickness

end module subgrid_ice_margin
