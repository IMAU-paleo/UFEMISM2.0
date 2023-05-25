MODULE ice_model_utilities

  ! Generally useful functions used by the ice model.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, distribute_from_master_int_1D
  USE math_utilities                                         , ONLY: is_floating

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

subroutine determine_masks( mesh, ice)
    ! Determine the different masks, on both the Aa and the Ac mesh

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'determine_masks'
    integer                             :: vi, ci, vc
    real(dp), dimension(:), allocatable :: Hi, Hb, SL
    integer, dimension(:), allocatable  :: mask

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate local full-array geometry for each process
    allocate( Hi   ( mesh%nV ))
    allocate( Hb   ( mesh%nV ))
    allocate( SL   ( mesh%nV ))
    allocate( mask ( mesh%nV ))

    ! Fill in the full arrays so all processes have them
    call gather_to_all_dp_1D( ice%Hi, Hi)
    call gather_to_all_dp_1D( ice%Hb, Hb)
    call gather_to_all_dp_1D( ice%SL, SL)

    mask = C%type_land

    ! === Basic masks ===
    ! ===================

    ! Land
    ! ====

    ! Start out with land everywhere, fill in the rest based on input.
    ice%mask_land   = .true.
    ice%mask_ocean  = .false.
    ice%mask_lake   = .false.
    ice%mask_ice    = .false.
    ice%mask_sheet  = .false.
    ice%mask_shelf  = .false.
    ice%mask_coast  = .false.
    ice%mask_margin = .false.
    ice%mask_gl_gr  = .false.
    ice%mask_gl_fl  = .false.
    ice%mask_cf_gr  = .false.
    ice%mask_cf_fl  = .false.

    do vi = 1, mesh%nV

      ! Ocean
      ! =====

      ! Both open and shelf-covered
      if (is_floating( Hi( vi), Hb( vi), SL( vi))) then
        ice%mask_ocean( vi) = .true.
        ice%mask_land(  vi) = .false.
        mask(           vi) = C%type_ocean
      end if

      ! Ice
      ! ===

      if (Hi( vi) > 0._dp) then
        ice%mask_ice( vi)  = .true.
      end if

      ! Ice sheet
      ! =========

      if (ice%mask_ice( vi) .and. ice%mask_land( vi)) then
        ice%mask_sheet( vi) = .true.
        mask(           vi) = C%type_sheet
      end if

      ! Ice shelf
      ! =========

      if (ice%mask_ice( vi) .and. ice%mask_ocean( vi)) then
        ice%mask_shelf( vi) = .true.
        mask(           vi) = C%type_shelf
      end if

    end do

    ! === Transitional masks ===
    ! ==========================

    do vi = 1, mesh%nV

      ! Coastline
      ! =========

      if (ice%mask_land( vi) .and. (.not. ice%mask_ice( vi))) then
        ! Ice-free land bordering ocean equals coastline
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ocean( vc)) then
            ice%mask_coast( vi) = .true.
            mask( vi) = C%type_coast
          end if
        end do
      end if

      ! Ice margin
      ! ==========

      if (ice%mask_ice( vi)) then
        ! Ice bordering non-ice equals margin
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (.not. ice%mask_ice( vc)) then
            ice%mask_margin( vi) =  .true.
            mask( vi) = C%type_margin
          end if
        end do
      end if

      ! Grounding line (ice sheet side)
      ! ===============================

      if (ice%mask_sheet( vi)) then
        ! Sheet bordering shelf equals grounding line
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_shelf( vc)) then
            ice%mask_gl_gr( vi) = .true.
            mask( vi) = C%type_groundingline_gr
          end if
        end do
      end if

      ! Grounding line (ice shelf side)
      ! ===============================

      if (ice%mask_shelf( vi)) then
        ! Shelf bordering sheet equals floating side of grounding line
        do ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_sheet( vc)) then
            ice%mask_gl_fl( vi) =  .true.
            mask( vi) = C%type_groundingline_fl
          end if
        end do
      end if

      ! Calving front (grounded)
      ! ========================

      if (ice%mask_sheet( vi)) then
        ! Ice sheet bordering open ocean equals calving front
        do ci = 1, mesh%nC(vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ocean( vc) .and. (.not. ice%mask_ice( vc))) then
            ice%mask_cf_gr( vi) = .true.
            mask( vi) = C%type_calvingfront_gr
          end if
        end do
      end if

      ! Calving front
      ! =============

      if (ice%mask_shelf( vi)) then
        ! Ice shelf bordering open ocean equals calving front
        do ci = 1, mesh%nC(vi)
          vc = mesh%C( vi,ci)
          if (ice%mask_ocean( vc) .and. (.not. ice%mask_ice( vc))) then
            ice%mask_cf_fl( vi) = .true.
            mask( vi) = C%type_calvingfront_fl
          end if
        end do
      end if

    end do ! vi = 1, mesh%nV

    ! === Diagnostic mask ===
    ! =======================

    call distribute_from_master_int_1D( mask, ice%mask)

    ! === Finalisation ===
    ! ====================

    deallocate( Hi  )
    deallocate( Hb  )
    deallocate( SL  )
    deallocate( mask)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_masks

END MODULE ice_model_utilities
