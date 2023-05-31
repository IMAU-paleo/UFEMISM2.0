MODULE ice_model_utilities

  ! Generally useful functions used by the ice model.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync, ierr
  USE control_resources_and_error_messaging                  , ONLY: init_routine, finalise_routine, crash
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D, distribute_from_master_int_1D
  USE math_utilities                                         , ONLY: is_floating
  USE mesh_remapping                                         , ONLY: Atlas, create_map_from_xy_grid_to_mesh
  USE petsc_basic                                            , ONLY: mat_petsc2CSR
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, deallocate_matrix_CSR_dist
  USE grid_basic                                             , ONLY: gather_gridded_data_to_master_dp_2D

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

    ! Fill in the full arrays so all processes have them
    call gather_to_all_dp_1D( ice%Hi, Hi)
    call gather_to_all_dp_1D( ice%Hb, Hb)
    call gather_to_all_dp_1D( ice%SL, SL)

    ! Allocate local full-array diagnostic mask for each process
    allocate( mask ( mesh%nV ))

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

    ! Distribute local mask to each process
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

  subroutine calc_bedrock_CDFs( mesh, refgeo, ice)
    ! Calculate the bedrock cumulative density functions

    implicit none

    ! In/output variables:
    type(type_mesh),               intent(in)    :: mesh
    type(type_reference_geometry), intent(in)    :: refgeo
    type(type_ice_model),          intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'calc_bedrock_CDFs'
    type(type_sparse_matrix_CSR_dp)              :: M_map
    logical                                      :: found_map, found_empty_page
    integer                                      :: mi, mi_valid
    real(dp), dimension(:,:  ), allocatable      :: Hb_grid_tot
    real(dp), dimension(:    ), allocatable      :: hb_list
    integer                                      :: vi_glob, vi_loc, k, n, i, j
    integer                                      :: n_grid_cells, ii0, ii1
    real(dp)                                     :: isc, wii0, wii1

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) write(0,*) '  Initialising sub-grid grounded-area fractions...'

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == refgeo%grid_raw%name .AND. Atlas( mi)%name_dst == mesh%name) THEN
        found_map = .TRUE.
        mi_valid  = mi
        EXIT
      END IF
    END DO

    ! If no appropriate mapping object could be found, create one.
    IF (.NOT. found_map) THEN
      found_empty_page = .FALSE.
      DO mi = 1, SIZE( Atlas,1)
        IF (.NOT. Atlas( mi)%is_in_use) THEN
          found_empty_page = .TRUE.
          CALL create_map_from_xy_grid_to_mesh( refgeo%grid_raw, mesh, Atlas( mi))
          mi_valid = mi
          EXIT
        END IF
      END DO
      ! Safety
      IF (.NOT. found_empty_page) CALL crash('No more room in Atlas - assign more memory!')
    END IF

    ! Turn map into CSR, so we can use it here
    CALL mat_petsc2CSR( Atlas( mi_valid)%M, M_map)

    ! Get complete gridded bedrock elevation on all processes
    ALLOCATE( Hb_grid_tot( refgeo%grid_raw%nx, refgeo%grid_raw%ny))
    CALL gather_gridded_data_to_master_dp_2D( refgeo%grid_raw, refgeo%Hb_grid_raw, Hb_grid_tot)
    CALL MPI_BCAST( Hb_grid_tot, refgeo%grid_raw%nx * refgeo%grid_raw%ny, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory for list of bedrock elevations
    allocate( hb_list( refgeo%grid_raw%nx * refgeo%grid_raw%ny ))
    hb_list = 0d0

    ! Initialise cumulative density function (CDF)
    ice%bedrock_cdf = 0._dp

    ! === Scan ===
    ! ============

    do vi_loc = 1, mesh%nV_loc

      ! Compute global index for this vertex
      vi_glob = vi_loc + mesh%vi1 - 1

      ! Clear the list
      hb_list = 0d0

      ! Skip vertices at edge of domain
      if (mesh%VBI( vi_glob) > 0) cycle

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      do k = M_map%ptr( vi_glob), M_map%ptr( vi_glob+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      end do

      ! Safety
      if (n_grid_cells == 0) then
        ! Use default mesh value
        ice%bedrock_cdf( vi_loc,:) = ice%Hb( vi_loc)
        ! And skip
        cycle
      end if

      ! === Cumulative density function ===
      ! ===================================

      ! Inefficient but easy sorting of hb_list
      do i = 1, n_grid_cells-1
      do j = i+1, n_grid_cells
        if (hb_list( i) > hb_list( j)) then
          hb_list( i) = hb_list( i) + hb_list( j)
          hb_list( j) = hb_list( i) - hb_list( j)
          hb_list( i) = hb_list( i) - hb_list( j)
        end if
      end do
      end do

      ! Set first (0%) and last bins (100%) of the CDF to the minimum
      ! and maximum bedrock elevations scanned, respectively

      ! NOTE: should the number of bins be configurable?

      ice%bedrock_cdf( vi_loc, 1) = hb_list( 1)
      ice%bedrock_cdf( vi_loc,11) = hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins,
      ! from the second (10%) to the tenth (90%)
      do i = 2, 10
        isc  = 1._dp + (real( n_grid_cells,dp) - 1._dp) * (real( i,dp) - 1._dp) / 10._dp
        ii0  = floor( isc)
        ii1  = ceiling( isc)
        wii0 = real( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf( vi_loc,i) = wii0 * hb_list( ii0) + wii1 * hb_list( ii1)
      end do

    end do
    CALL sync

    ! === Finalisation ===
    ! ====================

    ! Clean up after yourself
    deallocate( hb_list)
    deallocate( Hb_grid_tot)
    call deallocate_matrix_CSR_dist( M_map)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bedrock_CDFs

END MODULE ice_model_utilities
