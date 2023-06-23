MODULE ice_model_utilities

  ! Generally useful functions used by the ice model.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync, ierr
  USE control_resources_and_error_messaging                  , ONLY: init_routine, finalise_routine, crash
  USE model_configuration                                    , ONLY: C
  USE parameters                                             , ONLY: ice_density, seawater_density
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_logical_1D
  USE math_utilities                                         , ONLY: is_floating
  USE mesh_remapping                                         , ONLY: Atlas, create_map_from_xy_grid_to_mesh
  USE petsc_basic                                            , ONLY: mat_petsc2CSR
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, deallocate_matrix_CSR_dist
  USE grid_basic                                             , ONLY: gather_gridded_data_to_master_dp_2D
  USE mesh_operators                                         , ONLY: ddx_a_a_2D, ddy_a_a_2D, map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, &
                                                                     ddx_b_a_2D, ddy_b_a_2D

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE determine_masks( mesh, ice)
    ! Determine the different masks
    !
    ! mask_icefree_land       ! T: ice-free land , F: otherwise
    ! mask_icefree_ocean      ! T: ice-free ocean, F: otherwise
    ! mask_grounded_ice       ! T: grounded ice  , F: otherwise
    ! mask_floating_ice       ! T: floating ice  , F: otherwise
    ! mask_icefree_land_prev  ! T: ice-free land , F: otherwise (during previous time step)
    ! mask_icefree_ocean_prev ! T: ice-free ocean, F: otherwise (during previous time step)
    ! mask_grounded_ice_prev  ! T: grounded ice  , F: otherwise (during previous time step)
    ! mask_floating_ice_prev  ! T: floating ice  , F: otherwise (during previous time step)
    ! mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
    ! mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
    ! mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    ! mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_masks'
    INTEGER                                            :: vi, ci, vj
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_icefree_land_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_icefree_ocean_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_grounded_ice_tot
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_floating_ice_tot

    ! Add routine to path
    CALL init_routine( routine_name)

  ! === Basic masks ===
  ! ===================

    ! Store previous basic masks
    ice%mask_icefree_land_prev  = ice%mask_icefree_land
    ice%mask_icefree_ocean_prev = ice%mask_icefree_ocean
    ice%mask_grounded_ice_prev  = ice%mask_grounded_ice
    ice%mask_floating_ice_prev  = ice%mask_floating_ice

    ! Initialise basic masks
    ice%mask_icefree_land  = .FALSE.
    ice%mask_icefree_ocean = .FALSE.
    ice%mask_grounded_ice  = .FALSE.
    ice%mask_floating_ice  = .FALSE.
    ice%mask               = 0

    ! Calculate
    DO vi = mesh%vi1, mesh%vi2

      IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN
        ! Ice thickness is below the floatation thickness; either floating ice, or ice-free ocean

        IF (ice%Hi( vi) > 0._dp) THEN
          ! Floating ice

          ice%mask_floating_ice( vi) = .TRUE.
          ice%mask( vi) = C%type_floating_ice

        ELSE ! IF (ice%Hi( vi) > 0._dp) THEN
          ! Ice-free ocean

          ice%mask_icefree_ocean( vi) = .TRUE.
          ice%mask( vi) = C%type_icefree_ocean

        END IF ! IF (ice%Hi( vi) > 0._dp) THEN

      ELSE ! IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN
        ! Ice thickness is above the floatation thickness; either grounded ice, or ice-free land

        IF (ice%Hi( vi) > 0._dp) THEN
          ! Grounded ice

          ice%mask_grounded_ice( vi) = .TRUE.
          ice%mask( vi) = C%type_grounded_ice

        ELSE ! IF (ice%Hi( vi) > 0._dp) THEN
          ! Ice-free land

          ice%mask_icefree_land( vi) = .TRUE.
          ice%mask( vi) = C%type_icefree_land

        END IF ! IF (ice%Hi( vi) > 0._dp) THEN

      END IF ! IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! === Transitional masks ===
    ! ==========================

    ! Gather basic masks to all processes
    CALL gather_to_all_logical_1D( ice%mask_icefree_land, mask_icefree_land_tot)

    ! Initialise transitional masks
    ice%mask_gl_gr = .FALSE.
    ice%mask_gl_fl = .FALSE.
    ice%mask_cf_gr = .FALSE.
    ice%mask_cf_fl = .FALSE.

    DO vi = mesh%vi1, mesh%vi2

      ! Grounding line (grounded side)
      IF (mask_grounded_ice_tot( vi)) THEN
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          IF (mask_floating_ice_tot( vj)) THEN
            ice%mask_gl_gr( vi) = .TRUE.
            ice%mask( vi) = C%type_groundingline_gr
          END IF
        END DO
      END IF

      ! Grounding line (floating side)
      IF (mask_floating_ice_tot( vi)) THEN
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          IF (mask_grounded_ice_tot( vj)) THEN
            ice%mask_gl_fl( vi) = .TRUE.
            ice%mask( vi) = C%type_groundingline_fl
          END IF
        END DO
      END IF

      ! Calving front (grounded)
      IF (mask_grounded_ice_tot( vi)) THEN
        DO ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          IF (mask_icefree_ocean_tot( vj)) THEN
            ice%mask_cf_gr( vi) = .TRUE.
            ice%mask( vi) = C%type_calvingfront_gr
          END IF
        END DO
      END IF

      ! Calving front (floating)
      IF (mask_floating_ice_tot( vi)) THEN
        DO ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          IF (mask_icefree_ocean_tot( vj)) THEN
            ice%mask_cf_fl( vi) = .TRUE.
            ice%mask( vi) = C%type_calvingfront_fl
          END IF
        END DO
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_masks

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
    integer                                      :: vi, k, n, i, j
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

    do vi = mesh%vi1, mesh%vi2

      ! Clear the list
      hb_list = 0d0

      ! Skip vertices at edge of domain
      if (mesh%VBI( vi) > 0) cycle

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      do k = M_map%ptr( vi), M_map%ptr( vi+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      end do

      ! Safety
      if (n_grid_cells == 0) then
        ! Use default mesh value
        ice%bedrock_cdf( vi,:) = ice%Hb( vi)
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

      ice%bedrock_cdf( vi, 1) = hb_list( 1)
      ice%bedrock_cdf( vi,11) = hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins,
      ! from the second (10%) to the tenth (90%)
      do i = 2, 10
        isc  = 1._dp + (real( n_grid_cells,dp) - 1._dp) * (real( i,dp) - 1._dp) / 10._dp
        ii0  = floor( isc)
        ii1  = ceiling( isc)
        wii0 = real( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf( vi,i) = wii0 * hb_list( ii0) + wii1 * hb_list( ii1)
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

  subroutine calc_grounded_fractions( mesh, ice)
    ! Determine the sub-grid grounded-area fractions of all grid cells from a bedrock CDF

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_grounded_fractions'
    integer                             :: vi, il, iu, ti
    real(dp)                            :: hb_float, wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    ! === On the a-grid ===
    ! =====================

    do vi = mesh%vi1, mesh%vi2

      ! Edge vertices
      if (mesh%VBI( vi) > 0) then
        ! Either 0 or 1 depending on land mask
        if (ice%mask_icefree_land( vi) .or. ice%mask_grounded_ice( vi)) then
          ice%fraction_gr( vi) = 1._dp
        else
          ice%fraction_gr( vi) = 0._dp
        end if
        ! Then skip
        cycle
      end if

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      hb_float = ice%SL( vi) - ice%Hi( vi) * ice_density/seawater_density - ice%dHb( vi)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! hb_float as a linear interpolation of the numbers in the CDF.

      if     (hb_float <= minval( ice%bedrock_cdf( vi,:))) then
        ! All sub-grid points are above the floating bedrock elevation
        ice%fraction_gr( vi) = 1._dp
      elseif (hb_float >= maxval( ice%bedrock_cdf( vi,:))) then
        ! All sub-grid points are below the floating bedrock elevation
        ice%fraction_gr( vi) = 0._dp
      else
        ! Find the 2 elements in the CDF surrounding hb_float
        iu = 1
        do while (ice%bedrock_cdf( vi,iu) < hb_float)
          iu = iu+1
        end do
        il = iu-1

        ! Interpolate bedrock to get the weights for both CDF percentages
        wl = (ice%bedrock_cdf( vi,iu) - hb_float) / (ice%bedrock_cdf( vi,iu)-ice%bedrock_cdf( vi,il))
        wu = 1._dp - wl

        ! Interpolate between percentages, assuming that there are 11 CDF bins between
        ! 0% and 100%, at 10% intervals, i.e. element 1 is 0% =(1-1)*10%, element 2 is
        ! 10% = (2-1)*10%, and so on.
        ice%fraction_gr( vi) = 1._dp - (10._dp*(il-1)*wl + 10._dp*(iu-1)*wu) / 100._dp

        ! Safety
        ice%fraction_gr( vi) = min( 1._dp, max( 0._dp, ice%fraction_gr( vi)))

      end if

    end do

    ! === On the b-grid ===
    ! =====================

    ! Map from the a-grid to the b-grid
    call map_a_b_2D( mesh, ice%fraction_gr, ice%fraction_gr_b)

    ! Safety
    do ti = mesh%ti1, mesh%ti2
      ice%fraction_gr_b( ti) = min( 1._dp, max( 0._dp, ice%fraction_gr_b( ti)))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions

  SUBROUTINE calc_zeta_gradients( mesh, ice)
    ! Calculate all the gradients of zeta, needed to perform the scaled vertical coordinate transformation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_zeta_gradients'
!    REAL(dp), DIMENSION(:    ), POINTER                :: Hi_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dx2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dxdy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dy2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: Hi_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dx2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dxdy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dy2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: Hs_b
!    REAL(dp), DIMENSION(:    ), POINTER                :: Hs_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dx2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dxdy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dy2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dx2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dxdy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dy2_b
    INTEGER                                            :: vi,ti,k,ks
    REAL(dp)                                           :: Hi, zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
   !ALLOCATE( Hi_a(        mesh%vi1:mesh%vi2))
    ALLOCATE( dHi_dx_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( dHi_dy_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hi_dx2_a(  mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hi_dxdy_a( mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hi_dy2_a(  mesh%vi1:mesh%vi2))

    ALLOCATE( Hi_b(        mesh%ti1:mesh%ti2))
    ALLOCATE( dHi_dx_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( dHi_dy_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hi_dx2_b(  mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hi_dxdy_b( mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hi_dy2_b(  mesh%ti1:mesh%ti2))

   !ALLOCATE( Hs_a(        mesh%vi1:mesh%vi2))
    ALLOCATE( dHs_dx_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( dHs_dy_a(    mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hs_dx2_a(  mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hs_dxdy_a( mesh%vi1:mesh%vi2))
    ALLOCATE( d2Hs_dy2_a(  mesh%vi1:mesh%vi2))

    ALLOCATE( Hs_b(        mesh%ti1:mesh%ti2))
    ALLOCATE( dHs_dx_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( dHs_dy_b(    mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hs_dx2_b(  mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hs_dxdy_b( mesh%ti1:mesh%ti2))
    ALLOCATE( d2Hs_dy2_b(  mesh%ti1:mesh%ti2))

  ! Calculate gradients of Hi and Hs on both grids

    !CALL map_a_a_2D( mesh, ice%Hi_a, Hi_a       )
    CALL ddx_a_a_2D( mesh, ice%Hi  , dHi_dx_a   )
    CALL ddy_a_a_2D( mesh, ice%Hi  , dHi_dy_a   )
    CALL map_a_b_2D( mesh, ice%Hi  , Hi_b       )
    CALL ddx_a_b_2D( mesh, ice%Hi  , dHi_dx_b   )
    CALL ddy_a_b_2D( mesh, ice%Hi  , dHi_dy_b   )
    CALL ddx_b_a_2D( mesh, dHi_dx_b, d2Hi_dx2_a )
    CALL ddy_b_a_2D( mesh, dHi_dx_b, d2Hi_dxdy_a)
    CALL ddy_b_a_2D( mesh, dHi_dy_b, d2Hi_dy2_a )
    CALL ddx_a_b_2D( mesh, dHi_dx_a, d2Hi_dx2_b )
    CALL ddy_a_b_2D( mesh, dHi_dx_a, d2Hi_dxdy_b)
    CALL ddy_a_b_2D( mesh, dHi_dy_a, d2Hi_dy2_b )

    !CALL map_a_a_2D( mesh, ice%Hs_a, Hs_a       )
    CALL ddx_a_a_2D( mesh, ice%Hs  , dHs_dx_a   )
    CALL ddy_a_a_2D( mesh, ice%Hs  , dHs_dy_a   )
    CALL map_a_b_2D( mesh, ice%Hs  , Hs_b       )
    CALL ddx_a_b_2D( mesh, ice%Hs  , dHs_dx_b   )
    CALL ddy_a_b_2D( mesh, ice%Hs  , dHs_dy_b   )
    CALL ddx_b_a_2D( mesh, dHs_dx_b, d2Hs_dx2_a )
    CALL ddy_b_a_2D( mesh, dHs_dx_b, d2Hs_dxdy_a)
    CALL ddy_b_a_2D( mesh, dHs_dy_b, d2Hs_dy2_a )
    CALL ddx_a_b_2D( mesh, dHs_dx_a, d2Hs_dx2_b )
    CALL ddy_a_b_2D( mesh, dHs_dx_a, d2Hs_dxdy_b)
    CALL ddy_a_b_2D( mesh, dHs_dy_a, d2Hs_dy2_b )

    ! Calculate zeta gradients on all grids

    ! ak
    DO vi = mesh%vi1, mesh%vi2

      Hi = MAX( 10._dp, ice%Hi( vi))

      DO k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dt_ak(    vi,k) = ( 1._dp / Hi) * (ice%dHs_dt( vi) - zeta * ice%dHi_dt( vi))

        ice%dzeta_dx_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dx_a( vi) - zeta * dHi_dx_a( vi))
        ice%dzeta_dy_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dy_a( vi) - zeta * dHi_dy_a( vi))
        ice%dzeta_dz_ak(    vi,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_ak(  vi,k) = (dHi_dx_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dx2_a(  vi) - zeta * d2Hi_dx2_a(  vi))
        ice%d2zeta_dxdy_ak( vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dxdy_a( vi) - zeta * d2Hi_dxdy_a( vi))
        ice%d2zeta_dy2_ak(  vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dy_ak( vi,k) + (1._dp / Hi) * (d2Hs_dy2_a(  vi) - zeta * d2Hi_dy2_a(  vi))

      END DO
    END DO

    ! bk
    DO ti = mesh%ti1, mesh%ti2

      Hi = MAX( 10._dp, Hi_b( ti))

      DO k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dx_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bk(    ti,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_bk(  ti,k) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bk( ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bk(  ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bk( ti,k) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      END DO
    END DO

    ! bks
    DO ti = mesh%ti1, mesh%ti2

      Hi = MAX( 10._dp, Hi_b( ti))

      DO ks = 1, mesh%nz-1

        zeta = mesh%zeta_stag( ks)

        ice%dzeta_dx_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bks(    ti,ks) = (-1._dp / Hi)

        ice%d2zeta_dx2_bks(  ti,ks) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bks( ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bks(  ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_zeta_gradients

END MODULE ice_model_utilities
